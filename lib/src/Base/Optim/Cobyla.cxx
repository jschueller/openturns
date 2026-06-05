//                                               -*- C++ -*-
/**
 *  @brief Cobyla is an actual implementation for OptimizationAlgorithmImplementation using the cobyla library
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "openturns/Cobyla.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/Log.hxx"
#include "prima/prima.hpp"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(Cobyla)

static const Factory<Cobyla> Factory_Cobyla;

/* Default constructor */
Cobyla::Cobyla()
  : OptimizationAlgorithmImplementation()
  , rhoBeg_(ResourceMap::GetAsScalar("Cobyla-DefaultRhoBeg"))
{
  // Nothing to do
}

Cobyla::Cobyla(const OptimizationProblem & problem)
  : OptimizationAlgorithmImplementation(problem)
  , rhoBeg_(ResourceMap::GetAsScalar("Cobyla-DefaultRhoBeg"))
{
  checkProblem(problem);
}

Cobyla::Cobyla(const OptimizationProblem & problem,
               const Scalar rhoBeg)
  : OptimizationAlgorithmImplementation(problem)
  , rhoBeg_(rhoBeg)
{
  checkProblem(problem);
}


/* Virtual constructor */
Cobyla * Cobyla::clone() const
{
  return new Cobyla(*this);
}

/** Check whether this problem can be solved by this solver.  Must be overloaded by the actual optimisation algorithm */
void Cobyla::checkProblem(const OptimizationProblem & problem) const
{
  if (problem.hasMultipleObjective())
    throw InvalidArgumentException(HERE) << "Error: " << this->getClassName() << " does not support multi-objective optimization";
  if (problem.hasResidualFunction())
    throw InvalidArgumentException(HERE) << getClassName() << " does not support least-square problems";
  if (!problem.isContinuous())
    throw InvalidArgumentException(HERE) << "Error: " << getClassName() << " does not support non continuous problems";

}

/* Performs the actual computation by calling the Cobyla algorithm
 */
void Cobyla::run()
{
  result_ = OptimizationResult(getProblem());

  const UnsignedInteger dimension = getProblem().getDimension();

  const Point startingPoint(getStartingPoint());
  if (startingPoint.getDimension() != dimension)
    throw InvalidArgumentException(HERE) << "Invalid starting point dimension (" << startingPoint.getDimension() << "), expected " << dimension;
  if (dimension == 0)
    throw InvalidArgumentException(HERE) << "Dimension of the problem is zero";

  if (getProblem().hasBounds())
  {
    const Interval bounds(getProblem().getBounds());
    if (!bounds.contains(startingPoint))
    {
      LOGWARN(OSS() << "Starting point is not inside bounds startingPoint=" << startingPoint.__str__() << " bounds=" << bounds);
    }
  }

  const Scalar rhoEnd = getMaximumAbsoluteError();
  const int maxFun = getMaximumCallsNumber();

  // initialize history
  evaluationInputHistory_ = Sample(0, dimension);
  evaluationOutputHistory_ = Sample(0, 1);
  equalityConstraintHistory_ = Sample(0, getProblem().getEqualityConstraint().getOutputDimension());
  inequalityConstraintHistory_ = Sample(0, getProblem().getInequalityConstraint().getOutputDimension());
  result_ = OptimizationResult(getProblem());

  t0_ = std::chrono::steady_clock::now();

  // Build x0
  Eigen::VectorXd x0(dimension);
  for (UnsignedInteger i = 0; i < dimension; ++i)
    x0(i) = getStartingPoint()[i];

  // Wrapper for the objective function that tracks evaluations
  const Function probObj(getProblem().getObjective());
  const Bool minimization = getProblem().isMinimization();

  auto objective = [&](const Eigen::VectorXd& x) -> double {
    Point inP(x.size());
    std::copy(x.data(), x.data() + x.size(), inP.begin());
    const Point outP = probObj(inP);
    evaluationInputHistory_.add(inP);
    evaluationOutputHistory_.add(outP);
    return minimization ? outP[0] : -outP[0];
  };

  // Build nonlinear constraint function matching old COBYLA convention
  //   inequality g(x) <= 0  → constraint value = g(x)  (should be <= 0)
  //   equality   h(x) = 0   → constraint value = h(x) + eps and -h(x) + eps  (both <= 0)
  const Function probIneq(getProblem().getInequalityConstraint());
  const Function probEq(getProblem().getEqualityConstraint());
  const bool hasIneq = getProblem().hasInequalityConstraint();
  const bool hasEq = getProblem().hasEqualityConstraint();

  prima::NonlinearConstraintFunction nlcFunc;
  if (hasIneq || hasEq)
  {
    nlcFunc = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
      Point inP(x.size());
      std::copy(x.data(), x.data() + x.size(), inP.begin());

      int totalNlcon = 0;
      if (hasIneq) totalNlcon += probIneq.getOutputDimension();
      if (hasEq) totalNlcon += 2 * probEq.getOutputDimension();

      Eigen::VectorXd result(totalNlcon);
      int idx = 0;

      if (hasIneq)
      {
        const Point ineqVal = probIneq(inP);
        inequalityConstraintHistory_.add(ineqVal);
        for (UnsignedInteger i = 0; i < probIneq.getOutputDimension(); ++i)
          result(idx++) = ineqVal[i];
      }

      if (hasEq)
      {
        const Point eqVal = probEq(inP);
        equalityConstraintHistory_.add(eqVal);
        const Scalar eps = getMaximumConstraintError();
        for (UnsignedInteger i = 0; i < probEq.getOutputDimension(); ++i)
        {
          result(idx++) = eqVal[i] + eps;
          result(idx++) = -eqVal[i] + eps;
        }
      }

      return result;
    };
  }

  // Build bounds
  prima::Bounds* boundsPtr = nullptr;
  prima::Bounds primaBounds;
  if (getProblem().hasBounds())
  {
    const Interval bounds(getProblem().getBounds());
    Eigen::VectorXd lb(dimension);
    Eigen::VectorXd ub(dimension);
    for (UnsignedInteger i = 0; i < dimension; ++i)
    {
      lb(i) = bounds.getFiniteLowerBound()[i] ? bounds.getLowerBound()[i] : -std::numeric_limits<double>::infinity();
      ub(i) = bounds.getFiniteUpperBound()[i] ? bounds.getUpperBound()[i] : std::numeric_limits<double>::infinity();
    }
    primaBounds = prima::Bounds(lb, ub);
    boundsPtr = &primaBounds;
  }

  // Options
  prima::MinimizeOptions opts;
  opts.quiet = true;
  opts.rhobeg = rhoBeg_;
  opts.rhoend = rhoEnd;
  opts.maxfun = maxFun;
  opts.iprint = Log::HasDebug() ? 1 : 0;

  opts.callback = [&](const Eigen::VectorXd& /* x */, double /* f */, int nf, int /* info */, double /* cstrv */, const Eigen::VectorXd& /* nlconstr */) -> bool {
    if (stopCallback_.first && stopCallback_.first(stopCallback_.second))
      return true;

    if (progressCallback_.first)
      progressCallback_.first((100.0 * nf) / maxFun, progressCallback_.second);

    return false;
  };

  const auto primaResult = prima::minimize(
      objective, x0, "cobyla",
      boundsPtr, nullptr,
      (hasIneq || hasEq) ? &nlcFunc : nullptr,
      opts);

  // Map prima status to OT status
  result_.setCallsNumber(primaResult.nfev);
  if (primaResult.status == prima::SMALL_TR_RADIUS || primaResult.status == prima::FTARGET_ACHIEVED || primaResult.status == prima::CALLBACK_TERMINATE)
    result_.setStatus(OptimizationResult::SUCCESS);
  else if (primaResult.status == prima::MAXFUN_REACHED)
    result_.setStatus(OptimizationResult::MAXIMUMCALLS);
  else
    result_.setStatus(OptimizationResult::FAILURE);
  result_.setStatusMessage(primaResult.message);

  setResultFromEvaluationHistory(evaluationInputHistory_, evaluationOutputHistory_, inequalityConstraintHistory_, equalityConstraintHistory_);

  // check for timeout
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  const Scalar timeDuration = std::chrono::duration<Scalar>(t1 - t0_).count();
  result_.setTimeDuration(timeDuration);
  if ((getMaximumTimeDuration() > 0.0) && (timeDuration > getMaximumTimeDuration()))
  {
    result_.setStatus(OptimizationResult::TIMEOUT);
    result_.setStatusMessage(OSS() << "Cobyla optimization timeout after " << timeDuration << "s");
  }

  if (result_.getStatus() != OptimizationResult::SUCCESS)
  {
    if (getCheckStatus())
      throw InternalException(HERE) << "Solving problem by cobyla method failed (" << result_.getStatusMessage() << ")";
    else
      LOGWARN(OSS() << "The Cobyla algorithm failed. The error message is " << result_.getStatusMessage());
  }
}

/* RhoBeg accessor */
Scalar Cobyla::getRhoBeg() const
{
  return rhoBeg_;
}

void Cobyla::setRhoBeg(const Scalar rhoBeg)
{
  rhoBeg_ = rhoBeg;
}

/* String converter */
String Cobyla::__repr__() const
{
  OSS oss;
  oss << "class=" << Cobyla::GetClassName()
      << " " << OptimizationAlgorithmImplementation::__repr__()
      << " rhoBeg=" << rhoBeg_;
  return oss;
}

/* Method save() stores the object through the StorageManager */
void Cobyla::save(Advocate & adv) const
{
  OptimizationAlgorithmImplementation::save(adv);
  adv.saveAttribute("rhoBeg_", rhoBeg_);
}

/* Method load() reloads the object from the StorageManager */
void Cobyla::load(Advocate & adv)
{
  OptimizationAlgorithmImplementation::load(adv);
  adv.loadAttribute("rhoBeg_", rhoBeg_);
}

END_NAMESPACE_OPENTURNS
