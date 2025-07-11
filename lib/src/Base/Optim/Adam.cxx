//                                               -*- C++ -*-
/**
 *  @brief Adam optimization solver
 *
 *  Copyright 2005-2025 Airbus-EDF-IMACS-ONERA-Phimeca
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
#include <cmath>
#include "openturns/Adam.hxx"
#include "openturns/Log.hxx"
#include "openturns/PersistentObjectFactory.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(Adam)

static const Factory<Adam> Factory_Adam;

/* Default constructor */
Adam::Adam()
  : OptimizationAlgorithmImplementation()
{
}

Adam::Adam (const OptimizationProblem & problem)
  : OptimizationAlgorithmImplementation(problem)
{
  checkProblem(problem);
}

/* Virtual constructor */
Adam * Adam::clone() const
{
  return new Adam(*this);
}


/* Check whether this problem can be solved by this solver.  Must be overloaded by the actual optimisation algorithm */
void Adam::checkProblem(const OptimizationProblem & problem) const
{
  if (problem.hasMultipleObjective())
    throw InvalidArgumentException(HERE) << "Adam does not support multi-objective optimization";

  if (problem.hasResidualFunction())
    throw InvalidArgumentException(HERE) << "Adam does not support least squares problems";

  if (!problem.isContinuous())
    throw InvalidArgumentException(HERE) << "Adam does not support non continuous problems";
}

/* Performs the actual computation by using the Adam algorithm
 */
void Adam::run()
{
  // initialize();

  /* Get a local copy of the level function */
  const Function levelFunction(getProblem().getObjective());
  /* Get a local copy of the level value */
  // const Scalar levelValue = getProblem().getLevelValue();
  /* Current point -> u */
  // currentPoint_ = getStartingPoint();
  Bool stop = false;
  UnsignedInteger iterationNumber = 0;
  const UnsignedInteger initialEvaluationNumber = levelFunction.getEvaluationCallsNumber();
  const UnsignedInteger dimension = getProblem().getDimension();
  
  Scalar absoluteError = -1.0;
  Scalar constraintError = -1.0;
  Scalar relativeError = -1.0;
  Scalar residualError = -1.0;

  /* Compute the level function at the current point -> G */
  // currentLevelValue_ = levelFunction(currentPoint_)[0];

  UnsignedInteger evaluationNumber = levelFunction.getEvaluationCallsNumber() - initialEvaluationNumber;

  // reset result
  result_ = OptimizationResult(getProblem());
  // result_.store(currentPoint_, Point(1, currentLevelValue_), absoluteError, relativeError, residualError, constraintError);

  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

  Point m(dimension);
  Point v(dimension);
  Point x(getStartingPoint());
  Point xPrev(x);
  Point y;
  Point yPrev(levelFunction(x));
  
  if (getProblem().hasBounds() && !getProblem().getBounds().contains(x))
  {
    const Point lowerBound(getProblem().getBounds().getLowerBound());
    const Point upperBound(getProblem().getBounds().getUpperBound());
    const Interval::BoolCollection finiteLowerBound(getProblem().getBounds().getFiniteLowerBound());
    const Interval::BoolCollection finiteUpperBound(getProblem().getBounds().getFiniteUpperBound());
    for (UnsignedInteger i = 0; i < dimension; ++ i)
    {
      if (finiteLowerBound[i] && (x[i] < lowerBound[i]))
        x[i] = lowerBound[i];
      if (finiteUpperBound[i] && (x[i] > upperBound[i]))
        x[i] = upperBound[i];
    }
  }
  
  while ((!stop) && (iterationNumber <= getMaximumIterationNumber()))
  {
    /* Go to next iteration */
    ++ iterationNumber;

    /* Compute the level function gradient at the current point -> Grad(G) */
    Point g(levelFunction.gradient(x) * Point(1, 1.0));
    LOGDEBUG(OSS() << "current point=" << x << " current level value=" << y << " current gradient=" << g);
    /* Compute the current Lagrange multiplier */
    const Scalar normGradientSquared = g.normSquare();
    /* In case of a null gradient, throw an internal exception */
    if (!(normGradientSquared > 0.0))
      throw InternalException(HERE) << "Adam algorithm: the gradient of the level function is zero at point x=" << x;

    for (UnsignedInteger i = 0; i < dimension; ++ i)
    {
      m[i] = beta1_ * m[i] + (1.0 - beta1_) * g[i];
      v[i] = beta2_ * v[i] + (1.0 - beta2_) * g[i]*g[i];

      const Scalar mhat = m[i] / (1.0 - std::pow(beta1_, iterationNumber + 1.0));
      const Scalar vhat = v[i] / (1.0 - std::pow(beta2_, iterationNumber + 1.0));
      
      x[i] -= alpha_ * mhat / (std::sqrt(vhat) + epsilon_);
    }

    if (getProblem().hasBounds() && !getProblem().getBounds().contains(x))
    {
      const Point lowerBound(getProblem().getBounds().getLowerBound());
      const Point upperBound(getProblem().getBounds().getUpperBound());
      const Interval::BoolCollection finiteLowerBound(getProblem().getBounds().getFiniteLowerBound());
      const Interval::BoolCollection finiteUpperBound(getProblem().getBounds().getFiniteUpperBound());
      for (UnsignedInteger i = 0; i < dimension; ++ i)
      {
        if (finiteLowerBound[i] && (x[i] < lowerBound[i]))
          x[i] = lowerBound[i];
        if (finiteUpperBound[i] && (x[i] > upperBound[i]))
          x[i] = upperBound[i];
      }
    }
    
    y = levelFunction(x);

    // update number of evaluations
    evaluationNumber = levelFunction.getEvaluationCallsNumber() - initialEvaluationNumber;

    std::cout << "evaluationNumber=" << evaluationNumber << std::endl;
    
    /* Check if convergence has been achieved */
    absoluteError = (x - xPrev).normInf();
    // constraintError = std::abs(currentLevelValue_ - levelValue);
    const Scalar pointNorm = x.norm();
    relativeError = pointNorm > 0.0 ? absoluteError / pointNorm : -1.0;
    
    residualError = (y - yPrev).norm();
    residualError = y.norm() > 0.0 ? residualError / y.norm() : -1.0;

    constraintError = 0.0;
    // if (getProblem().hasBounds() && !getProblem().getBounds().contains(x))
    // {
    //   const Point lowerBound(getProblem().getBounds().getLowerBound());
    //   const Point upperBound(getProblem().getBounds().getUpperBound());
    //   const Interval::BoolCollection finiteLowerBound(getProblem().getBounds().getFiniteLowerBound());
    //   const Interval::BoolCollection finiteUpperBound(getProblem().getBounds().getFiniteUpperBound());
    //   for (UnsignedInteger i = 0; i < dimension; ++ i)
    //   {
    //     if (finiteLowerBound[i] && (x[i] < lowerBound[i]))
    //       constraintError = std::max(constraintError, lowerBound[i] - x[i]);
    //     if (finiteUpperBound[i] && (x[i] > upperBound[i]))
    //       constraintError = std::max(constraintError, x[i] - upperBound[i]);
    //   }
    // }

    stop = ((absoluteError < getMaximumAbsoluteError()) && (relativeError < getMaximumRelativeError())) || ((residualError < getMaximumResidualError()) && (constraintError < getMaximumConstraintError()));

    xPrev = x;
    yPrev = y;
    
    // update result
    result_.setCallsNumber(evaluationNumber);
    result_.setIterationNumber(iterationNumber);
    result_.store(x, y, absoluteError, relativeError, residualError, constraintError);

    if (evaluationNumber > getMaximumCallsNumber())
    {
      stop = true;
      result_.setStatus(OptimizationResult::MAXIMUMCALLS);
      result_.setStatusMessage(OSS() << "Adam reached the maximum calls number after " << evaluationNumber << " evaluations");
    }

    // callbacks
    if (progressCallback_.first)
    {
      progressCallback_.first((100.0 * evaluationNumber) / getMaximumCallsNumber(), progressCallback_.second);
    }
    if (stopCallback_.first && stopCallback_.first(stopCallback_.second))
    {
      stop = true;
      LOGTRACE(OSS() << "Adam was stopped by user");
      result_.setStatus(OptimizationResult::INTERRUPTION);
      result_.setStatusMessage(OSS() << "Adam was stopped by user");
    }

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    const Scalar timeDuration = std::chrono::duration<Scalar>(t1 - t0).count();
    if ((getMaximumTimeDuration() > 0.0) && (timeDuration > getMaximumTimeDuration()))
    {
      LOGTRACE("Optim timeout");
      stop = true;
      result_.setStatus(OptimizationResult::TIMEOUT);
      result_.setStatusMessage(OSS() << "Adam optimization timeout after " << timeDuration << "s");
    }
  }

  std::cout << "result" << result_ << std::endl;
  
  if (result_.getStatus() != OptimizationResult::SUCCESS)
  {
    if (getCheckStatus())
      throw InternalException(HERE) << "Adam raised an exception: " << result_.getStatusMessage();
    else
      LOGWARN(OSS() << "Adam algorithm failed: " << result_.getStatusMessage());
  }
} // run()


/* Alpha accessor */
Scalar Adam::getAlpha() const
{
  return alpha_;
}

void Adam::setAlpha(const Scalar alpha)
{
  alpha_ = alpha;
}

/* Beta1 accessor */
Scalar Adam::getBeta1() const
{
  return beta1_;
}

void Adam::setBeta1(const Scalar beta1)
{
  beta1_ = beta1;
}

/* Beta2 accessor */
Scalar Adam::getBeta2() const
{
  return beta2_;
}

void Adam::setBeta2(const Scalar beta2)
{
  beta2_ = beta2;
}

/* Epsilon accessor */
Scalar Adam::getEpsilon() const
{
  return epsilon_;
}

void Adam::setEpsilon(const Scalar epsilon)
{
  epsilon_ = epsilon;
}

/* String converter */
String Adam::__repr__() const
{
  OSS oss;
  oss << "class=" << Adam::GetClassName()
      << " " << OptimizationAlgorithmImplementation::__repr__()
      << " alpha_=" << alpha_
      << " beta1_=" << beta1_
      << " beta2_=" << beta2_;
  return oss;
}

/* Method save() stores the object through the StorageManager */
void Adam::save(Advocate & adv) const
{
  OptimizationAlgorithmImplementation::save(adv);
  adv.saveAttribute("alpha_", alpha_);
  adv.saveAttribute("beta1_", beta1_);
  adv.saveAttribute("beta2_", beta2_);
  adv.saveAttribute("epsilon_", epsilon_);
}

/* Method load() reloads the object from the StorageManager */
void Adam::load(Advocate & adv)
{
  OptimizationAlgorithmImplementation::load(adv);
  adv.loadAttribute("alpha_", alpha_);
  adv.loadAttribute("beta1_", beta1_);
  adv.loadAttribute("beta2_", beta2_);
  adv.loadAttribute("epsilon_", epsilon_);
}

END_NAMESPACE_OPENTURNS
