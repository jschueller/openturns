//                                               -*- C++ -*-
/**
 *  @brief HiGHS linear solver
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
#include "openturns/HiGHS.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/OTconfig.hxx"

#ifdef OPENTURNS_HAVE_HIGHS
#include "Highs.h"
#endif

#include <chrono>

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(HiGHS)

static const Factory<HiGHS> Factory_HiGHS;

/* Constructor with no parameters */
HiGHS::HiGHS()
  : OptimizationAlgorithmImplementation()
{
  // Nothing to do here
}


/* Constructor that sets starting sample */
HiGHS::HiGHS(const OptimizationProblem & problem)
  : OptimizationAlgorithmImplementation(problem)

{
}


/* Check whether this problem can be solved by this solver */
// Inherited method: never used by HiGHS because solver_ calls its own checkProblem method.
void HiGHS::checkProblem(const OptimizationProblem & ) const
{
  // Nothing to do
}


void HiGHS::run()
{
#ifdef OPENTURNS_HAVE_HIGHS
  const UnsignedInteger problemDimension = getProblem().getDimension();
  if (problemDimension == 0) throw InvalidArgumentException(HERE) << "No problem has been set.";

  HighsModel model;
  model.lp_.num_col_ = problemDimension;
  model.lp_.num_row_ = getProblem().getInequalityConstraint().getOutputDimension();
  model.lp_.sense_ = getProblem().isMinimization() ? ObjSense::kMinimize : ObjSense::kMaximize;
  const Scalar f0 = getProblem().getObjective()(Point(problemDimension, 0.0))[0];
  model.lp_.offset_ = f0;
  const Matrix grad0(getProblem().getObjective().gradient((Point(problemDimension, 0.0))));
  model.lp_.col_cost_.resize(model.lp_.num_col_);
  for (UnsignedInteger i = 0; i < problemDimension; ++ i)
    model.lp_.col_cost_[i] = grad0(i, 0);
  model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;

  model.lp_.integrality_.resize(model.lp_.num_col_);
  for (int col = 0; col < model.lp_.num_col_; ++ col)
    if (getProblem().getVariablesType()[col] == OptimizationProblemImplementation::CONTINUOUS)
      model.lp_.integrality_[col] = HighsVarType::kContinuous;
    else // INTEGER or BOOL
      model.lp_.integrality_[col] = HighsVarType::kInteger;

  // Create a Highs instance
  Highs highs;
  // HighsStatus return_status;
  //
  // Pass the model to HiGHS
  HighsStatus return_status = highs.passModel(model);
  if (return_status != HighsStatus::kOk)
    throw InvalidArgumentException(HERE) << "Cannot initialize highs model";
  // assert(return_status == HighsStatus::kOk);
  // If a user passes a model with entries in
  // model.lp_.a_matrix_.value_ less than (the option)
  // small_matrix_value in magnitude, they will be ignored. A logging
  // message will indicate this, and passModel will return
  // HighsStatus::kWarning
  //
  // Get a const reference to the LP data in HiGHS
  // const HighsLp& lp = highs.getLp();

  // Solve the model
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  return_status = highs.run();
  if (return_status != HighsStatus::kOk)
    throw InvalidArgumentException(HERE) << "Solve not ok";
  
  //
  // Get the model status
  const HighsModelStatus& model_status = highs.getModelStatus();
  if (model_status != HighsModelStatus::kOptimal)
    throw InvalidArgumentException(HERE) << "Solve not optimal";

  // cout << "Model status: " << highs.modelStatusToString(model_status) << endl;
  //
  // Get the solution information
  const HighsInfo& info = highs.getInfo();
  // cout << "Simplex iteration count: " << info.simplex_iteration_count << endl;
  // cout << "Objective function value: " << info.objective_function_value << endl;
  // cout << "Primal  solution status: "
  //      << highs.solutionStatusToString(info.primal_solution_status) << endl;
  // cout << "Dual    solution status: "
  //      << highs.solutionStatusToString(info.dual_solution_status) << endl;
  // cout << "Basis: " << highs.basisValidityToString(info.basis_validity) << endl;
  // const bool has_values = info.primal_solution_status;
  // const bool has_duals = info.dual_solution_status;
  // const bool has_basis = info.basis_validity;
  //
  // Get the solution values and basis
  // const HighsSolution& solution = highs.getSolution();
  // const HighsBasis& basis = highs.getBasis();
  //
  // Report the primal and solution values and basis
  // for (int col = 0; col < lp.num_col_; col++) {
  //   cout << "Column " << col;
  //   if (has_values) cout << "; value = " << solution.col_value[col];
  //   if (has_duals) cout << "; dual = " << solution.col_dual[col];
  //   if (has_basis)
  //     cout << "; status: " << highs.basisStatusToString(basis.col_status[col]);
  //   cout << endl;
  // }
  // for (int row = 0; row < lp.num_row_; row++) {
  //   cout << "Row    " << row;
  //   if (has_values) cout << "; value = " << solution.row_value[row];
  //   if (has_duals) cout << "; dual = " << solution.row_dual[row];
  //   if (has_basis)
  //     cout << "; status: " << highs.basisStatusToString(basis.row_status[row]);
  //   cout << endl;
  // }

  const UnsignedInteger initialCallsNumber = getProblem().getObjective().getCallsNumber();
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    Scalar timeDuration = std::chrono::duration<Scalar>(t1 - t0).count();
    result_.setTimeDuration(timeDuration);

    const UnsignedInteger callsNumber = getProblem().getObjective().getCallsNumber() - initialCallsNumber;
    result_.setCallsNumber(callsNumber);

    if (callsNumber > getMaximumCallsNumber())
    {
      result_.setStatus(OptimizationResult::MAXIMUMCALLS);
      result_.setStatusMessage(OSS() << "HiGHS reaches maximum calls number after " << callsNumber << " calls");
    }

    if ((getMaximumTimeDuration() > 0.0) && (timeDuration > getMaximumTimeDuration()))
    {
      result_.setStatus(OptimizationResult::TIMEOUT);
      result_.setStatusMessage(OSS() << "HiGHS optimization timeout after " << timeDuration << "s");
    }
/*
    // callbacks
    if (progressCallback_.first)
    {
      progressCallback_.first((100.0 * callsNumber) / getMaximumCallsNumber(), progressCallback_.second);
    }
    if (stopCallback_.first)
    {
      Bool stop = stopCallback_.first(stopCallback_.second);
      if (stop)
      {
        LOGWARN(OSS() << "HiGHS was stopped by user");
        break;
      }
    }
  }*/
#else
  throw NotYetImplementedException(HERE) << "No HiGHS support";
#endif
}


/* Virtual constructor */
HiGHS * HiGHS::clone() const
{
  return new HiGHS(*this);
}

/* String converter */
String HiGHS::__repr__() const
{
  OSS oss;
  oss << "class=" << getClassName()
      << " " << OptimizationAlgorithmImplementation::__repr__();
  return oss;
}



/* Method save() stores the object through the StorageManager */
void HiGHS::save(Advocate & adv) const
{
  OptimizationAlgorithmImplementation::save(adv);
  // adv.saveAttribute("solver_", solver_);
}

/* Method load() reloads the object from the StorageManager */
void HiGHS::load(Advocate & adv)
{
  OptimizationAlgorithmImplementation::load(adv);
  // adv.loadAttribute("solver_", solver_);
}

END_NAMESPACE_OPENTURNS
