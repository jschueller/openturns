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
#ifndef OPENTURNS_ADAM_HXX
#define OPENTURNS_ADAM_HXX

#include "openturns/OTprivate.hxx"
#include "openturns/OptimizationAlgorithmImplementation.hxx"
#include "openturns/OptimizationProblem.hxx"
#include "openturns/Point.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class Adam
 */
class OT_API Adam :
  public OptimizationAlgorithmImplementation
{
  CLASSNAME
public:


  /** Default constructor */
  Adam();

  /** Constructor with parameters */
  explicit Adam(const OptimizationProblem & problem);

  /** Virtual constructor */
  Adam * clone() const override;

  /** Performs the actual computation. Must be overloaded by the actual optimisation algorithm */
  void run() override;

  /** Alpha accessor */
  Scalar getAlpha() const;
  void setAlpha(const Scalar alpha);

  /** Beta1 accessor */
  Scalar getBeta1() const;
  void setBeta1(const Scalar beta1);

  /** Beta2 accessor */
  Scalar getBeta2() const;
  void setBeta2(const Scalar beta1);
  
  /** Smooth accessor */
  Scalar getEpsilon() const;
  void setEpsilon(const Scalar tau);

  /** String converter */
  String __repr__() const override;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;

protected:
  /** Check whether this problem can be solved by this solver.  Must be overloaded by the actual optimisation algorithm */
  void checkProblem(const OptimizationProblem & problem) const override;

private:
  /** Multiplicative decrease of linear step */
  Scalar alpha_ = 1e-3;

  /** Armijo factor */
  Scalar beta1_ = 0.9;

  /** Growing factor in penalization term */
  Scalar beta2_ = 0.999;

  /** Abdo Rackwitz current penalization factor */
  Scalar epsilon_ = 1e-8;

}; /* class Adam */

END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_ADAM_HXX */
