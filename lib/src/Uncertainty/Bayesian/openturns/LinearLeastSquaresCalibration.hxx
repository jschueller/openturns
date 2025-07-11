//                                               -*- C++ -*-
/**
 *  @brief LinearLeastSquaresCalibration algorithm
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
#ifndef OPENTURNS_LINEARLEASTSQUARESCALIBRATION_HXX
#define OPENTURNS_LINEARLEASTSQUARESCALIBRATION_HXX

#include "openturns/OTprivate.hxx"
#include "openturns/CalibrationAlgorithmImplementation.hxx"
#include "openturns/Sample.hxx"
#include "openturns/Function.hxx"
#include "openturns/LeastSquaresMethod.hxx"
#include "openturns/ResourceMap.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class LINEARLEASTSQUARESCALIBRATION
 *
 * @brief The class implements the linear least-squares calibration
 *
 */
class OT_API LinearLeastSquaresCalibration
  : public CalibrationAlgorithmImplementation
{
  CLASSNAME

  /** The NonLinearLeastSquaresCalibration class is closely related with the LinearLeastSquaresCalibration class */
  friend class NonLinearLeastSquaresCalibration;

public:

  /** Default constructor */
  LinearLeastSquaresCalibration();

  /** Parameter constructor */
  LinearLeastSquaresCalibration(const Function & model,
                                const Sample & inputObservations,
                                const Sample & outputObservations,
                                const Point & startingPoint,
                                const String & methodName = ResourceMap::GetAsString("LinearLeastSquaresCalibration-Method"));

  /** Parameter constructor */
  LinearLeastSquaresCalibration(const Sample & modelObservations,
                                const Matrix & gradientObservations,
                                const Sample & outputObservations,
                                const Point & startingPoint,
                                const String & methodName = ResourceMap::GetAsString("LinearLeastSquaresCalibration-Method"));

  /** String converter */
  String __repr__() const override;

  /** Performs the actual computation. Must be overloaded by the actual calibration algorithm */
  void run() override;

  /** Model observations accessor */
  Sample getModelObservations() const;

  /** Model gradient wrt the parameter accessor */
  Matrix getGradientObservations() const;

  /** StartingPoint accessor */
  Point getStartingPoint() const;

  /** Least squares method name accessor */
  String getMethodName() const;

  /* Here is the interface that all derived class must implement */

  /** Virtual constructor */
  LinearLeastSquaresCalibration * clone() const override;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;

private:
  /* Compute the posterior and error distributions */
  static Distribution ComputePosteriorAndErrorDistribution(const Point & thetaStar,
      const Point & r,
      const SquareMatrix & gramInverse,
      const UnsignedInteger outputDimension,
      Distribution & error);

  /* The model observations */
  Sample modelObservations_;

  /* The model gradient wrt the parameter */
  Matrix gradientObservations_;

  /* The least squares method name */
  String methodName_;

}; /* class LINEARLEASTSQUARESCALIBRATION */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_LINEARLEASTSQUARESCALIBRATION_HXX */
