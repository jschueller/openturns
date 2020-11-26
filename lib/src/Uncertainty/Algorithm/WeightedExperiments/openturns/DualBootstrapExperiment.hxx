//                                               -*- C++ -*-
/**
 *  @brief Bootstrap with two samples
 *
 *  Copyright 2005-2020 Airbus-EDF-IMACS-ONERA-Phimeca
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
#ifndef OPENTURNS_DUALBOOTSTRAPEXPERIMENT_HXX
#define OPENTURNS_DUALBOOTSTRAPEXPERIMENT_HXX

#include "openturns/WeightedExperimentImplementation.hxx"

BEGIN_NAMESPACE_OPENTURNS



/**
 * @class DualBootstrapExperiment
 *
 * The class describes the probabilistic concept of bootstrapExperiment plan
 */
class OT_API DualBootstrapExperiment
  : public PersistentObject
{
  CLASSNAME
public:


  /** Default constructor */
  DualBootstrapExperiment();

  /** Parameters constructor */
  DualBootstrapExperiment(const Sample & sample1, const Sample & sample2);

  /** Virtual constructor */
  DualBootstrapExperiment * clone() const override;

  /** String converter */
  String __repr__() const override;

  /** Generate new samples */
  virtual Sample generate(Sample & sample2Out) const;

  /** Selection generation */
  static Indices GenerateSelection(const UnsignedInteger size,
                                   const UnsignedInteger length);

protected:

private:
  /** The size of the sample to be generated */
  UnsignedInteger size_;

  /** The samples from which we resample */
  Sample sample1_;
  Sample sample2_;

}; /* class DualBootstrapExperiment */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_DUALBOOTSTRAPEXPERIMENT_HXX */
