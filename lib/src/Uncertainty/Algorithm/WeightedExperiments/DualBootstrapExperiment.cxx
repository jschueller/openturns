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
#include "openturns/DualBootstrapExperiment.hxx"
#include "openturns/UserDefined.hxx"
#include "openturns/RandomGenerator.hxx"
#include "openturns/PersistentObjectFactory.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(DualBootstrapExperiment)

static const Factory<DualBootstrapExperiment> Factory_DualBootstrapExperiment;

/* Default constructor */
DualBootstrapExperiment::DualBootstrapExperiment():
  PersistentObject()
{
  // Nothing to do
}

/* Constructor with parameters */
DualBootstrapExperiment::DualBootstrapExperiment(const Sample & sample1, const Sample & sample2)
  : PersistentObject()
  , size_(sample1.getSize())
  , sample1_(sample1)
  , sample2_(sample2)
{
  if (sample2.getSize() != sample1.getSize())
    throw InvalidArgumentException(HERE) << "sample sizes are not equal";
}

/* Virtual constructor */
DualBootstrapExperiment * DualBootstrapExperiment::clone() const
{
  return new DualBootstrapExperiment(*this);
}

/* String converter */
String DualBootstrapExperiment::__repr__() const
{
  OSS oss;
  oss << "class=" << GetClassName()
      << " name=" << getName()
      << " sample1=" << sample1_
      << " sample2=" << sample2_
      << " size=" << size_;
  return oss;
}

Sample DualBootstrapExperiment::generate(Sample & sample2) const
{
  const Indices selection(GenerateSelection(size_, size_));
  sample2 = sample2_.select(selection);
  return sample1_.select(selection);
}


/* Selection generation */
Indices DualBootstrapExperiment::GenerateSelection(const UnsignedInteger size,
    const UnsignedInteger length)
{
  const RandomGenerator::UnsignedIntegerCollection selection(RandomGenerator::IntegerGenerate(size, length));
  return Indices(selection.begin(), selection.end());
}

END_NAMESPACE_OPENTURNS
