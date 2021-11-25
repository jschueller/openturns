//                                               -*- C++ -*-
/**
 *  @brief Metropolis-Hastings algorithm
 *
 *  Copyright 2005-2021 Airbus-EDF-IMACS-ONERA-Phimeca
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

#include "openturns/RandomWalkMetropolisHastings.hxx"
#include "openturns/RandomGenerator.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/ConditionalDistribution.hxx"
#include "openturns/Log.hxx"
#include "openturns/SpecFunc.hxx"

BEGIN_NAMESPACE_OPENTURNS


CLASSNAMEINIT(RandomWalkMetropolisHastings)

static const Factory<RandomWalkMetropolisHastings> Factory_RandomWalkMetropolisHastings;

/* Default constructor */
RandomWalkMetropolisHastings::RandomWalkMetropolisHastings()
  : MetropolisHastingsImplementation()
{
  // Nothing to do
}


/* Parameters constructor */
RandomWalkMetropolisHastings::RandomWalkMetropolisHastings(const Distribution & targetDistribution,
                                                          const Point & initialState,
                                                          const Distribution & instrumental,
                                                           const Indices & marginalIndices)
  : MetropolisHastingsImplementation(targetDistribution, initialState, marginalIndices)
  , adaptationRange_(Interval(ResourceMap::GetAsScalar("RandomWalkMetropolisHastings-DefaultAdaptationLowerBound"), ResourceMap::GetAsScalar("RandomWalkMetropolisHastings-DefaultAdaptationUpperBound")))
  , adaptationExpansionFactor_(ResourceMap::GetAsScalar("RandomWalkMetropolisHastings-DefaultAdaptationExpansionFactor"))
  , adaptationShrinkFactor_(ResourceMap::GetAsScalar("RandomWalkMetropolisHastings-DefaultAdaptationShrinkFactor"))
  , adaptationPeriod_(ResourceMap::GetAsUnsignedInteger("RandomWalkMetropolisHastings-DefaultAdaptationPeriod"))
{
  setInstrumental(instrumental);
}

/* Parameters constructor */
RandomWalkMetropolisHastings::RandomWalkMetropolisHastings(const Function & targetLogPDF,
                                                          const Domain & support,
                                                          const Point & initialState,
                                                          const Distribution & instrumental,
                                                           const Indices & marginalIndices)
  : MetropolisHastingsImplementation(targetLogPDF, support, initialState, marginalIndices)
  , adaptationRange_(Interval(ResourceMap::GetAsScalar("RandomWalkMetropolisHastings-DefaultAdaptationLowerBound"), ResourceMap::GetAsScalar("RandomWalkMetropolisHastings-DefaultAdaptationUpperBound")))
  , adaptationExpansionFactor_(ResourceMap::GetAsScalar("RandomWalkMetropolisHastings-DefaultAdaptationExpansionFactor"))
  , adaptationShrinkFactor_(ResourceMap::GetAsScalar("RandomWalkMetropolisHastings-DefaultAdaptationShrinkFactor"))
  , adaptationPeriod_(ResourceMap::GetAsUnsignedInteger("RandomWalkMetropolisHastings-DefaultAdaptationPeriod"))
{
  setInstrumental(instrumental);
}


/* String converter */
String RandomWalkMetropolisHastings::__repr__() const
{
  return OSS() << "class=" << RandomWalkMetropolisHastings::GetClassName()
         << " name=" << getName()
         << " derived from " << MetropolisHastingsImplementation::__repr__()
         << " instrumental=" << instrumental_;
}


RandomWalkMetropolisHastings* RandomWalkMetropolisHastings::clone() const
{
  return new RandomWalkMetropolisHastings(*this);
}


Point RandomWalkMetropolisHastings::getCandidate() const
{
  // re-adapt if necessary
  if ((samplesNumber_ < getBurnIn()) && ((samplesNumber_ % adaptationPeriod_) == (adaptationPeriod_ - 1)))
  {
    // compute the current acceptation rate
    Scalar rho = 1.0 * acceptedNumberAdaptation_ / (1.0 * adaptationPeriod_);

    if (rho < adaptationRange_.getLowerBound()[0])
      // if the acceptance rate it too low, make smaller steps
      delta_ *= adaptationShrinkFactor_;
    else if (rho > adaptationRange_.getUpperBound()[0])
      // if the acceptance rate is too high, make larger steps
      delta_ *= adaptationExpansionFactor_;

    // reset the counter
    acceptedNumberAdaptation_ = 0;

    if (getVerbose())
      LOGTRACE(OSS() << "rho=" << rho << " delta=" << delta_);
  }

  Point prop(delta_ * instrumental_.getRealization());
  Point newState(currentState_);
  for (UnsignedInteger j = 0; j < marginalIndices_.getSize(); ++ j)
    newState[marginalIndices_[j]] += prop[j];
  return newState;
}


void RandomWalkMetropolisHastings::setInstrumental(const Distribution & instrumental)
{
  if (instrumental.getDimension() != marginalIndices_.getSize())
    throw InvalidArgumentException(HERE) << "The instrumental density dimension does not match the block size.";
  if (!(instrumental.getSkewness().norm() < ResourceMap::GetAsScalar("Distribution-DefaultQuantileEpsilon")))
    throw InvalidArgumentException(HERE) << "The instrumental density is not symmetric.";
  if (!(instrumental.getMean().norm() < ResourceMap::GetAsScalar("Distribution-DefaultQuantileEpsilon")))
    throw InvalidArgumentException(HERE) << "The instrumental density must have a null mean.";
  instrumental_ = instrumental;
}


Distribution RandomWalkMetropolisHastings::getInstrumental() const
{
  return instrumental_;
}


void RandomWalkMetropolisHastings::setAdaptationRange(const Interval & adaptationRange)
{
  if (adaptationRange.getDimension() != 1)
    throw InvalidDimensionException(HERE) << "Range should be 1-d. Got " << adaptationRange.getDimension();
  adaptationRange_ = adaptationRange;
}

Interval RandomWalkMetropolisHastings::getAdaptationRange() const
{
  return adaptationRange_;
}


void RandomWalkMetropolisHastings::setAdaptationExpansionFactor(Scalar adaptationExpansionFactor)
{
  if (!(adaptationExpansionFactor > 1.0))
    throw InvalidArgumentException(HERE) << "Expansion factor should be > 1. Got " << adaptationExpansionFactor;
  adaptationExpansionFactor_ = adaptationExpansionFactor;
}

Scalar RandomWalkMetropolisHastings::getAdaptationExpansionFactor() const
{
  return adaptationExpansionFactor_;
}


void RandomWalkMetropolisHastings::setAdaptationShrinkFactor(Scalar adaptationShrinkFactor)
{
  if (!(adaptationShrinkFactor > 0.0) || !(adaptationShrinkFactor < 1.0))
    throw InvalidArgumentException(HERE) << "Shrink factor should be in (0, 1). Got " << adaptationShrinkFactor;
  adaptationShrinkFactor_ = adaptationShrinkFactor;
}

Scalar RandomWalkMetropolisHastings::getAdaptationShrinkFactor() const
{
  return adaptationShrinkFactor_;
}

void RandomWalkMetropolisHastings::setAdaptationPeriod(const UnsignedInteger adaptationPeriod)
{
  if (!adaptationPeriod)
    throw InvalidArgumentException(HERE) << "The adaptation period must be positive.";
  adaptationPeriod_ = adaptationPeriod;
}


UnsignedInteger RandomWalkMetropolisHastings::getAdaptationPeriod() const
{
  return adaptationPeriod_;
}

/* Method save() stores the object through the StorageManager */
void RandomWalkMetropolisHastings::save(Advocate & adv) const
{
  MetropolisHastingsImplementation::save(adv);
  adv.saveAttribute("instrumental_", instrumental_);
  adv.saveAttribute("delta_", delta_);
  adv.saveAttribute("adaptationRange_", adaptationRange_);
  adv.saveAttribute("adaptationExpansionFactor_", adaptationExpansionFactor_);
  adv.saveAttribute("adaptationShrinkFactor_", adaptationShrinkFactor_);
  adv.saveAttribute("adaptationPeriod_", adaptationPeriod_);
}

/* Method load() reloads the object from the StorageManager */
void RandomWalkMetropolisHastings::load(Advocate & adv)
{
  MetropolisHastingsImplementation::load(adv);
  adv.loadAttribute("instrumental_", instrumental_);
  adv.loadAttribute("delta_", delta_);
  adv.loadAttribute("adaptationRange_", adaptationRange_);
  adv.loadAttribute("adaptationExpansionFactor_", adaptationExpansionFactor_);
  adv.loadAttribute("adaptationShrinkFactor_", adaptationShrinkFactor_);
  adv.loadAttribute("adaptationPeriod_", adaptationPeriod_);
}


END_NAMESPACE_OPENTURNS
