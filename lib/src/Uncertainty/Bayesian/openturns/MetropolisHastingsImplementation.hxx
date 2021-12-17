//                                               -*- C++ -*-
/**
 *  @brief MetropolisHastingsImplementation base class
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
#ifndef OPENTURNS_METROPOLISHASTINGSIMPLEMENTATION_HXX
#define OPENTURNS_METROPOLISHASTINGSIMPLEMENTATION_HXX

#include "openturns/OTprivate.hxx"
#include "openturns/RandomVectorImplementation.hxx"
#include "openturns/ResourceMap.hxx"
#include "openturns/Distribution.hxx"
#include "openturns/HistoryStrategy.hxx"

BEGIN_NAMESPACE_OPENTURNS

/**
 * @class MetropolisHastingsImplementation
 *
 * @brief MetropolisHastingsImplementation base class
 *
 * This class is abstract so it can not be instanciated. It must be derived.
 */
class OT_API MetropolisHastingsImplementation
  : public RandomVectorImplementation
{
  CLASSNAME
  friend class Gibbs;
public:

  /** Default constructor */
  MetropolisHastingsImplementation();

  /** Constructor with parameters*/
  MetropolisHastingsImplementation(const Distribution & targetDistribution,
                                   const Point & initialState,
                                   const Indices & marginalIndices = Indices());

  /** Constructor with parameters*/
  MetropolisHastingsImplementation(const Function & targetLogPDF,
                                   const Domain & support,
                                   const Point & initialState,
                                   const Indices & marginalIndices = Indices());

  /** Virtual constructor */
  MetropolisHastingsImplementation * clone() const override;

  /** String converter */
  String __repr__() const override;

  /** Get a realization */
  Point getRealization() const override;

  /** Compute the log-likelihood */
  virtual Scalar computeLogLikelihood(const Point & state) const;

  /** Compute the log-likelihood */
  virtual Scalar computeLogPosterior(const Point & state) const;

  /** Target distribution accessor */
  virtual Distribution getTargetDistribution() const;

  /** Target log-pdf accessor */
  virtual Function getTargetLogPDF() const;
  virtual Domain getTargetLogPDFSupport() const;

  /** Initial state accessor */
  virtual Point getInitialState() const;

  /** Marginal indices accessor */
  virtual Indices getMarginalIndices() const;

  /** Likelihood accessor */
  virtual void setLikelihood(const Distribution & conditional,
                            const Sample & observations,
                            const Function & linkFunction = Function(),
                            const Sample & covariates = Sample(0, 0));

  /** Conditional accessor */
  Distribution getConditional() const;

  /** Link function accessor */
  Function getLinkFunction() const;

  /** Obervations accessor */
  Sample getObservations() const;

  /** Covariates accessor */
  Sample getCovariates() const;

  /** Burnin accessor */
  void setBurnIn(UnsignedInteger burnIn);
  UnsignedInteger getBurnIn() const;

  /** Thinning accessor */
  void setThinning(UnsignedInteger thinning);
  UnsignedInteger getThinning() const;

  /** History accessor */
  void setHistory(const HistoryStrategy & strategy);
  HistoryStrategy getHistory() const;

  /** Dimension accessor */
  UnsignedInteger getDimension() const override;

  /** Verbosity accessor */
  void setVerbose(const Bool verbose);
  Bool getVerbose() const;

  /** Acceptance rate accessor*/
  Scalar getAcceptanceRate() const;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;

protected:
  Point initialState_;
  mutable Point currentState_;
  Indices marginalIndices_;
  mutable HistoryStrategy history_;

  // total number of samples
  mutable UnsignedInteger samplesNumber_ = 0;

  // number of accepted candidates
  mutable UnsignedInteger acceptedNumber_ = 0;

  // number of accepted candidates since adaptation
  mutable UnsignedInteger acceptedNumberAdaptation_ = 0;

  /** Propose a new point in the chain */
  virtual Point getCandidate() const;

  // Copy the attributes of an other instance
  void copyAttributes(const MetropolisHastingsImplementation * mhi);

private:
  // target distribution
  Distribution targetDistribution_;

  // ... when defined via its log-pdf
  Function targetLogPDF_;
  Domain support_;

  // likelihood term
  Distribution conditional_;
  Function linkFunction_;
  Sample covariates_;
  Sample observations_;

  // number of first samples discarded to reach stationary regime
  UnsignedInteger burnIn_ = 0;

  // number of samples skipped at each generation
  UnsignedInteger thinning_ = 0;

  // verbosity flag
  Bool verbose_ = false;

  // unnormalized log-posterior density of the current state
  mutable Scalar currentLogPosterior_ = 0.0;

  // prior log pdf
  Scalar computeLogPDFPrior(const Point & state) const;

  void setTargetDistribution(const Distribution & targetDistribution);
  void setTargetLogPDF(const Function & targetLogPDF, const Domain & support);

}; /* class MetropolisHastingsImplementation */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_METROPOLISHASTINGSIMPLEMENTATION_HXX */
