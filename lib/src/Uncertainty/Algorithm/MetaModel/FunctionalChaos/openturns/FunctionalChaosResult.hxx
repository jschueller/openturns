//                                               -*- C++ -*-
/**
 *  @brief The result of a chaos expansion
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
#ifndef OPENTURNS_FUNCTIONALCHAOSRESULT_HXX
#define OPENTURNS_FUNCTIONALCHAOSRESULT_HXX

#include "openturns/MetaModelResult.hxx"
#include "openturns/Point.hxx"
#include "openturns/Sample.hxx"
#include "openturns/Indices.hxx"
#include "openturns/Collection.hxx"
#include "openturns/PersistentCollection.hxx"
#include "openturns/Function.hxx"
#include "openturns/Distribution.hxx"
#include "openturns/OrthogonalBasis.hxx"

BEGIN_NAMESPACE_OPENTURNS



/**
 * @class FunctionalChaosResult
 *
 * The result of a chaos expansion
 */

class OT_API FunctionalChaosResult
  : public MetaModelResult
{
  CLASSNAME

public:

  typedef Collection<Function>           FunctionCollection;
  typedef PersistentCollection<Function> FunctionPersistentCollection;

  /** Default constructor */
  FunctionalChaosResult();

  /** Parameter constructor */
  FunctionalChaosResult(const Sample & inputSample,
                        const Sample & outputSample,
                        const Distribution & distribution,
                        const Function & transformation,
                        const Function & inverseTransformation,
                        const OrthogonalBasis & orthogonalBasis,
                        const Indices & I,
                        const Sample & alpha_k,
                        const FunctionCollection & Psi_k);

  /** Virtual constructor */
  FunctionalChaosResult * clone() const override;

  /** String converter */
  String __repr__() const override;
  String __str__(const String & offset = "") const override;
  String __repr_markdown__() const override;

  /** Distribution accessor */
  virtual Distribution getDistribution() const;

  /** IsoProbabilisticTransformation accessor */
  virtual Function getTransformation() const;

  /** InverseIsoProbabilisticTransformation accessor */
  virtual Function getInverseTransformation() const;

  /** Orthogonal basis accessor */
  virtual OrthogonalBasis getOrthogonalBasis() const;

  /** Indices accessor */
  virtual Indices getIndices() const;

  /** Coefficients accessor */
  virtual Sample getCoefficients() const;

  /** Reduced basis accessor */
  virtual FunctionCollection getReducedBasis() const;

  /** Composed meta model accessor */
  virtual Function getComposedMetaModel() const;

  /** Residuals accessor */
  virtual Sample getSampleResiduals() const;

  /** isLeastSquares_ accessor */
  virtual Bool isLeastSquares() const;

  /** involvesModelSelection accessor */
  virtual Bool involvesModelSelection() const;

  /** isLeastSquares_ accessor */
  virtual void setIsLeastSquares(const Bool isLeastSquares);

  /** involvesModelSelection accessor */
  virtual void setInvolvesModelSelection(const Bool involvesModelSelection);

  /** Conditional expectation accessor */
  virtual FunctionalChaosResult getConditionalExpectation(const Indices & conditioningIndices) const;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;

  /** Selection history accessor */
  IndicesCollection getIndicesHistory() const;
  Collection<Point> getCoefficientsHistory() const;
  void setSelectionHistory(Collection<Indices> & indicesHistory, Collection<Point> & coefficientsHistory);
  Graph drawSelectionHistory() const;

  /** Error history accessor */
  void setErrorHistory(const Point & errorHistory);
  Point getErrorHistory() const;
  Graph drawErrorHistory() const;

protected:

private:
  /** The input vector distribution */
  Distribution distribution_;

  /** The isoprobabilistic transformation maps the distribution into the orthogonal measure */
  Function transformation_;

  /** The inverse isoprobabilistic transformation */
  Function inverseTransformation_;

  /** The orthogonal basis */
  OrthogonalBasis orthogonalBasis_;

  /** The collection of indices */
  Indices I_;

  /** The collection of Alpha_k coefficients */
  Sample alpha_k_;

  /** The collection of vectors upon which we project the (composed) model */
  FunctionPersistentCollection Psi_k_;

  /** Composed meta model */
  Function composedMetaModel_;

  /** Selection history */
  PersistentCollection<Indices> indicesHistory_;
  PersistentCollection<Point> coefficientsHistory_;

  /** Error history */
  Point errorHistory_;

  /** Is regression? */
  Bool isLeastSquares_ = true;

  /** Is model selection? */
  Bool involvesModelSelection_ = false;

} ; /* class FunctionalChaosResult */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_FUNCTIONALCHAOSRESULT_HXX */
