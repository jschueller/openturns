//                                               -*- C++ -*-
/**
 *  @brief SlicedInverseRegression
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
#include "openturns/SlicedInverseRegression.hxx"
#include "openturns/PersistentObjectFactory.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(SlicedInverseRegression);

static Factory<SlicedInverseRegression> Factory_SlicedInverseRegression;


/* Default constructor */
SlicedInverseRegression::SlicedInverseRegression()
  : PersistentObject()
{
  // Nothing to do
}

SlicedInverseRegression::SlicedInverseRegression(const Sample & inputSample,
                                                 const Sample & outputSample)
  : PersistentObject()
  , inputSample_(inputSample)
  , outputSample_(outputSample)
  , modesNumber_(inputSample.getDimension())
{
  if (outputSample.getDimension() != 1)
    throw InvalidArgumentException(HERE) << "Supervision variable must be of dimension 1";
}

/* Virtual constructor method */
SlicedInverseRegression * SlicedInverseRegression::clone() const
{
  return new SlicedInverseRegression(*this);
}

void SlicedInverseRegression::run()
{ 
  const UnsignedInteger size = inputSample_.getSize();
  const UnsignedInteger inputDimension = inputSample_.getDimension();
  const Indices supervisionIndices = outputSample_.argsort();
  Collection<Indices> list_chunk;
  UnsignedInteger offset = 0;
  Indices chunk_population;
  for (UnsignedInteger i = 0; i < sliceNumber_; ++ i)
  {
    const UnsignedInteger localSize = std::min(size / sliceNumber_, size - offset);
    Indices chunk(localSize);
    std::copy(supervisionIndices.begin() + offset, supervisionIndices.begin() + offset + localSize, chunk.begin());
    list_chunk.add(chunk);
    chunk_population.add(localSize);
    offset += localSize;
  }
  const Point mean(inputSample_.computeMean());
  Sample X_centered(inputSample_);
  X_centered -= mean;
  Matrix Q(size, inputDimension);
  for (UnsignedInteger i = 0; i < size; ++ i)
    for(UnsignedInteger j = 0; j < inputDimension; ++ j)
      Q(i, j) = inputSample_(i, j);
  Matrix R;
  Q.computeQRInPlace(R);
  Sample Z(size, inputDimension);
  for (UnsignedInteger i = 0; i < size; ++ i)
    for(UnsignedInteger j = 0; j < inputDimension; ++ j)
      Z(i, j) = Q(i, j) * std::sqrt(1.0 * size);
  Matrix zMeans(sliceNumber_, inputDimension);
  for(UnsignedInteger j = 0; j < sliceNumber_; ++ j)
  {
    const Point zMean(Z.select(list_chunk[j]).computeMean() * std::sqrt(1.0 * chunk_population[j]));
    for(UnsignedInteger i = 0; i < inputDimension; ++ i)
      zMeans(j, i) = zMean[i];
  }
  SymmetricMatrix M((zMeans.transpose() * zMeans / size).getImplementation());
  SquareMatrix eigenVectors;
  const Point eigenValues = M.computeEVInPlace(eigenVectors);
  SquareMatrix eigenVectorsRev(inputDimension);
  Point eigenValuesRev(inputDimension);
  for(UnsignedInteger j = 0; j < inputDimension; ++ j)
  {
    eigenValuesRev[j] = eigenValues[inputDimension - 1 - j];
    for(UnsignedInteger i = 0; i < inputDimension; ++ i)
      eigenVectorsRev(i, j) = eigenVectors(i, inputDimension - 1 - j);
  }
  TriangularMatrix R2((R * std::sqrt(1.0 * size)).getImplementation(), false);
  Matrix directions = R2.solveLinearSystemInPlace(eigenVectorsRev);

  // prune directions and associated eigenvalues
  directions.getImplementation()->resize(inputDimension, modesNumber_);
  eigenValuesRev.resize(modesNumber_);

  // normalize directions
  for(UnsignedInteger j = 0; j < modesNumber_; ++ j)
  {
    Scalar normJ = 0.0;
    for(UnsignedInteger i = 0; i < inputDimension; ++ i)
      normJ += directions(i, j) * directions(i, j);
    normJ = 1.0 / std::sqrt(normJ);
    for(UnsignedInteger i = 0; i < inputDimension; ++ i)
      directions(i, j) *= normJ;
  }
  result_ = SlicedInverseRegressionResult(directions, mean, eigenValuesRev);
}

SlicedInverseRegressionResult SlicedInverseRegression::getResult() const
{
  return result_;
}

/* String converter */
String SlicedInverseRegression::__repr__() const
{
  OSS oss;
  oss << "class=" << SlicedInverseRegression::GetClassName();
  return oss;
}


/* Slice number accessor */
UnsignedInteger SlicedInverseRegression::getSliceNumber() const
{
  return sliceNumber_;
}

void SlicedInverseRegression::setSliceNumber(const UnsignedInteger sliceNumber)
{
  sliceNumber_ = sliceNumber;
}

/* Modes number accessor */
void SlicedInverseRegression::setModesNumber(const UnsignedInteger modesNumber)
{
  if (modesNumber_ > inputSample_.getDimension())
    throw InvalidArgumentException(HERE) << "Cannot use more than " <<inputSample_.getDimension() << " modes";
  throw NotYetImplementedException(HERE) << "Cannot do that yet";
  modesNumber_ = modesNumber;
}

UnsignedInteger SlicedInverseRegression::getModesNumber() const
{
  return modesNumber_;
}

/* Method save() stores the object through the StorageManager */
void SlicedInverseRegression::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute( "inputSample_", inputSample_ );
  adv.saveAttribute( "outputSample_", outputSample_ );
  adv.saveAttribute( "sliceNumber_", sliceNumber_ );
  adv.saveAttribute( "modesNumber", modesNumber_ );
  adv.saveAttribute( "result_", result_ );
}

/* Method load() reloads the object from the StorageManager */
void SlicedInverseRegression::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute( "inputSample_", inputSample_ );
  adv.loadAttribute( "outputSample_", outputSample_ );
  adv.loadAttribute( "sliceNumber_", sliceNumber_ );
  adv.loadAttribute( "modesNumber", modesNumber_ );
  adv.loadAttribute( "result_", result_ );
}

END_NAMESPACE_OPENTURNS
