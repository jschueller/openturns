//                                               -*- C++ -*-
/**
 *  @brief RandomGenerator implements methods to control the random generator
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
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
#include <random>

#include "pcg_random.hpp"

#include "openturns/RandomGenerator.hxx"
#include "openturns/RandomGeneratorState.hxx"
#include "openturns/ResourceMap.hxx"


BEGIN_NAMESPACE_OPENTURNS

static ::pcg64 & GetGenerator()
{
  // Function-local static defers construction (and the ResourceMap read)
  // until first use, avoiding the static-initialization-order fiasco that
  // would occur if Generator_ were initialized at dynamic-init time before
  // ResourceMap is ready.
  static ::pcg64 generator = ::pcg64(ResourceMap::GetAsUnsignedInteger("RandomGenerator-InitialSeed"));
  return generator;
}

Bool RandomGenerator::IsInitialized_ = false;

/* DefaultConstructor */
RandomGenerator::RandomGenerator()
{
  // Nothing to do
}

/* Seed accessor */
void RandomGenerator::SetSeed(const UnsignedInteger seed)
{
  GetGenerator().seed(seed);
  IsInitialized_ = true;
}

/* State accessor */
void RandomGenerator::SetState(const RandomGeneratorState & state)
{
  std::stringstream oss;
  oss << state.getBuffer();
  oss >> GetGenerator();
  IsInitialized_ = true;
  return;
}

/* Seed accessor */
RandomGeneratorState RandomGenerator::GetState()
{
  std::stringstream oss;
  oss << GetGenerator();
  return RandomGeneratorState(oss.str());
}

void RandomGenerator::Reset()
{
  SetSeed(ResourceMap::GetAsUnsignedInteger("RandomGenerator-InitialSeed"));
}

void RandomGenerator::Initialize()
{
  if (!IsInitialized_)
    Reset();
}

/* Generate a pseudo-random number uniformly distributed over ]0, 1[ */
Scalar RandomGenerator::Generate()
{
  Initialize();
  std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
  return uniform(GetGenerator());
}

/* Generate a pseudo-random integer uniformly distributed over [[0,...,n-1]] */
UnsignedInteger RandomGenerator::IntegerGenerate(const UnsignedInteger n)
{
  if (n == 0) throw InvalidArgumentException(HERE) << "Error: n must be positive, got n=" << n;
  Initialize();
  std::uniform_int_distribution<UnsignedInteger> uniform(0, n - 1);
  return uniform(GetGenerator());
}

/* Generate a pseudo-random vector of numbers uniformly distributed over ]0, 1[ */
Point RandomGenerator::Generate(const UnsignedInteger size)
{
  Point result(size);
  Initialize();
  std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
  for (UnsignedInteger i = 0; i < size; ++ i)
    result[i] = uniform(GetGenerator());
  return result;
}

/* Generate a pseudo-random vector of numbers uniformly distributed over [[0,...,n-1]] */
RandomGenerator::UnsignedIntegerCollection RandomGenerator::IntegerGenerate(const UnsignedInteger size, const UnsignedInteger n)
{
  if (n == 0) throw InvalidArgumentException(HERE) << "Error: n must be positive, got n=" << n;
  UnsignedIntegerCollection result(size);
  Initialize();
  std::uniform_int_distribution<UnsignedInteger> uniform(0, n - 1);
  for (UnsignedInteger i = 0; i < size; ++ i)
    result[i] = uniform(GetGenerator());
  return result;
}

END_NAMESPACE_OPENTURNS
