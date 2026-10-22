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
#include <mutex>
#include <random>

#include "pcg_random.hpp"

#include "openturns/RandomGenerator.hxx"
#include "openturns/RandomGeneratorState.hxx"
#include "openturns/ResourceMap.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/OSS.hxx"

BEGIN_NAMESPACE_OPENTURNS

CLASSNAMEINIT(RandomGenerator)

static const Factory<RandomGenerator> Factory_RandomGenerator;

struct RandomGenerator::GeneratorImpl
{
  ::pcg64 engine_;
  mutable std::mutex mutex_;

  explicit GeneratorImpl(const UnsignedInteger seed)
    : engine_(seed)
  {
    // Nothing to do
  }

  GeneratorImpl(const UnsignedInteger seed,
                const UnsignedInteger stream)
    : engine_(seed, stream)
  {
    // Nothing to do
  }

  GeneratorImpl(const GeneratorImpl & other)
  {
    std::lock_guard<std::mutex> lock(other.mutex_);
    engine_ = other.engine_;
  }
};

Bool RandomGenerator::IsInitialized_ = false;

/* Default constructor: seeded from ResourceMap */
RandomGenerator::RandomGenerator()
  : p_impl_(new GeneratorImpl(ResourceMap::GetAsUnsignedInteger("RandomGenerator-InitialSeed")))
{
  // Nothing to do
}

/* Constructor from seed */
RandomGenerator::RandomGenerator(const UnsignedInteger seed)
  : p_impl_(new GeneratorImpl(seed))
{
  // Nothing to do
}

/* Constructor from seed and stream */
RandomGenerator::RandomGenerator(const UnsignedInteger seed,
                                 const UnsignedInteger stream)
  : p_impl_(new GeneratorImpl(seed, stream))
{
  // Nothing to do
}

/* Copy constructor (deep copy of the generator state) */
RandomGenerator::RandomGenerator(const RandomGenerator & other)
  : PersistentObject(other)
  , p_impl_(new GeneratorImpl(*other.p_impl_))
{
  // Nothing to do
}

/* Assignment operator (deep copy of the generator state) */
RandomGenerator & RandomGenerator::operator=(const RandomGenerator & other)
{
  if (this != &other)
  {
    PersistentObject::operator=(other);
    p_impl_.reset(new GeneratorImpl(*other.p_impl_));
  }
  return *this;
}

/* Virtual constructor */
RandomGenerator * RandomGenerator::clone() const
{
  return new RandomGenerator(*this);
}

/* String converter */
String RandomGenerator::__repr__() const
{
  return OSS(true) << "class=" << getClassName()
         << " state=" << getState().getBuffer();
}

String RandomGenerator::__str__(const String & ) const
{
  return OSS(false) << getClassName() << "(" << getState().getBuffer() << ")";
}

/* Method save() stores the object through the StorageManager */
void RandomGenerator::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute("state_", getState().getBuffer());
}

/* Method load() reloads the object from the StorageManager */
void RandomGenerator::load(Advocate & adv)
{
  PersistentObject::load(adv);
  String buffer;
  adv.loadAttribute("state_", buffer);
  setState(RandomGeneratorState(buffer));
}

/* Accessor to the global generator (shared stream for backward compatibility) */
const RandomGenerator & RandomGenerator::GetGlobalGenerator()
{
  // Function-local static defers construction (and the ResourceMap read)
  // until first use, avoiding the static-initialization-order fiasco
  static const RandomGenerator generator(ResourceMap::GetAsUnsignedInteger("RandomGenerator-InitialSeed"));
  return generator;
}

/* Seed accessor */
void RandomGenerator::setSeed(const UnsignedInteger seed) const
{
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  p_impl_->engine_.seed(seed);
}

/* Stream accessor: select one of the 2^127 independent PCG streams */
void RandomGenerator::setStream(const UnsignedInteger stream) const
{
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  p_impl_->engine_.set_stream(stream);
}

/* State accessor */
void RandomGenerator::setState(const RandomGeneratorState & state) const
{
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  std::stringstream oss;
  oss << state.getBuffer();
  oss >> p_impl_->engine_;
}

/* State accessor */
RandomGeneratorState RandomGenerator::getState() const
{
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  std::stringstream oss;
  oss << p_impl_->engine_;
  return RandomGeneratorState(oss.str());
}

/* Advance the generator by delta steps without generating */
void RandomGenerator::advance(const UnsignedInteger delta) const
{
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  p_impl_->engine_.advance(delta);
}

/* Spawn n independent child generators derived from this generator */
RandomGenerator::RandomGeneratorCollection RandomGenerator::spawn(const UnsignedInteger n) const
{
  // Children derive both their initial state and stream from the parent,
  // so repeated calls produce distinct reproducible streams.
  // Child streams are drawn over the 127-bit PCG stream space, hence
  // collisions are negligible for practical purposes.
  RandomGeneratorCollection children(n);
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  for (UnsignedInteger i = 0; i < n; ++ i)
  {
    const pcg_extras::pcg128_t state = (static_cast< pcg_extras::pcg128_t >(p_impl_->engine_()) << 64) | p_impl_->engine_();
    const pcg_extras::pcg128_t stream = (static_cast< pcg_extras::pcg128_t >(p_impl_->engine_()) << 64) | p_impl_->engine_();
    RandomGenerator child;
    child.p_impl_->engine_.seed(state, stream);
    children[i] = child;
  }
  return children;
}

/* Generate a pseudo-random number uniformly distributed over [0, 1[ */
Scalar RandomGenerator::generate() const
{
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
  return uniform(p_impl_->engine_);
}

/* Generate a pseudo-random integer uniformly distributed over [[0,...,n-1]] */
UnsignedInteger RandomGenerator::integerGenerate(const UnsignedInteger n) const
{
  if (n == 0) throw InvalidArgumentException(HERE) << "Error: n must be positive, got n=" << n;
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  std::uniform_int_distribution<UnsignedInteger> uniform(0, n - 1);
  return uniform(p_impl_->engine_);
}

/* Generate a pseudo-random vector of numbers uniformly distributed over [0, 1[ */
Point RandomGenerator::generate(const UnsignedInteger size) const
{
  Point result(size);
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  std::uniform_real_distribution<Scalar> uniform(0.0, 1.0);
  for (UnsignedInteger i = 0; i < size; ++ i)
    result[i] = uniform(p_impl_->engine_);
  return result;
}

/* Generate a pseudo-random vector of integers uniformly distributed over [[0,...,n-1]] */
RandomGenerator::UnsignedIntegerCollection RandomGenerator::integerGenerate(const UnsignedInteger size,
    const UnsignedInteger n) const
{
  if (n == 0) throw InvalidArgumentException(HERE) << "Error: n must be positive, got n=" << n;
  UnsignedIntegerCollection result(size);
  std::lock_guard<std::mutex> lock(p_impl_->mutex_);
  std::uniform_int_distribution<UnsignedInteger> uniform(0, n - 1);
  for (UnsignedInteger i = 0; i < size; ++ i)
    result[i] = uniform(p_impl_->engine_);
  return result;
}

/* Seed accessor */
void RandomGenerator::SetSeed(const UnsignedInteger seed)
{
  GetGlobalGenerator().setSeed(seed);
  IsInitialized_ = true;
}

/* Stream accessor on the global generator */
void RandomGenerator::SetStream(const UnsignedInteger stream)
{
  GetGlobalGenerator().setStream(stream);
  IsInitialized_ = true;
}

/* Advance the global generator by delta steps without generating */
void RandomGenerator::Advance(const UnsignedInteger delta)
{
  Initialize();
  GetGlobalGenerator().advance(delta);
}

/* Spawn n independent child generators derived from the global generator */
RandomGenerator::RandomGeneratorCollection RandomGenerator::Spawn(const UnsignedInteger n)
{
  Initialize();
  return GetGlobalGenerator().spawn(n);
}

/* State accessor */
void RandomGenerator::SetState(const RandomGeneratorState & state)
{
  GetGlobalGenerator().setState(state);
  IsInitialized_ = true;
  return;
}

/* State accessor */
RandomGeneratorState RandomGenerator::GetState()
{
  return GetGlobalGenerator().getState();
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

/* Generate a pseudo-random number uniformly distributed over [0, 1[ */
Scalar RandomGenerator::Generate()
{
  Initialize();
  return GetGlobalGenerator().generate();
}

/* Generate a pseudo-random integer uniformly distributed over [[0,...,n-1]] */
UnsignedInteger RandomGenerator::IntegerGenerate(const UnsignedInteger n)
{
  Initialize();
  return GetGlobalGenerator().integerGenerate(n);
}

/* Generate a pseudo-random vector of numbers uniformly distributed over [0, 1[ */
Point RandomGenerator::Generate(const UnsignedInteger size)
{
  Initialize();
  return GetGlobalGenerator().generate(size);
}

/* Generate a pseudo-random vector of numbers uniformly distributed over [[0,...,n-1]] */
RandomGenerator::UnsignedIntegerCollection RandomGenerator::IntegerGenerate(const UnsignedInteger size, const UnsignedInteger n)
{
  Initialize();
  return GetGlobalGenerator().integerGenerate(size, n);
}

END_NAMESPACE_OPENTURNS
