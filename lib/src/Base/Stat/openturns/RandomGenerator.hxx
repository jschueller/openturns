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
#ifndef OPENTURNS_RANDOMGENERATOR_HXX
#define OPENTURNS_RANDOMGENERATOR_HXX

#include <memory>

#include "openturns/OTprivate.hxx"
#include "openturns/PersistentObject.hxx"
#include "openturns/Collection.hxx"
#include "openturns/Point.hxx"
#include "openturns/RandomGeneratorState.hxx"


BEGIN_NAMESPACE_OPENTURNS

/**
 * @class RandomGenerator
 *
 * RandomGenerator implements methods to control the random generator.
 * It can be used either through its static methods, which share a global
 * generator for backward compatibility, or as an instantiable object holding
 * its own independent random stream (see spawn()).
 */
class OT_API RandomGenerator
  : public PersistentObject
{
  CLASSNAME
public:

  typedef Collection<UnsignedInteger> UnsignedIntegerCollection;
  typedef Collection<RandomGenerator> RandomGeneratorCollection;

  /** Default constructor: seeded from ResourceMap */
  RandomGenerator();

  /** Constructor from seed */
  explicit RandomGenerator(const UnsignedInteger seed);

  /** Constructor from seed and stream */
  RandomGenerator(const UnsignedInteger seed,
                  const UnsignedInteger stream);

  /** Copy constructor (deep copy of the generator state) */
  RandomGenerator(const RandomGenerator & other);

  /** Assignment operator (deep copy of the generator state) */
  RandomGenerator & operator=(const RandomGenerator & other);

  /** Virtual constructor */
  RandomGenerator * clone() const override;

  /** String converter */
  String __repr__() const override;
  String __str__(const String & offset = "") const override;

  /** Method save() stores the object through the StorageManager */
  void save(Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(Advocate & adv) override;

  /** Seed accessor */
  void setSeed(const UnsignedInteger seed) const;

  /** Stream accessor: select one of the 2^127 independent PCG streams */
  void setStream(const UnsignedInteger stream) const;

  /** State accessor */
  void setState(const RandomGeneratorState & state) const;
  RandomGeneratorState getState() const;

  /** Advance the generator by delta steps without generating */
  void advance(const UnsignedInteger delta) const;

  /** Spawn n independent child generators derived from this generator */
  RandomGeneratorCollection spawn(const UnsignedInteger n) const;

  /** Generate a pseudo-random number uniformly distributed over [0, 1[ */
  Scalar generate() const;
  /** Generate a pseudo-random integer uniformly distributed over [[0,...,n-1]] */
  UnsignedInteger integerGenerate(const UnsignedInteger n) const;

  /** Generate a pseudo-random vector of numbers uniformly distributed over [0, 1[ */
  Point generate(const UnsignedInteger size) const;
  /** Generate a pseudo-random vector of integers uniformly distributed over [[0,...,n-1]] */
  UnsignedIntegerCollection integerGenerate(const UnsignedInteger size,
      const UnsignedInteger n) const;

  /** Seed accessor */
  static void SetSeed(const UnsignedInteger seed);

  /** Stream accessor on the global generator */
  static void SetStream(const UnsignedInteger stream);

  /** Advance the global generator by delta steps without generating */
  static void Advance(const UnsignedInteger delta);

  /** Spawn n independent child generators derived from the global generator */
  static RandomGeneratorCollection Spawn(const UnsignedInteger n);

  /** State accessor */
  static void SetState(const RandomGeneratorState & state);
  static RandomGeneratorState GetState();

  /** Generate a pseudo-random number uniformly distributed over [0, 1[ */
  static Scalar Generate();
  /** Generate a pseudo-random integer uniformly distributed over [[0,...,n-1]] */
  static UnsignedInteger IntegerGenerate(const UnsignedInteger n);

  /** Generate a pseudo-random vector of numbers uniformly distributed over [0, 1[ */
  static Point Generate(const UnsignedInteger size);
  /** Generate a pseudo-random vector of integers uniformly distributed over [[0,...,n-1]] */
  static UnsignedIntegerCollection IntegerGenerate(const UnsignedInteger size, const UnsignedInteger n);

  /* Reset with default seed */
  static void Reset();

#ifndef SWIG
  /** Accessor to the global generator (shared stream for backward compatibility) */
  static const RandomGenerator & GetGlobalGenerator();
#endif

private:

  struct GeneratorImpl;
  std::shared_ptr<GeneratorImpl> p_impl_;

  static void Initialize();

  static Bool IsInitialized_;

}; /* class RandomGenerator */


END_NAMESPACE_OPENTURNS

#endif /* OPENTURNS_RANDOMGENERATOR_HXX */
