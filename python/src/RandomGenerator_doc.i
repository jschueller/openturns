%feature("docstring") OT::RandomGenerator
"Uniform random generator.

Refer to :any:`uniform_random_generator`.

The random generator of uniform(0,1) samples is based on the PCG algorithm."

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::GetState
"Get the state of the random generator.

Returns
-------
particularState : :class:`~openturns.RandomGeneratorState`
    State of the random generator."

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::SetState
"Set the state of the random generator.

Parameters
----------
particularState : :class:`~openturns.RandomGeneratorState`
    State of the random generator.

Notes
-----
This method fixes the entire state of the random generator algorithm
thanks the specification of the entire state *particularState* usually
previously obtained thanks to the :meth:`GetState` method."

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::Generate
"Generate a pseudo-random vector.

Parameters
----------
size : positive int
    Number of realizations required. When not given, by default it is equal to
    1.

Returns
-------
numbers : float, :class:`~openturns.Point`
    Pseudo-random numbers uniformly distributed over :math:`[0, 1[`.

Examples
--------
>>> import openturns as ot
>>> ot.RandomGenerator.SetSeed(0)
>>> print('%.6f' % ot.RandomGenerator.Generate())
0.629877
>>> print(ot.RandomGenerator.Generate(2))
[0.882805,0.135276]"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::IntegerGenerate
"Generate a pseudo-random integer.

:Available usages:
    IntegerGenerate(*n*)

    IntegerGenerate(*size, n*)

Parameters
----------
n : positive int
    Upper bound of the interval where the pseudo-random integers are.
size : positive int
    Number of integers to generate. When not given, by default it is equal to
    1.

Returns
-------
integer : int, :class:`~openturns.UnsignedIntegerCollection`
    Pseudo-random integers uniformly distributed over :math:`[0,...,n-1]`.

Examples
--------
>>> import openturns as ot
>>> ot.RandomGenerator.SetSeed(0)
>>> print(ot.RandomGenerator.IntegerGenerate(30))
24
>>> print(ot.RandomGenerator.IntegerGenerate(5, 30))
[26,21,21,22,26]"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::SetSeed
R"RAW(Set the seed of the random generator.

Parameters
----------
n : int :math:`\in [0, 2^{32}-1]`

Notes
-----
This method fixes a particular state of the random generator algorithm
thanks to the seed :math:`n`. The seed of the random generator is automatically
initialized to 0 when a session is launched.)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::Reset
"Reset to the default state of the random generator.

Uses the default initial seed set in the :class:`~openturns.ResourceMap` entry `RandomGenerator-InitialSeed`.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::SetStream
"Select an independent stream of the global random generator.

Parameters
----------
stream : int
    Stream index, selects one of the 2^127 independent PCG streams.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::Advance
"Advance the global random generator without generating.

Parameters
----------
delta : int
    Number of steps to skip.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::Spawn
"Spawn independent child generators from the global generator.

Parameters
----------
n : positive int
    Number of child generators.

Returns
-------
children : collection of :class:`~openturns.RandomGenerator`
    Independent generators derived from the global generator state,
    suitable for parallel simulations.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::setSeed
"Set the seed of the generator.

Parameters
----------
seed : int
    Seed value.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::setStream
"Select an independent stream of the generator.

Parameters
----------
stream : int
    Stream index, selects one of the 2^127 independent PCG streams.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::setState
"Set the state of the generator.

Parameters
----------
state : :class:`~openturns.RandomGeneratorState`
    State of the generator.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::getState
"Get the state of the generator.

Returns
-------
state : :class:`~openturns.RandomGeneratorState`
    State of the generator.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::advance
"Advance the generator without generating.

Parameters
----------
delta : int
    Number of steps to skip.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::spawn
"Spawn independent child generators from this generator.

Parameters
----------
n : positive int
    Number of child generators.

Returns
-------
children : collection of :class:`~openturns.RandomGenerator`
    Independent generators derived from this generator state,
    suitable for parallel simulations.

Examples
--------
>>> import openturns as ot
>>> rng = ot.RandomGenerator(0)
>>> children = rng.spawn(2)
>>> len(children)
2
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::generate
"Generate pseudo-random numbers uniformly distributed over :math:`[0, 1[`.

Parameters
----------
size : positive int
    Number of realizations required. When not given, by default it is equal to
    1.

Returns
-------
numbers : float, :class:`~openturns.Point`
    Pseudo-random numbers.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::RandomGenerator::integerGenerate
"Generate pseudo-random integers uniformly distributed over :math:`[0,...,n-1]`.

Parameters
----------
n : positive int
    Upper bound of the interval.
size : positive int
    Number of integers to generate. When not given, by default it is equal to
    1.

Returns
-------
integers : int, :class:`~openturns.UnsignedIntegerCollection`
    Pseudo-random integers.
"
