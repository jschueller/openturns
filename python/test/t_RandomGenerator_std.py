#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import os

ot.TESTPREAMBLE()

# Test the initialization method
ot.RandomGenerator.SetSeed(0)

# Test the state accessor for reading
initialState = ot.RandomGenerator.GetState()
print("initialState=", repr(initialState))
# Test the generation of one realization
realization = ot.RandomGenerator.Generate()
print("one realization=%.6f" % realization)
state = ot.RandomGenerator.GetState()
print("state=", repr(state))

# Test the generation of several realizations
size = 10000
several = ot.RandomGenerator.Generate(size)
print(size, " realizations, first=%.6f" % several[0], ", last=%.6f" % several[size - 1])

# Test the state accessor for writing
ot.RandomGenerator.SetState(initialState)
print("restoring initial state, one realization=%.6f" % ot.RandomGenerator.Generate())
# Test the uniformity of the random generator
size = 100000
manyPoints = ot.RandomGenerator.Generate(size)
slice = 10
frequencies = ot.Point(slice, 0.0)
mean = 0.0
for i in range(size):
    frequencies[int(slice * manyPoints[i])] = (
        frequencies[int(slice * manyPoints[i])] + 1
    )
    mean = mean + manyPoints[i]

mean = mean / size
print("mean of", size, "realizations= %.12f" % mean)
for i in range(slice):
    frequencies[i] = frequencies[i] / size
    print(
        "frequency in slice [",
        (1.0 * i) / slice,
        ",",
        (i + 1.0) / slice,
        "]=",
        frequencies[i],
    )

print(
    "One integer generation in [[0, ",
    slice - 1,
    "]]=",
    ot.RandomGenerator.IntegerGenerate(slice),
)
print(
    "Another integer generation in [[0, ",
    slice - 1,
    "]]=",
    ot.RandomGenerator.IntegerGenerate(slice),
)
manyUnsignedInteger = ot.RandomGenerator.IntegerGenerate(size, slice)
frequencies = ot.Point(slice, 0.0)
mean = 0.0
for i in range(size):
    frequencies[manyUnsignedInteger[i]] = frequencies[manyUnsignedInteger[i]] + 1
    mean = mean + manyUnsignedInteger[i]

mean = mean / size
print("mean of ", size, " integer realizations=", mean)
for i in range(slice):
    frequencies[i] = frequencies[i] / size
    print("frequency for value ", i, "=", frequencies[i])

# IntegerGenerate(0) should throw (div by zero guard)
with ott.assert_raises(TypeError):
    ot.RandomGenerator.IntegerGenerate(0)

# Independent per-distribution streams: distributions with spawned
# generators must be reproducible and isolated from the global stream
ot.RandomGenerator.SetSeed(0)
children = ot.RandomGenerator.Spawn(2)

uniform = ot.Uniform(0.0, 1.0)
uniform.setRandomGenerator(children[0])
assert uniform.hasRandomGenerator()

normal = ot.Normal(0.0, 1.0)
normal.setRandomGenerator(children[1])
assert normal.hasRandomGenerator()

x1 = uniform.getSample(5)
y1 = normal.getSample(5)

# Reproducibility: same spawned generators give the same samples
ot.RandomGenerator.SetSeed(0)
children = ot.RandomGenerator.Spawn(2)
uniform2 = ot.Uniform(0.0, 1.0)
uniform2.setRandomGenerator(children[0])
normal2 = ot.Normal(0.0, 1.0)
normal2.setRandomGenerator(children[1])
ott.assert_almost_equal(x1, uniform2.getSample(5))
ott.assert_almost_equal(y1, normal2.getSample(5))

# Isolation: local draws do not advance the global stream
ot.RandomGenerator.SetSeed(0)
before = ot.RandomGenerator.Generate()
ot.RandomGenerator.SetSeed(0)
uniform.getSample(100)
after = ot.RandomGenerator.Generate()
ott.assert_almost_equal(before, after)

# Clearing falls back to the global generator
uniform.clearRandomGenerator()
assert not uniform.hasRandomGenerator()

# Persistence of the attached generator state
uniform.setRandomGenerator(ot.RandomGenerator(42))
buffer = uniform.getRandomGenerator().getState().getBuffer()
fileName = "study_local_rng.xml"
study = ot.Study(fileName)
study.add("uniform", uniform)
study.save()
loaded = ot.Uniform()
study2 = ot.Study(fileName)
study2.load()
study2.fillObject("uniform", loaded)
os.remove(fileName)
assert loaded.hasRandomGenerator()
assert loaded.getRandomGenerator().getState().getBuffer() == buffer

