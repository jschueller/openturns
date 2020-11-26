#! /usr/bin/env python

from __future__ import print_function
import openturns as ot
from openturns.usecases.ishigami_function import IshigamiModel

ot.TESTPREAMBLE()

ishigami = IshigamiModel()
X = ishigami.distributionX
f = ishigami.model

size = 10
x = X.getSample(size)
y = f(x)

print(x, y)

experiment = ot.DualBootstrapExperiment(x, y)
y1 = ot.Sample()
x1, y1 = experiment.generate()
print(x1, y1)
assert len(x1) == size, "unmatching size"
assert len(x1) == len(y1), "unmatching sizes"

# check that boostrap selection is consistent across the two samples
for i in range(size):
    i_x = x.find(x1[i])
    i_y = y.find(y1[i])
assert i_x == i_y, "unmatching index"
