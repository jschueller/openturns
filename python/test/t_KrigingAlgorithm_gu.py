#! /usr/bin/env python

from __future__ import print_function
from openturns import *
from openturns.testing import *
import openturns as ot

TESTPREAMBLE()

# Test 1

for prior in [ot.GeneralLinearModelAlgorithm.NONE, ot.GeneralLinearModelAlgorithm.JOINTLYROBUST, ot.GeneralLinearModelAlgorithm.REFERENCE]:
    for parametrization in [ot.CovarianceModelImplementation.STANDARD, ot.CovarianceModelImplementation.INVERSE, ot.CovarianceModelImplementation.LOGINVERSE]:

        ot.RandomGenerator.SetSeed(0)

        print('prior=', prior)
        print('parametrization=', parametrization)
        sampleSize = 6
        dimension = 1

        f = SymbolicFunction(['x0'], ['x0 * sin(x0)'])

        X = Sample(sampleSize, dimension)
        X2 = Sample(sampleSize, dimension)
        for i in range(sampleSize):
            X[i, 0] = 3.0 + i
            X2[i, 0] = 2.5 + i
        X[0, 0] = 1.0
        X[1, 0] = 3.0
        X2[0, 0] = 2.0
        X2[1, 0] = 4.0
        Y = f(X)
        Y2 = f(X2)

        # create algorithm
        basis = ConstantBasisFactory(dimension).build()
        covarianceModel = SquaredExponential([1e-05], [4.11749])
        covarianceModel.setScaleParametrization(parametrization)

        algo = KrigingAlgorithm(X, Y, covarianceModel, basis)
        algo.setScalePrior(prior)
        algo.run()

        # perform an evaluation
        result = algo.getResult()
        print('parameters=', result.getCovarianceModel().getParameter())
        #print("X=", X)
        #print("f(X)=", Y)

        assert_almost_equal(result.getMetaModel()(X), Y)
        assert_almost_equal(result.getResiduals(), [1.32804e-07], 1e-3, 1e-3)
        assert_almost_equal(result.getRelativeErrors(), [5.20873e-21])

        # Kriging variance is 0 on learning points
        var = result.getConditionalCovariance(X)

        # assert_almost_equal could not be applied to matrices
        # application to Point
        covariancePoint = Point(var.getImplementation())
        theoricalVariance = Point(sampleSize * sampleSize)
        assert_almost_equal(covariancePoint, theoricalVariance, 8.95e-7, 8.95e-7)
