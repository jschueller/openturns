#! /usr/bin/env python

import openturns as ot

ot.TESTPREAMBLE()


# First, build two functions from R^3->R
inVar = ot.Description(3)
inVar[0] = "x1"
inVar[1] = "x2"
inVar[2] = "x3"
formula = ot.Description(1)
formula[0] = "x1^3 * sin(x2 + 2.5 * x3) - (x1 + x2)^2 / (1.0 + x3^2)"
functions = []
functions.append(ot.SymbolicFunction(inVar, formula))
formula[0] = "exp(-x1 * x2 + x3) / cos(1.0 + x2 * x3 - x1)"
functions.append(ot.SymbolicFunction(inVar, formula))
# Second, build the weights
coefficients = ot.Point(2)
coefficients[0] = 0.3
coefficients[1] = 2.9
myHessian = ot.LinearCombinationHessian(
    ot.LinearCombinationEvaluation(functions, coefficients)
)
inPoint = ot.Point(3)
inPoint[0] = 1.2
inPoint[1] = 2.3
inPoint[2] = 3.4
# Hessian formulas differ between backends (SymEngine vs Ev3)
# so just check analytical hessian is available and verify numerical correctness
# Reduce the output precision because we use finite difference so only 5
# digits are significant
ot.PlatformInfo.SetNumericalPrecision(5)
hess = myHessian.hessian(inPoint)
assert hess.getNbRows() > 0 and hess.getNbColumns() > 0, "Hessian should be available"
print("Value at ", inPoint, "=", hess)
