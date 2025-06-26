#! /usr/bin/env python

import openturns as ot

ot.TESTPREAMBLE()

# from highs c++ example: https://github.com/ERGO-Code/HiGHS/blob/master/examples/call_highs_from_cpp.cpp
f = ot.SymbolicFunction(["x0", "x1"], ["x0+x1+3"])
bounds = ot.Interval([0.0, 1.0], [4.0, 1e30])
ineq = ot.SymbolicFunction(["x0", "x1"], ["7-x1", "x0+2*x1-5", "15-x0-2*x1", "3*x0+2*x1-6"])
problem = ot.OptimizationProblem(f)
problem.setBounds(bounds)
problem.setInequalityConstraint(ineq)
problem.setVariablesType(
    [ot.OptimizationProblemImplementation.CONTINUOUS] * 4
    + [ot.OptimizationProblemImplementation.INTEGER] * 2
)
