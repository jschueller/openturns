#! /usr/bin/env python

import openturns as ot
import openturns.experimental as otexp
import openturns.testing as ott

dim = 2
f = ot.SymbolicFunction(["x1", "x2"], ["1+100*(x2-x1^2)^2+(1-x1)^2"])
x0 = [1e-3] * dim
algo = otexp.Adam()
algo.setAlpha(1e-2)
for minimization in [True, False]:
    problem = ot.OptimizationProblem(f)
    problem.setMinimization(minimization)
    algo.setProblem(problem)
    algo.setStartingPoint(x0)
    try:
        algo.run()
    except Exception as e:
        print("-- ", e)
        continue
    result = algo.getResult()
    x = result.getOptimalPoint()
    y = result.getOptimalValue()[0]
    print("x^=", x, "y^=", y)
    # Adam does not read the actual minimum at [1, 1]
    if minimization:
        assert y < 1.1
    else:
        assert y > 1e2
