"""
Use the LOLA-Voronoi sequential design of experiment
====================================================
"""

# %%
# The LOLA-Voronoi sequential experiment helps to generate an optimized design allowing
# better approximations of functions taking into account empty regions and gradient values.
# It can be relevant to build a design of experiment for a metamodel: we will compare the
# LOLA-Voronoi design against a Sobol' design as learning points for a kriging metamodel.

# %%
import openturns as ot
import openturns.experimental as otexp
import openturns.viewer as otv
import math

# %%
# Lets use the Rastrigin standard optimization test-case


def rastriginPy(X):
    A = 10.0
    delta = [x**2 - A * math.cos(3.0 * math.pi * x) for x in X]
    y = A + sum(delta)
    return [y]


dim = 2
f1 = ot.PythonFunction(dim, 1, rastriginPy)
print(f1([1.0, 1.0]))
distribution = ot.JointDistribution([ot.Uniform(-0.5, 0.5)] * 2)

# %%
# Plot the gradient


def pyGradNorm(X):
    grad = f1.gradient(X)
    grad.squareElements()
    s = grad.computeSumElements()
    return [s**0.5]


gradNorm = ot.PythonFunction(f1.getInputDimension(), 1, pyGradNorm)

# %%
# Plot the function
graph = f1.draw(distribution.getRange().getLowerBound(), distribution.getRange().getUpperBound())
otv.View(graph)

# %%
# Plot the gradient
graph = gradNorm.draw(distribution.getRange().getLowerBound(), distribution.getRange().getUpperBound())
graph.setTitle("Gradient norm as a function of (x0, y0)")
otv.View(graph)

# %%
# Lets define an initial design of experiments
x0 = ot.LowDiscrepancyExperiment(ot.HaltonSequence(), distribution, 50).generate()
y0 = f1(x0)

# %%
# Plot the initial input sample
graph = ot.Graph("LOLA-Voronoi", "x1", "x2", True)
cloud1 = ot.Cloud(x0)
cloud1.setPointStyle("fcircle")
cloud1.setColor("blue")
graph.add(cloud1)
otv.View(graph)

# %%
# Instantiate the algorithm from the initial DOE and the distribution
algo = otexp.LOLAVoronoi(x0, y0, distribution)

# %%
# Iteratively generate new samples
for i in range(10):
    x = algo.generate(5)
    y = f1(x)
    algo.update(x, y)

    # Update the graph
    cloud2 = ot.Cloud(x)
    cloud2.setPointStyle("fcircle")
    cloud2.setColor("red")
    graph.add(cloud2)
    otv.View(graph)


# %%
# Lets compare metamodels from LOLA samples and from few other designs
xLola, yLola = algo.getX(), algo.getY()
learnSize = xLola.getSize()


def runMetaModel(x, y, tag):
    cov = ot.ExponentialModel([1.1] * 2, [10.0])
    kriging = ot.KrigingAlgorithm(x, y, cov)
    kriging.setOptimizeParameters(False)
    kriging.run()
    result = kriging.getResult()
    metamodel = result.getMetaModel()
    if kriging.getOptimizeParameters():
        print(f"cov={result.getCovarianceModel()}")
    yPred = metamodel(xRef)
    validation = ot.MetaModelValidation(yRef, yPred)
    mse = validation.computeMeanSquaredError()
    print(f"{tag} mse={mse} r2={validation.computeR2Score()}")


# %%
# Generate a large validation sample by Monte Carlo
nRef = int(1e6)
xRef = distribution.getSample(nRef)
yRef = f1(xRef)

# %%
# Build a metamodel from Sobol' samples
xSobol = ot.LowDiscrepancyExperiment(ot.SobolSequence(), distribution, learnSize).generate()
ySobol = f1(xSobol)
runMetaModel(xSobol, ySobol, "Sobol")

# %%
# Build a metamodel on the LOLA global samples
# We observe that the metamodel on the LOLA design is better than the Sobol' experiment
runMetaModel(xLola, yLola, "LOLA")

# %%
otv.View.ShowAll()
