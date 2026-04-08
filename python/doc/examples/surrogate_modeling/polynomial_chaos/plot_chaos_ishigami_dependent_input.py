"""
Create a FCE for dependent inputs: transformation vs domination
===============================================================
"""

# %%
#
# In this example, we create a functional chaos expansion for the
# :ref:`Ishigami function<use-case-ishigami>` when the input distribution
# has dependent marginals.
#
# refer to :ref:`functional_chaos` to learn more about functional chaos expansion.
#
# We provide one input sample and one output sample of the Ishigami function. We build two meta models:
#
# - Meta model 1: we use an isoprobabilistic transformation that maps the input distribution to another one with independent marginals.
# - Meta model 2: we use the domination method: the basis of the projection space is not orthonormal to the input distribution.
#


# %%
# Define the IshigamiModel model
# ------------------------------

# %%
from openturns.usecases import ishigami_function
import openturns as ot
import openturns.viewer as otv

# %%
# We load the Ishigami model. The `IshigamiModel` data class contains the input
# distribution :math:`\vect{X}=(X_1, X_2, X_3)` in `im.distribution` and the Ishigami
# function in `im.model`.
im = ishigami_function.IshigamiModel()
input_names = im.distribution.getDescription()

# %%
# We change the distribution in order to introduce some dependence. We use the normal copula.
copula = ot.BlockIndependentCopula([ot.ClaytonCopula(1.5), ot.IndependentCopula(1)])
input_dist = ot.JointDistribution([im.distribution.getMarginal(i) for i in range(im.dim)], copula)

# %%
# We can draw the function from a sample.
sampleSize = 1000
inputTrain = im.distribution.getSample(sampleSize)
outputTrain = im.model(inputTrain)

# %%
# Display relationships between the output and the inputs.
grid = ot.VisualTest.DrawPairsXY(inputTrain, outputTrain)
view = otv.View(grid, figure_kw={"figsize": (12.0, 4.0)})

# %%
# We draw the histogramm of the output values. The output disstribution has two modes.
graph = ot.HistogramFactory().build(outputTrain).drawPDF()
graph.setTitle("Ishigami outputs")
graph.setXTitle("y")
view = otv.View(graph)

# %%
# Meta model 1: Transformation method
# -----------------------------------
# A the input distribution has dependent marginals, we use an
# :ref:`isoprobabilistic_transformation` that maps the input distribution to
# another one with independent marginals.
#
# As the adaptive strategy is not specified, the following one is used:
#
# - the associated basis is built as the tensorization of the univariate polynomials family
#   orthonormal to the standard representative of the input marginals distribution. Then the basis is
#   orthonormal to the distribution denoted by :math:`\tilde{\mu}`;
# - the enumerate function is chosen according to
#   the `FunctionalChaosAlgorithm-QNorm` parameter of the :class:`~openturns.ResourceMap`:
#   if this parameter is equal to 1, then the :class:`~openturns.LinearEnumerateFunction` class is used,
#   otherwise, the :class:`~openturns.HyperbolicAnisotropicEnumerateFunction` class is used. The default
#   value of the key being 0.5, we use the Hyperbolic Anisotropic EnumerateFunction.
# - the first elements of the basis are used to build the approximation space. The number of elements
#   is computed from the total degree (using the enumerate function of the basis) specified as default value
#   in `FunctionalChaosAlgorithm-TotalDegree` of the :class:`~openturns.ResourceMap`.
#
# The projection stategy is not specified neither: we use the least-squares strategy  with no model selection
# if the key `FunctionalChaosAlgorithm-Sparse` of the :class:`~openturns.ResourceMap` is *False* and
# with model selection in the other case,  using the selection algorithm specified by the key
# `FunctionalChaosAlgorithm-FittingAlgorithm`. Considering the default values of the keys, we use a
# least-squares strategy with no model selection.
#
# The functional chaos algorithm uses an isoprobabilistic transformation that maps :math:`\mu` into
# :math:`\tilde{\mu}`. Note that the basis is not polynomials, due to the action of the
# isoprobabilistic transformation which is not affine.

# %%
# We create the functional chaos algorithm.
chaos_algo = ot.FunctionalChaosAlgorithm(inputTrain, outputTrain, input_dist)
chaos_algo.setUseDomination(False)
chaos_algo.run()

# %%
# We get the result and the resulting meta model.
chaos_result = chaos_algo.getResult()
metamodel = chaos_result.getMetaModel()

# %%
# In order to validate the metamodel, we generate a test sample.
n_valid = 1000
inputTest = input_dist.getSample(n_valid)
outputTest = im.model(inputTest)
metamodel_predictions = metamodel(inputTest)

# %%
# We draw the validation graph and we get the :math:`R^2` score: :math:`R^2 = 86.91\%`.
val = ot.MetaModelValidation(outputTest, metamodel_predictions)
r2Score = val.computeR2Score()[0]
print(f"r2Score with Transformation method = {r2Score:.6f}")
graph = val.drawValidation()
graph.setTitle(f"R2={r2Score * 100:.2f}%, use domination = false")
view = otv.View(graph)


# %%
# Meta model 2: Domination method
# -------------------------------
# Now, we want to use the domination method, which means that the basis created by the
# adaptive strategy is used to project the model. This basis is not orthonormal to :math:`\mu`.
#
# We use the :meth:`~openturns.FunctionalChaosAlgorithm.setUseDomination` method.
# We implement the same steps as before, until the validation graph.
chaos_algo.setUseDomination(True)
chaos_algo.run()

chaos_result = chaos_algo.getResult()
metamodel_dom = chaos_result.getMetaModel()
metamodel_dom_predictions = metamodel_dom(inputTest)
val = ot.MetaModelValidation(outputTest, metamodel_dom_predictions)
r2Score = val.computeR2Score()[0]
print(f"r2Score with Domination method = {r2Score:.6f}")
graph = val.drawValidation()
graph.setTitle(f"R2={r2Score * 100:.2f}%, use domination = true")
view = otv.View(graph)

# %%
# We can see that the meta model obtained with the domination metho is largely better than the previous one.
# %%
# Display all figures
otv.View.ShowAll()
