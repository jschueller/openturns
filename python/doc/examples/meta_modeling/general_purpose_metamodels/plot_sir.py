"""
Sliced inverse regression on 2-d data
=====================================
"""

# %%
import openturns as ot
import openturns.experimental as otexp
import matplotlib.pyplot as plt

# %%
# Create 2-d data X and 1-d feature Y
N = 100
X = ot.Normal([0.0] * 2, [0.1] * 2).getSample(N)
X += [[-i / (N - 1), 2 * i / (N - 1)] for i in range(N)]
X = X - X.computeMean()
f = ot.SymbolicFunction(["x1", "x2"], ["4*(x1+2*x2)+2"])
Y = f(X) + ot.Normal(0.0, 0.2).getSample(N)

# %%
# Plot data
plt.scatter(X.getMarginal(0), X.getMarginal(1), c=Y)
plt.xlabel('x1')
plt.ylabel('x2')
plt.axis('square')
plt.title("2D colored dataset")

# %%
# Run the SIR algorithm
algo = otexp.SlicedInverseRegression(X, Y)
algo.run()
transformation = algo.getResult().getTransformation()

# %%
# Show SIR direction
plt.scatter(X.getMarginal(0), X.getMarginal(1), c=Y)
plt.xlabel('x1')
plt.ylabel('x2')
plt.axis('square')
plt.title("First direction of SIR")
dir_sir = algo.getResult().getDirections()
plt.plot([-dir_sir[0, 0] * 30, dir_sir[0, 0] * 30], [-dir_sir[1, 0] * 30, dir_sir[1, 0] * 30], c='red')

# %%
plt.show()
