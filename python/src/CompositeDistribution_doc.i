%feature("docstring") OT::CompositeDistribution
R"RAW(Composite distribution.

Helper class for defining the push-forward distribution of a given univariate
distribution by a given scalar function.

We note :math:`X` a scalar random variable which distribution is :math:`distX`,
which probability density function is :math:`f_X`.

Then :math:`distY` is the distribution of the scalar random variable
:math:`Y=g(X)`, which probability density function :math:`f_Y` is defined as:

.. math::

    \displaystyle f_Y(y) = \sum_{k =0}^{k=N} \frac{f_X (g^{-1}(y))}{|g'\circ g^{-1}(y)|}1_{y \in g([a_k, a_{k+1}))}

with :math:`a_0=\inf \supp{f_X}`, :math:`a_N=\sup \supp{f_X}` and
:math:`(a_1, \dots, a_N)` such that :math:`g` is monotone over
:math:`[a_k, a_{k+1})` for :math:`0 \leq k \leq N`.

Both constructors assume :math:`g` is piecewise monotone on the support of
:math:`distX`, i.e. there exists a partition
:math:`a_0 \leq a_1 \leq \dots \leq a_N` (non-decreasing) with :math:`g` monotone on each
:math:`[a_k, a_{k+1}]` (or :math:`(a_k, a_{k+1})` when :math:`a_k=a_{k+1}`).
:math:`g` is assumed continuous on each open piece;
discontinuities and singularities (vertical asymptotes) inside the support
are not detected and lead to an incorrect range, PDF and CDF. A repeated
value :math:`a_k=a_{k+1}` is required when a singular point has distinct
left and right limits (e.g. :math:`1/x` at :math:`x=0`) so that the two
one-sided values can be given separately in :math:`v`.

Available constructors:
    CompositeDistribution(*g, distX*)

    CompositeDistribution(*g, distX, a, v*)

The first constructor (*auto-detection*) determines the partition
automatically by searching the zeros of :math:`g'` on the numerical support
of :math:`distX`. The interval :math:`[x_{\min}, x_{\max}]` (range shrunk by
*Distribution-DefaultQuantileEpsilon*) is split into :math:`N` regular
sub-intervals and a :class:`~openturns.Brent` solver brackets each root of
:math:`g'`. It requires :math:`g` to be differentiable with finite derivative
and finite values on the sampled points and roots, otherwise a
``NotDefinedException`` is raised. If :math:`g` oscillates more than
:math:`N` times some monotone pieces may be missed.

The second constructor (*user-provided partition*) does not evaluate
:math:`g'` and must be used to obtain correct results when :math:`g` is
discontinuous or has singularities, or to provide the exact partition for a
piecewise monotone :math:`g`. The bounds :math:`a` must include any
discontinuity or singular point — with a repeated value :math:`a_k=a_{k+1}`
when the singularity has distinct left and right limits — and :math:`v` may
contain :math:`\pm \infty` (``SpecFunc.Infinity``) to encode vertical
asymptotes (one entry per side of the repeated bound); the range is then the
interval :math:`[\min v, \max v]`. A :class:`~openturns.Mixture` of
composites on each continuous piece is an alternative workaround.

Parameters
----------
g : :class:`~openturns.Function`, :math:`\Rset \rightarrow \Rset`
    Must have input and output dimension 1.
distX : :class:`~openturns.Distribution`, univariate
a : sequence of float of dimension :math:`N+1`, :math:`a[0]=\inf \supp{distX}`, :math:`a[N]=\sup \supp{distX}`
    The bounds of the intervals on which :math:`g` is monotone, sorted in
    non-decreasing (ascending) order; repeated values are allowed.
    Must include any discontinuity or singularity; a repeated bound
    :math:`a_k=a_{k+1}` is required when the singularity has distinct
    one-sided limit values, so that :math:`v[k]` and :math:`v[k+1]` can
    encode the left and right limits separately (see example with
    :math:`1/x` at :math:`x=0`).
v : sequence of float of dimension :math:`N+1`,
    The values taken by :math:`g` on each bound: :math:`v[k]=g(a[k])`,
    or the one-sided limit at a repeated bound.
    May contain :math:`\pm \infty` (as ``SpecFunc.Infinity``) at singularities;
    a singularity with distinct one-sided values requires the two
    corresponding entries of :math:`v` (e.g. :math:`-\infty` and
    :math:`+\infty` for :math:`1/x` at :math:`0`).

Returns
-------
distY : :class:`~openturns.Distribution`, univariate
    :math:`distY` is the push-forward distribution of :math:`distX` by :math:`g`.

Notes
-----
Its first moments are obtained by numerical integration.
The mathematical support of a CompositeDistribution is defined as the image of the numerical support of the antecedent.
Consequently, its bounds are finite.

- *CompositeDistribution-StepNumber* (int, default 256): :math:`N`, number of
  sub-intervals for the automatic monotonicity detection in the first
  constructor.
- *CompositeDistribution-SolverEpsilon* (float, default 1e-14): accuracy of
  the :class:`~openturns.Brent` solver used to invert :math:`g` and to find
  zeros of :math:`g'`.
- *Distribution-DefaultQuantileEpsilon* (float): margin used to shrink the
  antecedent range to its numerical support.

Examples
--------
Create a distribution:

>>> import openturns as ot
>>> g = ot.SymbolicFunction(['x'], ['sin(x) + cos(x)'])
>>> distY = ot.CompositeDistribution(g, ot.Normal(1.0, 0.5))

>>> g = ot.SymbolicFunction(['x'], ['abs(x)'])
>>> a = [-1.0, 0.0, 2.0]
>>> v = [1.0, 0.0, 2.0]
>>> distZ = ot.CompositeDistribution(g, ot.Uniform(-1.0, 2.0), a, v)

>>> distX = ot.Normal(0.0, 1.0)
>>> a0 = distX.getRange().getLowerBound()
>>> aN = distX.getRange().getUpperBound()
>>> a = [a0[0], 0.0, 0.0, aN[0]]
>>> g = ot.SymbolicFunction(['x'], ['1.0/x'])
>>> v = [g(a0)[0], -ot.SpecFunc.Infinity, ot.SpecFunc.Infinity, g(aN)[0]]
>>> distT = ot.CompositeDistribution(g, distX, a, v)


Draw a sample:

>>> sample = distT.getSample(5)
)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::CompositeDistribution::getFunction
R"RAW(Accessor to the function.

Returns
-------
g :  :class:`~openturns.Function`, :math:`\Rset \rightarrow \Rset`
    the function :math:`g`.)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::CompositeDistribution::getAntecedent
"Accessor to the antecedent distribution.

Returns
-------
distX : :class:`~openturns.Distribution`, univariate
    Antecedent distribution :math:`distX`."

// ---------------------------------------------------------------------

%feature("docstring") OT::CompositeDistribution::setFunction
R"RAW(Fix the function through which the distribution is push-forwarded.

Parameters
----------
g :  :class:`~openturns.Function`, :math:`\Rset \rightarrow \Rset`
    the function :math:`g`.)RAW"
// ---------------------------------------------------------------------

%feature("docstring") OT::CompositeDistribution::setAntecedent
"Fix the antecedent distribution which is push-forwarded.

Parameters
----------
distX : :class:`~openturns.Distribution`, univariate
   Distribution of the antecedent :math:`distX`."

// ---------------------------------------------------------------------

%feature("docstring") OT::CompositeDistribution::setSolver
"Solver accessor.

Parameters
----------
solver : :class:`~openturns.Solver`
    The solver used for PDF/CDF computations.
"

// ---------------------------------------------------------------------

%feature("docstring") OT::CompositeDistribution::getSolver
"Solver accessor.

Returns
-------
solver : :class:`~openturns.Solver`
    The solver used for PDF/CDF computations.
"
