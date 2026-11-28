%feature("docstring") OT::CuhreIntegration
R"RAW(Multidimensional integration algorithm.

Parameters
----------
algorithmName : str
    The name of the integration routine to use. Only 'cuhre' is available.

Notes
-----
The Cuhre algorithm approximates the definite integral:

.. math::

    \int_{I_n} f(x_1, ..., x_n) \di{x_1}...\di{x_n}


which takes values in :math:`\Rset^p`, with :math:`f: \Rset^n \mapsto \Rset^p`,
:math:`I_n=[a_1, b_1]\times\cdots\times[a_n, b_n]\subset\Rset^n`,
having :math:`a_i<b_i` for all :math:`i=1,...,n`.

Cuhre is deterministic: it combines Genz-Malik cubature rules of degree 7, 9,
11 or 13 for subregion estimation within an adaptative subdivision scheme
[berntsen1991]_. The integration domain must be at least 2-dimensional.
The result is returned as a :class:`~openturns.Point` of dimension
:math:`p`, which is the number of components of the return value of the
:class:`~openturns.Function` that one will integrate.

The algorithm and its parameters are discussed in [hahn2005]_.

The following :class:`~openturns.ResourceMap` keys are used:

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Key
     - Description
   * - ``CuhreIntegration-DefaultMaximumAbsoluteError`` (``Scalar``, default: ``0.0001``)
     - Default absolute error.
   * - ``CuhreIntegration-DefaultMaximumRelativeError`` (``Scalar``, default: ``0.0001``)
     - Default relative error.
   * - ``CuhreIntegration-DefaultMaximumCallsNumber`` (``UnsignedInteger``, default: ``1000000``)
     - Default maximum number of integrand calls.
   * - ``CuhreIntegration-mineval`` (``UnsignedInteger``, default: ``0``)
     - Minimum number of integrand evaluations.
   * - ``CuhreIntegration-key`` (``UnsignedInteger``, default: ``0``)
     - Cuhre integration key: 7, 9, 11 (3D only) or 13 (2D only) select the
       cubature rule, 0 selects the default rule (13 in 2D, 11 in 3D, 9
       otherwise).
   * - ``CuhreIntegration-flags`` (``UnsignedInteger``, default: ``0``)
     - Global flags, only bit 2 is used: if set, only the last iteration
       contributes to the final result.

Examples
--------
Create a Cuhre algorithm:

>>> import openturns as ot
>>> algo = ot.experimental.CuhreIntegration('cuhre')
)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::integrate
R"RAW(Evaluation of the integral of :math:`f w` on an interval.

Parameters
----------
f : :class:`~openturns.Function`, :math:`f: \Rset^n \mapsto \Rset^p`
    The integrand function.
interval : :class:`~openturns.Interval`, :math:`interval \subset \Rset^n`
    The integration domain.

Returns
-------
value : :class:`~openturns.Point`
    Approximation of the integral. Its dimension is :math:`p`.

Examples
--------
>>> import openturns as ot
>>> f = ot.SymbolicFunction(['x', 'y', 'z'], ['sin(x) * cos(y) * exp(z)'])
>>> interval = ot.Interval([0.0] * 3, [1.0] * 3)
>>> algoC = ot.experimental.CuhreIntegration('cuhre')
>>> value = algoC.integrate(f, interval)  # doctest: +SKIP
)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::getMaximumRelativeError
"Accessor to the relative accuracy required in the integration.

Returns
-------
maximumRelativeError : float
    The relative accuracy that the *integrate* method attempts to achieve."

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::setMaximumRelativeError
"Accessor to the relative accuracy required in the integration.

Parameters
----------
maximumRelativeError : float
    The relative accuracy that the *integrate* method attempts to achieve."

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::getMaximumAbsoluteError
"Accessor to the absolute accuracy required in the integration.

Returns
-------
maximumAbsoluteError : float
    The absolute accuracy that the *integrate* method attempts to achieve."

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::setMaximumAbsoluteError
"Accessor to the absolute accuracy required in the integration.

Parameters
----------
maximumAbsoluteError : float
    The absolute accuracy that the *integrate* method attempts to achieve."

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::getAlgorithmName
"Accessor to the name of the integration routine to be used.

Returns
-------
algorithmName : str
    The name of the integration routine used by the algorithm."

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::GetAlgorithmNames
"Retrieves the names of the available integration algorithms.

Returns
-------
algoName : :class:`~openturns.Description`
    The names of the available integration algorithms."

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::setAlgorithmName
"Accessor to the name of the integration routine to be used.

Parameters
----------
algorithmName : str
    The name of the integration routine used by the algorithm. Only 'cuhre' is available."

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::getMaximumCallsNumber
"Accessor to the maximal number of integrand calls.

Returns
-------
N : int
    The maximal number of times the algorithm will attempt to call the integrand."

// ---------------------------------------------------------------------

%feature("docstring") OT::CuhreIntegration::setMaximumCallsNumber
"Accessor to the maximal number of integrand calls.

Parameters
----------
N : int, :math:`N>0`
    The maximal number of times the algorithm will attempt to call the integrand."
