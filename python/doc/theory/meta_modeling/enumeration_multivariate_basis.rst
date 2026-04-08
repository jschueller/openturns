.. _enumeration_multivariate_basis:

Tensorized multivariate basis enumeration functions
---------------------------------------------------

Enumeration functions (refer to :ref:`enumeration_strategy`) help to enumerate a multivariate basis
built as the tensorization of univariate basis, using the indexation of each marginal basis.

Such a  multivariate basis is used, for example, in the functional chaos expansion setting
(refer to  :ref:`functional_chaos`).

Let consider some :math:`\inputDim` univariate basis, denoted by
:math:`(\pi_{k}^{(i)})_{k \geq 0}` for :math:`1 \leq i \leq \inputDim`, where each
:math:`\pi_{k}^{(i)}: \Rset \rightarrow \Rset`.

Let denote by :math:`\{\psi_{\vect{\alpha}},\vect{\alpha} \in \Nset^\inputDim\}` a multivariate basis
built as the tensorization of the univariate basis. The multivariate basis term
:math:`\psi_{\vect{\alpha}}` is defined by the product:

.. math::

    \psi_{\vect{\alpha}} (\vect{\xi}) = \pi_{\alpha_1}^{(1)}(\xi_1) \times \dots \times
    \pi_{\alpha_{\inputDim}}^{({\inputDim})}(\xi_{\inputDim})

for any :math:`\vect{\xi} \in \mathbb{R}^{\inputDim}`.

When the univariate basis are polynomials such that :math:`\alpha_i`
is the degree of :math:`\pi_{\alpha_i}^{(i)}`, then the multi-index represents
the marginal degrees of the polynomial :math:`\psi_{\vect{\alpha}}`. In that case,
the length of the multi-index is the total degree of the polynomial.

Several enumeration functions can be used:

- the linear enumeration function,
- the hyperbolic enumeration function,
- the anisotropic hyperbolic enumeration function,
- the infinity norm enumeration function.

We detail the interest of each one in the functional chaos expansion setting.

Linear enumeration function
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The linear enumeration function (defined in :ref:`enumeration_strategy`, equation
:eq:`linearEnumFct`) explores the marginal degrees linearly, in order
of increasing total degree. The exploration is symmetric with respect to the components.

Hyperbolic enumeration function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The hyperbolic enumeration function is inspired by the so-called
*sparsity-of-effects principle*, which states that most models are
principally governed by main effects and low-order interactions.
Accordingly, one wishes to define an enumeration function which first
selects those multi-indices related to main effects, i.e. with a
reasonably small number of nonzero components, prior to selecting
those associated with higher-order interactions.

The hyperbolic enumeration functions (defined in :ref:`enumeration_strategy`, equation
:eq:`hyperBolEnumFct`) are based on the q-norm  (defined in :ref:`enumeration_strategy`,
equation :eq:`eq_q_set`).
They penalize the indices associated with high-order interactions all the more since :math:`q` is low.

The library implements the partition based on the length of the multi-index.

Anisotropic hyperbolic enumeration function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The anisotropic hyperbolic enumeration functions (defined in :ref:`enumeration_strategy`, equation
:eq:`anisotropEnumFct`) lead to first select the basis terms
depending on a specific subset of input variables. This function emphasizes multivariate polynomials
for whose components associated to small coefficients are polynomials of high degree.

Therefore, when the model is governed by main effects of specific inputs, the associated weights
should be small compared to the other ones.

Infinity norm enumeration function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The infinity norm enumeration function (defined in :ref:`enumeration_strategy`, equation
:eq:`infEnumFct`) is based on polynomials with fixed maximum
marginal degree in ascending order.

This function is used to define the polynomial space as large as possible with fixed maximum
marginal degree: when a tensorized Gaussian quadrature formula
(:class:`~openturns.GaussProductExperiment`) is obtained by tensorization of univariate formulas,
this space (with the appropriate marginal degrees) is the space of polynomials whose integrals are
computed exactly by the quadrature rule.

.. topic:: API:

    - See :class:`~openturns.LinearEnumerateFunction`
    - See :class:`~openturns.HyperbolicAnisotropicEnumerateFunction`
    - See :class:`~openturns.NormInfEnumerateFunction`

.. topic:: Examples:

    - See :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_functional_chaos`
    - See :doc:`/auto_surrogate_modeling/polynomial_chaos/plot_enumeratefunction`
    - See :doc:`/auto_surrogate_modeling/fields_surrogate_models/plot_fieldfunction_metamodel`

.. topic:: References:

    - [blatman2009]_
