.. _uniform_random_generator:

Uniform Random Generator
------------------------

Generating simulations according to a distribution is based on
generating simulations according to a Uniform distribution on
:math:`[0,1]` : several techniques exist then to transform a
realization according to a uniform distribution onto a realization
according to a distribution which cumulative distribution function is
:math:`F` (refer to for each distribution).

Thus, the quality of the random generation of simulation is entirely
based on the quality of the *deterministic* algorithm which simulates
realizations of the Uniform(0,1) distribution.

We use the `PGC <https://www.pcg-random.org/index.html>`_ algorithm,
which stands for Permuted Congruential Generator.
The specific member of the PCG family that we use is PCG XSL RR 128-bit state / 64-bit output
as described `here <https://www.cs.hmc.edu/tr/hmc-cs-2014-0905.pdf>`_

.. topic:: API:

    - See :class:`~openturns.RandomGenerator`
    - See :class:`~openturns.RandomGeneratorState` to save the generator state

.. topic:: Examples:

    - See :doc:`/auto_numerical_methods/general_methods/plot_random_generator`

