
.. _options:

Solver options and accessors
=============================

The solver consists of one main object with a constructor and a ``solve()`` method, and *accessors*. 
The options of the constructor and ``solve()`` method are used to give more information as to the desired run, such as number of Matsubara frequencies, number of Monte-Carlo cycles, or which observables have to be measured. The accessors are members of the solver class which are used both to set the inputs of the solver (like the non-interacting Green's function ``G0_iw`` or the dynamical interactions ``D0_iw``) and to read out the observables at the end of the computation (such as ``G_tau``, ``G_iw`` or ``nn_iw``).

The options of the constructor are:

+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+
| Parameter Name | Type        | Default | Documentation                                                                                            |
+================+=============+=========+==========================================================================================================+
| beta           | double      | --      | Inverse temperature                                                                                      |
+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+
| gf_struct      | gf_struct_t | --      | Structure of the GF (names, sizes of blocks)                                                             |
+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+
| n_tau          | int         | 10001   | Number of time slices for :math:`Delta(\tau)`/:math:`G(\tau)`/:math:`F(\tau)`                            |
+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+
| n_tau_k        | int         | 10001   | Number of time slices for :math:`K(\tau)`                                                                |
+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+
| n_tau_jperp    | int         | 10001   | Number of time slices for :math:`J_\perp(\tau)`                                                          |
+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+
| n_tau_nn       | int         | 101     | Number of Legendre coefficients for G(l)                                                                 |
+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+
| n_w_b_nn       | int         | 32      | Number of bosonic Matsub. freqs for :math:`nn(i\omega)`, :math:`\mathcal{D}_0(i\omega)`, :math:`J_\perp` |
+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+
| n_iw           | int         | 1025    | Number of fermionic Matsubara frequencies for :math:`G_0(i\omega)`, :math:`G`, :math:`F`, :math:`\Sigma` |
+----------------+-------------+---------+----------------------------------------------------------------------------------------------------------+


The options of the `solve()` method and the solver accessors are summarized below. They are identical in the C++ and Python interface.

C++ interface
-------------------

.. toctree::
   :maxdepth: 1

   ../reference/triqs_ctseg/solver_core



Python interface
----------------------


.. autoclass:: triqs_ctseg.SolverCore
   :members:
