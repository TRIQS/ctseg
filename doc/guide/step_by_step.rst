.. _step_by_step: 

Step-by-step guide
==================

Below is a step-by-step guide to setting up a CTQMC calculation with the CTSEG solver. 

Step 1 - Choose the construction parameters
*******************************************

We first need to specify the parameters that will define the structure of the solver object. 
They can be conveniently supplied as a Python dictionary::

    constr_params = {
        "gf_struct": gf_struct,
        "beta": 100,
        "n_tau": 1001,
        "n_tau_k": 1001
    }

* ``gf_struct`` specifies the block structure of the impurity Green's function (see below). 

* ``beta`` is the inverse temperature. 

* ``n_tau`` is the number of points of the imaginary time grid on which the fermionic two-point functions (input hybridization :math:`\Delta(\tau)`, measured Green's function :math:`G(\tau)`, etc.) will be sampled. 

* ``n_tau_k`` is the number of points of the imaginary time grid on which the bosonic two-point functions (input dynamical intercation :math:`D_0(\tau)`, measured density correlation function :math:`\chi(\tau)`, etc.) will be sampled. 

Green's function structure
--------------------------

The parameter ``gf_struct`` is a list of pairs of the form ``(block_name, block_size)``. 
The choice of the blocks is not unique, but choosing the smallest possible blocks 
ensures the best performance. The impurity problems solved by CTSEG are always spin-conserving, 
it is therfore always possible to define at least two blocks ``up`` and ``down``. For example, for a 
single orbital::

    gf_struct = [("up", 1), ("down", 1)]

It is suboptimal, but equivalent, to define a single-block structure for the same single-orbital problem::

    gf_struct = [("orb", 2)]

.. note::
    The choice of the block structure affects the Markov chain, unless ``move_move`` is disabled. 

.. warning::
    For a problem with :math:`\mathcal{J}_{\perp}` interactions, the structure **must** have
    two blocks of size 1. 

Step 2 - Construct the solver
*****************************

With the parameters set, we are ready to construct the solver object::

    from triqs_ctseg import Solver
    S = Solver(**constr_params)

Step 3 - Supply the solver inputs
*********************************

We now need to provide the solver with the parameters of the action that define the impurity problem. 
Here, we will give examples for a single-orbital impurity (``gf_struct = [("up", 1), ("down", 1)]``). 

Interaction Hamiltonian
-----------------------

The static part of the electron-electron interaction is supplied as a TRIQS operator. In the CTSEG case,
it is always a density-density interaction. Density operators are represented by ``n(block_name, orbital_number)``. 
For our single orbital problem:: 

    from triqs.operators import *
    h_int = U * n("up", 0) * n("down", 0)

For more complex Hamiltonians (in multi-orbital problems), it is possible to use the ``h_int_density`` function. 
For example, for a two-orbital system (``gf_struct = [("up", 2), ("down", 2)]``), one defines the interaction matrix for electrons of the same spin:: 

    Umat = matrix([[0, V],
                   [V, 0]])

And the interaction matrix for electrons of opposite spin:: 

    Upmat = matrix([[U, V],
                    [V, U]])

Then, the Hamiltonian is constructed as:: 

    from triqs.operators.util import h_int_density
    block_names = ["up", "down"]
    n_orb = 2
    h_int = h_int_density(block_names, n_orb, Umat, Upmat, off_diag = True)

The interaction Hamiltonian is passed to the solver as part of the solve parameters (see below). 

Hybridization function
----------------------

The CTSEG solver takes as an input the hybridization function :math:`\Delta(\tau)` that appears in the 
impurity action (see :doc:`CTSEG algorithm <../algorithm_implementation/ctseg>`). It is initialized as::

    from triqs.gf import *
    from triqs.gf.tools import *
    tau_mesh = MeshImTime(beta, 'Fermion', n_tau)
    Delta_tau = BlockGf(mesh = tau_mesh, gf_struct = gf_struct)

The data in ``Delta_tau`` can be specified manually (``Delta_tau["block_name"].data = ...``), or with one 
of the TRIQS built-in functions. For example::

    iw_mesh = MeshImFreq(beta, 'Fermion', n_tau//2)
    Delta_iw = BlockGf(mesh = iw_mesh, gf_struct = gf_struct)
    Delta_iw << Fourier(Delta_tau)
    Delta_iw << t**2 * SemiCircular(2 * t, chem_potential = mu)
    Delta_tau << Fourier(Delta_iw) 

The hybridization function is supplied to the solver via::

    S.Delta_tau << Delta_tau

This is different from the `CTHYB <https://triqs.github.io/cthyb/latest/>`_ solver, which takes as input the non-interacting impurity 
Green's function :math:`G_0(i\omega_n)`. It is defined as

.. math::
    G_0(i \omega_n) = \frac{1}{i \omega_n - \hat h_0 - \Delta(i \omega_n)}

where :math:`\hat h_0` is the non-interacting impurity Hamiltonian. Note that :math:`\Delta(i\omega_n)` vanishes in the high frequency limit. 
The following code extracts from a ``BlockGf`` ``G0_iw`` a ``BlockGf`` ``Delta_iw`` and a matrix ``h0`` for each block::

    def get_h0_Delta(G0_iw):
        h0_lst, Delta_iw = [], G0_iw.copy()
        for bl in G0_iw.indices:
            Delta_iw[bl] << iOmega_n - inverse(G0_iw[bl])
            tail, err = fit_hermitian_tail(Delta_iw[bl])
            Delta_iw[bl] << Delta_iw[bl] - tail[0]
            h0_lst.append(tail[0])
        return h0_lst, Delta_iw 

    h0_lst, Delta_iw = get_h0_Delta(G0_iw)

The eigenvalues of the ``h0`` provide the orbital dependent chemical potential (orbital energies)::

    import numpy as np
    from numpy import linalg
    mu_list = []
    for h0 in h0_lst:
        mu_list += [-l for l in linalg.eig(h0)[0].real]

If the ``h0`` are already expressed in the local eigenbasis, the obtained hybridization function can be 
directly used as the solver input. In the general case, they need to be rotated to the eigenbasis::

    rot_lst = [np.matrix(linalg.eig(h0_bl)[1]) for h0_bl in h0]
    for (bl, g0_bl), rot_bl in zip(Delta_iw, rot_lst):
        g0_bl << rot_bl.H * g0_bl * rot_bl

If a rotation is required, then the interactions should also be rotated accordingly. 

Chemical potential
------------------

The orbital-dependent chemical potential is passed to the solver as the solve parameter ``hartree_shift`` (see below), 
a list with one value per color. For example, for the single-orbital problem::

    hartree_shift = [mu, mu]

.. warning::

    In CTSEG, the meaning of the ``hartree_shift`` parameter is not the same as in CTHYB. In CTHYB, it represents 
    a shift of the chemical potential with respect to the one already contained in ``G0_iw``. In CTSEG, 
    ``hartree_shift`` is the full chemical potential with all shifts applied. 
    
If the list of chemical potentials is extracted from the CTHYB input ``G0_iw`` as above, then::

        hartree_shift_ctseg = [mu_list[i] + hartree_shift_cthyb[i] for i in range(n_colors)]

Dynamical density-density interaction
-------------------------------------

The dynamical density-density interaction :math:`D(\tau)` (see :doc:`CTSEG algorithm <../algorithm_implementation/ctseg>`) is initialized as:: 

    D_tau = GfImTime(indices = range(n_colors), beta = beta, statistic = "Boson", n_points = n_tau_k)

It is a matrix Green's function, for which no block structure is explicitly enforced. The data in 
``D_tau`` can be specified manually (``D_tau.data = ...``) or by using an analytical expression. 
For example:: 

    from triqs.gf.descriptors import Function
    wp = 1
    D_iw = GfImFreq(indices = range(n_colors), beta = beta, statistic = "Boson", n_points = n_tau_k//2)
    D_iw << Function(lambda w: wp**2 / (w**2 - wp**2))
    D_tau << Fourier(D_iw)

The dynamical interaction is supplied to the solver via::

    S.D0_tau << D_tau

Spin-spin interaction
---------------------

The perpendicular spin-spin interaction :math:`J_{\perp}(\tau)` (see :doc:`CTSEG algorithm <../algorithm_implementation/ctseg>`) is initialized as:: 

    Jperp_tau = GfImTime(indices = [0], beta = beta, statistic = "Boson", n_points = n_tau_k)

It is a :math:`1 \times 1` matrix Green's function. It is supplied to the solver via::

    S.Jperp_tau << Jperp_tau

If the impurity action contains a spin-spin interaction term of the form :math:`(1/2) \cdot Q(\tau - \tau') \sum_i s_i(\tau) s_i(\tau')`, 
it can be split into a density-density and a perpendicular spin-spin term. Indeed, making use of symmetry properties, we may replace in the action 

.. math::

    \frac{1}{2} Q(\tau - \tau') \sum_i s_i(\tau) \cdot s_i(\tau') \mapsto \frac{1}{2} Q(\tau - \tau') \left[s^+(\tau) s^-(\tau') + \frac{1}{4}\sum_{\sigma \sigma'} (-1)^{\sigma \sigma'} n_{\sigma}(\tau) n_{\sigma'}(\tau') \right]

The solver is then accordingly set up as:: 

    S.Jperp_tau << Q_tau
    S.D0_tau[0, 0] << D_tau + (1/4) * Q_tau[0, 0]
    S.D0_tau[0, 1] << D_tau - (1/4) * Q_tau[0, 0]
    S.D0_tau[1, 0] << D_tau - (1/4) * Q_tau[0, 0]
    S.D0_tau[1, 1] << D_tau + (1/4) * Q_tau[0, 0]

Solve parameters
----------------

The parameters required to perform a CTQMC run are conveniently supplied as a Python dictionary.
The following parameters need to be specified for every run. For example::

    solve_params = {
        "h_int": h_int, 
        "hartree_shift": hartree_shift,
        "length_cycle": 50,
        "n_warmup_cycles": 20000, 
        "n_cycles": 200000 
        }

* ``h_int`` is the local interaction Hamiltonian (see above).

* ``hartree_shift`` is the total orbital-dependent chemical potential (see above). 

* ``length_cycle`` is the length of a Monte Carlo cycle. Observables are sampled every ``length_cycle`` Monte Carlo moves (either accepted or rejected). 

* ``n_warmup_cycles`` is the number of cycles to do before any observables are samples, so as to "forget" the initial configuration. 

* ``n_cycles`` is the number of cycles used for the production run. 

Other parameters include: 

* **Measure control**. All the :doc:`measurements <measurements>` can be switched on and off. Some of the measurements (self-energy improved estimator,
  density correlation functions) can be time-consuming, and they are off by default. For example, to turn the improved estimator 
  measurement on, one should set ``solve_params["measure_ft"] = True``. 

* **Move control**. All the :doc:`Monte Carlo moves <moves>` can be switched on and off. This functionality exists to facilitate testing
  for developers. The solver chooses the relevant moves depending on its inputs, and regular users should not need move control. 

The complete list of parameters is available :doc:`here <../_ref/triqs_ctseg.solver.Solver.solve>`.

Step 4 - Run the solver 
***********************

The CTQMC run is triggered by::

    S.solve(**solve_params)

.. warning::
    The solver prints to the command line the interaction matrix ``U`` and chemical potential ``mu`` that are used internally. 
    In the presence of dynamical interactions, these are renormalized values, different from the input parameters contained 
    in ``h_int`` and ``hartree_shift`` (see :doc:`CTSEG algorithm <../algorithm_implementation/ctseg>`).

After it is done accumulating, the solver prints the average acceptance rates. Very low acceptance rates for all moves (below 0.01)
are generally a sign that something went wrong. However, some of the moves (``split_spin_segment``, ``regroup_spin_segment``)
often have low acceptance rates, even if the calculation runs as it should. 

Step 5 - Access the results
***************************

The results of the accumulation are stored in ``S.results``. For example, the impurity Green's function is accessed with::

    G_tau = S.results.G_tau

The results can be analyzed using the TRIQS plotting tools (``oplot``) or by directly extracting the data:: 

    G_up_data = G_tau["up"].data
    G_down_data = G_tau["down"].data
    times = np.array([tau.value for tau in G_tau.mesh])

For a rotationally-invariant impurity, the spin-spin correlation function :math:`\chi(\tau) = (1/3) \sum_{i = x, y, z} \langle s_i(\tau) s_i(0) \rangle` can be obtained from the density-density
correlation function::

    nn_tau = S.results.nn_tau
    chi_tau = 0.25 * (nn_tau[0, 0] + nn_tau[1, 1] - nn_tau[0, 1] - nn_tau[1, 0])

If rotational invariance is broken (for instance, in the presence of a Zeeman field), one needs to measure separately the 
perpendicular spin-spin correlation function:: 

    nn_tau = S.results.nn_tau
    sperp_tau = S.results.sperp_tau
    chi_tau = (0.25 * (nn_tau[0, 0] + nn_tau[1, 1] - nn_tau[0, 1] - nn_tau[1, 0]) + 2*sperp_tau[0, 0]) / 3

Step 6 - Save the results
*************************

The TRIQS ``h5`` module is convenient for saving the results to an hdf5 file. It is possible to save all the data at once 
by saving the solver object::

    from h5 import *
    with HDFArchive("results.h5", "a") as A:
        A["Solver"] = S

Running the solver in parallel
******************************

The CTSEG solver supports MPI parallelism. If the solver run is set up in a file ``script.py``, a parallel run 
is typically achieved with the command::

    mpirun -np <n_cores> python script.py

Each core then runs its own Markov chain of length ``n_cycles`` (starting from a different random number generator seed) 
and at the end the results from the different cores are averaged together. 

