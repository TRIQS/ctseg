#pragma once

#include <optional>
#include "params.hpp"
#include "work_data.hpp"
#include "inputs.hpp"
#include "results.hpp"

namespace triqs_ctseg {

  /// Main solver class
  class solver_core {

    // Inverse temperature
    double beta;

    // The set of inputs
    inputs_t inputs;

    // mpi communicator
    mpi::communicator c;

    public:
    /**Set of parameters used in the construction of the ``solver_core`` class.
  *
  *.. include:: ../../python/triqs_ctseg/parameters_constr_params_t.rst 
  */
    constr_params_t constr_params;

    /**Set of parameters used by the last call to ``solve()``.
  *
  *.. include:: ../../python/triqs_ctseg/parameters_solve_params_t.rst 
  */
    solve_params_t solve_params;

    /// The set of results. See :doc:`Measurements <../guide/measurements>`.
    // Will be passed to measures and initialized by them.
    results_t results;

    /// Initialize the solver
    CPP2PY_ARG_AS_DICT solver_core(constr_params_t const &p);

    /// Solve the impurity problem
    CPP2PY_ARG_AS_DICT void solve(solve_params_t const &p);

    // Green's function views for Python interface
    // do NOT add const here : python uses a non const object and non const view

    /// Hybridization function :math:`\Delta(\tau)`
    block_gf_view<imtime> Delta_tau() { return inputs.delta; }

    /// Dynamical spin-spin interaction :math:`\mathcal{J}_\perp(\tau)`
    gf_view<imtime> Jperp_tau() { return inputs.jperpt; }

    /// Dynamical density-density interaction :math:`D_0(\tau)`
    block2_gf_view<imtime> D0_tau() { return inputs.d0t; }

    // --------------- h5 -------------------------
    CPP2PY_IGNORE static std::string hdf5_format() { return "CTSEG_SolverCore"; }
    friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s);
    CPP2PY_IGNORE static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name);
  };

} // namespace triqs_ctseg
