#pragma once

#include <optional>
#include "params.hpp"
#include "work_data.hpp"
#include "inputs.hpp"
#include "results.hpp"

/// Main solver class
class solver_core {

  // Inverse temperature
  double beta;

  // Keep the construction parameter
  constr_params_t constr_params;

  // Keep the last solver parameters from the last call
  std::optional<solve_params_t> last_solve_params;

  // The set of inputs
  inputs_t inputs;

  // mpi communicator
  mpi::communicator c;

  public:
  // The set of results. Will be passed to measures and initialized by them.
  results_t results;

  CPP2PY_ARG_AS_DICT solver_core(constr_params_t const &p);

  CPP2PY_ARG_AS_DICT void solve(solve_params_t const &p);

  // Green's function views for Python interface
  // do NOT add const here : python use a non const object and non const view

  /// Hybridization function $\Delta^\sigma_{ab}(\tau)$
  block_gf_view<imtime> Delta_tau() { return inputs.delta; }

  /// Dynamical spin-spin interactions $\mathcal{J}_\perp(\tau)$
  gf_view<imtime> Jperp_tau() { return inputs.jperpt; }

  /// Dynamical density-density interactions $D_0(\tau)$
  gf_view<imtime> D0_tau() { return inputs.d0t; }

  // ---------------h5 -------------------------
  static std::string hdf5_format() { return "CTSEG_J2"; }
  friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s);
  CPP2PY_IGNORE static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name);
};
