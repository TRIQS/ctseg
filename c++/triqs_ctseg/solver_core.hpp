// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once

#include <optional>
#include "params.hpp"
#include "work_data.hpp"
#include "inputs.hpp"
#include "results.hpp"
#include <triqs/utility/macros.hpp>

namespace triqs_ctseg {

  /// Continuous-time hybridization-expansion quantum Monte Carlo solver.
  class solver_core {

    // Inverse temperature
    double beta;

    // The set of inputs
    inputs_t inputs;

    // mpi communicator
    mpi::communicator c;

    public:
    /// Parameters used for constructing the solver.
    constr_params_t constr_params;

    /// Parameters passed to the ``solve()`` method.
    solve_params_t solve_params;

    /// Container for all results accumulated by the CTQMC simulation.
    results_t results;

    /**
     * @brief Initialize the solver.
     * @param p Parameters used for constructing the solver class.
     */
    solver_core(constr_params_t const &p);

    /**
     * @brief Solve the impurity problem.
     * @param p Parameters controlling the MC simulation and measurements.
     */
    void solve(solve_params_t const &p);

    // Green's function views for Python interface
    // do NOT add const here : python uses a non const object and non const view

    /// Hybridization function \f$ \Delta(\tau) \f$.
    C2PY_PROPERTY_GET(Delta_tau) block_gf_view<imtime> Delta_tau() { return inputs.Delta; }

    /// Dynamical spin-spin interaction \f$ \mathcal{J}_\perp(\tau) \f$.
    C2PY_PROPERTY_GET(Jperp_tau) gf_view<imtime> Jperp_tau() { return inputs.Jperpt; }

    /// Dynamical density-density interaction \f$ D_0(\tau) \f$.
    C2PY_PROPERTY_GET(D0_tau) block2_gf_view<imtime> D0_tau() { return inputs.D0t; }

    // --------------- h5 -------------------------
    static std::string hdf5_format() { return "CTSEG_SolverCore"; }
    friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s);
    C2PY_IGNORE static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name);
  };

} // namespace triqs_ctseg
