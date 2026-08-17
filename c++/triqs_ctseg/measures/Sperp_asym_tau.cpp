// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./Sperp_asym_tau.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  Sperp_asym_tau::Sperp_asym_tau(params_t const &p, work_data_t const &wdata, configuration_t const &config,
                                 results_t &results)
     : wdata{wdata}, config{config}, results{results}, beta{p.beta} {

    sm_sp_tau   = gf<imtime>({beta, Boson, p.n_tau_chi2}, {1, 1});
    sp_sm_tau   = gf<imtime>({beta, Boson, p.n_tau_chi2}, {1, 1});
    sm_sp_tau() = 0;
    sp_sm_tau() = 0;
    Z           = 0;
  }

  // -------------------------------------

  void Sperp_asym_tau::accumulate(double s) {

    LOG("\n =================== MEASURE < S-(tau)S+(0) > and < S+(tau)S-(0) > ================ \n");

    Z += s;

    for (auto const &[k, line] : itertools::enumerate(config.Jperp_list)) {
      auto dtau_sm_sp = double(line.tau_Sminus - line.tau_Splus);
      auto dtau_sp_sm = double(line.tau_Splus - line.tau_Sminus);
      // Jperp lines carry a -J/2 expansion weight; the oriented estimator is
      // four times the corresponding S_x S_x contribution.
      sm_sp_tau[closest_mesh_pt(dtau_sm_sp)] += 2.0 / (real(wdata.Jperp(dtau_sm_sp)(0, 0)));
      sp_sm_tau[closest_mesh_pt(dtau_sp_sm)] += 2.0 / (real(wdata.Jperp(dtau_sp_sm)(0, 0)));
    }
  }

  // -------------------------------------

  void Sperp_asym_tau::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    sm_sp_tau = mpi::all_reduce(sm_sp_tau, c);
    sp_sm_tau = mpi::all_reduce(sp_sm_tau, c);
    sm_sp_tau = sm_sp_tau / (-beta * Z * sm_sp_tau.mesh().delta());
    sp_sm_tau = sp_sm_tau / (-beta * Z * sp_sm_tau.mesh().delta());

    // Fix the point at zero and beta.
    sm_sp_tau[0] *= 2;
    sm_sp_tau[sm_sp_tau.mesh().size() - 1] *= 2;
    sp_sm_tau[0] *= 2;
    sp_sm_tau[sp_sm_tau.mesh().size() - 1] *= 2;

    results.Sminus_Splus_tau = std::move(sm_sp_tau);
    results.Splus_Sminus_tau = std::move(sp_sm_tau);
  }

} // namespace triqs_ctseg::measures
