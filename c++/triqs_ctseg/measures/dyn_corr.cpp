// Copyright (c) 2026--present, The Simons Foundation
// Copyright (c) 2026--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./dyn_corr.hpp"
#include "./dyn_density_prefactors.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  dyn_corr::dyn_corr(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta          = p.beta;
    n_color       = config.n_color();
    n_tau_bosonic = p.n_tau_bosonic;

    phi_n   = nda::zeros<double>(n_color, n_color);
    phi_phi = nda::zeros<double>(n_color, n_color);
  }

  // -------------------------------------

  void dyn_corr::accumulate(double s) {

    LOG("\n =================== MEASURE RETARDED STATIC CORRELATIONS  ================ \n");

    Z += s;

    auto accumulate_interval = [&](int density_color, tau_t tau_left, tau_t tau_right) {
      for (int source_color = 0; source_color < n_color; ++source_color)
        phi_n(source_color, density_color) +=
           s * retarded_density_prefactor_integral(wdata, config, source_color, tau_left, tau_right);
    };

    for (int density_color = 0; density_color < n_color; ++density_color) {
      for (auto const &seg : config.seglists[density_color]) {
        if (is_full_line(seg)) {
          accumulate_interval(density_color, tau_t::beta(), tau_t::zero());
        } else if (is_cyclic(seg)) {
          accumulate_interval(density_color, tau_t::beta(), seg.tau_cdag);
          accumulate_interval(density_color, seg.tau_c, tau_t::zero());
        } else {
          accumulate_interval(density_color, seg.tau_c, seg.tau_cdag);
        }
      }
    }

    double const dtau = beta / (n_tau_bosonic - 1);
    auto phi          = nda::zeros<double>(n_color);
    for (int i = 0; i < n_tau_bosonic; ++i) {
      auto tau = tau_t(i == n_tau_bosonic - 1 ? beta : i * dtau);
      for (int source_color = 0; source_color < n_color; ++source_color)
        phi(source_color) = retarded_density_prefactor(wdata, config, source_color, tau);

      double weight = s * dtau * ((i == 0 or i == n_tau_bosonic - 1) ? 0.5 : 1.0);
      for (int a = 0; a < n_color; ++a)
        for (int b = 0; b < n_color; ++b) phi_phi(a, b) += weight * phi(a) * phi(b);
    }
  }

  // -------------------------------------

  void dyn_corr::collect_results(mpi::communicator const &c) {

    Z       = mpi::all_reduce(Z, c);
    phi_n   = mpi::all_reduce(phi_n, c);
    phi_phi = mpi::all_reduce(phi_phi, c);

    phi_n   = phi_n / (Z * beta);
    phi_phi = phi_phi / (Z * beta);

    results.dyn_phi_n   = std::move(phi_n);
    results.dyn_phi_phi = std::move(phi_phi);
  }

} // namespace triqs_ctseg::measures
