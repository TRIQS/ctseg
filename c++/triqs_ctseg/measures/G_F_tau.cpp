// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./G_F_tau.hpp"
#include "../logs.hpp"
#include <triqs/utility/legendre.hpp>
#include <cmath>

namespace triqs_ctseg::measures {

  G_F_tau::G_F_tau(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta          = p.beta;
    measure_G_tau = p.measure_G_tau;
    measure_F_tau = p.measure_F_tau and wdata.rot_inv;
    measure_G_l   = p.measure_G_l;
    measure_F_l   = p.measure_F_l and wdata.rot_inv;
    gf_struct     = p.gf_struct;

    if (measure_G_tau or measure_F_tau) {
      G_tau   = block_gf<imtime>{triqs::mesh::imtime{beta, Fermion, p.n_tau_G}, p.gf_struct};
      F_tau   = block_gf<imtime>{triqs::mesh::imtime{beta, Fermion, p.n_tau_G}, p.gf_struct};
      G_tau() = 0;
      F_tau() = 0;
    }
    if (measure_G_l or measure_F_l) {
      G_l   = block_gf<legendre>{triqs::mesh::legendre{beta, Fermion, p.n_l}, p.gf_struct};
      F_l   = block_gf<legendre>{triqs::mesh::legendre{beta, Fermion, p.n_l}, p.gf_struct};
      G_l() = 0;
      F_l() = 0;
    }
    Z = 0;
  }

  // -------------------------------------

  void G_F_tau::accumulate(double s) {

    LOG("\n =================== MEASURE G(tau) ================ \n");

    Z += s;

    for (auto [bl_idx, det] : itertools::enumerate(wdata.dets)) {
      long N      = det.size();
      auto *g_tau = measure_G_tau ? &G_tau[bl_idx] : nullptr;
      auto *f_tau = measure_F_tau ? &F_tau[bl_idx] : nullptr;
      auto *g_l   = measure_G_l ? &G_l[bl_idx] : nullptr;
      auto *f_l   = measure_F_l ? &F_l[bl_idx] : nullptr;
      for (long id_y : range(N)) {
        auto y        = det.get_y(id_y);
        double f_fact = 0;
        if (measure_F_tau or measure_F_l) f_fact = fprefactor(bl_idx, y);
        for (long id_x : range(N)) {
          auto x    = det.get_x(id_x);
          auto Minv = det.inverse_matrix(id_y, id_x);
          // beta-periodicity is implicit in the argument, just fix the sign properly
          auto val  = (y.first >= x.first ? s : -s) * Minv;
          auto dtau = double(y.first - x.first);
          if (measure_G_tau) (*g_tau)[closest_mesh_pt(dtau)](y.second, x.second) += val;
          if (measure_F_tau) (*f_tau)[closest_mesh_pt(dtau)](y.second, x.second) += val * f_fact;
          if (measure_G_l or measure_F_l) {
            double poly_arg = 2.0 * dtau / beta - 1.0;
            auto Tn         = triqs::utility::legendre_generator();
            Tn.reset(poly_arg);
            for (auto l : (measure_G_l ? g_l->mesh() : f_l->mesh())) {
              auto const p_l = Tn.next();
              if (measure_G_l) (*g_l)[l](y.second, x.second) += val * p_l;
              if (measure_F_l) (*f_l)[l](y.second, x.second) += val * f_fact * p_l;
            }
          }
        }
      }
    }
  }

  // -------------------------------------

  void G_F_tau::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    if (measure_G_tau) {
      G_tau = mpi::all_reduce(G_tau, c);
      G_tau = G_tau / (-beta * Z * G_tau[0].mesh().delta());

      // Fix the point at zero and beta, for each block
      for (auto &g : G_tau) {
        g[0] *= 2;
        g[g.mesh().size() - 1] *= 2;
      }
      // store the result (not reused later, hence we can move it).
      results.G_tau = std::move(G_tau);
    }

    if (measure_F_tau) {
      F_tau = mpi::all_reduce(F_tau, c);
      F_tau = F_tau / (-beta * Z * F_tau[0].mesh().delta());

      for (auto &f : F_tau) {
        f[0] *= 2;
        f[f.mesh().size() - 1] *= 2;
      }
      results.F_tau = std::move(F_tau);
    }

    if (measure_G_l) {
      G_l = mpi::all_reduce(G_l, c);
      for (auto &g : G_l)
        for (auto l : g.mesh()) g[l] *= -std::sqrt(2.0 * l.index() + 1.0) / (beta * Z);
      results.G_l = std::move(G_l);
    }

    if (measure_F_l) {
      F_l = mpi::all_reduce(F_l, c);
      for (auto &f : F_l)
        for (auto l : f.mesh()) f[l] *= -std::sqrt(2.0 * l.index() + 1.0) / (beta * Z);
      results.F_l = std::move(F_l);
    }
  }

  // -------------------------------------

  double G_F_tau::fprefactor(long const &block, std::pair<tau_t, long> const &y) {
    int color    = wdata.block_to_color(block, y.second);
    double I_tau = 0;
    for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
      auto ntau = n_tau(y.first, sl); // Density to the right of y.first in sl
      if (c != color) I_tau += wdata.U(c, color) * ntau;
      if (wdata.has_Dt) {
        I_tau -= K_overlap(sl, y.first, false, wdata.Kprime, c, color);
        if (c == color) I_tau -= 2 * real(wdata.Kprime(0)(c, c));
      }
      if (wdata.has_Jperp) {
        I_tau -= 4 * real(wdata.Kprime_spin(0)(c, color)) * ntau;
        I_tau -= 2 * K_overlap(sl, y.first, false, wdata.Kprime_spin, c, color);
      }
    }
    return I_tau;
  }

} // namespace triqs_ctseg::measures
