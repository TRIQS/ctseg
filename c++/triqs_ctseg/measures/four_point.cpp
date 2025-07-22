// Copyright (c) 2024--present, The Simons Foundation
// Copyright (c) 2024--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./four_point.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  four_point::four_point(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta           = p.beta;
    n_w_bosonic    = p.n_w_b_vertex;
    n_w_fermionic  = p.n_w_f_vertex;
    measure_g2w    = p.measure_g2w;
    measure_g3w    = p.measure_g3w;
    mesh_bosonic   = triqs::mesh::imfreq(beta, Boson  , n_w_bosonic  );
    mesh_fermionic = triqs::mesh::imfreq(beta, Fermion, n_w_fermionic);

    if (measure_g2w) g2w.resize(wdata.gf_struct.size());
    if (measure_g3w) g3w.resize(wdata.gf_struct.size());
    for (auto const &[bl1_idx, bl1] : itertools::enumerate(wdata.gf_struct)) {
      auto &[bl1_name, bl1_size] = bl1;
      block_names.push_back(bl1_name);
      if (measure_g2w) g2w[bl1_idx].resize(wdata.gf_struct.size());
      if (measure_g3w) g3w[bl1_idx].resize(wdata.gf_struct.size());
      for (auto const &[bl2_idx, bl2] : itertools::enumerate(wdata.gf_struct)) {
        auto &[bl2_name, bl2_size] = bl2;
        if (measure_g2w) {
          g2w[bl1_idx][bl2_idx] = gf<prod<imfreq, imfreq>, tensor_valued<4>>(
              { mesh_bosonic , mesh_fermionic },
              make_shape(bl1_size, bl1_size, bl2_size, bl2_size));
          g2w[bl1_idx][bl2_idx]() = 0;
        }
        if (measure_g3w) {
          g3w[bl1_idx][bl2_idx] = gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>(
              { mesh_bosonic , mesh_fermionic , mesh_fermionic },
              make_shape(bl1_size, bl1_size, bl2_size, bl2_size));
          g3w[bl1_idx][bl2_idx]() = 0;
        }
      }
    }

  }

  // -------------------------------------

  void four_point::accumulate(double s) {

    LOG("\n ============ MEASURE FOUR-POINT CORRELATION FUNCTION  ============ \n");

    /// Measure the four-point correlation function
    /**
    * The four-point correlation function is defined as:
    
    $$g_{abcd}^{(4)\,ph}(\omega, \nu, \nu') = \frac{1}{\beta} \int_0^\beta d\tau_2\, d\tau_3\, d\tau_4\, 
    \exp\left[i\omega(\tau_2 - \tau_3) + i\nu(\tau_2 - \tau_1) + i\nu'(\tau_4 - \tau_3)\right] \left\langle c_a^\dagger(\tau_1) c_b(\tau_2) c_c^\dagger(\tau_3) c_d(\tau_4) \right\rangle
    = \left\langle c_a^\dagger(-nu) c_b(nu + omega) c_c^\dagger(-nu' - omega) c_d(nu') \right\rangle
    $$
    
    * The number of fermionic (bosonic) frequencies is specified through the
    parameters ``n_w_f_vertex`` (``n_w_b_vertex``).
    */

    Z += s;

    auto Mw = compute_Mw();
    auto const &nb_blocks = wdata.gf_struct.size();

    if (measure_g2w) {
      auto nw = compute_nw();
      for (auto const &b1 : range(nb_blocks)) {
        for (auto const &b2 : range(nb_blocks)) {
          auto const &block_shape = g2w[b1][b2].target_shape();
          for (auto const &c : range(block_shape[2])) {
            auto col = wdata.block_to_color(b2, c);
            for (auto const &a : range(block_shape[0])) {
              for (auto const &b : range(block_shape[1])) {
                for (auto const &w : mesh_bosonic) {
                  for (auto const &nu1 : mesh_fermionic) {
                    g2w[b1][b2][w, nu1](a, b, c, c) -= s * Mw[b1][-nu1, nu1 + w](a, b) * nw[col][-w];
                  } // nu1
                } // w
              } // b
            } // a
          } // c
        } // b2
      } // b1
    } // measure_g2w

    if (measure_g3w) {
      for (auto const &b1 : range(nb_blocks)) {
        for (auto const &b2 : range(nb_blocks)) {
          auto const &block_shape = g3w[b1][b2].target_shape();
          for (auto const &a : range(block_shape[0])) {
            for (auto const &b : range(block_shape[1])) {
              for (auto const &c : range(block_shape[2])) {
                for (auto const &d : range(block_shape[3])) {
                  for (auto const &nu1 : mesh_fermionic) {
                    for (auto const &nu2 : mesh_fermionic) {
                      for (auto const &w : mesh_bosonic) {
                        g3w[b1][b2][w, nu1, nu2](a, b, c, d) += s * 
                        Mw[b1][-nu1, nu1 + w](a, b) * Mw[b2][-nu2 - w, nu2.value()](c, d);
                        if (b1 == b2)
                          g3w[b1][b2][w, nu1, nu2](a, b, c, d) -= s *
                          Mw[b1][-nu1, nu2.value()](a, d) * Mw[b2][-nu2 - w, nu1 + w](c, b);
                      } // w
                    } // nu2
                  } // nu1
                } // d
              } // c
            } // b
          } // a
        } // b2
      } // b1
    } // measure_g3w

  }

  // -------------------------------------

  void four_point::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    if (measure_g2w) g2w = mpi::all_reduce(g2w, c);
    if (measure_g3w) g3w = mpi::all_reduce(g3w, c);
    for (auto const &b1 : range(wdata.gf_struct.size())) {
      for (auto const &b2 : range(wdata.gf_struct.size())) {
        if (measure_g2w) g2w[b1][b2] = g2w[b1][b2] / (Z * beta);
        if (measure_g3w) g3w[b1][b2] = g3w[b1][b2] / (Z * beta);
      }
    }

    if (measure_g2w) {
      g2w_block = make_block2_gf(block_names, block_names, g2w);
      results.g2w = std::move(g2w_block);
    }
    if (measure_g3w) {
      g3w_block = make_block2_gf(block_names, block_names, g3w);
      results.g3w = std::move(g3w_block);
    }

  }

  // -------------------------------------

  block_gf<prod<imfreq, imfreq>> four_point::compute_Mw() {

    // Mw[bl][nu1, nu2](a, b) = < c^dagger_{bl,a} (nu1) c_{bl,b} (nu2) >
    
    int n_w_aux = n_w_fermionic + n_w_bosonic - 1;
    auto aux_mesh = triqs::mesh::imfreq(beta, Fermion, n_w_aux);
    auto Mw = block_gf(prod(aux_mesh, aux_mesh), wdata.gf_struct);
    Mw() = 0;

    auto w0 = aux_mesh[0].value();
    auto dw = aux_mesh[1].value() - aux_mesh[0].value();

    for (auto const &[bl, det] : itertools::enumerate(wdata.dets)) {
      long N = det.size();
      for (long i : range(N)) {
        auto [tau_i, a] = det.get_x(i);
        for (long j : range(N)) {
          auto [tau_j, b] = det.get_y(j);
          auto Mij = det.inverse_matrix(j, i);
          auto exp_i = std::exp(w0 * double(tau_i));
          auto exp_i_dw = std::exp(dw * double(tau_i));
          auto exp_j0 = std::exp(w0 * double(tau_j));
          auto exp_j_dw = std::exp(dw * double(tau_j));
          auto exp_j = exp_j0; 
          for (auto const &nu1 : aux_mesh) {
            exp_j = exp_j0;
            auto expMij = exp_i * Mij; 
            for (auto const &nu2 : aux_mesh) {
              Mw[bl][nu1, nu2](a, b) += expMij * exp_j;
              exp_j = exp_j * exp_j_dw;
            } // nu2
            exp_i = exp_i * exp_i_dw;
          } // nu1
        } // tau_j
      } // tau_i
    } // bl
    return Mw;

  }

  // -------------------------------------

  std::vector<gf<imfreq, scalar_valued>> four_point::compute_nw() {

    /* 
    $$ n(col)[\omega] = \sum_{\text{segments}} 
    \int_{\tau_\mathrm{start}}^{\tau_\mathrm{end}}e^{i\omega\tau}\mathrm{d}\tau $$
    */

    std::vector<gf<imfreq, scalar_valued>> nw;
    nw.resize(wdata.n_color);
    for (auto const &c : range(wdata.n_color)) {
      nw[c] = gf<imfreq, scalar_valued>(mesh_bosonic);
      nw[c]() = 0;

      for (auto const &s: config.seglists[c]) {

        double tau_c = double(s.tau_c);
        double tau_cdag = double(s.tau_cdag);

        // Zero frequency: Add up the all the segment length
        if (!is_cyclic(s))
          nw[c][0] += tau_c - tau_cdag;
        else
          nw[c][0] += beta - tau_cdag + tau_c;

        // Compute remaining frequencies
        for (auto const &w : mesh_bosonic) {
          if (w.n == 0) continue;
          nw[c][w] += (std::exp(w.value() * tau_c) - std::exp(w.value() * tau_cdag)) / w.value();
        } // w
        
      } // s
    } // c
    return nw;

  }

} // namespace triqs_ctseg::measures
