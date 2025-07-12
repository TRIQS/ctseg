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
    mesh_bosonic   = triqs::mesh::imfreq(beta, Boson  , n_w_bosonic  );
    mesh_fermionic = triqs::mesh::imfreq(beta, Fermion, n_w_fermionic);

    g3w.resize(wdata.gf_struct.size());
    for (auto const &[bl1_idx, bl1] : itertools::enumerate(wdata.gf_struct)) {
      auto &[bl1_name, bl1_size] = bl1;
      block_names.push_back(bl1_name);
      g3w[bl1_idx].resize(wdata.gf_struct.size());
      for (auto const &[bl2_idx, bl2] : itertools::enumerate(wdata.gf_struct)) {
        auto &[bl2_name, bl2_size] = bl2;
        g3w[bl1_idx][bl2_idx] = gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>(
            { mesh_bosonic , mesh_fermionic , mesh_fermionic },
            make_shape(bl1_size, bl1_size, bl2_size, bl2_size));
        g3w[bl1_idx][bl2_idx]() = 0;
      }
    }

  }

  // -------------------------------------

  void four_point::accumulate(double s) {

    LOG("\n ============ MEASURE FOUR-POINT CORRELATION FUNCTION  ============ \n");

    /// Measure the four-point correlation function
    /**
    * The four-point correlation function is defined as:
    
    $$\chi^{\sigma\sigma'}_{abcd}(i\omega, i\omega',i\Omega) = G^{2,\sigma,
    \sigma'}_{abcd}(i\omega, i\omega',i\Omega) = \langle c_{a\sigma}(i\omega)
    c^\dagger_{b\sigma}(i\omega+i\Omega) c_{c\sigma'}(i\omega'+i\Omega)
    c^\dagger_{d\sigma'}(i\omega') \rangle$$
    
    * The number of fermionic (bosonic) frequencies is specified through the
    parameters ``n_w_f_vertex`` (``n_w_b_vertex``).
    */

    Z += s;

    auto Mw = compute_Mw();

    for (auto const &b1 : range(wdata.gf_struct.size())) {
      for (auto const &b2 : range(wdata.gf_struct.size())) {
        for (auto const &a : range(g3w[b1][b2].target_shape()[0])) {
          for (auto const &b : range(g3w[b1][b2].target_shape()[1])) {
            for (auto const &c : range(g3w[b1][b2].target_shape()[2])) {
              for (auto const &d : range(g3w[b1][b2].target_shape()[3])) {
                for (auto const &n1 : mesh_fermionic) {
                  for (auto const &n4 : mesh_fermionic) {
                    for (auto const &m : mesh_bosonic) {
                      auto n2 = n1 + m;
                      auto n3 = n4 + m;
                      g3w[b1][b2][m, n1, n4](a, b, c, d) += s * 
                      Mw[b1](a, b, n1.n + n_w_fermionic + n_w_bosonic - 1, n2.n + n_w_fermionic + n_w_bosonic - 1) * 
                      Mw[b2](c, d, n3.n + n_w_fermionic + n_w_bosonic - 1, n4.n + n_w_fermionic + n_w_bosonic - 1);
                    } // m
                  } // n4
                } // n1
              } // d
            } // c
          } // b
        } // a
      } // b2
    } // b1
    
  }

  // -------------------------------------

  void four_point::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);
    g3w = mpi::all_reduce(g3w, c);
    for (auto const &b1 : range(wdata.gf_struct.size())) {
      for (auto const &b2 : range(wdata.gf_struct.size())) {
        g3w[b1][b2] = g3w[b1][b2] / (Z * beta);
      }
    }

    g3w_block = make_block2_gf(block_names, block_names, g3w);
    results.g3w = std::move(g3w_block);

  }

  // -------------------------------------

  std::vector<array<dcomplex, 4>> four_point::compute_Mw() {

    // Mw(a, b, c ,d) = < c^dagger_a (nu[c]) c_b (nu[d]) >

    std::vector<array<dcomplex, 4>> Mw(wdata.gf_struct.size());
    int n_w_aux = n_w_fermionic + n_w_bosonic - 1;
    auto aux_mesh = triqs::mesh::imfreq(beta, Fermion, n_w_aux);
    auto w0 = aux_mesh[0].value();
    auto dw = aux_mesh[1].value() - aux_mesh[0].value();

    for (auto const &[bl, bl_info] : itertools::enumerate(wdata.gf_struct)) {
      auto const &[bl_name, bl_size] = bl_info;
      Mw[bl].resize(make_shape(bl_size, bl_size, n_w_aux, n_w_aux));
      Mw[bl]() = 0;
    }

    for (auto const &[bl, det] : itertools::enumerate(wdata.dets)) {
      long N = det.size();
      for (long i : range(N)) {
        auto [tau_i, a] = det.get_x(i);
        for (long j : range(N)) {
          auto [tau_j, b] = det.get_y(j);
          auto Mij = det.inverse_matrix(i, j);
          auto exp_i = std::exp(w0 * double(tau_i));
          auto exp_i_dw = std::exp(dw * double(tau_i));
          auto exp_j0 = std::exp(w0 * double(tau_j));
          auto exp_j_dw = std::exp(dw * double(tau_j));
          auto exp_j = exp_j0; 
          for (int n : range(aux_mesh.size())) {
            exp_j = exp_j0;
            auto expMij = exp_i * Mij; 
            for (int m : range(aux_mesh.size())) {
              Mw[bl](a, b, n, m) += expMij * exp_j;
              exp_j = exp_j * exp_j_dw;
            }
            exp_i = exp_i * exp_i_dw;
          }
        }
      }
    }
    return Mw;
  }

} // namespace triqs_ctseg::measures
