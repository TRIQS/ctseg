// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#ifdef TRIQS_WITH_NFFT

#include "./chi4.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  chi4::chi4(params_t const &p, work_data_t const &wdata, configuration_t const & /*config*/, results_t &results)
     : wdata{wdata}, results{results}, buf_arr(wdata.gf_struct.size()) {

    beta = p.beta;

    // Meshes
    mesh_bosonic   = triqs::mesh::imfreq{p.beta, Boson, p.n_w_chi4_b};
    mesh_fermionic = triqs::mesh::imfreq{p.beta, Fermion, p.n_w_chi4_f};

    // Enlarged auxiliary fermionic mesh for M: needs to cover all arguments
    // -nu1, nu1+w, -nu2-w, nu2 where |nu1|,|nu2| < n_w_chi4_f and |w| < n_w_chi4_b
    int n_w_aux = p.n_w_chi4_f + p.n_w_chi4_b;
    aux_mesh    = triqs::mesh::imfreq{p.beta, Fermion, n_w_aux};

    // Initialize intermediate M on enlarged uniform mesh
    M = block_gf{prod{aux_mesh, aux_mesh}, wdata.gf_struct};

    // Create NFFT buffers: type1 for M (uniform grid)
    for (auto const &[bl, bl_struct] : itertools::enumerate(wdata.gf_struct)) {
      auto &[bl_name, bl_size] = bl_struct;
      buf_arr(bl) = nda::array<triqs::utility::nfft::buffer_t<2>, 2>(bl_size, bl_size);
      for (int a = 0; a < bl_size; ++a)
        for (int b = 0; b < bl_size; ++b)
          buf_arr(bl)(a, b) = triqs::utility::nfft::buffer_t<2>{
             slice_target_to_scalar(M[bl], a, b).data(), p.nfft_buf_size, p.nfft_tol};
    }

    // Initialize chi4 accumulator
    chi4_acc.resize(wdata.gf_struct.size());
    for (auto const &[bl1_idx, bl1] : itertools::enumerate(wdata.gf_struct)) {
      auto &[bl1_name, bl1_size] = bl1;
      block_names.push_back(bl1_name);
      chi4_acc[bl1_idx].resize(wdata.gf_struct.size());
      for (auto const &[bl2_idx, bl2] : itertools::enumerate(wdata.gf_struct)) {
        auto &[bl2_name, bl2_size] = bl2;
        chi4_acc[bl1_idx][bl2_idx] = gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>(
           {mesh_bosonic, mesh_fermionic, mesh_fermionic},
           make_shape(bl1_size, bl1_size, bl2_size, bl2_size));
      }
    }
  }

  // -------------------------------------

  void chi4::accumulate(double s) {

    LOG("\n ============ MEASURE CHI4 ============ \n");

    Z += s;

    // Reset M
    M() = 0;

    // Fill NFFT buffers with determinant matrix entries
    // Convention: Mw[nu1, nu2](a, b) = sum M_ji exp(i*nu1*tau_i) exp(i*nu2*tau_j)
    // Push {tau_i, tau_j} with M_ji to buffer for orbital pair (a, b)
    for (auto const &[bl, det] : itertools::enumerate(wdata.dets)) {
      long N = det.size();
      for (long i : range(N)) {
        auto [tau_i, a] = det.get_x(i); // c operator
        for (long j : range(N)) {
          auto [tau_j, b]  = det.get_y(j); // cdag operator
          auto Mij         = det.inverse_matrix(j, i);
          buf_arr(bl)(a, b).push_back({double(tau_i), double(tau_j)}, Mij);
        }
      }
    }

    // Flush all NFFT buffers -> M is now filled on enlarged uniform mesh
    for (auto &buf_block : buf_arr)
      for (auto &buf : buf_block) buf.flush();

    // Assembly: chi4[w, nu1, nu2](a,b,c,d) += s * M[-nu1, nu1+w](a,b) * M[-nu2-w, nu2](c,d)
    //           chi4[w, nu1, nu2](a,b,c,d) -= s * M[-nu1, nu2](a,d) * M[-nu2-w, nu1+w](c,b)  [same block]
    auto const &nb_blocks = wdata.gf_struct.size();
    for (auto b1 : range(nb_blocks)) {
      for (auto b2 : range(nb_blocks)) {
        auto const &block_shape = chi4_acc[b1][b2].target_shape();
        for (auto a : range(block_shape[0])) {
          for (auto b : range(block_shape[1])) {
            for (auto c : range(block_shape[2])) {
              for (auto d : range(block_shape[3])) {
                for (auto const &nu1 : mesh_fermionic) {
                  for (auto const &nu2 : mesh_fermionic) {
                    for (auto const &w : mesh_bosonic) {
                      // Direct term
                      chi4_acc[b1][b2][w, nu1, nu2](a, b, c, d) +=
                         s * M[b1][-nu1, nu1 + w](a, b) * M[b2][-nu2 - w, nu2.value()](c, d);
                      // Exchange term (same block only)
                      if (b1 == b2)
                        chi4_acc[b1][b2][w, nu1, nu2](a, b, c, d) -=
                           s * M[b1][-nu1, nu2.value()](a, d) * M[b2][-nu2 - w, nu1 + w](c, b);
                    } // w
                  } // nu2
                } // nu1
              } // d
            } // c
          } // b
        } // a
      } // b2
    } // b1
  }

  // -------------------------------------

  void chi4::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    for (auto b1 : range(wdata.gf_struct.size()))
      for (auto b2 : range(wdata.gf_struct.size())) {
        chi4_acc[b1][b2] = mpi::all_reduce(chi4_acc[b1][b2], c);
        chi4_acc[b1][b2] = chi4_acc[b1][b2] / (Z * beta);
      }

    auto chi4_block = make_block2_gf(block_names, block_names, chi4_acc);
    results.chi4    = std::move(chi4_block);
  }

} // namespace triqs_ctseg::measures

#endif // TRIQS_WITH_NFFT
