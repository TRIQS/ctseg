// Copyright (c) 2024--present, The Simons Foundation
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#ifdef TRIQS_WITH_NFFT

#include "./chi3.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  chi3::chi3(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results}, buf_arr(wdata.gf_struct.size()) {

    beta = p.beta;

    // Construct DLR2D Matsubara mesh (PH channel)
    chi3_mesh = triqs::mesh::dlr2d_imfreq{p.beta, p.dlr_wmax, p.dlr_eps, triqs::mesh::PH, p.dlr2d_compress_grid};

    // Initialize intermediate M on DLR2D mesh
    M = block_gf{chi3_mesh, wdata.gf_struct};

    // Build target frequencies for type3 NFFT
    // Convention from four_point.cpp compute_Mw():
    //   Mw[nu1, nu2](a, b) = sum_{i,j} M_ji * exp(i*nu1*tau_i) * exp(i*nu2*tau_j)
    //   where tau_i from get_x (c operator, orbital a), tau_j from get_y (cdag operator, orbital b)
    //
    // chi3 uses M[-nu1, nu2] where (nu1, nu2) are DLR2D mesh points.
    // NFFT dim0 pairs with tau_i, dim1 pairs with tau_j.
    // Target: {-nu1, nu2}
    target_mf.reserve(chi3_mesh.size());
    for (auto mp : chi3_mesh) {
      auto [nu1, nu2] = mp.value();
      target_mf.push_back({-nu1, nu2});
    }

    // Create nfft buffers: type3 for M (DLR2D targets)
    for (auto const &[bl, bl_struct] : itertools::enumerate(wdata.gf_struct)) {
      auto &[bl_name, bl_size] = bl_struct;
      buf_arr(bl) = nda::array<triqs::utility::nfft::buffer_t<2>, 2>(bl_size, bl_size);
      for (int a = 0; a < bl_size; ++a)
        for (int b = 0; b < bl_size; ++b)
          buf_arr(bl)(a, b) = triqs::utility::nfft::buffer_t<2>{
             slice_target_to_scalar(M[bl], a, b).data(), target_mf, p.nfft_buf_size, p.nfft_tol};
    }

    // Bosonic mesh for nw: need to cover all differences nu1 - nu2 from DLR2D mesh
    nw_mesh = triqs::mesh::imfreq{p.beta, Boson, chi3_mesh.max_n() + 1};

    // Initialize nw per color
    nw.resize(wdata.n_color, gf<imfreq, scalar_valued>{nw_mesh});

    // Initialize chi3 accumulator
    chi3_acc.resize(wdata.gf_struct.size());
    for (auto const &[bl1_idx, bl1] : itertools::enumerate(wdata.gf_struct)) {
      auto &[bl1_name, bl1_size] = bl1;
      block_names.push_back(bl1_name);
      chi3_acc[bl1_idx].resize(wdata.gf_struct.size());
      for (auto const &[bl2_idx, bl2] : itertools::enumerate(wdata.gf_struct)) {
        auto &[bl2_name, bl2_size] = bl2;
        chi3_acc[bl1_idx][bl2_idx] = gf<dlr2d_imfreq, tensor_valued<4>>(
           chi3_mesh, make_shape(bl1_size, bl1_size, bl2_size, bl2_size));
      }
    }
  }

  // -------------------------------------

  void chi3::accumulate(double s) {

    LOG("\n ============ MEASURE CHI3 ============ \n");

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

    // Flush all NFFT buffers -> M is now filled on DLR2D mesh
    for (auto &buf_block : buf_arr)
      for (auto &buf : buf_block) buf.flush();

    // Compute nw from segment lists (analytic integration)
    // nw[col][w] = sum_segments integral_{tau_cdag}^{tau_c} exp(i*w*tau) dtau
    for (auto &n : nw) n() = 0;

    for (int col = 0; col < wdata.n_color; ++col) {
      nw[col][nw_mesh(0)] = -beta; // Convention: nw[omega = 0] = beta * (density - 1), matching g2w

      for (auto const &seg : config.seglists[col]) {
        double tau_c    = double(seg.tau_c);
        double tau_cdag = double(seg.tau_cdag);

        // Zero frequency: add segment lengths -> nw(0) = beta * (n_col - 1)
        nw[col][nw_mesh(0)] += double(seg.length());

        // Non-zero frequencies
        for (auto const &w : nw_mesh) {
          if (w.n == 0) continue;
          nw[col][w] += (std::exp(w.value() * tau_c) - std::exp(w.value() * tau_cdag)) / w.value();
        }
      }
    }

    // Assembly: chi3[nu1, nu2](a,b,c,c) -= s * M[-nu1, nu2](a,b) * nw[col][nu1 - nu2]
    auto nb_blocks = wdata.gf_struct.size();
    for (auto b1 : range(nb_blocks)) {
      for (auto b2 : range(nb_blocks)) {
        auto const &block_shape = chi3_acc[b1][b2].target_shape();
        for (auto c : range(block_shape[2])) {
          auto col = wdata.block_to_color(b2, c);
          for (auto mp : chi3_mesh) {
            auto [nu1, nu2] = mp.value();
            auto Omega      = nu2 - nu1; // bosonic frequency
            auto nw_val     = nw[col][-Omega];
            for (auto a : range(block_shape[0]))
              for (auto b : range(block_shape[1]))
                chi3_acc[b1][b2][mp](a, b, c, c) -= s * M[b1][mp](a, b) * nw_val;
          }
        }
      }
    }
  }

  // -------------------------------------

  void chi3::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    for (auto b1 : range(wdata.gf_struct.size()))
      for (auto b2 : range(wdata.gf_struct.size())) {
        chi3_acc[b1][b2] = mpi::all_reduce(chi3_acc[b1][b2], c);
        chi3_acc[b1][b2] = chi3_acc[b1][b2] / (Z * beta);
      }

    auto chi3_block = make_block2_gf(block_names, block_names, chi3_acc);
    results.chi3    = std::move(chi3_block);
  }

} // namespace triqs_ctseg::measures

#endif // TRIQS_WITH_NFFT
