// Copyright (c) 2024--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE.txt in the root of this distribution for details.

#include "./four_point.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  four_point::four_point(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta           = p.beta;
    measure_g3w    = p.measure_g3w;
    measure_f3w    = p.measure_f3w;

    auto g3w_vec = [&]() {
      std::vector<gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>> green_v;
      for (auto const &[bl1_name, bl1_size] : wdata.gf_struct)
        for (auto const &[bl2_name, bl2_size] : wdata.gf_struct)
          green_v.emplace_back(gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>(
              {{beta, Boson,   p.n_w_b_vertex, imfreq::option::all_frequencies},
               {beta, Fermion, p.n_w_f_vertex, imfreq::option::all_frequencies},
               {beta, Fermion, p.n_w_f_vertex, imfreq::option::all_frequencies}},
              make_shape(bl1_size, bl1_size, bl2_size, bl2_size)));
      return green_v;
    };

    for (auto const &[bl_name, bl_size] : wdata.gf_struct) block_names.push_back(bl_name);
    auto bosonic_block_names = std::vector<std::string>{};
    for (auto const &str1 : block_names)
      for (auto const &str2 : block_names)
        bosonic_block_names.push_back(str1 + "|" + str2);

    g3w = make_block_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>(bosonic_block_names, g3w_vec());
    f3w = make_block_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>(bosonic_block_names, g3w_vec());
    g3w() = 0;
    f3w() = 0;

    n_w_fermionic = std::get<1>(g3w[0].mesh().components()).last_index() + 1;
    n_w_bosonic   = std::get<0>(g3w[0].mesh().components()).last_index() + 1;

    w_ini = (2 * (- n_w_fermionic - n_w_bosonic + 1) + 1) * M_PI / beta;
    w_inc = 2 * M_PI / beta;

  }

  // -------------------------------------

  void four_point::accumulate(double s) {

    LOG("\n ============ MEASURE FOUR-POINT CORRELATION FUNCTION  ============ \n");

    /// Measure the four-point correlation function
    /**
    *The four-point correlation function is defined as:
    *
    *$$\chi^{\sigma\sigma'}_{abcd}(i\omega, i\omega',i\Omega) = G^{2,\sigma,
    *\sigma'}_{abcd}(i\omega, i\omega',i\Omega) = \langle c_{a\sigma}(i\omega)
    *c^\dagger_{b\sigma}(i\omega+i\Omega) c_{c\sigma'}(i\omega'+i\Omega)
    *c^\dagger_{d\sigma'}(i\omega') \rangle$$
    *
    * Its improved estimator is the Fourier transform of
    *
    *$$F^{3,\sigma\sigma'}_{abcd}(\tau,\tau',\tau'') = \int \mathrm{d}\bar{\tau}
    *\sum_{e\bar{\sigma}} \mathcal{U}^{\sigma\bar{\sigma}}_{ae}(\bar{\tau}-\tau)
    *\langle c_{a\sigma}(\tau) c^\dagger_{b\sigma}(\tau') c_{c\sigma'}(\tau'')
    *c^\dagger_{d\sigma'}(0) \rangle$$
    *
    * The vertex corresponding to this correlation function is evaluated separately.
    */

    Z += s;

    auto const mesh_fermionic = std::get<1>(g3w[0].mesh());
    auto const mesh_bosonic   = std::get<0>(g3w[0].mesh());

    Mw_vector = compute_Mw(false);

    if (measure_g3w) {
      for (long bl = 0; bl < g3w.size(); bl++) { // bl : 'upup', 'updn', ...
        long b1 = bl / wdata.gf_struct.size();
        long b2 = bl % wdata.gf_struct.size();
        for (int a = 0; a < g3w[bl].target_shape()[0]; a++) {
          for (int b = 0; b < g3w[bl].target_shape()[1]; b++) {
            for (int c = 0; c < g3w[bl].target_shape()[2]; c++) {
              for (int d = 0; d < g3w[bl].target_shape()[3]; d++) {
                for (int n1 = -n_w_fermionic; n1 < n_w_fermionic; n1++) {
                  for (int n4 = -n_w_fermionic; n4 < n_w_fermionic; n4++) {
                    for (int m = -n_w_bosonic + 1; m < n_w_bosonic; m++) {
                      int n2 = n1 + m;
                      int n3 = n4 + m;
                      // This structure is ugly. Need someone who familiar with TRIQS to prune this part.
                      auto freq_1 = mesh_fermionic[n1 + n_w_fermionic];
                      auto freq_2 = mesh_fermionic[n4 + n_w_fermionic];
                      auto freq_3 = mesh_bosonic[m + n_w_bosonic - 1];
                      g3w[bl][freq_3, freq_1, freq_2](a, b, c, d) += s * Mw(b1, a, b, n1, n2) * Mw(b2, c, d, n3, n4);
                      if (b1 == b2)
                        g3w[bl][freq_3, freq_1, freq_2](a, b, c, d) -= s * Mw(b1, a, d, n1, n4) * Mw(b2, c, b, n3, n2);
                    } // m
                  } // n4
                } // n1
              } // d
            } // c
          } // b
        } // a
      } // bl
    } // measure_g3w

    if (measure_f3w) {
      nMw_vector = compute_Mw(true);
      for (long bl = 0; bl < f3w.size(); bl++) { // bl : 'upup', 'updn', ...
        long b1 = bl / wdata.gf_struct.size();
        long b2 = bl % wdata.gf_struct.size();
        for (int a = 0; a < f3w[bl].target_shape()[0]; a++) {
          for (int b = 0; b < f3w[bl].target_shape()[1]; b++) {
            for (int c = 0; c < f3w[bl].target_shape()[2]; c++) {
              for (int d = 0; d < f3w[bl].target_shape()[3]; d++) {
                for (int n1 = -n_w_fermionic; n1 < n_w_fermionic; n1++) {
                  for (int n4 = -n_w_fermionic; n4 < n_w_fermionic; n4++) {
                    for (int m = -n_w_bosonic + 1; m < n_w_bosonic; m++) {
                      int n2 = n1 + m;
                      int n3 = n4 + m;
                      // Please prune this
                      auto freq_1 = mesh_fermionic[n1 + n_w_fermionic];
                      auto freq_2 = mesh_fermionic[n4 + n_w_fermionic];
                      auto freq_3 = mesh_bosonic[m + n_w_bosonic - 1];
                      f3w[bl][freq_3, freq_1, freq_2](a, b, c, d) += s * nMw(b1, a, b, n1, n2) * Mw(b2, c, d, n3, n4);
                      if (b1 == b2)
                        f3w[bl][freq_3, freq_1, freq_2](a, b, c, d) -= s * nMw(b1, a, d, n1, n4) * Mw(b2, c, b, n3, n2);
                    } // m
                  } // n4
                } // n1
              } // d
            } // c
          } // b
        } // a
      } // bl
    } // measure_f3w
    
  }

  // -------------------------------------

  void four_point::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);
    if (measure_g3w) {
      g3w = mpi::all_reduce(g3w, c);
      g3w = g3w / (Z * beta);

      std::vector<std::vector<gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>>> g3w_vec(wdata.gf_struct.size(), std::vector<gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>>(wdata.gf_struct.size()));
      for (int b1 : range(wdata.gf_struct.size())) {
        for (int b2 : range(wdata.gf_struct.size())) {
          g3w_vec[b1][b2] = g3w[b1 * wdata.gf_struct.size() + b2];
        }
      }

      auto g3w_block = make_block2_gf(block_names, block_names, g3w_vec);
      results.g3w = std::move(g3w_block);
    }

    if (measure_f3w) {
      f3w = mpi::all_reduce(f3w, c);
      f3w = f3w / (Z * beta);

      std::vector<std::vector<gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>>> f3w_vec(wdata.gf_struct.size(), std::vector<gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>>(wdata.gf_struct.size()));
      for (int b1 : range(wdata.gf_struct.size())) {
        for (int b2 : range(wdata.gf_struct.size())) {
          f3w_vec[b1][b2] = f3w[b1 * wdata.gf_struct.size() + b2];
        }
      }

      auto f3w_block = make_block2_gf(block_names, block_names, f3w_vec);
      results.f3w = std::move(f3w_block);
    }

  }

  // -------------------------------------

  double four_point::fprefactor(long const &block, std::pair<tau_t, long> const &y) {

    // Copied from G_F_tau.cpp
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

  // -------------------------------------

  std::vector<array<dcomplex, 4>> four_point::compute_Mw(bool is_nMw) {

    int n_w_aux = 2 * (n_w_fermionic + n_w_bosonic - 1) > 0 ? 2 * (n_w_fermionic + n_w_bosonic - 1) : 0;
    std::vector<array<dcomplex, 4>> result;
    result.resize(wdata.gf_struct.size());

    for (auto const &[bl, bl_pair] : itertools::enumerate(wdata.gf_struct)) {
      auto const &[bl_name, bl_size] = bl_pair;
      result[bl].resize(make_shape(bl_size, bl_size, n_w_aux, n_w_aux));
      result[bl]() = 0;
    }

    for (auto const &[bl, det] : itertools::enumerate(wdata.dets)) {
      long N = det.size();
      y_exp_ini.resize(N);
      y_exp_inc.resize(N);
      x_exp_ini.resize(N);
      x_exp_inc.resize(N);
      y_inner_index.resize(N);
      x_inner_index.resize(N);

      for (long id : range(N)) {
        auto y = det.get_y(id);
        auto x = det.get_x(id);
        y_exp_ini(id) = std::exp(dcomplex(0, w_ini * double(std::get<0>(y))));
        y_exp_inc(id) = std::exp(dcomplex(0, w_inc * double(std::get<0>(y))));
        x_exp_ini(id) = std::exp(dcomplex(0, -w_ini * double(std::get<0>(x))));
        x_exp_inc(id) = std::exp(dcomplex(0, -w_inc * double(std::get<0>(x))));
        y_inner_index(id) = std::get<1>(y);
        x_inner_index(id) = std::get<1>(x);
      }

      for (long id_y : range(N)) {
        auto y = det.get_y(id_y);
        int yj = y_inner_index(id_y);
        double f_fact = is_nMw ? fprefactor(bl, y) : 1.0;

        for (long id_x : range(N)) {
          int xi = x_inner_index(id_x);
          dcomplex y_exp = y_exp_ini(id_y);
          dcomplex x_exp = x_exp_ini(id_x);
          auto Minv = det.inverse_matrix(id_y, id_x);

          for (int n_1 : range(n_w_aux)) {
            for (int n_2 : range(n_w_aux)) {
              auto val = Minv * y_exp * x_exp;
              result[bl](yj, xi, n_1, n_2) += val * f_fact;
              x_exp *= x_exp_inc(id_x);
            }
            x_exp = x_exp_ini(id_x);
            y_exp *= y_exp_inc(id_y);
          }
        }
      }
    }
    return result;

  }

} // namespace triqs_ctseg::measures
