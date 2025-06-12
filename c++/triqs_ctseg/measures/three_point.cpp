// Copyright (c) 2025--present, The Simons Foundation
// Copyright (c) 2025--present, Max Planck Institute for Polymer Research, Mainz, Germany
// Copyright (c) 2025--present, EPFL, Switzerland
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./three_point.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  three_point::three_point(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta           = p.beta;
    measure_g2w    = p.measure_g2w;
    measure_f2w    = p.measure_f2w;

    auto g2w_vec = [&]() {
      std::vector<gf<prod<imfreq, imfreq>, tensor_valued<3>>> green_v;
      for (auto const &[bl1_name, bl1_size] : wdata.gf_struct)
        for (auto const &[bl2_name, bl2_size] : wdata.gf_struct)
          green_v.emplace_back(gf<prod<imfreq, imfreq>, tensor_valued<3>>(
              {{beta, Boson,   p.n_w_b_vertex, imfreq::option::all_frequencies},
               {beta, Fermion, p.n_w_f_vertex, imfreq::option::all_frequencies}},
              make_shape(bl1_size, bl2_size, bl2_size)));
      return green_v;
    };

    for (auto const &[bl_name, bl_size] : wdata.gf_struct) block_names.push_back(bl_name);
    auto bosonic_block_names = std::vector<std::string>{};
    for (auto const &str1 : block_names)
      for (auto const &str2 : block_names)
        bosonic_block_names.push_back(str1 + "|" + str2);
    
    g2w = make_block_gf<prod<imfreq, imfreq>, tensor_valued<3>>(bosonic_block_names, g2w_vec());
    f2w = make_block_gf<prod<imfreq, imfreq>, tensor_valued<3>>(bosonic_block_names, g2w_vec());
    g2w() = 0;
    f2w() = 0;

    n_w_fermionic = std::get<1>(g2w[0].mesh().components()).last_index() + 1;
    n_w_bosonic   = std::get<0>(g2w[0].mesh().components()).last_index() + 1;

    w_ini = (2 * (- n_w_fermionic - n_w_bosonic + 1) + 1) * M_PI / beta;
    w_inc = 2 * M_PI / beta;

  }

  // -------------------------------------

  void three_point::accumulate(double s) {

    LOG("\n ============ MEASURE THREE-POINT CORRELATION FUNCTION  ============ \n");

    /// Measure the three-point (triangular) correlation function
    /**
    * The three-point correlation function is defined as:

    $$X_{ab}(i\omega_n,i\nu_m) = \int_0^\beta d\tau \int_0^\beta d\tau'
    e^{i\omega_n\tau} e^{i\nu_m\tau'} X_{ab}(\tau,\tau')$$

    with $X=G^{2},F^{2}$ defined as:

    $$G^{2,\sigma\sigma'}_{abc}(\tau,\tau') = -\langle T_\tau
    c_{a\sigma}(\tau)c_{b\sigma}^\dagger(0)n_{c\sigma'}(\tau') \rangle$$

    * Its improved estimator is the Fourier transform of

    $$F^{2,\sigma\sigma'}_{abc}(\tau,\tau') = -\int_0^\beta d\tilde{\tau}
    \sum_{d\bar{\sigma}} \langle T_\tau n_{d\bar{\sigma}}(\tilde{\tau})
    \mathcal{U}^{\sigma\bar{\sigma}}_{ad}(\tilde{\tau}-\tau)
    c_{a\sigma}(\tau)c_{b\sigma}^\dagger(0)n_{c\sigma'}(\tau') \rangle$$

    * The number of fermionic (bosonic) frequencies is specified through the
    parameters ``n_w_f_vertex`` (``n_w_b_vertex``).
    */

    Z += s;

    auto const mesh_fermionic = std::get<1>(g2w[0].mesh());
    auto const mesh_bosonic   = std::get<0>(g2w[0].mesh());

    nw_vector = compute_nw();

    if (measure_g2w) {
      Mw_vector = compute_Mw(false);
      for (long bl = 0; bl < g2w.size(); bl++) { // bl : 'upup', 'updn', ...
        long b1 = bl / wdata.gf_struct.size();
        long b2 = bl % wdata.gf_struct.size();
        for (int c = 0; c < g2w[bl].target_shape()[2]; c++) {
          auto col = wdata.block_to_color(b2, c);
          for (int a = 0; a < g2w[bl].target_shape()[0]; a++) {
            for (int b = 0; b < g2w[bl].target_shape()[1]; b++) {
              for (int m = -n_w_bosonic + 1; m < n_w_bosonic; m++) {
                for (int n1 = -n_w_fermionic; n1 < n_w_fermionic; n1++) {
                  int n2 = n1 + m; // set remaining frequency
                  // This structure is ugly. Need someone who familiar with TRIQS to prune this part.
                  auto freq_1 = mesh_bosonic[m + n_w_bosonic - 1];
                  auto freq_2 = mesh_fermionic[n1 + n_w_fermionic];
                  g2w[bl][freq_1, freq_2](a, b, c) -= s * Mw(b1, a, b, n1, n2) * nw(col, m);
                } // n1
              } // m
            } // b
          } // a
        } // c
      } // bl
    }

    if (measure_f2w) {
      nMw_vector = compute_Mw(true);
      for (long bl = 0; bl < f2w.size(); bl++) { // bl : 'upup', 'updn', ...
        long b1 = bl / wdata.gf_struct.size();
        long b2 = bl % wdata.gf_struct.size();
        for (int c = 0; c < f2w[bl].target_shape()[2]; c++) {
          auto col = wdata.block_to_color(b2, c);
          for (int a = 0; a < f2w[bl].target_shape()[0]; a++) {
            for (int b = 0; b < f2w[bl].target_shape()[1]; b++) {
              for (int m = -n_w_bosonic + 1; m < n_w_bosonic; m++) {
                for (int n1 = -n_w_fermionic; n1 < n_w_fermionic; n1++) {
                  int n2 = n1 + m; // set remaining frequency
                  // This structure is ugly. Need someone who familiar with TRIQS to prune this part.
                  auto freq_1 = mesh_bosonic[m + n_w_bosonic - 1];
                  auto freq_2 = mesh_fermionic[n1 + n_w_fermionic];
                  f2w[bl][freq_1, freq_2](a, b, c) -= s * nMw(b1, a, b, n1, n2) * nw(col, m);
                } // n1
              } // m
            } // b
          } // a
        } // c
      } // bl
    }

  }

  // -------------------------------------

  void three_point::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);
    if (measure_g2w) {
      g2w = mpi::all_reduce(g2w, c);
      g2w = g2w / (Z * beta);

      std::vector<std::vector<gf<prod<imfreq, imfreq>, tensor_valued<3>>>> g2w_vec(wdata.gf_struct.size(), std::vector<gf<prod<imfreq, imfreq>, tensor_valued<3>>>(wdata.gf_struct.size()));
      for (int b1 : range(wdata.gf_struct.size())) {
        for (int b2 : range(wdata.gf_struct.size())) {
          g2w_vec[b1][b2] = g2w[b1 * wdata.gf_struct.size() + b2];
        }
      }

      auto g2w_block = make_block2_gf(block_names, block_names, g2w_vec);
      results.g2w = std::move(g2w_block);
    }

    if (measure_f2w) {
      f2w = mpi::all_reduce(f2w, c);
      f2w = f2w / (Z * beta);

      std::vector<std::vector<gf<prod<imfreq, imfreq>, tensor_valued<3>>>> f2w_vec(wdata.gf_struct.size(), std::vector<gf<prod<imfreq, imfreq>, tensor_valued<3>>>(wdata.gf_struct.size()));
      for (int b1 : range(wdata.gf_struct.size())) {
        for (int b2 : range(wdata.gf_struct.size())) {
          f2w_vec[b1][b2] = f2w[b1 * wdata.gf_struct.size() + b2];
        }
      }

      auto f2w_block = make_block2_gf(block_names, block_names, f2w_vec);
      results.f2w = std::move(f2w_block);
    }

  }

  // -------------------------------------

  double three_point::fprefactor(long const &block, std::pair<tau_t, long> const &y) {

    // This function appears in G_F_tau.cpp, four_point.cpp and three_point.cpp
    // If furtherly used, consider to put it in a separate file
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

  std::vector<array<dcomplex, 4>> three_point::compute_Mw(bool is_nMw) {

    // This function appears in four_point.cpp and three_point.cpp
    // If furtherly used, consider to put it in a separate file
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

  // -------------------------------------

  array<dcomplex, 2> three_point::compute_nw() {

    // See Thomas Ayral's doctoral thesis (11.35)
    array<dcomplex, 2> result;
    result.resize(make_shape(wdata.n_color, n_w_bosonic * 2 - 1));
    result() = 0;

    for (int orb: range(wdata.n_color)) {

      for (auto s: config.seglists[orb]) {

        double tau_c = double(s.tau_c);
        double tau_cdag = double(s.tau_cdag);
        
        // Zero frequency: Add up the all the segment length
        if (!is_cyclic(s))
          result(orb, n_w_bosonic - 1) += tau_c - tau_cdag;
        else
          result(orb, n_w_bosonic - 1) += beta - tau_cdag + tau_c;

        // Compute remaining frequencies
        double wm = (-n_w_bosonic + 1) * 2 * M_PI / beta;
        double dw = 2 * M_PI / beta;

        for (auto m = 0; m < result.shape()[1]; ++m) {
          if (m == n_w_bosonic - 1) {
            wm += dw;
            continue;
          }
          dcomplex numerator = std::exp(dcomplex(0, wm * tau_c)) - std::exp(dcomplex(0, wm * tau_cdag));
          result(orb, m) += numerator / dcomplex(0., wm);
          wm += dw;
        }

      } // s

    } // orb

    return result;
  }

} // namespace triqs_ctseg::measures