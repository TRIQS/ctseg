// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./nn_static.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  nn_static::nn_static(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta    = p.beta;
    n_color = config.n_color();
    nn      = nda::zeros<double>(n_color, n_color);

    for (auto const &[bl1, bl1_size] : wdata.gf_struct)
      for (auto const &[bl2, bl2_size] : wdata.gf_struct)
        nn_bins_.emplace_back(nda::array<dcomplex, 2>(nda::zeros<dcomplex>(bl1_size, bl2_size)), 128, 1);
  }

  // -------------------------------------

  void nn_static::accumulate(double s) {

    LOG("\n =================== MEASURE <nn>  ================ \n");

    Z += s;
    ++N_;

    // Accumulate nn and build per-step matrix for error analysis
    nda::array<dcomplex, 2> nn_step = nda::zeros<dcomplex>(n_color, n_color);
    for (int a = 0; a < n_color; ++a)
      for (int b = 0; b < n_color; ++b) {
        double ov = 0;
        for (auto const &sa : config.seglists[a])
          for (auto const &sb : config.seglists[b]) ov += overlap(sa, sb);
        nn(a, b) += s * ov;
        nn_step(a, b) = dcomplex(s * ov / beta);
      }

    // Feed block-pair slices into lin_binning
    int bin_idx = 0;
    for (long x1 = 0; auto const &[bl1, bl1_size] : wdata.gf_struct) {
      for (long x2 = 0; auto const &[bl2, bl2_size] : wdata.gf_struct) {
        nn_bins_[bin_idx] << nda::array<dcomplex, 2>(nn_step(range(x1, x1 + bl1_size), range(x2, x2 + bl2_size)));
        ++bin_idx;
        x2 += bl2_size;
      }
      x1 += bl1_size;
    }
  }
  // -------------------------------------

  void nn_static::collect_results(mpi::communicator const &c) {

    Z  = mpi::all_reduce(Z, c);
    nn = mpi::all_reduce(nn, c);
    nn = nn / Z / beta;

    std::map<std::pair<std::string, std::string>, nda::matrix<double>> nn_block;
    for (long x1 = 0; auto &[bl1, bl1_size] : wdata.gf_struct) {
      for (long x2 = 0; auto &[bl2, bl2_size] : wdata.gf_struct) {
        nn_block[{bl1, bl2}] = nn(range(x1, x1 + bl1_size), range(x2, x2 + bl2_size));
        x2 += bl2_size;
      }
      x1 += bl1_size;
    }

    // store the result (not reused later, hence we can move it).
    results.nn_static = std::move(nn_block);

    // Compute error bars from linear binning
    N_        = mpi::all_reduce(N_, c);
    auto norm = std::abs(dcomplex(Z) / dcomplex(N_));

    std::map<std::pair<std::string, std::string>, nda::matrix<double>> nn_errors;
    int bin_idx = 0;
    for (auto &[bl1, bl1_size] : wdata.gf_struct) {
      for (auto &[bl2, bl2_size] : wdata.gf_struct) {
        auto [m, err, tau]    = nn_bins_[bin_idx].mean_error_and_tau(c);
        nn_errors[{bl1, bl2}] = nda::matrix<double>(nda::abs(err) / norm);
        ++bin_idx;
      }
    }
    results.nn_static_errors = std::move(nn_errors);
  }

} // namespace triqs_ctseg::measures
