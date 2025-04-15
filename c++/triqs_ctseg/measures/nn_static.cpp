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
  }

  // -------------------------------------

  void nn_static::accumulate(double s) {

    LOG("\n =================== MEASURE <nn>  ================ \n");

    Z += s;

    for (int a = 0; a < n_color; ++a)
      for (int b = 0; b < n_color; ++b) {
        for (auto const &sa : config.seglists[a]) {
          for (auto const &sb : config.seglists[b]) //
            nn(a, b) += s * overlap(sa, sb);
        }
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
  }

} // namespace triqs_ctseg::measures
