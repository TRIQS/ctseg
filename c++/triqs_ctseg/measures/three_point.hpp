// Copyright (c) 2025--present, The Simons Foundation
// Copyright (c) 2025--present, Max Planck Institute for Polymer Research, Mainz, Germany
// Copyright (c) 2025--present, EPFL, Switzerland
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace triqs_ctseg::measures {

  struct three_point {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;
    double w_ini, w_inc;
    bool measure_g2w;
    bool measure_f2w;
    int n_w_fermionic;
    int n_w_bosonic;
    std::vector<std::string> block_names;
    std::vector<array<dcomplex, 4>> Mw_vector, nMw_vector;
    array<dcomplex, 2> nw_vector;
    nda::vector<dcomplex> y_exp_ini, y_exp_inc, x_exp_ini, x_exp_inc;
    nda::vector<int> y_inner_index, x_inner_index;

    block_gf<prod<imfreq, imfreq>, tensor_valued<3>> g2w;
    block_gf<prod<imfreq, imfreq>, tensor_valued<3>> f2w;

    double Z = 0;

    three_point(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
    double fprefactor(long const &block, std::pair<tau_t, long> const &y);
    
    std::vector<array<dcomplex, 4>> compute_Mw(bool is_nMw);
    array<dcomplex, 2> compute_nw();

    dcomplex Mw(long const &block, int const &i, int const &j, int const &n1, int const &n2) {
      return Mw_vector[block](i, j, n1 + n_w_fermionic + n_w_bosonic - 1, n2 + n_w_fermionic + n_w_bosonic - 1);
    }

    dcomplex nMw(long const &block, int const &i, int const &j, int const &n1, int const &n2) {
      return nMw_vector[block](i, j, n1 + n_w_fermionic + n_w_bosonic - 1, n2 + n_w_fermionic + n_w_bosonic - 1);
    }

    dcomplex nw(int const &c, int const &m) {
      return nw_vector(c, m + n_w_bosonic - 1);
    }

  };

} // namespace triqs_ctseg::measures