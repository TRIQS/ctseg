// Copyright (c) 2024--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE.txt in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace triqs_ctseg::measures {

  struct four_point {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;
    double w_ini, w_inc;
    bool measure_g3w;
    bool measure_f3w;
    int n_w_fermionic;
    int n_w_bosonic;
    std::vector<std::string> block_names;
    std::vector<array<dcomplex, 4>> Mw_vector, nMw_vector;
    nda::vector<dcomplex> y_exp_ini, y_exp_inc, x_exp_ini, x_exp_inc;
    nda::vector<int> y_inner_index, x_inner_index;

    block_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> g3w;
    block_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> f3w;

    double Z = 0;

    four_point(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
    double fprefactor(long const &block, std::pair<tau_t, long> const &y);
    
    std::vector<array<dcomplex, 4>> compute_Mw(bool is_nMw);

    dcomplex Mw(long const &block, int const &i, int const &j, int const &n1, int const &n2) {
      return Mw_vector[block](i, j, n1 + n_w_fermionic + n_w_bosonic - 1, n2 + n_w_fermionic + n_w_bosonic - 1);
    }

    dcomplex nMw(long const &block, int const &i, int const &j, int const &n1, int const &n2) {
      return nMw_vector[block](i, j, n1 + n_w_fermionic + n_w_bosonic - 1, n2 + n_w_fermionic + n_w_bosonic - 1);
    }

  };

} // namespace triqs_ctseg::measures