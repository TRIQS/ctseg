// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace triqs_ctseg::measures {

  struct G_F_tau {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;
    bool measure_F_tau;
    gf_struct_t gf_struct;

    block_gf<imtime> G_tau;
    block_gf<imtime> F_tau;

    double Z;

    G_F_tau(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
    double fprefactor(long const &block, std::pair<tau_t, long> const &y);
  };

} // namespace triqs_ctseg::measures
