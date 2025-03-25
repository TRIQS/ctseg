// Copyright (c) 2024--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"
#include "../util.hpp"

namespace triqs_ctseg::measures {

  struct state_hist {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;

    nda::vector<double> H;

    double Z = 0;

    state_hist(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };

} // namespace triqs_ctseg::measures
