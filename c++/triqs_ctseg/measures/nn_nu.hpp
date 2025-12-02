// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace triqs_ctseg::measures {

  struct nn_nu {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;
    double dtau;
    int ntau;
    std::vector<long> block_number, index_in_block;

    gf<dlr_imfreq> q_nu;
    block2_gf<dlr_imfreq> q_nu_block;

    double Z = 0;
    int n_color;

    nn_nu(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };

} // namespace triqs_ctseg::measures
