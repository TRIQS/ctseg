// Copyright (c) 2023--present, The Simons Foundation
// Copyright (c) 2023--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../results.hpp"
#include "../work_data.hpp"
#include <triqs/stat/log_binning.hpp>
#include <triqs/stat/lin_binning.hpp>

namespace triqs_ctseg::measures {

  struct densities {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;

    bool measure_densities_;

    nda::array<double, 1> n;
    double Z = 0;
    long N_  = 0;

    // Log-binning for auto-correlation: [0] = perturbation order, [1..n_color] = sign * density
    std::vector<triqs::stat::log_binning<dcomplex>> log_accs_;

    // Linear binning for density errors (one per block)
    std::vector<triqs::stat::lin_binning<nda::array<dcomplex, 1>>> dens_bins_;

    densities(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };

} // namespace triqs_ctseg::measures
