// Copyright (c) 2023--present, The Simons Foundation
// Copyright (c) 2023--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../results.hpp"
#include "../work_data.hpp"
#include <triqs/stat/lin_binning.hpp>

namespace triqs_ctseg::measures {

  struct average_sign {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;

    double N = 0;
    double Z = 0;

    std::optional<triqs::stat::lin_binning<dcomplex>> sign_bins_;

    average_sign(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
    std::string report() const;
  };

} // namespace triqs_ctseg::measures
