// Copyright (c) 2025--present, The Simons Foundation
// Copyright (c) 2025--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"
#include "../invariants.hpp"

namespace triqs_ctseg::moves {

  class double_insert_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    std::vector<int> colors          = std::vector<int>(2);
    std::vector<tau_t> window_length = std::vector<tau_t>(2);
    std::vector<segment_t> prop_seg  = std::vector<segment_t>(2);
    double det_sign;
    bool is_same_block;

    public:
    double_insert_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_) {};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };

} // namespace triqs_ctseg::moves
