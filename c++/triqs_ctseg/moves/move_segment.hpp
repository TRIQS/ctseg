// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE.txt in the root of this distribution for details.

#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"
#include "../invariants.hpp"

namespace triqs_ctseg::moves {

  class move_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    bool flipped; // whether we flip an antisegment
    int origin_color, dest_color;
    segment_t origin_segment;
    long origin_index, dest_index;
    double det_sign;
    std::vector<segment_t> sl, dsl;

    public:
    move_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_) {};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };

} // namespace triqs_ctseg::moves
