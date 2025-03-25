// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"
#include "../invariants.hpp"

namespace triqs_ctseg::moves {

  class split_spin_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal dataseg;
    long line_idx, idx_c_up, idx_c_dn, idx_cdag_up, idx_cdag_dn;
    tau_t tau_up, tau_dn;
    double ln_trace_ratio, prop_ratio, det_sign;
    std::tuple<long, long, tau_t> propose(int color);

    public:
    split_spin_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_) {}

    // ------------------
    double attempt();
    double accept();
    void reject();
  };

} // namespace triqs_ctseg::moves
