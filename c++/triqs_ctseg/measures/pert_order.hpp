// Copyright (c) 2022--present, The Simons Foundation
// Copyright (c) 2022--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../results.hpp"
#include "../work_data.hpp"
#include <triqs/stat/lin_binning.hpp>

namespace triqs_ctseg::measures {

  struct pert_order {

    pert_order(std::function<int()> get_order, std::optional<std::vector<double>> &hist_opt,
               std::optional<double> &average_order_opt, std::optional<double> &average_order_error_opt);

    /// Accumulate pert order into histogram
    void accumulate(double s);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &c);

    /// Report current running average
    std::string report() const;

    private:
    // Function to get the pert order
    std::function<int()> get_order;

    // Histogram
    std::vector<double> &hist;

    // Average order
    double &average_order;
    std::optional<double> &average_order_error;

    // Linear binning for error estimation
    triqs::stat::lin_binning<dcomplex> order_bins_;

    // Running sum for report
    double order_sum_ = 0.0;

    // Accumulation counter
    long N = 0;
  };

} // namespace triqs_ctseg::measures
