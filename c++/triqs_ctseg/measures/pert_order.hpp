// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../results.hpp"
#include "../work_data.hpp"

namespace triqs_ctseg::measures {

  struct pert_order {

    pert_order(std::function<int()> get_order, std::optional<std::vector<double>> &hist_opt,
               std::optional<double> &average_order_opt);

    /// Accumulate pert order into histogram
    void accumulate(double s);

    /// Reduce and normalize
    void collect_results(mpi::communicator const &c);

    private:
    // Function to get the pert order
    std::function<int()> get_order;

    // Histogram
    std::vector<double> &hist;

    // Average order
    double &average_order;

    // Accumulation counter
    long N = 0;
  };

} // namespace triqs_ctseg::measures
