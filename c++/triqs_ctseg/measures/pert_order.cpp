// Copyright (c) 2024--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./pert_order.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  pert_order::pert_order(std::function<int()> get_order, std::optional<std::vector<double>> &hist_opt,
                         std::optional<double> &average_order_opt)
     : get_order{get_order}, hist{hist_opt.emplace(4, 0.0)}, average_order{average_order_opt.emplace(0.0)} {}

  // -------------------------------------

  void pert_order::accumulate(double) {
    auto order = get_order();
    while (order >= hist.size()) hist.resize(2 * hist.size());
    hist[order] += 1;
    ++N;
  }

  // -------------------------------------

  void pert_order::collect_results(mpi::communicator const &c) {
    N = mpi::all_reduce(N, c);

    // Make sure that all mpi threads have an equally sized hist
    auto max_size = mpi::all_reduce(hist.size(), c, MPI_MAX);
    hist.resize(max_size, 0.0);

    // Reduce hist over mpi threads
    hist = mpi::all_reduce(hist, c);

    // Normalize and Calculate average order
    for (int order : range(hist.size())) {
      hist[order] /= N;
      average_order += hist[order] * order;
    }
  }

} // namespace triqs_ctseg::measures
