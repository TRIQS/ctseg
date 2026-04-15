// Copyright (c) 2024--present, The Simons Foundation
// Copyright (c) 2024--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./pert_order.hpp"
#include <sstream>
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  pert_order::pert_order(std::function<int()> get_order, std::optional<std::vector<double>> &hist_opt,
                         std::optional<double> &average_order_opt, std::optional<double> &average_order_error_opt)
     : get_order{get_order},
       hist{hist_opt.emplace(4, 0.0)},
       average_order{average_order_opt.emplace(0.0)},
       average_order_error{average_order_error_opt},
       order_bins_(dcomplex{0.0}, 128, 1) {}

  // -------------------------------------

  void pert_order::accumulate(double) {
    auto order = get_order();
    while (order >= hist.size()) hist.resize(2 * hist.size());
    hist[order] += 1;
    order_sum_ += order;
    order_bins_ << dcomplex(double(order));
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

    auto [m, err, tau]  = order_bins_.mean_error_and_tau(c);
    average_order_error = std::abs(err);
  }

  // -------------------------------------

  std::string pert_order::report() const {
    std::ostringstream os;
    os << "Average perturbation order: " << (N > 0 ? order_sum_ / N : 0.0);
    return os.str();
  }

} // namespace triqs_ctseg::measures
