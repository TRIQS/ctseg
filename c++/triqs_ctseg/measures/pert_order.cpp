// Copyright (c) 2022-2024 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Nikita Kavokine, Olivier Parcollet, Nils Wentzell

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
