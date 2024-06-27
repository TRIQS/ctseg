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
