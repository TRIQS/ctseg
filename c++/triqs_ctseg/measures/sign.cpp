// Copyright (c) 2023-2024 Simons Foundation
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
// Authors: Nikita Kavokine, Nils Wentzell

#include "sign.hpp"
#include <itertools/itertools.hpp>
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  sign::sign(params_t const &, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {
    Z = 0.0;
    N = 0.0;
  }

  // -------------------------------------

  void sign::accumulate(double s) {
    Z += s;
    N += 1.0;
  }

  // -------------------------------------

  void sign::collect_results(mpi::communicator const &c) {
    Z            = mpi::all_reduce(Z, c);
    N            = mpi::all_reduce(N, c);
    results.sign = Z / N;
  }

} // namespace triqs_ctseg::measures
