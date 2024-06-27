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

#include "densities.hpp"
#include <itertools/itertools.hpp>
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  densities::densities(params_t const &, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    n = nda::zeros<double>(config.n_color());
  }

  // -------------------------------------

  void densities::accumulate(double s) {

    Z += s;
    for (auto const &[c, seglist] : itertools::enumerate(config.seglists)) {
      double sum = 0;
      for (auto &seg : seglist) sum += double(seg.length()); // accounts for cyclicity
      n[c] += s * sum;
    }
  }

  // -------------------------------------

  void densities::collect_results(mpi::communicator const &c) {
    Z = mpi::all_reduce(Z, c);
    n = mpi::all_reduce(n, c);
    n /= (Z * tau_t::beta());
    for (long offset = 0; auto [bl_name, bl_size] : wdata.gf_struct) {
      results.densities[bl_name] = n[range(offset, offset + bl_size)];
      offset += bl_size;
    }
    if (c.rank() == 0) {
      SPDLOG_INFO("Densities:");
      for (auto &[bl, dens] : results.densities) SPDLOG_INFO("  {}: {}", bl, dens);
    }
  }

} // namespace triqs_ctseg::measures
