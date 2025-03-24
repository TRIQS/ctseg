// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE.txt in the root of this distribution for details.

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

    std::map<std::string, nda::array<double, 1>> densities;
    for (long offset = 0; auto [bl_name, bl_size] : wdata.gf_struct) {
      densities[bl_name] = n[range(offset, offset + bl_size)];
      offset += bl_size;
    }
    if (c.rank() == 0) {
      SPDLOG_INFO("Densities:");
      for (auto &[bl, dens] : densities) SPDLOG_INFO("  {}: {}", bl, dens);
    }

    results.densities = std::move(densities);
  }

} // namespace triqs_ctseg::measures
