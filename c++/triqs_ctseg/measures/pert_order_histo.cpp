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

#include "./pert_order_histo.hpp"
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  pert_order_histo::pert_order_histo(params_t const &p, work_data_t const &wdata,
                                                     configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    histo_delta = triqs::stat::histogram{0, p.histogram_max_order};
    histo_Jperp = triqs::stat::histogram{0, p.histogram_max_order};
  }

  // -------------------------------------

  void pert_order_histo::accumulate(double) {

    // For the config, compute the number of segments
    long n_segments =
       std::accumulate(begin(config.seglists), end(config.seglists), 0, [](long r, auto &&v) { return r + v.size(); });
    long n_Splusminus = 2 * config.Jperp_list.size();
    long delta_order =
       n_segments - n_Splusminus; // half # of c, cdag operators  = (2 * n_segments - 2* n_Splusminus) /2

    histo_delta << delta_order;
    histo_Jperp << config.Jperp_list.size();
  }

  // -------------------------------------

  void pert_order_histo::collect_results(mpi::communicator const &c) {

    histo_delta = mpi::all_reduce(histo_delta, c);
    histo_Jperp = mpi::all_reduce(histo_Jperp, c);

    // Normalize
    results.pert_order_histo_Delta = pdf(histo_delta);
    results.pert_order_histo_Jperp = pdf(histo_Jperp);
  }

} // namespace triqs_ctseg::measures
