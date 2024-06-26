// Copyright (c) 2024 Simons Foundation
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
// Authors: Hao Lu, Nils Wentzell

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"
#include "../util.hpp"

namespace triqs_ctseg::measures {

  struct state_hist {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;

    nda::vector<double> H;

    double Z = 0;

    state_hist(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };

} // namespace triqs_ctseg::measures
