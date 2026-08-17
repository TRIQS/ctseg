// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace triqs_ctseg::measures {

  struct Sperp_asym_tau {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;

    gf<imtime> sm_sp_tau;
    gf<imtime> sp_sm_tau;

    double Z = 0;

    Sperp_asym_tau(params_t const &params, work_data_t const &wdata, configuration_t const &config,
                   results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
  };

} // namespace triqs_ctseg::measures
