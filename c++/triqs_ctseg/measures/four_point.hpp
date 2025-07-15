// Copyright (c) 2024--present, The Simons Foundation
// Copyright (c) 2024--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace triqs_ctseg::measures {

  struct four_point {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;
    int n_w_bosonic, n_w_fermionic;
    triqs::mesh::imfreq mesh_bosonic, mesh_fermionic;
    std::vector<std::string> block_names;

    std::vector<std::vector<gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>>>> g3w;
    block2_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> g3w_block;

    double Z = 0;

    four_point(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
    
    block_gf<prod<imfreq, imfreq>> compute_Mw();

  };

} // namespace triqs_ctseg::measures