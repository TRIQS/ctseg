// Copyright (c) 2025--present, The Simons Foundation
// Copyright (c) 2025--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include "../configuration.hpp"

namespace triqs_ctseg::measures {

  struct visualize_config {

    mpi::communicator _c;

    configuration_t const &config;

    double Z = 0;

    visualize_config(configuration_t const &config): config{config} {};

    void accumulate(double s) {

      Z += s;

      if (_c.rank() == 0) {

        std::ofstream outFile("configuration.txt", std::ios::app);

        outFile << "Configuration is " << config << std::endl;

        outFile.close();

      }

    }

    void collect_results(mpi::communicator const &c) {

      Z = mpi::all_reduce(Z, c);

    }

  };

} // namespace triqs_ctseg::measures