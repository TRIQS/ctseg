// Copyright (c) 2023--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE.txt in the root of this distribution for details.

#include "average_sign.hpp"
#include <itertools/itertools.hpp>
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  average_sign::average_sign(params_t const &, work_data_t const &wdata, configuration_t const &config,
                             results_t &results)
     : wdata{wdata}, config{config}, results{results} {
    Z = 0.0;
    N = 0.0;
  }

  // -------------------------------------

  void average_sign::accumulate(double s) {
    Z += s;
    N += 1.0;
  }

  // -------------------------------------

  void average_sign::collect_results(mpi::communicator const &c) {
    Z = mpi::all_reduce(Z, c);
    N = mpi::all_reduce(N, c);

    results.average_sign = Z / N;
  }

} // namespace triqs_ctseg::measures
