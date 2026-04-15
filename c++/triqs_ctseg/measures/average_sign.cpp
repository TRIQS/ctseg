// Copyright (c) 2023--present, The Simons Foundation
// Copyright (c) 2023--present, Max Planck Institute for Polymer Research, Mainz, Germany
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "average_sign.hpp"
#include <itertools/itertools.hpp>
#include <sstream>
#include "../logs.hpp"

namespace triqs_ctseg::measures {

  average_sign::average_sign(params_t const &, work_data_t const &wdata, configuration_t const &config,
                             results_t &results)
     : wdata{wdata}, config{config}, results{results} {
    sign_bins_.emplace(dcomplex{0.0}, 128, 1);
  }

  // -------------------------------------

  void average_sign::accumulate(double s) {
    Z += s;
    N += 1.0;
    *sign_bins_ << dcomplex(s);
  }

  // -------------------------------------

  void average_sign::collect_results(mpi::communicator const &c) {
    Z = mpi::all_reduce(Z, c);
    N = mpi::all_reduce(N, c);

    results.average_sign = Z / N;

    auto [m, err, tau]         = sign_bins_->mean_error_and_tau(c);
    results.average_sign_error = std::abs(err);
  }

  // -------------------------------------

  std::string average_sign::report() const {
    std::ostringstream os;
    os << "Average sign: " << Z / N;
    return os.str();
  }

} // namespace triqs_ctseg::measures
