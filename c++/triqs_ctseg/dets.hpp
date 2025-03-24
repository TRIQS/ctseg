// Copyright (c) 2022--present, The TRIQS/ctseg Authors and their Assignees
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE.txt in the root of this distribution for details.

#pragma once
#include <triqs/gfs.hpp>
#include <triqs/det_manip.hpp>

#include "./tau_t.hpp"

using namespace triqs::gfs;
using namespace triqs::mesh;

namespace triqs_ctseg {

  /// A lambda to adapt Delta(tau) for the call by det_manip.
  struct Delta_block_adaptor {
    gf<imtime, matrix_real_valued> Delta;

    double operator()(std::pair<tau_t, int> const &x, std::pair<tau_t, int> const &y) const {
      double res = Delta(double(x.first - y.first))(x.second, y.second);
      return (x.first >= y.first ? res : -res); // x,y first are tau_t, wrapping is automatic in
                                                // the - operation, but need to compute the sign
    }
  };

  using det_t = triqs::det_manip::det_manip<Delta_block_adaptor>;

} // namespace triqs_ctseg
