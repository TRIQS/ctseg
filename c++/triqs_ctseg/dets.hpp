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

#pragma once
#include <triqs/gfs.hpp>
#include <triqs/det_manip.hpp>

#include "./tau_t.hpp"

using namespace triqs::gfs;
using namespace triqs::mesh;

namespace triqs_ctseg {

  /// A lambda to adapt Delta(tau) for the call by det_manip.
  struct delta_block_adaptor {
    gf<imtime, matrix_real_valued> delta;

    double operator()(std::pair<tau_t, int> const &x, std::pair<tau_t, int> const &y) const {
      double res = delta(double(x.first - y.first))(x.second, y.second);
      return (x.first >= y.first ? res : -res); // x,y first are tau_t, wrapping is automatic in
                                                // the - operation, but need to compute the sign
    }
  };

  using det_t = triqs::det_manip::det_manip<delta_block_adaptor>;

} // namespace triqs_ctseg
