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
// Authors: Nikita Kavokine, Hao Lu, Olivier Parcollet, Nils Wentzell

#pragma once

#include <triqs/gfs.hpp>
using namespace triqs::gfs;

namespace triqs_ctseg {

  inline long modulo(long i, long N) { return (i + N) % N; }

  // Same as std::lower_bound, but i-th element of vector is returned by f[i]
  // f is called on 0:N strictly
  long lower_bound(auto f, long N, auto const &value) {
    long first = 0, count = N;

    while (count > 0) {
      long step = count / 2;
      assert(first + step < N);
      if (f(first + step) < value) {
        first += step + 1;
        count -= step + 1;
      } else
        count = step;
    }
    return first;
  }

  // If we have an ordered det_manip d, we get the lower_bound index for x
  long det_lower_bound_x(auto const &d, auto const &x) {
    return lower_bound([&d](long i) { return d.get_x(i).first; }, d.size(), x);
  }
  long det_lower_bound_y(auto const &d, auto const &y) {
    return lower_bound([&d](long i) { return d.get_y(i).first; }, d.size(), y);
  }

  // Integer power
  constexpr unsigned int ipow(unsigned int n, unsigned int m) { return m == 0 ? 1 : m == 1 ? n : n * ipow(n, m - 1); }

  // Block2Gf constructor
  template <typename M> block2_gf<M> make_block2_gf(M const &m, gf_struct_t const &gf_struct) {

    std::vector<std::vector<gf<M>>> gf_vecvec;
    std::vector<std::string> block_names;

    for (auto const &[bl1, bl1_size] : gf_struct) {
      block_names.push_back(bl1);
      std::vector<gf<M>> gf_vec;
      for (auto const &[bl2, bl2_size] : gf_struct) { gf_vec.emplace_back(m, make_shape(bl1_size, bl2_size)); }
      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
  }

} // namespace triqs_ctseg
