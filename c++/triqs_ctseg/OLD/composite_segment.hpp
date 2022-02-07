/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet
 *
 * CTSEG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * CTSEG is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * CTSEG. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#pragma once
#include "./segment.hpp"

namespace triqs_ctseg {

// a spin descriptor: creation op + annihilation op
// ordering of c and a is enforced
struct spin_desc {
  colored_const_iterator c, a;
  spin_desc() : c(), a() {}
  spin_desc(colored_const_iterator it1, colored_const_iterator it2)
      : c(it1), a(it2) {
    if (!it1->dagger)
      std::swap(c, a);
  }
  colored_const_iterator sup() { return c->tau > a->tau ? c : a; }
  colored_const_iterator inf() { return c->tau > a->tau ? a : c; }
};

/// A minimal struct describing a composite_segment
struct composite_segment {
  spin_desc plus, minus;
  composite_segment() : plus(), minus() {}
  composite_segment(spin_desc const &plus_, spin_desc const &minus_)
      : plus(plus_), minus(minus_) {}
  composite_segment(segment const &seg_1, segment const &seg_2)
      : plus((seg_1.l->color == 0) ? seg_1.c() : seg_2.c(),
             (seg_1.l->color == 0) ? seg_2.a() : seg_1.a()),
        minus((seg_1.l->color == 0) ? seg_2.c() : seg_1.c(),
              (seg_1.l->color == 0) ? seg_1.a() : seg_2.a()) {}
};
struct composite_segment_desc {
  qmc_time_t tl, tr;
};

} // namespace triqs_ctseg
