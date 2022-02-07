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

#include "./precompute_nw.hpp"

namespace triqs_ctseg {

// precomputes the Fourier transform of n(tau)

// call operator which performs the actual computation
void precompute_nw::operator()() {

  nw() = 0;

  for (int orb = 0; orb < config->n_colors(); ++orb) {

    nw(orb, 0) += config->trace.overlap_matrix(
        orb, orb); // length of all segments in orb1

    for (auto it = config->ops_map().begin(orb);
         it != config->ops_map().end(orb);
         ++it) {      // iterate over all operators in orb1
      double sgn = 1; // different signs for creator/annihilator
      if (it->dagger)
        sgn = -1;

      // fast update of exponentials
      double wm = 0;
      double dw = 2 * M_PI / beta;
      std::complex<double> exp = 1.0;
      std::complex<double> dexp =
          std::exp(std::complex<double>(0., dw * it->tau));

      for (size_t m = 1; m < nw.shape()[1]; ++m) {
        wm += dw;
        exp *= dexp;
        nw(orb, m) += sgn * exp / std::complex<double>(0., wm);
      }
    }
  }
}

} // namespace triqs_ctseg
