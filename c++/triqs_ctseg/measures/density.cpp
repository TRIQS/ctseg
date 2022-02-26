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
#include "density.hpp"

namespace density {

  density::density(configuration const &config, nda::array<double, 1> &n_per_color) : config{config}, n_per_color{n_per_color} { n_per_color = 0; }

  // -------------------------------------

  void density::accumulate(double s) {

    Z += s;
    for (auto const &[c, seglist] : itertools::enumerate(config.seglists)) {
      double sum = 0;
      for (auto &seg : seglist) sum += double(seg.tau_c - seg.tau_cdag); // NB : take the cyclic - and cast to double
      n_per_color[c] += s * sum;
    }
  }

  // -------------------------------------

  void density::collect_results(mpi::communicator const &c) {
    Z              = mpi::all_reduce(Z, c);
    n_per_color    = mpi::all_reduce(n_per_color, c);
    nnn_per_colort = n_per_color / Z;
  }
} // namespace density
