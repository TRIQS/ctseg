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
#include <itertools/itertools.hpp>

namespace measures {

  density::density(configuration_t const &config, results_t &results) : config{config}, densities{results.densities} {
    densities = nda::array<double, 1>(config.n_color());
  }

  // -------------------------------------

  void density::accumulate(double s) {

    Z += s;
    for (auto const &[c, seglist] : itertools::enumerate(config.seglists)) {
      double sum = 0;
      for (auto &seg : seglist) sum += double(seg.tau_c - seg.tau_cdag); // NB : take the cyclic - and cast to double
      densities[c] += s * sum;
    }
  }

  // -------------------------------------

  void density::collect_results(mpi::communicator const &c) {
    Z         = mpi::all_reduce(Z, c);
    densities = mpi::all_reduce(densities, c);
    densities /= Z;
  }
} // namespace measures
