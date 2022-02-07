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
#include "./measure_chipmt.hpp"
namespace triqs_ctseg {

measure_chipmt::measure_chipmt(const qmc_parameters *params_,
                               const configuration *config_,
                               gf<imtime> &chipmt_)
    : params(params_), config(config_), chipmt(chipmt_), beta(params->beta),
      Noverbeta(1.0 / chipmt.mesh().delta()) {
  Z = 0;
  // chipmt()=0; //tail problem when using mpi
  chipmt.data()() = 0;
}

void measure_chipmt::accumulate(double s) {
  Z += s;
  // loop times of operators
  for (int i = 0; i < config->boson_lines.size(); i++) {
    auto tau = double(config->boson_lines[i].plus.c->tau -
                      config->boson_lines[i].minus.c->tau);
    chipmt[closest_mesh_pt(tau)](0, 0) += -s / config->boson_lines.jperp(tau);
  }
}

void measure_chipmt::collect_results(mpi::communicator const &c) {
  Z = mpi::all_reduce(Z, c);
  chipmt = mpi::all_reduce(chipmt, c);
  chipmt = chipmt / (beta * Z / Noverbeta);
  chipmt[0] *= 2;
  chipmt[chipmt.mesh().size() - 1] *= 2;
}
} // namespace triqs_ctseg
