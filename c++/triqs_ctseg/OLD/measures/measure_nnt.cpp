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
#include "./measure_nnt.hpp"

namespace triqs_ctseg {

measure_nnt::measure_nnt(const qmc_parameters *params_,
                         const configuration *config_, block_gf<imtime> &nnt_)
    : params(params_), config(config_), nnt(nnt_), beta(params->beta),
      Noverbeta(1.0 / nnt[0].mesh().delta()), Z(0.0) {
  for (auto &g : nnt)
    g() = 0;
}

void measure_nnt::accumulate(double s) {

  Z += s;

  for (int orb2 = 0; orb2 < config->n_colors(); orb2++) {
    auto bl_ind_j = config->color_to_block_and_inner_index(orb2); // slow
    if (config->edge_state(orb2) == 1) { // if n_2(tau=0)=1
      for (int orb1 = 0; orb1 < config->n_colors(); orb1++) {
        auto bl_ind_i = config->color_to_block_and_inner_index(orb1); // slow
        int bl = bl_ind_i.first * config->n_blocks() + bl_ind_j.first;

        int tau_idx = 0;
        double local_density; //(0 or 1)
        double tau_lim;
        bool last = false;
        auto it = config->ops_map().end(orb1);

        do {
          if (it != config->ops_map().begin(orb1)) {
            it--;
            tau_lim = double(it->tau);
            local_density = it->right_occupation()[orb1];
          } else {
            last = true;
            tau_lim = double(params->beta);
            local_density = config->edge_state(orb1);
          }

          while (tau_idx <= std::floor(tau_lim * Noverbeta + 0.5)) {
            nnt[bl][tau_idx](bl_ind_i.second, bl_ind_j.second) +=
                s * local_density;
            tau_idx++;
          }

        } while (!last);

      } // orb1
    }   // if occupied
  }     // orb2
}

/// average results and normalize nnt
void measure_nnt::collect_results(mpi::communicator const &c) {
  Z = mpi::all_reduce(Z, c);
  nnt = mpi::all_reduce(nnt, c);
  nnt = nnt / Z;
}
} // namespace triqs_ctseg
