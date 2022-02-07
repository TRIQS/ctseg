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
#include "./measure_nnw.hpp"

namespace triqs_ctseg {

measure_nnw::measure_nnw(const qmc_parameters *params_,
                         const configuration *config_, block_gf<imfreq> &nnw_,
                         std::shared_ptr<precompute_nw> aux_measure_n_)
    : params(params_), config(config_), nnw(nnw_), nw(aux_measure_n_),
      beta(params->beta), n_imfreq(nnw[0].mesh().last_index() + 1), Z(0.0) {
  for (auto &g : nnw)
    g() = 0;
  nw->minimum_size(n_imfreq);
}

/// accumulate nnw
void measure_nnw::accumulate(double s) {

  Z += s;

  for (int bl = 0; bl < nnw.size(); bl++) {
    int b1 = bl / config->gf_struct().size();
    int b2 = bl % config->gf_struct().size();
    for (int a = 0; a < nnw[bl].target_shape()[0]; a++) {
      int i = config->block_and_inner_index_to_color(b1, a);
      for (int b = 0; b < nnw[bl].target_shape()[1]; b++) {
        int j = config->block_and_inner_index_to_color(b2, b);
        for (int m = 0; m < n_imfreq; ++m) {
          nnw[bl][m](a, b) +=
              s * nw->get(i, m) * std::conj(nw->get(j, m)) / beta;
        } // m
      }   // b
    }     // a
  }       // bl

} // accumulate

/// average results and normalize nnw
void measure_nnw::collect_results(mpi::communicator const &c) {
  Z = mpi::all_reduce(Z, c);
  nnw = mpi::all_reduce(nnw, c);
  nnw = nnw / Z;
  for (int bl = 0; bl < nnw.size(); bl++)
    for (auto const &w : nnw[bl].mesh())
      if (w.n < 0)
        nnw[bl][w] = nnw[bl](matsubara_freq(-w.n, beta, Boson));
}
} // namespace triqs_ctseg
