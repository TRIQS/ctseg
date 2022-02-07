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
#include "../configuration.hpp"
#include "../precompute_nw.hpp"

namespace triqs_ctseg {

/// measurement for the density-density correlation function on Matsubara
/// frequencies
/**
 *Measures the Fourier transform of
 *$$\chi_{a,b}^{\sigma\sigma'}(\tau)\equiv\langle
 *n_{a\sigma}(\tau)n_{b\sigma'}(0)\rangle$$
 */
struct measure_nnw {

  const qmc_parameters *params;
  const configuration *config;

  block_gf<imfreq> &nnw; //~ Norb*Norb matrix
  const double beta;
  const int n_imfreq;
  double Z;

  std::shared_ptr<precompute_nw> nw;

  measure_nnw(const qmc_parameters *params_, const configuration *config_,
              block_gf<imfreq> &nnw_,
              std::shared_ptr<precompute_nw> aux_measure_n_);

  /// accumulate nnw
  void accumulate(double s);

  /// average results and normalize nnw
  void collect_results(mpi::communicator const &c);
};

} // namespace triqs_ctseg
