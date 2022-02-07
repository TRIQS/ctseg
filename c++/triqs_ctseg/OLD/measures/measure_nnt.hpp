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

namespace triqs_ctseg {

/// measurement for the density-density correlation function in imaginary time
/**
 *Measures of the imaginary-time density-density correlation function

  $$\chi^{\sigma\sigma'}_{ab}(\tau) = \langle T_\tau
 n_{a\sigma}(\tau)n_{b\sigma'}(0) \rangle$$

 * It is triggered through the switch ``measure_nnt``.
 *In contrast to the Green's function, the measurement slows down considerably
 with the number of points ``n_tau_nn`` on the grid.

 *In most cases, it is more efficicient to measure it on Matsubara frequencies
 as specified through ``measure_nnw``:

  $$\chi^{\sigma\sigma'}_{ab}(i\nu_m) = \int_0^\beta d\tau\, e^{i\nu_m\tau}
 \chi^{\sigma\sigma'}_{ab}(\tau)$$

 * where the number of bosonic frequencies :math:`\nu_m=2m\pi/\beta` is
 specified through ``n_w_b_nn``.
  @note Does not exploit time translational invariance
 */
struct measure_nnt {

  const qmc_parameters *params;
  const configuration *config;

  block_gf<imtime> &nnt; //~ Norb*Norb matrix
  const double beta;
  const double Noverbeta;
  double Z;

  measure_nnt(const qmc_parameters *params_, const configuration *config_,
              block_gf<imtime> &nnt_);
  void accumulate(double s);

  /// average results and normalize nnt
  void collect_results(mpi::communicator const &c);
};
} // namespace triqs_ctseg
