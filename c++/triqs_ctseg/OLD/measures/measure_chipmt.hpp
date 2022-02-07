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
#include "../qmc_parameters.hpp"

namespace triqs_ctseg {
/// Measure for the transverse correlation function in imaginary time
/**
 *Measure $$\chi_\pm(\tau) = 0.5 \langle T s_+(\tau)s_-(0)\rangle$$
 *Binning meshes:
 * let N be the number of "time slices" of G, i.e. the physical number of
 *storage bins, N = Gt[k].numberTimeSlices with indices n=0,...,N-1; then for
 * integer mesh
 * the mapping from tau to the bin index n is
 * n = int(floor(tau (N-1)/beta + 0.5))
 * the first bin is [0,beta/(2(N-1)) )
 * the last bin is [beta-beta/(2(N-1)), beta ] (tau=beta goes to the last bin)
 * the n-th bin is [n beta/(N-1) - beta /(2(N-1)), n beta/(N-1) + beta/(2(N-1))
 *) width of inner bins is beta/(N-1), width of outer bins (indices 0, N-1) is
 *beta/(2(N-1)) center of bins is at _integer_ multiples of beta/(N-1)
 */

struct measure_chipmt {
  const qmc_parameters *params;
  const configuration *config;
  gf<imtime> &chipmt;
  double beta;
  double Noverbeta;
  double Z;
  /// constructor
  measure_chipmt(const qmc_parameters *params_, const configuration *config_,
                 gf<imtime> &chipmt_);
  /// accumulate the Green's function
  void accumulate(double s);
  /// reduce and normalize G
  void collect_results(mpi::communicator const &c);
};

} // namespace triqs_ctseg
