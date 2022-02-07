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
#include "../precompute_fprefactor.hpp"
#include "../qmc_parameters.hpp"

namespace triqs_ctseg {

using namespace triqs::gfs;
using namespace triqs::mesh;

/// Measure for the Green's function in imaginary time
/**
 *The imaginary-time Green's function is defined as
   $$G^\sigma_{ab}(\tau) = -\langle T_\tau
 c_{a\sigma}(\tau)c_{b\sigma}^\dagger(0) \rangle$$
 * and the corresponding improved estimator is given by
   $$F_{ab}^{\sigma}(\tau) = -\int_0^\beta d\tilde{\tau} \sum_{c\sigma'} \langle
 T_\tau n_{c\sigma'}(\tilde{\tau})
 \mathcal{U}^{\sigma\sigma'}_{ac}(\tilde{\tau}-\tau)
 c_{a\sigma}(\tau)c_{b\sigma'}^\dagger(0) \rangle$$

 *The imaginary-time measurement is most efficient. The performance of the
 algorithm does not scale with the number of points in the grid on which it is
 measured, so this number can and should be chosen large. By Nyquist's theorem,
 the Fourier transform will be correctly reproduce the function in the frequency
 domain on the first :math:`N_\omega\approx N_\tau/4\pi` frequencies.
 *
 *These measurements are turned on by setting ``measure_gt`` and ``measure_ft``
 to ``true``, respectively. *The number of time points on the grid is specified
 through ``n_tau`` and is the same for both observables.
 */
struct measure_gt {

  const qmc_parameters *params;
  const configuration *config;

  block_gf<imtime> &gt;
  block_gf<imtime> &ft;

  std::shared_ptr<precompute_fprefactor> fprefactor;

  double beta, Noverbeta, Z;
  accumulator<double> gt_stack = {0.0, -1, -1};
  accumulator<double> Z_stack  = {0.0, -1, -1};
  int counter;
  double accum;

  /// constructor
  measure_gt(const qmc_parameters *params_, const configuration *config_,
             std::shared_ptr<precompute_fprefactor> fprefactor_,
             block_gf<imtime> &gt_, block_gf<imtime> &ft_);

  /// accumulate the Green's function
  void accumulate(double s);

  /// reduce and normalize G
  void collect_results(mpi::communicator const &c);
};
} // namespace triqs_ctseg
