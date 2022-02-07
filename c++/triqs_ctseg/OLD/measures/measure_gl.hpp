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
#include "../types.hpp"
#include <triqs/utility/legendre.hpp>

namespace triqs_ctseg {

/// Measure for the Green's function in Legendre basis
/**
 *The Legendre coefficients of the Green's function and the improved estimator
 are defined as

 $$X(l)=\sqrt{2l+1}\int_0^\beta d\tau\,P_l[x(\tau)]X(\tau)$$

 * with $X=G^\sigma_{ab},F^\sigma_{ab}$, $x(\tau)=2\tau/\beta-1$ and $P_l(x)$
 are the Legendre polynomials, defined in the [-1,1] interval.
 *
 *These measurements are controlled through the switches and parameter
 ``measure_gl``, ``measure_fl`` and ``n_legendre_g``.

 *The Legendre Green's function may be transformed to the Matsubara basis
 through the unitary transformation

   $G_a(i\omega_n) = \sum_{l\geq 0}T_{nl} G_a(l)$

 * where

   $T_{nl} = (-1)^ni^{l+1}\sqrt{2l+1}j_l\left(\frac{(2n+1)\pi}{2}\right)$

 * with the spherical Bessel functions $j_l(z)$.
 */
struct measure_gl {

  const qmc_parameters *params;
  const configuration *config;

  g_l_t &gl;
  g_l_t &fl;

  std::shared_ptr<precompute_fprefactor> fprefactor;

  double beta;
  double Z;
  int n_l;

  triqs::utility::legendre_generator Tn;

  // constructor
  measure_gl(const qmc_parameters *params_, const configuration *config_,
             std::shared_ptr<precompute_fprefactor> fprefactor_, g_l_t &gl_,
             g_l_t &fl_);

  // accumulate the Green's function
  void accumulate(double s);

  // reduce and normalize G
  void collect_results(mpi::communicator const &c);
};

} // namespace triqs_ctseg
