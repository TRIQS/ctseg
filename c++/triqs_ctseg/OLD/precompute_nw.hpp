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
#include "./configuration.hpp"
#include "./qmc_parameters.hpp"

namespace triqs_ctseg {

// precomputes the Fourier transform of n(tau)

struct precompute_nw { // models concept of mc_generic auxiliary precomputation

  const configuration *config;
  const qmc_parameters *params;

  typedef std::complex<double> nw_t;
  array<nw_t, 2> nw;

  double beta;

  // constructor
  precompute_nw(const qmc_parameters *params_, const configuration *config_)
      : params(params_), config(config_), beta(params->beta) {}

  // resize
  void minimum_size(
      long size) { // max to make sure it's large enough if resized twice
    nw.resize(make_shape(config->n_colors(), std::max(size, nw.shape()[1])));
  }

  // accessor
  nw_t get(int color, int n) { return nw(color, n); }

  nw_t get(int color, triqs::gfs::matsubara_freq iom) {
    return nw(color, iom.n);
  }

  CLEF_IMPLEMENT_LAZY_METHOD(precompute_nw, get);

  // call operator which performs the actual computation
  void operator()();
};

} // namespace triqs_ctseg
