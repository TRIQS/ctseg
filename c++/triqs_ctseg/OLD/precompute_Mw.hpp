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
#include "./params.hpp"
#include "./precompute_fprefactor.hpp"
#include "./qmc_parameters.hpp"
#include <nda/nda.hpp>

namespace triqs_ctseg {

/// precomputes the Fourier transform of the M-matrix and nM-matrix
/**
  *Measures the M and nM matrices:

   $$M_\sigma(w_1,w_2) = \sum_{ij} M^{\sigma}_{ji} e^{i w_1 \tau_i} e^{-i w_2
  \tau_j}$$
   $$nM_\sigma(w_1,w_2) = \sum_{ij} I^\sigma_j M^{\sigma}_{ji} e^{i w_1 \tau_i}
  e^{-i w_2 \tau_j}$$

  * The measurement of $I^\sigma_j$ is done in [[precompute_fprefactor]].

  @note It is needed for the measurement of vertex functions. It models the
  concept of mc_generic auxiliary precomputation
  @warning If NFFT is used, N_w must be chosen such that the number of
  frequencies in the NFFT is >=8! Otherwise segfaults may occur.  (requires nfft
  version >= 3.2)
 */
struct precompute_Mw {

  const solve_params_t *p;
  const qmc_parameters *params;
  const configuration *config;

  std::vector<array<dcomplex, 4>> Mw;
  /// M-matrix left-multiplied by density operator; needed for improved
  /// estimator
  std::vector<array<dcomplex, 4>> nMw;

  std::shared_ptr<precompute_fprefactor> fprefactor;

  int n_w_fermionic; // number of fermionic frequencies for vertex functions
  int n_w_bosonic;   // number of bosonic frequencies for vertex functions
  int n_w_aux;       // number of fermionic frequencies for which to measure M
  int n_w_aux_1, n_w_aux_2;
  double beta;
  double w_ini, w_inc;

#ifdef HAVE_NFFT
  vector<dcomplex> c_exp;
  vector<dcomplex> cdag_exp;
#else
  // moved from #else because they are needed in any case
  vector<dcomplex> c_exp_ini;
  vector<dcomplex> c_exp_inc;
  vector<dcomplex> cdag_exp_ini;
  vector<dcomplex> cdag_exp_inc;
  vector<int> c_inner_index, cdag_inner_index;
#endif

  /// constructor
  /*
   * @param params_ QMC parameters
   * @param p_ solve parameters
   * @param fprefactor fprefactor (see [precompute_fprefactor])
   * @param n_w_fermionic max fermionic Matsubara index
   * @param n_w_bosonic max bosonic Matsubara index
   */
  precompute_Mw(const qmc_parameters *params_, const solve_params_t *p_,
                const configuration *config_,
                std::shared_ptr<precompute_fprefactor> fprefactor_,
                int n_w_fermionic_, int n_w_bosonic_);

  /// accessors
  /*
   * @param block block index
   * @param i orbital index
   * @param j orbital index
   * @param n1 matsubara index (fermionic)
   * @param n2 matsubara index (fermionic)
   * @return $M^\sigma_{ij}(i\omega_n_1,i\omega_n_2)$
   * @warning -n_w_fermionic <= n1 < n_w_fermionic
   * @warning -n_w_fermionic <= n2 < n_w_fermionic+n_w_bosonic
   */
  dcomplex getM(int block, int i, int j, int n1, int n2) {
    return Mw[block](i, j, n1 + n_w_fermionic, n2 + n_w_fermionic);
  }
  dcomplex getnM(int block, int i, int j, int n1, int n2) {
    return nMw[block](i, j, n1 + n_w_fermionic, n2 + n_w_fermionic);
  }

  dcomplex getM(int block, int i, int j, triqs::gfs::matsubara_freq nu1,
                triqs::gfs::matsubara_freq nu2) {
    return Mw[block](i, j, nu1.n + n_w_fermionic, nu2.n + n_w_fermionic);
  }

  dcomplex getnM(int block, int i, int j, triqs::gfs::matsubara_freq nu1,
                 triqs::gfs::matsubara_freq nu2) {
    return nMw[block](i, j, nu1.n + n_w_fermionic, nu2.n + n_w_fermionic);
  }

  CLEF_IMPLEMENT_LAZY_METHOD(precompute_Mw, getM);
  CLEF_IMPLEMENT_LAZY_METHOD(precompute_Mw, getnM);

  /// call operator which performs the actual computation
  void operator()();
};

} // namespace triqs_ctseg
