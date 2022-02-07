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
#include "./measure_gw.hpp"

// implements frequency measurement of Green's function
// if NFFT is used, N_w must be chosen such that the number of
// frequencies in the NFFT is >=8! Otherwise segfaults may occur.
// requires nfft version >= 3.2
#ifdef HAVE_NFFT
#include "nfft3.h"
#endif

namespace triqs_ctseg {

// constructor
measure_gw::measure_gw(const qmc_parameters *params_, const solve_params_t *p_,
                       const configuration *config_,
                       std::shared_ptr<precompute_fprefactor> fprefactor_,
                       block_gf<imfreq> &gw_, block_gf<imfreq> &fw_,
                       block_gf<imfreq> &sigmaw_,
                       block_gf_const_view<imfreq> g0w_)
    : params(params_), p(p_), config(config_), fprefactor(fprefactor_), gw(gw_),
      fw(fw_), sigmaw(sigmaw_), g0w(g0w_), beta(params->beta),
      n_w(gw[0].mesh().size()) {
  Z = 0;
  for (int i = 0; i < gw.size(); i++)
    gw[i]() = 0;
  for (int i = 0; i < fw.size(); i++)
    fw[i]() = 0;
}

/// accumulation of Green's function using NFFT
void measure_gw::accumulate(double s) {
  Z += s;

  if (p->use_nfft_for_gw) {
#ifdef HAVE_NFFT
    int n_w_meas = n_w + n_w % 2; // make sure n_w_meas is even

    for (int orb1 = 0; orb1 < config->n_colors(); orb1++) {
      int size = config->det(orb1).size();

      std::vector<dcomplex> coeff_f;
      if (params->measure_fw)
        coeff_f.resize(size * size); // buffer for coefficients of F

      // init a one dimensional plan
      nfft_plan plan;
      nfft_init_1d(&plan, n_w_meas, size * size);

      c_exp.resize(config->det(orb1).size());
      cdag_exp.resize(config->det(orb1).size());

      for (int i = 0; i < config->det(orb1).size();
           i++) { // precompute exponentials
        double c_time = double(config->det(orb1).get_y(i));
        double cdag_time = double(config->det(orb1).get_x(i));
        c_exp(i) = std::exp(dcomplex(0, M_PI * (1 + n_w_meas) * c_time / beta));
        cdag_exp(i) =
            std::exp(dcomplex(0, -M_PI * (1 + n_w_meas) * cdag_time / beta));
      } // i

      for (int j = 0; j < config->det(orb1).size(); j++) {
        double f_pref;
        if (params->measure_fw) {
          qmc_time_t tau_c =
              config->det(orb1).get_y(j); // time of the non-dagger
          f_pref = fprefactor->get(orb1, tau_c);
        }

        for (int i = 0; i < config->det(orb1).size(); i++) {

          qmc_time_t tau = config->det(orb1).get_y(j) -
                           config->det(orb1).get_x(i); // qmc_time wraps around
          (plan).x[i * size + j] =
              (double)tau / beta - 0.5; // init the nonequidistant times

          dcomplex coeff = s * config->det(orb1).inverse_matrix(j, i) *
                           c_exp(j) * cdag_exp(i);
          if (params->measure_fw)
            coeff_f[i * size + j] = coeff * f_pref;
          (plan).f[i * size + j][0] = coeff.real();
          (plan).f[i * size + j][1] = coeff.imag();
        }
      }

      /** precompute psi, the entries of the matrix B */
      if ((plan).nfft_flags &
          PRE_ONE_PSI) //(PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI| PRE_FULL_PSI)
        nfft_precompute_one_psi(&(plan));

      std::vector<dcomplex> prefactor(n_w_meas);
      for (int i = 0; i < n_w_meas; ++i)
        prefactor[i] = std::exp(dcomplex(0., M_PI * (i - n_w_meas / 2)));

      // approx. adjoint
      // nfft_check(&(plan));
      if (size < params->nfft_threshold)
        nfft_adjoint_direct(&(plan));
      else
        nfft_adjoint(&(plan));
      for (int n = 0; n < n_w; ++n) {
        dcomplex meas =
            dcomplex((plan).f_hat[n][0], (plan).f_hat[n][1]) * prefactor[n];
        gw[orb1][n](0, 0) += meas;
      }

      if (params->measure_fw) {
        for (int i = 0; i < config->det(orb1).size(); i++) {
          for (int j = 0; j < config->det(orb1).size(); j++) {
            (plan).f[i * size + j][0] = coeff_f[i * size + j].real();
            (plan).f[i * size + j][1] = coeff_f[i * size + j].imag();
          }
        }
        nfft_adjoint(&(plan));
        for (int n = 0; n < n_w; ++n) {
          dcomplex meas =
              dcomplex((plan).f_hat[n][0], (plan).f_hat[n][1]) * prefactor[n];
          fw[orb1][n](0, 0) += meas;
        }
      }
      // finalise the one dimensional plan
      nfft_finalize(&plan);

    } // orb1
#else
    TRIQS_RUNTIME_ERROR << "cthyb: measure_gw: NO NFFT INSTALLED! Exiting...";
#endif
  } else { // don't use NFFT
    for (int orb1 = 0; orb1 < gw.size(); orb1++) {
      c_exp.resize(config->det(orb1).size());
      cdag_exp.resize(config->det(orb1).size());
      c_inner_index.resize(config->det(orb1).size());
      cdag_inner_index.resize(config->det(orb1).size());

      for (int i = 0; i < config->det(orb1).size();
           i++) { // precompute exponentials
        auto c_time = config->det(orb1).get_y(i);
        auto cdag_time = config->det(orb1).get_x(i);
        c_exp(i) =
            std::exp(dcomplex(0, M_PI * double(std::get<0>(c_time)) / beta));
        c_inner_index(i) = std::get<1>(c_time);
        cdag_exp(i) = std::exp(
            dcomplex(0, -M_PI * double(std::get<0>(cdag_time)) / beta));
        cdag_inner_index(i) = std::get<1>(cdag_time);
      } // i

      for (int j = 0; j < config->det(orb1).size(); j++) {
        qmc_time_t tau;
        double f_pref;
        if (params->measure_fw) {
          auto y = config->det(orb1).get_y(j); // time of the non-dagger
          tau = std::get<0>(y);                // time of the non-dagger
          f_pref = fprefactor->get(std::get<2>(y), tau);
        }
        for (int i = 0; i < config->det(orb1).size(); i++) {
          dcomplex exp = c_exp(j) * cdag_exp(i);
          dcomplex incr = exp * exp;
          for (int n = 0; n < n_w / 2; ++n) { // fill only positive freqs
            dcomplex meas = s * config->det(orb1).inverse_matrix(j, i) * exp;
            gw[orb1][n](cdag_inner_index(i), c_inner_index(j)) += meas;
            if (params->measure_fw)
              fw[orb1][n](cdag_inner_index(i), c_inner_index(j)) +=
                  meas * f_pref;
            exp *= incr;
          }
        }
      }
    } // orb1
  }   // if use_nfft_for_gw
} // accumulate

// reduce and normalize G
void measure_gw::collect_results(mpi::communicator const &c) {
  Z = mpi::all_reduce(Z, c);
  gw = mpi::all_reduce(gw, c);
  gw = gw / (-beta * Z);

  for (int bl = 0; bl < gw.size(); bl++)
    for (auto const &w : gw[bl].mesh())
      if (w.n < 0)
        gw[bl][w] = conj(gw[bl](matsubara_freq(-w.n - 1, beta, Fermion)));

  if (params->measure_fw) {
    fw = mpi::all_reduce(fw, c);
    fw = fw / (-beta * Z);
    for (int bl = 0; bl < fw.size(); bl++)
      for (auto const &w : fw[bl].mesh())
        if (w.n < 0)
          fw[bl][w] = conj(fw[bl](matsubara_freq(-w.n - 1, beta, Fermion)));

    nda::clef::placeholder<0> orb_;
    nda::clef::placeholder<1> iom_;

    // sigmaw[orb_](iom_) << fw[orb_](iom_)/gw[orb_](iom_); //does not work
    for (int orb = 0; orb < gw.size(); orb++)
      for (auto &w : fw[orb].mesh())
        sigmaw[orb][w] = fw[orb](w) / gw[orb](w);
  } // end if measure_fw
  else {
    for (int orb = 0; orb < gw.size(); orb++)
      for (auto &w : fw[orb].mesh())
        sigmaw[orb][w] = 1 / g0w[orb](w) - 1 / gw[orb](w);
  }
}
} // namespace triqs_ctseg
