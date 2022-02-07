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
#include "./precompute_Mw.hpp"

#ifdef HAVE_NFFT
#include "nfft3.h"
#endif

namespace triqs_ctseg {

precompute_Mw::precompute_Mw(const qmc_parameters *params_,
                             const solve_params_t *p_,
                             const configuration *config_,
                             std::shared_ptr<precompute_fprefactor> fprefactor_,
                             int n_w_fermionic_, int n_w_bosonic_)
    : params(params_), p(p_), config(config_), fprefactor(fprefactor_),
      beta(params->beta), n_w_fermionic(n_w_fermionic_),
      n_w_bosonic(n_w_bosonic_),
      n_w_aux(2 * n_w_fermionic + n_w_bosonic > 1
                  ? 2 * n_w_fermionic + n_w_bosonic - 1
                  : 0),
      w_ini((2 * (-n_w_fermionic) + 1) * M_PI / beta), // starting frequency
      w_inc(2 * M_PI / beta)                           // frequency increment
{
  Mw.resize(config->n_blocks());
  nMw.resize(config->n_blocks());
  int bl = 0;
  for (auto const &[bl_name, bl_size] : config->gf_struct()) {
    Mw[bl].resize(make_shape(bl_size, bl_size, n_w_aux, n_w_aux));
    nMw[bl].resize(make_shape(bl_size, bl_size, n_w_aux, n_w_aux));
    bl++;
  }
}

void precompute_Mw::operator()() {
  if (p->use_nfft_for_Mw) { // use nfft ...
#ifdef HAVE_NFFT //... but only if its installed, otherwise triqs error
    for (int bl = 0; bl < config->n_blocks(); bl++) {
      Mw[bl]() = 0.;
      nMw[bl]() = 0.;
    }

    int n_w_meas = n_w_aux + n_w_aux % 2; // make sure n_w_meas is even
    int n_shift = (n_w_bosonic - (n_w_bosonic % 2)) / 2;

    // compute M for all orbitals
    for (int bl = 0; bl < config->n_blocks(); bl++) {
      int size = config->det(bl).size();

      // nfft calculation
      nfft_plan nfft_c;
      nfft_plan nfft_cdagger;

      // init a one dimensional plan
      nfft_init_1d(&nfft_c, n_w_meas, size);
      nfft_init_1d(&nfft_cdagger, n_w_meas, size);

      c_exp.resize(config->det(bl).size());
      cdag_exp.resize(config->det(bl).size());

      // precompute exponentials
      for (int i = 0; i < config->det(bl).size(); i++) {
        double c_time = double(config->det(bl).get_y(i));
        double cdag_time = double(config->det(bl).get_x(i));
        c_exp(i) =
            std::exp(dcomplex(0, M_PI * (1 + 2 * n_shift) * c_time / beta));
        cdag_exp(i) =
            std::exp(dcomplex(0, -M_PI * (1 + 2 * n_shift) * cdag_time / beta));

        nfft_c.x[i] = c_time / beta - 0.5;
        nfft_cdagger.x[i] = -cdag_time / beta + 0.5;
      }

      // precompute psi, the entries of the matrix B
      if (nfft_c.nfft_flags & PRE_ONE_PSI) {
        nfft_precompute_one_psi(&nfft_c);
        nfft_precompute_one_psi(&nfft_cdagger);
      }

      std::vector<dcomplex> coeff_jm(
          size * n_w_meas, 0.); // the coefficients for the second FT step
      std::vector<int> freq(n_w_meas);
      for (int i = 0; i < n_w_meas; ++i)
        freq[i] = i - (n_w_meas) / 2;

      // first ft step: for each j and m do NFFT on i
      for (int j = 0; j < config->det(bl).size(); j++) {   // j
        for (int i = 0; i < config->det(bl).size(); i++) { // i
          dcomplex coeff = config->det(bl).inverse_matrix(j, i) * cdag_exp(i);
          nfft_cdagger.f[i][0] = coeff.real();
          nfft_cdagger.f[i][1] = coeff.imag();
        }
        // nfft_check(&nfft_cdagger);
        if (size < params->nfft_threshold)
          nfft_adjoint_direct(&nfft_cdagger);
        else
          nfft_adjoint(&nfft_cdagger);

        for (int n2 = 0; n2 < n_w_meas; ++n2) {
          dcomplex prefactor = std::exp(dcomplex(0., -M_PI * freq[n2]));
          coeff_jm[j * n_w_meas + n2] =
              dcomplex(nfft_cdagger.f_hat[n2][0], nfft_cdagger.f_hat[n2][1]) *
              prefactor * c_exp(j);
        }
      }
      // second ft step: for all m, perform NFFT on j
      for (int n2 = 0; n2 < n_w_aux;
           ++n2) { // only needed for n_w_aux instead of n_w_meas
        for (int j = 0; j < size; ++j) {
          nfft_c.f[j][0] = coeff_jm[j * n_w_meas + n2].real();
          nfft_c.f[j][1] = coeff_jm[j * n_w_meas + n2].imag();
        }
        // nfft_check(&nfft_c);
        if (size < params->nfft_threshold)
          nfft_adjoint_direct(&nfft_c);
        else
          nfft_adjoint(&nfft_c);

        for (int n1 = 0; n1 < n_w_aux;
             ++n1) { // only needed for N_w_aux instead of n_omega_meas
          dcomplex prefactor = std::exp(dcomplex(0., M_PI * freq[n1]));
          Mw[bl](0, 0, n1, n2) =
              dcomplex(nfft_c.f_hat[n1][0], nfft_c.f_hat[n1][1]) * prefactor;
        }
      }

      // so far no check whether the improved estimator is used
      // first ft step for F: for each j and m do NFFT on i
      for (int j = 0; j < config->det(bl).size(); j++) { // j
        qmc_time_t tau = config->det(bl).get_y(j);       // time of the dagger
        double f_pref = fprefactor->get(bl, tau);

        for (int i = 0; i < config->det(bl).size(); i++) { // i
          dcomplex coeff =
              config->det(bl).inverse_matrix(j, i) * cdag_exp(i) * f_pref;
          nfft_cdagger.f[i][0] = coeff.real();
          nfft_cdagger.f[i][1] = coeff.imag();
        }
        // nfft_check(&nfft_cdagger);
        if (size < params->nfft_threshold)
          nfft_adjoint_direct(&nfft_cdagger);
        else
          nfft_adjoint(&nfft_cdagger);

        for (int n2 = 0; n2 < n_w_meas; ++n2) {
          dcomplex prefactor = std::exp(dcomplex(0., -M_PI * freq[n2]));
          coeff_jm[j * n_w_meas + n2] =
              dcomplex(nfft_cdagger.f_hat[n2][0], nfft_cdagger.f_hat[n2][1]) *
              prefactor * c_exp(j);
        }
      }
      // second ft step: for all m, perform NFFT on j
      for (int n2 = 0; n2 < n_w_aux;
           ++n2) { // only needed for n_w_aux instead of n_w_meas
        for (int j = 0; j < size; ++j) {
          nfft_c.f[j][0] = coeff_jm[j * n_w_meas + n2].real();
          nfft_c.f[j][1] = coeff_jm[j * n_w_meas + n2].imag();
        }
        // nfft_check(&nfft_c);
        if (size < params->nfft_threshold)
          nfft_adjoint_direct(&nfft_c);
        else
          nfft_adjoint(&nfft_c);

        for (int n1 = 0; n1 < n_w_aux;
             ++n1) { // only needed for N_w_aux instead of n_omega_meas
          dcomplex prefactor = std::exp(dcomplex(0., M_PI * freq[n1]));
          nMw[bl](0, 0, n1, n2) =
              dcomplex(nfft_c.f_hat[n1][0], nfft_c.f_hat[n1][1]) * prefactor;
        }
      }

      nfft_finalize(&nfft_c); // deallocation
      nfft_finalize(&nfft_cdagger);

    } // bl
#else
    TRIQS_RUNTIME_ERROR
        << "cthyb: precompute_Mw: NO NFFT INSTALLED! Exiting...";
#endif
  } else // don't use nfft
  {
    // call operator which performs the actual computation
    for (int bl = 0; bl < config->n_blocks(); bl++) {
      Mw[bl]() = 0.;
      nMw[bl]() = 0.;
    }

    // compute M for all orbitals
    for (int block = 0; block < config->n_blocks(); block++) {
      c_exp_ini.resize(config->det(block).size());
      c_exp_inc.resize(config->det(block).size());
      cdag_exp_ini.resize(config->det(block).size());
      cdag_exp_inc.resize(config->det(block).size());
      c_inner_index.resize(config->det(block).size());
      cdag_inner_index.resize(config->det(block).size());

      // precompute exponentials
      for (int i = 0; i < config->det(block).size(); i++) {
        auto c_time = config->det(block).get_y(i); // a triplet
        auto cdag_time = config->det(block).get_x(i);
        c_exp_ini(i) =
            std::exp(dcomplex(0, w_ini * double(std::get<0>(c_time))));
        c_exp_inc(i) =
            std::exp(dcomplex(0, w_inc * double(std::get<0>(c_time))));
        cdag_exp_ini(i) =
            std::exp(dcomplex(0, -w_ini * double(std::get<0>(cdag_time))));
        cdag_exp_inc(i) =
            std::exp(dcomplex(0, -w_inc * double(std::get<0>(cdag_time))));

        c_inner_index(i) = std::get<1>(c_time);
        cdag_inner_index(i) = std::get<1>(cdag_time);
      }

      // set values of Mw
      for (int j = 0; j < config->det(block).size(); j++) { // j
        auto y = config->det(block).get_y(j);               // annihilation
        // so far no check whether the improved estimator is used
        qmc_time_t tau = std::get<0>(y); // time of the annihilation
        int c = std::get<2>(y);          // color
        double f_pref = fprefactor->get(c, tau);
        int cj = c_inner_index(j);

        for (int i = 0; i < config->det(block).size(); i++) { // i
          int cdagi = cdag_inner_index(i);

          dcomplex c_exp = c_exp_ini(j); // do not move to outer loop
          dcomplex cdag_exp = cdag_exp_ini(i);

          dcomplex Mji = config->det(block).inverse_matrix(j, i);

          for (int n1 = 0; n1 < n_w_aux; ++n1) {
            for (int n2 = 0; n2 < n_w_aux; ++n2) {
              dcomplex meas = Mji * c_exp * cdag_exp;

              Mw[block](cj, cdagi, n1, n2) += meas; // put back after debugging
              nMw[block](cj, cdagi, n1, n2) += f_pref * meas;
              cdag_exp *= cdag_exp_inc(i); // update cdag_exp
            }                              // n2
            cdag_exp = cdag_exp_ini(i);    // reset cdag_exp
            c_exp *= c_exp_inc(j);         // update c_exp
          }
        } // i
      }   // j
    }     // block

  } // if use nfft
} // operator ()
} // namespace triqs_ctseg
