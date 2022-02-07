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
#include "./measure_gl.hpp"

namespace triqs_ctseg {

// constructor
measure_gl::measure_gl(const qmc_parameters *params_,
                       const configuration *config_,
                       std::shared_ptr<precompute_fprefactor> fprefactor_,
                       g_l_t &gl_, g_l_t &fl_)
    : params(params_), config(config_), fprefactor(fprefactor_), gl(gl_),
      fl(fl_), beta(params->beta), n_l(gl[0].mesh().size()), Tn() {
  Z = 0;
  for (int i = 0; i < gl.size(); i++)
    gl[i]() = 0;
  for (int i = 0; i < fl.size(); i++)
    fl[i]() = 0;
}

// accumulate the Green's function
void measure_gl::accumulate(double s) {

  Z += s;

  // this measures Legendre coefficients of +<T c_1(tau)c^*_1(tau')> ; overall
  // sign is taken care of below

  // 3 implementations to compare speed.
  // All 3 give the same results ...
  // on clang (OS X), this is comparable, with PLAIN, LAMBDA2 equal and about
  // 10-20% faster code than LAMBDA. But LAMBDA will always be correct, even if
  // the grid type (half_bin, full_bin) of gt changes Origin of this difference
  // is partly due to the fact that fac is a const computed out of the loop.

//#define TRIQS_CTSEG_MEASURE_GL_WITH_LAMBDA
//#define TRIQS_CTSEG_MEASURE_GL_WITH_LAMBDA2
#define TRIQS_CTSEG_MEASURE_GL_PLAIN

#ifdef TRIQS_CTSEG_MEASURE_GL_WITH_LAMBDA
  for (int k = 0; k < config->n_colors(); k++) {
    foreach (config->det(k),
             [this, k, s](qmc_time_t const &x, qmc_time_t const &y, double M) {
               // this is not optimal in terms of speed, since map is searched
               // every time (compare to last variant)
               double f_pref;
               if (params->measure_fl)
                 f_pref = this->fprefactor->get(k, y);
               auto tau = y - x; // folds back to [0,beta]
               double arg = 2 * (double)tau / beta - 1;
               Tn.reset(arg);
               for (auto l : gl[k].mesh()) {
                 auto tn = Tn.next();
                 this->gl[k][l](0, 0) += (y > x ? s * M * tn : -s * M * tn);
                 if (params->measure_fl)
                   this->fl[k][l](0, 0) +=
                       (y > x ? s * f_pref * M * tn : -s * f_pref * M * tn);
               }
             })
      ;
  }

#endif

#ifdef TRIQS_CTSEG_MEASURE_GL_WITH_LAMBDA2
  for (int k = 0; k < config->n_colors(); k++)
    foreach (config->det(k),
             [this, k, s](qmc_time_t const &x, qmc_time_t const &y, double M) {
               double f_pref;
               if (params->measure_fl)
                 f_pref = this->fprefactor->get(k, y);
               auto tau = y - x; // folds back to [0,beta]
               if (y <= x)
                 M = -M;
               double arg = 2 * (double)tau / beta - 1;
               Tn.reset(arg);
               for (auto l : gl[k].mesh()) {
                 auto tn = Tn.next();
                 this->gl[k][l](0, 0) += s * M * tn;
                 if (params->measure_fl)
                   this->fl[k][l](0, 0) += f_pref * s * M * tn;
               }
             })
      ;
#endif

#ifdef TRIQS_CTSEG_MEASURE_GL_PLAIN
  // loop over species and times of operators
  for (int k = 0; k < gl.size(); k++) {
    for (int i = 0; i < config->det(k).size(); i++) {
      auto y = config->det(k).get_y(i);
      double f_pref;
      if (params->measure_fl)
        f_pref = fprefactor->get(std::get<2>(y), std::get<0>(y));
      for (int j = 0; j < config->det(k).size(); j++) {
        auto [tau_x, x_inner, x_col] = config->det(k).get_x(i);
        auto [tau_y, y_inner, y_col] = config->det(k).get_y(j);
        double tau = double(tau_y - tau_x);
        double sym_sign = 1;
        if (tau < 0) {
          tau += beta;
          sym_sign = -1;
        }
        const double meas = s * sym_sign * config->det(k).inverse_matrix(j, i);
        const double x_l = 2 * tau / beta - 1;
        Tn.reset(x_l);

        for (auto l : gl[k].mesh()) {
          auto tn = Tn.next();
          gl[k][l](x_inner, y_inner) +=
              meas *
              tn; // replace with access on meshpoint and bin function later
          if (params->measure_fl)
            fl[k][l](x_inner, y_inner) +=
                f_pref * meas *
                tn; // replace with access on meshpoint and bin function later
        }
      }
    }
  }
#endif
}

// reduce and normalize G
void measure_gl::collect_results(mpi::communicator const &c) {

  Z = mpi::all_reduce(Z, c);

  gl = mpi::all_reduce(gl, c);

  using std::sqrt;
  for (int k = 0; k < config->n_colors(); k++)
    for (auto l : gl[k].mesh())
      gl[k][l](0, 0) = sqrt(2.0 * l + 1) / (-beta * Z) * gl[k][l](0, 0);

  // gl[orb_][l_](0,0) << sqrt(2.0*l_+1)/(-beta * Z) * gl[orb_][l_](0,0);

  if (params->measure_fl) {
    fl = mpi::all_reduce(fl, c);
    for (int k = 0; k < config->n_colors(); k++)
      for (auto l : fl[k].mesh())
        fl[k][l](0, 0) = sqrt(2.0 * l + 1) / (-beta * Z) * fl[k][l](0, 0);
  }
}

} // namespace triqs_ctseg
