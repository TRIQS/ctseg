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
#include "./measure_gt.hpp"

namespace triqs_ctseg {

measure_gt::measure_gt(const qmc_parameters *params_,
                       const configuration *config_,
                       std::shared_ptr<precompute_fprefactor> fprefactor_,
                       block_gf<imtime> &gt_, block_gf<imtime> &ft_)
    : params(params_), config(config_), fprefactor(fprefactor_), gt(gt_),
      ft(ft_), beta(params->beta), gt_stack(), Z_stack(), counter(0),
      accum(0.0), Noverbeta(1.0 / gt[0].mesh().delta()) {
  Z = 0;
  for (int i = 0; i < gt.size(); i++)
    gt[i]() = 0;
  for (int i = 0; i < ft.size(); i++)
    ft[i]() = 0;
}

void measure_gt::accumulate(double s) {
  Z += s;
  counter += 1;

  // 3 implementations to compare speed.
  // All 3 give the same results ...
  // on clang (OS X), this is comparable, with PLAIN, LAMBDA2 equal and about
  // 10-20% faster code than LAMBDA. But LAMBDA will always be correct, even if
  // the grid type (half_bin, full_bin) of gt changes Origin of this difference
  // is partly due to the fact that fac is a const computed out of the loop.

//#define TRIQS_CTSEG_MEASURE_GT_WITH_LAMBDA
//#define TRIQS_CTSEG_MEASURE_GT_WITH_LAMBDA2
#define TRIQS_CTSEG_MEASURE_GT_PLAIN

#ifdef TRIQS_CTSEG_MEASURE_GT_WITH_LAMBDA
  for (int k = 0; k < config->n_colors(); k++) {
    foreach (config->det(k),
             [this, k, s](qmc_time_t const &x, qmc_time_t const &y, double M) {
               this->gt[k](closest_mesh_pt(double(y - x)))(0, 0) +=
                   (y > x ? s * M : -s * M);
               // this is not optimal in terms of speed, since map is searched
               // every time (compare to last variant)
               if (params->measure_ft) {
                 auto f_pref = this->fprefactor->get(k, y);
                 this->ft[k](closest_mesh_pt(double(y - x)))(0, 0) +=
                     (y > x ? s * f_pref * M : -s * f_pref * M);
               }
             })
      ;
  }
#endif

#ifdef TRIQS_CTSEG_MEASURE_GT_WITH_LAMBDA2
  // const double fac = 1.0/gt[0].mesh().delta();
  const qmc_time_t delta_tau =
      div_by_int(params->beta, gt[0].mesh().size()) +
      qmc_time_t::epsilon(beta); // careful, only true for half bins

  for (int k = 0; k < config->n_colors(); k++)
    foreach (config->det(k),
             [this, k, s, delta_tau](qmc_time_t const &x, qmc_time_t const &y,
                                     double M) {
               auto tau = y - x;
               if (y <= x)
                 M = -M;
               this->gt[k][floor_div(tau, delta_tau)](0, 0) += s * M;
               if (params->measure_ft) {
                 auto f_pref = this->fprefactor->get(k, y);
                 this->ft[k][floor_div(tau, delta_tau)](0, 0) +=
                     s * f_pref * M; // change for f
               }
             })
      ;

#endif

#ifdef TRIQS_CTSEG_MEASURE_GT_PLAIN

  accum = 0.0;
  auto closest_index = [this](auto p) {
    double x = double(p.value) + 0.5 * gt[0].mesh().delta();
    int n = std::floor(x / gt[0].mesh().delta());
    return n;
  };
  // loop over species and times of operators
  for (int block = 0; block < gt.size(); block++) {
    for (int j = 0; j < config->det(block).size(); j++) {
      double f_pref;
      auto y = config->det(block).get_y(j);
      if (params->measure_ft)
        f_pref = fprefactor->get(std::get<2>(y), std::get<0>(y));
      for (int i = 0; i < config->det(block).size(); i++) {
        auto [tau_x, x_inner, x_col] = config->det(block).get_x(i);
        auto [tau_y, y_inner, y_col] = config->det(block).get_y(j);
        double tau = double(tau_y - tau_x);
        int st = (tau_y >= tau_x) ? 1 : -1;
        double meas = s * config->det(block).inverse_matrix(j, i);
        gt[block][closest_mesh_pt(tau)](x_inner, y_inner) += st * meas;
        if (block == 0 and x_inner == 0 and y_inner == 0 and
            closest_index(closest_mesh_pt(tau)) == 0)
          accum += st * meas;
        if (params->measure_ft)
          ft[block][closest_mesh_pt(tau)](x_inner, y_inner) +=
              st * meas * f_pref;
      }
    }
  }
  gt_stack << accum / (-beta * 1. / Noverbeta);
  Z_stack << s;
#endif
}

void measure_gt::collect_results(mpi::communicator const &c) {

  Z = mpi::all_reduce(Z, c);
  gt = mpi::all_reduce(gt, c);
  gt = gt / (-beta * Z / Noverbeta);
  for (int k = 0; k < gt.size(); k++) {
    gt[k][0] *= 2;
    gt[k][gt[k].mesh().size() - 1] *= 2;
  }

  if (params->measure_ft) {
    ft = mpi::all_reduce(ft, c);
    ft = ft / (-beta * Z / Noverbeta);
    for (int k = 0; k < ft.size(); k++) {
      ft[k][0] *= 2;
      ft[k][ft[k].mesh().size() - 1] *= 2;
    }
  }
}

} // namespace triqs_ctseg
