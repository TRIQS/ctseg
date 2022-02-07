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
#include "./measure_g2w.hpp"

namespace triqs_ctseg {

  measure_g2w::measure_g2w(const qmc_parameters *params_, const configuration *config_, std::shared_ptr<precompute_Mw> Mw_,
                           std::shared_ptr<precompute_nw> nw_, block_f_om_nu_tv_t &g2w_, block_f_om_nu_tv_t &f2w_)
     : params(params_),
       config(config_),
       Mw(Mw_),
       nw(nw_),
       g2w(g2w_),
       f2w(f2w_),
       beta(params->beta),
       n_w_fermionic(std::get<0>(g2w[0].mesh().components()).last_index() + 1),
       n_w_bosonic(std::get<1>(g2w[0].mesh().components()).last_index() + 1) {
    Z = 0;
    for (auto &g : g2w) g() = 0;
    for (auto &g : f2w) g() = 0;
    nw->minimum_size(n_w_bosonic);
  }

void measure_g2w::accumulate(double s) {
  Z += s;
  if (params->measure_g2w) {
    for (int bl = 0; bl < g2w.size(); bl++) { // bl : 'upup','updn',...
      int b1 = bl / config->gf_struct().size();
      int b2 = bl % config->gf_struct().size();
      for (int c = 0; c < g2w[bl].target_shape()[2]; c++) {
        auto col = config->block_and_inner_index_to_color(b2, c);
        for (int a = 0; a < g2w[bl].target_shape()[0]; a++) {
          for (int b = 0; b < g2w[bl].target_shape()[1]; b++) {
            for (int m = 0; m < n_w_bosonic; m++) {
              for (int n1 = -(m - 1) / 2 - 1; n1 < n_w_fermionic; n1++) {
                int n2 = n1 + m; // set remaining frequency
                g2w[bl][{n1, m}](a, b, c) -=
                    s * Mw->getM(b1, a, b, n1, n2) * nw->get(col, m);
              } // m
            }   // n1
          }     // c
        }       // b
      }         // a
    }           // bl
  }
  if (params->measure_f2w) {
    for (int bl = 0; bl < f2w.size(); bl++) { // bl : 'upup','updn',...
      int b1 = bl / config->gf_struct().size();
      int b2 = bl % config->gf_struct().size();
      for (int c = 0; c < f2w[bl].target_shape()[2]; c++) {
        auto col = config->block_and_inner_index_to_color(b2, c);
        for (int a = 0; a < f2w[bl].target_shape()[0]; a++) {
          for (int b = 0; b < f2w[bl].target_shape()[1]; b++) {
            for (int n1 = -n_w_fermionic; n1 < n_w_fermionic; n1++) {
              for (int m = 0; m < n_w_bosonic; m++) {
                int n2 = n1 + m; // set remaining frequency
                f2w[bl][{n1, m}](a, b, c) -=
                    s * Mw->getnM(b1, a, b, n1, n2) * nw->get(col, m);
              } // m
            }   // n1
          }     // c
        }       // b
      }         // a
    }           // bl
  }

  g2w_0_0_stack << (-s * Mw->getM(0, 0, 0, 0, 0) * nw->get(0, 0)).real() / beta;
  g2w_10_0_stack << (-s * Mw->getM(0, 0, 0, 10, 10) * nw->get(0, 0)).real() /
                        beta;
  g2w_m10_0_stack << (-s * Mw->getM(0, 0, 0, -10, -10) * nw->get(0, 0)).real() /
                         beta;
  /*
     if(params->measure_g2w) {
     g2w(orb1, orb2)(inu_, iom_) = g2w(orb1_, orb2_)(inu_, iom_) - s *
     Mw->getnM(orb1_, inu_, inu+iom_) * nw->get(orb2_, iom_);
     }

     if(params->measure_f2w) {
     f2w(orb1, orb2)(inu_, iom_) = f2w(orb1_, orb2_)(inu_, iom_) - s *
     Mw->getnM(orb1_, inu_, inu+iom_) * nw->get(orb2_, iom_);
     }
   */
}

void measure_g2w::collect_results(mpi::communicator const &c) {
  Z = mpi::all_reduce(Z, c);

  for (int bl = 0; bl < g2w.size(); bl++) {
    g2w[bl] = mpi::all_reduce(g2w[bl], c);
    f2w[bl] = mpi::all_reduce(f2w[bl], c);
    g2w[bl] = g2w[bl] / (beta * Z);
    f2w[bl] = f2w[bl] / (beta * Z);
    // fill the rest by symmetry
    for (int m = 0; m < n_w_bosonic; m++)
      for (int n1 = -n_w_fermionic; n1 < -(m - 1) / 2; n1++)
        g2w[bl][{n1, m}] = conj(g2w[bl][{-m - n1 - 1, m}]);
  } // for bl

  auto autocorrelation_time_from_binning = [&c](auto const &acc) {
    auto [errs, counts] = acc.log_bin_errors_all_reduce(c);
    if (errs.size() == 0) return double{NAN};
    return std::max(0.0, tau_estimate_from_errors(errs[int(0.7 * errs.size())], errs[0]));
  };

  /// print out error bars for some components
  try {
    std::cout << "rank \t n,m \t average +/- error \t autocorrelation_time"
              << std::endl;
    std::cout << c.rank() << "\t 0,0 \t "
              << mean_and_err_mpi(c, g2w_0_0_stack.linear_bins()) << "\t"
              << autocorrelation_time_from_binning(g2w_0_0_stack) << std::endl;
    std::cout << c.rank() << "\t 10,0 \t "
              << mean_and_err_mpi(c, g2w_10_0_stack.linear_bins()) << "\t"
              << autocorrelation_time_from_binning(g2w_10_0_stack) << std::endl;
    std::cout << c.rank() << "\t -10,0 \t "
              << mean_and_err_mpi(c, g2w_m10_0_stack.linear_bins()) << "\t"
              << autocorrelation_time_from_binning(g2w_m10_0_stack)
              << std::endl;
  } catch (triqs::exception const &e) {
    std::cerr << "Warning: " << e.what() << std::endl;
  }
}

} // namespace triqs_ctseg
