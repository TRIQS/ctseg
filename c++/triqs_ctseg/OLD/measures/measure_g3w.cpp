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
#include "./measure_g3w.hpp"

namespace triqs_ctseg {

/// constructor
measure_g3w::measure_g3w(const qmc_parameters *params_,
                         const configuration *config_,
                         std::shared_ptr<precompute_Mw> Mw_,
                         gf_3w_container_t &g3w_, gf_3w_container_t &f3w_)
    : params(params_), config(config_), Mw(Mw_), g3w(g3w_), f3w(f3w_),
      beta(params->beta),
      n_w_fermionic(std::get<0>(g3w[0].mesh().components()).last_index() + 1),
      n_w_bosonic(std::get<2>(g3w[0].mesh().components()).last_index() + 1) {
  Z = 0;
  for (auto &g : g3w)
    g() = 0;
  for (auto &g : f3w)
    g() = 0;
}

/// accumulate the Green's function
void measure_g3w::accumulate(double s) {
  Z += s;

  if (params->measure_g3w) {
    for (int bl = 0; bl < g3w.size(); bl++) { // bl : 'upup','updn',...
      int b1 = bl / config->gf_struct().size();
      int b2 = bl % config->gf_struct().size();
      for (int a = 0; a < g3w[bl].target_shape()[0]; a++) {
        for (int b = 0; b < g3w[bl].target_shape()[1]; b++) {
          for (int c = 0; c < g3w[bl].target_shape()[2]; c++) {
            for (int d = 0; d < g3w[bl].target_shape()[3]; d++) {
              for (int n1 = -n_w_fermionic; n1 < n_w_fermionic; n1++) {
                for (int n4 = -n_w_fermionic; n4 < n_w_fermionic; n4++) {
                  for (int m = 0; m < n_w_bosonic; m++) {
                    // set remaining frequencies
                    int n2 = n1 + m;
                    int n3 = n4 + m;
                    g3w[bl][{n1, n4, m}](a, b, c, d) +=
                        s * Mw->getM(b1, a, b, n1, n2) *
                        Mw->getM(b2, c, d, n3, n4);
                    if (b1 == b2)
                      g3w[bl][{n1, n4, m}](a, b, c, d) -=
                          s * Mw->getM(b1, a, d, n1, n4) *
                          Mw->getM(b2, c, b, n3, n2);
                  } // m
                }   // n4
              }     // n1
            }       // d
          }         // c
        }           // b
      }             // a
    }               // orb1
  }
  if (params->measure_f3w) {
    for (int bl = 0; bl < f3w.size(); bl++) { // bl : 'upup','updn',...
      int b1 = bl / config->gf_struct().size();
      int b2 = bl % config->gf_struct().size();
      for (int a = 0; a < f3w[bl].target_shape()[0]; a++) {
        for (int b = 0; b < f3w[bl].target_shape()[1]; b++) {
          for (int c = 0; c < f3w[bl].target_shape()[2]; c++) {
            for (int d = 0; d < f3w[bl].target_shape()[3]; d++) {
              for (int n1 = -n_w_fermionic; n1 < n_w_fermionic; n1++) {
                for (int n4 = -n_w_fermionic; n4 < n_w_fermionic; n4++) {
                  for (int m = 0; m < n_w_bosonic; m++) {
                    // set remaining frequencies
                    int n2 = n1 + m;
                    int n3 = n4 + m;
                    f3w[bl][{n1, n4, m}](a, b, c, d) +=
                        s * Mw->getnM(b1, a, b, n1, n2) *
                        Mw->getM(b2, c, d, n3, n4);
                    if (b1 == b2)
                      f3w[bl][{n1, n4, m}](a, b, c, d) -=
                          s * Mw->getnM(b1, a, d, n1, n4) *
                          Mw->getM(b2, c, b, n3, n2);
                  } // m
                }   // n4
              }     // n1
            }       // d
          }         // c
        }           // b
      }             // a
    }               // orb1
  }

  /*
  clef::placeholder<0> orb1_;
  clef::placeholder<1> orb2_;
  clef::placeholder<2> inu1_;
  clef::placeholder<3> inu2_;
  clef::placeholder<4> iom_;
  if(params->measure_g3w) {
    g3w(orb1_, orb2_)(inu1_,inu2_,iom_) << g3w(orb1_, orb2_)(inu1_,inu2_,iom_)
                                          + s * Mw->getM(orb1_, inu1_,
  inu1_+iom_) * Mw->getM(orb2_, inu2_+iom_, inu2_     )
                                          - s * Mw->getM(orb1_, inu1_, inu2_ ) *
  Mw->getM(orb2_, inu2_+iom_, inu1_+iom_) * kronecker(orb1_,orb2_);
  }
  if(params->measure_f3w){
    f3w(orb1_, orb2_)(inu1_,inu2_,iom_) << f3w(orb1_, orb2_)(inu1_,inu2_,iom_)
                                          + s * Mw->getnM(orb1_, inu1_,
  inu1_+iom_) * Mw->getM(orb2_, inu2_+iom_, inu2_     )
                                          - s * Mw->getnM(orb1_, inu1_, inu2_ )
  * Mw->getM(orb2_, inu2_+iom_, inu1_+iom_) * kronecker(orb1_,orb2_);
  }
 */
}

/// reduce and normalize G
void measure_g3w::collect_results(mpi::communicator const &c) {

  Z = mpi::reduce(Z, c);

  // has to be done later by bosonic frequencies, to minimize extra storage
  g3w = mpi::reduce(g3w, c);
  f3w = mpi::reduce(f3w, c);

  if (c.rank() == 0) {
    if (params->measure_g3w)
      g3w = g3w / (beta * Z);
    if (params->measure_f3w)
      f3w = f3w / (beta * Z);
  }
}

} // namespace triqs_ctseg
