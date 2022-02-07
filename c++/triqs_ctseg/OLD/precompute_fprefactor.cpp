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
#include "./precompute_fprefactor.hpp"

namespace triqs_ctseg {
precompute_fprefactor::precompute_fprefactor(const qmc_parameters *params_,
                                             const configuration *config_)
    : params(params_), config(config_), fprefactors(config->n_colors()) {
  if (params->jperp_interactions and !params->full_spin_rot_inv)
    TRIQS_RUNTIME_ERROR
        << "Improved estimator cannot be computed if Jperp non zero and the "
           "interactions in the spin sector are not rotationally invariant";
}

void precompute_fprefactor::operator()() {

  for (int orb1 = 0; orb1 < config->n_colors(); orb1++) {
    fprefactors[orb1].clear();
    auto hint = fprefactors[orb1].begin(); // hint for map insertion

    // set fprefactors for all annihilation operators in the fullop_map
    for (auto it = config->ops_map().begin(orb1);
         it != config->ops_map().end(orb1); ++it) {
      if (it->dagger)
        continue; // need fpref only for annihilation ops
      fprefactor_t fprefactor = 0.0;

      for (int orb2 = 0; orb2 < config->n_colors(); orb2++) {
        if (it->right_occupation()[orb2] and
            orb1 != orb2) // if occupied and different lines (U_ii=0)
          fprefactor +=
              params->U(orb1, orb2); // in the dyn case this is \tilde{U}

        if (params->dynamical_U) {
          double sgn2 = (orb1 == orb2) ? 1 : -1;
          if (it->right_occupation()[orb2]) { // if occupied
            if (orb1 == orb2)
              fprefactor +=
                  -2 * params->Kprime[closest_mesh_pt(0.0)](orb1, orb1)
                           .real(); // occupied anyway in the case orb1==orb2
            if (params->full_spin_rot_inv)
              fprefactor += -sgn2 * 4 *
                            params->Kperpprime[closest_mesh_pt(0.0)](0, 0)
                                .real(); // occupied anyway
          }

          // contribution from all operators on line orb2
          for (auto it2 = config->ops_map().begin(orb2);
               it2 != config->ops_map().end(orb2); ++it2) {
            qmc_time_t tau = it2->tau - it->tau;
            double sgn = (it2->dagger)
                             ? 1
                             : -1; // different sign for creator/annihilator
            fprefactor -= sgn * params->Kprime((double)tau)(orb1, orb2).real();
            if (params->full_spin_rot_inv) {
              fprefactor -=
                  sgn * 2 * sgn2 * params->Kperpprime((double)tau)(0, 0).real();
            }
          }
        } // end dynamical_U
      }   // end for orb2

      // insert and check
      int size = fprefactors[orb1].size(); // store 'old' size
      hint = fprefactors[orb1].insert(
          hint, fprefactormap_t::value_type(it->tau, fprefactor));
      // if size of map did not increase, element was already present; this
      // should not happen.
      if (fprefactors[orb1].size() == size)
        TRIQS_RUNTIME_ERROR
            << "insertion error while inserting f_prefactor element on line "
            << orb1;
    } // it
  }   // orb1
}

} // namespace triqs_ctseg
