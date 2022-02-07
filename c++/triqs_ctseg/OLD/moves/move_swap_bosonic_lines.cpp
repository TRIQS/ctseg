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
#include "./move_swap_bosonic_lines.hpp"
namespace triqs_ctseg {

move_swap_bosonic_lines::move_swap_bosonic_lines(
    const qmc_parameters *params_, configuration *config_,
    triqs::mc_tools::random_generator &RND_)
    : params(params_), config(config_), RND(RND_){};

//----------------------------------------------

double move_swap_bosonic_lines::attempt() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "\n =================== ATTEMPT SWAP BOSONIC LINES "
                 "================ \n"
              << std::endl;
#endif
  tried_swapping = false;
  if (config->boson_lines.size() < 2)
    return 0.0;
  int index1 = RND(config->boson_lines.size());
  int index2 = RND(config->boson_lines.size());
  if (index1 == index2)
    return 0.0;

  double r1, r2, r3, r4;
  if (index1 < index2)
    std::swap(index1, index2); // remove line with largest index first
  std::tie(r1, cseg1) = config->boson_lines.remove(index1);
  std::tie(r2, cseg2) = config->boson_lines.remove(index2);
  r3 = config->boson_lines.add(composite_segment{cseg1.plus, cseg2.minus});
  r4 = config->boson_lines.add(composite_segment{cseg2.plus, cseg1.minus});
  tried_swapping = true;
#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "attempting to swap bosonic lines #: " << index1 << "("
              << double(cseg1.plus.c->tau) << ", " << double(cseg1.minus.c->tau)
              << ") and " << index2 << "(" << double(cseg2.plus.c->tau) << ", "
              << double(cseg2.minus.c->tau) << ")" << std::endl;
    std::cerr << "\t r1: " << r1 << ", r2 = " << r2 << std::endl;
    std::cerr << "\t r3: " << r3 << ", r4 = " << r4 << std::endl;
    std::cerr << "\t SWAP proba : " << r3 * r4 * r1 * r2 << std::endl;
  }
#endif
  return r3 * r4 * r1 * r2;
}

//----------------------------------------------
double move_swap_bosonic_lines::accept() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> ACCEPT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
#endif
  return 1.0;
}

//----------------------------------------------
void move_swap_bosonic_lines::reject() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> REJECT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;
  if (tried_swapping) {
    config->boson_lines.pop_back();
    config->boson_lines.pop_back();
    config->boson_lines.add(cseg1);
    config->boson_lines.add(cseg2);
  }
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
#endif
}
} // namespace triqs_ctseg
