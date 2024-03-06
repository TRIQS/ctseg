#pragma once
#include <mpi/mpi.hpp>
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/det_manip.hpp>

#include "params.hpp"
#include "inputs.hpp"
#include "util.hpp"
#include "dets.hpp"

struct work_data_t {

  work_data_t(params_t const &p, inputs_t const &inputs, mpi::communicator c);

  nda::matrix<double> U;  // Uab n_a n_b
  gf_struct_t gf_struct;  // gf_struct of the Green function (input copied)
  int n_color;            // number of color
  nda::vector<double> mu; // chemical potential per color

  std::vector<long> gf_block_size_partial_sum; // data for block_to_color method

  bool has_delta = false; // There is a non-zero hybridization term
  bool has_Dt    = false; // Have D nn term
  bool has_jperp = false; // Have Jperp
  bool rot_inv   = true;  // Is the spin-spin interaction rotationally invariant? 
  bool minus_sign = false; // Has a move ever produced a negative sign?
  bool offdiag_delta = false; // Does Delta(tau) have blocks of size larger than 1? 

  // FIXME : Could be real_valued if profiling show some gain
  gf<imtime> K, Kprime, Jperp, Kprime_spin;

  // Hybridization function
  using delta_target_t = matrix_real_valued;
  block_gf<imtime, delta_target_t> delta;

  // The determinants
  std::vector<det_t> dets;

  // Map (block, idx) of gf to a color
  int block_to_color(int block, int idx) const { return gf_block_size_partial_sum[block] + idx; }

  // Color to (block, idx) conversion tables
  std::vector<long> block_number; // block numbers corresponding to colors
  std::vector<long> index_in_block; // index in block of a given color

  // Find block of color
  long find_block_number(int color) const {
    long bl = 0;
    long colors_so_far = 0;
    for (auto const &[s, l] : gf_struct) {
      colors_so_far += l;
      if (color < colors_so_far) {
        return bl;
      }
      bl++;
    }
    ALWAYS_EXPECTS((colors_so_far == n_color), "Error in color-to-block conversion.");
    return 0;
  }

  // Find index of color in block
  long find_index_in_block(int color) const {
    long colors_so_far = 0;
    for (auto const &[s, l] : gf_struct) {
      colors_so_far += l;
      if (color < colors_so_far) { return color - (colors_so_far - l); }
    }
    ALWAYS_EXPECTS((colors_so_far == n_color), "Error in color-to-block conversion.");
    return 0;
  }
  
};
