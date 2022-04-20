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

  nda::vector<double> mu; // chemical potential per color
  nda::matrix<double> U;  // Uab n_a n_b
  gf_struct_t gf_struct;  // gf_struct of the Green function (input copied)
  int n_color;            // number of color

  bool has_Dt    = false; // Have D nn term
  bool has_jperp = false; // Have Jperp
  bool rot_inv   = true;  // ???

  // FIXME : real_valued ?? : remove real everywhere ?
  // Jperp : should be scalar ? or matrix ?
  gf<imtime> K, Kprime, Jperp, Kprime_spin;

  // Hybridization function
  using delta_target_t = matrix_real_valued;
  block_gf<imtime, delta_target_t> delta;

  // The determinants
  std::vector<det_t> dets;

  //
  int block_to_color(long block, long idx) const {
    int color = 0;
    for (int i = 0; i <= block; i++) {
      if (i == block)
        color += int(idx); // index within the block
      else
        color += gf_struct[i].second; // block size
    }
    return color;
  }
};
