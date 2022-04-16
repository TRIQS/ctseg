#pragma once
#include "types.hpp"
#include "params.hpp"
#include "inputs.hpp"

#include <mpi/mpi.hpp>
#include <triqs/det_manip.hpp>

// ---------------------------------------------------
/// Working data
struct work_data_t {
  work_data_t(params_t const &p, inputs_t const &inputs, mpi::communicator c);

  nda::vector<double> mu;
  nda::matrix<double> U;
  gf_struct_t gf_struct;
  int n_color;

  bool has_Dt, has_jperp, rot_inv = true;
  // FIXME : real_valued ?? : remove real everywhere ?
  // Jperp : should be scalar ? or matrix ?
  gf<imtime> K, Kprime, Jperp, Kprime_spin;

  using delta_target_t = matrix_real_valued;
  block_gf<imtime, delta_target_t> delta; // Hybridization function
  std::vector<det_t> dets;                // The determinants

  int block_to_color(long const &block, long const &idx) const {
    int color = 0;
    for (int i = 0; i <= block; i++) {
      if (i == block)
        color += int(idx); // index within the block
      else
        color += this->gf_struct[i].second; // block size
    }
    return color;
  }
};
