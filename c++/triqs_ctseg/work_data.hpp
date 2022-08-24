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

  bool has_Dt    = false; // Have D nn term
  bool has_jperp = false; // Have Jperp
  bool rot_inv   = true;  // ???

  // FIXME : Could be real_valued if profiling show some gain
  gf<imtime> K, Kprime, Jperp, Kprime_spin;

  // Hybridization function
  block_gf<imtime, delta_target_t> delta;

  // The determinants
  std::vector<det_t> dets;

  // map (block, idx) of gf to a color
  int block_to_color(int block, int idx) const { return gf_block_size_partial_sum[block] + idx; }
};
