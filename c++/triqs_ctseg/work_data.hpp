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

  int n_color;
  double beta;

  nda::vector<double> mu;
  nda::matrix<double> U;

  bool has_Dt, has_jperp;
  // FIXME : real_valued ?? : remove real everywhere ?
  // Jperp : should be scalar ? or matrix ?
  gf<imtime> K, Kprime, Jperp;

  // FIXME off diagonal delta ??
  using delta_target_t = matrix_real_valued;
  block_gf<imtime, delta_target_t> delta; // Hybridization function
  std::vector<det_t> dets;                // The determinants
};
