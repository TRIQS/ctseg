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

  qmc_time_factory_t fac;
  qmc_time_t const qmc_beta = fac.get_upper_pt();
  qmc_time_t const qmc_zero = fac.get_lower_pt();

  nda::vector<double> mu;
  nda::matrix<double> U;

  bool has_Dt, has_jperp;
  gf<imtime> K, Kprime;

  // FIXME off diagonal delta ??
  using delta_target_t = matrix_real_valued;
  block_gf<imtime, delta_target_t> delta; // Hybridization function
  std::vector<det_t> dets;                // The determinants

  public:
  // Random time generation that excludes values at boundaries
  qmc_time_t make_random_time(triqs::mc_tools::random_generator &rng, qmc_time_t const &tau1, qmc_time_t const &tau2);
};
