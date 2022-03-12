#pragma once
#include "types.hpp"
#include "params.hpp"
#include "inputs.hpp"

#include <mpi/mpi.hpp>
#include <triqs/det_manip.hpp>

using qmc_time_factory_t = triqs::utility::time_segment;

/// A lambda to adapt the Delta function for the call of the det.
struct delta_block_adaptor {
  gf<imtime, matrix_real_valued> delta; // make a copy. Needed in the real case anyway.

  double operator()(std::pair<qmc_time_t, int> const &x, std::pair<qmc_time_t, int> const &y) const {
    //det_scalar_t res = delta[closest_mesh_pt(double(x.first - y.first))](x.second, y.second);
    double res = delta(double(x.first - y.first))(x.second, y.second);
    return (x.first >= y.first ? res :
                                 -res); // x,y first are time_pt, wrapping is automatic in 
                                        // the - operation, but need to compute the sign
  }
};

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
  using det_t = triqs::det_manip::det_manip<delta_block_adaptor>;
  std::vector<det_t> dets; // The determinants
};
