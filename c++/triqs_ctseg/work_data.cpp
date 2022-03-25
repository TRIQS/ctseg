#include "work_data.hpp"
#include <nda/basic_functions.hpp>
#include <triqs/gfs/functions/functions2.hpp>
#include <triqs/operators/util/extractors.hpp>
#include "logs.hpp"
#include "spdlog/common.h"
#include "spdlog/spdlog.h"

work_data_t::work_data_t(params_t const &p, inputs_t const &inputs, mpi::communicator c) : fac{p.beta} {

  // Set logger level
  spdlog::set_level(spdlog::level::info);
#ifdef EXT_DEBUG
  spdlog::set_level(spdlog::level::debug);
#endif
#ifdef PRINT_CONFIG
  spdlog::set_level(spdlog::level::trace);
#endif

  beta = p.beta;

  // Number of colors from Green's function structure
  n_color = 0;
  for (auto const &[bl_name, bl_size] : p.gf_struct) { n_color += bl_size; }

  // Color dependent chemical potential
  mu = nda::zeros<double>(n_color);

  if (p.hartree_shift.size() > 0) {
    ALWAYS_EXPECTS((p.hartree_shift.size() == n_color), "Hartree shift size is not {}", n_color);
    mu += p.hartree_shift;
  }

  // .............. Interactions .................
  // Extract the U from the operator
  auto U_full = triqs::operators::utils::dict_to_matrix(triqs::operators::utils::extract_U_dict2(p.h_int), p.gf_struct);
  U           = nda::matrix<double>{real(U_full)};

  // Do we have D(tau) and J_perp(tau)? Yes, unless the data is 0
  has_Dt    = max_element(abs(inputs.d0t.data())) > 1.e-13;
  has_jperp = max_element(abs(inputs.jperpt.data())) > 1.e-13;

  // Check: no J_perp or D(tau) implementation for more than 2 colors
  if (n_color != 2) {
    ALWAYS_EXPECTS((not has_jperp), "Error : has_jperp is true and we have {} colors instead of 2", n_color);
    ALWAYS_EXPECTS((not has_Dt), "Error : has_Dt is true and we have {} colors instead of 2", n_color);
  }
  // Report
  if (c.rank() == 0) {
    spdlog::info("mu = {}\n U = {}", mu, U);
    spdlog::info("dynamical_U = {}\n jperp_interactions = {}\n ", has_Dt, has_jperp);
    //spdlog::info("Spdlog level is {}", spdlog::get_level());
  }

  if (has_Dt) {
    // Compute interaction kernels K(tau), K'(tau) by integrating D(tau)
    K                    = gf<imtime>({beta, Boson, p.n_tau_k}, {n_color, n_color});
    Kprime               = gf<imtime>({beta, Boson, p.n_tau_k}, {n_color, n_color});
    auto D_data          = inputs.d0t.data()(0, 0, range());
    auto first_integral  = nda::zeros<dcomplex>(p.n_tau_k);
    auto second_integral = nda::zeros<dcomplex>(p.n_tau_k);
    std::partial_sum(D_data.begin(), D_data.end(), first_integral.begin());
    std::partial_sum(first_integral.begin(), first_integral.end(), second_integral.begin());
    first_integral *= beta / (p.n_tau_k - 1);
    second_integral *= (beta / (p.n_tau_k - 1)) * (beta / (p.n_tau_k - 1));
    Kprime.data()(0, 1, range()) = first_integral - second_integral(p.n_tau_k - 1) / beta;
    Kprime.data()(1, 0, range()) = Kprime.data()(0, 1, range());
    auto ramp                    = nda::zeros<double>(p.n_tau_k);
    for (auto n : range(p.n_tau_k)) { ramp(n) = n * beta / (p.n_tau_k - 1); }
    K.data()(0, 1, range()) = second_integral - second_integral(0) - ramp * second_integral(p.n_tau_k - 1) / beta;
    K.data()(1, 0, range()) = K.data()(0, 1, range());

    // Renormalize U and mu
    U -= real(2 * Kprime(0));
    mu += real(Kprime.data()(0, 1, 0)); // FIXME: true?
  }

  // .............  Determinants .....................

  delta = map([](gf_const_view<imtime> d) { return real(d); }, inputs.delta);
  if (c.rank() == 0) {
    for (auto const &bl : range(delta.size())) {
      if (!is_gf_real(delta[bl], 1e-10)) {
        spdlog::info("WARNING: The Delta(tau) block number {} is not real in tau space\n", bl);
        spdlog::info("WARNING: max(Im[Delta(tau)]) = {} \n", max_element(abs(imag(delta[bl].data()))));
        spdlog::info("WARNING: Disregarding the imaginary component in the calculation.\n");
      }
      dets.emplace_back(delta_block_adaptor{real(delta[bl])}, p.det_init_size);
      dets.back().set_singular_threshold(p.det_singular_threshold);
      dets.back().set_n_operations_before_check(p.det_n_operations_before_check);
      dets.back().set_precision_warning(p.det_precision_warning);
      dets.back().set_precision_error(p.det_precision_error);
    }
  }
} // work_data constructor

// Random time generation that excludes values at boundaries, accounts for cyclicity
qmc_time_t work_data_t::make_random_time(triqs::mc_tools::random_generator &rng, qmc_time_t const &tau1,
                                         qmc_time_t const &tau2) {
  auto l   = tau2 - tau1;
  auto tau = tau1 + fac.get_random_pt(rng, qmc_zero, l);
  auto dt  = fac.get_epsilon();
  if (tau == tau1) {
    if (tau1 < tau2)
      return tau1 + dt;
    else
      return tau1 - dt;
  }
  if (tau == tau2) {
    if (tau1 < tau2)
      return tau2 - dt;
    else
      return tau2 + dt;
  }
  return tau;
}
