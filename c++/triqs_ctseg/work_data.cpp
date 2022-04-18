#include "work_data.hpp"
#include <nda/basic_functions.hpp>
#include <nda/traits.hpp>
#include <triqs/gfs/functions/functions2.hpp>
#include <triqs/operators/util/extractors.hpp>
#include "logs.hpp"
#include "spdlog/common.h"
#include "spdlog/spdlog.h"

work_data_t::work_data_t(params_t const &p, inputs_t const &inputs, mpi::communicator c) {

  // Set logger level
  spdlog::set_pattern("%v");
  spdlog::set_level(spdlog::level::info);
  if constexpr (ctseg_debug) spdlog::set_level(spdlog::level::debug);

  double beta = p.beta;
  gf_struct   = p.gf_struct;
  n_color     = count_colors(gf_struct);

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

  // For numerical integration of the D0 and Jperp
  auto ramp = nda::zeros<double>(p.n_tau_k);
  for (auto n : range(p.n_tau_k)) { ramp(n) = n * beta / (p.n_tau_k - 1); }

  // Dynamical interactions
  if (has_Dt) {
    // Compute interaction kernels K(tau), K'(tau) by integrating D(tau)
    K      = gf<imtime>({beta, Boson, p.n_tau_k}, {n_color, n_color});
    Kprime = gf<imtime>({beta, Boson, p.n_tau_k}, {n_color, n_color});
    for (auto c1 : range(n_color)) {
      for (auto c2 : range(n_color)) {
        nda::array<dcomplex, 1> D_data = inputs.d0t.data()(range(), c1, c2);
        auto first_integral            = nda::zeros<dcomplex>(p.n_tau_k);
        auto second_integral           = nda::zeros<dcomplex>(p.n_tau_k);
        // Trapezoidal integration
        for (int i = 1; i < D_data.size(); ++i) {
          first_integral(i)  = first_integral(i - 1) + (D_data(i) + D_data(i - 1)) / 2;
          second_integral(i) = second_integral(i - 1) + (first_integral(i) + first_integral(i - 1)) / 2;
        }
        // Noramlize by bin size
        first_integral *= beta / (p.n_tau_k - 1);
        second_integral *= (beta / (p.n_tau_k - 1)) * (beta / (p.n_tau_k - 1));
        // Enforce K(0) = K(beta) = 0
        Kprime.data()(range(), c1, c2) = first_integral - second_integral(p.n_tau_k - 1) / beta;
        K.data()(range(), c1, c2)      = second_integral - ramp * second_integral(p.n_tau_k - 1) / beta;
        // Renormalize U and mu
        if (c1 != c2) U(c1, c2) -= real(2 * Kprime.data()(0, c1, c2));
      }
      mu(c1) += real(Kprime.data()(0, c1, c1));
    }
  }

  // J_perp interactions
  if (has_jperp) {
    Jperp = inputs.jperpt;
    if (not has_Dt)
      rot_inv = false;
    else {
      Kprime_spin = gf<imtime>({beta, Boson, p.n_tau_k}, {n_color, n_color}); // used in computation of F(tau)
      // Integrate Jperp to obtain the S_z.S_z part of K'(tau) (called Kprime_spin)
      auto Kprime_J                  = Jperp;
      nda::array<dcomplex, 1> J_data = Jperp.data()(range(), 0, 0);
      auto first_integral            = nda::zeros<dcomplex>(p.n_tau_k);
      // Trapezoidal integration
      for (int i = 1; i < J_data.size(); ++i) {
        first_integral(i) = first_integral(i - 1) + (J_data(i) + J_data(i - 1)) / 2;
      }
      // Noramlize by bin size
      first_integral *= beta / (p.n_tau_k - 1);
      // Enforce Kprime_J(beta/2) = 0
      Kprime_J.data()(range(), 0, 0) = first_integral - first_integral((p.n_tau_k - 1) / 2);
      // Kprime_spin = +/- Kprime_J depending on color
      for (auto c1 : range(n_color)) {
        for (auto c2 : range(n_color)) {
          Kprime_spin.data()(range(), c1, c2) = (c1 == c2 ? 1 : -1) * Kprime_J.data()(range(), 0, 0) / 4;
        }
      }
      auto Kprime_0 = gf<imtime>({beta, Boson, p.n_tau_k}, {n_color, n_color});
      Kprime_0      = Kprime - Kprime_spin;
      // The "remainder" Kprime_0 must be color-independent for there to be rotational invariance
      if (max_element(abs(Kprime_0.data()(range(), 0, 0) - Kprime_0.data()(range(), 0, 1))) > 1.e-13) rot_inv = false;
    }
  }

  // Report
  if (c.rank() == 0) {
    spdlog::info("mu = {} \nU = {}", mu, U);
    spdlog::info("Dynamical interactions = {}, J_perp interactions = {}", has_Dt, has_jperp);
    if (p.measure_ft and !rot_inv)
      spdlog::info("WARNING: Cannot measure F(tau) because spin-spin interaction is not rotationally invariant.");
  }

  // ................  Determinants .....................

  delta = map([](gf_const_view<imtime> d) { return real(d); }, inputs.delta);
  for (auto const &bl : range(delta.size())) {
    if (!is_gf_real(delta[bl], 1e-10)) {
      if (c.rank() == 0) {
        spdlog::info("WARNING: The Delta(tau) block number {} is not real in tau space", bl);
        spdlog::info("WARNING: max(Im[Delta(tau)]) = {}", max_element(abs(imag(delta[bl].data()))));
        spdlog::info("WARNING: Disregarding the imaginary component in the calculation.");
      }
    }
    dets.emplace_back(delta_block_adaptor{real(delta[bl])}, p.det_init_size);
    dets.back().set_singular_threshold(p.det_singular_threshold);
    dets.back().set_n_operations_before_check(p.det_n_operations_before_check);
    dets.back().set_precision_warning(p.det_precision_warning);
    dets.back().set_precision_error(p.det_precision_error);
  }
} // work_data constructor
