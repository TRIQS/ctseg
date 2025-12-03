// Copyright (c) 2022--present, The Simons Foundation
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./nn_nu.hpp"
#include "../logs.hpp"
#include <triqs/mesh.hpp>

namespace triqs_ctseg::measures {

  nn_nu::nn_nu(params_t const &p, work_data_t const &wdata, configuration_t const &config, results_t &results)
     : wdata{wdata}, config{config}, results{results} {

    beta           = p.beta;
    ntau           = p.n_tau_chi2;
    dtau           = p.beta / (ntau - 1);
    n_color        = config.n_color();
    block_number   = wdata.block_number;
    index_in_block = wdata.index_in_block;

    auto m = dlr_imfreq(p.beta, Boson, p.dlr_omega_max, p.dlr_epsilon);
    std::cout << "Number of DLR frequencies for nn_nu measurement : " << m.size() << std::endl;
    q_nu_block = make_block2_gf<dlr_imfreq>(m, p.gf_struct);
    q_nu       = gf<dlr_imfreq>(m, {n_color, n_color});
    q_nu()     = 0;
  }

  // -------------------------------------

  void nn_nu::accumulate(double s) {

    LOG("\n =================== MEASURE n_n(i nu)  ================ \n");

    Z += s;

    // Fourier transform of the characteristic function on [a,b], for n!=0
    auto ksi1 = [&](double a, double b, matsubara_freq const &nu) -> dcomplex {
      EXPECTS(b > a);
      if (nu.n == 0) return (b - a);
      dcomplex inu = nu;
      return (std::exp(inu * b) - std::exp(inu * a)) / inu;
    };

    // Fourier transform of the characteristic function on the segment seg = [tau_c, tau_cdag], accounting for cyclicity
    auto ksi = [&](segment_t const &seg, matsubara_freq const &nu) -> dcomplex {
      if (not is_cyclic(seg))
        return ksi1(double(seg.tau_cdag), double(seg.tau_c), nu);
      else
        return ksi1(double(seg.tau_cdag), beta, nu) + ksi1(0, double(seg.tau_c), nu);
    };

    // Compute n_a(nu)
    auto n_a_nu = [&](long a, matsubara_freq const &nu) -> dcomplex {
      dcomplex sum = 0.0;
      for (auto const &seg : config.seglists[a]) sum += ksi(seg, nu);
      return sum;
    };

    // <n_a(i nu) n_b(-i nu) >
    auto R = nda::range(n_color);
    nda::array<dcomplex, 1> n(n_color);
    for (auto const &nu : q_nu.mesh()) { // loop over bosonic frequencies
      for (auto a : R) { n(a) = n_a_nu(a, nu); }
      for (auto a : R)
        for (auto b : R) q_nu[nu](a, b) += s * n(a) * std::conj(n(b));
    }
  }

  // -------------------------------------

  void nn_nu::collect_results(mpi::communicator const &c) {

    Z = mpi::all_reduce(Z, c);

    q_nu = mpi::all_reduce(q_nu, c);
    q_nu = q_nu / Z / beta;

    // store the result
    for (int c1 : range(n_color)) {
      for (int c2 : range(n_color)) {
        q_nu_block(block_number[c1], block_number[c2]).data()(range::all, index_in_block[c1], index_in_block[c2]) =
           q_nu.data()(range::all, c1, c2);
      }
    }
    results.nn_nu = std::move(q_nu_block);
  }

} // namespace triqs_ctseg::measures
