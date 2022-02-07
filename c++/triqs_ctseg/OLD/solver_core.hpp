/*******************************************************************************
 * CTSeg: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet Copyright (C) 2019 Simons Foundation author: N. Wentzell
 *
 * CTSEG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * CTSEG is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * CTSEG. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#pragma once

#include "params.hpp"
#include "qmc_parameters.hpp"
#include "types.hpp"
#include "block_matrix.hpp"

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/operators/util/extractors.hpp>

#include <optional>

namespace triqs_ctseg {

/// Main solver class
/**
  *
  Worker which runs the quantum Monte-Carlo simulation.
*/
class solver_core {

  double percent_done_, average_sign_;
  double beta;
  int n_color;
  gf_struct_t gf_struct;
  int n_blocks;
  std::vector<std::string> block_names;
  std::vector<int> block_sizes;

  // Green's function containers

  block_gf<imtime> delta_raw;
  block_gf<imfreq> g0w;
  // blocks structure: (N_blocks)^2 (e.g if G is (up,dn), then (up|up, up|dn,
  // dn|up, dn|dn) )
  block_gf<imfreq> d0w; // retarded interactions
  block_gf<imtime> kt;  // retarded interaction kernel
  block_gf<imtime> kprimet;

  gf<imtime> kperpprimet; // derivative of retarded interaction kernel
  gf<imfreq> jperpw;      // perp spin retarded interaction kernel
  gf<imtime> jperpt;      // perp spin retarded interaction kernel
  gf<imtime> chipmt;      // chi-plusminus spinspin susceptibility

  // imaginary-time Green's functions
  block_gf<imtime> gt;
  block_gf<imtime> ft;
  block_gf<imtime> nnt;

  // Matsubara Green's functions
  block_gf<imfreq> gw;
  block_gf<imfreq> fw;
  block_gf<imfreq> sw;
  block_gf<imfreq> nnw;

  // Legendre Green's functions
  g_l_t gl;
  g_l_t fl;

  // two-frequency (three-leg) Green's functions
  block_f_om_nu_tv_t g2w;
  block_f_om_nu_tv_t f2w;
  //
  // three-frequency (four-leg) Green's functions
  gf_3w_container_t g3w;
  gf_3w_container_t f3w;
  //
  //  array<gf_2w_t, 1> lambda_spin, lambda_charge ;

  // histogram
  matrix<double> hist;
  matrix<double> hist_composite;

  // density-density matrix
  block_matrix<double> nn_matrix;

  // histogram of impurity states (sector statistics)
  vector<double> state_hist;

  // define the communicator, here MPI_COMM_WORLD
  mpi::communicator c;

public:
  TRIQS_WRAP_ARG_AS_DICT // Wrap the solver parameters as a dictionary in python
                         // with the clang tool
                         /// constructor
                         /**
                           Allocates the main observables.
                           */
                         solver_core(constr_params_t const &p);

  TRIQS_WRAP_ARG_AS_DICT // Wrap the solver parameters as a dictionary in python
                         // with the clang tool
      /// solve method: starts the Metropolis algorithm
      /**
       * Steps:
       *
       * - extract $\Delta^\sigma_{ab}(\tau)$ and $\mu^\sigma_a$ from
       * $\mathcal{G}^\sigma_{ab}(i\omega)$
       *
       * - if $\mathcal{D}^{\sigma\sigma'}_{0,ab}(i\Omega)\neq 0$,  extract
       * $K(^{\sigma\sigma'}_{ab}\tau)$ and $\partial_\tau
       * K^{\sigma\sigma'}_{ab}(\tau)$ from
       * $\mathcal{D}^{\sigma\sigma'}_{0,ab}(i\Omega)$
       *
       * - if $\mathcal{J}_{\perp,a}(i\Omega)\neq 0$,  extract $\partial_\tau
       * K_{\perp,a}(\tau)$ from $\mathcal{J}_{\perp,a}(i\Omega)$
       * - add the moves and measures according to the parameters supplied by
       * the user
       *
       * - start the Monte-Carlo simulation
       *
       * - finalize the Monte Carlo simulation
       */
      void
      solve(solve_params_t const &p);

  // Struct containing the parameters relevant for the solver construction
  constr_params_t constr_params;

  // Struct containing the parameters relevant for the solve process
  std::optional<solve_params_t> last_solve_params;

  void sanity_check(solve_params_t const &p, int n_w, int n_w_b);

  // Green's function views for Python interface
  // do NOT add const here : python use a non const object and non const view

  /// Weiss field $\mathcal{G}^{\sigma}_{0,ab}(i\omega)$
  block_gf_view<imfreq> G0_iw() { return g0w; }
  /// Density-density retarded interactions
  /// $\mathcal{D}^{\sigma\sigma'}_{0,ab}(i\Omega)$
  block_gf_view<imfreq> D0_iw() { return d0w; }
  /// Dynamical spin-spin interaction, perpendicual components:
  /// $\mathcal{J}_\perp(i\Omega)$
  gf_view<imfreq> Jperp_iw() { return jperpw; }

  /// Hybridization function $\Delta^\sigma_{ab}(\tau)$
  block_gf_view<imtime> Delta_tau() { return delta_raw; }
  /// Dynamical kernel $K(\tau)$
  block_gf_view<imtime> K_tau() { return kt; }
  /// Derivative of the dynamical kernel $\partial_\tau K(\tau)$
  block_gf_view<imtime> Kprime_tau() { return kprimet; }
  /// Derivative of the dynamical kernel $\partial_\tau K_\perp(\tau)$
  gf_view<imtime> Kperpprime_tau() { return kperpprimet; }
  /// Dynamical spin-spin interactions $\mathcal{J}_\perp(\tau)$
  gf_view<imtime> Jperp_tau() { return jperpt; }

  /// Impurity Green's function $G^\sigma_{ab}(\tau)$ (see [[measure_gt]])
  block_gf_view<imtime> G_tau() { return gt; }
  /// Improved estimator function $F^\sigma_{ab}(\tau)$ (see [[measure_gt]])
  block_gf_view<imtime> F_tau() { return ft; }
  /// Density-density correlation function $\langle n^\sigma_{a}(\tau)
  /// n^{\sigma'}_{b}(0)\rangle$ (see [[measure_nnt]])
  block_gf_view<imtime> nn_tau() { return nnt; }

  /// Impurity Green's function $G^\sigma_{ab}(i\omega)$ (see [[measure_gw]])
  block_gf_view<imfreq> G_iw() { return gw; }
  /// Improved estimator function $F^\sigma_{ab}(i\omega)$ (see [[measure_gw]])
  block_gf_view<imfreq> F_iw() { return fw; }
  /// Impurity self-energy $\Sigma^\sigma_{ab}(i\omega)$ (see [[measure_gw]])
  block_gf_view<imfreq> Sigma_iw() { return sw; }
  /// Density-density correlation function $\mathrm{FT}\left[\langle
  /// n^\sigma_{a}(\tau) n^{\sigma'}_{b}(0)\rangle\right]$ (see [[measure_nnw]])
  block_gf_view<imfreq> nn_iw() { return nnw; }

  /// 3-point correlation function $\chi^{\sigma\sigma'}_{abc}(i\omega,i\Omega)$
  /// (see [[measure_g2w]])
  block_f_om_nu_tv_vt G_2w() { return g2w; }
  /// 3-point improved estimator (see [[measure_g2w]])
  block_f_om_nu_tv_vt F_2w() { return f2w; }

  /// 4-point correlation function
  /// $\chi^{\sigma\sigma'}_{abcd}(i\omega,i\omega',i\Omega)$ (see
  /// [[measure_g3w]])
  gf_3w_container_t G_3w() { return g3w; }
  /// 4-point improved estimator (see [[measure_g3w]])
  gf_3w_container_t F_3w() { return f3w; }

  /// Impurity Green's function in Legendre basis $G^\sigma_{ab}(n)$ (see
  /// [[measure_gl]])
  g_l_vt G_l() { return gl; }
  /// Improved estimator function in Legendre basis $G^\sigma_{ab}(n)$ (see
  /// [[measure_gl]])
  g_l_vt F_l() { return fl; }
  /// Spin spin correlation function $\langle s_+(\tau) s_-(0)\rangle$ (see
  /// [[measure_chipmt]])
  gf_view<imtime> chipm_tau() { return chipmt; }

  /// density-density static correlation $\langle n^\sigma_a n^{\sigma'}_b
  /// \rangle$ (see [[measure_nn]])
  block_matrix<double> const &nn() const { return nn_matrix; }
  /// histogram of hybridization perturbation order (see [[measure_hist]])
  matrix_view<double> histogram() const { return hist; }
  /// histogram of $\mathcal{J}_\perp$ perturbation order (see
  /// [[measure_hist_composite]])
  matrix_view<double> histogram_composite() const { return hist_composite; }
  /// histogram of the boundary states of the trace (see [[measure_statehist]])
  vector_view<double> state_histogram() const { return state_hist; }

  /// Monte Carlo sign
  double average_sign() { return average_sign_; }

  double percent_done() const { return percent_done_; }
};

} // namespace triqs_ctseg
