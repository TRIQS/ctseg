/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
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
#include "solver_core.hpp"

#include "configuration.hpp"
#include "evaluate_vertex.hpp"
#include "precompute_Mw.hpp"
#include "precompute_fprefactor.hpp"
#include "precompute_nw.hpp"
#include "measures.hpp"
#include "moves.hpp"

#include <itertools/itertools.hpp>
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>

namespace triqs_ctseg {

solver_core::solver_core(constr_params_t const &p) {

  beta = p.beta;
  gf_struct = p.gf_struct;

  delta = block_gf<imtime>(triqs::mesh::imtime{beta, Fermion, p.n_tau}, p.gf_struct);






  n_color = 0;
  std::vector<std::string> block_names;
  for (auto const &[bl_name, bl_size] : gf_struct) {
    block_names.push_back(bl_name);
    n_color += bl_size;
  }



  // FIXME BlockGf

  // prepare frequency-dependent Green's functions
  auto gf_imfreq_vec = [&](int n_mat) {
    std::vector<gf<imfreq>> green_v;
    for (auto const &[bl_name, bl_size] : gf_struct)
      green_v.emplace_back(
          gf<imfreq>{{beta, Fermion, n_mat}, make_shape(bl_size, bl_size)});
    return green_v;
  };
  auto gf_imtime_vec = [&](int n_tau) {
    std::vector<gf<imtime>> green_v;
    for (auto const &[bl_name, bl_size] : gf_struct)
      green_v.emplace_back(
          gf<imtime>{{beta, Fermion, n_tau}, make_shape(bl_size, bl_size)});
    return green_v;
  };
  auto gf_legendre_vec = [&](size_t n_l) {
    std::vector<gf<triqs::gfs::legendre>> green_v;
    for (auto const &[bl_name, bl_size] : gf_struct)
      green_v.emplace_back(gf<triqs::gfs::legendre>{
          {beta, Fermion, n_l}, make_shape(bl_size, bl_size)});
    return green_v;
  };

  // the order of combined blocks is, e.g. 11 12 21 22
  auto bosonic_block_names = std::vector<std::string>{};
  for (auto const &str1 : block_names)
    for (auto const &str2 : block_names)
      bosonic_block_names.push_back(str1 + "|" + str2);

  // if the blocks are not of the same size, pay attention:
  // the "off-diagonal" blocks will not be square
  //------------ freq
  auto D_imfreq_vec = [&](int n_iw) {
    std::vector<gf<imfreq>> green_v;
    for (auto const &[bl1_name, bl1_size] : gf_struct)
      for (auto const &[bl2_name, bl2_size] : gf_struct)
        green_v.emplace_back(
            gf<imfreq>{{beta, Boson, n_iw}, make_shape(bl1_size, bl2_size)});
    return green_v;
  };
  //------------ time
  auto D_imtime_vec = [&](int n_tau) {
    std::vector<gf<imtime>> green_v;
    for (auto const &[bl1_name, bl1_size] : gf_struct)
      for (auto const &[bl2_name, bl2_size] : gf_struct)
        green_v.emplace_back(
            gf<imtime>{{beta, Boson, n_tau}, make_shape(bl1_size, bl2_size)});
    return green_v;
  };

  // input containers
  delta_raw = make_block_gf<imtime>(block_names, gf_imtime_vec(p.n_tau));
  g0w = make_block_gf<imfreq>(block_names, gf_imfreq_vec(p.n_iw));

  d0w = make_block_gf<imfreq>(bosonic_block_names, D_imfreq_vec(p.n_w_b_nn));
  kt = make_block_gf<imtime>(bosonic_block_names, D_imtime_vec(p.n_tau_k));
  kprimet = make_block_gf<imtime>(bosonic_block_names, D_imtime_vec(p.n_tau_k));

  // FIXME : scalar function 
  jperpt = gf<imtime>({beta, Boson, p.n_tau_jperp}, {1, 1});
  jperpw = gf<imfreq>({beta, Boson, p.n_w_b_nn}, {1, 1});
  chipmt = gf<imtime>({beta, Boson, p.n_tau_jperp}, {1, 1});

  // imaginary-time Green's functions
  gt = make_block_gf<imtime>(block_names, gf_imtime_vec(p.n_tau));
  ft = make_block_gf<imtime>(block_names, gf_imtime_vec(p.n_tau));
  nnt = make_block_gf<imtime>(bosonic_block_names, D_imtime_vec(p.n_tau_nn));

  // Matsubara Green's functions
  gw = make_block_gf<imfreq>(block_names, gf_imfreq_vec(p.n_iw));
  fw = make_block_gf<imfreq>(block_names, gf_imfreq_vec(p.n_iw));
  sw = make_block_gf<imfreq>(block_names, gf_imfreq_vec(p.n_iw));
  nnw = make_block_gf<imfreq>(bosonic_block_names, D_imfreq_vec(p.n_w_b_nn));


  // Prepare output containers
// FIXME : measure


  std::vector<matrix<double>> nn_v;
  for (auto const &[bl1_name, bl1_size] : gf_struct)
    for (auto const &[bl2_name, bl2_size] : gf_struct)
      nn_v.push_back(matrix<double>(make_shape(bl1_size, bl2_size)));
  nn_matrix = block_matrix<double>{bosonic_block_names, nn_v};


  // HIsto gram 
  int Nmax = 500;
  hist.resize(n_color, Nmax);
  hist() = 0;
  hist_composite.resize(n_color, Nmax);
  hist_composite() = 0;
  state_hist.resize(ipow(2, n_color));
  state_hist() = 0.0;
}

// ---------------------------------------------------------------------------

void solver_core::solve(solve_params_t const &solve_params) {

  last_solve_params = solve_params;

  // http://patorjk.com/software/taag/#p=display&f=Calvin%20S&t=TRIQS%20ctint
  if (c.rank() == 0)
    std::cout << "\n"
                 "╔╦╗╦═╗╦╔═╗ ╔═╗  ┌─┐┌┬┐┌─┐┌─┐┌─┐\n"
                 " ║ ╠╦╝║║═╬╗╚═╗  │   │ └─┐├┤ │ ┬\n"
                 " ╩ ╩╚═╩╚═╝╚╚═╝  └─┘ ┴ └─┘└─┘└─┘\n";

  // Merge constr_params and solve_params
  params_t p(constr_params, solve_params);

  vector<double> mu(n_color);
  mu() = 0.;

  if (p.hartree_shift.size() > 0) {
    if (p.hartree_shift.size() != mu.size())
      TRIQS_RUNTIME_ERROR
          << "Mismatch in the size of the Hartree shift and mu array";
    for (int i : range(mu.size()))
      mu[i] = p.hartree_shift[i];
  }

  bool dynamical_U = !is_zero(d0w) or !is_zero(jperpw);
  bool full_spin_rot_inv = false;
  bool jperp_interactions = !is_zero(jperpw);

  if (dynamical_U) {
    // extract kt and kprimet from d0w
    for (int bl = 0; bl < kt.size(); bl++) {
      for (auto const &t : kt[bl].mesh()) {
        for (auto i = 0; i < kt[bl].target_shape()[0]; i++) {
          for (auto j = 0; j < kt[bl].target_shape()[1]; j++) {
            kt[bl][t](i, j) = 0.0;
            kprimet[bl][t](i, j) = 0.0;

            for (auto const &w : d0w[bl].mesh()) {
              if (w.n > 0) {
                double mats = dcomplex(w).imag();
                kt[bl][t](i, j) += real(d0w[bl](w)(i, j)) / (mats * mats) *
                                   (1 - cos(mats * t));

                kprimet[bl][t](i, j) +=
                    real(d0w[bl](w)(i, j)) / mats * sin(mats * t);
              }
            }
            kt[bl][t](i, j) *= 2. / beta;
            kprimet[bl][t](i, j) *= 2. / beta;
            kt[bl][t](i, j) +=
                0.5 * real(d0w[bl](0)(i, j)) * t * (t / beta - 1.);
            kprimet[bl][t](i, j) += real(d0w[bl](0)(i, j)) * (t / beta - .5);
          }
        }
      }
    }
    // extract kt and kprimet from d0w
    for (auto &t : kperpprimet.mesh()) {
      kperpprimet[t] = 0.0;
      for (auto &w : jperpw.mesh()) {
        double mats = dcomplex(w).imag();
        if (w.n > 0)
          kperpprimet[t] += real(jperpw(w)) / mats * sin(mats * t);
      }
      kperpprimet[t] *= 2. / beta;
      kperpprimet[t] += real(jperpw(0)) * (t / beta - .5);
      kperpprimet[t] *= 0.25; // D^z = J_z/4.; K pertains to D, not J
    }

    // FIXME CHECK 
    gf<imfreq> D_sp_minus_one_fourth_Jperp(jperpw.mesh(), {1, 1});
    auto ind_UV = [&](std::string const &s) {
      return get_index(d0w.block_names(), s);
    };
    for (auto const &iom : D_sp_minus_one_fourth_Jperp.mesh())
      D_sp_minus_one_fourth_Jperp[iom] =
          0.5 *
              (d0w[ind_UV(block_names[0] + "|" + block_names[0])](iom)(0, 0) -
               d0w[ind_UV(block_names[0] + "|" + block_names[1])](iom)(0, 0)) -
          0.25 * jperpw(iom)(0, 0);
    full_spin_rot_inv =
        is_zero(D_sp_minus_one_fourth_Jperp) and jperp_interactions;
  } // end if dynamical_U

  // Extract the U from the operator
  auto U_full = triqs::operators::utils::dict_to_matrix(
      triqs::operators::utils::extract_U_dict2(p.h_int), gf_struct);
  matrix<double> U(U_full.shape());
  U() = real(U_full);

  // KEEP 
  // Get the new values for mu and Umatrix from Kprime
  if (dynamical_U) {
    for (size_t i = 0; i < U.shape()[0]; ++i) {
      auto bl_ind_i = color_to_block_and_inner_index_impl(i, gf_struct);
      auto bl_ii = bl_ind_i.first * gf_struct.size() +
                   bl_ind_i.first; // double index (e.g. up|up)
      mu(i) +=
          kprimet[bl_ii][closest_mesh_pt(0.0)](bl_ind_i.second, bl_ind_i.second)
              .real();
      for (size_t j = 0; j < U.shape()[1]; ++j) {
        auto bl_ind_j = color_to_block_and_inner_index_impl(j, gf_struct);
        auto bl_ij = bl_ind_i.first * gf_struct.size() + bl_ind_j.first;
        if (i != j)
          U(i, j) -= 2. * kprimet[bl_ij][closest_mesh_pt(0.0)](bl_ind_i.second,
                                                               bl_ind_j.second)
                              .real();
      }
    }
  }

  // FIXME Gather the checks
  if (jperp_interactions and g0w.size() != 2)
    TRIQS_RUNTIME_ERROR << "Jperp interactions implemented only for case with "
                           "one orbital with two spins";

  // ................   QMC  ...................

  auto CTQMC = triqs::mc_tools::mc_generic<double>(p.random_name, p.random_seed,
                                                   p.verbosity);
  // FIXME spdlog ...
  if (c.rank() == 0) {
    std::cout << "mu = " << mu << std::endl;
    std::cout << "U = " << U << std::endl;
    if (dynamical_U)
      std::cout << "dynamical_U" << endl;
    if (jperp_interactions)
      std::cout << "jperp_interactions";
    if (full_spin_rot_inv)
      cout << " with full spin rotational invariance";
    cout << endl;
  }

  // FIXME : params --> data qmc_data ...
  //
  qmc_parameters params(n_color, beta, U, mu, kt, kprimet, kperpprimet,
                        gf_struct, dynamical_U, jperp_interactions,
                        full_spin_rot_inv, p);
  
  // FIXME : constructire new config
  configuration config(&params, gf_struct, delta_raw, jperpt,
                       p.keep_Jperp_negative);
  // initialize moves
  if (p.move_insert_segment)
    CTQMC.add_move(move_insert_segment(&params, &config, CTQMC.get_rng()),
                   "segment insertion");
  if (p.move_remove_segment)
    CTQMC.add_move(move_remove_segment(&params, &config, CTQMC.get_rng()),
                   "segment removal");
  if (p.move_move)
    CTQMC.add_move(move_move_segment(&params, &config, CTQMC.get_rng()),
                   "segment move to another line");
  if (p.move_swap_empty_lines)
    CTQMC.add_move(move_swap_empty_lines(&params, &config, CTQMC.get_rng()),
                   "swap empty lines");
  if (jperp_interactions)
    CTQMC.add_move(move_insert_spin_segment(&params, &config, CTQMC.get_rng()),
                   "insert composite segment");
  if (jperp_interactions)
    CTQMC.add_move(move_remove_spin_segment(&params, &config, CTQMC.get_rng()),
                   "remove composite segment");
  if (jperp_interactions)
    CTQMC.add_move(move_swap_bosonic_lines(&params, &config, CTQMC.get_rng()),
                   "swap bosonic lines");
  if (p.move_group_into_spin_segment)
    CTQMC.add_move(
        move_group_into_spin_segment(&params, &config, CTQMC.get_rng()),
        "group into composite segment");
  if (p.move_split_spin_segment)
    CTQMC.add_move(move_split_spin_segment(&params, &config, CTQMC.get_rng()),
                   "split composite segment");
  if (p.move_group_into_spin_segment2)
    CTQMC.add_move(
        move_group_into_spin_segment2(&params, &config, CTQMC.get_rng()),
        "group into composite segment2");
  if (p.move_split_spin_segment2)
    CTQMC.add_move(move_split_spin_segment2(&params, &config, CTQMC.get_rng()),
                   "split composite segment2");

   // FIXME ? REMOVE ? 

  // initialize precomputations shared between measurements
  auto aux_measure_f =
      std::make_shared<precompute_fprefactor>(&params, &config);
  if (p.measure_fw || p.measure_ft || p.measure_fl || p.measure_f2w ||
      p.measure_f3w || p.measure_g2w || p.measure_g3w)
    CTQMC.add_measure_aux(
        aux_measure_f); // MUST BE ADDED BEFORE MEASURE_M SINCE IT IS USED BY IT

  auto aux_measure_M = std::make_shared<precompute_Mw>(
      &params, &p, &config, aux_measure_f, p.n_w_f_vertex, p.n_w_b_vertex);
  if (p.measure_g2w || p.measure_g3w)
    CTQMC.add_measure_aux(aux_measure_M);
  auto aux_measure_n = std::make_shared<precompute_nw>(&params, &config);
  if (p.measure_g2w || p.measure_nnw)
    CTQMC.add_measure_aux(aux_measure_n);

  // initialize measurements
  if (p.measure_gt)
    CTQMC.add_measure(measure_gt(&params, &config, aux_measure_f, gt, ft),
                      "G(tau) measurement");
  if (p.measure_nnt)
    CTQMC.add_measure(measure_nnt(&params, &config, nnt),
                      "nn(tau) measurement");
  if (p.measure_nn)
    CTQMC.add_measure(measure_nn(&params, &config, &nn_matrix),
                      "density measurement");
  if (p.measure_sign)
    CTQMC.add_measure(measure_sign(&average_sign_),
                      "measurement of the MC sign");
 
   // FIXME  ? Keep ? ?
 
  if (p.measure_hist)
    CTQMC.add_measure(measure_hist(&config, &hist), "histogram measurement");
 
  if (p.measure_hist_composite)
    CTQMC.add_measure(measure_hist_composite(&config, &hist_composite),
                      "histogram of composite order measurement");
  if (p.measure_statehist)
    CTQMC.add_measure(measure_statehist(&params, &config, &state_hist),
                      "impurity state histogram measurement");

  // Run and collect results
  auto _solve_status =
      CTQMC.warmup_and_accumulate(p.n_warmup_cycles, p.n_cycles, p.length_cycle,
                                  triqs::utility::clock_callback(p.max_time));
  CTQMC.collect_results(c);

} // solve

} // namespace triqs_ctseg
