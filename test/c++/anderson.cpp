/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet
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
#include <cmath>
#include <triqs/test_tools/gfs.hpp>
#include <triqs_ctseg/solver_core.hpp>

TEST(CtHybSpin, Anderson) {
  // start the mpi
  mpi::communicator world;

  double beta    = 20.0;
  double U       = 1.0;
  double mu      = 1.3;
  double epsilon = 0.2;
  //int n_cycles   = 10000 / world.size();
  int n_cycles = 1;
  //int n_warmup_cycles = 1000;
  int n_warmup_cycles = 1;
  int length_cycle    = 50;
  int random_seed     = 23488 + 28 * world.rank();

  // prepare the parameters
  constr_params_t param_constructor;
  solve_params_t param_solve;

  param_constructor.beta      = beta;
  param_constructor.gf_struct = {{"up", 1}, {"down", 1}};
  int n_iw                    = 200;

  // create solver instance
  solver_core ctqmc(param_constructor);

  param_solve.h_int                = U * n("up", 0) * n("down", 0);
  param_solve.hartree_shift        = {mu, mu};
  param_solve.n_cycles             = n_cycles;
  param_solve.n_warmup_cycles      = n_warmup_cycles;
  param_solve.length_cycle         = length_cycle;
  param_solve.random_seed          = random_seed;
  param_solve.measure_gt           = true;
  param_solve.move_insert_segment  = true;
  param_solve.move_remove_segment  = true;
  param_solve.move_split_segment   = true;
  param_solve.move_regroup_segment = true;

  // prepare G0
  nda::clef::placeholder<0> om_;
  auto delta_w   = gf<imfreq>({beta, Fermion, n_iw}, {1, 1});
  auto delta_tau = gf<imtime>({beta, Fermion, param_constructor.n_tau}, {1, 1});
  delta_w(om_) << 1.0 / (om_ - epsilon);
  delta_tau()          = fourier(delta_w);
  ctqmc.Delta_tau()[0] = delta_tau;
  ctqmc.Delta_tau()[1] = delta_tau;

  // solve!!
  ctqmc.solve(param_solve);

  // Save the results
  if (world.rank() == 0) {
    h5::file G_file("anderson.out.h5", 'w');
    h5_write(G_file, "(ctqmc.G_tau()[0])", ctqmc.G_tau()[0]);
  }
  if (world.rank() == 0) {
    h5::file G_file("anderson.ref.h5", 'r');
    gf<imtime> g;
    h5_read(G_file, "(ctqmc.G_tau()[0])", g);
    EXPECT_GF_NEAR(g, ctqmc.G_tau()[0]);
  }
}
MAKE_MAIN;
