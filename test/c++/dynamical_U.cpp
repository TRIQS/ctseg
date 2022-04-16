#include <cmath>
#include <triqs/test_tools/gfs.hpp>
#include <triqs_ctseg/solver_core.hpp>

TEST(CTSEG, Dynamical_U) {
  // Start the mpi
  mpi::communicator world;

  double beta         = 20.0;
  double U            = 2.0;
  double mu           = 1.0;
  double epsilon      = 0.2;
  int n_cycles        = 10000;
  int n_warmup_cycles = 1000;
  int length_cycle    = 50;
  int random_seed     = 23488;
  int n_iw            = 5000;

  // Prepare the parameters
  constr_params_t param_constructor;
  solve_params_t param_solve;

  param_constructor.beta      = beta;
  param_constructor.gf_struct = {{"up", 1}, {"down", 1}};
  param_constructor.n_tau     = 20001;
  param_constructor.n_tau_k   = 20001;

  // Create solver instance
  solver_core ctqmc(param_constructor);

  param_solve.h_int           = U * n("up", 0) * n("down", 0);
  param_solve.hartree_shift   = {mu, mu};
  param_solve.n_cycles        = n_cycles;
  param_solve.n_warmup_cycles = n_warmup_cycles;
  param_solve.length_cycle    = length_cycle;
  param_solve.random_seed     = random_seed;
  // Measures
  param_solve.measure_gt  = true;
  param_solve.measure_nnt = false;
  // Moves
  param_solve.move_insert_segment  = true;
  param_solve.move_remove_segment  = true;
  param_solve.move_split_segment   = true;
  param_solve.move_regroup_segment = true;
  param_solve.move_move_segment    = true;

  // Prepare delta
  nda::clef::placeholder<0> om_;
  auto delta_w   = gf<imfreq>({beta, Fermion, n_iw}, {1, 1});
  auto delta_tau = gf<imtime>({beta, Fermion, param_constructor.n_tau}, {1, 1});
  delta_w(om_) << 1.0 / (om_ - epsilon);
  delta_tau()          = fourier(delta_w);
  ctqmc.Delta_tau()[0] = delta_tau;
  ctqmc.Delta_tau()[1] = delta_tau;

  // Prepare dynamical interaction
  double l  = 1.0; // electron boson coupling
  double w0 = 1.0; // screening frequency
  auto D0w  = gf<imfreq>({beta, Boson, n_iw}, {1, 1});
  auto D0t  = gf<imtime>({beta, Boson, param_constructor.n_tau}, {1, 1});
  D0w(om_) << 2 * l * l * w0 / (om_ * om_ - w0 * w0);
  D0t()                                = fourier(D0w);
  ctqmc.D0_tau().data()(range(), 0, 0) = D0t.data()(range(), 0, 0);
  ctqmc.D0_tau().data()(range(), 0, 1) = D0t.data()(range(), 0, 0);
  ctqmc.D0_tau().data()(range(), 1, 0) = D0t.data()(range(), 0, 0);
  ctqmc.D0_tau().data()(range(), 1, 1) = D0t.data()(range(), 0, 0);

  // Solve!!
  ctqmc.solve(param_solve);

  // Save the results
  if (world.rank() == 0) {
    h5::file G_file("dynamical_U.out.h5", 'w');
    h5_write(G_file, "(ctqmc.G_tau()[0])", ctqmc.results.G_tau()[0]);
  }
  if (world.rank() == 0) {
    h5::file G_file("dynamical_U.ref.h5", 'r');
    gf<imtime> g;
    h5_read(G_file, "(ctqmc.G_tau()[0])", g);
    EXPECT_GF_NEAR(g, ctqmc.results.G_tau()[0]);
  }
}
MAKE_MAIN;
