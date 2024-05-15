#include <cmath>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/test_tools/arrays.hpp>
#include <triqs_ctseg-J/solver_core.hpp>

using triqs::operators::n;

TEST(CTSEGJ, Anderson) {
  // Start the mpi
  mpi::communicator c;

  double beta         = 20.0;
  double U            = 1.0;
  double mu           = 0.5;
  double epsilon      = 0.2;
  int n_cycles        = 10000;
  int n_warmup_cycles = 1000;
  int length_cycle    = 50;
  int random_seed     = 23488;
  int n_iw            = 1000;
  int n_tau           = 10001;
  double precision    = 1.e-13;

  // Prepare the parameters
  constr_params_t param_constructor;
  solve_params_t param_solve;

  param_constructor.beta      = beta;
  param_constructor.gf_struct = {{"up", 1}, {"down", 1}};
  param_constructor.n_tau     = n_tau;

  // Create solver instance
  solver_core Solver(param_constructor);

  param_solve.h_int           = U * n("up", 0) * n("down", 0);
  param_solve.hartree_shift   = {mu, mu};
  param_solve.n_cycles        = n_cycles;
  param_solve.n_warmup_cycles = n_warmup_cycles;
  param_solve.length_cycle    = length_cycle;
  param_solve.random_seed     = random_seed;
  // Measures
  param_solve.measure_ft  = true;
  param_solve.measure_nnt = true;
  param_solve.measure_nn  = true;

  // Prepare delta
  nda::clef::placeholder<0> om_;
  auto delta_w   = gf<imfreq>({beta, Fermion, n_iw}, {1, 1});
  auto delta_tau = gf<imtime>({beta, Fermion, param_constructor.n_tau}, {1, 1});
  delta_w(om_) << 1.0 / (om_ - epsilon);
  delta_tau()           = fourier(delta_w);
  Solver.Delta_tau()[0] = delta_tau;
  Solver.Delta_tau()[1] = delta_tau;

  // Solve
  Solver.solve(param_solve);

  // Save the results
  if (c.rank() == 0) {
    h5::file out_file("anderson.out.h5", 'w');
    h5_write(out_file, "G_tau", Solver.results.G_tau);
    h5_write(out_file, "F_tau", Solver.results.F_tau.value());
    h5_write(out_file, "nn_tau", Solver.results.nn_tau.value());
    h5_write(out_file, "nn_static", Solver.results.nn_static.value());
    h5_write(out_file, "densities", Solver.results.densities);
  }

  // Compare with reference
  if (c.rank() == 0) {
    h5::file ref_file("anderson.ref.h5", 'r');
    block_gf<imtime> G_tau, F_tau;
    gf<imtime> nn_tau;
    nda::matrix<double> nn_static;
    nda::array<double, 1> densities;
    h5_read(ref_file, "G_tau", G_tau);
    h5_read(ref_file, "F_tau", F_tau);
    h5_read(ref_file, "nn_tau", nn_tau);
    h5_read(ref_file, "nn_static", nn_static);
    h5_read(ref_file, "densities", densities);
    EXPECT_ARRAY_NEAR(densities, Solver.results.densities, precision);
    EXPECT_BLOCK_GF_NEAR(G_tau, Solver.results.G_tau, precision);
    EXPECT_BLOCK_GF_NEAR(F_tau, Solver.results.F_tau.value(), precision);
    EXPECT_GF_NEAR(nn_tau, Solver.results.nn_tau.value(), precision);
    EXPECT_ARRAY_NEAR(nn_static, Solver.results.nn_static.value(), precision);
  }
}
MAKE_MAIN;
