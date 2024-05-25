#include <cmath>
#include <triqs/test_tools/gfs.hpp>
#include <triqs_ctseg/solver_core.hpp>

using triqs::operators::n;
TEST(CTSEGJ, Dynamical_U) {
  
  mpi::communicator c; // Start the mpi

  double beta         = 20.0;
  double U            = 2.0;
  double mu           = 1.0;
  double epsilon      = 0.2;
  int n_cycles        = 10000;
  int n_warmup_cycles = 1000;
  int length_cycle    = 50;
  int random_seed     = 23488;
  int n_iw            = 5000;
  int n_tau           = 10001;
  double precision    = 1.e-13;

  // Prepare the parameters
  constr_params_t param_constructor;
  solve_params_t param_solve;

  // Construction parameters
  param_constructor.beta      = beta;
  param_constructor.gf_struct = {{"up", 1}, {"down", 1}};
  param_constructor.n_tau     = n_tau;
  param_constructor.n_tau_k   = n_tau;

  // Create solver instance
  solver_core Solver(param_constructor);

  // Solve parameters
  param_solve.h_int           = U * n("up", 0) * n("down", 0);
  param_solve.chemical_potential   = {mu, mu};
  param_solve.n_cycles        = n_cycles;
  param_solve.n_warmup_cycles = n_warmup_cycles;
  param_solve.length_cycle    = length_cycle;
  param_solve.random_seed     = random_seed;
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

  // Prepare dynamical interaction
  double l  = 1.0; // electron boson coupling
  double w0 = 1.0; // screening frequency
  auto D0w  = gf<imfreq>({beta, Boson, n_iw}, {1, 1});
  auto D0t  = gf<imtime>({beta, Boson, param_constructor.n_tau}, {1, 1});
  D0w(om_) << 2 * l * l * w0 / (om_ * om_ - w0 * w0);
  D0t()                                    = fourier(D0w);
  Solver.D0_tau().data()(range::all, 0, 0) = D0t.data()(range::all, 0, 0);
  Solver.D0_tau().data()(range::all, 0, 1) = D0t.data()(range::all, 0, 0);
  Solver.D0_tau().data()(range::all, 1, 0) = D0t.data()(range::all, 0, 0);
  Solver.D0_tau().data()(range::all, 1, 1) = D0t.data()(range::all, 0, 0);

  // Solve
  Solver.solve(param_solve);

  // Save the results
  if (c.rank() == 0) {
    h5::file out_file("dynamical_U.out.h5", 'w');
    h5_write(out_file, "G_tau", Solver.results.G_tau);
    h5_write(out_file, "F_tau", Solver.results.F_tau.value());
    h5_write(out_file, "nn_tau", Solver.results.nn_tau.value());
    h5_write(out_file, "nn_static", Solver.results.nn_static.value());
    h5_write(out_file, "densities", Solver.results.densities);
  }

  // Compare with reference
  if (c.rank() == 0) {
    h5::file ref_file("dynamical_U.ref.h5", 'r');
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
