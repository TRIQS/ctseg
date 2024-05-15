#include <iostream>

#include <cmath>
#include <triqs/test_tools/arrays.hpp>
#include <triqs_ctseg/tau_t.hpp>

TEST(tau, basic) {

  double beta = 40;
  tau_t::set_beta(beta);
  double precision = 1.e-13;

  EXPECT_EQ(tau_t::n_max + 1, 0);

  auto dzero = tau_t::zero();
  auto dbeta = tau_t::beta();
  auto deps  = tau_t::epsilon();

  EXPECT_EQ(dbeta - dzero, dbeta);
  EXPECT_EQ(dbeta + dzero, dbeta);

  EXPECT_EQ(dbeta + deps, dzero);
  EXPECT_EQ(dbeta - deps, (tau_t{tau_t::n_max - 1}));
  EXPECT_EQ((dbeta - deps) + deps, dbeta);

  auto t1 = tau_t{tau_t::n_max / 5};
  auto t2 = tau_t{(tau_t::n_max / 5) * 4};

  // almost cyclic
  EXPECT_EQ(t1 + dbeta, t1 - deps);
  EXPECT_EQ(t1 - dbeta, t1 + deps);

  // the eps is 0 for beta < 1000 or so, at machine precision double
  EXPECT_EQ(double(t1 + deps), double(t1));

  EXPECT_EQ(t1 + t2, dbeta);

  EXPECT_NEAR(t1 + t2 - (double(t1) + double(t2)), 0, precision);

  // no cyclicity : t2 > t1
  EXPECT_NEAR(t2 - t1 - (double(t2) - double(t1)), 0, precision);

  // using cyclicity !
  EXPECT_NEAR(t1 - t2 - (double(t1) - double(t2)), beta, precision);

  // just as double
  EXPECT_NEAR(t1 * t2 - (double(t1) * double(t2)), 0, precision);
  EXPECT_NEAR(t1 / t2 - (double(t1) / double(t2)), 0, precision);

  EXPECT_NEAR(double(-t1), 32, precision);
}
