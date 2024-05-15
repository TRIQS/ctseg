#include <iostream>

#include <cmath>
#include <triqs/test_tools/arrays.hpp>
#include <triqs_ctseg/tau_t.hpp>
#include <triqs_ctseg/configuration.hpp>

using vs_t = std::vector<segment_t>;

double beta      = 10;
double precision = 1.e-13;

tau_t make_tau(double x) { return tau_t{x}; } //uint64_t((x / tau_t::get_beta()) * double(tau_t::n_max))};}
tau_t make_tau(tau_t x) { return x; }

segment_t S(auto x, auto y) { return {make_tau(x), make_tau(y)}; }

// ------------------------------

TEST(segment, overlap) {
  tau_t::set_beta(beta);
  EXPECT_NEAR(overlap(S(3, 2), S(2.5, 1.5)), 0.5, precision);
  // noncyclic, cyclic
  EXPECT_NEAR(overlap(S(3, 2), S(2.5, 10)), 0.5, precision);
  // cyclic, noncyclic
  EXPECT_NEAR(overlap(S(2.5, 10), S(3, 2)), 0.5, precision);
}

// ------------------------------

TEST(segment, flip) {
  tau_t::set_beta(beta);
  auto v  = vs_t{S(4, 3), S(2, 1)};
  auto vf = vs_t{S(3, 2), S(1, 4)};
  //std::cout  << flip(v) << std::endl;
  //std::cout  << vf << std::endl;
  EXPECT_EQ(flip(v), vf);
  EXPECT_EQ(flip(flip(v)), v);
}

// ------------------------------

TEST(segment, lower_bound) {
  tau_t::set_beta(beta);

  auto tau1 = make_tau(3);
  auto v    = vs_t{S(tau1, 2.5), S(2, 1)};

  EXPECT_EQ(lower_bound(v, tau1) - v.begin(), 0);

  EXPECT_EQ(lower_bound(v, tau_t{2.0}) - v.begin(), 1);
  EXPECT_EQ(lower_bound(v, tau_t{2.1}) - v.begin(), 1);
  EXPECT_EQ(lower_bound(v, tau_t{1.9}) - v.begin(), 2);
}

// ------------------------------

// FIXME : more case, with cyclic ...
TEST(segment, is_insertable) {
  tau_t::set_beta(beta);

  auto v = vs_t{S(3, 2.5), S(2, 1)};

  EXPECT_TRUE(is_insertable_into(S(2.4, 2.3), v));
  EXPECT_FALSE(is_insertable_into(S(2.1, 2.3), v));
  EXPECT_FALSE(is_insertable_into(S(2.5, 2.3), v));
  EXPECT_FALSE(is_insertable_into(S(2.1, 2), v));

  EXPECT_TRUE(is_insertable_into(S(0.5, 5), v));
}

// TEST OVERLAP
//
