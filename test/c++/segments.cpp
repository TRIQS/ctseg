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

TEST(segment, find) {
  tau_t::set_beta(beta);

  auto tau1 = make_tau(3);
  auto v    = vs_t{S(tau1, 2.5), S(2, 1)};

  EXPECT_EQ(find_segment(v, tau1) - v.begin(), 0);

  EXPECT_EQ(find_segment(v, tau_t{2.0}) - v.begin(), 1);
  EXPECT_EQ(find_segment(v, tau_t{2.1}) - v.begin(), 1);
  EXPECT_EQ(find_segment(v, tau_t{1.9}) - v.begin(), 2);
}

// ------------------------------

// FIXME : more case, with cyclic ...
TEST(segment, is_insertable) {
  tau_t::set_beta(beta);

  auto v    = vs_t{S(3, 2.5), S(2, 1)};

  EXPECT_TRUE(is_insertable_into(S(2.4, 2.3), v));
  EXPECT_FALSE(is_insertable_into(S(2.1, 2.3), v));
  EXPECT_FALSE(is_insertable_into(S(2.5, 2.3), v));
  EXPECT_FALSE(is_insertable_into(S(2.1, 2), v));

  EXPECT_TRUE(is_insertable_into(S(0.5, 5), v));
  
}

// TEST OVERLAP 
// 
