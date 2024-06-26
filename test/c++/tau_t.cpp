// Copyright (c) 2022-2024 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Olivier Parcollet, Nils Wentzell

#include <iostream>

#include <cmath>
#include <triqs/test_tools/arrays.hpp>
#include <triqs_ctseg/tau_t.hpp>

using triqs_ctseg::tau_t;

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
