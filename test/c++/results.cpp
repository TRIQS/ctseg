// Copyright (c) 2026, The Simons Foundation
// This file is part of TRIQS/ctseg and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gtest/gtest.h>

#include <triqs_ctseg/results.hpp>

using namespace triqs_ctseg;

TEST(ResultsH5, MissingNewCorrelationFields) {
  results_t source;
  source.average_sign = 1.0;
  source.dyn_phi_n    = nda::matrix<double>{{1.0, 2.0}, {3.0, 4.0}};
  source.dyn_phi_phi  = nda::matrix<double>{{5.0, 6.0}, {7.0, 8.0}};

  h5::file file("results_missing_new_fields.h5", 'w');
  h5_write(file, "results", source);
  auto group = h5::group(file).open_group("results");
  group.unlink("dyn_phi_n", true);
  group.unlink("dyn_phi_phi", true);

  results_t restored;
  EXPECT_NO_THROW(h5_read(file, "results", restored));
  EXPECT_FALSE(restored.dyn_phi_n.has_value());
  EXPECT_FALSE(restored.dyn_phi_phi.has_value());
}
