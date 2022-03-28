// Copyright (c) 2013-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2013-2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2016 Igor Krivenko
// Copyright (c) 2018-2020 Simons Foundation
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
// Authors: Michel Ferrero, Igor Krivenko, Olivier Parcollet, Priyanka Seth, Hugo U. R. Strand, Nils Wentzell

#pragma once
#include <h5/h5.hpp>
#include <limits>
#include <iostream>
#include <cmath>

namespace triqs::utility {

  /**
  * Discretized Imaginary Time dimtime
  * A point in imaginary time, i.e. $\tau \in [0,\beta]$, but defined on a very thin grid.
  *
  * * Regular type.
  *
  * * **Rationale**: the position in the segment is given by an uint64_t,
  *   i.e. a very long integer.
  *   This allows exact comparisons, which notoriously dangerous on floating point number.
  *
  */
  struct dimtime_t {

    double const beta              = 0;
    double const val               = 0;
    uint64_t const n               = 0;
    static constexpr uint64_t Nmax = std::numeric_limits<uint64_t>::max();
   
    ///
    dimtime_t() = default;
    
    // Not for users. Use the factories
    dimtime_t(uint64_t n_, double beta_) : beta(beta_), val(beta_ * (double(n_) / double(Nmax))), n(n_) {}

    /// Comparisons (using integer, so it is safe)
    auto operator<=>(dimtime_t const &tau) const { return n <=> tau.n; }
    bool operator==(dimtime_t const &tau) const { return n == tau.n; }

    /// To cast to double, but it has to be done explicitly.
    explicit operator double() const { return val; } // cast to a double

    // ---- Factories ----

    /// Get a random point in $[0, tau[$
    static dimtime_t random(auto &rng, dimtime_t const &tau) { return {rng(tau.n), tau.beta}; }

    /// Get a random point in $[0,\beta[$
    static dimtime_t random(auto &rng, double beta) { return {rng(Nmax), beta}; }

    // Get a random point in $[tau1, tau2[$
    static dimtime_t random(auto &rng, dimtime_t const &tau1, dimtime_t const &tau2) {
      return {rng(tau2.n - tau1.n) + tau1.n, tau1.beta};
    }

    /// Get minimum point (i.e. $\tau =0$)
    static dimtime_t zero(double beta) { return {0, beta}; }

    /// Get epsilon, defined as $\epsilon = \beta /N_max$
    static dimtime_t epsilon(double beta) { return {1, beta}; }
  };

  // ----- arithmetic operations --------

  // NB adding and substracting is cyclic on [0, beta]

  inline dimtime_t operator+(dimtime_t const &a, dimtime_t const &b) {
    bool wrapped = ((dimtime_t::Nmax - std::max(a.n, b.n)) < std::min(a.n, b.n));
    if (!wrapped)
      return {(a.n + b.n) % dimtime_t::Nmax, a.beta};
    else
      return {((a.n + b.n) + 1) % dimtime_t::Nmax, a.beta};
  }

  inline dimtime_t operator-(dimtime_t const &a, dimtime_t const &b) {
    uint64_t n = (a.n >= b.n ? (a.n - b.n) % dimtime_t::Nmax : dimtime_t::Nmax - (b.n - a.n));
    return {n, a.beta};
  }

  /// unary -
  inline dimtime_t operator-(dimtime_t const &a) { return {dimtime_t::Nmax - a.n, a.beta}; }

  /// division by integer
  //inline dimtime_t div_by_int(dimtime_t const &a, size_t b) { return {a.n / b, a.beta}; }

  /// Multiplication by int
  //inline dimtime_t mult_by_int(dimtime_t const &a, size_t b) { return {a.n * b, a.beta}; }

  /// floor_div(x,y) = floor (x/y), but computed on the grid.
  //inline size_t floor_div(dimtime_t const &a, dimtime_t const &b) { return a.n / b.n; }

  /// other operations below decay to double
  inline double operator*(dimtime_t const &a, dimtime_t const &b) { return double(a) * double(b); }
  inline double operator/(dimtime_t const &a, dimtime_t const &b) { return double(a) / double(b); }

  inline double operator +(dimtime_t const &x, double y) { return double(x) + y; }
  inline double operator +(double y, dimtime_t const &x) { return y + double(x); }
  
  inline double operator -(dimtime_t const &x, double y) { return double(x) - y; }
  inline double operator -(double y, dimtime_t const &x) { return y - double(x); }
  
  inline double operator *(dimtime_t const &x, double y) { return double(x) * y; }
  inline double operator *(double y, dimtime_t const &x) { return y * double(x); }
  
  inline double operator /(dimtime_t const &x, double y) { return double(x) / y; }
  inline double operator /(double y, dimtime_t const &x) { return y / double(x); }

  // ----- IO --------

  /// Stream insertion
  inline std::ostream &operator<<(std::ostream &out, dimtime_t const &p) {
    return out << p.val << " [dimtime_t : beta = " << p.beta << " n = " << p.n << "]";
  }

  /// Write into HDF5
  inline void h5_write(h5::group fg, std::string const &subgroup_name, dimtime_t const &g) {
    auto gr = fg.create_group(subgroup_name);
    h5_write(gr, "beta", g.beta);
    h5_write(gr, "val", g.val);
    h5_write(gr, "n", g.n);
  }

  /// Read from HDF5
  inline void h5_read(h5::group fg, std::string const &subgroup_name, dimtime_t &g) {
    auto gr = fg.open_group(subgroup_name);
    h5_read(gr, "beta", g.beta);
    h5_read(gr, "val", g.val);
    h5_read(gr, "n", g.n);
  }
} // namespace triqs::utility
