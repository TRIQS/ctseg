#pragma once
#include <h5/h5.hpp>
#include <limits>
#include <iostream>
#include <cmath>
#include <assert.h>

#include <fmt/ostream.h>

/**
  * Discretized Imaginary Time type tau_t.
  * A point in imaginary time, i.e. $\tau \in [0,\beta]$, but defined on a very thin grid.
  *
  * * Regular type.
  *
  * * **Rationale**: the position in the segment is given by an uint64_t, i.e. a very long integer.
  *   This allows exact comparisons, which are notoriously dangerous on floating point number.
  *
  */

class tau_t {

  /// Inverse temperature associated with all $\tau$ points
  inline static double _beta;

  /// $\tau$ value, represented as an integer on a very fine grid
  uint64_t n = 0;

  public:
  /// Maximum value that can be stored inside a uint64_t
  static constexpr uint64_t n_max = std::numeric_limits<uint64_t>::max();

  /// Set the temperature
  static void set_beta(double beta) { _beta = beta; }

  ///
  tau_t() = default;

  /// Not for users. Use the factories
  tau_t(uint64_t n_) : n(n_) {}
  /// For test only, not for users. Use the factories
  tau_t(double x) {
    if ((x > _beta || x < 0)) {
      throw std::invalid_argument("Time tau must be in the range [0, beta]");
    } else if (x == _beta)
      n = n_max;
    else
      n = uint64_t((x / _beta) * double(n_max));
  }

  /// Comparisons (using integer, so it is safe)
  auto operator<=>(tau_t const &tau) const { return n <=> tau.n; }
  bool operator==(tau_t const &tau) const { return n == tau.n; }

  /// To cast to double, but it has to be done explicitly.
  explicit operator double() const { return _beta * (double(n) / double(n_max)); }

  /// tau_t at tau = beta
  static tau_t beta() { return {uint64_t{n_max}}; }

  /// $\tau = 0$
  static tau_t zero() { return {uint64_t{0}}; }

  /// Get epsilon, defined as $\epsilon = \beta /N_max$
  static tau_t epsilon() { return {uint64_t{1}}; }

  // Get a random point in $]tau1, tau2[$
  static tau_t random(auto &rng, tau_t const &tau1, tau_t const &tau2) {
    auto n1 = tau1.n + 1;
    return {rng(tau2.n - n1) + n1};
  }

  /// Get a random point in $]0, tau[$
  static tau_t random(auto &rng, tau_t const &tau) { return {rng(tau.n - 1) + 1}; }

  // ----- Some basic op, with cyclicity --------

  friend tau_t operator+(tau_t const &a, tau_t const &b) { return {a.n + b.n}; }

  friend tau_t operator-(tau_t const &a, tau_t const &b) { return {a.n - b.n}; }

  friend tau_t operator-(tau_t const &a) { return {tau_t::n_max - a.n}; }

  // ----- IO --------

  /// Stream insertion
  friend std::ostream &operator<<(std::ostream &out, tau_t const &p) {
    return out << double(p) << " [tau_t : n = " << p.n << "]";
  }
};

template <> struct fmt::formatter<tau_t> : ostream_formatter {};

// ----- arithmetic operations --------

/// other operations below decay to double
inline double operator*(tau_t const &a, tau_t const &b) { return double(a) * double(b); }
inline double operator/(tau_t const &a, tau_t const &b) { return double(a) / double(b); }

inline double operator+(tau_t const &x, double y) { return double(x) + y; }
inline double operator+(double y, tau_t const &x) { return y + double(x); }

inline double operator-(tau_t const &x, double y) { return double(x) - y; }
inline double operator-(double y, tau_t const &x) { return y - double(x); }

inline double operator*(tau_t const &x, double y) { return double(x) * y; }
inline double operator*(double y, tau_t const &x) { return y * double(x); }

inline double operator/(tau_t const &x, double y) { return double(x) / y; }
inline double operator/(double y, tau_t const &x) { return y / double(x); }
