#pragma once
#include <vector>
//#include "params.hpp"
#include "types.hpp"

// FIXME : EXPLAIN ....

// Segment: (time of c, time of c^dagger)
struct segment_t {

  tau_t tau_c, tau_cdag; // time of c and cdag
  bool J_c = false, J_cdag = false;

  // Length of segment (accounts for cyclicity)
  tau_t length() const { return tau_c - tau_cdag; };

  static segment_t full_line() { return {tau_t::beta(), tau_t::zero()}; }
};

// Store a couple of S+, S-
// FIXME : need to store the color, which is today just 0,1
struct jperp_line_t {
  tau_t tau_Sminus, tau_Splus; // times of the S-, S+
};

// ----------

struct configuration_t {
  // A list of segments for each color, ORDERED on tau_c of segments, with decreasing order.
  std::vector<std::vector<segment_t>> seglists;
  // List of Jperp lines, not ordered.
  std::vector<jperp_line_t> Jperp_list;

  configuration_t(int n_color) : seglists(n_color) {}

  [[nodiscard]] int n_color() const { return seglists.size(); }
};

// ------------------- Invariants ---------------------------

// FIXME : separate files
void check_invariant(configuration_t const &config, std::vector<det_t> const &dets);

// ------------------- Functions to manipulate config --------------------------

// Comparison of segments. Returns 1 if s1 is left of s2 (we order segments by decreasing time).
inline bool operator<(segment_t const &s1, segment_t const &s2) { return s1.tau_c > s2.tau_c; };

// Equality of segments.
inline bool operator==(segment_t const &s1, segment_t const &s2) {
  return s1.tau_c == s2.tau_c and s1.tau_cdag == s2.tau_cdag;
};

// Whether a segment is wrapped around beta/0
inline bool is_cyclic(segment_t const &seg) { return seg.tau_cdag > seg.tau_c; }

// Check whether segment is a full line.
inline bool is_full_line(segment_t const &seg) { return seg.tau_c - seg.tau_cdag == tau_t::beta(); }

// Make a list of time ordered (decreasing) operators
// vector of (time, color, is_dagger)
std::vector<std::tuple<tau_t, int, bool>> make_time_ordered_op_list(configuration_t const &config);

// Find index of first segment starting left of seg.tau_c.
std::vector<segment_t>::const_iterator find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg);

// Overlap between two non-cyclic segments.
double overlap_seg(segment_t const &seg1, segment_t const &seg2);

// Checks if two segments are completely disjoint (accounting for boundaries)
bool disjoint(segment_t seg1, segment_t seg2);

// Overlap between segment and a list of segments.
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg);

// Checks if segment is movable to a given color
bool is_insertable(std::vector<segment_t> const &seglist, segment_t const &seg);

// Contribution of the dynamical interaction kernel K to the overlap between a segment and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau_c, tau_t const &tau_cdag,
                 gf<imtime, matrix_valued> const &K, int c1, int c2);

// Contribution of the dynamical interaction kernel K to the overlap between an operator and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau, bool is_c, gf<imtime, matrix_valued> const &K,
                 int c1, int c2);

// Length occupied by all segments for a given color
double density(std::vector<segment_t> const &seglist);

// Find the state at tau = 0 or beta
std::vector<bool> boundary_state(configuration_t const &config);

// Find segments corresponding to bosonic line
std::pair<std::vector<segment_t>::const_iterator, std::vector<segment_t>::const_iterator>
find_spin_segments(int line_idx, configuration_t const &config);

// Flip config
std::vector<segment_t> flip(std::vector<segment_t> const &sl);

// Sign of a config
double config_sign(configuration_t const &config, std::vector<det_t> const &dets);

// Same as std::lower_bound, but i-th element of vector is returned by f[i]
// f is called on 0:N strictly
long lower_bound(auto f, long N, auto const &value) {
  long first = 0, count = N;

  while (count > 0) {
    long step = count / 2;
    assert(first + step < N);
    if (f(first + step) < value) {
      first += step + 1;
      count -= step + 1;
    } else
      count = step;
  }
  return first;
}

// Find the indices of the segments whose cdag are in ]wtau_left,wtau_right[
std::vector<long> cdag_in_window(tau_t const &wtau_left, tau_t const &wtau_right,
                                 std::vector<segment_t> const &seglist);

// ----------- DEBUG code --------------
// Print config
std::ostream &operator<<(std::ostream &out, configuration_t const &config);
