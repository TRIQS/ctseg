#pragma once
#include <vector>
#include "tau_t.hpp"
#include "logs.hpp"
#include "dets.hpp"
#include "work_data.hpp"

// The MC configuration and associated functions.
// NB : all the time ordering are in DECREASING order,
// in agreement with the usual convention

// ================= Data structures =================

// The configuration is made of n_color ordered lists of segments
// and (if needed) a list of J_perp lines. 

// --------------- Segment -------------------
//
// A segment represents a couple (c, cdag), at given tau_t times in [0, beta]
// Segment is cyclic if tau_c < tau_cdag, non-cyclic otherwise. 
//
// Some segments are aligned between 2 colors, corresponding to a S^+, S^- operators
// The S operators will not be linked to Delta, the hybridization, but to Jperp lines.
// The segment stores two booleans (J_c, J_cdag) that indicate whether its operators are 
// linked to J_perp lines. 
//
// Special case : full lines.
// When a line has no operator, we need to take into account 2 states : empty or full.
// The full line case is encoded as a [beta, 0] segment, as it naturally yields the correct overlap.
// This special case will however need to be treated separately in moves
//
struct segment_t {

  tau_t tau_c, tau_cdag;            // time of c and cdag
  bool J_c = false, J_cdag = false; // Whether c (resp cdag) is part of a S operator

  /// Length of segment (accounts for cyclicity : length from c to cdag possibly across 0/beta)
  tau_t length() const { return tau_c - tau_cdag; };

  /// A segment [beta, 0] represent a full line
  static segment_t full_line() { return {tau_t::beta(), tau_t::zero()}; }
};

// operator descriptor: time, color, is_cdag[creation/annihilation]
struct colored_ops_t {
  tau_t tau;
  int color;
  bool is_cdag;
};

// simple alias
using vec_seg_iter_t = std::vector<segment_t>::const_iterator;

// ----------------- J_perp line -------------------
// Stores the times of a couple (S+, S-)
// Note: J_perp expansion only implemented for single orbital 
// (two colors, spin up and spin down). 
struct jperp_line_t {
  tau_t tau_Sminus, tau_Splus; // times of the S-, S+
};

// --------------- Configuration ----------------------
// The configuration is a list of of segments for each color, 
// and a list of J_perp lines.
struct configuration_t {
  // A list of segments for each color.
  // NB: ordered in DECREASING time order of the tau_c.
  std::vector<std::vector<segment_t>> seglists;

  // List of Jperp lines, NOT ordered.
  std::vector<jperp_line_t> Jperp_list;

  // Construct from the number of colors
  configuration_t(int n_color) : seglists(n_color) {}

  // Accessor number of colors
  int n_color() const { return seglists.size(); }
};

// ===================  Functions to manipulate segments ===================

// Comparison of segments. s1 < s2 if s1 is left of s2
// we order segments by decreasing tau_c (independent of tau_cdag !).
inline bool operator<(segment_t const &s1, segment_t const &s2) { return s1.tau_c > s2.tau_c; };

// Equality of segments based on tau_c, tau_cdag.
// NB : J_c, J_cdag are not important here
inline bool operator==(segment_t const &s1, segment_t const &s2) {
  return s1.tau_c == s2.tau_c and s1.tau_cdag == s2.tau_cdag;
};

// Whether a segment is wrapped around beta/0
// [c cdag] is not while [cdag c] is. NB : time in decreasing order, left to right
inline bool is_cyclic(segment_t const &seg) { return seg.tau_cdag > seg.tau_c; }

// Check whether segment is a full line.
inline bool is_full_line(segment_t const &seg) { return seg == segment_t::full_line(); }

// Overlap between two non-cyclic segments.
double overlap(segment_t const &s1, segment_t const &s2);

// Flip a segment. J are set to default
inline segment_t flip(segment_t const &s) { return {s.tau_cdag, s.tau_c}; }

// =================== Functions to manipulate std::vector<segment_t> ========

// lower_bound : find segment at tau if present or the first after tau
vec_seg_iter_t lower_bound(std::vector<segment_t> const &seglist, tau_t const &tau);

// Value of n (= 0 or 1) at tau = beta = 0
// 1 iif there is a cyclic segment or a full line
int n_at_boundary(std::vector<segment_t> const &seglist);

// Find density (0 or 1)in seglist to the right of time tau.
int n_tau(tau_t const &tau, std::vector<segment_t> const &seglist);

// Flip config
std::vector<segment_t> flip(std::vector<segment_t> const &sl);

// Overlap between segment and a list of segments.
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg);

// Checks if segment seg can be inserted into the list, i.e. without
// overlap with other segment.
bool is_insertable_into(segment_t const &seg, std::vector<segment_t> const &seglist);

// Find the indices of the segments whose cdag are in ]wtau_left,wtau_right[
std::vector<long> cdag_in_window(tau_t const &wtau_left, tau_t const &wtau_right,
                                 std::vector<segment_t> const &seglist);

// Fix the list after a change of operator c time in some move
// to restore the invariants
// 1 -if first segment is cyclic (c has move beyond beta),  put it at end
// 2- if last segment is such that its tau is > tau of first (c has moved beyond 0), put it first
// Only useful when # segments > 1
inline void fix_ordering_first_last(std::vector<segment_t> &sl) {
  if (sl.size() <= 1) return;
  if (is_cyclic(sl[0])) std::rotate(begin(sl), begin(sl) + 1, end(sl));
  if (sl.back().tau_c > sl[0].tau_c) std::rotate(begin(sl), end(sl) - 1, end(sl));
}

// Contribution of the dynamical interaction kernel K to the overlap between a segment and a list of segments.
// Computes the sum of the s_a s_b K(tau_a - tau_b) where s_a is 1 for cdag and - 1 for c
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau_c, tau_t const &tau_cdag,
                 gf<imtime, matrix_valued> const &K, int c1, int c2);

// Contribution of the dynamical interaction kernel K to the overlap between an operator and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau, bool is_c, gf<imtime, matrix_valued> const &K,
                 int c1, int c2);

// List of operators containing all colors.
std::vector<colored_ops_t> colored_ordered_ops(std::vector<std::vector<segment_t>> const &seglists);

// ===================  PRINTING ========================

std::ostream &operator<<(std::ostream &out, std::vector<segment_t> const &sl);

std::ostream &operator<<(std::ostream &out, configuration_t const &config);
template <> struct fmt::formatter<configuration_t> : ostream_formatter {};

std::ostream &operator<<(std::ostream &out, std::vector<colored_ops_t> const &col);
template <> struct fmt::formatter<std::vector<colored_ops_t>> : ostream_formatter {};
