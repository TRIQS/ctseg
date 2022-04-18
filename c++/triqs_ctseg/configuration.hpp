#pragma once
#include <vector>
//#include "params.hpp"
#include "types.hpp"

#include "dets.hpp"
#include "tau_t.hpp"

// The MC configuration and associated functions

// ------------------- Data structures -------------------

// The configuration is made of n_color ordered list of segments

// NB : all the time ordering are in DECREASING order,
// in agreement with the usual convertion

// --------------- Segment ----------------------
//
// A segment represent a couple (c, cdag), at given times in [0, beta]
// cyclic segment if of the form  [ cdag c], ordinary non cyclic [c cdag]
//
// Some segments are aligned between 2 colors, corresponding to a S^+, S^- operators
// The S operators will not be linked to Delta, the hybridization, but to J lines
// This information is stored also in the segment (J_c, J_cdag field).
//
// Special case : full lines.
// When a line has no operator, we need to take into account 2 states : empty or full
// The full line case is coded as a [beta, 0] segment, as it naturally yields the correct overlap
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

// simple alias
using vec_seg_iter_t = std::vector<segment_t>::const_iterator;

// Stores a couple of S+, S-
// FIXME : need to store the color, which is today just 0,1
struct jperp_line_t {
  tau_t tau_Sminus, tau_Splus; // times of the S-, S+
};

// --------------- Configuration ----------------------
//
// The configuration is a list of of segments for each color.

struct configuration_t {
  // A list of segments for each color.
  // NB ORDERED on tau_c of segments, with decreasing order.
  std::vector<std::vector<segment_t>> seglists;

  // List of Jperp lines, NOT ordered.
  std::vector<jperp_line_t> Jperp_list;

  // Construct from the number of colors
  configuration_t(int n_color) : seglists(n_color) {}

  // Accessor number of colors
  int n_color() const { return seglists.size(); }
};

// ------------------- Functions to manipulate segments --------------------------

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

// Checks if two segments are completely disjoint (accounting for boundaries)
inline bool disjoint(segment_t const &s1, segment_t const &s2) {
  return s1.tau_cdag > s2.tau_c or s2.tau_cdag > s1.tau_c;
}

// Check whether time is in segment.
bool tau_in_seg(tau_t const &tau, segment_t const &seg);

// Overlap between two non-cyclic segments.
double overlap_seg(segment_t const &seg1, segment_t const &seg2);

// ------------------- Functions to manipulate std::vector<segment_t> --------------------------

// Total lenght of segment
// FIXME : rename total_length ???
double density(std::vector<segment_t> const &seglist);

// Find density in seglist at time tau.
double n_tau(tau_t const &tau, std::vector<segment_t> const &seglist);

// Find index of first segment starting left of seg.tau_c.
vec_seg_iter_t find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg);

// FIXME : reverse order ...
// Overlap between segment and a list of segments.
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg);

// Checks if segment seg can be inserted into the list, i.e. without
// overlap with other segment.
// FIXME : rename ??
bool is_insertable(std::vector<segment_t> const &seglist, segment_t const &seg);

// FIXME : Comment. Parameters ?
// Contribution of the dynamical interaction kernel K to the overlap between a segment and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau_c, tau_t const &tau_cdag,
                 gf<imtime, matrix_valued> const &K, int c1, int c2);

// FIXME : Comment. Parameters ?
// Contribution of the dynamical interaction kernel K to the overlap between an operator and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau, bool is_c, gf<imtime, matrix_valued> const &K,
                 int c1, int c2);

// ------------------- Functions to manipulate config --------------------------

// Return the value of n = 0 or 1, for a given color, at tau = beta = 0
int n_at_boundary(configuration_t const &config, int color);

// Find segments corresponding to bosonic line
std::pair<vec_seg_iter_t, vec_seg_iter_t> find_spin_segments(int line_idx, configuration_t const &config);

// Flip config
std::vector<segment_t> flip(std::vector<segment_t> const &sl);

// Sign of a config
double config_sign(configuration_t const &config, std::vector<det_t> const &dets);

// Find the indices of the segments whose cdag are in ]wtau_left,wtau_right[
std::vector<long> cdag_in_window(tau_t const &wtau_left, tau_t const &wtau_right,
                                 std::vector<segment_t> const &seglist);

// ----------- DEBUG code --------------

// Print config
std::ostream &operator<<(std::ostream &out, configuration_t const &config);
