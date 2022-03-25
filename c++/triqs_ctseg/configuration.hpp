#pragma once
#include <vector>
#include "params.hpp"
#include "types.hpp"

// Segment: (time of c, time of c^dagger)
struct segment_t {
  qmc_time_t tau_c, tau_cdag; // time of c and cdag
  bool J_at_start = false, J_at_end = false;
  [[nodiscard]] double length() const {
    auto length = double(tau_c - tau_cdag);
    return (length == 0) ? double(tau_c) : length;
  };
};

// Whether a segment is wrapped around beta/0
inline bool is_cyclic(segment_t const &seg) { return seg.tau_cdag > seg.tau_c; }

struct jperp_line_t {
  qmc_time_t tau_Splus, tau_Sminus; // times of the S+, S-
  //bool Splus_at_left;               /// USEFUL ???
};

// ----------

struct configuration_t {
  std::vector<std::vector<segment_t>>
     seglists; // list of segment per color : seglist[color] is ORDERED on tau, with decreasing order.
  std::vector<jperp_line_t> Jperp_list;

  configuration_t(int n_color) : seglists(n_color) {}

  [[nodiscard]] int n_color() const { return seglists.size(); }
};

// ------------------- Invariants ---------------------------

void check_invariant(configuration_t const &config, std::vector<det_t> const &dets);

// ------------------- Functions to manipulate config --------------------------

// Make a list of time ordered (decreasing) operators
// vector of (time, color, is_dagger)
std::vector<std::tuple<qmc_time_t, int, bool>> make_time_ordered_op_list(configuration_t const &config);

// Comparison of segments. Returns 1 if s1 is left of s2 (we order segments by decreasing time).
inline bool operator<(segment_t const &s1, segment_t const &s2) { return s1.tau_c > s2.tau_c; };

// Overlap between two non-cyclic segments.
double overlap_seg(segment_t const &seg1, segment_t const &seg2);

// Find index of first segment starting left of seg.tau_c.
inline auto find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg) {
  auto seg_iter = std::upper_bound(seglist.begin(), seglist.end(), seg);
  return (seg_iter == seglist.begin()) ? seg_iter : --seg_iter;
};

// Overlap between segment and a list of segments.
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg, qmc_time_factory_t const &fac);

// Contribution of the dynamical interaction kernel K to the overlap between a segment and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, segment_t const &seg,
                 gf_const_view<imtime, scalar_valued> const &K);

// Length occupied by all segments for a given color
double density(std::vector<segment_t> const &seglist);

// Check whether segment is a full line.
inline bool is_full_line(segment_t const &seg, qmc_time_factory_t const &fac) {
  auto qmc_zero = fac.get_lower_pt();
  auto qmc_beta = fac.get_upper_pt();
  return seg.tau_cdag == qmc_zero && seg.tau_c == qmc_beta;
}

/// Find the state at tau = 0 or beta
std::vector<bool> boundary_state(configuration_t const &config);

// --------- DEBUG code --------------
// print config + h5 config

std::ostream &operator<<(std::ostream &out, configuration_t const &config);
