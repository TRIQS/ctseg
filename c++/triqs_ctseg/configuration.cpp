#include "configuration.hpp"
#include "logs.hpp"

// Make a list of time ordered (decreasing) operators

std::vector<std::tuple<tau_t, int, bool>> make_time_ordered_op_list(configuration_t const &config) {

  std::vector<std::tuple<tau_t, int, bool>> result;

  result.reserve(config.seglists.size()
                 * config.seglists[0].size()); // optional : simple heuristics to reserve some space

  for (auto const &[color, seglist] : itertools::enumerate(config.seglists))
    for (auto const &seg : seglist) {
      result.emplace_back(seg.tau_c, color, false);
      result.emplace_back(seg.tau_cdag, color, true);
    }

  std::sort(std::begin(result), std::end(result), std::greater{});
  return result;
}

// ---------------------------

// Find index of first segment starting left of seg.tau_c.
std::vector<segment_t>::const_iterator find_segment_left(std::vector<segment_t> const &seglist, segment_t const &seg) {
  auto seg_iter = std::upper_bound(seglist.begin(), seglist.end(), seg);
  return (seg_iter == seglist.begin()) ? seg_iter : --seg_iter;
};

// ---------------------------

// Overlap between two non-cyclic segments.
double overlap_seg(segment_t const &seg1, segment_t const &seg2) {
  if (seg1.tau_cdag >= seg2.tau_c or seg2.tau_cdag >= seg1.tau_c)
    return 0;
  else {
    tau_t tau_start = std::min(seg1.tau_c, seg2.tau_c);
    tau_t tau_end   = std::max(seg1.tau_cdag, seg2.tau_cdag);
    return double(tau_start - tau_end);
  };
};

// ---------------------------

// Checks if two segments are completely disjoint (accounting for boundaries)
bool disjoint(segment_t seg1, segment_t seg2) {
  if (seg1.tau_cdag > seg2.tau_c or seg2.tau_cdag > seg1.tau_c) return true;
  return false;
}

// ---------------------------

// Overlap between segment and a list of segments.
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg) {
  if (seglist.empty()) return 0;

  auto beta = seg.tau_c.beta();
  auto zero = seg.tau_c.zero();

  // If seg is cyclic, split it
  if (is_cyclic(seg))
    return overlap(seglist, segment_t{beta, seg.tau_cdag}) + overlap(seglist, segment_t{seg.tau_c, zero});
  double result = 0;
  // Isolate last segment
  auto last_seg = seglist.back();
  // In case last segment is cyclic, split it and compute its overlap with seg
  if (is_cyclic(last_seg))
    result += overlap_seg(seg, segment_t{beta, last_seg.tau_cdag}) + overlap_seg(seg, segment_t{last_seg.tau_c, zero});
  else
    result += overlap_seg(seg, last_seg);

  // Compute overlap of seg with the remainder of seglist
  for (auto it = find_segment_left(seglist, seg); it->tau_c > seg.tau_cdag && it != --seglist.end(); ++it) //
    result += overlap_seg(*it, seg);
  return result;
};

// ---------------------------

// Checks if segment is movable to a given color
bool is_insertable(std::vector<segment_t> const &seglist, segment_t const &seg) {
  bool result = true;
  if (seglist.empty()) return result;
  auto beta = seg.tau_c.beta();
  auto zero = seg.tau_c.zero();
  // If seg is cyclic, split it
  if (is_cyclic(seg))
    return is_insertable(seglist, segment_t{beta, seg.tau_cdag}) and is_insertable(seglist, segment_t{seg.tau_c, zero});
  // In case last segment in list is cyclic, split it and check its overlap with seg
  if (is_cyclic(seglist.back())) {
    result = result and disjoint(seg, segment_t{beta, seglist.back().tau_cdag})
       and disjoint(seg, segment_t{seglist.back().tau_c, zero});
  } else
    result = result and disjoint(seg, seglist.back());
  // Check overlap of seg with the remainder of seglist
  for (auto it = find_segment_left(seglist, seg); it->tau_c >= seg.tau_cdag and it != --seglist.end(); ++it) {
    result = result and disjoint(*it, seg);
  }
  return result;
}

// ---------------------------

// Contribution of the dynamical interaction kernel K to the overlap between a segment and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau_c, tau_t const &tau_cdag,
                 gf<imtime, matrix_valued> const &K, int c1, int c2) {
  if (seglist.empty()) return 0;
  double result = 0;
  for (auto seg_in_list : seglist) {
    result += real(K(double(tau_c - seg_in_list.tau_c))(c1, c2) + K(double(tau_cdag - seg_in_list.tau_cdag))(c1, c2)
                   - K(double(tau_cdag - seg_in_list.tau_c))(c1, c2) - K(double(tau_c - seg_in_list.tau_cdag))(c1, c2));
  }
  return result;
}

// Contribution of the dynamical interaction kernel K to the overlap between an operator and a list of segments.
double K_overlap(std::vector<segment_t> const &seglist, tau_t const &tau, bool is_c, gf<imtime, matrix_valued> const &K,
                 int c1, int c2) {
  if (seglist.empty()) return 0;
  double result = 0;
  for (auto seg_in_list : seglist) {
    result += real(K(double(tau - seg_in_list.tau_c))(c1, c2) - K(double(tau - seg_in_list.tau_cdag))(c1, c2));
  }
  return is_c ? result : -result;
}
// ---------------------------

// Length occupied by all segments for a given color
double density(std::vector<segment_t> const &seglist) {
  if (seglist.empty()) return 0;
  double result = 0;
  for (auto const &seg : seglist) result += double(seg.tau_c - seg.tau_cdag);
  return result;
}
// ---------------------------

// Find the state at tau = 0 or beta
std::vector<bool> boundary_state(configuration_t const &config) {
  int N = config.n_color();
  std::vector<bool> res(N);
  for (auto const &[c, sl] : itertools::enumerate(config.seglists))
    res[c] = (sl.empty() ? false : is_cyclic(sl.back()) or is_full_line(sl.back()));
  return res;
}

// ---------------------------

// Find segments corresponding to bosonic line
std::pair<std::vector<segment_t>::const_iterator, std::vector<segment_t>::const_iterator>
find_spin_segments(int line_idx, configuration_t const &config) {
  auto const &line    = config.Jperp_list[line_idx];
  auto const &sl_up   = config.seglists[0];
  auto const &sl_down = config.seglists[1];
  // In spin up line, the c conneted to the J line is a at tau_Sminus
  auto c_up  = segment_t{line.tau_Sminus, line.tau_Sminus};
  auto it_up = std::lower_bound(sl_up.cbegin(), sl_up.cend(), c_up);
  // In spin down line, the c conneted to the J line is a at tau_Splus
  auto c_down  = segment_t{line.tau_Splus, line.tau_Splus};
  auto it_down = std::lower_bound(sl_down.cbegin(), sl_down.cend(), c_down);
  return std::make_pair(it_up, it_down);
}

// ---------------------------

// Flip seglist
std::vector<segment_t> flip(std::vector<segment_t> const &sl) {
  auto fsl = std::vector<segment_t>{};
  fsl.reserve(sl.size() + 1);
  if (sl.empty()) // Flipped seglist is full line
    fsl.emplace_back(segment_t{tau_t::beta(), tau_t::zero()});
  else if (sl.size() == 1 and is_full_line(sl[0])) { // Do nothing: flipped config empty
  } else {                                           // Swap c and cdag
    for (auto i : range(sl.size())) {
      if (is_cyclic(sl.back())) {
        long ind = (i == 0) ? long(sl.size()) - 1 : i - 1;
        fsl.emplace_back(segment_t{sl[ind].tau_cdag, sl[i].tau_c, sl[ind].J_cdag, sl[i].J_c});
      } else {
        long ind = (i == sl.size() - 1) ? 0 : i + 1;
        fsl.emplace_back(segment_t{sl[i].tau_cdag, sl[ind].tau_c, sl[i].J_cdag, sl[ind].J_c});
      }
    } // loop over segs
  }   // general case
  return fsl;
}

// ---------------------------

// Sign of a configuration
double config_sign(configuration_t const &config, std::vector<det_t> const &dets) {
  double sign = 1.0;
  // For every color we compute the sign of the permutation that takes
  // [(c_dag c) (c_dag c) (c_dag c) ...] with the cdag in increasing time order
  // to the completely time-ordered list of operators (with increasing time)
  for (auto c : range(config.n_color())) {
    auto s = long(dets[c].size());
    if (s != 0) {
      // We first compute the sign of the permutation that takes
      // [(c_dag c) (c_dag c) (c_dag c) ...] with the cdag time-ordered to
      // [(c_dag c_dag ... cdag)(c c ... c)] with the c and c_dag time-ordered
      if ((s * (s - 1) / 2) % 2 == 1) sign *= -1;
      // We then compute the sign of the permutation that takes
      // [(c_dag c_dag ... cdag)(c c ... c)] with the c and c_dag time-ordered
      // to the completely time-ordered list of operators
      int idx_c = 0, idx_cdag = 0;
      for (int n = 0; n < 2 * s - 1; ++n) {
        if (dets[c].get_x(idx_cdag).first < dets[c].get_y(idx_c).first) {
          if (idx_cdag < s - 1) ++idx_cdag;
        } else {
          // Count the number of transpositions
          if ((s - idx_cdag) % 2 == 1) sign *= -1;
          if (idx_c < s - 1) ++idx_c;
        }
      }
    }
  }
  return sign;

  /*   
    // Old sign computation, works only without J_perp
    for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
    if (not sl.empty()) {
      bool starts_with_dagger = false;
      auto s                  = dets[c].size();
      if (s != 0) starts_with_dagger = dets[c].get_x(s - 1) > dets[c].get_y(s - 1);
      if (starts_with_dagger) sign *= (s % 2 == 0) ? 1 : -1; 

    }
  } */
}

// ---------------------------

// Find the indices of the segments whose cdag are in ]wtau_left,wtau_right[
std::vector<long> cdag_in_window(tau_t const &wtau_left, tau_t const &wtau_right,
                                 std::vector<segment_t> const &seglist) {
  std::vector<long> found_indices;
  if (seglist.empty()) return found_indices; // should never happen, but protect
  found_indices.reserve(seglist.size());
  for (auto it = find_segment_left(seglist, segment_t{wtau_left, wtau_left});
       it->tau_cdag > wtau_right and it != --seglist.end(); ++it) {
    if (it->tau_cdag < wtau_left) found_indices.push_back(std::distance(seglist.cbegin(), it));
  }
  // Check separately for last segment (may be cyclic)
  if (seglist.back().tau_cdag < wtau_left and seglist.back().tau_cdag > wtau_right)
    found_indices.push_back(seglist.size() - 1);
  return found_indices;
}

// ---------------------------

// Print config
std::ostream &operator<<(std::ostream &out, configuration_t const &config) {
  for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
    out << '\n';
    for (auto const &[i, seg] : itertools::enumerate(sl))
      out << "Color " << c << ". Position " << i << " : [ J:" << seg.J_c << " " << seg.tau_c << ", " << seg.tau_cdag
          << " J:" << seg.J_cdag << "]\n";
  }
  out << "\nSpin lines : \n";
  for (auto const &[i, line] : itertools::enumerate(config.Jperp_list)) {
    out << "S_minus : [" << line.tau_Sminus << "] S_plus : [" << line.tau_Splus << "]\n";
  }
  return out;
}
