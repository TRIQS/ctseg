#include "configuration.hpp"
#include "logs.hpp"

// ------------------- Invariants ---------------------------

void check_invariant(configuration_t const &config, std::vector<det_t> const &dets) {
  for (auto const &[c, sl] : itertools::enumerate(config.seglists))
    if (not sl.empty()) {
      auto const &D = dets[c];
      for (int i = 0; i < sl.size() - 1; ++i) {
        ALWAYS_EXPECTS((sl[i].tau_cdag != sl[i].tau_c), "Degenerate segment in color {} at position {} in config \n{}",
                       c, i, config);
        ALWAYS_EXPECTS((sl[i].tau_cdag > sl[i + 1].tau_c), "Time order error in color {} at position {} in config \n{}",
                       c, i, config);
        ALWAYS_EXPECTS(not is_cyclic(sl[i]), "Segment in color {} at position {} should not by cyclic in config \n{}",
                       c, i, config); // only last segment can be cyclic
        ALWAYS_EXPECTS(
           (sl[i].tau_cdag == D.get_x(i).first),
           "Det error in color {}. The segment in position {} has tau_cdag = {}, while the det has {} in row {}. Config : \n{}",
           c, i, sl[i].tau_cdag, D.get_x(i).first, i, config);
        ALWAYS_EXPECTS(
           (sl[i].tau_c == D.get_y(i).first),
           "Det error in color {}. The segment in position {} has tau_c = {}, while the det has {} in column {}.  Config : \n{}",
           c, i, sl[i].tau_c, D.get_y(i).first, i, config);
      }
    }
  LOG("Invariants OK.");
}

// ---------------------------

// Make a list of time ordered (decreasing) operators

std::vector<std::tuple<qmc_time_t, int, bool>> make_time_ordered_op_list(configuration_t const &config) {

  std::vector<std::tuple<qmc_time_t, int, bool>> result;

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
    qmc_time_t tau_start = std::min(seg1.tau_c, seg2.tau_c);
    qmc_time_t tau_end   = std::max(seg1.tau_cdag, seg2.tau_cdag);
    return double(tau_start - tau_end); // FIXME: overlap of two full lines
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
double overlap(std::vector<segment_t> const &seglist, segment_t const &seg, qmc_time_factory_t const &fac) {
  if (seglist.empty()) return 0;

  auto zero = fac.get_lower_pt();
  auto beta = fac.get_upper_pt();

  // If seg is cyclic, split it
  if (is_cyclic(seg))
    return overlap(seglist, segment_t{beta, seg.tau_cdag}, fac) + overlap(seglist, segment_t{seg.tau_c, zero}, fac);
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
bool is_insertable(std::vector<segment_t> const &seglist, segment_t const &seg, qmc_time_factory_t fac) {
  bool result = true;
  if (seglist.empty()) return result;
  auto zero = fac.get_lower_pt();
  auto beta = fac.get_upper_pt();
  // If seg is cyclic, split it
  if (is_cyclic(seg))
    return is_insertable(seglist, segment_t{beta, seg.tau_cdag}, fac)
       and is_insertable(seglist, segment_t{seg.tau_c, zero}, fac);
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
double K_overlap(std::vector<segment_t> const &seglist, qmc_time_t const &tau_c, qmc_time_t const &tau_cdag,
                 gf<imtime, matrix_valued> const &K, int c1, int c2) {
  if (seglist.empty()) return 0;
  double result = 0;
  for (auto seg_in_list : seglist) {
    result += real(K(double(tau_c - seg_in_list.tau_c))(c1, c2) + K(double(tau_cdag - seg_in_list.tau_cdag))(c1, c2)
                   - K(double(tau_cdag - seg_in_list.tau_c))(c1, c2) - K(double(tau_c - seg_in_list.tau_cdag))(c1, c2));
  }
  return result;
}
// ---------------------------

// Length occupied by all segments for a given color
double density(std::vector<segment_t> const &seglist) {
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
    res[c] = (sl.empty() ? false : is_cyclic(sl.back())); // FIXME  BUG full line
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
  std::pair<std::vector<segment_t>::const_iterator, std::vector<segment_t>::const_iterator> result;
  result.first  = it_up;
  result.second = it_down;
  return result;
}

// ---------------------------

// Same as std::upper_bound, but i-th element of vector is returned by f[i]
long upper_bound(auto f, long N, auto const &value) {
  long first = 0, count = N;

  while (count > 0) {
    long step = count / 2;
    if (value < f(first + step))
      count = step;
    else {
      first += step + 1;
      count -= step + 1;
    }
  }
  return first;
}

// Print config
std::ostream &operator<<(std::ostream &out, configuration_t const &config) {
  for (auto const &[c, sl] : itertools::enumerate(config.seglists)) {
    out << '\n';
    for (auto const &[i, seg] : itertools::enumerate(sl))
      out << "Color " << c << ". Position " << i << " : [" << seg.tau_c << ", " << seg.tau_cdag << "]\n";
  }
  return out;
}
