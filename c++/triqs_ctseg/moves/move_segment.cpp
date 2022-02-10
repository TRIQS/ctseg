#include "move_segment.hpp"

namespace moves {

  // Checks if two segments overlap (even just at their boundaries) 
  bool move_segment::no_overlap(segment_t seg1, segment_t seg2) {
      if (seg1.tau_c >= seg2.tau_cdag && seg1.tau_c <= seg2.tau_c) return false; 
      if (seg1.tau_cdag >= seg2.tau_cdag && seg1.tau_c <= seg2.tau_c) return false; 
      return true; 
  }

  // Checks if a segment is movable to a color
  bool move_segment::is_movable(autstd::vector<segment_t> const &seglist,segment_t const &seg) {
    bool result = true; 
    if seglist.empty() return result;
    // If seg is cyclic, split it
    if (seg.tau_c < seg.tau_cdag) return is_movable(seglist,segment_t{params.beta,tau_cdag}) && is_movable(seglist,segment_t{tau_c,0});
    // Isolate last segment 
    segment_t last_seg = seglist.back();
    // In case last segment is cyclic, split it and check its overlap with seg
    if (last_seg.tau_c < last_seg.tau_cdag) {
        result = result && no_overlap(seg,segment_t{params.beta,last_seg.tau_cdag}) && no_overlap(seg,segment_t{last_seg.tau_c,0});
    }
    else result = result && no_overlap(seg,last_seg);
    // Check overlap of seg with the remainder of seglist
    auto ind = find_segment_left(seglist, seg);
    while (seglist[ind].tau_c < seg.tau_cdag && ind != --seglist.end()) {
        result = result && no_overlap(seglist[ind],seg);
        ++ind;
    }
    return result; 
    }

  double move_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT MOVE ================ \n");

    // ------------ Choice of segment and colors --------------
    // Select origin color
    origin_color = rng(data.n_color);
    auto &sl       = config.seglists[origin_color];
    SPDLOG_LOGGER_TRACE("Moving from color {}", origin_color);

    // If color has no segments, nothing to move
    if sl.empty() return 0; 

    // Select segment to move 
    origin_index = rng(sl.size());
    origin_segment = sl[origin_index];

    SPDLOG_LOGGER_TRACE("Moving c at {}, cdag at {}", origin_segment.tau_c, origin_segment.tau_cdag);

    // Select destination color 
    destination_color = rng(data.n_color - 1);
    if (destination_color >= origin_color) ++destination_color; 
    auto &dsl = config.seglists[destination_color];

    // Reject if chosen segment overlaps with destination color 
    if !is_movable(dsl,origin_segment) return 0; 
    destination_index = find_segment_left(dsl,origin_segment) + 1;

    // ------------  Trace ratio  -------------
        // FIXME : here we will need chemical potential, field, etc
    double ln_trace_ratio = 0;
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // FIXME

    // ------------  Proposition ratio ------------
    double prop_ratio = (dsl.size() + 1) / sl.size();

    SPDLOG_LOGGER_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double move_segment::accept() {

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    data.dets[color].complete_operation();
    // Regroup segments 
    auto &sl = config.seglists[color];
    auto &dsl = config.seglists[destination_color];
    dsl.insert(destination_index,origin_seg); 
    sl.erase(origin_index);
    

    // FIXME ??? SIGNE ???
    double sign_ratio = 1; 
    return sign_ratio;
  }

  //--------------------------------------------------
  void move_segment::reject() {
    SPDLOG_LOGGER_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    data.dets[color].reject_last_try();
  }
};
} // namespace moves
