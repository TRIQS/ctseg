#include "move_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double move_segment::attempt() {

    LOG("\n =================== ATTEMPT MOVE ================ \n");

    // ------------ Choice of segment and colors --------------
    // Select origin color
    origin_color = rng(wdata.n_color);
    auto &sl     = config.seglists[origin_color];
    LOG("Moving from color {}", origin_color);

    // If color has no segments, nothing to move
    if (sl.empty()) {
      LOG("Nothing to move!");
      return 0;
    }

    // Select segment to move
    origin_index   = rng(sl.size());
    origin_segment = sl[origin_index];
    LOG("Moving segment at position {}", origin_index);

    // Select destination color
    destination_color = rng(wdata.n_color - 1);
    if (destination_color >= origin_color) ++destination_color;
    auto &dsl = config.seglists[destination_color];
    LOG("Moving to color {}", destination_color);

    // Reject if chosen segment overlaps with destination color
    if (not is_insertable(dsl, origin_segment)) {
      LOG("Space is occupied in destination color.");
      return 0;
    }

    // Find where origin segment should be inserted in destination color
    destination_it  = std::upper_bound(dsl.begin(), dsl.end(), origin_segment);
    long dest_index = std::distance(dsl.cbegin(), destination_it);
    LOG("Moving to position {}", dest_index);

    // ------------  Trace ratio  -------------

    double ln_trace_ratio = (wdata.mu(destination_color) - wdata.mu(origin_color)) * origin_segment.length();
    double trace_ratio    = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
    // segment the tau_c/tau_cdag belongs to.
    double det_ratio      = 1.0;
    bool moving_full_line = is_full_line(origin_segment);
    if (not moving_full_line) {
      det_ratio = wdata.dets[destination_color].try_insert(dest_index, dest_index, {origin_segment.tau_cdag, 0},
                                                           {origin_segment.tau_c, 0})
         * wdata.dets[origin_color].try_remove(origin_index, origin_index);
    }

    // ------------  Proposition ratio ------------
    double prop_ratio = int(sl.size()) / (int(dsl.size()) + 1);

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double move_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Update the dets
    wdata.dets[origin_color].complete_operation();
    wdata.dets[destination_color].complete_operation();

    double sign_ratio = 1;

    // Add the segment at destination
    auto &dsl = config.seglists[destination_color];
    dsl.insert(destination_it, origin_segment);

    // Remove the segment at origin
    auto &sl = config.seglists[origin_color];
    sl.erase(sl.begin() + origin_index);

    // Check invariant
#ifdef CHECK_INVARIANTS
    check_invariant(config, wdata.dets);
#endif
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}.", det_sign,
                   sign_ratio);
    SPDLOG_TRACE("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void move_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[origin_color].reject_last_try();
    wdata.dets[destination_color].reject_last_try();
  }
}; // namespace moves
