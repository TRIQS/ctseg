#include "remove_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double remove_segment::attempt() {

    LOG("\n =================== ATTEMPT REMOVE ================ \n");

    // ------------ Choice of segment --------------

    // Select removal color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    LOG("Removing at color {}", color);

    // If color is empty, nothing to remove
    if (sl.empty()) {
      LOG("Color is empty.");
      return 0;
    }

    // Select segment to remove
    prop_seg_idx = rng(sl.size());
    prop_seg     = sl[prop_seg_idx];
    if (is_full_line(prop_seg, fac)) {
      LOG("Cannot remove full line.");
      return 0;
    }

    LOG("Removing segment at position {} : c at {}, cdag at {}", prop_seg_idx, prop_seg.tau_c, prop_seg.tau_cdag);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], prop_seg, fac);
        ln_trace_ratio -= wdata.mu(c) * prop_seg.length();
        if (wdata.has_Dt)
          ln_trace_ratio -= K_overlap(config.seglists[c], prop_seg, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    auto det_ratio = wdata.dets[color].try_remove(prop_seg_idx, prop_seg_idx);

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_intervals = std::max(1.0, sl.size() - 1.0);
    // Insertion window for reverse move, initialise at (beta,0)
    qmc_time_t wtau_left, wtau_right;
    if (current_number_segments != 1) {
      bool is_last_segment  = prop_seg_idx == sl.size() - 1;
      bool is_first_segment = prop_seg_idx == 0;
      // Look at segments on either side of prop_seg, accounting for cyclicity
      wtau_right = sl[is_last_segment ? 0 : prop_seg_idx + 1].tau_c;
      wtau_left  = sl[is_first_segment ? sl.size() - 1 : prop_seg_idx - 1].tau_cdag;
    }
    qmc_time_t window_length = (current_number_segments == 1) ? wdata.qmc_beta : wtau_left - wtau_right;
    double prop_ratio        = current_number_segments
       / (future_number_intervals * window_length * window_length / (current_number_segments == 1 ? 1 : 2));
    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    double prod = trace_ratio * det_ratio * prop_ratio;

    return (std::isfinite(prod)) ? prod : det_sign;
  }

  //--------------------------------------------------

  double remove_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Update the dets
    wdata.dets[color].complete_operation();

    auto &sl = config.seglists[color];
    // Compute the sign ratio
    double sign_ratio = is_cyclic(sl[prop_seg_idx]) ? -1 : 1;
    LOG("Sign ratio is {}", sign_ratio);

    // Remove the segment
    sl.erase(sl.begin() + prop_seg_idx);

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
  void remove_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
} // namespace moves
