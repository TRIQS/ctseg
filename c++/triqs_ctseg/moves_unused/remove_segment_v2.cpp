#include "remove_segment_v2.hpp"
#include "../logs.hpp"
#include <cmath>

namespace moves {

  double remove_segment_v2::attempt() {

    LOG("\n =================== ATTEMPT REMOVE ================ \n");

    // Select removal color
    color = rng(config.n_color());
    LOG("Removing at color {}", color);
    need_flip = false;

    current_density = density(config.seglists[color]);
    if (rng() < 1 - current_density / wdata.beta) {
      //if (rng() < 0.5) {
      need_flip = true;
      sl        = flip(config.seglists[color], wdata.beta);
      LOG("Removing antisegment.");
    } else {
      sl = config.seglists[color];
    }

    // If color is empty, nothing to remove
    if (sl.empty()) {
      LOG("Color is empty.");
      return 0;
    }

    // Select segment to remove
    prop_seg_idx = rng(sl.size());
    prop_seg     = sl[prop_seg_idx];
    if (is_full_line(prop_seg)) {
      LOG("Cannot remove full line.");
      return 0;
    }
    if (prop_seg.J_c or prop_seg.J_cdag) {
      LOG("Segement is connected to a spin line, cannot remove.");
      return 0;
    }

    if (need_flip)
      LOG("Removing antisegment at position {} : c at {}, cdag at {}", prop_seg_idx, prop_seg.tau_cdag, prop_seg.tau_c);
    else
      LOG("Removing segment at position {} : c at {}, cdag at {}", prop_seg_idx, prop_seg.tau_c, prop_seg.tau_cdag);

    // ------------  Trace ratio  -------------

    double trace_sign     = need_flip ? -1 : 1;
    double ln_trace_ratio = -trace_sign * wdata.mu(color) * prop_seg.length();
    for (auto c : range(config.n_color())) {
      if (c != color) { ln_trace_ratio -= -trace_sign * wdata.U(color, c) * overlap(config.seglists[c], prop_seg); }
      if (wdata.has_Dt)
        ln_trace_ratio -= (c == color ? 1 : trace_sign)
           * K_overlap(config.seglists[c], prop_seg.tau_c, prop_seg.tau_cdag, wdata.K, color, c);
    }
    if (wdata.has_Dt)
      ln_trace_ratio +=
         -real(wdata.K(double(prop_seg.tau_c - prop_seg.tau_cdag))(color, color)); // Correct double counting
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    auto det_c_time    = [&](long i) { return wdata.dets[color].get_y(i).first; };
    auto det_cdag_time = [&](long i) { return wdata.dets[color].get_x(i).first; };
    long det_index_c = 0, det_index_cdag = 0;
    // Times are ordered in det. We insert tau_cdag as a line (first index) and tau_c as a column.
    // c and cdag are inverted if we flip
    if (need_flip) {
      det_index_c    = lower_bound(det_c_time, wdata.dets[color].size(), prop_seg.tau_cdag);
      det_index_cdag = lower_bound(det_cdag_time, wdata.dets[color].size(), prop_seg.tau_c);
    } else {
      det_index_c    = lower_bound(det_c_time, wdata.dets[color].size(), prop_seg.tau_c);
      det_index_cdag = lower_bound(det_cdag_time, wdata.dets[color].size(), prop_seg.tau_cdag);
    }
    LOG("Det: removing c at position {}, cdag at position {}.", det_index_c, det_index_cdag);
    double det_ratio = wdata.dets[color].try_remove(det_index_cdag, det_index_c);

    // ------------  Proposition ratio ------------

    double future_density = current_density - prop_seg.length();
    double density_ratio  = (wdata.beta - future_density) / (current_density);
    if (need_flip) {
      future_density = wdata.beta - future_density;
      density_ratio  = (future_density) / (wdata.beta - current_density);
    }
    //density_ratio                  = 1;
    double current_number_segments = sl.size();
    double future_number_intervals = std::max(1.0, sl.size() - 1.0);
    // Insertion window for reverse move
    tau_t wtau_left, wtau_right;
    if (current_number_segments != 1) {
      bool is_last_segment  = prop_seg_idx == sl.size() - 1;
      bool is_first_segment = prop_seg_idx == 0;
      // Look at segments on either side of prop_seg, accounting for cyclicity
      wtau_right = sl[is_last_segment ? 0 : prop_seg_idx + 1].tau_c;
      wtau_left  = sl[is_first_segment ? sl.size() - 1 : prop_seg_idx - 1].tau_cdag;
    }
    tau_t window_length = (current_number_segments == 1) ? tau_t::beta() : wtau_left - wtau_right;
    double prop_ratio   = density_ratio * current_number_segments
       / (future_number_intervals * window_length * window_length / (current_number_segments == 1 ? 1 : 2));
    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    //det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    det_sign    = 1;
    double prod = trace_ratio * abs(det_ratio) * prop_ratio;

    return (std::isfinite(prod)) ? prod : det_sign;
  }

  //--------------------------------------------------

  double remove_segment_v2::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Insert the times into the det
    wdata.dets[color].complete_operation();

    // Compute the sign ratio
    double sign_ratio = 1;
    LOG("Sign ratio is {}", sign_ratio);

    // Remove the segment
    sl.erase(sl.begin() + prop_seg_idx);
    if (need_flip) {
      config.seglists[color] = flip(sl, wdata.beta);
    } else
      config.seglists[color] = sl;

      // Check invariant
#ifdef CHECK_INVARIANTS
    check_invariant(config, wdata.dets);
#endif
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}.", det_sign,
                   sign_ratio);
    //LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void remove_segment_v2::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
