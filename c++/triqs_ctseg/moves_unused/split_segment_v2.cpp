#include "split_segment_v2.hpp"
#include "../logs.hpp"
#include <cmath>

namespace moves {

  double split_segment_v2::attempt() {

    LOG("\n =================== ATTEMPT INSERT ================ \n");

    // Select insertion color
    color = rng(wdata.n_color);
    LOG("Inserting at color {}", color);
    need_flip = false;

    current_density = density(config.seglists[color]);
    if (rng() < 1.0) {
      need_flip = true;
      sl        = flip(config.seglists[color], wdata.beta);
      LOG("Inserting antisegment.");
    } else {
      sl = config.seglists[color];
    }

    LOG("Flip = {}", need_flip);
    // Select insertion window [wtau_left,wtau_right]
    tau_t wtau_left = tau_t::beta(), wtau_right = tau_t::zero();

    if (not sl.empty()) {
      if (is_full_line(sl.back())) {
        LOG("Full line, cannot insert.");
        return 0;
      }
      // Randomly choose one existing segment
      long seg_idx         = rng(sl.size());
      wtau_left            = sl[seg_idx].tau_cdag; // wtau_left is cdag of this segment
      bool is_last_segment = seg_idx == sl.size() - 1;
      wtau_right = sl[is_last_segment ? 0 : seg_idx + 1].tau_c; // wtau_right is c of next segment, possibly cyclic
    }

    LOG("Insertion window is wtau_left = {}, wtau_right = {}", wtau_left, wtau_right);
    tau_t window_length = sl.empty() ? tau_t::beta() : wtau_left - wtau_right;

    // Choose two random times in insertion window
    auto dt1 = tau_t::random(rng, window_length);
    auto dt2 = tau_t::random(rng, window_length);
    if (dt1 == dt2) {
      LOG("Generated equal times");
      return 0;
    }
    if (dt1 > dt2 and !sl.empty()) std::swap(dt1, dt2); // if inserting into an empty line, two ways to insert
    prop_seg = segment_t{wtau_left - dt1, wtau_left - dt2};
    // Iterator pointing to prop_seg if it is inserted in the list of segments.
    prop_seg_it = std::upper_bound(sl.begin(), sl.end(), prop_seg);

    if (need_flip)
      LOG("Inserting antisegment at position {}, with c at {}, cdag at {}.", std::distance(sl.begin(), prop_seg_it),
          prop_seg.tau_cdag, prop_seg.tau_c);
    else
      LOG("Inserting segment at position {}, with c at {}, cdag at {}.", std::distance(sl.begin(), prop_seg_it),
          prop_seg.tau_c, prop_seg.tau_cdag);

    // ------------  Trace ratio  -------------

    double trace_sign     = need_flip ? -1 : 1;
    double ln_trace_ratio = trace_sign * wdata.mu(color) * prop_seg.length();
    for (auto c : range(wdata.n_color)) {
      if (c != color) { ln_trace_ratio += -trace_sign * wdata.U(color, c) * overlap(config.seglists[c], prop_seg); }
      if (wdata.has_Dt)
        ln_trace_ratio += (c == color ? 1 : trace_sign)
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
    double det_ratio = 0;
    // Times are ordered in det. We insert tau_cdag as a line (first index) and tau_c as a column.
    // c and cdag are inverted if we flip
    if (need_flip) {
      det_index_c    = lower_bound(det_c_time, wdata.dets[color].size(), prop_seg.tau_cdag);
      det_index_cdag = lower_bound(det_cdag_time, wdata.dets[color].size(), prop_seg.tau_c);
      det_ratio =
         wdata.dets[color].try_insert(det_index_cdag, det_index_c, {prop_seg.tau_c, 0}, {prop_seg.tau_cdag, 0});
    } else {
      det_index_c    = lower_bound(det_c_time, wdata.dets[color].size(), prop_seg.tau_c);
      det_index_cdag = lower_bound(det_cdag_time, wdata.dets[color].size(), prop_seg.tau_cdag);
      det_ratio =
         wdata.dets[color].try_insert(det_index_cdag, det_index_c, {prop_seg.tau_cdag, 0}, {prop_seg.tau_c, 0});
    }
    LOG("Det: inserting c at position {}, cdag at position {}.", det_index_c, det_index_cdag);

    // ------------  Proposition ratio ------------

    /* double future_density = current_density + prop_seg.length();
    double density_ratio  = (wdata.beta - future_density) / (wdata.beta - current_density);
    if (need_flip) {
      future_density = wdata.beta - future_density;
      density_ratio  = future_density / current_density;
    } */
    double density_ratio            = 1.0;
    double current_number_intervals = std::max(1.0, double(sl.size()));
    double future_number_segments   = double(sl.size()) + 1;
    double prop_ratio               = density_ratio
       * (current_number_intervals * window_length * window_length / (sl.empty() ? 1 : 2)) / future_number_segments;
    // Account for absence of time swapping when inserting into empty line.

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * abs(det_ratio) * prop_ratio;
    //det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    det_sign = 1;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double split_segment_v2::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Insert the times into the det
    wdata.dets[color].complete_operation();

    // Compute the sign ratio
    double sign_ratio = 1;
    LOG("Sign ratio is {}", sign_ratio);

    // Insert the segment in an ordered list
    sl.insert(prop_seg_it, prop_seg);
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
  void split_segment_v2::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
