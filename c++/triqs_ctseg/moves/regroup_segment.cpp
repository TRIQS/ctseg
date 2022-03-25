#include "regroup_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double regroup_segment::attempt() {

    LOG("\n =================== ATTEMPT REGROUP ================ \n");

    // ------------ Choice of segment --------------
    // Select color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    LOG("Regrouping at color {}", color);

    // If no segments nothing to regroup
    if (sl.empty()) {
      LOG("Color is empty");
      return 0;
    }

    // Select pair of segments (or cyclic segment) to regroup
    making_full_line = sl.size() == 1;
    if (making_full_line) {
      if (is_full_line(sl[0], fac)) {
        LOG("Segment is full line");
        return 0; // If segment is full line nothing to regroup
      }
      left_segment_index  = 0;
      right_segment_index = 0;
    } else {
      left_segment_index  = rng(sl.size());
      right_segment_index = (left_segment_index == sl.size() - 1) ? 0 : left_segment_index + 1;
    }
    left_segment  = sl[left_segment_index];
    right_segment = sl[right_segment_index];

    auto inserted_segment =
       segment_t{left_segment.tau_cdag, right_segment.tau_c}; // "antisegment" : careful with order of c, cdag
    auto inserted_segment_length = double(inserted_segment.tau_c - inserted_segment.tau_cdag);

    LOG("Regroup: removing c at {}, cdag at {}", right_segment.tau_c, left_segment.tau_cdag);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(color, c) * overlap(config.seglists[c], inserted_segment, fac);
        ln_trace_ratio += wdata.mu(c) * inserted_segment_length;
        if (wdata.has_Dt)
          ln_trace_ratio += K_overlap(config.seglists[c], inserted_segment, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    // We remove a cdag (first index) from the left segment and a c (second index) from the right segment.
    auto det_ratio = wdata.dets[color].try_remove(left_segment_index, right_segment_index);
    // additional sign due to "roll" is not computed at this stage, cf accept

    // ------------  Proposition ratio ------------

    double future_number_segments   = making_full_line ? 1 : int(sl.size()) - 1;
    double current_number_intervals = sl.size();
    // Length of future segment
    qmc_time_t l = wdata.qmc_beta;
    if (not making_full_line) l = left_segment.tau_c - right_segment.tau_cdag;
    double prop_ratio = current_number_intervals / (future_number_segments * l * l / (making_full_line ? 1 : 2));

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    double prod = trace_ratio * det_ratio * prop_ratio;

    return (std::isfinite(prod)) ? prod : det_sign;
  }

  //--------------------------------------------------

  double regroup_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Update the dets
    wdata.dets[color].complete_operation();

    // Regroup segments and compute the sign ratio
    auto &sl          = config.seglists[color];
    double sign_ratio = 1;
    if (making_full_line) {
      if (is_cyclic(left_segment)) sign_ratio = -1;
      sl[left_segment_index] = segment_t{wdata.qmc_beta, wdata.qmc_zero};
    } else {
      auto new_segment           = segment_t{left_segment.tau_c, right_segment.tau_cdag};
      bool had_cyclic_segment    = is_cyclic(left_segment) or is_cyclic(right_segment);
      bool adding_cyclic_segment = !had_cyclic_segment and is_cyclic(new_segment);
      if (adding_cyclic_segment) sign_ratio = -1;
      bool need_roll = is_cyclic(left_segment) or adding_cyclic_segment;
      if (need_roll) {
        // In this case, a cdag that belonged to the first segment now belongs to the last segment:
        // we need to "roll up" the det and this produces an additional sign.
        sign_ratio *= wdata.dets[color].roll_matrix(det_t::Up);
      }
      // Update the left segment
      sl[left_segment_index] = new_segment;
      // Remove the "other" segment
      sl.erase(sl.begin() + right_segment_index);
    }
    LOG("Sign ratio is {}", sign_ratio);

    // Check invariant
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}.", det_sign,
                   sign_ratio);

#ifdef CHECK_INVARIANTS
    check_invariant(config, wdata.dets);
#endif
    SPDLOG_TRACE("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void regroup_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
} // namespace moves
