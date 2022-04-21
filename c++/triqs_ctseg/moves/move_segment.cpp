#include "move_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double move_segment::attempt() {

    LOG("\n =================== ATTEMPT MOVE ================ \n");

    // ------------ Choice of segment and colors --------------
    // Select origin color
    origin_color = rng(config.n_color());
    LOG("Moving from color {}", origin_color);

    // Select destination color, different from origin_color
    dest_color = rng(config.n_color() - 1);
    if (dest_color >= origin_color) ++dest_color; // little trick to select another color
    LOG("Moving to color {}", dest_color);

    // Do we want to move an antisegment ?
    flipped = (rng(2) == 0);

    if (flipped) {
      // if we want to move an antisegment, we simply flip the configuration
      sl  = flip(config.seglists[origin_color]);
      dsl = flip(config.seglists[dest_color]);
      LOG("Moving antisegment.");
    } else {
      sl  = config.seglists[origin_color];
      dsl = config.seglists[dest_color];
    }

    // If color has no segments, nothing to move
    if (sl.empty()) {
      LOG("Nothing to move!");
      return 0;
    }

    // Select segment to move
    origin_index   = rng(sl.size());
    origin_segment = sl[origin_index];
    LOG("Moving segment at position {}", origin_index);

    // Reject if the segment has spin lines attached
    if (origin_segment.J_c or origin_segment.J_cdag) {
      LOG("Segment has spin line attached: cannot move.");
      return 0;
    }

    // Reject if chosen segment overlaps with destination color
    if (not is_insertable_into(origin_segment, dsl)) {
      LOG("Space is occupied in destination color.");
      return 0;
    }

    // Find where origin segment should be inserted in destination color
    dest_index = std::upper_bound(dsl.begin(), dsl.end(), origin_segment) - dsl.cbegin();
    LOG("Moving to position {}", dest_index);

    // ------------  Trace ratio  -------------

    double ln_trace_ratio =
       (flipped ? -1 : 1) * (wdata.mu(dest_color) - wdata.mu(origin_color)) * double(origin_segment.length());

    if (wdata.has_Dt) {
      auto tau_c    = origin_segment.tau_c;
      auto tau_cdag = origin_segment.tau_cdag;
      if (flipped) std::swap(tau_c, tau_cdag);

      for (auto const &[c, slist] : itertools::enumerate(config.seglists)) {
        ln_trace_ratio += K_overlap(slist, tau_c, tau_cdag, wdata.K, dest_color, c);
        ln_trace_ratio -= K_overlap(slist, tau_c, tau_cdag, wdata.K, origin_color, c);
      }
      auto Kl = wdata.K(double(origin_segment.length())); // matrix
      ln_trace_ratio -= real(Kl(origin_color, origin_color));
      ln_trace_ratio -= real(Kl(dest_color, dest_color));
      ln_trace_ratio += 2 * real(Kl(origin_color, dest_color));
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    auto compute_det_ratio = [](auto &D_orig, auto &D_dest, segment_t const &seg_orig) {
      return D_dest.try_insert(det_lower_bound_x(D_dest, seg_orig.tau_cdag), det_lower_bound_y(D_dest, seg_orig.tau_c),
                               {seg_orig.tau_cdag, 0}, {seg_orig.tau_c, 0})
         * D_orig.try_remove(det_lower_bound_x(D_orig, seg_orig.tau_cdag), det_lower_bound_y(D_orig, seg_orig.tau_c));
    };

    // Times are ordered in det. We insert tau_cdag as a line (first index) and tau_c as a column.
    // c and cdag are inverted if we flip
    double det_ratio = 1;
    if (not is_full_line(origin_segment))
      det_ratio = compute_det_ratio(wdata.dets[origin_color], wdata.dets[dest_color],
                                    (flipped ? flip(origin_segment) : origin_segment));

    // ------------  Proposition ratio ------------

    double prop_ratio = double(sl.size()) / (dsl.size() + 1);

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double move_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = config_sign(config, wdata.dets);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Update the dets
    wdata.dets[origin_color].complete_operation();
    wdata.dets[dest_color].complete_operation();

    // Add the segment at destination
    dsl.insert(begin(dsl) + dest_index, origin_segment);

    // Remove the segment at origin
    sl.erase(begin(sl) + origin_index);

    if (flipped) {
      config.seglists[origin_color] = flip(sl);
      config.seglists[dest_color]   = flip(dsl);
    } else {
      config.seglists[origin_color] = std::move(sl);
      config.seglists[dest_color]   = std::move(dsl);
    }
    // WARNING : do not use sl, dsl AFTER !

    double final_sign = config_sign(config, wdata.dets);
    double sign_ratio = final_sign / initial_sign;
    LOG("Final sign is {}", final_sign);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata.dets);

    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}. Config: ",
                   det_sign, sign_ratio, config);
    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void move_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[origin_color].reject_last_try();
    wdata.dets[dest_color].reject_last_try();
  }
}; // namespace moves
