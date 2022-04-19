#include "move_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double move_segment::attempt() {

    LOG("\n =================== ATTEMPT MOVE ================ \n");

    // ------------ Choice of segment and colors --------------
    // Select origin color
    origin_color = rng(config.n_color());
    LOG("Moving from color {}", origin_color);

    // Select destination color
    destination_color = rng(config.n_color() - 1);
    if (destination_color >= origin_color) ++destination_color;
    LOG("Moving to color {}", destination_color);

    need_flip = false;

    //double current_density = density(config.seglists[origin_color]);
    //if (rng() < current_density / double(tau_t::beta())) {
    if (rng(2) == 0) {
      need_flip = true;
      sl        = flip(config.seglists[origin_color]);
      dsl       = flip(config.seglists[destination_color]);
      LOG("Moving antisegment.");
    } else {
      sl  = config.seglists[origin_color];
      dsl = config.seglists[destination_color];
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
    destination_it = std::upper_bound(dsl.begin(), dsl.end(), origin_segment);
    LOG("Moving to position {}", std::distance(dsl.cbegin(), destination_it));

    // ------------  Trace ratio  -------------

    double ln_trace_ratio =
       (need_flip ? -1 : 1) * (wdata.mu(destination_color) - wdata.mu(origin_color)) * double(origin_segment.length());

    if (wdata.has_Dt) {
      tau_t tau_c, tau_cdag;
      if (need_flip) {
        tau_c    = origin_segment.tau_cdag;
        tau_cdag = origin_segment.tau_c;
      } else {
        tau_c    = origin_segment.tau_c;
        tau_cdag = origin_segment.tau_cdag;
      }
      for (auto const &[c, slist] : itertools::enumerate(config.seglists)) {
        ln_trace_ratio += K_overlap(slist, tau_c, tau_cdag, wdata.K, destination_color, c);
        ln_trace_ratio -= K_overlap(slist, tau_c, tau_cdag, wdata.K, origin_color, c);
      }
      ln_trace_ratio -= real(wdata.K(double(origin_segment.length()))(origin_color, origin_color));
      ln_trace_ratio -= real(wdata.K(double(origin_segment.length()))(destination_color, destination_color));
      ln_trace_ratio += 2 * real(wdata.K(double(origin_segment.length()))(origin_color, destination_color));
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    // We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
    // segment the tau_c/tau_cdag belongs to.
    auto &D_orig            = wdata.dets[origin_color];
    auto &D_dest            = wdata.dets[destination_color];
    auto det_c_time_orig    = [&](long i) { return D_orig.get_y(i).first; };
    auto det_cdag_time_orig = [&](long i) { return D_orig.get_x(i).first; };
    auto det_c_time_dest    = [&](long i) { return D_dest.get_y(i).first; };
    auto det_cdag_time_dest = [&](long i) { return D_dest.get_x(i).first; };
    long det_index_c_orig = 0, det_index_cdag_orig = 0;
    long det_index_c_dest = 0, det_index_cdag_dest = 0;
    double det_ratio = 1;
    // Times are ordered in det. We insert tau_cdag as a line (first index) and tau_c as a column.
    // c and cdag are inverted if we flip
    if (not is_full_line(origin_segment)) {
      if (need_flip) {
        det_index_c_orig    = lower_bound(det_c_time_orig, D_orig.size(), origin_segment.tau_cdag);
        det_index_cdag_orig = lower_bound(det_cdag_time_orig, D_orig.size(), origin_segment.tau_c);
        det_index_c_dest    = lower_bound(det_c_time_dest, D_dest.size(), origin_segment.tau_cdag);
        det_index_cdag_dest = lower_bound(det_cdag_time_dest, D_dest.size(), origin_segment.tau_c);
        det_ratio           = D_dest.try_insert(det_index_cdag_dest, det_index_c_dest, {origin_segment.tau_c, 0},
                                                {origin_segment.tau_cdag, 0})
           * D_orig.try_remove(det_index_cdag_orig, det_index_c_orig);
      } else {
        det_index_c_orig    = lower_bound(det_c_time_orig, D_orig.size(), origin_segment.tau_c);
        det_index_cdag_orig = lower_bound(det_cdag_time_orig, D_orig.size(), origin_segment.tau_cdag);
        det_index_c_dest    = lower_bound(det_c_time_dest, D_dest.size(), origin_segment.tau_c);
        det_index_cdag_dest = lower_bound(det_cdag_time_dest, D_dest.size(), origin_segment.tau_cdag);
        det_ratio           = D_dest.try_insert(det_index_cdag_dest, det_index_c_dest, {origin_segment.tau_cdag, 0},
                                                {origin_segment.tau_c, 0})
           * D_orig.try_remove(det_index_cdag_orig, det_index_c_orig);
      }
    }

    // ------------  Proposition ratio ------------

    // double future_dest_density = density(dsl) + origin_segment.length();
    // double density_ratio       = (double(tau_t::beta()) - future_dest_density) / (double(tau_t::beta()) - current_density);
    // if (need_flip) {
    //   future_dest_density = double(tau_t::beta()) - future_dest_density;
    //   density_ratio       = future_dest_density / current_density;
    // }
    double density_ratio = 1;

    double prop_ratio = density_ratio * int(sl.size()) / (int(dsl.size()) + 1);

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
    wdata.dets[destination_color].complete_operation();

    // Add the segment at destination
    dsl.insert(destination_it, origin_segment);

    // Remove the segment at origin
    sl.erase(sl.begin() + origin_index);

    if (need_flip) {
      config.seglists[origin_color]      = flip(sl);
      config.seglists[destination_color] = flip(dsl);
    } else {
      config.seglists[origin_color]      = sl;
      config.seglists[destination_color] = dsl;
    }

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
    wdata.dets[destination_color].reject_last_try();
  }
}; // namespace moves
