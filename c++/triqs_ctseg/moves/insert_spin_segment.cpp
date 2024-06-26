#include "insert_spin_segment.hpp"
#include "../logs.hpp"
#include <cmath>

namespace triqs_ctseg::moves {

  insert_spin_segment::insert_spin_segment(work_data_t &data_, configuration_t &config_,
                                           triqs::mc_tools::random_generator &rng_)
     : wdata(data_), config(config_), rng(rng_) {
    ALWAYS_EXPECTS(config.n_color() == 2, "spin add/remove move only implemented for n_color == 2, got {}",
                   config.n_color());
  }

  // --------------------------------------------------

  double insert_spin_segment::attempt() {

    ALWAYS_EXPECTS((config.n_color() == 2),
                   "Insert spin segment only implemented for n_color = 2, but here n_color = {}", config.n_color());

    LOG("\n =================== ATTEMPT INSERT SPIN ================ \n");

    // ------------ Choice of segments --------------

    // Select spin insertion color
    orig_color = rng(config.n_color());
    auto &sl   = config.seglists[orig_color];
    LOG("Inserting spin: splitting at color {}", orig_color);

    if (sl.empty()) {
      LOG("ABORT : Empty line, cannot insert spin.");
      return 0;
    }

    dest_color = (orig_color == 0) ? 1 : 0;
    auto &dsl  = config.seglists[dest_color];

    // Randomly choose one existing segment
    prop_seg_idx = rng(sl.size());
    prop_seg     = sl[prop_seg_idx];

    // Is it a full line ?
    splitting_full_line = is_full_line(prop_seg);
    if (splitting_full_line) LOG("Inserting spin from full line.");

    auto dt1 = tau_t::random(rng, prop_seg.length());
    auto dt2 = tau_t::random(rng, prop_seg.length());

    if (dt1 == dt2) { // Almost never, but protect
      LOG("ABORT : Generated equal times");
      return 0;
    }

    // If splitting a full line, the order of tau_left and tau_right is not fixed
    if (dt1 > dt2 and not splitting_full_line) std::swap(dt1, dt2);

    tau_left  = prop_seg.tau_c - dt1; // dt1 < dt2
    tau_right = prop_seg.tau_c - dt2;
    spin_seg  = segment_t{tau_left, tau_right, true, true};
    LOG("Inserting spins in segment at position {}, at tau_left = {}, tau_right = {}", prop_seg_idx, tau_left,
        tau_right);

    if (not is_insertable_into(spin_seg, dsl)) {
      LOG("ABORT : Space is occupied on other line.");
      return 0;
    }

    // ------------  Trace ratio  -------------

    double ln_trace_ratio = (wdata.mu(dest_color) - wdata.mu(orig_color)) * spin_seg.length();
    if (wdata.has_Dt) {
      for (auto [c, slist] : itertools::enumerate(config.seglists)) {
        // "antisegment" - careful with order
        ln_trace_ratio += K_overlap(slist, spin_seg.tau_cdag, spin_seg.tau_c, wdata.K, orig_color, c);
        ln_trace_ratio += K_overlap(slist, spin_seg.tau_c, spin_seg.tau_cdag, wdata.K, dest_color, c);
      }
      // Add interactions of the inserted operators with themselves
      auto Kl = wdata.K(double(spin_seg.length())); // matrix
      ln_trace_ratio -= real(Kl(orig_color, orig_color));
      ln_trace_ratio -= real(Kl(dest_color, dest_color));
      ln_trace_ratio += 2 * real(Kl(orig_color, dest_color));
    }
    double trace_ratio = std::exp(ln_trace_ratio);
    trace_ratio *= -(real(wdata.Jperp(double(spin_seg.length()))(0, 0)) / 2);

    // ------------  Det ratio  ---------------

    double det_ratio = 1;

    // ------------  Proposition ratio ------------

    // T direct -> 1/ length(prop_seg) ^2  1/ ncolor 1/ sl.size  * (splitting_full_line ?  1 : 2)
    //                                   # The last because of the swap we have 2 ways to get to the same (dt1, dt2) except in full line case
    // T inverse -> 1/ Jperlist size in new config * (splitting_full_line ? 1/2 : 1)  # see reverse move
    // the (splitting_full_line ...) simplify to a ratio of 2
    double prop_ratio = (double(config.n_color()) * sl.size() * prop_seg.length() * prop_seg.length() / 2)
       / double(config.Jperp_list.size() + 1);

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double insert_spin_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    auto &sl  = config.seglists[orig_color];
    auto &dsl = config.seglists[dest_color];

    // Split segment at origin
    if (splitting_full_line) {
      sl[prop_seg_idx] = segment_t{tau_right, tau_left, true, true};
    } else {
      auto new_seg_left  = segment_t{prop_seg.tau_c, tau_left, prop_seg.J_c, true};
      auto new_seg_right = segment_t{tau_right, prop_seg.tau_cdag, true, prop_seg.J_cdag};
      // Update the proposed segment
      sl[prop_seg_idx] = new_seg_left;
      // Insert a new segment
      // cf split_segment. the new_seg_right need to be insert at front if it is part of a cyclic segment and is not cyclic itself
      bool insert_at_front = is_cyclic(prop_seg) and not is_cyclic(new_seg_right);
      sl.insert(begin(sl) + (insert_at_front ? 0 : prop_seg_idx + 1), new_seg_right);
    }

    // Insert segment at destination
    dsl.insert(std::upper_bound(begin(dsl), end(dsl), spin_seg), spin_seg);

    // Insert Jperp line
    auto &jl = config.Jperp_list;
    if (dest_color == 0)
      jl.push_back(jperp_line_t{spin_seg.tau_c, spin_seg.tau_cdag});
    else
      jl.push_back(jperp_line_t{spin_seg.tau_cdag, spin_seg.tau_c});

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata);
    LOG("Configuration is {}", config);

    return 1.0;
  }

  //--------------------------------------------------
  void insert_spin_segment::reject() { LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n"); }

} // namespace triqs_ctseg::moves
