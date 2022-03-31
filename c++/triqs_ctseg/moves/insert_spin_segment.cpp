#include "insert_spin_segment.hpp"
#include "../logs.hpp"
#include <cmath>

namespace moves {

  double insert_spin_segment::attempt() {

    ALWAYS_EXPECTS((wdata.n_color == 2), "Insert spin segment only implemented for n_color = 2, but here n_color = {}",
                   wdata.n_color);

    LOG("\n =================== ATTEMPT INSERT SPIN ================ \n");

    // ------------ Choice of segments --------------

    // Select spin insertion color
    orig_color = rng(wdata.n_color);
    auto &sl   = config.seglists[orig_color];
    LOG("Inserting spin: splitting at color {}", orig_color);

    if (sl.empty()) {
      LOG("ABORT : Empty line, cannot insert spin.");
      return 0;
    }

    // FIXME if n_color > 2
    dest_color = orig_color == 0 ? 1 : 0;
    auto &dsl  = config.seglists[dest_color];

    // Randomly choose one existing segment
    prop_seg_idx = rng(sl.size());
    prop_seg     = sl[prop_seg_idx];

    // Is it a full line ?
    splitting_full_line = is_full_line(prop_seg);
    if (splitting_full_line) LOG("Inserting spin from full line.");

    // Choose spin insertion points (tau_left,tau_right)
    // FIXME : S?IMPLy length(prop_seg) Should be a dimtime_t ...
    // IT does cover full line ...
    //dimtime_t prop_seg_length = prop_seg.tau_c - prop_seg.tau_cdag;

    dimtime_t prop_seg_length = splitting_full_line ? wdata.qmc_beta : prop_seg.tau_c - prop_seg.tau_cdag;
    dimtime_t dt1             = dimtime_t::random(rng, prop_seg_length);
    dimtime_t dt2             = dimtime_t::random(rng, prop_seg_length);

    if (dt1 == dt2) { // Almost never, but protect
      LOG("ABORT : Generated equal times");
      return 0;
    }

    // If splitting a full line, the order of tau_left and tau_right is not fixed
    if (dt1 > dt2 and !splitting_full_line) std::swap(dt1, dt2);

    // FIXME : bug? why tau_c ?  tau_c_dag ? why on the left ??
    tau_left  = prop_seg.tau_c - dt1; // dt1 < dt2
    tau_right = prop_seg.tau_c - dt2;
    spin_seg  = segment_t{tau_left, tau_right, true, true};
    LOG("Inserting spins at tau_left = {}, tau_right = {}", tau_left, tau_right);

    if (not is_insertable(dsl, spin_seg)) {
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
        if (splitting_full_line)
          ln_trace_ratio -= K_overlap(slist, wdata.qmc_beta, wdata.qmc_zero, wdata.K, orig_color, c);
      }
      // Add interactions of the inserted operators with themselves
      // Why real ?? K is not real ??
      ln_trace_ratio -= real(wdata.K(spin_seg.length())(orig_color, orig_color));
      ln_trace_ratio -= real(wdata.K(spin_seg.length())(dest_color, dest_color));
    }
    double trace_ratio = std::exp(ln_trace_ratio);
    trace_ratio *= -real(wdata.Jperp(double(spin_seg.tau_c - spin_seg.tau_cdag))(0, 0)) / 2;

    // ------------  Det ratio  ---------------

    double det_ratio = 1;

    // ------------  Proposition ratio ------------

    double prop_ratio = (double(wdata.n_color) * sl.size() * prop_seg_length * prop_seg_length / 2)
       / double(config.Jperp_list.size() + 1);
    if (splitting_full_line) prop_ratio *= 4; // Account for absence of time swapping when splitting full line,
    // and two possible choices for the recombination in reverse move

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
      // FIXME simply
      // auto new_segment   = segment_t{tau_right, tau_left, true, true};
      auto new_segment   = segment_t{tau_right, tau_left};
      new_segment.J_c    = true;
      new_segment.J_cdag = true;
      sl[prop_seg_idx]   = new_segment;
    } else {
      //auto new_seg_left      = segment_t{prop_seg.tau_c, tau_left, false, true}; // J_c, J_cdag
      //auto new_seg_right     = segment_t{tau_right, prop_seg.tau_cdag, true, false};
      auto new_seg_left   = segment_t{prop_seg.tau_c, tau_left};
      new_seg_left.J_cdag = true;
      auto new_seg_right  = segment_t{tau_right, prop_seg.tau_cdag};
      new_seg_right.J_c   = true;

      bool segment_overboard = is_cyclic(prop_seg) and !is_cyclic(new_seg_right);
      // Index of the rightmost of the two produced segments (the one to the left is always at prop_seg_idx)
      // FIXME splitting_full_line is false here ...
      int right_seg_idx = (splitting_full_line or segment_overboard) ? 0 : prop_seg_idx + 1;
      // Update the proposed segment
      sl[prop_seg_idx] = new_seg_left;
      // Insert a new segment
      sl.insert(sl.begin() + right_seg_idx, new_seg_right);
    }

    // Insert segment at destination
    auto spin_seg_it = std::upper_bound(dsl.begin(), dsl.end(), spin_seg);
    //auto spin_seg_it = std::upper_bound(dsl.begin(), dsl.end(), spin_seg);
    spin_seg.J_c    = true; // FIXME : move in decl of spin_seg
    spin_seg.J_cdag = true;
    dsl.insert(spin_seg_it, spin_seg);

    // Insert Jperp line
    auto &jl = config.Jperp_list;
    jperp_line_t jline;
    if (dest_color == 0)
      jline = jperp_line_t{spin_seg.tau_c, spin_seg.tau_cdag};
    else
      jline = jperp_line_t{spin_seg.tau_cdag, spin_seg.tau_c};

    // FIXME : push_back ?
    jl.insert(jl.end(), jline);

    // Check invariant
    // FIXME : simply
    //if (dest_color == 0)
    //jl.push_back(jperp_line_t{spin_seg.tau_c, spin_seg.tau_cdag});
    //else
    //jl.push_back(jperp_line_t{spin_seg.tau_cdag, spin_seg.tau_c});

    // Add the check of J_c flags in the invariants
#ifdef CHECK_INVARIANTS
    check_invariant(config, wdata.dets);
#endif
    SPDLOG_TRACE("Configuration is {}", config);

    return 1.0;
  }

  //--------------------------------------------------
  void insert_spin_segment::reject() { LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n"); }
}; // namespace moves
