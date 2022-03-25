#include "insert_segment.hpp"
#include "../logs.hpp"
#include <cmath>

namespace moves {

  double insert_segment::attempt() {

    LOG("\n =================== ATTEMPT INSERT ================ \n");

    // ------------ Choice of segment --------------

    // Select insertion color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    LOG("Inserting at color {}", color);

    // Select insertion window [tau1,tau2] (defaults to [beta,0])
    qmc_time_t tau1      = wdata.qmc_beta;
    qmc_time_t tau2      = wdata.qmc_zero;
    bool config_is_empty = sl.empty();
    if (not config_is_empty) {
      // Randomly choose one existing segment
      long ind_segment     = rng(sl.size());
      tau1                 = sl[ind_segment].tau_cdag; // tau1 is cdag of this segment
      bool is_last_segment = ind_segment == sl.size() - 1;
      tau2 = sl[is_last_segment ? 0 : ind_segment + 1].tau_c; // tau2 is c of next segment, possibly cyclic
      if (tau2 == wdata.qmc_beta and tau1 == wdata.qmc_zero) {
        LOG("Full line, cannot insert.");
        return 0; // If segment is a full line, cannot insert
      }
    }

    // Choose new segment within insertion window
    qmc_time_t l = config_is_empty ? wdata.qmc_beta : tau1 - tau2;
    LOG("Insertion window size is {}", l);
    auto dt1 = fac.get_random_pt(rng, wdata.qmc_zero, l);
    auto dt2 = fac.get_random_pt(rng, wdata.qmc_zero, l);
    if (dt1 == dt2 or dt1 == wdata.qmc_zero or dt2 == wdata.qmc_zero or dt1 == l or dt2 == l) {
      LOG("Generated time at boundary");
      return 0;
    }
    if (dt1 > dt2 and not config_is_empty) std::swap(dt1, dt2); // if inserting into an empty line, two ways to insert
    proposed_segment             = segment_t{tau1 - dt1, tau1 - dt2};
    auto proposed_segment_length = double(proposed_segment.tau_c - proposed_segment.tau_cdag); // can be cyclic.
    // The index of the segment if it is inserted in the list of segments.
    proposed_segment_insert_it = std::upper_bound(sl.begin(), sl.end(), proposed_segment);

    LOG("Inserting c at {}, cdag at {}", double(proposed_segment.tau_c), double(proposed_segment.tau_cdag));

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(color, c) * overlap(config.seglists[c], proposed_segment, fac);
        ln_trace_ratio += wdata.mu(c) * proposed_segment_length;
        if (wdata.has_Dt)
          ln_trace_ratio += K_overlap(config.seglists[c], proposed_segment, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    // pos is the position of the proposed segment if inserted, converted from iterator to int
    long pos = std::distance(sl.begin(), proposed_segment_insert_it);
    // We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
    // segment the tau_c/tau_cdag belongs to.
    auto det_ratio =
       wdata.dets[color].try_insert(pos, pos, {proposed_segment.tau_cdag, 0}, {proposed_segment.tau_c, 0});

    // ------------  Proposition ratio ------------

    double current_number_intervals = std::max(1.0, double(sl.size()));
    double future_number_segments   = double(sl.size()) + 1;
    double prop_ratio = (current_number_intervals * l * l / (config_is_empty ? 1 : 2)) / future_number_segments;
    // Account for absence of time swapping when inserting into empty line.

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double insert_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Insert the times into the det
    wdata.dets[color].complete_operation();

    // Compute the sign ratio
    double sign_ratio = 1;
    auto &sl          = config.seglists[color];
    if (is_cyclic(proposed_segment)) sign_ratio = -1.0;
    LOG("Sign ratio is {}", sign_ratio);

    // Insert the segment in an ordered list
    sl.insert(proposed_segment_insert_it, proposed_segment);

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
  void insert_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
}; // namespace moves
