#include "insert_segment.hpp"
//#include <triqs/gfs/functions/functions2.hpp>
namespace moves {

  double insert_segment::attempt() {

    SPDLOG_LOGGER_TRACE("\n =================== ATTEMPT INSERT ================ \n", void);

    auto qmc_beta = time_point_factory.get_upper_pt();
    auto qmc_zero = time_point_factory.get_lower_pt();

    // ------------ Choice of segment --------------

    // Select insertion color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    SPDLOG_LOGGER_TRACE("Inserting at color {}", color);

    // Select insertion window [tau1,tau2] (defaults to [beta,0])
    qmc_time_t tau1      = qmc_beta;
    qmc_time_t tau2      = qmc_zero;
    bool config_is_empty = sl.empty();
    if (not config_is_empty) {
      // Randomly choose one existing segment
      long ind_segment     = rng(sl.size());
      tau1                 = sl[ind_segment].tau_cdag; // tau1 is cdag of this segment
      bool is_last_segment = ind_segment == sl.size() - 1;
      tau2 = sl[is_last_segment ? 0 : ind_segment + 1].tau_c; // tau2 is c of next segment, possibly cyclic
      if (tau2 == qmc_beta and tau1 == qmc_zero) return 0;    // If segment is a full line, cannot insert
    }

    // Choose new segment within insertion window
    qmc_time_t l = tau1 - tau2;
    auto dt1     = time_point_factory.get_random_pt(rng, qmc_zero, l);
    auto dt2     = time_point_factory.get_random_pt(rng, qmc_zero, l);
    if (dt1 == dt2) return 0;
    if (dt1 > dt2 and not config_is_empty) std::swap(dt1, dt2); // if inserting into an empty line, two ways to insert
    proposed_segment             = segment_t{tau1 - dt1, tau1 - dt2};
    auto proposed_segment_length = double(proposed_segment.tau_c - proposed_segment.tau_cdag); // can be cyclic.
    // The index of the segment if it is inserted in the list of segments.
    proposed_segment_insert_it = std::upper_bound(sl.begin(), sl.end(), proposed_segment);

    SPDLOG_LOGGER_TRACE("Inserting c at{}, cdag at {}", proposed_segment.tau_c, proposed_segment.tau_cdag);

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(color, c) * overlap(config.seglists[c], proposed_segment, time_point_factory);
        ln_trace_ratio += wdata.mu(c) * proposed_segment_length;
        if (wdata.has_Dt)
          ln_trace_ratio += K_overlap(config.seglists[c], proposed_segment, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    // pos is the position of the proposed segment if inserted, converted from iterator to int
    long pos = std::distance(proposed_segment_insert_it, sl.begin());
    // We insert tau_cdag as a line (first index) and tau_c as a column (second index). The index always corresponds to the
    // segment the tau_c/tau_cdag belongs to.
    auto det_ratio =
       wdata.dets[color].try_insert(pos, pos, {proposed_segment.tau_cdag, 0}, {proposed_segment.tau_c, 0});

    // ------------  Proposition ratio ------------

    double current_number_intervals =
       config_is_empty ? 2 : double(sl.size()); // Account for absence of time swapping when inserting into empty line.
    double future_number_segments = double(sl.size()) + 1;
    double prop_ratio             = future_number_segments / (current_number_intervals * l * l / 2);

    SPDLOG_LOGGER_TRACE("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    return trace_ratio * det_ratio * prop_ratio;
  }

  //--------------------------------------------------

  double insert_segment::accept() {

    SPDLOG_LOGGER_TRACE("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n", void);

    // Insert the times into the det
    wdata.dets[color].complete_operation();
    // Insert the segment in an ordered list
    config.seglists[color].insert(proposed_segment_insert_it, proposed_segment);

    // FIXME ??? SIGNE ???
    double sign_ratio = 1; // ???config->trace.complete_insert_segment();

    // SPDLOG_LOGGER_TRACE("Configuration {}", config);
    // Check invariant ??
    // config->trace.check_overlap_matrix_from_scratch();

    return sign_ratio;
  }

  //--------------------------------------------------
  void insert_segment::reject() {
    SPDLOG_LOGGER_TRACE("\n - - - - - ====> REJECT - - - - - - - - - - -\n", void);
    wdata.dets[color].reject_last_try();

    // SPDLOG_LOGGER_TRACE("Configuration {}", config);
    // Check invariant ??
    // config->trace.check_overlap_matrix_from_scratch();
  }
}; // namespace moves
