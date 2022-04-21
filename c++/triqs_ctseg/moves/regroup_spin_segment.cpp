#include "regroup_spin_segment.hpp"
#include "../logs.hpp"
#include "triqs_ctseg/configuration.hpp"
#include <cmath>

namespace moves {

  double regroup_spin_segment::attempt() {

    LOG("\n =================== ATTEMPT REGROUP SPIN ================ \n");

    ln_trace_ratio = 0;
    prop_ratio     = 1 / (double(config.Jperp_list.size()) + 1);

    // ----------- Propose move in each color ----------
    bool prop_failed;
    std::tie(idx_c_up, idx_cdag_dn, tau_up, prop_failed) = propose(0); // spin up
    if (prop_failed) return 0;
    std::tie(idx_c_dn, idx_cdag_up, tau_dn, prop_failed) = propose(1); // spin down
    if (prop_failed) return 0;

    // ----------- Trace ratio -----------
    // Correct for the overlap between the two modified segments
    auto &sl_up     = config.seglists[0];
    auto &sl_dn     = config.seglists[1];
    auto old_seg_up = sl_up[idx_c_up];
    auto old_seg_dn = sl_dn[idx_c_dn];
    auto new_seg_up = segment_t{tau_up, old_seg_up.tau_cdag};
    auto new_seg_dn = segment_t{tau_dn, old_seg_dn.tau_cdag};

    ln_trace_ratio += -wdata.U(0, 1)
       * (overlap(new_seg_up, new_seg_dn) + overlap(old_seg_up, old_seg_dn) //
          - overlap(new_seg_up, old_seg_dn) - overlap(new_seg_dn, old_seg_up));

    // Correct for the dynamical interaction between the two operators that have been moved
    if (wdata.has_Dt) {
      ln_trace_ratio -= real(wdata.K(double(tau_up - old_seg_dn.tau_c))(0, 1));
      ln_trace_ratio -= real(wdata.K(double(tau_dn - old_seg_up.tau_c))(0, 1));
      ln_trace_ratio += real(wdata.K(double(tau_dn - tau_up))(0, 1));
      ln_trace_ratio += real(wdata.K(double(old_seg_up.tau_c - old_seg_dn.tau_c))(0, 1));
    }

    double trace_ratio = std::exp(ln_trace_ratio);
    trace_ratio *= -real(wdata.Jperp(double(tau_up - tau_dn))(0, 0)) / 2;

    // ----------- Det ratio -----------
    double det_ratio = 1;

    // Spin up
    auto &D_up = wdata.dets[0];
    det_ratio *= D_up.try_remove(det_lower_bound_x(D_up, sl_up[idx_cdag_up].tau_cdag),
                                 det_lower_bound_y(D_up, sl_up[idx_c_up].tau_c));

    // Spin down
    auto &D_dn = wdata.dets[1];
    det_ratio *= D_dn.try_remove(det_lower_bound_x(D_dn, sl_dn[idx_cdag_dn].tau_cdag),
                                 det_lower_bound_y(D_dn, sl_dn[idx_c_dn].tau_c));

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double regroup_spin_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = config_sign(config, wdata.dets);
    LOG("Initial sign is {}. Initial configuration: {}", initial_sign, config);

    // Update dets
    wdata.dets[0].complete_operation();
    wdata.dets[1].complete_operation();

    // Update the segments
    auto &sl_up = config.seglists[0];
    auto &sl_dn = config.seglists[1];

    // Update tau_c
    sl_up[idx_c_up].tau_c = tau_up;
    sl_dn[idx_c_dn].tau_c = tau_dn;

    // Update spin tags
    sl_up[idx_c_up].J_c       = true;
    sl_dn[idx_cdag_dn].J_cdag = true;
    sl_dn[idx_c_dn].J_c       = true;
    sl_up[idx_cdag_up].J_cdag = true;

    fix_ordering_first_last(sl_up);
    fix_ordering_first_last(sl_dn);

    // Add spin line
    config.Jperp_list.push_back(jperp_line_t{tau_up, tau_dn});

    // Compute sign
    double final_sign = config_sign(config, wdata.dets);
    double sign_ratio = final_sign / initial_sign;
    LOG("Final sign is {}", final_sign);

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata.dets);
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}. Config: {} ",
                   det_sign, sign_ratio, config);
    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------

  void regroup_spin_segment::reject() {

    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[0].reject_last_try();
    wdata.dets[1].reject_last_try();
  }

  //--------------------------------------------------

  std::tuple<long, long, tau_t, bool> regroup_spin_segment::propose(int color) {

    auto &sl        = config.seglists[color];
    int other_color = 1 - color;
    auto &dsl       = config.seglists[other_color];

    // --------- Eliminate cases where move is impossible ---------
    if (sl.empty() or dsl.empty()) {
      LOG("Line is empty, cannot regroup spin.");
      return {0, 0, tau_t::zero(), true};
    }
    if (is_full_line(sl[0]) or is_full_line(dsl[0])) {
      LOG("Line is full, cannot regroup spin.");
      return {0, 0, tau_t::zero(), true};
    }

    // --------- Randomly choose a c operator -----------

    auto idx_c = rng(sl.size());
    LOG("Spin {}: regrouping c at position {}.", (color == 0) ? "up" : "down", idx_c);
    if (sl[idx_c].J_c) {
      LOG("Spin {}: cannot regroup because c is connected to a spin line.", (color == 0) ? "up" : "down");
      return {0, 0, tau_t::zero(), true};
    }

    auto tau_c = sl[idx_c].tau_c;

    // -------- Propose new position for the c ---------

    // Determine window in which the c can be moved
    auto idx_left   = (idx_c == 0) ? sl.size() - 1 : idx_c - 1;
    auto wtau_left  = sl[idx_left].tau_cdag;
    auto wtau_right = sl[idx_c].tau_cdag;
    if (idx_c == idx_left) {
      wtau_left  = tau_t::beta();
      wtau_right = tau_t::zero();
    }
    auto window_length = wtau_left - wtau_right;

    // Find the cdag in opposite spin that are within the window
    auto cdag_list = cdag_in_window(wtau_left, wtau_right, dsl);
    if (cdag_list.empty()) {
      LOG("Spin {}: cannot regroup because there are no suitable cdag operators.", (color == 0) ? "up" : "down");
      return {0, 0, tau_t::zero(), true};
    }
    // Choose one of them randomly
    auto idx_cdag = cdag_list[rng(cdag_list.size())];
    if (dsl[idx_cdag].J_cdag) {
      LOG("Spin {}: cannot regroup because chosen cdag is connected to a spin line.", (color == 0) ? "up" : "down");
      return {0, 0, tau_t::zero(), true};
    }
    auto tau_c_new = dsl[idx_cdag].tau_cdag;
    auto new_seg   = segment_t{tau_c_new, sl[idx_c].tau_cdag};
    LOG("Spin {}: moving c from {} to {}.", (color == 0) ? "up" : "down", tau_c, tau_c_new);

    // FIXME FACTOR ?
    // -------- Trace ratio ---------
    ln_trace_ratio += wdata.mu(color) * (double(new_seg.length()) - double(sl[idx_c].length()));
    LOG("Spin {}: ln trace ratio = {}", (color == 0) ? "up" : "down", ln_trace_ratio);
    for (auto const &[c, slc] : itertools::enumerate(config.seglists)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(c, color) * overlap(slc, new_seg);
        ln_trace_ratio -= -wdata.U(c, color) * overlap(slc, sl[idx_c]);
      }
      if (wdata.has_Dt) {
        ln_trace_ratio += K_overlap(slc, tau_c_new, true, wdata.K, c, color);
        ln_trace_ratio -= K_overlap(slc, tau_c, true, wdata.K, c, color);
      }
    }
    if (wdata.has_Dt) ln_trace_ratio -= real(wdata.K(double(tau_c_new - tau_c))(color, color));

    // --------- Prop ratio ---------
    prop_ratio *= double(sl.size()) * cdag_list.size() / window_length;

    return {idx_c, idx_cdag, tau_c_new, false};
  }

}; // namespace moves
