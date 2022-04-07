#include "split_spin_segment.hpp"
#include "../logs.hpp"
#include <cmath>

namespace moves {

  double split_spin_segment::attempt() {

    LOG("\n =================== ATTEMPT SPLIT SPIN ================ \n");

    // ------------ Choose a spin line --------------

    auto &jl = config.Jperp_list;

    if (jl.empty()) {
      LOG("No bosonic lines!");
      return 0;
    }

    line_idx   = rng(jl.size());
    auto &line = jl[line_idx];
    LOG("Splitting S_plus at {}, S_minus at {}", line.tau_Splus, line.tau_Sminus);

    ln_trace_ratio = 0;
    det_ratio      = 1;
    prop_ratio     = 1;

    // ----------- Propose move in each color ----------

    std::tie(idx_up, tau_up)     = propose(0); // spin up
    std::tie(idx_down, tau_down) = propose(1); // spin down

    // Correct for the dynamical interaction between the two operators that have been moved
    if (wdata.has_Dt) {
      auto &sl_up   = config.seglists[0];
      auto &sl_down = config.seglists[1];
      ln_trace_ratio -= real(wdata.K(double(tau_up - sl_down[idx_down].tau_c))(0, 1));
      ln_trace_ratio -= real(wdata.K(double(tau_down - sl_up[idx_up].tau_c))(0, 1));
      ln_trace_ratio += real(wdata.K(double(tau_down - tau_up))(0, 1));
    }

    double trace_ratio = std::exp(ln_trace_ratio);
    trace_ratio /= -real(wdata.Jperp(double(line.tau_Splus - line.tau_Sminus))(0, 0)) / 2;
    prop_ratio *= jl.size();

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double split_spin_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = config_sign(config, wdata.dets);
    // Update the dets
    wdata.dets[0].complete_operation();
    wdata.dets[1].complete_operation();

    // Update the segments
    auto &sl_up             = config.seglists[0];
    auto &sl_down           = config.seglists[1];
    sl_up[idx_up].tau_c     = tau_up;
    sl_up[idx_up].J_c       = false;
    sl_down[idx_down].tau_c = tau_down;
    sl_down[idx_down].J_c   = false;

    // Remove Jperp line
    auto &jl = config.Jperp_list;
    jl.erase(jl.begin() + line_idx);

    // Compute sign
    double final_sign = config_sign(config, wdata.dets);
    double sign_ratio = final_sign / initial_sign;
    LOG("Sign ratio is {}", sign_ratio);

    // Check invariant
    if constexpr (check_invariants or ctseg_debug) check_invariant(config, wdata.dets);
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}.", det_sign,
                   sign_ratio);
    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------

  void split_spin_segment::reject() {

    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[0].reject_last_try();
    wdata.dets[1].reject_last_try();
  }

  //--------------------------------------------------

  std::pair<long, tau_t> split_spin_segment::propose(int color) {

    auto &sl         = config.seglists[color];
    auto &line       = config.Jperp_list[line_idx];
    int other_color  = 1 - color;
    std::string spin = (color == 0) ? "up" : "down";

    // --------- Find the c connected to the chosen spin line ----------
    segment_t seg_c;
    if (color == 0)
      // In spin up color, the c conneted to the J line is a at tau_Sminus
      seg_c = segment_t{line.tau_Sminus, line.tau_Sminus};
    else
      // In spin down color, the c conneted to the J line is a at tau_Splus
      seg_c = segment_t{line.tau_Splus, line.tau_Splus};
    auto it = std::lower_bound(sl.cbegin(), sl.cend(), seg_c);

    // -------- Propose new position for the c ---------
    // Determine window in which the c will be moved
    auto idx            = std::distance(sl.cbegin(), it);
    auto idx_left       = (idx == 0) ? sl.size() - 1 : idx - 1;
    tau_t wtau_left     = sl[idx_left].tau_cdag;
    tau_t wtau_right    = sl[idx].tau_cdag;
    tau_t window_length = (idx_left == idx) ? tau_t::beta() : wtau_left - wtau_right;
    // Choose random time in window
    tau_t dt        = tau_t::random(rng, window_length);
    tau_t tau_c_new = sl[idx_left].tau_cdag - dt;
    tau_t tau_c     = sl[idx].tau_c;

    LOG("Spin {}: moving c at position {} from {} to {}.", spin, idx, tau_c, tau_c_new);
    auto new_seg = segment_t{tau_c_new, sl[idx].tau_cdag};

    // -------- Trace ratio ---------
    for (auto const &[c, slc] : itertools::enumerate(config.seglists)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(c, color) * overlap(slc, new_seg);
        ln_trace_ratio -= -wdata.U(c, color) * overlap(slc, sl[idx]);
      }
      ln_trace_ratio +=
         K_overlap(slc, tau_c_new, true, wdata.K, c, color) - K_overlap(slc, tau_c, true, wdata.K, c, color);
    }

    // -------- Det ratio ----------
    auto &D          = wdata.dets[color];
    auto det_c_time  = [&](long i) { return D.get_y(i).first; };
    long det_index_c = lower_bound(det_c_time, D.size(), sl[idx].tau_c);
    det_ratio *= D.try_change_col(det_index_c, {tau_c_new, 0});

    // --------- Prop ratio ---------
    auto &dsl = config.seglists[other_color];
    prop_ratio *= window_length / (double(sl.size()) * cdag_in_window(wtau_left, wtau_right, dsl).size());

    return std::pair<long, tau_t>{idx, tau_c_new};
  }
}; // namespace moves
