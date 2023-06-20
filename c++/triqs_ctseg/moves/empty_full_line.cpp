#include "empty_full_line.hpp"
#include "../logs.hpp"
#include <cmath>

namespace moves {

  double empty_full_line::attempt() {

    LOG("\n =================== ATTEMPT EMPTY FULL LINE ================ \n");

    // Select color
    color    = rng(config.n_color());
    auto &sl = config.seglists[color];
    LOG("Emptying at color {}", color);

    if (sl.empty()) {
      LOG("Color is empty, nothing to remove.");
      return 0;
    }

    if (not is_full_line(sl.back())) {
      LOG("Color has no full line, cannot empty.");
      return 0;
    }

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = -wdata.mu(color) * tau_t::beta(); // chemical potential
    // Now the overlaps
    for (auto c : range(config.n_color())) {
      if (c != color) ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], segment_t::full_line());
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------
    auto det_ratio = 1.0;

    // ------------  Proposition ratio ------------
    double prop_ratio = 1.0;

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;

    return (std::isfinite(prod) ? prod : 1.0);
  }

  //--------------------------------------------------

  double empty_full_line::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Remove the full line
    auto &sl = config.seglists[color];
    sl.erase(sl.begin());

    // Check invariant
    if constexpr (print_logs or ctseg_debug) check_invariant(config, wdata.dets);

    LOG("Configuration is {}", config);

    return 1.0;
  }

  //--------------------------------------------------
  void empty_full_line::reject() { LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n"); }
}; // namespace moves
