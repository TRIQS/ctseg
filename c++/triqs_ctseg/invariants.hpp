#pragma once
#include <vector>
#include "configuration.hpp"
#include "work_data.hpp"

namespace triqs_ctseg {

  void check_invariant(configuration_t const &config, work_data_t const &wdata);

  void check_segments(configuration_t const &config);

  void check_dets(configuration_t const &config, work_data_t const &wdata);

  void check_jlines(configuration_t const &config);

} // namespace triqs_ctseg
