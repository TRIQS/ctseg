#include "inputs.hpp"

namespace triqs_ctseg {

  void h5_write(h5::group h5group, std::string subgroup_name, inputs_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "delta", c.delta);
    h5_write(grp, "jperpt", c.jperpt);
    h5_write(grp, "d0t", c.d0t);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, inputs_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "delta", c.delta);
    h5_read(grp, "jperpt", c.jperpt);
    h5_read(grp, "d0t", c.d0t);
  }

} // namespace triqs_ctseg
