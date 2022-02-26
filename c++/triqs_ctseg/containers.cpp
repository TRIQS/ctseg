#include "containers.hpp"

void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

  h5_write(grp, "G_tau", c.G_tau);
  h5_write(grp, "chi_tau", c.chi_tau);
  h5_write(grp, "densities", c.densities);
}

//------------------------------------

void h5_read(h5::group h5group, std::string subgroup_name, results_t &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

  h5_read(grp, "G_tau", c.G_tau);
  h5_read(grp, "chi_tau", c.chi_tau);
  h5_read(grp, "densities", c.densities);
}
