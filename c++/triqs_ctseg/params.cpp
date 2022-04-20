#include "params.hpp"

void h5_write(h5::group h5group, std::string subgroup_name, constr_params_t const &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

  h5_write(grp, "beta", c.beta);
  h5_write(grp, "gf_struct", c.gf_struct);
  h5_write(grp, "n_tau", c.n_tau);
  h5_write(grp, "n_tau_k", c.n_tau_k);
}

//------------------------------------

void h5_read(h5::group h5group, std::string subgroup_name, constr_params_t &c) {

  h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

  h5_read(grp, "beta", c.beta);
  h5_read(grp, "gf_struct", c.gf_struct);
  h5_read(grp, "n_tau", c.n_tau);
  h5_read(grp, "n_tau_k", c.n_tau_k);
}
