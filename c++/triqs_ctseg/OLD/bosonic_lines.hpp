#pragma once
#include "./composite_segment.hpp"
#include <boost/mpi/communicator.hpp>

namespace triqs_ctseg {
using namespace triqs::gfs;
using namespace triqs::mesh;

class bosonic_lines {
  std::vector<composite_segment> bosonic_links;
  gf<imtime> Jperp;
  bool keep_negative;

public:
  bosonic_lines(gf<imtime> const &Jperp_, bool keep_negative_)
      : bosonic_links(), Jperp(Jperp_), keep_negative(keep_negative_) {
    int positive_counter = 0;
    for (auto const &t : Jperp.mesh())
      if (Jperp[t](0, 0).real() > 0) {
        // std::stringstream fs; fs << "WARNING: ctseg : Jperp(tau) is positive
        // for some value of tau (causality violation?) : Jperp("<<t<< ")= "<<
        // Jperp[t](0,0);
        positive_counter++;
        // throw std::invalid_argument(fs.str());
        // std::cerr << fs.str() << std::endl;
        if (keep_negative)
          Jperp[t](0, 0) = -1e-7;
      }
    mpi::communicator c;
    if (positive_counter && keep_negative && c.rank() == 0)
      std::cerr << " WARNING: " << positive_counter
                << " points in Jperp_tau are positive and hence manually set "
                   "to -1e-7 "
                << std::endl;
    if (positive_counter && !keep_negative && c.rank() == 0)
      std::cerr << " WARNING: " << positive_counter
                << " points in Jperp_tau are positive " << std::endl;
  };

  size_t size() const { return bosonic_links.size(); }
  double jperp(composite_segment const &seg) const {
    return Jperp[closest_mesh_pt(seg.plus.a->tau - seg.minus.c->tau)](0, 0)
        .real();
  }
  double jperp(qmc_time_t const &tau) const {
    return Jperp[closest_mesh_pt(tau)](0, 0).real();
  }
  double jperp(double const &tau) const {
    return Jperp[closest_mesh_pt(double(tau))](0, 0).real();
  }

  double add(composite_segment const &cseg) {
    bosonic_links.push_back(cseg);
    return -jperp(cseg) / 2;
  }

  std::pair<double, composite_segment> remove(int n) {
    if (n >= size())
      TRIQS_RUNTIME_ERROR << " out of bounds";
    auto cseg = bosonic_links[n];
    bosonic_links.erase(bosonic_links.begin() + n);
    return {-2.0 / jperp(cseg), cseg};
  }
  /// returns true if tau is contained in bosonic links
  bool contains(qmc_time_t tau) {
    for (auto it = bosonic_links.begin(); it != bosonic_links.end(); it++)
      if (it->plus.c->tau == tau || it->plus.a->tau == tau ||
          it->minus.c->tau == tau || it->minus.a->tau == tau)
        return true;
    return false;
  }
  void pop_back() { bosonic_links.pop_back(); }

  composite_segment const &operator[](int i) const { return bosonic_links[i]; }

  // print
  friend std::ostream &operator<<(std::ostream &out, bosonic_lines const &b) {
    for (auto it = b.bosonic_links.begin(); it != b.bosonic_links.end(); it++)
      out << "(" << double(it->plus.c->tau) << ", " << double(it->minus.a->tau)
          << ") - ";
    return out;
  }
  friend void h5_write(h5::group fg, std::string subgroup_name,
                       bosonic_lines const &b) {
    h5::group gr = fg.create_group(subgroup_name);
    array<double, 2> a(b.size(), 2);
    int i = 0;
    for (auto it = b.bosonic_links.begin(); it != b.bosonic_links.end();
         it++, i++) {
      a(i, 0) = double(it->plus.c->tau);
      a(i, 1) = double(it->minus.a->tau);
    }
    if (b.size() > 0)
      h5_write(gr, "times_table", a);
  }
};
} // namespace triqs_ctseg
