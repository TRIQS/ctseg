#pragma once
#include "./types.hpp"
#include <boost/lexical_cast.hpp>

template <typename T> std::string str(T a) {
  return boost::lexical_cast<std::string>(a);
}
inline bool exists(std::string filename) {
  std::ifstream f(filename.c_str());
  bool e = f.good();
  if (e)
    f.close();
  return e;
}

namespace triqs::gfs {
bool is_zero(gf_const_view<imfreq> g, double epsilon = 1e-13);
bool is_zero(block_gf_const_view<imfreq> g, double epsilon = 1e-13);
} // namespace triqs::gfs

namespace triqs_ctseg {
/// returns the block number and the index within the block of a given color
std::pair<int, int>
color_to_block_and_inner_index_impl(int const color_number,
                                    gf_struct_t const &gf_struct);

/// from a block index and inner_index within this block, return color index
int block_and_inner_index_to_color_impl(int block, int inner_index,
                                        gf_struct_t const &gf_struct);
template <typename V> int get_index(std::vector<V> const &v, V const &y) {
  int i = 0;
  for (auto const &x : v) {
    if (x == y)
      return i;
    i++;
  }
  TRIQS_RUNTIME_ERROR << "[get_index] no element " << y
                      << " in this vector! (v[0] = " << v[0] << "...)";
}
std::tuple<int, std::vector<std::string>, std::vector<int>>
decode_gf_struct(const gf_struct_t &gf_struct);
} // namespace triqs_ctseg
