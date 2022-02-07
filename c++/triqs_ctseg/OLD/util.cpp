#include "./util.hpp"

namespace triqs::gfs {
bool is_zero(gf_const_view<imfreq> g, double epsilon) {
  for (auto const &x : g.data()) {
    if (std::abs(x) > epsilon)
      return false;
  }
  return true;
}
bool is_zero(block_gf_const_view<imfreq> G, double epsilon) {
  for (auto const &g : G)
    if (!is_zero(g, epsilon))
      return false;
  return true;
}
} // namespace triqs::gfs

namespace triqs_ctseg {
std::pair<int, int>
color_to_block_and_inner_index_impl(int const color_number,
                                    gf_struct_t const &gf_struct) {
  int bl = 0;
  int counter = 0;
  for (auto const &[bl_name, bl_size] : gf_struct) {
    for (auto idx : range(bl_size)) {
      if (counter == color_number)
        return {bl, idx};
      counter++;
    }
    bl++;
  }
  TRIQS_RUNTIME_ERROR << "color " << color_number << " is incorrect";
}

/// from a block index and inner_index within this block, return color index
int block_and_inner_index_to_color_impl(int block, int inner_index,
                                        gf_struct_t const &gf_struct) {
  int color = 0, bl = 0;
  for (auto const &[bl_name, bl_size] : gf_struct) {
    if (bl == block)
      for (auto i : range(bl_size)) {
        if (i == inner_index)
          return color;
        color++;
      }
    else
      color += bl_size;
    bl++;
  }
  // if (color>=m.n_colors()) TRIQS_RUNTIME_ERROR <<"Wrong block,inner_index
  // pair: " << block << ","<< inner_index ;
  return color;
}

std::tuple<int, std::vector<std::string>, std::vector<int>>
decode_gf_struct(const gf_struct_t &gf_struct) {
  int n_blocks = gf_struct.size();
  std::vector<std::string> block_names = std::vector<std::string>{};
  std::vector<int> block_sizes = std::vector<int>{};
  for (auto [bl_name, bl_size] : gf_struct) {
    block_names.push_back(bl_name);
    block_sizes.push_back(bl_size);
  }
  return std::make_tuple(n_blocks, block_names, block_sizes);
}
} // namespace triqs_ctseg
