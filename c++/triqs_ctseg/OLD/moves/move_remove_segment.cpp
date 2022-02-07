#include "./move_remove_segment.hpp"

namespace triqs_ctseg {
move_remove_segment::move_remove_segment(
    const qmc_parameters *params_, configuration *config_,
    triqs::mc_tools::random_generator &RND_)
    : params(params_), config(config_), RND(RND_){};

double move_remove_segment::attempt() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "\n =================== ATTEMPT REMOVE ================ \n"
              << std::endl;
#endif
  color = RND(params->n_color); // pick up color
  if (config->hyb_dets.seg_number(color) == 0)
    return 0.0; // no segment to remove
  int op_index =
      RND(config->hyb_dets.seg_number(color) * 2); // pick up an operator
  auto seg = config->hyb_dets.find(color, op_index);

  if (config->ops_map().cyclic_right(seg.l) != seg.r)
    return 0.0; // forbid removal if there are ops of the same color in between

  double trace_ratio = std::exp(config->trace.try_remove_segment(seg));

#ifdef CTSEG_DEBUG
  // config->hyb_dets.check_mat_inv();
#endif

  double det_ratio = config->hyb_dets.try_remove(seg.l, seg.r);
  auto space = config->ops_map().space_around(seg.l, seg.r);
  int n_ops_before = 2 * config->hyb_dets.seg_number(color);
  double prop_ratio =
      (config->ops_map().seg_number(color) == 1)
          ? 2.0 / (params->beta * params->beta)
          : 2.0 * n_ops_before /
                (params->beta * space); // compute proposal probability ratio

#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    std::cerr << "Trying to remove operator #" << op_index << " on line "
              << color << std::endl;
    std::cerr << "Seg# " << op_index << " is " << seg.l->tau << ", "
              << seg.r->tau << std::endl;
    std::cerr << "***RATIOS: TR= " << trace_ratio << "|DET= " << det_ratio
              << "|PROP = " << prop_ratio << std::endl;
  }
#endif

  double s = (det_ratio > 0.0 ? 1.0 : -1.0);
  double prod = trace_ratio * det_ratio * prop_ratio;
  return (std::isfinite(prod) ? prod : s);
}

//--------------------------------------------------

double move_remove_segment::accept() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> ACCEPT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;
  config->hyb_dets.complete_remove();
  double sign_ratio = config->trace.complete_remove_segment().first;
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    config->print_to_h5();
  std::cerr << "config: " << *config << std::endl;

  config->trace.check_overlap_matrix_from_scratch();
#endif
  return sign_ratio;
}

//--------------------------------------------------
void move_remove_segment::reject() {
#ifdef CTSEG_DEBUG
  if (config->print_condition())
    std::cout << "- - - - - ====> REJECT - - - - - - - - - - -" << std::endl;
#endif
  config->id++;

  config->hyb_dets.reject_remove();
  config->trace.reject_remove_segment();
#ifdef CTSEG_DEBUG
  if (config->print_condition()) {
    // std::cerr << "config: " << *config << std::endl;
    std::cerr << "config: " << *config << std::endl;
    config->print_to_h5();
  }
  config->trace.check_overlap_matrix_from_scratch();
#endif
}
} // namespace triqs_ctseg
