#pragma once
#include "../work_data.hpp"
#include "../configuration.hpp"

namespace moves {
  class insert_split_segment {
    work_data_t &wdata;
    configuration_t &config;
    triqs::mc_tools::random_generator &rng;

    // Internal data
    int color = 0;
    segment_t proposed_segment;
    std::vector<segment_t>::const_iterator seg_it;
    std::vector<segment_t>::const_iterator insert_it;
    tau_t tau1;
    tau_t tau2;
    long seg_idx;
    long insert_idx;
    long right_seg_idx;
    bool is_inside;
    bool splitting_full_line;
    bool insert_into_empty_line;
    double det_sign;

    void prep_insertion(std::vector<segment_t> const &seglist, tau_t tau);

    public:
    insert_split_segment(work_data_t &data_, configuration_t &config_, triqs::mc_tools::random_generator &rng_)
       : wdata(data_), config(config_), rng(rng_){};
    // ------------------
    double attempt();
    double accept();
    void reject();
  };
}; // namespace moves