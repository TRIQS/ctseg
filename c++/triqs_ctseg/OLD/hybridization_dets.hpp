/*******************************************************************************
 * CTSEG: TRIQS hybridization-expansion segment solver
 *
 * Copyright (C) 2013-2018 by T. Ayral, H. Hafermann, P. Delange, M. Ferrero, O.
 *Parcollet
 *
 * CTSEG is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * CTSEG is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * CTSEG. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#pragma once

#include "./hyb_opmap.hpp"
#include "./hyb_opmap_no_storage.hpp"
#include "./t_ordered_colored_c_ops.hpp"
#include <triqs/det_manip/det_manip.hpp>

namespace triqs_ctseg {

using namespace triqs::gfs;
using namespace triqs::mesh;

// use the first one in the master branch !
typedef hyb_opmap HybOpmapImpl;
// typedef hyb_opmap_no_storage HybOpmapImpl;

/**
 * This class handles the determinants of hybridisation and stores the c,cdag
 *which are hybridized
 *
 *- add : connects two operators to hybridisation lines
 *- remove : disconnects two operators
 *
 * Each operation has a try/complete/reject mode, hence the object has 3
 *possibles states. (cf status) with assert to prevent misuse.
 *
 * Part of the implementation is done in 2 implementations classes hyb_opmap and
 *hyb_opmap_no_storage The first one stores an additionnal map with the c,cdag
 *which are hybridized The second one is a simpler version, where all c, cdag
 *are always hybridized, hence there is no need to store an additional map.
 *
 */
class hybridization_dets : HybOpmapImpl {

  // This callable object adapts the hybridization function for the call of the
  // det.
  struct hybridization_adaptor {
    gf<imtime, matrix_real_valued> d;
    hybridization_adaptor(gf<imtime, matrix_real_valued> const &hybridization)
        : d(std::move(hybridization)) {
      int positive_counter = 0;
      for (auto const &t : d.mesh())
        if (d[t](0, 0) > 0) {
          std::stringstream fs;
          fs << "WARNING: ctseg : solve aborted\n Hybridization Delta(tau) is "
                "positive for some value of tau (causality violation) : Delta("
             << t << ")= " << d[t](0, 0);
          positive_counter++;
          // throw std::invalid_argument(fs.str());
          // std::cerr << fs.str() << std::endl;
          d[t](0, 0) = -1e-7;
        }
      mpi::communicator c;
      if (positive_counter && c.rank() == 0)
        std::cerr << " WARNING: " << positive_counter
                  << " points in the hybridization function are positive and "
                     "hence manually set to -1e-7 "
                  << std::endl;
    }
    /// Call operator
    /**
      The dets contain triplets made of
       0- the time of the operator,
       1- its inner_index (i.e. index of the color in its block and
       2- its color number (which is inserted too to simplify things)
      */
    double operator()(std::tuple<qmc_time_t, int, int> const &x,
                      std::tuple<qmc_time_t, int, int> const &y) const {
      double res = d[triqs::gfs::closest_mesh_pt(double(
          std::get<0>(x) - std::get<0>(y)))](std::get<1>(x), std::get<1>(y));
      return (std::get<0>(x) >= std::get<0>(y)
                  ? res
                  : -res); // get<0>(x),get<0>(y) are qmc_time_t, the wrapping
                           // is automatic in the - operation, but need to
                           // compute the sign
    }
  };

  typedef triqs::det_manip::det_manip<hybridization_adaptor> det_type;
  std::vector<det_type> dets;
  std::pair<int, int> current_color_insert, current_color_remove,
      current_color_shift;
  colored_const_iterator current_op1_insert, current_op2_insert;
  colored_const_iterator current_op1_remove, current_op2_remove;

  enum class status_t {
    neutral,
    tried_insert,
    tried_remove,
    tried_shift
  }; // state of the object
  status_t status_insert, status_remove, status_shift;

public:
  hybridization_dets(t_ordered_colored_c_ops const &ops_map,
                     gf_struct_t const &gf_struct_,
                     block_gf<imtime> const &hyb_raw)
      : HybOpmapImpl(ops_map, gf_struct_), status_insert(status_t::neutral),
        status_remove(status_t::neutral), status_shift(status_t::neutral) {
    mpi::communicator comm;
    for (int i = 0; i < hyb_raw.size(); i++) {
      if (!is_gf_real(hyb_raw[i])) {
        if (comm.rank() == 0)
          std::cerr << "Warning: The Delta(tau) block number " << i
                    << " is not real in tau space\n";
      }
      dets.emplace_back(hybridization_adaptor(real(hyb_raw[i])), 100);
    }
  };

  gf_struct_t const &gf_struct() const { return HybOpmapImpl::gf_struct; }

  // read only access to the det (for measure)
  det_type const &det(int k) const { return dets[k]; }

  // number of hybridization lines /segments on line 'color'
  int seg_number(int color) const { return HybOpmapImpl::seg_number(color); }

  int antiseg_number(int color) const {
    return HybOpmapImpl::antiseg_number(color);
  }

  // given a color and a int, returns the corresponding *hybridized* segment
  segment find(int color, size_t op_index) const {
    return HybOpmapImpl::find(color, op_index);
  }

  //*************** Addition ******************************

  /*
   * Insert two operators op1, op2 into the HybOpMapImpl, attempt and complete
   * insertion in determinant, update status Returns the determinant ratio NB:
   * operators do not need to be end of a segment
   */
  double add(colored_const_iterator const &op1,
             colored_const_iterator const &op2) {
    auto r = try_add(op1, op2);
    complete_add();
    return r;
  }

  /*
   * Insert two operators op1, op2 into the HybOpMapImpl, attempt insertion in
   * determinant and update status Returns the determinant ratio
   * @note operators do not need to be end of a segment
   * @note order of op1 op2 does not matter here
   */
  double try_add(colored_const_iterator const &op1,
                 colored_const_iterator const &op2) {
    assert(status_insert == status_t::neutral);
    assert(op1->dagger != op2->dagger);
    assert(op1->color == op2->color);
    current_color_insert =
        color_to_block_and_inner_index_impl(op1->color, gf_struct());
    this->insert_in_map(op1, op2);
    double det_ratio =
        (op1->dagger)
            ? dets[current_color_insert.first].try_insert(
                  find_index(op1), find_index(op2),
                  std::make_tuple(op1->tau, current_color_insert.second,
                                  op1->color),
                  std::make_tuple(op2->tau, current_color_insert.second,
                                  op2->color))
            : dets[current_color_insert.first].try_insert(
                  find_index(op2), find_index(op1),
                  std::make_tuple(op2->tau, current_color_insert.second,
                                  op2->color),
                  std::make_tuple(op1->tau, current_color_insert.second,
                                  op1->color));

    status_insert = status_t::tried_insert;
    return det_ratio;
  }

  /*
   * Complete insertion into determinant and update status
   */
  void complete_add() {
    assert(status_insert == status_t::tried_insert);
    dets[current_color_insert.first].complete_operation();
    status_insert = status_t::neutral;
  }

  /*
   * Reject insertion into determinant, removes inserted pair from HybOpMapImpl
   * and update status
   */
  void reject_add() {
    if (status_insert == status_t::neutral)
      return; // can be called also from neutral state
    assert(status_insert == status_t::tried_insert);
    this->reject_insert_in_map(block_and_inner_index_to_color_impl(
        current_color_insert.first, current_color_insert.second, gf_struct()));
    dets[current_color_insert.first].reject_last_try();
    status_insert = status_t::neutral;
  }

  //*************** Removal ******************************

  /*
   * Remove two operators op1, op2 at the two times from HybOpMapImpl, attempt
   * and complete removal from determinant Returns determinant ratio NB:
   * operators do not need to be end of a segment
   */
  double remove(colored_const_iterator const &op1,
                colored_const_iterator const &op2) {
    auto r = try_remove(op1, op2);
    complete_remove();
    return r;
  }

  /*
   * attempt removal from determinant
   * Returns determinant ratio
   * NB: operators do not need to be end of a segment
   */
  double try_remove(colored_const_iterator const &op1,
                    colored_const_iterator const &op2) {
    assert(status_remove == status_t::neutral);
    assert(op1->dagger != op2->dagger);
    assert(op1->color == op2->color);
    current_color_remove =
        color_to_block_and_inner_index_impl(op1->color, gf_struct());
    current_op1_remove = op1;
    current_op2_remove = op2;
    double det_ratio = (op1->dagger)
                           ? dets[current_color_remove.first].try_remove(
                                 find_index(op1), find_index(op2))
                           : dets[current_color_remove.first].try_remove(
                                 find_index(op2), find_index(op1));
    status_remove = status_t::tried_remove;
    return det_ratio;
  }

  /*
   * Complete removal from determinant and remove both operators from
   * HybOpMapImpl
   */
  void complete_remove() {
    assert(status_remove == status_t::tried_remove);
    dets[current_color_remove.first].complete_operation();
    this->remove_from_map(current_op1_remove);
    this->remove_from_map(current_op2_remove);
    status_remove = status_t::neutral;
  }

  /*
   * Reject removal from determinant and update status
   */
  void reject_remove() {
    if (status_remove == status_t::neutral)
      return; // can be called also from neutral state
    assert(status_remove == status_t::tried_remove);
    dets[current_color_remove.first].reject_last_try();
    status_remove = status_t::neutral;
  }

  /**  second version [special spin-spin+ move_move] where upon try_remove, the
   * operator is actually removed from hyb_opmap and put back at reject remove
   * **/

  /*
   * attempt removal from determinant
   * Returns determinant ratio
   * NB: operators do not need to be end of a segment
   */
  double try_remove2(colored_const_iterator const &op1,
                     colored_const_iterator const &op2) {
    assert(status_remove == status_t::neutral);
    assert(op1->dagger != op2->dagger);
    assert(op1->color == op2->color);
    current_color_remove =
        color_to_block_and_inner_index_impl(op1->color, gf_struct());
    current_op1_remove = op1;
    current_op2_remove = op2;
    double det_ratio = (op1->dagger)
                           ? dets[current_color_remove.first].try_remove(
                                 find_index(op1), find_index(op2))
                           : dets[current_color_remove.first].try_remove(
                                 find_index(op2), find_index(op1));
    this->remove_from_map(current_op1_remove);
    this->remove_from_map(current_op2_remove);
    status_remove = status_t::tried_remove;
    return det_ratio;
  }

  /*
   * Complete removal from determinant
   * @note contrary to complete_remove(), op1/op2 have already been removed
   */
  void complete_remove2() {
    assert(status_remove == status_t::tried_remove);
    dets[current_color_remove.first].complete_operation();
    status_remove = status_t::neutral;
  }

  /*
   * Reject removal from determinant and update status
   * insert op1, op2 back into hyb_opmap
   */
  void reject_remove2(colored_const_iterator const &op1,
                      colored_const_iterator const &op2) {
    if (status_remove == status_t::neutral)
      return; // can be called also from neutral state
    assert(status_remove == status_t::tried_remove);
    dets[current_color_remove.first].reject_last_try();
    this->insert_in_map(op1, op2);

    status_remove = status_t::neutral;
  }

  //*************** Shift operator******************************

  double shift(qmc_time_t const &old_tau,
               colored_const_iterator const &new_op) {
    std::cerr << "shifting..." << std::endl;
    auto r = try_shift(old_tau, new_op);
    std::cerr << "complete shifting..." << std::endl;
    complete_shift();
    return r;
  }

  double try_shift(qmc_time_t const &old_tau,
                   colored_const_iterator const &new_op) {
    assert(status_shift == status_t::neutral);
    current_color_shift =
        color_to_block_and_inner_index_impl(new_op->color, gf_struct());
    this->change_in_map(old_tau, new_op);
    std::cerr << "end change in map" << std::endl;
    double det_ratio =
        new_op->dagger
            ? dets[current_color_shift.first].try_change_row(
                  find_index(new_op),
                  std::make_tuple(new_op->tau, current_color_shift.second,
                                  new_op->color))
            : dets[current_color_shift.first].try_change_col(
                  find_index(new_op),
                  std::make_tuple(new_op->tau, current_color_shift.second,
                                  new_op->color));

    std::cerr << "soon done" << std::endl;
    status_shift = status_t::tried_shift;
    return det_ratio;
  }

  void complete_shift() {
    assert(status_shift == status_t::tried_shift);
    dets[current_color_shift.first].complete_operation();
    status_shift = status_t::neutral;
  }

  void reject_shift(qmc_time_t const &new_tau,
                    colored_const_iterator const &old_op) {
    if (status_shift == status_t::neutral)
      return; // can be called also from neutral state
    assert(status_shift == status_t::tried_shift);

    this->change_in_map(new_tau, old_op);
    dets[current_color_shift.first].reject_last_try();
    status_shift = status_t::neutral;
  }

#ifdef CTSEG_DEBUG
  //----------------------------------
  // returns delta(op1->tau - op2 -> tau2)
  // double check_mat_inv (){
  // for(int i=0;i<dets.size();i++)  dets[i].check_mat_inv() ;
  //}
#endif

  friend std::ostream &operator<<(std::ostream &out,
                                  hybridization_dets const &h) {
    return out << static_cast<HybOpmapImpl const &>(h);
  }
  friend void h5_write(h5::group fg, std::string subgroup_name,
                       hybridization_dets const &h) {
    h5::group gr = fg.create_group(subgroup_name);
    h5_write(gr, "hyb_dets", static_cast<HybOpmapImpl const &>(h));
  }
};

} // namespace triqs_ctseg
