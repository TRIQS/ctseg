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
#include "./types.hpp"
#include <triqs/utility/dressed_iterator.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/time_pt.hpp>

namespace triqs_ctseg {

typedef triqs::utility::time_pt
    qmc_time_t; // a double on a very thin grid to avoid collisions
class insertion_error : triqs::runtime_error {
}; // a specific exception to handle very exceptional insertion errors in maps.

/**
 * A storage for t-ordered colored canonical operators with :
 *
 *- a const_iterator iterating on all operators.
 *- a colored_const_iterator iterating on operators on a given color.
 *- cyclic_left/right method to advance these operators right and left,
 *cyclically.
 *- insert/erase method
 */
class t_ordered_colored_c_ops {

  // operator descriptor: color, dagger[creation/annihilation], occupations of
  // all lines to the right of the operator note: right_occupations is mutable :
  // it allows to have a const_iterator only and allow to modify it... simpler
  // code.
  typedef std::vector<bool> occupations_t;
  struct c_op_desc {
    int color;
    bool dagger;
    mutable occupations_t right_occupations;
  };

  // the time ordered list of operators containing all colors
  // time are ordered in decreasing order, in agreement with the whole physics
  // literature.
  typedef std::map<qmc_time_t, c_op_desc, std::greater<qmc_time_t>> fullopmap_t;

  // a time ordered list of operators for one given color
  // for each color, we store a "submap" of iterators of a the given color
  typedef std::map<qmc_time_t, fullopmap_t::const_iterator,
                   std::greater<qmc_time_t>>
      opmap_t;
  fullopmap_t fullopmap;
  std::vector<opmap_t> opmaps;

  // We use utility::dressed_iterator to dress iterators
  // _cdress is a simple struct of refs to dress the iterators (Cf doc)
  // It is used for both the opmap and the fullopmap iterators, hence 2
  // constructors. the mutable part is accessible only from this class and
  // trace_c_ops class below.
  struct _cdress {
    _cdress(fullopmap_t::const_iterator const &_it)
        : tau(_it->first), color(_it->second.color), dagger(_it->second.dagger),
          right_occupations(_it->second.right_occupations) {}
    _cdress(opmap_t::const_iterator const &_it)
        : // _cdress(_it->second){} // with delegating constructor in C++11
          tau(_it->second->first), color(_it->second->second.color),
          dagger(_it->second->second.dagger),
          right_occupations(_it->second->second.right_occupations) {}
    qmc_time_t const &tau;
    int const &color;
    bool const &dagger;
    std::vector<bool> const &right_occupation() const {
      return right_occupations;
    }

  private:
    std::vector<bool> &right_occupations; // dangerous, can be modified even for
                                          // const object ...
    friend class t_ordered_colored_c_ops;
    friend class trace_c_ops;
    friend std::ostream &operator<<(std::ostream &out,
                                    t_ordered_colored_c_ops const &c);
    friend void h5_write(h5::group fg, std::string subgroup_name,
                         t_ordered_colored_c_ops const &c);
  };

public:
  // constructor
  t_ordered_colored_c_ops(int n_color) : opmaps(n_color) {}

  /************************ custom iterators and associated functions
   * ********************************************/
  typedef triqs::utility::dressed_iterator<fullopmap_t::const_iterator, _cdress>
      const_iterator;
  typedef triqs::utility::dressed_iterator<opmap_t::const_iterator, _cdress>
      colored_const_iterator;

  const_iterator begin() const { return fullopmap.begin(); }
  const_iterator end() const { return fullopmap.end(); }

  colored_const_iterator begin(int color) const {
    return opmaps[color].begin();
  }
  colored_const_iterator end(int color) const { return opmaps[color].end(); }

  // return iterator to operator at a given time
  const_iterator find(qmc_time_t t) const { return fullopmap.find(t); }
  colored_const_iterator find_colored(int color, qmc_time_t t) const {
    return opmaps[color].find(t);
  }

  /* decolor an iterator.
   * NB : this is the ONLY function to change an iterator to another.
   * In particular there is no implicit cast between operators.
   * \note: cannot decolor end()!
   */
  friend const_iterator decolor(colored_const_iterator const &ci) {
    return ci.get()->second;
  }

  // returns the iterator to the operator to the right of 'it' (cyclically)
  colored_const_iterator cyclic_right(colored_const_iterator it) const {
    auto color = it->color; // store it before incrementing
    ++it;
    if (it == end(color))
      it = begin(color);
    return it;
  }

  // returns the iterator to the operator to the right of 'it' (cyclically)
  const_iterator cyclic_right(const_iterator it) const {
    ++it;
    if (it == end())
      it = begin();
    return it;
  }

  // returns the iterator to the operator to the left of 'it' (cyclically)
  colored_const_iterator cyclic_left(colored_const_iterator it) const {
    if (it == begin(it->color))
      it = end(it->color);
    --it;
    return it;
  }

  // returns the iterator to the operator to the left of 'it' (cyclically)
  const_iterator cyclic_left(const_iterator it) const {
    if (it == begin())
      it = end();
    --it;
    return it;
  }
  /************************end of
   * iterators******************************************************/

  /// Distance between the right and left neighbors of a given pair {l,r}
  double space_around_decolor(colored_const_iterator const &l,
                              colored_const_iterator const &r) const {
    return double(cyclic_left(decolor(l))->tau - cyclic_right(decolor(r))->tau);
  }
  double space_around(colored_const_iterator const &l,
                      colored_const_iterator const &r) const {
    return double(cyclic_left(l)->tau - cyclic_right(r)->tau);
  }

  // return number of operators of color 'color' between l  (excluded) and r
  // (excluded)
  int op_count_in_between(colored_const_iterator const &l,
                          colored_const_iterator const &r, int color) const {
    int counter = 0;
    auto it = decolor(l);
    auto it_r = decolor(r);
    while (it != it_r) {
      it = cyclic_right(it);
      if (it->color == color)
        counter++;
    }
    return counter;
  }

  // number of colors
  size_t n_colors() const { return opmaps.size(); }

  // number of operators on line 'color'
  size_t op_number(int color) const { return opmaps[color].size(); }

  // total number of operators
  size_t total_op_number() const { return fullopmap.size(); }

  // Number of segments of a given color
  int seg_number(int color) const { return op_number(color) / 2; }

  /**
   * insert operator
   * \arg tau imaginary time
   * \arg color line
   * \arg dagger: operator type (creation/annihilation)
   * \arg hint: hint of position in map (as for std::map)
   * \return colored_const_iterator to op
   * \exception : insert_error if the time tau is already occupied.
   */
  colored_const_iterator insert(qmc_time_t const &tau, int color, bool dagger,
                                colored_const_iterator const &hint) {
    auto decolored_hint = (hint == end(color)) ? end() : decolor(hint);
    // we check that insertion is successfull simply by controlling the size of
    // the container. Cf std::map doc : with an hint, the insert method does not
    // return a bool to indicate success of insertion
    auto s = fullopmap.size();
    auto it = fullopmap.insert(
        decolored_hint,
        std::make_pair(tau,
                       c_op_desc{color, dagger, occupations_t(opmaps.size())}));
    if (fullopmap.size() != s + 1)
      throw insertion_error();
    s = opmaps[color].size();
    auto cit = opmaps[color].insert(hint, std::make_pair(tau, it));
    if (opmaps[color].size() != s + 1) {
      fullopmap.erase(it);
      throw insertion_error();
    } // in fact can not happen
    return cit;
  }

  /**
   * insert operator
   * \arg tau imaginary time
   * \arg color line
   * \arg dagger: operator type (creation/annihilation)
   * \return colored_const_iterator to op
   * \exception : insert_error if the time tau is already occupied.
   */
  colored_const_iterator insert(qmc_time_t const &t, int color, bool dagger) {
    return insert(t, color, dagger, end(color));
  }

  ///
  void erase(colored_const_iterator const &it) {
    auto c = it->color;
    fullopmap.erase(decolor(it));
    opmaps[c].erase(it);
  }

  /// shift an operator to tau_new. No check on neighbors positions
  colored_const_iterator shift_operator(colored_const_iterator const &it,
                                        qmc_time_t const &tau_new) {
    auto color = it->color;
    auto dagger = it->dagger;
    auto hint = cyclic_right(it);
    if (hint == it)
      hint = end(color);
    erase(it);
    return insert(tau_new, color, dagger, hint);
  }

  /*
   * Returns iterator to right (cyclic) neighbor of the given color of time tau.
   * If line is empty, returns begin(color)(==end(color))
   * Takes care of beta-periodicity.
   * If line is non empty, never returns end(color).
   * If tau is exactly on an operator (not ambiguous because we use qmc_time_t,
   * hence a grid), it returns this operator. TO RETURN its neighbour, change to
   * upper_bound, cf STL doc.
   */
  colored_const_iterator right_neighbor(qmc_time_t const &tau,
                                        int color) const {
    if (op_number(color) == 0)
      return begin(color);
    colored_const_iterator it = opmaps[color].lower_bound(tau);
    if (it == end(color))
      return begin(color);
    return it;
  }

  // returns true if there are no operators of the same color between s1 and s2
  // or between s2 and s1
  bool are_colored_neighbors(colored_const_iterator const &t1,
                             colored_const_iterator const &t2) const {
    return (t2 == cyclic_right(t1) || t1 == cyclic_right(t2));
  }

  // print
  friend std::ostream &operator<<(std::ostream &out,
                                  t_ordered_colored_c_ops const &c) {
    for (auto op = c.begin(); op != c.end(); ++op) {
      out << "(" << double(op->tau) << "," << op->color << "," << op->dagger
          << ":";
      for (int i = 0; i < op->right_occupations.size(); i++)
        out << op->right_occupations[i] << ",";
      out << ") - ";
    }
    return out;
  }
  friend void h5_write(h5::group fg, std::string subgroup_name,
                       t_ordered_colored_c_ops const &c) {
    h5::group gr = fg.create_group(subgroup_name);
    array<double, 2> a(c.total_op_number(), 4 + c.n_colors());
    int i = 0;
    for (auto op = c.begin(); op != c.end(); ++op, ++i) {
      a(i, 0) = double(op->tau);
      // a(i, 1) = op->tau.get_n();
      a(i, 2) = op->color;
      a(i, 3) = op->dagger;
      for (int j = 0; j < c.n_colors(); j++)
        a(i, 4 + j) = op->right_occupations[j];
    }

    if (c.total_op_number() > 0)
      h5_write(gr, "op_table", a);
  }
};
} // namespace triqs_ctseg
