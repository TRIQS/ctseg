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
#include "./trace_c_ops.hpp"

namespace triqs_ctseg {

trace_c_ops::trace_c_ops(const qmc_parameters *params_, int *id_)
    : params(params_), id(id_), ops(params->n_color),
      _full_lines(params->n_color, false),
      _overlap_matrix(params->n_color, params->n_color),
      sign(params->n_color, 1.0), overlaps(params->n_color, 0.0) {
  assert(params->n_color == params->n_color);
  status = status_t::neutral;
  _overlap_matrix() = 0;
  for (int i = 0; i < params->n_color; i++) {
    _full_lines[i] = (i % 2 == 1); // alternated
  }
  for (int i = 0; i < params->n_color; i++) {
    for (int j = 0; j < params->n_color; j++) {
      _overlap_matrix(i, j) =
          (_full_lines[i] && _full_lines[j]) ? double(params->beta) : 0.0;
    }
  }
}

//***************   Segment insertion ******************************
// if empty/full line, modifies state of line from s.tl to s.tr
std::pair<double, segment>
trace_c_ops::try_insert_segment(segment_desc const &s) {
  assert(status == status_t::neutral);

  color = s.color;
  auto rn = ops.right_neighbor(s.tr, color);
  is_segment =
      (ops.op_number(color) == 0) ? !(full_lines(color)) : !(rn->dagger);
  // inserting in the operator list
  colored_const_iterator l, r;

  l = ops.insert(s.tl, color, !is_segment, rn);

  auto it_l = decolor(l);
  // if very first operator, right_occ initialized to false->update using
  // full_line
  // NB this right_occupations will be changed just below for color color
  it_l->right_occupations = (ops.total_op_number() == 1)
                                ? _full_lines
                                : ops.cyclic_left(it_l)->right_occupations;

  // set right_occupations before inserting 'right' operator
  try {
    r = ops.insert(s.tr, color, is_segment, l);
  } catch (insertion_error const &e) {
    ops.erase(l);
    throw;
  }

  // update "right_occupations": all operators between left(included) and
  // right(excluded)
  auto it_r = decolor(r);
  it_r->right_occupations = ops.cyclic_left(it_r)->right_occupations;
  auto it = it_l;
  do {
    it->right_occupations[color] = !(it->right_occupations[color]);
    it = ops.cyclic_right(it);
  } while (it != it_r);

#ifdef DEBUG_CTSEG_VERIF
  mpi::communicator w;
  for (auto const &op : ops) {
    if (op.dagger != !op.right_occupations[op.color]) {
      std::cerr << " [" << w.rank() << "][#" << *id << "] try insert " << s.tl
                << s.tr << "[full_lines[col=" << color
                << "]] = " << full_lines(color) << std::endl;
      std::cerr << "verif " << ops << std::endl;
      throw "stop";
    }
  }
#endif

  // update length and overlaps for future reuse in accept
  seg = segment{l, r};
  length = seg.length();
  compute_overlaps(seg);

  double ln_trace_ratio =
      params->mu(color) * length - dot(params->U(range(), color), overlaps);
  if (params->dynamical_U)
    ln_trace_ratio += dynamical_contribution(seg);

  status = status_t::tried_insert;
  return {ln_trace_ratio, seg};
}

//-------------------------------------

// return sign_ratio
double trace_c_ops::complete_insert_segment() {
  assert(status == status_t::tried_insert);

  // update full_lines if we just added an anti segment on a full_line
  if (ops.seg_number(color) == 1 && _full_lines[color])
    _full_lines[color] = false;

  // update overlaps
  for (int i = 0; i < params->n_color; i++) {
    _overlap_matrix(color, i) += overlaps[i];
    _overlap_matrix(i, color) += overlaps[i];
  }
  _overlap_matrix(color, color) += length;

#ifdef CTSEG_DEBUG_VERIF
  if (overlap_matrix(1, 0) > params->beta + 0.00001) {
    std::cerr << " error com ins " << _overlap_matrix << overlaps << "\n"
              << *this << std::endl;
    throw "stop";
  }
  if (overlap_matrix(0, 0) < -.0000001) {
    std::cerr << " error com rem " << _overlap_matrix << overlaps << "\n"
              << *this << std::endl;
    throw "stop";
  }
#endif

  status = status_t::neutral;
  return update_sign(color);
}

//-------------------------------------
void trace_c_ops::rm_segment_and_update() {

  // update "right_occupations" between the 2 operators to be removed
  auto it = ops.cyclic_right(decolor(seg.l));
  while (it != decolor(seg.r)) {
    it->right_occupations[color] = !it->right_occupations[color];
    it = ops.cyclic_right(it);
  }
  // erase operators from maps
  ops.erase(seg.l);
  ops.erase(seg.r);

  // update full_lines if we remove the last segment of a color
  if (ops.seg_number(color) == 0) {
    _full_lines[color] = !is_segment;
#ifdef CTSEG_DEBUG_VERIF
    for (auto const &op : ops) {
      if (op.right_occupations[color] != full_lines(color)) {
        std::cout << "rm_segment_and_update verif " << ops << std::endl;
        throw "stop";
      }
    }
#endif
  }
}

//------------------------------------

/* remove a segment */
double trace_c_ops::try_remove_segment(segment const &s) {
  assert(status == status_t::neutral);
  seg = s;
  color = s.l->color;
  is_segment = !(s.l->dagger);
  s_desc = segment_desc{s.l->tau, s.r->tau, color};

  length = seg.length();
  compute_overlaps(seg);

  double ln_trace_ratio =
      -params->mu(color) * length + dot(params->U(range(), color), overlaps);
  if (params->dynamical_U)
    ln_trace_ratio -= dynamical_contribution(seg);

  status = status_t::tried_remove;
  return ln_trace_ratio;
}

//-------------------------------------

std::pair<double, segment_desc> trace_c_ops::complete_remove_segment() {
  assert(status == status_t::tried_remove);
  //assert(ops.are_neighbors(seg.l, seg.r));
  rm_segment_and_update();
  // updating cached observables
  for (int i = 0; i < params->n_color; i++) {
    _overlap_matrix(color, i) -= overlaps[i];
    _overlap_matrix(i, color) -= overlaps[i];
  }
  _overlap_matrix(color, color) -= length;

#ifdef CTSEG_DEBUG_VERIF
  if (overlap_matrix(1, 0) > params->beta + 0.00001) {
    std::cerr << " error com rem " << _overlap_matrix << overlaps << *this
              << std::endl;
    throw "stop";
  }
  if (overlap_matrix(0, 0) < -.0000001) {
    std::cerr << " error com rem " << _overlap_matrix << overlaps << *this
              << std::endl;
    throw "stop";
  }
#endif
  status = status_t::neutral;
  return {update_sign(color), s_desc};
}

//-------------------------------------
std::pair<double, colored_const_iterator>
trace_c_ops::try_shift_operator(colored_const_iterator const &c_old,
                                qmc_time_t tau_new) {
  // assert(status==status_t::neutral || );
  double ln_trace_ratio = 0.0;
  if (params->dynamical_U)
    ln_trace_ratio -= dynamical_contribution(c_old);

  auto tau_old = c_old->tau;
  color = c_old->color;
  auto old_r_occ = decolor(c_old)->right_occupations;
  auto dagger = c_old->dagger;

  auto R_old = ops.cyclic_right(decolor(c_old));
  auto L_old_same_line = ops.cyclic_left(c_old);
  auto R_old_same_line = ops.cyclic_right(c_old);
  auto c_new = ops.shift_operator(c_old, tau_new); // erases c_old

  auto c_new_decolored = decolor(c_new);
  auto R_new = ops.cyclic_right(c_new_decolored);
  auto L_new_same_line = ops.cyclic_left(c_new);
  auto R_new_same_line = ops.cyclic_right(c_new);
  //assert(R_new_same_line->tau == R_old_same_line);
  //assert(L_new_same_line->tau == L_old_same_line);
  auto L = L_new_same_line;
  auto R = R_new_same_line;

  bool update_right_of_old = L->tau - tau_old < L->tau - c_new->tau;
  auto first = (update_right_of_old) ? R_old : R_new;
  auto last = (update_right_of_old) ? R_new : R_old;

  c_new_decolored->right_occupations =
      ops.cyclic_left(c_new_decolored)->right_occupations;
  if (!update_right_of_old)
    c_new_decolored->right_occupations[color] =
        !c_new_decolored->right_occupations[color];
  if (R_old->tau == R_new->tau)
    c_new_decolored->right_occupations = old_r_occ;

  auto it = first;
  while (it != last) {
    if (it->tau != c_new->tau)
      it->right_occupations[color] = !it->right_occupations[color];
    it = ops.cyclic_right(it);
  }

  // compute change in overlaps: insert back old op, measure overlap, remove
  auto c_old_temp = ops.insert(tau_old, color, !dagger); // no hint for now
  auto c_old_temp_decol = decolor(c_old_temp);
  c_old_temp_decol->right_occupations = old_r_occ;
  auto seg = (ops.cyclic_left(c_old_temp) == c_new)
                 ? segment{c_new, c_old_temp}
                 : segment{c_old_temp, c_new};
  compute_overlaps(seg);
  length = seg.length();
  ops.erase(c_old_temp);

  // compute overlap variation
  ln_trace_ratio =
      params->mu(color) * length - dot(params->U(range(), color), overlaps);
  if (params->dynamical_U)
    ln_trace_ratio += dynamical_contribution(c_new);

  // update _overlaps_matrix
  for (int i = 0; i < params->n_color; i++) {
    _overlap_matrix(color, i) += overlaps[i];
    _overlap_matrix(i, color) += overlaps[i];
  }
  _overlap_matrix(color, color) += length;

  // status = status_t::tried_shift;
  return {ln_trace_ratio, c_new};
}

double trace_c_ops::complete_shift_operator() {
  // assert(status==status_t::tried_shift);

  status = status_t::neutral;
  return update_sign(0) * update_sign(1);
}

//-------------------------------------
double trace_c_ops::swap_full_empty(int c1, int c2) {
  if (ops.seg_number(c1) > 0 || ops.seg_number(c2) > 0)
    return 0.0;
  _full_lines[c1] = !_full_lines[c1];
  _full_lines[c2] = !_full_lines[c2];
  // flip right_occupations for color c, for all operators
  for (auto const &op : ops) {
    op.right_occupations[c1] = !op.right_occupations[c1];
    op.right_occupations[c2] = !op.right_occupations[c2];
  }

  int fact =
      full_lines(c1)
          ? -1
          : 1; // if c1 is *now* full, it was like an antisegment(empty) before
  int full_line_before = (fact == 1) ? c2 : c1;
  // auto delta_U = data->U(all,c1)-data->U(all,c2); //does not compile: dot
  // complains
  range all;
  vector<double> delta_U = params->U(all, c1) - params->U(all, c2);
  double ln_trace_ratio =
      (params->beta * (params->mu(c1) - params->mu(c2)) -
       dot(_overlap_matrix(all, full_line_before), delta_U)) *
      fact;

  // swap indices c2 and c1
  auto temp_overlap = _overlap_matrix;
  for (size_t i = 0; i < ops.n_colors(); ++i) {
    if (i != c1 && i != c2) {
      _overlap_matrix(i, c2) = temp_overlap(i, c1);
      _overlap_matrix(c2, i) = temp_overlap(c1, i);
      _overlap_matrix(i, c1) = temp_overlap(i, c2);
      _overlap_matrix(c1, i) = temp_overlap(c2, i);
    }
  }
  _overlap_matrix(c1, c1) = temp_overlap(c2, c2);
  _overlap_matrix(c2, c2) = temp_overlap(c1, c1);
  return ln_trace_ratio;
}

//-------------------------------------

void trace_c_ops::compute_overlaps(segment const &s) {
  overlaps() = 0;
  int fact = s.l->dagger ? -1 : 1; // if antiseg = -1
  auto it1 = decolor(s.l), it2 = ops.cyclic_right(it1);
  auto const last = decolor(s.r);
  while (it1 != last) {
    const double d = fact * double(it1->tau - it2->tau);
    for (int i = 0; i < ops.n_colors(); ++i)
      overlaps(i) += it1->right_occupations[i] * d;
    it1 = it2;
    it2 = ops.cyclic_right(it1);
  }

  if (ops.total_op_number() == 2)
    for (int i = 0; i < ops.n_colors(); ++i)
      if ((i != s.l->color) && full_lines(i))
        overlaps(i) = s.length();

  overlaps(s.l->color) = 0;

#ifdef CTSEG_DEBUG_VERIF
  for (int i = 0; i < ops.n_colors(); ++i)
    if (std::abs(overlaps(i)) > params->beta + 0.00001) {
      std::cerr << " error " << overlaps << *this << std::endl;
      throw "stop";
    }
#endif
}

//-------------------------------------
double trace_c_ops::compute_length(int c) {
  double length = 0.0;
  if (ops.seg_number(c) > 0) { // at least one segment
    auto begin = ops.begin(c);
    if (begin->dagger)
      begin = ops.cyclic_right(begin);
    auto last = begin;
    auto it_a = begin;
    auto it_c = ops.cyclic_right(begin);

    do {
      length += segment{it_a, it_c}.length();
      it_a = ops.cyclic_right(it_c);
      it_c = ops.cyclic_right(it_a);
    } while (it_a != begin);

  } else { // empty/full_line
    length = _full_lines[c] ? double(params->beta) : 0.0;
  }

  return length;
}
//-------------------------------------
void trace_c_ops::check_overlap_matrix_from_scratch() {
  auto overlap_matrix = _overlap_matrix; // copy
  overlap_matrix() = 0.0;
  for (int c = 0; c < ops.n_colors(); c++) {
    // std::cerr << "\t looking at color " << c << std::endl;
    if (ops.seg_number(c) > 0) { // at least one segment
      auto begin = ops.begin(c);
      if (begin->dagger)
        begin = ops.cyclic_right(begin);
      auto last = begin;
      auto it_a = begin;
      auto it_c = ops.cyclic_right(begin);
      do {
        auto s = segment{it_a, it_c};
        // std::cerr << "\t\t looking at seg = "
        // <<double(s.l->tau)<<","<<double(s.r->tau)<<") on line "<<c<<std::endl;
        // std::cerr << "\t\t length  = " << s.length() << std::endl;
        this->compute_overlaps(s);

        // std::cerr << "\t\t overlap of seg of " << c << " with rest : " <<
        // overlaps() << std::endl;
        for (int c2 = 0; c2 < ops.n_colors(); c2++)
          if (c != c2)
            overlap_matrix(c, c2) += overlaps(c2);
        overlap_matrix(c, c) += s.length();
        it_a = ops.cyclic_right(it_c);
        it_c = ops.cyclic_right(it_a);
        assert(it_c->dagger);
      } while (it_a != begin);
      // std::cerr << "\t overlap of " << c << " with rest : " << overlaps() <<
      // std::endl; std::cerr << "\t overlap_matrix " << overlap_matrix <<
      // std::endl;

    } else { // empty/full_line
      // std::cerr << "\t full_lines "<<c<<" = "<<_full_lines[c] << std::endl;
      if (_full_lines[c]) {
        for (int c2 = 0; c2 < ops.n_colors(); c2++)
          if (c != c2)
            overlap_matrix(c, c2) = this->compute_length(c2);
        overlap_matrix(c, c) = double(params->beta);
      }
    }
  }
  // now compare to _overlap_matrix
  auto A = overlap_matrix - _overlap_matrix;
  auto sum_A = sum(A);
#ifdef CTSEG_DEBUG
  std::cerr << "******* sum of errors on overlaps = " << sum_A << std::endl;
#endif
  if (std::abs(sum_A) > 1e-10)
    TRIQS_RUNTIME_ERROR << "too large deviation in overlap_matrix : computed "
                        << overlap_matrix << " instead of " << _overlap_matrix;
}

//-------------------------------------
double trace_c_ops::dynamical_contribution(colored_const_iterator cit) {
  double contrib = 0.0;
  for (auto it = ops.begin(); it != ops.end(); ++it) {
    if (it != decolor(cit)) {
      auto c = params
                   ->K[closest_mesh_pt(double(it->tau - cit->tau))](it->color,
                                                                    cit->color)
                   .real();
      if (it->dagger == cit->dagger)
        contrib += c;
      else
        contrib -= c;
    }
  }
  return contrib;
}
//-------------------------------------

double trace_c_ops::dynamical_contribution(segment const &s) {
  double contrib = 0.0;
  for (auto it = ops.begin(); it != ops.end(); ++it) {
    if (it != decolor(s.l)) {
      auto c = params
                   ->K[closest_mesh_pt(double(it->tau - s.l->tau))](it->color,
                                                                    s.l->color)
                   .real();
      if (it->dagger == s.l->dagger)
        contrib += c;
      else
        contrib -= c;

      if (it != decolor(s.r)) {
        auto c2 = params
                      ->K[closest_mesh_pt(double(it->tau - s.r->tau))](
                          it->color, s.r->color)
                      .real();
        if (it->dagger == s.r->dagger)
          contrib += c2;
        else
          contrib -= c2;
      }
    }
  }
  return contrib;
}
//-------------------------------------

std::ostream &operator<<(std::ostream &out, trace_c_ops const &t) {
  out << t.ops;
  out << std::endl << "|FULL_LINES= ";
  for (auto it = t._full_lines.begin(); it != t._full_lines.end(); ++it)
    out << *it << ", ";
  return out;
}
//-------------------------------------
void h5_write(h5::group fg, std::string subgroup_name, trace_c_ops const &t) {
  h5::group gr = fg.create_group(subgroup_name);
  h5_write(gr, "ops", t.ops);
  h5_write(gr, "beta", double(t.params->beta));
  h5_write(gr, "overlap_matrix", t._overlap_matrix);
  array<double, 2> a(t._full_lines.size(), 2);
  for (int i = 0; i < t._full_lines.size(); i++) {
    a(i, 0) = i;
    a(i, 1) = t._full_lines[i];
  }
  h5_write(gr, "full_lines", a);
}

} // namespace triqs_ctseg
