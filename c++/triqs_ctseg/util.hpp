#pragma once

// Same as std::lower_bound, but i-th element of vector is returned by f[i]
// f is called on 0:N strictly
long lower_bound(auto f, long N, auto const &value) {
  long first = 0, count = N;

  while (count > 0) {
    long step = count / 2;
    assert(first + step < N);
    if (f(first + step) < value) {
      first += step + 1;
      count -= step + 1;
    } else
      count = step;
  }
  return first;
}

// If we have an ordered det_manip d, we get the lower_bound index for x
long det_lower_bound_x(auto const &d, auto const &x) {
  return lower_bound([&d](long i) { return d.get_x(i).first; }, d.size(), x);
}
long det_lower_bound_y(auto const &d, auto const &y) {
  return lower_bound([&d](long i) { return d.get_y(i).first; }, d.size(), y);
}
