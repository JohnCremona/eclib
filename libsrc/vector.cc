// vector.cc: manage implementations of integer vector classes
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
// 
// This file is part of the eclib package.
// 
// eclib is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// eclib is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with eclib; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include <eclib/marith.h>
#include <eclib/vector.h>
#include "random.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i

#include "vec.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l

#include "vec.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar bigint
#define vec vec_m
#define mat mat_m
#define subspace subspace_m
#define svec svec_m
#define smat smat_m
#define smat_elim smat_m_elim

#include "vec.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim

vec_m to_vec_m(const vec_i& v)
{
  const vector<int> & vi = v.get_entries();
  vector<bigint> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), [](const int& x) {return bigint(x);});
  return vec_m(w);
}

vec_m to_vec_m(const vec_l& v)
{
  const vector<long> & vi = v.get_entries();
  vector<bigint> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), [](const long& x) {return bigint(x);});
  return vec_m(w);
}

vec_i to_vec_i(const vec_m& v)
{
  const vector<bigint> & vi = v.get_entries();
  auto toint = [](const bigint& a) {return is_int(a)? I2int(a) : int(0);};
  vector<int> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), toint);
  return vec_i(w);
}

vec_i to_vec_i(const vec_l& v)
{
  const vector<long> & vi = v.get_entries();
  auto toint = [](const long& a) {return int(a);};
  vector<int> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), toint);
  return vec_i(w);
}

vec_l to_vec_l(const vec_m& v)
{
  const vector<bigint> & vi = v.get_entries();
  auto tolong = [](const bigint& a) {return is_long(a)? I2long(a) : long(0);};
  vector<long> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), tolong);
  return vec_l(w);
}

vec_l to_vec_l(const vec_i& v)
{
  const vector<int> & vi = v.get_entries();
  auto tolong = [](const int& a) {return long(a);};
  vector<long> w(vi.size());
  std::transform(vi.begin(), vi.end(), w.begin(), tolong);
  return vec_l(w);
}
