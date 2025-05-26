// matrix.cc: manage implementation of integer matrix classes
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
 
#include <eclib/matrix.h>

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef SCALAR_OPTION

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i
#define SCALAR_OPTION 1

#include "mat.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef SCALAR_OPTION

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l
#define SCALAR_OPTION 2

#include "mat.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef SCALAR_OPTION

#define scalar bigint
#define vec vec_m
#define mat mat_m
#define subspace subspace_m
#define SCALAR_OPTION 0

#include "mat.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef SCALAR_OPTION

mat_m to_mat_m(const mat_i& m)
{
  const vector<int> & mij = m.get_entries();
  vector<bigint> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), [](const int& x) {return bigint(x);});
  return mat_m(m.nrows(), m.ncols(), n);
}

mat_m to_mat_m(const mat_l& m)
{
  const vector<long> & mij = m.get_entries();
  vector<bigint> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), [](const long& x) {return bigint(x);});
  return mat_m(m.nrows(), m.ncols(), n);
}

mat_i to_mat_i(const mat_m& m)
{
  const vector<bigint> & mij = m.get_entries();
  auto toint = [](const bigint& a) {return is_int(a)? I2int(a) : int(0);};
  vector<int> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), toint);
  return mat_i(m.nrows(), m.ncols(), n);
}

mat_l to_mat_l(const mat_m& m)
{
  const vector<bigint> & mij = m.get_entries();
  auto tolong = [](const bigint& a) {return is_long(a)? I2long(a) : long(0);};
  vector<long> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), tolong);
  return mat_l(m.nrows(), m.ncols(), n);
}
