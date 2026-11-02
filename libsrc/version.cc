// version.cc: implementation of version functions
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

#include <iostream>
#include "eclib/version.h"

string eclib_version()
{
  string v = ECLIB_VERSION;
  return v;
}

vector<int> eclib_date()
{
  string y = ECLIB_YEAR,  m = ECLIB_MONTH, d = ECLIB_DAY;
  return {atoi(y.c_str()), atoi(m.c_str()), atoi(d.c_str())};
}

void show_version(ostream& os)
{
  os << "eclib version " << eclib_version() << ", ";
#ifdef NO_MPFP
  os << "using NTL ZZs but no multiprecision floating point";
#else
  os << "using NTL ZZs and NTL real and complex multiprecision floating point";
#endif
  os << endl;
}

int sgn(int a)   {return (a==0? 0: (a>0? 1: -1));}

int compare_eclib_version(int y, int m, int d)
{
  vector<int> date = eclib_date();
  int s;
  // compare years
  s = sgn(date[0] - y);
  if (s!=0) // different years, OK if y is smaller
    return s;
  // same year, compare months
  s = sgn(date[1] - m);
  if (s!=0) // different months, OK if m is smaller
    return s;
  // same year and month, compare days
  return sgn(date[2] - d);
}
