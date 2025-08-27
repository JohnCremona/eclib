// version.cc: implementation of functions show_version(), eclib_version()
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
 
#include <iostream>
#include "eclib/templates.h"

string eclib_version()
{
  string v = VERSION; // defined by autotools, of the form v<yyyy><mm><dd>
  return v;
}

vector<int> eclib_date()
{
  vector<int> date;
  string v = eclib_version(); // 8 chars long, yyyymmdd
  date.push_back(atoi(v.substr(0,4).c_str())); // chars 0-4
  date.push_back(atoi(v.substr(4,2).c_str())); // chars 5-6
  date.push_back(atoi(v.substr(6,2).c_str())); // chars 7-8
  return date;
}

void show_version(ostream& os)
{
  os << "eclib version " << VERSION << ", ";
#ifdef NO_MPFP
  os << "using NTL bigints but no multiprecision floating point";
#else
  os << "using NTL bigints and NTL real and complex multiprecision floating point";
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
