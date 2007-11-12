// version.cc: implementation of function show_version()
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include <iostream>

using namespace std;

void show_version()
{
  cerr << "Version compiled on " << __DATE__ << " at " << __TIME__ << " by GCC " << __VERSION__ << "\n";
  cerr << "using base arithmetic option ";
#ifdef LiDIA_ALL
    cerr << "LiDIA_ALL (LiDIA bigints and multiprecision floating point)";
#else
#ifdef LiDIA_INTS
  cerr << "LiDIA_INTS (LiDIA bigints and no multiprecision floating point)";
#else
#ifdef NTL_ALL
  cerr << "NTL_ALL (NTL bigints and multiprecision floating point)";
#else
#ifdef NTL_INTS
  cerr << "NTL_INTS (NTL bigints and no multiprecision floating point)";
#else
  cerr << "libg++, no multiprecision floating point";
#endif
#endif
#endif
#endif
  cerr<<endl;
}

