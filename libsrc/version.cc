// version.cc: implementation of function show_version()
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
#include <eclib/templates.h>


void show_version()
{
  cerr << "Version compiled on " << __DATE__ << " at " << __TIME__ << " by GCC " << __VERSION__ << "\n";
#ifdef MPFP
  cerr << "using NTL bigints and NTL real and complex multiprecision floating point";
#else
  cerr << "using NTL bigints but no multiprecision floating point";
#endif
  cerr<<endl;
}

