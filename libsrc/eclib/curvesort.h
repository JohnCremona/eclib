// FILE CURVESORT.H:  isogeny class id codes etc
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

#ifndef _ECLIB_CURVESORT_H
#define _ECLIB_CURVESORT_H      1
                           //flags that this file has been included

#include <string>
using namespace std;

#define USE_NEW_CODE_LETTERS


int booknumber0(int level, int form);  // permutes numbers starting from 0

const char alphabet[] = "abcdefghijklmnopqrstuvwxyz";

int booknumber(int level, int form);  // permutes numbers starting from 1

// new-new codes (from 01.08.05) are:
// a,b,...,z,ba,bb,...,bz,ca,cb,... etc., i.e. straight base 26 with
// digits a=0, b=1, ..., z=25

// Function to convert new code to integer (from 0) for any length of code.
// NB N=176400 has 516 < 26^2 newforms, with codes from a to vt!
int codeletter_to_int(string code);  // i counts from 0!

// Function to convert integer (from 0) to new code

string new_codeletter(int i);  // i counts from 0!

#ifdef USE_NEW_CODE_LETTERS
inline string codeletter(int i) {return new_codeletter(i);}
#else
inline string codeletter(int i) {return old_codeletter(i);}
#endif

#endif
