// method.h:  preprocessor definitions for linear algebra options
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
 
#if     !defined(_ECLIB_METHOD_H)
#define _ECLIB_METHOD_H      1       //flags that this file has been included

// Linear algebra options:  SCALAR_OPTION is 1 (int), 2 (long), or 3 (bigint)

#ifndef SCALAR_OPTION    // So you can override the setting at compile time
#define SCALAR_OPTION 1  // int
//#define SCALAR_OPTION 2  // long
//#define SCALAR_OPTION 3  // bigint
#endif

// types.h presets scalar, vec, mat, subspace, ssubspace, svec, smat, smat_elim
// to be int/long/bigint and *_i/*_l/*_m according to SCALAR
#include "types.h"

#define MODULUS DEFAULT_MODULUS  // (set in xmod.h) used for modular linear algebra

#endif
