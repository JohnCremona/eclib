// matrix.cc: manage implementation of integer matrix classes
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
