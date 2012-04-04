// matrix.cc: manage implementation of integer matrix classes
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
 
#include "arith.h"
#include "matrix.h"

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i
#include "mat.cc"
#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l
#include "mat.cc"
#undef scalar
#undef vec
#undef mat
#undef subspace
