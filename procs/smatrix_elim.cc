// smatrix_elim.cc: implementation of class smat_elim
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
 
// Original version by Luiz Figueiredo

#include "arith.h"
#include "smatrix_elim.h"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim
#undef ssubspace

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i
#define svec svec_i
#define smat smat_i
#define smat_elim smat_i_elim
#define ssubspace ssubspace_i
#include "smat_elim.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim
#undef ssubspace

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l
#define svec svec_l
#define smat smat_l
#define smat_elim smat_l_elim
#define ssubspace ssubspace_l
#include "smat_elim.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim
#undef ssubspace
