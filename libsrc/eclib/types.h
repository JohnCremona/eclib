// types.h: typedefs for scalar, vector, matrix types
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
 
#if     !defined(_ECLIB_TYPES_H)
#define _ECLIB_TYPES_H      1       //flags that this file has been included

#include <eclib/smatrix_elim.h>

// SCALAR_OPTION may be set to 1 or 2 or 3 by user

#if (SCALAR_OPTION==1) // scalar is int
typedef int scalar;
typedef vec_i vec;
typedef mat_i mat;
typedef subspace_i subspace;
typedef ssubspace_i ssubspace;
typedef svec_i svec;
typedef smat_i smat;
typedef smat_i_elim smat_elim;
#else
#if (SCALAR_OPTION==2) // scalar is long
typedef long scalar;
typedef vec_l vec;
typedef mat_l mat;
typedef subspace_l subspace;
typedef ssubspace_l ssubspace;
typedef svec_l svec;
typedef smat_l smat;
typedef smat_l_elim smat_elim;
#else
#if (SCALAR_OPTION==3) // scalar is bigint
typedef long scalar;
typedef vec_l vec;
typedef mat_l mat;
typedef subspace_l subspace;
typedef ssubspace_l ssubspace;
typedef svec_l svec;
typedef smat_l smat;
typedef smat_l_elim smat_elim;
#endif
#endif
#endif

#endif
