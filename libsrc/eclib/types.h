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
 
#include <eclib/smatrix_elim.h>

// SCALAR_OPTION may be set to 1 (int), 2 (long), or 3 (bigint) by user

#ifndef SCALAR_OPTION
#define SCALAR_OPTION 1 // int
#endif

#if (SCALAR_OPTION==2) // long
#if !defined(scalar_type_defined)
const string scalar_type = "long";
#define scalar_type_defined
#endif
typedef long scalar;
typedef vec_l vec;
typedef mat_l mat;
typedef subspace_l subspace;
typedef ssubspace_l ssubspace;
typedef svec_l svec;
typedef smat_l smat;
typedef smat_l_elim smat_elim;
typedef form_finder_l form_finder;
#define to_vec to_vec_l
#else
#if (SCALAR_OPTION==3) // bigint
#if !defined(scalar_type_defined)
const string scalar_type = "bigint";
#define scalar_type_defined
#endif
typedef bigint scalar;
typedef vec_m vec;
typedef mat_m mat;
typedef subspace_m subspace;
typedef ssubspace_m ssubspace;
typedef svec_m svec;
typedef smat_m smat;
typedef smat_m_elim smat_elim;
typedef form_finder_m form_finder;
#define to_vec to_vec_m
#else
#if (SCALAR_OPTION==1) // int
#if !defined(scalar_type_defined)
const string scalar_type = "int";
#define scalar_type_defined
#endif
typedef int scalar;
typedef vec_i vec;
typedef mat_i mat;
typedef subspace_i subspace;
typedef ssubspace_i ssubspace;
typedef svec_i svec;
typedef smat_i smat;
typedef smat_i_elim smat_elim;
typedef form_finder_i form_finder;
#define to_vec to_vec_i
#else
#if !defined(scalar_type_defined)
const string scalar_type = "undefined";
#define scalar_type_defined
#endif
#endif
#endif
#endif
