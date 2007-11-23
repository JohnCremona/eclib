// method.h:  preprocessor definitions for linear algebra options
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
 
#if     !defined(_METHOD_H)
#define _METHOD_H      1       //flags that this file has been included

// Linear algebra options:

#ifndef METHOD     // So you can override the setting at compile time
#define METHOD 2
#endif
//=0 for standard int arithmetic
//=1 for long-long arithmetic (obsolete) 
//=2 for ints with modular method (usually best in practice)
//=3 for multi-length method (slower)
//=4 for multi-length modular (not really used)
//=5 for standard long arithmetic
//=6 for longs with modular method

#if (METHOD==2)||(METHOD==4)||(METHOD==6)
#define MODULAR   // Causes linear algebra to be done modulo global MODULUS
#endif

#if (METHOD==3)||(METHOD==4)
#define MULTI
#endif

//The next two cause scalar, vector, matrix to be defined properly:
#if (METHOD==0)||(METHOD==2)
#define SCALAR_OPTION 1
#endif
#if (METHOD==5)||(METHOD==6)
#define SCALAR_OPTION 2
#endif

#include "arith.h"   // defines BIGPRIME
#include "vector.h"
#include "matrix.h"
#include "subspace.h"
#include "smatrix_elim.h"

#ifdef MULTI
#define SCALAR bigint
#define VEC vec_m
#define MAT mat_m
#define SUBSP msubspace
#include "marith.h"
#include "mvector.h"
#include "mmatrix.h"
#include "msubspace.h"
#define MODULUS atoI(string("6074000003").c_str())  // will convert
						    // this string to
						    // a bigint
#else
#define MODULUS BIGPRIME  // used for modular linear algebra
#define SCALAR scalar
#define VEC vec
#define MAT mat
#define SUBSP subspace
#endif

#ifdef MODULAR
#define EIGENSPACE(a,b) peigenspace(a,b,MODULUS)
#define SUBEIGENSPACE(a,b,c) psubeigenspace(a,b,c,MODULUS)
#define COMBINE(a,b) pcombine(a,b,MODULUS)
#define RESTRICT(a,b) prestrict(a,b,MODULUS)
#else
#define EIGENSPACE(a,b) eigenspace(a,b)
#define SUBEIGENSPACE(a,b,c) subeigenspace(a,b,c)
#define COMBINE(a,b) combine(a,b)
#define RESTRICT(a,b) restrict(a,b)
#endif

#if (METHOD==0)
#define form_finder form_finder0
#else
#if (METHOD==1)
#define form_finder form_finder1
#else
#if (METHOD==2)
#define form_finder form_finder2
#else
#if (METHOD==3)
#define form_finder form_finder3
#else
#if (METHOD==4)
#define form_finder form_finder4
#else
#if (METHOD==5)
#define form_finder form_finder5
#else
#if (METHOD==6)
#define form_finder form_finder6
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#endif
