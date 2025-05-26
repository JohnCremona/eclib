// flinterface.h: used to provide interface to flint's linear algebra mod p
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

// This interface was written with help from Fredrik Johansson

// This file will be included either with SCALAR_OPTION=1 or with SCALAR_OPTION=1

#if FLINT

// Possible FLINT support is governed by the macros FLINT and __FLINT_VERSION:
//
// FLINT==0: no FLINT support (or a version <2.3)
// FLINT==1: FLINT support using 64-bit nmod_mat (if __FLINT_VERSION<3)
//                     for both SCALAR_OPTION=1 (scalar=in) and 2 (scalar=long)
//                 support using 32-bit hmod_mat (if __FLINT_VERSION>=3) when scalar=int
//                         via Fredrik Stromberg's mini interface to gr_mat

//#define TRACE_FLINT_RREF

#include <gmp.h>
#include <flint/flint.h> // must include this first to set __FLINT_VERSION

#if (__FLINT_VERSION>2)
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/nmod.h>
#endif

#include <flint/nmod_mat.h>
#include <flint/profiler.h>

#if (__FLINT_VERSION>2)&&(SCALAR_OPTION==1) // using 32-bit int scalars

typedef unsigned int hlimb_t;

typedef struct
{
    hlimb_t * entries;
    slong r;
    slong c;
    hlimb_t ** rows;
    nmod_t mod;
}
hmod_mat_struct;

typedef hmod_mat_struct hmod_mat_t[1];

#define hmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])
#define hmod_mat_nrows(mat) ((mat)->r)
#define hmod_mat_ncols(mat) ((mat)->c)

void
hmod_mat_init(hmod_mat_t mat, slong rows, slong cols, hlimb_t n);

void
hmod_mat_clear(hmod_mat_t mat);

void
hmod_mat_mul(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B);

slong
hmod_mat_rref(hmod_mat_t mat);

#undef uscalar
#undef mod_mat
#undef mod_mat_init
#undef mod_mat_clear
#undef mod_mat_entry
#undef mod_mat_nrows
#undef mod_mat_ncols
#undef mod_mat_rref
#undef mod_mat_mul
#define uscalar hlimb_t // unsigned int
#define mod_mat hmod_mat_t // uses unsigned ints
#define mod_mat_init hmod_mat_init
#define mod_mat_clear hmod_mat_clear
#define mod_mat_entry hmod_mat_entry
#define mod_mat_nrows hmod_mat_nrows
#define mod_mat_ncols hmod_mat_ncols
#define mod_mat_rref hmod_mat_rref
#define mod_mat_mul hmod_mat_mul

#else // __FLINT_VERSION<3 or SCALAR_OPTION=2 -- using 64-bit int scalars

#undef uscalar
#undef mod_mat
#undef mod_mat_init
#undef mod_mat_clear
#undef mod_mat_entry
#undef mod_mat_nrows
#undef mod_mat_ncols
#undef mod_mat_rref
#undef mod_mat_mul
#define uscalar mp_limb_t // unsigned long
#define mod_mat nmod_mat_t // uses unsigned longs
#define mod_mat_init nmod_mat_init
#define mod_mat_clear nmod_mat_clear
#define mod_mat_entry nmod_mat_entry
#define mod_mat_nrows nmod_mat_nrows
#define mod_mat_ncols nmod_mat_ncols
#define mod_mat_rref nmod_mat_rref
#define mod_mat_mul nmod_mat_mul

#endif // __FLINT_VERSION >=3 or <3

#endif // FLINT
