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

#if     !defined(_ECLIB_FLINTERFACE_H)
#define _ECLIB_FLINTERFACE_H      1       //flags that this file has been included

// This interface was written with help from Fredrik Johansson

// Currently only used to provide ref_via_flint() for Zmat<int> and
// Zmat<long> via FLINT's rref()

// FLINT support using 32-bit hmod_mat for scalar=int via Fredrik
// Stromberg's mini interface to gr_mat, or using 64-bit nmod_mat for
// scalar=long, or using fmpx_mod_mat_y for scalar=bigint.

// This wrapper interface provides the necessary hmod_mat functions

//#define TRACE_FLINT_RREF

#include <gmp.h>
#include <flint/flint.h> // must include this first to set __FLINT_VERSION

#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/profiler.h>

typedef unsigned int hlimb_t;

typedef struct
{
    hlimb_t * entries;
    slong r;
    slong c;
#if (__FLINT_VERSION==3)&&(__FLINT_VERSION_MINOR<3)
    hlimb_t ** rows;
#else
    slong stride;
#endif
    nmod_t mod;
}
hmod_mat_struct;

typedef hmod_mat_struct hmod_mat_t[1];

#if (__FLINT_VERSION==3)&&(__FLINT_VERSION_MINOR<3)
#define hmod_mat_entry(mat,i,j) ((mat)->rows[(i)][(j)])
#else
#define hmod_mat_entry nmod_mat_entry
#endif
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

#endif // _ECLIB_FLINTERFACE_H

// uscalar: either hlimb_t (unsigned int) or mp_limb_t (unsigned long)
