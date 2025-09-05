// convert.h: declarations of integer conversion functions PARI/NTL/FLINT
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2025 John Cremona
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
 
#ifndef _ECLIB_CONVERT_H
#define _ECLIB_CONVERT_H

// Conversions between multipreciaion integers in NTL, FLINT, PARI:

#include <gmp.h>
#include <flint/flint.h> // must include this first to set __FLINT_VERSION
#include <flint/fmpz.h>
#include <NTL/ZZ.h>
#include "pari_init.h"

using NTL::ZZ;
using PARI::GEN;

ZZ FLINT_to_NTL  (const fmpz_t& a);   // from FLINT to NTL
ZZ PARI_to_NTL   (const GEN& a);      // from PARI to NTL
fmpz_t* NTL_to_FLINT  (const ZZ& a);  // from NTL to FLINT
fmpz_t* PARI_to_FLINT (const GEN& a); // from PARI to FLINT
GEN NTL_to_PARI   (const ZZ& a);      // from NTL to PARI
GEN FLINT_to_PARI (const fmpz_t& a);  // from FLINT to PARI

#endif

// end of file convert.h
