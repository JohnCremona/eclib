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
#include "interface.h"
#include "pari_init.h"

using NTL::ZZ;
using PARI::GEN;

// Warning: types fmpz_t and GEN are indistinguishable to the compiler

inline ZZ to_NTL (const int& a) {return to_ZZ(a);}    // from int to NTL integer
inline ZZ to_NTL (const long a) {return to_ZZ(a);}   // from long to NTL integer
ZZ FLINT_to_NTL (const fmpz_t& a); // from FLINT integer to NTL integer
ZZ PARI_to_NTL (const GEN& a);    // from PARI integer to NTL integer

fmpz_t* to_FLINT  (const int& a);  // from int to FLINT integer
fmpz_t* to_FLINT  (const long& a);  // from long to FLINT integer
fmpz_t* to_FLINT  (const ZZ& a);  // from NTL integer to FLINT integer
fmpz_t* to_FLINT  (const GEN& a); // from PARI integer to FLINT integer
fmpz_t* to_FLINT  (const INT& a); // from INT to FLINT integer

inline void set(int&a, const fmpz_t& z) {a = fmpz_get_si(z);}
inline void set(long&a, const fmpz_t& z) {a = fmpz_get_si(z);}
inline void set(ZZ&a, const fmpz_t& z) {a = FLINT_to_NTL(z);}
inline void set(INT&a, fmpz_t z) {a = INT(z);}

inline GEN to_PARI (const int& a) {return stoi(a);}  // from int to PARI integer
inline GEN to_PARI (const long& a) {return stoi(a);} // from long to PARI integer
GEN to_PARI (const ZZ& a);      // from NTL integer to PARI integer
GEN to_PARI (const fmpz_t& a);  // from FLINT integer to PARI integer

ZZ to_ZZ (const INT& a); // from INT to ZZ

inline INT to_INT(const int& x) {return INT(x);}
inline INT to_INT(const long& x) {return INT(x);}
inline INT to_INT(const ZZ& x) {return INT(*to_FLINT(x));}
inline INT to_INT(const GEN& x) {return INT(*to_FLINT(x));}

#endif

// end of file convert.h
