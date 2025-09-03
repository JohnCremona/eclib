// pari_interface.h: functions using libpari (integer factorization and ellap)
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

#ifndef _ECLIB_PARI_INTERFACE_H
#define _ECLIB_PARI_INTERFACE_H      1
                           //flags that this file has been included

#include "interface.h"

pair<vector<bigint>, vector<bigint>> factor_via_pari(const bigint& n);

int is_prime_via_pari(const bigint& p);

long ellap(long a1, long a2, long a3, long a4, long a6, long p);
bigint ellap(const bigint& a1, const bigint& a2, const bigint& a3, const bigint& a4, const bigint& a6, const bigint& p);

#if(0) // obsolete versions using string interface
#include <string>
string factor(const string n);
int is_prime(const string p);
#endif

#endif
