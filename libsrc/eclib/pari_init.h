// pari_init.h: initialization of libpari
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

#ifndef _ECLIB_PARI_INIT_H
#define _ECLIB_PARI_INIT_H      1
                           //flags that this file has been included

#include <iostream>
#include <pari/pari.h>       // must be included after any NTL

#define DEFAULT_PARI_SIZE 100000000
#define DEFAULT_PARI_MAX_PRIME 1000000

void eclib_pari_init(long max_prime=DEFAULT_PARI_MAX_PRIME);

#endif
