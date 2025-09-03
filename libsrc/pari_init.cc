// pari_init.cc: initialization of libpari
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

#include <iostream>

#include "eclib/interface.h"
#include "eclib/pari_init.h"       // must be included after any NTL

void eclib_pari_init(long max_prime)
{
  if (!avma) {
    long pari_size = strtol(getenv_with_default("PARI_SIZE", "DEFAULT_PARI_SIZE").c_str(), NULL, 0);
    if (pari_size==0) // e.g. syntax error in the environment variable PARI_SIZE
      pari_size = DEFAULT_PARI_SIZE;
    // the first parameter is the maximum stack size in bytes
    // the second parameter is the maximum precomputed prime
    pari_init(pari_size, max_prime);
  }
}
