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
namespace PARI{
  extern "C" {
#include <pari/pari.h>
  }
#undef coeff // since pari defines this as a macro with 3 args but NTL as a function with 2
#undef ulong // since pari defines this to be pari_ulong
}

using PARI::pari_sp;
using PARI::avma;
using PARI::itos;
using PARI::stoi;
using PARI::pari_printf;
using PARI::ellinit;
using PARI::ellap;
using PARI::Z_factor;
using PARI::nbrows;
using PARI::isprime;
using PARI::mkvecn;
using PARI::GEN;
using PARI::cgetg;
using PARI::itou;
using PARI::utoi;
using PARI::absi;
using PARI::negi;
using PARI::is_bigint;
using PARI::binary_2k;
using PARI::fromdigits_2k;
using PARI::pari_printf;
using PARI::t_VEC;
// NOT using PARI::signe, PARI::lg, PARI::gel which are macros

#define DEFAULT_PARI_SIZE 100000000
#define DEFAULT_PARI_MAX_PRIME 1000000

void eclib_pari_init(long max_prime=DEFAULT_PARI_MAX_PRIME);

#endif
