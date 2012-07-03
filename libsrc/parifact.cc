/* parifact.cc: integer factorization using libpari, interface via strings */
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////

#ifdef USE_PARI_FACTORING

#include <pari/pari.h>
#include <eclib/parifact.h>

//#define DEBUG_GPFACT
#include <iostream>
using namespace std;

string
factor(const string n)
{
  if (!bot) {
    pari_init(1000000, 1000000);
  }
#ifdef DEBUG_GPFACT
  std::cout<<"factor called with "<<n<<endl;
#endif

  pari_sp av=avma;  // store pari stack pointer
  GEN x = strtoi(n.c_str());
  setsigne(x,1);
  x = gel(Z_factor(x),1);
  settyp(x,t_VEC);
 
  string ans(GENtostr(x));
  avma=av;         // restore pari stackpointer
  return ans;
}

int 
is_prime(const string p)
{
  if (!bot) {
    pari_init(1000000, 1000000);
  }
  pari_sp av=avma;  // store pari stack pointer

#ifdef DEBUG_GPFACT
  std::cout<<"is_prime called with "<<p<<"..."<<flush;
#endif
  int ans = (isprime((GEN)strtoi(p.c_str()))==1);
#ifdef DEBUG_GPFACT
  std::cout<<"and returns "<<ans<<std::endl;
#endif
  
  avma=av;         // restore pari stackpointer
  return ans;
}

#endif
