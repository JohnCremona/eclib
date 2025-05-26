/* parifact.cc: integer factorization using libpari, interface via strings */
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

#include <eclib/parifact.h>
#include <eclib/interface.h> // for getenv_with_default
#include <pari/pari.h>


//#define DEBUG_GPFACT

#include <iostream>

#define DEFAULT_PARI_SIZE 100000000
#define DEFAULT_PARI_MAX_PRIME 1000000

void eclib_pari_init(long max_prime=DEFAULT_PARI_MAX_PRIME)
{
  if (!avma) {
    long pari_size = strtol(getenv_with_default("PARI_SIZE", "DEFAULT_PARI_SIZE").c_str(), NULL, 0);
    if (pari_size==0) // e.g. syntax error in the environment variable PARI_SIZE
      pari_size = DEFAULT_PARI_SIZE;
#ifdef DEBUG_GPFACT
    std::cout<<"calling pari_init with pari_size = "<<pari_size<<endl;
#endif
    // the first parameter is the maximum stack size in bytes
    // the second parameter is the maximum precomputed prime
    pari_init(pari_size, max_prime);
  }
}

string
factor(const string n)
{
  eclib_pari_init();
  pari_sp av=avma;  // store pari stack pointer
#ifdef DEBUG_GPFACT
  std::cout<<"factor called with "<<n<<endl;
#endif
  GEN x = strtoi(n.c_str());
  setsigne(x,1);
  x = gel(Z_factor(x),1);
  settyp(x,t_VEC);

  string ans(GENtostr(x));
#ifdef DEBUG_GPFACT
  std::cout<<"factor returns "<<ans<<endl;
#endif
  avma=av;         // restore pari stackpointer
  return ans;
}

int 
is_prime(const string p)
{
  eclib_pari_init();
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

//#define DEBUG_ELLAP 1
long
ellap(long a1, long a2, long a3, long a4, long a6, long p)
{
  eclib_pari_init();
  pari_sp av=avma;  // store pari stack pointer
#ifdef DEBUG_ELLAP
  std::cout<<"ellap called with ["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"], p="<<p<<endl;
#endif
  GEN ai = mkvecn(5, stoi(a1), stoi(a2), stoi(a3), stoi(a4), stoi(a6));
  GEN pp = stoi(p);
  long ap = itos(ellap(ellinit(ai, pp, 0), pp));
#ifdef DEBUG_ELLAP
  std::cout<<"ellap returns "<<ap<<endl;
#endif
  avma=av;         // restore pari stackpointer
  return ap;
}
