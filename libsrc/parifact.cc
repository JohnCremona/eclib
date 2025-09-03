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

#include <iostream>

#include "eclib/interface.h"
#include "eclib/parifact.h"
#include "eclib/pari_init.h"
#include "eclib/convert.h"

//#define DEBUG_GPFACT

vector<bigint> factor_via_pari(const bigint& n)
{
  eclib_pari_init();
  pari_sp av=avma;  // store pari stack pointer

  vector<bigint> plist;
  if (n<2) return plist;
#ifdef DEBUG_GPFACT
  cout << "In factor_via_pari("<<n<<")\n";
#endif
  GEN pn = NTL_to_PARI(n);
#ifdef DEBUG_GPFACT
  pari_printf(" - n as t_INT: %Ps\n", pn);
#endif
  pn = gel(Z_factor(pn), 1); // a t_VEC of primes
#ifdef DEBUG_GPFACT
  pari_printf(" - plist as t_VEC: %Ps\n", pn);
#endif
  int np = lg(pn)-1;
  plist.resize(np);
  for(int i=0; i<np; i++)
    plist[i] = PARI_to_NTL(gel(pn, i+1));
#ifdef DEBUG_GPFACT
  cout << " - plist as vector<bigint>: " << plist << endl;
#endif

  avma=av;         // restore pari stackpointer
  return plist;
}

int is_prime_via_pari(const bigint& p)
{
  eclib_pari_init();
  pari_sp av=avma;  // store pari stack pointer

  int ans = isprime(NTL_to_PARI(p));

  avma=av;         // restore pari stackpointer
  return ans;
}

//#define DEBUG_ELLAP

bigint
ellap(const bigint& a1, const bigint& a2, const bigint& a3, const bigint& a4, const bigint& a6, const bigint& p)
{
  eclib_pari_init();
  pari_sp av=avma;  // store pari stack pointer
#ifdef DEBUG_ELLAP
  std::cout<<"ellap called with ["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"], p="<<p<<endl;
#endif
  GEN ai = mkvecn(5, NTL_to_PARI(a1), NTL_to_PARI(a2), NTL_to_PARI(a3), NTL_to_PARI(a4), NTL_to_PARI(a6));
  GEN pp = NTL_to_PARI(p);
  bigint ap = PARI_to_NTL(ellap(ellinit(ai, pp, 0), pp));
#ifdef DEBUG_ELLAP
  std::cout<<"ellap returns "<<ap<<endl;
#endif
  avma=av;         // restore pari stackpointer
  return ap;
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
