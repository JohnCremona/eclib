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
 
// Conversions between multipreciaion integers in NTL, FLINT, PARI:

//#define DEBUG_CONVERT

#ifdef DEBUG_CONVERT
#include <iostream>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
#endif

#include "eclib/convert.h"
#include <NTL/ZZ_limbs.h>

using namespace NTL;

#if FLINT_LONG_LONG
typedef unsigned long long int ulong;
typedef long long int slong;
#else
typedef unsigned long int ulong;
typedef long int slong;
#endif

ZZ FLINT_to_NTL(const fmpz_t& a)  // from FLINT to NTL
{
  if (fmpz_fits_si(a))
    return ZZ(fmpz_get_si(a));

  // deal with sign (a is certainly nonzero now)
  int sign_a = fmpz_sgn(a);
  fmpz_t aa; fmpz_init(aa);
  fmpz_abs(aa, a);

  // extract limbs
  slong n = fmpz_size(aa);
  ulong* limbs = new ulong[n];
  fmpz_get_ui_array(limbs, n, aa);

  // construct ZZ from limbs
  ZZ b;
  ZZ_limbs_set(b, (ZZ_limb_t*)limbs, n);

  // return with sign
  return (sign_a<0? -b : b);
}

ZZ PARI_to_NTL(const GEN& a)  // from PARI to NTL
{
  // if a fits in a long int it is easy:
  if (!is_bigint(a))
    return ZZ(itos(a));

  // deal with sign (a is certainly nonzero now)
  int sign_a = signe(a);
  GEN aa = absi(a);

#ifdef DEBUG_CONVERT
  pari_printf("In NTL_to_PARI() with t_INT a =  %Ps, abs(a) = %Ps\n", a, aa);
#endif
  // extract digits as a t_VEC
  GEN digits = binary_2k(a, NTL_BITS_PER_LIMB_T);
#ifdef DEBUG_CONVERT
  cout << "NTL_BITS_PER_LIMB_T = " << NTL_BITS_PER_LIMB_T << endl;
  pari_printf("digits:  %Ps\n", digits);
#endif

  // convert this to an array of ZZ_limb_t
  int n = lg(digits)-1;
  ZZ_limb_t* limbs = new ZZ_limb_t[n];
  for (int i=1; i<lg(digits); i++)
    {
      limbs[n-i] = itou(gel(digits,i));
    }
#ifdef DEBUG_CONVERT
  cout << "Limbs for NTL: ";
  for (int i=0; i<n; i++)
    {
      if (i) cout << ", ";
      cout << limbs[i];
    }
  cout << endl;
#endif

  // construct ZZ from limbs
  ZZ b;
  ZZ_limbs_set(b, limbs, n);

  // return with sign
  return (sign_a<0? -b : b);
}

fmpz_t* NTL_to_FLINT(const ZZ& a)  // from NTL to FLINT
{
  if (IsZero(a))
    {
      fmpz_t* b = new fmpz_t[1];
      fmpz_init_set_si(b[0], 0);
      return b;
    }

  // deal with sign (now a is nonzero)
  int sign_a = sign(a);
  ZZ aa = abs(a);

#ifdef DEBUG_CONVERT
  cout<<"In NTL_to_FLINT() with ZZ a = " << a << ", abs(a) = " << aa <<endl;
#endif
  // extract limbs
  int n = aa.size();
  const ZZ_limb_t* limbs = ZZ_limbs_get(aa);
#ifdef DEBUG_CONVERT
  cout<<" #limbs = " << n << "\nlimbs: ";
  for(int i=0; i<n; i++)
    {
      if (i) cout << ", ";
      cout << limbs[i] << " ";
    }
  cout << endl;
#endif

  // construct an fmpz_t from these limbs:
  fmpz_t* b = new fmpz_t[1]; fmpz_init(b[0]);
  fmpz_set_ui_array(b[0], limbs, n);

  // adjust sign if necessary
  if (sign_a<0)
    fmpz_neg(b[0],b[0]);
  return b;
}

fmpz_t* PARI_to_FLINT(const GEN& a)  // from PARI to FLINT
{
  // if a fits in a long int it is easy:
  if (!is_bigint(a))
    {
      fmpz_t* b = new fmpz_t[1];
      fmpz_set_si(b[0], itos(a));
      return b;
    }
  // deal with sign (a is certainly nonzero here)
  int sign_a = signe(a);
  GEN aa = absi(a);

  // extract digits as a t_VEC
  GEN digits = binary_2k(aa, FLINT_BITS); // n's digits in base 2^64, as a t_VEC

  // convert this to an array of ulongs
  int n = lg(digits)-1;
  ulong* limbs = new ulong[n];
  for (int i=1; i<=n; i++)
    {
      limbs[n-i] = itou(gel(digits,i));
    }

  // construct an fmpz_t from these limbs:
  fmpz_t* b = new fmpz_t[1]; fmpz_init(b[0]);
  fmpz_set_ui_array(b[0], limbs, n);

  // adjust sign if necessary
  if (sign_a<0)
    fmpz_neg(b[0],b[0]);
  return b;
}

GEN NTL_to_PARI(const ZZ& a)  // from NTL to PARI
{
  if (IsZero(a))
    return stoi(0);

  // deal with sign (now a is nonzero)
  int sign_a = sign(a);
  ZZ aa = abs(a);

#ifdef DEBUG_CONVERT
  cout<<"In NTL_to_PARI() with ZZ a = " << a << ", abs(a) = " << aa <<endl;
#endif
  // extract limbs
  int n = aa.size();
  const ZZ_limb_t* limbs = ZZ_limbs_get(aa);

#ifdef DEBUG_CONVERT
  cout<<" #limbs = " << n << "\nlimbs:         ";
  for(int i=0; i<n; i++)
    {
      if (i) cout << ", ";
      cout << limbs[i] << " ";
    }
  cout << endl;
#endif

  // convert limbs array to a t_VEC
  GEN v = cgetg(n+1, t_VEC);
  for (int i=1; i<=n; i++)
    gel(v,i) = utoi(limbs[n-i]);

#ifdef DEBUG_CONVERT
    pari_printf("limbs as t_VEC: %Ps\n", v);
    cout << "NTL_BITS_PER_LIMB_T = " << NTL_BITS_PER_LIMB_T << endl;
#endif

  // construct t_INT from limbs
  GEN b = fromdigits_2k(v, NTL_BITS_PER_LIMB_T);
#ifdef DEBUG_CONVERT
  pari_printf("as t_INT: %Ps\n", b);
#endif

  // return with sign
  return (sign_a<0? negi(b): b);
}

GEN FLINT_to_PARI(const fmpz_t& a)  // from FLINT to PARI
{
  if (fmpz_fits_si(a))
    return stoi(fmpz_get_si(a));

  // deal with sign (a is certainly nonzero now)
  int sign_a = fmpz_sgn(a);

  fmpz_t aa; fmpz_init(aa);
  fmpz_abs(aa, a);

  // extract limbs
  slong n = fmpz_size(aa);
  ulong* limbs = new ulong[n];
  fmpz_get_ui_array(limbs, n, aa);

  // convert limbs array to a t_VEC
  GEN v = cgetg(n+1, t_VEC);
  for (int i=1; i<=n; i++)
    {
      gel(v,i) = utoi(limbs[n-i]);
    }

  // construct t_INT from limbs
  GEN b = fromdigits_2k(v, FLINT_BITS);

  // return with sign
  return (sign_a<0? negi(b): b);
}

