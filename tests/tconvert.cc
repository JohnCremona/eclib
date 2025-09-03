// tconvert.cc -- test program for integer conversion functions
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
 
#include <cassert>
#include <eclib/interface.h>
#include <eclib/convert.h>
#include <eclib/pari_init.h>

std::ostream& operator<<(std::ostream& s, const fmpz_t& a);

// Compiler cannot distinguish between fmpz_t (aka 'long int[1]') and GEN (aka 'long int*')
//std::ostream& operator<<(std::ostream& s, const GEN& a);

int test(const ZZ& a);

int main()
{
  eclib_pari_init();

  cout << "Testing conversions between types ZZ (NTL integers), fmpz_t (FLINT integers) and t_INT (PARI integers)" << endl;

  ZZ a;
  // cout << "Enter an integer a: "; cin >> a; cout << endl;
  // if (test(a))
  //   cout << "a = " << a << " PASSED" << endl << endl;
  // else
  //   cout << "a = " << a << " FAILED" << endl << endl;

  ZZ a0 = to_ZZ(0);
  ZZ a1 = to_ZZ("12345");
  ZZ a2 = to_ZZ("-98765");
  ZZ a3 = to_ZZ("-8472893749012374104710000000000000000000000000000000000001");
  ZZ a4 = to_ZZ("83475623875628564398568356325856876198561566179865781346578165843561854643856198564189564184651856148356184654187561856148561856143856431851");

  for (auto ai: {a0,a1,a2,a3,a4})
    if (test(ai))
      cout << "a = " << ai << " PASSED" << endl << endl;
    else
      cout << "a = " << ai << " FAILED" << endl << endl;
}

int test(const ZZ& a)
{
  cout << "Conversion tests" << endl;
  cout << "a = " << a << " (as ZZ)" << endl;

  fmpz_t* b = NTL_to_FLINT(a);
  cout << "a = " << b[0] << " (as fmpz_t)" << endl;
  ZZ a2 = FLINT_to_NTL(b[0]);
  cout << "converting back to ZZ...\na = " << a2 << endl;
  if (a==a2)
    cout << "OK" <<endl;
  else
    {
      cout << "******************* WRONG! ******************" <<endl;
      return 0;
    }
  GEN c = NTL_to_PARI(a);
  pari_printf("a = %Ps (as t_INT)\n", c);
  ZZ a3 = PARI_to_NTL(c);
  cout << "converting back to ZZ...\na = " << a3 << endl;
  if (a==a3)
    cout << "OK" <<endl;
  else
    {
      cout << "******************* WRONG! ******************" <<endl;
      return 0;
    }

  cout << "Converting from FLINT to PARI\n";
  GEN d = FLINT_to_PARI(b[0]);
  pari_printf("a = %Ps (as t_INT)\n", d);
  cout << "Converting from PARI back to FLINT\n";
  fmpz_t* e = PARI_to_FLINT(d);
  cout << "a = " << e[0] << " (as fmpz_t)" << endl;
  cout << "and then back from FLINT to NTL\n";
  ZZ f = FLINT_to_NTL(e[0]);
  cout << "a = " << f << " (as ZZ)" << endl;
  if (a==f)
    cout << "OK" <<endl;
  else
  {
      cout << "******************* WRONG! ******************" <<endl;
      return 0;
    }
  return 1;
}

std::ostream& operator<<(std::ostream& s, const fmpz_t& a)
{
  char* st = fmpz_get_str(NULL, 10, a);
  s << std::string(st);
  flint_free(st);
  return s;
}

// std::ostream& operator<<(std::ostream& s, const GEN& a)
// {
//   char* st = pari_sprintf("%Ps", a);
//   s << std::string(st);
//   free(st);
//   return s;
// }
