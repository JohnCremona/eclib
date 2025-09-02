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
 
#include <eclib/interface.h>
#include <eclib/convert.h>
#include <cassert>


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

int test(const ZZ& a)
{
  cout << "Conversion tests" << endl;
  cout << "a = " << a << " (as ZZ)" << endl;

  fmpz_t* b = NTL_to_FLINT(a);
  char* st = fmpz_get_str(NULL, 10, b[0]);
  cout << "a = " << std::string(st) << " (as fmpz_t)" << endl;
  flint_free(st);
  ZZ a2 = FLINT_to_NTL(*b);
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
  char* st1 = fmpz_get_str(NULL, 10, e[0]);
  cout << "a = " << std::string(st1) << " (as fmpz_t)" << endl;
  flint_free(st1);
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

  a = to_ZZ("12345");
  if (test(a))
    cout << "a = " << a << " PASSED" << endl << endl;
  else
    cout << "a = " << a << " FAILED" << endl << endl;

  a = to_ZZ("-98765");
  if (test(a))
    cout << "a = " << a << " PASSED" << endl << endl;
  else
    cout << "a = " << a << " FAILED" << endl << endl;

  a = to_ZZ("-8472893749012374104710000000000000000000000000000000000001");
  if (test(a))
    cout << "a = " << a << " PASSED" << endl << endl;
  else
    cout << "a = " << a << " FAILED" << endl << endl;

  a = to_ZZ(0);
  if (test(a))
    cout << "a = " << a << " PASSED" << endl << endl;
  else
    cout << "a = " << a << " FAILED" << endl << endl;

  a = to_ZZ("83475623875628564398568356325856876198561566179865781346578165843561854643856198564189564184651856148356184654187561856148561856143856431851");
  if (test(a))
    cout << "a = " << a << " PASSED" << endl << endl;
  else
    cout << "a = " << a << " FAILED" << endl << endl;
}
