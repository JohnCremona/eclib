// TLATCONST.CC -- test of lattice constant procedures
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
//

#include <eclib/htconst.h>

int main(){
#ifdef MPFP
  set_precision(100);
#endif
  initprimes("PRIMES",0);

  cout << "Table of Gamma values\n";
  cout << "---------------------\n\n";

  long n;
  cout << "n\tGamma(n)\tGamma(n+1/2)\n";
  for (n=1; n<=20; n++)
    cout << n << "\t" << Gamma_n(n) << "\t" << Gamma_n_plus_half(n) << endl;

  cout << "\n\n";

  cout << "Table of Lattice constants\n";
  cout << "--------------------------\n\n";

  cout << "n\tgamma_n\tgamma_n^n\n";
  for (n=1; n<=20; n++)
    {
      bigfloat gam = lattice_const(n);
      cout << n << "\t" << gam << "\t" << power(gam,n) << endl;
    }
  return 0;
}// end main()

