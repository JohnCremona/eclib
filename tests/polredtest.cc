// FILE POLREDTEST.CC: test program for functions for reducing ZZX
//                     polynomials (polredabs/polredbest) via libpari
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

#include "eclib/polred.h"

int main()
{
  cout << "Program polredtest: conversions between ZZX and t_POL and reduction of polynomials."
       << endl << endl;

  eclib_pari_init();
  int d;
  ZZ den;
  while (1)
    {
      cerr << "Enter degree d: ";
      cin >> d;
      if (d<1) exit(0);
      cerr << "Enter " << d+1 << " integers, starting with the leading coefficient: ";
      ZZX f;
      ZZ c;
      for (int i=0; i<=d; i++)
        {
          cin >> c;
          SetCoeff(f, d-i, c);
        }
      cout << "f = " << str(f) << endl;
      GEN P = ZZX_to_t_POL(f);
      pari_printf(" -as a t_POL: %Ps\n", P);
      ZZX g = t_POL_to_ZZX(P, den);
      assert (den==1);
      cout << " -back to ZZX: " << str(g) << endl;
      if (f==g)
        cout << " OK " << endl;
      else
        cout << " *** WRONG *** " << endl;

      if (IsIrreducible(f))
        {
          for (int canonical: {0,1})
            {
              if (canonical)
                cout << "Applying polredabs..." << endl;
              else
                cout << "Applying polredbest..." << endl;
              ZZX h;
              g = polred(f, h, den, canonical);
              cout << "... reduced polynomial is g = " << str(g);
              if (f==g)
                cout << " -- no change"<< endl;
              else
                cout << " -- polynomial has been reduced"<< endl;
              if (den==1)
                cout << "A root of f is a = " << str(h, "b") << " where g(b)=0" << endl;
              else
                cout << "A root of f is a = (" << str(h, "b") << ") / " << den << " where g(b)=0" << endl;
            }
        }
      else
        {
          cout << "f is reducible:\n";
          display_factors(f);
        }
      cout << endl;
    }
  exit(0);
}
