// cubics.cc: Program for listing integer cubics with given discriminant
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
 
#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>

int main()
{
  initprimes("PRIMES");
  bigint disc;
  int verbose=0, include_reducibles=1, gl2=0;
  cerr << "Verbosity level (0, 1 or 2): ";
  cin >> verbose;

  while(cerr << "Enter discriminant (positive or negative, 0 to stop): ", cin >> disc, !is_zero(disc))
    {
      cerr << "Include reducible cubics? (0 or 1): ";
      cin >> include_reducibles;
      cerr << "Use GL(2,Z)-equivalence instead of SL(2,Z)? (0 or 1): ";
      cin >> gl2;
      vector<cubic> glist = reduced_cubics(disc, include_reducibles, gl2, verbose);
      if (glist.size()==0)
        {
          cout<< "No ";
          if (!include_reducibles) cout << "irreducible";
          cout << " cubics with discriminant " << disc << endl;
        }
      else
        {
          cout << glist.size();
          if (!include_reducibles) cout << " irreducible";
          cout << " reduced cubics with discriminant " << disc;
          cout << " (up to " << (gl2?"GL":"SL") << "(2,Z)-equivalence)";
          if (glist.size()>0) cout<< " : " << glist;
          cout << endl;
        }
    }
}
