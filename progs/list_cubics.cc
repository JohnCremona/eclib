// list_cubics.cc: Program for listing integer cubics with given or bounded discriminant
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2022 John Cremona
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

void output_cubics(const bigint& disc, int include_reducibles=1, int gl2=0, int verbose=0)
{
  vector<cubic> glist = reduced_cubics(disc, include_reducibles, gl2, verbose);
  if (glist.size()==0)
    {
      if(verbose>1)
        {
          cout<< "No ";
          if (!include_reducibles) cout << "irreducible";
          cout << " cubics with discriminant " << disc << endl;
        }
    }
  else
    {
      if (verbose)
        cout << glist.size() << " with discriminant ";
      cout << disc << " " << glist << endl;
    }
}

int main()
{
  initprimes("PRIMES");
  bigint disc, absdisc, maxdisc;
  int verbose=0, include_reducibles=1, gl2=0, single=0;

  cerr << "Verbosity level (0, 1 or 2): ";
  cin >> verbose;

  cerr << "Single discriminants (1) or a range (0)? ";
  cin >> single;

  if (single)
    {
      cerr << "Include reducible cubics? (0 or 1): ";
      cin >> include_reducibles;

      cerr << "Use GL(2,Z)-equivalence instead of SL(2,Z)? (0 or 1): ";
      cin >> gl2;

      while(cerr << "Enter discriminant (positive or negative, 0 to stop): ",	cin >> disc, !is_zero(disc))
        {
          if (verbose)
            {
              cout << (include_reducibles? "Cubics": "Irreducible cubics") << " with discriminant ";
              cout << disc;
              cout << " up to " << (gl2?"GL":"SL") << "(2,Z)-equivalence";
              cout << endl;
            }
          output_cubics(disc, include_reducibles, gl2, verbose);
        }
    }
  else
    {
      while(cerr << "Enter discriminant bound (positive or negative, 0 to stop): ",	cin >> maxdisc, !is_zero(maxdisc))
        {
          cerr << "Include reducible cubics? (0 or 1): ";
          cin >> include_reducibles;

          cerr << "Use GL(2,Z)-equivalence instead of SL(2,Z)? (0 or 1): ";
          cin >> gl2;

          int neg=(maxdisc<0);

          cout << (include_reducibles? "Cubics with ": "Irreducible cubics with ");
          cout << (neg? "negative discriminant down to ": "positive discriminant up to ");
          cout << maxdisc;
          cout << ", up to " << (gl2? "GL": "SL") << "(2,Z)-equivalence";
          cout << endl;

          if (neg)
            ::negate(maxdisc);
          for(absdisc=1; absdisc<=maxdisc; absdisc++)
	    {
	      disc=absdisc;
	      if(neg) ::negate(disc);
              output_cubics(disc, include_reducibles, gl2, verbose);
	    }
          cout<<endl;
        }
    }
}
