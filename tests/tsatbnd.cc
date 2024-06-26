// tsatbnd.cc -- test saturation index bound
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

#include <eclib/interface.h>
#include <eclib/method.h>
#include <eclib/curve.h>
#include <eclib/getcurve.h>
#include <eclib/points.h>
#include <eclib/cperiods.h>
#include <eclib/polys.h>
#include <eclib/curvemod.h>
#include <eclib/pointsmod.h>
#include <eclib/htconst.h>
#include <eclib/saturate.h>
//#include <eclib/elog.h>
//#include <eclib/sieve_search.h>
//#include <eclib/mwprocs.h>


int main()
{
  //  set_precision("Enter precision in bits");
  set_precision(100);
  initprimes("PRIMES",0);
  int verbose = 1;
  cerr<<"verbose (0/1)? ";             cin >>verbose;

  Curvedata C;

  while (getcurve(C, verbose))
    {

      cout << "======================================================\n\n";
      cout << "E = " << (Curve)C <<endl;

      Point P(C);
      int j=0, npts;
      cerr<<"enter number of points: ";      cin >> npts;
      vector<Point> points; points.reserve(npts);

      while(j<npts)
        {
          cerr<<"\n  enter point "<<(j+1)<<" : ";
          cin >> P;
          if ( P.isvalid() )
            {
              points.push_back(P);
              j++;
            }
          else
            {
              cerr<<"point "<<P<<" not on curve.\n\n";
            }
        }
      cout<<npts<<" points entered:\n";
      for ( auto& Pi : points)
        cout << Pi << " (height "<<height(Pi)<<")"<<endl;

      bigfloat reg = regulator(points);
      cout << "Regulator of input points: " << reg << endl;
      bigint egri = egr_index(points);
      cout << "Subgroup of EGR points has index " << egri << endl;

      cout<< "\n Images of points in component groups at bad primes and R\n\n";
      cout<< "p   Group   Image(s)" <<endl;
      bigint u, r, s, t;
      Curvedata Cmin = C.minimalize(u,r,s,t);
      vector<Point> points_min(points.size());
      std::transform(points.begin(), points.end(), points_min.begin(),
                     [&Cmin, u, r, s, t] (const Point& Pi) {return transform(Pi, &Cmin, u, r, s, t);});
      CurveRed CR(Cmin);
      vector<bigint> plist = getbad_primes(CR);
      ComponentGroups CG(CR);
      vector<vector<int> >  ims;
      for( const auto& p : plist)
        {
          cout<<p<<"  "<<CG.ComponentGroup(p)<<"  ";
          ims = MapPointsToComponentGroup(CG, points_min,  p);
          for( const auto& im : ims)
            cout<<im;
          cout<<endl;
        }
      bigint p(0);
      cout<<"R "<< CG.ComponentGroup(p)<<"  ";
      ims = MapPointsToComponentGroup(CG, points_min,  p);
      for( const auto&  im : ims)
        cout<<im;
      cout<<endl<<endl;

      bigint tam_prod = global_Tamagawa_number(CR, 1); // include real place
      bigint tam_exp = global_Tamagawa_exponent(CR, 1); //
      cout << "Global Tamagawa number:   " << tam_prod <<endl;
      cout << "Global Tamagawa exponent: " << tam_exp <<endl;

      bigfloat lambda_egr = lower_height_bound(CR, 1);
      cout << "lower height bound (egr points): "<< lambda_egr <<endl;
      bigint egr_bound = index_bound(points, 1, verbose);
      cout << "bound on saturation index (egr points): "<< egr_bound <<endl;

      bigfloat lambda_all = lower_height_bound(CR, 0);
      cout << "lower height bound (all points): "<< lambda_all <<endl;
      bigint all_bound = index_bound(points, 0, verbose);
      cout << "bound on saturation index (all points): "<< all_bound <<endl;
    }
}

//end of file tsatbnd.cc





