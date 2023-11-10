// tdivpol.cc -- test for division poly functions in divpol.h/cc 
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

#include <eclib/curve.h>

#include <eclib/polys.h>
#include <eclib/divpol.h>
#include <eclib/points.h>

int main()
{
  initprimes("PRIMES",0);

  Curve E(BIGINT(0),BIGINT(0),BIGINT(1),BIGINT(-7),BIGINT(6));

  Curvedata C(E, 0);
  cout << "Curve " << E << endl;

  int i;
  cout << "\nDivision polynomials:\n\n";
  for(i=2; i<12; i+=1)
    cout<<"Division Poly ("<<i<<") = \t" << division_polynomial(&C,i) << endl;

  cout << "\nMultiplication-by-n maps (x coordinate):\n\n";

  bigint a1,a2,a3,a4,a6;
  E.getai(a1,a2,a3,a4,a6);
  vector<int> nn = {2,3,5,7};
  for(vector<int>::iterator n = nn.begin(); n!=nn.end(); n++)
    {
      i = *n;
      cout << "n = " << i << ":\n";
      cout<<"Numerator:  \t" << mul_by_n_num(a1,a2,a3,a4,a6,i) << endl;
      cout<<"Denominator:\t" << mul_by_n_den(a1,a2,a3,a4,a6,i) << endl;
    }

#if(1)
  cout << "\nTesting division of points:\n";

  Point P0(C, BIGINT(0),BIGINT(2)) ;
  Point P1(C, BIGINT(1),BIGINT(0)) ;
  Point P2(C, BIGINT(2),BIGINT(0)) ;

  vector<Point> Plist = {P0, P1, P2, P0+P1, P0+P2, P1+P2, P0+P1+P2};
  vector<int> mlist = {2,3,5,7,11,13};
  for (vector<Point>::iterator Pi = Plist.begin(); Pi!=Plist.end(); Pi++)
    {
      Point P = *Pi;
      cout << "\nP = " << P << endl;
      for (vector<int>::iterator mi = mlist.begin(); mi!=mlist.end(); mi++)
        {
          int m = *mi;
          Point Q = m*P;
          vector<Point> newP = Q.division_points(m);
          cout << m << "*P = " << Q << ", divided by " << m << " gives back ";
          if (newP.size()==1 && newP[0]==P)
            {
              cout << "P, OK" << endl;
            }
          else
            {
              cout << newP << "???" << endl;
            }
        }
    }
#endif
}


//end of file tdivpol.cc





