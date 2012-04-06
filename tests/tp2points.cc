// tp2points.cc:  test program for P2Point class
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include "p2points.h"

int main()
{
  cout << "Test program for P2Point class" << endl;
  P2Point P, Q, R;
  cout << "Point input formats are [x:y:z], [x,y], [x/z,y/z] with any type of brackets" << endl; 
  cout << "Enter a point P: ";
  cin >> P;
  cout << "P="<<P<<endl;
  cout << "output_pari(P) = "; output_pari(cout,P); cout<<endl;
  cout << "P==P? "<< (P==P) <<endl;
  P=P2Point(2,4,-6);
  cout << "After P=P2Point(2,4,-6), P="<<P<<endl;
  Q=P;
  cout << "After Q=P, Q="<<Q<<endl;
  Q=transform(P,BIGINT(3),BIGINT(4),BIGINT(5),BIGINT(6),0);
  cout << "After Q=transform(P,3,4,5,6,0), Q="<<Q<<endl;
  R=transform(Q,BIGINT(3),BIGINT(4),BIGINT(5),BIGINT(6),1);
  cout << "After R=transform(Q,3,4,5,6,1), R="<<R<<endl;
  cout << "R==P? "<< (R==P) <<endl;
  bigint x,y,z; bigrational X,Y;  bigfloat rx,ry;
  P.getcoordinates(x,y,z);
  cout<<"Projective coordinates of P are X="<<x<<", Y="<<y<<", Z="<<z<<endl;
  P.getaffinecoordinates(X,Y);
  cout<<"Affine coordinates of P are X="<<X<<", Y="<<Y<<endl;
  P.getrealcoordinates(rx,ry);
  cout<<"Real affine coordinates of P are X="<<rx<<", Y="<<ry<<endl;
  cout<<"P.isintegral()? "<<P.isintegral()<<endl;
}

// end of file: tp2points.cc
