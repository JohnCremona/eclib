// mspace.cc: test program for msubspace class
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
 
#include "msubspace.h"
const bigint MBIGPRIME=atoI("6074000003");
// will convert this string to an bigint
//This is nearly the largest p such that (p/2)^2 < 2^63.

int main()
{
cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)

cout << "\nM-Subspace package test program\n\n";

int i,r=0;

while (cout << "Enter size of square matrix M: ", cin >> r, r>0 )
{
  mat_m m(r,r);
  cout << "Enter entries of M: ";
  cin >> m;
  cout << " M = " << m;
  cout << "Trace(M) = " << trace(m) << endl;
  //
  mat_m mpower=m;
  for (i=2; i<=r; i++)
    {mpower=mpower*m;
     cout << "m^" << i << " = " << mpower;
     cout << "Trace(m^" << i << ") = " << trace(mpower) << endl;
   }
//
  {
    vector<bigint> cp = charpoly(m);
    cout << "char. poly. of m has coefficients " << cp << endl;
  }
  cout << "det(M) = " << determinant(m) << endl;
  cout << "rank(M) = " << rank(m) << endl;
  cout << "nullity(M) = " << nullity(m) << endl;
  {
    msubspace ker = kernel(m);
    mat_m kerbasis = basis(ker);
    cout << "kernel(m) has basis\n" << kerbasis;
    vec_i kerpivs = pivots(ker);
    cout << "pivots: " << kerpivs << "\n";
    bigint kerdenom = denom(ker);
    cout << "denom:  " << kerdenom  << "\n";
  }
  {
    msubspace im = image(m);
    cout << "image(m) has basis\n" << basis(im);
    cout << "pivots: " << pivots(im) << "\n";
    cout << "denom:  " << denom(im)  << "\n";
  }
  {
    bigint lambda;
    cout << "Enter lambda: "; cin >> lambda;
    msubspace elambda = eigenspace(m,lambda);
    cout << "eigenspace for lambda = " << lambda << " has basis\n" << basis(elambda);
    cout << "with dimension " << dim(elambda) << endl;
    cout << "\nNow repeating eigenspace calculation modulo " << MBIGPRIME << endl;
    msubspace elp = lift(peigenspace(m,lambda,MBIGPRIME),MBIGPRIME);
    cout << "eigenspace for lambda has basis\n" << basis(elp);
    cout << "with dimension " << dim(elp) << endl;
  }
}
 cout<<endl;
abort();

}
