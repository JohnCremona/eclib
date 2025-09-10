// mspace.cc: test program for subspace class with bigint scalars
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
 
#include <eclib/linalg.h>

const bigint modulus(default_modulus<bigint>());

int main()
{
cin.flags( cin.flags() | ios::dec );

cout << "\nM-Subspace package test program\n\n";

int i,r=0;

while (cout << "Enter size of square matrix M: ", cin >> r, r>0 )
{
  mat_m m(r,r);
  cout << "Enter entries of M: ";
  cin >> m;
  cout << " M = \n" << m << endl;
  cout << "Trace(M) = " << m.trace() << endl;
  //
  mat_m mpower=m;
  for (i=2; i<=r; i++)
    {mpower=mpower*m;
     cout << "m^" << i << " = \n" << mpower << endl;
     cout << "Trace(m^" << i << ") = " << mpower.trace() << endl;
   }
//
  {
    vector<bigint> cp = m.charpoly();
    cout << "char. poly. of m has coefficients " << cp << endl;
  }
  cout << "det(M) = " << m.determinant() << endl;
  cout << "rank(M) = " << m.rank() << endl;
  cout << "nullity(M) = " << m.nullity() << endl;
  {
    subspace_m ker = kernel(m);
    mat_m kerbasis = basis(ker);
    cout << "kernel(m) has basis\n" << kerbasis << endl;
    vec_i kerpivs = pivots(ker);
    cout << "pivots: " << kerpivs << "\n";
    bigint kerdenom = denom(ker);
    cout << "denom:  " << kerdenom  << "\n";
  }
  {
    subspace_m im = image(m);
    cout << "image(m) has basis\n" << basis(im) << endl;
    cout << "pivots: " << pivots(im) << "\n";
    cout << "denom:  " << denom(im)  << "\n";
  }
  {
    bigint lambda;
    cout << "Enter lambda: "; cin >> lambda;
    subspace_m elambda = eigenspace(m,lambda);
    cout << "eigenspace for lambda = " << lambda << " has basis\n" << basis(elambda) << endl;
    cout << "with dimension " << dim(elambda) << endl;
    cout << "\nNow repeating eigenspace calculation modulo " << modulus << endl;
    subspace_m elp;
    lift(peigenspace(m,lambda,modulus),modulus, elp);
    cout << "eigenspace for lambda has basis\n" << basis(elp) << endl;
    cout << "with dimension " << dim(elp) << endl;
  }
}
 cout<<endl;
}
