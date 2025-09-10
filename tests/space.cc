// space.cc: Test of subspace package
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

const scalar modulus(default_modulus<scalar>());

int main()
{
cout << "Subspace test program with scalar type " << scalar_type << "\n\n";

int times,ntimes,r=0;

while (cout << "Enter size of square matrix M: ", cin >> r, r>0 )
{
  mat m(r,r);
  cout << "Enter entries of M: ";
  cin >> m;
  cout << " M = \n" << m << endl;
  cout << "Trace(M) = " << m.trace() << endl;
/*
  int i;
  mat mpower=m;
  for (i=2; i<=r; i++)
    {mpower=mpower*m;
     cout << "m^" << i << " = \n" << mpower << endl;
     cout << "Trace(m^" << i << ") = " << mpower.trace() << endl;
   }
  {
    vector<long> cp = m.charpoly();
    cout << "char. poly. of m has coefficients " << cp << endl;
  }
  cout << "det(M) = " << m.determinant() << endl;
  cout << "rank(M) = " << m.rank() << endl;
  cout << "nullity(M) = " << m.nullity() << endl;
*/
  cout << endl << "Enter number of times to repeat kernel tests: ";
  cin >> ntimes;
  {
    subspace ker;
    times=ntimes;
    while(times--) ker = kernel(m);
    mat kerbasis = basis(ker);
    cout << "kernel(m) has basis\n" << kerbasis << endl;
    vec_i kerpivs = pivots(ker);
    cout << "pivots: " << kerpivs << "\n";
    scalar kerdenom = denom(ker);
    cout << "denom:  " << kerdenom  << "\n";
  }
  scalar p = modulus;
  cout << "Now compute kernel mod p, p = " << p << endl;
  {
    times=ntimes; subspace ker;
    while(times--) ker = pkernel(m, p);
    mat kerbasis = basis(ker);
    cout << "kernel(m) has basis\n" << kerbasis << endl;
    vec_i kerpivs = pivots(ker);
    cout << "pivots: " << kerpivs << "\n";
    subspace oldker;
    times=ntimes;
    while(times--) oldker = oldpkernel(m, p);
    if((kerbasis!=basis(oldker))||(kerpivs!=pivots(oldker)))
      {
	cout << "!!! Differs from old version !!!" << endl;
	cout << "Old basis = \n"<<basis(oldker)<<endl;
      }
    else
      cout << "!!! Agrees  with old version !!!" << endl;
  }
/*
  {
    subspace im = image(m);
    cout << "image(m) has basis\n" << basis(im);
    cout << "pivots: " << pivots(im) << "\n";
    cout << "denom:  " << denom(im)  << "\n";
  }
  {
    scalar lambda;
    cout << "Enter lambda: "; cin >> lambda;
    subspace elambda = eigenspace(m,lambda);
    cout << "eigenspace for lambda has basis\n" << basis(elambda);
    cout << "with dimension " << dim(elambda) << endl;
  }
*/
}
 cout<<endl;
}
