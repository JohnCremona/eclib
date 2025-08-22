// mattest.cc: Matrix package test program
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
 
#include <iostream>
#include <eclib/interface.h>
#include <eclib/timer.h>
#include <eclib/arith.h>
#include <eclib/types.h>

int main(void)
{
  init_time();
  start_time();
  cout << "Matrix test program with scalar type " << scalar_type << ".\n\n";

  long i; int r;
  vec_i pc(1),npc(1);
  vec poly(1);
  scalar two(2), three(3);
  cout << "Enter size of a square matrix A: "; cin >> r;
  mat a(r,r);
  cout << "Enter entries of A: "; cin >> a;
  cout << "A = " << a;
  cout << "Using A.output(cout): ";  a.output(cout);
  cout << "Using A.output_pari(cout): ";  a.output_pari(cout);
  cout << "Using A.output_pretty(cout): \n";  a.output_pretty(cout);

  cout << "Enter any number "; cin >> i;
  cout << "Creating an array of 3 matrices\n";
  vector<mat> matlist(3);
  matlist[0] = a;
  matlist[1] = two*a;
  matlist[2] = three*a;
  cout << " A=" << matlist[0];
  cout << "2A=" << matlist[1];
  cout << "3A=" << matlist[2];

  for (i=1; i<=r; i++)
    cout << "row(A,"<<i<<") = " << a.row(i) << endl;
  cout << "A = " << a;
  for (int j=1; j<=r; j++)
    cout << "col(A,"<<j<<") = " << a.col(j) << endl;
  cout << "A = " << a;
  cout << "directsum(A,A) = " << directsum(a,a);
  cout << "Enter any number "; cin >> i;

  mat b = a;
  cout << "B = A = " << b;
  cout << "Enter any number "; cin >> i;
  cout << "B==A?" << (b==a) << endl;
  cout << "B!=A?" << (b!=a) << endl;
  b+=a;
  cout << "after B+:=A, A = " << a << "and B = " << b;
  cout << "Enter any number "; cin >> i;
  b-=a;
  cout << "after B-:=A, A = " << a << "and B = " << b;
  cout << "Enter any number "; cin >> i;
  b*=two;
  cout << "after B*:=2, A = " << a << "and B = " << b;
  cout << "Enter any number "; cin >> i;
  b/=two;
  cout << "after B/:=2, A = " << a << "and B = " << b;
  cout << "Enter any number "; cin >> i;
  cout << "A+B=" << (a+b);
  cout << "Now A = " << a << "and B = " << b;
  cout << "Enter any number "; cin >> i;
  cout << "A-B=" << (a-b);
  cout << "Now A = " << a << "and B = " << b;
  cout << "Enter any number "; cin >> i;
  cout << "A*B=" << (a*b);
  cout << "Now A = " << a << "and B = " << b;
  cout << "Enter any number "; cin >> i;
  cout << "-A=" << (-a);
  cout << "Now A = " << a;
  cout << "-A=" << (-a);
  cout << "Now A = " << a;
  cout << "Enter any number "; cin >> i;
  vector<scalar> cp = a.charpoly();
  cout << "char. poly. of A has coefficients " << cp << endl;
  cout << "det(A) = " << a.determinant() << endl;
  mat aug = colcat(a,mat::identity_matrix(r));
  cout << "Augmented matrix = " << aug << endl;

  long rk, ny;
  scalar denom;
  int method=0;
  cout << "Which echelon method? (0=standard,1=longlong,2=modular) ";
  cin>>method;
  cout << "\nUsing method " << method;
  if(method==2) cout << " (modulus = " << DEFAULT_MODULUS << ")";
  cout << endl;
  mat ref = echelon(aug, pc, npc, rk, ny, denom, method);
  cout << "Echelon matrix = " << ref;
  cout << "pivotal columns: " << pc << endl;
  cout << "nonpivotal columns: " << npc << endl;
  cout << "Denom = " << denom << endl;

  for (i=1,rk=0; (i<=r)&&(pc[i]<=r); i++,rk++) ;
  ny = r-rk;
  cout << "Rank = " << rk << endl;
  cout << "Nullity = " << ny << endl;
  if (rk<r) // non-invertible!
    {
      cout << "A is not invertible; rk = " << rk << endl;
    }
  else
    {
      mat ainv = ref.slice(1,r,r+1,r+r);
      cout << "A has inverse ";
      if (denom>1)
        cout << "(1/" << denom << ")*";
      cout << ainv;
      cout << "Check: A.A^(-1) = I ?";
      if (a*ainv == mat::scalar_matrix(r,denom))
        cout << " True!";
      else
        cout << " False!";
      cout << endl;
    }

  stop_time();
  // commented out for sutomatice tests:
  //cout << "cpu time = "; show_time(); cout << endl;
}
