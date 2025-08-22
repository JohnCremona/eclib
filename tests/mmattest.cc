// mmattest.cc: Multiprecision matrix package test program
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
 
#include <time.h>
#include <eclib/method.h>
#include <eclib/matrix.h>

int main(void)
{
  time_t starttime,stoptime;
  time(&starttime);
  cout << "\nMultiprecision matrix package test program.\n\n";
  {
    long i, r;
    mat_m a,aug,ref;
    vec_i pc(1),npc(1); vec_m poly(1);
    {
      cout << "Enter size of a square matrix A: "; cin >> r;
      a.init(r,r);
      cout << "Enter entries of A: "; cin >> a;
      cout << "A = " << a;
    }

    {
      cout << "Creating an array of 3 matrices\n";
      vector<mat_m> matlist(3);
      matlist[0] = a;
      matlist[1] = bigint(2)*a;
      matlist[2] = bigint(3)*a;
      cout << " A=" << matlist[0];
      cout << "2A=" << matlist[1];
      cout << "3A=" << matlist[2];
    }

{
for (i=1; i<=r; i++)
 cout << "row(A,"<<i<<") = " << a.row(i) << endl;
cout << "A = " << a;
for (int j=1; j<=r; j++)
 cout << "col(A,"<<j<<") = " << a.col(j) << endl;
cout << "A = " << a;
cout << "directsum(A,A) = " << directsum(a,a);
cout << "Enter any number "; cin >> i;
}
{
  mat_i sa = to_mat_i(a);
  cout << "After shortening to a matrix of ints, A = " << sa;
}
mat_m b = a;
{
cout << "B = A = " << b;
cout << "Enter any number "; cin >> i;
cout << "B==A?" << (b==a) << endl;
cout << "B!=A?" << (b!=a) << endl;
b+=a;
cout << "after B+:=A, A = " << a << "and B = " << b;
cout << "Enter any number "; cin >> i;
}
{
b-=a;
cout << "after B-:=A, A = " << a << "and B = " << b;
cout << "Enter any number "; cin >> i;
bigint two; two=2;
b*=two;
cout << "after B*:=2, A = " << a << "and B = " << b;
cout << "Enter any number "; cin >> i;
b/=two;
cout << "after B/:=2, A = " << a << "and B = " << b;
cout << "Enter any number "; cin >> i;
}
{
cout << "A+B=" << (a+b);
cout << "Now A = " << a << "and B = " << b;
cout << "Enter any number "; cin >> i;
cout << "A-B=" << (a-b);
cout << "Now A = " << a << "and B = " << b;
cout << "Enter any number "; cin >> i;
cout << "A*B=" << (a*b);
cout << "Now A = " << a << "and B = " << b;
}
{
cout << "Enter any number "; cin >> i;
cout << "-A=" << (-a);
cout << "Now A = " << a;
cout << "-A=" << (-a);
cout << "Now A = " << a;
cout << "Enter any number "; cin >> i;
}
{
vector<bigint> cp = a.charpoly();
cout << "char. poly. of A has coefficients " << cp << endl;
cout << "det(A) = " << a.determinant() << endl;
}
{
  aug = colcat(a, mat_m::identity_matrix(r));
  cout << "Augmented matrix = " << aug << endl;
}

long rk, ny;
bigint denom;
{
int method;
cout << "Which echelon method? (0=standard,1=longlong,2=modular) ";
cin>>method;
ref = echelon(aug, pc, npc, rk, ny, denom, method);
cout << "Echelon matrix = " << ref;
cout << "pivotal columns: " << pc << endl;
cout << "nonpivotal columns: " << npc << endl;
cout << "Denom = " << denom << endl;
}

for (i=1,rk=0; (i<=r)&&(pc[i]<=r); i++,rk++) ;
ny = r-rk;
cout << "Rank = " << rk << endl;
cout << "Nullity = " << ny << endl;
if (rk<r) // non-invertible!
 {cout << "A is not invertible; rk = " << rk << endl;
 }
else
 {mat_m ainv = ref.slice(1,r,r+1,r+r);
  cout << "A has inverse ";
  if (!is_one(denom)) cout << "(1/" << denom << ")*";
  cout << ainv;
  cout << "Check: A.A^(-1) = I ?";
  if (a*ainv == mat_m::scalar_matrix(r,denom)) cout << " True!";
  else cout << " False!";
  cout << endl;
 }
// Clear up:
 ref.init(); 
 pc.init(); npc.init(); 
 aug.init(); 
 a.init(); 
}
time(&stoptime);
//cout << "cpu time = " << (stoptime-starttime) << " seconds\n";
}
