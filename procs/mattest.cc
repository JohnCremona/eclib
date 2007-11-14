// mattest.cc: Matrix package test program
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
 
#include "timer.h"
#include "arith.h"
#include "matrix.h"

int main(void)
{
  init_time();
  start_time();
  cout << "\nMatrix package test program.\n\n";
  {
    long i,j; scalar r;
    mat a,aug,ref;
    vec pc(1),npc(1),poly(1);
{
cout << "Enter size of a square matrix A: "; cin >> r;
a.init(r,r);
cout << "Enter entries of A: "; cin >> a;
cout << "A = " << a;
cout << "Using A.output(cout): ";  a.output(cout);
cout << "Using A.output_pari(cout): ";  a.output_pari(cout);
cout << "Using A.output_pretty(cout): \n";  a.output_pretty(cout);
}

 char* filename = new char[20];
 cout<< "Enter a filename for matrix binary output: "; cin>>filename;
 a.dump_to_file(filename);
 cout<< "Matrix dumped to file " << filename << endl;
 mat ax;
 ax.read_from_file(filename);
 cout<< "Matrix reread from file " << filename << endl;
 delete[] filename;
 cout << "B = " << ax;
 if(a==ax) cout<<"agree"; else cout << "WRONG!";
 cout<<endl;
 
 cout << "Enter any number "; cin >> i;

//  The following should work, doesn't in TURBO C++, but DOES with GCC!!

{
cout << "Creating an array of 3 matrices\n";
mat* matlist = new mat[3];
matlist[0] = a;
matlist[1] = 2*a;
matlist[2] = 3*a;
cout << " A=" << matlist[0];
cout << "2A=" << matlist[1];
cout << "3A=" << matlist[2];
delete[] matlist;
}

{
for (i=1; i<=r; i++)
 cout << "row(A,"<<i<<") = " << a.row(i) << endl;
cout << "A = " << a;
for (j=1; j<=r; j++)
 cout << "col(A,"<<j<<") = " << a.col(j) << endl;
cout << "A = " << a;
cout << "directsum(A,A) = " << directsum(a,a);
cout << "Enter any number "; cin >> i;
}
mat b = a;
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
b*=2;
cout << "after B*:=2, A = " << a << "and B = " << b;
cout << "Enter any number "; cin >> i;
b/=2;
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
vector<long> cp = charpoly(a);
cout << "char. poly. of A has coefficients " << cp << endl;
cout << "det(A) = " << determinant(a) << endl;
}
{
aug = colcat(a,idmat(r));
cout << "Augmented matrix = " << aug << endl;
}

long rk, ny;
scalar denom;
{
int method=0;
cout << "Which echelon method? (0=standard,1=longlong,2=modular) ";
cin>>method;
cout << "\nUsing method " << method;
if(method==2) cout << " (modulus = " << BIGPRIME << ")";
cout << endl;
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
 {mat ainv = ref.slice(1,r,r+1,r+r);
  cout << "A has inverse ";
  if (denom>1) cout << "(1/" << denom << ")*";
  cout << ainv;
  cout << "Check: A.A^(-1) = I ?";
  if (a*ainv == denom*idmat(r)) cout << " True!";
  else cout << " False!";
  cout << endl;
 }
// Clear up:
 ref.init(); 
 pc.init(); npc.init(); 
 aug.init(); 
 a.init(); 
}
stop_time();
//cerr << "cpu time = "; show_time(); cerr << endl;
}
