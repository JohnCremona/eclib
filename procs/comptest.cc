// comptest.cc: test program for complex functions
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
 
#include "compproc.h"
#include "marith.h"

int main(void)
{ 
  set_precision("Enter number of decimal places");
  bigfloat x=to_bigfloat(3.125), y=to_bigfloat(4.25);

   bigcomplex z = bigcomplex(x,y);
   bigcomplex a = bigcomplex(to_bigfloat(1),to_bigfloat(1));
   bigcomplex b = bigcomplex(to_bigfloat(2),to_bigfloat(1));
   bigcomplex c;
   cout << "z = " << z << "\n";
   cout << " has real part = " << real(z) << "\n";
   cout << " and imaginary part = " << imag(z) << "i\n";
   cout << "z has complex conjugate = " << conj(z) << "\n";
   c = cagm(z,to_bigfloat(1));
   cout << "AGM(" << z << "," << bigcomplex(to_bigfloat(1)) << ") = ";
   cout << c << "\n\n";
   cout << "AGM(" << a << "," << b << ") = "<<flush;
   c = cagm(a,b);
   cout << c << "\n\n";
   bigint ia,ib,ic,id,ie; long nr;
   vector<bigint> roots;

   cout << "Enter Integer coefficients of a (monic) cubic:";
   cin>>ia>>ib>>ic;
   roots = Introotscubic(ia,ib,ic);
   nr=roots.size();
   if (nr==0) cout << "No integer roots"<<endl;
   else       cout << "The "<<nr<<" root(s) are:\n"<<roots<<endl;

// Test of complex cube root

  bigcomplex rootz;
  bigfloat three(to_bigfloat(3));
  bigcomplex w =  bigcomplex(to_bigfloat(-1), 
			     sqrt(three))/to_bigfloat(2);

  cout << "Enter a real or complex: ";  cin >> z;
  rootz=exp(log(z)/three);
  cout << "Main cube root = " << rootz << endl;
  cout << "whose cube is    " << pow(rootz,3) << endl;
  rootz*=w;
  cout << "Next cube root = " << rootz << endl;
  cout << "whose cube is    " << pow(rootz,3) << endl;
  rootz*=w;
  cout << "Next cube root = " << rootz << endl;
  cout << "whose cube is    " << pow(rootz,3) << endl;

// Test for quartic root-finding: 
  
  bigfloat xa,xb,xc,xd,xe;
  cout << "Enter real coefficients a b c d e of a quartic:";
  cin>>xa>>xb>>xc>>xd>>xe;

  int iroot;
  vector<bigcomplex> croots = solverealquartic((xa), (xb),(xc),(xd),(xe));
  cout<<"Quartic [" <<xa << "," << xb << "," << xc << "," << xd << "," << xe 
      << "] has roots:\n";
  for(iroot=0; iroot<4; iroot++) cout<<croots[iroot]<<"\n";
  cout<<endl;

   abort();
}
