// rattest.cc: test of rational number package
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
 
#include <iostream>
#include "rat.h"

int main(void)
{
 cout << "Test run of rational number procs.\n\n";
 long n = 4, d = 6, g;
 g = gcd(n,d);
 cout << "gcd(" << n << "," << d << ") = " << g << "\n";
 rational q = rational(n,d);
 cout << "q = " << q << " has denominator " << den(q) << \
         " and numerator " <<  num(q) << "\n";
 while (0) {
 cout << "Enter numerator and denominator\n";
 cin>>ws;  if(cin.eof()) {cout<<endl; exit(0);}
 cin >> n;
 cin >> d;
 q=rational(n,d);
 cout << "q = " << q << " has denominator " << den(q) << \
         " and numerator " <<  num(q) << "\n";
 }
 cout << "Enter a rational (either n or n/d): ";
 cin>>q;
 cout << "q = " << q << " has denominator " << den(q) <<	\
   " and numerator " <<  num(q) << "\n";
 
 rational q1,q2,q3;
 cout << "Enter three rationals separated by whitespace: ";
 cin>>q1>>q2>>q3;
 cout << "q1 = " << q1 << "\t";
 cout << "q2 = " << q2 << "\t";
 cout << "q3 = " << q3 << "\n";

 cout << "Enter three rationals in the format [q1:q2:q3]: ";
 char c;
 cin>>c>>q1>>c>>q2>>c>>q3>>c;
 cout << "[q1:q2:q3] = [" << q1 << ":";
 cout << q2 << ":";
 cout << q3 << "]\n";
}
