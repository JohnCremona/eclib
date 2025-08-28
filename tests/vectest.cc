// vectest.cc: Test of vector package
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

int main(void)
{
  cout << "Vector test program with scalar type " << scalar_type << ".\n\n";
 int i,j,k,n;
 scalar one(1), two(2), three(3);

 cout << "iota(10) = " << vec::iota(10) << endl;
 cout << "Enter n : "; cin >> n;
 vec v;
 cout << "Uninitialized new vec v = " << v << endl;
 vec v2(v);
 cout << "Copy of v = " << v2 << endl;
 v.init(n);
 cout << "Initialized new vec v = " << v << endl;
 vec v3(v);
 cout << "Copy of v = " << v3 << endl;
 cout << "Enter new entries of v: ";
 cin >> v;
 cout << "Now v = " << v << endl;

 vec w(3);
 cout << "w = " << w <<  endl;
 w = v;
 cout << "After w=v, " << endl;
 cout << "Now v = " << v << endl;
 cout << "Now w = " << w << endl;
 cout << "w==v: " << (w==v) << endl;
 cout << "w!=v: " << (w!=v) << endl;
 cout << "Enter i : "; cin >> i;
 w*=two;
 cout << "After w*=2, w = " << w << endl;
 cout << "3*v = " << (three*v) << endl;
 cout << "Now v = " << v << endl;
 cout << "v+w = " << v+w << endl;
 cout << "Now v = " << v << endl;
 cout << "v-w = " << v-w << endl;
 cout << "Now v = " << v << endl;
 cout << "v*w = " << v*w << endl;
 cout << "Now v = " << v << endl;
 cout << "w/2 = " << w/two << endl;
 cout << "Now w = " << w << endl;
 cout << "-v  = " << -v  << endl;
 cout << "+v  = " << +v  << endl;
 cout << "+w  = " << +w  << endl;
 cout << "v = " << v << "; w = " << w << endl;

 cout << "Elements of v: \n";
 for (i=1; i<=n; i++) cout << "v[" << i << "] = " << v[i] << endl;

 cout << "Member test: Enter a test number: " ;
 scalar vi;
 cin >> vi; cout << vi;
 if (member(vi,v)) cout << " IS "; else cout << " IS NOT ";
 cout << "a member of v." << endl;

 cout << "Subscript test\n";
 cout << "Enter length of subscript vec:";
 int m; cin >> m; vec_i index(m);
 cout << "Enter subscript vector:";
 cin >> index;
 vec vv = v[index];
 cout << "The sub-vector is " << vv << endl;

 cout << "Change one entry of v.  Index?"; cin >> i;
 cout << "New entry?"; scalar x; cin >> x;
 v[i]=x;
 cout << "New entry: v[" << i << "] = " << v[i] << endl;
 cout << "Now v = " << v << endl;

 cout << "Initial slice; length? "; cin >> j;
 cout << "Slice = " << v.slice(j) << endl;
 cout << "Now v = " << v << endl;
 cout << "General slice; beginning, end? "; cin >> j >> k;
 cout << "Slice = " << v.slice(j,k) << endl;
 cout << "Now v = " << v << endl;
 cout << "w = " << w << "; content(w) = " << content(w) << endl;
 make_primitive(w);
 cout << "After make_primitive(w), w = " << w << endl;

 vec u(n);
 cout << "u = "<< u << endl;
 swapvec(u,v);
 cout << "After swapvec(u,v):\nu = " << u << "\nv = " << v << endl;
}
