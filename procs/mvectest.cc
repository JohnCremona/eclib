// Test of vec_m package
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
 
#include "mvector.h"

int main(void)
{
 cout << "Test run of vector package.\n\n"; 
 int i,j,k;
 scalar n;
 cout << "Enter n : "; cin >> n;
 vec_m v;
 cout << "Uninitialized new vector v = " << v << endl;
 vec_m v2(v);
 cout << "Copy of v = " << v2 << endl;
 v.init(n);
 cout << "Initialized new vector v = " << v << endl;
 vec_m v3(v);
 cout << "Copy of v = " << v3 << endl;
 cout << "Enter new entries of v: ";
 cin >> v;
 cout << "Now v = " << v << endl;

 vec_m w(3);
 cout << "w = " << w <<  endl;
 w = v;
 cout << "After w=v, " << endl;
 cout << "Now v = " << v << endl;
 cout << "Now w = " << w << endl;
 cout << "w==v: " << (w==v) << endl;
 cout << "w!=v: " << (w!=v) << endl;
 cout << "Enter i : "; cin >> i;
 bigint two; two=2;
 w*=two;
 cout << "After w*=2, w = " << w << endl;
 cout << "3*v = " << (3*v) << endl;
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
 bigint vi;
 cin >> vi; cout << vi;
 if (member(vi,v)) cout << " IS "; else cout << " IS NOT ";
 cout << "a member of v." << endl;

 cout << "Subscript test\n";
 cout << "Enter length of subscript vector:";
 int m; cin >> m; vec index(m);
 cout << "Enter subscript vector:";
 cin >> index;
 vec_m vv = v[index];
 cout << "The sub-vector is " << vv << endl;

 cout << "Change one entry of v.  Index?"; cin >> i;
 cout << "New entry?"; bigint x; cin >> x;
 v[i]=x;
 cout << "New entry: v[" << i << "] = " << v[i] << endl;
 cout << "Now v = " << v << endl;

 cout << "Initial slice; length? "; cin >> j;
 cout << "Slice = " << v.slice(j) << endl;
 cout << "Now v = " << v << endl;
 cout << "General slice; beginning, end? "; cin >> j >> k;
 cout << "Slice = " << v.slice(j,k) << endl;
 cout << "Now v = " << v << endl;
 cout << "w = " << w << "; mvecgcd(w) = " << mvecgcd(w) << endl;
 makeprimitive(w);
 cout << "After makeprimitive(w), w = " << w << endl;

 vec_l sv = v.shorten(n);
 cout << "v shortened to a vector of longs: " << sv << endl;

 vec_m u(n);
 cout << "u = "<< u << endl;
 swapvec(u,v);
 cout << "After swapvec(u,v):\nu = " << u << "\nv = " << v << endl;
}
