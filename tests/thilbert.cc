// thilbert.cc: test of Hilbert symbol functions
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include <eclib/marith.h>
#include <eclib/quadratic.h>
#include <eclib/conic.h>
#include <eclib/hilbert.h>

//#define AUTO

int main()
{
  bigint a,b,p,x0,y0,z0; int resp, res, checkres;
  bigint zero, one;  zero=0; one=1;
#ifdef AUTO
  long la,lb,abmax;
  cout<<"Enter max for a,b: "; cin>>abmax;
  for(la=-abmax; la<=abmax; la++)
  for(lb=-abmax; lb<=abmax; lb++)
    {
      a=la; b=lb;
      if(a*b==0) continue;
#else
  while(1) {
    cout<<"Enter nonzero a and b: "; cin>>a>>b;
    if(a*b==zero) break;
#endif
    res=checkres=0;
    cout<<"("<<a<<","<<b<<")_p\n";
    resp = local_hilbert(a,b,0);
    cout<<0<<"\t"<<resp<<endl;
    res = res|resp;
    checkres = checkres^resp;

    resp = local_hilbert(a,b,2);
    cout<<2<<"\t"<<resp<<endl;
    res = res|resp;
    checkres = checkres^resp;

    vector<bigint> plist = vector_union(pdivs(a),pdivs(b));
    for (vector<bigint>::iterator pr = plist.begin(); pr!=plist.end(); pr++)
      {
	p=*pr; if(p==2) continue;
	resp=local_hilbert(a,b,p);
	cout << p << "\t" << resp << endl;
	res = res|resp;
	checkres = checkres^resp;
      }
    cout<<"\nGlobal symbol = " << res << endl;
    cout<<"Check (should be 0) = " << checkres << endl;

    int gres = global_hilbert(a,b,plist,p);
    if(res==gres)
      cout<<"--agrees with single call to global_hilbert()\n";
    else
      cout<<"--DISAGREES with single call to global_hilbert()\n";

    quadratic q(a,zero,b);
    int oldres = !solve_conic(q,one,x0,y0,z0,4);
    if(oldres==res)
      {
	cout<<"--agrees with solve_conic()\n";
      }
    else cout<<"--DISAGREES with solve_conic()\n";
  }
  cout<<endl;
}
