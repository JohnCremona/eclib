// ratpoint.cc:  program to search for points on quartics and map to curve
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
 
#include "marith.h"
#include "points.h"
#include "mquartic.h"
#include "msoluble.h"
#include "qc.h"
#include "version.h"


int getquartic(quartic& g)
{
  bigint a, b, c, d, e;
  
  cout << "Enter quartic coefficients a,b,c,d,e ?" << endl;
  char ch; cin>>ch;
  if(ch=='(') cin>>a>>ch>>b>>ch>>c>>ch>>d>>ch>>e>>ch;
     else 
     {
       cin.putback(ch);
       cin >> a >> b >> c >> d >> e;
     }
     
  if (is_zero(a)&&is_zero(b)&&is_zero(c)&&is_zero(d)&&is_zero(e))
     return 0;
    
  g=quartic(a,b,c,d,e);  // will set its own invariants, roots and type
  return 1;
}

int main()
{
  show_version();
#ifdef LiDIA
  long lidia_precision=40;
  cout<<"Enter number of decimal places: "; cin>>lidia_precision;
  bigfloat::precision(lidia_precision);
#else
  cout.precision(15);
#endif
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  
  int verb; cout << "Verbose? "; cin >> verb;
  initprimes(string("PRIMES").c_str(),0);
  int modopt=0;
  //  cout<<"moduli option (0 (Stoll)/ 1/2/3)?";      cin >> modopt;

  quartic g;

  while (getquartic(g))
    {
      double hlim;
      cout << "Limit on height? "; cin >> hlim;

      bigint I = g.getI(), J=g.getJ(), zero=BIGINT(0);
      Curvedata IJ_curve(zero,zero,zero,-27*I,-27*J,0);
      bigint tr_u,tr_r,tr_s,tr_t;
      Curvedata E = IJ_curve.minimalize(tr_u,tr_r,tr_s,tr_t);

      cout << "I = " << I << ", J = " << J << "\n";
      cout << "Minimal model for Jacobian: " << (Curve)E << endl;

      bigint badp;
      vector<bigint> plist = pdivs(6*g.getdisc());
      unsigned int i, els, els1;

      cout << "Checking local solublity in R:\n";
      els = ((g.gettype()>1)||(g.geta()>BIGINT(0)));
      if(!els) cout << "Not locally soluble over R\n";

      cout << "Checking local solublity at primes " << plist << ":\n";
      els1 = qpsoluble(g,BIGINT(2));
      if(!els1) cout << "Not locally soluble at p = 2\n";
      els = els&els1;

      for (i=1; i<plist.size(); i++)
	{
	  els1=new_qpsoluble(g,plist[i],verb);
	  if(!els1) cout << "Not locally soluble at p = "<<plist[i]<<"\n";
	  els = els&els1;
	}
      if(!els) continue;
      
      cout << "Everywhere locally soluble.\n";
      
      quartic_sieve qs(&g,modopt,verb);
      cout << "Searching for points on "<<g<<" up to height "<<hlim<<endl;

      if(qs.search(hlim))
	{
	  bigint x, y, z;
	  qs.getpoint(x,y,z);
	  cout << "(x:y:z) = (" << x << ":" << y << ":" << z << ")\n";
	  Point P;
	  qc(g,x,y,z,&E,&IJ_curve,tr_u,tr_r,tr_s,tr_t,P,1);
	  cout << "Curve = " << (Curve)(E) << "\n";
	  cout << "Point = " << P << "\n";
	  cout << "height = " << height(P)<< "\n";
	}
      else cout << "No point found!\n";
    }
}
