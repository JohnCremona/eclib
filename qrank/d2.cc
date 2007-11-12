// d2.cc:  program for 2nd descent from a phi-descent homogeneous space
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
#include "unimod.h"
#include "points.h"
#include "mquartic.h"
#include "transform.h"
#include "msoluble.h"
//#include "samir.h"
#include "qc.h"
#include "quadratic.h"
#include "conic.h"
#include "minim.h"
#include "reduce.h"
#include "sqfdiv.h"
#include "desc2.h"


bigcomplex roots[4];
int getquartic(quartic& g);  // special version for a*x^4+c*x^2+e quartics

int main()
{
  set_precision("Enter number of decimal places");

  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  
  bigint zero; zero=0;
  int alldesc=0, verb=0, selmer_only=0; 
  cout << "Verbose? "; cin >> verb;
  initprimes("PRIMES",0);
  double hlim=8;
  cout << "Limit on height? "; cin >> hlim;
  cout << "Stop after first point found? "; cin >> alldesc; alldesc=1-alldesc;
  cout << "Selmer only (0/1: if 1, just tests whether second descent possible)? "; 
  cin >> selmer_only;

  quartic g;

  while (getquartic(g))
    {

      bigint I = g.getI(), J=g.getJ();
      Curvedata IJ_curve(zero,zero,zero,-27*I,-27*J,0);
      bigint tr_u,tr_r,tr_s,tr_t;
      Curvedata E = IJ_curve.minimalize(tr_u,tr_r,tr_s,tr_t);

      bigint d1=g.geta(), c=g.getcc(), d2=g.gete();
      bigint d = d1*d2, cdash = -2*c, ddash = sqr(c)-4*d;
      vector<bigint> plist = pdivs(6*d*ddash);
      vector<bigint> supp  = support(ddash);

      long mask=0;

      Curvedata E_cd(zero,c,zero,d,zero);
      cout << "cd-curve       (nearer):\t" << (Curve)E_cd << endl;
      Curvedata E_cd_dash(zero,cdash,zero,ddash,zero);
      cout << "cd-dash-curve (further):\t" << (Curve)E_cd_dash << endl;
      Point P(&E_cd);
      Point Q(&E_cd_dash);
      cout << "I = " << I << ", J = " << J << "\n";
      cout << "Minimal model for Jacobian: " << (Curve)E << endl;

      cout << "Checking local solublity at primes " << plist << ":\n";
      int i, els, els1; bigint two; two=2;
      els = els1 = qpsoluble(g,two);
      if(!els1) cout << "Not locally soluble at p = 2\n";
      for (i=1; i<plist.size(); i++)
	{
	  els1=new_qpsoluble(g,plist[i]);
	  if(!els1) cout << "Not locally soluble at p = "<<plist[i]<<"\n";
	  els = els&els1;
	}
      if(!els) continue;
      
      cout << "Everywhere locally soluble.\n";
      cout<<"------------------------------------------\n";

      bigint x,y,z,x2,xz,z2;

      int res = desc2(c,d1,d2,plist,supp,supp,mask,hlim,x,y,z,verb,selmer_only,alldesc);

      cout << "\nRESULTS\n\n";

      switch(res)
	{
	case 1:
	  cout<<"Quartic has rational point "; show_xyz(x,y,z);cout<<endl;
	  xz=x*z; x2=x*x; z2=z*z;
	  P = Point(&E_cd, d1*x2*z, d1*x*y, z*z2);

	  cout << "Point on  (c,d)  curve = " << P << "\n";
	  cout << "height = " << height(P)<< "\n";
	  Q = Point(&E_cd_dash,y*y*xz,y*(d1*x2*x2-d2*z2*z2),pow(xz,3));
	  cout << "Point on (c',d') curve = " << Q << "\n";
	  cout << "height = " << height(Q)<< "\n\n";
	  break;
	case -1:
	  cout<<"Quartic has no rational point (no ELS descendents)\n\n";
	  break;
	case 0:
	default:
	  cout<<"Undecided: quartic has ELS descendents but no rational\n";
	  cout<<"points were found on any.\n\n";
	}
    }
}


int getquartic(quartic& g)  // special version for a*x^4+c*x^2+e quartics
{
  bigint a, b, c, d, e;
  
  cout << "Enter quartic coefficients (a,0,c,0,e) or just a c e " << endl;
  char ch; cin>>ch;
  if(ch=='(') cin>>a>>ch>>b>>ch>>c>>ch>>d>>ch>>e>>ch;
  else 
    {
      cin.putback(ch);
      cin >> a >> c >> e;
      b=0; d=0;
    }

     
  if (sign(a)==0)  return 0;
    
  g=quartic(a,b,c,d,e);
  return 1;
}
