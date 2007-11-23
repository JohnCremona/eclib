// TRED.CC:  Test program for reduction of quartics to minimal I,J
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
//    
#include "marith.h"
#include "unimod.h"
#include "points.h"
#include "mquartic.h"
#include "transform.h"
#include "msoluble.h"
#include "minim.h"
#include "reduce.h"

int getquartic(quartic& g);

int main()
{
  set_precision(string("Enter number of decimal places").c_str());
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  
  int verb=1; //0; 
  //   cout << "Verbose? "; cin >> verb;
  initprimes(string("PRIMES").c_str(),verb);

  quartic g;

  while (getquartic(g))
    {
      cout<<"Quartic is "<<g<<endl;
      bigint I = g.getI(), J=g.getJ();
      cout<<"I = "<<I<<"\nJ = "<<J<<endl;
      bigint ga=g.geta(), gb=g.getb(), gc=g.getcc(), gd=g.getd(), ge=g.gete();
      bigint p, badp;
      vector<bigint> plist = pdivs(g.getdisc());
      cout<<"Bad primes: "<<plist<<endl;
      scaled_unimod m;

      cout << "Attempting to minimize the quartic.\n";
      bigint newa(ga), newb(gb), newc(gc), newd(gd), newe(ge);
      cout << "First partial minimization without assuming local solvability:\n";
      minim_all(newa,newb,newc,newd,newe,I,J,plist,m,0,1);
      quartic newg(newa,newb,newc,newd,newe);
      cout<<"Result has coefficients: "<<newg<<endl;
      cout<<"I = "<<I<<"\nJ = "<<J<<endl;
      plist = pdivs(newg.getdisc());
      cout<<"Bad primes: "<<plist<<endl;

      cout << "Now check local solvability:\n";
      int locsol = locallysoluble(newg, plist, badp);
      if(!locsol)
	{
	  cout << "Not locally soluble at p = " << badp << endl;
	  continue;
	}
      cout << "Everywhere locally soluble\n";
      cout << "Final minimization of I, J:\n";
      minim_all(newa,newb,newc,newd,newe,I,J,plist,m,1,1);
      newg.assign(newa,newb,newc,newd,newe);
      cout<<"Result has coefficients: "<<newg<<endl;
      cout<<"I = "<<I<<"\nJ = "<<J<<endl;
      plist = pdivs(newg.getdisc());
      cout<<"Bad primes: "<<plist<<endl;
      cout<<"transform = "<<m<<endl;
      if(check_transform(ga,gb,gc,gd,ge,m,newa,newb,newc,newd,newe))
	{
	  cout << "Attempting to reduce the quartic.\n";
	  unimod m1;
	  reduce(newa,newb,newc,newd,newe,m1);
	  newg.assign(newa,newb,newc,newd,newe,newg.getroots(),0,I, J, 4*pow(I,3)-J*J);
	  cout<<"Finished reducing g, new coefficients: "<<newg<<endl;
	  cout<<"I = "<<I<<"\nJ = "<<J<<endl;
	  cout<<"extra reducing transform = "<<m1<<endl;
	  m *= m1;
	  cout<<"total transform = "<<m<<endl;
	  if(check_transform(ga,gb,gc,gd,ge,m,newa,newb,newc,newd,newe))
	    {
	      cout << "OK\n";
	    }
	  else
	    {
	      cout << "check_transform fails after reduction!\n";
	    }
	}
      else
	{
	  cout << "check_transform fails after minimalization!\n";
	}
    }
}

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
     
  if (sign(a)==0&&sign(b)==0&&sign(c)==0&&sign(d)==0&&sign(e)==0)
     return 0;
    
  g=quartic(a,b,c,d,e);  // will set its own invariants, roots and type

  return 1;
}
