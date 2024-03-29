// REDUCE_QUARTICS.CC:  Program for minimisation and reduction of quartics
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
//    
#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/points.h>
#include <eclib/mquartic.h>
#include <eclib/transform.h>
#include <eclib/msoluble.h>
#include <eclib/minim.h>
#include <eclib/reduce.h>

int getquartic(quartic& g);

int main()
{
  set_precision(600);
  cin.flags( cin.flags() | ios::dec );
  
  int verb=1; //0; 
  //   cout << "Verbose? "; cin >> verb;
  initprimes("PRIMES",verb);

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
