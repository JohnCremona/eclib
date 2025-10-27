// tequiv.cc: test program for quartic equivalence
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
#include <eclib/mquartic.h>
#include <eclib/mequiv.h>
#define NEQPLIST 0        // Number of primes for equiv-test sieving

vector<long> eqplist;

int getquartic(quartic& g, int verbose)
{
  ZZ a, b, c, d, e;
  
  if(verbose)  cout << "Enter quartic coefficients a,b,c,d,e ?" << endl;
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
  g.set_equiv_code(eqplist);
  return 1;
}

int main()
{
  initprimes("PRIMES",0);
  cout.precision(50);
  cin.flags( cin.flags() | ios::dec );
  
  int verb; cout << "Verbose? "; cin >> verb;
  while (1)
    {
  long nq;  cout << "How many quartics to check (0 to quit)? ";  cin >> nq;
  cout<<endl<<endl;
  if (nq<1) {return 0;}
  vector<quartic> glist(nq);
  vector<ZZ> dlist;
  int i,j;

  for(i=0; i<nq; i++)
    {
      getquartic(glist[i], verb);
      cout<<(i+1)<<": "<<glist[i];
      if(verb) cout<<endl;
      int eq=0;
      for(j=0; (j<i)&&(!eq); j++)
	{
	  eq = new_equiv(glist[i],glist[j],verb);
	  if(eq) cout<<" equivalent to #"<<(j+1)<<endl;
	}
      if(!eq) cout<<" new"<<endl;
    }
    }
  return(0);
}
