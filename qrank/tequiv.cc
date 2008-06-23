// tequiv.cc: test program for quartic equivalence
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
#include "mquartic.h"
#include "mequiv.h"
#define NEQPLIST 0        // Number of primes for equiv-test sieving

vector<long> eqplist;

int getquartic(quartic& g, int verbose)
{
  bigint a, b, c, d, e;
  
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
  initprimes(string("PRIMES").c_str(),0);
#ifdef LiDIA
  long lidia_precision=40;
  cout<<"Enter number of decimal places: "; cin>>lidia_precision;
  bigfloat::precision(lidia_precision);
#else
  cout.precision(15);
#endif
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  
  int verb; cout << "Verbose? "; cin >> verb;
  long nq;  cout << "How many quartics to check? ";  cin >> nq;
  cout<<endl<<endl;
  quartic* glist = new quartic[nq];
  vector<bigint> dlist;
  int i,j;

  for(i=0; i<nq; i++)
    {
      getquartic(glist[i], verb);
      cout<<(i+1)<<": "<<glist[i];
      if(verb) cout<<endl;
#ifndef NEW_EQUIV
      if(i==0) dlist = sqdivs(glist[i].getdisc());
#endif
      int eq=0;
      for(j=0; (j<i)&&(!eq); j++)
	{
#ifdef NEW_EQUIV
	  eq = new_equiv(glist+i,glist+j,verb);
#else
	  eq = equiv(glist+i,glist+j,dlist,verb);
#endif
	  if(eq) cout<<" equivalent to #"<<(j+1)<<endl;
	}
      if(!eq) cout<<" new"<<endl;
    }
  delete[] glist;
  return(0);
}
