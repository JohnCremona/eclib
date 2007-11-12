// fixforms.cc -- special prog for recomputing specific newforms
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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

#include <time.h>
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"

//#define AUTOLOOP

int main(void)
{
 cout << "Program fixforms.  Using METHOD = " << METHOD << " to find newforms" << endl;
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 long starttime,stoptime,cputime;
 time(&starttime);
 long n; 
 cout << "Enter level: "; cin >> n;
 int plus=1;
 int verbose=0;
 cout << "See the hecke matrices (0/1)? "; cin >> verbose;

 cout << ">>>Level " << n << "\t" << flush;
 homspace hplus(n,plus,0);
 int genus = hplus.h1dim();
 cout << "Dimension = " << genus << "\n";

 subspace s(genus);
 long p, ap, newdim=genus;
 while(cout<<"Enter prime p and eigenvalue ap: ", cin>>p>>ap, (p>1)&&(newdim>1))
   {
     cout<<"Using p = " << p << ", ap = " << ap << endl;
     cout<<"Computing matrix..."<<flush;
     matrix t0 = transpose(hplus.newheckeop(p,verbose));
     matrix t  = restrict(t0,s);
     cout<<"done."<<endl;
     SCALAR lambda = hplus.h1denom()*ap*denom(s);
#ifdef MODULAR
     SUBSP temp = peigenspace(t,lambda,MODULUS);
     SUBSP newspace = pcombine(s,temp,MODULUS);
#else
     SUBSP newspace = combine(s,eigenspace(t,lambda,METHOD));
#endif
     newdim = dim(newspace);
     cout << "Current dimension = " <<  newdim << "\n";
     s = newspace;
   }
 if(newdim==1)
   {
     VEC bas1 = basis(s).col(1);
#ifdef MODULAR
// Lifts to Z from mod MODULUS
  if(verbose) cout<<"Before lifting: v="<<bas1<<endl;
  if(!liftok(bas1,MODULUS))
    cout << "Problem lifting eigenvector!"<<endl; 
//if(verbose) cout<<"After  lifting:\nv="<<bas1<<endl;
#else
  makeprimitive(bas1);
#endif
#ifdef MULTI
  vector bas =  bas1.shorten();
#else
  vector bas =  bas1;
#endif
     cout << "Final basis vector = \n" << bas << endl;
   }
 else
   {
     cout << "Final dimension = " <<newdim << " != 1" << endl;
   }
}       // end of main()
