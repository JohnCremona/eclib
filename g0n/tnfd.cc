// FILE TNFD.CC:  test program for nfd (d-dimensional newform) class
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
//
//
#include <iostream>
#include "timer.h"
#include "msubspace.h"
#include "moddata.h"
#include "symb.h"
#include "homspace.h"
#include "nfd.h"

#define OUTPUT_PARI_STYLE
//#define DEBUG
//#define COMPARE_OLD

int main()
{
 init_time();
 cout << "Program tnfd." << endl;
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 long n=1; 
 int plus=1;
 int verbose=1;
 int w_split=0;
 int mult_one=0;
 int one_p=0;
 cout << "Verbose output? (0/1) "; cin >> verbose;
 // cout << "Plus space (0/1)? "; cin >> plus;
 while (cout<<"Enter level: ", cin>>n, n>0)
   {
     cout << ">>>Level " << n << "\t";
     homspace hplus(n,plus,0);
     int dimh = hplus.h1dim();
     cout << "dimension = " << dimh << endl;

     cout << "Split into W-eigenspaces (0/1)? "; 
     cin >> w_split;
     cout << "Multiplicity 1 eigenspaces only? (0/1)? "; 
     cin >> mult_one;
     cout << "Use just one T_p (1) or a linear combination (0)? "; 
     cin >> one_p;
     nfd form = nfd(&hplus, one_p, w_split, mult_one, verbose);
     long dims = dim(form.S);
     if(dims==0) continue;

     bigint den=form.dHS;
     int i, ip, nap=5, bad;
     cout<<"Number of ap? ";  cin>>nap;
     primevar pr; long p;
     start_time();
     for(ip=0; ip<nap; ip++, pr++)
       {
	 p=pr; bad = ::div(p,n);
	 if(verbose)
	   {
	     mat_m tp = form.oldheckeop(p);
	     if(den>1) cout<<den<<"*";
	     cout<<"Matrix of ";
	     if(bad) cout<<"W("; else cout<<"T(";
	     cout <<p<<") = "; 
	     showmatrix(tp);
	     vector<bigint> cptp = charpoly(tp);
	     for(i=0; i<dims; i++)
	       {
		 bigint temp = cptp[i];
		 divide_exact(temp,form.Hscales[dims-i],temp);
		 divide_exact(temp,form.Sscales[dims-i],temp);
		 cptp[i]=temp;
	       }
	     cout<<"char poly = "<<cptp<<endl;
	   }

	 vec_m apvec = form.ap(p);

	 if(bad) cout<<"w_"; else cout<<"a_";
	 cout<<p<<" = ";
	 
	 if(dims==1)
	   {
	     cout<<apvec[1]<<endl; 
	   }
	 else
	   {
	     cout<<apvec<<endl;
	   }
       }
     stop_time();
     show_time();
     cout<<endl;
   }
 abort();
}       // end of main()

