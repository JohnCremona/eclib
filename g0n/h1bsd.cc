// FILE H1BSD.CC: Program to compute L^(r)(f,1) for newforms 
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
#include <fstream>
#include "moddata.h"
#include "symb.h"
#include "oldforms.h"
#include "homspace.h"
#include "cperiods.h"     //from qcurves, for computing conductors
#include "newforms.h"
#include "periods.h"
#include "curvesort.cc"

#define AUTOLOOP

int main(void)
{
  set_precision("Enter number of decimal places");
 int limit,n=1; 
 char* code = new char[20];
#ifdef AUTOLOOP
 cout<<"Enter first and last N: ";cin>>n>>limit; 
 n--; cout<<endl;
 while (n<limit) { n++;
#else
 while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
 newforms nf(n,1,0,0);
 int noldap=25;
 nf.createfromdata(noldap,0); // do not create from scratch if data absent
 for(int i=0; i<nf.n1ds; i++)
   {
     codeletter(i,code);
     //     i=booknumber0(n,i);
     newform& nfi = nf.nflist[i];
     ldash1 x(&nf, &nfi);  
     x.compute();
     bigfloat lf1 = abs(x.value());
     long r = x.rank();
     cout << n << "\t" << code 
	  << "\tRank = " << r << "\tL^(r)(f,1)/r! = " << lf1 << endl;
   }
}       // end of if(n)
}       // end of while()
}       // end of main()
