// FILE LF1PER.CC:   -- program to compute L(f_chi,1) for newforms
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
// (needs upgrading)

#include <time.h>
#include <fstream.h>
#include "marith.h"
#include "moddata.h"
#include "symb.h"
#include "oldforms.h"
#include "homspace.h"
#include "newforms.h"
#include "h1newforms.h"
#include "periods.h"

//#define AUTOLOOP

int main(void)
{
 cout.precision(15);
 int limit,n=1; 
#ifdef AUTOLOOP
 cout<<"Enter first and last N: ";cin>>n>>limit; 
 n--; cout<<endl;
 while (n<limit) { n++;
#else
 while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
 int usedata=1;
 newforms nf(n,usedata);
 int inf; cout << "Enter form number: "; cin>>inf; inf--;
 for(int i=inf; (i==inf)&&(i<nf.n1ds); i++)
   {
     for(int c=0; c<80; c++) cout<<"-"; 
     cout<<endl<<"Form number "<<i+1<<endl;
     lfchi lx(nf.nflist+i);
     int s = nf.nflist[i].sfe;
     rational loverp = nf.nflist[i].loverp;  // Converting from rational to double 
     double x0=10, y0=10;
     int first1=1, first3=1;
     if(loverp!=0) 
       {
	 double y=lx.value(1);
	 if(abs(y)>0.001)
	   {
	     x0=y/double(loverp); first1=0;
	     cout << "Real period = " << x0 << "(Ratio+ = 1)" << endl;
	   }
       }
     for(primevar l; l<200; l++)
       {
	 int ell = l;
	 int threemod4=(l%4)==3;
	 if(legendre(-n,l)!=s) continue;  // i.e. skip to next l
	 double y = lx.scaled_value(ell);
	 cout << "l = " << l ;
	 cout << "\tPer(l) = " ;
	 if(threemod4) cout << "i*";
	 cout << y;
	 if(first3&&threemod4&&abs(y)>0.001){y0=y; first3=0;}
	 if(first1&&!threemod4&&abs(y)>0.001){x0=y; first1=0;}
	 if(threemod4) cout << "\t\tRatio- = " << y/y0 << endl;
	 else cout << " Ratio+ = " << y/x0 << endl;
       }
   }
}       // end of if(n)
}       // end of while()
}       // end of main()
