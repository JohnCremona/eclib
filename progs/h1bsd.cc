// FILE H1BSD.CC: Program to compute L^(r)(f,1) for newforms 
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
#include <fstream>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/oldforms.h>
#include <eclib/homspace.h>
#include <eclib/cperiods.h>     //from qcurves, for computing conductors
#include <eclib/newforms.h>
#include <eclib/periods.h>
#include <eclib/curvesort.h>

#define AUTOLOOP

int main(void)
{
  set_precision(60);
 int limit,n=1; 
#ifdef AUTOLOOP
 cout<<"Enter first and last N: ";cin>>n>>limit; 
 n--; cout<<endl;
 while (n<limit) { n++;
#else
 while (n>1) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
 newforms nf(n,0);
 int noldap=25;
 nf.createfromdata(1,noldap,0); // do not create from scratch if data absent
 for(int xi=0; xi<nf.n1ds; xi++)
   {
     int i = xi;
     string code = codeletter(xi);
     i=booknumber0(n,i);
     newform& nfi = nf.nflist[i];
     bigfloat lf1 = nfi.special_value();
     long r = nfi.rank();
     rational loverp = nfi.loverp;
     // loverp = L(f,1)/x where the period lattice is [x,yi] or [2x,x+yi], so in BSD we want L(f,1)/2x
     loverp /= 2;
     cout << n << " ";
     cout << codeletter(i) << "\t";
     cout << "Rank = " << r << "\t";
     cout << "L^(r)(f,1)/r! = " << lf1 << " \t";
     cout << "L(f,1)/period = " << loverp << endl;
   }
}       // end of if(n)
}       // end of while()
}       // end of main()
