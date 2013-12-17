// FILE PCURVE.CC
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
// computes a possible equation from periods for newforms 
// ONLY needs eigs data, not full intdata 
// BUT only finds curves "up to isogeny" as it does not 
// determine the full Gamma_0(N) period lattice.

#include <iostream>
#include <eclib/marith.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/oldforms.h>
#include <eclib/homspace.h>
#include <eclib/cperiods.h>     //from qcurves, for computing conductors
#include <eclib/newforms.h>
#include <eclib/periods.h>
#include <eclib/pcprocs.h>

//#define SINGLE
#ifndef SINGLE
#define AUTOLOOP
#endif

#define DMAX 50  // Upper bound for d = (2,2) matrix entry used
#define LMAX 200  // Upper bound twisting prime l-


int main(void)
{
  initprimes("PRIMES",0);
  set_precision("Enter number of decimal places");
 long limit,n=1; 
 int dump=1, detail; 
 long maxn, dmax=DMAX;
 //#ifdef SINGLE
 detail=2;
 //#else
 cout << "See details? "; cin>>detail;
 //#endif
 cout << "Enter max d: ";  cin>>dmax;
 cout << "Enter max scaling factor for periods: "; cin>>maxn;
#ifdef AUTOLOOP
 cout<<"Enter first and last N: ";cin>>n>>limit; 
 n--; cout<<"\n";
 while (n<limit) { n++;
#else
 while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
  cout << "N = " << n << endl;

  newforms nf(n,detail);
  int noldap=25;
  nf.createfromdata(1,noldap);
  int squarelevel=nf.squarelevel;
  long fac=nf.sqfac;
  long nnf = nf.n1ds; 
  long inf = 1; 
#ifdef SINGLE
  cout << "Enter form number: "; cin>>inf;
  nnf=inf;
#endif
  primevar pr; long p0;                  // First "good" prime
  while (p0=(long)pr, ::divides(p0,n)) pr++; 
  
  int anyfound=0;
 for(long i=inf-1; i<nnf; i++)
   {
     if(detail) cout<<"\n"<<"Form number "<<i+1<<"\n";
     else cout<<(i+1)<<" ";
     newform* nfi = &((nf.nflist)[i]);
     if(detail>1) {cout<<"\n Newform details:\n";nfi->display();cout<<endl;}
     bigfloat x0=to_bigfloat(10), y0=to_bigfloat(10), ratio;
// flags to say that we have real/imag period 
//          and that we have a cycle giving both together
     int rp_known=0, have_both=0;
     long lplus=nfi->lplus, mplus=nfi->mplus;
     long s = nfi->sfe;
     //nfi->dotplus=1;     
     nfi->dotminus=1;  // In case wrong values were on the file!

// STEP 1: find the real period using L(f,1) is L/P nonzero,
//         or via L(f,chi,1)

     // Now done within get_matrix_periods()

// STEP 2: compute periods of lots of matriced in Gamma_0(N) over all
// symbols {0,b/d} for d<dmax.  We are looking for (i) a symbol whose
// imaginary part is nonzero; (ii) a symbol whose real and imaginary
// parts are both non-zero, which we store.

     have_both = nf.find_matrix( i, dmax, rp_known, x0, y0);
     if(!have_both)
       cout<<"Problem!  find_matrix() returns 0!"<<endl;

     if(detail) 
       {
	 cout << "Minimal periods found: x0 = "<<x0<<", y0 = "<<y0<<"\n";
	 cout << "Matrix ("<<nfi->a<<","<<nfi->b<<";"<<n*nfi->c<<","<<nfi->d<<"):\t";
	 cout << "dotplus = "<< nfi->dotplus << ", dotminus = "<< nfi->dotminus<< "\n";
	 cout << "Searching for scaling factors.\n";
       }

// STEP 3: Now we have nonzero integer multiples of real and imag
// periods, we search for a sub-multiple of both which gives a genuine
// curve.  
     
     long nx, ny; int type;  
     long maxnx=maxn; 
     if(rp_known) maxnx=1;
     int found = get_curve(n, fac, maxnx, maxn, x0, y0, nx, ny, type, detail);

     if(found)
       {
         anyfound=1;
	 // cout<<"before rescaling, dotplus="<<nfi->dotplus<<", dotminus="<<nfi->dotminus<<endl;
	 // cout<<"rescaling by nx = "<<nx<<", ny="<<ny<<endl;
	 nfi->dotplus *= nx;
	 nfi->dotminus *= ny;
	 nfi->type=type;
	 cout << "[(" <<nfi->a<<","<<nfi->b<<";"<<nfi->c
	      <<","<<nfi->d<<"),"<<nfi->dotplus<<","<<nfi->dotminus
	      <<";"<<nfi->type<<"]"<<endl;

// STEP 4:  We find a suitable twisting prime for the imaginary period

	 if(!squarelevel)
	   {
	     bigfloat y1 = y0/to_bigfloat(ny);
	     int ok = nf.find_lminus(i, 0, y1);
	     if(!ok)
	       {
		 cout<<"No suitable twisting prime "
		     <<" found for imaginary period!"<<endl;
	       }
	   }
	 if(detail) 
	   {
	     cout<<"Updated newform data:\n"; nfi->display();
	     cout<<"test reconstruction of curve from updated data:"<<endl;
	   }
	 // test reconstruction of the curve
	 bigfloat rperiod;
	 Curve C = nf.getcurve(i, -1, rperiod, detail);
	 Curvedata CD(C,1);  // The 1 causes minimalization
	 cout << "Curve = \t";
	 cout << (Curve)CD << "\t";
	 CurveRed CR(CD);
	 bigint nc = getconductor(CR);
	 cout << "N = " << getconductor(CR) << endl;
       }
     else cout<<"No curve found\n";

// dump complete data file into NF_DIR/x$N.

   }       // end of forms loop
 if (anyfound)
   {
     cout<<"Store newform data? "; cin>>dump;
     if(dump)
       {
         nf.output_to_file();
       }
   }
}       // end of if(n)
}       // end of while()
}       // end of main()
