// FILE HPMCURVE.CC: 
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
// computes a possible equation from L(f_chi,1) for newforms 
// ONLY needs eigs data, not full intdata 
// BUT only finds curves "up to isogeny" as it cannot determine 
// the full Gamma_0(N) period lattice.

#include <fstream>
#include "marith.h"
#include "moddata.h"
#include "symb.h"
#include "oldforms.h"
#include "homspace.h"
#include "newforms.h"
#include "cperiods.h"     //from qcurves, for computing conductors
#include "h1newforms.h"
#include "periods.h"

#ifndef SINGLE
#define AUTOLOOP
#endif

void ratapprox(bigfloat x, long& a, long& b);

int main(void)
{
  set_precision("Enter number of decimal places");
 long limit,n=1; 
 int detail, maxn, pmax=100;  // Upper bound for twisting primes used
 cout << "See details? "; cin>>detail;
 cout << "Enter max prime for twisting: ";  cin>>pmax;
 cout << "Enter max scaling factor for real and imaginary periods: "; cin>>maxn;
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
  newforms nf(n,1);
  long fac=nf.sqfac;
  long fac6=(odd(fac)?fac:2*fac);
  if(detail&&(fac>1)) cout<<"c6 factor " << fac6 << endl;
  for(int i=0; (i<nf.n1ds); i++)
    //int i=4;
   {
     if(detail) 
       cout<<"\n"<<"Form number "<<i+1<<"\n";
     else cout<<(i+1)<<" ";
     newform& nfi = nf.nflist[i];
     lfchi lx(&nf, &nfi);
     int s = nfi.sfe;
     bigfloat x0, y0, ratio; long a,b;
     x0=y0=10;
     int first1=1, first3=1;
     long lplus, lminus, mplus=1, mminus=1;
     lminus=lplus=0; int ell=0;
     for(primevar l; (ell<pmax); l++)
       {
	 ell = l; 
	 int ell1mod4 = (ell%4==1);
	 if(legendre(-n,l)!=s) continue;  // i.e. skip to next l
	 lx.compute(ell);
	 bigfloat y = lx.scaled_value();
	 if(first3&&(!ell1mod4)&&abs(y)>0.001) {y0=y; first3=0; lminus=ell;}
	 if(first1&&(ell1mod4)&&abs(y)>0.001) {x0=y; first1=0; lplus=ell;}
	 ratio = (ell1mod4? y/x0: y/y0);
	 ratapprox(ratio, a, b);
	 if(b<501)   // avoid spurious values
	   {
	     if(ell1mod4) mplus=lcm(mplus,b);
	     else       mminus=lcm(mminus,b);
	   }
       	 if(detail) 
	   cout << "l = "<<ell<< "\tL(chi,1)="<<y<<"\tRatio = " << ratio << " =~= "<<a<<"/"<<b<<"\n";
       } // end of primes loop
     if(first3||first1)
       {
	 cout<<"No suitable twisting primes found under " << pmax << ".\n";
	 cout<<"Rerun with a bigger bound, possibly increasing the number of eigs also."<<endl;
       }
     if(detail) cout << "lplus = " << lplus << "; mplus = " << mplus << "\n";
     if(detail) cout << "lminus = "<< lminus <<"; mminus = "<< mminus<< "\n";
     x0/=to_bigfloat(mplus);
     y0/=to_bigfloat(mminus);
     if(detail) cout << "Minimal periods found: x0 = " << x0 << ", y0 = " << y0 << "\n";
     int maxny=maxn, maxnx=maxn;
     int type, nx, ny; bigcomplex w1, w2; bigcomplex c4, c6;
     bigfloat x1=x0, y1=y0;
     int notfound=1;

     for(nx=1; notfound&&(nx<=maxnx); nx++)
       {x1=x0/to_bigfloat(nx);
     for(ny=1; notfound&&(ny<=maxny); ny++)
       {
	 y1=y0/to_bigfloat(ny);
	 for(type=2; (type>0); type--) 
	   {
	     if(type==2){w1=bigcomplex(x1,to_bigfloat(0)); w2=bigcomplex(to_bigfloat(0),y1);}
	     else {w1=bigcomplex(2*x1,to_bigfloat(0)); w2=bigcomplex(x1,y1);}
	     bigcomplex tau=normalize(w1,w2);
	     getc4c6(w1,w2,c4,c6);
	     bigint ic4 = fac*Iround(real(c4)/fac);
	     bigint ic6 = fac6*Iround(real(c6)/fac6);
	     int close=abs(I2bigfloat(ic4)-real(c4))<0.001;
	     int validc4c6 = valid_invariants(ic4,ic6);
	     if((validc4c6||close)&&detail)
	       {
		 cout << "type = " << type << ", nx,ny = " << nx<<","<<ny << "\n"; 
		 cout << "w1 = " << w1 << ", w2 = " << w2 << "\n";
		 cout << "c4 = " << real(c4) << ", c6 = " << real(c6) << "\n";
		 cout << "ic4 = " << ic4 << ", ic6 = " << ic6 << "\n";
	       }
	     if(validc4c6)
	       {
		 Curve C(ic4,ic6);
		 Curvedata CD(C,1);
		 CurveRed CR(CD);
		 bigint cond = getconductor(CR);
		 if(cond==n)
		   {
		     if(detail)cout<<"Curve is ";
		     else if(!notfound)cout<<"   ";
		     cout << (Curve)CD << "  ";
		     if(detail) cout << "N = " << cond << "  ";
		     cout << "(l- = "<<lminus<<", m- = "<<mminus*ny<<", l+ = "<<lplus<<", m+ = "<<mplus*nx<<", type = "<<type<<", index = "<<mplus*mminus*nx*ny<<")\n";
		     notfound = 0;
		   }
	       }
	   } // end of type loop
       } // end of ny loop
       } // end of nx loop
     if(notfound) 
       cout << "No curve found: try again with higher limits for maxny" << "\n";
}       // end of forms loop
}       // end of if(n)
}       // end of while()
}       // end of main()


void ratapprox(bigfloat x, long& a, long& b)
{
  long c, x0, x1, x2, y0, y1, y2;
  bigfloat xx, diff, eps = to_bigfloat(1.0e-6);
  xx = x; x0 = 0; x1 = 1; y0 = 1; y1 = 0;
//cout<<"ratapprox("<<x<<"): \n";
  diff = 1; c=x2=y2=0;
  while ( diff > eps )
    { c = longify( xx + ( (xx>0) ? 0.5 : -0.5 ) ); // ie round(xx)
      x2 = x0 + c*x1; x0 = x1; x1 = x2;
      y2 = y0 + c*y1; y0 = y1; y1 = y2;
      diff = abs( x - (to_bigfloat(x2)/to_bigfloat(y2)) );
//cout<<"a = "<<x2<<", b = "<<y2<<", diff = "<<diff<<endl;
      if ( abs(xx - c) < eps ) diff = 0;
      else xx = 1/(xx - c);
    }
  a = x2; b = y2;
  if ( b < 0 ) {a=-a; b=-b; }
  if ( x < 0 ) {a=-a;}
}

