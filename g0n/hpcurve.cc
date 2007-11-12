// FILE HPCURVE.CC: 
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
#include "compproc.h"
#include "moddata.h"
#include "symb.h"
#include "oldforms.h"
#include "homspace.h"
#include "cperiods.h"     //from qcurves, for computing conductors
#include "newforms.h"
#include "periods.h"

//#define SINGLE
#ifndef SINGLE
#define AUTOLOOP
#endif

void ratapprox(bigfloat x, long& a, long& b);

int main(void)
{
  set_precision("Enter number of decimal places");
 initprimes("PRIMES",0);
 long limit,n=1; 
 int dump, detail; 
 long maxn, pmax=100;  // Upper bound for twisting primes used
#ifdef SINGLE
 detail=1;
#else
 cout << "See details? "; cin>>detail;
#endif
 cout << "Enter max prime for twisting: ";  cin>>pmax;
 cout << "Enter max scaling factor for imaginary period: "; cin>>maxn;
#ifdef SINGLE
 dump=0;
#else
 cout << "Dump parameters to file? "; cin >> dump;
#endif
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
  long rootn=(long)(sqrt((double)n)+0.1); int squarelevel=(n==rootn*rootn);
  if(squarelevel)
    {
      cout << "--square level, not implemented.\n";
      continue;
    }
  int usedata=1;
  newforms nf(n,usedata,5,0);
  long fac=nf.sqfac;
  long fac6=(odd(fac)?fac:2*fac);
  if(detail&&(fac>1)) cout<<"c6 factor " << fac6 << endl;
  long nnf = nf.n1ds; 
  long inf = 1; 
#ifdef SINGLE
  cout << "Enter form number: "; cin>>inf;
  nnf=inf;
#endif
  primevar pr; long p0;                  // First "good" prime
  while (p0=(long)pr, ::div(p0,n)) pr++; 
  
  ofstream data;
  if(dump&&(nnf>0))
    {
      char *datafilename = new char[20];
      sprintf(datafilename,"intdata/f%ld",n);
      data.open(datafilename);;
      delete datafilename;
    }
 for(long i=inf-1; i<nnf; i++)
   {
     if(detail) 
       cout<<"\n"<<"Form number "<<i+1<<"\n";
     else cout<<(i+1)<<" ";
     newform& nfi = nf.nflist[i];
     lfchi lx(&nf, &nfi);
     long s = nfi.sfe;
     long dp0 = 1+p0-nfi.aplist[nf.npdivs];
     rational loverp = nf.nflist[i].loverp;  if(num(loverp)<0) loverp=-loverp;
     bigfloat x0=to_bigfloat(10), y0=to_bigfloat(10), ratio; long a,b;
     int first3=1;
     long lplus=1, lminus=1, mplus=1, mminus=1;
     if(num(loverp)!=0) 
       {
	 lx.compute(1);
	 bigfloat y=lx.value();
	 if(abs(y)>0.001)
	   {
	     x0=y/to_bigfloat(loverp); lplus=1; mplus=num(dp0*loverp);
	     loverp/=2;
//	     cout << "Real period = " << x0 << "\n";
	   }
       }
     else
       {
	 lplus=nf.nflist[i].lplus;
	 mplus=nf.nflist[i].mplus;
	 lx.compute(lplus);
	 x0 = lx.scaled_value()/to_bigfloat(mplus);
       }
     lminus=0;
     for(primevar l; (l<pmax); l++)
       {
	 long ell = l; 
	 if(ell%4!=3) continue;
	 if(legendre(-n,l)!=s) continue;  // i.e. skip to next l
	 lx.compute(ell);
	 bigfloat y = lx.scaled_value();
	 if(first3&&abs(y)>0.001) {y0=y; first3=0; lminus=ell;}
	 if(lminus==ell) {ratio=1; a=b=1;}
	 else
	   {
	     if(first3) ratio = 0; else ratio = y/y0;
	     ratapprox(ratio, a, b);
	     if(b<2501) mminus=lcm(mminus,b);  // avoid spurious values
	   }
	 if(detail) 
	   cout << "l = "<<ell<< "\tRatio- = " << ratio << " =~= "<<a<<"/"<<b<<"\n";
       } // end of primes loop
     if(first3)
       {
	 cout<<"No suitable twisting primes found under " << pmax << ".\n";
	 cout<<"Rerun with a bigger bound, possibly increasing the number of eigs also."<<endl;
       }
     if(detail) cout << "lplus = " << lplus << "; mplus = " << mplus << "\n";
     if(detail) cout << "lminus = "<< lminus <<"; mminus = "<< mminus<< "\n";
     y0/=to_bigfloat(mminus);
     if(detail) cout << "Minimal periods found: x0 = " << x0 << ", y0 = " << y0 << "\n";
     long maxny=maxn;
     long type, ny; bigcomplex w1, w2; bigcomplex c4, c6;
     int notfound=1;
     for(ny=1; notfound&&(ny<=maxny); ny++)
       {
	 for(type=2; (type>0); type--) 
	   {
	     if(type==2){w1=bigcomplex(x0,to_bigfloat(0)); w2=bigcomplex(to_bigfloat(0),y0/to_bigfloat(ny));}
	     else {w1=bigcomplex(2*x0,to_bigfloat(0)); w2=bigcomplex(x0,y0/to_bigfloat(ny));}
	     bigcomplex tau=normalize(w1,w2);
	     getc4c6(w1,w2,c4,c6);
	     bigint ic4 = fac*Iround(real(c4)/fac);
	     bigint ic6 = fac6*Iround(real(c6)/fac6);
	     int close=abs(I2double(ic4)-real(c4))<0.001;
	     int validc4c6 = valid_invariants(ic4,ic6);
	     if((validc4c6||close)&&detail)
	       {
		 cout << "type = " << type << ", ny = " << ny << "\n"; 
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
		     cout << "N = " << cond << "  ";
		     cout << "(lminus = "<<lminus<<", mminus = "<<mminus*ny<<", type = "<<type<<")\n";
		     if(notfound&&dump)
		       {
			 long xmminus = mminus*ny;
			 data << "0\t" << s << "\t" << num(loverp) << "\t" << den(loverp) << "\t" << type << "\n";
			 data << "\t0\t" << dp0 << "\t" << lplus << "\t" << mplus << "\t" << lminus << "\t" << xmminus << "\n";
			 data << "\t0\t0\t0\t0\t0\t0\n";
		       }
		     notfound = 0;
		   }
	       }
	   } // end of type loop
       } // end of ny loop
     if(notfound) 
       {
	 cout << "No curve found: try again with higher limits for maxny" << "\n";
	 if(dump)
	   {
	     data << "0\t" << s << "\t" << num(loverp) << "\t" << den(loverp) << "\t1\n";
	     data << "\t0\t" << dp0 << "\t" << lplus << "\t" << mplus << "\t" << lminus << "\t" << mminus << "\n";
	     data << "\t0\t0\t0\t0\t0\t0\n";
	   }
       }
}       // end of forms loop
if(dump) data.close();
}       // end of if(n)
}       // end of while()
}       // end of main()


void ratapprox(bigfloat x, long& a, long& b)
{
  long c, x0, x1, x2, y0, y1, y2;
  bigfloat xx, diff, eps = to_bigfloat(1.0e-6);
  xx = x; x0 = 0; x1 = 1; y0 = 1; y1 = 0;
//cout<<"ratapprox("<<x<<"): \n";
  diff = 1; x2=y2=c=1;
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

