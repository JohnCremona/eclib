// G2.CC  -- Genus 2 eigenvalues (preliminary version)
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

#include "msubspace.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"


void sqfsplit(long d, long& d0, long& d1);

int main(void)
{
 cout << "Program g2.  Using METHOD = " << METHOD << " to find newforms" << endl;
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 int n=1; 
 int plus=1;
 int verbose=1;
 cout << "Long output? (0/1) "; cin >> verbose;
 // cout << "Plus space (0/1)? "; cin >> plus;
 while (cout<<"Enter level: ", cin>>n, n>0)
   {
     cout << ">>>Level " << n << "\t";
     homspace hplus(n,plus,0);
     int genus = hplus.h1dim();
     cout << "Dimension = " << genus << "\n";
     longlist badprimes = hplus.plist;
     int nq = badprimes.length; int firstq=0;  // =0 for all W's
     if (genus>0)
       {
	 mmatrix id = idmat(genus), tp, tp2, tp0, m;
	 msubspace s;
	 bigint a, b; long d, dims, np=10,ip=0; 
	 primevar pr;
	 while (n%pr==0) pr++;
	 long p=pr;
	 cout << "Computing T_p for p = " << p << "\n";
	 tp = hplus.newheckeop(p,0);
	 tp2= tp*tp;
	 if(verbose) {cout<<"Tp         = \n"; cout<<tp;}
	 if(verbose) {cout<<"Tp squared = \n"; cout<<tp2;}
	 int looking=1;
	 while(looking)
	   {
	     cout<<"Enter coeffs a, b for quadratic polynomial X^2+a*X+b: ";
	     cin>>a>>b;
	     m = tp2+a*tp+b*id;
	     if(verbose) 
	       {cout<<"Matrix of T_p^2+a*T_p+b = \n"; cout<<m;}
	     s = kernel(m);
	     dims = dim(s);
	     looking = (dims!=2);
	     if(looking)
	       {
		 cout<<"Eigendimension = " << dims << ", not 2!  Try again.\n";
	       }
	   }
	 d=I2long(denom(s));
	 tp0 = restrict(tp,s);
	 if(verbose)
	   {
	     cout<<"Eigenspace S has dimension = 2";
	     if(d!=1) cout<<", denom = "<<d;
	     cout<<endl;
	     cout<<"Basis of S = \n"; cout<<basis(s);
	     cout<<"Restriction of T("<<p<<") to S = ";
	     if(d>1) cout<<"1/"<<d<<" * ";
	     cout<<"\n";
	     cout<<tp0;
	     cout<<endl;
	   }
// Compute W_q:
	 long iq;
	 for(iq=0; iq<nq; iq++)
	   {
	     long q = badprimes[iq];
	     mmatrix wq = restrict(hplus.heckeop(q,0),s);
	     if(verbose)
	       {
		 cout<<"Restriction of W("<<q<<") to S = ";
		 if(d>1) cout<<"1/"<<d<<" * ";
		 cout<<"\n";
		 cout<<wq;
	       }
	     long eq = I2long(wq(1,1))/d;
	     cout << "Eigenvalue of W("<<q<<") = " << eq << endl;
	   }

	 cout<<"How many ap? "; cin>>np;
	 long am=I2long(tp0(1,1)), bm=I2long(tp0(1,2)), cm=I2long(tp0(2,1)), dm=I2long(tp0(2,2));
	 long del = ((am-dm)*(am-dm)+4*bm*cm)/(d*d);
	 long del0, del1;
	 sqfsplit(del, del0, del1);
	 int onemod4 = (del0%4)==1;
	 cout << "\nEigenvalues in Q(sqrt("<<del0<<")), ";
	 if(onemod4)
	   cout << "with respect to basis [1, sqrt("<<del0<<")]/2:\n\n";
	 else
	   cout << "with respect to basis [1, sqrt("<<del0<<")]:\n\n";
	 cout << "p:\ts\tt\n";
	 for (; pr.ok()&&ip<np; pr++, ip++)
	   {
	     while (n%pr==0) pr++;
	     p=pr;
	     if(p<50) tp = hplus.newheckeop(p,0);
	     else tp = hplus.heckeop(p,0);
	     tp0 = restrict(tp,s);
	     if(verbose)
	       {
		 cout<<"Restriction of T("<<p<<") to S = ";
		 if(d>1) cout<<"1/"<<d<<" * ";
		 cout<<"\n";
		 cout<<tp0;
	       }
	     long s = I2long(tp0(1,1)+tp0(2,2))/d;
	     long t = (I2long(tp0(2,1))*del1)/cm;
	     if(!onemod4) {s/=2; t/=2;}
	     cout << p << ":\t" << s << "\t" << t << "\n";
	   } 
       }      // end of if(genus>0)
     
   }       // end of while()
 abort();
}       // end of main()

// Function to write (positive) integer d as d0*d1^2 with d0 square-free

void sqfsplit(long d, long& d0, long& d1)
{
  longlist plist = pdivs(d);
  long j,p,ep,ip,np=plist.length;
  d0=d1=1;
  for(ip=0; ip<np; ip++)
    {
      p = plist[ip];
      ep=val(p,d);
      if(odd(ep)) {d0*=p; ep--;}
      ep/=2;
      for(j=0; j<ep; j++) {d1*=p;}
    }
  if(d!=d0*d1*d1)
    cout<<"Error in sqfsplit("<<d<<"): d0="<<d0<<", d1="<<d1<<endl;
}
