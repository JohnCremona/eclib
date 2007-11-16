// FILE H1BSDCURISOG.CC: program to compute isogenous curves & BSD data
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
#include <iomanip>
#include "moddata.h"
#include "symb.h"
#include "oldforms.h"
#include "homspace.h"
#include "cperiods.h"     //from qcurves
#include "newforms.h"
#include "periods.h"

#include "points.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "ffmod.h"
#include "divpol.h"
#include "tlss.h"
#include "elog.h"
#include "saturate.h"
#include "sieve_search.h"
#include "mwprocs.h"
#include "isogs.h"

#include "../qrank/descent.h"

//#define AUTOLOOP
//#define POS_RANK_ONLY  1  // for testing
//#define RANK_ZERO_ONLY 1
//#define RANK_TWO_ONLY 1
//#define SQUAREFREE_ONLY

//#define DEBUG_BSD

#include "curvesort.cc"

int main(void)
{
  set_precision("Enter number of decimal places");
 
 long limit, n=1, hlim1=10, hlim2=15; 
 bigint nn;
 int verbose=0;
 char genfile[30];
 char* code = new char[20];
#ifndef RANK_ZERO_ONLY
 cerr << "See detail (0/1)? "; cin >> verbose;
 cerr << "Limits on naive height in point search? "; 
 cerr << "(searches up to first limit on all curves, up to second limit only if rank is deficient after that): ";
 cin>>hlim1>>hlim2;
 cerr << "filename for dumping generators? "; cin >> genfile;
 ofstream genout;
 genout.open(genfile);
 if(!genout.is_open()) {cerr<<"Unable to open file " << genfile << "\n"; abort();}
#endif

#ifdef AUTOLOOP
 cerr<<"Enter first and last N: ";cin>>n>>limit; 
 n--; cerr<<endl;
 while (n<limit) { n++;
#else
 while (n>0) { cerr<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
 nn=n;
 // Temporary code: check that n is square-free
#ifdef SQUAREFREE_ONLY
 longlist plist=pdivs(n);
 int sqfree=1; long i;
 for(i=0; (i<plist.length)&&sqfree; i++) sqfree=(val(plist(i),n)==1);
 if(!sqfree) continue;
#endif

 int plus=1, cuspidal=0;
 newforms nf(n,plus,cuspidal,0); 
 int noldap=25;
 nf.createfromdata(noldap,0); // do not create from scratch if data absent
 long nclasses = nf.n1ds;
 for(int xi=0; xi<nclasses; xi++)
   { int i=xi;
     codeletter(xi,code);
     i=booknumber0(n,i);
     newform& nfi = nf.nflist[i];
     ldash1 x(&nf, &nfi);  
     long r = x.rank();
     int type = nfi.type;
     rational loverp = nfi.loverp;
     bigint nloverp, dloverp;
     nloverp=abs(num(loverp)); 
     dloverp=abs(den(loverp));
     x.compute();
     bigfloat lf1 = abs(x.value()); // = r!*L^{(r)}(f,1) -- note the r! factor!
#ifdef DEBUG_BSD
	 cout<<"\nL^{(r)}(f,1)/r!: " << lf1 << "\n";
#endif
     /*
     switch(r) 
       {
       case 2: lf1/=2; break; 
       case 3: lf1/=6; break; 
       case 4: lf1/=24; break; 
       case 5: lf1/=120; break; 
       }
     */

#ifdef POS_RANK_ONLY
     if(r==0) continue;    // skip this newform/curve
#endif
#ifdef RANK_TWO_ONLY
     if(r!=2) continue;    // skip this newform/curve
#endif
#ifdef RANK_ZERO_ONLY
     if(r!=0) continue;    // skip this newform/curve
#endif
     if(verbose) cout << "Class " << n << code << ": " << "r = " << r << "\n"; 

     bigfloat rperiod;
     Curve C = nf.getcurve(i, -1, rperiod);
     rperiod = abs(rperiod*(type));

     Curvedata CD(C,1);  // The 1 causes minimalization; else we get 
                         // [0,0,0,-27c4,-54c6]

     IsogenyClass icl(CD,verbose);
     icl.grow();
     if(verbose) icl.displaycurves(cout);
     vector<CurveRed> clist = icl.getcurves();
     long ic, nc = clist.size();

     for(ic=0; ic<nc; ic++)
       {
	 Curvedata & CDi = clist[ic];
	 if(verbose)cout << "CDi = " << CDi << "\t";
	 long nt = CDi.get_ntorsion();
//	 cout.form("%4d %s %2d ",n,code,ic+1);
#ifdef DEBUG_BSD
	 cout<<"\n===========================\n";
#endif
	 cout << setw(6) << n << " " << code << " " << setw(3) << (ic+1) << " ";
	 cout << (Curve)CDi;
	 cout << "\t" << r << "\t" << nt << "\t";

	 CurveRed CR(CDi);
#ifdef DEBUG_BSD
	 cout<<endl;
	 CR.display(cout);
#endif
	 long pcp = prodcp(CR);
#ifdef DEBUG_BSD
	 cout<<"\nProduct of cp: ";
#endif
	 cout << pcp << "\t";

	 bigint cond = getconductor(CR);

	 Cperiods CP(CDi); bigcomplex wR,wRI;
	 CP.getwRI(wR,wRI);
#ifdef DEBUG_BSD
	 cout<<"\nPeriods:"<<CP<<endl;
#endif
	 bigfloat rperiod_i = getconncomp(CDi)*real(wR);
#ifdef DEBUG_BSD
	 cout<<"\nReal Period*: ";
#endif
	 cout << rperiod_i << "\t";
#ifdef DEBUG_BSD
	 cout<<"\nL(f,1): ";
#endif
	 cout << lf1 << "\t";
	 bigfloat loverp_i = lf1/rperiod_i;
#ifdef DEBUG_BSD
	 cout<<"\nL(f,1)/RP*: ";
	 cout << loverp_i << "\t";
#endif
	 bigfloat RS = (loverp_i*to_bigfloat(nt*nt)) / (to_bigfloat(pcp));
#ifdef DEBUG_BSD
	 cout<<"\nRS: ";
#endif
	 //	 cout << RS << "\t";

	 if(nn!=cond) {cout<<"Bad curve ***!!!***\n";}
	 else {
	   if(r==0)
	     {
#ifdef DEBUG_BSD
// N.B. The following code does give some errors!  The "algebraic"
// code for L/P has a bug which has not been fixed, so we use the
// analytic value which is ok.
	       if(ic==0) 
		 {
		   bigfloat diff = I2bigfloat(dloverp)*loverp_i-I2bigfloat(nloverp);
		   if(abs(diff)>0.01)
		     {
		       cout<<"analytic  L/P = "<<loverp_i<<endl;
		       cout<<"algebraic L/P = "<<loverp<<endl;
		     }
		 }
#endif
	       cout << "1\t";               // regulator
	       bigfloat S = RS;
	       long roundS = I2long(Iround(S));
	       if(abs(S-roundS)>0.1) cout << S <<  " ***!!!***";
	       else
		 {
		   cout << roundS;
		   long rootS=(long)(sqrt((double)roundS)+0.1);
		   int squareS=(roundS==rootS*rootS);
		   if(!squareS) cout << " ***!!!***";
		 }
	       cout << endl;
	     }
	   else
	     {
	       if(verbose) cout << "\nRS = " << RS << endl;
	       
	       mw mwbasis(&CDi,verbose,1,r); // stop when rank r is reached
	       mwbasis.search(to_bigfloat(hlim1));
	       long mwr = mwbasis.getrank();
//	       if((mwr<r) && (((nt%2)==0)||(r>1)) ) // 2-descent for rank 2
	       if((mwr<r) && (((nt%2)==0)) ) // no 2-descent for rank 2
		 {
		   if(verbose)
		     {
		       cout<<"Shortfall in rank ("<<mwr<<"<"<<r<<"), doing 2-descent on "<<(Curve)CDi<<"..."<<endl;
		     }
		   long lim1 = 20; // quartic search bound 1 (naive, non-log)
		   long lim2 = 10; // quartic search bound 2 (sieve, log)
		   two_descent two_d(&CDi, 0, 0, lim1, lim2);
		   two_d.saturate(0); // processes 2-desc points only
		   vector<Point> shortlist = two_d.getbasis();
		   mwbasis.process(shortlist);
		   mwr = mwbasis.getrank();
		   if(verbose)
		     cout<<"After 2-descent, rank of points is "<<mwr<<endl;
		 }
	       if(mwr<r)
		 {
		   if(verbose)
		     {
		       cout<<"Shortfall in rank ("<<mwr<<"<"<<r<<"), doing second search..."<<endl;
		     }
		   mwbasis.search(to_bigfloat(hlim2));
		   mwr = mwbasis.getrank();
		 }
       
	       bigint index; vector<long> unsat;
	       int sat_ok = mwbasis.saturate(index,unsat);
	       if(!sat_ok) 
		 {
		   cerr << "saturation possibly incomplete at primes " << unsat << "\n";
		 }
	       bigfloat reg = mwbasis.regulator();
	       vector<Point> gens = mwbasis.getbasis();
	       if(r!=mwr)
		 {
		   cout << "RS = " << RS <<"\t";
		   cout<<"Warning: points found have rank " << mwr 
		       << ", not " << r << ", with hlim2="<<hlim2<<endl;
		 }
	       else {
		 if(verbose)
		   {
		     cout << "Basis of points found: " << gens << endl;
		     cout << "Regulator = " << reg << endl;
		   }
		 cout << reg << "\t";
		 bigfloat S = RS/reg;
		 cout << S;
		 long roundS = I2long(Iround(S));
		 if (abs(S-roundS)>0.1) cout << " ***!!!***";
		 else
		 {
		     long rootS=(long)(sqrt((double)roundS)+0.1);
		     int squareS=(roundS==rootS*rootS);
		     if(!squareS) cout << " ***!!!***";
		 }
		 cout << endl;
	       }
// Dump to genfile:
	       genout << n << "\t" << code << "\t" << ic+1;
	       genout << "\t" << (Curve)CDi ;
	       genout << "\t" << r ;
	       for(int ip=0; ip<mwr; ip++) genout << "\t" << gens[ip];
	       genout << endl;
	     }
	 }
       } // end of curves loop
   }     // end of classes loop

}       // end of if(n)
}       // end of while()
#ifndef RANK_ZERO_ONLY
genout.close();
#endif 
}       // end of main()
