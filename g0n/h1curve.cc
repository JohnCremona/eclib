// FILE H1CURVE.CC: Program to list curves
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
#include "compproc.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "curve.h"     //from qcurves
#include "cperiods.h"     //from qcurves, for computing conductors
#include "newforms.h"
#include "periods.h"

#ifndef SINGLE   // so Makefile can override
#define AUTOLOOP
#endif

// If this is defined, the fisc6 output is placed in subdirectory
// fixc6 in file fixc6/fixc6.N (one per level); otherwise it is all
// put in ./fixc6.extra
//#define FIXC6_OUTPUT_TO_SUBDIR

#define BOOKORDER       // if defined, sorts newforms/curves into order
                        // in the Book (relevant up to 500 only)
#include "curvesort.cc"

vector<pair<int,int> > bad_ones; // holds bad (n,i) list

int checkap(const level* iN, const newform& nf, CurveRed& CR, long pmax=100);

int main(void)
{
  set_precision("Enter number of decimal places");
 int verb=0; 
#ifdef SINGLE
 verb=0;
#else
 cout<<"See detail? "; cin>>verb;
#endif
 int limit,n=1; 
 char* code = new char[20];
#ifdef AUTOLOOP
 cout<<"Enter first and last N: ";cin>>n>>limit; 
 n--; cout<<endl;
 cout<<endl<<"Table of curves computed from newforms via periods"<<endl;
#ifdef BOOKORDER
 cout<<"(reordered to agree with Book for levels up to 1000)"<<endl;
#endif
 if(!verb)
   {
     cout << "\nN \t# \t";
     cout<<"[a1,a2,a3,a4,a6]";
     cout<<"\t\tConductor\n";
   }
 while (n<limit) { n++;
#else
 while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
 int usedata=1;
 int plus=1,cuspidal=0;
 newforms nf(n,plus,cuspidal,verb);
 int noldap=25;
 nf.createfromdata(noldap,0); // do not create from scratch if data absent
 // nf.createfromolddata();
 //nf.output_to_file();
 int nnf = nf.n1ds; 
 int inf = 1; 
#ifndef SINGLE
 if(verb>1) nf.display();
#else
// if(nnf>1) 
   {
     cout << "Enter form number (between 1 and "<<nnf<<"): "; cin>>inf;
     if((inf<1)||(inf>nnf)) 
       {
	 cout << "Not in range!\n"; inf=1; nnf=0;
       }
     else nnf=inf;
   }
#endif

 for(int xi=inf-1; xi<nnf; xi++)
   { int i = xi;
#ifdef BOOKORDER
     i=booknumber0(n,i);
#endif
     codeletter(xi,code);
     if(verb) cout << "\nForm number " << i+1 << ": " << endl;
     else     cout << n << "\t" << code << "\t";
     //#ifdef SINGLE
     if(verb) nf.nflist[i].display();
     //#endif

     bigfloat rperiod;
     Curve C = nf.getcurve(i, -1, rperiod, verb);
     Curvedata CD(C,1);  // The 1 causes minimalization
     if(getdiscr(Curvedata(CD,0))!=getdiscr(CD))
       {
	 cout << "Non-minimal curve = \t" << C << ", minimal curve = \t";
       }
     else if(verb) cout << "Curve = \t";
     cout << (Curve)CD << "\t";
     CurveRed CR(CD);
     bigint nc = getconductor(CR);
     cout << "N = " << getconductor(CR);
     if(n!=nc) 
       {
	 cout<<" ------WRONG CONDUCTOR!";
	 bad_ones.push_back(pair<int,int>(n,i+1));
       }
     else
       {
	 if(!checkap(&nf, nf.nflist[i],  CR))
	   cout<<" ----- a_p do not agree!";
	 //     else
	 //     cout<<" ----- a_p agree for p<100";
       }
     cout<<endl;
#ifdef SINGLE
     bigint c6=getc6(CD),  c4=getc4(CD);
     char* f = new char[20];
#ifdef FIXC6_OUTPUT_TO_SUBDIR
     sprintf(f,"fixc6/fixc6.%ld",n);
#else
     sprintf(f,"fixc6.extra");
#endif
     ofstream xout(f,ios::app);
     delete[] f;
     xout<<" "<<n<<" "<<(i+1)<<" "<<c6<<endl;
     xout.close();
     cout<<"c4: "<<n<<" "<<(i+1)<<" "<<c4<<endl;
     cout<<"c6: "<<n<<" "<<(i+1)<<" "<<c6<<endl;
#endif
     if(verb) cout<<endl;
   }
}       // end of if(n)
}       // end of while()

 if (bad_ones.size()>0) 
   {
     cout<<"\nNumber of bad curves: "<<bad_ones.size()<<endl;
     cout<<"List of bad curves\n";
     for(unsigned int i=0; i<bad_ones.size(); i++)
       cout<<bad_ones[i].first<<" "<<bad_ones[i].second<<"\n";
     cout<<endl;
   }
}       // end of main()

int checkap(const level* iN, const newform& nf, CurveRed& CR, long pmax)
{
  vector<long> aplist = nf.aplist;
  vector<long> primelist = primes(aplist.size());
  unsigned int i;
  bigint ap, p=BIGINT(0);
  int ok=1, ok1;
  for(i=0; (i<aplist.size())&&(p<=pmax); i++)
    {
      p=primelist[i];
      ap=Trace_Frob(CR,p);
      ok1 =  (ap==BIGINT(aplist[i]));
      if(!ok1) cout<<"p="<<p<<": ap(E)="<<ap<<" but ap(f)="<<aplist[i]<<endl;
      ok = ok && ok1;
    }
  return ok;
}
