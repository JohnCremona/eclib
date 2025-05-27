// FILE H1CURVE.CC: Program to list curves
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
#include <eclib/compproc.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/oldforms.h>
#include <eclib/curve.h>     //from qcurves
#include <eclib/cperiods.h>     //from qcurves, for computing conductors
#include <eclib/newforms.h>
#include <eclib/periods.h>

#ifndef SINGLE   // so Makefile can override
#define AUTOLOOP
#endif

#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output
                          // otherwise, sorts newforms into Cremona order before output

// If this is defined, the fisc6 output is placed in subdirectory
// fixc6 in file fixc6/fixc6.N (one per level); otherwise it is all
// put in ./fixc6.extra
//#define FIXC6_OUTPUT_TO_SUBDIR

#include <eclib/curvesort.h>

vector<pair<int,int> > bad_ones; // holds bad (n,i) list

int checkap(const level* iN, const newform& nf, CurveRed& CR, long pmax=100);

int main(void)
{
  int prec0 = 100;
  int maxprec = 500;
  int prec = prec0;
  int delta_prec = 50;
  set_precision(prec);
 int verb=0;
#ifdef SINGLE
 verb=1;
#else
 cout<<"See detail? "; cin>>verb;
#endif
 long n=1;
#ifdef AUTOLOOP
 int limit;
 cout<<"Enter first and last N: ";cin>>n>>limit;
 n--; cout<<endl;
 cout<<endl<<"Table of curves computed from newforms via periods"<<endl;
#ifdef LMFDB_ORDER
 cout<<"(in LMFDB order)"<<endl;
#else
 cout<<"(in Cremona order)"<<endl;
#endif
 if(!verb)
   {
     cout << "\nN \t# \t";
     cout<<"[a1,a2,a3,a4,a6]";
     cout<<"\t\tConductor\n";
   }
 while (n<limit) { n++;
#else
 while (n>1) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
 newforms nf(n,verb);
 int noldap=25;
 nf.createfromdata(1,noldap,0); // do not create from scratch if data absent
#ifdef LMFDB_ORDER
  nf.sort_into_LMFDB_label_order();
#else
  nf.sort_into_Cremona_label_order();
#endif
 int nnf = nf.n1ds;
 int inf = 0;
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
     else
       {
         nnf=inf;
         inf--;
       }
   }
#endif

 for(int xi=inf; xi<nnf; xi++)
   { int i = xi;
     if(verb) cout << "\nForm number " << i+1 << ": " << endl;
     else     cout << n << "\t" << codeletter(xi) << "\t";
     //#ifdef SINGLE
     if(verb) nf.nflist[i].display();
     //#endif

     bigfloat rperiod;
     Curve C;
     Curvedata CD;
     CurveRed CR;
     bigint nc;
     prec = prec0-delta_prec;
     set_precision(prec);
     C = Curve();
     while (C.isnull() && (prec<maxprec))
       {
	 prec += delta_prec;
	 set_precision(prec);
	 C = nf.getcurve(i, -1, rperiod, verb);
         if (!C.isnull())
           {
             CD = Curvedata(C,1);  // The 1 causes minimalization
             CR = CurveRed(CD);
             nc = getconductor(CR);
             if(n!=nc) C = Curve();  // wrong conductor; reset to null curve
           }
       }
     if(C.isnull())
       {
         cout << "bad curve, even after trying precision up to "<<maxprec<<endl;
	 bad_ones.push_back(pair<int,int>(n,i+1));
       }
     else
       {
         if(getdiscr(Curvedata(C,0))!=getdiscr(CD))
           {
             cout << "Non-minimal curve = \t" << C << ", minimal curve = \t";
           }
         else if(verb) cout << "Curve = \t";
         cout << (Curve)CD << "\t";
         cout << "N = " << nc;
         if(n!=nc)
           {
             cout<<" ------WRONG CONDUCTOR!";
           }
         else
           {
             if(!checkap(&nf, nf.nflist[i],  CR))
               cout<<" ----- a_p do not agree!";
             //     else
             //     cout<<" ----- a_p agree for p<100";
           }
         cout<<endl;
       }
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
  long p=1;
  int ok=1;
  for(i=0; (i<aplist.size())&&(p<=pmax); i++)
    {
      p = primelist[i];
      long ap = CR.ap(p);
      if (ap!=aplist[i])
        {
          cout<<"p="<<p<<": ap(E)="<<ap<<" but ap(f)="<<aplist[i]<<endl;
          ok = 0;
        }
    }
  return ok;
}
