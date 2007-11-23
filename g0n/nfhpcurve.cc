// FILE nfhpcurve.cc main newform- and curve-finding program
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

#include "interface.h"
#include "timer.h"
#include "xsplit.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "cperiods.h"     //from qcurves, for computing conductors
#include "newforms.h"
#include "periods.h"
#include "pcprocs.h"

//#define AUTOLOOP
#define MAXNY 100
#define MAXD 10

int main(void)
{
 init_time();
 start_time();
 long n=110, stopp; 
 int output, verbose;
 long maxn = MAXNY;
 long dmax = MAXD;

 cout << "Program nfhpcurve.  Using METHOD = " << METHOD 
      << " to find newforms" << endl;
 set_precision(string("Enter number of decimal places").c_str());
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 cout << "Verbose output? "; cin>>verbose;
 cout << "How many primes for Hecke eigenvalues? ";
 cin  >> stopp; cout << endl;
 output=1; 
 cout << "Output newforms to file? (0/1) ";  cin >> output;
#ifdef AUTOLOOP
 int limit;
 cout<<"Enter first and last N: ";cin>>n>>limit; n--;
 while (n<limit) { n++;
#else
     while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
  cout << "\n>>>Level " << n;
  // Temporary code to skip non-square-free levels
  //
  //  if(!is_squarefree(n)) {cout<<" --not square-free, skipping!\n"; continue;}
  //
  if(verbose)cout<<endl; else cout<< ":\t"<<flush;
  if(!is_valid_conductor(n))
    {
      cout<<"Not a valid conductor!"<<endl;
      output_to_file_no_newforms(n);
      cout << "Finished level "<<n<<endl;
      continue;
    }
  int plus=1, cuspidal=0;
  newforms nf(n,plus,cuspidal,verbose); 
  int noldap=25;
  nf.createfromscratch(noldap);
  if(verbose) nf.display();
  else          cout << nf.n1ds << " newform(s) found."<<endl;
  nf.addap(stopp);

  long nnf = nf.n1ds, inf; 
  if(nnf==0)
    {
      if(output) nf.output_to_file();
      cout << "No newforms.\n";
      cout << "Finished level "<<n<<endl;
      continue;
    }

  // Thus far, as in tmanin

  // Now we search for curves as in pcurve.cc


  long rootn=(long)(sqrt((double)n)+0.1); 
  int squarelevel=(n==rootn*rootn);
  cout<<"Computing "<<nnf<<" curves...\n";

  long fac=nf.sqfac;
  if(verbose&&(fac>1)) cout<<"c4 factor " << fac << endl;

  int* success=new int[nnf];
  long nsucc=0;
  for(inf=0; inf<nnf; inf++) success[inf]=0;

  while(nsucc<nnf){

  for(inf=0; inf<nnf; inf++)
   {
     if(success[inf]) continue;
     if(verbose) 
       cout<<"\n"<<"Form number "<<inf+1<<"\n";
     else cout<<(inf+1)<<" ";
     newform* nfi = &((nf.nflist)[inf]);

     int rp_known=0;
     nfi->dotplus=1; nfi->dotminus=1;  // need to reset these in case this is not 1st attempt!
     bigfloat x0=to_bigfloat(10), y0=to_bigfloat(10);
     int have_both = nf.find_matrix( inf, dmax, rp_known, x0, y0);
     if(!have_both)
       cout<<"Problem!  find_matrix() returns 0!"<<endl;

     if(verbose) 
       {
	 cout << "Minimal periods found: x0 = "<<x0<<", y0 = "<<y0<<"\n";
	 cout << "Matrix ("<<nfi->a<<","<<nfi->b<<";"<<n*nfi->c<<","<<nfi->d<<"):\t";
	 cout << "dotplus = "<< nfi->dotplus << ", dotminus = "<< nfi->dotminus<< "\n";
	 cout << "Searching for scaling factors.\n";
       }
    
     long nx, ny;
     long maxnx=maxn; if(rp_known) maxnx=1;
     int found = get_curve(n, fac, maxnx, maxn, x0, y0, nx, ny, nfi->type, verbose);

     if(found)
       {
	 success[inf]=1; nsucc++;
	 nfi->dotplus *= nx;
	 nfi->dotminus *= ny;
	 cout << "[(" <<nfi->a<<","<<nfi->b<<";"<<nfi->c
	      <<","<<nfi->d<<"),"<<nfi->dotplus<<","<<nfi->dotminus
	      <<";"<<nfi->type<<"]"<<endl;

// STEP 4:  We find a suitable twisting prime for the imaginary period

	 if(!squarelevel)
	   {
	     bigfloat y1 = y0/to_bigfloat(ny);
	     int ok = nf.find_lminus(inf, 0, y1);
	     if(!ok)
	       {
		 cout<<"No suitable twisting prime "
		     <<" found for imaginary period!"<<endl;
	       }
	   }
       }
     else cout<<"No curve found"<<endl;

   } // ends loop over newforms

  if(nsucc==nnf)
    {
      cout<<"All curves found successfully!"<<endl;
      cout << "Finished level "<<n<<endl;
      if(output) nf.output_to_file();
    }
  else
    {
      cout<<(nnf-nsucc)<<" curve(s) missing."<<endl;
      int newstopp=stopp+500;
      if(newstopp>32000)
      {
	  cout<<"Cannot compute more ap, something must be wrong in newform data"<<endl;
	  if(output) nf.output_to_file();
	  abort();
      }
      cout<<"Computing some more ap: from "<<stopp+1<<" to "
	  <<newstopp<<"..."<<endl;
      nf.addap(newstopp);
      stopp=newstopp;
      if(maxn<10000) 
	{
	  maxn+=200;
	  cout << "Now working with maxny =  " <<maxn<< endl;
	}
#ifdef MPFP
      if(decimal_precision()<50) 
	{
	  set_precision(decimal_precision()+10);
	  cout << "Now working to "<<decimal_precision()<<" decimal places" << endl;
	}
#endif
    }
  delete[]success;
  }   // ends while(nsucc<nnf)
  
}       // end of if(n>0)
     }  // end of while(n>0) or while(n<limit)
 }       // end of main()

