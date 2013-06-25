// FILE nfhpcurve_prof.cc profiling nfhpcurve program
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

#include <eclib/interface.h>
#include <eclib/timer.h>
#include <eclib/xsplit.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/oldforms.h>
#include <eclib/cperiods.h>     //from qcurves, for computing conductors
#include <eclib/newforms.h>
#include <eclib/periods.h>
#include <eclib/pcprocs.h>

#define AUTOLOOP
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output
#define PROFILE

#define MAXNY 10000
#define MAXD 10
#define MAXNAP 20000

int main(void) {

#if PROFILE
  // Initiate global timers
  init_time();
  
  // Open timer log file
  ofstream profile.open('times.dat',std::ios::out);
  if( profile == NULL ) {
    cout << "Cannot open output file ... exiting" << endl;
    exit( EXIT_FAILURE );
  }

  // Initialise start timer
  start_time();
#endif

  long n=1, stopp; // have a dud (but positive) value of n here to avoid mishaps
  int output, verbose;
  long maxn = MAXNY;
  long dmax = MAXD;

  cout << "Program nfhpcurve.  Using METHOD = " << METHOD 
       << " to find newforms" << endl;
  set_precision(25);
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
  if (n>0) {
    cout << "\n>>>Level " << n;
    // Temporary code to skip non-square-free levels
    //
    //  if(!is_squarefree(n)) {cout<<" --not square-free, skipping!\n"; continue;}
    //
    if(verbose)cout<<endl; else cout<< ":\t"<<flush;
    if(!is_valid_conductor(n)) {
      cout<<"Not a valid conductor!"<<endl;
      output_to_file_no_newforms(n);
      cout << "Finished level "<<n<<endl;
      continue;
    }
  int plus=1, cuspidal=0;
  newforms nf(n,verbose); 
  int noldap=25;
  nf.createfromscratch(plus,noldap);
#ifdef LMFDB_ORDER
  nf.sort();
  nf.make_projcoord(); // needed for when we add more ap
#endif
  if( verbose ) nf.display();
  else          cout << nf.n1ds << " newform(s) found."<<endl;
  nf.addap(stopp);

  long nnf = nf.n1ds, inf; 
  if( nnf == 0 ) {
    if(output) nf.output_to_file();
      cout << "No newforms.\n";
      cout << "Finished level "<<n<<endl;
      continue;
    }

  // Thus far, as in tmanin

  // Now we search for curves as in pcurve.cc

  if(output) nf.output_to_file();

  long rootn=(long)(sqrt((double)n)+0.1); 
  int squarelevel=(n==rootn*rootn);
  cout<<"Computing "<<nnf<<" curves...\n";

  long fac=nf.sqfac;
  if(verbose&&(fac>1)) cout<<"c4 factor " << fac << endl;

  int* success=new int[nnf];
  long nsucc=0;
  for(inf=0; inf<nnf; inf++) success[inf]=0;

  while(nsucc<nnf){

  for(inf=0; inf<nnf; inf++) {
     if(success[inf]) continue;
     if(verbose) 
       cout<<"\n"<<"Form number "<<inf+1<<"\n";
     else cout<<(inf+1)<<" ";
     newform* nfi = &((nf.nflist)[inf]);

     int rp_known=0;
     // nfi->dotplus=1;  // This will have been set correctly
     nfi->dotminus=1;
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

         bigfloat rperiod;
         Curve C = nf.getcurve(inf, -1, rperiod);
         if (C.isnull())
           {
             success[inf]=0; nsucc--;
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
      if(output) nf.output_to_file();
      int newstopp=stopp;
      if (newstopp<1000)
        newstopp+=500;
      else
        newstopp+=1000;
      if(newstopp>MAXNAP)
      {
	  cout<<"Cannot compute more ap, something must be wrong in newform data"<<endl;
	  abort();
      }
      cout<<"Computing some more ap: from "<<stopp+1<<" to "
	  <<newstopp<<"..."<<endl;
      nf.addap(newstopp);
      if(output) nf.output_to_file();
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
  }   // ends while(nsucc<nnf)
  delete[]success;
  
}       // end of if(n>0)
     
#if PROFILE
// Stop timer and write to file
stop_timer();
profile << "Level " << n << ": " << get_time() << "seconds" << end;'
#endif

     }  // end of while(n>0) or while(n<limit)

#if PROFILE
// Close output file
profile.close()
#endif

 }       // end of main()

