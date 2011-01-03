// FILE nfhpmcurve.cc main newform- and curve-finding program
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

#define AUTOLOOP

int main(void)
{
 init_time();
 start_time();
 long n=110, stopp, stopp0; 
 long prec0=20;
 int output, verbose;

 cout << "Program nfhpcurve.  Using METHOD = " << METHOD 
      << " to find newforms" << endl;
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 cout << "Verbose output? "; cin>>verbose;
 cout << "How many primes for Hecke eigenvalues? ";
 cin  >> stopp0; cout << endl;
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
  newforms nf(n,verbose); 
  int noldap=25; // stopp0 must be at least this big!
  if (stopp0<noldap) stopp0=noldap;
  nf.createfromscratch(1,noldap);
  if(verbose) nf.display();
  else          cout << nf.n1ds << " newform(s) found."<<endl;
  stopp=stopp0;
  // Must make sure that this is big enough to catch the bad primes...?
  nf.addap(stopp);

  long nnf = nf.n1ds, inf; 
  if(nnf==0)
    {
      cout << "No newforms.\n";
      cout << "Finished level "<<n<<endl;
      if(output) nf.output_to_file();
      continue;
    }

  // Thus far, as in tmanin

  if(verbose)
    {
      cout << endl;
      cout << "Creating minus space..."<<endl;
    }
  nf.set_sign(-1);
  if(verbose)
    {
      cout<<"sign = "<<nf.get_sign()<<endl;
    }
  nf.makebases(1);
  if(verbose)
    {
      cout << "...finished creating minus data."<<endl<<endl;
      nf.display();
      cout << "Filling in data for newforms..."<<endl;
    }
  nf.set_sign(0);
  nf.merge();
  if(verbose) 
    {
      cout << "...finished filling in full data."<<endl;
      cout << "Updated newforms: ";
      nf.display();
    }
  if(output) nf.output_to_file();

  // Now we compute the curves
  cout<<"Computing "<<nnf<<" curves...\n";

  long fac=nf.sqfac;
  if(verbose&&(fac>1)) cout<<"c4 factor " << fac << endl;

  stopp = stopp0; // will be increased if necessary
  set_precision(prec0);
  bigfloat rperiod;
  vector<int> forms;
  for(inf=0; inf<nnf; inf++) forms.push_back(inf);

  while(forms.size()>0)
    {
      forms = nf.showcurves(forms,verbose);
      if(forms.size()==0)
	{
	  cout<<"All curves found successfully!"<<endl;
	  cout << "Finished level "<<n<<endl;
	}
      else
	{      
	  cout<<forms.size()<<" curve(s) missing: ";
	  for(vector<int>::const_iterator inf=forms.begin(); inf!=forms.end(); inf++)
	    cout<<(*inf+1)<<" ";
	  cout<<endl;
	  int newstopp;
	  if(stopp<500) 
	    newstopp=2*stopp;
	  else
	    newstopp=stopp+500;
	  if(newstopp>32000)
	    {
	      cout<<"Cannot compute more ap, something must be wrong in newform data"<<endl;
	      abort();
	    }
	  cout<<"Computing some more ap: from "<<stopp+1<<" to "
	      <<newstopp<<"..."<<endl;
	  nf.set_sign(1);
	  nf.makeh1(1);
	  nf.addap(newstopp);
	  stopp=newstopp;
	  if(output) nf.output_to_file();
#ifdef MPFP
	  if(decimal_precision()<50) 
	    {
	      set_precision(decimal_precision()+5);
	      cout << "Now working to "<<decimal_precision()<<" decimal places" << endl;
	    }
#endif
	} 
    }
 }       // end of if(n>0)
     }  // end of while(n>0) or while(n<limit)
 }       // end of main()

