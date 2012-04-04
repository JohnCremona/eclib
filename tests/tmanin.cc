// TMANIN.CC: Program for finding newforms & computing ap
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
#include "interface.h"
#include "timer.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "cperiods.h"
#include "newforms.h"

#define AUTOLOOP

int main(void)
{
 init_time();
 start_time();
 long n=110, stopp; 
 int output, verbose, sign=1, cuspidal=0;

 cout << "Program tmanin.  Using METHOD = " << METHOD << " to find newforms" << endl;
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 cout << "Verbose output? "; cin>>verbose;
 cout << "How many primes for Hecke eigenvalues? ";
 cin  >> stopp; cout << endl;
 output=1; 
 cout << "Output Hecke eigenvalues to file? (0/1) ";  cin >> output;
 cout << "Sign? (-1/0/1) ";  cin >> sign;
#ifdef AUTOLOOP
 long limit;
 cout<<"Enter first and last N: ";cin>>n>>limit; n--;
 while (n<limit) { n++;
#else
     while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
  cout << "\n>>>Level " << n;
  if(verbose)cout<<endl; else cout<< ":\t";
  newforms nf(n,verbose); 
  int noldap=25;
  nf.createfromscratch(sign,noldap);
  if(verbose>1) nf.display();
  else          cout << nf.n1ds << " newform(s) found."<<endl;
  if(verbose&&nf.n1ds>0) 
    cout<<"Computing ap for primes up to "<<prime_number(stopp)<<endl;
  nf.addap(stopp);
  if(output) nf.output_to_file();
  //  }
}       // end of if(n)
}       // end of while()
stop_time();
show_time(cerr); cerr<<endl;
}       // end of main()
