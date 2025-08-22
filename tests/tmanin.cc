// TMANIN.CC: Program for finding newforms & computing ap
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
//
#include <eclib/interface.h>
#include <eclib/timer.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/oldforms.h>
#include <eclib/cperiods.h>
#include <eclib/newforms.h>

#define AUTOLOOP
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output

int main(void)
{
 init_time();
 start_time();
 long n=110, stopp;
 int output, verbose, sign=1;

 cout << "Program tmanin." << endl;
 cerr << "Verbose output? "; cin>>verbose;
 cerr << "How many primes for Hecke eigenvalues? ";
 cin  >> stopp; cerr << endl;
 output=1;
 cerr << "Output Hecke eigenvalues to file? (0/1) ";  cin >> output;
 cerr << "Sign? (-1/0/1) ";  cin >> sign;
#ifdef AUTOLOOP
 long limit;
 cerr<<"Enter first and last N: ";cin>>n>>limit; n--;
 while (n<limit) { n++;
#else
     while (n>1) { cerr<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
  cout << ">>>Level " << n;
  if(verbose)cout<<endl; else cout<< ":\t";
  newforms nf(n,verbose); 
  int noldap=25;
  nf.createfromscratch(sign,noldap);
#ifdef LMFDB_ORDER
  nf.sort_into_LMFDB_label_order();
#else
  nf.sort_into_Cremona_label_order();
#endif
  nf.make_projcoord(); // needed for when we add more ap
  if(verbose>1) nf.display();
  else          cout << nf.n1ds << " newform(s) found.";
  if(verbose&&nf.n1ds>0) 
    cout<<"\nComputing ap for primes up to "<<prime_number(stopp)<<endl;
  nf.addap(stopp);
  if(output)
    {
      nf.output_to_file(1,0); // full nf data
      nf.output_to_file(1,1); // small nf data
    }
  //  }
}       // end of if(n)
}       // end of while()
stop_time(); cout<<endl;
show_time(cerr); cerr<<endl;
}       // end of main()
