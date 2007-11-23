// tmrank.cc:  test program for mwrank: read from (e.g.) tmrank.in
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
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
                                                   
#include "mwprocs.h"
#include "descent.h"
#include "version.h"
#include "timer.h"

bigint a1,a2,a3,a4,a6;
Curve C;
Curvedata CD;
int verbose;
long hlimq=5; 
long naux=-1;

//#define SELMER_ONLY

int getcurve(void)
{
  if (verbose) cout << "Enter curve coefficients a1,a2,a3,a4,a6 ?" << "\n";
  cin >> a1 >> a2 >> a3 >> a4 >> a6;
  C = Curve(a1,a2,a3,a4,a6);
  CD = Curvedata(C,1);     // "1" means minimise
  return (a1!=0||a2!=0||a3!=0||a4!=0||a6!=0);
}

int main()
{
  show_version();
  set_precision(20);
  init_time();
  int second_descent=1;
  int selmer_only=0;
#ifdef SELMER_ONLY
  selmer_only=1;
#endif  
  cout << "Verbose mode? (0/1)\n"; cin >> verbose;
  initprimes(string("PRIMES").c_str(),verbose);
//  cout << "Number of sieving primes?\n"; cin >> naux;
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  while (getcurve())
    {
      int filerank;  // for testing
      cin >> filerank;
      if(verbose) cout<<"\n\n";
      cout << "Curve "<< C << " :\t";
      if (verbose) cout << endl; else cout<<flush;
      long rank;
 
      start_time();
      two_descent two_d(&CD, verbose, selmer_only, 
			20, hlimq, 
			naux, second_descent);
      stop_time();
 
      if (two_d.ok())
        {   
#ifdef SELMER_ONLY
	  rank = two_d.getselmer();
	  cout << "(r2) Selmer rank = " << rank << "\t" << flush;
#else
	  rank=two_d.getrank();
          cout << "Rank = " << rank << "\t" << flush;
#endif
	  show_time();cout<<endl;
          if (rank!=filerank) cout << "Wrong! rank of "<<C
				   <<" should be " << filerank 
				   << " not "<<rank<<endl;
        }
      else cout << "Failed to compute rank\n";
    }
}
