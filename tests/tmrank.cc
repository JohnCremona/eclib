// tmrank.cc:  test program for mwrank: read from (e.g.) tmrank.in
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
                                                   
#include <eclib/mwprocs.h>
#include <eclib/descent.h>
#include <eclib/version.h>
#include <eclib/timer.h>

bigint a1,a2,a3,a4,a6;
Curve C;
Curvedata CD;
int verbose=0;
long hlimq=5; 
long naux=-1;

//#define SELMER_ONLY
//#define TIMINGS

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
  show_version(cerr);
  set_precision(50);
#ifdef TIMINGS
  init_time();
#endif
  int second_descent=1;
  int selmer_only=0;
#ifdef SELMER_ONLY
  selmer_only=1;
#endif  
  //  cout << "Verbose mode? (0/1)\n"; cin >> verbose;
  initprimes("PRIMES",verbose);
//  cout << "Number of sieving primes?\n"; cin >> naux;
  cin.flags( cin.flags() | ios::dec );
  while (getcurve())
    {
      int filerank;  // for testing
      cin >> filerank;
      if(verbose) cout<<"\n\n";
      cout << "Curve "<< C << " :\t";
      if (verbose) cout << endl; else cout<<flush;
      long rank;

#ifdef TIMINGS
      start_time();
#endif
      two_descent two_d(&CD, verbose, selmer_only,
			20, hlimq,
			naux, second_descent);
#ifdef TIMINGS
      stop_time();
#endif

      if (two_d.ok())
        {
#ifdef SELMER_ONLY
	  rank = two_d.getselmer();
	  cout << "(r2) Selmer rank = " << rank << "\t" << flush;
#else
	  rank=two_d.getrank();
          cout << "Rank = " << rank << "  ";
#endif
#ifdef TIMINGS
	  show_time();
#endif
          cout<<endl;
          if (rank!=filerank) cout << "Wrong! rank of "<<C
                                   <<" should be " << filerank
				   << " not "<<rank<<endl;
        }
      else cout << "Failed to compute rank\n";
    }
}
