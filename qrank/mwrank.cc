// mwrank.cc:  program to find Mordell-Weil group via 2-descent and saturation
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
#include "getcurve.h"
#include "timer.h"
#include "options.h"
#include "descent.h"
#include "version.h"

#define MAX_HEIGHT 20 // will never search for points on curve of naive 
                      // logarithmic height greater than this 
#define SEARCH_METHOD 0  // 0 for Stoll sieve (fastest)
                         // 1, 2, 3 for JC's sieve

int main(int argc, char **argv)
{
  mrank_options opt;
  opt.set(argc, argv);
  int verbose = (opt.get_verbose());
  //opt.show();
  if(!opt.get_quiet()) {opt.banner(1);  show_version();}
#if defined(LiDIA_ALL) || defined(NTL_ALL)
  long lidia_precision=opt.get_precision();
  set_precision(lidia_precision);
  if(verbose) cerr<<"Using "<<
#if defined(LiDIA_ALL)
		"LiDIA"
#else
		"NTL"
#endif
		  <<" multiprecision floating point with "
		  <<lidia_precision <<" decimal places.\n";
#else
  cout.precision(15);
  cerr.precision(15);
#endif
  initprimes(string("PRIMES").c_str(),0);
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)

  init_time();

  long sat_bd = opt.get_saturation_bound();
  int selmer_only = (opt.get_selmer_only());
  int second_descent = (opt.get_second_descent());
  long hlimq = (opt.get_hlimq());
  long hlimq0 = hlimq; if(hlimq0!=0) hlimq0=10;
  
  Curvedata CD;
  vector<bigrational> ai;

  //  while ( getcurve(CD,!opt.get_quiet())	 )
  while ( getcurve(ai,!opt.get_quiet())	 )
    {
      cout << "Curve "<< 
	"["<<ai[0]<<","<<ai[1]<<","<<ai[2]<<","<<ai[3]<<","<<ai[4]<<"] :\t";
      start_time();

//  Step 1: 2-descent

      two_descent two_d(ai, verbose, selmer_only, 
			hlimq0, hlimq, 
			opt.get_naux(), second_descent);
      two_d.report_rank();
     

      if((!two_d.ok())||selmer_only) continue;  // to next curve
      if(!verbose&&!opt.get_ptl()&&!opt.get_output_pari()) continue;
             // no point in saturating since no points are output

//  Step 2: saturation

      two_d.saturate(sat_bd);  

//  Step 3: output

      if(verbose||opt.get_ptl())
	two_d.show_gens();
      if(verbose) 
	two_d.show_result_status();

//  Optional Pari interface output:

      if(opt.get_output_pari())
	{
	  two_d.pari_output();
	  cout<<endl;
	}

      stop_time();
      if(verbose) {show_time(cerr); cerr<<endl;}
    }
  cout<<endl;
}


