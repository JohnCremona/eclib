// options.h:   declaration & implementation of class to handle mwrank options
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
 
#include "GetOpt.h"

#define DEFAULT_QUIET 0
#define DEFAULT_VERBOSE 1
#define DEFAULT_PRECISION 15
#define DEFAULT_HLIMQ 10
#define DEFAULT_NAUX 8
#define DEFAULT_HLIMC 0
#define DEFAULT_TRACEEQUIV 0
#define DEFAULT_PTL -99  // if not set manually will be set to same as verbose
#define DEFAULT_OUTPUT_PARI 0
#define DEFAULT_SELMER 0  // if set will do only local tests to get Selmer rank
#define DEFAULT_2ND_DESCENT 1 // if set will do 2nd descent (2-torsion curves)
#define DEFAULT_SAT_BOUND 1000 // bound on saturation primes

class mrank_options {
private:
  int quiet;            // 0/1, controls header output
  int verbose;          // 0-3, controls output verbosity
  long precision;       // 1-\infty, controls LiDIA precision
  long hlimq;           // 1-20, height limit for quartic search
  long naux;            // number of primes used in syzygy sieve
  long hlimc;           // 1-15, height limit for curve search
  int ptl;              // 0/1, controls whether points are output
  int traceequiv;       // 0/1, controls equivalence tracing (debugging only)
  int output_pari;      // 0/1, controls pari-abbreviated output
  int selmer_only;      // 0/1, if set only computes Selmer rank
  int second_descent;   // 0/1, if set does 2nd descent
  int saturation_bound; // 0-infty, controls saturation

public:

  mrank_options(void)
//set default values
    :quiet(DEFAULT_QUIET),   
     verbose(DEFAULT_VERBOSE), 
     precision(DEFAULT_PRECISION),   
     hlimq(DEFAULT_HLIMQ), 
     naux(DEFAULT_NAUX),   
     hlimc(DEFAULT_HLIMC), 
     ptl(DEFAULT_PTL), 
     traceequiv(DEFAULT_TRACEEQUIV),
     output_pari(DEFAULT_OUTPUT_PARI), 
     selmer_only(DEFAULT_SELMER), 
     second_descent(DEFAULT_2ND_DESCENT),
     saturation_bound(DEFAULT_SAT_BOUND)
  { ; }

  mrank_options(int q, int v, long p, long hq, long nx, long hc, 
	       int pl, int teq, int o, int sel, int d2, long sat)
    :quiet(q), verbose(v), precision(p), hlimq(hq), naux(nx), hlimc(hc), 
     ptl(pl), traceequiv(teq),  output_pari(o), selmer_only(sel), 
     second_descent(d2), saturation_bound(sat) {;}

  void set(GetOpt& getopt)
    {
      int option_char;
      while ((option_char = getopt ()) != EOF)
	switch (option_char)
	  {
	  case 'h': help();  abort(); break;
	  case 'q': quiet = 1; break;
	  case 'p': precision = atoi (getopt.optarg); break;
	  case 'v': verbose = atoi (getopt.optarg); break;
	  case 'b': hlimq = atoi (getopt.optarg); break;
	  case 'x': naux = atoi (getopt.optarg); break;
	  case 'l': ptl = 1; break;
	  case 't': traceequiv = 1; break;
	  case 'o': output_pari = 1; break;
	  case 's': selmer_only = 1; break;
	  case 'd': second_descent = 0;  break;
	  case 'c': hlimc = atoi (getopt.optarg);
	    if(hlimc>15)
	      {
		cout << "NB: reducing hlimc to 15\n";
	      }
	    break;
	  case 'S': saturation_bound = atoi (getopt.optarg); break;
	  case '?': cerr<< "usage: mwrank"<<
	    " [q p<precision> v<verbosity> b<hlim_q> x<naux>  c<hlim_c> l t o s d>]\n";
	  }
      if(ptl==-99) ptl=(verbose>0);
      if(naux<1) naux=1;  // syzygy sieving MUST have p=3 in it.
    }

  void set(int argc, char **argv)
    { 
      GetOpt getopt (argc, argv, "hqp:v:b:x:ltosdc:S:");
      set(getopt);
    }

  void help(void)
    {	    
      cerr << "mrank/mwrank command line options (can be in any order):\n\n";
      cerr << "-h\t""help""\t\tprints this info and quits\n";
      cerr << "-q\t""quiet""\t\tturns OFF banner display\n";
      cerr << "-v n\t""verbosity""\tsets verbosity to n (default=1)\n";
      cerr << "-o\t""PARI/GP output""\tturns ON extra PARI/GP short output (default is OFF)\n";
      cerr << "-p n\t""precision""\tsets precision to n decimals (default=15)\n";
      cerr << "\t\t\t(irrelevant unless compiled with multiprecision option)\n";
      cerr << "-b n\t""quartic bound""\tbound on quartic point search (default=10)\n";
      cerr << "-x n\t""n aux""\t\tnumber of aux primes used for sieving (default=6)\n";
      cerr << "-l\t""list""\t\tturns ON listing of points (default ON unless v=0)\n";
      cerr << "-t\t""trace""\t\tturns ON trace of quartic equivalence testing (debugging only)\n";
      //      cerr << "-c n\t""curve bound""\tbound on curve point search\n";
      //      cerr << "\t\t\t(default="<<DEFAULT_HLIMC<<", 0 for no point search)\n";
      //      cerr << "\t\t\t(use -1 for automatic)\n";
      cerr << "-s\t""selmer_only""\tif set, computes Selmer rank only (default: not set)\n";
      cerr << "-d\t""skip_2nd_descent""\tif set, skips the second descent for curves with 2-torsion (default: not set)\n";
      cerr << "-S n\t""sat_bd""\t\tupper bound on saturation primes (default="<<DEFAULT_SAT_BOUND<<", -1 for automatic)\n";
      cerr << endl;
    }
  
  void banner(int which) // 0 for mrank, 1 for mwrank
    {	    
      if(which) {cerr << "Program mwrank: ";}
      else      {cerr << "Program mrank: ";}

      cerr << "uses 2-descent (via 2-isogeny if possible) to\n"; 
      cerr << "determine the rank of an elliptic curve E over Q, and list a\n"; 
      cerr << "set of points which generate E(Q) modulo 2E(Q).\n";		    
      if(which)
	cerr << "and finally saturate to obtain generating points on the curve.\n";

      cerr << "For more details see the file mwrank.doc.\n";
      cerr << "For details of algorithms see the author's book.\n\n";
      cerr << "Please acknowledge use of this program in published work, \n";
      cerr << "and send problems to john.cremona@gmail.com.\n\n";
    }

  void show()
    {
      cerr << "List of current options:\n";
      cerr << "Quiet mode "; 
      if(quiet)cerr<<"ON"; else cerr<<"OFF"; cerr<<"\n";
      cerr << "PARI/GP output "; 
      if(output_pari)cerr<<"ON"; else cerr<<"OFF"; cerr<<"\n";
      cerr << "Precision = " << precision << " decimal places (only relevant for LiDIA version)\n";
      cerr << "Verbosity level = " << verbose << "\n";
      cerr << "Limit on height for point search on quartics: "<<hlimq<<"\n";
      cerr << "Number of auxiliary primes for syzygy sieving: "<<naux<<"\n";
      cerr << "Tracing of quartic equivalence testing is ";
      if(traceequiv)cerr<<"ON"; else cerr<<"OFF"; cerr<<"\n";
      cerr << "Point listing is ";
      if(ptl)cerr<<"ON"; else cerr<<"OFF"; cerr<<"\n";
      cerr << "Selmer-rank-only flag is ";
      if(selmer_only)cerr<<"ON"; else cerr<<"OFF"; cerr<<"\n";
      cerr << "do-second-descent flag is ";
      if(second_descent)cerr<<"ON"; else cerr<<"OFF"; cerr<<"\n";
      //      cerr << "Limit on height for point search on curve:    "<<hlimq<<" (only relevant for mwrank)\n";
      cerr << "Saturation bound = "<<saturation_bound<<"\n";
      cerr<<endl;
    }

  int get_quiet() {return quiet;}
  int get_verbose() {return verbose;}
  int get_output_pari() {return output_pari;}
  long get_precision() {return precision;}
  long get_hlimq() {return hlimq;}
  long get_naux() {return naux;}
  long get_hlimc() {return hlimc;}
  int get_ptl() {return ptl;}
  int get_traceequiv() {return traceequiv;}
  int get_selmer_only() {return selmer_only;}
  int get_second_descent() {return second_descent;}
  long get_saturation_bound() {return saturation_bound;}

  void set_quiet(int q) {quiet=q;}
  void set_verbose(int v) {verbose=v;}
  void set_output_pari(int o) {output_pari=o;}
  void set_precision(long p) {precision=p;}
  void set_hlimq(long h) {hlimq=h;}
  void set_naux(long n) {if(n>0) naux=n;  else naux=1;}
  void set_hlimc(long h) {hlimc=h;}
  void set_ptl(int p) {ptl=p;}
  void set_traceequiv(int t) {traceequiv=t;}
  void set_selmer_only(int s) {selmer_only=s;}
  void set_second_descent(int d) {second_descent=d;}
  void set_saturation_bound(long sat) {saturation_bound=sat;}
};
