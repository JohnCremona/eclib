// BLANK.CC:  systematic mwrank
//
#include "points.h"  // from qcurves library
#include "mwprocs.h" //  ""     "       "
#include "getcurve.h"
#include "mquartic.h"
#include "mrank1.h"
#include "mrank2.h"
#include "htconst.h"
#include "version.h"
#include "timer.h"
#include "options.h"

#define MAX_HEIGHT 20 // will never search for points on curve of naive 
                      // logarithmic height greater than this 
#define SEARCH_METHOD 0  // 0 for Stoll sieve (fastest)
                         // 1, 2, 3 for JC's sieve
void process(Curvedata& CD,
	     Curvedata& CD_orig,
	     const bigint& u, const bigint& r, 
	     const bigint& s, const bigint& t,
	     int change, mrank_options opt);

void process_ai(long a1, long a2, long a3, long a4, long a6, 
		mrank_options opt);
void process_a4a6(long a4, long a6, 
		  mrank_options opt);
void process_a4a6max(long a4a6max, 
		     mrank_options opt);
void process_a4a6minmax(long a4a6min, long a4a6max, 
			mrank_options opt);

int main(int argc, char **argv)
{
  mrank_options opt;
  opt.set(argc, argv);
//  opt.show();
  if(!opt.get_quiet()) {opt.banner(1);  show_version();}
#ifdef LiDIA
  long lidia_precision=opt.get_precision();
  bigfloat::precision(lidia_precision);
  cerr<<"Using LiDIA multiprecision floating point with "<<lidia_precision
    <<" decimal places.\n";
#else
  cout.precision(15);
  cerr.precision(15);
#endif
  initprimes("PRIMES",0);
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  if(opt.get_selmer_only()) cout<<"Only computing Selmer Groups!\n";

  long i, rank, srank;
  int change, certain, success, fullmw;
  
  Curvedata CD, CD_orig;
  bigint u,r,s,t;             // transform CD_orig -> CD

  long a4a6min, a4a6max;
  long a, b, c, a1, a2, a3, a4, a6;
  cout<<"Enter first and last max(|a_4|,|a_6|): "; 
  cin>>a4a6min>>a4a6max;
  process_a4a6minmax(a4a6min, a4a6max, opt);
}

void process_a4a6max(long a4a6max, mrank_options opt)
{
  process_a4a6minmax(0,a4a6max,opt);
}

void process_a4a6minmax(long a4a6min, long a4a6max, mrank_options opt)
{
  if(a4a6min==0) 
    {
      process_a4a6(0,0, opt);
      a4a6min=1;
    }
  long a, b;
  for(a=a4a6min; a<=a4a6max; a++)
    {
      //b=0
      process_a4a6( a,0, opt);
      process_a4a6(-a,0, opt);
      process_a4a6(0, a, opt);
      process_a4a6(0,-a, opt);

      //0<b<a
      for(b=1; b<a; b++)
	{
	  process_a4a6( a, b, opt);
	  process_a4a6( a,-b, opt);
	  process_a4a6(-a, b, opt);
	  process_a4a6(-a,-b, opt);
	  process_a4a6( b, a, opt);
	  process_a4a6( b,-a, opt);
	  process_a4a6(-b, a, opt);
	  process_a4a6(-b,-a, opt);
 	}
      //b=a
      process_a4a6( a, a, opt);
      process_a4a6(-a, a, opt);
      process_a4a6( a,-a, opt);
      process_a4a6(-a,-a, opt);
    }
}

void process_a4a6(long a4, long a6, mrank_options opt)
{
  //  cout<<"("<<a4<<","<<a6<<")\n";
  long a1, a2, a3;
  for(a1=0; a1<=1; a1++)
    for(a2=-1; a2<=1; a2++)
      for(a3=0; a3<=1; a3++)
	process_ai(a1,a2,a3,a4,a6, opt);
}

void process_ai(long a1, long a2, long a3, long a4, long a6, mrank_options opt)
{
  bigint A1, A2, A3, A4, A6; 
  A1=a1, A2=a2, A3=a3, A4=a4, A6=a6; 
  Curvedata CD_orig(Curve(A1,A2,A3,A4,A6));
  if(CD_orig.isnull()) return;
  bigint u, r, s, t;
  Curvedata CD=CD_orig.minimalize(u,r,s,t);
  if(CD.isnull()) return;
  int change=(((Curve)CD)!=((Curve)CD_orig));
  if(change) return;
  cout //<< "Curve "
       << (Curve)CD
       //<<" :" 
       <<"\t";
  long ntor = CD.get_ntorsion();
  CurveRed CR(CD);
  bigint N = getconductor(CR);
  cout << N << "\t" << ntor << "\t";

  process(CD, CD_orig, u, r, s, t, change, opt);
}

void process(Curvedata& CD,
	     Curvedata& CD_orig,
	     const bigint& u, const bigint& r, 
	     const bigint& s, const bigint& t,
	     int change, mrank_options opt)
{
  int do_infinite_descent = (opt.get_hlimc()!=0);
  int verbose = (opt.get_verbose());
  int selmer_only = (opt.get_selmer_only());
  int second_descent = (opt.get_second_descent());
  int i, certain, fullmw;
  long rank=0, srank=0;
  vector<Point> plist, shortlist;
  bigfloat ht_c, ht, maxht, hlim;

  rank2 r2(&CD, verbose, 20, opt.get_hlimq(), second_descent, selmer_only);
 
  if (r2.ok())
    {   
      srank =  r2.getrank_ub();
      if(!selmer_only) rank = r2.getrank();
      cout << rank << "\t" << srank << endl;
    }
  else
    {
      rank1 r1(&CD, verbose, opt.get_traceequiv(), selmer_only, 20,
	       opt.get_hlimq(), opt.get_naux());
      
      if (r1.ok())
	{    
	  srank = r1.getselmer();
	  if(!selmer_only) rank = r1.getrank();
	  cout << rank << "\t" << srank << endl;
	}
      else cout << "Failed to compute rank\n";
    }

  /*
//cout<<"plist = "<<plist<<endl;
//cout<<"shortlist = "<<shortlist<<endl;

  mw mwbasis(&CD, (verbose>1));
  mwbasis.process(shortlist);
  rank = mwbasis.getrank();
  fullmw=(rank==0);
  
  if(do_infinite_descent&&(rank>0))
    {
      fullmw=1;  // but may be reset to 0 below...
      // Compute the height up to which to search
      ht_c = height_constant(CD);
      if(verbose) 
	cout << "Height Constant = " << ht_c << "\n";
      maxht=0;
      for (i=0; i<plist.size(); i++)
	{ 
	  ht  = height(plist[i]);
	  if(ht>maxht) maxht=ht;
	}
      if(verbose)
	cout<<"\nMax height = "<<maxht<<endl;
      if(rank==1) maxht/=9;  // in this case the index is at least 3
      hlim = maxht+ht_c;
      if(verbose)
	cout<<"Bound on naive height of extra generators = "<<hlim<<endl;
      long hlimc = opt.get_hlimc();
      if((hlimc>=0) && (hlimc<hlim))
	{
	  fullmw=0;
	  hlim = hlimc;
	  if(verbose)
	    cout << "Only searching up to height " << hlim << endl;
	}
      else
	if((hlimc<0) && (MAX_HEIGHT<hlim))
	  {
	    fullmw=0;
	    hlim = MAX_HEIGHT;
	    if(verbose)
	      cout << "Only searching up to height "<<MAX_HEIGHT<<"\n";
	  }
// Do the extra search
      mwbasis.search(hlim, SEARCH_METHOD,0);
      rank = mwbasis.getrank();
      if(verbose)
	cout<<"After point search, rank of points found is "<<rank<<endl;
    }
  else
    if(verbose)
      cout<<"After descent, rank of points found is "<<rank<<endl;
  
// End  of extra search for points
//
// Final output of points:

  plist = mwbasis.getbasis();
  if(change&&verbose&&(rank>0)) 
    cout<<"Transferring points back to original curve "
	<<(Curve)(CD_orig)<<"\n";
  for (i=0; i<rank; i++)
    { 
      Point P = plist[i];
      if(change) 
	{
	  P = shift(P,&CD_orig,u,r,s,t,1);
	  plist[i] = P;
	}
      if(verbose||opt.get_ptl())
	{
	  cout << "\nGenerator "<<(i+1)<<" is "<<P<<"; ";
	  cout << "height "<<height(P);
	  if(!P.isvalid()) cout<<" --warning: NOT on curve!";
	}
    }

// Final report

  if(verbose||opt.get_ptl())
    {
      cout<<endl<<endl;
      if(certain)
	{
	  if(fullmw)
	    {
	      cout << "The rank and full Mordell-Weil basis have been determined unconditionally.\n";
	      cout << "Regulator = "<<mwbasis.regulator()<<endl<<endl;
	    }
	  else
	    {		
	      cout << "The rank has been determined unconditionally.\n";
	      if(rank>0)
		{
		  cout << "The basis given is for a subgroup of full rank of the Mordell-Weil group\n";
		  cout << " (modulo torsion), possibly of index greater than 1.\n";
		  cout << "Regulator (of this subgroup) = "<<mwbasis.regulator()<<endl;
		}
	      else
		{
		  cout << "Regulator = 1\n";
		}
	      cout<<endl;
	    }
	}
      else // not certain of the rank
	{
	  cout << "The rank has not been completely determined, \n";
	  cout << "only a lower bound of "<<rank<<" and an upper bound of "<<srank<<".\n";
	  long diff=srank-rank;
	  if(odd(diff))
	    {
	      cout<<"Standard parity conjectures would increase the lower bound by 1 to "<<(rank+1);
	      if(diff==1)
		{
		  cout<<",\n implying that the rank was exactly "<<srank;
		}
	      cout<<".";
	    }
	  cout<<endl;
	  if(fullmw)
	    {
	      if(rank>0)
		{
		  cout << "If the lower bound is strict, the basis given ";
		  cout << "is for the full Mordell-Weil group (modulo torsion).\n";
		  cout << "Regulator (of this subgroup) = "<<mwbasis.regulator()<<endl;
		}
	    }
	  else
	    {
	      if(rank>0)
		{
		  cout << "Even if the lower bound is strict, ";
		  cout << "the basis given is for a subgroup of the Mordell-Weil group\n ";
		  cout << " (modulo torsion), possibly of index greater than 1.\n";
		  cout << "Regulator (of this subgroup) = "<<mwbasis.regulator()<<endl;
		}
	      cout<<endl;
	    }
	}
    }
  
//Pari interface output:
  if(opt.get_output_pari())
    {
      cout<<"[["<<rank;
      if(rank<srank)  cout<<","<<srank;
      cout<<"],[";
      for(i=0; i<rank; i++)
	{
	  if(i) cout<<",";
	  output_pari(cout,plist[i]);
	}
      cout<<"]]\n";
    }
  */  
  if(verbose) cout<<endl;
}
