// CONGNO.CC:  Congruent number curves
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
//#define AUTOLOOP

int main(int argc, char **argv)
{
  mrank_options opt;
  opt.set(argc, argv);
//  opt.show();
  if(!opt.get_quiet()) {show_version();}
#ifdef LiDIA_ALL
  long lidia_precision=opt.get_precision();
  bigfloat::set_precision(lidia_precision);
  cerr<<"Using LiDIA multiprecision floating point with "<<lidia_precision
    <<" decimal places.\n";
#else
  cout.precision(15);
  cerr.precision(15);
#endif
  initprimes("PRIMES",0);
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)

  init_time();
  if(TIME_CONICS) init_conic_timer();
  long i, rank, srank, rankb;
  int change, certain=0, fullmw;
  int do_infinite_descent = (opt.get_hlimc()!=0);
  int selmer_only = (opt.get_selmer_only());
  int second_descent = (opt.get_second_descent());
  int verbose = (opt.get_verbose());
  bigfloat ht_c, ht, maxht, hlim;
  
  Curvedata CD, CD_orig;
  bigint u,r,s,t;             // transform CD_orig -> CD
  vector<Point> plist, shortlist;

  bigint d, dmin, dmax, zero;
  zero=0;
#ifdef AUTOLOOP
  cout<<"dmin: ";cin>>dmin;
  cout<<"dmax: ";cin>>dmax;
  for(d=dmin; d<=dmax; d+=1)
#else
    while(cout<<"d: ",cin>>d, d>0)
#endif
    {
      CD_orig=Curvedata(Curve(zero,zero,zero,-d*d,zero),0);
      CD=CD_orig.minimalize(u,r,s,t);
      change=(((Curve)CD)!=((Curve)CD_orig));
      cout<<"d="<<d<<": ";
      if(verbose)
	{
	  cout<<endl;
      if(change)
	{
	  cout<<"Original curve "<<(Curve)CD_orig<<"\n";
	  cout<<"Working with minimal curve "<<(Curve)CD<<"\n";
	  cout<<"\t[u,r,s,t] = ["<<u<<","<<r<<","<<s<<","<<t<<"]\n";
	}
      else 
	{
	  cout << "Curve "<< (Curve)CD<<" :\t";
	  change=0;
	}
	}      

      start_time();
      rank=srank=rankb=0;

      rank2 r2(&CD, verbose, 20, opt.get_hlimq(), second_descent, selmer_only);
 
      if (r2.ok())
        {   
	  certain = r2.getcertain();
	  srank = rankb = rank = r2.getrank();
	  if(certain) 
	    cout << "Rank = " << rank;
	  else 
	    {
	      srank = rankb = r2.getrank_ub();
	      cout<<"d="<<d<<": ";
	      cout << rank << " <= rank <= " << rankb;
	    }
	  if(verbose||opt.get_ptl()) cout << endl;
	  else cout<<"\n";

	  shortlist = r2.getgens();
	  if(do_infinite_descent)
	    {
	      r2.makepoints();  // only done automatically for small ranks
	      plist = r2.getpoints();
	    }
        }
      else
        {
          rank1 r1(&CD, verbose, opt.get_traceequiv(), selmer_only, 20,
		   opt.get_hlimq(), opt.get_naux());
 
          if (r1.ok())
            {    
	      certain = r1.getcertain();
	      srank = rankb = rank = r1.getrank();
	      if(certain) 
		cout << "Rank = " << rank;
	      else 
		{
		  srank = r1.getselmer();
		  cout << rank << " <= rank <= selmer-rank = " << srank;
		}
	      cout<<"\n";
	      if(opt.get_ptl()&&(rank>0))
		{
		  cout<<"\n";
		  if(change) r1.listpoints(&CD_orig, u, r, s, t);
		  else       r1.listpoints();
		}
	      shortlist = r1.getgens();
	      if(do_infinite_descent)  plist = r1.getpoints();
	    }
          else cout << "Failed to compute rank\n";
        }

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
	      long diff=rankb-rank;
	      if(odd(diff))
		{
		  cout<<"Standard parity conjectures would increase the lower bound by 1 to "<<(rank+1);
		  if(diff==1)
		    {
		      cout<<",\n implying that the rank was exactly "<<rankb;
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
	  cout<<"\n[["<<rank;
	  if(rank<rankb) cout<<","<<rankb;
	  cout<<"],[";
	  for(i=0; i<rank; i++)
	    {
	      if(i) cout<<",";
	      output_pari(cout,plist[i]);
	    }
	  cout<<"]]\n";
	}

      stop_time();
      if(verbose) 
	show_time();
      if(verbose) cout<<endl;
    }
  cout<<endl;
  if(TIME_CONICS) {show_conic_timer(); cout<<endl;}
}
