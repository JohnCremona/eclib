// MRANK.CC:  Mordell-Weil Rank  without infinite descent step
//
#include "points.h"  // from qcurves library
#include "mwprocs.h" //  ""     "       "
#include "getcurve.h"
#include "mquartic.h"
#include "mrank1.h"
#include "mrank2.h"
#include "version.h"
#include "timer.h"
#include "options.h"

int main(int argc, char **argv)
{
  mrank_options opt;
  opt.set(argc, argv);
  if(!opt.get_quiet()) {opt.banner(0);  show_version();}
#ifdef LiDIA
  long lidia_precision=opt.get_precision();
  bigfloat::precision(lidia_precision);
  cerr<<"Using LiDIA multiprecision floating point with "<<lidia_precision
    <<" decimal places.\n";
#else
  cout.precision(15);
  cerr.precision(15);
#endif

  init_time();
  long i, rank, srank;
  int change, certain, success;
  initprimes("PRIMES",0);
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  
  Curvedata CD, CD_orig;
  bigint u,r,s,t;             // transform CD_orig -> CD
  PointArray plist, shortlist;

  while ( getcurve(CD,CD_orig,u,r,s,t,change,!opt.get_quiet())	 )
    {
      start_time();
      rank2 r2(&CD, opt.get_verbose(), 20, opt.get_hlimq());
 
      if (r2.ok())
        {   
	  certain = r2.getcertain();
	  srank = rank = r2.getrank();
	  if(certain) 
	    cout << "Rank = " << rank;
	  else 
	    {
	      srank = r2.getselmer();
	      cout << rank << " <= rank <= selmer-rank = " << srank;
	    }
	  if(opt.get_verbose()||opt.get_ptl()) 
	    cout << endl;
	  else cout<<"\t";

	  shortlist = r2.getgens();
	  if(opt.get_ptl()) 
	    {
	      if(rank<=5)
		{
		  if(change) r2.listpoints(&CD_orig, u, r, s, t);
		  else       r2.listpoints();
		}
	      else
		{
		  if(change) r2.listgens(&CD_orig, u, r, s, t);
		  else       r2.listgens();
		}
	    }
        }
      else
        {
          rank1 r1(&CD, opt.get_verbose(), opt.get_traceequiv(), 20,
		   opt.get_hlimq(), opt.get_naux());
 
          if (r1.ok())
            {    
	      certain = r1.getcertain();
	      srank = rank = r1.getrank();
	      if(certain) 
		cout << "Rank = " << rank;
	      else 
		{
		  srank = r1.getselmer();
		  cout << rank << " <= rank <= selmer-rank = " << srank;
		}

	      if(opt.get_verbose()||opt.get_ptl()) cout << endl;
	      if(opt.get_ptl()&&(rank>0))
		{
		  if(change) r1.listpoints(&CD_orig, u, r, s, t);
		  else       r1.listpoints();
		}
	      plist = r1.getgens();
            }
          else cout << "Failed to compute rank\n";
        }
//Pari interface output:
      if(opt.get_output_pari())
	{
	  if(!opt.get_verbose()&&!opt.get_ptl()) cout<<"\n";
	  mw mwbasis(&CD, opt.get_verbose());
	  mwbasis.process(plist);
	  plist = mwbasis.getbasis();
	  cout<<"[["<<rank;
	  if(rank<srank) cout<<","<<srank;
	  cout<<"],[";
	  for(i=0; i<plist.getlength(); i++)
	    {
	      if(i) cout<<",";
	      cout<<"[";
	      Point P = plist[i];
	      bigrational R(getX(P),getZ(P));
	      cout<<R<<",";
	      bigrational S(getY(P),getZ(P));
	      cout<<S<<"]";
	    }
	  cout<<"]]\n";
	}
      stop_time();
      if(opt.get_verbose>=0) show_time();
    }
  cout<<endl;
}
