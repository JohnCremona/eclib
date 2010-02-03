// descent.cc: implementation of classes rank12 and two_descent
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
 

// rank12 is a common base for separate classes rank1 and rank2 (for
// computing rank via general 2-descent and descent via 2-isogeny
// repectively); class two_descent is a user interface to these

#include "compproc.h"
#include "points.h"
#include "mwprocs.h"
#include "mquartic.h"
#include "descent.h"
#include "mrank1.h"
#include "mrank2.h"

#define PRE_SATURATION_SEARCH_LIMIT 8

// Constructor:
//
// sel is selmer_only switch
// firstlim is bound on |x|+|z|
// secondlim is bound on log max {|x|,|z| }, i.e. logarithmic
// n_aux only relevant for general 2-descent when 2-torsion trivial
// n_aux=-1 causes default to be used (depends on method)
// second_descent only relevant for descent via 2-isogeny

two_descent::two_descent(Curvedata* ec, 
	      int verb, int sel, 
	      long firstlim, long secondlim, 
	      long n_aux, int second_descent)
  :verbose(verb), selmer_only(sel), e_orig(*ec) 
  {
    qai.resize(5);
    bigint a1,a2,a3,a4,a6;
    ec->getai(a1,a2,a3,a4,a6);
    v=BIGINT(1); 
    qai[0]=a1; qai[1]=a2; qai[2]=a3; qai[3]=a4; qai[4]=a6;
    do_the_descent(firstlim,secondlim,n_aux,second_descent);
  }

//#define DEBUG_Q_INPUT

two_descent::two_descent(vector<bigrational> ai, 
	      int verb, int sel, 
	      long firstlim, long secondlim, 
	      long n_aux, int second_descent)
  :verbose(verb), selmer_only(sel)
  {
    // Construct Curvedata from rational coeffs & get the scaling
    // factor v used 
#ifdef DEBUG_Q_INPUT
    cout<<"two_descent constructor called with "<<ai<<endl;
#endif
    qai=ai;
    e_orig = Curvedata(ai, v);
    if(e_orig.isnull()) exit(1);
    if(verbose&&(v!=1)) 
      cout<<"integral model = "<<(Curve)e_orig<<" with scale factor "<<v<<endl;
    // Do the work as in the previous constructor
    do_the_descent(firstlim,secondlim,n_aux,second_descent);
  }

void two_descent::do_the_descent(long firstlim, long secondlim, 
				 long n_aux, int second_descent)
{
    e_min=e_orig.minimalize(u,r,s,t);
    if(verbose&&(e_min!=e_orig))
      {
	cout<<"Working with minimal curve "<<(Curve)e_min
	    <<" via [u,r,s,t] = ["<<u<<","<<r<<","<<s<<","<<t<<"]\n";
      }
    two_torsion_exists = (two_torsion(e_min).size()>1) ;
    if(two_torsion_exists) 
      r12=new rank2(&e_min,verbose,selmer_only,firstlim,secondlim,second_descent);
    else
      r12=new rank1(&e_min,verbose,selmer_only,firstlim,secondlim,n_aux);

    success=r12->ok();
    rank = r12->getrank();
    rank_bound = r12->getrankbound();
    selmer_rank = r12->getselmer();
    certain=r12->getcertain();
    // the last parameter here is set to a known upper bound on the rank
    mwbasis = new mw(&e_min,verbose>2,1,rank_bound);
}

void two_descent::report_rank() const
{
  if(!success) {cout << "Failed to compute rank\n"; return;}
  if(selmer_only)
    {
      cout << "selmer-rank = " << selmer_rank << endl;
      cout << "upper bound on rank = " << rank_bound << endl;
    }
  else
    {
      if(certain) 
	cout << "Rank = " << rank << endl;
      else 
	{
	  if(two_torsion_exists)
	    cout<< rank << " <= rank <= " << rank_bound << endl;
	  else
	    cout<< rank << " <= rank <= selmer-rank = " << selmer_rank << endl;
	}   
    } 
}

void two_descent::saturate(long sat_bd)
{

// Do a quick search for points on the curve before processing points
  bigfloat hlim=to_bigfloat(PRE_SATURATION_SEARCH_LIMIT);
  bigfloat oldreg=mwbasis->regulator();
  if(verbose) cout <<"Searching for points (bound = "<<hlim<<")..." << flush;
  mwbasis->search(hlim);
  if(verbose) cout<<"done:"<<endl;
  long search_rank=mwbasis->getrank();
  bigfloat newreg=mwbasis->regulator();
  if(verbose) cout<<"  found points of rank "<<search_rank
		  <<"\n  and regulator "<<newreg<<endl;
  
  if(verbose) cout <<"Processing points found during 2-descent..." << flush;
  mwbasis->process(r12->getgens(),0); // no saturation yet
  if(verbose) cout <<"done:"<<endl;
  rank = mwbasis->getrank();
  if(verbose) 
    {
      if(rank>search_rank) cout << "2-descent increases rank to "<<rank<<", ";
      if(rank<search_rank) cout << "2-descent only finds rank lower bound of "<<rank<<", ";
      cout <<"  now regulator = "<<mwbasis->regulator()<<endl;
    }
  sat_bound=sat_bd; // store for reporting later
  if(sat_bd==0) 
    {
      fullmw=0;
      if(verbose) cout <<"No saturation being done" << endl;
    }
  else
    {
//  Saturate
      if(verbose) cout <<"Saturating (bound = "<<sat_bd<<")..." << flush;
      bigint index; vector<long> unsat;
      int sat_ok = mwbasis->saturate(index,unsat,sat_bd,1);
      // The last parameter 1 says not to bother with 2-saturation!
      if(verbose) cout <<"done:"<<endl;

// Report outcome
      if(verbose)
	{
	  if(index>1) 
	    {
	      cout <<"  *** saturation increased group by index "<<index<<endl;
	      cout <<"  *** regulator is now "<<mwbasis->regulator()<<endl;
	    }
	  else cout << "  points were already saturated."<<endl;
	}
      if(!sat_ok) 
	{
	  cout << "*** saturation possibly incomplete at primes " << unsat << "\n";
	}
      rank = mwbasis->getrank();
      fullmw=sat_ok; // (rank==0);
    }
}

vector<P2Point> two_descent::getbasis() // returns points on original model
{
  vector<Point>plist=mwbasis->getbasis();
  vector<P2Point>qlist(rank);
  for (int i=0; i<rank; i++)
    qlist[i] = scale(transform(plist[i],&e_orig,u,r,s,t,1),v,0);
  return qlist;
}

vector<Point> two_descent::getpbasis() // returns points on integral model
{
  vector<Point>plist=mwbasis->getbasis();
  for (int i=0; i<rank; i++)
    plist[i] = transform(plist[i],&e_orig,u,r,s,t,1);
  return plist;
}

void two_descent::show_gens()  // display points on original model
{
  if(verbose&&(rank>0)) 
    cout<<"Transferring points from minimal curve "<<(Curve)e_min 
	<<" back to original curve " 
	<<"["<<qai[0]<<","<<qai[1]<<","<<qai[2]<<","<<qai[3]<<","<<qai[4]<<"]"
	<<endl;
  cout<<endl;
  vector<Point>plist=mwbasis->getbasis();
  for (int i=0; i<rank; i++)
    { 
      Point P = plist[i];
      bigfloat h=height(P); // must compute height on minimal model!
      P = transform(P,&e_orig,u,r,s,t,1);
      cout << "Generator "<<(i+1)<<" is "<<scale(P,v,0)<<"; " 
	   << "height "<<h;
      if(!P.isvalid()) cout<<" --warning: NOT on curve!";
      cout<<endl;
    }
  cout<<endl;
  cout << "Regulator = "<<mwbasis->regulator()<<endl<<endl;
}

void two_descent::show_result_status()
{
  if(certain)
    {
      if(fullmw)
	{
	  cout << "The rank and full Mordell-Weil basis have been determined unconditionally.\n";
	}
      else
	{		
	  cout << "The rank has been determined unconditionally.\n";
	  if(rank>0)
	    {
	      cout << "The basis given is for a subgroup of full rank of the Mordell-Weil group\n";
	      cout << " (modulo torsion), possibly of index greater than 1\n";
	      if(sat_bound>0)
		cout << " (but not divisible by any prime less than "
		     <<sat_bound<<" unless listed above).\n";
	    }
	  cout<<endl;
	}
    }
  else // not certain of the rank
    {
      cout << "The rank has not been completely determined, \n";
      cout << "only a lower bound of "<<rank
	   <<" and an upper bound of "<<rank_bound<<".\n";
      cout<<endl;
      if(fullmw)
	{
	  if(rank>0)
	    {
	      cout << "If the rank is equal to the lower bound, the basis given ";
	      cout << "is for the full Mordell-Weil group (modulo torsion).\n";
	    }
	}
      else
	{
	  if(rank>0)
	    {
	      cout << "Even if the lower bound is strict, ";
	      cout << "the basis given is for a subgroup of the Mordell-Weil group\n ";
	      cout << " (modulo torsion), possibly of index greater than 1.\n";
	    }
	  cout<<endl;
	}
    }
}

void two_descent::pari_output()
{
  vector<P2Point>plist=getbasis();
  cout<<"[["<<rank;
  if(rank<rank_bound)  cout<<","<<rank_bound;
  cout<<"],[";
  for(int i=0; i<rank; i++)
    {
      if(i) cout<<",";
      output_pari(cout,plist[i]);
    }
  cout<<"]]\n";
}

