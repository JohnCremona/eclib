// point_search.cc:  program to find points up to given naive height
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
 
#include <eclib/interface.h>
#include <eclib/compproc.h>
#include <eclib/points.h>
#include <eclib/polys.h>
#include <eclib/curvemod.h>
#include <eclib/pointsmod.h>
#include <eclib/ffmod.h>
#include <eclib/divpol.h>
#include <eclib/tlss.h>
#include <eclib/elog.h>
#include <eclib/saturate.h>
#include <eclib/getcurve.h>
#include <eclib/sieve_search.h>
#include <eclib/mwprocs.h>

int main()
{
#ifdef MPFP
  set_precision(70);
#endif
  initprimes("PRIMES",0);

  bigfloat ht_limit;
  bigint u,r,s,t;
  int verbose = 1, modopt=0, pp=1;
  long blength, rank, maxrank;
  cerr<<"\nenter search limit: ";      cin>>ht_limit;
  cerr<<"verbose (0/1)? ";             cin >>verbose;
  cerr<<"process points found (1) or just list them (0): "; cin>>pp;
  //  cout<<"moduli option (0 (Stoll)/ 1/2/3)?";      cin >> modopt;
  int verb=1;
  bigint v;
  vector<bigrational> ai(5);

  while (getcurve(ai,verb))
    {
      Curvedata C(ai,v);
      cout << "Input curve ";
      cout <<"["<<ai[0]<<","<<ai[1]<<","<<ai[2]<<","<<ai[3]<<","<<ai[4]<<"]" << endl;
      Curvedata C_min = C.minimalize(u,r,s,t);
      int change_flag = (v!=1) || (((Curve)C_min) != (Curve)C);
      if(change_flag)
	{
	  cout<<"Searching on standard minimal model "<<(Curve)C_min<<endl;
	  cout<<"(points found will be transferred back at end)"<<endl;
	  cout<<"Transformation: \t[u,r,s,t] = ["<<u<<","<<r<<","<<s<<","<<t<<"]";
	  if(v!=1) cout<<" with scale factor "<<v;
	  cout<<endl;
	}
      
      P2Point P0;
      cerr<<"enter number of known points: ";      cin >>blength;
      cerr<<"enter max rank to stop when this is reached (-ve for none): ";
      cin >>maxrank; if(maxrank<0) maxrank=999;
      vector<Point> known_points(blength);
      if (blength)
	{ 
	  for (long j=0; j<blength; j++)
	    { cerr<<"\n  enter point "<<(j+1); 
	      if(change_flag) cerr<<" on original curve";
	      cerr<<": ";
	      cin >> P0;
	      Point P(C,scale(P0,v,1));
	      if ( P.isvalid() ) 
		{
		  Point Q = transform(P,&C_min,u,r,s,t,0);
		  known_points[j] = Q;
		  if(change_flag) 
		    {
		      cout<<P0<<" maps to "<<Q<<" on "<<(Curve)C_min<<", with height "<<height(Q)<<endl;
		    }
		}
	      else 
		{
		  cout<<"Bad point -- not on curve. Start all over again.\n";
		  break; 
		}
	    }
	}

      mw mwbasis(&C_min, verbose, pp, maxrank);
      if(blength)
	{
	  mwbasis.process(known_points);
	  rank = mwbasis.getrank();
	  cout<<"Rank of known points is "<<rank
	      <<" with regulator "<<mwbasis.regulator()<<endl<<endl;
	}

      mwbasis.search(ht_limit, modopt, verbose);
      
      if(pp)
	{
	  long index;
	  vector<long> unsat;
	  if(verbose) 
	    {
	      cout <<"Regulator (before saturation) = "<<mwbasis.regulator()<<"\n";
	      cout <<"Saturating..." << flush;
	    }
	  int sat_ok = mwbasis.saturate(index,unsat);
	  if(verbose)
	    {
	      cout <<"finished saturation (index was " << index << ")\n";
	      cout <<"Regulator (after saturation) = "<<mwbasis.regulator()<<"\n";
	    }
	  if(!sat_ok) 
	    {
	      cout << "saturation possibly incomplete at primes " << unsat << "\n";
	    }
	  rank = mwbasis.getrank();
	  cout<<"Rank of points found is "<<rank<<endl;
	  vector<Point> b = mwbasis.getbasis();
	  vector<P2Point> bb(rank);
	  for (long i=0; i<rank; i++)
	    { 
	      Point Q = b[i];
	      bigfloat hQ = height(Q);
	      cout << "Generator "<<(i+1)<<" is "<<Q<<"; ";
	      cout << "height "<<hQ<<endl;
	      P0 = scale(transform(Q,&C,u,r,s,t,1),v,0);
	      bb[i] = P0;
	      if(change_flag)
		cout<<"\t--maps back to "<<P0<<" on input curve"<<endl;
	    }
	  if(rank>1) cout<<"Regulator =  "<<mwbasis.regulator()<<endl;
	  if(sat_ok)
	    {
	      cout<<"Points have been successfully saturated"<<endl;
	    }
	  else
	    {
	      cout<<"Saturation incomplete "<<endl;
	    }
	  cout<<endl;
// Pari output:
	  cout<<"[";
	  for(long i=0; i<rank; i++)
	    {
	      if(i) cout<<",";
	      output_pari(cout,bb[i]);
	    }
	  cout<<"]\n"<<endl;
	}
    }
  cout<<endl;
}

//end of file point_search.cc
