// findinf.cc:  program to find points up to given naive height
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
 
#include "interface.h"
#include "compproc.h"

#include "matrix.h"
#include "subspace.h"

#include "points.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "ffmod.h"
#include "divpol.h"
#include "tlss.h"
#include "elog.h"
#include "saturate.h"

#include "sieve_search.h"

#include "mwprocs.h"

int main()
{
  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);

  bigfloat ht_limit;
  bigint u,r,s,t;
  int verbose = 1, modopt=0, pp=1, change_flag; 
  long blength, rank, maxrank;
  cout<<"\nenter search limit: ";      cin>>ht_limit;
  cout<<"verbose (0/1)? ";             cin >>verbose;
  cout<<"process points found (1) or just list them (0): "; cin>>pp;
  //  cout<<"moduli option (0 (Stoll)/ 1/2/3)?";      cin >> modopt;
  Curve E;

  while (1)
    {
      cout<<"\nInput a curve: ";      cin >> E;
      if ( E.isnull() ) break;
      Curvedata C(E);
      cout << "Curve " << E << endl;
      Curvedata C_min = C.minimalize(u,r,s,t);
      change_flag = ((Curve)C_min) != E;
      if(change_flag)
	{
	  cout<<"Searching on standard minimal curve "<<(Curve)C_min<<endl;
	  cout<<"(points found will be transferred back at end)"<<endl;
	  cout<<"Transformation: \t[u,r,s,t] = ["<<u<<","<<r<<","<<s<<","<<t<<"]\n";
	}
      
      Point P(C);
      Point Q(C_min);
      cout<<"enter number of known points: ";      cin >>blength;
      cout<<"enter max rank to stop when this is reached (-ve for none): ";
      cin >>maxrank; if(maxrank<0) maxrank=999;
      vector<Point> known_points(blength);
      if (blength)
	{ 
	  for (long j=0; j<blength; j++)
	    { cout<<"\n  enter point "<<(j+1); 
	    if(change_flag) cout<<" on original curve";
	      cout<<": ";
	      cin >> P;
	      if ( P.isvalid() ) 
		{
		  cout<<"Height on original curve = "<<height(P)<<endl;
		  Q = shift(P,&C_min,u,r,s,t,0);
		  known_points[j] = Q;
		  if(change_flag) 
		    {
		      cout<<P<<" on "<<E<<"\nmaps to "<<
			Q<<" on "<<(Curve)C_min<<", with height "<<height(Q)<<endl;
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
	  cout<<"Rank of known points is "<<rank<<endl;
	  cout<<"\nregulator is "<<mwbasis.regulator()<<endl<<endl;
	}

      mwbasis.search(ht_limit, modopt, verbose);
      
      if(pp)
	{
	  bigint index;
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
	  for (long i=0; i<rank; i++)
	    { 
	      Q = b[i];
	      bigfloat hQ = height(Q);
	      cout << "\tGenerator "<<(i+1)<<" is "<<Q<<"; ";
	      cout << "height "<<hQ<<endl;
	      if(change_flag)
		{
		  P = shift(Q,&C,u,r,s,t,1);
		  cout<<"\t--maps back to "<<P<<" on curve "<<E<<endl;
		  b[i]=P;
		}
	    }
	  if(rank>1) cout<<"\t\tregulator is "<<mwbasis.regulator()<<endl;
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
	      output_pari(cout,b[i]);
	    }
	  cout<<"]\n";
	}
    }
  cout<<endl;
}

//end of file findinf.cc





