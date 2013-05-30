// FILE CHECKGENS.CC -- Program to check input gens are Mordell-Weil basis
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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


#include <fstream>
#include <eclib/points.h>
#include <eclib/mwprocs.h>
#include <eclib/sifter.h>
#include <eclib/compproc.h>
#include <eclib/egr.h>
#include <eclib/htconst.h>
#include <eclib/elog.h>

#define INPUT_CLASS_IS_LETTER // we only use letters now!

Curve C;
Curvedata CD;
vector<Point> plist;
int verbose;
string genfile;
string ccode;

#include <eclib/curvesort.h> // for codeletter() function

int main()
{  
//  cout<<"precision = "<<decimal_precision()<<endl;
  set_precision(15);
//  cout<<"precision = "<<decimal_precision()<<endl;
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)

  int rank, rank2, i;
  long cond, ncurve, nclass;
  bigfloat reg, hmax;;
  cerr<<"verbose (0/1)? ";             cin >>verbose;
  verbose=0;
//  cerr<<"\nLimit on height of search (-1 to just check points are on curve)? ";
//  cin >> hmax;
  hmax=-1;
  cerr <<"\n";
  cerr << "input filename for curves and generators? "; 
  cin >> genfile;
  ifstream genin;
  genin.open(genfile.c_str());
  if(!genin.is_open()) {cerr<<"Unable to open file " << genfile << "\n"; abort();}
  cerr<<endl;

  while (!genin.eof())
    {
      genin >> cond;
      //      cout<<"input conductor="<<cond<<endl;
      if(cond==0) break;
#ifdef INPUT_CLASS_IS_LETTER
      genin >> ccode;
#else
      genin >> nclass;
      ccode = codeletter(nclass-1);
      //      cout<<"After input, nclass="<<nclass<<endl;
#endif
      genin >> ncurve;
      //      cout<<"After input, ncurve="<<ncurve<<endl;
      genin>>C;
      //      cout<<"After input, C="<<C<<endl;
      cout << cond<<ccode<<ncurve<<"\t";
      cout << C << ":\t";
      CD = Curvedata(C);
      if(verbose) cout <<endl;
      genin >> rank;
      cout << "r = "<<rank<<"\t";
      if(verbose) cout <<endl;
      if (rank==0) {cout<<" OK\n"; continue;}
      Point P(CD);
      plist.clear();
      while(plist.size()<rank)
	{
	  genin >> P; 
	  cout << P;
	  if(P.isvalid()) 
	    {
	      plist.push_back(P);
	      cout << " OK ";
	    }
	  else 
	    {
	      cout << " not on curve! ";  abort();
	    }
	}

      // A test of our division-by-2 code:

      for(vector<Point>::const_iterator Pi=plist.begin(); Pi!=plist.end(); Pi++)
        {
          Point P = *Pi;
          Point Q = 2*P;
          vector<Point> Qhalves = division_points(CD,Q,2);
          if(Qhalves.size()==0)
            {
              cout<<"ERROR: Q="<<Q<<"=2*P for P="<<P<<" but halving Q failed!"<<endl;
            }
        }


      // Do nothing at all apart from checking points are on curve
      // (since we have a separate program for saturation now):

      if(hmax<0) {cout<<endl;continue;}

      bigfloat ht, maxht = to_bigfloat(0);
      for (i=0; i<rank; i++)
	{ 
	  ht  = height(plist[i]);
	  if(ht>maxht) maxht=ht;
	}

      mw mwbasis(&CD,verbose,1);
      mwbasis.process(plist);
      rank2 = mwbasis.getrank();
      int rank_under = (rank2<rank);
      if(rank_under)
	{
	  cout<<"\nInput points only have rank "<<rank2<<", not "<<rank<<"\n";
//	  continue;
	}
      else cout<<" INDEP ";
      // Now we check that the points are independent in E(Q)/2E(Q):

      long naux=rank;
      long r2=0;
      while((naux<100)&&(r2<rank))
	{
	  naux+=5;
	  sifter box(&CD, naux, verbose);
	  box.process(plist);
	  r2 = box.getrank();
	}
      if(r2==rank)
	{
	  cout<<" INDEP_MOD_2 ";
	  if(verbose)
	    cout<<"(Points are independent mod 2)"<<endl;
	}
      else
	cout<<"\nPoints may be dependent mod 2, r2 = "<<r2<<endl;



      if(hmax<1) {cout<<endl;}
      else
	{
      reg = mwbasis.regulator();
	
      if(verbose) cout<<"\nMax height = "<<maxht<<endl;
      if(rank==1) 
        {
          if(r2==rank) maxht/=9;	else maxht/=4;
        }
      // since in this case we know that the index is at least 2 (or 3).
      double ht_bound = height_constant(CD);
      if(verbose) cout << "height bound = " << ht_bound << "\n";
      bigfloat hlim = maxht+ht_bound;
      if(verbose||(hlim>hmax)) 
	{
	  if(!verbose) cout<<"\n";
	  cout<<"Bound on naive height of extra generators = "<<hlim<<endl;
	}      
      if(hlim>hmax)
	{
	  cout<<"Only searching up to height "<<hmax<<endl;
	  hlim=hmax;
	}
      else if(!verbose) cout<<endl;
      
      sieve s(&CD, &mwbasis, 2, 0);
      s.search(hlim);

      rank2 = mwbasis.getrank();
      if(rank2!=rank)
	{
	  cout << "\nPoints found have rank " << rank2<<", not "<<rank<<"\n";
	  continue;
	}
      if(rank_under)
	cout << "New basis: "<<mwbasis.getbasis()<<endl;

      bigfloat reg2 = mwbasis.regulator();

      long index = I2long(Iround((sqrt(reg/reg2)+0.01)));

      if(index>1)
	{
	  if(!verbose) cout<<"\n";
	  cout << "Points found contain original with index " << index <<"\n";

	  vector<Point> b = mwbasis.getbasis();
	  for (i=0; i<rank; i++)
	    { Point P = b[i];
	      cout << "\nGenerator "<<(i+1)<<" is "<<P<<"; ";
	      cout << "height "<<height(P);
	    }
	  cout<<"\nNew regulator is "<<reg2<<" (old was "<<reg<<")"<<endl<<endl;
	}
	}
//      else cout << "\n";
    }
}         // end main()


