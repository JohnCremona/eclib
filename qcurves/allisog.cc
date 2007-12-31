// allisog.cc: program to find curves isogenous to input curves
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
 
//#define TEST
//#define DUMP_DATA

#include "matrix.h"
#include "isogs.h"
#include "getcurve.h"

int main(){
  set_precision(string("Enter number of decimal places").c_str());
  int verbose=0;
  cout << "Verbose? (0/1) " << flush;  cin >> verbose;

  initprimes(string("PRIMES").c_str(),0);
	
  bigint v;
  vector<bigrational> ai(5);

  while (getcurve(ai,verbose))
    {
      Curvedata D(Curvedata(ai,v),1);
      CurveRed C(D);
      IsogenyClass cl(C ,verbose);
      cl.grow();
      //      cl.display(cout);
      vector<CurveRed> crs=cl.getcurves();
      vector<Curve> cs;
      for(unsigned int i=0; i<crs.size(); i++) cs.push_back((Curve)(crs[i]));
      cout<<cs<<"\t";
      cl.getmatrix().output_pari(cout); 
      cout<<endl;
#ifdef TEST
      vector<CurveRed> curves = cl.getcurves();
      cout << "Curve list via getcurves(): " << curves << endl;
      vector<long> m = cl.getmat();
      long ncurves = curves.size();
      long i,j;
      cout << "Isogeny matrix via getmat():\n";
      cout << "\t"; for(j=0; j<ncurves; j++) cout<<(j+1)<<"\t";  cout<<"\n";
      for(i=0; i<ncurves; i++)
	{
	  cout<<(i+1)<<"\t"; 
	  for(j=0; j<ncurves; j++) 
	    cout<<m[i*ncurves+j]<<"\t";  
	  cout<<"\n";
	}
      cout << "Isogeny matrix via mat():\n";
      cout << "\t"; for(j=0; j<ncurves; j++) cout<<(j+1)<<"\t";  cout<<"\n";
      for(i=0; i<ncurves; i++)
	{
	  cout<<(i+1)<<"\t"; 
	  for(j=0; j<ncurves; j++) 
	    cout<<cl.mat(i,j)<<"\t";  
	  cout<<"\n";
	}
      cout<<endl;
#endif
#ifdef DUMP_DATA
      cout << "Output from dumpdata (rank set arbitrarily to 99):\n\n";
      cl.dumpdata(cout,99);
#endif
    }
}
