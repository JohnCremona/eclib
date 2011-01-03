// FILE moreap.cc: computes more ap for given level(s)
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2010 John Cremona
// 
// This file is part of the mwrank/g0n package.
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
#include <fstream>
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "cperiods.h"     //from qcurves, for computing conductors
#include "newforms.h"

//#define AUTOLOOP

int main(void)
{
 int n=1; 
 int lastn, stopp, output, showeigs, showforms, findcurves;
 int nnf, nap;
 cout << "Program moreap\n";
 cout << "---------------\n\n";
 cout << "For each N, assumes that the file newforms/xN exists, and computes more\n";
 cout << "Hecke eigenvalues.\n";
 cout << "Output new eigs to file (1/0)? ";  cin>>output;
 cout << "Output new eigs to screen (1/0)? "; cin>>showeigs;
 cout << "Display newforms (1/0)? "; cin>>showforms;
 cout << "Attempt curve construction (1/0)? "; cin>>findcurves;
 if(findcurves) set_precision("Enter number of decimal places");

#ifdef AUTOLOOP
 cout << "How many primes for Hecke eigenvalues? ";
 cin  >> stopp; cout << endl;
 cout<<"Enter first and last N: ";cin>>n>>lastn; n--;
 while (n<lastn) { n++;
#else
 while (cout<<"Enter level: ", cin>>n, n>0) { 
 cout << "How many primes for Hecke eigenvalues? ";
 cin  >> stopp; cout << endl;
#endif
  cout << "\n>>>Level " << n << "\t";
  newforms nf(n,showforms); 
  nf.createfromdata(1,0,0);
  if (showforms) nf.newforms::display();

  nnf = nf.n1ds;
  if(nnf==0)
    {
      cout<<"No newforms."<<endl;
      continue;
    }
  nap = nf.nflist[0].aplist.size();
  if(nap>=stopp)
    {
      cout<<"Already have "<<nap<<" eigenvalues on file, no need to compute more."<<endl;
    }
  else
    {
      nf.makebases(1);
      cout << "About to start computing ap..."<<flush;
      nf.make_projcoord();
      nf.find_jlist();
      nf.addap(stopp);
      cout << "...done."<<endl;
      if(output) nf.output_to_file();
    }
  if(findcurves)
    {
      int inf, nnf = nf.n1ds;
      // Now we compute the curves
      cout<<"Computing "<<nnf<<" curves...\n";
      vector<int> forms;
      for(inf=0; inf<nnf; inf++) forms.push_back(inf);
      vector<int> failures = nf.showcurves(forms,0);
      if(failures.size()>0)
	{
	  cout<<"No curve found for "<<failures.size()<<" forms: "<<failures<<endl;
	}
      else
	{
	  cout<<"All curves found OK"<<endl;
	}
    }
}       // end of while()
}       // end of main()
