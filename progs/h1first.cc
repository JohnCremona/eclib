// FILE h1first.cc :  h1 (full space)
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

#include <eclib/marith.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/oldforms.h>
#include <eclib/cperiods.h>     //from qcurves, for computing conductors
#include <eclib/newforms.h>
#include <eclib/periods.h>
#include <eclib/pcprocs.h>

#define AUTOLOOP
#define SHOWCURVES
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output

int main(void)
{
  long prec0=25;
  int verbose,output,limit=210,n=130;
  cout << "Program h1first.  Using METHOD = " << METHOD << endl;
  cerr << "Verbose output? "; cin>>verbose;
  cerr << "Output updated newform data? "; cin>>output;
#ifdef AUTOLOOP
  cerr<<"Enter first and last N: ";cin>>n>>limit;  n--;
  while (n<limit) { n++;
#else
  while (n>0) { cout<<"\n\nEnter level: "; cin>>n;
#endif
  if (n>0)
    {
      if(verbose) cout << "\n\n";
      cout << ">>>Level " << n << "<<<" << endl;
      newforms nf(n,verbose>1);
      if(verbose)
        {
          cout << "Reading newform data from file..." << flush;
        }
      nf.createfromdata(1,0,0);
      long inf, nnf = nf.n1ds;
      if(verbose)
        {
          cout << "done: " << nnf << " newforms." << endl;
        }
      if(nnf==0)
        {
          cout << "No newforms.\n";
          cout << "Finished level "<<n<<endl;
          continue;
        }
#ifdef LMFDB_ORDER
      nf.sort();
#endif
      if(verbose)
        {
          cout << "Finding +1 eigenvectors..." << flush;
        }
      nf.makebases(1);
      if(verbose)
        {
          cout << "done.\nNow finding -1 eigenvectors..." << flush;
        }
      nf.set_sign(-1);
      nf.makebases(1);
      if(verbose)
        {
          cout << "done.\nNow filling in data for newforms..."<<flush;
        }
      nf.set_sign(0);
      nf.merge();

      if(verbose)
        {
          cout << "done.\nUpdated newforms: ";
          nf.display();
        }
      if(output)
        {
          nf.output_to_file();
          if(verbose)
            cout << "saved updated newform data to file" << endl;
        }

#ifdef SHOWCURVES
      // Now we compute the curves
      if(verbose) cout<<"Computing "<<nnf<<" curves...\n";
      vector<int> forms;
      for(inf=0; inf<nnf; inf++) forms.push_back(inf);
      vector<int> failures = nf.showcurves(forms,0);
      if(failures.size()>0)
	{
	  cout<<"No curve found for "<<failures.size()<<" forms: "<<failures<<endl;
	}
      else
	{
	  if(verbose) cout<<"All curves found OK"<<endl;
	}
#endif
    }       // end of if(n)
  }       // end of while()
  if(verbose) cout<<"Finished"<<endl;
  }       // end of main()

