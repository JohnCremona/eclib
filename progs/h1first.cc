// FILE h1first.cc :  h1 (full space)
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

#include <eclib/cperiods.h>     //from qcurves, for computing conductors
#include <eclib/periods.h>
#include <eclib/pcprocs.h>

#ifndef SINGLE   // so Makefile can override
#define AUTOLOOP
#endif
#define SHOWCURVES
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output
                          // otherwise, sorts newforms into Cremona order before output

const scalar modulus(default_modulus<scalar>());

int main(void)
{
  int verbose,output,curve_output,n=130;
  cout << "Program h1first." << endl;
  cerr << "Verbose output? "; cin>>verbose;
  cerr << "Output updated newform data? "; cin>>output;
  cerr << "Output updated curve file? (0/1) ";  cin >> curve_output;
  string curve_filename;
#ifdef AUTOLOOP
  int limit=210;
  cerr<<"Enter first and last N: ";cin>>n>>limit;  n--;
  while (n<limit) { n++;
#else
#ifndef SINGLE
  while (n>1)
#endif
    { cout<<"\n\nEnter level: "; cin>>n;
#endif
  if (n>1)
    {
      if (curve_output)
        curve_filename=single_curve_filename(n);
      else
        curve_filename="no";
      if(verbose) cout << "\n\n";
      cout << ">>>Level " << n << "<<<" << endl;
      newforms nf(n, modulus, verbose>1);
      if(verbose)
        {
          cout << "Reading newform data from file..." << flush;
        }
      nf.createfromdata(1,0,0);
      long nnf = nf.n1ds;
      if(verbose)
        {
          cout << "done: " << nnf << " newforms." << endl;
        }
      if(nnf==0)
        {
          cout << "No newforms.\n";
          cout << "Finished level "<<n<<endl;
#ifdef SINGLE
          exit(0);
#else
          continue;
#endif
        }
#ifdef LMFDB_ORDER
      nf.sort_into_LMFDB_label_order();
#else
      nf.sort_into_Cremona_label_order();
#endif
      int all_nf = 1; // default; means do all
#ifdef SINGLE
      all_nf = 0;
      nf.nf_subset.clear();
      cout << "Enter list of form numbers (between 1 and "<<nnf<<"), ending with 0: ";
      int inf = 1;
      cin>>inf;
      while((inf>0)&&(inf<=nnf))
        {
          nf.nf_subset.push_back(inf-1); // actually count from 0 despite UI using 1 as base
          cin>>inf;
        }
      nnf = nf.nf_subset.size();
      cout << endl;
#endif
      if(verbose)
        {
          cout << "Working on " << nnf << " newforms: " << nf.nf_subset << " ..." << endl;
          cout << "Finding +1 eigenvectors..." << flush;
        }
      nf.makebases(1, all_nf);
      if(verbose)
        {
          cout << "done.\nNow finding -1 eigenvectors..." << flush;
        }
      nf.set_sign(-1);
      nf.makebases(1, all_nf);
      if(verbose)
        {
          cout << "done.\nNow filling in data for newforms..."<<flush;
          cout << "about to call merge" << endl;
        }
      nf.merge(all_nf);

      if(verbose)
        {
          cout << "done.\nUpdated newform(s): ";
          if (all_nf)
            {
              nf.display();
            }
          else
            {
              cout<<"updated newforms "<<nf.nf_subset<<":"<<endl;
              for( const auto& nfi : nf.nf_subset)
                {
                  cout<<"# "<<nfi<<":\t"<<endl;
                  nf.nflist[nfi].display();
                }
            }
        }
      if(output)
        {
          nf.output_to_file();
          if(verbose)
            cout << "saved updated newform data to file" << endl;
        }

#ifdef SHOWCURVES
      // Now we compute the curves (all, unless one_only)
      if(verbose) cout<<"Computing "<<nnf<<" curves...\n";
      vector<int> forms;
      if (all_nf)
        {
          for(int i=0; i<nf.n1ds; i++) forms.push_back(i);
        }
      else
        {
          for(int i=0; i<nnf; i++) forms.push_back(nf.nf_subset[i]);
        }
      vector<int> failures = nf.showcurves(forms,0,curve_filename);
      if(failures.size()>0)
        {
          cout<<"No curve found for "<<failures.size()<<" forms: "<<failures<<endl;
        }
      else
        {
          if(verbose) cout<<"All curves found OK"<<endl;
          if (!all_nf && curve_output)
            {
              if(verbose)
                cout<<"Recomputing all curves and outputting to "<<curve_filename<<endl;
              forms.clear();
              for(int i=0; i<nf.n1ds; i++) forms.push_back(i);
              nf.showcurves(forms,verbose,curve_filename);
            }
        }

#endif
    }       // end of if(n)
    }       // end of while()
  if(verbose) cout<<"Finished"<<endl;
  }       // end of main()

