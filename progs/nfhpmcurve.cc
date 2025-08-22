// FILE nfhpmcurve.cc main newform- and curve-finding program
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
#include <eclib/timer.h>
#include <eclib/newforms.h>

#define AUTOLOOP
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output
                          // otherwise, sorts newforms into Cremona order before output

#define MAXNAP 20000
#define BITPREC0 100  // initial bit precision
#define BITPRECX  20  // step-size for increasing bit precision
#define BITPRECMAX 300  // maximum bit precision

int main(void)
{
 init_time();
 start_time();
 long n=110, stopp, stopp0; 
 long prec0=BITPREC0;
 int output, verbose;

 cout << "Program nfhpmcurve.  Using METHOD = " << METHOD 
      << " to find newforms" << endl;
 cout << "Verbose output? "; cin>>verbose;
 cout << "How many primes for Hecke eigenvalues? ";
 cin  >> stopp0; cout << endl;
 output=1; 
 cout << "Output newforms to file? (0/1) ";  cin >> output;
#ifdef AUTOLOOP
 int limit;
 cout<<"Enter first and last N: ";cin>>n>>limit; n--;
 while (n<limit) { n++;
#else
     while (n>1) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
  cout << "\n>>>Level " << n;
  // Temporary code to skip non-square-free levels
  //
  //  if(!is_squarefree(n)) {cout<<" --not square-free, skipping!\n"; continue;}
  //
  if(verbose)cout<<endl; else cout<< ":\t"<<flush;
  if(!is_valid_conductor(n))
    {
      cout<<"Not a valid conductor!"<<endl;
      if (output) // output extended full nf data and small nf data
	{
	  output_to_file_no_newforms(n,1,0);
	  output_to_file_no_newforms(n,1,1);
	}
      cout << "Finished level "<<n<<endl;
      continue;
    }
  newforms nf(n,verbose); 
  int noldap=25; // stopp0 must be at least this big!
  if (stopp0<noldap) stopp0=noldap;
  nf.createfromscratch(1,noldap);
#ifdef LMFDB_ORDER
  nf.sort_into_LMFDB_label_order();
#else
  nf.sort_into_Cremona_label_order();
#endif
  nf.make_projcoord(); // needed for when we add more ap
  if(verbose) nf.display();
  else          cout << nf.n1ds << " newform(s) found."<<endl;
  stopp=stopp0;
  nf.addap(stopp);

  long nnf = nf.n1ds, inf; 
  if(nnf==0)
    {
      cout << "No newforms.\n";
      cout << "Finished level "<<n<<endl;
      if(output)  // output extended full nf data and small nf data
	{
	  nf.output_to_file(1,0);
	  nf.output_to_file(1,1);
	}
      continue;
    }

  // Thus far, as in tmanin

  if(verbose)
    {
      cout << endl;
      cout << "Creating minus space..."<<endl;
    }
  nf.set_sign(-1);
  if(verbose)
    {
      cout<<"sign = "<<nf.get_sign()<<endl;
    }
  nf.makebases(1);
  if(verbose)
    {
      cout << "...finished creating minus data."<<endl<<endl;
      nf.display();
      cout << "Filling in data for newforms..."<<endl;
    }
  nf.set_sign(0);
  nf.merge();
  if(verbose) 
    {
      cout << "...finished filling in full data."<<endl;
      cout << "Updated newforms: ";
      nf.display();
    }
  if(output)  // output extended full nf data and small nf data
    {
      nf.output_to_file(1,0);
      nf.output_to_file(1,1);
    }
  // Now we compute the curves
  cout<<"Computing "<<nnf<<" curves...\n";

  long fac=nf.sqfac;
  if(verbose&&(fac>1)) cout<<"c4 factor " << fac << endl;

  stopp = stopp0; // will be increased if necessary
  set_precision(prec0);
  vector<int> forms;
  for(inf=0; inf<nnf; inf++) forms.push_back(inf);

  while(forms.size()>0)
    {
      forms = nf.showcurves(forms,verbose,"no");
      if(forms.size()==0)
	{
	  cout<<"All curves found successfully!"<<endl;
	  cout << "Finished level "<<n<<endl;
	  if(output)  // output extended full nf data and small nf data
	    {
	      nf.output_to_file(1,0);
	      nf.output_to_file(1,1);
	    }
	}
      else
	{
	  cout<<forms.size()<<" curve(s) missing: ";
	  for( int i : forms)
	    cout<<i+1<<" ";
	  cout<<endl;
	  int newstopp;
	  if(stopp<500)
	    newstopp=2*stopp;
	  else
	    newstopp=stopp+500;
	  if(newstopp>MAXNAP)
	    {
	      cout<<"Cannot compute more ap, something must be wrong in newform data"<<endl;
	    }
          else
            {
              cout<<"Computing some more ap: from "<<stopp+1<<" to "
                  <<newstopp<<"..."<<endl;
              nf.set_sign(1);
              nf.makeh1(1);
              nf.addap(newstopp);
              stopp=newstopp;
              if(output)  // output extended full nf data and small nf data
                {
                  nf.output_to_file(1,0);
                  nf.output_to_file(1,1);
                }
#ifdef MPFP
              if(bit_precision()<BITPRECMAX)
                {
                  set_precision(bit_precision()+BITPRECX);
                  cout << "Now working with bit precision "<<bit_precision()<< endl;
                }
#endif
            }
        }
    }
 }       // end of if(n>1)
     }  // end of while(n>1) or while(n<limit)
 }       // end of main()

