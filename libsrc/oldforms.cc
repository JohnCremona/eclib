// FILE OLDFORMS.CC: implementation of class oldforms
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

#include <eclib/cperiods.h>
#include <eclib/oldforms.h>
#include <eclib/newforms.h>

inline int testbit(long a, long i) {return (a& (1<<i));}

bool file_exists(string filename)
{
  return ifstream(filename.c_str()).good();
}

// returns a filename of the form "$ECDB/curves.a0000-a9999" where "a"
// is a single digit 0..9 or two digits 10..99, chosen so that N is in
// the range, where 10<N<10^6.  The environment variable ECDB should
// be set to point to the curves subdirectory of a clone of
// https://github.com/JohnCremona/ecdata; if it is not set then a
// subdirectory ./curves of the current directory is used.
string ecdb_filename(long N)
{
  long a = N/10000; // rounded down
  stringstream s;
  s << getenv_with_default("ECDB","./curves");
  s << "/curves." << a << "0000-" << a << "9999";
  //  cout << "Filename for curves of conductor " << d << " is " << s.str() << endl;
  return s.str();
}

// returns a filename of the form "$TCURVES/curves.N" where N is a
// conductor.  The environment variable TCURVES should be set to point
// to a suitable ubdirectory, for example one where curves are output
// by progs/nfhpcurve; if it is not set then a subdirectory ./tcurves
// of the current directory is used.
string single_curve_filename(long N)
{
  stringstream s;
  s << getenv_with_default("TCURVES_DIR","./tcurves");
  s << "/curves." << N;
  //  cout << "Filename for curves of conductor " << N << " is " << s.str() << endl;
  return s.str();
}

// Combination of the above: return the name of a single conductor's
// filename if it exists, else that of a database curve file
string curve_filename(long N)
{
  string s = single_curve_filename(N);
  if (!file_exists(s))
    s = ecdb_filename(N);
  return s;
}

// Read from a curves file and store the list of curve whose conductor
// is N.  Each line of the file should have the format like
// 11 a 1 [0,-1,1,-10,-20] 0 5 0
// with N first and 1 curve per isogeny class.  Fields 2,5,6,7 are
// ignored.  The curves must be sorted by conductor.

vector<CurveRed> get_curves(string filename, long N)
{
  vector<CurveRed> curves;
  ifstream infile(filename.c_str());
  if (!infile.is_open())
    {
      cerr<<"Unable to open file "<<filename<<" for curve input"<<endl;
      return curves;
    }
  Curve C;
  long iN, dum;
  string code;
  infile >> iN;
  // skip over curves whose conductors are too small
  while ((iN<N)&&!infile.eof())
    {
      infile >> code >> dum >> C >> dum >> dum >> dum;
      // check for end of file:
      infile >> ws;
      if (!infile.eof())
	infile >> iN;
    }
  if (infile.eof()) return curves;
  // read curves whose conductors are just right
  while ((iN==N)&&!infile.eof())
    {
      infile >> code >> dum >> C;
      CurveRed CR(Curvedata(C,0));
      if (getconductor(CR)!=N)
	cerr<<"Wrong conductor "<<getconductor(CR)<<" for "<<C<<": should be "<<N<<endl;
      curves.push_back(CR);
      infile >> dum >> dum >> dum;
      // check for end of file:
      infile >> ws;
      if (!infile.eof()) 
	infile >> iN;
    }
  // cout << "Read "<<curves.size()<<" curve(s) of conductor "<<N<<" from "<<filename<<endl;
  // for (dum=0; dum<curves.size(); dum++) cout<<(Curve)curves[dum]<<endl;
  return curves;
}

// Implementation of oldforms member functions

oldforms::oldforms(long intp, const level* iN, scalar mod, int verbose, int plus)
  : noldclasses(0), nap(intp), ntp(intp), totalolddim(0), N(iN), modulus(mod), plusflag(plus)
{
   for (const auto& M : N->dlist)
     {
       if(M<11) continue;
       if(M==(N->N)) continue;
       getoldclasses(M,verbose);
     }
   if(verbose) cout<<"Finished getting oldclasses "<<endl;
   for (long i=0; i<noldclasses; i++) totalolddim+=oldclassdims[i];
}

void oldforms::getoldclasses(long d, int verbose) 
{
  long n = N->N;
  if ((d>10) && (n>d))
    {
      if(verbose) cout << "Getting oldclasses for divisor M = " << d << "\n";
      newforms olddata(d, modulus, verbose);
      string curve_file = curve_filename(d);
      if (file_exists(curve_file))
	{
	  vector<CurveRed> curves = get_curves(curve_file,d);
	  olddata.createfromcurves_mini(curves);
	}
      else
	{
	  // if a newforms file exists for level d it will read data
	  // from it, otherwise it will create the newforms at level d
	  // (and store them)
	  olddata.createfromdata(1,25,
				 1 /* create_from_scratch_if_absent */,
				 1 /* small_data_ok */);
	}
      long nforms = olddata.n1ds;
      long oldnap = olddata.nap;
      if(nforms==0) return;
      if(verbose>1) cout << "Computing W multiplicities." << "\n";
      long m = n/d;
      long k=0, xmult, mult, j, beta;
      vector<long> betalist; // =new long[N->npdivs];
      for (const auto& qj : N->plist)
	{
	  beta=val(qj,m);
	  if(beta>0) k++;
	  betalist.push_back(beta);
	}
      if(verbose>1) cout<<"betas: "<<betalist<<endl;
      vector<long> nextoldformap(nap);
      primevar pr; long iform, c, ip, aq; int bit;
      for(iform=0; iform<nforms; iform++)
	{
	  vector<long>& aqlist=olddata.nflist[iform].aqlist;
	  nextoldformap = olddata.nflist[iform].aplist;
	  if(verbose>1)
	    {
	      cout<<"form #"<<(iform+1)<<": "<<"aqlist="<<aqlist<<endl;
	      cout<<"aplist before adjusting (size "<<oldnap<<") ="<<nextoldformap<<endl;
	    }
	  for (c=0; c<(1<<k); c++) // 2^k different oldclasses
	    {
	      if(verbose>1) cout<<"c="<<c<<endl;
	      mult=1; j=0;
	      auto betai=betalist.begin();
	      auto aqj=aqlist.begin();
              for ( const auto& p : N->plist)
		{
                  ip = prime_pi(p);
		  if(verbose>1) cout<<"p="<<p<<" (ip="<<ip<<")"<<endl;
                  beta = *betai++;
                  if(::divides(p,d))
                    aq = *aqj++;
                  else
                    aq=1;
                  if (beta==0)
                    {
                      if(ip<oldnap)
                        {
                          nextoldformap[ip] = aq;
                          if(verbose>1)
                            cout<<"setting entry #"<<ip<<" to "<<aq<<endl;
                        }
                    }
                  else
                    {
                      bit = testbit(c,j++);
                      if(ip<oldnap)
                        {
                          nextoldformap[ip] = bit?1:-1;
                          if(verbose>1)
                            cout<<"setting entry #"<<ip<<" to "<<aq<<endl;
                        }
                      if (odd(beta))
                        xmult =  (beta+1)/2;
                      else
                        {
                          xmult=beta/2 +1;
                          if(::divides(p,d) && (aq==-1))
                            xmult--;
                        }
                      if (!bit) xmult=1+beta-xmult;
                      mult*=xmult;
                    }
		}
	      if (mult>0)
		{
		  oldformap.push_back(nextoldformap);
		  oldclassdims.push_back(mult);
		  oldlevels.push_back(d);
		  noldclasses++;
		}
	    }
	}
    }
}
 
long oldforms::dimoldpart(vector<long> aplist) const
{ long ans = 0;
 if (aplist.size()==0) return 0;   // all lists "start with" a null list!
 for (long i=0; i<noldclasses; i++)
   if (startswith(oldformap[i] , aplist, aplist.size())) 
     ans += oldclassdims[i];
 if(!plusflag) ans*=2;
 return ans;
}

void oldforms::display(void) const
{
if (noldclasses>0)
  {
    long nap0=nap; if(nap0>20) nap0=20;
    cout << "\nOld classes\n~~~~~~~~~~~\n";
    cout << "Level   Dimension " << primes(nap0) << "\n";
    for (long i=0; i<noldclasses; i++)
      { 
	cout << oldlevels[i] << "       " << oldclassdims[i] << "       ";
	cout << vector<long>(oldformap[i].begin(),oldformap[i].begin()+nap0); 
	cout << "\n";
    }
  }
 cout<<"Total number of oldclasses = "<<noldclasses<<"\n";
 cout<<"Total dimension of oldclasses = "<<totalolddim<<"\n";
}
