// FILE OLDFORMS.CC: implementation of class oldforms
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

#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/cperiods.h>
#include <eclib/oldforms.h>
#include <eclib/newforms.h>

inline int testbit(long a, long i) {return (a& (1<<i));}

// Implementation of oldforms member functions

oldforms::oldforms(long intp, const level* iN, int verbose, int plus)
{  
   N=iN;
   ntp=intp;
   nap = ntp;
   noldclasses = totalolddim = 0;
   plusflag=plus;
   vector<long>::const_iterator d;
   for(d=(N->dlist).begin();d!=(N->dlist).end();d++)
     {
       long M=*d;
       if(M<11) continue;
       if(M==(N->modulus)) continue;
       getoldclasses(M,verbose);
     }
   if(verbose) cout<<"Finished getting oldclasses "<<endl;
   for (long i=0; i<noldclasses; i++) totalolddim+=oldclassdims[i];
}

void oldforms::getoldclasses(long d, int verbose) 
{
  long n = N->modulus;
  if ((d>10) && (n>d))
    {
      if(verbose) cout << "Getting oldclasses for divisor M = " << d << "\n";
      newforms olddata(d,verbose);
      olddata.createfromdata(1,25,1);
      long nforms=olddata.n1ds;
      if(nforms==0) return;
      if(olddata.nap<25)
	{
	  olddata.addap(25);
	  olddata.output_to_file();
	}
      if(verbose>1) cout << "Computing W multiplicities." << "\n";
      long m = n/d;
      long k=0, xmult, mult, j, beta, q; 
      vector<long> betalist; // =new long[N->npdivs];
      vector<long>::const_iterator qj=(N->plist).begin();
      while(qj!=(N->plist).end())
	{
	  beta=val(*qj++,m);
	  if(beta>0) k++;
	  betalist.push_back(beta);
	}
      if(verbose>1) cout<<"betas: "<<betalist<<endl;
      vector<long> nextoldformap(nap);
      vector<long>::iterator betai;
      vector<long>::const_iterator aqj;
      primevar pr; long iform, c, ip, p, aq; int bit;
      for(iform=0; iform<nforms; iform++)
	{ 
	  vector<long>& aqlist=olddata.nflist[iform].aqlist;
	  nextoldformap = olddata.nflist[iform].aplist;
	  if(verbose>1) 
	    {
	      cout<<"form #"<<(iform+1)<<": "<<"aqlist="<<aqlist<<endl;
	      cout<<"aplist before adjusting="<<nextoldformap<<endl;
	    }
	  for (c=0; c<(1<<k); c++) // 2^k different oldclasses
	    {  
	      if(verbose>1) cout<<"c="<<c<<endl;
	      mult=1; j=0; 
	      betai=betalist.begin();
	      aqj=aqlist.begin();
              for (qj=(N->plist).begin(); qj!=(N->plist).end(); qj++)
		{
		  p = *qj;
                  ip = prime_pi(p);
		  if(verbose>1) cout<<"p="<<p<<" (ip="<<ip<<")"<<endl;
                  beta = *betai++;
                  if(::divides(p,d))
                    aq = *aqj++;
                  else
                    aq=1;
                  if (beta==0)
                    {
                      nextoldformap[ip] = aq;
                      if(verbose>1)
                        cout<<"setting entry #"<<ip<<" to "<<aq<<endl;
                    }
                  else
                    {
                      bit = testbit(c,j++);
                      nextoldformap[ip] = bit?1:-1;
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
