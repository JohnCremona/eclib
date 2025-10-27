// sqfdiv.cc : implementation of class sqfdiv for managing square-free divisors
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
 
#include "eclib/bitspace.h"
#include "eclib/sqfdiv.h"

//#define DEBUG

sqfdiv::sqfdiv(const ZZ& dd, int posd, vector<ZZ>* plist)
  :primebase(plist), d(1), np(0), positive(posd), factor(0)
{
  ZZ p;
  for(unsigned long i=0; i<plist->size(); i++)
    if(p=(*primebase)[i],div(p,dd)) {d*=p; np++;}
  maxnsub=2<<np;    // = 2^(np+1) 
  maxngens=np+1; ngens=0;  nsub=1;
  subgp.resize(maxnsub);
  gens.resize(maxngens);
  pivs.resize(maxngens);
  subgp[0]=1;
  if(positive) {nsub=2; subgp[1]=-1; ngens=1; gens[0]=-1; pivs[0]=-1;}
}

void sqfdiv::usediv(const ZZ& ee)
{
  ZZ e = sqfred(ee,*primebase);
#ifdef DEBUG
  cout << "usediv called with e = " << ee 
    << " reduced mod squares to "<<e<<"\n";
#endif
  int triv=(e==1); long i;
  for(i=0; (i<ngens)&&!triv; i++)
    {
      long pivi = pivs[i];
      if(pivi==-1) {e=abs(e);}
      else { const ZZ& p = (*primebase)[pivi];
	     if(div(p,e)) e = sqfmul(e,gens[i]);
	   }
      triv = (e==1);
#ifdef DEBUG
      cout << "After using gen "<<gens[i]<<", reduced to " << e <<"\n";
#endif
    }
  if(triv) 
    {
#ifdef DEBUG
      cout << ee << " is in subgroup already...\n";
#endif
      return;
    }
#ifdef DEBUG
  cout << ee << " not in subgroup; reduced to " << e << "; adding...\n";
#endif
// update subgroup and gens:
  gens[ngens] = e;
  for(i=0; i<nsub; i++)
    {
      subgp[i+nsub]=sqfmul(subgp[i],e);
    }
  nsub*=2;

// find new pivotal prime:

  ZZ p; long valp=0, newpiv=0;
  for(i=primebase->size(); (i>0)&&(!valp); i--)
    {
      p=(*primebase)[i-1];
      if(div(p,d)) {valp=val(p,e); newpiv=i-1;}
    }
  // now if valp!=0, p|d and p|e to an odd power
  if(valp) 
    {d/=p; np--; factor++; pivs[ngens++]=newpiv;}
  else 
    if((e<0)&& !positive) 
      {positive=1; factor++; pivs[ngens++]=-1;}
  //else e is square over the support of d so should have returned earlier!
#ifdef DEBUG
  cout << "New gens:     " << vector<ZZ>(gens.begin(),gens.begin()+ngens) << endl;
  cout << "New pivs:     " << vector<long>(pivs.begin(),pivs.begin()+ngens) << endl;
  cout << "New subgroup: " << vector<ZZ>(subgp.begin(),subgp.begin()+nsub) << endl;
#endif
}

vector<ZZ> sqfdiv::getdivs() const
{
 long nd=1<<np;
 if(!positive) nd*=2;
// cout << "Constructing divisor list for d = " << d;
// cout << ", positive = " << positive;
// cout << ", number of divisors = " << nd << endl;
 vector<ZZ> dlist(nd);
 dlist[0]=1; 
 nd=1;
 if(!positive) {dlist[nd++]=-1;}
 for(unsigned long i=0; i<primebase->size(); i++)
   {
     const ZZ& p = (*primebase)[i];
     if(ndiv(p,d)) continue;
     for (long k=0; k<nd; k++)
       dlist[nd+k] = p*dlist[k];
     nd*=2;
   }
 return dlist;
}

vector<ZZ> sqfdiv::getsupp(int bothsigns) const
{
 int use_minus_one = (!positive)||bothsigns;
 vector<ZZ> supp;
 if(use_minus_one) {supp.push_back(ZZ(-1));}
 for(unsigned long i=0; i<primebase->size(); i++)
   {
     const ZZ& p = (*primebase)[i];
     if(ndiv(p,d)) continue;
     supp.push_back(p);
   }
 return supp;
}

void sqfdiv::display()
{
  cout << "Current reduced d = " << d << "\n";
  cout << "np = " << np << ", positive = " << positive << ", log_2(factor) = ";
  cout << factor << "\n";
  cout << "Subgroup gens     = " << vector<ZZ>(gens.begin(),gens.begin()+ngens) << endl;
  cout << "Subgroup elements = " << vector<ZZ>(subgp.begin(),subgp.begin()+nsub) << endl;
}

ZZ sqfred(const ZZ& a, const vector<ZZ>& plist)
{
  ZZ ans; ans=sign(a);
  for(unsigned long i=0; i<plist.size(); i++)
    {
     const ZZ& p = plist[i];
     if(odd(val(p,a))) ans*=p;
    }
  return ans;
}

ZZ sqfmul(const ZZ& a, const ZZ& b)
{ const ZZ& g = gcd(a,b);
  const ZZ& ans = (a/g)*(b/g);
  return ans;
}

ZZ makenum(const vector<ZZ>& supp, long mask)
{
  ZZ ans; ans=1; 
  long i, ns=supp.size();
  for(i=0; i<ns; i++) if(testbit(mask,i)) ans=sqfmul(ans,supp[i]);
  return ans;
}

long makeindex(const vector<ZZ>& supp, const ZZ& n, ZZ& n0)
{
  if(is_zero(n)) return 0;
  long i, ns = supp.size(), index=0;  n0=1;
  for(i=0; i<ns; i++)
    {
      ZZ pi = supp[i];
      if(sign(pi)<0) // special case, supp might have -1 as well as primes
	{
	  if(sign(n)<0) {setbit(index,i); n0=-n0;}
	}
      else // pi is really a prime
	{
	  if(odd(val(pi,n))) {setbit(index,i); n0*=pi;}
	}
    }
  return index;
}

// support(n) is like pdivs(n) but includes -1 always (except n=0, 
//but it should never be called with 0)
vector<ZZ> support(const ZZ& n)
{
  vector<ZZ> supp;
  if(is_zero(n)) { return supp;}
  vector<ZZ> supp_pos  = pdivs(n);
  supp.push_back(ZZ(-1));
  supp.insert(supp.end(),supp_pos.begin(),supp_pos.end());
  return supp;
}
