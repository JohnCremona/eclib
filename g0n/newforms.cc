// FILE  NEWFORMS.CC: implementation of newforms class
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
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

#include <iomanip>
#ifdef LiDIA_INTS
#include <LiDIA/bigint_matrix.h>
#include <LiDIA/bigmod_matrix.h>
#endif
#include "interface.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "cperiods.h"
#include "newforms.h"

// Functions for ordering newforms
// (1) Old ordering (first aq, then ap for good p);
// with the order for eigenvalues being
// 1,-1 or 0,1,-1,2,-2,...
// (2) New ordering (ap for all p in natural order)
// with the order for eigenvalues being
// -1,1 or ...,-2,-1,0,1,2,... (plain numerical order)

// less_ap(a,b) returns +1,0,-1 according to whether a is
// before/equal/after b in the above ordering
int less_ap(long a, long b, int old=0) 
{ 
  if(!old) return sign(b-a); // the simple new ordering!
  if(a==b) return 0;
  int s = sign(abs(b)-abs(a));
  if(s) return s;
  else return sign(a-b); // this way round! +1 before -1 etc
}

int less_apvec(const vector<long>& v, const vector<long>& w, int old=0); 
int less_apvec(const vector<long>& v, const vector<long>& w, int old) 
{
  vector<long>::const_iterator vi=v.begin(), wi=w.begin();
  while(vi!=v.end())
    {
      int s = less_ap(*vi++,*wi++,old); 
      if(s) return s;
    }
  return 0;
}

struct less_newform_old : public binary_function<newform, newform, bool> { 
  bool operator()(const newform& f, const newform& g) 
  {
    int s = less_apvec(f.aqlist,g.aqlist,1);
    if(s==0) s = less_apvec(f.aplist,g.aplist,1);
    return (s==1);
  }
};

struct less_apvec_function : public binary_function<const vector<long>&, const vector<long>&, bool> { 
  bool operator()(const vector<long>& f, const vector<long>& g) 
  {
    return 1==less_apvec(f,g);
  }
};

vector<long> eiglist(const newform& f, int oldorder)
{
  /*
  cout<<"Entering eiglist with f.aqlist="<<f.aqlist<<"\nand f.aplist=";
  vec_out(cout,f.aplist,10);
  cout<<endl;
  */
  vector<long> eigs;
  primevar pr;
  long N = (f.nf)->modulus;
  vector<long>::const_iterator aqi=f.aqlist.begin();
  vector<long>::const_iterator api=f.aplist.begin();
  vector<long>::iterator eigsi;
  if(oldorder)
    {
      eigs.resize(f.aplist.size());
      eigsi=eigs.begin();
      while(aqi!=f.aqlist.end()) 
	*eigsi++ = *aqi++;
      while(api!=f.aplist.end())
	{
	  if(!::div(pr,N)) *eigsi++ = *api;     
	  api++; pr++;
	}
    }
  else
    {
      eigs=f.aplist; // copy; now adjust the aq:
      eigsi=eigs.begin();
      while((aqi!=f.aqlist.end())&&(eigsi!=eigs.end()))
	{
	  if(::div(pr,N)) *eigsi = (*aqi++);
	  eigsi++; pr++;
	}
    }
  /*
  cout<<"Leaving eiglist with eigs=";
  vec_out(cout,eigs,10);
  cout<<endl;
  */
  return eigs;
}

struct less_newform_new : public binary_function<newform, newform, bool> { 
  bool operator()(const newform& f, const newform& g) 
  {
    //    return less_apvec(eiglist(f),eiglist(g),0)==1;    
    return less_apvec(f.aplist,g.aplist,0)==1;    
  }
};

// Newform constructor given the ap and aq lists and extra data (but
// no homology basis), e.g. after reading from newforms file

newform::newform(const vector<int>& data, const vector<long>& aq, const vector<long>& ap, newforms* nfs) : nf(nfs)
{
  sfe=data[0];
  ap0=data[1];
  np0=data[2];
  dp0=data[3];
  loverp=rational(dp0,np0);
  lplus=data[4];
  mplus=data[5];
  lminus=data[6];
  mminus=data[7];
  a=data[8];
  b=data[9];
  c=data[10];
  d=data[11];
  dotplus=data[12];
  dotminus=data[13];
  type=data[14];
  degphi=data[15];
  aqlist=aq;
  aplist=ap;
  index=-1;
  pdot=0;
}

// Newform constructor, given the homology basis vector(s) and
// Hecke eigenvalues

newform::newform(const vec& vplus, const vec& vminus, const vector<long>& ap, newforms* nfs,long ind)
   :nf(nfs), sign(nfs->sign), bplus(vplus),bminus(vminus),aplist(ap),index(ind) 
{
  int verbose=(nf->verbose);
  long n = nf->modulus;
  
  if(verbose) 
    {
      cout<<"Creating H1"; 
      if(sign==+1) cout<<"+"; 
      if(sign==-1) cout<<"-"; 
      cout<<" newform from aplist..."<<endl;
      if(verbose>2)
        {
          if(sign!=-1) cout<<"bplus = "<<bplus<<endl;
          if(sign!=+1) cout<<"bminus = "<<bminus<<endl;
        }
    }
  
  // check_expand_contract();

  // Fixing the eigenvalue lists: ap is indexed by primes in natural
  // order we need to extract aq (computing any not yet there).
  
  // At the same time we change the entries in aplist for bad primes q
  // from the Wq-eigenvalue to the newform coefficient.
  fixup_eigs();

  // Compute cuspidalfactors and type (only if sign=0):

  type = 0;
  find_cuspidal_factors();

  // Compute projected coordsplus/minus and denomplus/minus

  find_coords_plus_minus();
  
  // Compute pdot, dp0, loverp (unless sign is -1)

  find_bsd_ratio();
  
  // Find deg(phi) (only if sign is 0)
  degphi = 0;
  find_degphi();

  // Find twisting primes if N non-square

  lplus=mplus=0;
  lminus=mminus=0;
  find_twisting_primes();

  // find a,b,c,d,dotplus,dotminus

  a=b=c=d=0;
  dotplus=dotminus=0;
  find_matrix();
}

int newform::check_expand_contract()
{
  int success=1;
  long denom = nf->h1->h1denom();
  vec bplusx, bminusx, tvec;
  if (sign!=-1) 
    {
      bplusx= nf->h1->extend_coords(bplus);
      tvec = nf->h1->contract_coords(bplusx);
      tvec /= denom;
      if (tvec!=bplus)
	{
	  success=0;
	  cout<<"! bplus ="<<bplus<<" extends to "<<bplusx<<" which contracts to "<<tvec<<endl;
	}
    }
  if (sign!=+1) 
    {
      bminusx= nf->h1->extend_coords(bminus);
      tvec = nf->h1->contract_coords(bminusx);
      tvec /= denom;
      if (tvec!=bminus)
	{
	  success=0;
	  cout<<"! bminus="<<bminus<<"  extends to "<<bminusx<<" which contracts to "<<tvec<<endl;
	}
    }
  return success;
}

void newform::fixup_eigs()
{ 
  long denom = nf->h1->h1denom();
  aqlist.resize(nf->npdivs);
  vector<long>::iterator api=aplist.begin(), pi=nf->plist.begin();
  vector<long>::iterator aqi=aqlist.begin();
  primevar pr;   long q, i, j;
  long n = nf->modulus;
  while((api!=aplist.end())&&(aqi!=aqlist.end()))
    {
      q=pr.value(); pr++;
      if(::div(q,n)) 
	{
	  *aqi++=*api; 
	  *api=(::div(q*q,n)? 0: -*api);
	  pi++;
	}
      api++;
    }
  if(aqi!=aqlist.end()) // compute missing aq
    {
      long piv;
      ssubspace espace;
      if(sign==-1)
        espace=make1d(bminus,piv); 
      else
        espace=make1d(bplus,piv); 
      piv*=denom;
      while(aqi!=aqlist.end()) // compute missing aq
	{
	  q=*pi++;
	  if(nf->verbose) cout<<"Computing Wq for q="<<q<<"..."<<flush;
	  smat Wq = nf->h1->s_heckeop_restricted(q,espace,1,0);
	  long aq = Wq.elem(1,1) / piv;
	  if(nf->verbose) cout<<"aq ="<<aq<<endl;
	  *aqi++=aq;
	}
    }
  if(nf->verbose) cout<<"aqlist = "<<aqlist<<endl;
  
  //Compute sfe:

  sfe=-1;
  for(i=0; i<(nf->npdivs); i++) sfe*=aqlist[i];
  if(nf->verbose) cout<<"sfe = "<<sfe<<endl;
}

// Before recovering eigenbases, we need to put back the aq into the
// aplist (and resort, for efficiency).
void newform::unfix_eigs()
{ 
  vector<long>::iterator api=aplist.begin();
  vector<long>::iterator aqi=aqlist.begin();
  primevar pr;
  long n = nf->modulus;
  while((api!=aplist.end())&&(aqi!=aqlist.end()))
    {
      if(::div(pr.value(),n)) *api=*aqi++; 
      api++;
      pr++;
    }
}

// After recovering eigenbases, we need to replace the ap for bad p
void newform::refix_eigs()
{ 
  vector<long>::iterator api=aplist.begin();
  primevar pr;
  long n = nf->modulus, np = nf->npdivs, ip=0, q;
  while((api!=aplist.end())&&(ip<np))
    {
      q=pr.value();
      if(::div(q,n)) 
	{
	  *api=(::div(q*q,n)? 0: -*api);
	  ip++;
	}
      api++;
      pr++;
    }
}

void newform::find_bsd_ratio()
{
  // get ap for p=p0:

  primevar pr; 
  vector<long>::const_iterator api=aplist.begin();
  while(pr.value()!=nf->p0) {pr++; api++;}
  ap0=*api;
  np0 = 1 + (nf->p0) - ap0;
  if(nf->verbose) cout<<"ap0 = "<<ap0<<"\tnp0 = "<<np0<<endl;

  if(sign==-1) return;

  // DO NOT scale pdot by denom: factor will cancel when used to compute ap
  pdot = abs((nf->mvp)*bplus); 
  dp0=pdot;
  // DO scale dp0 since it is used to compute L/P
  if(dp0!=0)
    {
      if(denomplus>1)
	{
	  if(::div(denomplus,dp0))  dp0/=denomplus;
	  else 
	    cout<<"newform constructor error: dp0 not divisible by denomplus!"
		<<endl; 
	}
    }
  loverp = rational(dp0,np0);   
  if(nf->verbose) 
    {
      cout<<"pdot = "<<pdot<<endl;
      cout<<"dp0 = "<<dp0<<endl;
      cout<<"np0 = "<<np0<<endl;
      cout<<"loverp = "<<loverp<<endl;
    }
}

void newform::find_coords_plus_minus()
{
  int verbose = nf->verbose;
  int i, ncoords=nf->h1->coord_vecs.size()-1;
  svec cvi;
  if(sign!=-1)
    coordsplus=vec(ncoords);
  if(sign!=+1)  
    coordsminus=vec(ncoords);
  //  if(verbose) cout<<"About to compute coordsplus/minus"<<endl;
  for(i=1; i<=ncoords; i++)
    {
      cvi = nf->h1->coord_vecs[i];
      if(sign!=-1)
        coordsplus[i]=dotmodp(cvi,bplus,92681);
      if(sign!=+1) 
        coordsminus[i]=dotmodp(cvi,bminus,92681);
    }

  if(sign!=+1) 
    {
      denomminus=vecgcd(coordsminus)*cuspidalfactorminus;
      if(verbose) cout<<"denomminus   = "<<denomminus<<endl;
    }
  if(sign!=-1) 
    {
      denomplus=vecgcd(coordsplus)*cuspidalfactorplus;
      if(verbose) cout<<"denomplus   = "<<denomplus<<endl;
    }  
}

void newform::find_cuspidal_factors()
{
  vec bplusc, bminusc;
  int verbose = nf->verbose;

  if(!(nf->h1->cuspidal))
    {
      if(sign!=-1) // do this if sign = 0,1
        {
          bplusc=(nf->h1->tkernbas)*bplus;
          cuspidalfactorplus = vecgcd(bplusc);
          bplusc /= cuspidalfactorplus;
        }
      if(sign!=+1) // do this if sign = 0,-1 
	{
	  bminusc=(nf->h1->tkernbas)*bminus;
	  cuspidalfactorminus = vecgcd(bminusc);
	  bminusc/= cuspidalfactorminus;
        }
      if(sign==0)  // do this only if sign = 0
        {
	  type=3-vecgcd(bplusc-bminusc);
	  if(verbose) cout<<"Lattice type = "<<type<<endl;
        }
   
      if(verbose&&(cuspidalfactorplus*cuspidalfactorminus>1))
	{
          if(sign!=-1)
            {
              cout<<"cuspidalfactorplus  = "<<cuspidalfactorplus<<endl; 
              if(verbose>2) cout<<"bplusc = "<<bplusc<<endl;
            }
	  if(sign!=+1) 
            {
              cout<<"cuspidalfactorminus = "<<cuspidalfactorminus<<endl;
              if(verbose>2) cout<<"bminusc = "<<bminusc<<endl;
            }
	}
    } 
}
  
void newform::find_degphi()
{
  if(sign!=0) return;
#ifdef DEG_PHI
    if(nf->verbose) cout<<"computing deg(phi)..."<<flush;
    degphi=jumpinfo->degphi(bplusc,bminusc,type);
    if(nf->verbose) cout<<"done..."<<flush;
#else
    degphi=0;
#endif
}

void newform::find_twisting_primes()
{
  int verbose=(nf->verbose);
  if(verbose) cout<<"computing twisting primes (sign="<<sign<<")..."<<flush;
  if(sign!=-1) 
    if(dp0!=0)
      {
	lplus=1; // so we need not search for a prime 1(mod 4) below
	mplus=1; // dummy value, not used
      }
    else
      {
	lplus=0;
	mplus =0;
      }
  if(sign!=+1) 
    {
      lminus=0;
      mminus=0;
    }  
  if(nf->squarelevel) return;

  long n = nf->modulus;
  for (primevar lvar; lvar.ok() && 
	 (((sign!=-1)&&(mplus==0)) ||
	  ((sign!=+1)&&(mminus==0))); lvar++)
    {
      //cout << "Trying l = " << lvar << endl;
      while (n%lvar==0) {lvar++;}
      long l = lvar;
      //cout << "Trying l = " << l << endl;
      if (legendre(-n,l)!=sfe) continue;
      //cout << "Legendre condition passed... " << endl;

      if((sign!=-1)&&(mplus==0)&&(l%4==1))
	{
	  lplus = l;   
	  //cout << "Trying lplus = " << l << "\n";
	  map<long,vec>::const_iterator vi = nf->mvlplusvecs.find(l);
	  if(vi==nf->mvlplusvecs.end())
	    mplus = abs((nf->mvlplusvecs[l]=nf->h1->manintwist(l))*bplus);
	  else
	    mplus = abs((vi->second)*bplus);
	  if((denomplus>1)&&(mplus!=0))
	    {
	      if(::div(denomplus,mplus))  mplus/=denomplus;
	      else 
		cout<<"Warning in newform constructor: mplus not divisible by denomplus!"
		    <<endl; 
	    }
	}
      if((sign!=+1)&&(mminus==0)&&(l%4==3))
	{
	  lminus = l;   
	  //cout << "Trying lminus = " << l << "\n";
	  map<long,vec>::const_iterator vi = nf->mvlminusvecs.find(l);
	  if(vi==nf->mvlminusvecs.end())
	    mminus = abs((nf->mvlminusvecs[l]=nf->h1->manintwist(l))*bminus);
	  else
	    mminus = abs((vi->second)*bminus);
	  if((denomminus>1)&&(mminus!=0))
	    {
	      if(::div(denomminus,mminus))  mminus/=denomminus;
	      else 
		cout<<"Warning in newform constructor: mminus="<<mminus<<" is not divisible by denomminus="<<denomminus<<"!"
		    <<endl; 
	    }
	}
    }

  if(verbose) 
    {
      cout<<"done..."<<flush;
      cout<<"lplus = "<<lplus<<endl;
      cout<<"mplus = "<<mplus<<endl;
      cout<<"lminus = "<<lminus<<endl;
      cout<<"mminus = "<<mminus<<endl;
    }
}

void newform::find_matrix()
{
  int verbose=(nf->verbose);
  if(verbose) cout<<"computing a,b,c,d..."<<flush;
  long n = nf->modulus;
  int found=0;
  vec v;
  for(d=2; !found; d++)
    {
      if(1==gcd(d,n))
        {
          for(b=1; (b<d) && !found; b++)
            {
              if(1==bezout(d,-n*b,a,c))
                {
		  //                  cout<<"b/d = "<<b<<"/"<<d<<": ";
                  v = nf->h1->chain(b,d).as_vec();
		  //                  cout<<"v="<<v<<endl;
//                if(!(nf->h1->cuspidal)) v = nf->h1->cuspidalpart(v);
                  if(sign!=-1)
                    {
                      dotplus=v*bplus;
                      if(::div(denomplus,dotplus))  
                        dotplus/=denomplus;
                      else 
                        cout<<"Warning in find_matrix: dotplus not divisible by denomplus!"<<endl; 
                      dotplus=abs(dotplus);
                    }
                  if(sign!=+1)
                    {
                      dotminus=v*bminus;
                      if(::div(denomminus,dotminus))  
                        dotminus/=denomminus;
                      else 
                        cout<<"Warning in find_matrix: dotminus not divisible by denomminus!"<<endl; 
                      dotminus=abs(dotminus);
                    }
                  found=(((dotplus!=0)||(sign==-1))&&
                         ((dotminus!=0)||(sign==+1)));
                }
            }
        }
    }
  b--; d--;  //because they get incremented BEFORE the loop end-test
  if(d<0) {a=-a; b=-b; c=-c; d=-d;} // because we need d>0 for integration

  if(verbose) 
    {
      cout<<"done: ";
      cout << "[(" <<a<<","<<b<<";"<<c
	   <<","<<d<<"),"<<dotplus<<","<<dotminus
	   <<";"<<type<<"]"<<endl;
    }
}

void newform::add_more_ap(int nap)
{
  if(aplist.size()>=nap) return;
  int verbose=(nf->verbose);
  long piv, p, ap;
  // Do not make the espace right away, as it is possible that the
  // only ap we are missing are aq which we already have...
  ssubspace espace;
  int have_espace=0;

  primevar pr(nap,aplist.size()+1);
  while(aplist.size()<nap)
    {
      p=pr;
      if(::div(p,nf->modulus))
	{
	  if(::div(p*p,nf->modulus)) 
	    ap=0; 
	  else
	    ap=-aqlist[find(nf->plist.begin(),nf->plist.end(),p)-nf->plist.begin()];
	}
      else
	{	  
	  if(verbose>1) cout<<"Computing Tp for p="<<p<<endl;
	  if(!have_espace)
	    {
              if(sign==-1)
                espace=make1d(bminus,piv);
              else
                espace=make1d(bplus,piv);
	      piv*=nf->h1->h1denom();
	      have_espace=1;
	    }
	  ap = (nf->h1->s_heckeop_restricted(p,espace,1,0)).elem(1,1) / piv;
	}
      aplist.push_back(ap);
      pr++;
    }
  if(verbose>1) cout<<"aplist = "<<aplist<<endl;
}

newforms::~newforms(void)
{
  delete of;  
  delete h1plus;
  delete h1minus;
  delete h1full;
}

void newforms::makeh1(int s)
{
  if(s==1)
    {
      if(!h1plus)
	{
	  if(verbose) cout<<"Constructing H1 (with sign=+1) ..."<<flush;
	  h1plus = new homspace(modulus,1,0,0 /*verbose*/);
	  if(verbose) cout<<"done"<<endl;
	}
      h1 = h1plus;
      return;
    }
  if(s==-1)
    {
      if(!h1minus)
	{
	  if(verbose) cout<<"Constructing H1 (with sign=-1) ..."<<flush;
	  h1minus = new homspace(modulus,-1,0,0 /*verbose*/);
	  if(verbose) cout<<"done"<<endl;
	}
      h1 = h1minus;
      return;
    }
  if(s==0)
    {
      if(!h1full)
	{
	  if(verbose) cout<<"Constructing H1 (with sign=0) ..."<<flush;
	  h1full = new homspace(modulus,0,0,0 /*verbose*/);
	  if(verbose) cout<<"done"<<endl;
	}
      h1 = h1full;
      return;
    }
  cout<<"Error in makeh1(s): s = "<<s<<" should be one of 0,1,-1"<<endl;
  return;
}

void newforms::createfromscratch(int s, long ntp)
{
  sign = s;
  makeh1(s);
  //  cout<<"Constructing oldforms with sign="<<sign<<endl; 
  of = new oldforms(ntp,h1,(verbose>1),sign); // h1 provides the level*
  if(verbose>1) of->display();
  maxdepth = of->nap;
  long mindepth = npdivs; // must include at least one good p, and it
			  // is cheap to continue recursing after
			  // reaching dimension 1.
  n1ds = 0;
  int upperbound = h1->dimension-(of->totalolddim);
  if(upperbound>0)  // Else no newforms certainly so do no work!
    {  
       mvp=h1->maninvector(p0); 
       //       cout<<"mvp                 = "<<mvp<<endl;
       if(verbose>1) cout<<"h1 denom = "<<h1->h1denom()<<endl;
       long totalmult=upperbound;
       if (totalmult==0) n1ds=0;      
       else 
	 {
	   form_finder ff(this,(sign!=0),maxdepth,mindepth,1,0,verbose);
	   basisflag=0;
	   ff.find();
	 }
    }
  if(verbose)
    {
      cout << "Total dimension = " << h1->dimension << endl;
      cout << "Number of rational newforms = " << n1ds <<endl;
      if(h1->dimension==of->totalolddim+n1ds)
	cout<<"The whole space splits over Q" << endl;
    }

  if(n1ds==0) return;

  int i,j,nap,maxnap=0;

  if((n1ds>1)&&(modulus<130000)) // reorder into old order
    {
      if(verbose) cout<<"Reordering newforms into old order as N<130000"<<endl;
      // if(verbose) cout<<"Before sorting:\n"; display();
      sort(1);
      // if(verbose) cout<<"After sorting:\n"; display();
    }
  // At this point the newforms may contain different numbers of ap,
  // so we need to even these up, which we do by computing more ap for
  // those which need it.
   if(n1ds>1)
    {
      for(i=0; i<n1ds; i++) 
	if((nap=nflist[i].aplist.size())>maxnap) maxnap=nap;
      if(verbose) 
	cout<<"Max number of ap in newforms so far is "<<maxnap<<endl;
      for(i=0; i<n1ds; i++) 
	if((nap=nflist[i].aplist.size())<maxnap) 
	  {
	    if(verbose)
	      cout<<"Newform #"<<(i+1)<<" has only "<<nap
		  <<" ap so we need to compute more..."<<endl;
	    nflist[i].add_more_ap(maxnap);
	  }  
    }

   // Compute homspace::projcoord, so projchain can be used
   // Replaces coord_vecs of homspace with projections onto eigenspaces
   make_projcoord();

   // Look for a j0 such that nflist[i].bplus/bminus[j0]!=0 for all i, or a set of such j
   find_jlist();
}  


void newforms::find_jlist()
{
  int i, j, ok=0; j0=0;
  for(j=1; (!ok)&&(j<=h1->h1dim()); j++)
    {
      ok=1;
      for (i=0; (i<n1ds)&&ok; i++)
        if(sign==-1)
          ok=(nflist[i].bminus[j]!=0);
        else
          ok=(nflist[i].bplus[j]!=0);
      if(ok) j0=j;
    }
  if(ok)
    {
      if(verbose>1)  cout<<"j0="<<j0<<endl;
      jlist.insert(j0);
      for (i=0; i<n1ds; i++) 
	{
	  nflist[i].j0 = j0;
	  if(sign==-1)
	    nflist[i].fac = nflist[i].bminus[j0];
	  else
	    nflist[i].fac = nflist[i].bplus[j0];
	}
    }
  else
    {
      if(verbose) 
	cout<<"Failed to find j0 such that nflist[i].bplus/bminus[j0]!=0 for all i"
	    <<endl;
      // Find out which pivots we'll be using:
      for (i=0; i<n1ds; i++)
	{
	  vec& bas = nflist[i].bplus;
	  j=1; while(bas[j]==0) j++;
	  jlist.insert(j);
	  nflist[i].j0 = j;
	  nflist[i].fac = nflist[i].bplus[j];
	}
      if(verbose)  cout<<"jlist="<<jlist<<endl;
    }
}

// Compute homspace::projcoord, so projchain can be used
// Replaces coord_vecs of homspace with projections onto eigenspaces
void newforms::make_projcoord()
{
  h1->projcoord.init(h1->coord_vecs.size()-1,n1ds);
  int j;
  if(sign==-1)
    for (j=1; j<=n1ds; j++)
      h1->projcoord.setcol(j, nflist[j-1].coordsminus);
  else
    for (j=1; j<=n1ds; j++)
      h1->projcoord.setcol(j, nflist[j-1].coordsplus);
}

long newforms::dimoldpart(const vector<long> l) 
{
  return of->dimoldpart(l);
}

// if(!cuspidal) we really should check here that the basis vector b1
// is in ker(delta), by checking that b1*h1->deltamat == 0

void newforms::use(const vec& b1, const vec& b2, const vector<long> aplist)
{
  if(basisflag) // we already have all the data except the
                // basis vector, so not much needs doing
    {
      if(verbose) 
	cout<<"Filling in data for for newform #"<<(j1ds+1)<<": bases..."<<flush;
      nflist[j1ds].sign=sign;
      if(sign==+1)
        nflist[j1ds].bplus=b1;
      if(sign==-1)
        nflist[j1ds].bminus=b1; // formfinder puts the basis vector in b1
      if(sign==0) 
	{
	  nflist[j1ds].bplus=b1;
	  nflist[j1ds].bminus=b2;
	}
      if(verbose) 
	cout<<"type and cuspidal factors..."<<flush;
      nflist[j1ds].find_cuspidal_factors();
      if(verbose) 
	cout<<"coords..."<<flush;
      nflist[j1ds].find_coords_plus_minus();
      if(sign==0)
	{
	  if(verbose) 
	    cout<<"twisting primes..."<<flush;
	  nflist[j1ds].find_twisting_primes();
	  if(verbose) 
	    cout<<"matrix..."<<flush;
	  nflist[j1ds].find_matrix();
	}
      if(verbose) 
	cout<<"done."<<endl;
      j1ds++;
      if(verbose) 
	cout<<"Finished filling in data for newform #"<<j1ds<<endl;
      return;
    }

  // Code for initial newform construction

  // We use the newform constructor to do all the work, given the basis vector(s) and aplist:

  n1ds++;
  if(verbose) 
    {
      cout<<"Constructing newform #"<<n1ds<<" with eigs ";
      vec_out(cout,aplist,10);
      cout<<endl;
    }
  if(sign==-1)
    nflist.push_back(newform(b1,b1,aplist,this)); // only 2nd vector used
  else
    nflist.push_back(newform(b1,b2,aplist,this));
  if(verbose) 
    cout<<"Finished constructing newform #"<<n1ds<<" with sign = "<<sign<<endl;
}

// Sort newforms 
void newforms::sort(int oldorder)
{
  if(oldorder)
    ::sort(nflist.begin(),nflist.end(),less_newform_old());
  else
    ::sort(nflist.begin(),nflist.end(),less_newform_new());
}
  
// Before recovering eigenbases, we need to put back the aq into the
// aplist (and resort, for efficiency).
void newforms::unfix_eigs()
{
  for(int i=0; i<n1ds; i++)
    nflist[i].unfix_eigs();
}

// After recovering eigenbases, we need to refix the aplist
void newforms::refix_eigs()
{
  for(int i=0; i<n1ds; i++)
    nflist[i].refix_eigs();
}

void newforms::display(void) const
{
 if (n1ds==0) {cout<<"No newforms."<<endl; return;}
 cout << "\n"<<n1ds<<" newform(s) at level " << modulus << ":" << endl;
 cout<<"p0="<<p0<<endl;
 // if(dim(mvp)!=0) cout<<"mvp="<<mvp<<endl;
 cout<<"#ap=\t"<<nflist[0].aplist.size()<<endl;
 long i;
 for(i=0; i<n1ds; i++)
   {cout<<i+1<<":\t";
    nflist[i].display();
  }
}

void newforms::display_modular_symbol_map(void) const
{
 long i,j,k,m;
 rational rplus, rminus;
 for(i=0; i<h1->nsymb; i++)
   {
     symb s = h1->symbol(i);
     cout<<s<<" = "<<modsym(s)<<" -> ";
     j=h1->coordindex[i];
     int sg=::sign(j); j=abs(j);
     //     cout<<"j="<<j<<"("<<sg<<")"<<endl;
     if(j==0) 
       for(k=0; k<n1ds; k++) 
	 if(sign!=0) 
	   cout<<"0 "; 
	 else 
	   cout<<"(0,0) ";
     else 
       for(k=0; k<n1ds; k++) 
         {
           //           cout<<"coordsplus = "<<nflist[k].coordsplus<<endl;
           //           cout<<"coordsminus = "<<nflist[k].coordsminus<<endl;
           if(sign!=-1)
             rplus = rational(sg*nflist[k].coordsplus[j],nflist[k].cuspidalfactorplus);
           if(sign!=+1)
             rminus = rational(sg*nflist[k].coordsminus[j],nflist[k].cuspidalfactorminus);
           if(sign==+1) 
             cout<<rplus<<" ";
           else if(sign==-1) 
             cout<<rminus<<" ";
           else
             cout<<"("<<rplus<<","<<rminus<<") ";
         }
     cout<<endl;
   }
}

void newform::display(void) const
{
  cout << "aplist = ";  
  vec_out(cout,aplist,20);  // outputs at most 20 eigs.
  cout<< endl;
  //  cout << "basis = " << bplus << endl;
  cout << "aq = " << aqlist<<endl;
  cout << "ap0 = " << ap0
       <<", dp0 = " << dp0
       <<", np0 = " << np0;
  if(pdot!=0) cout <<", pdot = " << pdot;
  cout <<endl;
  cout << "SFE = " << sfe << ",\tL/P = " << loverp << endl;
  if(lplus>0) cout << "lplus = " << lplus << ", mplus = " << mplus << endl;
  if(lminus>0) cout << "lminus = " << lminus << ", mminus = " << mminus << endl;
  if(a!=0) 
    {
      cout << "[(" <<a<<","<<b<<";"<<c
           <<","<<d<<"),"<<dotplus<<","<<dotminus
           <<";";
      if(type) 
        cout<<type;
      else 
        cout<<"?";
      cout<<"]"<<endl;
    }
    
  if(index!=-1)cout << "Splitting index = " << index << endl;
}

void putout(ofstream& of, short a, int binflag)
{
  if(binflag)
    of.write((char*)&a,sizeof(short));
  else
    of<<setw(5)<<a;
}
void putout(ofstream& of, int a, int binflag)
{
  if(binflag)
    of.write((char*)&a,sizeof(int));
  else
    of<<setw(10)<<a;
}
void putout(ofstream& of, long a, int binflag)
{
  if(binflag)
    of.write((char*)&a,sizeof(long));
  else
    of<<setw(15)<<a;
}
void nl(ofstream& of, int binflag) 
{if(!binflag) of<<"\n";}

void newforms::output_to_file(int binflag) const
{
  long i,j;
  char* name;
  if(binflag) 
    name = nf_filename(modulus,'x');
  else 
    name = nf_filename(modulus,'e');
  ofstream out(name);
  if(!out)
    {
      cout<<"Unable to open file "<<name<<" for newform output"<<endl;
      delete[] name;
      abort();
    }
  delete[] name;
  
  if(n1ds==0) 
    {  
      putout(out,(int)0,binflag);
      putout(out,(int)0,binflag);
      putout(out,(int)0,binflag);
      out.close(); 
      return;
    }
    
  // Line 1:  #newforms, #aq, #ap
  putout(out,(int)n1ds,binflag);
  putout(out,(int)nflist[0].aqlist.size(),binflag);
  putout(out,(int)nflist[0].aplist.size(),binflag);
  nl(out,binflag);
  // Line 2:  blank line
  nl(out,binflag);
  // Line 3:  sign of f.e. for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].sfe,binflag); 
  nl(out,binflag);
  // Line 4:  ap0 for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].ap0,binflag); 
  nl(out,binflag);
  // Line 5:  np0 for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].np0,binflag); 
  nl(out,binflag);
  // Line 6:  dp0 for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].dp0,binflag); 
  nl(out,binflag);
  // Line 7:  lplus for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].lplus,binflag);	    
  nl(out,binflag);
  // Line 8:  mplus each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].mplus,binflag);	
  nl(out,binflag);
  // Line 9:  lminus each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].lminus,binflag);	
  nl(out,binflag);
  // Line 10:  mminus for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].mminus,binflag);	
  nl(out,binflag);
  // Line 11:  matrix entry a for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].a,binflag); 
  nl(out,binflag);
  // Line 12:  matrix entry b for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].b,binflag);	
  nl(out,binflag);
  // Line 13:  matrix entry c for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].c,binflag);	
  nl(out,binflag);
  // Line 14:  matrix entry d for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].d,binflag);	
  nl(out,binflag);
  // Line 15:  dotplus for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].dotplus,binflag); 
  nl(out,binflag);
  // Line 16:  dotminus for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].dotminus,binflag); 
  nl(out,binflag);
  // Line 17:  lattice type for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].type,binflag); 
  nl(out,binflag);
  // Line 18:  deg(phi) for each newform
  for(i=0; i<n1ds; i++) putout(out,(int)nflist[i].degphi,binflag);	
  nl(out,binflag);
  // Line 19:  blank line
  nl(out,binflag);
  
  // Lines 20-(20+#aq):  aq for each newform;  then blank line
  for(j=0; j<nflist[0].aqlist.size(); j++)
    {
      for(i=0; i<n1ds; i++) putout(out,(short)nflist[i].aqlist[j],binflag); 
      nl(out,binflag);
    }
  nl(out,binflag);
  
  // Lines (21+#aq)-(20+#aq+#ap):  ap for each newform
  for(j=0; j<nflist[0].aplist.size(); j++)
    {
      for(i=0; i<n1ds; i++) putout(out,(short)nflist[i].aplist[j],binflag); 
      nl(out,binflag);
    }
  out.close();
}

// Read in newform data from file NF_DIR/xN 

void newforms::createfromdata(int s, long ntp, int create_from_scratch_if_absent)
{
  sign = s;
  long i, j, n = modulus;
  if(verbose) cout << "Retrieving newform data for N = " << n << endl;

  char* name = nf_filename(modulus,'x');
  ifstream datafile(name);
  if(!datafile.is_open())
    {
      if(verbose) cout<<"Unable to open file "<<name<<" for newform input"<<endl;
      delete[] name;
      if(create_from_scratch_if_absent)
	{
	  if(verbose) cout<<"Creating from scratch instead"<<endl;
	  createfromscratch(sign, ntp);
	  output_to_file();
	  if(verbose) cout << "Finished creating newform data for N = " << n << endl;
	  if(verbose) display();
	  return;
	}
      else
	{
	  cout<<"Quitting"<<endl;
	  abort();
	}
    }
  delete[] name;

  short temp_short;
  int temp_int;
  datafile.read((char*)&temp_int,sizeof(int));   // = number of newforms
  n1ds=temp_int;
  datafile.read((char*)&temp_int,sizeof(int));   // = number of bad primes
  datafile.read((char*)&temp_int,sizeof(int));   // = number of eigs
  nap=temp_int;
  if(n1ds==0) 
    {
      if(verbose) cout << "No newforms at level " << n << endl;
      datafile.close();
      return;
    }

  vector<int> * data = new vector<int>[n1ds];
  vector<long> * aq = new vector<long>[n1ds];
  vector<long> * ap = new vector<long>[n1ds];

  // read extra data for each newform
  for(i=0; i<n1ds; i++) data[i].resize(16);
  long ntotal = 16*n1ds;
  int* batch_i = new int[ntotal];
  datafile.read((char*)batch_i,ntotal*sizeof(int));
  int *batch_i_ptr = batch_i;
  for(j=0; j<16; j++) 
    for(i=0; i<n1ds; i++) 
      data[i][j]=*batch_i_ptr++;
  delete[] batch_i;
  //  cout<<"Raw  data:\n";  for(i=0; i<n1ds; i++) cout<<data[i]<<endl;

  // read aq for each newform
  for(i=0; i<n1ds; i++) aq[i].resize(npdivs);
  ntotal = npdivs*n1ds;
  short* batch = new short[ntotal];
  datafile.read((char*)batch,ntotal*sizeof(short));
  short *batchptr = batch;
  for(j=0; j<npdivs; j++) 
    for(i=0; i<n1ds; i++) 
      aq[i][j]=*batchptr++;
  //  cout<<"Raw  aq:\n";  for(i=0; i<n1ds; i++) cout<<aq[i]<<endl;

  // read ap for each newform
  for(i=0; i<n1ds; i++) ap[i].resize(nap);
  ntotal = nap*n1ds;
  delete[] batch;
  batch = new short[ntotal];
  datafile.read((char*)batch,ntotal*sizeof(short));
  batchptr = batch;
  for(j=0; j<nap; j++) 
    for(i=0; i<n1ds; i++) 
      ap[i][j]=*batchptr++;
  //  cout<<"Raw  ap:\n";  for(i=0; i<n1ds; i++) cout<<ap[i]<<endl;
  delete[] batch;
  datafile.close();

  // construct the newforms from this data
  nflist.reserve(n1ds);
  for(i=0; i<n1ds; i++) 
    nflist.push_back(newform(data[i],aq[i],ap[i],this));

  delete[] ap; delete[] aq; delete[] data;

  if(verbose) 
    {
      cout << "Finished reading newform data for N = " << n << endl;
      if(verbose>1) display();
    }
}


// Create from a list of Hecke eigenvalues from an elliptic curve

vector<long> eiglist(CurveRed& C, int nap)
{
  long N = I2long(getconductor(C));
  long p; bigint pp;
  vector<long> ans;
  for(primevar pr(nap); pr.ok(); pr++)
    {
      p=pr; pp=BIGINT(p);
      if(N%p==0)
	ans.push_back(LocalRootNumber(C,pp));
      else
	ans.push_back(I2long(Trace_Frob(C,pp)));
    }
  //  cout<<"eiglist("<<(Curve)C<<") = "<<ans<<endl;
  return ans;
}

// Create from a list of elliptic curves of the right conductor:

void newforms::createfromcurve(int s, CurveRed C, int nap)
{
  vector<CurveRed> Clist; Clist.push_back(C);
  return createfromcurves(s,Clist,nap);
}
void newforms::createfromcurves(int s, vector<CurveRed> Clist, int nap)
{
  if(verbose) cout << "In newforms::createfromcurves()..."<<endl;
  sign=s;
  int ncurves = Clist.size();
  if(ncurves==0) return;
  if(verbose) cout << "Making homspace..."<<flush;
  makeh1(sign);
  if(verbose) cout << "done." << endl;
  mvp=h1->maninvector(p0); 
  if(verbose) cout << "Making form_finder (nap="<<nap<<")..."<<flush;
  form_finder splitspace(this, (sign!=0), nap, 0, 1, 0, verbose);
  if(verbose) cout << "Recovering eigenspace bases with form_finder..."<<endl;
  // j1ds counts through the newforms as they are found
  basisflag=0; j1ds=0;
  vector< vector<long> > eigs(ncurves);
  int i,j;

  for(i=0; i<ncurves; i++) 
    eigs[i]=eiglist(Clist[i],nap);
  n1ds=0; nflist.resize(0);
  splitspace.recover(eigs);  // NB newforms::use() determines what is
			     // done with each one as it is found;
			     // this depends on basisflag and sign
  if(verbose) cout << "...done."<<endl;
}

// Read in newform data from old-style data files eigs/xN and intdata/[ep]N

void newforms::createfromolddata()
{
  long i, j, n = modulus;
  if(verbose) 
    cout << "Retrieving old-style newform data for N = " << n << endl;

  char* eigsname = new char[20];
  sprintf(eigsname,"eigs/x%ld",n);
  ifstream eigsfile(eigsname);

  if(!eigsfile.is_open())
    {
      cout<<"Unable to open file "<<eigsname<<" for eigs input"<<endl;
      delete[] eigsname;
      abort();
      return;
    }
  delete[] eigsname;

  short temp;
  eigsfile.read((char*)&temp,sizeof(short));   // # newforms
  n1ds=temp;
  eigsfile.read((char*)&temp,sizeof(short));   // # irrational newforms
  eigsfile.read((char*)&temp,sizeof(short));   // # ap
  nap=temp;
  if(n1ds==0) 
    {
      if(verbose) cout << "No newforms at level " << n << endl;
      eigsfile.close();
      return;
    }

  // read ap for each newform
  vector<long> * ap = new vector<long>[n1ds];
  for(i=0; i<n1ds; i++) ap[i].resize(nap);
  long ntotal = nap*n1ds;
  short* batch = new short[ntotal];
  eigsfile.read((char*)batch,ntotal*sizeof(short));
  eigsfile.close();
  short* batchptr = batch;
  for(j=0; j<nap; j++) 
    for(i=0; i<n1ds; i++) 
      ap[i][j]=*batchptr++;
  //  cout<<"Raw  ap:\n";  for(i=0; i<n1ds; i++) cout<<ap[i]<<endl;
  delete batch;

  // extract aq for each newform
  vector<long> * aq = new vector<long>[n1ds];
  for(i=0; i<n1ds; i++) aq[i].resize(npdivs);
  primevar pr; long q, k, a;
  for(k=0, j=0; j<plist.size(); j++)
    {
      q = plist[j];
      int q2divN = ::div(q*q,modulus);
      while((long)pr!=q) {pr++; k++;}
      for(i=0; i<n1ds; i++)
	{
	  a = ap[i][k];
	  aq[i][j] = a;
	  ap[i][k] = (q2divN? 0: -a);
	}
    }
  //  cout<<"Raw  aq:\n";  for(i=0; i<n1ds; i++) cout<<aq[i]<<endl;

  // read extra data for each newform
  vector<int> * data = new vector<int>[n1ds];
  for(i=0; i<n1ds; i++) data[i].resize(16);
  char* intdataname = new char[20];
  sprintf(intdataname,"intdata/e%ld",n);
  ifstream intdatafile(intdataname);
  if(!intdatafile.is_open())
    {
     intdataname[8] = 'p';
     intdatafile.clear();
     intdatafile.open(intdataname);
     if(!intdatafile.is_open())
       {
	 cout<<"Unable to open data file "<<intdataname<<" for data input"<<endl;
	 delete[] intdataname;
	 abort();
	 return;
       }
    }
  delete[] intdataname;
  
  long nloverp, dloverp, dp0, np0;
  for(i=0; i<n1ds; i++)
    {
      //      cout<<"Reading intdata for form #"<<(i+1)<<endl;
      intdatafile >> data[i][15]; // degphi
      //      cout<<"degphi = "<<data[i][15]<<endl;
      intdatafile >> data[i][0];  // sfe
      //      cout<<"sfe = "<<data[i][0]<<endl;
      intdatafile >> nloverp;     // num(L/P)
      intdatafile >> dloverp;     // den(L/P)
      intdatafile >> data[i][14]; // type
      //      cout<<"type = "<<data[i][14]<<endl;
      intdatafile >> dp0; data[i][3]=dp0;  // dp0
      //      cout<<"dp0 = "<<data[i]=dp0[3]<<endl;
      intdatafile >> np0; data[i][2]=np0;  // np0
      //      cout<<"np0 = "<<data[i][2]<<endl;
      if(dp0==0) data[i][3]=(2*nloverp*np0)/dloverp;
      data[i][1]=1+p0-data[i][2]; // ap0 (not in intdata file)
      //      cout<<"ap0 = "<<data[i][1]<<endl;
      intdatafile >> data[i][4];  // lplus
      intdatafile >> data[i][5];  // mplus
      intdatafile >> data[i][6];  // lminus
      intdatafile >> data[i][7];  // mminus
      intdatafile >> data[i][8];  // a
      intdatafile >> data[i][9];  // b
      intdatafile >> data[i][10]; // c
      intdatafile >> data[i][11]; // d
      intdatafile >> data[i][12]; // dotplus
      intdatafile >> data[i][13]; // dotminus
    }  
  intdatafile.close();
  //  cout<<"Raw  data:\n";  for(i=0; i<n1ds; i++) cout<<data[i]<<endl;


  // construct the newforms from this data
  nflist.reserve(n1ds);
  for(i=0; i<n1ds; i++) 
    nflist.push_back(newform(data[i],aq[i],ap[i],this));
  //  delete ap; delete aq; delete data;

  if(verbose) 
    {
      cout << "Finished reading oldstyle newform data for N = " << n << endl;
      display();
    }
}

// Construct bases (homology eigenvectors) from eigenvalue lists:
void newforms::makebases(int flag)
{
  if(n1ds==0) return;
  
  if(((sign==-1)||(dim(nflist[0].bplus)>0)) && 
     ((sign==+1)||(dim(nflist[0].bminus)>0))) 
    return;
  if(verbose) cout << "Making homspace..."<<flush;
  makeh1(sign);
  if(verbose) cout << "done." << endl;
  mvp=h1->maninvector(p0); 
  if(verbose) cout << "Making form_finder (nap="<<nap<<")..."<<flush;
  form_finder splitspace(this, (sign!=0), nap, 0, 1, 0, verbose);
  if(verbose) cout << "Recovering eigenspace bases with form_finder..."<<endl;
  // basisflag controls what ::use() does with the nfs when found
  // j1ds counts through the newforms as they are found
  basisflag=flag; j1ds=0;
  vector< vector<long> > eigs(n1ds);
  int i,j;

  unfix_eigs();
  if((n1ds>1)&&(modulus<130000)) // reorder into old order
    {
      if(verbose) cout<<"Reordering newforms before recovery as N<130000"<<endl;
      if(verbose>1) {cout<<"Before sorting:\n"; display();}
      sort();
      if(verbose>1) {cout<<"After sorting:\n"; display();}
    }
  for(i=0; i<n1ds; i++) 
    {
      eigs[i] = nflist[i].aplist;
      //      cout<<i+1<<": ";vec_out(cout,eigs[i],10); cout<<endl;
    }

  splitspace.recover(eigs);  // NB newforms::use() determines what is
			     // done with each one as it is found;
			     // this depends on basisflag and sign
  if(verbose) cout << "...done."<<endl;
  refix_eigs();
  if((n1ds>1)&&(modulus<130000)) // reorder into old order
    {
      if(verbose) cout<<"Reordering newforms after recovery as N<130000"<<endl;
      if(verbose>1) {cout<<"Before sorting:\n"; display();}
      sort(1);
      if(verbose>1) {cout<<"After sorting:\n"; display();}
    }
}

void newforms::merge()
{
  if(n1ds==0) return;
  if(verbose) cout << "Making homspace..."<<flush;
  makeh1(0);
  if(verbose) cout << "done." << endl;
  vec bplus, bminus;
  j1ds = 0; basisflag = 1;
  mvlplusvecs.clear();
  mvlminusvecs.clear();
  for(int inf=0; inf<n1ds; inf++)
   {
     if(verbose) cout << "Newform #"<<(inf+1)<<":"<<endl;
     if(verbose) cout<< "-about to extend bplus,bminus..."<<flush;
     bplus.init(h1->nsymb);
     bminus.init(h1->nsymb);
     int i,j;
     for(i=1; i<=h1->nsymb; i++)
       {      
	 j = h1plus->coordindex[i-1];
	 if (j==0) bplus[i] = 0;
	 else if (j>0) bplus[i] =  nflist[inf].coordsplus[j];
	 else if (j<0) bplus[i] = -nflist[inf].coordsplus[-j];
	 j = h1minus->coordindex[i-1];
	 if (j==0) bminus[i] = 0;
	 else if (j>0) bminus[i] =  nflist[inf].coordsminus[j];
	 else if (j<0) bminus[i] = -nflist[inf].coordsminus[-j];
       }
     if(verbose) cout<< "done, about to contract bplus,bminus..."<<flush;
     bplus = h1->contract_coords(bplus);
     bplus /= vecgcd(bplus);
     bminus = h1->contract_coords(bminus);
     bminus /= vecgcd(bminus);
     if(verbose) cout<< "done."<<endl;
     if(verbose>1)
       {
	 cout << " new bplus  = "<<bplus <<":"<<endl;
	 cout << " new bminus = "<<bminus<<":"<<endl;
       } 
     // These new dual eigenvectors are used to compute all
     // additional data needed for curve and modular symbol
     // computation (scaling and cuspidal factors and type)
     use(bplus, bminus, nflist[inf].aplist);
   }
  
}

vector<long> newforms::apvec(long p) //  computes a[p] for each newform
{
  //cout<<"In apvec with p = "<<p<<endl;
  vector<long> apv(n1ds);
  vec v;
  long i,j,iq,ap; 
  if(::div(p,modulus)) // we already have all the aq
    { 
      if(::div(p*p,modulus))
	for (i=0; i<n1ds; i++) apv[i] = 0;
      else
	{
	  iq = find(plist.begin(),plist.end(),p)-plist.begin();
	  for (i=0; i<n1ds; i++) apv[i] = -nflist[i].aqlist[iq];
	}
      return apv;
    }
  
  // now p is a good prime

  long maxap=(long)(2*sqrt((double)p)); // for validity check

  map<long,vec> images; // [j,v] stores image of j'th M-symbol in v
                        // (so we don't compute any more than once)
  vec bas, imagej;
  long fac;
  long p2=(p-1)>>1; // (p-1)/2
  long sl, sg, x1, x2, x3, y1, y2, y3, a, b, c, q, r;
  long u1,u2,u3;
  long ind;
  
  // Compute the image of the necessary M-symbols (hopefully only one)
  //cout<<"Computing images of M-symbols"<<endl<<flush;
  //cout<<"jlist = "<<jlist<<endl;
  for(std::set<long>::const_iterator jj=jlist.begin(); jj!=jlist.end(); jj++)
    {
      imagej=vec(n1ds); // initialised to 0
      j=*jj;
      //cout<<"Computing image of "<<j<<"'th M-symbol "<<endl;
      symb s = h1->symbol(h1->freegens[j-1]);
      //cout<<" = "<<s<<"..."<<flush;
      long u=s.cee(),v=s.dee(); 
      mat& pcd = h1->projcoord;
// Matrix [1,0;0,p]
      ind = h1->coordindex[h1->index2(u,p*v)];
      if(ind>0) imagej+=pcd.row(ind);
      else if(ind<0) imagej-=pcd.row(-ind);
// Matrix [p,0;0,1]
      ind = h1->coordindex[h1->index2(p*u,v)];
      if(ind>0) imagej+=pcd.row(ind);
      else if(ind<0) imagej-=pcd.row(-ind);
// Other matrices
      for(sg=0; sg<2; sg++) // signs
	for(r=1; r<=p2; r++)
	  {
	    a = -p; 
	    b = sg ? -r : r ;
	    u1=u*p; u2=v-u*b;
	    ind = h1->coordindex[h1->index2(u1,u2)];
	    if(ind>0) imagej+=pcd.row(ind);
	    else if(ind<0) imagej-=pcd.row(-ind);
	    while(b!=0)
	      {
		c=mod(a,b); q=(a-c)/b;
		if(q==1) {u3=  u2-u1;} else {u3=q*u2-u1;}
		a=-b; b=c; u1=u2; u2=u3;
		ind = h1->coordindex[h1->index2(u1,u2)];
		if(ind>0) imagej+=pcd.row(ind);
		else if(ind<0) imagej-=pcd.row(-ind);
	      }
	  }    
      images[j]=imagej/(h1->h1denom());
      //cout<<" image is "<<imagej<<endl;
    }

  for (i=0; i<n1ds; i++)
    {
// recover eigenvalue:
      apv[i]=ap=images[nflist[i].j0][i+1]/nflist[i].fac;
// check it is in range:
      if((ap>maxap)||(-ap>maxap))
	{
	  cout<<"Error:  eigenvalue "<<ap<<" for p="<<p
	      <<" for form # "<<(i+1)<<" is outside valid range "
	      <<-maxap<<"..."<<maxap<<endl;
	  abort();
	}
    }  
  return apv;
}

void newforms::addap(long last) // adds ap for primes up to the last'th prime
{
  if(n1ds==0) return;
  long i, j, p;

  j=0;
  if(verbose>1)  // output the ap already known...
    for(primevar pr(nflist[0].aplist.size()); pr.ok(); pr++, j++) 
      {
	p=(long)pr;
	if(ndiv(p,modulus)) cout<<"p="; else cout<<"q="; 
	cout<<p<<":\t";
	{
	  for (i=0; i<n1ds; i++) cout<<nflist[i].aplist[j]<<"\t";
	  cout<<endl;
	}
      }
  // Now compute and output the rest of the ap...
  for(primevar pr(last,1+nflist[0].aplist.size()); pr.ok(); pr++) 
    {
      p=(long)pr;
      vector<long> apv=apvec(p);
      if(verbose>1) 
	{
	  if(ndiv(p,modulus)) cout<<"p="; else cout<<"q="; 
	  cout<<p<<":\t";
	  for (i=0; i<n1ds; i++) cout<<apv[i]<<"\t";
	  cout<<endl;
	}
      for (long i=0; i<n1ds; i++) nflist[i].aplist.push_back(apv[i]);
    }
}

void output_to_file_no_newforms(long n, int binflag)
{
  char* name;
  if(binflag)
    name = nf_filename(n,'x');
  else 
    name = nf_filename(n,'e');
  ofstream out(name);
  delete[] name;
  if(binflag)
    {
      int a=0;
      out.write((char*)&a,sizeof(int));
      out.write((char*)&a,sizeof(int));
      out.write((char*)&a,sizeof(int));
    }
  else
    {
      out<<"0 0 0\n";
    }
  out.close();
   
}
char* nf_filename(long n, char c)
{
  char* nf_dir = getenv("NF_DIR"); 
  string nf_file;
  if (nf_dir==NULL) 
    nf_file = string("./newforms");
  else
    nf_file = string(nf_dir);
  char* filename=new char[20];
  sprintf(filename,"%s/%c%ld",nf_file.c_str(),c,n);
  return filename;
}

  // for the i'th newform return the value of the modular symbol {0,r}
rational newforms::plus_modular_symbol(const rational& r, long i) const
{
  return rational(h1->nfprojchain(num(r),den(r),nflist[i].coordsplus), 
		  nflist[i].cuspidalfactorplus);  
}

rational newforms::minus_modular_symbol(const rational& r, long i) const
{
  return rational(h1->nfprojchain(num(r),den(r),nflist[i].coordsminus), 
		  nflist[i].cuspidalfactorminus);  
}

pair<rational,rational> newforms::full_modular_symbol(const rational& r, long i) const
{
  mat m(h1->coord_vecs.size()-1,2);
  m.setcol(1,nflist[i].coordsplus);
  m.setcol(2,nflist[i].coordsminus);
  vec a = h1->projchain(num(r),den(r),m);
  rational a1(a[1],nflist[i].cuspidalfactorplus);
  rational a2(a[2],nflist[i].cuspidalfactorminus);
  return pair<rational,rational> ( a1, a2 );
}

  // Attempt to compute and display the elliptic curve for each
  // newform; return a list of newform indices where this failed.
vector<int> newforms::showcurves(vector<int> forms, int verbose)
{
  if((verbose>1)&&(sqfac>1)) cout<<"c4 factor " << sqfac << endl;

  bigfloat rperiod;
  vector<int> badcurves; // will hold the indices of forms for which we fail to find a curve
  vector<int>::const_iterator inf; // will iterate through the forms to be used

  for(inf=forms.begin(); inf!=forms.end(); inf++)
   {
     if(verbose) 
       cout<<"\n"<<"Form number "<<*inf+1<<"\n";
     else cout<<(*inf+1)<<" ";

     newform* nfi = &(nflist[*inf]);
     Curve C = getcurve(*inf,-1,rperiod,verbose);
     Curvedata CD(C,1);  // The 1 causes minimalization
     if(verbose) cout << "\nCurve = \t";
     cout << (Curve)CD << "\t";
     CurveRed CR(CD);
     cout << "N = " << getconductor(CR) << endl;
     if(verbose) cout<<endl;

     if(getconductor(CR)!=modulus)
       {
	 cout<<"No curve found"<<endl;
	 badcurves.push_back(*inf);
       }
   }

  return badcurves;
}
