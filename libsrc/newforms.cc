// FILE  NEWFORMS.CC: implementation of newforms class
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

#include <iomanip>
#include <eclib/newforms.h>
#include <eclib/periods.h>
#include <eclib/curvesort.h>

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

// Compare two ap-vectors lexicographically, using less_ap(.,.,old):
int less_apvec(const vector<long>& v, const vector<long>& w, int old=0);
int less_apvec(const vector<long>& v, const vector<long>& w, int old)
{
  auto wi=w.begin();
  for ( const auto& vi : v)
    {
      int s = less_ap(vi,*wi++,old);
      if(s) return s;
    }
  return 0;
}

// Old newform sorting comparison function: first by aq, then by ap,
// each lexicographically, with individual eigs compared first by
// absolute value, then positive first:

struct old_newform_comparer {
  bool operator()(const newform& f, const newform& g)
  {
    int s = less_apvec(f.aqlist,g.aqlist,1);
    if(s==0) s = less_apvec(f.aplist,g.aplist,1);
    return (s==1);
  }
}
  less_newform_old;

// New newform sorting comparison function: only ap,
// lexicographically, with individual eigs compared by
// numerical value:

struct new_newform_comparer {
  bool operator()(const newform& f, const newform& g)
  {
    return less_apvec(f.aplist,g.aplist,0)==1;
  }
}
  less_newform_new;

vector<long> eiglist(const newform& f, int oldorder)
{
  /*
  cout<<"Entering eiglist with f.aqlist="<<f.aqlist<<"\nand f.aplist=";
  vec_out(cout,f.aplist,10);
  cout<<endl;
  */
  vector<long> eigs;
  primevar pr;
  long N = (f.nf)->N;
  auto aqi=f.aqlist.begin();
  auto api=f.aplist.begin();
  if(oldorder)
    {
      eigs.resize(f.aplist.size());
      auto eigsi=eigs.begin();
      while(aqi!=f.aqlist.end())
	*eigsi++ = *aqi++;
      while(api!=f.aplist.end())
	{
	  if(ndivides(pr,N)) *eigsi++ = *api;
	  api++; pr++;
	}
    }
  else
    {
      eigs=f.aplist; // copy; now adjust the aq:
      auto eigsi=eigs.begin();
      while((aqi!=f.aqlist.end())&&(eigsi!=eigs.end()))
	{
	  if(::divides(pr.value(),N)) *eigsi = (*aqi++);
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
  rk=-1;
}

// Newform constructor, given the homology basis vector(s) and
// Hecke eigenvalues

newform::newform(const vec& vplus, const vec& vminus, const vector<long>& ap, newforms* nfs,long ind)
  :nf(nfs), sign(nfs->sign), bplus(vplus),bminus(vminus),index(ind),aplist(ap),rk(-1)
{
  int verbose=(nf->verbose);

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

  // Compute coordsplus/minus and denomplus/minus

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

  // Set default values for optimalityfactrplus/minus (will be reset if constructing from curves)
  optimalityfactorplus  = 1;
  optimalityfactorminus = 1;
}

int newform::check_expand_contract()
{
  int success=1;
  scalar denom = nf->h1->h1denom();
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
  scalar denom = nf->h1->h1denom();
  aqlist.resize(nf->npdivs);
  auto api=aplist.begin();
  auto pi=nf->plist.begin();
  auto aqi=aqlist.begin();
  primevar pr;   long q, i;
  long N = nf->N;
  while((api!=aplist.end())&&(aqi!=aqlist.end()))
    {
      q=pr.value(); pr++;
      if(::divides(q,N))
	{
	  *aqi++=*api;
	  *api=(::divides(q*q,N)? 0: -*api);
	  pi++;
	}
      api++;
    }
  if(aqi!=aqlist.end()) // compute missing aq
    {
      scalar piv;
      ssubspace espace(1, nf->modulus);
      if(sign==-1)
        espace=make1d(bminus,piv, nf->modulus);
      else
        espace=make1d(bplus,piv, nf->modulus);
      piv*=denom;
      while(aqi!=aqlist.end()) // compute missing aq
	{
	  q=*pi++;
	  if(nf->verbose) cout<<"Computing Wq for q="<<q<<"..."<<flush;
	  smat Wq = nf->h1->s_heckeop_restricted(q,espace,1,0);
	  long aq = I2long(Wq.elem(1,1) / piv);
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
  auto api=aplist.begin();
  auto aqi=aqlist.begin();
  primevar pr;
  long N = nf->N;
  while((api!=aplist.end())&&(aqi!=aqlist.end()))
    {
      if(::divides(pr.value(),N)) *api=*aqi++;
      api++;
      pr++;
    }
}

// After recovering eigenbases, we need to replace the ap for bad p
void newform::refix_eigs()
{
  auto api=aplist.begin();
  primevar pr;
  long N = nf->N, np = nf->npdivs, ip=0;
  while((api!=aplist.end())&&(ip<np))
    {
      long q=pr.value();
      if(::divides(q,N))
	{
	  *api=(::divides(q*q,N)? 0: -*api);
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
  auto api = aplist.begin();
  while(pr.value()!=nf->p0) {pr++; api++;}
  ap0=*api;
  np0 = 1 + (nf->p0) - ap0;
  if(nf->verbose) cout<<"ap0 = "<<ap0<<"\tnp0 = "<<np0<<endl;

  if(sign==-1) return;

  pdot = I2long((nf->mvp)*bplus); // should be negative since L(f,1)>=0
  if (pdot>0)
    // NB This will ensure that plus modular symbols have the right
    // sign for curves where L(E,1) is nonzero, but more work is
    // necessary for the plus symbols when L(Em1)=0, and for minus
    // symbols.  The additional work is done in find_matrix().
    {
      coordsplus = -coordsplus;
      bplus = -bplus;
      pdot  = -pdot;
    }
  dp0=abs(pdot);
  // DO NOT scale pdot by denom: factor will cancel when used to compute ap
  // DO scale dp0 since it is used to compute L/P
  if(dp0!=0)
    {
      if(denomplus>1)
	{
	  if(::divides(denomplus,dp0))  dp0/=denomplus;
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
  // ncoords is the same as ngens in homspace, i.e. the number of symbols aftre 2-term relations
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
        coordsplus[i]=dotmodp(cvi,bplus, nf->modulus);
      if(sign!=+1)
        coordsminus[i]=dotmodp(cvi,bminus, nf->modulus);
    }
  scalar contp = content(coordsplus);
  contplus = I2long(contp);
  if (contplus>1) coordsplus/=contp;
  scalar contm = content(coordsminus);
  contminus = I2long(contm);
  if (contminus>1) coordsminus/=contm;

  if(sign!=+1)
    {
      denomminus=contminus*cuspidalfactorminus;
      if(verbose>1) cout<<"coordsminus   = "<<coordsminus<<endl;
      if(verbose) cout<<"denomminus   = "<<denomminus<<endl;
    }
  if(sign!=-1)
    {
      denomplus=contplus*cuspidalfactorplus;
      if(verbose>1) cout<<"coordsplus   = "<<coordsplus<<endl;
      if(verbose) cout<<"denomplus   = "<<denomplus<<endl;
    }
}

void newform::find_cuspidal_factors()
{
  vec bplusc, bminusc;
  int verbose = nf->verbose;

  cuspidalfactorplus=1;
  cuspidalfactorminus=1;

  if(!(nf->h1->cuspidal))
    {
      if(sign!=-1) // do this if sign = 0,1
        {
          bplusc=(nf->h1->tkernbas)*bplus;
          scalar cfp = content(bplusc);
          cuspidalfactorplus = I2long(cfp);
          bplusc /= cfp;
        }
      if(sign!=+1) // do this if sign = 0,-1
	{
	  bminusc=(nf->h1->tkernbas)*bminus;
          scalar cfm = content(bminusc);
	  cuspidalfactorminus = I2long(cfm);
	  bminusc/= cfm;
        }
      if(sign==0)  // do this only if sign = 0
        {
	  type = (int)I2long((3-content(bplusc-bminusc)));
	  if(verbose) cout<<"Lattice type = "<<type<<endl;
          if((type!=1)&&(type!=2))
            {
              cerr<<"Error: lattice type computed to be "<<type<<", should be 1 or 2!"<<endl;
            }
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
    {
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
    }
  if(sign!=+1)
    {
      lminus=0;
      mminus=0;
    }
  if(nf->squarelevel) return;

  long N = nf->N;
  for (primevar lvar; lvar.ok() &&
	 (((sign!=-1)&&(mplus==0)) ||
	  ((sign!=+1)&&(mminus==0))); lvar++)
    {
      //cout << "Trying l = " << lvar << endl;
      while (N%lvar==0) {lvar++;}
      long l = lvar;
      //cout << "Trying l = " << l << endl;
      if (legendre(-N,l)!=sfe) continue;
      //cout << "Legendre condition passed... " << endl;

      if((sign!=-1)&&(mplus==0)&&(l%4==1))
	{
	  lplus = l;
	  //cout << "Trying lplus = " << l << "\n";
	  auto  vi = nf->mvlplusvecs.find(l);
	  if(vi==nf->mvlplusvecs.end())
	    mplus = I2long((nf->mvlplusvecs[l]=nf->h1->manintwist(l))*bplus);
	  else
	    mplus = I2long((vi->second)*bplus);
          // We force mplus>0 to fix the sign of the modular symbol to
          // agree with L(f*chi,1)>0, since L(f*chi,1) is real and a
          // positive multiple of mplus.  This uses the fact that the
          // Gauus sum is +sqrt(l).
          if (mplus<0)
            {
              mplus = -mplus;
              bplus = -bplus;
              coordsplus = -coordsplus;
            }
	  if((denomplus>1)&&(mplus!=0))
	    {
	      if(::divides(denomplus,mplus))  mplus/=denomplus;
	      else
		cout<<"Warning in newform constructor: mplus not divisible by denomplus!"
		    <<endl;
	    }
	}
      if((sign!=+1)&&(mminus==0)&&(l%4==3))
	{
	  lminus = l;
	  //cout << "Trying lminus = " << l << "\n";
	  auto vi = nf->mvlminusvecs.find(l);
	  if(vi==nf->mvlminusvecs.end())
	    mminus = I2long((nf->mvlminusvecs[l]=nf->h1->manintwist(l))*bminus);
	  else
	    mminus = I2long((vi->second)*bminus);
          // We force mminus<0 to fix the sign of the modular symbol
          // to agree with L(f*chi,1)>0, since L(f*chi,1) is real and
          // a negative multiple of mminus.  This uses the fact that
          // the Gauus sum is +i*sqrt(l).
          if (mminus>0)
            {
              mminus = -mminus;
              bminus = -bminus;
              coordsminus = -coordsminus;
            }
	  if((denomminus>1)&&(mminus!=0))
	    {
	      if(::divides(denomminus,mminus))  mminus/=denomminus;
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
  long N = nf->N;
  int found=0;
  vec v;
  for(d=2; !found; d++)
    {
      if(1==gcd(d,N))
        {
          for(b=1; (b<d) && !found; b++)
            {
              if(1==bezout(d,-N*b,a,c))
                {
		  //  cout<<"b/d = "<<b<<"/"<<d<<": ";
                  v = nf->h1->coords(b,d).as_vec();
		  //  cout<<"v="<<v<<endl;
                  if(sign!=-1)
                    {
                      dotplus = I2long(v*bplus);
                      if(::divides(denomplus,dotplus))
                        dotplus/=denomplus;
                      else
                        cout<<"Warning in find_matrix: dotplus not divisible by denomplus!"<<endl;
                    }
                  if(sign!=+1)
                    {
                      dotminus = I2long(v*bminus);
                      if(::divides(denomminus,dotminus))
                        dotminus/=denomminus;
                      else
                        cout<<"Warning in find_matrix: dotminus not divisible by denomminus!"<<endl;
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
  if((int)aplist.size()>=nap) return;
  int verbose=(nf->verbose);
  // Do not make the espace right away, as it is possible that the
  // only ap we are missing are aq which we already have...
  ssubspace espace(1, nf->modulus);
  int have_espace=0;
  scalar piv;

  primevar pr(nap,aplist.size()+1);
  while((int)aplist.size()<nap)
    {
      long p=pr, ap;
      if(::divides(p,nf->N))
	{
	  if(::divides(p*p,nf->N))
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
                espace=make1d(bminus,piv, nf->modulus);
              else
                espace=make1d(bplus,piv, nf->modulus);
	      piv*=nf->h1->h1denom();
	      have_espace=1;
	    }
	  ap = I2long((nf->h1->s_heckeop_restricted(p,espace,1,0)).elem(1,1) / piv);
	}
      aplist.push_back(ap);
      pr++;
    }
  if(verbose>1) cout<<"aplist = "<<aplist<<endl;
}

// Compute analytic rank and special value if not set
void newform::compute_rank()
{
  if (rk==-1) // not yet computed
    {
      ldash1 x(nf, this);
      Lvalue = abs(x.value()); // = L^{(r)}(f,1) -- note the r! factor!
      rk = 0;
      if (num(loverp)==0) // else trivially 0
        rk = x.rank();
    }
}

long newform::rank()
{
  compute_rank();  return rk;
}

bigfloat newform::special_value()
{
  compute_rank();  return Lvalue;
}

// To find optimality factors when created from a curve E.  Here
// CP_opt is the periods of this newform, and of the optimal curve E0,
// which we compute earlier in the newforms class since that is where
// getperiods() is implemented.

void newform::find_optimality_factors(const CurveRed& E, int i)
{
  int verbose=(nf->verbose);
  bigcomplex w1,w2;
  bigfloat x0, y0, xE, yE;

  // Definitions: for the period lattice of the optimal curve or
  // newform, x0 = (type)*(least real period) and y0 = (type)*(least
  // imag period)/i.  Here type = #components = 1 (Delta<0) or 2
  // (Deta>0).  Note that if we are in the plus or minus quotient then
  // we can still compute these even though we do not now the type
  // since they generate the projections of the period lattice to the
  // real (resp. imaginary) axis.

  // Similarly xE, yE for the input curve E, though here we do know
  // the type.

  // Compute the real and/or imginary periods of the newform, which
  // are those of the optimal curve in the isogeny class:
  if (sign==+1)
    {
      nf->get_real_period(i,x0);
      x0 = 2*abs(x0);
    }
  else if (sign==-1)
    {
      // NB it is impossible to get the scaling right in this case
      nf->get_imag_period(i,y0);
      y0 = abs(y0);
    }
  else
    {
      Cperiods CP_opt = nf->getperiods(i);
      int opt_type = CP_opt.getwRI(w1,w2);
      x0 = opt_type * abs(w1.real());
      y0 = abs(w2.imag());
    }

  // Compute the real and/or imginary periods of the input curve, which may not be optimal:
  Cperiods CP(E);
  int Etype = CP.getwRI(w1,w2);
  xE = Etype * abs(w1.real()); // least real period
  yE = abs(w2.imag()); // least imag period
  // now xE, yE are twice the least real/imag part of a period of E in both cases

  // Now we find rational approximations to the rations x0/xE and
  // y0/yE.  These will have very small numerators and denominators.
  long n,d;
  if (sign!=-1)
    {
      ratapprox(x0/xE,n,d,163);
      optimalityfactorplus = rational(n,d);
      if (verbose) cout << "x-ratio (optimalityfactorplus) = " << (x0/xE) << " = " <<n<<"/"<<d<<" = "<<optimalityfactorplus << endl;
    }
  if (sign!=+1)
    {
      ratapprox(y0/yE,n,d,163);
      optimalityfactorminus = rational(n,d);
      if (verbose) cout << "y-ratio (optimalityfactorminus) = " << (y0/yE) << " = " <<n<<"/"<<d<<" = "<< optimalityfactorminus << endl;
    }
}

// Adjust the sign of dotplus/dotminus, mplus/mminus, pdot and the
// associated primitive coordinate vectors by computing a period
// numerically.  Here we only need sufficient precision to determine
// the sign, knowing the values to be nonzero rationals with small
// denominator.

void newform::sign_normalize()
{
  int verbose = (nf->verbose);
  int sign = (nf->sign);

  periods_direct integrator(nf, this);
  integrator.compute();
  bigfloat x0 = integrator.rper();
  bigfloat y0 = integrator.iper();
  if(verbose>1)
    cout<<"integral over {0,"<<b<<"/"<<d<<"} = ("<<x0<<")+i*("<<y0<<")"<<endl;
  if (sign!=-1)
    {
      if (dotplus*x0<0)
        {
          if (verbose)
            cout<<"flipping sign for plus symbols"<<endl;
          coordsplus = -coordsplus;
          bplus = -bplus;
          dotplus = -dotplus;
          pdot = -pdot;
          mplus = -mplus;
        }
    }
  if (sign!=+1)
    {
      if (dotminus*y0<0)
        {
          if (verbose)
            cout<<"flipping sign for minus symbols"<<endl;
          coordsminus = -coordsminus;
          dotminus = -dotminus;
          bminus = -bminus;
          mminus = -mminus;
        }
    }
  if (verbose>1)
    {
      cout<<"Approximate numerical values:"<<endl;
      if (sign==0)
        {
          cout<<"x = "<<(x0/dotplus)<<endl;
          cout<<"y = "<<(y0/dotminus)<<endl;
          cout<<"integral over {0,"<<b<<"/"<<d<<"} = ("<<x0<<")+i*("<<y0<<")"<<endl;
        }
      if (sign==+1)
        {
          cout<<"x = "<<(x0/dotplus)<<endl;
          cout<<"integral over {0,"<<b<<"/"<<d<<"} has real      part "<<x0<<endl;
        }
      if (sign==-1)
        {
          cout<<"y = "<<(y0/dotminus)<<endl;
          cout<<"integral over {0,"<<b<<"/"<<d<<"} has imaginary part "<<x0<<endl;
        }
    }
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
	  h1plus = new homspace(N, modulus, 1 /*plusflag*/,  0 /*cuspidal*/,  0 /*verbose*/);
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
	  h1minus = new homspace(N,modulus, -1, /*plusflag*/ 0 /*cuspidal*/, 0 /*verbose*/);
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
	  h1full = new homspace(N,modulus, 0 /*plusflag*/, 0 /*cuspidal*/, 0 /*verbose*/);
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
  of = new oldforms(ntp,h1, modulus, (verbose>1),sign); // h1 provides the level*
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
       form_finder ff(this, modulus, (sign!=0),maxdepth,mindepth,1,0,verbose);
       basisflag=0;
       ff.find();
    }
  if(verbose)
    {
      cout << "Total dimension = " << h1->dimension << endl;
      cout << "Number of rational newforms = " << n1ds <<endl;
      if(h1->dimension==of->totalolddim+n1ds)
	cout<<"The whole space splits over Q" << endl;
    }

  if(n1ds==0) return;

  if((n1ds>1)&&(N<130000)) // reorder into old order
    {
      if(verbose) cout<<"Reordering newforms into old order as N<130000"<<endl;
      // if(verbose) cout<<"Before sorting:\n"; display();
      sort(1);
      // if(verbose) cout<<"After sorting:\n"; display();
    }
  // At this point the newforms may contain different numbers of ap,
  // so we need to even these up, which we do by computing more ap for
  // those which need it.
  if(n1ds>0)
    {
      int nap, maxnap=0;
      for(int i=0; i<n1ds; i++)
	if((nap=nflist[i].aplist.size())>maxnap)
          maxnap=nap;
      if(verbose)
	cout<<"Max number of ap in newforms so far is "<<maxnap
	    <<", increasing to " << DEFAULT_SMALL_NAP << endl;
      if (maxnap < DEFAULT_SMALL_NAP)
        maxnap = DEFAULT_SMALL_NAP;
      for(int i=0; i<n1ds; i++)
        {
          if((nap=nflist[i].aplist.size())<maxnap)
            {
              if(verbose)
                cout<<"Newform #"<<(i+1)<<" has only "<<nap
                    <<" ap so we need to compute more..."<<endl;
              nflist[i].add_more_ap(maxnap);
            }

          // Now if necessary we adjust the sign of dotplus/dotminus and the
          // associated primitive coordinate vectors by computing a period
          // numerically.  Here we only need sufficient precision to determine
          // the sign, knowing the values to be nonzero rationals with small
          // denominator.

          if(verbose)
            cout<<"Newform #"<<(i+1)<<": fixing sign normalization using approximate periods"<<endl;
          nflist[i].sign_normalize();
        }
    }
   // Compute homspace::projcoord, so proj_coords can be used
   // Replaces coord_vecs of homspace with projections onto eigenspaces
   // NB if #newforms>1 this MUST be re-called after any sorting of newforms
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
          const vec& bas = (sign==-1? nflist[i].bminus: nflist[i].bplus);
          nflist[i].fac = bas[j0];
          if (verbose>1)
            {
              cout<<"Newform #"<<(i+1)<<": bplus = "<<bas<<endl;
              cout<<"   fac = "<<nflist[i].fac<<endl;
            }
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
	  const vec& bas = nflist[i].bplus;
	  j=1; while(bas[j]==0) j++;
	  jlist.insert(j);
	  nflist[i].j0 = j;
	  nflist[i].fac = nflist[i].bplus[j];
          if (verbose>1)
            {
              cout<<"Newform #"<<(i+1)<<": bplus = "<<bas<<endl;
              cout<<"   fac = "<<nflist[i].fac<<endl;
            }
	}
      if(verbose)  cout<<"jlist="<<jlist<<endl;
    }
}

// Compute homspace::projcoord, so proj_coords can be used
// Replaces coord_vecs of homspace with projections onto eigenspaces
// NB if #newforms>1 this MUST be re-called after any sorting of newforms
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
      int n = nf_subset[j1ds++];
      newform& nf = nflist[n];
      if(verbose)
	cout<<"Filling in data for for newform #"<<(n+1)<<": bases..."<<flush;
      nf.sign=sign;
      if(sign==+1)
        nf.bplus=b1;
      if(sign==-1)
        nf.bminus=b1; // form_finder puts the basis vector in b1
      if(sign==0)
	{
	  nf.bplus=b1;
	  nf.bminus=b2;
	}
      if(verbose)
	cout<<"type and cuspidal factors..."<<flush;
      nf.find_cuspidal_factors();
      if(verbose)
	cout<<"coords..."<<flush;
      nf.find_coords_plus_minus();
      if(sign==0)
	{
	  if(verbose)
	    cout<<"twisting primes..."<<flush;
	  nf.find_twisting_primes();
	  if(verbose)
	    cout<<"matrix..."<<flush;
	  nf.find_matrix();
	}
      if(verbose)
        {
          cout<<"done."<<endl;
          cout<<"Finished filling in data for newform #"<<(n+1)<<endl;
        }
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
    ::sort(nflist.begin(),nflist.end(),less_newform_old);
  else
    ::sort(nflist.begin(),nflist.end(),less_newform_new);
}

void newforms::sort_into_Cremona_label_order()
{
  int sort_system = level_range(N);
  switch (sort_system)
    {
    case 0: case 1:
      sort(1);
      break;
    case 2: default:
      sort(0);
      break;
    case 3:
      unfix_eigs();
      sort(0);
      refix_eigs();
    }
  if (sort_system>0)
    return;
  // apply the permutation defined in curvesort.cc
  vector<newform> sorted_nflist(n1ds);
  for (int i=0; i<n1ds; i++)
    sorted_nflist[i] = nflist[booknumber0(N,i)];
  nflist = sorted_nflist;
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
 cout << "\n"<<n1ds<<" newform(s) at level " << N << ":" << endl;
 cout<<"p0="<<p0<<endl;
 // if(dim(mvp)!=0) cout<<"mvp="<<mvp<<endl;
 cout<<"#ap=\t"<<nflist[0].aplist.size()<<endl;
 long i;
 for(i=0; i<n1ds; i++)
   {cout<<i+1<<":\t";
    nflist[i].display();
  }
}

void newforms::display_modular_symbol_map(int check) const
{
 rational rplus, rminus;
 vector<bigfloat> x0s;
 vector<bigfloat> y0s;
 bigfloat x0,y0;
 if (check)
   for(long k=0; k<n1ds; k++)
     {
       cout<<"getting period(s) for newform # "<<(k+1)<<endl;
       get_both_periods(k,x0,y0);
       cout<<"x0="<<x0<<", y0="<<y0<<endl;
       x0s.push_back(x0);
       y0s.push_back(y0);
     }

 for(long i=0; i<h1->nsymb; i++)
   {
     symb s = h1->symbol(i);
     modsym ms = modsym(s);
     cout<<s<<" = "<<ms<<" -> ";
     rational alpha = ms.alpha(), beta = ms.beta();
     long a=0, b=0, c=0, d=0, g=0;
     // intialise all, otherwise the compiler worries
     if (num(alpha)==0)
       {
         b = num(beta);
         d = den(beta);
         g = bezout(-N*b,d,c,a); // so g=a*d-b*N*c
       }
     long j=h1->coordindex[i];
     long sg=::sign(j); j=abs(j);
     // cout<<"g="<<g<<", j="<<j<<"("<<sg<<")"<<endl;
     if(j==0)
       for(long k=0; k<n1ds; k++)
	 if(sign!=0)
	   cout<<"0 ";
	 else
	   cout<<"(0,0) ";
     else
       for(long k=0; k<n1ds; k++)
         {
           if (check && (g==1))
             {
               // NB g==1 implies that a,b,c,d are defined and define
               // a matrix in SL(2,Z)
               periods_direct integrator(this,&(nflist[k]));
               integrator.compute(a,b,c,d);
               x0 = integrator.rper();
               y0 = integrator.iper();
             }
           if(sign!=-1)
             {
               long nrplus = I2long(sg*nflist[k].coordsplus[j]);
               long drplus = nflist[k].cuspidalfactorplus;
               rplus = rational(nrplus,drplus);
               if (check && (g==1))
                 {
                   long n = I2long(Iround(drplus *x0/x0s[k])); // should = nrplus
                   if (n != nrplus)
                     {
                       cout << "plus check fails: rplus = "<<nrplus<<"/"<<drplus<<" = "<<rplus<<endl;
                       cout << "real part = "<<x0<<endl;
                       cout << "x0        = "<<x0s[k]<<endl;
                       cout << "ratio     = "<<x0/x0s[k]<<endl;
                       cout << "scaled ratio = "<<n<<endl;
                     }
                 }
               rplus *= nflist[k].optimalityfactorplus;
             }
           if(sign!=+1)
             {
               long nrminus = I2long(sg*nflist[k].coordsminus[j]);
               long drminus = nflist[k].cuspidalfactorminus;
               rminus = rational(nrminus,drminus);
               if (check && (g==1))
                 {
                   long n = I2long(Iround(drminus *y0/y0s[k])); // should = nrminus
                   if (n != nrminus)
                     {
                       cout << "minus check fails: rminus = "<<rminus<<endl;
                       cout << "imag part = "<<y0<<endl;
                       cout << "y0        = "<<y0s[k]<<endl;
                       cout << "ratio     = "<<y0/y0s[k]<<endl;
                       cout << "scaled ratio = "<<n<<endl;
                     }
                 }
               rminus *= nflist[k].optimalityfactorminus;
             }
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

void newforms::output_to_file(int binflag, int smallflag) const
{
  long i,j;
  char prefix = 'e';
  if(binflag) prefix = 'x';
  string name = smallflag ? small_nf_filename(N,prefix)
                          : nf_filename(N,prefix);
  ofstream out(name.c_str());
  if(!out)
    {
      cerr<<"Unable to open file "<<name<<" for newform output"<<endl;
      return;
    }
  // else
  //   {
  //     cout<<"--outputting newforms data to "<<name<<" (smallflag="<<smallflag<<")"<<endl;
  //     display();
  //   }

  if(n1ds==0)
    {
      putout(out,(int)0,binflag);
      putout(out,(int)0,binflag);
      putout(out,(int)0,binflag);
      out.close();
      return;
    }

  // Line 1:  #newforms, #aq, #ap
  int nap = nflist[0].aplist.size();
  if (smallflag)
    {
      if (nap>=DEFAULT_SMALL_NAP) nap=DEFAULT_SMALL_NAP;
      else
	{
	  cout<<"Warning: small newforms output will only have" << nap
	      << "a_p (at least " << DEFAULT_SMALL_NAP <<"required" << endl;
	}
    }
  putout(out,(int)n1ds,binflag);
  putout(out,(int)nflist[0].aqlist.size(),binflag);
  putout(out, nap,binflag);
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
  for(j=0; j<int(nflist[0].aqlist.size()); j++)
    {
      for(i=0; i<n1ds; i++) putout(out,(short)nflist[i].aqlist[j],binflag);
      nl(out,binflag);
    }
  nl(out,binflag);

  // Lines (21+#aq)-(20+#aq+#ap):  ap for each newform
  for(j=0; j<nap; j++)
    {
      for(i=0; i<n1ds; i++) putout(out,(short)nflist[i].aplist[j],binflag);
      nl(out,binflag);
    }
  out.close();
}

// Read in newform data from file NF_DIR/xN

void newforms::createfromdata(int s, long ntp,
			      int create_from_scratch_if_absent,
			      int small_data_ok)
{
  sign = s;
  long i, j;
  if(verbose) cout << "Retrieving newform data for N = " << N << endl;

  string name = nf_filename(N,'x');
  ifstream datafile(name.c_str());
  if(small_data_ok && !datafile.is_open())
    {
      if (verbose)
	cout << "No file "<<name<<" exists, trying ";
      name = small_nf_filename(N,'x');
      if (verbose)
	cout << name << "instead..."<<endl;
      datafile.open(name.c_str());
    }
  if(!datafile.is_open())
    {
      if(verbose) cout<<"Unable to open file "<<name<<" for oldform input"<<endl;
      if(create_from_scratch_if_absent)
	{
	  if(verbose) cout<<"Creating from scratch instead"<<endl;
	  createfromscratch(sign, ntp);
	  output_to_file(1,0); // full newform data
	  output_to_file(1,1); // short newform data (only DEFAULT_SMALL_NAP ap)
	  if(verbose)
            {
              cout << "Finished creating newform data for N = " << N << endl;
              display();
            }
	}
      return;
    }

  int temp_int;
  datafile.read((char*)&temp_int,sizeof(int));   // = number of newforms
  n1ds=temp_int;
  datafile.read((char*)&temp_int,sizeof(int));   // = number of bad primes
  datafile.read((char*)&temp_int,sizeof(int));   // = number of eigs
  nap=temp_int;
  if(n1ds==0)
    {
      if(verbose) cout << "No newforms at level " << N << endl;
      datafile.close();
      return;
    }

  vector<vector<int>> data(n1ds, vector<int>(16));
  vector<vector<long>> aq(n1ds, vector<long>(npdivs));
  vector<vector<long>> ap(n1ds, vector<long>(nap));

  // read extra data for each newform (16 ints each)
  long ntotal = 16*n1ds;
  int* batch_i = new int[ntotal];
  datafile.read((char*)batch_i,ntotal*sizeof(int));
  int *batch_i_ptr = batch_i;
  for(j=0; j<16; j++)
    for(i=0; i<n1ds; i++)
      data[i][j]=*batch_i_ptr++;
  delete[] batch_i;
  //  cout<<"Raw  data:\n";  for(i=0; i<n1ds; i++) cout<<data[i]<<endl;

  // read aq for each newform (npdivs longs each)
  ntotal = npdivs*n1ds;
  short* batch = new short[ntotal];
  datafile.read((char*)batch,ntotal*sizeof(short));
  short *batchptr = batch;
  for(j=0; j<npdivs; j++)
    for(i=0; i<n1ds; i++)
      aq[i][j]=*batchptr++;
  //  cout<<"Raw  aq:\n";  for(i=0; i<n1ds; i++) cout<<aq[i]<<endl;

  // read ap for each newform (nap longs each)
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

  if(verbose)
    {
      cout << "Finished reading newform data for N = " << N << endl;
      if(verbose>1) display();
    }
}


// Create from a list of Hecke eigenvalues from an elliptic curve

vector<long> eiglist(CurveRed& C, int nap)
{
  long N = I2long(getconductor(C));
  vector<long> ans;
  for(primevar pr(nap); pr.ok(); pr++)
    {
      long p=pr; bigint pp(p);
      if(N%p==0)
	ans.push_back(LocalRootNumber(C,pp));
      else
	ans.push_back(C.ap(p));
    }
  return ans;
}

// extract the eigenvalues for bad primes
vector<long> aqlist(vector<long> aplist, long N)
{
  long iq=0, naq = pdivs(N).size();
  //cout << "Setting aq of size "<<naq<<endl;
  auto api = aplist.begin();
  vector<long> aq(naq);
  for(primevar pr; (iq<naq)&&pr.ok(); pr++)
    {
      long p=pr;
      if (N%p==0)
	{
	  //cout << "Setting aq["<<p<<"] = "<<*api<<endl;
	  aq[iq++] = *api;
	}
      api++;
    }
  return aq;
}

// Create from a list of elliptic curves of the right conductor:

void newforms::createfromcurve(int s, const CurveRed& C, int nap)
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
  long nap_default = long(100*sqrt(N));
  if (nap<nap_default)
    {
      if(verbose)
        {
          cout << "--increasing nap from "<< nap << " to " << nap_default << endl;
        }
      nap=nap_default;
    }
  if(verbose) cout << "Making form_finder (nap="<<nap<<")..."<<flush;
  form_finder splitspace(this, modulus, (sign!=0), nap, 0, 1, 0, verbose);
  if(verbose) cout << "Recovering eigenspace bases with form_finder..."<<endl;
  // j1ds counts through the newforms as they are found
  basisflag=0; j1ds=0;
  vector< vector<long> > eigs(ncurves);
  int i;

  for(i=0; i<ncurves; i++)
    eigs[i]=eiglist(Clist[i],nap);
  n1ds=0; nflist.resize(0);
  splitspace.recover(eigs);  // NB newforms::use() determines what is
			     // done with each one as it is found;
			     // this depends on basisflag and sign

  // Get period lattice of each newform (and hence of optimal curves)
  for(i=0; i<ncurves; i++)
    {
      if(verbose)
        cout<<"Finding optimality scaling factors using approximate periods"<<endl;
      nflist[i].find_optimality_factors(Clist[i], i);
      if(verbose)
        cout<<"Fixing sign normalization using approximate periods"<<endl;
      nflist[i].sign_normalize();
    }
  if(verbose) cout << "...done."<<endl;
}

void newforms::createfromcurves_mini(vector<CurveRed> Clist, int nap)
{
  if(verbose) cout << "In newforms::createfromcurves_mini()..."<<endl;
  n1ds = Clist.size();
  nflist.reserve(n1ds);

  if (n1ds>0) // construct the ap and aq vectors from the curves
    {
      long N = I2long(getconductor(Clist[0]));
      for(long i=0; i<n1ds; i++)
	{
	  vector<long> ap=eiglist(Clist[i],nap);
	  vector<long> aq=aqlist(ap,N);
	  // dummy data -- these fields are not set and will not be used
	  vector<int> data(16,0);
	  newform nf(data,aq,ap,this);
	  if (verbose)
	    {
	      cout<<"adding this newform: "<<endl;
	      nf.display();
	    }
	  nflist.push_back(newform(data,aq,ap,this));
	}
    }
  if(verbose) cout << "...done."<<endl;
}

// Read in newform data from old-style data files eigs/xN and intdata/[ep]N

void newforms::createfromolddata()
{
  long i, j;
  if(verbose)
    cout << "Retrieving old-style newform data for N = " << N << endl;

  stringstream eigsname;
  eigsname << "eigs/x" << N;
  ifstream eigsfile(eigsname.str().c_str());

  if(!eigsfile.is_open())
    {
      cout<<"Unable to open file "<<eigsname.str()<<" for eigs input"<<endl;
      return;
    }

  short temp;
  eigsfile.read((char*)&temp,sizeof(short));   // # newforms
  n1ds=temp;
  eigsfile.read((char*)&temp,sizeof(short));   // # irrational newforms
  eigsfile.read((char*)&temp,sizeof(short));   // # ap
  nap=temp;
  if(n1ds==0)
    {
      if(verbose) cout << "No newforms at level " << N << endl;
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
  delete[] batch;

  // extract aq for each newform
  vector<long> * aq = new vector<long>[n1ds];
  for(i=0; i<n1ds; i++) aq[i].resize(npdivs);
  primevar pr; long k, a;
  for(k=0, j=0; j<int(plist.size()); j++)
    {
      long q = plist[j];
      int q2divN = ::divides(q*q,N);
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
  stringstream intdataname;
  intdataname << "intdata/e" << N;
  ifstream intdatafile(intdataname.str().c_str());
  if(!intdatafile.is_open())
    {
      cout<<"Unable to open data file "<<intdataname.str()<<" for data input"<<endl;
      return;
    }

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
      cout << "Finished reading oldstyle newform data for N = " << N << endl;
      display();
    }
}

// Construct bases (homology eigenvectors) from eigenvalue lists:
void newforms::makebases(int flag, int all_nf)
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
  form_finder splitspace(this, modulus, (sign!=0), nap, 0, 1, 0, verbose);
  if(verbose) cout << "Recovering eigenspace bases with form_finder..."<<endl;
  // basisflag controls what ::use() does with the nfs when found
  // j1ds counts through the newforms as they are found
  basisflag=flag;
  j1ds = 0; // counts through newforms as they are recovered

  if (all_nf)
    {
      nf_subset.resize(n1ds);
      std::iota(nf_subset.begin(), nf_subset.end(), 0);
    }

  unfix_eigs();
  sort();
  vector< vector<long> > eigs(nf_subset.size());
  std::transform(nf_subset.begin(), nf_subset.end(), eigs.begin(),
                 [this](const int& nfi) {return nflist[nfi].aplist;});

  splitspace.recover(eigs);  // NB newforms::use() determines what is
			     // done with each one as it is found;
			     // this depends on basisflag and sign
  if(verbose)
    cout << "...done."<<endl;
  refix_eigs();
  if(verbose>1)
    {
      cout<<"Reordering newforms after recovery"<<endl;
      cout<<"Before sorting:\n";
      display();
    }
  sort(int(N<130000)); // old order for N<130000, else new order
  if(verbose>1)
    {
      cout<<"After sorting:\n";
      display();
    }
}

void newforms::merge(int all_nf)
{
  if(n1ds==0) return;
  if(verbose) cout << "Making homspace..."<<flush;
  makeh1(0);
  if(verbose) cout << "done." << endl;
  vec bplus, bminus;
  j1ds = 0;
  basisflag = 1;
  mvlplusvecs.clear();
  mvlminusvecs.clear();
  if (verbose>1)
    cout<<"merging newforms " << nf_subset << endl;
  unfix_eigs();
  sort();
  for (const auto nfi : nf_subset)
   {
     int inf = nf_subset[nfi];
     if(verbose)
       {
         cout << "Newform #"<<(inf+1)<<":"<<endl;
         cout << "-about to extend bplus,bminus..."<<flush;
       }
     bplus.init(h1->nsymb);
     bminus.init(h1->nsymb);
     for(int i=1; i<=h1->nsymb; i++)
       {
	 int j = h1plus->coordindex[i-1];
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
     bplus /= content(bplus);
     bminus = h1->contract_coords(bminus);
     bminus /= content(bminus);
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
  refix_eigs();
  sort(int(N<130000)); // old order for N<130000, else new order
}


void update(const mat& pcd, vec& imagej, long ind)
{
  if(ind>0) imagej.add_row(pcd,ind);
  else
    if(ind<0) imagej.sub_row(pcd,-ind);
}

vector<long> newforms::apvec(long p) //  computes a[p] for each newform
{
  //cout<<"In apvec with p = "<<p<<endl;
  vector<long> apv(n1ds);
  if(::divides(p,N)) // we already have all the aq
    {
      if(::divides(p*p,N))
	for (long i=0; i<n1ds; i++)
          apv[i] = 0;
      else
	{
	  long iq = find(plist.begin(),plist.end(),p)-plist.begin();
	  for (long i=0; i<n1ds; i++)
            apv[i] = -nflist[i].aqlist[iq];
	}
      return apv;
    }

  // now p is a good prime

  long maxap=(long)(2*sqrt((double)p)); // for validity check

  map<long,vec> images; // [j,v] stores image of j'th M-symbol in v
                        // (so we don't compute any more than once)
  vec bas, imagej;
  long p2=(p-1)>>1; // (p-1)/2
  long sg, a, b, c, q, r;
  long u1,u2,u3;

  // Compute the image of the necessary M-symbols (hopefully only one)
  //cout<<"Computing images of M-symbols"<<endl<<flush;
  //cout<<"jlist = "<<jlist<<endl;

  for( const auto& j : jlist)
    {
      imagej=vec(n1ds); // initialised to 0
      symb s = h1->symbol(h1->freegens[j-1]);
      //cout<<"Computing image of "<<j<<"'th M-symbol "<<s<<endl;
      //cout<<" = "<<s<<"..."<<flush;
      long u=s.cee(),v=s.dee();
      const mat& pcd = h1->projcoord;
      //cout<<"projcoord = "<<pcd;
// Matrix [1,0;0,p]
      long ind = h1->coordindex[h1->index2(u,p*v)];
      update(pcd,imagej,ind);
      //cout<<"(1) (u1,u2)=("<<u<<","<<p*v<<") partial image index is "<<ind<<", subtotal="<<imagej<<endl;
// Matrix [p,0;0,1]
      ind = h1->coordindex[h1->index2(p*u,v)];
      update(pcd,imagej,ind);
      //cout<<"(2) (u1,u2)=("<<p*u<<","<<v<<") partial image index is "<<ind<<", subtotal="<<imagej<<endl;
// Other matrices
      for(sg=0; sg<2; sg++) // signs
	for(r=1; r<=p2; r++)
	  {
	    a = -p;
	    b = sg ? -r : r ;
	    u1=u*p; u2=v-u*b;
	    ind = h1->coordindex[h1->index2(u1,u2)];
            update(pcd,imagej,ind);
            //cout<<"(3) (u1,u2)=("<<u1<<","<<u2<<") partial image index is "<<ind<<", subtotal="<<imagej<<endl;
	    while(b!=0)
	      {
		c=mod(a,b); q=(a-c)/b;
		if(q==1) {u3=  u2-u1;} else {u3=q*u2-u1;}
		a=-b; b=c; u1=u2; u2=u3;
		ind = h1->coordindex[h1->index2(u1,u2)];
                update(pcd,imagej,ind);
                //cout<<"(4) (u1,u2)=("<<u1<<","<<u2<<") partial image index is "<<ind<<", subtotal="<<imagej<<endl;
	      }
	  }
      images[j]=imagej;
      //cout<<" image is "<<imagej<<endl;
    }

  for (long i=0; i<n1ds; i++)
    {
// recover eigenvalue:
      //cout<<"Numerator =   "<< images[nflist[i].j0][i+1] <<endl;
      //cout<<"Denominator = "<< nflist[i].fac << endl;
      long ap = I2long(images[nflist[i].j0][i+1]/nflist[i].fac);
      ap *= (sign==-1? nflist[i].contminus: nflist[i].contplus);
      ap /= I2long(h1->h1denom());
      apv[i]=ap;
      // check it is in range:
      if((ap>maxap)||(-ap>maxap))
	{
	  cout<<"Error:  eigenvalue "<<ap<<" for p="<<p
	      <<" for form # "<<(i+1)<<" is outside valid range "
	      <<-maxap<<"..."<<maxap<<endl;
          break; // no point in trying to compute any more.
	}
    }
  return apv;
}

void newforms::addap(long last) // adds ap for primes up to the last'th prime
{
  if(n1ds==0) return;
  if(verbose>1)  // output the ap already known...
    {
      long j=0;
      for(primevar pr(nflist[0].aplist.size()); pr.ok(); pr++, j++)
        {
          long p = pr;
          if(ndivides(p,N)) cout<<"p="; else cout<<"q=";
          cout<<p<<":\t";
          {
            for (long i=0; i<n1ds; i++) cout<<nflist[i].aplist[j]<<"\t";
            cout<<endl;
          }
        }
    }
  // Now compute and output the rest of the ap...
  for(primevar pr(last,1+nflist[0].aplist.size()); pr.ok(); pr++)
    {
      long p = pr;
      vector<long> apv=apvec(p);
      if(verbose>1)
	{
	  if(ndivides(p,N)) cout<<"p="; else cout<<"q=";
	  cout<<p<<":\t";
	  for (long i=0; i<n1ds; i++) cout<<apv[i]<<"\t";
	  cout<<endl;
	}
      for (long i=0; i<n1ds; i++) nflist[i].aplist.push_back(apv[i]);
    }
}

void output_to_file_no_newforms(long n, int binflag, int smallflag)
{
  char prefix = 'e';
  if(binflag)  prefix = 'x';
  string name = smallflag ? small_nf_filename(n,prefix)
                          : nf_filename(n,prefix);
  ofstream out(name.c_str());
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

// for the i'th newform return the value of the modular symbol {0,r} (default) or {oo,r}
rational newforms::plus_modular_symbol(const rational& r, long i, int base_at_infinity) const
{
  rational a(I2long(h1->nfproj_coords(num(r),den(r),nflist[i].coordsplus)),
		  nflist[i].cuspidalfactorplus);
  // {oo,r} = {0,r}+{oo,0} and loverp={oo,0} (not {0,oo}!)
  if (base_at_infinity) a+=nflist[i].loverp;
  a *= nflist[i].optimalityfactorplus;
  return a;
}

rational newforms::minus_modular_symbol(const rational& r, long i, int base_at_infinity) const
{
  // Ignore the value of base_at_infinity as it does not affect the minus symbol
  rational a(I2long(h1->nfproj_coords(num(r),den(r),nflist[i].coordsminus)),
		  nflist[i].cuspidalfactorminus);
  a *= nflist[i].optimalityfactorminus;
  return a;
}

pair<rational,rational> newforms::full_modular_symbol(const rational& r, long i, int base_at_infinity) const
{
  mat m(h1->coord_vecs.size()-1,2);
  m.setcol(1,nflist[i].coordsplus);
  m.setcol(2,nflist[i].coordsminus);
  vec a = h1->proj_coords(num(r),den(r),m);
  rational a1(I2long(a[1]),nflist[i].cuspidalfactorplus);
  // {oo,r} = {0,r}+{oo,0} and loverp={oo,0} (not {0,oo}!)
  if (base_at_infinity) a1 += nflist[i].loverp;
  a1 *= nflist[i].optimalityfactorplus;
  rational a2(I2long(a[2]),nflist[i].cuspidalfactorminus);
  a2 *= nflist[i].optimalityfactorminus;
  return pair<rational,rational> ( a1, a2 );
}

  // Attempt to compute and display the elliptic curve for each
  // newform; return a list of newform indices where this failed.
vector<int> newforms::showcurves(vector<int> forms, int verbose, string filename)
{
  if((verbose>1)&&(sqfac>1)) cout<<"c4 factor " << sqfac << endl;

  ofstream curve_out;
  int output_curves = (filename!="no");
  if (output_curves) curve_out.open(filename.c_str());
  bigfloat rperiod;
  bigint a1,a2,a3,a4,a6, NC;
  vector<int> badcurves; // will hold the indices of forms for which we fail to find a curve

  for( const auto& inf : forms)
   {
     if(verbose)
       cout<<"\n"<<"Form number "<<inf+1<<"\n";
     else cout<<(inf+1)<<" ";

     if (output_curves)
       curve_out << N << " "<< codeletter(inf) << " 1 ";

     Curve C = getcurve(inf,-1,rperiod,verbose);
     Curvedata CD(C,1);  // The 1 causes minimalization
     if(verbose) cout << "\nCurve = \t";
     cout << (Curve)CD << "\t";
     CurveRed CR(CD);
     NC = getconductor(CR);
     cout << "N = " << NC << endl;
     if(verbose) cout<<endl;

     if(NC!=N)
       {
	 cout<<"No curve found"<<endl;
	 badcurves.push_back(inf);
         if (output_curves)
           curve_out<<endl;
       }
     else
       if (output_curves)
         {
           C.getai(a1,a2,a3,a4,a6);
           curve_out<<"["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"]";
           int nt = CD.get_ntorsion();
           int r = nflist[inf].rank(); // analytic rank
           curve_out<<" "<<r<<" "<<nt<<" 0"<<endl;
         }
   }
  if (output_curves)
    curve_out.close();

  return badcurves;
}
