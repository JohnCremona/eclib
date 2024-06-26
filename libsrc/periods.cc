// FILE PERIODS.CC : implementation of classes for integrating newforms
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
//
#include <eclib/periods.h>
#ifndef MPFP // Multi-Precision Floating Point
#include <eclib/fixc6.h>
#endif

#define CHECK_PERIODS // check that curves constructed from periods
                      // have the same periods...

//#define TRACE_CACHE
//#define TRACE_USE

/////////////////////////////////////
//  functions for character class  //
/////////////////////////////////////

character::character(long m)
{
  modul=m;
  chartable.resize(m);
  init();
}

void character::init()
{
  if ( modul==1 ) chartable[0] = 1;
  else { long i=modul; while(i--) chartable[i] =  legendre(i,modul);}
}

void character::reset(long m)
{
  modul=m;
  chartable.resize(m);
  init();
}

/////////////////////////////////
//  functions for summer class //
/////////////////////////////////

#ifndef NEW_OP_ORDER
vector<long> resort_aplist(const level* iN,
			   const vector<long>& primelist,
			   const vector<long>& apl)
{
  long N = iN->modulus;
//Sort out the aplist, since apl has the W_q eigs first:
  long nap=apl.size();
  //  cout<<"resort_aplist(), using "<<nap<<" primes"<<endl;
  vector<long> aplist;
  aplist.reserve(nap);
  long ip = iN->npdivs;
  for(long i=0; i<nap; i++)
    { long p = primelist[i];
    auto pi = find(iN->plist.begin(),iN->plist.end(),p);
    if(pi==iN->plist.end()) // then p is good
      {
	aplist.push_back( apl[ip++]);
// 	cout << "i = "<<i<<",\tp = " << p << "\ta_p = " << aplist[i] << endl;
      }
    else // p is bad
      {
	if(::divides(p*p,N))
	  {
	    aplist.push_back(0);
// 	    cout << "i = "<<i<<",\tp = " << p << "\ta_p = " << aplist[i] << endl;
	  }
	else
	  {
	    long j = pi-(iN->plist.begin());  // p is j'th bad prime
	    aplist.push_back(- apl[j]);
// 	    cout << "i = "<<i<<",\tp = " << p << "\ta_p = " << aplist[i] << endl;
	  }
      }
    }
//   cout << "After rearranging, aplist = " << aplist << endl;
  return aplist;
}
#endif

void summer::initaplist(const level* iN, const vector<long>& apl)
{
  N = iN->modulus;
  primelist = primes(apl.size());   //First #apl primes, indexed from 0
#ifdef NEW_OP_ORDER
  aplist=apl;
#else
  aplist=resort_aplist(iN,primelist,apl);
#endif
}


void summer::use1(long n, long an)
{
 bigfloat can(to_bigfloat(-an));
 can /=  to_bigfloat(n);
 if (n<rootlimit) {an_cache[n]=an;
#ifdef TRACE_CACHE
 cout<<"Caching an["<<n<<"]="<<an<<endl;
#endif
 }
#ifdef TRACE_USE
 cout<<"use1(("<<n<<","<<an<<")"<<endl;
#endif
 if (n<limit1) { sum1 += func1(n) * can;}
}

void summer::use2(long n, long an)
{
 bigfloat can(to_bigfloat(-an)); can /=  to_bigfloat(n);
#ifdef TRACE_USE
 cout<<"use2("<<n<<","<<an<<")"<<endl;
 cout<<"(n,an,can)=("<<n<<","<<an<<","<<can<<"):\t";
#endif
  if (n<rootlimit) {an_cache[n]=an;
#ifdef TRACE_CACHE
  cout<<"Caching an["<<n<<"]="<<an<<endl;
#endif
  }
  if (n<limit1) { bigfloat term=func1(n) * can; sum1 += term;
#ifdef TRACE_USE
  cout<<"term="<<term<<",\tsum1="<<sum1<<"\n";
#endif
  }
  if (n<limit2) { bigfloat term = func2(n) * can; sum2+=term;
#ifdef TRACE_USE
  cout<<term<<",\t"<<sum2<<"\n";
#endif
  }
}

// This calls use(m,am) for all m=n*2^a*3^b*5^c*7^d<limit
void summer::use2357(long n, long an)
{
  long m2=n;
  for(long i2=0; (i2<=n2p)&&(m2<limit); i2++, m2*=2)
    {
      long am2=an*a2p_cache[i2];
      long m23=m2;
      for(long i3=0; (i3<=n3p)&&(m23<limit); i3++, m23*=3)
	{
	  long am23=am2*a3p_cache[i3];
	  long m235=m23;
	  for(long i5=0; (i5<=n5p)&&(m235<limit); i5++, m235*=5)
	    {
	      long am235=am23*a5p_cache[i5];
	      long m2357=m235;
	      for(long i7=0; (i7<=n7p)&&(m2357<limit); i7++, m2357*=7)
		{
		  long am2357=am235*a7p_cache[i7];
		  use(m2357,am2357);
		}
	    }
	}
    }
}

void summer::add(long n, long pindex, long y, long z)
{
  // y = a[n], z=a[n/p]
  // where p is the pindex'th prime, which is the largest prime dividing n

  //    cout<<"In add() with n="<<n<<",\tpindex="<<pindex<<",\ty="<<y
  //        <<",\tz="<<z<<endl;

 long x,p,pn,ip,istart=pindex;
  int triv = (y==0);
  if ( ! triv ) { use(n,y); istart=0; }

  ip=istart;
  p=primelist[ip];
  if(triv&&(p>rootlimit)) return;

  for(pn=n*p; (ip<=pindex) && (pn<=limit); ip++)
    {
      p=primelist[ip];
      pn=p*n;
      if(pn<=limit)
	{
	  x = y * aplist[ip];
	  if ( (ip==pindex)  && (N%p)) { x -=  p*z; }   // x = a[pn]
	  add(pn,ip,x,y);
	}
    }
}

void summer::add2357(long n, long pindex, long y, long z)
{
  // y = a[n], z=a[n/p]
  // where p is the pindex'th prime, which is the largest prime dividing n

#ifdef TRACE_USE
  cout<<"In add() with n="<<n<<",\tpindex="<<pindex<<",\ty="<<y
      <<",\tz="<<z<<endl;
#endif

 long x,p,pn,ip,istart=pindex;
  int triv = (y==0);
  if ( ! triv ) { use2357(n,y); istart=4; } // skip 2,3,5,7

  ip=istart;
  p=primelist[ip];
  if(triv&&(p>rootlimit)) return;

  for(pn=n*p; (ip<=pindex) && (pn<=limit); ip++)
    {
      p=primelist[ip];
      pn=p*n;
      if(pn<=limit)
	{
	  x = y * aplist[ip];
	  if ( (ip==pindex)  && (N%p)) { x -=  p*z; }   // x = a[pn]
	  add2357(pn,ip,x,y);
	}
    }
}

void summer::sumit()
{
  static double log2 = log(2.0);
  static double log3 = log(3.0);
  static double log5 = log(5.0);
  static double log7 = log(7.0);
  double loglimit = log((double)limit);
#ifdef TRACE_USE
  cout<<"In sumit(), limit="<<limit<<", rootlimit="<<rootlimit<<endl;
#endif
  sum1=sum2=to_bigfloat(0);
  long j;
  an_cache[1]=1;

  n2p=(long)(floor(loglimit/log2));
  a2p_cache.resize(n2p+1);
  a2p_cache[0]=1;
  a2p_cache[1]=aplist[0];
  for(j=2; j<=n2p; j++)
    {
      a2p_cache[j]=a2p_cache[j-1]*aplist[0];
      if(ndivides(2,N)) a2p_cache[j]-=2*a2p_cache[j-2];
    }
#ifdef TRACE_CACHE
  cout<<"a2p_cache = "<<a2p_cache<<endl;
#endif

  n3p=(long)(floor(loglimit/log3));
  a3p_cache.resize(n3p+1);
  a3p_cache[0]=1;
  a3p_cache[1]=aplist[1];
  for(j=2; j<=n3p; j++)
    {
      a3p_cache[j]=a3p_cache[j-1]*aplist[1];
      if(ndivides(3,N)) a3p_cache[j]-=3*a3p_cache[j-2];
    }
#ifdef TRACE_CACHE
  cout<<"a3p_cache = "<<a3p_cache<<endl;
#endif

  n5p=(long)(floor(loglimit/log5));
  a5p_cache.resize(n5p+1);
  a5p_cache[0]=1;
  a5p_cache[1]=aplist[2];
  for(j=2; j<=n5p; j++)
    {
      a5p_cache[j]=a5p_cache[j-1]*aplist[2];
      if(ndivides(5,N)) a5p_cache[j]-=5*a5p_cache[j-2];
    }
#ifdef TRACE_CACHE
  cout<<"a5p_cache = "<<a5p_cache<<endl;
#endif

  n7p=(long)(floor(loglimit/log7));
  a7p_cache.resize(n7p+1);
  a7p_cache[0]=1;
  a7p_cache[1]=aplist[3];
  for(j=2; j<=n7p; j++)
    {
      a7p_cache[j]=a7p_cache[j-1]*aplist[3];
      if(ndivides(7,N)) a7p_cache[j]-=7*a7p_cache[j-2];
    }
#ifdef TRACE_CACHE
  cout<<"a7p_cache = "<<a7p_cache<<endl;
#endif

  use2357(1,1);
  unsigned long ip=4; long p=11,ap;
#ifdef TRACE_USE
  cout<<"Using primes with indices "<<ip<<" to "<<primelist.size()
      <<" and less than "<<rootlimit<<endl;
#endif
  while( (ip+1<primelist.size()) && (p<rootlimit)  )
  {
    ap=aplist[ip];
#ifdef TRACE_USE
    cout<<"Using prime("<<ip<<") = "<<p<<", ap = "<<ap<<endl;
#endif
    add2357(p,ip,ap,1);
    p=primelist[++ip];
  }
#ifdef TRACE_USE
  cout<<endl;
#endif
#ifdef TRACE_USE
  cout<<"Using primes with indices "<<ip<<" to "<<primelist.size()
      <<" and less than "<<limit<<endl;
#endif
  while( (ip+1<primelist.size()) && (p<=limit)  )
  {
    ap=aplist[ip];
#ifdef TRACE_USE
    cout<<"Using prime("<<ip<<") = "<<p<<", ap = "<<ap<<endl;
#endif
    if(ap!=0)
      {
	use(p,ap);
	long m=2, n=2*p;
	while(n<limit)
	  {
#ifdef TRACE_CACHE
	    cout<<"Getting cached an["<<m<<"]="<<an_cache[m]<<endl;
#endif
	    use(n,ap*an_cache[m]);
	    m++; n+=p;
	  }
      }
    p=primelist[++ip];
  }
}

/////////////////////////////////////////////
//  functions for periods_via_lfchi class  //
/////////////////////////////////////////////

periods_via_lfchi::periods_via_lfchi (const level* iN, const newform* f)
  : chi1(f->lplus), chi2(f->lminus)
{
  type = f->type; dp0 = f->dp0; mplus = f->mplus; mminus = f->mminus;
  if(f->lplus==1) mplus=f->np0; // L/P = dp0/np0 so now P=L*mplus/dp0
  initaplist(iN,f->aplist);
  rootmod =  sqrt(to_bigfloat(N));
  factor1 =  exp(-(TWOPI)/ (to_bigfloat(chi1.modulus()) * rootmod));
  factor2 =  exp(-(TWOPI)/ (to_bigfloat(chi2.modulus()) * rootmod));
  long bp = bit_precision();
 // MUST keep brackets around TWOPI!:
  bigfloat rootlimit1=(bp-log(1-factor1))*rootmod/(TWOPI) ;
  bigfloat rootlimit2=(bp-log(1-factor1))*rootmod/(TWOPI) ;
  Iasb(limit1,rootlimit1);
  Iasb(limit2,rootlimit2);
#ifdef TRACE_USE
  cout << "Bit precision = "<<bp<<endl;
  cout << "Basic limits on n in sums = " << limit << endl;
#endif
  limit1=chi1.modulus()*limit1;
  limit2=chi2.modulus()*limit2;
  limit = ( limit2 > limit1 ? limit2 : limit1);
#ifdef TRACE_USE
  //  cout<<"periods_via_lfchi.limit = "<<limit<<endl;
  //  cout<<"Largest prime = "<<primelist[primelist.size()-1]<<endl;
#endif
/*
  if(primelist[primelist.size()-1]<limit)
    cout<<"\nLARGEST PRIME "<<primelist[primelist.size()-1]<<" IS BELOW LIMIT "<< limit <<" required for bit precision "<<bp<<endl;
*/
  rootlimit=sqrt(to_bigfloat(limit));
  an_cache.resize(I2long(Ifloor(rootlimit+1)),0);
#ifdef TRACE_USE
  cout << "factor1 = "<<factor1<<", factor2 = "<<factor2<<endl;
  cout << "Limits on n in sums = " << limit1 <<", "<<limit2 << endl;
#endif
#ifdef TRACE_CACHE
  cout << "Initial an_cache = "<<an_cache<<endl;
#endif
}

void periods_via_lfchi::compute(void)
{
  sumit();

  long q = chi1.modulus(), l = chi2.modulus();
  if (q==1) {rp = 2*sum1*mplus/to_bigfloat(dp0);}
  else      {rp = 2*sum1*sqrt(to_bigfloat(q))/to_bigfloat(mplus);}
  ip =            2*sum2*sqrt(to_bigfloat(l))/to_bigfloat(mminus);
}

//////////////////////////////////////////
//  functions for periods_direct class  //
//////////////////////////////////////////

periods_direct::periods_direct(const level* iN, const newform*f)
{
  eps_N = -(f->sfe);
  initaplist(iN, f->aplist);
  factor1 = -(TWOPI)  / sqrt(to_bigfloat(N));
  type=f->type;
  a=f->a,b=f->b,c=f->c,d=f->d;
}

void periods_direct::compute(long ta, long tb, long tc, long td)
{
  a=ta; b=tb; c=tc; d=td;
  compute();
}

void periods_direct::compute(void)
{
  if(d==0)
   {
    cout<<"Problem: cannot compute periods for matrix with d=0!"<<endl;
    rp=ip=0;
    return;
   }

  if (d<0) { a=-a;b=-b;c=-c;d=-d;}

  bigfloat drecip =  to_bigfloat(1) / to_bigfloat(d);
  if (int(ctab.size())!=d) // else same d as before, no need to recompute
    {
      int j; bigfloat x;
      ctab.clear();
      stab.clear();
      for(j=0; j<d; j++)
        {
          x = (TWOPI)*to_bigfloat(j)*drecip;
          ctab.push_back(cos(x));
          stab.push_back(sin(x));
        }
    }
  // cout<<"d = "<<d<<endl;
  // cout<<"ctab = "<<ctab<<endl;
  // cout<<"stab = "<<stab<<endl;
  theta1 = to_bigfloat(b) * drecip;
  theta2 = to_bigfloat(c) * drecip;
  b = posmod(b,d);
  c = posmod(c,d);
  factor2 = factor1 * drecip;
  long bp = bit_precision();
  Iasb(limit1,(-bp-log((1-exp(factor1))/3))/(factor1)) ;
  Iasb(limit2,(-bp-log((1-exp(factor2))/3))/(factor2)) ;

  limit = limit2;
  rootlimit=sqrt(to_bigfloat(limit));
  an_cache.resize(I2long(Ifloor(rootlimit+1)),0);
#ifdef TRACE_USE
  cout << "Limits on n in sums = " << limit1 <<", "<<limit2 << endl;
  cout<<"periods_direct.limit = "<<limit<<endl;
  cout<<"Largest prime = "<<primelist[primelist.size()-1]<<endl;
#endif
/*
  if(primelist[primelist.size()-1]<limit)
    cout<<"\nLARGEST PRIME "<<primelist[primelist.size()-1]<<" IS BELOW LIMIT "<< limit <<" required for bit precision "<<bp<<endl;
*/
#ifdef TRACE_CACHE
  cout << "Initial an_cache = "<<an_cache<<endl;
  cout << "Calling sumit()" << endl;
#endif

  sumit();

#ifdef TRACE_CACHE
  cout << "Returned from sumit()" << endl;
  cout<<"rp = "<<rp<<endl;
  cout<<"ip = "<<ip<<endl;
#endif
  rp = sum1;
  ip = sum2;
}

void periods_direct::use(long n, long an)
// specially defined, as does not fit simply into prototype provided by
// summer class
{
  if(n>limit) return;
#ifdef TRACE_USE
  cout<<"periods_direct::use("<<n<<","<<an<<")"<<endl;
#endif
  if (n<rootlimit) {an_cache[n]=an;
#ifdef TRACE_CACHE
  cout<<"Caching an["<<n<<"]="<<an<<endl;
#endif
  }
  bigfloat dn = to_bigfloat(n);
  bigfloat coeff = to_bigfloat(an)/dn;
  bigfloat ef2 = coeff * exp(dn*factor2);
  int nbd = (n*b)%d, ncd = (n*c)%d;
  if(eps_N==-1)
    {
      if(n<limit1)
	{
	  bigfloat ef1 = coeff * exp(dn*factor1);
#ifdef TRACE_USE
	  cout<<"n="<<n<<", ef1="<<ef1<<", ef2="<<ef2<<"\n";
#endif
	  sum1 -= (2*ef1);
	}
      sum1 += ef2*(ctab[nbd]+ctab[ncd]);
      sum2 += ef2*(stab[nbd]+stab[ncd]);
    }
  else
    {
      sum1 += ef2*(ctab[nbd]-ctab[ncd]);
      sum2 += ef2*(stab[nbd]-stab[ncd]);
    }
#ifdef TRACE_USE
  cout<<"done: sum1 = "<<sum1<<", sum2 = "<<sum2<<endl;
#endif
}

///////////////////////////////////////
//  functions for part_period class  //
///////////////////////////////////////

part_period::part_period(const level* iN, const newform*f)
{
  initaplist(iN, f->aplist);
}

void part_period::compute(const bigcomplex& z0)
{
  x0=(TWOPI)*real(z0);
  y0=(TWOPI)*imag(z0);
  compute();
}

void part_period::compute()
{
  long bp = bit_precision();
  Iasb(limit,bp/y0);
  limit1=limit2=limit;
  rootlimit=sqrt(to_bigfloat(limit));
  an_cache.resize(I2long(Ifloor(rootlimit+1)),0);
#ifdef TRACE_USE
  cout<<"Using limit = "<<limit<<endl;
#endif
#ifdef TRACE_CACHE
  cout << "Initial an_cache = "<<an_cache<<endl;
#endif
  sumit();

  rp = sum1;
  ip = sum2;
}


//////////////////////////////////
//  functions for ldash1 class  //
//////////////////////////////////

ldash1::ldash1(const level* iN, const newform* f)
{
  init(iN, f->aplist, f->sfe, f->loverp);
}

ldash1::ldash1(const  newforms* nf, long i)
{
  const newform* nfi = &((nf->nflist)[i]);
  init(nf, nfi->aplist, nfi->sfe, nfi->loverp);
}

void ldash1::init(const level* iN, const vector<long>& f_aplist, long f_sfe, const rational& f_loverp)
{
  initaplist(iN, f_aplist);

  rootmod=sqrt(to_bigfloat(N));
  factor1 = (TWOPI)/rootmod;
  long maxp = prime_number(aplist.size());
  limit  =  I2long(Ifloor((30+bit_precision())/factor1));
  //cout<<"ldash1::init() with N="<<N<<", bit_precision="<<bit_precision()<<endl;
  //cout<<"number of terms to use = "<<limit<<endl;
  if(limit>maxp) limit=maxp;
  limit1 = limit;
  rootlimit=sqrt(to_bigfloat(limit));
  an_cache.resize(I2long(Ifloor(rootlimit+1)),0);
#ifdef TRACE_CACHE
  cout << "Initial an_cache = "<<an_cache<<endl;
#endif
//decide on r (first estimate)
  computed=0;
  r=0;

  if(f_sfe==-1) { r=1;}
  else if(f_loverp==0) {r=2; }
}

void ldash1::compute(void)
{
  static const bigfloat two=to_bigfloat(2);
  if(computed) return;
  sumit(); ld1=two*sum1; computed=1;
  if(r==0) return;

  while(abs(ld1)<0.0001) // ?? What's a sensible value??
    {
      //      cout<<"L^(r)(1) small for r="<<r<<", increasing to ";
      r+=2;
      //      cout<<r<<endl;
      sumit(); ld1=two*sum1;
    }
}

bigfloat ldash1::func1(long n)
{
#ifdef MPFP
  long l = bit_precision();
  set_precision(l+20);
#endif
  bigfloat ans = -G(factor1*to_bigfloat(n));
#ifdef MPFP
  set_precision(l);
#endif
  return ans;
}

/////////////////////////////////
//  functions for lfchi class  //
/////////////////////////////////

lfchi::lfchi (const level* iN, const newform* f)
{
  initaplist(iN, f->aplist);
  rootmod =  sqrt(to_bigfloat(N));
  long bp = bit_precision();
  Iasb(limit0,bp*rootmod/(TWOPI)) ; // brackets essential
}

void lfchi::compute(long ell)
{
  chi.reset(ell);
  limit = limit1 = ell*limit0;
  factor1 =  exp(-(TWOPI)/ (to_bigfloat(ell) * rootmod));
  rootlimit=sqrt(to_bigfloat(limit));
  an_cache.resize(I2long(Ifloor(rootlimit+1)),0);
#ifdef TRACE_USE
  cout<<"lfchi.factor1 = "<<factor1<<" for ell="<<ell<<endl;
  cout<<"lfchi.limit = "<<limit<<" for ell="<<ell<<endl;
#endif
#ifdef TRACE_CACHE
  cout << "Initial an_cache = "<<an_cache<<endl;
#endif
  sumit();

  val = -2*sum1; // - since sumit() uses -a_n
}

////////////////////////////////////////////////////////

int newforms::get_real_period(long i, bigfloat& x, int verbose) const
{
  const newform* nfi = &(nflist[i]);
  int rank0 = (num(nfi->loverp)!=0);
  lfchi lx(this,nfi);

  if(rank0)
    {
      if(verbose)
	cout << "Computing real period via L(f,1): ";
      lx.compute(1);
      if(verbose) cout<<"L(f,1) = "<<lx.value()<<"; ";
      rational lop=nfi->loverp;
      x = abs(lx.value())*to_bigfloat(den(lop))/to_bigfloat(num(lop));
      if(verbose) cout<<"real period = "<<x<<endl;
      return 1;
    }

  long lplus=nfi->lplus;
  long mplus=nfi->mplus;
  if(mplus!=0)
    {
      if(verbose)
	cout << "Computing real period via L(f,chi,1) with chi mod "
	     <<lplus<<"...";
      lx.compute(lplus);
      if(verbose) cout<<"L(f,chi,1) = "<<lx.scaled_value()<<"; ";
      x = abs(lx.scaled_value()/to_bigfloat(mplus));
      if(verbose) cout<<"real period = "<<x<<endl;
      return 1;
    }

  // we only reach here if sfe=-1 and level is square

  periods_direct pd(this,nfi);
  if(verbose) cout<<"...computing directly...";
  pd.compute();
  x = pd.rper();
  long dotplus = (nfi->dotplus);
  if (!dotplus) return 0;
  x /= dotplus;
  if(verbose) cout<<"real period (after scaling by "<<dotplus<<") = "<<x<<endl;
  return 1;
}

int newforms::get_imag_period(long i, bigfloat& y, int verbose) const
{
  const newform* nfi = &(nflist[i]);
  lfchi lx(this,nfi);

  long lminus=(nfi->lminus);
  long mminus=(nfi->mminus);
  if(mminus!=0)
    {
      if(verbose)
	cout << "Computing imaginary period via L(f,chi,1) with chi mod "
	     <<lminus<<"...";
      lx.compute(lminus);
      if(verbose) cout<<"L(f,chi,1) = "<<(lx.scaled_value())<<"; ";
      y = (lx.scaled_value()/to_bigfloat(mminus));
      if(verbose) cout<<"imaginary period = "<<y<<endl;
      return 1;
    }
  return 0;
}

////////////////////////////////////////////////////////

Cperiods newforms::getperiods(long i, int method, int verbose)
  // method=1 to use periods_direct
  // method=0 to use periods_via_lfchi
  // method=-1 (default) to use whichever is best
{
  newform* nfi = &(nflist[i]);
  if(method==-1) // find and use best method
    {
      if((this->squarelevel) || ((nfi->lplus)==0) || ((nfi->lminus)==0))
	method=1;
      else
	{
	  long d = (nfi->d);
	  if(d>0)
	    method = ( d < (nfi->lplus)) || ( d < (nfi->lminus));
	  else method=0;
	}
    }
  //  method=1;
  if(method==1)
       {
         if(verbose)
	   {
	     cout<<"Finding periods -- direct method "<<endl;
	     cout << "using matrix ("<<nfi->a<<","<<nfi->b<<";"<<nfi->c
		  <<","<<nfi->d<<"), dotplus="<<nfi->dotplus
		  <<", dotminus="<<nfi->dotminus<<"; type="<<nfi->type<<endl;
	   }
         periods_direct per(this,nfi);
         per.compute();
         return Cperiods(per.rper()/(nfi->dotplus),
                         per.iper()/(nfi->dotminus),
                         (nfi->type));;
       }
     else
       {
         if(verbose)
	   cout<<"Finding periods -- via L(f_chi) using twists by "
	       <<nfi->lplus<<" and "<<nfi->lminus<<endl;
         periods_via_lfchi per(this,nfi);
         per.compute();
         return per.getperiods();
       }
}

Curve newforms::getcurve(long i, int method, bigfloat& rperiod, int verbose)
// Cannot just use trans_to_curve() as we need to fixc6...
// which DOES use n and i
{
  long n = modulus;  // NOT redundant, used in fixc6
  long fac = sqfac;
  long fac6=(odd(fac)?fac:2*fac);
//  if(fac>1) cout<<"factor = "<<fac<<endl;
//
// p^2 | n => p | c4, c6 so we take advantage of the fact that c4, c6
// must be multiples of fac when rounding.
//
  bigcomplex wR, wRI, w1, w2, c4, c6;
  Cperiods cp = getperiods(i, method, verbose);
  if(verbose)
    cout<<cp<<endl;
  cp.getwRI(wR, wRI);
  rperiod = real(wR);
  cp.getwi(w1, w2);
  getc4c6(w2,w1,c4,c6);  // from compproc.h
  bigfloat rc4 = real(c4), rc6 = real(c6);
  if(verbose)
    cout << "c4 = " << rc4 << "\nc6 = " << rc6 << endl;
  bigint ic4 = fac*Iround(rc4/fac);
  bigint ic6 = fac6*Iround(rc6/fac6);
  if(verbose)
    {
      cout << "After rounding";
      if(fac>1)
        cout << ", using factors " << fac << " for c4 and " << fac6 << " for c6";
      cout<<":\n";
      cout << "ic4 = " << ic4 << "\nic6 = " << ic6 << endl;
    }
  // To fix the c4 or c6 values insert data in files fixc4.data and
  // fixc6.data; NB the index here (i) starts at 0, but class fixc6
  // adjusts so in the data files, start at 1

#ifndef MPFP // Multi-Precision Floating Point
  c4c6fixer(n,i,ic4,ic6);
  if(verbose)
    {
      cout << "After fixing: \n";
      cout << "ic4 = " << ic4 << "\nic6 = " << ic6 << endl;
    }
#endif
  Curve C(ic4,ic6);
  if(C.isnull()) return C;
  Curvedata CD(C,1);  // The 1 causes minimalization
  CurveRed CR(CD);
  if (getconductor(CR) != n)
    {
      if (verbose)
        cout << "Constructed curve "<<C<<" has wrong conductor "<<getconductor(CR)<<endl;
      C = Curve(); // reset to null curve
      return C;
    }
#ifdef CHECK_PERIODS
// Check periods were correct:
//  verbose=1;
  Cperiods cpC(CD);
  bigcomplex wRC, wRIC;
  cpC.getwRI(wRC, wRIC);
  wR=abs(wR); wRI=abs(imag(wRI));
  wRC=abs(wRC); wRIC=abs(imag(wRIC));
  //cout<<"C = "<<CD<<endl;
  //  cout<<"type(C)="<<get_lattice_type(cpC)<<endl;
  if((get_lattice_type(cp)!=get_lattice_type(cpC)))
    {
      cout<<"Period lattice type of constructed curve does not match that"
	  <<" of the newform"<<endl;
      cout<<"Lattice type of C: "<<get_lattice_type(cpC)<<endl;
      cout<<"Lattice type of f: "<<get_lattice_type(cp)<<endl;
    }
  else if(verbose) cout<<"Lattice type checks OK"<<endl;
  if((abs((wR-wRC)/wRC)>0.0001))
    {
      cout<<"Real period of constructed curve does not match that"
	  <<" of the newform (using bit precision "<<bit_precision()<<")"<<endl;
      cout<<"Real period of C: "<<real(wRC)<<endl;
      cout<<"Real period of f: "<<real(wR)<<endl;
      cout<<"Ratio = "<<real(wR)/real(wRC)<<endl;
    }
  else if(verbose) cout<<"Real period checks OK"<<endl;
  if((abs((wRI-wRIC)/wRIC)>0.0001))
    {
      cout<<"Second period of constructed curve does not match that"
	  <<" of the newform (using bit precision "<< bit_precision()<<")"<<endl;
      cout<<"Imag part of second period of C: "<<real(wRIC)<<endl;
      cout<<"Imag part of second period of f: "<<real(wRI)<<endl;
      cout<<"Ratio of imaginary parts = "<<real(wRI)/real(wRIC)<<endl;
    }
  else if(verbose) cout<<"Imag period checks OK"<<endl;
#endif // CHECK_PERIODS
  return C;
}

bigfloat CG(int r, bigfloat x) // Cohen's G_r(x)
{
  // cout<<"CG("<<r<<","<<x<<") = ..."<<flush;
  static const bigfloat one = to_bigfloat(1);
  bigfloat n=one, emx=exp(-x), ans=x, term=x, y;
  vector<bigfloat> Av(r+1, one);  // indexed from 0 to r
  while(!is_approx_zero(emx*term*Av[r]))
    {
      n++;
      // update A-vector, term and sum
      for(int j=1;j<=r;j++)
        Av[j]+=(Av[j-1]/n);
      term*=x;
      term/=n;
      y = Av[r]*term;
      ans+=y;
      // cout<<"\n"<<y<<"\t"<<ans<<"\t"<<emx*ans;
      if (is_approx_zero(y/ans)) break;
    }
  // cout<<"\n"<<emx*ans<<" (using "<<n<<" terms)"<<endl;
  return emx*ans;
}

bigfloat Glarge(int r, bigfloat x) // Cohen's Gamma_r(x) for large x
{
  static const bigfloat zero = to_bigfloat(0);
  static const bigfloat one = to_bigfloat(1);
  static const bigfloat two = to_bigfloat(2);
  bigfloat emx2=2*exp(-x), ans=zero, mxinv=-one/x;
  bigfloat y, term = mxinv, n=zero;
  if (is_approx_zero(abs(emx2*term)))
    return zero;
  // cout<<"Glarge("<<r<<","<<x<<") = ..."<<flush;
  vector<bigfloat> Cv(r+1, zero);  // indexed from 0 to r
  Cv[0]=one;
  while((n<1000) && !is_approx_zero(abs(emx2*term)))
    {
      n++;
      // update C-vector, term and sum
      for(int j=r;j>0;j--) Cv[j]+=(Cv[j-1]/n);
      term*=n;
      term*=mxinv;
      y = Cv[r]*term;
      ans+=y;
      // cout<<"\n"<<y<<"\t"<<ans<<"\t"<<emx2*ans;
      if (is_approx_zero(y/ans)) break;
    }
  // cout<<"\n"<<emx2*ans<<" (using "<<n<<" terms)"<<endl;
  return emx2*ans;
}

bigfloat Q(int r, bigfloat x)  // Q_r(x) polynomial from AMEC p.44
{
#ifdef MPFP // Multi-Precision Floating Point
  static const bigint nz2=to_ZZ("3772654005711327105320428325461179161744295822071095339706353540767904529098322739007189721774317982928833");
  bigfloat zeta2; MakeRR(zeta2,nz2,-350);
  static const bigint nz3=to_ZZ("2756915843737840912679655856873838262816890298077497105924627168570887325226967786076589016002130138897164");
  bigfloat zeta3; MakeRR(zeta3,nz3,-350);
  static const bigint nz4=to_ZZ("2482306838570394152367457601777793352247775704274910416102594171643891396599068147834147756326957412925856");
  bigfloat zeta4; MakeRR(zeta4,nz4,-350);
#else
 static const bigfloat zeta2 = 1.6449340668482264364724151666460251892189499;
 static const bigfloat zeta3 = 1.20205690315959428539973816151144999076498629;
 static const bigfloat zeta4 = 1.08232323371113819151600369654116790277475095;
#endif
 static const bigfloat two = to_bigfloat(2);
 static const bigfloat three = to_bigfloat(3);
 static const bigfloat four = to_bigfloat(4);
 static const bigfloat nine = to_bigfloat(9);
 static const bigfloat sixteen = to_bigfloat(16);
 static const bigfloat twentyfour = to_bigfloat(24);
 static const bigfloat const1 = nine*zeta4/sixteen;
 static const bigfloat const2 = zeta3/three;
 static const bigfloat const3 = zeta4/four;
 static const bigfloat half = to_bigfloat(1)/two;
 static const bigfloat third = to_bigfloat(1)/three;
 static const bigfloat twenty4th = to_bigfloat(1)/twentyfour;
  switch(r) {
  case 1: default: return x;
  case 2: return (x*x+zeta2)*half;
  case 3: return x*(x*x*third+zeta2)*half-const2;
  case 4: return const1 +x*(-const2+x*(const3+x*x*twenty4th));
  }
}

bigfloat P(int r, bigfloat x)  // P_r(x) polynomial from AMEC p.44
{
  static bigfloat gamma=Euler();
  return Q(r,x-gamma);
}

bigfloat Gsmall(int r, bigfloat x) // Cohen's Gamma_r(x) for small x
{
  bigfloat a=P(r,-log(x)), b = CG(r,x);
  return (r%2? a+b : a-b);
}

bigfloat G(int r, bigfloat x)  // G_r(x)
{
#ifndef MPFP // Multi-Precision Floating Point
  static const  bigfloat x0 = to_bigfloat(14);
#else
  // log(2) = 0.693147180599..
  // this value only determines whether to use Gsmall or Glarge:
  static const  bigfloat x0 = to_bigfloat(0.693147*bit_precision());
#endif
  //  cout<<"switch point = "<<x0<<endl;
  //  cout<<"x="<<x<<endl;

  return x<x0? Gsmall(r,x) : Glarge(r,x);
}

bigfloat ldash1::G(bigfloat x)  // G_r(x)
{
#ifndef MPFP // if without Multi-Precision Floating Point, avoid underflow: log(MINDOUBLE)=-708.396
  if((r<2) && (x>708))
    return to_bigfloat(0);
  if(r==1)
    return myg1(x);
#endif
  return r==0? exp(-x) : ::G(r,x);
}

// myg1(x) is the same as G_1(x) = G(1,x)

//#define TRACEG1

bigfloat myg1(bigfloat x)
{
  static bigfloat gamma=Euler();
#ifndef MPFP // Multi-Precision Floating Point
  if(x>708) return to_bigfloat(0);
#endif
  if(x<2)
    {
      bigfloat ans = -log(x) - gamma;
      bigfloat p = to_bigfloat(-1);
#ifdef TRACEG1
      cout<<"Computing myg1 for x = "<<x<<" using series"<<endl;
#endif
      int ok=0;
      bigfloat term;
//The following does not work for large x!  (The terms get too big before
//they get small, and overflow destroys the result.)
      for (long n=1; (n<5000) && !ok; n++)
	{
	  p /=n; p*= -x;
	  term = p/n;
	  ans += term;
	  ok=is_approx_zero(term/ans);
#ifdef TRACEG1
	  cout<<"   n="<<n<<": term = " << term << ", ans = "<<ans<<endl;
#endif
	}
#ifdef TRACEG1
      cout << "returning g1 = " << ans << endl;
#endif
      return ans;
    }
  //  else  x>2, use continued fraction form from B-G-Z p.478
#ifdef TRACEG1
  cout<<"Computing myg1 for x = "<<x<<" using continued fraction"<<endl;
#endif
  bigfloat a0, a1, a2, b0, b1=x, b2; a0=0; b0=1;
  bigfloat ans, newans; ans=0;
  bigfloat ca;
  a1=exp(-x);
  for(long k=2; k<10000; k++)
    {
      if (k&1) //then k is odd
	{ ca = (k-1)/2;
	  a2=x*a1+ca*a0; a0=a1; a1=a2;
	  b2=x*b1+ca*b0; b0=b1; b1=b2;
	}
      else //k is even
	{
	  ca = k/2;
	  a2=a1+ca*a0; a0=a1; a1=a2;
	  b2=b1+ca*b0; b0=b1; b1=b2;
	}
      newans=a2/b2;
#ifdef TRACEG1
      cout<<"   k="<<k<<": approx = " << newans << endl;
#endif
      if (is_approx_zero(ans-newans))
	{
//	  cout << "returning g1 = " << newans << endl;
	  return newans;
	}
      ans=newans;
    }
  cout << "In function g1, continued fraction method, reached end of loop!\n";
  ans=0;
  return ans;
}
