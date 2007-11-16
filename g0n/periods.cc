// FILE PERIODS.CC : implementation of classes for integrating newforms
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
//
#include "compproc.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "cperiods.h"  // from qcurves
#include "newforms.h"
#include "periods.h"
#include "fixc6.h"

//#define CHECK_PERIODS // check that curves constructed from periods 
                      // have the same periods...

//#define TRACE_CACHE
//#define TRACE_USE

/////////////////////////////////////
//  functions for character class  //
/////////////////////////////////////
 
character::character(long m)
{
  modul=m;
  chartable = new int[m];
  init();
}

character::~character() 
{
  delete [] chartable;
}

void character::init()
{
  if ( modul==1 ) chartable[0] = 1;
  else { long i=modul; while(i--) chartable[i] =  legendre(i,modul);}
}

void character::reset(long m)
{
  delete [] chartable;
  modul=m;
  chartable = new int[m];
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
  long i, j, p, ip = iN->npdivs;
  for(i=0; i<nap; i++)
    { p = primelist[i];
    vector<long>::const_iterator pi = find(iN->plist.begin(),iN->plist.end(),p);
    if(pi==iN->plist.end()) // then p is good
      {
	aplist.push_back( apl[ip++]);
// 	cout << "i = "<<i<<",\tp = " << p << "\ta_p = " << aplist[i] << endl;
      }
    else // p is bad
      {
	if(::div(p*p,N)) 
	  {
	    aplist.push_back(0);
// 	    cout << "i = "<<i<<",\tp = " << p << "\ta_p = " << aplist[i] << endl;
	  }
	else 
	  {
	    j = pi-(iN->plist.begin());  // p is j'th bad prime
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
  nap = apl.size();
  primelist = primes(nap);   //First nap primes, indexed from 0
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
  long m2=n,m23,m235,m2357;
  long am2,am23,am235,am2357;
  long i2,i3,i5,i7;
  for(i2=0; (i2<=n2p)&&(m2<limit); i2++, m2*=2)
    {
      am2=an*a2p_cache[i2];
      m23=m2;
      for(i3=0; (i3<=n3p)&&(m23<limit); i3++, m23*=3)
	{
	  am23=am2*a3p_cache[i3];
	  m235=m23;
	  for(i5=0; (i5<=n5p)&&(m235<limit); i5++, m235*=5)
	    {
	      am235=am23*a5p_cache[i5];
	      m2357=m235;
	      for(i7=0; (i7<=n7p)&&(m2357<limit); i7++, m2357*=7)
		{
		  am2357=am235*a7p_cache[i7];
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
      if(ndiv(2,N)) a2p_cache[j]-=2*a2p_cache[j-2];
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
      if(ndiv(3,N)) a3p_cache[j]-=3*a3p_cache[j-2];
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
      if(ndiv(5,N)) a5p_cache[j]-=5*a5p_cache[j-2];
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
      if(ndiv(7,N)) a7p_cache[j]-=7*a7p_cache[j-2];
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
  while( (ip<primelist.size()) && (p<rootlimit)  )
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
  while( (ip<primelist.size()) && (p<=limit)  )
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
  type = f->type;
  dp0=f->dp0, mplus=f->mplus, mminus=f->mminus;
  if(f->lplus==1) mplus=f->np0; // L/P = dp0/np0 so now P=L*mplus/dp0
  initaplist(iN,f->aplist);
  rootmod =  sqrt(to_bigfloat(N));  
  factor1 =  exp(-(TWOPI)/ (to_bigfloat(chi1.modulus()) * rootmod));
  factor2 =  exp(-(TWOPI)/ (to_bigfloat(chi2.modulus()) * rootmod));
  long dp = decimal_precision();
 // MUST keep brackets around TWOPI!:
  bigfloat rootlimit1=(dp*LOG10-log(1-factor1))*rootmod/(TWOPI) ;
  bigfloat rootlimit2=(dp*LOG10-log(1-factor1))*rootmod/(TWOPI) ; 
  Iasb(limit1,rootlimit1);
  Iasb(limit2,rootlimit2);
#ifdef TRACE_USE 
  cout << "Decimal precision = "<<dp<<endl;
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
    cout<<"\nLARGEST PRIME "<<primelist[primelist.size()-1]<<" IS BELOW LIMIT "<< limit <<" required for decimal precision "<<dp<<endl;
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
  dotplus=f->dotplus, dotminus=f->dotminus;
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
    cerr<<"Problem: cannot compute periods for matrix with d=0!"<<endl;
    rp=ip=0;
    return;
   }

  if (d<0) { a=-a;b=-b;c=-c;d=-d;}

  bigfloat drecip =  to_bigfloat(1) / to_bigfloat(d);
  theta1 = to_bigfloat(b) * drecip;
  theta2 = to_bigfloat(c) * drecip;
  factor2 = factor1 * drecip;
  long dp = decimal_precision();
  Iasb(limit1,(-dp*LOG10-log((1-exp(factor1))/3))/(factor1)) ;
  Iasb(limit2,(-dp*LOG10-log((1-exp(factor2))/3))/(factor2)) ;

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
    cout<<"\nLARGEST PRIME "<<primelist[primelist.size()-1]<<" IS BELOW LIMIT "<< limit <<" required for decimal precision "<<dp<<endl;
*/
#ifdef TRACE_CACHE
  cout << "Initial an_cache = "<<an_cache<<endl;
  cout << "Calling sumit()" << endl;
#endif

  sumit();

#ifdef TRACE_CACHE
  cout << "Returned from sumit()" << endl;
#endif
  //  cout<<"rp = "<<rp<<endl;
  //  cout<<"ip = "<<ip<<endl;
  rp = sum1; if(dotplus!=1)  rp/=to_bigfloat(dotplus);
  ip = sum2; if(dotminus!=1) ip/=to_bigfloat(dotminus);
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
  bigfloat dan = to_bigfloat(an);
  bigfloat coeff = -dan/dn;
  bigfloat ef2 = coeff * exp(dn*factor2);
  bigfloat dn2pi = dn*(TWOPI);
  bigfloat c1 = ef2*cos(dn2pi*theta1),
           s1 = ef2*sin(dn2pi*theta1),
           c2 = ef2*cos(dn2pi*theta2),
           s2 = ef2*sin(dn2pi*theta2);
  if(eps_N==-1) 
    {
      if(n<limit1)
	{
	  bigfloat ef1 = coeff * exp(dn*factor1);
#ifdef TRACE_USE 
	  cout<<"n="<<n<<", ef1="<<ef1<<", ef2="<<ef2<<"\n";
#endif
	  sum1 += (2*ef1); 
	}
      sum1 -= (c1+c2);
      sum2 -= (s1+s2);
    }
  else
    {
      sum1 += (c1-c2);
      sum2 += (s1-s2);
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
  long dp = decimal_precision();
  Iasb(limit,dp*LOG10/y0);
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
  long maxp = prime_number(nap);
  limit  =  I2long(Ifloor((15+decimal_precision())*log(10)/factor1));
  if(limit>maxp) limit=maxp;
  limit1 = limit;
  rootlimit=sqrt(to_bigfloat(limit));
  an_cache.resize(I2long(Ifloor(rootlimit+1)),0);
#ifdef TRACE_CACHE
  cout << "Initial an_cache = "<<an_cache<<endl;
#endif
//decide on r (first estimate)
  computed=0;
  r=0;  g = &myg0;

  if(f_sfe==-1) 
    {
	{r=1; g = &myg1;}
    }
  else if(f_loverp==0) {r=2; g = &myg2;}
}

void ldash1::compute(void) 
{
  if(computed) return;
  sumit(); ld1=to_bigfloat(2)*sum1; computed=1;
  if(r==0) return;
  if(abs(ld1)>0.0001) return;
  if(r==1) // suspect rank is in fact 3
    {
//    cout<<"\nSwitching rank from 1 to 3 since L'(1) = "<<ld1<<endl;
      r=3;  g=&myg3;  sumit(); ld1=2*sum1;
//    cout<<"Now L^(3)(1) = "<<ld1<<endl;
      return;
    }
  cout<<"\n!!! L^("<<r<<")(1) = "<<ld1<<" !!!"<<endl;
  r+=2;
  cout<<"Do we have a rank "<<r<<" curve?"<<endl;
  cout<<"If so we need to implement G_"<<r
      <<"(x) to compute L^("<<r<<")(1)"<<endl;
}

/////////////////////////////////
//  functions for lfchi class  //
/////////////////////////////////
 
lfchi::lfchi (const level* iN, const newform* f)
{
  initaplist(iN, f->aplist);
  rootmod =  sqrt(to_bigfloat(N));
  long dp = decimal_precision();
  Iasb(limit0,dp*LOG10*rootmod/(TWOPI)) ; // brackets essential
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

  val = 2*sum1;
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
      if(verbose) cout<<"L(f,1) = "<<abs(lx.value())<<"; ";
      rational lop=nfi->loverp;
      x = abs(lx.value()/to_bigfloat(lop)); 
      if(verbose) cout<<"real period = "<<x<<endl;
      return 1;
    }
  
  long lplus=abs(nfi->lplus);
  long mplus=abs(nfi->mplus);
  if(mplus!=0)
    {
      if(verbose) 
	cout << "Computing real period via L(f,chi,1) with chi mod "
	     <<lplus<<"...";
      lx.compute(lplus);
      if(verbose) cout<<"L(f,chi,1) = "<<abs(lx.scaled_value())<<"; ";
      x = abs(lx.scaled_value()/to_bigfloat(mplus));
      if(verbose) cout<<"real period = "<<x<<endl;
      return 1;
    }
  
  // we only reach here if sfe=-1 and level is square  
  if(verbose) 
    cout << "Unable to compute real period via L(f,chi,1)..."<<endl;
  return 0;
}

Cperiods newforms::getperiods(long i, int method, int verbose)
  // method=1 to use periods_direct
  // method=0 to use periods_via_lfchi
  // method=-1 (default) to use whichever is best
{
  newform* nfi = &(nflist[i]);
  if(method==-1) // find and use best method
    {
      if(this->squarelevel) 
	method=1;
      else
	{
	  long d = abs(nfi->d);
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
	 Cperiods cp = per.getperiods();
         return cp;
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
  if(verbose) cout<<cp<<endl;
  cp.getwRI(wR, wRI);
  rperiod = real(wR);
  cp.getwi(w1, w2);
  getc4c6(w2,w1,c4,c6);  // from compproc.h
  bigfloat rc4 = real(c4), rc6 = real(c6);
  if(verbose) cout << "c4 = " << rc4 << "\nc6 = " << rc6 << endl;
  bigint ic4 = fac*Iround(rc4/fac);
  bigint ic6 = fac6*Iround(rc6/fac6);
  if(verbose) cout << "After rounding";
  if(verbose&&(fac>1)) 
    cout << ", using factors " << fac << " for c4 and " << fac6 << " for c6";
  if(verbose) cout<<":\n";
  if(verbose) cout << "ic4 = " << ic4 << "\nic6 = " << ic6 << endl;
  // To fix the c4 or c6 values insert data in files fixc4.data and
  // fixc6.data; NB the index here (i) starts at 0, but class fixc6
  // adjusts so in the data files, start at 1

#ifndef MPFP // Multi-Precision Floating Point
  c4c6fixer(n,i,ic4,ic6);
  if(verbose) cout << "After fixing: \n";
  if(verbose) cout << "ic4 = " << ic4 << "\nic6 = " << ic6 << endl;
#endif
  Curve C(ic4,ic6);
#ifdef CHECK_PERIODS
  if(C.isnull()) return C;
// Check periods were correct:
  Curvedata CD(C,1);  // The 1 causes minimalization
  Cperiods cpC(CD);
  bigcomplex wRC, wRIC;
  cpC.getwRI(wRC, wRIC);
  wR=abs(wR); wRI=abs(imag(wRI));
  wRC=abs(wRC); wRIC=abs(imag(wRIC));
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
	  <<" of the newform"<<endl;
      cout<<"Real period of C: "<<real(wRC)<<endl;
      cout<<"Real period of f: "<<real(wR)<<endl;
      cout<<"Ratio = "<<real(wR)/real(wRC)<<endl;
    }
  else if(verbose) cout<<"Real period checks OK"<<endl;
  if((abs((wRI-wRIC)/wRIC)>0.0001))
    {
      cout<<"Second period of constructed curve does not match that"
	  <<" of the newform"<<endl;
      cout<<"Imag part of second period of C: "<<real(wRIC)<<endl;
      cout<<"Imag part of second period of f: "<<real(wRI)<<endl;
      cout<<"Ratio of imaginary parts = "<<real(wRI)/real(wRIC)<<endl;
    }
  else if(verbose) cout<<"Imag period checks OK"<<endl;
#endif // CHECK_PERIODS
  return C;
}

// avoid underflow: log(MINDOUBLE)=-708.396

bigfloat myg0(bigfloat x)
{
#ifndef MPFP // Multi-Precision Floating Point
  if(x>708) return to_bigfloat(0);
#endif
  return exp(-x);
}

//#define TRACEG1

bigfloat myg1(bigfloat x)
{
#ifndef MPFP // Multi-Precision Floating Point
  if(x>708) return to_bigfloat(0); 
#endif
  if(x<2)
    {
      bigfloat ans = -log(x) - Euler();
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

bigfloat CG(int r, bigfloat x) // Cohen's G_r(x)
{
  bigfloat emx=exp(-x), ans=x, term=x;
  vector<bigfloat> Av(r+1);  // indexed from 0 to r
  int j, n=1;
  for(j=0;j<=r;j++) Av[j]=to_bigfloat(1);
  while(!is_approx_zero(emx*term*Av[r]))
    //while(!is_approx_zero(term*Av[r]))
    {
      n++;
      // update A-vector, term and sum
      for(j=1;j<=r;j++) Av[j]+=(Av[j-1]/to_bigfloat(n)); 
      term*=(x/to_bigfloat(n));
      ans+=(Av[r]*term);
    }
  return emx*ans;
}

bigfloat Glarge(int r, bigfloat x) // Cohen's Gamma_r(x) for large x
{
  bigfloat emx=exp(-x), ans=to_bigfloat(0), term=-to_bigfloat(1)/x;
  vector<bigfloat> Cv(r+1);  // indexed from 0 to r
  int j, n=0;
  Cv[0]=to_bigfloat(1);
  for(j=1;j<=r;j++) Cv[j]=to_bigfloat(0); 
  //  cout<<"emx*term="<<emx*term<<endl;
  while(!is_approx_zero(abs(emx*term)))
    //while(!is_approx_zero(abs(term)))
    {
      n++;
      // update C-vector, term and sum
      for(j=r;j>0;j--) Cv[j]+=(Cv[j-1]/to_bigfloat(n)); 
      term*=(-to_bigfloat(n)/x);
      ans+=(Cv[r]*term);
      //      cout<<"term="<<term<<"; ans="<<ans<<endl;
    }
  return to_bigfloat(2)*emx*ans;
}

bigfloat Q(int r, bigfloat x)  // Q_r(x) polynomial from AMEC p.44
{
#ifdef MPFP // Multi-Precision Floating Point
  static const bigint nz2=atoI("3772654005711327105320428325461179161744295822071095339706353540767904529098322739007189721774317982928833");
  bigfloat zeta2; MakeRR(zeta2,nz2,-350);
  static const bigint nz3=atoI("2756915843737840912679655856873838262816890298077497105924627168570887325226967786076589016002130138897164");
  bigfloat zeta3; MakeRR(zeta3,nz3,-350);
  static const bigint nz4=atoI("2482306838570394152367457601777793352247775704274910416102594171643891396599068147834147756326957412925856");
  bigfloat zeta4; MakeRR(zeta4,nz4,-350);
#else
 const bigfloat zeta2 = 1.6449340668482264364724151666460251892189499;
 const bigfloat zeta3 = 1.20205690315959428539973816151144999076498629;
 const bigfloat zeta4 = 1.08232323371113819151600369654116790277475095;
#endif
  switch(r) {
  case 1: default: return x;
  case 2: return (x*x+zeta2)/to_bigfloat(2);
  case 3: return x*(x*x/to_bigfloat(3)+zeta2)/to_bigfloat(2)-zeta3/to_bigfloat(3);
  case 4: return to_bigfloat(9)*zeta4/to_bigfloat(16)+x*(-zeta3/to_bigfloat(3)+x*(zeta2/to_bigfloat(4)+x*x/to_bigfloat(24)));
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
  bigfloat x0 = 
#ifndef MPFP // Multi-Precision Floating Point
    to_bigfloat(14);
#else
    to_bigfloat(log(10)*decimal_precision());
#endif
    //  cout<<"switch point = "<<x0<<endl;
    //  cout<<"x="<<x<<endl;
  if(x<x0) 
    {
      return Gsmall(r,x);
    }
  else
    {  
      return Glarge(r,x);
    }
}

bigfloat ldash1::G(bigfloat x)  // G_r(x)
{
  switch(r) {
  case 0: return myg0(x);
  case 1: return myg1(x);
  default: return ::G(r,x);
  }
}
 
bigfloat myg2(bigfloat x)
{
  static bigfloat zero=to_bigfloat(0);
  static bigfloat one=to_bigfloat(1);
  static bigfloat two=to_bigfloat(2);
  static bigfloat twelve=to_bigfloat(12);
  static bigfloat gamma=Euler();
  if (x>20) return zero;
//#define TRACEG2
  bigfloat ans = -log(x) - gamma;
#ifdef TRACEG2
  cout<<"Computing myg2 for x = "<<x<<endl;
#endif
  ans=ans*ans/two + PI*PI/twelve;
  bigfloat p = one;
  int ok=0;
  bigfloat term;
  for (long n=1; (n<500) && !ok; n++)
    {
      p /=n;  p*= -x;
      term = (p/n)/n;
      ans += term;
      ok=is_approx_zero(term/ans);
#ifdef TRACEG2
      // cout<<"\tn="<<n<<" \tans = "<<ans<<endl;
#endif
    }
#ifdef TRACEG2
  cout<<"...returning ans = "<<ans<<endl;
#endif
   return ans;
}
 
bigfloat myg3(bigfloat x)
{
  static bigfloat zero=to_bigfloat(0);
  static bigfloat one=to_bigfloat(1);
  static bigfloat three=to_bigfloat(3);
  static bigfloat twelve=to_bigfloat(12);
#ifdef MPFP // Multi-Precision Floating Point
  static const bigint nz3=atoI("2756915843737840912679655856873838262816890298077497105924627168570887325226967786076589016002130138897164");
  bigfloat zeta3; MakeRR(zeta3,nz3,-350);
#else
 const bigfloat zeta3 = 1.20205690315959428539973816151144999076498629;
#endif
  if (x>20) return zero;
  bigfloat ans = -log(x) - Euler();
//cout<<"Computing myg3 for x = "<<x<<endl;
  ans = (PI*PI + 2*ans*ans) * ans/twelve - zeta3/three;
  bigfloat p = -one;
  bigfloat term;
  int ok = 0;
  for (long n=1; (n<500) && !ok; n++)
    {
      p *= -x/n;
      term = ((p/n)/n)/n;
      ans += term;
      ok=is_approx_zero(term/ans);
      //cout<<"   N="<<n<<" ans = "<<ans<<endl;
    }
   return ans;
}

long* getai(long c4, long c6)
{
//cout << "In getai with c4 = " << c4 << " and c6 = " << c6 << endl;
  long b2,b4,b6;
  b2 = (
             (c4%3 == 0 ? (odd(c4) ?   3 :   0)
              : (odd(c4) ?  -1 :  -4))
             *c6)
    % 12;
  if ( b2 < -5 ) { b2 += 12; }
  if ( b2 > 6 ) { b2 -= 12; }
//cout << "b2 = " << b2 << endl;
  b4 = (b2*b2-c4) / 24;
//cout << "b4 = " << b4 << endl;
  b6 = (-b2*b2*b2 + 36 * b2*b4 - c6) / 216;
//cout << "b6 = " << b6 << endl;
  long* ans = new long[5];
  ans[0] = (odd(b2) ? 1 : 0);
  ans[2] = (odd(b6) ? 1 : 0);
  ans[1] = (b2-ans[0]) / 4;
  ans[3] = (b4-ans[0]*ans[2]) / 2;
  ans[4] = (b6-ans[2]) / 4;
//cout << "a1,a2,a3,a4,a6 = " << ans[0] <<","<< ans[1] <<","<< ans[2] <<","<< ans[3] <<","<< ans[4] << endl;
  return ans;
}
