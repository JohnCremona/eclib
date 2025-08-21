// marith.cc: implementations of integer arithmetic functions (multiprecision)
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
 
#include <NTL/ZZXFactoring.h>

// We use the pari C library for factoring (via a string interface
// defined in parifact.h/cc)

#include <sstream>
#include <eclib/marith.h>

// Utilities for debugging output -- for example, from gdb you can give the command
//   p show(a)
// to see a where a is a bigint or vector<bigint>

bigint show(const bigint& a) {cout<<a<<endl; return a;}
vector<bigint> show(const vector<bigint>& a) {cout<<a<<endl; return a;}

// integers and rationals

bigint bezout(const bigint& aa, const bigint& bb, bigint& xx, bigint& yy)
{bigint ans; XGCD(ans,xx,yy,aa,bb); return ans;}
int divides(const bigint& a, const bigint& b, bigint& q, bigint& r)
  { DivRem(q,r,a,b); return IsZero(r);}
int divides(const bigint& a, long b, bigint& q, long& r)
  { r=DivRem(q,a,b); return (r==0);}

// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
bigint rounded_division(const bigint& a, const bigint& b)
{
  bigint q, r;
  DivRem(q,r,a,b);
  bigint r2 = r<<1;
  return (r2<-b? q-1: (r2>=b? q+1: q));
}

// oddsqrt works on odd n, called by isqrt
//
//#define oddsqrt(root,n) sqrt(root,n)   // builtin
//#define oddsqrt(root,n) sqrtq2(root,n) // 2-adic
#define oddsqrt(root,n) sqrtnr(root,n) // JC's Newton

// 2-adic version of isqrt:

int sqrtq2(bigint& root, const bigint& n)
{
  bigint a,r; long a0;
  ::divides(n,(long)8,r,a0);
  if(a0!=1) return 0;              // odd squares must be 1 mod 8
  if(r==1) {a0=3; r=0;}           // special case
  a=a0;
//  cout<<"odd part 1 mod 8 with quotient r = " << r << endl;
  bigint twok(8), twok3(1);
  long kminus1=2;
  while(r>0)
    {
      if(even(r))
	{
	  rshift(r,1,r);
	}
      else
	{ 
	  subx(r,a,r); rshift(r,1,r); subx(r,twok3,r);
	  setbit(a,kminus1);
	}
      lshift(twok,1,twok);  lshift(twok3,1,twok3);
      kminus1++;
//      cout<<"a="<<a<<", r="<<r<<endl;
    }
  if(is_zero(r))     {root=a; return 1;}
  if(r+a==(twok>>2)) {root=(twok>>1)-a; return 1;}
  return 0;
}

// Newton-R-type iteration as in Henri Cohen's book pp 38-39

  // some arrays borrowed from pari:
static int carresmod64[]={
    1,1,0,0,1,0,0,0,0,1, 0,0,0,0,0,0,1,1,0,0,
    0,0,0,0,0,1,0,0,0,0, 0,0,0,1,0,0,1,0,0,0,
    0,1,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,1,0,0,
    0,0,0,0};
static int carresmod63[]={
    1,1,0,0,1,0,0,1,0,1, 0,0,0,0,0,0,1,0,1,0,
    0,0,1,0,0,1,0,0,1,0, 0,0,0,0,0,0,1,1,0,0,
    0,0,0,1,0,0,1,0,0,1, 0,0,0,0,0,0,0,0,1,0,
    0,0,0};
static int carresmod65[]={
    1,1,0,0,1,0,0,0,0,1, 1,0,0,0,1,0,1,0,0,0,
    0,0,0,0,0,1,1,0,0,1, 1,0,0,0,0,1,1,0,0,1,
    1,0,0,0,0,0,0,0,0,1, 0,1,0,0,0,1,1,0,0,0,
    0,1,0,0,1};
static int carresmod11[]={1,1,0,1,1,1,0,0,0,1,0};
static int carresmod17[]={1,1,1,0,1,0,0,0,1,1,0,0,0,1,0,1,1};
static int carresmod19[]={1,1,0,0,1,1,1,1,0,1,0,1,0,0,0,0,1,1,0};
/*
static int carresmod23[]={1,1,1,1,1,0,1,0,1,1,0,0,1,1,0,0,1,0,1,0,0};
static int carresmod29[]={1,1,0,0,1,1,1,1,0,1,0,0,0,1,0,0,1,0,0,0,1,
			  0,1,1,1,1,0,0,1};
static int carresmod31[]={1,1,1,0,1,1,0,1,1,1,1,0,0,0,1,0,1,0,1,1,1,
			  0,0,0,0,1,0,0,1,0,0};
static int carresmod37[]={1,1,0,1,1,0,0,1,0,1,1,1,1,0,0,0,1,0,0,0,0,
			  1,0,0,0,1,1,1,1,0,1,0,0,1,1,0,1};
static int carresmod41[]={1,1,1,0,1,1,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,
			  1,0,1,0,1,0,0,0,0,0,1,1,1,0,0,1,1,0,1,1};
static int carresmod43[]={1,1,0,0,1,0,1,0,0,1,1,1,0,1,1,1,1,1,0,0,0,
			  1,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,0,1,0,1,1,0};
*/

int sqrtnr(bigint& root, const bigint& n)
{
  long l = 1+(lg(n)>>1);
  //  cout << "n = " << n << ", l = " << l << endl;
  bigint y;
  root=1; root<<=l;  // first approx, > sqrt(n)
  //  cout << "first approx root = " << root << endl;
  while(1)
    {
      divx(n,root,y); addx(y,root,y); rshift(y,1,y);
      if(y<root) {root=y;} else break;
      //      cout << "root = " << root << endl;
    }
  return (sqr(root)==n);
}

//NB The code here proved faster than any of the version builtin to NTL

int isqrt(const bigint& in, bigint& root)
{
//  cout<<"In isqrt with n = " << in << endl;
  root=0;
  if(sign(in)<0) return 0;
  if(sign(in)==0) return 1;
  long m,twopow=0;
  bigint n(in);
  while(even(n)) {n>>=1; twopow++;}
//  cout << "2-power = " << twopow << endl;
  if(odd(twopow)) return 0;            // 2 | a to an odd power
  twopow>>=1;                          // = power of 2 in root
  m = bigint_mod_long(n,931170240); // 931170240=64*63*65*11*17*19 < 2^30

  if (!carresmod64[m&63]) return 0;
  if (!carresmod63[m%63]) return 0;
  if (!carresmod65[m%65]) return 0;
  if (!carresmod11[m%11]) return 0;
  if (!carresmod17[m%17]) return 0;
  if (!carresmod19[m%19]) return 0;

  if(oddsqrt(root,n)) {lshift(root,twopow,root); return 1;}
  return 0;
}


int sqrt_mod_2_power(bigint& x, const bigint& a, int e)
{
  if(e==0) {x=0; return 1;}
  long a8 = posmod(a,8); // I2long(a%8);
                         //  if(a8<0) {a8+=8;}
  if(!(a8%2)) return 0;  // odd a only
  x=1;
  if(e==1) return 1;
  if(e==2) return ((a8%4)==1);
  if(a8!=1) return 0;
  if(e==3)  return 1;
  // Now e>=4 and a=1 (mod 8)  
  int k; bigint q1, q, q2; q1=4, q=8, q2=16;
  for(k=3; k<e; k++)
    {
      if(ndiv(q2,(sqr(x)-a))) x+=q1;
      q1=q; q=q2; q2*=2;
    }
#ifdef CHECK_SQRT_MOD
  if(ndiv(q,sqr(x)-a))
    cout<<"Error in sqrt_mod_2_power with a="<<a<<", a mod 8="<<a8<<", e="<<e
	<<": returns "<<x<<endl;
#endif  
  return 1;
}

int sqrt_mod_p_power(bigint& x, const bigint& a, const bigint& p, int e)
{
  if(p==2) {return sqrt_mod_2_power(x,a,e);}
  if(e==0) {x=0; return 1;}
  bigint a1 = a%p;
  if(a1==0) return 0;  // p ndiv a only
  if(legendre(a1,p)==-1) return 0;
  if(a1<0) a1+=p;  // since sqrt_mod_p wants it between 0 and p-1
  sqrt_mod_p(x,a1,p);  
  //  cout<<"sqrt("<<a1<<" mod "<<p<<") = "<<x<<endl;
  if(e==1) {return 1;}
  bigint s = invmod(2*x,p);

  int k; bigint q=p;
  for(k=1; k<e; k++)
    {
      q*=p;
      x-= s*(sqr(x)-a)%q;
      x%=q;
    }
#ifdef CHECK_SQRT_MOD
  if(ndiv(q,sqr(x)-a))
    cout<<"Error in sqrt_mod_p_power with a="<<a<<", p="<<p<<", e="<<e
	<<": returns "<<x<<endl;
#endif  
  return 1;
}

int sqrt_mod_m(bigint& x, const bigint& a, const bigint& m)
{
  static const bigint zero(0), one(1);
  // Some trivial cases require no work:
  if(is_one(m))  {x=zero; return 1;}
  if(is_zero(a)) {x=zero; return 1;}
  if(is_one(a))  {x=one; return 1;}
#ifdef CHECK_SQRT_MOD
  cout<<"Factorizing "<<m<<"..."<<flush;
#endif  
  vector<bigint> mpdivs = pdivs(m);
#ifdef CHECK_SQRT_MOD
  cout<<"prime factors are "<<mpdivs<<endl;
#endif  
  return sqrt_mod_m(x,a,m,mpdivs);
}

int sqrt_mod_m(bigint& x, const bigint& a, const bigint& m, const vector<bigint>& mpdivs)
{
  // Some trivial cases require no work:
  if(is_one(m))  {x=0; return 1;}
  if(is_zero(a)) {x=0; return 1;}
  if(is_one(a))  {x=1; return 1;}
  bigint mm, xp, q; int e;
  x=0;  mm=1;
  
  for( const auto& p : mpdivs)
    {
      e = val(p,m);
      if(e==0) continue;
      if(p==2)
	{if(!sqrt_mod_2_power(xp,a,e)) return 0;}
      else
	{if(!sqrt_mod_p_power(xp,a,p,e)) return 0;}
      q=pow(p,e);
      if(p==mpdivs.front())
	x=xp;
      else
	x=chrem(x,xp,mm,q);
      mm*=q;
    }
#ifdef CHECK_SQRT_MOD
  if(ndiv(m,sqr(x)-a))
    cout<<"Error in sqrt_mod_m with a="<<a<<", m="<<m
	<<": returns "<<x<<endl;
#endif  
  return 1;
}

int modsqrt(const bigint& a, const vector<bigint>& bplist, bigint& x)
     // Solves x^2=a mod b, returns success/fail
{
  // Assumes b square-free, primes factors in bplist
  bigint u, v, amodp, xmodp, m;
  x=0; m=1;
  for( const auto& p : bplist)
    {
      if(p==2)
	{
	  xmodp=odd(a);
	}
      else // odd p
	{
	  amodp = a%p;
	  if(is_zero(amodp)) xmodp=0;
	  else
	    {
	      if(legendre(amodp,p)==-1) return 0;
	      if(amodp<0) amodp+=p;  // sqrt_mod_p wants it between 0 and p-1
	      sqrt_mod_p(xmodp,amodp,p);
	    }
	}
      // Now Chinese xmodp with previous (x mod m)
      bezout(m,p,u,v);
      x = x*v*p+xmodp*u*m;
      m*=p;
      x = mod(x,m);
    }
  return 1;
}

//
// bigint divisor lists etc
//

extra_prime_class the_extra_primes;  // The one and only instance

void initprimes(const string pfilename, int verb)
{
  if(verb) 
    {
      cout<<"Computed " << nprimes() << " primes, ";
      cout << "largest is " << maxprime() << "\n";
    }
  the_extra_primes.read_from_file(pfilename,verb);
  if(verb) the_extra_primes.show();
}

extra_prime_class::~extra_prime_class()
{
  //  write_to_file(string("PRIMES").c_str());
}

void extra_prime_class::write_to_file(const string pfilename, int verb)
{
  if(the_primes.size()==0) return;
  if(verb) cout << "writing primes to file " << pfilename << endl;
  ofstream pfile(pfilename.c_str());
  copy(the_primes.begin(),the_primes.end(), ostream_iterator<bigint>(pfile, "\n"));  
  if(verb) cout << "finished writing primes to file " << pfilename << endl;
}

void extra_prime_class::read_from_file(const string pfilename, int verb)
{
  ifstream pfile(pfilename.c_str());
  if(!pfile)  // default: no primes file exists
    {
      return;
    } 
  pfile>>ws;
  if(pfile.eof())  //  primes file exists but is empty
    {
      return;
    } 
  if(verb) cout << "reading primes from file " << pfilename << endl;
  bigint xp;
  while(pfile>>xp>>ws, (xp!=0) )
	{
	  if(verb) cout << "read extra prime " << xp << endl;
	  the_extra_primes.add(xp);
	  if(pfile.eof()) break;
	}
  if(verb) cout << "finished reading primes from file " << pfilename << endl;
}


// n>0 will be changed;  returns prime factors from factor base and divides out from n

vector<bigint> pdivs_use_factorbase(bigint& n, const std::set<bigint>& factor_base)
{
  vector<bigint> plist;
  if(n<2) return plist;
  for (const auto& p : factor_base)
    {
      if (n==1) break;
      if(divide_out(n,p))
	plist.push_back(p);
    }
  return plist;
}

// n>0 will be changed;  returns prime factors p<pmax and divides out from n

vector<bigint> pdivs_trial_div(bigint& n, const bigint& pmax)
{
  vector<bigint> plist;
  if(n<2) return plist;
  primevar pr;
  long p=2, r; bigint mp, q; mp=2;
  while ( (n>1) && (pr.ok()) && (mp<=pmax))
    { 
      if (::divides(n,p,q,r))   // then found a prime factor
	{
	  plist.push_back(mp); // add it to the list
	  n=q;
	  divide_out(n,mp);     // divide it out from n
	}
      // Now we might be able to conclude that the cofactor is prime:
      if(n>1) if (sqr(mp)>n) 
	{
	  plist.push_back(n); 
	  the_extra_primes.add(n); 
	  n=1; 
	}
      pr++; p = pr.value(); mp=p;
    }
  return plist;
}

vector<bigint> pdivs_trial(const bigint& number, int trace)
{
  if(trace) cout<<"In pdivs_trial() with number = " << number << endl;

  vector<bigint> plist;
  bigint n = abs(number), q, mp, mr;
  if(n<2) return plist;
  // use prime base first...

  plist=pdivs_use_factorbase(n,the_extra_primes.the_primes);
  if(n<2) return plist;
  if(trace) cout<< "After using factor base, n= " <<n<<", plist = "<< plist << endl;

  plist = vector_union(plist,pdivs_trial_div(n, bigint(maxprime())));
  if(trace) cout<< "After using trial division, n= " <<n<<", plist = "<< plist << endl;

  if(n>1) if(ProbPrime(n)) 
    {plist.push_back(n); the_extra_primes.add(n); n=1; }

  if (n>1) // failed to factor it 
    {
      cout<<"\n***Failed to find prime factor for composite "<<n<<" using trial division factorization of "<<number<<endl;
      cout<<"*** --appending "<<n<<" to its list of prime divisors"<<endl;
      plist.push_back(n);
    }
  if(trace) cout<< "pdivs_trial() returns " << plist << endl;
  return plist;
}

#include <eclib/parifact.h>

int
is_prime(const bigint& n)
{
  ostringstream oss;
  oss<<n;
  return is_prime(oss.str().c_str());
}

vector<bigint>
read_vec_from_string(string vecstr)
{
  //  cout<<"parsing output string "<<vecstr<<endl;
  vector<bigint> plist;
  istringstream vecin(vecstr);
  bigint p;
  char c;
  vecin>>skipws>>c; // swallow leading "["
  while(c!=']')
    {
      vecin>>p;
      //      cout<<"Reading p="<<p<<" from string"<<endl;
      plist.push_back(p);
      vecin>>skipws>>c; // swallow ",", but it might turn out to be "]"
    }
  //  cout<<"Finished reading from string"<<endl;
  return plist;
}

vector<bigint>
factor(const bigint& n, int proof=1)
{
  ostringstream oss;
  oss<<n;
  vector<bigint> plist =  read_vec_from_string(factor(oss.str()));
  if(proof)
    for( const auto& p : plist)
      {
	if(!is_prime(p))
	  {
	    cout<<"WARNING:  pari's factor() returned p="<<p
		<<" for which pari's isprime(p) FAILS!! Please report.";
	  }
      }
  return plist;
}

// The following uses pari's factorization function.
// However, numbers less than
#define TRIAL_DIV_BOUND bigint(100000000)
// will be handled by trial division, and the libpari function will
// only be called once primes factors less than
#define TRIAL_DIV_PRIME_BOUND bigint(10000)
// have been divided out,  to reduce the overheads involved.

vector<bigint> pdivs_pari(const bigint& number, int trace)
{
  vector<bigint> plist;
  bigint n=abs(number);
  if(n<2) return plist; // empty!

  // for small n just use trial division...

  if(n<TRIAL_DIV_BOUND) 
    {
      return pdivs_trial(n,trace); 
    }
  if(trace) cout<<"pdivs_pari factoring "<<n<<endl;

  // use prime base first...

  plist=pdivs_use_factorbase(n,the_extra_primes.the_primes);
  if(trace&&plist.size()>0) 
    cout<<"after using factorbase, have factors "<<plist
	<<", and cofactor = "<<n<< endl;
  if(n<2) 
    {
      sort(plist.begin(),plist.end());
      return plist;
    }

  // now use small primes...

  plist = vector_union(plist,pdivs_trial_div(n,TRIAL_DIV_PRIME_BOUND));
  if(trace&&plist.size()>0) 
    cout<<"after using trial division up to "<<TRIAL_DIV_PRIME_BOUND<<", have factors "<<plist
	<<", and cofactor = "<<n<< endl;

  if(n<2) 
    {
      sort(plist.begin(),plist.end());
      return plist;
    }

  // finally call factor() which interfaces to libpari...

  plist = vector_union(plist,::factor(n));
  sort(plist.begin(),plist.end());
  if(trace) cout<<"pdivs_pari returns "<<plist<<endl;
  return plist;
}

vector<bigint> pdivs(const bigint& number, int trace)
{
  return pdivs_pari(number);
}

vector<bigint> posdivs(const bigint& number)
{
 const vector<bigint>& plist=pdivs(number);
 return posdivs(number, plist);
}

vector<bigint> posdivs(const bigint& number, const vector<bigint>& plist)
{
  static const bigint one(1);
  vector<bigint> dlist = {one};
  for (auto p: plist)
    dlist = multiply_lists(powers(p, val(p,number)), dlist);
  return dlist;
}

vector<bigint> alldivs(const bigint& number)
{
  const vector<bigint>& plist=pdivs(number);
  return alldivs(number, plist);
}

vector<bigint> alldivs(const bigint& number, const vector<bigint>& plist)
{
  static const bigint one(1);
  return multiply_lists(posdivs(number, plist), {one, -one});
}

vector<bigint> sqdivs(const bigint& number)
{
  const vector<bigint>& plist=pdivs(number);
  return sqdivs(number, plist);
}

vector<bigint> sqdivs(const bigint& number, const vector<bigint>& plist)
{
  static const bigint one(1);
  vector<bigint> dlist = {one};
  for (auto p: plist)
    dlist = multiply_lists(powers(p, val(p,number)/2), dlist);
  return dlist;
}

vector<bigint> sqfreedivs(const bigint& number)
{
 const vector<bigint>& plist=pdivs(number);
 return sqfreedivs(number, plist);
}

vector<bigint> sqfreedivs(const bigint& number, const vector<bigint>& plist)
{
 static const bigint one(1);
  vector<bigint> dlist = {one};
  for (auto p: plist)
    dlist = multiply_lists(powers(p, 1), dlist);
 return dlist;
}

void sqfdecomp(const bigint& a, bigint& a1, bigint& a2, vector<bigint>& plist, int trace_fact)
     // a must be non-zero, computes square-free a1 and a2>0 such that a=a1*a2^2
     // plist will hold prime factors of a1
{
  plist = pdivs(a, trace_fact);
  sqfdecomp(a,plist,a1,a2);
}

void sqfdecomp(const bigint& a, vector<bigint>& plist, bigint& a1, bigint& a2)
     // a must be non-zero, computes square-free a1 and a2>0 such that a=a1*a2^2
     // plist already holds prime factors of a
     // plist will hold prime factors of a1
{
  long j;
  vector<bigint> aplist;
  a1=1;  a2=1;
  for (const auto& p : plist)
    {
      long e = val(p,a);
      if(e==0) continue;
      if(e&1) {a1*=p; aplist.push_back(p);}
      e >>= 1;
      for(j=0; j<e; j++) a2*=p;
    }
  if(is_negative(a)) a1=-a1;
  plist=aplist;
}

// Given a, b, lem3 returns m1 etc so that a=c1^2*m1*m12, b=c2^2*m2*m12 
// with m1, m2, m12 pairwise coprime.   At all  times these equations hold, 
// and at each step the product m1*m2*m12 is decreased by a factor d, 
// so the process terminates when the coprimality condition is satisfied. 

void rusin_lem3(const bigint& a, const bigint& b,
	  bigint& m1, bigint& m2, bigint& m3, bigint& c1, bigint& c2)
{
  m1=a; m2=b; m3=1; c1=1; c2=1;
  if((a==0)||(b==0)) return;  // shouldn't happen
  bigint d;
  int check=3;
  while(check)
    {
      //      cout<<m1<<", "<<m2<<", "<<m3<<endl;
      d=abs(gcd(m1,m2));
      if(d>1) {check=3; m1/=d; m2/=d; m3*=d;}
      else check-=1;
      if(check)
	{
	  d=abs(gcd(m1,m3));
	  if(d>1) {check=3; m1/=d; m2*=d; m3/=d; c1*=d;}
	  else check-=1;
	}
      if(check)
	{
	  d=abs(gcd(m2,m3));
	  if(d>1) {check=3; m1*=d; m2/=d; m3/=d; c2*=d;}
	  else check-=1;
	}
    }
#ifdef CHECK_LEM3
  if( (a==sqr(c1)*m1*m3) && (b==sqr(c2)*m2*m3) 
      && (gcd(m1,m3)==1) && (gcd(m2,m3)==1) && (gcd(m1,m2)==1) )
    {;}
  else
    {
      cout<<"Error in rusin_lem3("<<a<<","<<b<<"), returning\n"
	  <<"c1="<<c1<<", c2="<<c2<<", m1="<<m1<<", m2="<<m2<<", m3="
	  <<m3<<endl;
    }
#endif
}

bigint chrem(const bigint& a1, const bigint& a2, 
	     const bigint& m1, const bigint& m2)
{
  bigint u,v,q,r,ans;
  bigint g = bezout(m1,m2,u,v);
  bigint l = m1*(m2/g);
  if(::divides(a2-a1,g,q,r))
    {
      ans= (a1+u*m1*q)%l;
#ifdef CHECK_CHREM
      if(div(m1,ans-a1)&&div(m2,ans-a2)) {;}
      else cout<<"Error in chrem("<<a1<<","<<a2<<","<<m1<<","<<m2
	       <<"): returning wrong value "<<ans<<endl;
#endif
      return ans;
    }
  cout<<"No solution in chrem to "<<a1<<" mod "<<m1
      <<", "<<a2<<" mod "<<m2<<endl; 
  ans = 0;
  return ans;
}


//
// general functions
//

#ifdef MPFP
bigint Iround(bigfloat x) {return RoundToZZ(x);}
bigint Ifloor(bigfloat x) {return FloorToZZ(x);}
bigint Iceil (bigfloat x) {return CeilToZZ(x);}
#else
bigint Iceil(double x)  {return -Ifloor(-x);}
bigint Iround(double x) {return (x>0?Ifloor(x+0.5):Iceil(x-0.5));}

#define BIG 100000   // not used in new version

bigint Ifloor(double x)  // bigfloats are just doubles in this case
{
  bigint ans; ans =0;
  int s=1;
  if(x==0.0) return ans;
  if(x<0) {x=-x; s=-1;}
//#define DEBUG_IFLOOR
 int e;
 frexp(x,&e);
#ifdef DEBUG_IFLOOR
 cout<<"x="<<x<<", e="<<e<<endl;
#endif
 if(e>53)  // precision loss -- issue warning
   {
     bigint err; err=1; lshift(err,e-53,err); err-=1;
#ifdef DEBUG_IFLOOR
     cout<<"Warning in Ifloor("<<s*x<<"): possible precision loss in "
         <<"converting to bigint; maximium rounding error "
	 << "2^"<<(e-53)<<"-1 = "<<err<<endl;
#endif
   }
 while((x>0)&&(e>0))   // This was >=0 in double->Integer original version
   {
#ifdef DEBUG_IFLOOR
     cout<<"x="<<x<<", e="<<e;
#endif
     setbit(ans,e-1);
#ifdef DEBUG_IFLOOR
     cout<<", ans="<<ans;
#endif
     x=fmod(x,ldexp((double)1,e-1)); 
#ifdef DEBUG_IFLOOR
     cout<<", new x="<<x<<endl;
#endif
     frexp(x,&e);
   }
 if((x>0)&&(s<0)) ++ans;        // adjust if fractional part non-zero
 if(s<0) ans=-ans;  // adjust if negative
 return ans; 
}

#endif

bigint mod(const bigint& a, const bigint& b)
{
  bigint bb(abs(b));
  bigint c=a%bb;
  bigint c2=c<<1;
  if (c2>  bb) return c-bb; 
  if (c2<=-bb) return c+bb;
  return c;
}

long mod(const bigint& a, long b)
{
  long bb=abs(b);  
  long c = bigint_mod_long(a,bb); 
  long c2=c<<1;
  if (c2>  bb) return c-bb; 
  if (c2<=-bb) return c+bb;
  return c;
}

bigint posmod(const bigint& a, const bigint& b)
{
  bigint bb(abs(b));
  bigint c=a%bb;
  if (c<0) return c+bb;
  return c;
}

long posmod(const bigint& a, long b)
{
  long bb = abs(b);
  long c = bigint_mod_long(a,bb); 
  if (c<0) return c+bb;
  return c;
}

int divide_exact(const bigint& aa, const bigint& bb, bigint& c)
     // c = a/b with error message if remainder is non-zero
{
  bigint a(aa), b(bb), r;  // BECAUSE DIVIDE() WAS CHANGING A!!!!
  //  cout<<"In divide_exact with a = " << a << ", b = " << b << endl;
  int ok = ::divides(a,b,c,r);
  //  cout<<"After divides(),  a = " << a << ", b = " << b << ", q = " << c << ", r = " << r << endl;
  if(!ok)
    {
      cout<<"Error in dividing "<<a<<" by "<<b<<": not exact, remainder = "<<r<<endl;
    }
  return ok;
}

long divide_out(bigint& a, const bigint& d)
// divides a by d as many times as possible returning number of times (but none if a=0!)
{
  if(is_zero(a)) return 0;
  bigint q, r;
  long count=0;
  while(::divides(a,d,q,r)) {a=q; count++;}
  return count;
}

long divide_out(bigint& a, long d)
// divides a by d as many times as possible returning number of times (but none if a=0!)
{
  if(is_zero(a)) return 0;
  bigint q;
  long r, count=0;
  while(::divides(a,d,q,r)) {a=q; count++;}
  return count;
}

#define VALUATION_OF_ZERO 99999

long val(const bigint& factor, const bigint& number)
{
 if (is_zero(number)) return VALUATION_OF_ZERO;
 bigint f = abs(factor);
 if ((f<2)) return VALUATION_OF_ZERO;  // error condition! N.B. This value 
 bigint n = number;                        // must be unlikely and POSITIVE.
 long e = divide_out(n,f);
 return e;
}

long val(long factor, const bigint& number)
{
 if (is_zero(number)) return VALUATION_OF_ZERO;
 long f = abs(factor);
 if ((f<2)) return VALUATION_OF_ZERO;  // error condition! N.B. This value 
 bigint n = number;                    // must be unlikely and POSITIVE.
 long e = divide_out(n,f);
 return e;
}

vector<int> valuations(const bigint& n, const vector<bigint>& primes)
{
  vector<int> ee(primes.size());
  auto p = primes.begin();
  std::generate(ee.begin(), ee.end(),
                [n, &p](){return val(*p++, n);});
  return ee;
}

vector<int> valuations(const bigint& n, const vector<int>& primes)
{
  vector<int> ee(primes.size());
  auto p = primes.begin();
  std::generate(ee.begin(), ee.end(),
                [n, &p](){return val(*p++, n);});
  return ee;
}

int div(const bigint& factor, const bigint& number) 
{ if (is_zero(factor)) return is_zero(number);
  else return (is_zero(number%factor));
}

int div(long factor, const bigint& number) 
{ if (factor==0) return is_zero(number);
  else return is_zero(number%factor);
}

long bezout(const bigint& aa, long bb, bigint& xx, bigint& yy)
{bigint a,b,c,x,oldx,newx,y,oldy,newy,q;
 oldx = 1; oldy = 0; x = 0; y = 1; a = aa; b = bb;
 while (sign(b)!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
   newy = oldy - q*y; oldy = y; y = newy;
  }
 if (sign(a)<0) {xx=-oldx; yy=-oldy; return -I2long(a);}
 else     {xx= oldx; yy= oldy; return  I2long(a);}
}

 
bigint invmod(const bigint& a, const bigint& p)
{bigint g,x,y;
 g=bezout(a,p,x,y);
 if (!is_one(g))
   {
     x=0;
     cerr << "invmod called with " << a << " and " << p << " -- not coprime!" << endl;
   }
 return x;
}

long invmod(const bigint& a, long p)
{bigint g,x,y;
 g=bezout(a,p,x,y);
 if (!is_one(g)) 
   {
     x=0;
     cerr << "invmod called with " << a << " and " << p << " -- not coprime!" << endl;
   }
 return I2long(x);
}


int m1pow(const bigint& a)
{ return (odd(a) ?  -1 : +1);
}
 
static int table8[8] = {0,1,0,-1,0,-1,0,1};
static int table4[4] = {0,1,0,-1};
static int table44[4][4] = {{0,0,0,0}, {0,1,0,1}, {0,0,0,0}, {0,1,0,-1}};
  
int chi2(const bigint& a)
{ 
  return table8[posmod(a,8)];
}
 
int chi4(const bigint& a)
{ 
  return table4[posmod(a,4)];
} 


int hilbert2(const bigint& a, const bigint& b)
{ 
  return table44[posmod(a,4)][posmod(b,4)];
}

int hilbert2(const bigint& a, long b)
{ 
  return table44[posmod(a,4)][posmod(b,4)];
}

int hilbert2(long a, const bigint& b)
{ 
  return table44[posmod(a,4)][posmod(b,4)];
}


static int leg(const bigint& a, const bigint& b) 
//nb this function is not intended for public use!
{ 
  bigint aa = a;
  bigint bb = b;
//  cout<<"leg("<<a<<","<<b<<") = "<<flush;
  bigint c;
  int ans = 1;
  while (bb>1) 
  {     aa = aa % bb;
        if (sign(aa)<0)       {aa=-aa; ans*=chi4(bb);}
        while (is_zero(aa%4)) {aa/=4;}
        if (is_zero(aa%2))    {aa/=2; ans *= chi2(bb);}
        ans*=hilbert2(aa,bb);
        c=bb; bb=aa; aa=c;
  }
//  cout << ans << endl;
  return ans;
}
 
static int leg(const long& a, const long& b) 
//nb this function is not intended for public use!
{ long aa = a, bb = b;
//  cout<<"leg("<<a<<","<<b<<") = "<<flush;
  int ans = 1;
  while (bb>1) 
  {     aa = aa % bb;
        if (aa<0)       {aa=-aa; ans*=chi4(bb);}
        while (!(aa&3)) {aa/=4;}
        if    (!(aa&1)) {aa/=2; ans *= chi2(bb);}
        ans*=hilbert2(aa,bb);
        long c=bb; bb=aa; aa=c;
  }
//  cout << ans << endl;
  return ans;
}
 

int legendre(const bigint& a, const bigint& b)
{ 
  return ((is_one(gcd(a,b)) && (odd(b))) ? leg(a,b) : 0);
}

int legendre(const bigint& aa, long b)
{ 
  if(!(b%2)) return 0;  // b was even
  long a=I2long(aa%b);
  long g=::gcd(a,b);
  if(g!=1)  return 0;
  return leg(a,b);
}

int kronecker(const bigint& x, const bigint& y)
{ 
  long r; bigint x1=x,y1=y,z;
  int s=1;

  if (is_zero(y1)) return (abs(x1)==1);

  if (is_negative(y)) 
    {
      y1= -y1; 
      if (is_negative(x1)) s = -1;
    }

  r=divide_out(y1,2);
  if (r)
  {
    if (odd(x1))
    {
      if (odd(r) && labs(posmod(x1,8)-4) == 1) s = -s;
    }
    else return 0;
  }
  x1=posmod(x1,y1);

  while (!is_zero(x1))
  {
    r=divide_out(x1,2);
    if (odd(r) && labs(posmod(y1,8)-4) == 1) s= -s;
    if ((posmod(y1,4)==3) && (posmod(x1,4)==3)) s= -s;
    z=y1%x1; y1=x1; x1=z;
  }
  return (y1==1)? s: 0;
}
 
int kronecker(const bigint& d, long n)
{ 
  return kronecker(mod(d,n),n);
}
 
long gcd(const bigint& a, long b)
{
  bigint bb(b);
  return I2long(gcd( a, bb ));
}

// Assuming a*d-b*c!=0, computes a reduced Z-basis for <(a,b),(c,d)>
void gauss_reduce(const bigint& a0, const bigint& b0, const bigint& c0, const bigint& d0,
                  bigint& a, bigint& b, bigint& c, bigint& d)
{
  a=a0; b=b0; c=c0; d=d0;
  // cout<<"Initial (a,b) = ("<<a<<","<<b<<")"<<"; (c,d) = ("<<c<<","<<d<<")"<<endl;
  bigint P = a*a+b*b, Q = a*c+b*d, R = c*c+d*d, one(1), t;
  t = one; // any nonzero will do
  while (!is_zero(t))
    {
      // cout<<"(a,b) = ("<<a<<","<<b<<")"<<"; (c,d) = ("<<c<<","<<d<<")"<<endl;
      // cout<<"(P,Q,R) = ("<<P<<","<<Q<<","<<R<<")"<<endl;
      t = rounded_division(Q,P);
      if (!is_zero(t))
        {
          // cout<<"Shift by "<<t<<endl;
          c -= t*a;
          d -= t*b;
          Q -= t*P;
          R = c*c+d*d;
        }
      if (R<P)
        {
          // cout<<"Invert"<<endl;
          t = -a; a = c; c = t;
          t = -b; b = d; d = t;
          t = P; P = R; R = t;
          t = one;
        }
    }
  // cout<<"Final (a,b) = ("<<a<<","<<b<<")"<<"; (c,d) = ("<<c<<","<<d<<")"<<endl;
}

int modrat(const bigint& n, const bigint& m,
           /* return values: */ bigint& a, bigint& b)
{
  static const bigint zero(0), one(1);
  bigint c,d, n1 = mod(n,m);
  gauss_reduce(n1,one,m,zero,a,b,c,d);
  bigint lim = sqrt(m>>1);
  return (abs(a) <= lim) && (abs(b) <= lim);
}

// Find the number of roots of X^3 + bX^2 + cX + d = 0 (mod p) and
// assign roots to a list of these. A stupid search is good enough
// since we only use this for very small p! Also it is tacitly assumed
// that the roots are distinct.
int nrootscubic(long b, long c, long d, long p, vector<long>& roots)
{
  int nr=0;
  long r=0;
  roots.clear();
  for (r = 0; r<p; r++)
    {
      if (((((r+b)*r+c)*r+d)%p)==0)
        {
          roots.push_back(r);
          nr++;
          if (nr==2) // find 3rd root the easy way
            {
              roots.push_back(posmod(-b-roots[0]-roots[1],p));
              std::sort(roots.begin(), roots.end());
              return 3;
            }
        }
    }
  return nr;
}

void ratapprox(bigfloat x, bigint& a, bigint& b, const bigint& maxd)
{
  bigint c, x0, x1, x2, y0, y1, y2;
  bigfloat rc, xx, diff, eps = to_bigfloat(1.0e-6);
  xx = x; x0 = 0; x1 = 1; y0 = 1; y1 = 0;
  diff = 1; c=x2=y2=0;
  while (!is_approx_zero(diff)) // ( diff > eps )
    { c = Iround( xx ); rc=I2bigfloat(c);
      x2 = x0 + c*x1; x0 = x1; x1 = x2;
      y2 = y0 + c*y1; y0 = y1; y1 = y2;
      diff = abs( x - I2bigfloat(x2)/I2bigfloat(y2) );
      //      cout<<"x2 = "<<x2<<",\ty2 = "<<y2<<",\tdiff = "<<diff<<endl;
      if ( abs(xx - rc) < eps )
        diff = 0;
      else
        {
          if ( (maxd>0) && (abs(y2)>maxd) ) // go back to previous
            {
              diff = 0;
              x2 = x0;
              y2 = y0;
            }
          else
            xx = 1/(xx - rc);
        }
    }
  a = x2; b = y2;
  if ( b < 0 )
    {::negate(a); ::negate(b); }
}

void ratapprox(bigfloat x, long& a, long& b, long maxd)
{
  long c, x0, x1, x2, y0, y1, y2;
  bigfloat xx, diff, eps = to_bigfloat(1.0e-7);
  xx = x; x0 = 0; x1 = 1; y0 = 1; y1 = 0;
  diff = 1; c=x2=y2=0;
  //cout<<"ratapprox("<<x<<") with maxd="<<maxd<<endl;
  while ( diff > eps )
    { int ok = longify( xx, c, 0); // ie round(xx)
      if (!ok)
        {
          cerr<<"failed to round "<<x<<" to a long int in ratapprox"<<endl;
          return;
        }
      x2 = x0 + c*x1; x0 = x1; x1 = x2;
      y2 = y0 + c*y1; y0 = y1; y1 = y2;
      //cout<<"(x2,y2)=("<<x2<<","<<y2<<")"<<endl;
      diff = abs( x - (to_bigfloat(x2)/to_bigfloat(y2)) );
      if ( abs(xx - c) < eps )
        diff = 0;
      else
        {
          if ( (maxd>0) && (abs(y2)>maxd) ) // go back to previous
            {
              diff = 0;
              x2 = x0;
              y2 = y0;
            }
          else
            xx = 1/(xx - c);
        }
    }
  a = x2; b = y2;
  if ( b < 0 )    {a=-a; b=-b; }
  //  if ( x < 0 ) {a=-a;}
}

int is_nth_power(const bigint& x, int n)
{
  if (is_zero(x))
    return 1;
  if ((x<0) && (n%2==0))
    return 0;
  vector<bigint> plist = pdivs(x);
  auto local_condition = [n,x](const bigint& p){return val(p,x)%n==0;};
  return std::all_of(plist.begin(), plist.end(), local_condition);
}

bigint prime_to_S_part(const bigint& x,  const vector<bigint>& S)
{
  if (is_zero(x))
    return x;
  bigint y(abs(x));
  for (auto p: S)
    {
      divide_out(y, p);
      if (y==1)
        break;
    }
  return y;
}

int is_S_unit(const bigint& x,  const vector<bigint>& S)
{
  return prime_to_S_part(x, S)==1;
}

// implementation of divisor_iterator class to iterate through divisors of a factored positive integer

// constructors
divisor_iterator::divisor_iterator(const vector<bigint>& P, const vector<long>& E)
  :PP(P), EE(E)
{
  np = PP.size();
  rewind();
  nd = 1;
  std::for_each(EE.cbegin(), EE.cend(), [this](auto e){nd *= (1+e);});
}

divisor_iterator::divisor_iterator(const bigint& N)
{
  PP = pdivs(N);
  np = PP.size();
  nd = 1;
  std::for_each(PP.cbegin(), PP.cend(),
                [this, N](auto p){
                  int e = val(p,N);
                  EE.push_back(e);
                  nd *= (1+e);
                });
  rewind();
}

divisor_iterator::divisor_iterator()
  :ok(1), np(0), nd(1)
{
  NN.resize(1, bigint(1));
}

void divisor_iterator::increment()
{
  if (!ok) return;
  // find first exponent which can be incremented, increment it and
  // then reset earlier ones to 0 and update partial products
  for (int ip=0; ip<np; ip++)
    {
      if (ee[ip]<EE[ip])
        {
          ee[ip]++;
          NN[ip]*=PP[ip];
          for (int jp=0; jp<ip; jp++)
            {
              ee[jp]=0;
              NN[jp]=NN[ip];
            }
          return;
        }
    }
  // we only get here when all ee[ip]==EE[ip] so the iteration stops
  ok = 0;
}

// report on current status
void divisor_iterator::report()
{
  cout<<"Divisor iterator status:"<<endl;
  cout<<"Primes:    "<<PP<<endl;
  cout<<"Exponents: "<<EE<<endl;
  cout<<"Number of divisors: "<<nd<<endl;
  cout<<"current exponents:  "<<ee<<endl;
}

// [n^e for 0 <= e <= maxexp]
vector<bigint> powers(const bigint& n, int maxexp)
{
  vector<bigint> npowers(1+maxexp);
  bigint np(1);
  npowers[0] = np;
  int e = 0;
  auto next_power = [n, &np, &e](){np*=n; e++; return np;};
  std::generate(npowers.begin()+1, npowers.end(), next_power);
  return npowers;
}

// [n^e for e in exponents]
vector<bigint> powers(const bigint& n, const vector<int>& exponents)
{
  vector<bigint> npowers(exponents.size());
  auto npower = [n](int e){return pow(n,e);};
  std::transform(exponents.cbegin(), exponents.cend(), npowers.begin(), npower);
  return npowers;
}

// Compute N from its factorization (lists of primes and exponents) --
// (name taken from gp)
bigint factorback(const vector<bigint>&PP, const vector<int>& EE)
{
  static const bigint one(1);
  auto power = [](const bigint& p, int e){return pow(p,e);};
  return std::transform_reduce(PP.cbegin(), PP.cend(), EE.cbegin(),
                               one,                  // initial value
                               std::multiplies<>(),  // how to combine terms,
                               power);               // how to form terms
}

// Maximum conductor for a given list of primes
bigint MaxN(const vector<bigint>&PP)
{
  static const bigint one(1);
  bigint N(one);
  std::for_each(PP.cbegin(), PP.cend(),
               [&N](const bigint& p){N *= pow(p, (p==2?8:p==3?5:2));});
  return N;
}

// Radical of N
bigint radical(const bigint& N)
{
  vector<bigint> PP = pdivs(N);
  return std::accumulate(PP.begin(), PP.end(), bigint(1), std::multiplies<bigint>());
}

// multiply all integers in a list by a constant:
vector<bigint> multiply_list(const bigint& a, const vector<bigint>& L)
{
  vector<bigint> aL(L.size());
  auto times_a = [a](const bigint& x){return a*x;};
  std::transform(L.begin(), L.end(), aL.begin(), times_a);
  return aL;
}

// multiply all integers in a list by all in a second list:
vector<bigint> multiply_lists(const vector<bigint>& L1, const vector<bigint>& L2)
{
  vector<bigint> L3;
  L3.reserve(L1.size()*L2.size());
  std::for_each(L1.cbegin(), L1.cend(),
                [L2, &L3](const bigint& x)
                {
                  auto xL2 = multiply_list(x,L2);
                  L3.insert(L3.end(), xL2.begin(), xL2.end());
                });
  return L3;
}

// multiply all integers in L by p^e for e in exponents:
vector<bigint> multiply_list_by_powers(const bigint& p, const vector<int>& exponents, const vector<bigint>& L)
{
  return multiply_lists(powers(p, exponents), L);
}

// convert a list of longs to a list of bigints:
vector<bigint> bigintify(const vector<long>& L)
{
  vector<bigint> LL(L.size());
  auto long2big = [](long x){return bigint(x);};
  std::transform(L.cbegin(), L.cend(), LL.begin(), long2big);
  return LL;
}

// end of file marith.cc
