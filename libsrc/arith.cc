// arith.cc: definitions of arithmetic functions (single precision)
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
 
#include "eclib/arith.h"

/* Prime number procs; adapted from Pari  */

// These have not been stl-ized at all since they work just fine...

primeclass the_primes;    // The one and only instance

primeclass::primeclass() 
{
  pdiffptr=0; // will be allocated in init()
  ifstream pfile("MAXPRIME");
  if(!pfile)
    {
      init(1000000);  // default value
    } 
  else
    {
      long maxnum;
      pfile>>maxnum;
      init(maxnum);
    }
}

primeclass::primeclass(long maxnum) 
{
  pdiffptr=0; // will be allocated in init()
  init(maxnum);
}

void primeclass::init(long maxnum)  /* initializes variable pdiffptr */
                                    /* to give primes up to maxnum  */
{
  long k,size=(maxnum+257)>>1;
  if(pdiffptr) delete [] pdiffptr;
  byteptr p= new unsigned char[size+1];
  if (!p) {cerr<<"Out of memory in primeclass::init!"<<endl; return;}
  memset(p, 0, size + 1); 
  byteptr q,r,s,fin=p+size;
  for(r=q=p,k=1;r<fin;)
    {
      do {r+=k; k+=2; r+=k;} while (*++q);
      for(s=r;s<fin;s+=k) *s=1;
    }
  r=p; *r++=2; *r++=1;  /* 2 and 3 */
  for(s=q=r-1;; s=q)
    {
      do q++; while (*q);
      if (q>=fin) break;
      *r++=(q-s)<<1;
      s=q;
    }
  *r++=0;
  NPRIMES=r-p-1;
  BIGGESTPRIME=((s - p) << 1) + 1;
  //    cout<<"Near end of init, NPRIMES = "<<NPRIMES<<endl;
  //    cout<<"Significant elements of pdiffptr: ";
  //    for(k=0; k<NPRIMES+1; k++) cout<<(int)p[k]<<" "; cout<<endl;
//  cout<<"BIGGESTPRIME = "<< BIGGESTPRIME << endl;
  pdiffptr = new unsigned char[NPRIMES+1];
  q=p; r=pdiffptr; k=NPRIMES+1;
  while(k--) {*r = *q; r++; q++;}
  delete [] p;
  reset();
//  cout<<"At end of init, NPRIMES = "<<NPRIMES<<endl;
  //    cout<<"First few elements of pdiffptr: ";
  //    for(k=0; k<10; k++) cout<<(int)pdiffptr[k]<<" "; cout<<endl;
}

primeclass::~primeclass()
{
  delete [] pdiffptr;
}

void primeclass::reset(void) {p_ind=0; p_val=0; p_aptr=pdiffptr;}
int primeclass::at_end(void) {return *p_aptr==0;}
int primeclass::advance(void)
{
    unsigned char d=*p_aptr;
    if(d) {p_ind++; p_val+=d; p_aptr++; return 1;}
    else  {return 0;}
}

long primeclass::number(long n)
// returns n'th prime from list, starting at n=1 for p=2
{
//  cout << "In primeclass::number("<<n<<")"<<endl;
  if(n<p_ind) reset();
  int ok=1;
//  cout << "Advancing to the "<<n<<"th prime...\n";
  while((p_ind<n)&&ok)
    {
      ok=advance();
    }
  if(!ok)
    {
      cerr<<"Not enough primes in primeclass.number("<<n<<") !"<<endl;
    }
  return p_val;
}

vector<long> primeclass::getfirst (long n)  /* returns list of first n primes */
{
//  cout << "In primeclass::getfirst("<<n<<")"<<endl;
  vector<long> ans;
  reset();
  int ok=1;
  for (long i=0; (i<n)&&ok; i++)
    {
      ok=advance();
      ans.push_back(p_val);
    }
  if(!ok)
    {
      cerr<<"Not enough primes in primeclass.getfirst("<<n<<") !"<<endl;
    }
  return ans;
}

// returns i>=0 such that p is the i'th prime
long prime_pi(long p)
{
  primevar pr;
  int ip=0;
  while ((long)pr<p) {++pr; ip++;}
  return ip;
}


long primdiv(long aa)
{
 primevar pr;
 long p=0;
 long a = labs(aa);
 while (pr.ok() && p==0)
   {
     long q=pr; ++pr;
     if (a%q==0) p = q;
     else if (q*q>a) p=a;   // N.B. this causes a=1 to return 1.  Beware!
   }
 if (p==0) {p=a;
            cout<<"No prime divisor found for "<<aa<<" so assuming prime!\n";
           }
 return p;
}


vector<long> pdivs(long aa)
{vector<long> plist;
 primevar pr;
 long a = abs(aa);
 while ( (a>1) && (pr.ok()))
 { long p = pr; ++pr;
 if (a%p==0) 
   {
     plist.push_back(p);
     while (a%p==0) a/=p;      //divide out by all the p's in a
   }
 else if (p*p>a) 
   {
     plist.push_back(a); a=1;
   }
 }
 if (a>1) {plist.push_back(a);}  //In case of p-factors outside range, assume
                                //the cofactor is prime.
 return plist;
}


vector<long> posdivs(long a, const vector<long>& plist)
{
// cout << "In posdivs with a = " << a << endl;
// cout << plist.size() << " primes: " <<endl; cout << plist << endl;
 long nd = 1;
 vector<long> dlist(1,1);  // cout << "Divisor 0 = 1" << endl;
 for (const auto& p : plist)
   {
     long e = val(p,a);
     dlist.resize((e+1)*dlist.size());
     for (long j=0; j<e; j++)
       for (long k=0; k<nd; k++)
	 dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   }
 return dlist;
}

vector<long> alldivs(long a, const vector<long>& plist)
{//cout << "In alldivs with a = " << a << endl;
// cout << plist.size() << " primes: " <<endl; cout << plist << endl;
 long nd = 2;
 vector<long> dlist(1,1);  // cout << "Divisor 0 =  1" << endl;
 dlist.push_back(-1);      // cout << "Divisor 1 = -1" << endl;
 for (const auto& p : plist)
 {
   long e = val(p,a);
   dlist.resize((e+1)*dlist.size());
   for (long j=0; j<e; j++)
     for (long k=0; k<nd; k++)
       dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
   nd*=(e+1);
 }
 return dlist;
}

vector<long> sqdivs(long a, const vector<long>& plist)
{
 long nd = 1;
 vector<long> dlist(1,1);
 for (const auto& p : plist)
   {
     long e = val(p,a)/2;
     dlist.resize((e+1)*dlist.size());
     for (long j=0; j<e; j++)
       for (long k=0; k<nd; k++)
	 dlist[nd*(j+1)+k] =  p*dlist[nd*j+k];
     nd*=(e+1);
   }
 return dlist;
}

vector<long> sqfreedivs(long a, const vector<long>& plist)
{
 long nd = 1;
 vector<long> dlist(1,1);
 for (const auto& p : plist)
 {
  long e = 1;
  dlist.resize((e+1)*dlist.size());
  for (long j=0; j<e; j++)
    for (long k=0; k<nd; k++)
      dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
  nd*=(e+1);
 }
 return dlist;
}

// gcc division truncates towards 0, while we need rounding, with a
// consistent behaviour for halves (they go up here).
//
// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
long rounded_division(long a, long b)
{
  std::ldiv_t qr = ldiv(a, b);
  long r = qr.rem, q = qr.quot;
  long r2 = r<<1;
  return (r2<-b? q-1: (r2>=b? q+1: q));
}

long mod(long a, long b)
{long c;
 if (b<0) b=-b;
 if (a>=0) c=a%b; else c=b-((-a)%b);
 if (c>(b>>1)) c-=b;
 return(c);
}

long mod(int a, long b)
{long c;
 if (b<0) b=-b;
 if (a>=0) c=a%b; else c=b-((-a)%b);
 if (c>(b>>1)) c-=b;
 return(c);
}

int mod(int a, int b)
{int c;
 if (b<0) b=-b;
 if (a>=0) c=a%b; else c=b-((-a)%b);
 if (c>(b>>1)) c-=b;
 return(c);
}

int mod(long a, int b)
{
  return (int)mod(a,(long)b);
}

long posmod(long a, long b)
{
  long c=a%b;
  if (c<0) return(c+b);
  return(c);
}

long posmod(int a, long b)
{
  long c=(long(a))%b;
  if (c<0) return(c+b);
  return(c);
}

int posmod(int a, int b)
{
  int c=a%b;
  if (c<0) return(c+b);
  return(c);
}

int posmod(long a, int b)
{
  return (int)posmod(a, (long)b);
}

long gcd(long a, long b)
{
  if ((a==1)||(b==1)) return 1;
  if (a==0) return abs(b);
  while (b!=0) {long c=a%b; a=b; b=c;}
  return abs(a);
}

int gcd(int a, int b)
{
  if ((a==1)||(b==1)) return 1;
  if (a==0) return abs(b);
  while (b!=0) {int c=a%b; a=b; b=c;}
  return abs(a);
}

long lcm(long a, long b)
{
  long g=gcd(a,b);
  if(g==0) return 0;
  return a*(b/g);
}

long bezout(long aa, long bb, long& xx, long& yy)
{long a = aa, b = bb, x = 0, oldx = 1, y = 1, oldy = 0;
 while (b!=0)
 { long q = a/b;
   long c    = a    - q*b; a    = b; b = c;
   long newx = oldx - q*x; oldx = x; x = newx;
   long newy = oldy - q*y; oldy = y; y = newy;
  }
 if (a<0) {xx=-oldx; yy=-oldy; return(-a);}
 else     {xx= oldx; yy= oldy; return( a);}
}

long invmod(long a, long p)
{long g,x,y;
 g=bezout(a,p,x,y);
 if (g==1) return x;
 else
   {
     cout << "invmod called with " << a << " and " << p << " -- not coprime!"<<endl;
     return 0;
   }
}

//#define DEBUG_MODRAT

// Assuming a*d-b*c!=0, computes a reduced Z-basis for <(a,b),(c,d)>
void gauss_reduce(long a0, long b0, long c0, long d0,
                  long& a, long& b, long& c, long& d)
{
  a=a0; b=b0; c=c0; d=d0;
#ifdef DEBUG_MODRAT
  cout<<"Initial (a,b) = ("<<a<<","<<b<<")"<<"; (c,d) = ("<<c<<","<<d<<")"<<endl;
#endif
  long P = a*a+b*b, Q = a*c+b*d, R = c*c+d*d, t=1;
  while (t)
    {
#ifdef DEBUG_MODRAT
      cout<<"(a,b) = ("<<a<<","<<b<<")"<<"; (c,d) = ("<<c<<","<<d<<")"<<endl;
      cout<<"(P,Q,R) = ("<<P<<","<<Q<<","<<R<<")"<<endl;
#endif
      t = rounded_division(Q,P);
      if (t)
        {
#ifdef DEBUG_MODRAT
          cout<<"Shift by "<<t<<endl;
#endif
          c -= t*a;
          d -= t*b;
          Q -= t*P;
          R = c*c+d*d;
        }
      if (R<P)
        {
#ifdef DEBUG_MODRAT
          cout<<"Invert"<<endl;
#endif
          t = -a; a = c; c = t;
          t = -b; b = d; d = t;
          t = P; P = R; R = t; Q=-Q;
          t = 1;
        }
    }
#ifdef DEBUG_MODRAT
  cout<<"Final (a,b) = ("<<a<<","<<b<<")"<<"; (c,d) = ("<<c<<","<<d<<")"<<endl;
#endif
}

// Set a, b so that a/b=n (mod m) with |a|, |b| minimal; return success if a^2, b^2 <= m/2
int modrat(int n, int m, int& a, int& b)
{
  long la,lb,ln=n,lm=m;
  int ok = modrat(ln,lm,la,lb);
  a=la; b=lb;
  return ok;
}

int old_modrat(long n, long m, long& a, long& b);
int new_modrat(long n, long m, long& a, long& b);

int modrat(long n, long m, long& a, long& b)
{
  return old_modrat(n, m, a, b);
  //return new_modrat(n, m, a, b);  // NB new version has problems on 32-bit
}

int old_modrat(long n, long m, long& a, long& b)
{long q=m, r=posmod(n,m), qq=0, rr=1;
#ifdef DEBUG_MODRAT
 cout<<"modrat("<<n<<","<<m<<")\n";
#endif
 float lim = sqrt(float(m)/2.0);
 a=r; b=1;
 if (r<lim)
   {
#ifdef DEBUG_MODRAT
     cout<<" = "<<a<<"/"<<b<<"\n";
#endif
     return 1;
   }
 while (r)
 {
   long quot = q/r, t, tt;
#ifdef DEBUG_MODRAT
   cout<<"q,r,qq,rr = "<<q<<" "<<r<<" "<<qq<<" "<<rr<<"\n";
#endif
   t  =  q-quot*r;   q = r;   r = t;
   tt = qq-quot*rr; qq = rr; rr = tt;
   if (r<lim)
     {
       if (abs(rr)<lim)
         {
           a=r; b=rr;
#ifdef DEBUG_MODRAT
           cout<<" success:  "<<a<<"/"<<b<<"\n";
#endif
           return 1;
         }
#ifdef DEBUG_MODRAT
       cerr << "***modrat failure: no reconstruction for " << n << " mod " << m << "\n";
#endif
       return 0;
     }
 }
 cerr << "***modrat error: common factor with " << n << " mod " << m << "\n";
 return 0;
}

int new_modrat(long n, long m, long& a, long& b)
{
#ifdef DEBUG_MODRAT
  cout<<"modrat("<<n<<","<<m<<")\n";
#endif
  long c,d, n1 = mod(n,m);
  gauss_reduce(n1,1,m,0,a,b,c,d);
#ifdef DEBUG_MODRAT
  cout<<" = "<<a<<"/"<<b<<"\n";
#endif
 float lim = sqrt(float(m)/2.0);
 return (abs(a) <= lim) && (abs(b) <= lim);
}

long val(long factor, long number)
{
 long n = abs(number), f = abs(factor);
 if ((n==0) || (f<2)) return 99999;  // error condition! N.B. This value 
                                     // must be unlikely and POSITIVE.
 long e = 0;
 while (n%f==0) {e++; n/=f;}
 return e;
}

int bezout(int aa, int bb, int& xx, int& yy)
{int a = aa, b = bb, x = 0, oldx = 1, y = 1, oldy = 0;
 while (b!=0)
 { long q = a/b;
   long c    = a    - q*b; a    = b; b = c;
   long newx = oldx - q*x; oldx = x; x = newx;
   long newy = oldy - q*y; oldy = y; y = newy;
  }
 if (a<0) {xx=-oldx; yy=-oldy; return(-a);}
 else     {xx= oldx; yy= oldy; return( a);}
}

int chi2(long a)
{ static const int table8[8] = {0,1,0,-1,0,-1,0,1};
  return table8[posmod(a,8)];
}

// set root to rounded sqrt(a) if a>=0, return 1 iff exact
int isqrt(long a, long& root)
{
  if (a<0) {return 0;}
  root = round(sqrt(a));
  return a==root*root;
}

// return rounded sqrt(a) (undefined for a<0)
long isqrt(const long a)
{
  long r;
  isqrt(a,r);
  return r;
}

long squarefree_part(long d)
{
  if (d==0) return d;
  long maxd = sqdivs(d).back();
  return (d/maxd)/maxd;
}

int chi4(long a)
{ static const int table4[4] = {0,1,0,-1};
  return table4[posmod(a,4)];
}

int hilbert2(long a, long b)
{ static int table44[4][4] = {{0,0,0,0},
                              {0,1,0,1},
                              {0,0,0,0},
                              {0,1,0,-1}};
  return table44[posmod(a,4)][posmod(b,4)];
}

int leg(long a, long b)  //nb this function is just a helper for legendre()
{ long aa = a;
  long bb = b;
  int ans = 1;
  while (bb>1)
  {     aa = aa % bb;
        if (aa<0)       {aa=-aa; ans*=chi4(bb);}
        while (!(aa%4)) {aa/=4;}
        if (!(aa%2))    {aa/=2; ans *= chi2(bb);}
        ans*=hilbert2(aa,bb);
        long c=bb; bb=aa; aa=c;
  }
  return ans;
}

int legendre(long a, long b)
{
  return (((gcd(a,b)==1) && (b%2)) ? leg(a,b) : 0);
}

// Function which returns 1 and sets e such that 2**e=n if n is a power of 2.
// If the "roundup" flag is set and n is not a power of 2 it increases n to
// the next power of 2 (and returns 0)

int intlog2(long& n, long& e, int roundup)
{ 
  e = 0;
  if (n<1) {if(roundup) n=1; return 0;}
  long m=n;
  while (m) { m >>= 1; e++; }
  e--;
  m=long(1)<<e;
//                     at this point m=2^e <= n < 2^(e+1)
  if(m==n) return 1;
  if(roundup) {n=m<<1; e++;}
  return 0;
}

// stolen from pari: base_math/arith1.cc

static const int longis64bit = sizeof(long)==4;

// The following function returns valuation(z,2) for a long int z:

int val2(unsigned long z);

int kronecker(long x, long y)
{
  long r,s=1,x1;

  if (y<=0)
  {
    if (y) { y= -y; if (x<0) s = -1; }
    else  return (labs(x)==1);
  }
  r=val2(y);              // = valuation(y,2)
  if (r)
  {
    if (odd(x))
    {
      if (odd(r) && labs((x&7)-4) == 1) s = -s;
      y>>=r;
    }
    else return 0;
  }
  x1=x%y; if (x1<0) x1+=y;
  while (x1)
  {
    r=val2(x1);
    if (r)
    {
      if (odd(r) && labs((y&7)-4) == 1) s= -s;
      x1>>=r;
    }
    if (y&2 && x1&2) s= -s;
    long z=y%x1; y=x1; x1=z;
  }
  return (y==1)? s: 0;
}

int val2(unsigned long z)
{
  int v=0;
  while (!(z&1)) {v++; z>>=1;}
  return v;
}

int is_squarefree(long n)
{
  if(n%4==0) return 0;
  if(n%9==0) return 0;
  if(n%25==0) return 0;
  if(n%49==0) return 0;
  auto plist = pdivs(n);
  return std::all_of(plist.begin(), plist.end(), [n] (const long& p) {return val(p,n)==1;});
}

int is_valid_conductor(long n)
{
  if (n<11) return 0;
  long m=n, e;
  e=0; while(!(m&1)) {e++; m>>=1;}  if(e>8) return 0;
  e=0; while(!(m%3)) {e++; m/=3;}   if(e>5) return 0;
  auto plist = pdivs(m);
  return std::all_of(plist.begin(), plist.end(), [m] (const long& p) {return val(p,m)<=2;});
}


// a=b*q+r, return 1 iff r==0
int divrem(long a, long b, long& q, long& r)
{
  std::ldiv_t qr = ldiv(a, b);
  r = qr.rem;
  q = qr.quot;
  return (r==0);
}

// a=b*q+r, return 1 iff r==0
int divrem(int a, int b, int& q, int& r)
{
  std::div_t qr = div(a, b);
  r = qr.rem;
  q = qr.quot;
  return (r==0);
}

// return list of integers from first to last inclusive
vector<long> range(long first, long last)
{
  vector<long> ans(last-first+1);
  std::iota(ans.begin(), ans.end(), first);
  return ans;
}

/* END OF FILE */
