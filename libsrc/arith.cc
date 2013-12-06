// arith.cc: definitions of arithmetic functions (single precision)
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
 
#include <eclib/arith.h>


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
  register long k,size=(maxnum+257)>>1;
  if(pdiffptr) delete [] pdiffptr;
  byteptr p= new unsigned char[size+1];
  if (!p) {cout<<"Out of memory in primeclass::init!"<<endl;abort();}
  memset(p, 0, size + 1); 
  register byteptr q,r,s,fin=p+size;
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
      //      cout<<"ind="<<ind<<"\tval="<<val<<"\td="<<(int)(*aptr)<<endl;
      ok=advance();
    }
  if(!ok)
    {
	cout<<"Not enough primes in primeclass.number("<<n<<") !"<<endl;
	abort();
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
	cout<<"Not enough primes in primeclass.getfirst("<<n<<") !"<<endl;
	abort();
    }
  return ans;
}

// returns i>=0 such that p is the i'th prime
long prime_pi(long p)
{
  primevar pr;
  int ip=0;
  while ((long)pr<p) {pr++; ip++;}
  return ip;
}


long primdiv(long aa)
{
 primevar pr;
 long p=0,q;
 long a = labs(aa);
 while (pr.ok() && p==0)
   {q=pr; pr++;
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
 primevar pr; long a = abs(aa);long p;
 while ( (a>1) && (pr.ok()))
 { p = pr; pr++;
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
 long j,k,p,e,nd = 1;
 vector<long> dlist(1,1);  // cout << "Divisor 0 = 1" << endl;
 vector<long>::const_iterator pr = plist.begin();
 while(pr!=plist.end())
   {
     p=*pr++; e = val(p,a);
     dlist.resize((e+1)*dlist.size());
     for (j=0; j<e; j++)
       for (k=0; k<nd; k++)
	 dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

vector<long> alldivs(long a, const vector<long>& plist)
{//cout << "In alldivs with a = " << a << endl;
// cout << plist.size() << " primes: " <<endl; cout << plist << endl;
 long j,k,p,e,nd = 2;
 vector<long> dlist(1,1);  // cout << "Divisor 0 =  1" << endl;
 dlist.push_back(-1);      // cout << "Divisor 1 = -1" << endl;
 vector<long>::const_iterator pr = plist.begin();
 while(pr!=plist.end())
 {
   p = *pr++; e = val(p,a);
   dlist.resize((e+1)*dlist.size());
   for (j=0; j<e; j++)
     for (k=0; k<nd; k++)
       dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
   nd*=(e+1);
 }
 return dlist;
}

vector<long> sqdivs(long a, const vector<long>& plist)
{
 long j,k,p,e,nd = 1;
 vector<long> dlist(1,1);
 vector<long>::const_iterator pr = plist.begin();
 while(pr!=plist.end())
   {
     p = *pr++; e = val(p,a)/2;
     dlist.resize((e+1)*dlist.size());
     for (j=0; j<e; j++)
       for (k=0; k<nd; k++)
	 dlist[nd*(j+1)+k] =  p*dlist[nd*j+k];
     nd*=(e+1);
   }
 return dlist;
}

vector<long> sqfreedivs(long a, const vector<long>& plist)
{
 long j,k,p,e,nd = 1;
 vector<long> dlist(1,1);
 vector<long>::const_iterator pr = plist.begin();
 while(pr!=plist.end())
 {
  p = *pr++; e = 1;
  dlist.resize((e+1)*dlist.size());
  for (j=0; j<e; j++)
    for (k=0; k<nd; k++)
      dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
  nd*=(e+1);
 }
 return dlist;
}

long mod(long a, long b)
{long c;
 if (b<0) b=-b;
 if (a>=0) c=a%b; else c=b-((-a)%b);
 if (c>(b>>1)) c-=b;
 return(c);
}

long posmod(long a, long b)
{
  long c=a%b;  
  if (c<0) return(c+b);
  return(c);
}

long gcd(long a, long b)
{
  long c;
  while (b!=0) {c=a%b; a=b; b=c;}
  return abs(a);
}

int gcd(int a, int b)
{
  int c;
  while (b!=0) {c=a%b; a=b; b=c;}
  return abs(a);
}

long lcm(long a, long b)
{
  long g=gcd(a,b);
  if(g==0) return 0;
  return a*(b/g);
}

long bezout(long aa, long bb, long& xx, long& yy)
{long a,b,c,x,oldx,newx,y,oldy,newy,q;
 oldx = 1; oldy = 0; x = 0; y = 1; a = aa; b = bb;
 while (b!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
   newy = oldy - q*y; oldy = y; y = newy;
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
     cout << "invmod called with " << a << " and " << p << " -- not coprime!\n"; 
     abort();
     return 0;
   }
}

int modrat(int n, int m, float lim, int& a, int& b)
{
  long la,lb,ln=n,lm=m;
  int ans = modrat(ln,lm,lim,la,lb);
  a=la; b=lb;
  return ans;
}

//#define DEBUG_MODRAT

int modrat(long n, long m, float lim, long& a, long& b)
{long q,r,t,qq,rr,tt,quot;
#ifdef DEBUG_MODRAT
 cout<<"modrat("<<n<<","<<m<<")\n";
#endif
 q=m; r=posmod(n,m); qq=0; rr=1; t=0; tt=0; a=r; b=1; 
 if (r<lim) 
   { 
#ifdef DEBUG_MODRAT
     cout<<" = "<<a<<"/"<<b<<"\n";
#endif
     return 1;
   }
 while (r!=0) 
 { 
   quot = q/r;
#ifdef DEBUG_MODRAT
   cout<<"q,r,t = "<<q<<" "<<r<<" "<<t<<"\n";
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
       cout << "\nmodrat error: no reconstruction for " << n << " mod " << m << "\n";
#endif
       return 0;
     }
 }
#ifdef DEBUG_MODRAT
 cout << "\nmodrat error: common factor with " << n << " mod " << m << "\n";
#endif
 return 0;
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

int intbezout(int aa, int bb, int& xx, int& yy)
{int a,b,c,x,oldx,newx,y,oldy,newy,q;
 oldx = 1; oldy = 0; x = 0; y = 1; a = aa; b = bb;
 while (b!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
   newy = oldy - q*y; oldy = y; y = newy;
  }
 if (a<0) {xx=-oldx; yy=-oldy; return(-a);}
 else     {xx= oldx; yy= oldy; return( a);}
}

long chi2(long a)
{ static long table8[8] = {0,1,0,-1,0,-1,0,1};
  return table8[posmod(a,8)];
}
 
long chi4(long a)
{ static long table4[4] = {0,1,0,-1};
  return table4[posmod(a,4)];
} 

long hilbert2(long a, long b)
{ static long table44[4][4] = {{0,0,0,0},
                              {0,1,0,1},
                              {0,0,0,0},
                              {0,1,0,-1}};
  return table44[posmod(a,4)][posmod(b,4)];
}

long leg(long a, long b)  //nb this function is not intended for public use
{ long aa = a;
  long bb = b;
  long ans = 1;
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
 
long legendre(long a, long b)
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
  m=1<<e;
//                     at this point m=2^e <= n < 2^(e+1)
  if(m==n) return 1;
  if(roundup) {n=m<<1; e++;}
  return 0;
}

// stolen from pari: base_math/arith1.cc

static const int longis64bit = sizeof(long)==4;

// The following function returns valuation(z,2) for a long int z:

long val2(unsigned long z);

long kronecker(long x, long y)
{
  long r,s=1,x1,z;

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
    z=y%x1; y=x1; x1=z;
  }
  return (y==1)? s: 0;
}

long val2(unsigned long z)
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
  vector<long>plist = pdivs(n);
  for(unsigned int i=0; i<plist.size(); i++) 
    if(val(plist[i],n)>1) return 0;
  return 1;
}

int is_valid_conductor(long n)
{
  long m=n, p, e;
  e=0; while(!(m&1)) {e++; m>>=1;}  if(e>8) return 0;
  e=0; while(!(m%3)) {e++; m/=3;}   if(e>5) return 0;
  vector<long>plist = pdivs(m);
  for(unsigned int i=0; i<plist.size(); i++) 
    {
      p = plist[i]; 
      e=0; while(!(m%p)) {e++; m/=p;} if(e>2) return 0;
    }
  return 1;
}




// The following is no longer used

long old_kronecker(long d, long n)
{ long ans=0; long m=n, d4=d%4; if(d4<0)d4+=4;
  if ((gcd(d,n)==1) && ((d4==0)||(d4==1)))
  { ans=1;
    while (!(m%4)) m/=4;
    if (!(m%2)) {m/=2; ans*=((d-1)%8 ? -1 : 1);}
    ans *= legendre(d,m);
  }
  return ans;
}

/* END OF FILE */
