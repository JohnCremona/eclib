// tsqrt.cc: test program for integer sqrt

#include "marith.h"

  // some arrays borrowed from pari:
static int carresmod64[]={
    1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,
    0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,
    0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,
    0,0,0,0};
static int carresmod63[]={
    1,1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,
    0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,
    0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,
    0,0,0};
static int carresmod65[]={
    1,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,
    0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,
    1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,
    0,1,0,0,1};
static int carresmod11[]={1,1,0,1,1,1,0,0,0,1,0};
static int carresmod17[]={1,1,1,0,1,0,0,0,1,1,0,0,0,1,0,1,1};
static int carresmod19[]={1,1,0,0,1,1,1,1,0,1,0,1,0,0,0,0,1,1,0};
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

static int carresmod101[]={1,1,0,0,1,1,1,0,0,1,0,0,0,1,1,0,1,1,0,1,1,1,1,1,1,1,0,0,0,0,1,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,0,1,1,1,1,1,1,1,0,1,1,0,1,1,0,0,0,1,0,0,1,1,1,0,0,1};

#define bigint_mod_long(a,m) I2long(a%m)

int jsqrtnr(bigint& root, const bigint& n)
{
  //  cout<<"jsqrtnr("<<n<<")\n";
  long l = 1+(lg(n)>>1);
  bigint y;
  root=1<<l;  // first approx, > sqrt(n)
//  cout << "first approx root = " << root << endl;
  while(1)
    {
      div(n,root,y); addx(y,root,y); rshift(y,1,y);
      if(y<root) {root=y;} else break;
//      cout << "root = " << root << endl;
    }
  return (sqr(root)==n);
}
#define oddsqrt(root,n) jsqrtnr(root,n) // JC's Newton

int jsqrt(const bigint& in, bigint& root, long np)
{
//  cout<<"In isqrt with n = " << in << endl;
  root=0;
  if(sign(in)<0) return 0;
  if(sign(in)==0) return 1;
  long l,m, twopow=0;
  bigint n(in);
  while(even(n)) {n>>=1; twopow++;}
//  cout << "2-power = " << twopow << endl;
  if(odd(twopow)) return 0;            // 2 | a to an odd power
  twopow>>=1;                     // power of 2 in root
  m = bigint_mod_long(n,931170240); // 931170240=64*63*65*11*17*19 < 2^30
  if (!carresmod64[m&63]) return 0;
  if (!carresmod63[m%63]) return 0;
  if (!carresmod65[m%65]) return 0;
  if (!carresmod11[m%11]) return 0;
  if (!carresmod17[m%17]) return 0;
  if (!carresmod19[m%19]) return 0;

  /*
  m = bigint_mod_long(n,1348781387); //1348781387= 23*29*31*37*41*43 < 2^31
  if (!carresmod23[m%23]) return 0;
  if (!carresmod29[m%29]) return 0;
  if (!carresmod31[m%31]) return 0;
  if (!carresmod37[m%37]) return 0;
  if (!carresmod41[m%41]) return 0;
  if (!carresmod43[m%43]) return 0;
  */
  /*
  l=bigint_mod_long(n,101);  if (!carresmod101[l]) return 0;
  long ip=60, p, lp;
  primevar pr(np+ip);
  while(ip--) {pr++;} // skip to 17
  while(np--) 
    {
      p=pr; pr++; 
      lp=legendre(n,p);
      if(lp==-1) return 0;
    }
  */
  if(oddsqrt(root,n)) {lshift(root,twopow,root); return 1;}
  return 0;
}

int main()
{
  bigint a,r;
  long np;
  //  cout<<"np = "<<endl; cin>>np;
  int is;
  while(1)
    {
      cin>>a;
      if(a<0) break;
      //      cout<<"a = "<<a<<", lg(a) = "<<lg(a)<<endl;
      is=isqrt(a,r);
      if(is) cout<<a<<" is the square of "<<r<<"\n";
    }
  return 0;
}
