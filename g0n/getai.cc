// FILE geti.cc: redundant

#include "marith.h"
#include "getai.h"

int c4c6ok(Integer c4,Integer c6)
{
  Integer disc = c4*c4*c4 - c6*c6;
  if (disc==0) return 0;          //
  if (disc%1728!=0) return 0;     // need c4^3-c6^2=1728D, with D|=0
  long x6= mod(c6,27);
  if((x6==9)||(x6==-9)) return 0; // need c6 != +-9 (mod 27)
  if((c6+1)%4==0) return 1;       // OK if c6 = -1 (mod 4)
  if(c4%16!=0) return 0;          // else need c4=0 (mod 16)
  x6=mod(c6,32);                  //      and
  return ((x6==0) || (x6==8));    //      c6 = 0,8 (mod 32).
}

IntegerArray getai(Integer c4,Integer c6)
{
  IntegerArray a(5);
  a[0]=0; a[1]=0; a[2]=0; a[3]=0; a[4]=0;
  if(!c4c6ok(c4,c6)) {cerr<<"Bad c4, c6: "<<c4<<", "<<c6<<endl; return a;}
  Integer b2 = mod(-c6,12);
  Integer b22 = b2*b2;
  Integer b4 = (b22-c4)/24;
  Integer b6 = (-b2*b22+36*b2*b4-c6)/216;
  a[0] /* =a1 */ = (odd(b2) ? 1 : 0);
  a[2] /* =a3 */ = (odd(b6) ? 1 : 0);
  a[1] /* =a2 */ = (b2-a[0]*a[0])/4;
  a[3] /* =a4 */ = (b4-a[0]*a[2])/2;
  a[4] /* =a6 */ = (b6-a[2]*a[2])/4;
  return a;
}
