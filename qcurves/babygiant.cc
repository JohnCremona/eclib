//
// babygiant.cc -- test for my_isomorphism_type() and related functions
//

#include "curve.h"
#include "points.h"
#include "polys.h"
#include "gf.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "ffmod.h"

#if defined(LiDIA_INTS) || defined(LiDIA_ALL)
#define NextPrime next_prime
#endif

int main()
{

  //  Curve E(BIGINT(0),BIGINT(0),BIGINT(1),BIGINT(-7),BIGINT(6));
  //  Curve E(BIGINT(0),BIGINT(-1),BIGINT(1),BIGINT(0),BIGINT(0));
  Curve E(BIGINT(1),BIGINT(0),BIGINT(0),atoI("-7630810"),atoI("8112786291"));
  Curvedata C(E);
  cout << "Curve " << E << endl;
  bigint disc = getdiscr(C);
  long iq = NextPrime(20050);
  bigint q = BIGINT(iq);
  bigint qmin, qmax;
  cout<<"Enter qmin and qmax: "; cin>>qmin>>qmax;
  q = NextPrime(qmin-1);
  while(q<=qmax) {
  while(disc%q==0)   q=NextPrime(q+1);
  cout<<"q = "<<q<<endl;

  galois_field Fq(q);
  curvemodq Cq = reduce_curve(C,q);
  cout<<"Curve mod q = "<<Cq<<endl;
  cout<<"Group order = "<<Cq.group_order()<<endl;

  pointmodq P1, P2;
  bigint n1, n2, n;
  my_isomorphism_type(Cq,n1,n2,P1,P2);
  //  Cq.isomorphism_type(n1,n2,P1,P2);
  n=n1*n2;
  cout<<"Group structure of "<<(Cq)<<" mod "<<q<<": \n";
  cout<<" gen 1 = "<<P1<<" (order "<<n1<<")\n";
  cout<<"Check: order is "<<order_point(P1)<<endl;
  if(n2>1) 
    {
      cout<<" gen 2 = "<<P2<<flush<<" (order "<<n2<<")"<<endl;
      cout<<"Check: order is "<<order_point(P2)<<endl;
    }

  cout<<endl;
  q=NextPrime(q+1);
  }
}


