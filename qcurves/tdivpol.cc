//
// tdivpol.cc -- test for division poly functions in divpol.h/cc 
//

#include "curve.h"

#include "polys.h"
#include "divpol.h"

int main()
{
  //  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);

  Curve E(BIGINT(0),BIGINT(0),BIGINT(1),BIGINT(-7),BIGINT(6));

  Curvedata C(E);
  cout << "Curve " << E << endl;

  cout<<"Division Poly (2) = \t"<<makepdivpol(&C,2)<<endl;
  int i;
  for(i=3; i<12; i+=2)
    cout<<"Division Poly ("<<i<<") = \t"<<makepdivpol(&C,i)<<endl;


}


//end of file tdivpol.cc





