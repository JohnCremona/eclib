//
// tperiods.cc -- input a curve, find its periods
//

//#define TEST

#include "compproc.h"
#include "cperiods.h"
#include "reader.h"

int main(){
  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);
	
  CurveReader input;
  Curve E;
  Curvedata D;
  CurveRed C;

  while (input>>E)
    {
      cout << "Curve input:        " << E << endl;
      Cperiods cp(E);
      cout << "periods: " << cp << endl; 

      Curve EE = cp.trans_to_curve();
      cout << "Curve from periods: " << EE << endl;
    }
}
