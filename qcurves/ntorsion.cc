//
// ntorsion.cc -- count torsion points 

#include "points.h"
#include "reader.h"

int main(){
  cerr<<"Program to count torsion on a curve.\n\n";
  set_precision("Enter number of decimal places");
  initprimes("PRIMES",1);
  CurveReader input;
  Curve F;
  Curvedata E(F);
  int ntor;
  while (input>>F) 
    {
      E = Curvedata(F,1);
      cout<<"Curve "<<F;
      ntor = ntorsion(E);
      cout<<" \t has " << ntor << " torsion point(s)" << endl;
    }
} //ends main

	


