//
// TCURVE.CC, adapted from TestElliptic.cc by JEC
//

#include "curve.h"

int main(void)
{
  set_precision(30);
  initprimes("PRIMES",0);
        
  Curve E;
  
  cout << "\nEnter a curve: " << endl ;
  cin >> E ;
  cout << "The curve is " << E << endl ; 
  /*
    cout << "\nEnter a curve: " << endl ;
    E.input(cin);
    cout << "The curve is " << E << endl ; 
    
    cout << "To test out different constructors: Using all specified:\n" ;
    E = Curve(0,0,1,-7,6) ;
    
    cout << "the curve is " << E << "\n" ;
    cout << "Using just a4 and a6 specified:\n" ;
    E = Curve(78, 89) ;
    cout << "the curve is " << E << "\n" ; 
    
    E = Curve(0, 0, 1, -7, 6);
    cout << "the curve is " << E << "\n" ; 
  */
  cout << "A test of invariants:\n" ;
  Curvedata cd(E,0) ;  // the 0 means no minimalization
  cout << "The curve is " << cd << endl ;
  cd = Curvedata(E,1) ;  // the 1 forces minimalization
  cout << "The minimal curve is "; 
  cout << cd;
  cout << endl ;
  
  /*
    cout << "A test of extended invariants:\n" ;
    CurvedataExtra cdx(cd) ;
    cout << "The extra curve data is "; 
    cout << cdx;
    cout << endl ;
  */
  
  cout <<"A test of Tate's algorithm:\n";
  CurveRed cdr(cd);
  cout << cdr << endl;
  cout <<"Full display:\n";
  cdr.display(cout);
  
  cout <<"Traces of Frobenius:\n";
  for(primevar pr(25); pr.ok(); pr++)
    {
      long p=pr;
      cout<<"p="<<p<<": ap="<<Trace_Frob(cdr,BIGINT(p));
      if(div(p,getdiscr(cdr))) cout<<" (bad reduction)";
      cout<<endl;
    }
  return 0;
}
