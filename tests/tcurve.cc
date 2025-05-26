// TCURVE.CC, test program for curve classes
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
//

#include <eclib/curve.h>

int main(void)
{
  set_precision(30);
  initprimes("PRIMES",0);
  //the_primes.init(25000000);
  Curve E;
  
  cout << "\nEnter a curve: " << endl ;
  cin>>ws;  if(cin.eof()) {cout<<endl; exit(0);}
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
      long p = pr;
      bigint ap = Trace_Frob(cdr,bigint(p));
      cout<<"p="<<p<<": ap="<<ap;
      if(div(p,getdiscr(cdr))) cout<<" (bad reduction)";
      cout<<endl;
    }
  // for(primevar pr(1270000); pr.ok(); pr++)
  //   {
  //     p=pr;
  //     ap = Trace_Frob(cdr,bigint(p));
  //   }

  cout <<"Testing construction from a non-integral model:\n";
  bigint a1,a2,a3,a4,a6;
  E.getai(a1,a2,a3,a4,a6);
  bigrational qa1(a1),qa2(a2),qa3(a3),qa4(a4),qa6(a6);
  bigint s(60), scale;
  bigint si=s;
  qa1/=si; si*=s;
  qa2/=si; si*=s;
  qa3/=si; si*=s;
  qa4/=si; si*=s; si*=s;
  qa6/=si;
  cout<<"After scaling down by "<<s<<", coeffs are ";
  cout<<"["<<qa1<<","<<qa2<<","<<qa3<<","<<qa4<<","<<qa6<<"]"<<endl;
  vector<bigrational> qai;
  qai.push_back(qa1);
  qai.push_back(qa2);
  qai.push_back(qa3);
  qai.push_back(qa4);
  qai.push_back(qa6);
  Curvedata Es(qai,scale);
  cout<<"Constructed curve is "<<(Curve)Es<<" with scale = "<<scale<<endl;
  return 0;
}
