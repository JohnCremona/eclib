// tround.cc: Test program for Iround etc

#include "marith.h"

int main()
{
  cout.precision(15);
  double x;
/*
  double d = 100000;
  x =-3.0063103472586036e+28;
  cout << "x = " << x << endl;
  double ax = fabs(x);
  cout << "ax = fabs(x) = " << ax << endl;
  double hx = floor(ax/d);
  cout << "hx = floor(ax/d) = " << hx << endl;
  double lx = ax-hx*d;
  cout << "ax-hx*"<<d<<" = " << lx << endl;
  double flx = fabs(lx);
  cout << "fabs("") = " << flx << endl;
*/
  cout << "Enter a double x: ";  cin >> x;

  bigint a = Iround(x);
  cout << "Iround(x) = " << a << "\n";
  bigint b = Ifloor(x);
  cout << "Ifloor(x) = " << b << "\n";
  bigint c = Iceil(x);
  cout << "Iceil(x)  = " << c << "\n";

  // Test of roundover.
  // With b>0, c=roundover(a,b) is the nearest int to a/b
  b=5;
  for(a=-b; a<=3*b; a+=1)
    {
      c=roundover(a,b);
      cout<<"roundover("<<a<<","<<b<<") = "<<c<<endl;
    }

}
