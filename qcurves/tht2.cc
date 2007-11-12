// THT2.CC:  computing height bounds  

#include "points.h"  // from qcurves library
#include "cubic.h"
#include "sieve_search.h"
#include "htconst.h"


bigint a1,a2,a3,a4,a6;
Curve C;
Curvedata CD;
vector<Point> plist;

int getcurve(void)
{
  //  cout << "Enter curve coefficients a1,a2,a3,a4,a6: " << endl;
  cin >> C;
  //  cout<<C<<endl;
  if(C.isnull()) return 0;  // quitting condition
  //  cin >> a1 >> a2 >> a3 >> a4 >> a6;
  //  C = Curve(a1,a2,a3,a4,a6);
  CD = Curvedata(C,0);     // "1" means minimise
  //  cout<<CD<<endl;
  return 1;
}


int main()
{
  set_precision(30);

  initprimes("PRIMES",0);
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)

  long nsilverbetter=0, ncpsbetter=0, n=0;
  double silvertotal=0, cpstotal=0, cps2total=0;
  while (getcurve())
    {
      if(C!=(Curve)CD)
	{
	  cout<<"Input curve "<<C<<", ";
	  cout<<"Working with minimal curve "<<CD<<":\t";
	}
      else 
	cout << "Curve "<< C <<":\t";

      double htc = height_constant(CD);
      cout << "Height bound = " << htc << "\n";
    }
}
