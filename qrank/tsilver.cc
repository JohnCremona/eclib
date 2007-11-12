/* TSILVER.CC:  Test program for Silverman & Siksek height bounds  */
/*                                */
#include "points.h"  // from qcurves library
#include "silver.h"


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
  initprimes("PRIMES",0);
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)

  long nsilverbetter=0, nsiksekbetter=0, n=0;
  
  while (getcurve())
    {
      if(C!=(Curve)CD)
	{
	  cout<<"Input curve "<<C<<", ";
	  cout<<"Working with minimal curve "<<CD<<":\t";
	}
      else 
	cout << "Curve "<< C <<":\t";

      double silver = silverman_bound(CD);
      cout << "Silverman bound = " << silver << ",\t";
      double siksek = siksek_bound(CD);
      cout << "Siksek bound = " << siksek << "\n";
      if (silver<siksek) nsilverbetter++; else nsiksekbetter++;
      n++;
      if(n%100==0) 
	cout<<"So far, siksek is better for "<<nsiksekbetter
	    <<" curves out of "<<n<<endl;
//       double best = height_constant(CD);
//       cout << "Best = " << best << "\n";
      
    }
  cout<<"Overall, siksek is better for "<<nsiksekbetter
      <<" curves out of "<<n<<endl;
  double rate = (100.0*nsiksekbetter)/n;
  cout << "For " << rate << "% of curves the new bound is better."<<endl;
}
