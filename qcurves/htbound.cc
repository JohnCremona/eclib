/* HTBOUND.CC:  Testing Siksek/Cremona/Prickett lower height bounds  */
/*                                */
#include "compproc.h"
#include "matrix.h"
#include "subspace.h"
#include "points.h"  // from qcurves library
#include "cperiods.h"
#include "points.h"
#include "sieve_search.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "saturate.h"
#include "mwprocs.h"
#include "htconst.h"

//#define HTB_DEBUG


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

  cout<<"Program to compute a lower bound on the height of all points with\n"; 
  cout<<"good reduction everywhere (and on the identity real component)\n";
  cout<<endl;

  long npositive=0, nnegative=0, n=0;
  
  while (getcurve())
    {
      if(C!=(Curve)CD)
	{
	  cout<<"Input curve "<<C<<", ";
	  cout<<"Working with minimal curve "<<CD<<":\t";
	}
      else 
	cout << "Curve "<< C <<":\t";

      //      int nt = CD.get_ntorsion();
      bigfloat lb = lower_height_bound(CD);
      cout << "Lower bound = "<<lb<<"\n";
      n++;
      if(sign(lb)>0) npositive++; else nnegative++;
      

      if((n%100)==0)
	{
	  cout<<"\nOut of "<<n<<" curves so far, "<<npositive
	      <<" ("<<(100.0*npositive)/n<<"%)"
	      <<" have a positive lower bound.\n";
	  cout<<endl;
	}
    }
  cout<<"\nOut of "<<n<<" curves, "<<npositive
      <<" ("<<(100.0*npositive)/n<<"%)"
      <<" have a positive lower bound.";
  cout<<endl;
}

