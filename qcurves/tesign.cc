// TESIGN.CC: test program for esign.cc, computing the sign of the
// functional equation over Q.

// Designed to work with same data file input as tmrank, so each input
// line should have the 5 coefficients followed by the rank (mod 2)
                                                  
#include "curve.h"
#include "esign.h"

bigint a1,a2,a3,a4,a6;
Curve C;
Curvedata CD;
int verbose;

int getcurve(void)
{
  //  if (verbose) cout << "Enter curve coefficients a1,a2,a3,a4,a6 ?" << "\n";
  //  cin >> a1 >> a2 >> a3 >> a4 >> a6;
  //  C = Curve(a1,a2,a3,a4,a6);
  cin>>C;
  CD = Curvedata(C,1);     // "1" means minimise
  return !C.isnull();
  //  return (a1!=0||a2!=0||a3!=0||a4!=0||a6!=0);
}

int main()
{
  set_precision(20);
  initprimes("PRIMES",verbose);
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)

  int allok=1, verbose=1, ok, sfe, filerank;
  long count=0, badcount=0;
  cout<<"verbose? "; cin>>verbose;

  while (getcurve())
    {
      count++;
      cin >> filerank;
      if(verbose) cout << "Curve "<< C << " :\t"<<flush;
      sfe = GlobalRootNumber(CD);
      ok = (sfe==m1pow(filerank));
      
      if (ok)
	{
	  if(verbose) 
	    cout<<"sign = "<<sfe<<", rank = "<<filerank<<": OK!"<<endl;
	  else
	    {
	      cout<<"."<<flush;
	      if(count%50==0)cout<<" "<<count<<endl;
	    }
	}
      else
	{
	  allok=0; badcount++;
	  if(!verbose) cout << "Curve "<< C << " :\t";
	  if(sfe)
	    cout << "Wrong! rank of "<<C<<" is " << filerank
		 << " but sign =  "<<sfe<<endl;
	  else 
	    cout << "computed sign =  0"<<endl;
	    
	}
    }
  if(allok) cout<<"All "<<count<<" curves ok"<<endl;
  else cout<<badcount<<" curves out of "<<count<<" were wrong"<<endl;
}


