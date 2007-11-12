//
// paricurves.cc -- read allgens and curves file, output for pari database
//

#include "interface.h"
#include "matrix.h"
#include "curve.h"
#include "points.h"
#include "../g0n/curvesort.cc"

#define CLASS_IS_LETTER // we only use letters now!

int main()
{
  //  set_precision("Enter number of decimal places");
  set_precision(100);
  initprimes("PRIMES",0);
  int verbose = 0;
  //  cout<<"verbose (0/1)? ";             cin >>verbose;
  int j, npts, npts2, ntor;

  long N, N2, ncurve, ncurve2;
  char code[20], code2[20];
  Curve E, E2;

  char filename[20],filename2[20];
  cerr<<"Enter allcurves filename : ";
  cin>>filename;
  cerr<<"Opening "<<filename<<endl;
  ifstream allcurves(filename);
  cerr<<"Enter allgens filename : ";
  cin>>filename2;
  cerr<<"Opening "<<filename2<<endl;
  ifstream allgens(filename2);

  while(!allcurves.eof()) {
  //  while(1) {

  // Input the curve's ID and the curve:

  cin >> N; if(N==0) abort();
#ifdef CLASS_IS_LETTER
  allcurves >> code;
#else 
  long nclass;
  allcurves >> nclass;
  codeletter((nclass-1),code);
#endif
  allcurves >> ncurve;
  allcurves >> E;
  Curvedata C(E);
  //  cout<<endl;
  //  cout << N<<code<<ncurve<<" = "<< E << endl;

  // Input the number of points and the points:

  Point P(C);
  allcurves >> npts >> ntor;
  vector<Point> points; points.reserve(npts);

  if(npts>0) // get points from allgens file
    {
      allgens>>N2>>code2>>ncurve2>>E2>>npts2;
      if(E==E2)
	{
	  j=0; 
	  while(j<npts)
	    { 
	      allgens >> P;
	      if ( !P.isvalid() ) 
		{
		  cout<<"point "<<P<<" not on curve.\n\n"; 
		  exit(1);
		}
	      points.push_back(P); 
	      j++;
	    }
	  //  cout<<"Input points: "<<points<<endl;
	}
      else
	{
	  cerr << "Curve mismatch!" <<endl;
	  cerr << "Curve from "<<filename<<":"<< E<<endl;
	  cerr << "Curve from "<<filename2<<": "<<E2<<endl;
	  abort(1);
	}
    }
  // Output a line suitable for the pari database:

  cout<<"[\""<<N<<code<<ncurve<<"\", "<<E<<", [";
  for(j=0; j<npts; j++) 
    {
      if(j) cout<<", ";
      output_pari(cout,points[j]);
    }
  cout<<"]]"<<endl;
  }
}


//end of file paricurves.cc

