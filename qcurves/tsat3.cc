//
// tsat3.cc -- test for saturate.h/cc reading from gens files directly
//

#include "interface.h"
#include "matrix.h"
#include "curve.h"
#include "points.h"
#include "cperiods.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "saturate.h"
#include "elog.h"
#include "sieve_search.h"
#include "mwprocs.h"
#include "../g0n/curvesort.cc"

#define USE_EGR
#define SAT_BND 2000 // saturation bound:  use -1 for global default

#ifdef USE_EGR
const int use_egr=1;
#else
const int use_egr=0;
#endif

#define INPUT_CLASS_IS_LETTER // we only use letters now!

//void codeletter(int i, char* code, int width=0);

int main()
{
  //  set_precision("Enter number of decimal places");
  set_precision(100);
  initprimes("PRIMES",0);
  int verbose = 0;
  //  cout<<"verbose (0/1)? ";             cin >>verbose;
  int j, npts;

  long N, ncurve;
  char code[20];
  Curve E;

  while(!(feof(stdin))) {

  // Input the curve's ID and the curve:

  cin >> N; if(N==0) break;
#ifdef INPUT_CLASS_IS_LETTER
  cin >> code;
#else
  long nclass;
  cin >> nclass; 
  codeletter((nclass-1),code);
#endif
  cin >> ncurve;
  cin >> E;
  Curvedata C(E);
  cout<<endl;
  cout << N<<code<<ncurve<<" = "<< E << endl;

  // Input the number of points and the points:

  Point P(C);
  cin >> npts;
  vector<Point> points; points.reserve(npts);
  j=0; 
  while(j<npts)
    { 
      cin >> P;
      if ( !P.isvalid() ) 
	{
	  cout<<"point "<<P<<" not on curve.\n\n"; 
	  abort();
	}
      points.push_back(P); 
      j++;
    }
  cout<<"Input points: "<<points<<endl;

  bigint index;
  vector<long> unsatprimes;
  int success = saturate_points(C, points, index, unsatprimes, SAT_BND, use_egr, (verbose));

  if(success)
    {
      cout<<"Saturation complete --";
      if(index==1) 
	cout<<"input points were saturated"<<endl;
      else
	{
	  cout<<"input points had index "<<index
	      <<" in their saturation."<<endl;
	  cout<<"Basis for saturation:\t"<<points<<endl;
	}
    }
  else
    cout<<"Saturation failed at "<<unsatprimes<<endl;

  for (int i=0; i<npts; i++)
    { 
      Point P = points[i];
      cout << "Generator "<<(i+1)<<" is "<<P<<"; ";
      cout << "height "<<height(P);
      if(!P.isvalid()) cout<<" --warning: NOT on curve!";
      cout<<endl;
    }
  cout<<endl;
  cout << "Regulator = "<<regulator(points)<<endl<<endl;

  // Finally output a line similar to the input line:

#ifdef INPUT_CLASS_IS_LETTER
  cout<<N<<"\t"<<code<<"\t"<<ncurve<<"\t"<<E<<"\t"<<npts;
#else
  cout<<N<<"\t"<<nclass<<"\t"<<ncurve<<"\t"<<E<<"\t"<<npts;
#endif
  for(j=0; j<npts; j++) cout<<"\t"<<points[j];
  cout<<endl;
  cout<<"=============================================================="<<endl;
  }
}

#if(0)
void codeletter(int i, char* code, int width)
{
  int n=width;    // pads string to this width with blanks
  code[n]='\0';
  while (n) code[--n]=' ';

  int nc = i%26;
  char c = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[nc];
  n = 1 + (i-nc)/26;
  if(width==0) code[n]='\0';
  while (n) code[--n]=c;
}
#endif

//end of file tsat3.cc

