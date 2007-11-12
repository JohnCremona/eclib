//
// tegr.cc -- test for finding egr subgroup from a set of points
//

#include "matrix.h"
#include "cperiods.h"
#include "points.h"
#include "sieve_search.h"
#include "cperiods.h"
#include "elog.h"
#include "egr.h"
#include "htconst.h"

void codeletter(int i, char* code, int width=0);

int main()
{
  //  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);
  int j, npts;

  long N, nclass, ncurve;
  char code[20];
  Curve E;

  while(1) {
    cin >> N; if(N==0) abort();
  cin >> nclass >> ncurve;
  cin >> E;
  Curvedata C(E);
  codeletter((nclass-1),code);
  cout<<endl;
  cout<<"==============================================================="<<endl;
  cout<<endl;
  cout << N<<code<<ncurve<<" = "<< E << endl;
  Point P(C);
  cin >> npts;
  vector<Point> points; points.reserve(npts);
  j=0; 
  while(j<npts)
    { 
      cin >> P;
      if ( P.isvalid() ) {points.push_back(P); j++;}
      else {cout<<"point "<<P<<" not on curve.\n\n"; }
    }
  cout<<npts<<" points entered:"<<points<<endl;
  //  bigfloat reg = regulator(points);
  //  cout<<"Regulator = "<<reg<<endl;

  bigint tam = Tamagawa_exponent(C);
  cout<<"Tamagawa exponent = "<<tam<<endl;

  vector<Point> egr_points=points;
  bigint egri = egr_index(egr_points);
  cout<<"Index of egr subgroup = "<<egri<<endl;
  }
}

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

//end of file tsatsieve.cc
