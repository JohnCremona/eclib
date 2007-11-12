//
// tsat2.cc -- test for saturate.h/cc reading from gens files directly
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

#define PMIN 2
#define PMAX 100

void codeletter(int i, char* code, int width=0);


int main()
{
  set_precision(100);
  initprimes("PRIMES",0);
  int verbose = 1;
  //  cout<<"verbose (0/1)? ";             cin >>verbose;
  int i, j, npts;

  long N, nclass, ncurve;
  char code[20];
  Curve E;

  long curvecount=0;
  long okcount=0;
  long upcount=0;
  vector<vector<long> > keeplist;  // list of curves which were not saturated, for report at end

  while(1) {
    cin >> N; if(N==0) break;
  cin >> nclass >> ncurve;
  cin >> E;
  codeletter((nclass-1),code);
  curvecount++;

  Curvedata C(E);
  saturator sieve(&C,verbose);

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
      else {cout<<"point "<<P<<" not on curve.\n\n"; abort();}
    }
  cout<<npts<<" points entered:"<<points<<endl;

  int pmax=PMAX;

  sieve.set_points(points);
  int index = sieve.do_saturation_upto(pmax);

  cout<<"Finished p-saturation for p up to "<<pmax;
  if(index>1) 
    {
      cout<<", index gain = "<<index<<endl;
      vector<Point> newpoints = sieve.getgens();
      cout<<"New generators:\n"<<newpoints<<endl;
      upcount++;
      vector<long> keep;
      keep.push_back(N);
      keep.push_back(nclass);
      keep.push_back(ncurve);
      keep.push_back(index);
      keeplist.push_back(keep);
    }
  else
    {
      cout<<", points were saturated"<<endl;
      okcount++;
    }
  }

  cout<<endl;
  cout<<"==============================================================="<<endl;
  cout<<endl;
  cout<<"Number of curves entered:    "<<curvecount<<endl;
  cout<<"Number saturated already:    "<<okcount<<endl;
  cout<<"Number which were saturated: "<<upcount<<endl;
  //  keeper.close();
  if(upcount>0)
    {
      cout<<"Curves which needed saturation: "<<endl;
      for(i=0; i<upcount; i++)
	{
	  vector<long> keep = keeplist[i];
	  codeletter((keep[1]-1),code);
	  cout<<keep[0]<<code<<keep[2]<<": index gained = "<<keep[3]<<endl;
	}
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

//end of file tsat2.cc





