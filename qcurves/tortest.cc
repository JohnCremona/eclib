//
// tortest.cc -- to test torsion routines
//

#include "points.h"
#include <fstream>

#define TORSION_ONLY // else computes and checks conductors too.

//#define TEST_3_TORSION

int main(){
  initprimes("PRIMES",0);
  cout<<"Enter filename: ";
  char filename[20]; cin>>filename; ifstream in(filename);
  in.flags( in.flags() | ios::dec );
  bigint a1,a2,a3,a4,a6;
  long T, nt; bigint N, nn;
  long ncurves=0,nbadnt=0;
#ifndef TORSION_ONLY    
  long nbadnn=0;
#endif

  cout<<"Starting, output . every 100 curves\n";
  int more=1;
  while (more) {
    Curve C;   
#ifdef TORSION_ONLY    
    in >> C >> T;
#else
    in >> C >> N >> T;
    CurveRed CR(C);
#endif
    more=!C.isnull();
    if(more){
    ncurves++;
    if((ncurves%100)==0) cout<<"."<<flush;
    Curvedata E(C);

    vector<Point>Tlist = torsion_points(E);
    nt=Tlist.size();
    //cout<<"Curve "<<(Curve)(E)<<" has torsion order "<<nt<<endl;
    if (T != nt)
    {cout<<"\n Bad ntorsion for "<<C<<": ";
     cout<<"T as read: "<<T<<"  T as calculated: "<<nt<<endl;
     nbadnt++;
    }
#ifdef TEST_3_TORSION
    vector<Point>T3list = three_torsion(E);
    unsigned int nt3=T3list.size();
    vector<Point>T3list2;
    for(unsigned int j=0; j<Tlist.size(); j++)
      if((3*Tlist[j]).iszero()) T3list2.push_back(Tlist[j]);
    if(nt3!=T3list.size())
      {
	cout<<"Error:  three_torsion gives "<<nt3<<" points "<<endl;
      }
    //    else cout<<"3-torsion ok"<<endl;
#endif
#ifndef TORSION_ONLY
    nn = getconductor(CR);
    if (N != nn)
    {cout<<"\n Bad conductor for "<<C<<": ";
     cout<<"N as read: "<<N<<"  N as calculated: "<<nn<<endl;
     nbadnn++;
    }
#endif
  }
  }
  cout<<"\n\n"<<ncurves<<" curves tested\n";
#ifndef TORSION_ONLY
  if(nbadnn) cout << nbadnn << " wrong conductors\n";
  else cout << "All conductors correct\n";
#endif
  if(nbadnt) cout << nbadnt << " wrong torsion orders\n";
  else cout << "All torsion orders correct\n";

  in.close();
} //ends main
