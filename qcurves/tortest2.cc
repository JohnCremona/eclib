//
// tortest2.cc -- to test analytic torsion routines
//

#include "points.h"
#include "elog.h"
#include <fstream>

#define TORSION_ONLY // else computes and checks conductors too.

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
    //    cout<<endl<<endl;
    Curvedata E(C);
    //    cout<<"Curve "<<(Curve)(E)<<" has input torsion order "<<T<<endl;

    //    nt = E.get_ntorsion();
    //cout<<"Curve "<<(Curve)(E)<<" has torsion order "<<nt<<endl;
    vector<Point>Tlist = torsion_points(E,T);
    nt=Tlist.size();
    //    cout<<"Computed "<<nt<<" torsion points "<<Tlist<<endl;
    if (T != nt)
    {cout<<"\n Bad ntorsion for "<<C<<": ";
     cout<<"T as read: "<<T<<"  T as calculated: "<<nt<<endl;
     cout<<Tlist<<endl;
     nbadnt++;
    }
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
