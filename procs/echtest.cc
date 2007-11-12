// echtest.cc: Echelon test program

#include <iostream.h>
#include <time.h>
#include <builtin.h>
#include "matrix.h"

int main(void)
{
  long starttime,stoptime,cputime;
  time(&starttime);
  cout << "\nEchelon test program.\n\n";
  
  long nr,nc,rk, ny,denom;
  cout << "Enter number of rows and number of columns: ";
  cin>>nr>>nc;
  matrix a(nr,nc);
  cout << "Enter entries of matrix: ";
  cin>>a;

  vector pc,npc;
  long method;
  cout << "Which echelon method? (0=standard,2=modular) ";
  cin>>method;
  cout << "\nUsing method " << method;
  if(method==2) cout << " (modulus = " << BIGPRIME << ")";
  cout << endl;
  matrix ref = echelon(a, pc, npc, rk, ny, denom, method);
  cout << "Echelon matrix = " << ref;
  cout << "pivotal columns: " << pc << endl;
  cout << "nonpivotal columns: " << npc << endl;
  cout << "Denom = " << denom << endl;

  ny = nr-rk;
  cout << "Rank = " << rk << endl;
  cout << "Nullity = " << ny << endl;

  time(&stoptime);
  cout << "cpu time = " << (stoptime-starttime) << " seconds\n";
}
