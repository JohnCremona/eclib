// factor.cc: factorization program

#include <iostream>
using namespace std;

#include "LiDIA/bigint.h"
#include <LiDIA/rational_factorization.h>
#include <LiDIA/timer.h>
using namespace LiDIA;

void fact(const bigint& n, int verb);

int main(int argc, char *argv[])
{
  bigint n;

  int verb=0;

  if(argc>1) 
    string_to_bigint(argv[1],n);
  else
    {
      cout << "\nPlease enter a number "; 
      cin >> n;
    }
  cout << "n = " << n << endl;
  fact(n,verb);
  cout<<endl;
}

void fact(const bigint& n, int verb)
{      
  timer T;
  T.set_print_mode(HMS_MODE);
  
  rational_factorization f;
  f.assign(n);
  f.verbose(verb);
  
  T.start_timer();
  f.factor();
  T.stop_timer();
  
  cout << "Factorization : " << f ;
  
  if (f.is_prime_factorization())
    cout<<"  is prime factorization\n";
  else 
    cout<<"  is not a prime factorization\n";
  cout << "cpu time = " << T << "\n";
}
