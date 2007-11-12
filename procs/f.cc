// f.cc: factorization using LiDIA

#include <LiDIA/rational_factorization.h>

using namespace LiDIA;

int main()
{
  rational_factorization f;
  bigint n;
  cout<<"Enter positive integer to be factored: "; cin>>n;
  f.assign(n);
  f.factor();
  if(f.is_prime_factorization())
    cout<<"Prime factorization = "<<f<<endl;
  else
    cout<<"Factorization = "<<f<<endl;
  return 0;
}

