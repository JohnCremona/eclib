#include <LiDIA/Fp_polynomial.h>
#include <LiDIA/factorization.h>

using namespace LiDIA;
using namespace std;

factorization<Fp_polynomial> fact_c(const bigint& p, long deg, bigint *c);

int main()
{
  bigint p;
  cout<<"Enter modulus p: "; cin>>p;
  cout<<"p = "<<p<<endl;

  long i,d;
  cout<<"Enter degree d: "; cin>>d;

  bigint * coeffs = new bigint[d+1];
  for(i=0; i<=d; i++)
    {
      cout<<"Enter coefficient of x^"<<i<<": "; cin>>coeffs[i];
    }
  
  factorization<Fp_polynomial> factors = fact_c(p,d,coeffs);

  cout<<"The factorization is:\n"<<factors<<endl;


}


factorization<Fp_polynomial> fact_c(const bigint& p, long deg, bigint *c)
{
  Fp_polynomial f;
  f.set_modulus(p);
  f.set_max_degree(deg);
  for (long i=0; i<=deg; i++) { f[i].assign(c[i]); }
  f.make_monic();
  cout<<" About to factor "<<f<<endl;
  cout<<"Leading coeff = "<<f.lead_coeff()<<endl;
  cout<<"Const Term    = "<<f.const_term()<<endl;
  cout<<"Modulus       = "<<f.modulus()<<endl;
  cout<<"Degree        = "<<f.degree()<<endl;
  factorization<Fp_polynomial> fact_f=factor(f);
  return fact_f;
}

