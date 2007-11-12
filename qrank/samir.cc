#include <LiDIA/Fp_polynomial.h>
#include <LiDIA/factorization.h>
#include "mquartic.h"
#include "mlocsol.h"
#include "samir.h"

/* Samir's Local Solubility Test FOR ODD PRIMES*/

int new_qpsoluble(const quartic& g, const bigint& p, int verbose)
{
  bigint a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
  return new_qpsoluble(a,b,c,d,e,p,verbose);
}

int new_qpsoluble_ace(const bigint& a, const bigint& c, const bigint& e, 
		      const bigint& p, int verbose)
{
  bigint b(0);
  return new_qpsoluble(a,b,c,b,e,p,verbose);
}

int new_qpsoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d, 
			const bigint& e, const bigint& p, int verbose)
{
  if(p<5000) 
    {
      if(verbose) 
	cout<<"new_qpsoluble with p<100 passing to old qpsoluble.\n";
      return qpsoluble(a,b,c,d,e,p);
    }
  if (new_zpsol(a,b,c,d,e,p,verbose)) return 1;
  else return new_zpsol(e,d,c,b,a,p,verbose);
} /* end of new_qpsoluble */

int new_zpsol(const bigint& a,const bigint& b,const bigint& c,const bigint& d,const bigint& e, const bigint& p, int verbose)
{
  bigint *coeff =new bigint[5];
  coeff[0]=a;  coeff[1]=b;  coeff[2]=c;  coeff[3]=d;  coeff[4]=e;
  int res=local_sol(p,coeff,verbose);
  delete[] coeff;
  return res;
}

factorization<Fp_polynomial> fact_c(const bigint& p,bigint *c)
{
  Fp_polynomial f;
  f.set_modulus(p);
  for (long i=0; i<5; i++) { f.set_coefficient(c[i],i); }
  factorization<Fp_polynomial> fact_f=factor(f);
  return fact_f;
}

/* Samir's Local Solubility Test for odd p */
int local_sol(const bigint& p,bigint *c, int verbose)
{

  if (verbose)
    {  cout << "---------------------------------------------\n";
       cout << "LOCAL_SOL \n";
       cout << c[4] << " " << c[3] << " " << c[2] << " " << c[1] << " ";
       cout << c[0] << "      p=" << p << "\n" << flush;
    }

  bigint r[2],t;
  long e,fl,i;
  bigint p2=p*p;
  int xxx=1;
  for (i=0; i<5 && xxx; i++)
     { xxx=(c[i]%p).is_zero(); }
  single_factor<Fp_polynomial> term;
  Fp_polynomial F;
  if (xxx)
    { // Step II  (s is zero on entry)
      if (verbose) { cout << "Step II\n"; }
      xxx=1;
      for (i=0; i<5 && xxx; i++) { xxx=(c[i]%p2).is_zero(); }
      bigint *dd=new bigint[5];
      if (xxx)
         { for (i=0; i<5; i++)   { dd[i]=c[i]/p2; }
           int res=local_sol(p,dd);
           delete[] dd;
           return res;
         }
      for (i=0; i<5; i++)        { dd[i]=c[i]/p; }
      factorization<Fp_polynomial> fact_f=fact_c(p,dd);
      // Is there a non-repeated root
      if (verbose) { cout << fact_f << "\n" << flush; }
      for (i=0; i<fact_f.no_of_components(); i++)
        { e=fact_f.prime_exponent(i);
          F=fact_f.prime_base(i).base();
          if (F.degree()==1 && e==1)
             { if (verbose) { cout << "Non-Repeated Root\n " << flush; }
               delete[] dd;
               return 1;
             }
        }

      /* Go thru each repeated root and make the
         required transformation
      */
      bigint *d=new bigint[5];
      fl=0;
      for (i=0; i<fact_f.no_of_components() && !fl; i++)
        { term=fact_f.prime_base(i);
          F=term.base();
          e=fact_f.prime_exponent(i);
          if (F.degree()==1 && e!=1)
            { bigint r=-F[0];
              if (verbose) { cout << "Repeated Root=" << r << "\n" << flush; }
              d[4]=dd[4]*p2*p;
              d[3]=p2*(dd[3]+4*dd[4]*r);
              d[2]=p*(dd[2]+6*dd[4]*r*r+3*dd[3]*r);
              d[1]=(2*dd[2]*r+4*dd[4]*r*r*r+3*dd[3]*r*r+dd[1]);
              d[0]=((((dd[4]*r+dd[3])*r+dd[2])*r+dd[1])*r+dd[0])/p;
              fl=local_sol(p,d);
            }
        }
      delete[] d;
      delete [] dd;
      return fl;
    }
// STEP I
  bigint un;
  xxx=1;
  for (i=4; i>=0 && xxx; i--)
     { un=c[i];
       xxx=(c[i]%p).is_zero();
     }
  // If leading non-zero term is a square return 1
  if (jacobi(un,p)==1)      { return 1; }
  // If f is a constant mod p and constant not a square return 0
  if (i==-1) { return 0; }
  // Factorize f

  factorization<Fp_polynomial> fact_f=fact_c(p,c);
  long nc=fact_f.no_of_components();
  // Check if of the form un g^2
  if (verbose) { cout << fact_f << "\n" << flush; }
  for (i=0; i<nc; i++)
    { if (fact_f.prime_exponent(i)%2!=0) { return 1; } }
  if (verbose) { cout << "Roots mod p\n"; }
  // Compute roots mod p
  long num_r=0;
  for (i=0; i<nc; i++)
     { term=fact_f.prime_base(i);
       F=term.base();
       if (F.degree()==1)
         { r[num_r]=-F[0];
           if (verbose) { cout << r[num_r] << "\n" << flush; }
           num_r=num_r+1;
         }
     }
  if (num_r==0) { return 0; }

  // Now need to determine g
  if (nc==1)
    { term=fact_f.prime_base(0);
      F=term.base();
      if (fact_f.prime_exponent(0)%4==0) { F=F*F; }
    }
  else
    { term=fact_f.prime_base(0);
      F=term.base();
      term=fact_f.prime_base(1);
      F=F*term.base();
    }
  bigint te,g0,g1,g2;
  g2=F[0];
  g1=F[1];
  g0=0;
  if (F.degree()==2) { g0=1; }
  if (verbose) { cout << "g = " << g0 << " " << g1 << " " << g2 << "\n" <<
flush; }
  // Now determine h
  bigint h0,h1,h2,h3,h4;
  h4=(c[4]-un*g0*g0)/p;
  h3=(c[3]-2*un*g1*g0)/p;
  h2=(c[2]-2*un*g2*g0-un*g1*g1)/p;
  h1=(c[1]-2*un*g1*g2)/p;
  h0=(c[0]-un*g2*g2)/p;
  if (verbose)
    { cout << "h =" << h4 << " " << h3 << " " << h2 << " " << h1 << " ";
      cout << h0 << "\n" << flush;
    }
  /* For each root which is also a root of h
     transform the equations and call again
  */
  bigint *d=new bigint[5];
  fl=0;
  for (i=0; i<num_r && !fl; i++)
    { te=(((h4*r[i]+h3)*r[i]+h2)*r[i])%p;
      te=((te+h1)*r[i]+h0)%p;
      if (te.is_zero())
        { if (verbose) { cout << "Using " << r[i] << "\n" << flush; }
          d[4]=c[4]*p2;
          d[3]=p*(c[3]+4*c[4]*r[i]);
          d[2]=(c[2]+6*c[4]*r[i]*r[i]+3*c[3]*r[i]);
          d[1]=(2*c[2]*r[i]+4*c[4]*r[i]*r[i]*r[i]+3*c[3]*r[i]*r[i]+c[1])/p;
          d[0]=((((c[4]*r[i]+c[3])*r[i]+c[2])*r[i]+c[1])*r[i]+c[0])/(p2);
          fl=local_sol(p,d);
        }
    }

  delete[] d;
  return fl;
}
