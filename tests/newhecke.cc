// NEWHECKE.CC  -- Hecke operators: factored characteristic polynomials on the new cuspidal subspace

#include <eclib/convert.h>
#include <eclib/timer.h>
#include <eclib/homspace.h>
#include <eclib/polys.h>

#define AUTOLOOP

int main(void)
{
 cout << "Program newhecke: test of new and/or cuspidal char polys of Hecke operators." << endl;
  scalar modulus = default_modulus<scalar>();

  int show_only_new_pols=1;
  cerr << "See only the new char polys (0/1)? ";
  cin >> show_only_new_pols;

  int show_only_cuspidal_pols=1;
  cerr << "See only the cuspidal char polys (0/1)? ";
  cin >> show_only_cuspidal_pols;

  int np;
  cerr << "How many T_p? ";
  cin >> np;

  long N;
#ifdef AUTOLOOP
  int firstn, lastn;
  cerr<<"Enter first and last levels: ";
  cin >> firstn >> lastn;
  firstn = max(firstn,2);
  for (N=firstn; N<=lastn; N++)
    {
#else
  while(cerr<<"Enter level: ", cin>>N, N>1)
    {
#endif
  cout << ">>>Level " << N << endl;
  const homspace* h = get_homspace(N, modulus);
  int dim = h->h1dim();
  int cdim = h->h1cuspdim();
  cout << "Dimension = " << dim << endl
       << "Cuspidal dimension = " << cdim
       << endl << endl;

  ZZX charpol;

  if (dim==0) continue; // to next level

  int ntp = 0;
  for (primevar pr(np); pr.ok()&&ntp<np; pr++)
    {
      long p = pr;
      if (divides(p,N))
        continue;
      ntp++;
      matop T(p, N);
      cout << "p = " << p << " T = "<<T.name()<<endl;
      if (!show_only_new_pols)
        {
          if (!show_only_cuspidal_pols)
            {
              charpol = get_poly(N, T, 0, modulus); // cuspidal=0, triv_char=0
              cout << "Full characteristic polynomial of " << T.name() << ": "
                   << str(charpol)
                   << endl;
              if (deg(charpol)>0)
                {
                  cout <<"Factors:"<<endl;
                  display_factors(charpol);
                  cout << endl;
                }
            }
          charpol = get_poly(N, T, 1, modulus); // cuspidal=1, triv_char=0
          cout << "Full cuspidal characteristic polynomial of " << T.name() << ": "
               << str(charpol)
               << endl;
          if (deg(charpol)>0)
            {
              cout <<"Factors:"<<endl;
              display_factors(charpol);
              cout << endl;
            }
        }
      charpol = get_new_poly(N, T, 1, modulus); // cuspidal=1, triv_char=0
      int dimnewcusp = deg(charpol);
      if (ntp==1)
        {
          cout << "Cuspidal newspace has dimension " << dimnewcusp << endl;
        }
      cout << "New cuspidal characteristic polynomial of " << T.name() << ": "
           << str(charpol)
           << endl;
      if (deg(charpol)>0)
        {
          cout << "Factors:" <<endl;
          display_factors(charpol);
          cout << endl;
        }
      cout << endl;
    }
}       // end of while()
exit(0);
}       // end of main()
