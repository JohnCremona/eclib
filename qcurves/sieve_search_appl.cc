// File        : sieve_search_appl.c
// Author      : Sophie Labour
// Last change : SL, Jun 30 1999, first definitions
// Last change : JC, Dec 2 1999, adapted for direct use in mwrank et al.

#include "interface.h"
#include "sieve_search.h"

int main()
{
  point_printer printer;
  point_counter counter;
  int deg,verbose,pflag;
  vector<bigint> coeff;
  double hlim;
  lidia_size_t i;

  cout<<"Height limit?"<<endl; cin>>hlim;
  cout<<"Verbose? (0/1)"<<endl; cin>>verbose;
  cout<<"Count points (0) or print them (1)?"<<endl; cin>>pflag;
  point_processor* processor = (pflag? (point_processor*)&printer : (point_processor*)&counter);

  while (1)
    {
      cout<<"Degree?"<<endl; cin>>deg;
      if (!deg) abort();
      coeff.resize(deg+1);
      cout<<"Coeffs?"<<endl;
      for (i=0;i<=deg;i++)
	cin>>coeff[i];
      cout<<"Degree = "<<deg<<"\nCoeffs = "<<coeff<<endl;
      qsieve s(processor,deg,coeff,hlim,verbose);
      s.search();
      if(!pflag) cout << counter.get_tally() << " points found";
      cout << endl;
    }
}
