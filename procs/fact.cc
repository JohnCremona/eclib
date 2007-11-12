// fact.cc: factorization program

#include "interface.h"
#include "marith.h"

int main(int argc, char *argv[])
{
  int verb=0;
  initprimes("PRIMES",verb);
  bigint n;

  //  cout << "Verbose? "; cin >> verb;
  if(argc>1)
    {
      if((argv[1][0]=='-')&&(argv[1][1]=='v')) 
	{
	  verb=1;
	  argv++;
	  argc--;
	}
    }

  if(argc>1) 
    n=atoI(argv[1]);
  else
    {
      cout << "\nPlease enter a number "; 
      cin >> n;
    }

  cout << "n = " << n << endl;
  vector<bigint> plist = pdivs(n,verb);
  cout<<plist<<endl << endl;
  if(verb) the_extra_primes.show();
}
