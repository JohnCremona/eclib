#include <iostream>
#include "log2.h"

using namespace std;

int main()
{
  long n,e,rup;
  while(1)
    {
      cout<<"Enter an integer n: ";
      cin >> n;
      if(n<0) break;
      if(intlog2(n,e,0))
	cout<<n<<" = 2^"<<e<<endl;
      else
	{
	  cout<<n<<" is not a power of 2\n";
	  intlog2(n,e,1);
	  cout<<"rounds up to "<<n<<" = 2^"<<e<<endl;	  
	}
    }
}
