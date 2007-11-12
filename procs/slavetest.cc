// slavetest.cc:  test file for class parislave

#include "interface.h"
#include "marith.h"
#include "gp.h"

int main()
{
  bigint n; n=1;
  while (n>0)
    {
      cout<<"Enter an integer: "<<flush;
      cin>>n;
      cout<<"n = "<<n<<endl;
      if(n>0) cout<<the_pari_slave.factor(n)<<endl;
    }
}
