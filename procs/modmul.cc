// modmul.cc:  test program for various modular multiplication methods

#include <iostream>
using namespace std;

#include "xmod.h"

//#define NTRIES 1000
#define NTRIES 10000

void all_sizes()
{
  cout << "char =        " << sizeof( char        ) << endl;
  cout << "short =       " << sizeof( short       ) << endl;
  cout << "int =         " << sizeof( int         ) << endl;
  cout << "long =        " << sizeof( long        ) << endl;
  cout << "long long =   " << sizeof( long long   ) << endl;
  cout << "float =       " << sizeof( float       ) << endl;
  cout << "double =      " << sizeof( double      ) << endl;
  cout << "long double = " << sizeof( long double ) << endl;
  cout << "int8_t =      " << sizeof( int8_t      ) << endl;
  cout << "int16_t =     " << sizeof( int16_t     ) << endl;
  cout << "int32_t =     " << sizeof( int32_t     ) << endl;
  cout << "int64_t =     " << sizeof( int64_t     ) << endl;
  cout << "uint8_t =     " << sizeof( uint8_t     ) << endl;
  cout << "uint16_t =    " << sizeof( uint16_t    ) << endl;
  cout << "uint32_t =    " << sizeof( uint32_t    ) << endl;
  cout << "uint64_t =    " << sizeof( uint64_t    ) << endl;
  cout << "size_t =      " << sizeof( size_t      ) << endl;
  cout << "void* =       " << sizeof( void*       ) << endl;
  cout << "intptr_t =    " << sizeof( intptr_t    ) << endl;
}

void int_test()
{
  int i,ntries=NTRIES;
  int b,c,ib;

  cout<<"int_test"<<endl;
  cout<<"Method: "<<XMOD_METHOD<<endl;
  cout<<"ntries: "<<ntries<<endl;

  i=ntries; 
  while(i--) 
    {
      for(b=-1000; b<1000; b++)
	{
	  if(b)
	    {
	      ib=invmod0(b);
	      c = xmod0(xmodmul0(b,ib)-1);
	      if(c!=0) 
		{
		  cout<<"Error inverting "<<b<<": invmod0 returned "<<ib<<endl;
		  cout<<"product-1 = "<<c<<" mod "<<BIGPRIME<<endl;
		  return;
		}
	      //	      else if(i==0) cout<<"Inverse of "<<b<<" is "<<ib<<"\n";
	    }
	}
    }
}

void long_test()
{
  int i,ntries=NTRIES;
  long b,c,ib;

  cout<<"long_test"<<endl;
  cout<<"Method: "<<XMOD_METHOD<<endl;
  cout<<"ntries: "<<ntries<<endl;

  i=ntries; 
  while(i--) 
    {
      for(b=-1000; b<1000; b++)
	{
	  if(b)
	    {
	      ib=invmod0(b);
	      c = xmod0(xmodmul0(b,ib)-1);
	      if(c!=0) 
		{
		  cout<<"Error inverting "<<b<<": invmod0 returned "<<ib<<endl;
		  cout<<"product-1 = "<<c<<" mod "<<BIGPRIME<<endl;
		  return;
		}
	      //	      else if(i==0) cout<<"Inverse of "<<b<<" is "<<ib<<"\n";
	    }
	}
    }
}

int main()
{
#ifdef LONG_IS_64_BIT
  cout<<"Using 64-bit arithmetic"<<endl;
#else
  cout<<"Using 32-bit arithmetic"<<endl;
#endif
  cout<<"Modulus = "<<BIGPRIME<<endl;
  cout<<"size of int = "<<sizeof(int)<<endl;
  cout<<"size of long = "<<sizeof(long)<<endl;
  cout<<"size of long long = "<<sizeof(long long)<<endl;

  all_sizes();

  //int_test();
    long_test();
}
