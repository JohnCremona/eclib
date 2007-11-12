// test.cc: matrix test program

#include <LiDIA/bigint_matrix.h>
#include <iostream>
#include "interface.h"
#include "method.h"
#include "mmatrix.h"

int main()
{
  matrix m;
  m.read_from_file("mat");
  int i,j,n = nrows(m);
  cout<<"Read matrix of size "<<n<<"*"<<ncols(m)<<" from file mat"<<endl;
  matrix m2 = m.slice(1,10,1,10);
  cout<<"Top left 10x10 block:"<<m2<<endl;
  
  long biggest=0;
  for(i=1; i<=n; i++) 
    for(j=1; j<=n; j++) 
      if(abs(m(i,j))>biggest) 
	{biggest=m(i,j);}
  cout<<"Largest entry = "<<biggest<<endl;

  m2 = m;
  if(n>100) 
    {
      n=100;
      m2 = m.slice(1,n,1,n);
      cout<<"From now on we work with the top left 100x100 block only"<<endl;
    }

#ifdef LiDIA_INTS
  bigint_matrix M(n,n);
  for(i=1; i<=n; i++) 
    for(j=1; j<=n; j++) 
      M.sto(i-1,j-1,m2(i,j));
  cout<<"Created LiDIA copy of M"<<endl;
  bigint d = M.det();
  cout<<"det(M) = "<<d<<endl;
#else
#endif

//    cout<<"Computing kernel..."<<flush;
//    subspace ker = kernel(m2);
//    cout<<"done.  Dimension = "<<dim(ker)<<endl;

  int a;
  bigint A;
  
  while(0) {
  cout<<"Enter an integer: "; cin>>A;
  cout<<"A = "<<A<<endl;
  if(is_int(A))
    {
      a = I2int(A);
      cout<<" as a int, A = "<<a<<endl;
    }
  else
    {
      cout<<"A does not fit into a int"<<endl;
    }
  }

  cout<<"m2.slice(1,10,1,10)=\n"<<m2.slice(1,10,1,10)<<endl;
  mmatrix mm(n,n);
  for(i=1; i<=n; i++) 
    for(j=1; j<=n; j++) 
      mm(i,j)=m2(i,j);
  cout<<"mm.slice(1,10,1,10)=\n"<<mm.slice(1,10,1,10)<<endl;
  
  SCALAR dummy;
  cout<<"MININT = "<<MININT<<endl;
  cout<<"MAXINT = "<<MAXINT<<endl;
  matrix m3 = mm.shorten(dummy);
  cout<<"m3.slice(1,10,1,10)=\n"<<m3.slice(1,10,1,10)<<endl;
  
  int ok = (m2==m3);
  cout<<"After lengthening and shortening, matrices ";
  if(!ok) cout<<" DO NOT ";
  cout<<"agree"<<endl;

}


