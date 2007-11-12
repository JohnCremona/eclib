// tilll_lidia: test program for illl integer lll reduction

#include <LiDIA/bigint_lattice.h>

#define bigint_vector math_vector<bigint>

bigint dot(const bigint_vector& a, const bigint_vector& b);

int main()
{
  int i, j, k, n=3;
  bigint_vector* b = new bigint_vector[n];
  bigint_vector* new_b = new bigint_vector[n];
  cout<<"Enter "<<n<<" vectors: ";
  for(i=0; i<n; i++)
    {
      b[i].set_capacity(n);
      cin>>b[i];
    }

  cout<<"Before reduction, vectors are:\n";
  for(i=0; i<n; i++)
    {
      cout<<b[i]<<endl;
    }
  cout<<endl;

  cout << "METHOD: using LiDIA's bigint_lattice class\n";

  bigint_lattice M(3,3);
  bigint_matrix T(3,3);
  bigint mij;
  for(i=0; i<3; i++) 
    for(j=0; j<=i; j++)
      {
	multiply(mij,b[i],b[j]);
	M.sto(i,j,mij);
	if(j<i) M.sto(j,i,mij);
      }
  M.set_gram_flag();
  cout<<"Gram matrix = "<<M<<endl;

  M.lll(T,0.7,2);
  cout<<"After LLL, Gram matrix = "<<M<<endl;
  cout<<"transform matrix = "<<T<<endl;

  new_b[0] = T(0,0)*b[0] + T(0,1)*b[1] + T(0,2)*b[2];
  new_b[1] = T(1,0)*b[0] + T(1,1)*b[1] + T(1,2)*b[2];
  new_b[2] = T(2,0)*b[0] + T(2,1)*b[1] + T(2,2)*b[2];

  cout<<"After reduction, vectors are:\n";
  for(i=0; i<n; i++)
    {
      cout<<new_b[i]<<endl;
    }
  cout<<endl;

  delete[] b;
  delete[] new_b;
}

bigint dot(const bigint_vector& a, const bigint_vector& b)
{
  bigint ans;
  multiply(ans,a,b);
}
