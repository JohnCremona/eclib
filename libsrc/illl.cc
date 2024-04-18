// illl.cc: implementations of functions for integer LLL
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
// 
// This file is part of the eclib package.
// 
// eclib is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// eclib is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with eclib; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////

#include <eclib/illl.h>

//#define DEBUG_LLL
//#define TRACE_LLL

bigint sdot(const vector<vec_m>& b, int i, int j)
{
  bigint ans(0);
  const vec_m& g=b[0], bi=b[i], bj=b[j];
  int n=dim(g);
  for(int k=1; k<=n; k++)
    ans += (g[k]*bi[k]*bj[k]);
  return ans;
}

void show(const int n, const vector<vec_m>& b, const vector<vector<bigint>>& lambda, const vector<bigint>& d)
{
  int i=0, j;
  cout<<"Vectors:\n";
  for (const auto& bi : b)
    {
      if (i) cout<<bi<<endl;
      i++;
    }
  cout<<endl;
  cout<<"d: ";
  for ( const auto& di : d) cout<<di<<"\t";
  cout<<endl;
  cout<<"Lambda matrix:\n";
  for(i=1; i<=n; i++)
    {
      for(j=1; j<=i; j++)
	if(j==i) cout<<d[i]<<"\t";
	else cout<<lambda[i-1][j-1]<<"\t";
      cout<<endl;
    }
  cout<<endl;
}

void redi(const int n, const int k, const int l,
          vector<vec_m>& b, vector<vector<bigint>>& lambda, const vector<bigint>& d)
{
#ifdef TRACE_LLL
  cout<<"R"<<k;
#endif
#ifdef DEBUG_LLL
  cout<<"In redi with k="<<k<<", l="<<l<<endl;
  show(n,b,lambda,d);
#endif
  int i;
  bigint lkl=lambda[k-1][l-1], dl=d[l], q;
  nearest(q,lkl,dl); // nearest integer to lkl/dl
  if(is_zero(q))
    return;
  b[k]-=q*b[l];
  lambda[k-1][l-1]-=q*dl;
  for(i=1; i<l; i++) lambda[k-1][i-1]-=q*lambda[l-1][i-1];
#ifdef DEBUG_LLL
  cout<<"Leaving redi with k="<<k<<", l="<<l<<endl;
  show(n,b,lambda,d);
#endif
}

void swapi(const int n, const int k, const int kmax,
	   vector<vec_m>& b, vector<vector<bigint>>& lambda, vector<bigint>& d)
{
#ifdef TRACE_LLL
  cout<<"S"<<k<<endl;
#endif
#ifdef DEBUG_LLL
  cout<<"In swapi with k="<<k<<", kmax="<<kmax<<endl;
  show(n,b,lambda,d);
#endif
  bigint t, lam, bb, dk=d[k], dk1=d[k-1];
  int i, j;
  std::swap(b[k-1],b[k]);
  for(j=1; j<=k-2; j++)
    {
      t = lambda[k-1][j-1];
      lambda[k-1][j-1]=lambda[k-2][j-1];
      lambda[k-2][j-1]=t;
    }
  lam = lambda[k-1][k-2];
  bb   = (d[k-2]*dk+sqr(lam))/dk1;
  for(i=k+1; i<=kmax; i++)
    {
      t = lambda[i-1][k-1];
      lambda[i-1][k-1] = (dk*lambda[i-1][k-2]-lam*t)/dk1;
      lambda[i-1][k-2] = (bb*t+lam*lambda[i-1][k-1])/dk;
    }
  d[k-1] = bb;
#ifdef DEBUG_LLL
  cout<<"Leaving swapi with k="<<k<<", kmax="<<kmax<<endl;
  show(n,b,lambda,d);
#endif
}


void step3(const int n, int& k, const int kmax,
	   vector<vec_m>& b, vector<vector<bigint>>& lambda, vector<bigint>& d)
{
  redi(n,k,k-1,b,lambda,d);
  bigint lhs = 4*(d[k]*d[k-2]+sqr(lambda[k-1][k-2]));
  bigint rhs = 3*sqr(d[k-1]);
#ifdef DEBUG_LLL
  cout<<"lhs="<<lhs<<", rhs="<<rhs<<endl;
#endif
  if( lhs<rhs )
    {
      swapi(n,k,kmax,b,lambda,d);
      k--; if(k<2) k=2;
      step3(n,k,kmax,b,lambda,d);
    }
  else
    {
      int l;
      for(l=k-2; l>0; l--) redi(n,k,l,b,lambda,d);
      k++;
    }
}

// b is an array of n+1 vectors indexed from 0 to n.
// b[1]...b[n] are the lattice basis, while b[0] holds the coefficients
// of the (diagonal) Gram matrix, so the inner product of b[i] and b[j]
// is sum(k,b[0][k]*b[i][k]*b[j][k]).
//

void lll_reduce(const int n, vector<vec_m>& b)
{
  int i, j, k, kmax;
  bigint u;
  vector<bigint> d(n+1);
  vector<vector<bigint>> lambda(n, vector<bigint>(n));

  k=2; kmax=1;
  d[0]=1; d[1]=sdot(b,1,1);

  while(k<=n)
    {
      vector<bigint>& lambda_k = lambda[k-1];
      if(k>kmax)
	{
	  kmax=k;
	  for(j=1; j<=k; j++)
	    {
	      vector<bigint>& lambda_j = lambda[j-1];
	      u=sdot(b,k,j);
	      for(i=1; i<j; i++)
		{
		  u=(d[i]*u-lambda_k[i-1]*lambda_j[i-1])/d[i-1];
		}
	      if(j<k) lambda_k[j-1]=u;
	      else
		{
		  if(u==0)
		    {
		      cout<<"lll_reduce(): input vectors dependent!\n";
		      return;
		    }
		  d[k]=u;
		}
	    }
	}
#ifdef DEBUG_LLL
	  cout<<"After step 2 with k="<<k<<", kmax="<<kmax<<endl;
	  show(n,b,lambda,d);
#endif
      step3(n,k,kmax,b,lambda,d);
    }
#ifdef TRACE_LLL
  cout<<endl;
#endif
}
