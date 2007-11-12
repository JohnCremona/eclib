// illl.cc: implementations of functions for integer LLL
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include "illl.h"

//#define DEBUG_LLL
//#define TRACE_LLL
#define DEBUG_LIST_SHORT_VECS

void show(const int n, const vec_m* b, const bigint** lambda, const bigint* d);

void redi(const int n, const int k, const int l, 
	  vec_m* b, bigint** lambda, const bigint* d);

void swapi(const int n, const int k, const int kmax, 
	   vec_m* b, bigint** lambda, bigint* d);

void step3(const int n, int& k, const int kmax, 
	   vec_m* b, bigint** lambda, bigint* d);

// b is an array of n+1 vectors indexed from 0 to n.  
// b[1]...b[n] are the lattice basis, while b[0] holds the coefficients 
// of the (diagonal) Gram matrix, so the inner product of b[i] and b[j] 
// is sum(k,b[0][k]b[i][k]*b[j][k]). 
//

void lll_reduce(const int n, vec_m* b)
{
  int i, j, k, kmax;
  bigint u;
  bigint* d = new bigint[n+1];
  bigint ** lambda = new bigint*[n];
  for(i=0; i<n; i++) lambda[i] = new bigint[n];

  k=2; kmax=1;
  d[0]=1; d[1]=sdot(b,1,1);

  while(k<=n)
    {
      bigint* lambda_k = lambda[k-1];
      if(k>kmax)
	{
	  kmax=k;
	  for(j=1; j<=k; j++)
	    {
	      bigint* lambda_j = lambda[j-1];
	      u=sdot(b,k,j);
	      for(i=1; i<j; i++)
		{
		  u=(d[i]*u-lambda_k[i-1]*lambda_j[i-1])/d[i-1];
//		  divide_exact(d[i]*u-lambda_k[i-1]*lambda_j[i-1],d[i-1],u);
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
  for(i=0; i<n; i++) delete [] lambda[i];
  delete [] lambda;
  delete [] d;
#ifdef TRACE_LLL
  cout<<endl;
#endif
}

void step3(const int n, int& k, const int kmax, 
	   vec_m* b, bigint** lambda, bigint* d)
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

void redi(const int n, const int k, const int l, 
	  vec_m* b, bigint** lambda, const bigint* d)
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
  if(is_zero(q)) return;
  b[k]-=q*b[l];
  lambda[k-1][l-1]-=q*dl;
  for(i=1; i<l; i++) lambda[k-1][i-1]-=q*lambda[l-1][i-1];
#ifdef DEBUG_LLL
  cout<<"Leaving redi with k="<<k<<", l="<<l<<endl;
  show(n,b,lambda,d);
#endif
}

void swapi(const int n, const int k, const int kmax, 
	   vec_m* b, bigint** lambda, bigint* d)
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
  swapvec(b[k-1],b[k]);
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

void show(const int n, const vec_m* b, const bigint** lambda, const bigint* d)
{
  int i, j;
  cout<<"Vectors:\n";
  for(i=1; i<=n; i++) cout<<b[i]<<endl;
  cout<<endl;
  cout<<"d: ";
  for(i=0; i<=n; i++)  cout<<d[i]<<"\t";
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

bigint sdot(const vec_m* b, int i, int j)
{
  bigint ans;
  const vec_m& g=b[0];
  const vec_m& bi=b[i];
  const vec_m& bj=b[j];
  int n=dim(g), k;
  for(k=1; k<=n; k++) ans += (g[k]*bi[k]*bj[k]);
  return ans;
}

//
// Uses Pohst-Zassenhaus Algorithm (page 190) to find all vectors of
// length < c, where the quadratic form is again given by b[0].
//
// NB The following DOES NOT WORK: I blindly implemented P-Z without
// doing the necessary preliminary completing of the square.
//
#if(0)
int dig(const int n, const vec_m* b, const int i, vec_m& x, 
	const mat_m& Q, vec_m& T, vec_m& U);


void list_short_vecs(const int n, vec_m* b, const bigint& c)
{
  mat_m Q(n,n);
  vec_m T(n), U(n), x(n);
  int i,j;
  for(i=1; i<=n; i++) for(j=1; j<=n; j++) Q(i,j)=sdot(b,i,j);
  T[n]=c; U[n]=0;
  dig(n, b, n, x, Q, T, U);
}

vec_m comb(const int n, const vec_m* b, const vec_m& x);

int dig(const int n, const vec_m* b, const int i, vec_m& x, 
	const mat_m& Q, vec_m& T, vec_m& U)
{
#ifdef DEBUG_LIST_SHORT_VECS
  cout<<"dig("<<i<<")"<<endl;
#endif
  bigint z=Ifloor(sqrt(I2bigfloat(T[i])/I2bigfloat(Q(i,i)))+0.1);
  bigint xmax =  z-U[i], xmin = -z-U[i], xi;
  int j, ok=1;
#ifdef DEBUG_LIST_SHORT_VECS
  cout<<"range for x["<<i<<"]: from "<<xmin<<" to "<<xmax<<endl;
#endif
  for(xi=xmin; ok&&(xi<=xmax); xi+=1)
    {
      x[i]=xi;
#ifdef DEBUG_LIST_SHORT_VECS
      cout<<"Setting x["<<i<<"] to "<<xi<<endl;
#endif
      if(i==1) 
	{
	  if(trivial(x)) return 0;
	  else 
	    {
	      vec_m v = comb(n,b,x);
#ifdef DEBUG_LIST_SHORT_VECS
	      cout<<"x="<<x<<"\t"; 
#endif
	      cout<<"v="<<v<<endl; 
	    }
	}
      else
	{
	  int i_minus_1=i-1;
	  U[i_minus_1]=0;
	  for(j=i; j<=n; j++) U[i_minus_1]+=Q(i_minus_1,j)*x[j];
	  T[i_minus_1]=T[i]-Q(i,i)*sqr(xi+U[i]);
	  ok=dig(n,b, i_minus_1,x,Q,T,U);
	}
    }
  return 1;
}

vec_m comb(const int n, const vec_m* b, const vec_m& x)
{
  vec_m v(n);
  int i;
  for(i=1; i<=n; i++) v+=x[i]*b[i];
  return v;
}
#endif

