// subspace.cc: implementations of subspace class
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
 
#include <eclib/subspace.h>

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i

#include "sub.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l

#include "sub.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar bigint
#define vec vec_m
#define mat mat_m
#define subspace subspace_m

#include "sub.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

///////////////////////////////////////////////////////////////////////////

// definitions of member operators and functions:

// assignment
template<class T>
void subspaceT<T>::operator=(const subspaceT<T>& s)
{
  pivots=s.pivots;
  basis=s.basis;
  denom=s.denom;
}

// Definitions of nonmember, nonfriend operators and functions:

template<class T>
subspaceT<T> combine(const subspaceT<T>& s1, const subspaceT<T>& s2)
{
  T d = s1.denom * s2.denom;
  const matT<T>& b1=s1.basis;
  const matT<T>& b2=s2.basis;
  matT<T> b = b1*b2;
  T g = b.content();
  if(g>1)
    {
      d/=g; b/=g;
    }
  vecT<int> p = s1.pivots[s2.pivots];
  return subspaceT<T>(b,p,d);
}

template<class T>
matT<T> expressvectors(const matT<T>& m, const subspaceT<T>& s)
{ vecT<int> p = pivots(s);
  long   n = dim(s);
  matT<T> ans(n,m.ncols());
  for (int i=1; i<=n; i++) ans.setrow(i, m.row(p[i]));
  return ans;
}

//This one is used a lot:
// M is nxn;
// S is a subspace of dim d<=n with nxd basis B, whose pivotal rows are a multiple den of the dxd identity
// return A such that den*M*B = B*A
//
// Algorithm: restricting to pivotal rows on both sides: B restricts to den*I and M to M1 (say), so
// den*M1*B=den*A, so A=M1*B

template<class T>
matT<T> restrict_mat(const matT<T>& M, const subspaceT<T>& S, int cr)
{
  if(dim(S)==M.nro) return M; // trivial special case, s is whole space
  const matT<T>& B = S.basis;
  matT<T> A = rowsubmat(M, S.pivots) * B;

  if(cr) // optional check that S is invariant under M
    {
      T m(DEFAULT_MODULUS);
      int check = (S.denom*matmulmodp(M,B,m) == matmulmodp(B,A,m));
      if (!check)
        cerr<<"Error in restrict_mat: subspace not invariant!"<<endl;
    }
  return A;
}

template<class T>
subspaceT<T> kernel(const matT<T>& m1, int method)
{
   long rank, nullity;
   T d;
   vecT<int> pcols,npcols;
   matT<T> m = echelon(m1,pcols,npcols, rank, nullity, d, method);
   matT<T> basis(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     basis.set(npcols[n],n,d);
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       basis.set(i,j, -m(r,npcols[j]));
   }
   return subspaceT<T>(basis, npcols, d);
}

template<class T>
subspaceT<T> image(const matT<T>& m, int method)
{
  vecT<int> p,np;
  long rank, nullity;
  T d;
  matT<T> b = transpose(echelon(transpose(m),p,np,rank,nullity,d,method));
  return subspaceT<T>(b,p,d);
}

template<class T>
subspaceT<T> eigenspace(const matT<T>& m1, const T& lambda, int method)
{
  matT<T> m = addscalar(m1,-lambda);
  return kernel(m,method);
}

template<class T>
subspaceT<T> subeigenspace(const matT<T>& m1, const T& l, const subspaceT<T>& s, int method)
{
  matT<T> m = restrict_mat(m1,s);
  subspaceT<T> ss = eigenspace(m, l*(denom(s)),method);
  return combine(s,ss );
}

template<class T>
subspaceT<T> pcombine(const subspaceT<T>& s1, const subspaceT<T>& s2, const T& pr)
{
  T   d = s1.denom * s2.denom;  // redundant since both should be 1
  const matT<T>& b1=s1.basis,  b2=s2.basis;
  const matT<T>& b = matmulmodp(b1,b2,pr);
  const vecT<int>& p = s1.pivots[s2.pivots];
  return subspaceT<T>(b,p,d);
}

// Same as restrict_mat, but modulo pr
template<class T>
matT<T> prestrict(const matT<T>& M, const subspaceT<T>& S, const T& pr, int cr)
{
  if(dim(S)==M.nro) return M; // trivial special case, s is whole space
  const matT<T>& B = S.basis;
  matT<T> A = matmulmodp(rowsubmat(M, S.pivots), B, pr);

  if(cr) // optional check that S is invariant under M
    {
      int check = (S.denom*matmulmodp(M,B,pr) == matmulmodp(B,A,pr));
      if (!check)
        cerr<<"Error in prestrict: subspace not invariant!"<<endl;
    }
  return A;
}

template<class T>
subspaceT<T> oldpkernel(const matT<T>& m1, const T& pr)   // using full echmodp
{
   long rank, nullity;
   vecT<int> pcols,npcols;
   matT<T> m = echmodp(m1,pcols,npcols, rank, nullity, pr);
   matT<T> basis(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     basis.set(npcols[n],n,T(1));
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       basis.set(i,j, mod(-m(r,npcols[j]),pr));
   }
   return subspaceT<T>(basis, npcols, T(1));
}

// using echmodp_uptri, with no back-substitution
template<class T>
subspaceT<T> pkernel(const matT<T>& m1, const T& pr)
{
  long rank, nullity;
  vecT<int> pcols,npcols;
  matT<T> m = echmodp_uptri(m1,pcols,npcols, rank, nullity, pr);
  matT<T> basis(m.ncols(),nullity);
  for(int j=nullity; j>0; j--)
    {
      int jj = npcols[j];
      basis(jj,j) = 1;
      for(int i=rank; i>0; i--)
        {
          T temp = -m(i,jj);
          for(int t=rank; t>i; t--)
            {
              int tt=pcols[t];
              temp -= xmodmul(m(i,tt),basis(tt,j),pr);
              temp = xmod(temp,pr);
            }
          basis(pcols[i],j) = mod(temp,pr);
        }
    }
  return subspaceT<T>(basis, npcols, T(1));
}

template<class T>
subspaceT<T> pimage(const matT<T>& m, const T& pr)
{
  vecT<int> p,np;
  long rank, nullity;
  const matT<T>& b = transpose(echmodp(transpose(m),p,np,rank,nullity,pr));
  return subspaceT<T>(b,p,T(1));
}

template<class T>
subspaceT<T> peigenspace(const matT<T>& m1, const T& lambda, const T& pr)
{
  const matT<T>& m = addscalar(m1,-lambda);
  return pkernel(m,pr);
}

template<class T>
subspaceT<T> psubeigenspace(const matT<T>& m1, const T& l, const subspaceT<T>& s, const T& pr)
{
  const matT<T>& m = prestrict(m1,s,pr);
  const subspaceT<T>& ss = peigenspace(m, l*(denom(s)),pr);
  return pcombine(s,ss,pr);
}

//Attempts to lift from a mod-p subspace to a normal Q-subspace by expressing
//basis as rational using modrat and clearing denominators

template<class T>
int lift(const subspaceT<T>& s, const T& pr, subspaceT<T>& ans)
{
  T dd;
  matT<T> m;
  int ok = liftmat(s.basis,pr,m,dd);
  if (!ok)
    cerr << "Failed to lift subspace from mod "<<pr<<endl;
  ans = subspaceT<T>(m, pivots(s), dd);
  return ok;
}
