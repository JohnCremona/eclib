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

// Instantiate subspaceT template classes for T=int, long, bigint

template class subspaceT<int>;
template class subspaceT<long>;
template class subspaceT<bigint>;

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
   matT<T> bas(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     bas.set(npcols[n],n,d);
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       bas.set(i,j, -m(r,npcols[j]));
   }
   return subspaceT<T>(bas, npcols, d);
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
   matT<T> bas(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     bas.set(npcols[n],n,T(1));
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       bas.set(i,j, mod(-m(r,npcols[j]),pr));
   }
   return subspaceT<T>(bas, npcols, T(1));
}

// using echmodp_uptri, with no back-substitution
template<class T>
subspaceT<T> pkernel(const matT<T>& m1, const T& pr)
{
  long rank, nullity;
  vecT<int> pcols,npcols;
  matT<T> m = echmodp_uptri(m1,pcols,npcols, rank, nullity, pr);
  matT<T> bas(m.ncols(),nullity);
  for(int j=nullity; j>0; j--)
    {
      int jj = npcols[j];
      bas(jj,j) = 1;
      for(int i=rank; i>0; i--)
        {
          T temp = -m(i,jj);
          for(int t=rank; t>i; t--)
            {
              int tt=pcols[t];
              temp -= xmodmul(m(i,tt),bas(tt,j),pr);
              temp = xmod(temp,pr);
            }
          bas(pcols[i],j) = mod(temp,pr);
        }
    }
  return subspaceT<T>(bas, npcols, T(1));
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

// Instantiate template functions for T=int

template int dim<int>(const subspaceT<int>& s);
template int denom<int>(const subspaceT<int>& s);
template vecT<int> pivots<int>(const subspaceT<int>& s);
template matT<int> basis<int>(const subspaceT<int>& s);
template subspaceT<int> combine<int>(const subspaceT<int>& s1, const subspaceT<int>& s2);
template matT<int> restrict_mat<int>(const matT<int>& m, const subspaceT<int>& s, int cr);
template subspaceT<int> pcombine<int>(const subspaceT<int>& s1, const subspaceT<int>& s2, const int& pr);
template matT<int> prestrict<int>(const matT<int>& m, const subspaceT<int>& s, const int& pr, int cr);
template int lift<int>(const subspaceT<int>& s, const int& pr, subspaceT<int>& ans);
template matT<int> expressvectors<int>(const matT<int>& m, const subspaceT<int>& s);
template subspaceT<int> kernel<int>(const matT<int>& m, int method=0);
template subspaceT<int> image<int>(const matT<int>& m, int method=0);
template subspaceT<int> eigenspace<int>(const matT<int>& m, const int& lambda, int method=0);
template subspaceT<int> subeigenspace<int>(const matT<int>& m, const int& l, const subspaceT<int>& s, int method=0);
template subspaceT<int> oldpkernel<int>(const matT<int>& m, const int& pr);
template subspaceT<int> pkernel<int>(const matT<int>& m, const int& pr);
template subspaceT<int> pimage<int>(const matT<int>& m, const int& pr);
template subspaceT<int> peigenspace<int>(const matT<int>& m, const int& lambda, const int& pr);
template subspaceT<int> psubeigenspace<int>(const matT<int>& m, const int& l, const subspaceT<int>& s, const int& pr);

// Instantiate template functions for T=long

template int dim<long>(const subspaceT<long>& s);
template long denom<long>(const subspaceT<long>& s);
template vecT<int> pivots<long>(const subspaceT<long>& s);
template matT<long> basis<long>(const subspaceT<long>& s);
template subspaceT<long> combine<long>(const subspaceT<long>& s1, const subspaceT<long>& s2);
template matT<long> restrict_mat<long>(const matT<long>& m, const subspaceT<long>& s, int cr);
template subspaceT<long> pcombine<long>(const subspaceT<long>& s1, const subspaceT<long>& s2, const long& pr);
template matT<long> prestrict<long>(const matT<long>& m, const subspaceT<long>& s, const long& pr, int cr);
template int lift<long>(const subspaceT<long>& s, const long& pr, subspaceT<long>& ans);
template matT<long> expressvectors<long>(const matT<long>& m, const subspaceT<long>& s);
template subspaceT<long> kernel<long>(const matT<long>& m, int method=0);
template subspaceT<long> image<long>(const matT<long>& m, int method=0);
template subspaceT<long> eigenspace<long>(const matT<long>& m, const long& lambda, int method=0);
template subspaceT<long> subeigenspace<long>(const matT<long>& m, const long& l, const subspaceT<long>& s, int method=0);
template subspaceT<long> oldpkernel<long>(const matT<long>& m, const long& pr);
template subspaceT<long> pkernel<long>(const matT<long>& m, const long& pr);
template subspaceT<long> pimage<long>(const matT<long>& m, const long& pr);
template subspaceT<long> peigenspace<long>(const matT<long>& m, const long& lambda, const long& pr);
template subspaceT<long> psubeigenspace<long>(const matT<long>& m, const long& l, const subspaceT<long>& s, const long& pr);

// Instantiate template functions for T=bigint

template int dim<bigint>(const subspaceT<bigint>& s);
template bigint denom<bigint>(const subspaceT<bigint>& s);
template vecT<int> pivots<bigint>(const subspaceT<bigint>& s);
template matT<bigint> basis<bigint>(const subspaceT<bigint>& s);
template subspaceT<bigint> combine<bigint>(const subspaceT<bigint>& s1, const subspaceT<bigint>& s2);
template matT<bigint> restrict_mat<bigint>(const matT<bigint>& m, const subspaceT<bigint>& s, int cr);
template subspaceT<bigint> pcombine<bigint>(const subspaceT<bigint>& s1, const subspaceT<bigint>& s2, const bigint& pr);
template matT<bigint> prestrict<bigint>(const matT<bigint>& m, const subspaceT<bigint>& s, const bigint& pr, int cr);
template int lift<bigint>(const subspaceT<bigint>& s, const bigint& pr, subspaceT<bigint>& ans);
template matT<bigint> expressvectors<bigint>(const matT<bigint>& m, const subspaceT<bigint>& s);
template subspaceT<bigint> kernel<bigint>(const matT<bigint>& m, int method=0);
template subspaceT<bigint> image<bigint>(const matT<bigint>& m, int method=0);
template subspaceT<bigint> eigenspace<bigint>(const matT<bigint>& m, const bigint& lambda, int method=0);
template subspaceT<bigint> subeigenspace<bigint>(const matT<bigint>& m, const bigint& l, const subspaceT<bigint>& s, int method=0);
template subspaceT<bigint> oldpkernel<bigint>(const matT<bigint>& m, const bigint& pr);
template subspaceT<bigint> pkernel<bigint>(const matT<bigint>& m, const bigint& pr);
template subspaceT<bigint> pimage<bigint>(const matT<bigint>& m, const bigint& pr);
template subspaceT<bigint> peigenspace<bigint>(const matT<bigint>& m, const bigint& lambda, const bigint& pr);
template subspaceT<bigint> psubeigenspace<bigint>(const matT<bigint>& m, const bigint& l, const subspaceT<bigint>& s, const bigint& pr);
