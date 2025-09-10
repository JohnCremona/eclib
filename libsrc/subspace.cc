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
 
#include "eclib/linalg.h"

// Instantiate subZspace template classes for T=int, long, bigint

template class subZspace<int>;
template class subZspace<long>;
template class subZspace<bigint>;

// definitions of member operators and functions:

// assignment
template<class T>
void subZspace<T>::operator=(const subZspace<T>& s)
{
  pivots=s.pivots;
  basis=s.basis;
  denom=s.denom;
}

// Definitions of nonmember, nonfriend operators and functions:

template<class T>
subZspace<T> combine(const subZspace<T>& s1, const subZspace<T>& s2)
{
  T d = s1.denom * s2.denom;
  const Zmat<T>& b1=s1.basis;
  const Zmat<T>& b2=s2.basis;
  Zmat<T> b = b1*b2;
  T g = b.content();
  if(g>1)
    {
      d/=g; b/=g;
    }
  Zvec<int> p = s1.pivots[s2.pivots];
  return subZspace<T>(b,p,d);
}

template<class T>
Zmat<T> expressvectors(const Zmat<T>& m, const subZspace<T>& s)
{ Zvec<int> p = pivots(s);
  long   n = dim(s);
  Zmat<T> ans(n,m.ncols());
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
Zmat<T> restrict_mat(const Zmat<T>& M, const subZspace<T>& S, int cr)
{
  if(dim(S)==M.nro) return M; // trivial special case, s is whole space
  const Zmat<T>& B = S.basis;
  // cerr<<"In restrict_mat() with M =\n" << M <<"\n and basis(S) = \n"<<B<<endl;
  // cerr<<"rowsubmat(M, S.pivots) =\n"<<rowsubmat(M, S.pivots)<<endl;
  Zmat<T> A = rowsubmat(M, S.pivots) * B;
  // cerr<<"A            =\n"<< A << endl;

  if(cr) // optional check that S is invariant under M
    {
      T modulus = default_modulus<T>();
      Zmat<T>
        RHS = matmulmodp(B,A, modulus),
        LHS = matmulmodp(M,B, modulus);
      LHS *= S.denom;
      LHS.reduce_mod_p(modulus);
      if (LHS!=RHS)
        {
          cerr<<"Error in restrict_mat: subspace not invariant!"<<endl;
          cerr<<"denom(S)*M*B =\n"<< LHS << endl;
          cerr<<"B*A          =\n"<< RHS << endl;
          cerr<<"(using modulus "<<modulus<<")"<<endl;
        }
    }
  return A;
}

template<class T>
subZspace<T> kernel(const Zmat<T>& m1, int method)
{
   long rank, nullity;
   T d;
   Zvec<int> pcols,npcols;
   Zmat<T> m = echelon(m1,pcols,npcols, rank, nullity, d, method);
   Zmat<T> bas(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     bas.set(npcols[n],n,d);
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       bas.set(i,j, -m(r,npcols[j]));
   }
   return subZspace<T>(bas, npcols, d);
}

template<class T>
subZspace<T> image(const Zmat<T>& m, int method)
{
  Zvec<int> p,np;
  long rank, nullity;
  T d;
  Zmat<T> b = transpose(echelon(transpose(m),p,np,rank,nullity,d,method));
  return subZspace<T>(b,p,d);
}

template<class T>
subZspace<T> eigenspace(const Zmat<T>& m1, const T& lambda, int method)
{
  Zmat<T> m = addscalar(m1,-lambda);
  return kernel(m,method);
}

template<class T>
subZspace<T> subeigenspace(const Zmat<T>& m1, const T& l, const subZspace<T>& s, int method)
{
  Zmat<T> m = restrict_mat(m1,s);
  subZspace<T> ss = eigenspace(m, l*(denom(s)),method);
  return combine(s,ss );
}

template<class T>
subZspace<T> pcombine(const subZspace<T>& s1, const subZspace<T>& s2, const T& pr)
{
  T   d = s1.denom * s2.denom;  // redundant since both should be 1
  const Zmat<T>& b1=s1.basis,  b2=s2.basis;
  const Zmat<T>& b = matmulmodp(b1,b2,pr);
  const Zvec<int>& p = s1.pivots[s2.pivots];
  return subZspace<T>(b,p,d);
}

// Same as restrict_mat, but modulo pr
template<class T>
Zmat<T> prestrict(const Zmat<T>& M, const subZspace<T>& S, const T& pr, int cr)
{
  if(dim(S)==M.nro) return M; // trivial special case, s is whole space
  const Zmat<T>& B = S.basis;
  Zmat<T> A = matmulmodp(rowsubmat(M, S.pivots), B, pr);

  if(cr) // optional check that S is invariant under M
    {
      int check = (S.denom*matmulmodp(M,B,pr) == matmulmodp(B,A,pr));
      if (!check)
        cerr<<"Error in prestrict: subspace not invariant!"<<endl;
    }
  return A;
}

template<class T>
subZspace<T> oldpkernel(const Zmat<T>& m1, const T& pr)   // using full echmodp
{
   long rank, nullity;
   Zvec<int> pcols,npcols;
   Zmat<T> m = echmodp(m1,pcols,npcols, rank, nullity, pr);
   Zmat<T> bas(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     bas.set(npcols[n],n,T(1));
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       bas.set(i,j, mod(-m(r,npcols[j]),pr));
   }
   return subZspace<T>(bas, npcols, T(1));
}

// using echmodp_uptri, with no back-substitution
template<class T>
subZspace<T> pkernel(const Zmat<T>& m1, const T& pr)
{
  long rank, nullity;
  Zvec<int> pcols,npcols;
  Zmat<T> m = echmodp_uptri(m1,pcols,npcols, rank, nullity, pr);
  Zmat<T> bas(m.ncols(),nullity);
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
  return subZspace<T>(bas, npcols, T(1));
}

template<class T>
subZspace<T> pimage(const Zmat<T>& m, const T& pr)
{
  Zvec<int> p,np;
  long rank, nullity;
  const Zmat<T>& b = transpose(echmodp(transpose(m),p,np,rank,nullity,pr));
  return subZspace<T>(b,p,T(1));
}

template<class T>
subZspace<T> peigenspace(const Zmat<T>& m1, const T& lambda, const T& pr)
{
  const Zmat<T>& m = addscalar(m1,-lambda);
  return pkernel(m,pr);
}

template<class T>
subZspace<T> psubeigenspace(const Zmat<T>& m1, const T& l, const subZspace<T>& s, const T& pr)
{
  const Zmat<T>& m = prestrict(m1,s,pr);
  const subZspace<T>& ss = peigenspace(m, l*(denom(s)),pr);
  return pcombine(s,ss,pr);
}

//Attempts to lift from a mod-p subspace to a normal Q-subspace by expressing
//basis as rational using modrat and clearing denominators

template<class T>
int lift(const subZspace<T>& s, const T& pr, subZspace<T>& ans)
{
  T dd;
  Zmat<T> m;
  int ok = liftmat(s.basis,pr,m,dd);
  if (!ok)
    cerr << "Failed to lift subspace from mod "<<pr<<endl;
  ans = subZspace<T>(m, pivots(s), dd);
  return ok;
}

// return a basis for the orthogonal complement of a<2^r (viewed as a bit vector of length r)
vector<long> dotperp(long a, int r)
{
  if (a==0) // trivial special case
    {
      vector<long> ans(r);
      for (int i=0; i<r; i++)
        ans[i] = 1<<i;
      return ans;
    }
  else
    {
      mat_i m(1,r);
      for (int j=1; j<=r; j++)
        m.set(1,j,bit(a,j-1));
      subspace_i ker = pkernel(m,2); // right kernel mod 2
      // assert (dim(ker)==r-1);
      mat_i bas = basis(ker);
      vector<long> ans(r-1, 0);
      for (int i=0; i<r-1; i++)
        {
          vec_i coli = bas.col(i+1);
          for (int j=0; j<r; j++)
            if (coli[j+1])
              ans[i] |= 1<<j;
        }
      return ans;
    }
}

// return a basis for the orthogonal complement of the span of a in alist (viewed as bit vectors of length r)
vector<long> dotperp(const vector<long>& alist, int r)
{
  int s = alist.size();
  mat_i m(s,r);
  for (int i=1; i<=s; i++)
    for (int j=1; j<=r; j++)
      m.set(i,j,bit(alist[i-1],j-1));
  subspace_i ker = pkernel(m,2); // right kernel mod 2
  // assert (dim(ker)==r-s);
  mat_i bas = basis(ker);
  vector<long> ans(r-s, 0);
  for (int i=0; i<r-s; i++)
    {
      vec_i coli = bas.col(i+1);
      for (int j=0; j<r; j++)
        if (coli[j+1])
          ans[i] |= 1<<j;
    }
  return ans;
}


// Instantiate template functions for T=int

template int dim<int>(const subZspace<int>& s);
template int denom<int>(const subZspace<int>& s);
template Zvec<int> pivots<int>(const subZspace<int>& s);
template Zmat<int> basis<int>(const subZspace<int>& s);
template subZspace<int> combine<int>(const subZspace<int>& s1, const subZspace<int>& s2);
template Zmat<int> restrict_mat<int>(const Zmat<int>& m, const subZspace<int>& s, int cr);
template subZspace<int> pcombine<int>(const subZspace<int>& s1, const subZspace<int>& s2, const int& pr);
template Zmat<int> prestrict<int>(const Zmat<int>& m, const subZspace<int>& s, const int& pr, int cr);
template int lift<int>(const subZspace<int>& s, const int& pr, subZspace<int>& ans);
template Zmat<int> expressvectors<int>(const Zmat<int>& m, const subZspace<int>& s);
template subZspace<int> kernel<int>(const Zmat<int>& m, int method=0);
template subZspace<int> image<int>(const Zmat<int>& m, int method=0);
template subZspace<int> eigenspace<int>(const Zmat<int>& m, const int& lambda, int method=0);
template subZspace<int> subeigenspace<int>(const Zmat<int>& m, const int& l, const subZspace<int>& s, int method=0);
template subZspace<int> oldpkernel<int>(const Zmat<int>& m, const int& pr);
template subZspace<int> pkernel<int>(const Zmat<int>& m, const int& pr);
template subZspace<int> pimage<int>(const Zmat<int>& m, const int& pr);
template subZspace<int> peigenspace<int>(const Zmat<int>& m, const int& lambda, const int& pr);
template subZspace<int> psubeigenspace<int>(const Zmat<int>& m, const int& l, const subZspace<int>& s, const int& pr);

// Instantiate template functions for T=long

template int dim<long>(const subZspace<long>& s);
template long denom<long>(const subZspace<long>& s);
template Zvec<int> pivots<long>(const subZspace<long>& s);
template Zmat<long> basis<long>(const subZspace<long>& s);
template subZspace<long> combine<long>(const subZspace<long>& s1, const subZspace<long>& s2);
template Zmat<long> restrict_mat<long>(const Zmat<long>& m, const subZspace<long>& s, int cr);
template subZspace<long> pcombine<long>(const subZspace<long>& s1, const subZspace<long>& s2, const long& pr);
template Zmat<long> prestrict<long>(const Zmat<long>& m, const subZspace<long>& s, const long& pr, int cr);
template int lift<long>(const subZspace<long>& s, const long& pr, subZspace<long>& ans);
template Zmat<long> expressvectors<long>(const Zmat<long>& m, const subZspace<long>& s);
template subZspace<long> kernel<long>(const Zmat<long>& m, int method=0);
template subZspace<long> image<long>(const Zmat<long>& m, int method=0);
template subZspace<long> eigenspace<long>(const Zmat<long>& m, const long& lambda, int method=0);
template subZspace<long> subeigenspace<long>(const Zmat<long>& m, const long& l, const subZspace<long>& s, int method=0);
template subZspace<long> oldpkernel<long>(const Zmat<long>& m, const long& pr);
template subZspace<long> pkernel<long>(const Zmat<long>& m, const long& pr);
template subZspace<long> pimage<long>(const Zmat<long>& m, const long& pr);
template subZspace<long> peigenspace<long>(const Zmat<long>& m, const long& lambda, const long& pr);
template subZspace<long> psubeigenspace<long>(const Zmat<long>& m, const long& l, const subZspace<long>& s, const long& pr);

// Instantiate template functions for T=bigint

template int dim<bigint>(const subZspace<bigint>& s);
template bigint denom<bigint>(const subZspace<bigint>& s);
template Zvec<int> pivots<bigint>(const subZspace<bigint>& s);
template Zmat<bigint> basis<bigint>(const subZspace<bigint>& s);
template subZspace<bigint> combine<bigint>(const subZspace<bigint>& s1, const subZspace<bigint>& s2);
template Zmat<bigint> restrict_mat<bigint>(const Zmat<bigint>& m, const subZspace<bigint>& s, int cr);
template subZspace<bigint> pcombine<bigint>(const subZspace<bigint>& s1, const subZspace<bigint>& s2, const bigint& pr);
template Zmat<bigint> prestrict<bigint>(const Zmat<bigint>& m, const subZspace<bigint>& s, const bigint& pr, int cr);
template int lift<bigint>(const subZspace<bigint>& s, const bigint& pr, subZspace<bigint>& ans);
template Zmat<bigint> expressvectors<bigint>(const Zmat<bigint>& m, const subZspace<bigint>& s);
template subZspace<bigint> kernel<bigint>(const Zmat<bigint>& m, int method=0);
template subZspace<bigint> image<bigint>(const Zmat<bigint>& m, int method=0);
template subZspace<bigint> eigenspace<bigint>(const Zmat<bigint>& m, const bigint& lambda, int method=0);
template subZspace<bigint> subeigenspace<bigint>(const Zmat<bigint>& m, const bigint& l, const subZspace<bigint>& s, int method=0);
template subZspace<bigint> oldpkernel<bigint>(const Zmat<bigint>& m, const bigint& pr);
template subZspace<bigint> pkernel<bigint>(const Zmat<bigint>& m, const bigint& pr);
template subZspace<bigint> pimage<bigint>(const Zmat<bigint>& m, const bigint& pr);
template subZspace<bigint> peigenspace<bigint>(const Zmat<bigint>& m, const bigint& lambda, const bigint& pr);
template subZspace<bigint> psubeigenspace<bigint>(const Zmat<bigint>& m, const bigint& l, const subZspace<bigint>& s, const bigint& pr);
