// subspace.cc: implementations of subspace class
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

// definitions of member operators and functions:

// assignment
template<class T>
void subZspace<T>::operator=(const subZspace<T>& s)
{
  pivots=s.pivots;
  basis=s.basis;
  denom=s.denom;
}

template<class T>
int subZspace<T>::contains(const Zvec<T>& v) const           // does this contain v?
{
  if (v.dim() != basis.nrows()) return 0;
  return basis*v[pivots]==denom*v;
}

template<class T>
int subZspace<T>::contains(const subZspace<T>& s) const      // does this contain s?
{
  // check ambient spaces are the same:
  if (basis.nrows()!=s.basis.nrows())
    return 0;
  // check dim(s) <= dim:
  if (basis.ncols()<s.basis.ncols())
    return 0;
  for(int j=1; j<=s.basis.ncols(); j++)
    if (!contains(s.basis.col(j)))
      return 0;
  return 1;
}

// Definitions of nonmember, nonfriend operators and functions:

template<class T>
subZspace<T> combine(const subZspace<T>& s1, const subZspace<T>& s2)
{
  T d = s1.den() * s2.den();
  const Zmat<T>& b1=s1.bas();
  const Zmat<T>& b2=s2.bas();
  Zmat<T> b = b1*b2;
  T g = b.content();
  if(g>1)
    {
      d/=g; b/=g;
    }
  Zvec<int> p = s1.pivs()[s2.pivs()];
  return subZspace<T>(b,p,d);
}
template subZspace<int> combine<int>(const subZspace<int>&, const subZspace<int>&);
template subZspace<long> combine<long>(const subZspace<long>&, const subZspace<long>&);
template subZspace<ZZ> combine<ZZ>(const subZspace<ZZ>&, const subZspace<ZZ>&);
template subZspace<INT> combine<INT>(const subZspace<INT>& , const subZspace<INT>&);

template<class T>
Zmat<T> expressvectors(const Zmat<T>& m, const subZspace<T>& s)
{ Zvec<int> p = s.pivs();
  long   n = s.dim();
  Zmat<T> ans(n,m.ncols());
  for (int i=1; i<=n; i++) ans.setrow(i, m.row(p[i]));
  return ans;
}
template Zmat<int> expressvectors<int>(const Zmat<int>&, const subZspace<int>&);
template Zmat<long> expressvectors<long>(const Zmat<long>&, const subZspace<long>&);
template Zmat<ZZ> expressvectors<ZZ>(const Zmat<ZZ>&, const subZspace<ZZ>&);
template Zmat<INT> expressvectors<INT>(const Zmat<INT>&, const subZspace<INT>&);

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
  if(S.dim()==M.nro) return M; // trivial special case, s is whole space
  const Zmat<T>& B = S.bas();
  // cerr<<"In restrict_mat() with M =\n" << M <<"\n and basis(S) = \n"<<B<<endl;
  // cerr<<"rowsubmat(M, S.pivs()) =\n"<<rowsubmat(M, S.pivs())<<endl;
  Zmat<T> A = rowsubmat(M, S.pivs()) * B;
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
template Zmat<int> restrict_mat<int>(const Zmat<int>&, const subZspace<int>&, int);
template Zmat<long> restrict_mat<long>(const Zmat<long>&, const subZspace<long>&, int);
template Zmat<ZZ> restrict_mat<ZZ>(const Zmat<ZZ>&, const subZspace<ZZ>&, int);
template Zmat<INT> restrict_mat<INT>(const Zmat<INT>&, const subZspace<INT>&, int);

template<class T>
subZspace<T> kernel(const Zmat<T>& m1)
{
   long rk, ny;
   T d;
   Zvec<int> pcols,npcols;
   Zmat<T> m = ref(m1, d, pcols, npcols, rk, ny);
   Zmat<T> bas(m.ncols(),ny);
   for (int n=1; n<=ny; n++)
     bas.set(npcols[n],n,d);
   for (int r=1; r<=rk; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=ny; j++)
       bas.set(i,j, -m(r,npcols[j]));
   }
   return subZspace<T>(bas, npcols, d);
}
template subZspace<int> kernel<int>(const Zmat<int>&);
template subZspace<long> kernel<long>(const Zmat<long>&);
template subZspace<ZZ> kernel<ZZ>(const Zmat<ZZ>&);
template subZspace<INT> kernel<INT>(const Zmat<INT>&);

template<class T>
subZspace<T> image(const Zmat<T>& m)
{
  Zvec<int> p,np;
  long rk, ny;
  T d;
  Zmat<T> b = transpose(ref(transpose(m),d,p,np,rk,ny));
  return subZspace<T>(b,p,d);
}
template subZspace<int> image<int>(const Zmat<int>& m);
template subZspace<long> image<long>(const Zmat<long>& m);
template subZspace<ZZ> image<ZZ>(const Zmat<ZZ>& m);
template subZspace<INT> image<INT>(const Zmat<INT>& m);

template<class T>
subZspace<T> eigenspace(const Zmat<T>& m1, const T& lambda)
{
  Zmat<T> m = addscalar(m1,-lambda);
  return kernel(m);
}
template subZspace<int> eigenspace<int>(const Zmat<int>& m, const int& lambda);
template subZspace<long> eigenspace<long>(const Zmat<long>& m, const long& lambda);
template subZspace<ZZ> eigenspace<ZZ>(const Zmat<ZZ>& m, const ZZ& lambda);
template subZspace<INT> eigenspace<INT>(const Zmat<INT>& m, const INT& lambda);

template<class T>
subZspace<T> subeigenspace(const Zmat<T>& m1, const T& l, const subZspace<T>& s)
{
  Zmat<T> m = restrict_mat(m1,s);
  subZspace<T> ss = eigenspace(m, l*(s.den()));
  return combine(s,ss);
}
template subZspace<int> subeigenspace<int>(const Zmat<int>&, const int&, const subZspace<int>&);
template subZspace<long> subeigenspace<long>(const Zmat<long>&, const long&, const subZspace<long>&);
template subZspace<ZZ> subeigenspace<ZZ>(const Zmat<ZZ>&, const ZZ&, const subZspace<ZZ>&);
template subZspace<INT> subeigenspace<INT>(const Zmat<INT>&, const INT&, const subZspace<INT>&);

template<class T>
subZspace<T> pcombine(const subZspace<T>& s1, const subZspace<T>& s2, const T& pr)
{
  T   d = s1.denom * s2.denom;  // redundant since both should be 1
  const Zmat<T>& b1=s1.basis,  b2=s2.basis;
  const Zmat<T>& b = matmulmodp(b1,b2,pr);
  const Zvec<int>& p = s1.pivots[s2.pivots];
  return subZspace<T>(b,p,d);
}
template subZspace<int> pcombine<int>(const subZspace<int>&, const subZspace<int>&, const int&);
template subZspace<long> pcombine<long>(const subZspace<long>&, const subZspace<long>&, const long&);
template subZspace<ZZ> pcombine<ZZ>(const subZspace<ZZ>&, const subZspace<ZZ>&, const ZZ&);
template subZspace<INT> pcombine<INT>(const subZspace<INT>&, const subZspace<INT>&, const INT&);

// Same as restrict_mat, but modulo pr
template<class T>
Zmat<T> prestrict(const Zmat<T>& M, const subZspace<T>& S, const T& pr, int cr)
{
  if(S.dim()==M.nrows()) return M; // trivial special case, s is whole space
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
template Zmat<int> prestrict<int>(const Zmat<int>&, const subZspace<int>&, const int&, int);
template Zmat<long> prestrict<long>(const Zmat<long>&, const subZspace<long>&, const long&, int);
template Zmat<ZZ> prestrict<ZZ>(const Zmat<ZZ>&, const subZspace<ZZ>&, const ZZ&, int);
template Zmat<INT> prestrict<INT>(const Zmat<INT>&, const subZspace<INT>&, const INT&, int);

template<class T>
subZspace<T> oldpkernel(const Zmat<T>& m1, const T& pr)
{
   long rk, ny;
   Zvec<int> pcols,npcols;
   Zmat<T> m = ref_mod_p(m1, pr, pcols, npcols, rk, ny);
   Zmat<T> bas(m.ncols(),ny);
   for (int n=1; n<=ny; n++)
     bas.set(npcols[n],n,T(1));
   for (int r=1; r<=rk; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=ny; j++)
       bas.set(i,j, mod(-m(r,npcols[j]),pr));
   }
   return subZspace<T>(bas, npcols, T(1));
}
template subZspace<int> oldpkernel<int>(const Zmat<int>&, const int&);
template subZspace<long> oldpkernel<long>(const Zmat<long>&, const long&);
template subZspace<ZZ> oldpkernel<ZZ>(const Zmat<ZZ>&, const ZZ&);
template subZspace<INT> oldpkernel<INT>(const Zmat<INT>&, const INT&);

template<class T>
subZspace<T> pkernel(const Zmat<T>& m1, const T& pr)
{
  long rk, ny;
  Zvec<int> pcols,npcols;
  Zmat<T> m = ref_mod_p(m1, pr, pcols, npcols, rk, ny);
  Zmat<T> bas(m.ncols(),ny);
  for(int j=ny; j>0; j--)
    {
      int jj = npcols[j];
      bas(jj,j) = 1;
      for(int i=rk; i>0; i--)
        {
          T temp = -m(i,jj);
          for(int t = rk; t>i; t--)
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
template subZspace<int> pkernel<int>(const Zmat<int>&, const int&);
template subZspace<long> pkernel<long>(const Zmat<long>&, const long&);
template subZspace<ZZ> pkernel<ZZ>(const Zmat<ZZ>&, const ZZ&);
template subZspace<INT> pkernel<INT>(const Zmat<INT>&, const INT&);

template<class T>
subZspace<T> pimage(const Zmat<T>& m, const T& pr)
{
  Zvec<int> p,np;
  long rk, ny;
  const Zmat<T>& b = transpose(ref_mod_p(transpose(m),pr, p,np, rk,ny));
  return subZspace<T>(b,p,T(1));
}
template subZspace<int> pimage<int>(const Zmat<int>&, const int&);
template subZspace<long> pimage<long>(const Zmat<long>&, const long&);
template subZspace<ZZ> pimage<ZZ>(const Zmat<ZZ>&, const ZZ&);
template subZspace<INT> pimage<INT>(const Zmat<INT>&, const INT&);

template<class T>
subZspace<T> peigenspace(const Zmat<T>& m1, const T& lambda, const T& pr)
{
  const Zmat<T>& m = addscalar(m1,-lambda);
  return pkernel(m,pr);
}
template subZspace<int> peigenspace<int>(const Zmat<int>&, const int&, const int&);
template subZspace<long> peigenspace<long>(const Zmat<long>&, const long&, const long&);
template subZspace<ZZ> peigenspace<ZZ>(const Zmat<ZZ>&, const ZZ&, const ZZ&);
template subZspace<INT> peigenspace<INT>(const Zmat<INT>&, const INT&, const INT&);

template<class T>
subZspace<T> psubeigenspace(const Zmat<T>& m1, const T& l, const subZspace<T>& s, const T& pr)
{
  const Zmat<T>& m = prestrict(m1,s,pr);
  const subZspace<T>& ss = peigenspace(m, l*(s.den()),pr);
  return pcombine(s,ss,pr);
}
template subZspace<int> psubeigenspace<int>(const Zmat<int>&, const int&, const subZspace<int>&, const int&);
template subZspace<long> psubeigenspace<long>(const Zmat<long>&, const long&, const subZspace<long>&, const long&);
template subZspace<ZZ> psubeigenspace<ZZ>(const Zmat<ZZ>&, const ZZ&, const subZspace<ZZ>&, const ZZ&);
template subZspace<INT> psubeigenspace<INT>(const Zmat<INT>&, const INT&, const subZspace<INT>&, const INT&);

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
  ans = subZspace<T>(m, s.pivs(), dd);
  return ok;
}
template int lift<int>(const subZspace<int>&, const int&, subZspace<int>&);
template int lift<long>(const subZspace<long>&, const long&, subZspace<long>&);
template int lift<ZZ>(const subZspace<ZZ>&, const ZZ&, subZspace<ZZ>&);
template int lift<INT>(const subZspace<INT>&, const INT&, subZspace<INT>&);

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
        m.set(1,j, NTL::bit(a,j-1));
      subspace_i ker = pkernel(m,2); // right kernel mod 2
      // assert (dim(ker)==r-1);
      mat_i bas = ker.bas();
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
      m.set(i,j, NTL::bit(alist[i-1],j-1));
  subspace_i ker = pkernel(m,2); // right kernel mod 2
  // assert (ker.dim()==r-s);
  mat_i bas = ker.bas();
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

// Instantiate subZspace template classes for T=int, long, ZZ, INT

template class subZspace<int>;
template class subZspace<long>;
template class subZspace<ZZ>;
template class subZspace<INT>;
