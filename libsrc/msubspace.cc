// msubspace.cc: implementations of multiprecision subspace class
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
 
//#define CHECK_RESTRICT   // define this to make restrict_mat and prestrict check
                           //  invariance of subspaces.
#include <eclib/marith.h>
#include <eclib/msubspace.h>

// Definitions of nonmember, nonfriend operators and functions:

msubspace combine(const msubspace& s1, const msubspace& s2)
{
  bigint d = s1.denom * s2.denom;
  const mat_m& b1=s1.basis;
  const mat_m& b2=s2.basis;
  mat_m b = b1*b2;
  bigint g = b.content();
  if(g>1)
    {
      d/=g; b/=g;
    }
  vec_i p = s1.pivots[s2.pivots];
  return msubspace(b,p,d);
}

mat_m restrict_mat(const mat_m& M, const msubspace& S)
{
  if(dim(S)==M.nro) return M; // trivial special case, s is whole space
  return rowsubmat(M, S.pivots) * S.basis;
}

msubspace kernel(const mat_m& mat, int method)
{
   long rank, nullity;
   bigint d;
   vec_i pcols,npcols;
   mat_m m = echelon(mat,pcols,npcols, rank, nullity, d, method);
   mat_m basis(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     basis.set(npcols[n],n,d);
   for (int r=1; r<=rank; r++)
     {
       int i = pcols[r];
       for (int j=1; j<=nullity; j++)
         basis.set(i,j, -m(r,npcols[j]));
     }
   msubspace ans(basis, npcols, d);
   return ans;
}

msubspace image(const mat_m& mat, int method)
{
  vec_i p,np;
  bigint d;
  long rank, nullity;
  mat_m b = transpose(echelon(transpose(mat),p,np,rank,nullity,d,method));
  msubspace ans(b,p,d);
  return ans;
}
 
msubspace eigenspace(const mat_m& mat, const bigint& lambda, int method)
{
  mat_m m = addscalar(mat,-lambda);
  msubspace ans = kernel(m,method);
  return ans;
}
 
msubspace subeigenspace(const mat_m& mat, const bigint& l, const msubspace& s, int method)
{
  mat_m m = restrict_mat(mat,s);
  msubspace ss = eigenspace(m, l*(denom(s)),method);
  msubspace ans = combine(s,ss );
  return ans;
}

msubspace pcombine(const msubspace& s1, const msubspace& s2, const bigint& pr)
{
  bigint   d = s1.denom * s2.denom;  // redundant since both should be 1
  mat_m b1=s1.basis,  b2=s2.basis;
  mat_m b = matmulmodp(b1,b2,pr);
  vec_i p = s1.pivots[s2.pivots];
  return msubspace(b,p,d);
}

mat_m prestrict(const mat_m& M, const msubspace& S, const bigint& pr)
{
  if(dim(S)==M.nro) return M; // trivial special case, s is whole space
  return matmulmodp(rowsubmat(M, S.pivots), S.basis, pr);
}

msubspace pkernel(const mat_m& mat, const bigint& pr)
{
   long rank, nullity;
   vec_i pcols,npcols;
   mat_m m = echmodp(mat,pcols,npcols, rank, nullity, pr);
   mat_m basis(m.ncols(),nullity);
   bigint one; one=1;
   for (int n=1; n<=nullity; n++)
     basis.set(npcols[n],n,one);
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       basis.set(i,j, -m(r,npcols[j]));
   }
   msubspace ans(basis, npcols, one);
   return ans;
}

msubspace pimage(const mat_m& mat, const bigint& pr)
{
  vec_i p,np;
  long rank, nullity;
  mat_m b = transpose(echmodp(transpose(mat),p,np,rank,nullity,pr));
  msubspace ans(b,p,BIGINT(1));
  return ans;
}

msubspace peigenspace(const mat_m& mat, const bigint& lambda, const bigint& pr)
{
  mat_m m = addscalar(mat,-lambda);
  msubspace ans = pkernel(m,pr);
  return ans;
}

msubspace psubeigenspace(const mat_m& mat, const bigint& l, const msubspace& s, const bigint& pr)
{
  mat_m m = prestrict(mat,s,pr);
  msubspace ss = peigenspace(m, l*(denom(s)),pr);
  msubspace ans = pcombine(s,ss,pr);
  return ans;
}


//Attempts to lift from a mod-p msubspace to a normal Q-msubspace by expressing
//basis as rational using modrat and clearing denominators
//
int lift(const msubspace& s, const bigint& pr, msubspace& ans, int trace)
{
  bigint dd;
  mat_m m;
  if (liftmat(s.basis,pr,m,dd,trace))
    {
      ans = msubspace(m, pivots(s), dd);
      return 1;
    }
  return 0;
}
