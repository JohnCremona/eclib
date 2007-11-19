// msubspace.cc: implementations of multiprecision subspace class
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
 
//#define CHECK_RESTRICT   // define this to make restrict and prestrict check
                           //  invariance of subspaces.
#include "marith.h"
#include "msubspace.h"

// Definitions of nonmember, nonfriend operators and functions:

msubspace combine(const msubspace& s1, const msubspace& s2)
{ 
  bigint d = s1.denom * s2.denom;
  mat_m b1=s1.basis, b2=s2.basis;
  long nr = b1.nro, nc = b2.nco;
  mat_m b = b1*b2;
  bigint g; long n=nr*nc; bigint* bp=b.entries;
  while ((n--)&&(!is_one(g))) g=gcd(g,*bp++);
  if(!(is_zero(g)||is_one(g)))
    {
      d/=g; bp=b.entries; n=nr*nc; while(n--) (*bp++)/=g;
    }
  vec_i p = s1.pivots[s2.pivots];
  return msubspace(b,p,d);
}
 
//This one is used a LOT
mat_m restrict(const mat_m& m, const msubspace& s)
{ long i,j,k,d = dim(s), n=m.nro;
  bigint dd = s.denom;
  mat_m ans(d,d);
  const mat_m& sb = s.basis;
  bigint *ap, *a=m.entries, *b=sb.entries, *bp, *c=ans.entries, *cp;
  int *pv=s.pivots.entries;
  for(i=0; i<d; i++)
    {
      bp=b; k=n; ap=a+n*(pv[i]-1);
      while(k--)
        {
          cp=c; j=d;
          while(j--)
            {
              *cp++ += *ap * *bp++;
            }
          ap++;
        }
      c += d;
    }
// N.B. The following check is strictly unnecessary and slows it down, 
// but is advisable! 
#ifdef CHECK_RESTRICT
  int check = 1; n = nrows(sb);
  for (i=1; (i<=n) && check; i++)
  for (j=1; (j<=d) && check; j++)
    check = (dd*m.row(i)*sb.col(j) == sb.row(i)*ans.col(j));
//#define BIGPRIME 92681
//int check = (dd*matmulmodp(m,sb,BIGPRIME) == matmulmodp(sb,ans,BIGPRIME));
  if (!check) 
    {
      cout<<"Error in restrict: msubspace not invariant!\n";
      abort();
    }
#endif
  return ans;
}
 
msubspace kernel(const mat_m& mat, int method)
{
   long rank, nullity, n, r, i, j;
   bigint d, zero; zero=0;
   vec_i pcols,npcols;
   mat_m m = echelon(mat,pcols,npcols, rank, nullity, d, method);
   long dim = ncols(m);
   mat_m basis(dim,nullity);
   for (n=1; n<=nullity; n++) basis.set(npcols[n],n,d);
   for (r=1; r<=rank; r++)
   { i = pcols[r];
     for (j=1; j<=nullity; j++) basis.set(i,j, -m(r,npcols[j]));
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
  mat_m m = restrict(mat,s);
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

mat_m prestrict(const mat_m& m, const msubspace& s, const bigint& pr)
{ long i,j,k,d = dim(s), n=m.nro;
  bigint dd = s.denom;  // will be 1 if s is a mod-p msubspace
  mat_m ans(d,d);
  const mat_m& sb = s.basis;
  bigint *ap, *a=m.entries, *b=sb.entries, *bp, *c=ans.entries, *cp;
  int *pv=s.pivots.entries;
  for(i=0; i<d; i++)
    {
      bp=b; k=n; ap=a+n*(pv[i]-1);
      while(k--)
        {
          cp=c; j=d;
          while(j--)
            {
              *cp += mod(*ap * *bp++, pr);
              *cp = mod(*cp, pr);
              cp++;
            }
          ap++;
        }
      c += d;
    }
#ifdef CHECK_RESTRICT
  mat_m& left = matmulmodp(m,sb,pr), right = matmulmodp(sb,ans,pr);
  if(dd!=1) left*=dd;
  int check = (left==right);
  if (!check) 
    {
      cout<<"Error in prestrict: msubspace not invariant!\n";
      abort();
    }
#endif
  return ans;
}
 
msubspace pkernel(const mat_m& mat, const bigint& pr)
{
   long rank, nullity, n, r, i, j;
   vec_i pcols,npcols;
   mat_m m = echmodp(mat,pcols,npcols, rank, nullity, pr);
   long dim = ncols(m);
   mat_m basis(dim,nullity);
   bigint one; one=1;
   for (n=1; n<=nullity; n++) basis.set(npcols[n],n,one);
   for (r=1; r<=rank; r++)
   { i = pcols[r];
     for (j=1; j<=nullity; j++) basis.set(i,j, -m(r,npcols[j]));
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
msubspace lift(const msubspace& s, const bigint& pr, int trace)
{
  bigint modulus=pr,dd,n,d; long nr,nc,nrc;
  int succ,success=1;
  bigint lim=sqrt(pr>>1);
  mat_m m = s.basis; bigint *mp;
  if(trace)
    {
      cout << "Lifting mod-p msubspace.\n basis mat_m mod "<<pr<<" is:\n";
      cout << m;
      cout << "Now lifting back to Q.\n";
      cout << "lim = " << lim << "\n";
    }
  nr = m.nro; nc = m.nco;  nrc = nr*nc; mp=m.entries;
  dd=1;
  while(nrc--)
    {  
      succ = modrat(*mp++,modulus,lim,n,d);
      success = success && succ;
      dd=lcm(d,dd);
    }
  if(!success) 
    cout << "Problems encountered with modrat lifting of msubspace." << endl;
  dd=abs(dd);
  if(trace) cout << "Common denominator = " << dd << "\n";
  nrc=nr*nc; mp=m.entries;
  while(nrc--)
      {
        *mp=mod(dd*(*mp),pr); 
        mp++;
      }
  msubspace ans(m, pivots(s), dd);
  return ans;
}
