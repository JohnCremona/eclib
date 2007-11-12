// sub.cc: implementation of subspace class
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

// Only to be included by subspace.cc

// Inline definitions of member operators and functions:  

subspace::subspace(int n) 
:denom(1),pivots(iota((scalar)n)),basis(idmat((scalar)n))
{}

subspace::subspace(const mat& b, const vec& p, scalar d)
:denom(d),pivots(p),basis(b)
{}

subspace::subspace(const subspace& s) 
:denom(s.denom),pivots(s.pivots),basis(s.basis) 
{}

// destructor -- no need to do anything as componenets have their own
subspace::~subspace() 
{}

// assignment
void subspace::operator=(const subspace& s) 
{
  pivots=s.pivots; 
  basis=s.basis; 
  denom=s.denom;
}

// Definitions of nonmember, nonfriend operators and functions:

subspace combine(const subspace& s1, const subspace& s2)
{ 
  scalar d = s1.denom * s2.denom;
  const mat& b1=s1.basis, b2=s2.basis;
  int nr = b1.nro, nc = b2.nco;
  mat b = b1*b2;
  scalar g=0; long n=nr*nc; scalar* bp=b.entries;
  while ((n--)&&(g!=1)) g=gcd(g,*bp++);
  if(g>1)
    {
      d/=g; bp=b.entries; n=nr*nc; while(n--) (*bp++)/=g;
    }
  vec p = s1.pivots[s2.pivots];
  return subspace(b,p,d);
}
 
//Don't think the following is ever actually used...
mat expressvectors(const mat& m, const subspace& s)
{ vec p = pivots(s);
  long   n = dim(s);
  mat ans(n,ncols(m));
  for (int i=1; i<=n; i++) ans.setrow(i, m.row(p[i]));
  return ans;
}
 
//This one is used a LOT
mat restrict(const mat& m, const subspace& s, int cr)
{ long i,j,k,d = dim(s), n=m.nro;
  if(d==n) return m; // trivial special case, s is whole space
  scalar dd = s.denom;
  mat ans(d,d);
  const mat& sb = s.basis;
  scalar *ap, *a=m.entries, *b=sb.entries, *bp, *c=ans.entries, *cp, *pv=s.pivots.entries;
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
  if(cr) {
//  int check = 1, n = nrows(b);
//  for (i=1; (i<=n) && check; i++)
//  for (j=1; (j<=d) && check; j++)
//   check = (dd*m.row(i)*b.col(j) == b.row(i)*ans.col(j));
    int check = (dd*matmulmodp(m,sb,BIGPRIME) == matmulmodp(sb,ans,BIGPRIME));
    if (!check) cerr<<"Warning from restrict: subspace not invariant!\n";
  }
  return ans;
}
 
subspace kernel(const mat& m1, int method)
{
   long rank, nullity, n, r, i, j;
   scalar d;
   vec pcols,npcols;
   mat m = echelon(m1,pcols,npcols, rank, nullity, d, method);
   int dim = ncols(m);
   mat basis(dim,nullity);
   for (n=1; n<=nullity; n++) basis.set(npcols[n],n,d);
   for (r=1; r<=rank; r++)
   { i = pcols[r];
     for (j=1; j<=nullity; j++) basis.set(i,j, -m(r,npcols[j]));
   }
   subspace ans(basis, npcols, d);
   return ans;
}
 
subspace image(const mat& m, int method)
{
  vec p,np;
  long rank, nullity;
  scalar d;
  mat b = transpose(echelon(transpose(m),p,np,rank,nullity,d,method));
  subspace ans(b,p,d);
  return ans;
}
 
subspace eigenspace(const mat& m1, scalar lambda, int method)
{
  mat m = addscalar(m1,-lambda);
  subspace ans = kernel(m,method);
  return ans;
}
 
subspace subeigenspace(const mat& m1, scalar l, const subspace& s, int method)
{
  mat m = restrict(m1,s);
  subspace ss = eigenspace(m, l*(denom(s)),method);
  subspace ans = combine(s,ss );
  return ans;
}

subspace pcombine(const subspace& s1, const subspace& s2, scalar pr)
{
  scalar   d = s1.denom * s2.denom;  // redundant since both should be 1
  const mat& b1=s1.basis,  b2=s2.basis;
  const mat& b = matmulmodp(b1,b2,pr);
  const vec& p = s1.pivots[s2.pivots];
  return subspace(b,p,d);
}

mat prestrict(const mat& m, const subspace& s, scalar pr, int cr)
{ int i,j,k,d = dim(s), n=m.nro;
  if(d==n) return m; // trivial special case, s is whole space
  scalar dd = s.denom;  // will be 1 if s is a mod-p subspace
  mat ans(d,d);
  const mat& sb = s.basis;
  scalar *ap, *a=m.entries, *b=sb.entries, *bp, *c=ans.entries, *cp, *pv=s.pivots.entries;
  for(i=0; i<d; i++)
    {
      bp=b; k=n; ap=a+n*(pv[i]-1);
      while(k--)
	{
	  cp=c; j=d;
	  while(j--)
	    {
	      *cp += xmodmul(*ap , *bp++, pr);
	      *cp = xmod(*cp, pr);
	      cp++;
	    }
	  ap++;
	}
      cp=c; j=d;
      while(j--)
	{
	  *cp = mod(*cp,pr);
	  cp++;
	}
      c += d;
    }
  if(cr) {
    const mat& left = dd*matmulmodp(m,sb,pr);
    const mat& right = matmulmodp(sb,ans,pr);
    int check = (left==right);
    if (!check) cerr<<"Warning from prestrict: subspace not invariant!\n";
  }
  return ans;
}
 
subspace oldpkernel(const mat& m1, scalar pr)   // using full echmodp
{
   long rank, nullity, n, r, i, j;
   vec pcols,npcols;
   mat m = echmodp(m1,pcols,npcols, rank, nullity, pr);
   int dim = ncols(m);
   mat basis(dim,nullity);
   for (n=1; n<=nullity; n++) basis.set(npcols[n],n,1);
   for (r=1; r<=rank; r++)
   { i = pcols[r];
     for (j=1; j<=nullity; j++) basis.set(i,j, mod(-m(r,npcols[j]),pr));
   }
   subspace ans(basis, npcols, 1);
   return ans;
}

// using echmodp_uptri, with no back-substitution
subspace pkernel(const mat& m1, scalar pr)
{
  long rank, nullity, i, j, jj, t, tt;
  vec pcols,npcols;
  mat m = echmodp_uptri(m1,pcols,npcols, rank, nullity, pr);
  int dim = ncols(m);
  mat basis(dim,nullity);
  for(j=nullity; j>0; j--)
    {
      jj = npcols[j];
      basis(jj,j) = 1;
      for(i=rank; i>0; i--)
        {
          scalar temp = -m(i,jj);
          for(t=rank; t>i; t--)
            {
              tt=pcols[t];
              temp -= xmodmul(m(i,tt),basis(tt,j),pr);
              temp = xmod(temp,pr);
            }
          basis(pcols[i],j) = mod(temp,pr);
        }
    }
  subspace ans(basis, npcols, 1);
  return ans;
}
 
subspace pimage(const mat& m, scalar pr)
{
  vec p,np;
  long rank, nullity;
  const mat& b = transpose(echmodp(transpose(m),p,np,rank,nullity,pr));
  subspace ans(b,p,1);
  return ans;
}
 
subspace peigenspace(const mat& m1, scalar lambda, scalar pr)
{
  const mat& m = addscalar(m1,-lambda);
  subspace ans = pkernel(m,pr);
  return ans;
}

subspace psubeigenspace(const mat& m1, scalar l, const subspace& s, scalar pr)
{
  const mat& m = prestrict(m1,s,pr);
  const subspace& ss = peigenspace(m, l*(denom(s)),pr);
  subspace ans = pcombine(s,ss,pr);
  return ans;
}


//Attempts to lift from a mod-p subspace to a normal Q-subspace by expressing
//basis as rational using modrat and clearing denominators
//
subspace lift(const subspace& s, scalar pr, int trace)
{
  scalar dd;
  const mat& m = liftmat(s.basis,pr,dd,trace);
  subspace ans(m, pivots(s), dd);
  return ans;
}
