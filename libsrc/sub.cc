// sub.cc: implementation of subspace class
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

// Only to be included by subspace.cc

// definitions of member operators and functions:

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
  const mat& b1=s1.basis;
  const mat& b2=s2.basis;
  mat b = b1*b2;
  scalar g = b.content();
  if(g>1)
    {
      d/=g; b/=g;
    }
  vec p = s1.pivots[s2.pivots];
  return subspace(b,p,d);
}

//Don't think the following is ever actually used...
mat expressvectors(const mat& m, const subspace& s)
{ vec p = pivots(s);
  long   n = dim(s);
  mat ans(n,m.ncols());
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

mat restrict_mat(const mat& M, const subspace& S, int cr)
{
  if(dim(S)==M.nro) return M; // trivial special case, s is whole space
  const mat& B = S.basis;
  mat A = rowsubmat(M, S.pivots) * B;

  if(cr) // optional check that S is invariant under M
    {
      int check = (S.denom*matmulmodp(M,B,DEFAULT_MODULUS) == matmulmodp(B,A,DEFAULT_MODULUS));
      if (!check)
        cerr<<"Error in restrict_mat: subspace not invariant!"<<endl;
    }
  return A;
}

subspace kernel(const mat& m1, int method)
{
   long rank, nullity;
   scalar d;
   vec pcols,npcols;
   mat m = echelon(m1,pcols,npcols, rank, nullity, d, method);
   mat basis(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     basis.set(npcols[n],n,d);
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       basis.set(i,j, -m(r,npcols[j]));
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
  mat m = restrict_mat(m1,s);
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

// Same as restrict_mat, but modulo pr
mat prestrict(const mat& M, const subspace& S, scalar pr, int cr)
{
  if(dim(S)==M.nro) return M; // trivial special case, s is whole space
  const mat& B = S.basis;
  mat A = matmulmodp(rowsubmat(M, S.pivots), B, pr);

  if(cr) // optional check that S is invariant under M
    {
      int check = (S.denom*matmulmodp(M,B,pr) == matmulmodp(B,A,pr));
      if (!check)
        cerr<<"Error in prestrict: subspace not invariant!"<<endl;
    }
  return A;
}

subspace oldpkernel(const mat& m1, scalar pr)   // using full echmodp
{
   long rank, nullity;
   vec pcols,npcols;
   mat m = echmodp(m1,pcols,npcols, rank, nullity, pr);
   mat basis(m.ncols(),nullity);
   for (int n=1; n<=nullity; n++)
     basis.set(npcols[n],n,1);
   for (int r=1; r<=rank; r++)
   {
     int i = pcols[r];
     for (int j=1; j<=nullity; j++)
       basis.set(i,j, mod(-m(r,npcols[j]),pr));
   }
   subspace ans(basis, npcols, 1);
   return ans;
}

// using echmodp_uptri, with no back-substitution
subspace pkernel(const mat& m1, scalar pr)
{
  long rank, nullity;
  vec pcols,npcols;
  mat m = echmodp_uptri(m1,pcols,npcols, rank, nullity, pr);
  mat basis(m.ncols(),nullity);
  for(int j=nullity; j>0; j--)
    {
      int jj = npcols[j];
      basis(jj,j) = 1;
      for(int i=rank; i>0; i--)
        {
          scalar temp = -m(i,jj);
          for(int t=rank; t>i; t--)
            {
              int tt=pcols[t];
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
int lift(const subspace& s, scalar pr, subspace& ans, int trace)
{
  scalar dd;
  mat m;
  if (liftmat(s.basis,pr,m,dd,trace))
    {
      ans = subspace(m, pivots(s), dd);
      return 1;
    }
  return 0;
}
