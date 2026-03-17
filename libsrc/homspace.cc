// FILE HOMSPACE.CC: Implemention of class homspace
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

#define USE_SMATS // Warning:  no longer testing without this switched on!

#include "eclib/cusp.h"
#include "eclib/homspace.h"
#include "eclib/timer.h"
#include "eclib/polys.h"

svec mat22::operator()(const symb& s, const homspace* h)const
{
  long u=s.ceered(),v=s.deered();
  return h->coords_cd(a*u+c*v,b*u+d*v);
}

vec mat22::operator()(const symb& s, const homspace* h, const mat& bas)const
{
  long u=s.cee(),v=s.dee();
  return h->proj_coords_cd(a*u+c*v,b*u+d*v, bas);
}

string opname(const long& p, const long& n=1)
{
  ostringstream ans;
  ans << (n%p ? "T" : "W") << "_" << p;
  return ans.str();
}

string gmatop::name() const
{
  ostringstream s;
  auto ci=coeffs.begin();
  int first = 1;
  for (auto T: ops)
    {
      scalar c = *ci++;
      if (c!=0)
        {
          if (!first)
            {
              s << "+";
            }
          first = 0;
          if (c!=1)
            s << "[" << c << "]";
          s<< T.name();
        }
    }
  return s.str();
}

matop::matop(long a, long b, long c, long d)
{
  mats.push_back(mat22(a,b,c,d));
  the_name = "generic";
}

matop::matop(long p, long n)
{
  the_name = opname(p, n); // W_p or T_p
  if (p==n)
    {
      mats.push_back(mat22(0,-1,n,0));
      return;
    }
  if ((n%p)==0)   // W involution, 1 term
    {
      long u,v,a,b;
      for (u=1, v=n; v%p==0; v/=p, u*=p) ;
      bezout(u,v,a,b);
      mats.push_back(mat22(u*a,-b,n,u));
      return;
    }
  // else  Hecke operator, p+1 terms
  {
    mats.resize(p+1);
    long j, p2 = p>>1;
    for (j=0; j<p; j++) mats[j] = mat22(1,j-p2,0,p);
    mats[p] = mat22(p,0,0,1);
  }
}

homspace::homspace(long n, scalar mod, int hp, int verbose)
  :symbdata(n), modulus(mod)
{
  plusflag = hp; // plusflag is a member of class level
  // cout<<"In homspace constructor, level="<<n<<", modulus="
  //     <<mod<<", plusflag="<<hp<<endl;
  init_time();
  long i,j,k, ngens=0;
   coordindex.resize(nsymb);
   vector<int> check(nsymb, 0);
   vector<int> gens(nsymb+1);    // N.B. Start of gens array is at 1 not 0
   vector<long> rel(4);

// 2-term relations:

if (plusflag!=0)
  for (j=0; j<nsymb; j++)
  {if (check[j]==0)
   { rel[0]=j;
     if (plusflag==-1)
       rel[1]=rof(j);
     else
       rel[1]=rsof(j);
     rel[2]=sof(j);
     rel[3]=sof(rel[1]);
     if (verbose>1)
       cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" "<<rel[3]<<endl;
     for (k=0; k<4; k++) check[rel[k]]=1;
     if ( (j==rel[2]) || (j==rel[3]) )
         for (k=0; k<4; k++) coordindex[rel[k]]=0;
     else
     {   ngens++;
         gens[ngens] = j;
         if (verbose>1)  cout << "gens["<<ngens<<"]="<<j<<endl;
         coordindex[rel[0]] =  ngens;
         coordindex[rel[1]] =  ngens;
         coordindex[rel[2]] = -ngens;
         coordindex[rel[3]] = -ngens;
     }
     }
   }
if (plusflag==0)
  {for (j=0; j<nsymb; j++)
   {if (check[j]==0)
    {rel[0]=j;
     rel[1]=sof(j);
     check[rel[0]] = check[rel[1]] = 1;
     if (j==rel[1])
         for (k=0; k<2; k++) coordindex[rel[k]]=0;
     else
     {   ngens++;
         gens[ngens] = j;
         coordindex[rel[0]] =  ngens;
         coordindex[rel[1]] = -ngens;
     }
    }
   }
 }


// end of 2-term relations
if (verbose)
{
  cout << "After 2-term relations, ngens = "<<ngens<<"\n";
// Compare with predicted value:
/*
 int nu2=(::divides((long)4,N)?0:1);
 static int nu2table[4] = {0,2,1,0};
 for(i=0; nu2&&(i<npdivs); i++)  nu2 *= nu2table[plist[i]%4];
 int ngens0=(nsymb-nu2)/2;
 cout<<"predicted value of ngens = "<<ngens0;
 if(!plusflag) if(ngens!=ngens0) cout<<" --WRONG!";
 cout<<endl;
*/
if (verbose>1)
  {
 cout << "gens = ";
 for (i=1; i<=ngens; i++) cout << gens[i] << " ";
 cout << "\n";
 cout << "coordindex = " << coordindex << "\n";
}}
//
// 3-term relations

//   long maxnumrel = 20+(2*ngens)/3;
 long maxnumrel = ngens+10;

   if (verbose)
     {
       cout << "ngens     = "<<ngens<<endl;
       cout << "maxnumrel = "<<maxnumrel<<endl;
       cout << "relation matrix has = "<<(maxnumrel*ngens)<<" entries..."<<endl;
     }
   {
#ifdef USE_SMATS
   smat relmat(maxnumrel,ngens);
   svec newrel(ngens);
#else
   mat relmat(maxnumrel,ngens);
   vec newrel(ngens);
#endif
   int numrel = 0;
   long ij; int fix;

   std::fill(check.begin(), check.end(), 0);
   for (k=0; k<nsymb; k++)
     {
       if (check[k]) continue;
       newrel.clear();
       rel[2]=tof(rel[1]=tof(rel[0]=k));
       if (verbose>1)
         cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
       for (j=0; j<3; j++)
         {
           ij = rel[j];
           check[ij] = 1;
           if (plusflag) check[rof(ij)] = 1;
           fix = coordindex[ij];
           if(fix!=0)
             newrel.add(abs(fix),scalar(fix>0?1:-1));
         }
       if(verbose>1)
         cout<<newrel<<"\n";
       if (content(newrel)!=0)
         {
           numrel++;
           make_primitive(newrel);
           if(numrel<=maxnumrel)
             relmat.setrow(numrel,newrel);
           else
             cout<<"Too many 3-term relations (numrel="<<numrel
                 <<", maxnumrel="<<maxnumrel<<")"<<endl;
         }
     }

   if (verbose)
     {
       cout << "Finished 3-term relations: numrel = "<<numrel<<" ( maxnumrel = "<<maxnumrel<<")"<<endl;

// end of 3-term relations

#ifdef USE_SMATS
       cout << "relmat has "<< get_population(relmat)<<" nonzero entries (density = "<<density(relmat)<<")"<<endl;
#endif

       if(verbose>1)
         cout << "relmat = " << relmat.as_mat().slice(numrel,ngens) << endl;
       cout << "Computing kernel..."<<endl;
     }

   vec_i pivs, npivs;
   coord_vecs.resize(ngens+1); // 0'th is unused
#ifdef USE_SMATS
   smat_elim sme(relmat, modulus);
   smat sp;
   liftmat(sme.kernel(npivs,pivs), modulus, sp, denom1);
   dimension = sp.ncols();
   for(i=1; i<=ngens; i++)
     coord_vecs[i]=sp.row(i);
   coord = sp.as_mat();
#else
   subspace sp = kernel(relmat);
   coord = basis(sp);
   dimension = dim(sp);
   pivs = pivots(sp);
   denom1 = denom(sp);
   for(i=1; i<=ngens; i++)
     coord_vecs[i]=svec(coord.row(i));
#endif

   //   cout<<"ngens = "<<ngens<<endl;
   if (verbose)
     {
       cout << "dimension = " << dimension << endl;
       if (verbose>1)
	 {
	   //       cout << "coord:\n" ; coord.output_pretty();
	   cout << "coord_vecs:\n";
           for(i=1; i<=ngens; i++)
             cout<<i<<": "<<coord_vecs[i].as_vec()<<"\n";
	   cout << "pivots = " << pivs <<endl;
	 }
     }
   freegens.resize(dimension);
   if (dimension>0)
   {
     for (i=0; i<dimension; i++) freegens[i] = gens[pivs[i+1]];
     if (verbose>1)
       cout << "freegens: " << freegens << "\n";
   }
  pivs.init();  npivs.init();
   }

   if (verbose) cout << "Number of cusps is " << ncusps << endl;

   cusplist cusps(ncusps, this);

   for (i=0; i<dimension; i++)
     {
       modsym m(symbol(freegens[i]));
       for (j=1; j>-3; j-=2)
	 {
	   rational c = (j==1 ? m.beta() : m.alpha());
           if (plusflag==-1)
	     k = cusps.index_1(c);   //adds automatically if new, ignores if [c]=[-c]
           else
             k = cusps.index(c);   //adds automatically if new
	 }
     }
   long ncusps_seen = cusps.count();
   if (verbose&&(ncusps!=ncusps_seen))
     cerr << "image of delta only sees " << ncusps_seen
          << " cusps out of " << ncusps << endl;
   if(verbose)
     {
       cout << "Number of cusps in image(delta) = " << ncusps_seen << endl;
       if(verbose>1) {cusps.display(); cout<<endl;}
       cout << "About to compute matrix of delta"<<endl;
     }
   mat deltamat=mat(ncusps,dimension);  // could make this sparse

   for (i=0; i<dimension; i++)
     {
       modsym m(symbol(freegens[i]));
       for (j=1; j>-3; j-=2)
	 {
	   rational c = (j==1 ? m.beta() : m.alpha());
           if (plusflag==-1)
             k = cusps.index_1(c);
           else
             k = cusps.index(c);
           if (k>0)
             deltamat(k,i+1) += j;
           if (k<0)
             deltamat(-k,i+1) -= j;
	 }
     }
   if (verbose)
     {
       cout << "delta matrix done: size "<<deltamat.nrows()<<"x"<<deltamat.ncols()<<". "<<endl;
       if(verbose>1)
	 cout<<"deltamat = "<<deltamat<<endl;
       cout << "About to compute kernel of delta"<<endl;
     }

   smat sdeltamat(deltamat);
   kern=kernel(sdeltamat, modulus);
   vec_i pivs, npivs;
   smat sk;
   liftmat(smat_elim(sdeltamat, modulus).kernel(npivs,pivs), modulus,sk,denom2);
   tkernbas = transpose(kern.bas());         // dim(kern) x rank
   deltamat.init(0); // clear space.
   if(verbose>1)
      cout<<"tkernbas = "<<tkernbas.as_mat()<<endl;

   if (verbose) cout << "done "<<endl;

   cuspidal_dimension = dim(kern);
   denom3 = denom1*denom2;

   freemods.resize(dimension);

   if (dimension>0)
   {
        if (verbose>1)  cout << "Freemods:\n";
        for (i=0; i<dimension; i++)
	  {
	    freemods[i] = modsym(symbol(freegens[i])) ;
	    if (verbose>1)
	      {
		cout << i << ": " << freemods[i] << "\n";
	      }
	  }
        if (verbose>1)
        {
	  cout << "Basis of ker(delta):\n";
	  cout << kern.bas().as_mat()<<endl;
	  cout << "pivots: " << pivots(kern) << endl;
        }
   }
   if (verbose) cout << "Finished constructing homspace." << endl;
}

// Extend a dual vector of length dimension to one of length nsymb:
vec homspace::extend_coords(const vec& v)
{
  //  cout<<"Extending vector "<<v<<endl;
  vec ans(nsymb);
  for(int i=1; i<=nsymb; i++)
    {
      int j = coordindex[i-1];
      if (j==0) ans[i] = 0;
      else if (j>0) ans[i] =  v*coord_vecs[j];
      else if (j<0) ans[i] = -v*coord_vecs[-j];
    }
  //  cout<<"returning "<<ans<<endl;
  return ans;
}

// Contract a dual vector of length nsymb to one of length dimension:
vec homspace::contract_coords(const vec& v)
{
  //  cout<<"Contracting vector "<<v<<endl;
  vec ans(dimension);
  int i;
  for(i=1; i<=dimension; i++)
    ans[i] = v[1+freegens[i-1]];
  //  cout<<"returning "<<ans<<endl;
  return ans;
}

svec homspace::coords_from_index(int ind) const
{
 long i= coordindex[ind];
 if (i>0) return  coord_vecs[i];
 if (i<0) return -coord_vecs[-i];
 return zero_coords();
}

vec homspace::proj_coords_from_index(int ind, const mat& bas) const
{
 long i= coordindex[ind];
 if (i>0) return  bas.row(i);
 if (i<0) return -bas.row(-i);
 return vec(bas.ncols());
}

scalar homspace::nfproj_coords_from_index(int ind, const vec& bas) const
{
 long i= coordindex[ind];
 if (i>0) return  bas[i];
 if (i<0) return -bas[-i];
 return scalar(0);
}

svec homspace::coords(const symb& s) const
{
  return coords_from_index(index(s));
}

void homspace::add_coords(svec& v, const symb& s) const
{
  v += coords_from_index(index(s));
}

svec homspace::coords_cd(long c, long d) const
{
  return coords_from_index(index2(c,d));
}

void homspace::add_coords_cd(svec& v, long c, long d) const
{
  v += coords_from_index(index2(c,d));
}

vec homspace::proj_coords_cd(long c, long d, const mat& bas) const
{
  return proj_coords_from_index(index2(c,d), bas);
}

void homspace::add_proj_coords_cd(vec& v, long c, long d, const mat& bas) const
{
  // cout << "In add_proj_coords_cd(v,c,d,bas) with v = "<<v<<", (c:d)=("<<c<<":"<<d<<")"<<endl;
  long n = coordindex[index2(c,d)];
  // cout << "(c:d) --> " << n <<endl;
  if (n>0) v.add_row(bas,n);
  else if (n<0) v.sub_row(bas,-n);
  // cout << "New v = " << v << endl;
}

scalar homspace::nfproj_coords_cd(long c, long d, const vec& bas) const
{
  return nfproj_coords_from_index(index2(c,d), bas);
}

void homspace::add_nfproj_coords_cd(scalar& a, long c, long d, const vec& bas) const
{
  a += nfproj_coords_from_index(index2(c,d), bas);
}

svec homspace::coords(long nn, long dd) const
{
   svec ans = zero_coords();
   add_coords(ans, nn, dd);
   return ans;
}

svec homspace::coords(const modsym& m) const
{
  svec ans = zero_coords();
  add_coords(ans, m);
  return ans;
}

void homspace::add_coords(svec& vv, const modsym& m) const
{
  rational al = m.alpha(), be=m.beta();
  long a=num(be), b=num(al), c=den(be), d=den(al), u, v;
  long de = a*d-b*c;
  if (de<0) {de=-de; b=-b; d=-d;}
  if (de==1)
    {
      add_coords_cd(vv, c, d);
      return;
    }
  // now de>1
  bezout(a,c, u,v); // =1
  long nu = b*u+v*d;
  //
  // now m = M{0,infinity} = U.{nu/de,infinity} where U=[a,-v;c,u], so
  // we find the CF expansion of nu/de.
  //
  long C=c, D=u, r=nu, s=de;
  while (s)
    {
      long t=s; s=mod(r,s); long q=(r-s)/t;  r=-t;
      long e=D; D=-C; C= q*C+e;
      add_coords_cd(vv, -D, C);
   }
}


vec homspace::proj_coords(long nn, long dd, const mat& bas) const
{
   vec ans = vec(bas.ncols());
   add_proj_coords(ans, nn, dd, bas);
   return ans;
}

scalar homspace::nfproj_coords(long nn, long dd, const vec& bas) const
{
  scalar ans(0);
  add_nfproj_coords(ans, nn, dd, bas);
  return ans;
}

void homspace::add_coords(svec& v, long nn, long dd) const
{
   add_coords_cd(v,0,1);
   long c=0, d=1, a=nn, b=dd;
   while (b)
     {
       long f=b; b=mod(a,b); long q=(a-b)/f; a= -f;
       long e=d; d=-c; c=q*c+e;
       add_coords_cd(v,c,d);
     }
}

void homspace::add_proj_coords(vec& v, long nn, long dd, const mat& bas) const
{
  // cout << "In add_proj_coords(v,n,d,bas) with v = " << v << ", n/d = "<<nn<<"/"<<dd<<endl;
  add_proj_coords_cd(v,0,1,bas);
  // cout << "After one step, v = " << v <<endl;
  long c=0, d=1, a=nn, b=dd;
   while (b)
   {
     long f=b; b=mod(a,b); long q=(a-b)/f; a= -f;
     long e=d; d=-c; c=q*c+e;
     add_proj_coords_cd(v,c,d,bas);
     // cout << "After another step, v = " << v <<endl;
   }
}

void homspace::add_nfproj_coords(scalar& aa, long nn, long dd, const vec& bas) const
{
   add_nfproj_coords_cd(aa,0,1,bas);
   long c=0, d=1, a=nn, b=dd;
   while (b)
   {
     long f=b; b=mod(a,b); long q=(a-b)/f; a= -f;
     long e=d; d=-c; c=q*c+e;
     add_nfproj_coords_cd(aa,c,d,bas);
   }
}

svec homspace::applyop(const matop& T, const rational& q) const
{ svec ans(dimension);
  long i=T.size();
  while (i--) add_coords(ans,T[i](q));
  return ans;
}

svec homspace::applyop(const matop& T, const modsym& m) const
{ svec ans(dimension);
  long i=T.size();
  while (i--)  ans += coords(T[i](m));
  return ans;
}

vec homspace::applyop_proj(const matop& T, const rational& q, const mat& bas) const
{ vec ans(bas.ncols());
  long i=T.size();
  while (i--) add_proj_coords(ans,T[i](q), bas);
  return ans;
}

mat homspace::calcop_restricted(const matop& T,	const subspace& s, int dual, int display) const
{
  long d=dim(s);
  mat m(d,dimension);
  for (long j=0; j<d; j++)
     {
       long jj = pivots(s)[j+1]-1;
       svec colj = applyop(T,freemods[jj]);
       m.setrow(j+1,colj.as_vec());
     }
  m = (smat(m)*smat(basis(s))).as_mat();
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display)
    {
      cout << "Matrix of " << T.name() << " = ";
      if (dimension>1) cout << "\n";
      m.output_pretty();
    }
  return m;
}

smat homspace::s_calcop_restricted(const matop& T, const ssubspace& s, int dual, int display) const
{
  long d=dim(s);
  smat m(d,dimension);
  for (long j=1; j<=d; j++)
     {
       long jj = pivots(s)[j];
       svec colj = applyop(T,freemods[jj-1]);
       m.setrow(j,colj);
     }
  //  m = m*basis(s);
  m = mult_mod_p(m,basis(s), modulus);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display)
    {
      cout << "Matrix of " << T.name() << " = ";
      if (dimension>1) cout << "\n";
      cout<<m<<endl;
    }
  return m;
}

mat homspace::calcop(const matop& T, int cuspidal, int dual, int display) const
{
  mat m(dimension,dimension);
  for (long j=0; j<dimension; j++)
     { svec colj = applyop(T,freemods[j]);
       m.setcol(j+1,colj.as_vec());
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) {  m=transpose(m);}
  if (display)
    {
      cout << "Matrix of " << T.name() << " = ";
      if (dimension>1) cout << "\n";
      m.output_pretty();
    }
  return m;
}

mat homspace::calcop(const gmatop& T, int cuspidal, int dual, int display) const
{
  if(display)
    cout<<"Computing " << T.name() <<"...";
  mat m(dimension,dimension);
  auto ci = T.coeffs.begin();
  auto Ti = T.ops.begin();
  while (ci!=T.coeffs.end())
    {
      scalar c = *ci++;
      if (c !=0 )
        {
          mat mi = calcop(*Ti, cuspidal, dual, display);
          if (c!=1)
            mi *= c;
          m += mi;
        }
      ++Ti;
    }
  if (display)
    {
      cout<<"done."<<endl;
      cout << "Matrix of " << T.name() << " = " << m;
      if (dimension>1) cout << endl;
    }
  return m;
}

//#define DEBUG_CHARPOL
ZZX homspace::charpoly(const matop& T, int cuspidal) const
{
  ZZ den = to_ZZ(cuspidal? denom3: denom1);
#ifdef DEBUG_CHARPOL
  cout << "Computing charpoly of " << T.name();
  if (cuspidal) cout << " on cuspidal subspace";
  else  cout << " on full space";
  mat M = calcop(T,cuspidal,0);
  cout << "Matrix: "; output_flat_matrix(M); cout << endl;
  cout << "denominator: " << den << endl;
  cout << "scaled char poly: " << scaled_charpoly(mat_to_mat_ZZ(M),den) << endl;
#endif
  return scaled_charpoly(mat_to_mat_ZZ(calcop(T,cuspidal,0)), den);
}

ZZX homspace::charpoly(const gmatop& T, int cuspidal) const
{
  ZZ den = to_ZZ(cuspidal? denom3: denom1);
  return scaled_charpoly(mat_to_mat_ZZ(calcop(T,cuspidal,0)), den);
}

vec homspace::calcop_col(int j, const matop& T, int display) const
{
  // j counts from 1
  vec colj = applyop(T,freemods[j-1]).as_vec();
  if (display)
    {
      cout << "Image of "<<j<<"-th generator under " << T.name() << " = "
           << colj << endl;
    }
  return colj;
}

mat homspace::calcop_cols(const vec_i& jlist, const matop& T, int display) const
{
  int d = dim(jlist);
  mat m(d,dimension);
  for (int i=1; i<=d; i++)
    m.setrow(i, applyop(T,freemods[jlist[i]-1]).as_vec());
  return m;
}

svec homspace::s_calcop_col(int j, const matop& T,
                            int display) const
{
  // j counts from 1
  svec colj = applyop(T,freemods[j-1]);
  if (display)
    {
      cout << "Image of "<<j<<"-th generator under " << T.name() << " = "
           << colj.as_vec() << endl;
    }
  return colj;
}

smat homspace::s_calcop_cols(const vec_i& jlist, const matop& T, int display) const
{
  int d = dim(jlist);
  smat m(d,dimension);
  for (int i=1; i<=d; i++)
    m.setrow(i, applyop(T,freemods[jlist[i]-1]));
  return m;
}

smat homspace::s_calcop(const matop& T, int cuspidal, int dual, int display) const
{
  smat m(dimension,dimension);
  for (long j=0; j<dimension; j++)
     { svec colj = applyop(T,freemods[j]);
       m.setrow(j+1,colj);
     }
  if(cuspidal)
    {
      m = restrict_mat(transpose(m),kern);
      if(dual) m = transpose(m);
    }
  else if(!dual) {m=transpose(m);}
  if (display)
    {
      cout << "Matrix of " << T.name() << " = ";
      if (dimension>1) cout << "\n";
      cout<<m.as_mat();
    }
  return m;
}

/* NOTE: heckeop actually only returns the Hecke operator for p not dividing
   the level.  For p that divide the level it actually computes the Atkin-Lehner
   operator.
*/

mat homspace::newheckeop(long p, int cuspidal, int dual, int display) const
{
  if(::divides(p,N)) return wop(p,cuspidal,dual,display);
  matop hmats(p); // constructs H-matrices
  long j, nmats=hmats.size();
  svec colj(dimension);    mat m(dimension,dimension);
  for (j=0; j<dimension; j++)
    {  symb s = symbol(freegens[j]);
       colj = hmats[0](s,this);
       colj+= hmats[1](s,this);
       for(long i=2; i<nmats; i++) colj+= (hmats[i](s,this));
       m.setcol(j+1,colj.as_vec());
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m=transpose(m);
  if (display)
    {
      cout << opname(p,N) << " = ";
      if (dimension>1) cout << "\n";
      m.output_pretty();
    }
  return m;
}

mat homspace::conj(int cuspidal, int dual, int display) const
{
 mat m(dimension,dimension);
 for (long j=1; j<=dimension; j++)
 {  symb s = symbol(freegens[j-1]);
    svec colj   =  coords_cd(-s.cee(),s.dee());
    m.setcol(j,colj.as_vec());
 }
 if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
 if(dual) m=transpose(m);
 if (display) cout << "Matrix of conjugation = " << m;
 return m;
}

vec homspace::conj_col(int j, int display) const
{
  // j counts from 1
  symb s = symbol(freegens[j-1]);
  vec colj   =  coords_cd(-s.cee(),s.dee()).as_vec();
  if (display) cout << "Column "<<j<<" of matrix of conjugation = " << colj << endl;
  return colj;
}

mat homspace::conj_cols(const vec_i& jlist, int display) const
{
  int d = dim(jlist);
  mat m(d,dimension);
  for (int i=1; i<=d; i++)
    {
      symb s = symbol(freegens[jlist[i]-1]);
      m.setrow(i, coords_cd(-s.cee(),s.dee()).as_vec());
    }
  return m;
}

svec homspace::s_conj_col(int j, int display) const
{
  // j counts from 1
  symb s = symbol(freegens[j-1]);
  svec colj   =  coords_cd(-s.cee(),s.dee());
  if (display) cout << "Column "<<j<<" of matrix of conjugation = " << colj.as_vec() << endl;
  return colj;
}

smat homspace::s_conj_cols(const vec_i& jlist, int display) const
{
  int d = dim(jlist);
  smat m(d,dimension);
  int i;
  for (i=1; i<=d; i++)
    {
      symb s = symbol(freegens[jlist[i]-1]);
      m.setrow(i, coords_cd(-s.cee(),s.dee()));
    }
  return m;
}

smat homspace::s_conj(int cuspidal, int dual, int display) const
{
 smat m(dimension,dimension);
 for (long j=1; j<=dimension; j++)
 {  symb s = symbol(freegens[j-1]);
    svec colj   =  coords_cd(-s.cee(),s.dee());
    m.setrow(j,colj);
 }
 if(cuspidal)
   {
     m = restrict_mat(transpose(m),kern);
     if(dual) m = transpose(m);
   }
 else if(!dual) {m=transpose(m);}
 if (display) cout << "Matrix of conjugation = " << m;
 return m;
}

// Computes matrix of conjugation restricted to a DUAL subspace; here
// the ambient space of s must be H_1(-;cusps) of dimension dimension
mat homspace::conj_restricted(const subspace& s,
			      int dual, int display) const
{
  long d = dim(s);
  mat m(d,dimension);
  for (long j=1; j<=d; j++)
    {
      long jj=pivots(s)[j];
      symb sy = symbol(freegens[jj-1]);
      svec colj   =  coords_cd(-sy.cee(),sy.dee());
      m.setrow(j,colj.as_vec());
    }
  m = matmulmodp(m,basis(s), modulus);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) cout << "Matrix of conjugation = " << m;
  return m;
}

smat homspace::s_conj_restricted(const ssubspace& s,
				 int dual, int display) const
{
  long d = dim(s);
  smat m(d,dimension);
  for (long j=1; j<=d; j++)
    {
      long jj=pivots(s)[j];
      symb sy = symbol(freegens[jj-1]);
      svec colj   =  coords_cd(-sy.cee(),sy.dee());
      m.setrow(j,colj);
    }
  //  cout<<"m = "<<m<<" = "<<m.as_mat()<<endl;
  m = mult_mod_p(m,basis(s), modulus);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) cout << "Matrix of conjugation = " << m.as_mat();
  return m;
}

mat homspace::fricke(int cuspidal, int dual, int display) const
{
  return calcop(matop(N,N),cuspidal,dual,display);
}

long homspace::op_prime(int i)  // the i'th operator prime for Tp or Wq
{
#ifdef NEW_OP_ORDER
  long p = prime_number(i+1);
  //  cout<<"opmat("<<i<<") using p="<<p<<endl;
#else
  long p = primelist[i];
#endif
  return p;
}

mat homspace::opmat(int i, int dual, int v)
{
  if(i==-1) return conj(0,dual,v);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::opmat(): called with i = " << i << endl;
      return mat(dimension);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    cout<<"Computing " << opname(p,N) << "..."<<flush;
  mat ans = heckeop(p,0,dual,v); // Automatically chooses W or T
  if(v)
    cout<<"done."<<endl;
  return ans;
}

vec homspace::opmat_col(int i, int j, int v)
{
  if(i==-1) return conj_col(j,v);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::opmat_col(): called with i = " << i << endl;
      return vec(dimension);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    {
      cout<<"Computing col "<<j<<" of " << opname(p,N) << "..."<<flush;
    }
  vec ans = heckeop_col(p,j,0); // Automatically chooses W or T
  if(v)
    {
      cout<<"done."<<endl;
    }
  return ans;
}

mat homspace::opmat_cols(int i, const vec_i& jlist, int v)
{
  if(i==-1) return conj_cols(jlist,v);
  int d = dim(jlist);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::opmat_cols(): called with i = " << i << endl;
      return mat(d,dimension);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    {
      cout<<"Computing "<<d<<" cols of " << opname(p,N) << "..."<<flush;
    }
  mat ans = heckeop_cols(p,jlist,0); // Automatically chooses W or T
  if(v)
    {
      cout<<"done."<<endl;
    }
  return ans;
}

svec homspace::s_opmat_col(int i, int j, int v)
{
  if(i==-1) return s_conj_col(j,v);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::opmat_col(): called with i = " << i << endl;
      return svec(dimension);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    {
      cout<<"Computing col "<<j<<" of " << opname(p,N) << "..."<<flush;
    }
  svec ans = s_heckeop_col(p,j,0); // Automatically chooses W or T
  if(v)
    {
      cout<<"done."<<endl;
    }
  return ans;
}

smat homspace::s_opmat_cols(int i, const vec_i& jlist, int v)
{
  if(i==-1) return s_conj_cols(jlist,v);
  int d = dim(jlist);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::opmat_col(): called with i = " << i << endl;
      return smat(d,dimension);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    cout<<"Computing "<<d<<" cols of " << opname(p,N) << "..."<<flush;
  smat ans = s_heckeop_cols(p,jlist,0); // Automatically chooses W or T
  if(v)
    cout<<"done."<<endl;
  return ans;
}

smat homspace::s_opmat(int i, int dual, int v)
{
  if(i==-1) return s_conj(0,dual,v);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::s_opmat(): called with i = " << i << endl;
      return smat(dimension);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    {
      cout<<"Computing " << opname(p,N) << "..."<<flush;
      smat ans = s_heckeop(p,0,dual,v>1); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop(p,0,dual,0); // Automatically chooses W or T
}

mat homspace::opmat_restricted(int i, const subspace& s, int dual, int v)
{
  if(i==-1) return conj_restricted(s,dual,v);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::opmat_restricted(): called with i = "
	  << i << endl;
      return mat(dim(s));  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    {
      cout<<"Computing " << opname(p,N) <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
      mat ans = heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
}

smat homspace::s_opmat_restricted(int i, const ssubspace& s, int dual, int v)
{
  if(i==-1) return s_conj_restricted(s,dual,v);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::s_opmat_restricted(): called with i = "
	  << i << endl;
      return smat(dim(s));  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    {
      cout<<"Computing " << opname(p,N) <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
      smat ans = s_heckeop_restricted(p,s,dual,v); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
}

#ifdef OLD_EIG_ORDER

static long*pm1={1,-1};

vector<long> T_eigrange(long p)
{
  vector<long> ans;
  ans.push_back(0);
  long aplim=1;
  while (aplim*aplim<=4*p) aplim++; aplim--;
  for(long ap=1; ap<=aplim; ap++)
    {
      ans.push_back(ap);
      ans.push_back(-ap);
    }
  return ans;
}

#else  // new eig order, in strict numerical order

static long pm1[]={-1,1};

vector<long> T_eigrange(long p)
{
  long aplim=3, four_p=p<<2;
  while (aplim*aplim<=four_p) aplim++;
  aplim--;
  vector<long> ans(1+2*aplim);
  iota(ans.begin(),ans.end(),-aplim);
  return ans;
}
#endif

vector<long> homspace::eigrange(long i)
{
  if((i<0)||(i>=nap)) return vector<long>(0);  // shouldn't happen
  long p = op_prime(i);
  if(::divides(p,N))  return vector<long>(pm1,pm1+2);
  return T_eigrange(p);
}

vec homspace::maninvector(long p) const
{
  svec tvec = coords(0,p);             // =0, but sets the right length.
  if (plusflag!=-1)
    {
      if (p==2)
	add_coords(tvec,1,2);
      else
	{
	  long p2=(p-1)>>1;
	  for (int i=1; i<=p2; i++) { add_coords(tvec,i,p); }
	  if(plusflag)
	    tvec *= scalar(2);
	  else
	    for (int i=1; i<=p2; i++) { add_coords(tvec,-i,p); }
	}
    }
  return tvec.as_vec();
}

vec homspace::manintwist(long p) const
{
 svec sum = coords(0,p);                   // =0, but sets the right length.
 for (long i=1; i<p; i++) sum += scalar(legendre(i,p))*coords(i,p);
 return sum.as_vec();
}

matop::matop(long p)
{
  //    case 31: nmats = 106; break;
  //    case 37: nmats = 128; break;
  //    case 41: nmats = 146; break;
  //    case 43: nmats = 154; break;
  //    case 47: nmats = 170; break;

  the_name = opname(p);
  switch (p) {
  case 2:
    mats.resize(4);
    mats[0]=mat22(2,0,0,1); mats[1]=mat22(2,1,0,1);
    mats[2]=mat22(1,0,1,2); mats[3]=mat22(1,0,0,2);
    return;
  case 3:
    mats.resize(6);
    mats[0]=mat22(1,0,0,3);
    mats[1]=mat22(3,1,0,1);
    mats[2]=mat22(1,0,1,3);
    mats[3]=mat22(3,0,0,1);
    mats[4]=mat22(3,-1,0,1);
    mats[5]=mat22(-1,0,1,-3);
    return;
  case 5:
    mats.resize(12);
    mats[0]=mat22(1,0,0,5);
    mats[1]=mat22(5,2,0,1);
    mats[2]=mat22(2,1,1,3);
    mats[3]=mat22(1,0,3,5);
    mats[4]=mat22(5,1,0,1);
    mats[5]=mat22(1,0,1,5);
    mats[6]=mat22(5,0,0,1);
    mats[7]=mat22(5,-1,0,1);
    mats[8]=mat22(-1,0,1,-5);
    mats[9]=mat22(5,-2,0,1);
    mats[10]=mat22(-2,1,1,-3);
    mats[11]=mat22(1,0,-3,5);
    return;
  case 7:
    mats.resize(18);
    mats[0]=mat22(1,0,0,7);
    mats[1]=mat22(7,3,0,1);
    mats[2]=mat22(3,-1,1,2);
    mats[3]=mat22(-1,0,2,-7);
    mats[4]=mat22(7,2,0,1);
    mats[5]=mat22(2,1,1,4);
    mats[6]=mat22(1,0,4,7);
    mats[7]=mat22(7,1,0,1);
    mats[8]=mat22(1,0,1,7);
    mats[9]=mat22(7,0,0,1);
    mats[10]=mat22(7,-1,0,1);
    mats[11]=mat22(-1,0,1,-7);
    mats[12]=mat22(7,-2,0,1);
    mats[13]=mat22(-2,1,1,-4);
    mats[14]=mat22(1,0,-4,7);
    mats[15]=mat22(7,-3,0,1);
    mats[16]=mat22(-3,-1,1,-2);
    mats[17]=mat22(-1,0,-2,-7);
    return;
  case 11:
    mats.resize(30);
    mats[0]=mat22(1,0,0,11);
    mats[1]=mat22(11,5,0,1);
    mats[2]=mat22(5,-1,1,2);
    mats[3]=mat22(-1,0,2,-11);
    mats[4]=mat22(11,4,0,1);
    mats[5]=mat22(4,1,1,3);
    mats[6]=mat22(1,0,3,11);
    mats[7]=mat22(11,3,0,1);
    mats[8]=mat22(3,1,1,4);
    mats[9]=mat22(1,0,4,11);
    mats[10]=mat22(11,2,0,1);
    mats[11]=mat22(2,1,1,6);
    mats[12]=mat22(1,0,6,11);
    mats[13]=mat22(11,1,0,1);
    mats[14]=mat22(1,0,1,11);
    mats[15]=mat22(11,0,0,1);
    mats[16]=mat22(11,-1,0,1);
    mats[17]=mat22(-1,0,1,-11);
    mats[18]=mat22(11,-2,0,1);
    mats[19]=mat22(-2,1,1,-6);
    mats[20]=mat22(1,0,-6,11);
    mats[21]=mat22(11,-3,0,1);
    mats[22]=mat22(-3,1,1,-4);
    mats[23]=mat22(1,0,-4,11);
    mats[24]=mat22(11,-4,0,1);
    mats[25]=mat22(-4,1,1,-3);
    mats[26]=mat22(1,0,-3,11);
    mats[27]=mat22(11,-5,0,1);
    mats[28]=mat22(-5,-1,1,-2);
    mats[29]=mat22(-1,0,-2,-11);
    return;
  case 13:
    mats.resize(38);
    mats[0]=mat22(1,0,0,13);
    mats[1]=mat22(13,6,0,1);
    mats[2]=mat22(6,-1,1,2);
    mats[3]=mat22(-1,0,2,-13);
    mats[4]=mat22(13,5,0,1);
    mats[5]=mat22(5,2,1,3);
    mats[6]=mat22(2,-1,3,5);
    mats[7]=mat22(-1,0,5,-13);
    mats[8]=mat22(13,4,0,1);
    mats[9]=mat22(4,-1,1,3);
    mats[10]=mat22(-1,0,3,-13);
    mats[11]=mat22(13,3,0,1);
    mats[12]=mat22(3,-1,1,4);
    mats[13]=mat22(-1,0,4,-13);
    mats[14]=mat22(13,2,0,1);
    mats[15]=mat22(2,1,1,7);
    mats[16]=mat22(1,0,7,13);
    mats[17]=mat22(13,1,0,1);
    mats[18]=mat22(1,0,1,13);
    mats[19]=mat22(13,0,0,1);
    mats[20]=mat22(13,-1,0,1);
    mats[21]=mat22(-1,0,1,-13);
    mats[22]=mat22(13,-2,0,1);
    mats[23]=mat22(-2,1,1,-7);
    mats[24]=mat22(1,0,-7,13);
    mats[25]=mat22(13,-3,0,1);
    mats[26]=mat22(-3,-1,1,-4);
    mats[27]=mat22(-1,0,-4,-13);
    mats[28]=mat22(13,-4,0,1);
    mats[29]=mat22(-4,-1,1,-3);
    mats[30]=mat22(-1,0,-3,-13);
    mats[31]=mat22(13,-5,0,1);
    mats[32]=mat22(-5,2,1,-3);
    mats[33]=mat22(2,-1,-3,8);
    mats[34]=mat22(-1,0,8,-13);
    mats[35]=mat22(13,-6,0,1);
    mats[36]=mat22(-6,-1,1,-2);
    mats[37]=mat22(-1,0,-2,-13);
    return;
  case 17:
    mats.resize(52);
    mats[0]=mat22(1,0,0,17);
    mats[1]=mat22(17,8,0,1);
    mats[2]=mat22(8,-1,1,2);
    mats[3]=mat22(-1,0,2,-17);
    mats[4]=mat22(17,7,0,1);
    mats[5]=mat22(7,-3,1,2);
    mats[6]=mat22(-3,-1,2,-5);
    mats[7]=mat22(-1,0,-5,-17);
    mats[8]=mat22(17,6,0,1);
    mats[9]=mat22(6,1,1,3);
    mats[10]=mat22(1,0,3,17);
    mats[11]=mat22(17,5,0,1);
    mats[12]=mat22(5,-2,1,3);
    mats[13]=mat22(-2,-1,3,-7);
    mats[14]=mat22(-1,0,-7,-17);
    mats[15]=mat22(17,4,0,1);
    mats[16]=mat22(4,-1,1,4);
    mats[17]=mat22(-1,0,4,-17);
    mats[18]=mat22(17,3,0,1);
    mats[19]=mat22(3,1,1,6);
    mats[20]=mat22(1,0,6,17);
    mats[21]=mat22(17,2,0,1);
    mats[22]=mat22(2,1,1,9);
    mats[23]=mat22(1,0,9,17);
    mats[24]=mat22(17,1,0,1);
    mats[25]=mat22(1,0,1,17);
    mats[26]=mat22(17,0,0,1);
    mats[27]=mat22(17,-1,0,1);
    mats[28]=mat22(-1,0,1,-17);
    mats[29]=mat22(17,-2,0,1);
    mats[30]=mat22(-2,1,1,-9);
    mats[31]=mat22(1,0,-9,17);
    mats[32]=mat22(17,-3,0,1);
    mats[33]=mat22(-3,1,1,-6);
    mats[34]=mat22(1,0,-6,17);
    mats[35]=mat22(17,-4,0,1);
    mats[36]=mat22(-4,-1,1,-4);
    mats[37]=mat22(-1,0,-4,-17);
    mats[38]=mat22(17,-5,0,1);
    mats[39]=mat22(-5,-2,1,-3);
    mats[40]=mat22(-2,-1,-3,-10);
    mats[41]=mat22(-1,0,-10,-17);
    mats[42]=mat22(17,-6,0,1);
    mats[43]=mat22(-6,1,1,-3);
    mats[44]=mat22(1,0,-3,17);
    mats[45]=mat22(17,-7,0,1);
    mats[46]=mat22(-7,-3,1,-2);
    mats[47]=mat22(-3,1,-2,-5);
    mats[48]=mat22(1,0,-5,17);
    mats[49]=mat22(17,-8,0,1);
    mats[50]=mat22(-8,-1,1,-2);
    mats[51]=mat22(-1,0,-2,-17);
    return;
  case 19:
    mats.resize(58);
    mats[0]=mat22(1,0,0,19);
    mats[1]=mat22(19,9,0,1);
    mats[2]=mat22(9,-1,1,2);
    mats[3]=mat22(-1,0,2,-19);
    mats[4]=mat22(19,8,0,1);
    mats[5]=mat22(8,-3,1,2);
    mats[6]=mat22(-3,1,2,-7);
    mats[7]=mat22(1,0,-7,19);
    mats[8]=mat22(19,7,0,1);
    mats[9]=mat22(7,2,1,3);
    mats[10]=mat22(2,-1,3,8);
    mats[11]=mat22(-1,0,8,-19);
    mats[12]=mat22(19,6,0,1);
    mats[13]=mat22(6,-1,1,3);
    mats[14]=mat22(-1,0,3,-19);
    mats[15]=mat22(19,5,0,1);
    mats[16]=mat22(5,1,1,4);
    mats[17]=mat22(1,0,4,19);
    mats[18]=mat22(19,4,0,1);
    mats[19]=mat22(4,1,1,5);
    mats[20]=mat22(1,0,5,19);
    mats[21]=mat22(19,3,0,1);
    mats[22]=mat22(3,-1,1,6);
    mats[23]=mat22(-1,0,6,-19);
    mats[24]=mat22(19,2,0,1);
    mats[25]=mat22(2,1,1,10);
    mats[26]=mat22(1,0,10,19);
    mats[27]=mat22(19,1,0,1);
    mats[28]=mat22(1,0,1,19);
    mats[29]=mat22(19,0,0,1);
    mats[30]=mat22(19,-1,0,1);
    mats[31]=mat22(-1,0,1,-19);
    mats[32]=mat22(19,-2,0,1);
    mats[33]=mat22(-2,1,1,-10);
    mats[34]=mat22(1,0,-10,19);
    mats[35]=mat22(19,-3,0,1);
    mats[36]=mat22(-3,-1,1,-6);
    mats[37]=mat22(-1,0,-6,-19);
    mats[38]=mat22(19,-4,0,1);
    mats[39]=mat22(-4,1,1,-5);
    mats[40]=mat22(1,0,-5,19);
    mats[41]=mat22(19,-5,0,1);
    mats[42]=mat22(-5,1,1,-4);
    mats[43]=mat22(1,0,-4,19);
    mats[44]=mat22(19,-6,0,1);
    mats[45]=mat22(-6,-1,1,-3);
    mats[46]=mat22(-1,0,-3,-19);
    mats[47]=mat22(19,-7,0,1);
    mats[48]=mat22(-7,2,1,-3);
    mats[49]=mat22(2,-1,-3,11);
    mats[50]=mat22(-1,0,11,-19);
    mats[51]=mat22(19,-8,0,1);
    mats[52]=mat22(-8,-3,1,-2);
    mats[53]=mat22(-3,-1,-2,-7);
    mats[54]=mat22(-1,0,-7,-19);
    mats[55]=mat22(19,-9,0,1);
    mats[56]=mat22(-9,-1,1,-2);
    mats[57]=mat22(-1,0,-2,-19);
    return;
  case 23:
    mats.resize(74);
    mats[0]=mat22(1,0,0,23);
    mats[1]=mat22(23,11,0,1);
    mats[2]=mat22(11,-1,1,2);
    mats[3]=mat22(-1,0,2,-23);
    mats[4]=mat22(23,10,0,1);
    mats[5]=mat22(10,-3,1,2);
    mats[6]=mat22(-3,-1,2,-7);
    mats[7]=mat22(-1,0,-7,-23);
    mats[8]=mat22(23,9,0,1);
    mats[9]=mat22(9,4,1,3);
    mats[10]=mat22(4,-1,3,5);
    mats[11]=mat22(-1,0,5,-23);
    mats[12]=mat22(23,8,0,1);
    mats[13]=mat22(8,1,1,3);
    mats[14]=mat22(1,0,3,23);
    mats[15]=mat22(23,7,0,1);
    mats[16]=mat22(7,-2,1,3);
    mats[17]=mat22(-2,-1,3,-10);
    mats[18]=mat22(-1,0,-10,-23);
    mats[19]=mat22(23,6,0,1);
    mats[20]=mat22(6,1,1,4);
    mats[21]=mat22(1,0,4,23);
    mats[22]=mat22(23,5,0,1);
    mats[23]=mat22(5,2,1,5);
    mats[24]=mat22(2,-1,5,9);
    mats[25]=mat22(-1,0,9,-23);
    mats[26]=mat22(23,4,0,1);
    mats[27]=mat22(4,1,1,6);
    mats[28]=mat22(1,0,6,23);
    mats[29]=mat22(23,3,0,1);
    mats[30]=mat22(3,1,1,8);
    mats[31]=mat22(1,0,8,23);
    mats[32]=mat22(23,2,0,1);
    mats[33]=mat22(2,1,1,12);
    mats[34]=mat22(1,0,12,23);
    mats[35]=mat22(23,1,0,1);
    mats[36]=mat22(1,0,1,23);
    mats[37]=mat22(23,0,0,1);
    mats[38]=mat22(23,-1,0,1);
    mats[39]=mat22(-1,0,1,-23);
    mats[40]=mat22(23,-2,0,1);
    mats[41]=mat22(-2,1,1,-12);
    mats[42]=mat22(1,0,-12,23);
    mats[43]=mat22(23,-3,0,1);
    mats[44]=mat22(-3,1,1,-8);
    mats[45]=mat22(1,0,-8,23);
    mats[46]=mat22(23,-4,0,1);
    mats[47]=mat22(-4,1,1,-6);
    mats[48]=mat22(1,0,-6,23);
    mats[49]=mat22(23,-5,0,1);
    mats[50]=mat22(-5,2,1,-5);
    mats[51]=mat22(2,-1,-5,14);
    mats[52]=mat22(-1,0,14,-23);
    mats[53]=mat22(23,-6,0,1);
    mats[54]=mat22(-6,1,1,-4);
    mats[55]=mat22(1,0,-4,23);
    mats[56]=mat22(23,-7,0,1);
    mats[57]=mat22(-7,-2,1,-3);
    mats[58]=mat22(-2,-1,-3,-13);
    mats[59]=mat22(-1,0,-13,-23);
    mats[60]=mat22(23,-8,0,1);
    mats[61]=mat22(-8,1,1,-3);
    mats[62]=mat22(1,0,-3,23);
    mats[63]=mat22(23,-9,0,1);
    mats[64]=mat22(-9,4,1,-3);
    mats[65]=mat22(4,1,-3,5);
    mats[66]=mat22(1,0,5,23);
    mats[67]=mat22(23,-10,0,1);
    mats[68]=mat22(-10,-3,1,-2);
    mats[69]=mat22(-3,1,-2,-7);
    mats[70]=mat22(1,0,-7,23);
    mats[71]=mat22(23,-11,0,1);
    mats[72]=mat22(-11,-1,1,-2);
    mats[73]=mat22(-1,0,-2,-23);
    return;
  case 29:
    mats.resize(96);
    mats[0]=mat22(1,0,0,29);
    mats[1]=mat22(29,14,0,1);
    mats[2]=mat22(14,-1,1,2);
    mats[3]=mat22(-1,0,2,-29);
    mats[4]=mat22(29,13,0,1);
    mats[5]=mat22(13,-3,1,2);
    mats[6]=mat22(-3,-1,2,-9);
    mats[7]=mat22(-1,0,-9,-29);
    mats[8]=mat22(29,12,0,1);
    mats[9]=mat22(12,-5,1,2);
    mats[10]=mat22(-5,-2,2,-5);
    mats[11]=mat22(-2,1,-5,-12);
    mats[12]=mat22(1,0,-12,29);
    mats[13]=mat22(29,11,0,1);
    mats[14]=mat22(11,4,1,3);
    mats[15]=mat22(4,1,3,8);
    mats[16]=mat22(1,0,8,29);
    mats[17]=mat22(29,10,0,1);
    mats[18]=mat22(10,1,1,3);
    mats[19]=mat22(1,0,3,29);
    mats[20]=mat22(29,9,0,1);
    mats[21]=mat22(9,-2,1,3);
    mats[22]=mat22(-2,-1,3,-13);
    mats[23]=mat22(-1,0,-13,-29);
    mats[24]=mat22(29,8,0,1);
    mats[25]=mat22(8,3,1,4);
    mats[26]=mat22(3,1,4,11);
    mats[27]=mat22(1,0,11,29);
    mats[28]=mat22(29,7,0,1);
    mats[29]=mat22(7,-1,1,4);
    mats[30]=mat22(-1,0,4,-29);
    mats[31]=mat22(29,6,0,1);
    mats[32]=mat22(6,1,1,5);
    mats[33]=mat22(1,0,5,29);
    mats[34]=mat22(29,5,0,1);
    mats[35]=mat22(5,1,1,6);
    mats[36]=mat22(1,0,6,29);
    mats[37]=mat22(29,4,0,1);
    mats[38]=mat22(4,-1,1,7);
    mats[39]=mat22(-1,0,7,-29);
    mats[40]=mat22(29,3,0,1);
    mats[41]=mat22(3,1,1,10);
    mats[42]=mat22(1,0,10,29);
    mats[43]=mat22(29,2,0,1);
    mats[44]=mat22(2,1,1,15);
    mats[45]=mat22(1,0,15,29);
    mats[46]=mat22(29,1,0,1);
    mats[47]=mat22(1,0,1,29);
    mats[48]=mat22(29,0,0,1);
    mats[49]=mat22(29,-1,0,1);
    mats[50]=mat22(-1,0,1,-29);
    mats[51]=mat22(29,-2,0,1);
    mats[52]=mat22(-2,1,1,-15);
    mats[53]=mat22(1,0,-15,29);
    mats[54]=mat22(29,-3,0,1);
    mats[55]=mat22(-3,1,1,-10);
    mats[56]=mat22(1,0,-10,29);
    mats[57]=mat22(29,-4,0,1);
    mats[58]=mat22(-4,-1,1,-7);
    mats[59]=mat22(-1,0,-7,-29);
    mats[60]=mat22(29,-5,0,1);
    mats[61]=mat22(-5,1,1,-6);
    mats[62]=mat22(1,0,-6,29);
    mats[63]=mat22(29,-6,0,1);
    mats[64]=mat22(-6,1,1,-5);
    mats[65]=mat22(1,0,-5,29);
    mats[66]=mat22(29,-7,0,1);
    mats[67]=mat22(-7,-1,1,-4);
    mats[68]=mat22(-1,0,-4,-29);
    mats[69]=mat22(29,-8,0,1);
    mats[70]=mat22(-8,3,1,-4);
    mats[71]=mat22(3,-1,-4,11);
    mats[72]=mat22(-1,0,11,-29);
    mats[73]=mat22(29,-9,0,1);
    mats[74]=mat22(-9,-2,1,-3);
    mats[75]=mat22(-2,-1,-3,-16);
    mats[76]=mat22(-1,0,-16,-29);
    mats[77]=mat22(29,-10,0,1);
    mats[78]=mat22(-10,1,1,-3);
    mats[79]=mat22(1,0,-3,29);
    mats[80]=mat22(29,-11,0,1);
    mats[81]=mat22(-11,4,1,-3);
    mats[82]=mat22(4,-1,-3,8);
    mats[83]=mat22(-1,0,8,-29);
    mats[84]=mat22(29,-12,0,1);
    mats[85]=mat22(-12,-5,1,-2);
    mats[86]=mat22(-5,2,-2,-5);
    mats[87]=mat22(2,1,-5,12);
    mats[88]=mat22(1,0,12,29);
    mats[89]=mat22(29,-13,0,1);
    mats[90]=mat22(-13,-3,1,-2);
    mats[91]=mat22(-3,1,-2,-9);
    mats[92]=mat22(1,0,-9,29);
    mats[93]=mat22(29,-14,0,1);
    mats[94]=mat22(-14,-1,1,-2);
    mats[95]=mat22(-1,0,-2,-29);
    return;
  default:
    long p2=(p-1)>>1;
    long sl, sg, x1, x2, x3, y1, y2, y3, a, b, c, q, r;
    mats.push_back(mat22(1,0,0,p));
    mats.push_back(mat22(p,0,0,1));
    sl = -2;
    for(sg=1; sg>sl; sg-=2) // i.e. sg = +1 then -1
      for(r=1; r<=p2; r++)
	{
	  x1=p; x2=-sg*r; y1=0; y2=1; a=-p; b=sg*r;
	  mats.push_back(mat22(x1,x2,y1,y2));
	  while(b!=0)
	    {
	      c=mod(a,b); q=(a-c)/b;
	      x3=q*x2-x1; y3=q*y2-y1;
	      a=-b; b=c; x1=x2; x2=x3; y1=y2; y2=y3;
	      mats.push_back(mat22(x1,x2,y1,y2));
	    }
	}
  }
}

// Functions for caching homspaces, full Hecke polynomials and new Hecke polynomials
// Keys are strings of the form Nlabel (for homspace) or Nlabel-Plabel (for Hecke polynomials)

string Nkey(const long& N)
{
  return to_string(N);
}

string NPkey(const long& N, const long& p)
{
  return to_string(N) + ":" + to_string(p);
}

string NTkey(const long& N, const matop& T)
{
  return to_string(N) + ":" + T.name();
}

// identical code to previous
string NTkey(const long& N, const gmatop& T)
{
  return to_string(N) + ":" + T.name();
}

// cache of homspaces keyed by level
map<string,homspace*> H1_dict;

// cache of operator matrices keyed by level and opname, not restricted to cuspidal
map<string, mat> full_mat_dict;

// cache of operator charpolys keyed by level and opname, not restricted to cuspidal
map<string, ZZX> poly_dict;
// cache of operator charpolys keyed by level and opname, restricted to cuspidal
map<string, ZZX> cuspidal_poly_dict;
// cache of operator new charpolys keyed by level and opname, not restricted to cuspidal
map<string, ZZX> new_poly_dict;
// cache of operator new charpolys keyed by level and opname, restricted to cuspidal
map<string, ZZX> new_cuspidal_poly_dict;

void output_poly_dict(ostream& os, map<string, ZZX> D)
{
  for (auto key_pol: D)
    os<<key_pol.first<<" "<< str(key_pol.second) << endl;
}

map<string, ZZX> input_poly_dict(istream& is)
{
  map<string, ZZX> D;
  string key;
  ZZX poly;
  while (!is.eof())
    {
      is >> key >> poly;
      D[key] = poly;
    }
  return D;
}

//#define DEBUG_GET_HOMSPACE
// Return cached (pointer to) homspace at level N, computing and and
// caching it first, for N>0
homspace* get_homspace(const long& N, scalar mod)
{
  if (N<1)
    {
#ifdef DEBUG_GET_HOMSPACE
      cout << "Not computing homspace at level " << N << " -- should be positive" << endl;
#endif
      return NULL;
    }

  string Nlabel = to_string(N);
  if (H1_dict.find(Nlabel) != H1_dict.end())
    {
#ifdef DEBUG_GET_HOMSPACE
      cout << "Getting homspace at level " << N << " from cache " << endl;
#endif
      return H1_dict[Nlabel];
    }

#ifdef DEBUG_GET_HOMSPACE
  cout << "Computing homspace at level " << N << endl;
#endif
  homspace* H = new homspace(N, mod, 1, 0); // sign=+1, verbose=0
#ifdef DEBUG_GET_HOMSPACE
  cout << "Caching homspace at level " << N << endl;
  cout << "dim = " << H->h1dim() << ", denom = " << H-> h1denom() << endl;
  cout << "cdim = " << H->h1cuspdim() << ", cdenom = " << H-> h1cdenom() << endl;
#endif
  H1_dict[Nlabel] = H;
  return H;
}

//#define DEBUG_GET_FULL_MAT

// Key is label(N)-T.name()
// Value is matrix of T on the full space (not restricted to cuspidal subspace)
mat get_full_mat(const long& N,  const matop& T, const scalar& mod)
{
  if (N<2)
    return mat(0,0);
  string NT = NTkey(N,T);
#ifdef DEBUG_GET_FULL_MAT
  cout << "In get_full_mat() with matop key " << NT << endl;
#endif
  if (full_mat_dict.find(NT) != full_mat_dict.end())
    {
#ifdef DEBUG_GET_FULL_MAT
      cout << "key " << NT << " is in full_mat_dict, returning cached matrix" << endl;
#endif
      return full_mat_dict[NT];
    }
#ifdef DEBUG_GET_FULL_MAT
  cout << "key " << NT << " not in full_mat_dict, computing matrix" << endl;
#endif
  homspace* H = get_homspace(N, mod);
  mat M = H->calcop(T,0,0); // cuspidal=0, dual=0
  full_mat_dict[NT] = M;
#ifdef DEBUG_GET_FULL_MAT
  cout << "caching and returning matrix of " << NT << endl;
  output_flat_matrix(M);
  cout << endl;
#endif
  return M;
}

// Key is label(N)-T.name()
// Value is matrix of T on the full space (not restricted to cuspidal subspace)
mat get_full_mat(const long& N,  const gmatop& T, const scalar& mod)
{
  if (N<2)
    return mat(0,0);
  string NT = NTkey(N,T);
#ifdef DEBUG_GET_FULL_MAT
  cout << "In get_full_mat() with gmatop key " << NT << endl;
#endif
  if (full_mat_dict.find(NT) != full_mat_dict.end())
    {
#ifdef DEBUG_GET_FULL_MAT
      cout << "key " << NT << " is in full_mat_dict, returning cached matrix" << endl;
#endif
      return full_mat_dict[NT];
    }
#ifdef DEBUG_GET_FULL_MAT
  cout << "key " << NT << " not in full_mat_dict, computing matrix" << endl;
#endif
  int d = get_homspace(N, mod)->h1dim();
  mat M(d,d);
  if (d)
    {
      auto ci = T.coeffs.begin();
      auto Ti = T.ops.begin();
      while (ci!=T.coeffs.end())
        {
          scalar c = *ci++;
          if (c !=0 )
            {
              mat Mi = get_full_mat(N, *Ti, mod);
              if (c!=1)
                Mi *= c;
              M += Mi;
            }
          ++Ti;
        }
    }
  // We do not need to add to the dict if this gmatop consist of a
  // single matop, since that will have been done in the loop.
  if (full_mat_dict.find(NT) == full_mat_dict.end())
    full_mat_dict[NT] = M;
#ifdef DEBUG_GET_FULL_MAT
  cout << "caching and returning matrix of " << NT << endl;
#endif
  return M;
}

//#define DEBUG_GET_POLY

// from either poly_dict, cuspidal_poly_dict depending on cuspidal flag
ZZX get_poly(const long& N,  const gmatop& T, int cuspidal, const scalar& mod)
{
  if (N<2)
    {
      ZZX pol;
      set(pol); // sets to 1
      return pol;
    }

  string NT = NTkey(N,T);
#ifdef DEBUG_GET_POLY
  cout << "In get_poly(), N = " << N << ", T = " << T.name()
       << ", cuspidal = " << cuspidal
       << ", mod = " << mod
       <<endl;
  cout << "key = " << NT << endl;
#endif
  auto poly_cache = (cuspidal?
                     cuspidal_poly_dict:
                     poly_dict);
  auto new_poly_cache = (cuspidal?
                         new_cuspidal_poly_dict:
                         new_poly_dict);
  if (poly_cache.find(NT) != poly_cache.end())
    {
#ifdef DEBUG_GET_POLY
      cout << "key is in cache, returning " << str(poly_cache[NT]) << endl;
#endif
      return poly_cache[NT];
    }

  homspace* H = get_homspace(N, mod);
  scalar den = (cuspidal? H->h1cdenom() :H->h1denom());
#ifdef DEBUG_GET_POLY
  cout << "Homspace for level " << N << " obtained from get_homspace()" << endl;
  cout << "den =  " << den << endl;
  cout << "About to call get_full_mat()"<< endl;
#endif

  mat M = get_full_mat(N, T, mod);  // dimension x dimension
#ifdef DEBUG_GET_POLY
  cout << "Full matrix M of size " << M.nrows() << " obtained from get_full_mat()" << endl;
  cout << "den =  " << den << endl;
  output_flat_matrix(M);  cout << endl;
#endif
  if (cuspidal)
    {
      // H->kern is the subspace ker(delta) of the full space
      M = restrict_mat(smat(M),H->kern).as_mat();  // cuspidal_dimension x cuspidal_dimension
#ifdef DEBUG_GET_POLY
      cout << "Cuspidal case" << endl;
      cout << "Restricted M to cuspidal subspace"
           << " (which has denominator " << den << "):" << endl;
      output_flat_matrix(M);
      cout << endl;
#endif
    }
  ZZX full_poly =  scaled_charpoly(mat_to_mat_ZZ(M), to_ZZ(den));
  if (T.is_simple())
    {
      poly_cache[NT] = full_poly;
      if (deg(full_poly)==0)
        new_poly_cache[NT] = full_poly;
    }
  return full_poly;
}

// from either new_poly_dict or new_cuspidal_poly_dict depending on cuspidal flag
ZZX get_new_poly(const long& N, const gmatop& T, int cuspidal, const scalar& mod)
{
  string NT = NTkey(N,T);
  auto new_poly_cache = (cuspidal?
                         new_cuspidal_poly_dict:
                         new_poly_dict);
  if (new_poly_cache.find(NT) != new_poly_cache.end())
    return new_poly_cache[NT];

  ZZX new_poly = get_poly(N, T, cuspidal, mod);
  if (deg(new_poly)==0)
    {
      new_poly_cache[NT] = new_poly;
      return new_poly;
    }
  vector<long> DD = posdivs(N);
  for( auto D : DD)
    {
      if (D==N)
        continue;
      ZZX new_poly_D = get_new_poly(D, T, cuspidal, mod);
      if (deg(new_poly_D)==0)
        continue;
      long M = N/D;
      int m = posdivs(M).size();
      for (int i=0; i<m; i++)
        {
          //essentially new_poly /= new_poly_D // but checking divisibility
          ZZX quo, rem;
          DivRem(quo, rem, new_poly, new_poly_D);
          if (IsZero(rem))
            new_poly = quo;
          else
            {
              cout << "Problem in get_new_poly("<<NT<<"), D="<<D<<endl;
              cout << "Dividing " << str(new_poly) << " by " << str(new_poly_D)
                   << " gives quotient " << str(quo) <<", remainder "<< str(rem) << endl;
              cout << "Old multiplicities are smaller than expected."<<endl;
              exit(1);
            }
        }
      if (deg(new_poly)==0) // nothing left, new dimension must be 0
        break;
    } // end of loop over divisors
  if (T.is_simple())
    new_poly_cache[NT] = new_poly;
  return new_poly;
}

// Return true iff T's new poly is squarefree and coprime to its old poly
int test_splitting_operator(const long& N, const gmatop& T, const scalar& mod, int verbose)
{
  if (verbose)
    cout << "Testing " << T.name() << "..." << flush;
  ZZX f_new = get_new_poly(N, T, 1, mod); // cuspidal=1, triv_char=0
  if (!IsSquareFree(f_new))
    {
      if (verbose>1)
        cout << "\n NO: new Hecke polynomial "<<str(f_new)<<" is not squarefree" << endl;
      return 0;
    }
  ZZX f_full = get_poly(N, T, 0, mod); // cuspidal=0, triv_char=0
  ZZX f_old = f_full / f_new;
  if (!AreCoprime(f_new, f_old))
    {
      if (verbose>1)
        cout << "\n NO: new Hecke polynomial "<<str(f_new)
             <<" is not coprime to old Hecke polynomial "<<str(f_old)<<endl
             <<" (full polynomial is "<<str(f_full)<<")"<<endl;
      return 0;
    }
  if (verbose>1)
    cout << "\n YES: new Hecke polynomial is squarefree and coprime to old Hecke polynomial" << endl;
  return 1;
}
