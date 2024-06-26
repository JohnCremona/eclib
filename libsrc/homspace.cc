// FILE HOMSPACE.CC: Implemention of class homspace
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

#define USE_SMATS // Warning:  no longer testing without this switched on!

#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/timer.h>

#ifdef USE_SMATS
#include <eclib/smatrix_elim.h>
#endif

const string W_opname("W");
const string T_opname("T");

matop::matop(long a, long b, long c, long d)
{
  mats.push_back(mat22(a,b,c,d));
}

matop::matop(long p, long n)
{
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
 
homspace::homspace(long n, int hp, int hcusp, int verbose) :symbdata(n)
{ 
  init_time();
   plusflag=hp;
   cuspidal=hcusp;
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
{cout << "After 2-term relations, ngens = "<<ngens<<"\n";
// Compare with predicted value:
/*
 int nu2=(::divides((long)4,modulus)?0:1);
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
   smat relmat((int)maxnumrel,(int)ngens);
   svec newrel(ngens);
#else
   mat relmat(maxnumrel,ngens);
   vec newrel(ngens);
   if (verbose) cout << "successfully allocated "<<endl;
#endif
   int numrel = 0;
   long ij; int fix;

   std::fill(check.begin(), check.end(), 0);
   for (k=0; k<nsymb; k++)
     {
     if (check[k]==0)
   {
#ifdef USE_SMATS
     newrel.clear();
#else
     for (j=1; j<=ngens; j++) newrel[j]=0;
#endif
     rel[2]=tof(rel[1]=tof(rel[0]=k));
     if (verbose>1)   cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
     for (j=0; j<3; j++)
     {ij = rel[j];
      check[ij] = 1;
      if (plusflag) check[rof(ij)] = 1;
      fix = coordindex[ij];
#ifdef USE_SMATS
      if(fix) newrel.add(abs(fix),(fix>0?1:-1));
#else
      if (verbose>1)  cout << fix << " ";
      if (fix!=0) newrel[abs(fix)] += sign(fix);
#endif
     }
#ifdef USE_SMATS
     if(verbose>1) cout<<newrel<<"\n";
#else
     if(verbose>1) cout<<"\n";
#endif
#ifdef USE_SMATS
     if(newrel.size()!=0) 
       {
	 numrel++;
	 make_primitive(newrel);
	 if(numrel<=maxnumrel)
	   relmat.setrow(numrel,newrel);
	 else 
	   cout<<"Too many 3-term relations (numrel="<<numrel
	       <<", maxnumrel="<<maxnumrel<<")"<<endl;
       }
#else
     if (verbose>1)  cout << newrel << "\n";
     long h = content(newrel);
     if (h!=0)
     {  if (h>1) newrel/=h;
	if(numrel<maxnumrel)
	  relmat.setrow(++numrel,newrel);
	else cout<<"Too many 3-term relations!"<<endl;
     }
#endif
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

   vec pivs, npivs;
#ifdef USE_SMATS
   smat_elim sme(relmat,MODULUS);
   //relmat=smat(0,0); // clear space
   int d1;
   start_time();
   smat sp;
   liftmat(sme.kernel(npivs,pivs),MODULUS,sp,d1);
   stop_time();
   if (verbose) {cout<<"time to compute kernel = "; show_time(); cout<<endl;}
   denom1=d1;
   if(verbose>1)
     {
       cout << "kernel of relmat = " << sp.as_mat() << endl;
       cout << "pivots = "<<pivs << endl;
       cout << "denom = "<<d1 << endl;
     }
   rk = sp.ncols();
   coord_vecs.resize(ngens+1); // 0'th is unused
   for(i=1; i<=ngens; i++)
     coord_vecs[i]=sp.row(i);
   //sp=smat(0,0); // clear space
#else
   subspace sp = kernel(relmat);
   rk = dim(sp);
   coord = basis(sp);
   pivs = pivots(sp);
   denom1 = denom(sp);
   for(i=1; i<=ngens; i++) 
     coord_vecs[i]=svec(coord.row(i));
   sp.clear(); 
   relmat.init(); newrel.init();
#endif

   //   cout<<"ngens = "<<ngens<<endl;
   if (verbose) 
     {
       cout << "rk = " << rk << endl;
       if (verbose>1) 
	 {
	   //       cout << "coord:\n" ; coord.output_pretty();
	   cout << "coord_vecs:\n";
           for(i=1; i<=ngens; i++) 
             cout<<i<<": "<<coord_vecs[i].as_vec()<<"\n";
	   cout << "pivots = " << pivs <<endl;
	 }
     }
   freegens.resize(rk);
   if (rk>0)
   {
     for (i=0; i<rk; i++) freegens[i] = gens[pivs[i+1]];
     if (verbose>1)
       cout << "freegens: " << freegens << "\n";
   }
  pivs.init();  npivs.init();
   }

   // Compute the number of cusps
   long maxncusps = 0;
   std::for_each(dlist.begin(), dlist.end(),
                 [&maxncusps, this] (const long& d)
                 {
                   long dd = ::gcd(d,modulus/d); // compute phi(dd):
                   std::for_each(plist.begin(), plist.end(),
                                 [&dd] (const long& p) {if ((dd%p)==0) dd=dd*(p-1)/p;});
                   maxncusps += dd;
                 });
   if (verbose) cout << "Number of cusps is "<<maxncusps<<endl;

   cusplist cusps(maxncusps, this);

   for (i=0; i<rk; i++)
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
   ncusps=cusps.count();

   if(verbose)
     {
       cout << "ncusps = " << ncusps << endl;
       if(verbose>1) {cusps.display(); cout<<endl;}
       cout << "About to compute matrix of delta"<<endl;
     }
   mat deltamat=mat(ncusps,rk);  // should make this sparse

   for (i=0; i<rk; i++)
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
   kern=kernel(sdeltamat, MODULUS);
   vec pivs, npivs;
   int d2;
   smat sk;
   liftmat(smat_elim(sdeltamat, MODULUS).kernel(npivs,pivs),MODULUS,sk,d2);
   denom2=d2;
   tkernbas = transpose(kern.bas());         // dim(kern) x rank
   deltamat.init(0); // clear space.
   if(verbose>1)
      cout<<"tkernbas = "<<tkernbas.as_mat()<<endl;

   if (verbose) cout << "done "<<endl;

   if(cuspidal)
     dimension = dim(kern);
   else
     dimension = rk;
   //   denom2 = 1;
   denom3 = denom1*denom2;

   freemods.resize(rk);
   needed.resize(rk);

   if (dimension>0)
   {
        if (verbose>1)  cout << "Freemods:\n";
        for (i=0; i<rk; i++)
	  {
	    freemods[i] = modsym(symbol(freegens[i])) ;
	    needed[i]   =  (cuspidal? ! trivial(kern.bas().row(i+1).as_vec())
			              : 1);
	    if (verbose>1)
	      {
		cout << i << ": " << freemods[i];
		if (!needed[i]) cout << " (not needed)";
		cout << "\n";
	      }
	  }
        if ((verbose>1)&&cuspidal)
        {
	  cout << "Basis of ker(delta):\n";
	  cout << kern.bas().as_mat()<<endl;
	  cout << "pivots: " << pivots(kern) << endl;
        }
   }
   if (verbose) cout << "Finished constructing homspace." << endl;
}

// Extend a dual vector of length rk to one of length nsymb:
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

  // Contract a dual vector of length nsymb to one of length rk: 
vec homspace::contract_coords(const vec& v)
{
  //  cout<<"Contracting vector "<<v<<endl;
  vec ans(rk);
  int i;
  for(i=1; i<=rk; i++)
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

vec homspace::proj_coords_from_index(int ind, const mat& m) const
{
 long i= coordindex[ind];
 if (i>0) return  m.row(i);
 if (i<0) return -m.row(-i);
 return vec(m.ncols());
}

long homspace::nfproj_coords_from_index(int ind, const vec& bas) const
{
 long i= coordindex[ind];
 if (i>0) return  bas[i];
 if (i<0) return -bas[-i];
 return 0;
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

vec homspace::proj_coords_cd(long c, long d, const mat& m) const 
{
  return proj_coords_from_index(index2(c,d), m);
}

void homspace::add_proj_coords_cd(vec& v, long c, long d, const mat& m) const 
{
  long n = coordindex[index2(c,d)];
  if (n>0) v.add_row(m,n);
  else if (n<0) v.sub_row(m,-n);
}

long homspace::nfproj_coords_cd(long c, long d, const vec& bas) const 
{
  return nfproj_coords_from_index(index2(c,d), bas);
}

void homspace::add_nfproj_coords_cd(long& a, long c, long d, const vec& bas) const 
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


vec homspace::proj_coords(long nn, long dd, const mat& m) const
{
   vec ans = vec(m.ncols());
   add_proj_coords(ans, nn, dd, m);
   return ans;
}

long homspace::nfproj_coords(long nn, long dd, const vec& bas) const
{
  long ans = 0;
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

void homspace::add_proj_coords(vec& v, long nn, long dd, const mat& m) const
{
   add_proj_coords_cd(v,0,1,m);
   long c=0, d=1, a=nn, b=dd;
   while (b)
   {
     long f=b; b=mod(a,b); long q=(a-b)/f; a= -f;
     long e=d; d=-c; c=q*c+e;
     add_proj_coords_cd(v,c,d,m);
   }
}

void homspace::add_nfproj_coords(long& aa, long nn, long dd, const vec& bas) const
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

svec homspace::applyop(const matop& mlist, const rational& q) const
{ svec ans(rk);
  long i=mlist.size();
  while (i--) add_coords(ans,mlist[i](q));
  return ans;
}

svec homspace::applyop(const matop& mlist, const modsym& m) const
{ svec ans(rk);
  long i=mlist.size();
  while (i--)  ans += coords((mlist[i])(m));
  return ans;
}

// copy of routine from ../procs/xsplit.cc:
mat sparse_restrict(const mat& m, const subspace& s);
smat restrict_mat(const smat& m, const subspace& s);

mat homspace::calcop_restricted(string opname, long p, const matop& mlist, 
				const subspace& s, int dual, int display) const
{
  long d=dim(s);
  mat m(d,rk);
  for (long j=0; j<d; j++) 
     { 
       long jj = pivots(s)[j+1]-1;
       svec colj = applyop(mlist,freemods[jj]);
       m.setrow(j+1,colj.as_vec());
     }
  m = (smat(m)*smat(basis(s))).as_mat();
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) 
    {
      cout << "Matrix of " << opname << "(" << p << ") = ";
      if (dimension>1) cout << "\n";
      m.output_pretty();
    }
  return m;
}

smat homspace::s_calcop_restricted(string opname, long p, const matop& mlist, 
				   const ssubspace& s, int dual, int display) const
{
  long d=dim(s);
  smat m(d,rk);
  for (long j=1; j<=d; j++) 
     { 
       long jj = pivots(s)[j];
       svec colj = applyop(mlist,freemods[jj-1]);
       m.setrow(j,colj);
     }
  //  m = m*basis(s);
  m = mult_mod_p(m,basis(s),MODULUS);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) 
    {
      cout << "Matrix of " << opname << "(" << p << ") = ";
      if (dimension>1) cout << "\n";
      cout<<m<<endl;
    }
  return m;
}

mat homspace::calcop(string opname, long p, const matop& mlist, 
		     int dual, int display) const
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { svec colj = applyop(mlist,freemods[j]);
       m.setcol(j+1,colj.as_vec());
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) {  m=transpose(m);}
  if (display) 
    {
      cout << "Matrix of " << opname << "(" << p << ") = ";
      if (dimension>1) cout << "\n";
      m.output_pretty();
    }
  return m;
}

vec homspace::calcop_col(string opname, long p, int j, const matop& mlist,
                         int display) const
{
  // j counts from 1
  vec colj = applyop(mlist,freemods[j-1]).as_vec();
  if (display)
    {
      cout << "Image of "<<j<<"-th generator under " << opname << "(" << p << ") = "
           << colj << endl;
    }
  return colj;
}

mat homspace::calcop_cols(string opname, long p, const vec& jlist, const matop& mlist,
                         int display) const
{
  int i, d = dim(jlist);
  mat m(d,rk);
  for (i=1; i<=d; i++)
    {
      m.setrow(i, applyop(mlist,freemods[jlist[i]-1]).as_vec());
    }
  return m;
}

svec homspace::s_calcop_col(string opname, long p, int j, const matop& mlist,
                         int display) const
{
  // j counts from 1
  svec colj = applyop(mlist,freemods[j-1]);
  if (display)
    {
      cout << "Image of "<<j<<"-th generator under " << opname << "(" << p << ") = "
           << colj.as_vec() << endl;
    }
  return colj;
}

smat homspace::s_calcop_cols(string opname, long p, const vec& jlist, const matop& mlist,
                         int display) const
{
  int i, d = dim(jlist);
  smat m(d,rk);
  for (i=1; i<=d; i++)
    {
      m.setrow(i, applyop(mlist,freemods[jlist[i]-1]));
    }
  return m;
}

smat homspace::s_calcop(string opname, long p, const matop& mlist, 
			int dual, int display) const
{
  smat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { svec colj = applyop(mlist,freemods[j]);
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
      cout << "Matrix of " << opname << "(" << p << ") = ";
      if (dimension>1) cout << "\n";
      cout<<m.as_mat();
    }
  return m;
}

/* NOTE: heckeop actually only returns the Hecke operator for p not dividing
   the level.  For p that divide the level it actually computes the Atkin-Lehner
   operator. 
*/

mat homspace::heckeop(long p, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = ((modulus%p) ? T_opname : W_opname);
 return calcop(name,p,matlist,dual,display);
}
 
smat homspace::s_heckeop(long p, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = ((modulus%p) ? T_opname : W_opname);
 return s_calcop(name,p,matlist,dual,display);
}
 
vec homspace::heckeop_col(long p, int j, int display) const
{
 matop matlist(p,modulus);
 string name = ((modulus%p) ? T_opname : W_opname);
 return calcop_col(name,p,j,matlist,display);
}

mat homspace::heckeop_cols(long p, const vec& jlist, int display) const
{
 matop matlist(p,modulus);
 string name = ((modulus%p) ? T_opname : W_opname);
 return calcop_cols(name,p,jlist,matlist,display);
}

svec homspace::s_heckeop_col(long p, int j, int display) const
{
 matop matlist(p,modulus);
 string name = ((modulus%p) ? T_opname : W_opname);
 return s_calcop_col(name,p,j,matlist,display);
}

smat homspace::s_heckeop_cols(long p, const vec& jlist, int display) const
{
 matop matlist(p,modulus);
 string name = ((modulus%p) ? T_opname : W_opname);
 return s_calcop_cols(name,p,jlist,matlist,display);
}

mat homspace::heckeop_restricted(long p,const subspace& s, 
				 int dual, int display) const
{
 matop matlist(p,modulus);
 string name = ((modulus%p) ? T_opname : W_opname);
 return calcop_restricted(name,p,matlist,s,dual,display);
}
 
smat homspace::s_heckeop_restricted(long p,const ssubspace& s, 
				    int dual, int display) const
{
 matop matlist(p,modulus);
 string name = ((modulus%p) ? T_opname : W_opname);
 return s_calcop_restricted(name,p,matlist,s,dual,display);
}
 
mat homspace::newheckeop(long p, int dual, int display) const
{
  if(::divides(p,modulus)) return wop(p,display);
  matop hmats(p); // constructs H-matrices
  long j, nmats=hmats.size();
  svec colj(rk);    mat m(rk,rk); 
  for (j=0; j<rk; j++) if (needed[j])
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
      cout << "Matrix of T(" << p << ") = ";
      if (dimension>1) cout << "\n";
      m.output_pretty();
    }
  return m;
}

mat homspace::wop(long q, int dual, int display) const
{
 matop matlist(q,modulus);
 return calcop(W_opname,q,matlist,dual,display);
}
 
smat homspace::s_wop(long q, int dual, int display) const
{
 matop matlist(q,modulus);
 return s_calcop(W_opname,q,matlist,dual,display);
}
 
mat homspace::conj(int dual, int display) const
{
 mat m(rk,rk);
 for (long j=1; j<=rk; j++) if (needed[j-1])
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

mat homspace::conj_cols(const vec& jlist, int display) const
{
  int d = dim(jlist);
  mat m(d,rk);
  int i;
  for (i=1; i<=d; i++)
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

smat homspace::s_conj_cols(const vec& jlist, int display) const
{
  int d = dim(jlist);
  smat m(d,rk);
  int i;
  for (i=1; i<=d; i++)
    {
      symb s = symbol(freegens[jlist[i]-1]);
      m.setrow(i, coords_cd(-s.cee(),s.dee()));
    }
  return m;
}

smat homspace::s_conj(int dual, int display) const
{
 smat m(rk,rk);
 for (long j=1; j<=rk; j++) if (needed[j-1])
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
// the ambient space of s must be H_1(-;cusps) of dimension rk
mat homspace::conj_restricted(const subspace& s, 
			      int dual, int display) const
{
  long d = dim(s);
  mat m(d,rk);
  for (long j=1; j<=d; j++) 
    {  
      long jj=pivots(s)[j];
      symb sy = symbol(freegens[jj-1]);
      svec colj   =  coords_cd(-sy.cee(),sy.dee());
      m.setrow(j,colj.as_vec());
    }
  m = matmulmodp(m,basis(s),MODULUS);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) cout << "Matrix of conjugation = " << m;
  return m;
}

smat homspace::s_conj_restricted(const ssubspace& s, 
				 int dual, int display) const
{
  long d = dim(s);
  smat m(d,rk);
  for (long j=1; j<=d; j++) 
    {  
      long jj=pivots(s)[j];
      symb sy = symbol(freegens[jj-1]);
      svec colj   =  coords_cd(-sy.cee(),sy.dee());
      m.setrow(j,colj);
    }
  //  cout<<"m = "<<m<<" = "<<m.as_mat()<<endl;
  m = mult_mod_p(m,basis(s),MODULUS);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) cout << "Matrix of conjugation = " << m.as_mat();
  return m;
}

mat homspace::fricke(int dual, int display) const
{
 matop frickelist(modulus,modulus);
 return calcop(W_opname,modulus,frickelist,dual,display);
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
  if(i==-1) return conj(dual,v);
  if((i<0)||(i>=nap)) 
    {
      cerr<<"Error in homspace::opmat(): called with i = " << i << endl;
      return mat(dimension);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v) 
    {
      cout<<"Computing " << ((::divides(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")..."<<flush;
      mat ans = heckeop(p,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return heckeop(p,dual,0); // Automatically chooses W or T
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
      cout<<"Computing col "<<j<<" of " << ((::divides(p,modulus)) ? W_opname : T_opname)
          <<"("<<p<<")..."<<flush;
    }
  vec ans = heckeop_col(p,j,0); // Automatically chooses W or T
  if(v)
    {
      cout<<"done."<<endl;
    }
  return ans;
}

mat homspace::opmat_cols(int i, const vec& jlist, int v)
{
  if(i==-1) return conj_cols(jlist,v);
  int d = dim(jlist);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::opmat_cols(): called with i = " << i << endl;
      return mat(d,rk);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    {
      cout<<"Computing "<<d<<" cols of " << ((::divides(p,modulus)) ? W_opname : T_opname)
          <<"("<<p<<")..."<<flush;
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
      cout<<"Computing col "<<j<<" of " << ((::divides(p,modulus)) ? W_opname : T_opname)
          <<"("<<p<<")..."<<flush;
    }
  svec ans = s_heckeop_col(p,j,0); // Automatically chooses W or T
  if(v)
    {
      cout<<"done."<<endl;
    }
  return ans;
}

smat homspace::s_opmat_cols(int i, const vec& jlist, int v)
{
  if(i==-1) return s_conj_cols(jlist,v);
  int d = dim(jlist);
  if((i<0)||(i>=nap))
    {
      cerr<<"Error in homspace::opmat_col(): called with i = " << i << endl;
      return smat(d,rk);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v)
    {
      cout<<"Computing "<<d<<" cols of " << ((::divides(p,modulus)) ? W_opname : T_opname)
          <<"("<<p<<")..."<<flush;
    }
  smat ans = s_heckeop_cols(p,jlist,0); // Automatically chooses W or T
  if(v)
    {
      cout<<"done."<<endl;
    }
  return ans;
}

smat homspace::s_opmat(int i, int dual, int v)
{
  if(i==-1) return s_conj(dual,v);
  if((i<0)||(i>=nap)) 
    {
      cerr<<"Error in homspace::s_opmat(): called with i = " << i << endl;
      return smat(dimension);  // shouldn't happen
    }
  long p = op_prime(i);
  if(v) 
    {
      cout<<"Computing " << ((::divides(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")..."<<flush;
      smat ans = s_heckeop(p,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop(p,dual,0); // Automatically chooses W or T
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
      cout<<"Computing " << ((::divides(p,modulus)) ? W_opname : T_opname) <<"("<<p
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
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
      cout<<"Computing " << ((::divides(p,modulus)) ? W_opname : T_opname) <<"("<<p
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
      smat ans = s_heckeop_restricted(p,s,dual,(v>2)); // Automatically chooses W or T
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
  if(::divides(p,modulus))  return vector<long>(pm1,pm1+2);
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
	    tvec *=2;
	  else
	    for (int i=1; i<=p2; i++) { add_coords(tvec,-i,p); }
	}
    }
  if(cuspidal)
    return cuspidalpart(tvec.as_vec());
  else
    return tvec.as_vec();
}

vec homspace::manintwist(long p) const
{
 svec sum = coords(0,p);                   // =0, but sets the right length.
 for (long i=1; i<p; i++) sum += legendre(i,p)*coords(i,p);
 if(cuspidal) return cuspidalpart(sum.as_vec()); 
 else return sum.as_vec();
}

#if(0) // no longer used
vec homspace::projmaninvector(long p) const  // Will only work after "proj"
{
  long i,p2;
  vec tvec = proj_coords(0,p);             // =0, but sets the right length.
  if (plusflag==-1) return tvec;
  if (p==2) 
    add_proj_coords(tvec,1,2);
  else
    { 
      p2=(p-1)>>1;
      for (i=1; i<=p2; i++) { add_proj_coords(tvec,i,p); }
      if(plusflag)   
	tvec *=2;
      else
	for (i=1; i<=p2; i++) { add_proj_coords(tvec,-i,p); }
    }
  return tvec;
}

vec homspace::projmaninvector(long p, const mat& m) const
{
  long i,p2;
  vec tvec = proj_coords(0,p,m);       // =0, but sets the right length.
  if (plusflag==-1) return tvec;
  if (p==2) 
    add_proj_coords(tvec,1,2,m);
  else
    { 
      p2=(p-1)>>1;
      for (i=1; i<=p2; i++) { add_proj_coords(tvec,i,p,m); }
      if(plusflag)   
	tvec *=2;
      else
	for (i=1; i<=p2; i++) { add_proj_coords(tvec,-i,p,m); }
    }
  return tvec;
}

vec homspace::newhecke(long p, long n, long d) const
                                     // Will only work after "proj"
{ vec tvec = proj_coords(p*n,d);
// cout<<"newhecke starts: tvec = "<<tvec<<endl; 
 long p2 = p>>1, dp = d*p, k;
  for (k=1+p2-p; k<=p2; k++) 
    {
      add_proj_coords(tvec,n+d*k,dp);
      // cout<<"tvec = "<<tvec<<endl; 
    }
  // cout<<"newhecke returns: tvec = "<<tvec<<endl; 
  return tvec;
}
#endif // no longer used

matop::matop(long p)
{

  //    case 31: nmats = 106; break;
  //    case 37: nmats = 128; break;
  //    case 41: nmats = 146; break;
  //    case 43: nmats = 154; break;
  //    case 47: nmats = 170; break;

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
