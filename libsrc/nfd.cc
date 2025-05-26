// FILE nfd.cc: implementation of class nfd (higher-dimensional newforms)
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

#include <eclib/marith.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/ZZXFactoring.h>
#include <eclib/nfd.h>

#define OUTPUT_PARI_STYLE

nfd::nfd(homspace* in_h1, int one_p, int w_split, int mult_one, int verbose)
{
  h1=in_h1;
  long n = h1->modulus;
  long dimh = h1->h1dim();
  long denh = h1->h1denom(); dH=denh;
  vector<long> badprimes = h1->plist;
  mat K = basis(h1->kern).as_mat();
  mat_m tp, tp1, m;
  long d, i,j, p;
  bigint ap1;

  Hscales.resize(dimh+1);
  Hscales[0]=1;
  for(i=1; i<=dimh; i++) Hscales[i]=Hscales[i-1]*denh;

// Compute the desired linear combination of Tp:

  if(one_p) // Compute one Tp:
    {
      primevar pr;
      while (n%pr==0) pr++;
      p=pr;
      cout << "Computing T_p for p = " << p << "..." << flush;
      tp = to_mat_m(transpose(h1->newheckeop(p,0)));
      cout<<"done."<<endl;
    }
  else
    {
      tp.init(dimh,dimh); // zero matrix
      while(1)
	{
	  cout<<"Enter a (good) prime (or 1 for id, 0 to stop): "; cin>>p;
	  if(p==0) break;
	  if(p==1)
	    {
	      cout<<"Coefficient of identity: "; cin>>ap1;
	      if(ap1!=0) tp = addscalar(tp,ap1);
	    }
	  else
	    {
	      cout << "Computing T_p for p = " << p << "..." << flush;
	      tp1 = to_mat_m(transpose(h1->newheckeop(p,0)));
	      cout<<"done."<<endl;
	      cout<<"coefficient of T_"<<p<<": "; cin>>ap1;
	      if(ap1!=1) tp1*=ap1;
	      if(ap1!=0) tp+=tp1;
	    }
	}
    }

// Compute the appropriate W-eigenspace and restrict to it

  subspace_m SW(dimh);
  int dimsw=dimh;
  if(w_split)
    {
      int nq = badprimes.size();
      for(i=0; (i<nq)&&(dimsw>0); i++)
	{
	  long q = badprimes[i];
	  bigint eq;
	  cout<<"Enter eigenvalue of W("<<q<<"): ";
	  cin>>eq;
	  eq *=dH;
	  mat_m wq = to_mat_m(transpose(h1->heckeop(q,0)));
	  if(dimsw<dimh)
	    {
	      SW=subeigenspace(wq,eq,SW);
	    }
	  else // still at top level
	    {
	      SW=eigenspace(wq,eq);
	    }
	  dimsw=dim(SW);
	  cout<<"eigenspace now has dimension "<<dimsw<<endl;
	}
      if(dimsw<dimh) tp = restrict_mat(tp,SW);
    }
  if(dimsw==0)
    {
      cout<<"This W-eigenspace is trivial!"<<endl;
      return;
    }

  mat_ZZ ntl_tp;
  ntl_tp.SetDims(dimsw,dimsw);
  for(i=1; i<=dimsw; i++)
    for(j=1; j<=dimsw; j++)
      ntl_tp(i,j)=tp(i,j);

  bigint swden=denom(SW);
  Sscales.resize(dimsw+1);
  Sscales[0]=1;
  for(i=1; i<=dimsw; i++) Sscales[i]=Sscales[i-1]*swden;

// Compute char poly of restriction of tp to this subspace:

  ZZX ntl_cptp; ZZ cont;
  CharPoly(ntl_cptp, ntl_tp);
  vec_pair_ZZX_long factors;
  //  SetCoeff(ntl_cptp,dimsw,1);
  for(i=0; i<dimsw; i++)
    {
      bigint temp = coeff(ntl_cptp,i);
      divide_exact(temp,Hscales[dimsw-i]*Sscales[dimsw-i],temp);
      SetCoeff(ntl_cptp,i,temp);
    }
  cout<<"char poly = "<<ntl_cptp<<endl;

  if(mult_one)
    {

// factor the charpoly:

      SquareFreeDecomp(factors,ntl_cptp);
      if(verbose)  cout<<"NTL char poly square-free factors = "<<factors<<endl;

      if(factors[0].b>1)
	{
	  cout<<"No factors of multiplicity 1"<<endl;
	  return;
	}
      else
	{
	  cout<<"Factors of multiplicity 1 are:"<<endl;
	}
      factor(cont,factors,factors[0].a);
    }
  else
    {
      factor(cont,factors,ntl_cptp);
      cout<<"Factors are:"<<endl;
    }

  long nf = factors.length();
  for(i=0; i<nf; i++)
    {
      cout<<(i+1)<<":\t"<<factors[i].a
	  <<"\t(degree "<<deg(factors[i].a)<<")";
      if(!mult_one) cout<<"\t to power "<<factors[i].b;
      cout<<endl;
    }

// select subspace:

  int looking=1;
  vector<bigint> coeffs;
  while(looking)
    {
      cout<<"Enter factor number: "; cin>>j;
      if((j<1)||(j>nf)) 
	{cout<<"Must be between 1 and "<<nf<<endl; continue;}
      j--;
      if(factors[j].b!=1) {cout<<"Multiplicity>1!\n"; continue;}
      d = deg(factors[j].a);
      cout<<"Degree = "<<d<<endl;
      m = tp;
      minpol.resize(d);
      coeffs.resize(d);
      for(i=d-1; i>=0; i--)
	{
	  minpol[i]=coeff(factors[j].a,i);
	  coeffs[i]=coeff(factors[j].a,i)*Hscales[d-i]*Sscales[d-i];
	  m = addscalar(m,coeffs[i]);
	  if(i) m = m*tp;
	}
      if(verbose) 
	{
	  cout<<"(unscaled) min poly = [1 ";
	  for(i=d-1; i>=0; i--) cout<<coeffs[i]<<" ";
	  cout<<"]"<<endl;
	}
      cout<<"(rescaled) min poly = "<<factors[j].a<<endl;
      S = kernel(m);
      if(dim(S)!=d) 
	{
	  cout<<"Problem: eigenspace has wrong dimension ("<<dim(S)<<")"
	      <<endl;
	}
      else looking=0;
    }
  
  //    if(verbose) 
      cout<<"finished constructing S, now restricting T_p to S"<<endl;

      tp0 = restrict_mat(tp,S);

  //  if(verbose) 
      cout<<"done.  now combining S and SW"<<endl;

  if(w_split)// make S a subspace of H_1, not of the W-eigenspace
    {
      mat_m SWbasis=basis(SW);
      bigint  SWden; SWden=denom(SW);
      subspace_m mSW(SWbasis,pivots(SW),SWden);
      S=combine(mSW,S);  
    }

  long dims=dim(S);
  dS=denom(S);
  long sden=I2long(dS);
  long sden2=sden*denh;
  dHS=dH*dS;

  //  if(verbose) 
    {
      if(sden2>1) cout<<sden2<<"*";
      cout<<"Matrix of T("<<p<<") restricted to S is ";
      showmatrix(tp0); cout<<endl;
    }
  // = matrix of T_p on irreducible subspace of dual space

    //  if(verbose) 
    {
      cout<<"The former poly is the min poly of alpha_1 = "
	  <<sden2<<"*alpha"<<endl;
    }
  cout<<"The latter is the min poly of alpha, ";
  cout<<"which is the eigenvalue of T("<<p<<")"<<endl;

  if(verbose) cout<<"Finished computing (dual) subspace S"<<endl;	 
  if(verbose>1||(sden2>1))
    {
      cout<<"S has denom "<<sden<<", cumulative denom = "<<sden2<<endl;
    }
  V = transpose(basis(S)); // so V is dims x dimh
  Sscales.resize(dims+1);
  Sscales[0]=1;
  for(i=1; i<=dims; i++) Sscales[i]=Sscales[i-1]*sden;
  
  mat_m A=transpose(tp0);
  W.init(dims,dims); Winv.init(dims,dims);
  vec_m v(dims);  v[1]=1; // so v=[1,0,...,0]
  W.setcol(1,v);
  for(i=2; i<=dims; i++) {v = A*v; W.setcol(i,v);}
  Wdetnum = inverse(W,Winv);
  WinvV = Winv*V;
  if(verbose)
    {
      cout<<"W     = ";showmatrix(W); cout<<endl;
      cout<<"W^(-1)= (1/"<<Wdetnum<<")*";showmatrix(Winv);cout<<endl;
      if(verbose>1)
	{
	  cout<<"WinvV = ";showmatrix(WinvV);cout<<endl;
	}
    }

// compute projcoord, precomputed projections of the modular symbol basis

  long ncoord = h1->coord_vecs.size()-1;
  projcoord.init(ncoord,dims);
  coord_fac=0;
  vec_m mrowi(dims), coordi(dimh);
  vec rowi(dims);
  for (i=1; i<=ncoord; i++)
    {
      coordi = to_vec_m((h1->coord_vecs[i]).as_vec());
      if(h1->cuspidal) coordi = to_vec_m(h1->cuspidalpart(to_vec_i(coordi)));
      mrowi = V*coordi;
      rowi = to_vec_i(mrowi);
      projcoord.setrow(i,rowi);
      coord_fac=gcd(coord_fac,(long)content(rowi));
    }
  if(verbose>1) cout<<"content of projccord = "<<coord_fac<<endl;
  if(coord_fac>1)  projcoord /= coord_fac;

  Wdetdenom = coord_fac;
  Wdetnum*=dHS;

  Winv_scaled=Winv;
  bigint g; g=content(Winv_scaled.row(1));
  for(i=2; i<=dims; i++)
    {
      Winv_scaled.multrow(i,Hscales[i-1]*Sscales[i-1]);
      g=gcd(g,content(Winv_scaled.row(i)));
    }
  // now g is the content of Winv_scaled
  Winv_scaled/=g;
  Wdetdenom*=g;
  g = gcd(Wdetnum,Wdetdenom);
  if(g>1) {Wdetnum/=g; Wdetdenom/=g;}
  cout<<"Basis for Hecke eigenvalues, in terms of powers of alpha:"<<endl;
  for(i=1; i<=dims; i++)
    {
      cout<<"("<<Wdetdenom<<"/"<<Wdetnum<<")*";
      cout<<Winv_scaled.col(i)<<endl;
    }
}

// ap_vec has length dim(S); last entries hold numerator and
// denominator of content
vec_m nfd::ap(long p)
{
  mat K = basis(h1->kern).as_mat();
  long rk = K.nrows();
  long k,l,n = h1->modulus, dims=dim(S);
  vec_m apvec(dims);
  int bad = ::divides(p,n);
  if(bad) return apvec; // temporary fix!
  matop matlist(p);
  // if(bad)
  //   matlist=matop(p,n);

  for(k=0; k<rk; k++)
    {
      long Kkj = K(k+1,pivots(S)[1]);
      if(Kkj!=0)
	{
	  bigint mKkj; mKkj = Kkj;
	  // if(bad)
	  //   {
          //     continue;  // not yet implemented
	  //   }
	  // else
	    {
	      symb s = h1->symbol(h1->freegens[k]);
	      for(l=0; l<matlist.size(); l++)
                apvec += mKkj*to_vec_m(matlist[l](s,h1,projcoord));
	    }
	}
    }
  return apvec;
}

mat_m nfd::oldheckeop(long p)
{
  return restrict_mat(transpose(to_mat_m(h1->newheckeop(p,0))),S);
}

mat_m nfd::heckeop(long p)
{
  mat K = basis(h1->kern).as_mat();
  long rk = K.nrows();
  long j,k,l,n = h1->modulus, dimh=h1->h1dim(), dims=dim(S);
  int bad = ::divides(p,n);
  matop matlist(p);
  if(bad)
    {
      cout<<"q = "<<p<<"\t";
      matlist=matop(p,n);
    }
  else
    {
      cout<<"p = "<<p<<"\t";
    }
  mat_m TE(dimh,dims);
  vec_m colj(dimh);
  for (j=0; j<dims; j++)
    { 
      colj.init(dimh);
      for(k=0; k<rk; k++)
	{
	  long Kkj = K(k+1,pivots(S)[j+1]);
	  if(Kkj!=0)
	    {
	      bigint mKkj; mKkj = Kkj;
	      if(bad)
		{
		  vec vt = (h1->applyop(matlist,h1->freemods[k])).as_vec();
		  if(h1->cuspidal) vt=h1->cuspidalpart(vt);
		  colj += (mKkj*to_vec_m(vt));
		}
	      else
		{
		  symb s = h1->symbol(h1->freegens[k]);
		  for(l=0; l<matlist.size(); l++)
		    {
		      vec vt = (matlist[l](s,h1)).as_vec();
		      if(h1->cuspidal) vt=h1->cuspidalpart(vt);
		      colj += mKkj*to_vec_m(vt);
		    }
		}
	    }
	}
      TE.setcol(j+1,colj);
    }
  return transpose(V*TE);
}

bigint inverse(const mat_m& a, mat_m& ainv)
{
  long d = a.nrows();
  mat_m aug=colcat(a, mat_m::identity_matrix(d));
  long rk, ny; vec pc,npc; bigint denom;
  mat_m ref = echelon(aug, pc, npc, rk, ny, denom);
  ainv = ref.slice(1,d,d+1,2*d);
  //  cout<<"Inverse = "<<denom<<"*"<<ainv<<endl;
  return denom;
}

void showmatrix(const mat_m& m)
{
#ifdef OUTPUT_PARI_STYLE
  long i,j, nc=m.ncols(),nr=m.nrows();
  cout << "[";
  for(i=0; i<nr; i++)
    {
      if(i) cout<<";";
      for(j=0; j<nc; j++) 
	{ 
	  if(j) cout<<","; 
	  cout<<m(i+1,j+1);
	}
    }
  cout << "]\n";
#else
  cout<<m;
#endif
}

void showmatrix(const mat& m)
{
#ifdef OUTPUT_PARI_STYLE
  m.output_pari(cout);
#else
  cout<<m;
#endif
}

