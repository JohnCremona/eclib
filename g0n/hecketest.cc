// FILE HECKETEST.CC  -- Test program for Hecke operators
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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
//
#include "interface.h"
#include "timer.h"
#ifdef LiDIA_INTS
#include <LiDIA/bigint_matrix.h>
#else
#ifdef NTL_INTS
#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/LLL.h>
#endif
#endif
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "smatrix_elim.h"
//#include "mmatrix.h"
//#include "msubspace.h"

//#define AUTOLOOP
//#define COMPARE_OLD
//#define CHECK_COMMUTE
//#define TEST_EIGS

double sparsity(const mat_m& m);
double sparsity(const mat& m);

int main(void)
{
 cout << "Program hecketest.  Using METHOD = " << METHOD << " to find newforms" << endl;
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 init_time();
 start_time();
 bigint P = BIGINT(BIGPRIME);
 int n=1; 
 int plus=1;
 int verbose=0;
 cout << "See the hecke matrices (0/1)? "; cin >> verbose;
 cout << "Plus space (0/1)? "; cin >> plus;
 int limit; 
#ifdef AUTOLOOP
     cout<<"Enter limit on level: ";cin>>limit;
     while (n<limit) { n++;
#else
     while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
 cout << ">>>Level " << n << "\t";
 homspace hplus(n,plus,0,0);
 int genus = hplus.h1dim();
 long den = hplus.h1denom();
 bigint den2; den2 = den*den;
 cout << "Dimension = " << genus << "\n";
 cout << "denominator = " << den << "\n";
 vector<long> badprimes = hplus.plist;
 int nq = badprimes.size(); int firstq=0;  // =0 for all W's
 if (genus>0)
   {
   mat_m id = idmat(genus);
   mat_m id2 = den2*id;
   mat_m* wqlist = new mat_m[nq];
   cout << "Computing conjmat...  " << flush;
   smat conjmat = hplus.s_conj(1);
   cout<<" done."<<endl;
   cout << "Computing +1 eigenspace...  " << flush;
   ssubspace h1plus = eigenspace(conjmat,den);
   cout<<" done, dimension = "<<dim(h1plus)<<endl;
   cout << "Computing -1 eigenspace...  " << flush;
   ssubspace h1minus = eigenspace(conjmat,-den);
   cout<<" done, dimension = "<<dim(h1minus)<<endl;

#if(1)
   for (int i=0; i<nq; i++)
     {long q=badprimes[i]; if(i<firstq) continue;
      cout << "Computing W("<<q<<")...  " << flush;
      mat wq = hplus.heckeop(q,verbose);
      cout << "done, sparsity = "<<sparsity(wq)<<". " << endl;
      if(verbose)
	{
	  cout<<"Computed matrix "; 
	  if(den>1) cout<<" (scaled by "<<den<<") ";
	  cout<<" = "<<wq<<endl;
	  //	  wq.output_pretty();
	}
      smat swq(wq);
      int e; long mult;
      for(e=1; e>-2; e-=2)
	{
	  /*
	  cout<<"Using modular matrix code..."<<flush;
	  start_time();
	  mult=dim(peigenspace(wq,e*den,MODULUS));
	  stop_time();
	  show_time(); 
	  cout<<"\nDimension of "<<e<<"-eigenspace="<<mult<<endl;
	  */
	  cout<<"Using sparse matrix code..."<<flush;
	  start_time();
	  mult=dim(eigenspace(swq,e*den));
	  stop_time();
	  show_time(); 
	  cout<<"\nDimension of "<<e<<"-eigenspace="<<mult<<endl;
	}
      wqlist[i]=wq;
#ifdef CHECK_COMMUTE
      if (mult_mod_p(swq,swq,BIGPRIME)==den*den*sidmat(genus)) 
	cout << "Involution!" << "\n";
      else
	cout << "NOT an involution...." << "\n";
#else
      cout<<endl;
#endif
    }

   int np=5,ip=0; 
   cout<<"How many T_p? "; cin>>np;
   mat_m* tplist = new mat_m[np];
   for (primevar pr(np+nq); pr.ok()&&ip<np; pr++, ip++)
     {while (n%pr==0) pr++;
      int p=pr;
      cout << "\nComputing T_p for p = " << p << "..." << flush;
#ifdef COMPARE_OLD
      cout<<endl;
      start_time();
      mat_m temp = hplus.heckeop(p,verbose);
      stop_time();
      cout<<"Time for old method: "; show_time();
#endif
      start_time();
      mat tp = hplus.newheckeop(p,verbose);
      if(verbose)
	{
	  cout<<"Computed matrix "; 
	  if(den>1) cout<<" (scaled by "<<den<<") ";
	  cout<<" = "<<tp<<endl;
	  //	  tp.output_pretty();
	}
      tplist[ip] = tp;
      //      cout<<"Copied  mat_m = "<<tplist[ip]<<endl;
      stop_time();
#ifdef COMPARE_OLD
      cout<<"Time for new method: "; show_time();
      if(temp!=tplist[ip]) cout<<"Matrices differ!\n";
#else
      cout << "done, sparsity = "<<sparsity(tplist[ip])<<". " << endl;
#endif
#ifdef TEST_EIGS
      vector<long> eigs=hplus.eigrange(nq+ip);
      cout<<"\nChecking for eigenvalues from "<<eigs<<endl;
      long i,j,k,n=genus,r;
      long nulty, nulty1, totalmult=0;
      SCALAR dummy;
      mat m = tplist[ip].shorten(dummy);
      //      if(verbose) cout<<"shortened matrix: \n"<<m<<endl;
      smat sm=smat(m);
#ifdef LiDIA_INTS
      bigint_matrix M(n,n);
      for(i=1; i<=n; i++) 
	for(j=1; j<=n; j++) 
	  M.sto(i-1,j-1,tplist[ip](i,j));
      //      cout<<"LiDIA matrix = "<<M<<endl;

      for(k=0; k<(eigs.size()); k++)
	{
	  long e = eigs[k];
	  cout<<"\nTrying eigenvalue e = "<<e<<endl;
	  for(i=0;i<n;i++) M.sto(i,i,tplist[ip](i+1,i+1)-e*den);
	  //	  cout<<"Adjusted LiDIA matrix = "<<M<<endl;
	  //	  nulty = (n-rank(M));
/*
	  start_time();
	  bigint detM = det(M);
	  nulty = (detM==0);
	  stop_time();
	  cout<<"det(M-e) = "<<detM<<endl;
	  if(nulty) cout<<" IS "; else cout<<" is NOT ";
	  cout<<"an eigenvalue "; show_time();
*/
	  start_time();
	  long rankM = rank(M);
	  stop_time();
	  cout << "Nullity(M-e) = "<<n-rankM; show_time();
	  totalmult+=(rankM<n);
	}
#else
#ifdef NTL_INTS
      for(k=0; k<eigs.size(); k++)
	{
	  long e = eigs[k];
	  cout<<"\nTrying eigenvalue e = "<<e<<" ("<<e<<")"<<endl;
	  e *= den;


	  cout<<"Computing nullity, using my (modular) matrix code..."<<flush;
	  start_time();
	  nulty=dim(peigenspace(m,e,MODULUS));
	  stop_time();
	  show_time(); 
	  cout<<"\n nullity="<<nulty<<endl;

	  cout<<"Computing nullity, using my sparse matrix code..."<<flush;
	  start_time();
	  nulty=dim(eigenspace(sm,e));
	  stop_time();
	  show_time(); 
	  cout<<"\n nullity="<<nulty<<endl;

	  cout<<"\n"<<e<<" (scaled) ";
	  if(nulty>0) cout<<" IS "; else cout<<" is NOT ";
	  cout<<"an eigenvalue"; 
	  if(nulty>0) cout<<", of multiplicity "<<nulty;
	  cout<<endl;

	}

      mat_ZZ M;
      M.SetDims(n,n);
      for(i=1; i<=n; i++) 
	for(j=1; j<=n; j++) 
	  M(i,j)=m(i,j);
      //      cout<<"NTL matrix = "<<M<<endl;

      //      Evaluate "rational" charpoly of T_p:
      cout<<"Computing product of (T-a*I) over possible eigs a..."<<flush;
      start_time();
      mat mm=m;
      mat MT=mm, m2=m*m;
      for(k=1; k<(eigs.size()); k++)
	{
	  long e = eigs[k]*den;
	  if(e>0) 
	    MT = MT * addscalar(m2,(-e*e));
	}
      stop_time();
      cout<<"...done, sparsity =  "<<sparsity(MT); show_time(); cout<<endl;
      cout<<"Computing kernel, using my (modular) matrix code..."<<flush;
      start_time();
      subspace  ker=pkernel(MT,MODULUS);
      nulty=dim(ker);
      int denker=denom(ker); // =1 for modular method, but set below
      stop_time();
      cout<<"done, nulty = "<<nulty; show_time();
      cout<<endl;
      vector<long> eigs1;
      long n1f = 0;

      if(nulty>0)
	{
	  cout<<"lifting kernel..."<<flush;
	  start_time();
	  mat MTR = liftmat(prestrict(m,ker,MODULUS),MODULUS,denker);
	  stop_time();
	  cout<<"done, denom(ker)="<<denker; show_time(); cout<<endl;
	  if(nulty<21) 
	    cout<<"Restriction of Tp to relevant subspace = \n" << MTR << endl;
	  
	  mat_ZZ Msub;
	  Msub.SetDims(nulty,nulty);
	  for(i=1; i<=nulty; i++) 
	    for(j=1; j<=nulty; j++) 
	      Msub(i,j)=MTR(i,j);
	  
	  ZZX cptp; ZZ cont;
	  cout<<"computing char poly..."<<flush;
	  start_time();
	  CharPoly(cptp, Msub);
	  stop_time();
	  cout<<"done ";  show_time(); cout<<endl;
	  vec_pair_ZZX_long factors;
	  cout<<"\nfactorizing char poly..."<<flush;
	  start_time();
	  factor(cont,factors,cptp);
	  stop_time();
	  cout<<"done ";  show_time(); cout<<endl;
	  cout<<"\nFactors are:"<<endl;
	  long nf = factors.length();
	  for(i=0; i<nf; i++)
	    {
	      cout<<(i+1)<<":\t"<<factors[i].a
		  <<"\t(degree "<<deg(factors[i].a)<<")"
		  <<"\t to power "<<factors[i].b;
	      cout<<endl;
	      if(deg(factors[i].a)==1) 
		{
		  long ap = -I2long(coeff(factors[i].a,0))/denker;
		  cout<<"Adding eigenvalue "<<ap<<endl;
		  eigs1.push_back(ap); n1f++;
		}
	    }

      cout<<"Rational eigenvalues (scaled by "<<den<<") are "<<eigs1<<endl;

      for(k=0; k<n1f; k++)
	{
	  long e = eigs1[k];
	  cout<<"\nTrying eigenvalue e = "<<e<<" ("<<e<<")"<<endl;
	  e *= denker;
	  ZZ ee; ee = e;
	  for(i=1;i<=nulty;i++) Msub(i,i)=MTR(i,i)-ee;
	  if(nulty<21) cout<<"Adjusted NTL matrix = "<<Msub<<endl;

	  ZZ detM; 

	  /*
	  start_time();
	  cout<<"Computing determinant, using randomized strategy..."<<flush;
	  detM = determinant(Msub);
	  stop_time();
	  show_time(); 
	  cout<<"\ndet="<<detM<<endl;
	  */

	  cout<<"Computing determinant, using deterministic strategy..."<<flush;
	  start_time();
	  detM = determinant(Msub,1);
	  stop_time();
	  show_time(); 
	  cout<<"\ndet="<<detM<<endl;

	  /*
	  cout<<"Computing rank, using NTL's LLL..."<<flush;
	  start_time();
	  mat_ZZ M2; M2=Msub;
	  r=image(detM,M2); // NB image() changes its 2nd arg!
	  stop_time();
	  show_time(); 
	  nulty1 = nulty-r;
	  cout<<"\n nullity="<<nulty1<<endl;
	  */

	  cout<<"Computing nullity, using my (modular) matrix code..."<<flush;
	  start_time();
	  nulty1=dim(peigenspace(MTR,e,MODULUS));
	  stop_time();
	  show_time(); 
	  cout<<"\n nullity="<<nulty1<<endl;

	  cout<<"\n"<<e<<" (scaled) ";
	  if(nulty>0) cout<<" IS "; else cout<<" is NOT ";
	  cout<<"an eigenvalue"; 
	  if(nulty>0) cout<<", of multiplicity "<<nulty1;
	  cout<<endl;

	  totalmult+=nulty1;
	}
	}
#endif
#endif
      cout<<"Total multiplicity of rational eigenvalues = "<<totalmult<<endl;
#endif // TEST_EIGS
#ifdef CHECK_COMMUTE
      for (int kp=firstq; kp<nq; kp++)
	{if (matmulmodp(wqlist[kp],tplist[ip],P)!=matmulmodp(tplist[ip],wqlist[kp],P))
	   {
	     cout << "Problem: T_p matrix "<<ip<<" and W_q matrix "<<kp<<" do not commute!" << "\n";
	   }
       }
      for (int jp=0; jp<ip; jp++)
	{if (matmulmodp(tplist[ip],tplist[jp],P)!=matmulmodp(tplist[jp],tplist[ip],P))
	   {
	     cout << "Problem: T_p matrices "<<ip<<" and "<<jp<<" do not commute!" << "\n";
	   }
	}
#endif
     }
   delete[] wqlist; delete[] tplist;
 }      // end of if(genus>0)

}       // end of if(n)
}       // end of while()
     //stop_time();
     //show_time();
 cout<<endl;
abort();
}       // end of main()
#endif

double sparsity(const mat_m& m)
  {
    double count=0;
    long i,j,nr=nrows(m), nc=ncols(m);
    for(i=0; i<nr; i++)
      for(j=0; j<nc; j++)
	if(!is_zero(m(i+1,j+1))) count=count+1;
    return count/(nr*nc);
  }
