// FILE HECKETEST.CC  -- Test program for Hecke operators
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
//

#include <eclib/timer.h>
#include <eclib/homspace.h>

//#define AUTOLOOP
//#define COMPARE_OLD
//#define CHECK_COMMUTE
//#define TEST_EIGS

vector<long> eigrange(long p)
{
  long aplim=3, four_p=p<<2;
  while (aplim*aplim<=four_p) aplim++;
  aplim--;
  vector<long> ans(1+2*aplim);
  iota(ans.begin(),ans.end(),-aplim);
  return ans;
}


int main(void)
{
 cout << "Program hecketest." << endl;
 init_time();
 start_time();
 scalar modulus = default_modulus<scalar>();
 int n=2;
 int plus=1;
 int verbose=0;
 cerr << "See the hecke matrices (0/1)? "; cin >> verbose;
 cerr << "Plus space (0/1)? "; cin >> plus;
#ifdef AUTOLOOP
     int limit; 
     cerr<<"Enter limit on level: ";cin>>limit;
     while (n<limit) { n++;
#else
     while (n>1) { cerr<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
 cout << ">>>Level " << n << "\t";
 homspace hplus(n, modulus, plus,0,0);
 int genus = hplus.h1dim();
 scalar den = hplus.h1denom();
 ZZ den2; den2 = den*den;
 cout << "Dimension = " << genus << "\n";
 cout << "denominator = " << den << "\n";
 vector<long> badprimes = hplus.plist;
 int nq = badprimes.size(); int firstq=0;  // =0 for all W's
 if (genus>0)
   {
   mat_m id = mat_m::identity_matrix(genus);
   mat_m id2 = den2*id;
   vector<mat_m> wqlist(nq);
   cout << "Computing conjmat...  " << flush;
   smat conjmat = hplus.s_conj(1);
   cout<<" done."<<endl;
   cout << "Computing +1 eigenspace...  " << flush;
   ssubspace h1plus = eigenspace(conjmat,den, modulus);
   cout<<" done, dimension = "<<dim(h1plus)<<endl;
   cout << "Computing -1 eigenspace...  " << flush;
   ssubspace h1minus = eigenspace(conjmat,-den, modulus);
   cout<<" done, dimension = "<<dim(h1minus)<<endl;

   int w_eigs=0;
   cout<<"Compute W-eigenspaces? "; cin>>w_eigs;
   for (int i=0; i<nq; i++)
     {long q=badprimes[i]; if(i<firstq) continue;
      cout << "Computing W("<<q<<")...  " << flush;
      mat wq = hplus.heckeop(q,verbose);
      cout << "done, sparsity = "<<sparsity(wq)<<". " << endl;
      if(verbose)
	{
	  cout<<"Computed matrix "; 
	  if(den>1) cout<<" (scaled by "<<den<<") ";
	  cout<<" = \n"<<wq<<endl<<endl;
	  //	  wq.output_pretty();
	}
      if(w_eigs) {
      smat swq(wq);
      int e;
      for(e=1; e>-2; e-=2)
	{
          long mult;
          /*
	  cout<<"Using modular matrix code..."<<flush;
	  start_time();
	  mult=dim(peigenspace(wq,e*den, modulus));
	  stop_time();
	  show_time(cerr); cerr<<endl;
	  cout<<"Dimension of "<<e<<"-eigenspace="<<mult<<endl;
	  */
	  cout<<"Using sparse matrix code..."<<endl;
	  start_time();
	  mult=dim(eigenspace(swq,e*den, modulus));
	  stop_time();
	  show_time(cerr); cerr<<endl;
	  cout<<"Dimension of "<<e<<"-eigenspace="<<mult<<endl;
	}
      }
      wqlist[i] = to_mat_m(wq);
#ifdef CHECK_COMMUTE
      if (mult_mod_p(swq,swq,modulus) == smat::scalar_matrix(genus, den*den))
	cout << "Involution!" << "\n";
      else
	cout << " *** Not an involution...." << "\n";
#else
      cout<<endl;
#endif
    }

   int np=5,ip=0; 
   cout<<"How many T_p? "; cin>>np;
   vector<mat_m> tplist(np);
   for (primevar pr(np+nq); pr.ok()&&ip<np; pr++, ip++)
     {while (n%pr==0) pr++;
      int p=pr;
      cout << "\nComputing T_p for p = " << p << "..." << flush;
#ifdef COMPARE_OLD
      cout<<endl;
      start_time();
      mat_m temp = mat_m(hplus.heckeop(p,verbose));
      stop_time();
      cout<<"Time for old method: "; show_time(cerr); cerr<<endl;
#endif // COMPARE_OLD
      start_time();
      mat_m tp = to_mat_m(hplus.newheckeop(p,verbose));
      if(verbose)
	{
	  cout<<"Computed matrix "; 
	  if(den>1) cout<<" (scaled by "<<den<<") ";
	  cout<<" = \n"<<tp<<endl<<endl;
	  //	  tp.output_pretty();
	}
      tplist[ip] = tp;
      //      cout<<"Copied  mat_m = "<<tplist[ip]<<endl;
      stop_time();
#ifdef COMPARE_OLD
      cout<<"Time for new method: "; show_time(cerr); cerr<<endl;
      if(temp!=tplist[ip]) cout<<"Matrices differ!\n";
#else
      cout << "done, sparsity = "<<sparsity(tplist[ip])<<". " << endl;
#endif // COMPARE_OLD
#ifdef TEST_EIGS
      vector<long> eigs = eigrange(p); // hplus.eigrange(nq+ip);
      cout<<"\nChecking for eigenvalues from "<<eigs<<endl;
      long i,j,k;
      long nulty, nulty1, totalmult=0;
      scalar dummy;
      mat m = tplist[ip].shorten(dummy);
      //      if(verbose) cout<<"shortened matrix: \n"<<m<<endl;
      smat sm=smat(m);

      for(k=0; k<eigs.size(); k++)
	{
	  long e = eigs[k];
	  cout<<"\nTrying eigenvalue e = "<<e<<" ("<<e<<")"<<endl;
	  e *= den;

	  /*
	  cout<<"Computing nullity, using my (modular) matrix code..."<<flush;
	  start_time();
	  nulty=dim(peigenspace(m,e,modulus));
	  stop_time();
	  show_time(cerr); cerr<<endl;
	  cout<<" nullity="<<nulty<<endl;
	  */
	  cout<<"Computing nullity, using my sparse matrix code..."<<flush;
	  start_time();
	  nulty=dim(eigenspace(sm,e));
	  stop_time();
	  show_time(cerr); cerr<<endl;
	  cout<<" nullity="<<nulty<<endl;

	  cout<<"\n"<<e<<" (scaled) ";
	  if(nulty>0) cout<<" IS "; else cout<<" is NOT ";
	  cout<<"an eigenvalue"; 
	  if(nulty>0) cout<<", of multiplicity "<<nulty;
	  cout<<endl;

	}

      mat_ZZ M;
      M.SetDims(genus,genus);
      for(i=1; i<=genus; i++) 
	for(j=1; j<=genus; j++) 
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
      cout<<"...done, sparsity =  "<<sparsity(MT); show_time(cerr); cerr<<endl;
      cout<<"Computing kernel, using my (modular) matrix code..."<<flush;
      start_time();
      subspace  ker = pkernel(MT,modulus);
      nulty=dim(ker);
      int denker=denom(ker); // =1 for modular method, but set below
      stop_time();
      cout<<"done, nulty = "<<nulty; show_time(cerr); cerr<<endl;
      cout<<endl;
      vector<long> eigs1;
      long n1f = 0;

      if(nulty>0)
	{
	  cout<<"lifting kernel..."<<flush;
	  start_time();
	  mat MTR;
          int ok = liftmat(prestrict(m,ker,modulus),modulus,MTR,denker);
          if (!ok)
            cout << "**!!!** failed to lift modular kernel\n" << endl;
	  stop_time();
	  cout<<"done, denom(ker)="<<denker; show_time(cerr); cout<<endl;
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
	  cout<<"done ";  show_time(cerr); cerr<<endl;
	  vec_pair_ZZX_long factors;
	  cout<<"\nfactorizing char poly..."<<flush;
	  start_time();
	  factor(cont,factors,cptp);
	  stop_time();
	  cout<<"done ";  show_time(cerr); cerr<<endl;
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
	  show_time(cerr); cerr<<endl;
	  cout<<"det="<<detM<<endl;
	  */

	  cout<<"Computing determinant, using deterministic strategy..."<<flush;
	  start_time();
	  detM = determinant(Msub,1);
	  stop_time();
	  show_time(cerr); cerr<<endl;
	  cout<<"det="<<detM<<endl;

	  cout<<"Computing nullity, using my (modular) matrix code..."<<flush;
	  start_time();
	  nulty1=dim(peigenspace(MTR,e,modulus));
	  stop_time();
	  show_time(cerr); cerr<<endl;
	  cout<<" nullity="<<nulty1<<endl;

	  cout<<"\n"<<e<<" (scaled) ";
	  if(nulty>0) cout<<" IS "; else cout<<" is NOT ";
	  cout<<"an eigenvalue"; 
	  if(nulty>0) cout<<", of multiplicity "<<nulty1;
	  cout<<endl;

	  totalmult+=nulty1;
	}
	}
      cout<<"Total multiplicity of rational eigenvalues = "<<totalmult<<endl;
#endif // TEST_EIGS
#ifdef CHECK_COMMUTE
      ZZ P(modulus);
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
#endif //CHECK_COMMUTE
     } // loop on p
   }      // end of if(genus>0)
 }       // end of if(n)
     }       // end of while(n>1) or while(n<limit)
 cerr<<endl;
exit(0);
     }       // end of main()
