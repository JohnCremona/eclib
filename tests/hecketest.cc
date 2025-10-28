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
#include <eclib/polys.h>

//#define AUTOLOOP

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
 ZZ Modulus = to_ZZ(modulus);
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
 ZZ Den = to_ZZ(den);
 cout << "Dimension = " << genus << "\n";
 cout << "denominator = " << den << "\n";
 vector<long> badprimes = hplus.plist;
 int nq = badprimes.size();
 if (genus==0)
   continue; // to next level

 cout << "Computing conjmat...  " << flush;
 smat conjmat = hplus.s_conj(1);
 cout<<" done."<<endl;
 cout << "Computing +1 eigenspace...  " << flush;
 ssubspace h1plus = eigenspace(conjmat,den, modulus);
 cout<<" done, dimension = "<<dim(h1plus)<<endl;
 cout << "Computing -1 eigenspace...  " << flush;
 ssubspace h1minus = eigenspace(conjmat,-den, modulus);
 cout<<" done, dimension = "<<dim(h1minus)<<endl;

 vector<mat_ZZ> Wqlist;
 int w_eigs=0;
 cout<<"Compute W-eigenspaces? "; cin>>w_eigs;
 for ( auto q: badprimes)
   {
     cout << "Computing W("<<q<<")...  " << flush;
     mat wq = hplus.heckeop(q,verbose);
     mat_ZZ Wq = mat_to_mat_ZZ(wq);
     cout << "done, sparsity = "<<sparsity(wq)<<". " << endl;
     if(verbose)
       {
         cout<<"Computed matrix ";
         if(den>1) cout<<" (scaled by "<<den<<") = ";
         output_flat_matrix(wq);
         cout<<endl;
       }
     smat swq(wq);
     if(w_eigs)
       {
         for( auto e: {1, -1})
           {
             long mult;
             cout<<"Using sparse matrix code..."<<endl;
             start_time();
             mult=dim(eigenspace(swq,e*den, modulus));
             stop_time();
             show_time(cerr); cerr<<endl;
             cout<<"Dimension of "<<e<<"-eigenspace="<<mult<<endl;
           }
       }
     if (check_involution(Wq, Den, Modulus))
       cout << "Involution!" << "\n";
     else
       cout << " *** Not an involution...." << "\n";
     if (check_commute(Wq, Wqlist, Modulus))
       cout << "Commutes with all previous Wq!" << endl;
     else
       cout << "*** Does NOT commute with all previous Wq." << endl;
     cout << "Factored char poly:"<<endl;
     display_factors(scaled_charpoly(Wq, Den));
     Wqlist.push_back(Wq);
   }

 int np=5,ip=0;
 cout<<"How many T_p? "; cin>>np;
 vector<mat_ZZ> Tplist;
 for (primevar pr(np+nq); pr.ok()&&ip<np; pr++, ip++)
   {
     while (n%pr==0) pr++;
     int p=pr;
     cout << "\nComputing T_p for p = " << p << "..." << flush;
     start_time();
     mat tp = hplus.newheckeop(p,verbose);
     mat_ZZ Tp = mat_to_mat_ZZ(tp);
     if(verbose)
       {
         cout<<"Computed matrix ";
         if(den>1) cout<<" (scaled by "<<den<<") = ";
         output_flat_matrix(tp);
         cout<<endl;
       }
     stop_time();
     cout << "done, sparsity = "<<sparsity(tp)<<". " << endl;
     if (check_commute(Tp, Wqlist, Modulus))
       cout << "Commutes with all Wq!" << endl;
     else
       cout << "*** Does NOT commute with all Wq." << endl;
     if (check_commute(Tp, Tplist, Modulus))
       cout << "Commutes with all previous Tp!" << endl;
     else
       cout << "*** Does NOT commute with all Tp." << endl;
     cout << "Factored char poly:"<<endl;
     display_factors(scaled_charpoly(Tp, Den));
     Tplist.push_back(Tp);
   } // end of primes loop
 }       // end of if(n)
}       // end of while(n>1) or while(n<limit)
  cerr<<endl;
  exit(0);
}       // end of main()
