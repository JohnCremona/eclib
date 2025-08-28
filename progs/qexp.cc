// FILE QEXP.CC: program for listing coefficients of q-expansions
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

#include <eclib/newforms.h>
#include <eclib/curvesort.h> // for codeletter

#define AUTOLOOP
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output

#define NAP 25      // number of ap to output
#define SEPCHAR "," // char to separate ap in output

const scalar modulus(default_modulus<scalar>());

int main(void)
{
  cerr<<"q-expansions of rational newforms";
#ifdef LMFDB_ORDER
  cerr<<" in LMFDB order (simple lexicographic)";
#else
  cerr<<" in Cremona order";
#endif
  cerr<<endl;
  int limit,firstn,n=1;
  int verbose=0;
  unsigned int nap;
#ifdef AUTOLOOP
  cerr<<"Enter first and last N: ";cin>>firstn>>limit; cerr<<endl;
  n=firstn-1;
  while (n<limit) { n++;
#else
  while (n>1) { cerr<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
 newforms nf(n, modulus, verbose);
 nf.createfromdata(1,25);
 int i, num = nf.n1ds;
 if(num>0){
 nap = nf.nflist[0].aplist.size();
 if (nap>NAP) nap=NAP;
 if(verbose)
   {
     cout << ">>> Level " << n << " <<<\t";
     cout << num << " newform(s) "<<endl;
   }
#ifdef LMFDB_ORDER
  nf.sort_into_LMFDB_label_order();
#else
  nf.sort_into_Cremona_label_order();
#endif
 for (i=0; i<num; i++)
   {
     cout<<n<<codeletter(i)<<": ";
     vector<long>v = nf.nflist[i].aplist;
     copy(v.begin(),v.begin()+nap, ostream_iterator<long>(cout, SEPCHAR));
     //vec_out(cout,v,20);
     cout<<"..."<<endl;
   }
 }
}       // end of if(n)
}       // end of while()
}       // end of main()
