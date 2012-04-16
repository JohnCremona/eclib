// FILE QEXP.CC: program for listing coefficients of q-expansions
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
//
#include <eclib/interface.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/oldforms.h>
#include <eclib/cperiods.h>
#include <eclib/newforms.h>
#include "curvesort.cc"

#define AUTOLOOP
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output

#define NAP 25      // number of ap to output
#define SEPCHAR "," // char to separate ap in output 

int main(void)
{
  cerr<<"q-expansions of rational newforms";
#ifdef LMFDB_ORDER
  cerr<<" in LMFDB order (simple lexicographic)";
#endif
  cerr<<endl;
  int limit,firstn,n=1; 
  int verbose=0;
  unsigned int nap;
  char* code = new char[20];
#ifdef AUTOLOOP
  cerr<<"Enter first and last N: ";cin>>firstn>>limit; cerr<<endl;
  n=firstn-1;
  while (n<limit) { n++;
#else
  while (n>0) { cerr<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
 newforms nf(n,verbose);
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
 nf.sort();
 // cout<<"After sort():"<<endl;
 for (i=0; i<num; i++)
   {
     codeletter(i,code);
     cout<<n<<code<<": ";
     vector<long>v = nf.nflist[i].aplist;
     copy(v.begin(),v.begin()+nap, ostream_iterator<long>(cout, SEPCHAR));
     //vec_out(cout,v,20);
     cout<<"..."<<endl;
   }
 }
}       // end of if(n)
}       // end of while()
}       // end of main()
