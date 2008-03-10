// FILE NFTEST.CC: test program for newform constructor
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
//
#include "interface.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "cperiods.h"
#include "newforms.h"

//#define AUTOLOOP

int main(void)
{
  int limit,firstn,n=1,count=0; 
  int verbose=0;
#ifdef AUTOLOOP
  cerr<<"Enter first and last N: ";cin>>firstn>>limit; 
  n=firstn-1; cout<<endl;
  while (n<limit) { n++;
#else
  verbose=1;
  while (n>0) { cerr<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
 cout << ">>> Level " << n << " <<<\t";
 int plusspace=1;
 int cuspidal=0;
 newforms nf(n,plusspace,cuspidal,verbose);
 cout << "\nAfter constructor, about to createfromdata() \n";
 nf.createfromdata(25);
 int num = nf.n1ds; count+=num;
 cout << num << " newform(s), "<<nf.nap<<" eigs on file." << endl;
 // nf.makebases();
 /*
 nf.sort();
 cout<<"After sort():"<<endl;
 nf.display();
 nf.sort(1);
 cout<<"After sort(1):"<<endl;
 nf.display();
 */
}       // end of if(n)
}       // end of while()
#ifdef AUTOLOOP
cout << "\n" << count << " newform(s) in range " << firstn << "..." << limit << endl;
#endif
}       // end of main()
