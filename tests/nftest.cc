// FILE NFTEST.CC: test program for newform constructor
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
//
#include <eclib/interface.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/oldforms.h>
#include <eclib/cperiods.h>
#include <eclib/newforms.h>

//#define AUTOLOOP
#define LMFDB_ORDER       // if defined, sorts newforms into LMFDB order before output

int main(void)
{
  int limit,firstn,n=2,count=0; 
  int verbose=1;
#ifdef AUTOLOOP
  cerr<<"Enter first and last N: ";cin>>firstn>>limit; 
  n=firstn-1; cout<<endl;
  while (n<limit) { n++;
#else
  verbose=1;
  while (n>1) { cerr<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
 cout << ">>> Level " << n << " <<<\t";
 newforms nf(n,verbose);
 cout << "\nAfter constructor, about to createfromdata() \n";
 nf.createfromdata(1,25);
 int num = nf.n1ds; count+=num;
 cout << num << " newform(s), "<<nf.nap<<" eigs on file." << endl;
#ifdef LMFDB_ORDER
 nf.sort_into_LMFDB_label_order();
#else
 nf.sort_into_Cremona_label_order();
#endif
 nf.display();
}       // end of if(n)
}       // end of while()
#ifdef AUTOLOOP
cout << "\n" << count << " newform(s) in range " << firstn << "..." << limit << endl;
#endif
}       // end of main()
