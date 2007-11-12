// FILE CHECKNF.CC: compute #newforms and compare with file -- not upgraded
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
#include <fstream>
#include "xsplit.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "newforms.h"

#define AUTOLOOP

int main(void)
{
 long n=110, limit=210, startp, stopp; 
 int output, verbose=0;
  char* name = new char[20];

 cout << "Program checknf.  Using METHOD = " << METHOD << " to find newforms" << endl;
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 // cout << "Verbose output? "; cin>>verbose;
#ifdef AUTOLOOP
     cout<<"Enter first and last N: ";cin>>n>>limit; n--;
     while (n<limit) { n++;
#else
     while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
  cout << ">>>Level " << n << "\t";
  newforms nf(n,0,10,verbose); // (level,use_old,depth for splitting,verbose)
  cout << nf.n1ds << " newform(s) found, \t";
  short num=0, nap=0;
  sprintf(name,"eigs/x%d\0",n);
  ifstream in(name);
  if(!in.is_open())
    {
      cout<<"Eigs file "<<name<<" does not exist!"<<endl;
    }
  else
    {
      in.read((char*)&num,sizeof(short));
      in.close();
    }
  cout << num << " newform(s) on file; ";
  if(num!=nf.n1ds) cout<<" numbers disagree!!!";
  else cout<<" OK";
  cout<< endl;
}       // end of if(n)
}       // end of while()
abort();
}       // end of main()
