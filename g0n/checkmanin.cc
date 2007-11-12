// FILE CHECKMANIN.CC: not upgraded
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

#include "xsplit.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "newforms.h"
#include "manin.h"

#define AUTOLOOP
#define VERB 0

int main(void)
{
 long firstn, lastn, n, startp, stopp; 
 long nbad=0;
 char* name = new char[20];
 int output=0, verbose=VERB;
 cout << "Program checkmanin.  Using METHOD = " << METHOD << " to find newforms" << endl;
#ifdef MODULAR
 cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
 cout << "Verbose output? "; cin>>verbose;
#ifdef AUTOLOOP
     cout<<"Enter first and last N: ";cin>>firstn>>lastn; n=firstn-1;
     while (n<lastn) { n++;
#else
     while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
  output=1; 
  if(verbose)cout<< "\n>>>Level " << n<<endl;

  int filenum, prognum;
  sprintf(name,"eigs/x%d\0",n);
  ifstream in(name);
  if(!in.is_open())
    {
      cout<<"Unable to open file "<<name<<" for newform input"<<endl;
    }
  else
    {
      short temp;
      in.read((char*)&temp,sizeof(short)); filenum=temp;
      in.close();
    }
  if(verbose) cout << filenum << " newform(s) on file." << endl;
  newforms nf(n,0,10,0,0); // (level,use_old,depth,cuspidal,verbose)
  prognum=nf.n1ds;
  if(verbose) cout << prognum << " newform(s) found."<<endl;

  if(prognum==filenum)
    {
      if(verbose) cout << "Numbers AGREE\n" << endl;
    }
  else
    {
      nbad++;
      if(verbose) cout << "Numbers DISAGREE\n" << endl;
      else
	{
	  cout << ">>>Level " << n << " -- ";
	  cout << "Numbers DISAGREE\n" << endl;
	  cout << filenum << " newform(s) on file." << endl;
	  cout << prognum << " newform(s) found."<<endl;
	}
    }
  if(n%10==0) cout << nbad << " disagreements found up to " <<n<<endl;
}       // end of if(n)
}       // end of while()
cout << nbad << " disagreements found in the range "
     <<firstn<<".."<<lastn<<endl;
abort();
}       // end of main()
