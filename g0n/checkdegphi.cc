// FILE CHECKDEGPHI.CC: not upgraded
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

#include <time.h>
#include <fstream.h>
#include <builtin.h>
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "splitter.h"
#include "oldforms.h"
#include "h1newforms.h"

#define AUTOLOOP

int main(void)
{
  int n=110; 
  int limit=210, startp, stopp, output;
  int disp=0; // cout << "Display newform info?"; cin>>disp;

#ifdef AUTOLOOP
  cout<<"Enter first and last N: ";cin>>n>>limit; n--;
  while (n<limit) { n++;
#else
  while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
  if (n>0)
    {
      output=1; 
//      cout << "\n>>>Level " << n << "\n";
      h1newforms m(n); // from eigs
      h1newforms nf(n, /* display */ 0, /* usefiledata */1);
      if(disp) m.display();
      for(int i=0; i<m.n1ds; i++)
	{
	  long degphinew = m.nflist[i].degphi;
	  long degphiold = nf.nflist[i].degphi;
	  int agree = (degphiold==degphinew);
	  if(!agree)
	    {
	      cout<<"N = "<<n<<"\tForm # "<<(i+1);
	      cout << "\tOld deg(phi) = " << degphiold;
	      cout << "\tNew deg(phi) = " << degphinew << endl;
	    }
	}
    }       // end of if(n)
}       // end of while()
}       // end of main()
