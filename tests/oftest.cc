// FILE OFTEST.CC  -- Test program for oldform class
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
#include "moddata.h"
#include "oldforms.h"

#define AUTOLOOP

int main(void)
{
 int n=11, limit; 
 int verbose; 
 cout << "Verbose details of oldform constructor? ";  
 cin >> verbose;

#ifdef AUTOLOOP
     cout<<"Enter first and last N: ";cin>>n>>limit; n--;
     while (n<limit) 
       { n++;
#else
     while (cout<<"Enter level: ", cin>>n, n>0) 
       {
#endif
	 cout << ">>>Level " << n << "\t";
	 moddata symbols(n);   // (not really needed except 
	                       //that the level data gets initialized properly 
                               //which IS needed for oldforms)
	 oldforms of(10,&symbols,verbose);       // default args ntp=5, verbose=0
	 of.display();
       }       // end of while()
}       // end of main()
