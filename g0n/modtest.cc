// MODTEST.CC  -- Test for modular symbols
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

#include "timer.h"
#include "moddata.h"
#include "symb.h"

#define AUTOLOOP

int main(void)
{
 init_time();
 start_time();
 int n=1; 
 int plus=1;
 int verbose=0;
 cout << "Display symbol details (0/1)? " << flush; cin >> verbose;
 int limit; 
#ifdef AUTOLOOP
 cout<<"Enter first and last levels: ";cin>>n>>limit; n--;
 while (n<limit) { n++;
#else
     while (n>0) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>0)
{
 cout << ">>>Level " << n << "\t";
 symbdata symbols(n);
 cout<<"("<<symbols.nsymb<<" symbols)\t";
 if(verbose) symbols.display();
 symbols.check();
}       // end of if(n)
}       // end of while()
stop_time();
show_time();
 cout<<endl; 
}       // end of main()
