// MODTEST.CC  -- Test for modular symbols
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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

#include <eclib/interface.h>
#include <eclib/timer.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>

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
     while (n>1) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
 cout << ">>>Level " << n << "\t";
 symbdata symbols(n);
 cout<<"("<<symbols.nsymb<<" symbols)\t";
 if(verbose) symbols.display();
 symbols.check();
}       // end of if(n)
}       // end of while()
stop_time();
show_time(cerr); cerr<<endl; 
}       // end of main()
