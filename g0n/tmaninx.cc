// FILE TMANINX.CC: Program for computing more ap (needs upgrading)
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
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "newforms.h"
#include "manin.h"

//#define AUTOLOOP

int main(void)
{
 int n=1; 
 int lastn, startp, stopp, output, showeigs, showforms;
 cout << "Program tmaninx\n";
 cout << "---------------\n\n";
 cout << "For each N, assumes that the file eigs/eN exists, and computes more\n";
 cout << "Hecke eigenvalues.\n";
 cout << "Output new eigs to file (1/0)? ";  cin>>output;
 cout << "Output new eigs to screen (1/0)? "; cin>>showeigs;
 cout << "Display newforms (1/0)? "; cin>>showforms;
#ifdef AUTOLOOP
 cout << "Which primes for Hecke eigenvalues (first#, last#)? ";
 cin >> startp >> stopp; cout << endl;
 cout<<"Enter first and last N: ";cin>>n>>lastn; n--;
 while (n<lastn) { n++;
#else
 while (cout<<"Enter level: ", cin>>n, n>0) { 
 cout << "Which primes for Hecke eigenvalues (first#, last#)? ";
 cin >> startp >> stopp; cout << endl;
#endif
  cout << "\n>>>Level " << n << "\t";
  manin machine(n,1,10,0,0);   
// (level, use_old, depth, cuspidal, verbose)
  if (showforms) machine.newforms::display();
  char* eigfilename = eigfile(n);
  machine.getap(startp,stopp,output,eigfilename,showeigs);
  delete[] eigfilename;
}       // end of while()
}       // end of main()
