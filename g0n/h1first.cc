// FILE h1first.cc :  h1 (full space)
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

#include "marith.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "cperiods.h"     //from qcurves, for computing conductors
#include "newforms.h"
#include "periods.h"
#include "pcprocs.h"

#define AUTOLOOP
#define SHOWCURVES

int main(void)
{
  set_precision("Enter number of decimal places");
  int verbose,limit=210,n=130; 
  cout << "Program h1first.  Using METHOD = " << METHOD << endl;
  cout << "Verbose output? "; cin>>verbose;
#ifdef AUTOLOOP
  cout<<"Enter first and last N: ";cin>>n>>limit;  n--;
  while (n<limit) { n++;
#else
  while (n>0) { cout<<"\n\nEnter level: "; cin>>n;
#endif
  if (n>0)
    {
      int usedata=0;
      if(verbose) cout << "\n\n";
      cout << ">>>Level " << n << "<<<" << endl;
      newforms nf(n,verbose); 
      nf.createfromdata(0,25,0);
      if(verbose) 
	cout << "finished reading newform data from file" << endl;     
      nf.makebases(1);
      if(verbose) 
	cout << "finished constructing h1newforms" << endl;     
      nf.output_to_file();
      if(verbose) 
	cout << "saved updated newform data to file" << endl;     
      
#ifdef SHOWCURVES
      for(int i=0; i<nf.n1ds; i++)
	{
	  if(verbose) cout << "\nForm number ";
	  cout << i+1 << ": ";
	  bigfloat rperiod;
	  Curve C = nf.getcurve(i, -1, rperiod, verbose);
	  Curvedata CD(C,1);  // The 1 causes minimalization
	  if(verbose) cout << "\nCurve = \t";
	  cout << (Curve)CD << "\t";
	  CurveRed CR(CD);
	  cout << "N = " << getconductor(CR) << endl;
	  if(verbose) cout<<endl;
	}
#endif
      
    }       // end of if(n)
  }       // end of while()
  cout<<"Finished"<<endl;
  }       // end of main()

