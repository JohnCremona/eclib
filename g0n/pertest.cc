// FILE PERTEST.CC: Program to find curve from input periods (needs upgrade)
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
#include "marith.h"
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "newforms.h"
#include "h1newforms.h"
#include "periods.h"
#include "curve.h"     //from qcurves, for computing conductors

//#define AUTOLOOP

int main(void)
{
 cout.precision(15);
 int limit,n=1; 
 double x,y; int type;
 Complex w1, w2, c4, c6;
//  cout << "Enter type: "; cin >> type;
//  if(type==2) cout << "Lattice [x,yi]\n";
//  else if(type==1) cout << "Lattice [2x,x+yi]\n";
//  else {cout << "Type must be 1 or 2, assuming 2" << endl; type=2;}
//  cout << "Enter x: "; cin >> x;
//  cout << "Enter y: "; cin >> y;
//
//1290#7
// type = 2;
// x = 3.40948866706572;
// y = 37.755995750309/12;
//1290#8
// type = 2;
// x = 8.86405604424407/106;
// y = 55.0347994600134/424;
//1518#14
// type = 2;
// x = 1.81141045097419;
// y = 9.26930968357432/20;
//1560#4
// type = 1;
// x = 8.63724186942595/11;
// y = 7.38662908312853/11;
// int scale;  cout << "Enter factor to scale up by: "; cin >> scale;
// y*=scale;
// 1590#18
// type=1;
// x = 0.087985643044491;
// y = 2.20157842288697/3;
// 1599#3
// type=1;
 x = 0.0214941423671314;
 y = 6.12762332794969/3.0;
 
 cout << "x = " << x << "; y = " << y << endl << endl;
 cout << "Enter  max. scaling factors for x and y: ";
 int maxnx, maxny;
 cin >> maxnx >> maxny;
 int nx, ny;
 for(type=1; type<=2; type++) 
   {
     cout << "Trying type = " << type << endl;
     for(ny=1; ny<=maxny; ny++)
       {
	 for(nx=1; nx<=maxnx; nx++)
	   {
	     if(type==2){w1=Complex(x/nx,0); w2=Complex(0,y/ny);}
	     else {w1=Complex(2*x/nx,0); w2=Complex(x/nx,y/ny);}
	     Complex tau=normalize(w1,w2);
	     getc4c6(w1,w2,c4,c6);
	     Integer ic4 = Iround(real(c4));
	     Integer ic6 = Iround(real(c6));
	     if(valid_invariants(ic4,ic6)) 
	       {
		 cout << "\n nx = " << nx << ", ny = " << ny << "\n"; 
		 cout << "w1 = " << w1 << "; w2 = " << w2 << endl;
		 cout << "c4 = " << c4 << "\nc6 = " << c6 << endl;
		 cout << "ic4 = " << ic4 << "\nic6 = " << ic6 << endl;
		 Curve C(ic4,ic6);
		 Curvedata CD(C,1);
		 cout << (Curve)CD << "\t";
		 CurveRed CR(CD);
		 Integer cond = getconductor(CR);
		 cout << "N = " << cond << endl;
	       }
	   }
       }
   }
}       // end of main()
