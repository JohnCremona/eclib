// FILE TESTG1.CC: test program for G1 function (needs upgrade)
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
//
#include <builtin.h>
#include <iostream.h>

#include <values.h>
#include "moddata.h"
#include "symb.h"
#include "oldforms.h"
#include "homspace.h"
#include "newforms.h"
#include "cperiods.h"     //from qcurves, for computing conductors
#include "h1newforms.h"
#include "periods.h"

int main(void)
{
  cout.precision(15);
  int which;
  cout << "Which function (0,1,2,3): "; cin>>which;
  double xmin, xmax, xstep;
  cout<<"Enter min and max x and step size: "; cin >> xmin >> xmax >> xstep;
  if(xmin<=0) xmin+=xstep;
  for(double x = xmin; x<xmax; x+=xstep)
    {
      switch(which) {
      default:
      case 0:      cout << "x = " << x << ", g0(x) = " << myg0(x) << endl; 
	break;
      case 1:      cout << "x = " << x << ", g1(x) = " << myg1(x) << endl; 
	break;
      case 2:      cout << "x = " << x << ", g2(x) = " << myg2(x) << endl; 
	break;
      case 3:      cout << "x = " << x << ", g3(x) = " << myg3(x) << endl; 
	break;
      }
    }
}
