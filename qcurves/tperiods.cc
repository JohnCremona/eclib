// tperiods.cc -- input a curve, find its periods
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
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
//#define TEST

#include "compproc.h"
#include "cperiods.h"
#include "reader.h"

int main(){
  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);
	
  CurveReader input;
  Curve E;
  Curvedata D;
  CurveRed C;

  while (input>>E)
    {
      cout << "Curve input:        " << E << endl;
      Cperiods cp(E);
      cout << "periods: " << cp << endl; 

      Curve EE = cp.trans_to_curve();
      cout << "Curve from periods: " << EE << endl;
    }
}
