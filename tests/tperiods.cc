// tperiods.cc -- input a curve, find its periods
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
//
//#define TEST

#include <eclib/compproc.h>
#include <eclib/cperiods.h>
#include <eclib/getcurve.h>

int main(){
  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);
	
  int verb=1;
  bigint v;
  vector<bigrational> ai(5);

  while (getcurve(ai,verb))
    {
      cout << "Input curve = ";
      cout <<"["<<ai[0]<<","<<ai[1]<<","<<ai[2]<<","<<ai[3]<<","<<ai[4]<<"]:" << endl;
      Curvedata CDin(ai,v); // v holds scaling factor
      Curvedata CD(CDin,1); // 1 for minimise
      cout << "Minimal model = "<<(Curve)CD<<endl;
      Cperiods cp(CD);
      cout << "Periods: " << cp << endl; 

      Curve EE = cp.trans_to_curve();
      cout << "Curve from periods: " << EE << endl;
    }
}
