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
#include "getcurve.h"

int main(){
  set_precision(string("Enter number of decimal places").c_str());
  initprimes(string("PRIMES").c_str(),0);
	
  int verb=1;
  bigint v;
  vector<bigrational> ai(5);

  while (getcurve(ai,verb))
    {
      cout << "Input curve = ";
      cout <<"["<<ai[0]<<","<<ai[1]<<","<<ai[2]<<","<<ai[3]<<","<<ai[4]<<"]:" << endl;
      Curvedata CD(Curvedata(ai,v),1);
      cout << "Minimal model = "<<(Curve)CD<<endl;
      Cperiods cp(CD);
      cout << "Periods: " << cp << endl; 

      Curve EE = cp.trans_to_curve();
      cout << "Curve from periods: " << EE << endl;
    }
}
