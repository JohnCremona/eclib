// conductor.cc: program to call Tate's algorithm and display conductors
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
 
#include "curve.h"
#include "getcurve.h"

int main(void)
{
  initprimes(string("PRIMES").c_str(),0);
        
  int verb=1;
  bigint v;
  vector<bigrational> ai(5);

  while (getcurve(ai,verb))
    {
      cout <<"["<<ai[0]<<","<<ai[1]<<","<<ai[2]<<","<<ai[3]<<","<<ai[4]<<"]:\t\t" << flush;
      Curvedata CD(Curvedata(ai,v),1);
      CurveRed CR(CD);
      cout << "N = " << getconductor(CR) << endl;
    }
}

