// bsdprocs.h : class for computing L^(r)(f,1) for a newform
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

#include <values.h>

const double twopi = 2.0*PI;
const double eps = 1.0e-20;   // ?? mindouble;
const double eulergamma = 0.577215664901532860606512;
const double zeta3 = 1.20205690315959428540;

double myg0(double x);
double myg1(double x);
double myg2(double x);
double myg3(double x);


class ldash1 {
private:
  int N,limit,nap,r;
  double rootmod, factor, sum;
  longlist aplist;  
  longlist primelist;
  double (*g)(double x);

  void use(int n, int an);
  void add(int n, int pindex, int y, int z);

public:
  ldash1 (newform* f); 
  ~ldash1 () {;}
  double value() const {return 2*sum;}
  int rank() const {return r;}
};

