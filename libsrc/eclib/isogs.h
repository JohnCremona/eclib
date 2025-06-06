// isogs.h:  declaration of class IsogenyClass and related functions
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
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
 
#ifndef _ECLIB_ISOGS_H_
#define _ECLIB_ISOGS_H_

#include "cperiods.h"  // which itself includes everything relevant
#include "matrix.h"

// isogeny functions:

vector<CurveRed> twoisog(const CurveRed& CR, int verbose);
vector<CurveRed> threeisog(const CurveRed& CR, int verbose);

vector<CurveRed> lisog(const CurveRed& CR, Cperiods& cp,
			  long ell, int verbose=0);
			  // cp must be normalized for lattice at call time
			  // returned array is of length up to 3
			  // of curves ell-isogenous to CR


int semistable(CurveRed& CR, const bigint& p);
int semistable(CurveRed& CR);

vector<long> getelllist(CurveRed& CR); 
// returns list of possible primes l ("ell") for which CR might have an l-isogeny.

#define MAXNCURVES 26    // max number of curves in an isogeny class.

class IsogenyClass {
private: 
  vector<CurveRed> curves;
  vector<long> llist; 
  long nell, ncurves, ndone;
  int ss, verb;
  Cperiods cp;
  vector<long> fromlist; // fromlist[i]=j if curve i first constructed from 
  vector<long> isoglist; // curve j via an isoglist[i]-isogeny
  vector<long> matij;
  void matset(long i, long j, long ell) { matij[i*MAXNCURVES+j]=ell;}
  void process(long i);  // process i'th curve
public:
  IsogenyClass(CurveRed& C, int verbose);
  void grow(void);       // does the work
  void display_llist(ostream& os)const {os<<llist;}
  void displaycurves(ostream& os)const;
  void displaymat(ostream& os)const;
  void display(ostream& os)const {displaycurves(cout); displaymat(cout);}
  void dumpdata(ostream& os, long rank);  // output for textab to input
  vector<long> getmat() const;
  mat_i getmatrix() const;
  long mat_entry(long i, long j)const {return matij[i*MAXNCURVES+j];}
  vector<CurveRed> getcurves() const {return curves;}
};


#endif
