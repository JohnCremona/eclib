// descent.h: declaration of classes rank12 and two_descent
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
#include "mwprocs.h"

// rank12 is a common base for separate classes rank1 and rank2 (for
// computing rank via general 2-descent and descent via 2-isogeny
// repectively); class two_descent is a user interface to these

class rank1;
class rank2;

class rank12 {
protected: 
  Curvedata* the_curve;
  int verbose, certain, success, selmer_only, do_second_descent;
  long num_aux;
  long rank, selmer_rank;
  Curvedata IJ_curve;  // [0,0,0,-27*I,-27*J]
  bigint tr_u,tr_r,tr_s,tr_t;  // transformation from latter to minimal curve
  long lim1, lim2;
public:

// Constructor:
//
// sel is selmer_only switch
// firstlim is bound on |x|+|z|
// secondlim is bound on log max {|x|,|z| }, i.e. logarithmic
// n_aux only relevant for general 2-descent when 2-torsion trivial
// n_aux=-1 causes default to be used (depends on method)
// second_descent only relevant for descent via 2-isogeny

  rank12(){;}
  rank12(Curvedata* ec, 
         int verb=0, int sel=0, 
         long firstlim=20, long secondlim=5, 
	 long n_aux=-1, int second_descent=1)
    :the_curve(ec), verbose(verb), selmer_only(sel), 
     do_second_descent(second_descent),
     num_aux(n_aux), lim1(firstlim), lim2(secondlim) {;}
  virtual ~rank12() {;}
  long getrank() const {return rank;}
  long getselmer() const {return selmer_rank;}
  long getcertain()  const {return certain;}
  int ok()      const {return success;}
  virtual void listpoints()=0;
  virtual void listpoints(Curvedata* CD_orig, 
			  const bigint& u, const bigint& r, 
			  const bigint& s, const bigint& t)=0;
  virtual vector<Point> getgens() const =0;
  virtual vector<Point> getpoints() =0;
};


class two_descent {
private: rank12 * r12;  // does all the work
  int verbose, two_torsion_exists, selmer_only;
  int success, certain, fullmw;
  long rank, selmer_rank, sat_bound;
  mw mwbasis;
public:
// Constructor:
//
// sel is selmer_only switch
// firstlim is bound on |x|+|z|
// secondlim is bound on log max {|x|,|z| }, i.e. logarithmic
// n_aux only relevant for general 2-descent when 2-torsion trivial
// n_aux=-1 causes default to be used (depends on method)
// second_descent only relevant for descent via 2-isogeny

  two_descent(Curvedata* ec, 
	      int verb=0, int sel=0, 
	      long firstlim=20, long secondlim=5, 
	      long n_aux=-1, int second_descent=1);
  ~two_descent() {delete r12;}
  long getrank() const {return rank;}
  long getselmer() const {return selmer_rank;}
  long getcertain()  const {return certain;}
  int ok()      const {return success;}
  int get2t() const {return two_torsion_exists;}
  int getfullmw() const {return fullmw;}
  bigfloat regulator() {return mwbasis.regulator();}
  vector<Point> getbasis() {return mwbasis.getbasis();}
  vector<Point> getbasis(Curvedata* CD_orig, 
			  const bigint& u, const bigint& r, 
			  const bigint& s, const bigint& t);
  void report_rank() const;
  void saturate(long sat_bd); // =0 for none
  void show_gens(int change,
		 Curvedata* CD_orig, 
		 const bigint& u, const bigint& r, 
		 const bigint& s, const bigint& t);
  void show_result_status();
  void pari_output(int change,
		   Curvedata* CD_orig, 
		   const bigint& u, const bigint& r, 
		   const bigint& s, const bigint& t);
};

