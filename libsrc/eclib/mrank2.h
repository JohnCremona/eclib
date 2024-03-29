// mrank2.h: declaration of class rank2 for descent via 2-isogeny
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

#if     !defined(_ECLIB_MRANK2_H)
#define _ECLIB_MRANK2_H      1       //flags that this file has been included

#include <eclib/descent.h>

class rank2 : public rank12 { // class to do 2-descent via 2-isogeny
private:
  vector<bigint> badprimes, supp0, supp1;  
  vector<bigint> elsgens0, elsgens1, els2gens0, els2gens1, glsgens0, glsgens1;  
  long selmer_rank_phi_Eprime, selmer_rank_phiprime_E, selmer_rank_Eprime;
  long best_rank_bound, best_isogeny, index2;
  long nt2gens0, nt2gens1, mask0, mask1;
  long els0, els1, gls0, gls1;   // after first descent
  long els20, els21, gls20, gls21;   // after second descent
  int d_is_sq, ddash_is_sq;
  bigint e2, e3, e2dash, e3dash, s2, s4, s6;
  Curvedata ee, eedash, Eprime;
  vector<Point> pointlist, fullpointlist, two_torsion;  
  int npoints, npoints1, fullnpoints, ntwo_torsion;
public: 
  rank2(Curvedata* ec, int verb,int sel=0, long l1=20, long l2=5,int second=1);
// lim1  is (bigint) bound on |x|+|z| in naive search
// lim2 is (double) bound on log max {|x|,|z| }, i.e. logarithmic
// sel is selmer_only switch
// second is do-second-descent switch
  Curvedata getEprime() const {return Eprime;}
  long getselmerprime() const {return selmer_rank_Eprime;}
  long getselmerphi() const {return selmer_rank_phi_Eprime;}
  long getselmerphiprime() const {return selmer_rank_phiprime_E;}
  void listgens();
  void listgens(Curvedata* CD_orig, const bigint& u, const bigint& r, 
		const bigint& s, const bigint& t);
  void listpoints();
  void listpoints(Curvedata* CD_orig, const bigint& u, const bigint& r, 
		  const bigint& s, const bigint& t);
  void makepoints();
  vector<Point> getpoints() 
  {
    if(fullnpoints==0) makepoints(); 
    return fullpointlist;
  }
  vector<Point> getgens() const {return pointlist;}
private:
  int testquartic(const bigint& c, const bigint& d1, const bigint& d2, int which);
  int second_descent(const bigint& c, const bigint& d1, const bigint& d2, int which);
  void makepoint(const bigint& c, const bigint& d1, const bigint& d2, 
		 const bigint& x, const bigint& y, const bigint& z, 
		 int which);
  void local_descent(const bigint& x0);
  void find_elsgens(int which, const bigint& c, const bigint& d);
  void find_els2gens(int which, const bigint& c, const bigint& d);
  void find_glsgens(int which, const bigint& c, const bigint& d);
  void makegens();
};

#endif
