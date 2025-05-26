// mrank1.h -- declaration of class rank1 for general 2-descent
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

#if     !defined(_ECLIB_MRANK1_H)
#define _ECLIB_MRANK1_H      1       //flags that this file has been included

#include <eclib/descent.h>
#include <eclib/mquartic.h>

class rank1 : public rank12 {
private:
  long nquarticsa, nquarticsb, nfirstlota, nfirstlotb;
  long sha_rank, sha2;
  int traceequiv, posdisc, disc_is_square, npairs, extra2, threediv, type;
  vector<quartic> qlista, qlistb;
  vector<int> qlistbflag;
  vector<bigcomplex> croots, cphi;
  vector<Point> pointlist1, pointlist2;
  long npoints1, npoints2;
  int have_eggpoint, have_large_quartics;
  long twoadic_index, global_index;
  long bsd_npairs; // only for testing
             // 1, 2 or 4: local/global index of "small" quartics
  bigint c4, c6, d1728, ii, jj, disc;
  long Imod2, Jmod2;
  bigfloat xii, xjj;
  vector<bigint> plist, dlist;
  vector<long> eqplist;  // primes used for equiv-sieving
  long n0, n1, n2, rank_B;
//
//
// Sieving stuff:
  int ipivot, pivflag;
  vector<long> auxs, amod, astepmod, hmod, hstepmod, hscalemod;
  vector<vector<long>> phimod;
  vector<int> aux_flags, aux_types;
  vector<vector<int>> squares, flaga;
  vector<vector<vector<int>>> flags;
  long ah_count, ah_sieve_0, ah_sieve_1, ah_sieve_2;
  long ah_rfail, ah_dfail, ah_efail, ah_extra2fail, ah_pass;

  void aux_init();  // define  auxiliary moduli and squares
  void flag_init(); // set up flag array
  vector<long> qeps(const quartic&, int x2); //computes eps of quartic
  void show_eps_vec(const vector<long>& vec);

  // process latest quartic found
  void addquartic(const bigint& a, const bigint& b, const bigint& c,
		  const bigint& d, const bigint& e);
  void getquartics();
  void getquartics1();
  void gettype(int t);
public:
  explicit rank1(Curvedata* ec,
                 int verb=0, int sel=0,
                 long lim1=20, long lim2=5, long n_aux=-1);
// lim1 is bound on |x|+|z| in naive search
// lim2 is bound on log max {|x|,|z| }, i.e. logarithmic
// sel is selmer_only switch
// n_aux is # sieving primes in quartic search
// n_aux=-1 causes default to be used (depends on method)
//
  void sortpoints();
  void listpoints();
  void listpoints(Curvedata* CD_orig, const bigint& u, const bigint& r,
		                      const bigint& s, const bigint& t);
  vector<Point> getgens() const;
  vector<Point> getpoints();
  long getselmerprime() const {return selmer_rank;}
  long getselmerphi() const {return selmer_rank;}
  long getselmerphiprime() const {return selmer_rank;}
  Curvedata getEprime() const {return *the_curve;}
};

#endif
