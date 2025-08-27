// File NEWFORMS.H
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

#ifndef _ECLIB_NEWFORMS_H
#define _ECLIB_NEWFORMS_H      1
                           //flags that this file has been included

#include <eclib/cperiods.h>
#include <eclib/xsplit.h> // for splitter_base class
#include <eclib/curve.h>
#include <eclib/oldforms.h>
#include <eclib/homspace.h>

class newforms;
class jumps;

#define DEFAULT_SMALL_NAP 25

/* Data stored in a newform (and in data files newforms/x$N and smallnf/x$N):
   (Numbers refer to lines of data file)
Items 1-18 are "int" while the ap and aq are "short"
3.  sfe : sign of functional equation (=-product of aq)
4.  ap0 : p0'th Hecke eigenvalue, p0=smallest good prime
5.  np0 : np0=1+p0-ap0
6.  dp0 : dp0/np0=L/P=L(f,1)/2x
7.  lplus : prime =1 (mod 4) with L(f,lplus,1) nonzero
8.  mplus : L(f,lplus,1)*sqrt(l)=mplus*x
9.  lminus : prime =3 (mod 4) with L(f,lminus,1) nonzero
10. mminus : L(f,lminus,1)*sqrt(-l)=mminus*yi
11-14. a, b, c, d : entries of a matrix M=[a,b;N*c,d] in Gamma_0(N) s.t.
15. dotplus       : the integral of f over {0,M(0)} is
16. dotminus      : dotplus*x+dotminus*yi
17. type : type 1 if period lattice = [2x,x+yi], type 2 if [x,yi]
18. degphi : degree of modular parametrization
aq : list of Wq-eigenvalues at bad primes
ap : list of Tp- & Wq-eigenvalues at all primes
*/

class newform {
  friend class newforms;
public:
  newforms *nf;  // the "parent"
  int sign;   // 1/-1 for old-style newform, 0 for old-style h1newform
  vec bplus,bminus; // DUAL eigenvectors
  int type;            // 2 for rectangular, 1 for non-rectangular period lattice
  long index;  // splitting index, -1 if not known
  vector<long> aplist, aqlist;
  long ap0;     // Eigenvalue of first "good" p
  long sfe;     // sign of functional equation
  long cuspidalfactorplus, cuspidalfactorminus;  // pdot =cuspidalfactor*np0
  long pdot,np0,dp0;  // np0=1+p0-ap0, pdot = maninvector(p0).bplus,
                      //                    = cuspidalfactor*dp0

  rational loverp;  // L(f,1)/x where x = least real part of a period
                    // =np0/dp0
  long lplus, lminus;  // primes = +1, -1 mod 4
  long mplus, mminus;  // mplus*x=sqrt(lplus)*L(fxlplus,1)
                       // mminus*yi=sqrt(-lminus)*L(fxlminus,1)
  long a,b,c,d,dotplus,dotminus; // matrix for period integration
  // Either type=1, lattice=[2x,x+yi]
  // Or     type=2, lattice=[x,yi]
  // & integral over [a,b;Nc,d] is dotplus*x+dotminus*yi
  long degphi;             // degree of Weil parametrization
  long rk;                 // analytic rank (-1 if not set)
  bigfloat Lvalue;         // L^(r)(f,1)/r!
  vec coordsplus, coordsminus;  // vector components of each freegen
  long denomplus, denomminus, contplus, contminus;
  int j0; scalar fac;
  rational optimalityfactorplus, optimalityfactorminus;

  newform(void) {;}
  //  newform(const vec& v, const vector<long>& ap, newforms* nfs,long ind=-1);
  newform(const vector<int>& data, const vector<long>& aq, const vector<long>& ap, newforms* nfs);
  newform(const vec& vplus, const vec& vminus, const vector<long>& ap, newforms* nfs,long ind=-1);

  void add_more_ap(int nap);
  void display(void) const;
  // Testing function
  int check_expand_contract();

  // Explanation of the following three utilities: after newform
  // searching, each newform's aplist contains eigenvalues not Fourier
  // coefficients: these are the same for good primes p but for bad
  // primes q the eigenvalue is for thw AL-operator W_q.  Before
  // sorting and outputting this needs "fixing up" as in fixup_eigs().
  // After reading in from file (e.g. to compute more ap) we need to
  // go back to the eigenvalue list using unfix_eigs() before
  // recreating eigenspaces and then reverse this afterwards using
  // refix_eigs().  Each of these functions has a version in the
  // newforms class too, which applies the operation to every newform.

  // To fix eigenvalues lists after finding a newform: use when aplist
  // contains AL-eigenvalues w_q for bad primes q.  This extracts
  // those into the list aqlist and replaces them with the Fourier
  // coefficients a_q (=0 if q^2|N else -w_q).
  void fixup_eigs();
  // To fix eigenvalues lists before/after recovering bases: use when
  // aplist contains Fourier coefficients for bad primes q.  This
  // replaces those with AL-eigenvalues from aqlist.
  void unfix_eigs();
  // Same as fixup_eigs except that aqlist is not (re)created. It
  // replaces AL-eigenvalues in aplist with Fourier coefficients.
  void refix_eigs();

  // To find BSD ratio:
  void find_bsd_ratio();
  // To find projected coords:
  void find_coords_plus_minus();
  // To find cuspidal factors:
  void find_cuspidal_factors();
  // To find twisting primes:
  void find_twisting_primes();
  // To find optimality factors when created from a curve:
  void find_optimality_factors(const CurveRed& E, int i=0);
  // To find deg(phi):
  void find_degphi();
   // To get matrix and scale factors
  void find_matrix();
   // To normalize signs:
  void sign_normalize();
  // Compute analytic rank (if not already done)
  void compute_rank();
  // Return analytic rank (compute if not already done)
  long rank();
  // L^(r)(f,1)/r!
  bigfloat special_value();
};


class newforms :public level, splitter_base<scalar>   {
  friend class newform;
private:
  scalar modulus;
  int verbose; long maxdepth, cuspidal, sign;
  int basisflag;  // is set, then use() only sets bases for newforms
		  // already defined.
  mat opmat(int i, int d, int v=0)
  {return h1->opmat(i,d,v);}
  vec opmat_col(int i, int j, int v=0)
  {return h1->opmat_col(i,j,v);}
  mat opmat_cols(int i, const vec_i& jlist, int v=0)
  {return h1->opmat_cols(i,jlist,v);}
  mat opmat_restricted(int i, const subspace& s, int d, int v=0)
  {return h1->opmat_restricted(i,s,d,v);}
  smat s_opmat(int i, int d, int v=0)
  {return h1->s_opmat(i,d,v);}
  svec s_opmat_col(int i, int j, int v=0)
  {return h1->s_opmat_col(i,j,v);}
  smat s_opmat_cols(int i, const vec_i& jlist, int v=0)
  {return h1->s_opmat_cols(i,jlist,v);}
  smat s_opmat_restricted(int i, const ssubspace& s, int d, int v=0)
  {return h1->s_opmat_restricted(i,s,d,v);}
  long matdim(void)  {return h1->dimension;}
  scalar matden(void)  {return h1->denom1;}
  vector<long> eigrange(int i) {return h1->eigrange(i);}
  long dimoldpart(const vector<long> l);
protected:
  vec mvp;
  map<long,vec> mvlplusvecs, mvlminusvecs;
  oldforms* of;
  homspace *h1, *h1plus, *h1minus, *h1full;
  int j0;                   // data used for ap computation
  std::set<long> jlist;
public:
  long n1ds, j1ds;
  vector<newform> nflist;
  vector<int> nf_subset;
  newforms(long n, scalar mod, int disp)
    :level(n), modulus(mod), verbose(disp), of(0), h1(0), h1plus(0), h1minus(0), h1full(0)
  {
    ; //cout<<"In newforms constructor, level="<<n<<", modulus="<<mod<<", verbose="<<verbose<<endl;
  }
  ~newforms(void);
  void display(void) const;
  void display_modular_symbol_map(int check=0) const;
  void output_to_file(int binflag=1, int smallflag=0) const;
  void set_sign(int s) {sign=s;}
  int  get_sign() {return sign;}
  void makeh1(int s);
// add newform with basis b1, eiglist l to current list (b2 not used):
  void use(const vec& b1, const vec& b2, const vector<long> l);

  // find newforms using homology; ntp is number of eigenvalues to use
  // for oldforms, *not* the number computed via homology (use addap()
  // for that):
  void createfromscratch(int s, long ntp);

  // read newforms from file, if it exists, otherwise (perhaps) revert
  // to createfromscratch; if small_data_ok is set then we only need
  // the eigenvalues so are content with a small version of the file
  // (fewer ap and no additional data):
  void createfromdata(int s, long ntp, int create_from_scratch_if_absent=1,
		      int small_data_ok=0);

  // Compute homspace::projcoord, so projchain can be used
  // Replaces coord_vecs of homspace with projections onto eigenspaces
  // NB if #newforms>1 this MUST be re-called after any sorting of newforms
  void make_projcoord();

  // Look for a j0 such that nflist[i].bplus/bminus[j0]!=0 for all i, or a set of such j
  void find_jlist();

  // Create from one or a list of elliptic curves of the right conductor:
  void createfromcurve(int s, const CurveRed& C, int nap=500);
  void createfromcurves(int s, vector<CurveRed> Clist, int nap=500);
  // Lazy version which does not get the homology eigenvectors, only
  // sets the aq & ap (for use as oldforms):
  void createfromcurves_mini(vector<CurveRed> Clist, int nap=25);

  // read newforms from old-style data files (eigs/x$N and intdata/e$N):
  void createfromolddata();

  // Construct bases (homology eigenvectors) from eigenvalue lists:
  // flag controls what ::use() does with the nfs when found
  void makebases(int flag, int all_nf=1);

  // Construct H1 newforms, given H1+ and H1- newforms
  void merge(int all_nf=1);

  vector<long> apvec(long p);  // computes a[p] for each newform
  void addap(long last); // adds ap for primes up to the last'th prime

  // Sort newforms
  void sort(int oldorder=0);
  void sort_into_Cremona_label_order();
  void sort_into_LMFDB_label_order() {sort(0);}

  // To fix eigenvalues lists before/after recovering bases.  See
  // comments for the same named methods in the newform class for
  // details.
  void unfix_eigs();
  void refix_eigs();

  // for the i'th newform return the value of the modular symbol {0,r} (default) or {oo,r}
  rational plus_modular_symbol(const rational& r, long i=0, int base_at_infinity=0) const;
  rational minus_modular_symbol(const rational& r, long i=0, int base_at_infinity=0) const;
  pair<rational,rational> full_modular_symbol(const rational& r, long i=0, int base_at_infinity=0) const;

  // next three implemented in periods.cc

  // Given newform with no intdata, compute least real (part of)
  // period -- unless sfe=-1 and n=square, in which case return 0
  int get_real_period(long i, bigfloat& x, int verbose=0) const;
  // Given newform with no intdata, compute least imag (part of) period
  int get_imag_period(long i, bigfloat& y, int verbose=0) const;
  // Given all data, compute the periods as a Cperiods
  Cperiods getperiods(long i, int method=-1, int verbose=0);
  // Given all data & Cperiods, compute the curve (using fixc6 etc)
  Curve getcurve(long i, int method, bigfloat& rperiod, int verbose=0);

  // Attempt to compute and display the elliptic curve for each
  // newform in forms; return a (sub)list of newform indices where this failed.
  // If filename is not "no", output to file also
  vector<int> showcurves(vector<int> forms, int verbose, string filename);

  // next three implemented in pcprocs.cc

  // Computes x0, y0 (real & imag parts of periods) & a matrix which
  // gives these scaled by dotplus & dotminus.  rp_known is set if we
  // know x0 to be the least real part of a period (usually true).

  int find_matrix(long i, long dmax, int&rp_known, bigfloat&x0, bigfloat&y0);

  // Compute both periods from the (known) matrix [a,b;Nc,d] and
  // scaling factors dotplus, dotminus.  Return success flag.
  int get_both_periods(long i, bigfloat&x0, bigfloat&y0) const;

  // Given an imaginary period y1, finds a prime lminus =3(mod 4) and
  // <=lmax for which L(f,lminus,1) is nonzero and hence a multiple
  // mminus of y1.
  // if lmax==0 it carries on until a suitable lminus is found
  int find_lminus(long i, long lmax, const bigfloat& y1);
};

void output_to_file_no_newforms(long n, int binflag=1, int smallflag=0);
vector<long> eiglist(const newform& f, int oldorder=0);

/******************************************************************************
 To sort the newforms of level N from the order in which they are
 stored in newforms/x<N> into the correct order to match the "Cremona
 labels" of isogeny classes, we apply one of the following
 procedures. Note that to sort into LMFDB order, sort(0) suffices.

0: nf.sort(1) and then permute according to booknumber(N,i)
1: nf.sort(1)
2: nf.sort(0)
3: nf.unfix_eigs(); nf.sort(0); nf.refix_eigs();

Here:

 -  unfix_eigs() replaces the coefficient aq in aplist for q|N with the AL-eigenvalue wq;
 -  refix_eigs() reverse this;
 -  sort(1) sorts first by lexicographically sorting aqlist (AL-eigenvalues) in order +1,-1,
                  next  by lexicographically sorting aplist in order 0, +1, -1, +2, -2, ...;
 -  sort(0) sorts by lexicographically sorting aplist in order ...,-2,-1,0,1,2,...

The difference between cases 0 and 1 is that for N<=450 the order of
the newform files is essentially random, being the order in which the
newforms were found at a time when the strategy used in the code was
evolving steadily -- this was in the late 1980s, running batch jobs on
a remote mainframe, so rerunning those levels was not an easy option.
Instead, for N<=450 we have hard-wired the permutation taking the
order produced by sort(1) -- which is correct without further
adjustment for 450<N<130000 -- to the published order.

The required permutation is defined in eclib/libsrc/curvesort.cc and
accessed via the function i -> booknumber0(N,i), which is not the
identity for exactly 146 levels N between 56 and 450 inclusive.  The
correct i'th newform is number booknumber(N,i) in the stored list.

The difference between cases 2 and 3 is that in case 2, aplist contains
the p'th coefficient for all p, while in case 3 the q'th coefficient
for q|N is replaced by the AL-eigenvale wq.

Case 2 is the LMFDB ordering and is correct for all N>230000, as well
as some (but not all!) N between 130100 and 130200.  It should have
been used for all N>130000 but was not (my mistake): case 3 is
actually the order in which the newforms are found using the strategy
in place since level 130000, but the line sort(0) was omitted for the
code in error.

*******************************************************************************/

// utility to determine which sort method should be used, depending on
// the level, to recreate the "Cremona label" order of newforms.

inline int level_range(long N)
{
  if (N<=450)
    return 0;
  if (N<130000)
    return 1;
  if (N>230000)
    return 2;
  if ((N>130100)&&(N<130200)&&(N!=130144)&&(N!=130146)&&(N!=130150)&&(N!=130190)&&(N!=130192))
    return 2;
  return 3;
}


#endif
