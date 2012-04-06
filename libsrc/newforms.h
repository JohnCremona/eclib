// File NEWFORMS.H
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

#include "xsplit.h"   // which includes method.h

class newforms;
class jumps;

/* Data stored in a newform (and in data files newforms/x$N):
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
  scalar type;            // 2 for rectangular, 1 for triangular 
			  //  period lattice
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
  vec coordsplus, coordsminus;  // vector components of each freegen
  long denomplus, denomminus; 
  int j0; long fac;

  newform(void) {;}
  //  newform(const vec& v, const vector<long>& ap, newforms* nfs,long ind=-1);
  newform(const vector<int>& data, const vector<long>& aq, const vector<long>& ap, newforms* nfs);
  newform(const vec& vplus, const vec& vminus, const vector<long>& ap, newforms* nfs,long ind=-1);

  void add_more_ap(int nap);
  void display(void) const;
  // Testing function
  int check_expand_contract();
  // To fix eigenvalues lists after finding a newform
  void fixup_eigs();
  // To fix eigenvalues lists before/after recovering bases
  void unfix_eigs();
  void refix_eigs();
  // To find BSD ratio:
  void find_bsd_ratio();
  // To find projected coords:
  void find_coords_plus_minus();
  // To find cuspidal factors:
  void find_cuspidal_factors();
  // To find twisting primes:
  void find_twisting_primes();
  // To find deg(phi):
  void find_degphi();
  // To get matrix and scale factors when sign==0 ...
  void find_matrix();
};


class newforms :public level, splitter_base   {
  friend class newform;
private:
  int verbose; long maxdepth, cuspidal, sign;
  int basisflag;  // is set, then use() only sets bases for newforms
		  // already defined.
  mat opmat(int i, int d, int v=0) 
  {return h1->opmat(i,d,v);}
  mat opmat_restricted(int i, const subspace& s, int d, int v=0) 
  {return h1->opmat_restricted(i,s,d,v);}
  smat s_opmat(int i, int d, int v=0) 
  {return h1->s_opmat(i,d,v);}
  smat s_opmat_restricted(int i, const ssubspace& s, int d, int v=0) 
  {return h1->s_opmat_restricted(i,s,d,v);}
  long matdim(void)  {return h1->dimension;} 
  long matden(void)  {return h1->denom3;}
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
  newforms(long n, int disp) 
    :level(n), verbose(disp), of(0), h1(0), h1plus(0), h1minus(0), h1full(0) {;}
  ~newforms(void);
  void display(void) const;
  void display_modular_symbol_map(void) const;
  void output_to_file(int binflag=1) const;
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
  // to createfromscratch:
  void createfromdata(int s, long ntp, int create_from_scratch_if_absent=1);

  // Compute homspace::projcoord, so projchain can be used
  // Replaces coord_vecs of homspace with projections onto eigenspaces
  void make_projcoord();

  // Look for a j0 such that nflist[i].bplus/bminus[j0]!=0 for all i, or a set of such j
  void find_jlist();

  // Create from one or a list of elliptic curves of the right conductor:
  void createfromcurve(int s, CurveRed C, int nap=25);
  void createfromcurves(int s, vector<CurveRed> Clist, int nap=25);

  // read newforms from old-style data files (eigs/x$N and intdata/e$N):
  void createfromolddata();

  // Construct bases (homology eigenvectors) from eigenvalue lists:
  // flag controls what ::use() does with the nfs when found
  void makebases(int flag);

  // Construct H1 newforms, given H1+ and H1- newforms
  void merge();

  vector<long> apvec(long p);  // computes a[p] for each newform
  void addap(long last); // adds ap for primes up to the last'th prime

  // Sort newforms 
  void sort(int oldorder=0);
  // To fix eigenvalues lists before/after recovering bases
  void unfix_eigs();
  void refix_eigs();
  
  // for the i'th newform return the value of the modular symbol {0,r}
  rational plus_modular_symbol(const rational& r, long i=0) const;
  rational minus_modular_symbol(const rational& r, long i=0) const;
  pair<rational,rational> full_modular_symbol(const rational& r, long i=0) const;

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
  vector<int> showcurves(vector<int> forms, int verbose);

  // next three implemented in pcprocs.cc

  // Computes x0, y0 (real & imag parts of periods) & a matrix which
  // gives these scaled by dotplus & dotminus.  rp_known is set if we
  // know x0 to be the least real part of a period (usually true).

  int find_matrix(long i, long dmax, int&rp_known, bigfloat&x0, bigfloat&y0);

  // Compute both periods from the (known) matrix [a,b;Nc,d] and
  // scaling factors dotplus, dotminus.  Return success flag.
  int get_both_periods(long i, bigfloat&x0, bigfloat&y0);

  // Given an imaginary period y1, finds a prime lminus =3(mod 4) and
  // <=lmax for which L(f,lminus,1) is nonzero and hence a multiple
  // mminus of y1.
  // if lmax==0 it carries on until a suitable lminus is found
  int find_lminus(long i, long lmax, const bigfloat& y1);
};

void output_to_file_no_newforms(long n, int binflag=1);
vector<long> eiglist(const newform& f, int oldorder=0);
char* nf_filename(long n, char c);

