// File NEWSPACE.H: classes Newspace and Newform for newforms of any dimension
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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


#ifndef _NEWSPACE_H
#define _NEWSPACE_H      1

#include "homspace.h"
#include "curvesort.h" // for codeletter
#include "polys.h"     // for poly_cmp
#include "field.h"

class Newform;
class Newspace;

// class for a d-dimensional newform, defined by an irreducible factor
// of the characteristic polynomial of some splitting operator T
class Newform {
  friend class Newspace;
private:
  long N; // the level
  Newspace* nsp;    // pointer to "parent" class holding global info
  int index;       // index (starting from 1) of this newforms in the list of all
  string lab;      // label suffix (a,b,c,...)
  int d;      // dim(S)
  Field F0;   // Hecke field (original)
  Field F;    // Hecke field (reduced)
  FieldIso Fiso; // isomorphism from F0 to F (possibly identity)

  subspace_m S; // irreducible subspace of modular symbol space
  ZZ denom_abs; // absolute denominator of S
  mat projcoord; // used to computed eigenvalues of any operator
  modsym key_symbol; // nsp->H1->freemods[pivots(S)[1] -1]
  mat_m basis_change_matrix;
  ZZ basis_change_denominator;
  // To construct a field element from a 'raw' coord vector c and
  // denominator d, replace c by basis_change_matrix*c and d by
  // basis_change_denominator*d and cancel.

  int sfe; // Sign of functional equation (minus product of AL eigs), 0 if not known

  // Dict of T(P) eigenvalues of good primes P:
  map<long, FieldElement> aPmap;
  // max(p) for p in aPmap:
  long maxP;
  // Dict of W(Q) eigenvalues in {+1,-1} of bad primes Q, triv char only:
  map<long, int> eQmap;
  // List of coefficients in F, indexed by integers M.
  vector<FieldElement> aMlist;
  // List of traces of coefficients in F, indexed by integers M.
  vector<ZZ> trace_list;

  // Fill dict eQmap, if triv_char; if aPmap not already filled, first
  // compute ntp aP
  void compute_AL_eigs(int ntp=10, int verbose=0);
  // Compute aM for one M, assuming value for M'<M are all known and in Mlist
  FieldElement compute_aM(const long& M);

public:
  // Fill dict aPmap of eigenvalues of first ntp good primes; put max(P) into maxP
  void compute_eigs(int ntp=10, int verbose=0);
  // Assuming aPmap filled, fill aMlist (Fourier coefficients), and
  // trace_list (their traces) using all P in aPmap if aPmap
  // already filled, else first compute ntp aP.
  void compute_coefficients(int ntp=10, int verbose=0);
  // Compute aP and AL-eigs and Fourier coeffs
  void compute_eigs_and_coefficients(int ntp=10, int verbose=0);

  // constructor from ambient Newspace using one irreducible factor of char
  // poly of Newspace's T_mat
  Newform(Newspace* x, int ind, const ZZX& f, int verbose=0);
  // constructor from ambient Newspace (read from file)
  Newform(Newspace* x, int i, int verbose=0);

  // NB We do not use automatic copy constructor and assignment since
  // when data items which have field pointers are copied these must be
  // changed to pointers to fields in the new Newform not the old.

  // copy constructor
  Newform(const Newform& x);
  // assignment
  Newform& operator=(const Newform& x);

  // Return the number of this newform (counting from 1)
  int get_index() const { return index;}
  // Use after sorting to reset the numbers and variable names
  void set_index(int i);

  // Functions for computing eigenvalues of Hecke operators:

  // eigenvalue of a general operator:
  FieldElement eig(const matop& T);
  // eigenvalue of T_p or Q_p:
  FieldElement ap(const long& p)  {return eig(matop(p,N));}

  // eigenvalue +-1 of a scalar involution operator
  int eps(const matop& T);

  // coefficient in F of integral M from aMlist or computed (and stored
  // in aMmap) using multiplicative relations.
  FieldElement aM(const long& M);

  // eigenvalue of a prime from eQmap or aPmap if P is in there;
  // otherwise compute it.
  FieldElement eig_P(const long& P);
  // Principal eigenvalue of a linear combination of the above:
  FieldElement eig_lin_comb(const vector<long>& Plist, const vector<scalar>& coeffs, int verb=0);
  // Characteristic polynomial of such a linear combination:
  ZZX char_pol_lin_comb(const vector<long>& Plist, const vector<scalar>& coeffs, int verb=0);

  // Return the Hecke field
  Field field(int original=0) const {return (original? F0: F);}
  // Return the degree of the Hecke field
  int dimension() const {return d;}
  ZZX poly(int original=0) const {return (original? F0.poly(): F.poly());}
  string label_suffix() const {return lab;}
  string label() const;

  // Output basis for the Hecke field
  // Optionally aP and AL data too
  void display(int aP=0, int AL=0, int traces=0) const;
  // Display aP data
  void display_aP() const;
  // Display AL eigenvalues
  void display_AL() const;

  // Compute aPmap for first ntp primes if empty, and return it
  map<long, FieldElement> TP_eigs(int ntp, int verbose=0);
  // Compute eQmap if empty and return it, first computing aPmap for
  // first ntp primes if necessary
  map<long, int> AL_eigs(int ntp=10, int verbose=0);
  // return the list of traces
  vector<ZZ> traces() const {return trace_list;}

  // Filename for this Newform:
  string filename() const;
  // Output newform data:
  void output_to_file() const;
  // Input newform data. Returns 0 if data not available, else 1.
  int input_from_file(int verb=0);
};

// function to sort newforms of the same level, by (1) traces, (2)
// dimension, (3) min poly.  Before computing any eigenvalues (1) will
// not apply so only (2) and (3) are relevant; resorted after
// computing eigenvalues and traces, (2) and (3) are irrelevant.

struct newform_comparison {
  bool operator()(const Newform& f1, const Newform& f2)
  {
    // (1) Sort by lists of traces; if these have not been computed
    // yet, or if the character is nontrivial, the trace lists will be
    // empty and we will fall through to the final tests.
    vector<ZZ> traces1 = f1.traces(), traces2 = f2.traces();
    int t = traces1<traces2; // true e.g. if f1 has smaller absolute Hecke field degree
    if(t) return 1;
    t = traces1>traces2; // true e.g. if f1 has larger absolute Hecke field degree
    if(t) return 0;

    // We only get here when the traces have not been computed.

    // (2) Sort by dimension, i.e. the degree of the Hecke field,

    int s = f1.dimension() - f2.dimension();
    if(s) return (s<0); // true if f1 has smaller dimension

    // (3) Sort by min poly (smaller degree comes before larger)
    ZZX pol1 = f1.poly(), pol2 = f2.poly();
    return poly_cmp(pol1, pol2);
  }
};

extern newform_comparison newform_cmp;

// class for the collection of all d-dimensional newforms
class Newspace {
  friend class Newform;
private:
  int verbose;

  long N; // the level
  string level_label; // the level's label
  vector<long> Ndivs; // divisors of N
  vector<long> badprimes; // prime divisrs of N

  homspace* H1;    // the ambient modular symbol space at level N
  int dimH, cdimH; // full dimension, cuspidal dimension
  scalar dH;       // denominator of H1

  mat_m T_mat;  // matrix of splitting operator
  string T_name;  // name of splitting operator
  vector<ZZX> factors; // list of multiplicity-1 irreducible factor of charpoly(T)

  // Return the char poly of T on the new cuspidal subspace using the
  // oldspaces to obtain the old factors with correct multiplicities.

  ZZX new_cuspidal_poly(const vector<long>& Plist, const vector<scalar>& coeffs,
                        const gmatop &T);

  // Return true iff this combo of ops T has char poly on the new
  // cuspidal subspace which is squarefree and coprime to both the old
  // cuspidal poly and the full Eisenstein poly. f_new is set to the new
  // cuspidal poly.
  int valid_splitting_combo(const vector<long>& Plist, const vector<scalar>& coeffs,
                            const gmatop &T, ZZX& f_new);

  // Find a linear combination T of up to maxnp operators T_P with
  // coefficients up to maxc, for good p>=minp, whose char poly on the newspace is
  // squarefree and coprime to its char poly on the oldspace.  Set
  // split_ok=1 if successful else 0.
  void find_T(int maxnp, int maxc, int minp);

public:
  vector<Newform> newforms; // the newforms
  Newspace(void) :verbose(0) {;}
  // constructor from a homspace, looking for a splitting operator
  // using linear combinations of up to maxnp primes >=minp,
  // coefficients up to maxc
  Newspace(homspace* h1, int maxnp, int maxc, int minp=2, int verb=0);
  int split_ok; // records whether the constructor was able to find a splitting operator

  // constructor from file
  explicit Newspace(const long& level, int verb=0);

  mat_m heckeop(const long& P, int cuspidal=0, int dual=0);
  mat_m heckeop(const matop& T, int cuspidal=0, int dual=0) const;
  mat_m heckeop(const gmatop& T, int cuspidal=0, int dual=0) const;

  int ok() const {return split_ok;}
  int nforms() const {return newforms.size();}
  long level() const {return N;}
  string label() {return level_label;}

  string splitopname() const {return T_name;}

  // Return a list of the degrees of the Hecke fields
  vector<int> dimensions() const;

  // output all newforms: Dimension, Hecke field; optionally aP and AL data and traces
  void display_newforms(int aP=0, int AL=0, int traces=0) const;
  // sort the list of newforms using newform_cmp
  void sort_newforms();
  // return the list of newforms
  vector<Newform> the_newforms() const {return newforms;}

  // filename for Newspace
  string filename();
  // output data for this Newspace and each Newform
  void output_to_file();
  // Input Newspace data and newform data for each newform. Returns 0 if data missing, else 1.
  int input_from_file(const long& level, int verb=0);
};

// dict of Newspaces read from file
extern map<string,Newspace*> Newspace_dict;  // Key: label(N)
Newspace* get_Newspace(const long& N, int verb=0);

string newspaces_directory(int create_if_necessary=1);

#endif
