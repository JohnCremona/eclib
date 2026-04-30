// FILE NEWSPACE.CC: implementation of class Newspace for newforms of any dimension
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

#include <assert.h>
#include <random>
#include "eclib/newspace.h"

// Degree bound: fields of degree up to this will be reduced via
// polredabs (giving a canonical defining polynomial); above this only
// polredbest will be used.
const int POLREDABS_DEGREE_UPPER_BOUND = 15;

newform_comparison newform_cmp;

// When this is set, use the old spaces info read in from file, with
// the correct multiplicities, to compute the new char poly of any
// linear combination of Hecke operators.

// Otherwise, use recursion to find the new polys.

#define USE_OLD_SPACES

Newform::Newform(Newspace* x, int ind, const ZZX& f, int verbose)
  :N(x->N), nsp(x), index(ind), lab(codeletter(ind-1))
{
  d = deg(f);
  string fstring = str(f);
  if (verbose)
    {
      cout << "Constructing Newform from factor f = "<< fstring <<" of degree "<<d<<endl;

      // Compute f(T); since T is scaled by dH and f(X) is not, we
      // evaluate dH^d*f(X/dH) at T; that is, we scale the coefficient of
      // X^i by dH^(d-i):

      if (verbose>1)
        {
          cout << "Finding kernel of f(T)..."<<endl;
          cout << " f = " << fstring <<endl;
          cout << "denH = " << nsp->dH << endl;
        }
    }
  mat_m fT = evaluate(scale_poly_up(f, to_ZZ(nsp->dH)), nsp->T_mat);
  if (verbose>1)
    cout << "Computing ker(f(T)) where f(T) = \n" << fT << endl;
  S = kernel(fT);
  if(dim(S)!=d)
    {
      cout<<"Problem: eigenspace has wrong dimension "<<dim(S)<<", not "<<d<<endl;
      exit(1);
    }
  denom_abs=to_ZZ(nsp->dH)*denom(S);
  key_symbol = nsp->H1->freemods[pivots(S)[1] -1];
  if (verbose)
    {
      cout<<"Finished constructing subspace S of dimension "<<d
          <<", absolute denom = "<<denom_abs<<endl;
      cout<<"Computing A, the restriction of T..." <<flush;
    }
  Qmat A(transpose(restrict_mat(nsp->T_mat,S)), denom_abs);
  if(verbose)
    cout<<"done. Checking its char poly..."<<flush;

  // Check that (scaled) charpoly(A) = f
  ZZX cpA = A.charpoly();
  if (cpA != f)
    {
      cout<<endl;
      cout<<"Error: f(X) =            "<<fstring<<endl;
      cout<<"but charpoly(A) = " << str(cpA) << endl;
      exit(1);
    }

  if(verbose)
    {
      cout<<"done."<<endl;
      if (verbose>1)
        cout<<"A (the matrix of T restricted to S) = \n"
            << A << endl;
      cout<<"f(X) is the min poly of A"<<endl;
    }

  // compute projcoord, precomputed projections the basis of S
  projcoord = nsp->H1->coord * to_mat(basis(S));
  if(verbose>1)
    {
      cout << "H1->coord = " << nsp->H1->coord << endl;
      cout << "basis(S)  = " << basis(S) << endl;
      cout << "projcoord = " << projcoord << endl;
    }

  if (d==1)
    {
      F = F0 = new Field();
      Fiso = FieldIso(*F);
    }
  else
    // Compute Hecke field from matrix A if degree d>1
    {
      string var = codeletter(index-1);
      F0 = new Field(A, basis_change_inverse, basis_change_matrix, var, verbose>1);
      F = new Field();
      int canonical = (d<=POLREDABS_DEGREE_UPPER_BOUND);
      if (verbose)
        {
          cout << "Applying "
               << (canonical? "polredabs": "polredbest");
          if (verbose>1)
            cout << " to field "
                 // << F0 << " --> "
                 << *F0;
          cout  << endl;
          if (!canonical)
            {
              cout << "Not applying polredabs as degree is greater than "
                   << POLREDABS_DEGREE_UPPER_BOUND <<endl;
            }
        }
      Fiso = F0->reduction_isomorphism(var, *F, canonical);
      if (verbose)
        cout << "Reduced field is " << *F << endl;
      if (Fiso.is_nontrivial() && verbose)
        {
          cout << "[replacing original Hecke field with polynomial " << ::str(F0->poly())
               << " with polredabs reduced field with polynomial " << ::str(F->poly()) << "]" << endl;
        }
    }
  if (verbose)
    {
      cout <<"Hecke field data:" << endl;
      cout << *F << endl;
    }

  sfe = 0;   // means unknown
}

// Constructor which will read from file
Newform::Newform(Newspace* x, int i, int verbose)
  :N(x->N), nsp(x), index(i), lab(codeletter(i-1))
{
  if (verbose)
    cout << "Constructing Newform from file data "<< endl;
  if (!input_from_file(verbose))
    cerr << "Unable to read Newform " << lab << endl;
}


string Newform::label() const
{
  return nsp->level_label + lab;
}

// eigenvalue of a general principal operator:
FieldElement Newform::eig(const matop& T)
{
  // cout << "Matrix of "<<T.name()<<" is\n" << nsp->H1->calcop_restricted(T, S, 0, 0) << endl;
  // cout << "key_symbol = " << key_symbol <<endl;
  // cout << "applyop_proj("<<T.name()<<") to this with projcoord gives \n";

  vec_m apv = to_vec_m(nsp->H1->applyop_proj(T, key_symbol, projcoord));
  // cout << "ap vector = " << apv <<endl;
  static const ZZ one(1);
  if (d==1)
    {
      return FieldElement(*F0, bigrational(apv[1], to_ZZ(denom_abs)));
    }
  else
    {
      FieldElement a(*F0, basis_change_matrix * Qvec(apv, denom_abs));
      return (Fiso.is_identity()? a : Fiso(a));
    }
}

// eigenvalue +-1 of a scalar involution operator
int Newform::eps(const matop& T) // T should be a scalar +- identity
{
  FieldElement e = eig(T);
  if (e.is_one())
    return 1;
  if (e.is_minus_one())
    return -1;
  cerr << "eps(" << T.name() << ") returns " << e << ", not +1 or -1" << endl;
  exit(1);
}

// eigenvalue of a prime from eQmap or aPmap if P is in there;
// otherwise compute it.
FieldElement Newform::eig_P(const long& P)
{
  // If P is bad, return AL eigenvalue
  auto it1 = eQmap.find(P);
  if (it1!=eQmap.end())
    return (*F)(it1->second);
  // If P is good and in aPmap, return stored Tp eigenvalue
  auto it2 = aPmap.find(P);
  if (it2!=aPmap.end())
    return it2->second;
  // If P is good and not in aPmap, compute, store and return Tp eigenvalue
  FieldElement aP = ap(P);
  aPmap[P] = aP;
  return aP;
}

//#define DEBUG_COEFFS

// Compute aM for one M, assuming value for M'<M are all known and in aMlist
FieldElement Newform::compute_aM(const long& M)
{
#ifdef DEBUG_COEFFS
  cout<<"Computing aM for M = " << M << endl;
#endif
  // largest M for which we already have the coefficient:
  long maxM = aMlist.size()-1;

  // Just return the ctored value if possible:
  if (M<=maxM)
    return aMlist[M];

  // Otherwise compute it
  if (M<1)
    return (*F)(0);
  if (M==1)
    return (*F)(1);

  // Now we assume that aMlist contains aM[M'] for proper divisors
  long p = primdiv(M); // smallest prime divisor
  long M1 = M;
  long e = divide_out(M1, p); // now M = M1*p^e

  if (M1>1) // use easy multiplicativity
    {
#ifdef DEBUG_COEFFS
      cout << " - using multiplicativity with factors " << M1 << " * " << M/M1 <<endl;
#endif
      return aMlist[M1] * aMlist[M/M1];
    }

#ifdef DEBUG_COEFFS
  cout << " - prime power case, p = " << p << ", e = " << e << endl;
#endif

  // Now M = p^e is a prime power.  aPmap has the aP for bad p as well
  // as good, from compute_AL_eigs()
  FieldElement aP = aPmap[p];
#ifdef DEBUG_COEFFS
  cout << " - a_p = " << aP << endl;
#endif
  if (e==1)
    return aP;

  // Now  e>=2, use recursion
  M1 = M/p;
  FieldElement a = aP*aMlist[M1];
  return (divides(p,N)? a : a - (*F)(p) * aMlist[M1/p]);
}

//#define DEBUG_COEFFS
// coefficient in F of integral ideal M from aMmap or computed (and
// stored in aMmap) using multiplicative relations.
FieldElement Newform::aM(const long& M)
{
  // largest M for which we already have the coefficient:
  long maxM = aMlist.size()-1;
#ifdef DEBUG_COEFFS
  cout<<"Fetching aM for M = " << M << endl;
#endif

  // Compute and store a_M' for M'<=M if necessary,(doing nothing
  // unless M>maxM:

  for (long M1=maxM+1; M1<=M; M1++)
    {
      FieldElement a = compute_aM(M1);
#ifdef DEBUG_COEFFS
      cout<<" computed a_"<<M1<<" = " << a << endl;
#endif
      aMlist.push_back(a);
      trace_list.push_back(a.trace().num());
    }
#ifdef DEBUG_COEFFS
  cout<<" returning a_"<<M<<" = " << aMlist[M] << endl;
#endif
  return aMlist[M];
}

// Assuming aPmap filled, fill aMmap (Fourier coefficients)
void Newform::compute_coefficients(int ntp, int verbose)
{
  if (aPmap.empty())
    {
      cout << "In compute_coefficients("<<ntp<<"), aPmap is empty so computing eigs..."<<endl;
      compute_eigs(ntp, verbose);
      cout << "Done\n";
      display(1,1,0);
      cout << "----------------------------------------------"<<endl;
    }

  // Compute a_m for m up to maxP (the max p for which we have aP).
  // This will trigger computation and storing of all coefficients up
  // to maxP and their traces.
  aMlist.clear(); aMlist.reserve(maxP+1);
  trace_list.clear(); trace_list.reserve(maxP+1);
  aM(maxP); // value discarded, just done for the side-effect
}

// Principal eigenvalue of a linear combination of the above:
FieldElement Newform::eig_lin_comb(const vector<long>& Plist, const vector<scalar>& coeffs,
                                   int verb)
{
  FieldElement a(*F);
  auto Pj = Plist.begin();
  auto cj = coeffs.begin();
  while (Pj!=Plist.end())
    {
      scalar c = *cj++;
      long P = *Pj++;
      if (c!=0)
        a += eig_P(P) * to_ZZ(c);
    }
  return a;
}

// Characteristic polynomial of such a linear combination:
ZZX Newform::char_pol_lin_comb(const vector<long>& Plist, const vector<scalar>& coeffs,
                               int verb)
{
  return eig_lin_comb(Plist, coeffs, verb).charpoly();
}

Newspace::Newspace(homspace* h1, int maxnp, int maxc, int minp, int verb)
  : verbose(verb), H1(h1)
{
  N = H1->N;
  level_label = to_string(N);
  if(verbose)
    cout << "In Newspace constructor at level " << N << endl;
  dimH = H1->h1dim();
  cdimH = H1->h1cuspdim();
  dH = H1->h1cdenom();
  badprimes = pdivs(N);
  Ndivs = posdivs(N, badprimes);
  Ndivs.pop_back(); // excude N itself

  if(verbose && dH>1)
    cout<<"H has dimension "<<dimH<<", cuspidal dimension "<<cdimH<<", denominator "<<dH<<endl;

  if (cdimH==0)
    {
      if (verbose)
        cout << "Cuspidal dimension is 0, so newspace is trivial" << endl;
      return;
    }

  // Find the splitting operator
  find_T(maxnp, maxc, minp);
  // Construct the newforms if that succeeded
  if (split_ok)
    {
      if (verbose)
        {
          cout << "Using splitting operator T = " << T_name << endl;
          cout << "Number of irreducible components is " << factors.size() << endl;
        }
      int i=1;
      newforms.reserve(factors.size());
      for (const auto& f: factors)
        {
          // The index i is stored in the newform on construction, but
          // these will change when we sort them
          newforms.push_back(Newform(this, i, f, verbose));
          ++i;
        }
    }
  else
    // abort if not
    cout << "Unable to find a suitable splitting operator!" << endl;

  // Sort the newforms (by character, dimension, polynomial; we have no traces yet)
  sort_newforms();
}

// sort the list of newforms using newform_cmp
void Newspace::sort_newforms()
{
  // Sort the newforms (by character, dimension, polynomial)
  std::sort(newforms.begin(), newforms.end(), newform_cmp);
  int i=1;
  for (Newform& f: newforms)
    f.set_index(i++);
}

// constructor from file
Newspace::Newspace(const long& level, int verb)
{
  input_from_file(level, verb);
}

// Compute the char poly of T on the new cuspidal subspace using the
// oldspaces to obtain the old factors with correct multiplicities.

// If triv_char=0: requires oldspace data for forms with all genus
// characters, so will only work over fields where this is
// implemented.

// If triv_char=1: only requires oldspace data for forms with
// trivial genus characters, so works over all fields.
ZZX Newspace::new_cuspidal_poly(const vector<long>& Plist, const vector<scalar>& coeffs,
                                const gmatop &T)
{
  // This will use the caches
  if (verbose>1)
    cout << "In new_cuspidal_poly() with T = " <<T.name() << endl;

  // Get/compute the char poly of T on the cuspidal subspace /
  // cuspidal, trivial character, subspace
  ZZX f_all = get_poly(N, T, 1, H1->modulus); // cuspidal=1
  if (verbose>1)
    {
      cout << "f_all = " << str(f_all) << endl;
      display_factors(f_all);
    }
  ZZX f_old;
  set(f_old); // set = 1
  for (auto D: Ndivs) // Ndivs contains all *proper* D|N
    {
      const Newspace* NSD = get_Newspace(D, verbose>1);
      if (NSD->nforms()==0)
        continue;
      long M = N/D;
      if (verbose>1)
        cout << "Divisor D = " << D << " has " << NSD->nforms() << " newforms" <<endl;
      vector<long> divs = posdivs(M);
      int mult = divs.size();
      if (verbose>1)
        cout << "# divisors of N/D is " << mult <<endl;
      for (auto form: NSD->newforms)
        {
          ZZX f_D = form.char_pol_lin_comb(Plist, coeffs, verbose>1);
          if (verbose>1)
            {
              cout << "T's poly for " << form.label_suffix() << " is f_D = " << str(f_D) << endl;
              display_factors(f_D);
            }
          for (int i=0; i<mult; i++)
            f_old *= f_D;
        }
    }
  if (verbose>1)
    {
      cout << "The product of all old polys is " << str(f_old) << endl;
      display_factors(f_old);
    }
  //essentially f_new = f_all / f_old; but checking divisibility
  ZZX f_new, rem;
  DivRem(f_new, rem, f_all, f_old);
  if (!IsZero(rem))
    {
      cout << "Problem in Newspace::new_cuspidal_poly() at level "<<level_label<<endl;
      cout << "Operator " << T.name() << " has full poly" << endl;
      display_factors(f_all);
      cout << "and old poly" << endl;
      display_factors(f_old);
      cout << " which does not divide the full poly. " << endl;
      exit(1);
    }

  if (T.is_simple()) // cache this new poly
    {
      string NT = NTkey(N,T);
      if (verbose)
        cout<<"Caching new cuspidal poly for " << NT <<endl; // << ": " << str(f_new) << endl;
      new_cuspidal_poly_dict[NT] = f_new;
    }
  return f_new;
}

// Return true iff this combo of ops T has char poly on the new
// cuspidal subspace which is squarefree and coprime to both the old
// cuspidal poly and the full Eisenstein poly. f_new is set to the new
// cuspidal poly. If triv_char=1 then same for the char poly of T on
// the new cuspidal trivial-character subspace, in which case f_new
// must also be coprime to the char poly of T on the new cuspidal
// nontrivial char subspace.
int Newspace::valid_splitting_combo(const vector<long>& Plist, const vector<scalar>& coeffs,
                                    const gmatop &T, ZZX& f_new)
{
  // f_new is the char poly of T on the new cuspidal subspace
  f_new = new_cuspidal_poly(Plist, coeffs, T);
  if (!IsSquareFree(f_new))
    {
      if (verbose>0)
        {
          cout << "\n NO: new Hecke polynomial " // <<str(f_new)
               << " for " << T.name() << " is not squarefree" << endl;
          // display_factors(f_new);
        }
      return 0;
    }
  // f_full is the char poly of T on the full space (not just the
  // cuspidal subspace, or the new subspace, or (if triv_char==1) just
  // the trivial character subspace:
  ZZX f_full = get_poly(N, T, 0, H1->modulus); // cuspidal=0

  // Hence the quotient f_other is the char poly of T on the
  // complement of the new cuspidal (or new cuspidal trivial char0
  // subspace.  We want f_new to be coprime to this:
  ZZX f_other = f_full / f_new;
  if (!AreCoprime(f_new, f_other))
    {
      if (verbose>0)
        {
          cout << "\n NO: new Hecke polynomial " // <<str(f_new)
               << " for " << T.name()
               << " is not coprime to old Hecke polynomial * Eisenstein polynomial ";
          if (verbose>1)
            {
              cout << str(f_other)<<endl
                   << " (full polynomial is "<<str(f_full)<<")";
            }
          cout <<endl;
        }
      return 0;
    }
  if (verbose>0)
    cout << "\n YES: new Hecke polynomial for " << T.name()
         << " is squarefree and coprime to old * Eisenstein Hecke polynomial" << endl;
  return 1;
}

// Find a linear combination T of up to maxnp operators T_P with
// coefficients up to maxc, for good p>=minp, whose char poly on the new cuspidal
// subspace is squarefree and coprime to its char polys on the
// oldspace and non-cuspidal subspace.
//
// Set split_ok=1 if successful else 0.

void Newspace::find_T(int maxnp, int maxc, int minp)
{
  split_ok = 0;
  // fill Plist with the first maxnp good primes >= minp
  vector<long> Plist = the_primes.getfirst(maxnp, N, minp);
  if (verbose)
    cout << "Trying linear combinations of T_p for p in " << Plist << endl;
  vector<matop> ops(maxnp);
  // insert maxnp T_P into ops:
  std::transform(Plist.begin(), Plist.end(), ops.begin(),
                 [this](long p){return matop(p,N);});

  BoundedWeightVectorGenerator combo_generator(maxnp, maxc);

  if (verbose)
    {
      cout << "Trying linear combinations with coefficients up to "<<maxc;
      if (verbose)
        {
          cout << " of (";
          for (auto T: ops)
            cout << " " << T.name();
          cout << ")";
        }
      cout << endl;
    }
  gmatop T_op;
  ZZX f;

  // local function to test a linear combination (lc) of ops,
  // returning success flag.  If successful, assigns T_op, T_name, f
  auto test_lc = [this, ops, Plist, &T_op, &f](const vector<scalar>& lc)
  {
    T_op = gmatop(ops, lc);
    T_name = T_op.name();
    if (verbose)
      cout << "Trying "<<lc<<": "<<T_name<<"..."<<flush;
#ifdef USE_OLD_SPACES
    return valid_splitting_combo(Plist, lc, T_op, f);
#else
    return test_splitting_operator(N, T_op, H1->modulus, verbose>1);
#endif
  };

  // First test individual ops
  vector<scalar> lc(maxnp); // all 0
  for (int i=0; i<maxnp; i++)
    {
      lc[i] = 1;
      split_ok = test_lc(lc);
      if (split_ok)
        {
          if (verbose)
            cout<<"OK!"<<endl;
          break;
        }
      lc[i] = 0;
      if (verbose)
        cout << " no good, continuing..." << endl;
    }

  // Now test 100 random linear combos
  for (int i=0; (!split_ok) && (i<100); i++)
  {
    split_ok = test_lc(random_vector(maxnp, -maxc, maxc));
    if (split_ok)
      {
        if (verbose)
          cout<<"OK!"<<endl;
        break;
      }
    if (verbose)
      cout << " no good, continuing..." << endl;
  }

  if (split_ok)
    {
      if (verbose)
        cout << " OK: using operator " << T_name << " to split off newforms" << endl;
      T_mat = heckeop(T_op, 0, 1); // not cuspidal,  dual
      if (verbose)
        cout << " Getting new cuspidal poly for " << T_name << " from cache" << endl;
      f = get_new_poly(N, T_op, 1, H1->modulus); // cuspidal=1 (cached)
      if (verbose)
        cout << "  New cuspidal poly for " << T_name << " from cache is " << str(f) << endl;
    }
  else
    return;

  factors.clear();
  if (verbose)
    {
      cout << " New cuspidal Hecke polynomial for operator " << T_name
           <<" is "<<str(f)<<endl;
    }
  NTL::vec_ZZX NTL_factors= SFFactor(f);
  ::sort(NTL_factors.begin(), NTL_factors.end(), poly_cmp);
  int nfactors = NTL_factors.length();
  if (verbose)
    cout<<"Irreducible factors:"<<endl;
  for(int i=0; i<nfactors; i++)
    {
      ZZX fi = NTL_factors[i];
      if (verbose)
        cout<<(i+1)<<":\t"<<str(fi)<<"\t(degree "<<deg(fi)<<")"<<endl; 
     factors.push_back(fi);
    }
  return;
}

// Compute aP and AL-eigs and Fourier coeffs
void Newform::compute_eigs_and_coefficients(int ntp, int verbose)
{
  compute_eigs(ntp, verbose);         // a(P) for good P
  compute_AL_eigs(ntp, verbose);      // e(Q) and a(Q) for bad Q
  compute_coefficients(ntp, verbose); // a(M) and traces for N(M)<=max N(P)
}

// Fill aPmap, dict of eigenvalues of first ntp good primes
void Newform::compute_eigs(int ntp, int verbose)
{
  long p=2;
  primevar pr(ntp); // iterator over first ntp primes
  while(pr.ok())
    {
      p = pr;
      ++pr;
      if (!divides(p,N)) // compute AL eigs separately, later
        {
          if (verbose) cout << "Computing a_p for p = " << p << "..." << flush;
          aPmap[p] = ap(p);
          if (verbose) cout << "done, a_p = " << aPmap[p] << endl;
        }
    }
  maxP = p; // record last prime
} // end of compute_eigs

// Fill dict eQmap *after* aPmap, and set sfe
void Newform::compute_AL_eigs(int ntp, int verbose)
{
  // Compute AL-eigs to fill map<long, FieldElement> eQmap and also
  // the entries in aPmap indexed by bad primes.

  if (aPmap.empty())
    compute_eigs(ntp, verbose);

  sfe = -1; // since the sign is *minus* the product of the AL eigs

  if (verbose)
    cout << "Computing W(Q) eigenvalues for Q in " << nsp->badprimes << endl;

  FieldElement zero(*F);
  for (auto q: nsp->badprimes)
    {
      if (verbose)
        cout << "Computing W(q) eigenvalue for q = " << q << "..." << flush;

      int  eq = eps(matop(q,N));

      if (verbose)
        cout << " eigenvalue " << eq << endl;

      // Store the AL eigenvalue in eQmap, update sfe and set aq to be
      // minus the AL eigenvalue if q||N else 0
      eQmap[q] = eq;
      sfe *= eq;
      aPmap[q] = (divides(q*q,N)? zero: (*F)(-eq));

    } // end of loop over bad primes q
} // end of Newform::compute_AL_eigs()

// output basis for the Principal Hecke field and character of one newform
// If full, also output multiplicative basis for the full Hecke field
// Optionally aP and AL (if trivial char) data too
// Optionally traces too
// Optionally principal eigs too
void Newform::display(int aP, int AL, int traces) const
{
  cout << "Newform #" << index << " (" << label() << ")" << endl;

  // Information about the dimension (degree of Hecke field):

  cout << " - Dimension: "<<d<<endl;

  // Information about sign of functional equation:

  if (sfe!=0)
    cout << " - Sign of functional equation = " << (sfe>0? "+1" : "-1") << endl;

  // Information about Hecke field:

  cout << " - Hecke field k_f = " << *F << endl;

  if (AL)
    {
      cout << endl;
      display_AL();
    }
  if (aP)
    {
      cout << endl;
      display_aP();
    }
  if (traces)
    {
      cout << endl;
      cout << "First " << trace_list.size()-1 << " traces: "; // << trace_list << endl;
      copy(trace_list.begin()+1, trace_list.end(), ostream_iterator<ZZ>(cout, ", "));
      cout << "..." << endl;
    }
}

// Display aP data (trivial char or C4 fields)
//#define CHECK_TRACES
//#define BASES
void Newform::display_aP() const
{
  if (aPmap.empty())
    {
      cout << "No aP known" << endl;
      return;
    }
#ifdef BASES
  Qmat bc;
  if (d>1)
    {
      bc = denom_abs * basis_change_inverse * Fiso.matrix().inverse();
      FieldElement gen = F->gen();
      FieldElement gen_power = F->operator()(1);
      for (int i=0; i<d; i++)
        {
          cout << gen_power << "\t --> "
               << bc * gen_power.coords() << endl;
          if (i < d-1) gen_power *= gen;
        }
      cout << endl;
    }
#endif
  cout << "a_p for first " << aPmap.size() << " primes:" << endl;
  for (auto x : aPmap)
    {
      FieldElement aP = x.second;
      cout << x.first << ":\t";
#ifdef BASES
      if (d>1) cout << bc * aP.coords() << "\t\t";
#endif
      cout << aP;
#ifdef CHECK_TRACES
      bigrational taP = aP.trace();
      cout << " (trace = " << taP << ")";
#endif
      auto it = eQmap.find(x.first);
      if (it != eQmap.end())
        {
          cout << "\t(AL eigenvalue = " << it->second << ")";
        }
      cout << endl;
    }
}

// Display A-L eigenvalues (trivial char or C4 fields)
void Newform::display_AL() const
{
  cout << "Atkin-Lehner eigenvalues:" << endl;
  for (auto x: eQmap)
    cout << x.first << ":\t" << x.second << endl;
}

map<long, FieldElement> Newform::TP_eigs(int ntp, int verbose)
{
  if (aPmap.empty())
    compute_eigs(ntp, verbose);
  return aPmap;
}

map<long, int> Newform::AL_eigs(int ntp, int verbose)
{
  if (eQmap.empty())
    compute_AL_eigs(ntp, verbose);
  return eQmap;
}

// return name of Newspace directory; if create_if_necessary, creates
// the directory if it does not yet exist
string newspaces_directory(int create_if_necessary)
{
  string dirname = getenv_with_default("NSP_DIR", "newspaces");
  if (create_if_necessary)
    {
      int res = std::system(("mkdir -p "+dirname).c_str());
      if (res)
        {
          cerr << "mkdir -p "<<dirname<<" failed, writing newspace files in current directory"<<endl;
          dirname = ".";
        }
    }
  return dirname;
}

// filename
string Newform::filename() const
{
  return newspaces_directory() + "/" + label();
}

string Newspace::filename()
{
  return newspaces_directory() + "/" + label();
}

// Use after sorting to reset the numbers and variable names
void Newform::set_index(int i)
  {
    index = i;
    lab = codeletter(i-1);
    // On creation from scratch, F0 exists and F is a
    // polredabs/polredbest isomorphic field, but after reading from a
    // file only F is set.
    F0->set_var(lab+string("0"));
    F->set_var(lab);
  }

// newform file output
void Newform::output_to_file() const
{
  ofstream out;
  out.open(filename().c_str());

  // Level, letter-code:
  out << nsp->label() << " " << lab << endl;

  // Dimension
  out << dimension() << endl;

  // Hecke field:
  out << F->str(1) << endl;  // raw=1
  // cout << "Principal Hecke field output:\n" << F->str(1) << endl;  // raw=1

  out << endl;

  // Output A-L eigenvalues if computed, else output nothing
  vector<long> bads = nsp->badprimes;
  for (auto& Q: bads)
    {
      int eQ=0;
      if (!eQmap.empty())
        eQ = eQmap.at(Q);
      out << Q << " " << eQ << endl;
    }
  out << endl;

  // Output aP
  for (auto x: aPmap)
    {
      FieldElement aP = x.second;
      out << x.first << " " << aP.str(1) << endl;
    }
}

// Input newform data (needs lab to be set to construct the filename).
// Returns 0 if data file missing, else 1
int Newform::input_from_file(int verb)
{
  string fname = filename();
  if (verb)
    cout << "Reading newform " << lab << " from " << fname << " (verb="<<verb<<")"<<endl;
  ifstream fdata(fname.c_str());
  if (!fdata.is_open())
    {
      cerr << "Newform file " << fname << " missing" << endl;
      return 0;
    }
  // Check field, level, letter-code:
  string dat;
  fdata >> dat;
  assert (dat==nsp->label());
  fdata >> dat;
  assert (dat==lab);

  // Dimension:
  fdata >> d;
  if (verb>1)
    cout << "--> dim = " << d << endl;

  // Hecke field (we only read F; F0 and Fiso are not defined):
  F = new Field();
  F->read(fdata);

 if (verb>1)
    {
      cout << "After reading Hecke field, before reading eigenvalues:" << endl;
      display();
    }

  fdata >> ws;
  sfe = -1;
  //cout << "Now reading eQ for Q in " << nsp->badprimes <<endl;
  // A-L eigenvalues (will read nothing if level is (1))
  for (auto Q: nsp->badprimes)
    {
      if (verb>1)
        cout << "reading AL eig for Q = " << Q << endl;
      long iQ;
      fdata >> iQ;
      if (verb>1)
        cout << "--> Q = " << Q << ", file prime is " << iQ << endl;
      if (iQ!=Q)
        cerr << "!!! Q = " << Q << " but file prime is " << iQ << endl;
      assert (iQ==Q);

      int eQ;
      fdata >> eQ;
      eQmap[Q] = eQ;
      sfe *= eQ;
      if (verb>1)
        cout << "AL eigenvalue for " << Q << " is " << eQ << endl;
    }
  if (verb>1)
    {
      cout << "After reading eQ, before reading aP:" << endl;
      display(0, 1, 0);
    }

  // aP
  long P; maxP=0;
  FieldElement aP(*F);
  // read whitespace, so if there are no aP on file it does not try to read any
  fdata >> ws;
  // keep reading lines until end of file
  while (!fdata.eof())
    {
      fdata >> P // prime label
            >> aP   // eigenvalue data
            >> ws;  // eat whitespace, including newline
      aPmap[P] = aP;
      if (verb>1)
        cout << "--> P = " << P
             << ": a_P = " << aP
             << endl;
      if ((P>maxP) && !divides(P,N)) maxP = P;
    }
  if (verb>1)
    {
      cout << "After reading aP from " << fname <<":" << endl;
      display(1,1,0);
    }
  // compute a_m and traces from a_p
  compute_coefficients();
  if (verb)
    {
      cout << "After reading everything from " << fname <<":" << endl;
      display(1,1,1);
    }
  return 1;
}

void Newspace::output_to_file()
{
  ofstream out;
  out.open(filename().c_str());

  // Field, level, number of newforms:
  out << label() << " " << newforms.size() << endl;
  if (newforms.empty())
    {
      out.close();
      return;
    }

  // Dimensions:
  vector<int> dims = dimensions();
  for (auto d: dims) out << " " << d;
  out << endl;
  out.close();

  // Output each newform data to a separate file:
  for (const auto& f: newforms)
    f.output_to_file();

}

int Newspace::input_from_file(const long& level, int verb)
{
  N = level;
  level_label = to_string(N);
  if (verb)
    cout << "In Newspace::input_from_file() with N="<< N << endl;
  badprimes = pdivs(N);
  Ndivs = posdivs(N, badprimes);
  Ndivs.pop_back(); // excude N itself

  string fname = filename();
  ifstream fdata(fname.c_str());
  if (!fdata.is_open())
    {
      cerr << "File " << fname << " not available for Newspace input" << endl;
      return 0;
    }
  string dat;
  fdata >> dat;
  assert (dat==level_label);
  int nnf;
  fdata >> nnf;
  if (verb)
    cout << "-> Level "<<level_label<<" has " << nnf << " newforms in "<<fname<<endl;
  if (nnf==0)
    {
      fdata.close();
      return 1;
    }

  vector<int> dims(nnf);
  fdata >> dims;
  if (verb>1)
    cout << "-> dims: " << dims <<endl;

  newforms.reserve(nnf);
  for (int i=1; i<=nnf; i++)
    {
      if (verb)
        cout << "About to read newform #" << i << " from file" << endl;
      Newform nf(this, i, verb);
      if (verb>1)
        {
          cout << "After reading newform #" << i << " from file:" << endl;
          nf.display(1,1,1);
        }
      newforms.push_back(nf);
    }

  if (verb>1)
    {
      cout << "Finished reading Newspace data with " << nnf << " newforms from " << fname << endl;
      for (int i=0; i<nnf; i++)
        cout << "#" << (i+1) << ": " << "dim = " << dims[i] << endl;
    }
  fdata.close();
  return 1;
}

// output basis for the Hecke field and character of all newforms
void Newspace::display_newforms(int aP, int AL, int traces) const
{
  for ( auto& nf : newforms)
    {
      nf.display(aP, AL, traces);
      cout<<endl;
    }
}

// Return a list of the degrees of the Hecke fields
vector<int> Newspace::dimensions() const
{
  vector<int> dims(newforms.size());
  std::transform(newforms.begin(), newforms.end(), dims.begin(),
                 [](const Newform& nf){return nf.dimension();});
  return dims;
}

mat_m Newspace::heckeop(const gmatop& T, int cuspidal, int dual) const
{
  return to_mat_m(H1->calcop(T, cuspidal, dual, 0)); // 0 display
}

mat_m Newspace::heckeop(const matop& T, int cuspidal, int dual) const
{
  return to_mat_m(H1->calcop(T, cuspidal, dual, 0)); // 0 display
}

mat_m Newspace::heckeop(const long& P, int cuspidal, int dual)
{
  return heckeop(matop(P, N), cuspidal, dual);
}

// dict of Newspaces read from file
map<string,Newspace*> Newspace_dict;  // Key: label(N)

Newspace* get_Newspace(const long& N, int verb)
{
  string Nlabel = to_string(N);
  if (Newspace_dict.find(Nlabel) != Newspace_dict.end())
    {
      if (verb)
        cout << "Newspace at level " << Nlabel << " retrieved from cache" << endl;
      return Newspace_dict[Nlabel];
    }
  if (verb)
    cout << "Newspace at level " << Nlabel << " not in cache, reading from file..." << endl;
  Newspace* NSP = new Newspace(N, verb);
  Newspace_dict[Nlabel] = NSP;
  if (verb)
    cout << "Newspace at level " << Nlabel << " read from file and cached" << endl;
  return NSP;
}

// Function to generate a random integer vector of a given size wit
// entries taken uniformly from [minv..maxv], optionally repeating
// until the vector is primitive

vector<scalar> random_vector(size_t size, scalar minv, scalar maxv, int primitive)
{
  if (size==1 && primitive==1 && minv<=1 and maxv>=1)
    return {1};
  // We use static in order to instantiate the random engine and the
  // distribution once only.  It may provoke some thread-safety
  // issues.
  std::uniform_int_distribution<int> distribution(minv, maxv);
  static std::default_random_engine generator;

  std::vector<int> v(size);
  int cont = 0;
  while (cont==0 || (primitive && cont>1))
    {
      std::generate(v.begin(), v.end(),
                    [&distribution]() { return distribution(generator); });
      cont = std::accumulate(v.begin(), v.end(), 0,
                             [](const int& x, const int& y) {return gcd(x,y);});
    }
  return v;
}
