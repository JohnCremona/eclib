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
#include "eclib/newspace.h"

// Degree bound: fields of degree up to this will be reduced via
// polredabs (giving a canonical defining polynomial); above this only
// polredbest will be used.
const int POLREDABS_DEGREE_UPPER_BOUND = 10;

newform_comparison newform_cmp;

// When this is set, use the old spaces info read in from file, with
// the correct multiplicities, to compute the new char poly of any
// linear combination of Hecke operators.

// Otherwise, use recursion to find the new polys.

#define USE_OLD_SPACES

// Set this to use sparse matrix reduction

//#define SPARSE_KERNEL

Newform::Newform(Newspace* x, int ind, const ZZX& f, int verbose)
  :N(x->N), nsp(x), index(ind), lab(codeletter(ind-1)), d(deg(f)), sfe(0)
{
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
          cout << "denH = " << nsp->dH << endl;
        }
    }
  mat_m fT = evaluate(scale_poly_up(f, to_ZZ(nsp->dH)), nsp->T_mat);
  if (verbose>1)
    cout << "Computing ker(f(T)) where f(T) = \n" << fT << endl;

#ifdef SPARSE_KERNEL
  const ZZ modulus(default_modulus<ZZ>());
  ssubspace_m sS = kernel(smat_m(fT), modulus);
  mat_m bas = sS.bas().as_mat();
  ZZ den;
  int ok = liftmat(bas, modulus, bas, den);
  assert(ok);
  S = subspace_m(bas, sS.pivs(), den);
#else
  S = kernel(fT, 3); // method 3: echelon via FLINT
#endif

  if(dim(S)!=d)
    {
      cout<<"Problem: eigenspace has wrong dimension "<<dim(S)<<", not "<<d<<endl;
      exit(1);
    }
  denom_abs=to_ZZ(nsp->dH)*denom(S);
  key_symbol = nsp->H1->freemods[pivots(S)[1] -1];
  if (verbose)
    {
      cout << "Finished constructing subspace S of dimension " << d
           << ", relative denom = "<< denom(S)
           << ", absolute denom = "<< denom_abs << endl;
    }

  // compute projcoord, precomputed projections the basis of S
  projcoord = nsp->H1->coord * to_mat(basis(S));
  if(verbose>1)
    {
      cout << "H1->coord = " << nsp->H1->coord << endl;
      cout << "basis(S)  = " << basis(S) << endl;
      cout << "projcoord = " << projcoord << endl;
    }

  // To define the Hecke field we could use f itself and the
  // restriction of T to S, but if T is a complicated linear
  // combination of T_p we do better to use a single T_p if possible,
  // provided that it acts irreducibly on S.

  // If the newform is ratioinal there is nothing more to do

  if (d==1)
    {
      F = F0 = new Field();
      Fiso = FieldIso(*F);
      if (verbose)
        cout <<"Hecke field is Q" << endl;
      return;
    }

  Qmat A; ZZX cpA; string opname;
  int found=0;
  for (auto op: nsp->ops)
    {
      A = Qmat(transpose(restrict_mat(nsp->heckeop(op),S)), denom_abs);
      cpA = A.charpoly();
      if (IsIrreducible(cpA))
        {
          found=1;
          opname = op.name();
          break;
        }
    }
  if (found)
    {
      if (verbose && (opname != nsp->T_name))
        cout << "Using " << opname << " to define the Hecke field with polynomial\n" << ::str(cpA)
             << "\ninstead of " << nsp->T_name << " polynomial\n" << ::str(f) << endl;
    }
  else
    {
      A = Qmat(transpose(restrict_mat(nsp->T_mat,S)), denom_abs);
      cpA = A.charpoly();

      // Check that (scaled) charpoly(A) = f
      if (cpA != f)
        {
          cout<<endl;
          cout<<"Error: f(X) =            "<<fstring<<endl;
          cout<<"but charpoly(A) = " << str(cpA) << endl;
          exit(1);
        }
      if (verbose)
        cout << "No single T_p for p in " << nsp->Plist << " acts irreducibly on S, so using "
             << nsp->T_name << " to define the Hecke field with polynomial " << ::str(cpA) << endl;
    }

  // Now define the field, using A and its char poly

  F0 = new Field(A, basis_change_inverse, basis_change_matrix, lab, verbose>1);
  F = new Field();
  int canonical = (d<=POLREDABS_DEGREE_UPPER_BOUND);
  string F0pol = ::str(F0->poly());
  if (verbose)
    {
      cout << "Applying "
           << (canonical? "polredabs": "polredbest")
           << " to " << F0pol;
      if (!canonical)
        cout << " (not applying polredabs as degree is greater than "
             << POLREDABS_DEGREE_UPPER_BOUND << ")";
      cout  << endl;
    }
  Fiso = F0->reduction_isomorphism(lab, *F, canonical);
  if (verbose)
    {
      cout << "Reduced field is " << *F << endl;
      if (Fiso.is_nontrivial())
        {
          cout << " -- replacing original Hecke field with polynomial\n" << F0pol
               << "\n   with polredabs reduced field with polynomial\n" << ::str(F->poly()) << endl;
        }
      cout <<"Hecke field data:" << endl;
      cout << *F << endl;
    }
}

// Constructor which will read from file
Newform::Newform(Newspace* x, int i, int fill_aPmap, int fill_aMlist, int verbose)
  :N(x->N), nsp(x), index(i), lab(codeletter(i-1))
{
  if (verbose)
    cout << "Constructing Newform from file data "<< endl;
  if (!input_from_file(fill_aPmap, fill_aMlist, verbose))
    cerr << "Unable to read Newform " << lab << endl;
}


string Newform::label() const
{
  return nsp->level_label + lab;
}

// eigenvalue of a general principal operator:
FieldElement Newform::eig(const matop& T)
{
  vec_m apv = to_vec_m(nsp->H1->applyop_proj(T, key_symbol, projcoord));
  static const ZZ one(1);
  if (d==1)
    {
      return FieldElement(*F0, bigrational(apv[1], to_ZZ(denom_abs)));
    }
  else
    {
      FieldElement a(*F0, basis_change_matrix * Qvec(apv, denom_abs));
      assert(a.is_integral());
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
      if (verbose)
        cout << "In compute_coefficients("<<ntp<<"), aPmap is empty so computing eigs..."<<endl;
      compute_eigs(ntp, verbose);
      if (verbose)
        {
          cout << "Done\n";
          display(1,1,0);
          cout << "----------------------------------------------"<<endl;
        }
    }

  // Compute a_m for m up to maxP (the max p for which we have aP).
  // This will trigger computation and storing of all coefficients up
  // to maxP and their traces.
  aMlist.clear(); aMlist.reserve(maxP+1);
  trace_list.clear(); trace_list.reserve(maxP+1);
  if (verbose)
    cout << "Computing a_m for m <= " << maxP << endl;
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

  // Sort the newforms (by dimension and polynomial; we have no traces yet)
  sort_newforms();
}

// sort the list of newforms using newform_cmp
void Newspace::sort_newforms()
{
  // Sort the newforms (by traces, dimension, polynomial)
  std::sort(newforms.begin(), newforms.end(), newform_cmp);
  int i=1;
  for (Newform& f: newforms)
    f.set_index(i++);
}

// constructor from file
Newspace::Newspace(const long& level, int fill_aPmaps, int fill_aMlists, int verb)
{
  input_from_file(level, fill_aPmaps, fill_aMlists, verb);
}

// Compute the char poly of T on the new cuspidal subspace using the
// oldspaces to obtain the old factors with correct multiplicities.

ZZX Newspace::new_cuspidal_poly(const vector<long>& Plist, const vector<scalar>& coeffs,
                                const gmatop &T)
{
  // This will use the caches
  if (verbose>1)
    cout << "In new_cuspidal_poly() with T = " <<T.name() << endl;

  // Get/compute the char poly of T on the cuspidal subspace

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
      const Newspace* NSD = get_oldspace(D, verbose>1);
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
// cuspidal poly.
int Newspace::valid_splitting_combo(const vector<long>& Plist, const vector<scalar>& coeffs,
                                    const gmatop &T, ZZX& f_new)
{
  // f_new is the char poly of T on the new cuspidal subspace
  f_new = new_cuspidal_poly(Plist, coeffs, T);
  if (!IsSquareFree(f_new))
    {
      if (verbose)
        cout << "NO: not squarefree";
      return 0;
    }
  // f_full is the char poly of T on the full space (not just the
  // cuspidal subspace, or the new subspace:
  ZZX f_full = get_poly(N, T, 0, H1->modulus); // cuspidal=0

  // Hence the quotient f_other is the char poly of T on the
  // complement of the new cuspidal (or new cuspidal trivial char0
  // subspace.  We want f_new to be coprime to this:
  ZZX f_other = f_full / f_new;
  if (!AreCoprime(f_new, f_other))
    {
      if (verbose)
        cout << "NO: not coprime to old Hecke polynomial * Eisenstein polynomial";
      return 0;
    }
  if (verbose)
    cout << "YES: squarefree and coprime to old * Eisenstein Hecke polynomial";
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
  Plist = the_primes.getfirst(maxnp, N, minp);
  if (verbose)
    cout << "Trying linear combinations of T_p for p in " << Plist << endl;

  // insert maxnp T_P into ops:
  ops.resize(maxnp);
  std::transform(Plist.begin(), Plist.end(), ops.begin(),
                 [this](long p){return matop(p,N);});

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
  auto test_lc = [this, &T_op, &f](const vector<scalar>& lc)
  {
    vector<scalar> xlc(lc);
    xlc.insert(xlc.end(), ops.size()-lc.size(), scalar(0));
    T_op = gmatop(ops, xlc);
    T_name = T_op.name();
    if (verbose)
      cout << "Trying "<<lc<<" ("<<T_name<<")..."<<flush;
#ifdef USE_OLD_SPACES
    return valid_splitting_combo(Plist, xlc, T_op, f);
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

  // Now test random linear combos of length 2...maxnp
  for (int np=2; np<=maxnp; np++)
    {
      for (int i=0; (!split_ok) && (i<5*np); i++)
        {
          vector<int> lci = random_vector(np, -maxc, maxc);
          vector<scalar> lc(lci.size());
          std::transform(lci.begin(), lci.end(), lc.begin(), [](const int& x) {return scalar(x);});
          split_ok = test_lc(lc);
          if (split_ok)
            {
              if (verbose)
                cout<<" -- OK!"<<endl;
              break;
            }
          if (verbose)
            cout << " -- no good, continuing..." << endl;
        }
      if (split_ok)
        break;
    }

  if (split_ok)
    {
      if (verbose)
        cout << " Using operator " << T_name << " to split off newforms" << endl;
      T_mat = heckeop(T_op);
      if (verbose>1)
        cout << " Getting new cuspidal poly for " << T_name << " from cache" << endl;
      f = get_new_poly(N, T_op, 1, H1->modulus); // cuspidal=1 (cached)
      if (verbose>1)
        cout << "  New cuspidal poly for " << T_name << " from cache is " << str(f) << endl;
    }
  else
    return;

  factors.clear();
  NTL::vec_ZZX NTL_factors= SFFactor(f);
  ::sort(NTL_factors.begin(), NTL_factors.end(), poly_cmp);
  int nfactors = NTL_factors.length();
  if (verbose)
    cout << " New cuspidal Hecke polynomial for operator " << T_name
         << " has " << nfactors << " irreducible factors:"<<endl;
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
  if (d>1)
    {
      if (verbose) cout << "a_p computed, now LLL-reducing Hecke order..." << flush;
      // (1) just reduce the integral basis
      // HO.LLL_reduce();
      // (2) LLL-reduce the ap coefficients
      // // Get list of ap
      vector<FieldElement> aplist(aPmap.size());
      std::transform(aPmap.begin(), aPmap.end(), aplist.begin(),
                     [](const auto& pap){return pap.second;});
      HO.LLL_reduce(aplist);
      if (verbose) cout << "done." << endl;
    }

  compute_AL_eigs(ntp, verbose);      // e(Q) and a(Q) for bad Q

  // Store integral coeffs of ap (good and bad p):
  for (auto x: aPmap)
    aPmap_int_coords[x.first] = HO.integral_coords(x.second);

  // Compute and store a(M) and traces for M <= maxP
  compute_coefficients(ntp, verbose);
}

// Fill aPmap, dict of eigenvalues of first ntp good primes
void Newform::compute_eigs(int ntp, int verbose)
{
  // HO = Order(*F); // to initialise with equation order instead of:
  ZZ bound(0);
  if (d>POLREDABS_DEGREE_UPPER_BOUND) bound = ZZ(100000000);
  HO = MaximalOrder(F, bound); // initialise with maximal order (or approximation)
  ZZ index_orig = HO.get_order_index();
  if (d>1 && verbose)
      cout << "Hecke order initialised, contains the equation order with index "
           << index_orig <<endl;
  long p=2;
  primevar pr(ntp); // iterator over first ntp primes
  while(pr.ok())
    {
      p = pr;
      ++pr;
      if (!divides(p,N)) // compute AL eigs separately, later
        {
          if (verbose) cout << "Computing a_p for p = " << p << "..." << flush;
          FieldElement a_p = ap(p);
          if (verbose)
            {
              cout << "done";
              if (verbose>1) cout << ", a_p = " << a_p;
              cout<< endl;
            }
          assert(a_p.is_integral());
          aPmap[p] = a_p;
          if (d>1)
            {
              if (!HO.contains(a_p))
                {
                  if (verbose)
                    cout << "a_"<<p<<" not in current Hecke order (denominator " << HO.denom(a_p)
                         << "), extending..." << endl;
                  ZZ rel_index = HO.extend_by(a_p, 1); // 0: no need to check that a_p is integral
                  if (verbose)
                    {
                      cout << "Hecke order grows by index " << rel_index << endl;
                      if (verbose>1)
                        {
                          cout << "New basis is\n";
                          for (auto b: HO.get_basis()) cout << b << " = " << b.coords() << endl;
                        }
                    }
                }
            }
        }
    }
  maxP = p; // record last prime
  ZZ ind = HO.get_order_index();
  if (d>1 && verbose)
    {
      cout << "After computing a_p for primes up to " << maxP << ", Hecke order ";
      if (is_one(ind))
        cout <<"is the equation order Z[" << F->gen() << "]" << endl;
      else
        {
          cout << "contains the equation order Z[" << F->gen() << "] with index " << ind;
          if (ind!=index_orig)
            cout << " (after enlarging by index " << ind/index_orig <<")";
          cout<< endl;
        }
    }
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

// output info about the newform: dimension, sign, Hecke field and
// order; optionally aP and AL data; optionally also traces.

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
  ZZ ind = (d==1? ZZ(1) : HO.get_order_index());
  int is_eqn_order = HO.is_equation_order(); // not just index=1: basis must be standard
  FieldElement gen = F->gen();
  if (d>1)
    {
      if (is_eqn_order)
        {
          cout << " - Hecke order is equation order Z[" << gen <<"] with standard basis.";
        }
      else
        {
          if (is_one(ind))
            cout << " - Hecke order is equation order Z[" << gen <<"] but with non-standard basis.";
          else
            cout << " - Hecke order contains equation order Z[" << gen <<"] with index " << ind <<".";
        }
      cout << endl;
    }

  if (AL)
    {
      cout << endl;
      display_AL();
    }
  if (aP)
    {
      cout << endl;
      if (aPmap_int_coords.empty())
        cout << "No aP known" << endl;
      else
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

// Display aP data.

void Newform::display_aP() const
{
  ZZ ind = (d==1? ZZ(1) : HO.get_order_index());
  int is_eqn_order = HO.is_equation_order(); // not just index=1: basis must be standard
  if (!is_eqn_order)
    {
      string g = F->gen().str();
      const mat_m& pcm = HO.get_pcm();
      cout << "Coordinates of power basis in Hecke order:\n";
      for (int i=0; i<d; i++)
        {
          if (i==0)
            cout << "1";
          else
            {
              cout << g;
              if (i>1) cout << "^" << i;
            }
          cout << "\t = " << pcm.col(i+1) << endl;
        }
      cout << endl;
    }

  cout << "a_p for first " << aPmap_int_coords.size() << " primes:" << endl;
  for (auto x : aPmap_int_coords)
    {
      const long& p = x.first;
      const vec_m& aPcoords = x.second;
      cout << p << ":\t";
      if (d==1)
        cout << aPcoords[1];
      else
        cout << aPcoords;
      auto it = eQmap.find(p);
      if (it != eQmap.end())
        cout << "\t(AL eigenvalue = " << it->second << ")";
      cout << endl;
    }
}

// Display A-L eigenvalues
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
  out << nsp->label() << " " << lab << "\n";

  // Dimension
  out << dimension() << "\n";

  // Hecke field:
  out << F->str(1) << "\n";  // raw=1
  // cout << "Principal Hecke field output:\n" << F->str(1) << "\n";  // raw=1

  // Hecke order:
  out << HO.get_order_index() << "\n";
  vec_out(out, HO.get_pcm().get_entries(), " ", " ", " ");
  out << "\n\n";

  // Output A-L eigenvalues if computed, else output nothing
  vector<long> bads = nsp->badprimes;
  for (auto& Q: bads)
    {
      int eQ=0;
      if (!eQmap.empty())
        eQ = eQmap.at(Q);
      out << Q << " " << eQ << "\n";
    }
  out << "\n";

  // Output aP
  for (auto x: aPmap_int_coords)
    {
      const long& p = x.first;
      const vec_m& apcoords = x.second;
      out << p << " ";
      if (d==1)
        out << apcoords[1]; // no brackets
      else
        out << apcoords;
      out << "\n";
    }
}

// Input newform data (needs lab to be set to construct the filename).
// Returns 0 if data file missing, else 1
int Newform::input_from_file(int fill_aPmap, int fill_aMlist, int verb)
{
  // verb = 2;
  string fname = filename();
  if (verb)
    cout << "Reading newform " << lab << " from " << fname << " (verb="<<verb<<")"<<endl;
  ifstream fdata(fname.c_str());
  if (!fdata.is_open())
    {
      cerr << "Newform file " << fname << " missing" << endl;
      return 0;
    }

  if (fill_aMlist && !fill_aPmap)
    {
      if (verb>1) cout << "(switching on fill_aPmap since fill_aMlist is on)" << endl;
      fill_aPmap = 1;
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

  // read integral basis info

  ZZ ib;
  fdata >> ib; // index of the equation order
  if (verb>1)
    cout << "Read order index = " << ib << endl;
  mat_m bcm(d,d);
  fdata >> bcm;
  if (verb>1)
    cout << "Read bcm = \n" << bcm << endl;
  HO = Order(*F,ib,bcm);
  if (verb>1)
    cout << "Hecke order: " << HO << endl;

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
  maxP=0;
  long P; ZZ aP1;
  FieldElement aP(*F);
  vec_m aPcoords(d), un = vec_m::unit_vector(1,1);

  // read whitespace, so if there are no aP on file it does not try to read any
  fdata >> ws;
  // keep reading lines until end of file
  while (!fdata.eof())
    {
      fdata >> P;    // prime
      if (d==1)
        {
          fdata >> aP1;
          if (fill_aPmap)
            aP = (*F)(aP1);
          aPcoords = aP1 * un;
        }
      else
        {
          fdata >> aPcoords;
          if (fill_aPmap)
            aP = HO(aPcoords);
        }
      if (fill_aPmap)
        aPmap[P] = aP;
      aPmap_int_coords[P] = aPcoords;
      if (verb>1)
        {
          cout << "--> P = " << P << ": a_P = ";
          if (fill_aPmap)
            cout << aP << " = ";
          cout << aPcoords << endl;
        }
      if ((P>maxP) && !divides(P,N))
        maxP = P;
      // eat whitespace, including newline, so we can detect eof
      fdata >> ws;
    }
  if (verb>1)
    {
      cout << "After reading aP from " << fname <<":" << endl;
      display(1,1,0);
    }
  if (fill_aMlist)
    {
      // compute a_m and traces from a_p
      compute_coefficients();
    }
  if (verb)
    {
      cout << "After reading everything from " << fname <<":" << endl;
      display(1,1, fill_aMlist); // no traces unless aMlist has been filled.
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

int Newspace::input_from_file(const long& level, int fill_aPmaps, int fill_aMlists, int verb)
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
      Newform nf(this, i, fill_aPmaps, fill_aMlists, verb);
      if (verb>1)
        {
          cout << "After reading newform #" << i << " from file:" << endl;
          nf.display(1,1, fill_aMlists); // only display traces if each aMlist have been filled
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

mat_m Newspace::heckeop(const gmatop& T) const
{
  return to_mat_m(transpose(get_full_mat(N,T,H1->modulus)));
}

mat_m Newspace::heckeop(const matop& T) const
{
  return to_mat_m(transpose(get_full_mat(N,T,H1->modulus)));
}

mat_m Newspace::heckeop(const long& P)
{
  return heckeop(matop(P, N));
}

// dict of Newspaces read from file
map<string,Newspace*> Newspace_dict;  // Key: label(N)

Newspace* get_oldspace(const long& N, int verb)
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
  Newspace* NSP = new Newspace(N,
                               1, // fill_aPmap
                               0, // do not fill aMlist, trace_list
                               verb);
  Newspace_dict[Nlabel] = NSP;
  if (verb)
    cout << "Newspace at level " << Nlabel << " read from file and cached" << endl;
  return NSP;
}
