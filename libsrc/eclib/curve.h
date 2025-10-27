// curve.h: declarations of elliptic curve classes
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
 
// originally adapted from Elliptic.h by Oisin McGuiness

// allow for multiple includes
#ifndef _ECLIB_ELLIPTIC_
#define _ECLIB_ELLIPTIC_

#include <map>
#include "bigrat.h"

class Curve; class Curvedata; class CurveRed;
class Point;
class IsogenyClass;

//general test:
int valid_invariants(const ZZ& c4, const ZZ& c6); //true if valid
void c4c6_to_ai(const ZZ& c4, const ZZ& c6, 
		ZZ& a1, ZZ& a2, ZZ& a3, ZZ& a4, 
		ZZ& a6, 
		ZZ& b2, ZZ& b4, ZZ& b6, ZZ& b8);
void c4c6_to_ai(const ZZ& c4, const ZZ& c6, 
		ZZ& a1, ZZ& a2, ZZ& a3, ZZ& a4, 
		ZZ& a6);
void minimise_c4c6(const ZZ& c4, const ZZ& c6, const ZZ& discr, 
                   ZZ& newc4, ZZ& newc6, ZZ& newdiscr, ZZ& u);

//
// base class for bare elliptic curve
//

class Curve{ 
friend class Point;
friend class IsogenyClass;
public:
  int isnull() const
    {return ((a1==0)&&(a2==0)&&(a3==0)&&(a4==0)&&(a6==0));}
  void getai(ZZ& aa1, ZZ& aa2, ZZ& aa3,
             ZZ& aa4, ZZ& aa6) const
    {aa1=a1; aa2=a2; aa3=a3; aa4=a4; aa6=a6; }
  vector<ZZ> a_invariants() const
  { return vector<ZZ>{a1,a2,a3,a4,a6}; }

// input and output
  void output(ostream& os) const
    { os<<"["<<a1<<","<<a2<<","<<a3<<","<<a4<<","<<a6<<"]";}
//N.B. No spaces in output for ease of input to GP

  void input(istream& is);
  void tex_print(ostream &) const ;
  // puts out TeX-ed equation of curve; never been used

// constructors 
  Curve(void)  :a1(0),a2(0),a3(0),a4(0),a6(0)  {;}
  Curve(const ZZ& c4, const ZZ& c6); //init by invariants
                                 //check valid for elliptic curve
                                 //if not create null curve
  Curve(const ZZ& aa1, const ZZ& aa2, const ZZ& aa3,
        const ZZ& aa4, const ZZ& aa6)
    :a1(aa1),a2(aa2),a3(aa3),a4(aa4),a6(aa6) {}
  Curve(long aa1, long aa2, long aa3, long aa4, long aa6)
    :a1(aa1),a2(aa2),a3(aa3),a4(aa4),a6(aa6) {}
  Curve(const Curve& c)
    : a1(c.a1), a2(c.a2), a3(c.a3), a4(c.a4), a6(c.a6)
      {}
  Curve(const bigrational& j); // one curve with this j-invariant
  void operator=(const Curve& c)
    { a1=c.a1; a2=c.a2; a3=c.a3; a4=c.a4; a6=c.a6; }
// no destructor is needed

// equality tests -- inherited, so you can compare Curvedata with Curve
  int operator==(const Curve& f) const
    {return ((a1==f.a1)&&(a2==f.a2)&&
             (a3==f.a3)&&(a4==f.a4)&&(a6==f.a6));}
  int operator!=(const Curve& f) const
    {return ((a1!=f.a1)||(a2!=f.a2)||
             (a3!=f.a3)||(a4!=f.a4)||(a6!=f.a6));}

protected:
  ZZ a1 ;
  ZZ a2 ;
  ZZ a3 ;
  ZZ a4 ;
  ZZ a6 ;
} ;


//
// derived class for curve with computed invariants
//  NB minimalize() is a member function, flagged if ever called; by
//  default not called, but see constructors for how to force.
//

class Curvedata : public Curve{
friend class CurveRed;  // bug fix; ridiculous, is a derived class
friend class IsogenyClass;
public:
  Curvedata() {discr_factored=0;}
  Curvedata(const ZZ& aa1, const ZZ& aa2, const ZZ& aa3,
        const ZZ& aa4, const ZZ& aa6, int min_on_init);
  /*
  Curvedata(const bigrational& qa1, const bigrational& qa2, 
	    const bigrational& qa3, const bigrational& qa4, 
	    const bigrational& qa6, ZZ& scale);
  */
  Curvedata(const vector<bigrational>& qai, ZZ& scale);
  Curvedata(const Curve& c, int min_on_init);
  Curvedata(const ZZ& cc4, const ZZ& cc6, int min_on_init);
  Curvedata(const Curvedata& c);
  Curvedata(const Curvedata& c, int min_on_init);
       // nb compiler cannot generate copy because constructor from
       // curve overrides.
       // But default assign, destruct suffice
  void operator=(const Curvedata& c);
  void minimalize();  // Changes self in situ
  void factor_discr()
  {if(!discr_factored){ the_bad_primes=pdivs(discr); discr_factored=1; }}
  Curvedata minimalize(ZZ& u, ZZ& r, ZZ& s, ZZ& t) const;
     // Self unchanged and returns transformation
  void transform(const ZZ& r, const ZZ& s, const ZZ& t);   // NB  u = 1;
                            // the more general case is not implemented here
  void output(ostream& os) const;
  void input(istream& is);
  long get_ntorsion();     // implemented in points.cc

  void getbi(ZZ& bb2, ZZ& bb4, ZZ& bb6, ZZ& bb8) const
    {bb2=b2; bb4=b4; bb6=b6; bb8=b8; }
  void getci(ZZ& cc4, ZZ& cc6) const
    {cc4=c4; cc6=c6; }
  friend inline ZZ getb2(const Curvedata& c) {return c.b2; }
  friend inline ZZ getb4(const Curvedata& c) {return c.b4; }
  friend inline ZZ getb6(const Curvedata& c) {return c.b6; }
  friend inline ZZ getb8(const Curvedata& c) {return c.b8; }
  friend inline ZZ getc4(const Curvedata& c) {return c.c4; }
  friend inline ZZ getc6(const Curvedata& c) {return c.c6; }
  friend inline ZZ getdiscr(const Curvedata& c) {return c.discr; }
  friend inline bigrational j_invariant(const Curvedata& c) {return bigrational(power(c.c4,3),c.discr); }
  friend inline int getconncomp(const Curvedata& c) {return c.conncomp; }
  friend inline vector<ZZ> getbad_primes(Curvedata& c)
    {
      if(!c.discr_factored) c.factor_discr();
      return c.the_bad_primes;
    }
  // NB the is_minimal function returns 0 when minimization has not
  // been done; the curve may still be minimal
  friend int is_minimal(const Curvedata& c) {return c.minimal_flag;}
protected:
  ZZ b2 ;
  ZZ b4 ;
  ZZ b6 ;
  ZZ b8 ;
  ZZ c4 ;
  ZZ c6 ;
  ZZ discr ;
  int minimal_flag;  // 0 if .minimalize() has not been called
  int discr_factored; // 0 if discr has not yet been factored
  vector<ZZ> the_bad_primes; //prime divisors of discriminant
  int conncomp ;    // number of components (1 or 2)
  long ntorsion; // 0 if .gettorsion() not called
} ;

// function to find "optimal x shift" of a given curve
Curvedata opt_x_shift(const Curvedata& C, ZZ& k);

// CurveRed class call Tates algorithm as constructor,
// stores the info as member variables

// class Kodaira_code just holds an int which "codes" the type as follows:
// (this coding originally from R.G.E.Pinch)
//
// Im                 -> 10*m
// I*m                -> 10*m+1
// I, II, III, IV     -> 1, 2, 3, 4
// I*, II*. III*, IV* -> 5, 6, 7, 8
//
class  Kodaira_code {
public:
  int code;
//
  Kodaira_code(int k = 0) : code(k) {;}
  Kodaira_code(const Kodaira_code& c) : code(c.code) {;}
  void operator=(int k) {code=k;}
  void operator=(const Kodaira_code& c) {code=c.code;}
  friend ostream& operator<<(ostream& os, const Kodaira_code& c);
};

// utility function for converting Kodaira codes to the Pari coding

// Kodaira Symbol        My coding    Pari Coding

// I0                    0              1
// I*0                   1             -1
// Im  (m>0)             10*m           m+4
// I*m (m>0)             10*m+1        -(m+4)
// II, III, IV           2, 3, 4        m
// II*. III*, IV*        7, 6, 5       -m

int PariKodairaCode(Kodaira_code Kod);

class Reduction_type {
public:
  int ord_p_discr;
  int ord_p_N;
  int ord_p_j_denom;
  Kodaira_code Kcode;  // NB the constructor makes this from an int
  int c_p;
  int local_root_number;
//
  Reduction_type(int opd=0, int opN=0, int opj=0, int kc=1, int cp=1, int rno=0)
    : ord_p_discr(opd), ord_p_N(opN), ord_p_j_denom(opj), Kcode(kc), c_p(cp), local_root_number(rno)
    {}
};

inline ostream& operator<<(ostream& os, const Reduction_type& R);


class CurveRed : public Curvedata {
friend class IsogenyClass;
protected:
  map<ZZ,Reduction_type> reduct_array;  // maps p -> its reduction type
  ZZ N;                      //the conductor
public:
  CurveRed() : Curvedata(), N(0) {}
  explicit CurveRed(const Curvedata& E);  // construct by Tate's algorithm
             // arg E need not be minimal, but the reduced form will be
  ~CurveRed();
  CurveRed(const CurveRed& E);
  void operator=(const CurveRed& E);
  ZZ conductor() {return N;}
  // The full display function is not const, since if called it will
  // compute and set the local root numbers if necessary
  void display(ostream& os); // full output
  void output(ostream& os) const;  // just the curvedata and conductor

  // Sort key for sorting lists of curves (LMFDB ordering):
  // (1) conductor
  // (2) list of ap for good p < NP_SORT
  // (3) a1,a2,a3,4,a6
  vector<ZZ> sort_key(const int NP_SORT=25) const;
  int operator<(const CurveRed& E) const {return sort_key()<E.sort_key();}
private:
  // functions for setting local root numbers:
  int neron(long p, int kod); // p = 2 or 3
  void setLocalRootNumber(const ZZ& p);
  void setLocalRootNumber2();
  void setLocalRootNumber3();
  void setLocalRootNumber_not_2_or_3(const ZZ& p);

public:
  // member access functions:
  friend inline vector<ZZ> getbad_primes(const CurveRed& c)
  {return c.the_bad_primes; }
  friend inline ZZ getconductor(const CurveRed& c) {return c.N; }
  int ord_p_discr(const ZZ& p) const;
  int ord_p_N(const ZZ& p) const;
  int ord_p_j_denom(const ZZ& p) const;
  int c_p(const ZZ& p) const;
  vector<ZZ> all_cp() const;
  ZZ prodcp() const;
  int LocalRootNumber(const ZZ& p) const;
  int GlobalRootNumber() const;

  friend Kodaira_code getKodaira_code(const CurveRed& c, const ZZ& p);
  // the returned value casts as a character array; to use coded as int,
  // say declared Kodaira_code Kc, just use public member Kc.code

  // Trace of Frobenius (via pari if p is good for long p)
  // (or 0 for additive reduction, +1 for split multiplicative, -1 for nonsplit)
  // These are not constant methods as a call to LocalRootNumber may be needed.
  long ap(long p) const;
  ZZ ap(const ZZ& p) const;

  // The local Tamagawa number.  Use p=0 for reals
  friend ZZ local_Tamagawa_number(const CurveRed& c, const ZZ& p);
  // The local Tamagawa exponent -- same as Tamagawa number unless the
  // component group is (2,2).  Use p=0 for reals
  friend ZZ local_Tamagawa_exponent(const CurveRed& c, const ZZ& p);
  // The global Tamagawa exponent, i.e. the lcm of the exponents of
  // the component groups at all bad primes (including infinity if
  // real_too is 1), which is the lcm of the local Tamagawa exponents.
  // So (with no further knowledge of the MW group) we know that m*P
  // is in the good-reduction subgroup for all P, with this m.
  friend ZZ global_Tamagawa_exponent(const CurveRed& c, int real_too);

  int has_good_reduction_outside_S(const vector<ZZ>& S)
  {
    return is_S_unit(N, S);
  }
};

// Here the CurveRed parameter is not const since a call to
// LocalRootNumber may have to compute and store it
inline ZZ Trace_Frob(CurveRed& c, const ZZ& p) {return c.ap(p);}
inline long Trace_Frob(CurveRed& c, const long& p) {return c.ap(p);}

inline int getord_p_discr(const CurveRed& c, const ZZ& p) {return c.ord_p_discr(p);}
inline int getord_p_N(const CurveRed& c, const ZZ& p) {return c.ord_p_N(p);}
inline int getord_p_j_denom(const CurveRed& c, const ZZ& p) {return c.ord_p_j_denom(p);}
inline int getc_p(const CurveRed& c, const ZZ& p) {return c.c_p(p);}
inline vector<ZZ> all_cp(const CurveRed& c) {return c.all_cp();}
inline ZZ prodcp(const CurveRed& c) {return c.prodcp();}
inline int LocalRootNumber(const CurveRed& c, const ZZ& p) {return c.LocalRootNumber(p);}
inline int GlobalRootNumber(const CurveRed& c) {return c.GlobalRootNumber();}

// The global Tamagawa number, = product of local ones.
ZZ global_Tamagawa_number(const CurveRed& c, int real_too);

// Tamagawa primes: primes dividing any Tamagawa number
vector<long> tamagawa_primes(const CurveRed& C, int real_too);

inline ostream& operator<<(ostream& os, const Curve& c)
{
  c.output(os);
  return os;
}

inline ostream& operator<<(ostream& os, const Curvedata& c)
{
  c.output(os);
  return os;
}

inline ostream& operator<<(ostream& os, const CurveRed& c)
{
  c.output(os);
  return os;
}

inline istream& operator>>(istream& is, Curve& c)
{
   c.input(is);
   return is ;
}

inline istream& operator>>(istream& is, Curvedata& c)
{
   c.input(is);
   return is ;
}
   //Reads the curve (ai) and computes the rest.

inline int GlobalRootNumber(const Curvedata& E)
{
  CurveRed C(E);
  return GlobalRootNumber(C);
}

// Quadratic twist of an elliptic curve (returns minimal model)
CurveRed QuadraticTwist(const CurveRed& E, const ZZ& D);
inline CurveRed QuadraticTwist(const CurveRed& E, long D)
{ return QuadraticTwist(E, ZZ(D)); }

// Given a list of elliptic curves E, and one discriminant D, return the
// list of quadratic twists of the curves by D
vector<CurveRed> QuadraticTwists(const vector<CurveRed>& EE, const ZZ& D);
inline vector<CurveRed> QuadraticTwists(const vector<CurveRed>& EE, long D)
{ return QuadraticTwists(EE, ZZ(D)); }

// Given a list of elliptic curves E, and one prime p, return the
// list of quadratic twists of the curves by:
// +p if p=1 (mod 4)
// -p if p=3 (mod 4)
// -4, 8 and -8 if p=2

vector<CurveRed> PrimeTwists(const vector<CurveRed>& EE, const ZZ& p);
inline vector<CurveRed> PrimeTwists(const vector<CurveRed>& EE, long p)
{ return PrimeTwists(EE, ZZ(p)); }

// Given a list of elliptic curves, and a list of primes, return a
// list of all quadratic twists of the curves by discriminants supported on
// those primes (including the original curves)

vector<CurveRed> AllTwists(const vector<CurveRed>& EE, const vector<ZZ>& PP);
inline vector<CurveRed> AllTwists(const vector<CurveRed>& EE, vector<long> PP)
{ return AllTwists(EE, ZZify(PP)); }

// end of file: curve.h

#endif
