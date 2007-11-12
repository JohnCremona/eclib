// vec.h: declaration of integer vector classes
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
 
// Should not be included directly by user: include vector.h instead

// We cannot use one macro to prevent this being included twice since
// we want it to be included twice by vector.h!  We must use two
// different ones

class svec;
class subspace;

class vec;
vec iota(scalar n);                      // (1,2,...,n)

class vec {
friend class svec;
friend class smat;
friend class smat_elim;
friend class vec_m;
friend class mat;
friend class subspace;
public:
    // constructors
        vec(long n=0);
        vec(long n, scalar* arr);
        vec(const vec&);                       // copy constructor
        ~vec();                                   // destructor
     // member functions & operators
        void init(long n=0);                 // (re)-initializes 
        vec& operator=(const vec&);         // assignment
        scalar& operator[](long i) const;            // the i'th component 
        vec& operator+=(const vec&);
        void addmodp(const vec&, scalar);
        vec& operator-=(const vec&);
        vec& operator*=(scalar);
        vec& operator/=(scalar);
        vec slice(long,long=-1) const;           // returns subvec
        vec operator[](const vec&) const;   // subscript composition
        void set(long i, scalar x);                  // sets v[i]=x
        void add(long i, scalar x);                  // v[i]+=x
        scalar sub(long i) const;                    // same as v[i] (no ref)
	const scalar* get_entries()const {return entries;}
     // non-member (friend) functions and operators
        friend long dim(const vec&);                  // the dimension
        friend scalar operator*(const vec&, const vec&);   // dot product
        friend scalar operator*(const svec&, const vec&); 
        friend vec operator*(const mat& m, const vec& v);
        friend int operator==(const vec&, const vec&);
        friend int operator!=(const vec&, const vec&);
        friend int trivial(const vec&);                  // v==zerovec?
  // add/sub row i of mat to v (implemented in mat.cc)
        friend void add_row_to_vec(vec& v, const mat& m, long i);
        friend void sub_row_to_vec(vec& v, const mat& m, long i);
        friend ostream& operator<< (ostream&s, const vec&);
        friend istream& operator>> (istream&s, vec&);
        friend vec iota(scalar n);                      // (1,2,...,n)
        friend scalar vecgcd(const vec&);
        friend void swapvec(vec& v, vec& w);
        friend int member(scalar a, const vec& v);//tests if a=v[i] for some i
        friend mat restrict(const mat& m, const subspace& s, int cr=0);
        friend mat_m restrict(const mat_m& m, const msubspace& s);
        friend mat prestrict(const mat& m, const subspace& s, scalar pr, int cr=0);
        friend mat_m prestrict(const mat_m& m, const msubspace& s, const bigint& pr);

// Implementation
private:
       long d;
       scalar * entries;
};

// Declaration of non-member, non-friend functions

vec operator+(const vec&);                   // unary
vec operator-(const vec&);                   // unary
vec operator+(const vec&, const vec&);
vec addmodp(const vec&, const vec&, scalar);
vec operator-(const vec&, const vec&);
inline vec operator*(scalar, const vec&);       // componentwise
vec operator/(const vec&, scalar);       // componentwise
void makeprimitive(vec& v);
void elim(const vec& a, vec& b, long pos);
void elim1(const vec& a, vec& b, long pos);
void elim2(const vec& a, vec& b, long pos, scalar lastpivot);
vec reverse(vec& order);
vec express(const vec& v, const vec& v1, const vec& v2);
vec lift(const vec& v, scalar pr);  //lifts a mod-p vector to a rational
                                        //and scales to a primitive vec in Z.
int liftok(vec& v, scalar pr);  //lifts a mod-p vector in place to a rational
                               //and scales to a primitive vec in Z; returns success
scalar dotmodp(const vec& v1, const vec& v2, scalar pr);

// inline function definitions

inline long dim(const vec& v) {return v.d;}

inline int operator!=(const vec& v, const vec& w) { return !(v==w);}

inline vec operator+(const vec& v) { return v;}

inline vec operator-(const vec& v) { return (-1)*v;}

inline vec operator+(const vec& v1, const vec& v2)
{ vec ans(v1); ans+=v2; return ans;}

inline vec addmodp(const vec& v1, const vec& v2, scalar pr)
{ vec ans(v1); ans.addmodp(v2,pr); return ans;}

inline vec operator-(const vec& v1, const vec& v2)
{ vec ans(v1); ans-=v2; return ans;}

inline vec operator*(scalar scal, const vec& v)
{ vec ans(v); ans*=scal; return ans;}

inline vec operator/(const vec& v, scalar scal)
{ vec ans(v); ans/=scal; return ans;}

inline void makeprimitive(vec& v)
{ scalar g=vecgcd(v); if (g>1) v/=g;}

inline void elim(const vec& a, vec& b, long pos)
{ (b*=a[pos])-=(b[pos]*a);}

inline void elim1(const vec& a, vec& b, long pos)
{ (b*=a[pos])-=(b[pos]*a); makeprimitive(b);}

inline void elim2(const vec& a, vec& b, long pos, scalar lastpivot)
{ ((b*=a[pos])-=(b[pos]*a))/=lastpivot;}

