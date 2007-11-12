// mmatrix.h: declarations of multiprecision integer matrix class
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
 
#if     !defined(_MMATRIX_H)
#define _MMATRIX_H      1       //flags that this file has been included

#include "mvector.h"
#include "matrix.h"

class mat_m {
friend class msubspace;
public:
     // constructors
        mat_m(long nr=0, long nc=0);
        mat_m(const mat_m&);                    // copy constructor
        mat_m(const mat_i&);
        mat_m(const mat_l&);
     // destructor
        ~mat_m();

     // member functions & operators
        void init(long nr=0, long nc=0);
        mat_m& operator=(const mat_m&);          // assignment with copy
        bigint& operator()(long i, long j) const;   // returns ref to (i,j) entry
        mat_m slice(long,long,long=-1,long=-1) const;// returns submatrix
        bigint sub(long i, long j) const;             // returns the (i,j) entry
        vec_m row(long i) const;                // returns row i (as a vector)
        vec_m col(long j) const;                // returns col j (as a vector)
        void set(long i, long j, const bigint& x);     // sets the (i,j) entry to x
        void add(long i, long j, const bigint& x);  // adds x to the (i,j) entry  
        void setrow(long i, const vec_m& v);
        void setcol(long i, const vec_m& v);
        void swaprows(long r1, long r2);
        void multrow(long r, const bigint& scal);
        void divrow(long r, const bigint& scal);
        void clearrow(long r);
        mat_m& operator+=(const mat_m&);
        mat_m& operator-=(const mat_m&);
        mat_m& operator*=(const bigint&);
        mat_m& operator/=(const bigint&);
  // shortens to a matrix o ints or longs if possible
  //the parameter here is a dummy just to distinguish these
        mat_i  shorten(int) const;
        mat_l shorten(long) const;

     // non-member (friend) functions and operators
        friend long nrows(const mat_m&);      
        friend long ncols(const mat_m&);      
        friend mat_m operator*(const mat_m&, const mat_m&);
        friend vec_m operator*(const mat_m&, const vec_m&);
        friend int operator==(const mat_m&, const mat_m&);
        friend ostream& operator<< (ostream&s, const mat_m&);
        friend istream& operator>> (istream&s, mat_m&);
        friend mat_m colcat(const mat_m& a, const mat_m& b);
        friend mat_m rowcat(const mat_m& a, const mat_m& b);
        friend mat_m directsum(const mat_m& a, const mat_m& b);
        friend void elimrows(mat_m& m, long r1, long r2, long pos); //plain elimination, no clearing
        friend void elimrows1(mat_m& m, long r1, long r2, long pos); //elimination + clearing
        friend void elimrows2(mat_m& m, long r1, long r2, long pos, const bigint& last); //elimination + divide by last pivot
        friend mat_m echelon0(const mat_m& m, vec_i& pcols, vec_i& npcols,
                                  long& rk, long& ny, bigint& d);
        friend void elimp(const mat_m& m, long r1, long r2, long pos, const bigint& pr);
        friend mat_m echmodp(const mat_m& m, vec_i& pcols, vec_i& npcols,
                                  long& rk, long& ny, const bigint& pr);
        friend msubspace combine(const msubspace& s1, const msubspace& s2);
        friend mat_m restrict(const mat_m& m, const msubspace& s);
        friend msubspace lift(const msubspace& s, const bigint& pr, int =0);
        friend msubspace pcombine(const msubspace& s1, const msubspace& s2, const bigint& pr);
        friend mat_m prestrict(const mat_m& m, const msubspace& s, const bigint& pr);
        friend mat_m matmulmodp(const mat_m&, const mat_m&, const bigint& pr);

// Implementation
private:
       long nro,nco;
       bigint * entries;  // stored in one array, by rows
};



// Declaration of non-friend functions

mat_m operator+(const mat_m&);                   // unary
mat_m operator-(const mat_m&);                   // unary
mat_m operator+(const mat_m& m1, const mat_m& m2);
mat_m operator-(const mat_m& m1, const mat_m& m2);
mat_m operator*(const bigint& scal, const mat_m& m);
mat_m operator/(const mat_m& m, const bigint& scal);
int operator!=(const mat_m& m1, const mat_m& m2);
mat_m midmat(long n);  // = multi-idmat
mat_m transpose(const mat_m& m);
mat_m submatrix(const mat_m& m, const vec_i& iv, const vec_i& jv);
mat_m echelon(const mat_m& m, vec_i& pcols, vec_i& npcols,
                          long& rk, long& ny, bigint& d, int method=0);
mat_m echelon(const mat_m& m, vec_l& pcols, vec_l& npcols,
                          long& rk, long& ny, bigint& d, int method=0);
long rank(const mat_m&);
long nullity(const mat_m&);
bigint trace(const mat_m&);
vector<bigint> charpoly(const mat_m&);
bigint determinant(const mat_m&);
mat_m addscalar(const mat_m&, const bigint&);
vec_m apply(const mat_m&, const vec_m&);

#endif
