// mat.h: declarations for integer matrix classes
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
 
// Not to be included directly by user: use matrix.h
//
// SCALAR_OPTION must be set to 1 or 2 by including file

#ifndef LONG_MIN
#define LONG_MIN (-LONG_MAX-1)
#endif
#ifndef INT_MIN
#define INT_MIN (-INT_MAX-1)
#endif

class mat {
friend class subspace;
friend class mat_m;
friend class smat;
friend class svec;
friend class smat_elim;
public:
     // constructors
        mat(long nr=0, long nc=0);
        mat(const mat&);                    // copy constructor
     // destructor
        ~mat();

     // member functions & operators
        void init(long nr=0, long nc=0);
        mat& operator=(const mat&);          // assignment with copy
        scalar& operator()(long i, long j) const;   // returns ref to (i,j) entry
        mat slice(long,long,long=-1,long=-1) const;// returns submatrix
        scalar sub(long i, long j) const;             // returns the (i,j) entry
        vec row(long i) const;                // returns row i (as a vec)
        vec col(long j) const;                // returns col j (as a vec)
        void set(long i, long j, scalar x);     // sets the (i,j) entry to x
        void add(long i, long j, scalar x);  // adds x to the (i,j) entry  
        void setrow(long i, const vec& v);
        void setcol(long i, const vec& v);
        void swaprows(long r1, long r2);
        void multrow(long r, scalar scal);
        void divrow(long r, scalar scal);
        void clearrow(long r);
        mat& operator+=(const mat&);
        mat& operator-=(const mat&);
        mat& operator*=(scalar);
        mat& operator/=(scalar);
	const scalar* get_entries()const{return entries;}
	void output(ostream&s=cout) const;
	void output_pari(ostream&s=cout)   const;
        void output_pretty(ostream&s=cout)   const;
        void dump_to_file(char* filename) const; // binary output
        void read_from_file(char* filename);     // binary input

     // non-member (friend) functions and operators
        friend long nrows(const mat&);      
        friend long ncols(const mat&);      
        friend void add_row_to_vec(vec& v, const mat& m, long i);
        friend void sub_row_to_vec(vec& v, const mat& m, long i);
        friend mat operator*(const mat&, const mat&);
	friend vec operator*(const mat&, const vec&);
        friend int operator==(const mat&, const mat&);
        friend istream& operator>> (istream&s, mat&);
        friend mat colcat(const mat& a, const mat& b);
        friend mat rowcat(const mat& a, const mat& b);
        friend mat directsum(const mat& a, const mat& b);
        friend void elimrows(mat& m, long r1, long r2, long pos); //plain elimination, no clearing
        friend void elimrows1(mat& m, long r1, long r2, long pos); //elimination + clearing
        friend void elimrows2(mat& m, long r1, long r2, long pos, scalar last); //elimination + divide by last pivot
	friend mat echelon0(const mat& m, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar& d);
	friend mat echelonl(const mat& m, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar& d);
	friend void elimp(const mat& m, long r1, long r2, long pos, scalar pr);
	friend void elimp1(const mat& m, long r1, long r2, long pos, scalar pr);
	friend mat echelonp(const mat& m, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar& d, scalar pr);
	friend mat echmodp(const mat& m, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar pr);
	friend mat echmodp_uptri(const mat& m, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar pr);
	friend subspace combine(const subspace& s1, const subspace& s2);
	friend mat restrict(const mat& m, const subspace& s, int cr);
	friend mat liftmat(const mat& mm, scalar pr, scalar& dd, int trace=0);
	friend subspace lift(const subspace& s, scalar pr, int trace=0);
	friend subspace pcombine(const subspace& s1, const subspace& s2, scalar pr);
	friend mat prestrict(const mat& m, const subspace& s, scalar pr, int cr);
	friend mat matmulmodp(const mat&, const mat&, scalar pr);
	friend mat echmodp_d(const mat& mat, vec& pcols, vec& npcols, long& rk, long& ny, double pr);
        friend double sparsity(const mat& m);
// Implementation
private:
       long nro,nco;
       scalar * entries;  // stored in one array, by rows
};

// Declaration of non-friend functions

inline ostream& operator<< (ostream&s, const mat&m) 
{m.output(s); return s;}

mat operator+(const mat&);                   // unary
mat operator-(const mat&);                   // unary
mat operator+(const mat& m1, const mat& m2);
mat operator-(const mat& m1, const mat& m2);
mat operator*(scalar scal, const mat& m);
mat operator/(const mat& m, scalar scal);
int operator!=(const mat& m1, const mat& m2);
mat idmat(scalar n);
mat transpose(const mat& m);
mat submat(const mat& m, const vec& iv, const vec& jv);
mat echelon(const mat& m, vec& pcols, vec& npcols,
                          long& rk, long& ny, scalar& d, int method=0);  // default method 0: scalars
long rank(const mat&);
long nullity(const mat&);
long trace(const mat&);
vector<long> charpoly(const mat&);
long determinant(const mat&);
mat addscalar(const mat&, scalar);
vec apply(const mat&, const vec&);
