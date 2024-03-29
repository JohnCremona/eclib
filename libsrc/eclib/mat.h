// mat.h: declarations for integer matrix classes
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
 
// Not to be included directly by user: use matrix.h which defines
// _ECLIB_MATRIX_H and includes this twice
//

#ifndef LONG_MIN
#define LONG_MIN (-LONG_MAX-1)
#endif
#ifndef INT_MIN
#define INT_MIN (-INT_MAX-1)
#endif

int liftmat(const mat& mm, scalar pr, mat& m, scalar& dd, int trace=0);
int lift(const subspace& s, scalar pr, subspace& ans, int trace=0);

class mat {
friend class subspace;
friend class mat_m;
friend class smat;
friend class svec;
friend class vec;
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
        long nrows() const {return nro;}
        long ncols() const {return nco;}
        long rank() const;
        long nullity() const;
        long trace() const;
        vector<long> charpoly() const;
        long determinant() const;
	void output(ostream&s=cout) const;
	void output_pari(ostream&s=cout)   const;
        void output_pretty(ostream&s=cout)   const;
        void dump_to_file(string filename) const; // binary output
        void read_from_file(string filename);     // binary input

     // non-member (friend) functions and operators
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
        friend mat ref_via_flint(const mat& M, vec& pcols, vec& npcols,
                                 long& rk, long& ny, scalar pr);
        friend mat ref_via_ntl(const mat& M, vec& pcols, vec& npcols,
                                 long& rk, long& ny, scalar pr);
        friend mat rref(const mat& M, vec& pcols, vec& npcols,
                                 long& rk, long& ny, scalar pr)
  {
#if FLINT
    return ref_via_flint(M, pcols, npcols, rk, ny, pr);
#else
    return ref_via_ntl(M, pcols, npcols, rk, ny, pr);
#endif
  }
        friend long rank_via_ntl(const mat& M, scalar pr);
        friend long det_via_ntl(const mat& M, scalar pr);
	friend subspace combine(const subspace& s1, const subspace& s2);
        friend mat restrict_mat(const mat& m, const subspace& s, int cr);
        friend int liftmat(const mat& mm, scalar pr, mat& m, scalar& dd, int trace);
        friend int lift(const subspace& s, scalar pr, subspace& ans, int trace);
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
mat addscalar(const mat&, scalar);
vec apply(const mat&, const vec&);

// Construct an NTL mat_lzz_p (matrix mod p) from a mat mod pr
mat_zz_p mat_zz_p_from_mat(const mat& M, scalar pr);

// Construct a mat (scalar type same as pr) from an NTL mat_lzz_p

mat mat_from_mat_zz_p(const mat_zz_p& A, scalar pr); // type of scalar fixes return type
