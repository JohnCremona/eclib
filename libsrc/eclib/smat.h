// smat.h: declarations for sparse integer matrix class smat
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
 
// Should not be included directly by user: include smatrix.h instead

#include "eclib/flinterface.h"

// We cannot use one macro to prevent this being included twice since
// we want it to be included twice by smatrix.h which defined _ECLIB_SMATRIX_H

// Original version by Luiz Figueiredo

int eqmodp(const smat&, const smat&, const scalar& p=DEFAULT_MODULUS);
int liftmat(const smat& mm, scalar pr, smat& m, scalar& dd, int trace=0);

class smat {

friend class smat_elim;

protected:
         int nco;            // number of columns
         int nro;            // number of rows
	 int **col;          // holds cols of entries
	 scalar **val;       // holds values of entries

public:
     // constructors

     smat (int nr=0, int nc = 0);
     smat (const smat&);                  // copy constructor
     smat (const mat &);                  // conversion constructor
     ~smat();                             // destructor

     // member functions & operators

     mat as_mat( ) const;
     smat& operator=(const smat&);       // assignment with copy
     scalar  elem(int i, int j) const;   // returns value of (i,j) entry
     smat& operator+= (const smat&);
     smat& operator+= (const scalar&); // adds scalar*identity
     smat& operator-= (const smat&);
     smat& operator-= (const scalar& s)   // subtracts scalar*identity
      {this->operator+=(-s); return *this;}
     smat& operator*= (scalar);
     void sub_mod_p(const scalar& lambdal, const scalar& p=DEFAULT_MODULUS); 
      // subtracts scalar*identity mod p
     void reduce_mod_p(const scalar& p=DEFAULT_MODULUS);
     smat& mult_by_scalar_mod_p (scalar scal, const scalar& p=DEFAULT_MODULUS);
     smat& operator/= (scalar);
     mat operator*( const mat& );
     void set_row ( int, int, int*, scalar* );
     smat select_rows(const vec& rows) const;
     void setrow ( int i, const svec& v); // i counts from 1
     void setrow ( int i, const vec& v); // i counts from 1
     svec row(int) const; // extract row i as an svec
     int nrows() const {return nro;}
     int ncols() const {return nco;}
     long rank(scalar mod=DEFAULT_MODULUS);
     long nullity(const scalar& lambda, scalar mod=DEFAULT_MODULUS); // nullity of this-lambda*I

     // non-member (friend) functions and operators

     friend inline vector<int> dim(const smat& A) 
     {vector<int>d; d.push_back(A.nro);d.push_back(A.nco);return d;}
     friend vec operator*  (smat& m, const vec& v);
     friend svec operator* ( const smat& A, const svec& v );
     friend svec operator* ( const svec& v, const smat& A );
     friend svec mult_mod_p( const smat& A, const svec& v, const scalar& p  );
     friend vec mult_mod_p( const smat& A, const vec& v, const scalar& p  );
     friend svec mult_mod_p( const svec& v, const smat& A, const scalar& p  );
     friend smat operator* ( const smat& A, const smat& B );
     friend smat mult_mod_p ( const smat& A, const smat& B, const scalar& p );
     friend smat mult_mod_p_flint ( const smat& A, const smat& B, const scalar& p );
     friend smat transpose(const smat&);
     friend int operator==(const smat&, const smat&);
  // Equality mod p:
     friend int eqmodp(const smat&, const smat&, const scalar& p);
     friend ostream& operator<< (ostream&s, const smat&);
     friend istream& operator>> (istream&s, smat&);
     friend int get_population (const smat& );      //mainly used for testing
     friend inline double density (const smat& m)
     {return (((double)(get_population(m)))/m.nro)/m.nco;}
     friend void random_fill_in( smat&, int, scalar ); //the elimination program
     friend smat sidmat(scalar);  // identity matrix
     friend int liftmat(const smat& mm, scalar pr, smat& m, scalar& dd, int trace);
     friend int liftmats_chinese(const smat& mm1, scalar pr1, const smat& mm2, scalar pr2,
                                 smat& m, scalar& dd);
 };

// Declaration of non-friend functions

smat operator+(const smat&);                   // unary
smat operator-(const smat&);                   // unary
smat operator+(const smat& m1, const smat& m2);
smat operator-(const smat& m1, const smat& m2);
smat operator*(scalar scal, const smat& m);
smat operator/(const smat& m, scalar scal);
int operator!=(const smat& m1, const smat& m2);
smat sidmat(scalar);  // identity matrix

inline void display_population(const smat& A)
{
  cout << " number of non-zero entries: " << get_population(A) << endl;
}
