// smat.h: declarations for sparse integer matrix class smat
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
 
// Original version by Luiz Figueiredo

class smat {

friend class smat_elim;

protected:
         int nco;            // number of columns
         int nro;            // number of rows
	 int **col;          // holds cols of entries
	 scalar **val;          // holds values of entries

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
     void sub_mod_p(const scalar& lambdal, const scalar& p=BIGPRIME); 
      // subtracts scalar*identity mod p
     void reduce_mod_p(const scalar& p=BIGPRIME);
     smat& mult_by_scalar_mod_p (scalar scal, const scalar& p=BIGPRIME);
     smat& operator/= (scalar);
     mat operator*( const mat& );
     void set_row ( int, int, int*, scalar* );
     smat select_rows(const vec& rows) const;
     void setrow ( int i, const svec& v); // i counts from 1
     svec row(int) const; // extract row i as an svec
     
     // non-member (friend) functions and operators

     friend inline int nrows(const smat& A) {return A.nro;}
     friend inline int ncols(const smat& A) {return A.nco;}
     friend inline vector<int> dim(const smat& A) 
     {vector<int>d; d.push_back(A.nro);d.push_back(A.nco);return d;}
     friend vec operator*  (smat& m, const vec& v);
     friend svec operator* ( const smat& A, const svec& v );
     friend svec operator* ( const svec& v, const smat& A );
     friend svec mult_mod_p( const smat& A, const svec& v, const scalar& p  );
     friend svec mult_mod_p( const svec& v, const smat& A, const scalar& p  );
     friend smat operator* ( const smat& A, const smat& B );
     friend smat mult_mod_p ( const smat& A, const smat& B, const scalar& p );
     friend smat transpose(const smat&);
     friend int operator==(const smat&, const smat&);
  // Equality mod p:
     friend int eqmodp(const smat&, const smat&, const scalar& p=BIGPRIME);
     friend ostream& operator<< (ostream&s, const smat&);
     friend istream& operator>> (istream&s, smat&);
     friend int get_population (const smat& );      //mainly used for testing
     friend inline double density (const smat& m)
     {return (((double)(get_population(m)))/m.nro)/m.nco;}
     friend void random_fill_in( smat&, int, scalar ); //the elimination program
     friend smat sidmat(scalar);  // identity matrix
     friend smat liftmat(const smat& mm, scalar pr, scalar& dd, int trace=0);
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
