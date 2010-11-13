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
 
class smat {

friend class smat_elim;

protected:
  int nco;            // number of columns
  int nro;            // number of rows
  vector<svec> rows;  // holds rows as svecs (length = nro+1, we count from 1)
  
public:
  // constructors
  
  smat (int nr=0, int nc = 0)
    :nco(nc), nro(nr) 
  {
    rows.resize(nro+1); 
    for(int i=1; i<=nro; i++) {rows[i]=svec(nco);}
  }
  smat (const smat & m)                  // copy constructor
    :nco(m.nco), nro(m.nro), rows(m.rows) {;}
  smat (const mat &);                  // conversion constructor
  
  // member functions & operators
  
  smat& operator= (const smat&);  // assignment with copy
  mat as_mat( ) const;            // conversion to ordinary matrix
  
  void setrow ( int i, const svec& v) {rows[i]=v;}
  svec row(int i) {return rows[i];}
  smat select_rows(const vec& v) const;
  scalar elem(int i, int j) const  {return rows[i].elem(j); }
  
  smat& operator+= (const smat&);
  smat& operator+= (const scalar&); // adds scalar*identity
  smat& operator-= (const smat&);
  smat& operator-= (const scalar& s);   // subtracts scalar*identity
  smat& operator*= (scalar);
  void sub_mod_p(const scalar& lambda); // subtracts scalar*identity
					// mod p
  void reduce_mod_p(const scalar& p=BIGPRIME);
  smat& mult_by_scalar_mod_p (scalar scal, const scalar& p=BIGPRIME);
  smat& operator/= (scalar);
  mat operator*( const mat& );
     
  // non-member (friend) functions and operators
  
  friend inline int nrows(const smat& A) {return A.nro;}
  friend inline int ncols(const smat& A) {return A.nco;}
  friend inline vector<int> dim(const smat& A) 
  {vector<int>d; d.push_back(A.nro);d.push_back(A.nco);return d;}
  friend vector<std::set<int> > row_supports(const smat& A);
  friend vec operator*  (smat& m, const vec& v);
  friend svec operator* ( const smat& A, const svec& v );
  friend svec operator* ( const svec& v, const smat& A );
  friend svec mult_mod_p( const smat& A, const svec& v, const scalar& p  );
  friend svec mult_mod_p( const svec& v, const smat& A, const scalar& p  );
  friend smat operator* ( const smat& A, const smat& B );
  friend smat mult_mod_p ( const smat& A, const smat& B, const scalar& p );
  friend smat transpose(const smat&);
  friend int operator==(const smat&, const smat&);
  friend int operator==(const smat&, const mat&);
  friend int operator==(const mat& m, const smat& sm) {return sm==m;}
  // Equality mod p:
  friend int eqmodp(const smat&, const smat&, const scalar& p=BIGPRIME);
  friend ostream& operator<< (ostream&s, const smat&);
  //     friend istream& operator>> (istream&s, smat&); // not implemented
  friend int get_population (const smat& );          // used for testing
  friend inline double density (const smat& m)
  {return (((double)(get_population(m)))/m.nro)/m.nco;}
  friend void random_fill_in( smat&, int, scalar& ); // elimination
  friend void random_fill_in( svec&, int, scalar& ); 
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
