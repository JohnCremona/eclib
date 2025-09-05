// smatrix.h: manage declarations for sparse integer matrix classes
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
 
// Original version by Luiz Figueiredo
 
#if     !defined(_ECLIB_SMATRIX_H)
#define _ECLIB_SMATRIX_H      1       //flags that this file has been included

#include "matrix.h"
#include "svector.h"

template<class T>
class sZmat {

  friend class sZmat_elim<T>;
  friend class ssubZspace<T>;

protected:
  int nco;            // number of columns
  int nro;            // number of rows
  int **col;          // holds cols of entries
  T **val;       // holds values of entries

public:
  // constructors

  sZmat<T> (int nr=0, int nc=0);
  sZmat<T> (const sZmat<T>&);                  // copy constructor
  explicit sZmat<T> (const Zmat<T> &);         // conversion constructor
  ~sZmat<T>();                             // destructor

  // member functions & operators

  Zmat<T> as_mat( ) const;
  sZmat<T>& operator=(const sZmat<T>&);       // assignment with copy
  T elem(int i, int j) const;   // returns value of (i,j) entry
  sZmat<T>& operator+= (const sZmat<T>&);
  sZmat<T>& operator+= (const T&); // adds T*identity
  sZmat<T>& operator-= (const sZmat<T>&);
  sZmat<T>& operator-= (const T& s)   // subtracts scalar*identity
  {this->operator+=(-s); return *this;}
  sZmat<T>& operator*= (T);
  void sub_mod_p(const T& lambdal, const T& p);
  // subtracts scalar*identity mod p
  void reduce_mod_p(const T& p);
  sZmat<T>& mult_by_scalar_mod_p (T scal, const T& p);
  sZmat<T>& operator/= (T);
  Zmat<T> operator*( const Zmat<T>& );
  void set_row ( int, int, int*, T* );
  sZmat<T> select_rows(const Zvec<int>& rows) const;
  void setrow ( int i, const sZvec<T>& v); // i counts from 1
  void setrow ( int i, const Zvec<T>& v); // i counts from 1
  sZvec<T> row(int) const; // extract row i as an svec
  int nrows() const {return nro;}
  int ncols() const {return nco;}
  int rank(T mod); // implemented in smat_elim.cc
  int nullity(const T& lambda, T mod); // nullity of this-lambda*I

  static sZmat<T> scalar_matrix(int n, const T& a);  // nxn scalar matrix a*I
  static sZmat<T> identity_matrix(int n) {return scalar_matrix(n, T(1));}  // nxn identity matrix I

  // non-member (friend) functions and operators

  friend vector<int> dim<>(const sZmat<T>& A);
  template<class T1> friend Zvec<T1> operator*(const sZmat<T1>& m, const Zvec<T1>& v);
  template<class T1> friend sZvec<T1> operator* ( const sZmat<T1>& A, const sZvec<T1>& v );
  template<class T1> friend sZvec<T1> operator* ( const sZvec<T1>& v, const sZmat<T1>& A );
  friend sZvec<T> mult_mod_p<>( const sZmat<T>& A, const sZvec<T>& v, const T& p  );
  friend Zvec<T> mult_mod_p<>( const sZmat<T>& A, const Zvec<T>& v, const T& p  );
  friend sZvec<T> mult_mod_p<>( const sZvec<T>& v, const sZmat<T>& A, const T& p  );
  template<class T1>  friend sZmat<T1> operator* ( const sZmat<T1>& A, const sZmat<T1>& B );
  friend sZmat<T> mult_mod_p<> ( const sZmat<T>& A, const sZmat<T>& B, const T& p );
  friend T maxabs<>( const sZmat<T>& A);
  friend sZmat<T> transpose<>(const sZmat<T>&);
  friend int operator==<>(const sZmat<T>&, const sZmat<T>&);
  friend int eqmodp<>(const sZmat<T>&, const sZmat<T>&, const T& p);
  friend ostream& operator<<<> (ostream&s, const sZmat<T>&);
  friend istream& operator>><> (istream&s, sZmat<T>&);
  friend int get_population<> (const sZmat<T>& );
  friend double density<> (const sZmat<T>& m);
  friend void random_fill_in<>( sZmat<T>&, int, int );
  friend int liftmat<>(const sZmat<T>& mm, T pr, sZmat<T>& m, T& dd);
  friend int liftmats_chinese<>(const sZmat<T>& mm1, T pr1, const sZmat<T>& mm2, T pr2,
                                sZmat<T>& m, T& dd);
};

// Declaration of non-friend functions

template<class T>
sZmat<T> operator+(const sZmat<T>&);                   // unary
template<class T>
sZmat<T> operator-(const sZmat<T>&);                   // unary
template<class T>
sZmat<T> operator+(const sZmat<T>& m1, const sZmat<T>& m2);
template<class T>
sZmat<T> operator-(const sZmat<T>& m1, const sZmat<T>& m2);
template<class T>
sZmat<T> operator*(T scal, const sZmat<T>& m);
template<class T>
sZmat<T> operator/(const sZmat<T>& m, T scal);
template<class T>
int operator!=(const sZmat<T>& m1, const sZmat<T>& m2);
template<class T>
inline void display_population(const sZmat<T>& A)
{cout << " number of non-zero entries: " << get_population(A) << endl;}

#endif
