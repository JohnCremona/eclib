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
class smatT {

  friend class smatT_elim<T>;
  friend class ssubspaceT<T>;

protected:
  int nco;            // number of columns
  int nro;            // number of rows
  int **col;          // holds cols of entries
  T **val;       // holds values of entries

public:
  // constructors

  smatT<T> (int nr=0, int nc=0);
  smatT<T> (const smatT<T>&);                  // copy constructor
  explicit smatT<T> (const matT<T> &);         // conversion constructor
  ~smatT<T>();                             // destructor

  // member functions & operators

  matT<T> as_mat( ) const;
  smatT<T>& operator=(const smatT<T>&);       // assignment with copy
  T elem(int i, int j) const;   // returns value of (i,j) entry
  smatT<T>& operator+= (const smatT<T>&);
  smatT<T>& operator+= (const T&); // adds T*identity
  smatT<T>& operator-= (const smatT<T>&);
  smatT<T>& operator-= (const T& s)   // subtracts scalar*identity
  {this->operator+=(-s); return *this;}
  smatT<T>& operator*= (T);
  void sub_mod_p(const T& lambdal, const T& p);
  // subtracts scalar*identity mod p
  void reduce_mod_p(const T& p);
  smatT<T>& mult_by_scalar_mod_p (T scal, const T& p);
  smatT<T>& operator/= (T);
  matT<T> operator*( const matT<T>& );
  void set_row ( int, int, int*, T* );
  smatT<T> select_rows(const vecT<int>& rows) const;
  void setrow ( int i, const svecT<T>& v); // i counts from 1
  void setrow ( int i, const vecT<T>& v); // i counts from 1
  svecT<T> row(int) const; // extract row i as an svec
  int nrows() const {return nro;}
  int ncols() const {return nco;}
  int rank(T mod); // implemented in smat_elim.cc
  int nullity(const T& lambda, T mod); // nullity of this-lambda*I

  static smatT<T> scalar_matrix(int n, const T& a);  // nxn scalar matrix a*I
  static smatT<T> identity_matrix(int n) {return scalar_matrix(n, T(1));}  // nxn identity matrix I

  // non-member (friend) functions and operators

  friend vector<int> dim<>(const smatT<T>& A);
  template<class T1> friend vecT<T1> operator*(const smatT<T1>& m, const vecT<T1>& v);
  template<class T1> friend svecT<T1> operator* ( const smatT<T1>& A, const svecT<T1>& v );
  template<class T1> friend svecT<T1> operator* ( const svecT<T1>& v, const smatT<T1>& A );
  friend svecT<T> mult_mod_p<>( const smatT<T>& A, const svecT<T>& v, const T& p  );
  friend vecT<T> mult_mod_p<>( const smatT<T>& A, const vecT<T>& v, const T& p  );
  friend svecT<T> mult_mod_p<>( const svecT<T>& v, const smatT<T>& A, const T& p  );
  template<class T1>  friend smatT<T1> operator* ( const smatT<T1>& A, const smatT<T1>& B );
  friend smatT<T> mult_mod_p<> ( const smatT<T>& A, const smatT<T>& B, const T& p );
  friend smatT<T> mult_mod_p_flint<> ( const smatT<T>& A, const smatT<T>& B, const T& p );
  friend T maxabs<>( const smatT<T>& A);
  friend smatT<T> transpose<>(const smatT<T>&);
  friend int operator==<>(const smatT<T>&, const smatT<T>&);
  friend int eqmodp<>(const smatT<T>&, const smatT<T>&, const T& p);
  friend ostream& operator<<<> (ostream&s, const smatT<T>&);
  friend istream& operator>><> (istream&s, smatT<T>&);
  friend int get_population<> (const smatT<T>& );
  friend double density<> (const smatT<T>& m);
  friend void random_fill_in<>( smatT<T>&, int, int );
  friend int liftmat<>(const smatT<T>& mm, T pr, smatT<T>& m, T& dd);
  friend int liftmats_chinese<>(const smatT<T>& mm1, T pr1, const smatT<T>& mm2, T pr2,
                                smatT<T>& m, T& dd);
};

// Declaration of non-friend functions

template<class T>
smatT<T> operator+(const smatT<T>&);                   // unary
template<class T>
smatT<T> operator-(const smatT<T>&);                   // unary
template<class T>
smatT<T> operator+(const smatT<T>& m1, const smatT<T>& m2);
template<class T>
smatT<T> operator-(const smatT<T>& m1, const smatT<T>& m2);
template<class T>
smatT<T> operator*(T scal, const smatT<T>& m);
template<class T>
smatT<T> operator/(const smatT<T>& m, T scal);
template<class T>
int operator!=(const smatT<T>& m1, const smatT<T>& m2);
template<class T>
inline void display_population(const smatT<T>& A)
{cout << " number of non-zero entries: " << get_population(A) << endl;}

#endif
