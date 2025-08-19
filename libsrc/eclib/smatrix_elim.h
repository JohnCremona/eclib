// smatrix_elim.h: manages declarations of sparse integer matrix classes
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
 
#ifndef _ECLIB_SMATRIX_ELIM_H
#define _ECLIB_SMATRIX_ELIM_H 1       //flags that this file has been included

#include "smatrix.h"
#include "subspace.h"

inline int
find( int X, const int* ptr, int ub, int lb = 0 ) {
  if( ptr[ub] < X ) return ub;
  while( ptr[lb] < X ) {
    int i = (ub + lb)/2;
    ptr[i] < X ? (lb = i+1) : (ub = i);
  }
  return lb;
}

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim
#undef ssubspace

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i
#define svec svec_i
#define smat smat_i
#define smat_elim smat_i_elim
#define ssubspace ssubspace_i

#include "smat_elim.h"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim
#undef ssubspace

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l
#define svec svec_l
#define smat smat_l
#define smat_elim smat_l_elim
#define ssubspace ssubspace_l

#include "smat_elim.h"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim
#undef ssubspace

#define scalar bigint
#define vec vec_m
#define mat mat_m
#define subspace subspace_m
#define svec svec_m
#define smat smat_m
#define smat_elim smat_m_elim
#define ssubspace ssubspace_m

#include "smat_elim.h"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim
#undef ssubspace

///////////////////////////////////////////////////////////////////////////

template<class T> class vecT;
template<class T> class svecT;
template<class T> class smatT;
template<class T> class smatT_elim;
template<class T> class matT;
template<class T> class subspaceT;
template<class T> class ssubspaceT;
template<class T> class smat_elim_list;
template<class T> class smat_elim_ordlist;

template<class T>
class smat_elim_list {
public:
  static int listsize;
  int maxsize;
  T *list_array;
  int num;
  int index;
  void put( const T& X)
  {
    if( num >= maxsize ) grow();
    list_array[ num++ ] = X;
  }
  int find( const T& X, int ub, int lb = 0 );
  void grow ();
  T next()
  {
    if( index < num ) return( list_array[index++] ); else return -1;
  }
  smat_elim_list( int m = 10);
  ~smat_elim_list( );
  void clear( int m = 0);
};

template<class T>
class smat_elim_ordlist : public smat_elim_list<T> {
public:
  void put( T& X);
  void put( smat_elim_list<T>& L);     // L must be ordered
  void remove( T& X );
  void remove( smat_elim_list<T>& L );     // L must be ordered
  smat_elim_ordlist( int m = 10) : smat_elim_list<T>(m) {;}
};

template<class T> ostream& operator<< (ostream&s, const smat_elim_list<T>&);

template<class T>
class smatT_elim : public smatT<T>{

private:

  T modulus;
  int rank;

  smat_elim_ordlist<T>* column; // an array of lists, oner per col, of row nos
		   // which have nonzero entries in that col
  int *position;
  int *elim_col;     // in which row was col eliminated;
  int *elim_row;       //order of elimination of rows;
  void clear_col(int,int, smat_elim_list<T>&, int fr = 0, int fc = 0,int M = 0,int *li =0);
  void check_col( int col, smat_elim_list<T>& L );
  void check_row (int d2, int row2, smat_elim_list<T>& L );
  int get_weight( int, const int* );
  int has_weight_one( int, const int* );
  int n_active_cols(); // number of active columns
  int n_active_rows(); // number of active rows
  long n_active_entries(); // number of active entries
  double active_density(); // density of non-eliminated part
  void report(); // report current state (density etc)

public:
  int get_rank( ) { return rank; }
  void init( );
  void step0();
  void step1();
  void step2();
  void step3();
  void step4();
  void standard( );
  void back_sub( );
  void sparse_elimination( );
  smatT<T> old_kernel( vecT<int>&, vecT<int>& );
  smatT<T> new_kernel( vecT<int>&, vecT<int>& );
  smatT<T> kernel( vecT<int>&, vecT<int>& );
  void normalize( int, int );
  void eliminate( const int&, const int& );
  void step5dense();
  void free_space( int col );
  void elim( int row1, int row2, T v2 );
  // constructor:
  explicit smatT_elim( const smatT<T>& sm, T mod) : smatT<T>( sm ), modulus(mod) { init(); };
  smatT_elim( int r = 0,int c = 0 );
  // destructor:
  ~smatT_elim();

  friend ostream& operator<<<> (ostream&s, const smat_elim_list<T>&);

};

template<class T>
inline ostream& operator<< (ostream&s, const smat_elim_list<T>& L)
{
  s<<"[";
  int n=L.num;
  int* x=L.list_array;
  while(n--) s<<(*x++)<<" ";
  s<<"]";
  return s;
}

template<class T> int dim(const ssubspaceT<T>& s)     {return s.basis.ncols();}
template<class T> vecT<int> pivots(const ssubspaceT<T>& s)  {return s.pivots;}
template<class T> smatT<T> basis(const ssubspaceT<T>& s)  {return s.basis;}
template<class T> ssubspaceT<T> combine(const ssubspaceT<T>& s1, const ssubspaceT<T>& s2);
template<class T> smatT<T> restrict_mat(const smatT<T>& m, const ssubspaceT<T>& s);

template<class T>
class ssubspaceT {

public:
  // constructors
  ssubspaceT<T>(int n=0);
  ssubspaceT<T>(const smatT<T>& b, const vecT<int>& p, T mod);
  ssubspaceT<T>(const ssubspaceT<T>& s);
  // assignment
  void operator=(const ssubspaceT<T>& s);

  // member functions & operators
  inline void clear() { pivots.init(); basis=smatT<T>(0,0);}
  inline vecT<int> pivs() const {return pivots;}  // the pivot vector
  inline smatT<T> bas() const {return basis;}   // the basis matrix
  inline T mod() const {return modulus;}   // the (prime) modulus

  // non-member (friend) functions and operators
  friend int dim(const ssubspaceT<T>& s)     {return s.basis.ncols();}
  friend vecT<int> pivots(const ssubspaceT<T>& s)  {return s.pivots;}
  friend smatT<T> basis(const ssubspaceT<T>& s)  {return s.basis;}
  friend ssubspaceT<T> combine<>(const ssubspaceT<T>& s1, const ssubspaceT<T>& s2);
  friend smatT<T> restrict_mat<>(const smatT<T>& m, const ssubspaceT<T>& s);

  // Implementation
private:
  T modulus;
  vecT<int> pivots;
  smatT<T> basis;
};


// Declarations of nonmember, nonfriend operators and functions:

template<class T>
ssubspaceT<T> kernel(const smatT<T>& sm, T m);
template<class T>
ssubspaceT<T> eigenspace(const smatT<T>& sm, T lambda, T m);
template<class T>
ssubspaceT<T> subeigenspace(const smatT<T>& sm, T l, const ssubspaceT<T>& s, T m);

// construction of a 1-dimensional sparse subspace from a vector:
template<class T>
ssubspaceT<T> make1d(const vecT<T>& bas, int piv, T m);


#endif
