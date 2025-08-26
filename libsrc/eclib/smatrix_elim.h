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

#include "ssubspace.h"
#include "subspace.h"

class smat_elim_list;
class smat_elim_ordlist;

template<class T>
class sZmat_elim : public sZmat<T>{

private:

  T modulus;
  int rank;

  smat_elim_ordlist* column; // an array of lists, oner per col, of row nos
		   // which have nonzero entries in that col
  int *position;
  int *elim_col;     // in which row was col eliminated;
  int *elim_row;       //order of elimination of rows;
  void clear_col(int,int, smat_elim_list&, int fr = 0, int fc = 0,int M = 0,int *li =0);
  void check_col( int col, smat_elim_list& L );
  void check_row (int d2, int row2, smat_elim_list& L );
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
  sZmat<T> old_kernel( Zvec<int>&, Zvec<int>& );
  sZmat<T> new_kernel( Zvec<int>&, Zvec<int>& );
  sZmat<T> kernel( Zvec<int>&, Zvec<int>& );
  void normalize( int, int );
  void eliminate( const int&, const int& );
  void step5dense();
  void free_space( int col );
  void elim( int row1, int row2, T v2 );
  // constructor:
  explicit sZmat_elim( const sZmat<T>& sm, T mod) : sZmat<T>( sm ), modulus(mod) { init(); };
  sZmat_elim( int r = 0,int c = 0 );
  // destructor:
  ~sZmat_elim();

};

class smat_elim_list {
public:
  static int listsize;
  int maxsize;
  int *list_array;
  int num;
  int index;
  void put( int X)
  {
    if( num >= maxsize ) grow();
    list_array[ num++ ] = X;
  }
  int find( int X, int ub, int lb = 0 );
  void grow ();
  int next()
  {
    if( index < num ) return( list_array[index++] ); else return -1;
  }
  smat_elim_list( int m = 10);
  ~smat_elim_list( );
  void clear( int m = 0);
};

class smat_elim_ordlist : public smat_elim_list {
public:
  void put( int X);
  void put( smat_elim_list& L);     // L must be ordered
  void remove( int X );
  void remove( smat_elim_list& L );     // L must be ordered
  smat_elim_ordlist( int m = 10) : smat_elim_list(m) {;}
};

ostream& operator<< (ostream&s, const smat_elim_list&);

inline ostream& operator<< (ostream&s, const smat_elim_list& L)
{
  s<<"[";
  int n=L.num;
  int* x=L.list_array;
  while(n--) s<<(*x++)<<" ";
  s<<"]";
  return s;
}

inline int
find( int X, const int* ptr, int ub, int lb = 0 ) {
  if( ptr[ub] < X ) return ub;
  while( ptr[lb] < X ) {
    int i = (ub + lb)/2;
    ptr[i] < X ? (lb = i+1) : (ub = i);
  }
  return lb;
}

#endif
