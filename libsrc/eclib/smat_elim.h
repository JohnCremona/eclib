// smat_elim.h: declaration of class smat_elim, for sparse elimination
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
 
// Not to be included directly by user: use smatrix_elim.h which
// defines _ECLIB_SMATRIX_ELIM_H and includes this twice

class smat_elim : public smat{

 private:
  
  class list {
  public:
    typedef int type;    //use this till smat is made a template class.
    static int listsize;
    int maxsize;
    type *list_array;
    int num;
    int index;
    void put( const type& X) { 
      if( num >= maxsize ) 
	{
// 	  cout<<"About to grow list from size "<<maxsize<<endl;
	  grow();
// 	  cout<<"Finished growing list, new size "<<maxsize<<endl;
	}
      list_array[ num++ ] = X; 
    }
    int find( const type& X, int ub, int lb = 0 );
    void grow ();
    type next() { 
      if( index < num ) return( list_array[index++] ); else return(-1); 
    }
    list( int m = 10);
    ~list( );
    void clear( int m=0); 
    
  };
  friend ostream& operator<< (ostream&s, const smat_elim::list&);

  class ordlist : public list {
  public:
    void put( type& X);
    void put( list& L);     // L must be ordered
    void remove( type& X );
    void remove( list& L );     // L must be ordered
    ordlist( int m = 10) : list(m) {;}
  };

  scalar modulus;
  int rank;

  ordlist* column; // an array of lists, oner per col, of row nos
		   // which have nonzero entries in that col
  int *position;
  int *elim_col;     // in which row was col eliminated;
  int *elim_row;       //order of elimination of rows;
  void clear_col(int,int,list&, int fr = 0, int fc = 0,int M = 0,int *li =0);
  void check_col( int col, list& L );
  void check_row (int d2, int row2, list& L ); 
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
  smat old_kernel( vec_i&, vec_i& );
  smat new_kernel( vec_i&, vec_i& );
  smat kernel( vec_i&, vec_i& );
  void normalize( int, int );
  void eliminate( const int&, const int& );
  void step5dense();
  void free_space( int col );
  void elim( int row1, int row2, scalar v2 );
  // constructor:
  explicit smat_elim( const smat& sm, scalar mod) : smat( sm ), modulus(mod) { init(); };
  smat_elim( int r = 0,int c = 0 );
  // destructor:
  ~smat_elim();
};

inline ostream& operator<< (ostream&s, const smat_elim::list& L)
{
  s<<"[";
  int n=L.num;
  int* x=L.list_array;
  while(n--) s<<(*x++)<<" ";
  s<<"]";
  return s;
}

class ssubspace {

public:
     // constructors
        ssubspace(int n=0);
        ssubspace(const smat& b, const vec_i& p, scalar mod);
	ssubspace(const ssubspace& s);
     // assignment
	void operator=(const ssubspace& s);

     // member functions & operators
        inline void clear() { pivots.init(); basis=smat(0,0);}
        inline vec_i pivs() const {return pivots;}  // the pivot vector
        inline smat bas() const {return basis;}   // the basis matrix
        inline scalar mod() const {return modulus;}   // the (prime) modulus

     // non-member (friend) functions and operators
        friend int dim(const ssubspace& s)     {return s.basis.ncols();}
        friend vec_i pivots(const ssubspace& s)  {return s.pivots;}
        friend smat basis(const ssubspace& s)  {return s.basis;}
	friend ssubspace combine(const ssubspace& s1, const ssubspace& s2);
	friend smat restrict_mat(const smat& m, const ssubspace& s);

// Implementation
private:
       scalar modulus;
       vec_i pivots;
       smat basis;
};


// Declarations of nonmember, nonfriend operators and functions:

ssubspace kernel(const smat& sm, scalar m);
ssubspace eigenspace(const smat& sm, scalar lambda, scalar m);
ssubspace subeigenspace(const smat& sm, scalar l, const ssubspace& s, scalar m);

// construction of a 1-dimensional sparse subspace from a vector:
ssubspace make1d(const vec& bas, scalar&piv, scalar m);
