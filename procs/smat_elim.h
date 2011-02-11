// smat_elim.h: declaration of class smat_elim, for sparse elimination
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
 
// Not to be included directly by user: use smatrix_elim.h

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
    void put( type& X) { 
      if( num >= maxsize ) 
	{
// 	  cout<<"About to grow list from size "<<maxsize<<endl;
	  grow();
// 	  cout<<"Finished growing list, new size "<<maxsize<<endl;
	}
      list_array[ num++ ] = X; 
    }
    int find( type& X, int ub, int lb = 0 );
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
  
  int rank;
  ordlist* column; // an array of lists, oner per col, of row nos
		   // which have nonzero entries in that col
  int *position;     
  int *elim_col;     // in which row was col eliminated;
  int *elim_row;       //order of elimination of rows;
  void clear_col(int,int,list&, int fr = 0, int fc = 0,int M = 0,int *li =0);
  void check_col( int col, list& L );
  void check_row (int d2, int row2, list& L ); 
  int get_weight( int, int* ); 

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
  smat kernel( vec&, vec& );
  void normalize( int, int );
  void eliminate( int&, int& );
  void step5dense();
  void free_space( int col );
  void elim( int row1, int row2, scalar v2 );
  // constructor:
  smat_elim( const smat& sm) : smat( sm ) { init(); };
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

long rank(smat& sm);

inline long nullity(const smat& sm, const scalar& lambda) // nullity of sm-lambda*I
{
  smat sma(sm); sma-=lambda;  return ncols(sm)-rank(sma);
}

class ssubspace {

public:
     // constructors
        ssubspace(int n=0);
        ssubspace(const smat& b, const vec& p);
	ssubspace(const ssubspace& s);
     // destructor
        ~ssubspace();
     // assignment
	void operator=(const ssubspace& s);

     // member functions & operators
        inline void clear() { pivots.init(); basis=smat(0,0);}
        inline vec pivs() const {return pivots;}  // the pivot vector
        inline smat bas() const {return basis;}   // the basis matrix

     // non-member (friend) functions and operators
        friend int dim(const ssubspace& s)     {return ncols(s.basis);}
        friend vec pivots(const ssubspace& s)  {return s.pivots;}
        friend smat basis(const ssubspace& s)  {return s.basis;}  
	friend ssubspace combine(const ssubspace& s1, const ssubspace& s2);
	friend smat restrict(const smat& m, const ssubspace& s);

// Implementation
private:
       vec pivots;
       smat basis;
};


// Declarations of nonmember, nonfriend operators and functions:

ssubspace kernel(const smat& m);
ssubspace eigenspace(const smat& m, scalar lambda);
ssubspace subeigenspace(const smat& m, scalar l, const ssubspace& s);

// construction of a 1-dimensional sparse subspace from a vector:
ssubspace make1d(const vec& bas, long&piv);
