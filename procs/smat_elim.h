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

#include <queue>
 
class smat_elim : public smat{

 private:
  
  int rank;                  // initially 0; incremented whenever a
			     // row is eliminated;
  int nrows_left, ncols_left; // # rows/cols not yet eliminated
  int ech_form, red_ech_form; // flags to show that reduction has taken place
  vector<std::set<int> > column;  // an array of ordered lists, one per
			     // col, of row nos which have nonzero
			     // entries in that col;
  vector<int> position;     // array indexed by row numbers r;
			    // initially -1; 0 iff row is zero; else c
			    // where row r was eliminated using pivot
			    // at (r,c) after clearing rest of column c
  vector<int> elim_col;     // array indexed by column numbers c;
			    // initially -1; else r where row r was
			    // eliminated using pivot at (r,c) after
			    // clearing rest of column c;
  vector<int> elim_row;     // array indexed by 1..rank;
			    // initially 0; else gives the order of
			    // elimination of rows;
  vector<int> light_col_flag;   // 0/1 flag columns which are `light' : not yet
			    // eliminated but no more than M entries
  queue<int> light_cols, light_rows; // lists of (very) light rows/cols
  void clear_col(int,int, int fr = 0, int fc = 0, int M = 0, int frl=0, int fcl=0);
  int get_weight( int row ); // compute the number of "light" columns
			     // intersecting row# row
  void eliminate( int&, int& );
  int finished() {return (nrows_left==0)&&(ncols_left==0);}
  int step4finished();

public:
  smat oldkernel( vec&, vec& );
  smat kernel( vec&, vec& );
  // constructor:
  smat_elim( const smat& sm) : smat( sm ) { ; };
  //  smat_elim( int r = 0,int c = 0 );
  // destructor:
  ~smat_elim() {;}
  int get_rank( ) { return rank; }
  void init_elim( );
  void elim_light_rows(int fr); // eliminate (row,col) for rows of wt <=fr
  void elim_light_cols(int fc); // eliminate (row,col) for columns of wt <=fc
  void step0() {elim_light_rows(1);}
  void step1() {elim_light_cols(1);}
  void step2() {elim_light_rows(2);}
  void step3() {elim_light_cols(2);}
  void step4(); // eliminate (row,col) for "light" columns 
  void step4new(); // eliminate (row,col) for "light" columns 
  void step5(); // all remaining elimination
  void step6(); // back-substitution
  void echelon_form( ); // steps 0,1,2,3,4,5 in turn
  void reduced_echelon_form( ); // steps 0,1,2,3,4,5,6 in turn
  int check_echelon();       // check that we do have echelon form
  int check_red_echelon();   // check that we do have reduced echelon form
  void show_progress() 
  {
    cout<<"#rows left = "<<nrows_left<<"\n#cols left = "<<ncols_left<<endl;
    cout<<"Population = "<<get_population(*this)<<endl;
  }
};

int rank(smat& sm);

// nullity of sm-lambda*I
inline int nullity(const smat& sm, const scalar& lambda) 
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
