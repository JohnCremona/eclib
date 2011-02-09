// smattest.cc: test of sparse matrix package
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

//#include <sys/types.h>
#include <sys/times.h>
//#include <sys/param.h>
//#include <iostream>
#include "random.cc"
#include "arith.h"

//#define SCALAR_OPTION 1 // ints
#define SCALAR_OPTION 2   // longs

#include "smatrix_elim.h"

long starttime,stoptime;

void random_fill_in( smat& sm, int max, scalar seed )
{
  int *intpos = new int [sm.nco];
  scalar *scalarval = new scalar [sm.nco];
  for( int r = 0; r < sm.nro; r++ )
    {
      int *lp = intpos; scalar *lv =scalarval;
      for( int i = 0; i < sm.nco; i++ ) { *lp++ = 0; *lv++ = 0; }
      int count = 0;
      int N = int( (max+1) * ran0( seed ) ); //number of entries in row i
      if( N == (max+1) ) N--;  // could occur !
      for( int s = 0; s < N; s++ )
	{
	  int v = int( 20 * ran0( seed ) ) - 10; // value between -10 & 9
	  if( v != 0 ) 
	    {  
	      int p = int( sm.nco * ran0( seed ) ); //position in matrix
	      if( p == sm.nco ) p--;
	      if( intpos[ p ] == 0 ) count++;
	      intpos[ p ] = 1;
	      scalarval[ p ] = v;
	    }
	}
      delete[] sm.col[r]; delete[] sm.val[r];
      int *ptr = sm.col[r] = new int [ count + 1 ];
      scalar *vptr = sm.val[r] = new scalar [ count ];
      *ptr++ = count;
      for( int l = 0; l < sm.nco; l++ )
	{
	  if( intpos[l] != 0 ) 
	    {
	      *ptr++ = l+1;
	      *vptr++ = scalarval[ l ];
	    }
	}
    }
  delete[] intpos;
  delete[] scalarval;
}

int main(void)
{ 
  cout << "enter 0 to exit\n";
  cout << "enter 1 to do all tests \n";
  cout << "enter 2 to operations \n"; 
  cout << "enter 3 to special speed test \n" << "enter> "<<endl;
  cout << "enter 4 to elimination test \n" << "enter> "<<endl;
  cout << "enter 5 to kernel test \n" << "enter> "<<endl;
  cout << "enter 6 to multiply smat by smat \n" << "enter> "<<endl;
  cout << "enter 7 to find eigenspaces \n" << "enter> "<<endl;
  int t,i;
  cin >> t;
  
  while( t != 0 ) {
    if( t == 1 )
      {
	int nr, nc;
	cout << "Test run of sparse matrix package.\n\n"; 
	cout << "enter # of rows:\n";
	cin >> nr;
	cout << "enter # of colums:\n";
	cin >> nc;
	smat sm(nr,nc);
	cout << "Unitialized new smat sm = " << sm << endl;
	cout << "Enter any number "; cin >> i;
//******** test of constructor from matrix  ***************
	//mat a; 
	int r;
	cout << "Enter size of a square matrix A: "; 
	cin >> r;
	mat a(r,r);
	cout << "Enter entries of A: "; 
	cin >> a;
	cout << "A = " << a;
	smat sm1(a);
	cout << "smat from matrix A = " << endl << sm1 << endl;
	cout << "Enter any number "; cin >> i;

 //*******  test of copy constructor  *******************
	smat sm2 = sm1;
	cout << "Copy of smat = " << endl << sm2 << endl;
	cout << "Enter any number "; cin >> i;
  
//********  test of assignment   **************
        smat sm3;
	sm3 = sm1;
	cout << "Copy using assignment= " << endl << sm3 << endl;
	cout << "Enter any number "; cin >> i;
	
//********* testing function as_matrix *****
	cout << "enter number of rows of smat for conversion " << endl;
	int numRow;
	cin >> numRow;
	cout << "Now enter number of columns" << endl;
	int numCol;
	cin >> numCol;
	smat T( numRow, numCol );
	// cout << " Now enter smat for conversion" << endl;
	// cout << " for each row, enter first value then position for each of the entries" << endl;
	// cout << " for each row terminate input by typing zero as a value." << endl;
	// cin >> T;
        int max=10, seed = 10;
        random_fill_in( T, max, seed );
	cout << endl << "smat = " << endl << T;
	cout << "smat as a matrix = "<< endl << T.as_mat( ) << endl;
	cout << "Enter any number "; cin >> i;

	//*********  testing ref to (i,j) entry **********
	cout << "Enter position: (row,col)";
	int j,k;
	cin >> j;
	cin >> k;
	cout << "T(" << j <<" , " << k << ") = " << T.elem(j,k) <<endl;
	
	// *** testing set_row *** //
	cout << "testing set_row" << endl;
	cout << "enter new row : which row ? (starting from zero)" << endl;
	int i, d;
	cin >> i;
	cout << "number of non-zero elements ?" << endl;
	cin >> d;
	int n = d;
	int *pos = new int [d+1];
	scalar *val = new scalar [d];
	cout << "values ? " << endl;
	while( n-- ) cin >> *val++;
	cout << " positions ? " << endl;
	n = d;
	while( n-- ) cin >> *pos++;
	pos -= d; val -= d;
	T.set_row(i, d, pos, val);
	cout << " new matrix : " << endl;
	cout << T << endl;
	delete[] pos; delete[] val;
	}
   
//*********testing operations ***********

    if( t == 1 || t == 2 )
      {
	cout << "testing operations " << endl;
 
	int row,col;
	cout << "Enter size of smat A row,col: "<< endl; 
	cin >> row >> col;
	smat A( row, col );
        random_fill_in( A, 10, 10 );
	// cout << "Enter entries of A: "<< endl; 
	// cin >> A;
	cout << "Enter size of smat B: row,col: "<< endl; 
	cin >> row >> col;
	smat B( row, col );
        random_fill_in( B, 10, 10 );
	// cout << "Enter entries of B: "<< endl; 
	// cin >> B;
	cout << "matrix A" << endl << A << endl;
	cout << "matrix B" << endl << B << endl;
	cout << "Enter any number "; cin >> i;
	smat C = A;
	cout << "C = A = " << endl << C;
	cout << "Enter any number "; cin >> i;
	cout << "B==A?" << (B==A) << endl;
	cout << "B!=A?" << (B!=A) << endl;
	cout << "C==A?" << (C==A) << endl;
	cout << "C!=A?" << (C!=A) << endl;
	cout << "Enter any number "; cin >> i;
	B+=A;
	cout << "after B+=A, A = " << endl << A <<endl <<"and B = "<< endl << B<<endl;
	cout << "Enter any number "; cin >> i;
	B-=A;
	cout << "after B-=A, A = " << A <<endl <<"and B = "<< B<<endl;
	cout << "Enter any number "; cin >> i;
	B*=2;
	cout << "after B*=2, B = " << B << endl;
	cout << "Enter any number "; cin >> i;
	B/=2;
	cout << "after B/=2, B = " << B << endl;
	cout << "Enter any number "; cin >> i;
	cout << "A+B=" << (A+B) << endl;
	cout << "Now A = " << A << "and B = " << B << endl;
	cout << "Enter any number "; cin >> i;
	cout << "A-B=" << (A-B) << endl;
	cout << "Now A = " << A << "and B = " << B;


	cout << "test addition of scalar to smat" << endl;
	B = A;
	scalar sc = 17;
	B += sc;
	cout<<"After adding 17 to A it is now:\n"<<B<<endl;
	B -= sc;
	cout<<"After subtracting 17 again it is now:\n"<<B<<endl;
	if(B==A) 
	  cout<<"scalar addition OK\n"; 
	else 
	  cout<<"scalar addition WRONG\n";
	A = smat(5);
	B = A;
	cout << "A = 5x5 zero smat:\n"<<B<<endl;
	B += sc;
	cout<<"After adding 17 to A it is now:\n"<<B<<endl;
	B -= sc;
	cout<<"After subtracting 17 again it is now:\n"<<B<<endl;
	if(B==A) 
	  cout<<"scalar addition OK\n"; 
	else 
	  cout<<"scalar addition WRONG\n";

	cout << "test multiplication of smat by matrix" << endl;
	cout << "enter smat: first row, col and then the rows" << endl;
	cin >> row >> col;
	smat sm(row,col);
        random_fill_in( sm, 10, 10 );
	// cin >> sm;
	cout << " now enter matrix: first row, col and then the entries"<<endl;
	cin >> row >> col;
	mat m(row, col);
	cin >> m;
	mat sm_as_m = sm.as_mat();
	cout << "the smat is (as a matrix) \n"<<sm_as_m<<endl;
	cout << "the matrix is \n"<<m<<endl;
	mat sm_times_m = sm*m;
	cout << " the product matrix is: \n" << sm_times_m << endl;
	if(sm_times_m==sm_as_m*m)
	  cout<<"Correct\n";
	else
	  cout<<"Wrong! Correct is \n"<<sm_as_m*m<<endl;
	
	A = smat(100,100);
        random_fill_in( A, 10, 10 );
	cout<<"Testing transpose function"<<endl;
	// cout << "A = "<<A<<" = "<<A.as_mat()<<endl;
	smat At = transpose(A);
	// cout << "A^t = "<<At<<" = "<<At.as_mat()<<endl;
	smat Att = transpose(At);
	cout<<"A=(A^t)^t ? ";
	if(A==Att) cout<<"yes!"<<endl; else cout<<"no!"<<endl;
      } 
    if( t == 6 ) {
	/* test of multiplication of smat by smat */
	int row1, row2, col1, col2;
	cout << "test multiplication of smat by smat" << endl;
	cout << "enter dimension first smat;" << endl;
	cin >> row1 >> col1;
	smat sm1(row1,col1);
	cout << "now enter first smat;" << endl;
	cin >> sm1;
	cout << "enter dimension second smat;" << endl;
	cin >> row2 >> col2;
	smat sm2(row2, col2);
	cin >> sm2;
	cout << " the product is: " << endl;
	cout << sm1*sm2 << endl;
	}

       if ( t == 3 )
	{
          int row,col;
	  cout << "enter size of matrices for speed test (row,col) " << endl;
	  cin >> row; cin >> col;
	  mat m1( row,col);
	  mat m2( row,col);
	  mat m3(row,col);
	  mat m(row,col);
	  for( int r=1; r<=row; r++ ) {
	    for( int c = 1; c<=5 && r+c < col; c++ ) {
	      m1.set( r, r+c,1 );
	      m2.set( r, r+c-1, 2 );
	      m3.set( r, r+c+1, -1 );
	    }
	  }
	  smat sm1(m1); smat sm2 (m2); smat sm3 (m3); smat sm(row,col);
	  cout << "want to see matrices ? ";
	  int see; cin >> see;
	  if( see )
	  {
	    cout << "sm1 = " << sm1 << endl;
	    cout << "sm2 = " << sm2 << endl;
	    cout << "sm3 = " << sm3 << endl;
	    }
	  cout << "Loop how many times ?" << endl;
	  int j, loop; cin >> loop;

	  starttime = clock();
	  for(j = 0; j < loop; j++) 
	    {  
	      m += m1;
	      m -= 2*m2; 
	      m += 3*m3; 
	      if(see>1) cout<<m<<endl; else cout<<"."<<flush;
	    }  
	  stoptime = clock();

	  if( see ) cout << " resulting matrix = " << m << endl;
	  cout << "matrix cpu time = " 
	       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
	       << " seconds"<<endl;
	  m1.init(0,0); m2.init(0,0); m3.init(0,0);
	  
	  starttime = clock();
	  for(j = 0; j < loop; j++) 
	    { 
	      sm += sm1;
	      sm -= 2*sm2; 
	      sm += 3*sm3; 
	      if(see>1) cout<<sm<<endl; else cout<<"."<<flush;
	    } 		  
	  stoptime = clock();

	  if( see ) cout << " resulting smat = " << sm << endl;
	  cout << "smat cpu time = " 
	       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
	       << " seconds"<<endl;
	  sm1 = smat(); sm2 = smat(); sm3 = smat();
	  if( sm == smat(m) ) cout << "results are equal" << endl;
	  else cout << "problem: results are not equal ! " << endl;
	  }

	  if( t == 4 )
	  {
	    cout << "enter size of matrix for elimination (row,col) "<< endl;
	    int nro,nc;
	    cin >> nro; cin >> nc;
	    smat sm(nro,nc);
	    cout << "enter matrix as an smat" << endl;
	    int max=10, seed = 10;
	    random_fill_in( sm, max, seed );
	    //	    cin >> sm;
	    cout << "display matrices? ( 0 = no; 1 = yes )" << endl;
	    int flag, flag2;
	    cin >> flag;
	    cout << "display fill-in information?" << endl;
	    cin >> flag2;
	    if( flag ) cout << "matrix A is " << sm << endl;
	    smat_elim A (sm );
	    if( flag2 ) 
	      { cout << "initial population: "; display_population(A); }
	    A.step0();
	    cout << "after step0 " << endl;
	    cout << "# of rows eliminated so far: " << A.get_rank( ) << endl;
	    if( flag ) cout << " matrix: " << A << endl;
	    if( flag2 ) display_population( A );
	    A.step1();
	    cout << "after step1" << endl;
	    cout << "# of rows eliminated so far: " << A.get_rank( ) << endl;
	    if( flag ) cout << " matrix: " << A << endl;
	    if( flag2 ) display_population( A );
	    A.step2();
	    cout << "after step2" << endl;
	    cout << "# of rows eliminated so far: " << A.get_rank( ) << endl;
	    if( flag ) cout << " matrix: " << A << endl;
	    if( flag2 ) display_population( A );
	    A.step3();
	    cout << "after step3" << endl; 
	    cout << "# of rows eliminated so far: " << A.get_rank( ) << endl;
	    if( flag ) cout << " matrix: " << A << endl;
	    if( flag2 ) display_population( A );
	    A.step4();
	    cout << "after step4" << endl;
	    cout << "# of rows eliminated so far: " << A.get_rank( ) << endl;
	    if( flag ) cout << " matrix: " << A << endl;
	    if( flag2 ) display_population( A );
	    A.standard( );
	    cout << "after standard " << endl;
	    cout << "rank of matrix: " << A.get_rank( ) << endl;
	    if( flag ) cout << " matrix: " << A << endl;
	    if( flag2 ) display_population( A );
	  }
    if( t == 5 )
      {
	cout << "test of kernel function" << endl;
	cout << "enter size of matrix for elimination (row,col) "<< endl;
	int nro,nc;
	cin >> nro; cin >> nc;
	smat sm(nro,nc);
	int rand;
	cout << "Do you want to input the matrix for elimination or do you want\n";
	cout << "a matrix with random entries? (1 for random and zero otherwise\n";
	cin >> rand;
	int flag;
	cout << "want to determine the rank using matrix? (1 = yes; 0 = no )";
	cin >> flag;
	if( rand ) 
	  {
	    cout << "How many maximum number of non-zero entries per row?\n";
	    int max, seed;
	    cin >> max;
	    cout << "enter seed for random number generator\n";
	    cin >> seed;
	    cout << "calculating matrix,\n";
	    random_fill_in( sm, max, seed );
	    cout << "done\n";
	  }
	else
	  {
	    cout << "enter matrix as an smat" << endl;
	    cin >> sm;
	  }
	smat_elim A( sm );
	vec pc, npc;
	
	if( flag ) {
	  long rk, ny; scalar pr = 92861;
	  mat m = sm.as_mat ();  
	  mat ker_mat = echmodp( m, pc, npc, rk, ny, pr);
	  cout << " rank using echmodp : " << rk;
	  int pop = 0;
	  int nro = nrows(ker_mat);
	  int nco = ncols(ker_mat);
	  for( int r = 1; r <= nro; r++ ) {
	    for( int c = 1; c <= nco; c++ ) {
	      pop += ( ker_mat( r, c ) != 0 );
	    }
	  }
	  cout << " number of non-zero entries: " << pop << endl;
	}
	/********A.step0 ();
	  cout << "step0 : A * ker_mat is : ";
	  if( (A*ker_mat) ==  smat(nro)) cout << "0";
	  else cout << "problem in step0!!!";
	  A.step1 ();
	  cout << "step1 : A * ker_mat is : ";
	  if( (A*ker_mat) ==  smat(nro)) cout << "0";
	  else cout << "problem in step1!!!";
	  A.step2();
	  cout << "step2 : A * ker_mat is : ";
	  if( (A*ker_mat) ==  smat(nro)) cout << "0";
	  else cout << "problem in step2!!!";
	  A.step3 ();
	  cout << "step3 : A * ker_mat is : ";
	  if( (A*ker_mat) ==  smat(nro)) cout << "0";
	  else cout << "problem in step3!!!";
	  A.step4 ();
	  cout << "step4 : A * ker_mat is : ";
	  if( (A*ker_mat) ==  smat(nro)) cout << "0";
	  else cout << "problem in step4!!!";
	  A.standard();
	  cout << "standard : A * ker_mat is : ";
	  if( (A*ker_mat) == smat(nro)) cout << "0";
	  else cout << "problem in standard!!!";*****/
	smat ker = A.kernel(pc,npc);
	cout << "rank is:" << dim( pc ) << endl;
	display_population(ker);
	smat result = mult_mod_p(sm,ker,BIGPRIME);
	// cout << "sm  is:\n" << sm.as_mat() << endl;
	// cout << "ker is:\n" << ker.as_mat() << endl;
	// cout << "result is:\n" << result.as_mat() << endl;
	if( result == smat(nro) ) cout << "kernel correct\n";
	else cout << "PROBLEM : kernel not correct!\n";
      }
    if(t==7)
      {
	cout << "Test of kernel and eigenspace operations\n";
	int row;
	cout << "Enter size of (square) smat A: "<< endl; 
	cin >> row;
	smat A( row, row );
	cout << "Enter entries of A: "<< endl; 
	cin >> A;
	cout << "A = \n"<<A <<endl;
	cout << "A (as matrix) = "; A.as_mat().output_pari(); cout <<endl;
	ssubspace ker = kernel(A);
	cout << "ker(A) has dimension " << dim(ker) << endl;
	cout << "basis =  " << basis(ker) << endl;

	cout << "Enter a possible eigenvalue: ";
	scalar lambda;
	cin >> lambda;
	cout<<"lambda = "<<lambda<<endl;
	ssubspace e = eigenspace(A,lambda);
	cout << "Eigenspace for lambda = "<<lambda<<" has dimension " << dim(e) << endl;
      }
    cout << "enter new value of t  ";
    cin >> t;
  }  
  cout<<endl;
  return(0);
}

