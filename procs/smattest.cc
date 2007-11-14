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
#include <time.h>
#include <sys/times.h>
//#include <sys/param.h>
//#include <iostream>
#include "arith.h"

//#define SCALAR_OPTION 1 // ints
#define SCALAR_OPTION 2   // longs

const int MAX_FILL=25;

#include "subspace.h"
#include "smatrix.h"
#include "smatrix_elim.h"

long starttime,stoptime;

int main(void)
{ 
  cout<<"BIGPRIME = "<<BIGPRIME<<endl;
  starttime = clock();
  //  cout<<"raw start time = "<<starttime<<endl;
  //  cout<<"CLOCKS_PER_SEC = "<<CLOCKS_PER_SEC<<endl;

  cout << "enter 0 to exit\n";
  cout << "enter 1 to do all tests \n";
  cout << "enter 2 to operations \n"; 
  cout << "enter 3 to special speed test \n" << "enter> "<<endl;
  cout << "enter 4 to elimination test \n" << "enter> "<<endl;
  cout << "enter 5 to kernel test \n" << "enter> "<<endl;
  cout << "enter 6 to multiply smat by smat \n" << "enter> "<<endl;
  cout << "enter 7 to find eigenspaces \n" << "enter> "<<endl;
  int t,i;
  scalar seed=10;
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
	cout << "Unitialized new smat sm = " << flush;
	smat sm(nr,nc);
	cout << "(constructor finished)"<< endl;
	cout << sm << endl;
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
	cout << "Copy of smat = " << endl;
	smat sm2 = sm1;
	cout << sm2 << endl;
	cout << "Enter any number "; cin >> i;
  
//********  test of assignment   **************
        smat sm3;
	cout << "Copy using assignment= " << endl;
	sm3 = sm1;
	cout << sm3 << endl;
	cout << "Enter any number "; cin >> i;
	
//********* testing function as_matrix *****
	cout << "enter number of rows of smat for conversion " << endl;
	int numRow;
	cin >> numRow;
	cout << "Now enter number of columns" << endl;
	int numCol;
	cin >> numCol;
	smat T( numRow, numCol );
	/*
	cout << " Now enter smat for conversion" << endl;
	cout << " for each row, enter first value then position for each of the entries" << endl;
	cout << " for each row terminate input by typing zero as a value." << endl;
	cin >> T;
	*/
	random_fill_in(T,MAX_FILL,seed);
	cout << endl << "smat = " << endl << T << endl;
	cout << "smat as a matrix = "<< endl << T.as_mat( ) << endl;
	cout << "Row support sets:\n";
	vector<std::set<int> > supps=row_supports(T);
	for(i=1; i<=numRow; i++) cout << i << ": "<< supps[i] << endl;
	cout << "Enter any number "; cin >> i;

	//*********  testing ref to (i,j) entry **********
	cout << "Enter position (row,col): ";
	int j,k;
	cin >> j;
	cin >> k;
	cout << "T(" << j <<" , " << k << ") = " << T.elem(j,k) <<endl;
	
	}
   
//*********testing operations ***********

    if( t == 1 || t == 2 )
      {
	cout << "testing operations " << endl;
 
	int row,col;
	cout << "Enter size of smat A row,col: "<< endl; 
	cin >> row >> col;
	smat A( row, col );
	/*
	cout << "Enter entries of A: "<< endl; 
	cin >> A;
	*/
	random_fill_in(A,MAX_FILL,seed);
	cout << "Enter size of smat B: row,col: "<< endl; 
	cin >> row >> col;
	smat B( row, col );
	/*
	cout << "Enter entries of B: "<< endl; 
	cin >> B;
	*/
	random_fill_in(B,MAX_FILL,seed);
	cout << "matrix A" << endl << A << endl;
	cout << "matrix B" << endl << B << endl;
	cout << "Enter any number "; cin >> i;
	smat C = A;
	cout << "C = A = " << endl << C << endl;
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
	cout << "Now A = " << A << "\nand B = " << B << endl;
	cout << "Enter any number "; cin >> i;
	cout << "A-B=" << (A-B) << endl;
	cout << "Now A = " << A << "\nand B = " << B;


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
	cout << "enter smat number of rows, cols: " << endl;
	cin >> row >> col;
	smat sm(row,col);
	/*
	cin >> sm;
	*/
	random_fill_in(sm,MAX_FILL,seed);

	cout << " now enter matrix: nrows, ncols and then the entries"<<endl;
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
	
	A = sm;//at(m);
	cout<<"Testing transpose function"<<endl;
	cout << "A = "<<A<<" = "<<A.as_mat()<<endl;
	smat At = transpose(A);
	cout << "A^t = "<<At<<" = "<<At.as_mat()<<endl;
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
	/*
	cout << "now enter first smat;" << endl;
	cin >> sm1;
	*/
	random_fill_in(sm1,MAX_FILL,seed);
	cout<<"sm1 = "<<sm1<<endl;
	cout << "enter dimension second smat;" << endl;
	cin >> row2 >> col2;
	smat sm2(row2, col2);
	/*
	cin >> sm2;
	*/
	random_fill_in(sm2,MAX_FILL,seed);
	cout<<"sm2 = "<<sm2<<endl;
	cout << " the product is: " << endl;
	cout << "product = " << flush; 
	smat prod = sm1*sm2;
	cout << prod << endl;
	cout << "product mod " << BIGPRIME << " = " << flush;
	smat prod_mod_p = mult_mod_p(sm1,sm2,BIGPRIME);
	cout << prod_mod_p << endl;
	if(prod==prod_mod_p)
	  cout<<"--agree"<<endl;
	else
	  cout<<"--DO NOT agree"<<endl;

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
	  cout << "finished creating matrices"<<endl;
	  smat sm1(m1); smat sm2 (m2); smat sm3 (m3); smat sm(row,col);
	  cout << "finished creating smats"<<endl;
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

	  cout <<"Starting loop using smats..."<<flush;
	  starttime=clock();
	  for(j = 0; j < loop; j++) 
	    { 
	      sm += sm1;
	      sm -= 2*sm2; 
	      sm += 3*sm3; 
	      if(see>1) cout<<sm<<endl; else cout<<"."<<flush;
	    }
	  stoptime=clock();

	  cout << "...done"<<endl;
	  if( see ) cout << " resulting smat = " << sm << endl;
	  cout << "smat cpu time = " 
	    //	       << ((stoptime-starttime)) << " ticks = "
	       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
	       << " seconds"<<endl;
	  sm1 = smat(); sm2 = smat(); sm3 = smat();

	  cout <<"Starting loop using mats..."<<flush;
	  starttime=clock();  
	  for(j = 0; j < loop; j++) 
	    {  
	      m += m1;
	      m -= 2*m2; 
	      m += 3*m3; 
	      if(see>1) cout<<m<<endl; else cout<<"."<<flush;
	    }
	  stoptime=clock();

	  cout << "...done"<<endl;
	  if( see ) cout << " resulting matrix = " << m << endl;
	  cout << "matrix cpu time = " 
	    //	       << ((stoptime-starttime)) << " ticks = "
	       << ((double)(stoptime-starttime)/CLOCKS_PER_SEC) 
	       << " seconds"<<endl;
	  m1.init(0,0); m2.init(0,0); m3.init(0,0);
	  
	  if( sm == m ) cout << "results are equal" << endl;
	  else cout << "problem: results are not equal ! " << endl;
	  }

	  if( t == 4 )
	  {
	    cout << "enter size of matrix for elimination (row,col) "<< endl;
	    int nro,nc,max;
	    cin >> nro; cin >> nc;
	    smat sm(nro,nc);
	    //	    cout << "enter matrix as an smat" << endl;
	    //	    cin >> sm;
	    cout << "How many maximum number of non-zero entries per row?\n";
	    cin >> max;
	    cout << "enter seed for random number generator\n";
	    cin >> seed;
	    random_fill_in(sm,max,seed);
	    cout << "display matrices? ( 0 = no; 1 = yes )" << endl;
	    int flag, flag2;
	    cin >> flag;
	    cout << "display fill-in information?" << endl;
	    cin >> flag2;
	    if( flag ) cout << "matrix A is " << sm;
	    if(flag>1) cout << " = " << sm.as_mat();
	    if( flag ) cout << endl;
	    smat_elim A (sm );
	    A.init_elim();
	    if( flag2 ) 
	      { cout << "initial population: "<< get_population(A)<<endl; }
	    A.show_progress();
	    for(int pass=1; pass<=3; pass++){cout<<"Pass "<<pass<<endl;
	    for(int fc=1; fc<=2; fc++)
	      {
		A.elim_light_rows(fc);
		cout << "after elim_light_rows("<<fc<<")" << endl; 
		cout << "# of rows eliminated so far: " << A.get_rank( ) << endl;
		A.show_progress();
		if( flag ) cout << " matrix: " << A << endl;
		if(flag>1) cout << " = " << A.as_mat();
		if( flag2 ) cout<<"("<<get_population( A )<<" entries)"<<endl;
		A.elim_light_cols(fc);
		cout << "after elim_light_cols("<<fc<<")" << endl; 
		cout << "# of rows eliminated so far: " << A.get_rank( ) << endl;
		A.show_progress();
		if( flag ) cout << " matrix: " << A << endl;
		if(flag>1) cout << " = " << A.as_mat();
		if( flag2 ) cout<<"("<<get_population( A )<<" entries)"<<endl;
	      }
	    }
	    A.step4new();
	    A.reduce_mod_p();
	    cout << "after step4" << endl;
	    cout << "# of rows eliminated so far: " << A.get_rank( ) << endl;
	    A.show_progress();
	    if( flag ) cout << " matrix: " << A << endl;
	    if(flag>1) cout << " = " << A.as_mat();
	    if( flag2 ) cout<<"("<<get_population( A )<<" entries)"<<endl;
	    A.step5( );
	    A.reduce_mod_p();
	    cout << "after step5 " << endl;
	    cout << "rank of matrix: " << A.get_rank( ) << endl;
	    A.show_progress();
	    if( flag ) cout << " matrix: " << A << endl;
	    if(flag>1) cout << " = " << A.as_mat();
	    if( flag2 ) cout<<"("<<get_population( A )<<" entries)"<<endl;
	    cout<<"Echelon form?         "<<A.check_echelon()<<endl;
	    cout<<"Reduced echelon form? "<<A.check_red_echelon()<<endl;
	    A.step6( );
	    A.reduce_mod_p();
	    cout << "after step6 " << endl;
	    cout << "rank of matrix: " << A.get_rank( ) << endl;
	    A.show_progress();
	    if( flag ) cout << " matrix: " << A << endl;
	    if(flag>1) cout << " = " << A.as_mat();
	    if( flag2 ) cout<<"("<<get_population( A )<<" entries)"<<endl;
	    cout<<"Echelon form?         "<<A.check_echelon()<<endl;
	    cout<<"Reduced echelon form? "<<A.check_red_echelon()<<endl;
	  }
    if( t == 5 )
      {
	cout << "test of kernel function" << endl;
	cout << "enter size of matrix for elimination (row,col) "<< endl;
	int nro,nc;
	cin >> nro; cin >> nc;
	smat sm(nro,nc);
	int rand=1;
	//	cout << "Do you want to input the matrix for elimination or do you want\n"; 
	//	cout << "a matrix with random entries? (1 for random and zero otherwise\n";
	//	cin >> rand;
	int flag=1, max;
	cout << "want to determine the rank using matrix? (1 = yes; 0 = no )\n";
	cin >> flag;
	if( rand ) 
	  {
	    cout << "How many maximum number of non-zero entries per row?\n";
	    cin >> max;
	    cout << "enter seed for random number generator\n";
	    cin >> seed;
	    cout << "calculating matrix,\n";
	    random_fill_in( sm, max, seed );
	    cout << "done\n";
	    //	    cout<<"Matrix = "<<sm<<"\n=\n";sm.as_mat().output_pari(cout);
	  }
	//	else
	//	  {
	//	    cout << "enter matrix as an smat" << endl;
	//	    cin >> sm;
	//	  }

	vec pc, npc;
	mat ker_bas;	
	if( flag ) {
	  long rk, ny; scalar pr = BIGPRIME;
	  mat ker_mat = echmodp( sm.as_mat(), pc, npc, rk, ny, pr);
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
	  //	  cout << "mat kernel = " << ker_mat<<endl;
	  ker_bas = basis(pkernel(sm.as_mat(),BIGPRIME));
	  //	  cout << "Basis = "<<ker_bas<<endl;
	  //	  cout << "A*kernel = " << matmulmodp(sm.as_mat(),ker_bas,BIGPRIME)<<endl;
	}
	smat_elim A( sm );
	smat result;
#define COMPARE_KERNEL
#ifdef COMPARE_KERNEL
	cout << "Computing kernel (old method)"<<endl;
	smat ker = A.oldkernel(pc,npc);
	cout << "rank is:   " << dim( pc ) << endl;
	cout << "nullity is:" << dim( npc ) << endl;
	cout<<"A now has "<<get_population( A )<<" entries"<<endl;
	cout<<"kernel has "<<get_population( ker )<<" entries"<<endl;
	//	cout<<"old kernel = "<<ker.as_mat()<<endl;
	/*
	result = mult_mod_p(sm,ker,BIGPRIME);
	result.reduce_mod_p();
	if( result == smat(nc,dim(npc)) ) cout << "old kernel correct\n";
	else cout << "PROBLEM : old kernel not correct!\n";
	*/
	A=smat_elim(sm);
#endif
	cout << "Computing kernel (new method)"<<endl;
	smat newker = A.kernel(pc,npc);
	cout << "rank is:   " << dim( pc ) << endl;
	cout << "nullity is:" << dim( npc ) << endl;
	cout<<"A now has "<<get_population( A )<<" entries"<<endl;
	cout<<"kernel has "<<get_population( newker )<<" entries"<<endl;
	//	cout<<"new kernel = "<<newker.as_mat()<<endl;
#ifdef COMPARE_KERNEL
	if(ker==newker) cout<<"old and new kernels are equal"<<endl;
	else 
	  {
	    cout<<"old and new kernels DISAGREE"<<endl;
	    //	    cout<<"old-new="<<(ker-newker).as_mat()<<endl;
	  }
#endif
	/*
	result = mult_mod_p(sm,newker,BIGPRIME);
	result.reduce_mod_p();
	if( result == smat(nc,dim(npc)) ) cout << "new kernel correct\n";
	else cout << "PROBLEM : new kernel not correct!\n";
	*/
	if(flag)
	  {
	    smat ker1(ker_bas); ker1.reduce_mod_p();
	    result = mult_mod_p(sm,ker1,BIGPRIME);
	    result.reduce_mod_p();
	    if( result == smat(nc,dim(npc)) ) cout << "non-sparse kernel correct\n";
	    else cout << "PROBLEM : new kernel not correct!\n";
	//	    if(ker1==ker) cout<<"Agrees with non-sparse result"<<endl;
	//	    else {
	//	      cout<<"Does NOT agree with non-sparse result"<<endl;
	//	      cout<<"Non-sparse kernel basis = "<<ker1.as_mat()<<endl;
	//	      cout<<"Sparse kernel basis     = "<<ker<<endl;
	    }
      }      
    if(t==7)
      {
	cout << "Test of kernel and eigenspace operations\n";
	int row;
	cout << "Enter size of (square) smat A: "<< endl; 
	cin >> row;
	smat A( row, row );
	random_fill_in(A,MAX_FILL,seed);
	//	cout << "Enter entries of A: "<< endl; 
	//	cin >> A;
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

