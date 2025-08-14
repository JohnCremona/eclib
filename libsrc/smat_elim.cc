// smat_elim.cc: implementation of class smat_elim
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
 
//  implements structured modular elimination
//  original written by Luiz Figueiredo

#include <eclib/timer.h>

#if(0)
// This special version only works modulo BIGPRIME, not a general modulus:

inline scalar xmm0(scalar a, scalar b)
{
  if (a==1) return b;
  if (a==-1) return -b;
  if (b==1) return a;
  if (b==-1) return -a;

  //return xmodmul(a,b,m);
  //return (a*b) % m;
  //return (scalar)(((long)a*(long)b) % (long)m);

  /* this works fine:
  return (a*(int64_t)b) % BIGPRIME;
  */
  // the following should work faster (no divisions!  Thanks to David Harvey)
  if(a<0) a+=BIGPRIME;
  if(b<0) b+=BIGPRIME;
 int64_t ab = a*(int64_t)b;
 int64_t r = ab-((INV_BIGPRIME*(ab>>30))>>32)*BIGPRIME;
 r -= ( ((r>=TWO_BIGPRIME)?BIGPRIME:0) + ((r>=BIGPRIME)?BIGPRIME:0) );
 /*
 // check:
 scalar r2 = (a*(int64_t)b) % BIGPRIME;
 if (r!=r2)
 {
 cout << "Problem with "<<a<<"*"<<b<<" (mod "<<BIGPRIME<<"): computed "<<r<<", not "<<r2<<endl;
 return r2;
}
 */
 return (scalar)r;
}

#endif

//#define TRACE_LISTS
//#define TRACE_FIND

int smat_elim::list::listsize;

void smat_elim::list::clear( int m) 
{ 
  delete [] list_array;
  list_array = new type [m]; num = 0; maxsize = m; index = 0;
}

smat_elim::list::list( int m) 
{ 
  list_array = new type [m]; num = 0; maxsize = m; index = 0;
}

smat_elim::list::~list( ) { delete [] list_array; }

void
smat_elim::list::grow()
{
  int growth = (maxsize==0? listsize : maxsize/2 + 1);
  type *new_array = new type [ maxsize + growth]; 
  if( !new_array )
    {
      cerr << "memory exhausted in elim::list::grow"<< endl;
      return;
    }
  type* newi = new_array;
  type *P = list_array;
  int s = maxsize;
  // while(s--) *newi++ = *P++;
  size_t n = s*sizeof(type);
  memmove(newi,P,n);
  maxsize += growth;
  delete [] list_array;  list_array=new_array;
}

int
smat_elim::list::find( const type& X, int ub, int lb )
{
  // returns highest  number i, lb <= i <= ub, such that list_array[i] <= X 
  // or returns ub+1 if list_array[ub]<X
  // or returns lb   if list_array[lb]>X

  int i;
#ifdef TRACE_FIND
  cout<<"\n\t\tfinding "<<X<<" in list "<<(*this)<<" from "<<lb<<" to "<<ub<<endl;
#endif
  if( list_array[ub] <  X ) 
    {
#ifdef TRACE_FIND
      cout<<"\t\tfind returns "<<(ub+1)<<endl;
#endif
      return ub+1;
    }
  while( list_array[lb] < X ) {
    i = (ub + lb)/2;
    list_array[i] < X ? (lb = i+1) : (ub = i);
  }
#ifdef TRACE_FIND
  cout<<"\t\tfind returns "<<(lb)<<endl;
#endif
  return lb;
}

void
smat_elim::ordlist::put( type& X )
{
#ifdef TRACE_LISTS
  cout<<"\tputting "<<X<<" into ordlist "<<(*this);
#endif
  if( num == maxsize ) grow();
  if( num == 0 ) {
    list_array[0] = X;
    num++;
  }
  else {
    int ind = find( X, num-1 );
    if( (ind==num)||(list_array[ind] != X )) {  // if X is not already in there
      // type *array = list_array + num -1;
      // for( int r = num; r > ind; r-- ) array[1] = *array--;
      // array[1] = X;
      type *source = list_array+num-1;
      type *dest = list_array+num;
      size_t n = sizeof(type)*(num-ind);
      memmove(dest, source, n);
      list_array[ind]=X;
      num++;
    }
  }
#ifdef TRACE_LISTS
  cout<<", result is "<<(*this)<<endl;
#endif
}

void
smat_elim::ordlist::put( list& L )   // L must be ordered
{
  if( L.num == 0 ) return;
#ifdef TRACE_LISTS
  cout<<"Inserting list "<<L<<" into ordlist "<<(*this)<<endl;
#endif
  L.index = index = 0;  //need to reset in case next() was used before.
  if( num == 0 ) 
    {
      for( int r = 0; r < L.num; r++ ) 
	{
	  type X = L.next();
	  this->put( X );
	}
      num = L.num;
    }
  else 
    {
      type *new_array = new type [ maxsize + L.num ];
      type *na = new_array;
      for( int r = 0, ind = 0; r < L.num; r++ ) 
	{
	  type X = L.next();
	  ind = find( X, num-1, ind );
	  if( list_array[ind] != X ) 
	    {
	      while( index < ind ) *na++ = next();
	      *na++ = X;
	    }	
	}	
      while( index < num ) *na++ = next();
      delete [] list_array;
      list_array = new_array;
      maxsize += L.num;
      L.index = index = 0;
      num = na - new_array;
    }
#ifdef TRACE_LISTS
  cout<<"Result is "<<(*this)<<endl;
#endif
}	

void
smat_elim::ordlist::remove( type& X )
{
#ifdef TRACE_LISTS
  cout<<"\tremoving "<<X<<" from ordlist "<<(*this);
#endif
  int ind = find( X, num-1 );
  if( list_array[ind] != X ) 
    { 
      cout<<endl;  
      cerr << "error in remove(1)\n"; 
      cerr<<"while removing "<<X<<" from "<<(*this)<<endl;
      return;
    }
  // type *array = list_array + ind;
  // for( int s = ind + 1; s < num; s++, array++ ) *array = array[1];
  type *source = list_array + ind +1;
  type *dest   = list_array + ind ;
  size_t n = sizeof(type)*(num-1-ind);
  memmove(dest,source,n);

  num--;
#ifdef TRACE_LISTS
  cout<<", result is "<<(*this)<<endl;
#endif
}

void
smat_elim::ordlist::remove( list& L )  // L must be ordered
{
  if( L.num == 0 ) return;
#ifdef TRACE_LISTS
  cout<<"Removing list "<<L<<" from ordlist "<<(*this)<<endl;
#endif
  L.index = 0;
  type X = L.next(); 
  int ind1 = find(X, num-1);
  int ind2 = ind1;
  if( list_array[ind1] != X )  
    { 
      cout<<endl;
      cerr << "error in remove(2)\n"; 
      cerr<<"while removing "<<L<<" from "<<(*this)<<endl;
      return; 
    }
  type *ar = list_array + ind1;
  index = ind1+1;
  for( int r = 1; r < L.num; r++ ) {
    X = L.next();
    ind2 = find( X, num-1, ind2 );
    if( list_array[ind2] != X )  
      { 
	cout<<endl;
	cerr << "error in remove(3)\n"; 
	cerr<<"while removing "<<L<<" from "<<(*this)<<endl;
	return; 
      }
    while( index < ind2 ) *ar++ = next();
    index++;
  }
  while( index < num ) *ar++ = next();
  L.index = index = 0;
  num = ar - list_array;
#ifdef TRACE_LISTS
  cerr<<"Result is "<<(*this)<<endl;
#endif
}

void smat_elim::init( )
{
  //  cout<<"smat_elim::init()  with smat:\n"<<(smat)(*this)<<endl;
  //  this->reduce_mod_p(modulus);
  //  cout<<"smat_elim::init()  after reducing:\n"<<(smat)(*this)<<endl;
  list::listsize = 10;
  rank = 0;
  position = new int[nro];
  int *p = position;
  elim_col = new int[nco];
  int* el = elim_col;
  elim_row = new int[nro];
  int* er = elim_row;
  column = new ordlist [nco];
  if( !column ) { cerr << "memory exhausted in smat_elim::init"<<endl; return; }
//   else {cout<<"Successfully created column array of length "<<nco<<endl;}
  int l,r;
  for( l = 0; l < nco; l++ ) *el++ = -1;
  for( r = 0; r < nro; r++ ) { *er++ = 0; *p++ = -1; }
  
  for( r = 0; r < nro; r++ ) {
    int d = *col[r];
    p = col[r] + 1;
    while( d-- ) (column + (*p++) - 1)->list::put(r);
  }
//   cout<<"At end of init(), columns are: \n";
//   for( l = 0; l < nco; l++ ) 
//     cout<<(l+1)<<": "<<column[l]<<"\n";
}

smat_elim::~smat_elim()
{
  delete [] position;
  delete [] elim_col;
  delete [] elim_row;
  delete [] column;
}

//#define TRACE_ELIM 1
//#define TRACE_DENSE 1

void smat_elim::sparse_elimination( )
{
#if TRACE_ELIM || TRACE_DENSE
  int pop=get_population(*this);
  double dens = density(*this); //(double)pop/(nco*nro);
  cout<<"Starting sparse elimination: "<<nro<<" rows, "<<nco<<" columns, "<<pop<<" entries (density = "<<dens<<")\n";
  //  cout<<"row weights:\n";
  //  for(int i=0; i<nro; i++) cout<<col[i][0]<<" ";
  //  cout<<endl;
  //  cout<<(*this)<<endl;
#endif
#if TRACE_ELIM
  cout<<"Starting step 0..."<<flush;
#endif
  step0();
#if TRACE_ELIM
  cout<<"finished\n"; 
  report();
  cout<<"Starting step 1..."<<flush;
#endif
  step1();
#if TRACE_ELIM
  cout<<"finished step 1\n";
  report();
  cout<<"Starting step 2..."<<flush;
#endif
  step2();
#if TRACE_ELIM
  cout<<"finished step 2\n";
  report();
  cout<<"Starting step 3..."<<flush;
#endif
  step3();
#if TRACE_ELIM
  cout<<"finished step 3\n";
  report();
  cout<<"Starting step 4..."<<flush;
#endif
  step4();
#if TRACE_ELIM
  cout<<"finished step 4\n";
  report();
#endif
#if(0)  // use dense method for final elimination
#if TRACE_ELIM || TRACE_DENSE
  cout << "Switching to dense mode..."<<endl;
#endif
  step5dense();
#if TRACE_ELIM || TRACE_DENSE
  cout<<"finished, ";
  report();
#endif
#else  // use sparse method for final elimination
#if TRACE_ELIM
    cout<<"Starting step 5 (remaining elimination)..."<<flush;
#endif
     standard( );
#if TRACE_ELIM
  cout<<"finished step 5 ";
  report();
#endif
#endif
}

smat smat_elim::kernel( vec_i& pc, vec_i& npc)
{
  return old_kernel(pc, npc);
}

// New version of kernel, not using back_sub() but constructing the
// kernel directly fro the "upper triangular" result of elim().

//#define TRACE_ELIM 1

smat smat_elim::new_kernel( vec_i& pc, vec_i& npc)
{
  int i,ir, j, jj, t, r, c;
  scalar v;
  static const scalar zero(0);
  static const scalar one(1);

#if TRACE_ELIM
  cout<<"Starting sparse_elimination()..."<<flush;
#endif
  sparse_elimination( );
#if TRACE_ELIM
  cout<<"finished sparse_elimination()"<<endl;
#endif

  int nullity = nco - rank;

  /* pc and npc hold the pivotal and non-pivotal column numbers, each
     is a vec, so indexed from 1, and the values are indexed from 1.

     The pivotal positions *in the order of elimination* are

     (elim_row[i]+1, position[elim_row[i]])  for 0<=i<rank

     where these row/col indices start at 1.
  */

  pc.init( rank );
  npc.init( nullity );

#if TRACE_ELIM
  cout<<"rank =    "<<rank<<endl;
  cout<<"nullity = "<<nullity<<endl;

  float dense = get_population(*this);
  dense /= (nro*nco);
  cout<<"density = "<<dense<<endl;
#endif

  /* set-up vecs pc & npc */
#if TRACE_ELIM
  cout<<"Finding pivotal and non-pivotal columns..."<<endl;
#endif

  /* ny is just a dummy index; after the loop it will equal nullity
     and k will equal rank: */
  int ny = 0, k = 0;

  /* find the pivotal and non-pivotal columns and the pivotal rows.

  */
  for( c = 1; c <= nco; c++ )  // loop through all columns
    {
      r = elim_col[c-1]+1;
      if( r > 0 )            // this is a pivotal column for row r
        {
          k++;               // the k'th pivot is (r,c)
          pc[k] = c;         // pc[k] is its column, c
        }
      else                     // non-pivotal column
        {
          npc[++ny] = c;       // record c as a non-pivotal column
        }
    }
#if TRACE_ELIM
  cout << "pivotal columns: "<<pc<<endl;
  cout << "non-pivotal columns: "<<npc<<endl;

  // density of the non-upper-triangular part:
  dense = 0;
  for (i=0; i<rank; i++)
    for (j=0; j<nullity; j++)
      if (elem(elim_row[i]+1, npc[j+1]) !=0)
        dense += 1;
  dense /= (rank*nullity);
  cout<<"density of block = "<<dense<<endl;
  start_time();
#endif


#if TRACE_ELIM
  cout<<"Constructing basis for kernel..."<<endl;
#endif

  /* We construct the basis matrix by rows for efficiency; the j'th
     column is the j'th basis vector. */

  /* There is a nullity x nullity identity matrix in rows indexed by
     npc, and the remaining entries, in rows position[elim_row[i]] for
     0<=i<rank, are given as follows, in column j for 1<=j<=nullity:

     set jj = npc[j];
     for i from rank-1 down to 0
         let ir = elim_row[i] and set

         basis[position[ir], j] =

         -M[ir+1,jj] - sum_{t=i+1}^{rank-1} M[ir+1,position[elim_row[t]]*basis[position[elim_row[t]], j]

     For fixed ir we put the entries M[ir+1,position[elim_row[t]] into array R.

   */

  smat bas( nco, nullity );

  /* First set the identity block */

  int *co;
  scalar *va;

  for(j=1; j<=nullity; j++)
    {
      jj = npc[j]-1;  // NB constructor gives this much
      bas.col[jj][0] = 1; // 1 entry in this row
      bas.col[jj][1] = j; // in column 1
      bas.val[jj][0] = one; // with value 1
    }

  /* set the other entries in order */

  // array to hold the rank*nullity dense entries
  // B[i] has length nullity for 0<=i<rank, and holds the entries
  // in row position(elim_row[i])-1 of basis

  scalar **B = new scalar*[rank];
  scalar **b = B;
  i = rank;
  while(i--)
    *b++ = new scalar[nullity];

  scalar *bij, *bij_nz;
  scalar *bi_nz = new scalar[nullity];
  scalar *R = new scalar[rank];
  scalar *Rt;
  scalar **Bt;
  scalar rr, ss;
  int *ij_nz;
  int *i_nz = new int[nullity];

  for(i=rank-1; i>=0; i--) // set B[i]
    {
      bij = B[i];
      bij_nz = bi_nz;
      ij_nz = i_nz;

      ir = elim_row[i];
      int nv=0; // counts # non-zero v

      for(t=0; t<rank; t++)
        R[t] = (t<i? zero : elem(ir+1, position[elim_row[t]]));

      for(j=0; j<nullity; j++) // set B[i][j], using B[t][j] for t>i
        {
          v = -elem(ir+1, npc[j+1]);
          Rt = R + rank-1;
          Bt = B + rank-1;
          t = rank-1-i;
          while(t--)
            {
              rr = *Rt--;
              if (is_nonzero(rr))
                {
                  ss = (*Bt)[j];
                  if (is_nonzero(ss))
                    {
                      v = mod(v - xmodmul(rr, ss, modulus), modulus);
                    }
                }
              Bt--; // must be outside the if(rr)
            }
          *bij++ = v;
          if (is_nonzero(v))
            {
              nv++;
              *bij_nz++ = v;
              *ij_nz++ = j+1;
            }
        }
#if TRACE_ELIM
      cout<<" setting row "<< position[ir]-1 <<" (from 0) of basis: "<< nv <<" non-zero entries out of "<<nullity<<endl;
#endif
      ir = position[ir]-1;
      co = bas.col[ir];
      va = bas.val[ir];
      if (nv > co[0]) // there are more entries in this row than the
                      // constructor gave us
        {
          delete [] co;
          delete [] va;
          co = bas.col[ir] = new int[nv+1];
          va = bas.val[ir] = new scalar[nv];
        }
      *co++ = nv;
      size_t nbytes = nv*sizeof(int);
      memmove(co,i_nz,nbytes);
      nbytes = nv*sizeof(scalar);
      memmove(va,bi_nz,nbytes);

#if TRACE_ELIM
      cout<<" finished setting row "<< ir << endl;
#endif
    }

  /* release memory dynamically allocated in this function using new */

  b = B;
  i = rank;
  while(i--)
    delete [] *b++;
  delete [] B;
  //  delete [] elim_row_inv;
  delete [] R;
  delete [] bi_nz;
  delete [] i_nz;

  #if TRACE_ELIM
  stop_time();
  cout<<"time for computing basis: ";
  show_time();
  cout<<endl;

  cout<<"Finished constructing basis for kernel"<<endl;
  cout<<" basis = "<<bas.as_mat()<<endl;
#endif

  bas.reduce_mod_p(modulus);
  return bas;
}

// old version of kernel which uses back_sub()

smat smat_elim::old_kernel( vec_i& pc, vec_i& npc)
{
  int i,n,r;
#if TRACE_ELIM
  cout<<"Starting sparse_elimination()..."<<flush;
#endif
  sparse_elimination( );
#if TRACE_ELIM
  cout<<"finished sparse_elimination()"<<endl;
#endif


  int nullity = nco - rank;
  if (nullity>0)
    {
#if TRACE_ELIM
  cout<<"Starting back-substitution..."<<flush;
#endif
      back_sub();
#if TRACE_ELIM
  cout<<"finished back-substitution"<<endl;
#endif

    }
  smat bas( nco, nullity );
  pc.init( rank );
  npc.init( nullity );

  /* set-up vecs pc & npc */
#if TRACE_ELIM
  cout<<"Setting up pc and npc..."<<flush;
#endif
  int ny = 0, k = 0;
  long *new_row = new long [ rank ];
  for( i = 1; i <= nco; i++ )
    {
      if( elim_col[i-1] > -1 )
        {
          pc[++k] = i;
          new_row[k-1] = elim_col[i-1];
        }
      else
        {
          npc[++ny] = i;
        }
    }

  /* write basis for kernel */
#if TRACE_ELIM
  cout<<"Constructing basis for kernel..."<<flush;
#endif
  for( n = 1; n <= nullity; n++ )
    { 
      i = npc[n]-1;
      bas.col[i][0] = 1;      //this much storage was granted in the
      bas.col[i][1] = n;      // in the constructor.
      bas.val[i][0] = 1;
    }

  scalar *aux_val = new scalar [nco];
  int *aux_col = new int [nco];
  for ( r=1; r<=rank; r++)
    { 
      i = pc[r]-1;
      int count = 0;
      int *axp = aux_col; scalar *axv = aux_val;
      int *posB = col[new_row[r-1]];
      int d = *posB++-1;
      scalar *valB = val[new_row[r-1]];
      for (int j = 1, h = 0; j<=nullity; j++) {
	while( *posB < npc[j] && h < d ) { posB++; h++; }
	if( *posB == npc[j] ) {	*axp++ = j; *axv++ = -valB[h]; count++; }
      }
      delete [] bas.col[i];
      delete [] bas.val[i];
      bas.col[i] = new int [count + 1];
      bas.val[i] = new scalar [count];
      int *pos = bas.col[i];
      scalar *vali = bas.val[i];
      axp = aux_col;
      axv = aux_val;
      *pos++ = count;

      // Using memmove only works when scalar is int or long, not bigint
      for( n = 0; n < count; n++ ) { *pos++ = *axp++; *vali++ = *axv++; }
      // size_t nbytes = count*sizeof(int);
      // memmove(pos,axp,nbytes);
      // nbytes = count*sizeof(scalar);
      // memmove(vali,axv,nbytes);
    }
  delete[]new_row;
  delete[]aux_val;
  delete[]aux_col;
  bas.reduce_mod_p(modulus);
#if TRACE_ELIM
  cout<<"Finished constructing basis for kernel"<<endl;
  cout<<"Basis = "<<bas.as_mat()<<endl;
#endif
  return bas;
}

void smat_elim::step0()
{
  /*This step eliminates all rows with zero or only one entry, 
   *  system is supposed to be homogeneous */ 

  list L(nro);
  int row,i,j,n;
  for( row = 0; row < nro; row++ )
    if( *col[row] < 2 ) L.put( row );

  while( (row = L.next()) != -1 ) { 
   if( *col[row] == 0 ) { position[ row ] = 0; continue; }
    else {               // only one entry in that row
      val[row][0] = 1;   // trivial normalization

      /* clear other rows with entry in that column */
      
      int colr = col[row][1];
      int N = (column + colr - 1)->num;   // # of rows in column col;
      for( j = 0; j < N; j++ ) {
	i = (column + colr - 1)->next();   // row to be cleared of col;
	if( i == row ) continue;
	int d = col[i][0]--;
	if( d == 2 ) L.put( i );
	int ind = find( colr, col[i]+1, d-1 );
	int *pos = col[i] + ind + 1;
	if( *pos != colr ) { cerr << "error in step0!"<<endl; return;}
	scalar *values = val[i] + ind;
	for( n = ind+1; n < d; n++, pos++, values++ ) 
	  { *pos = pos[1]; *values = values[1]; }
      }
      
      eliminate( row, colr );
      free_space( colr );
    }
  }
}

void smat_elim::step1 ()
{
  /* eliminates all rows which cut a column which has only one entry */
  
  list L(nco);
  int col0,col1;
#if TRACE_ELIM
  cout<<"Step 1, column weights:"<<endl;  
  //  for( col0 = 0; col0 < nco; col0++ ) cout<<(column+col0)->num<<" ";
  //  cout<<endl;
#endif
  for( col0 = 0; col0 < nco; col0++ )
    if( (column+col0)->num == 1 ) {col1=col0+1; L.put(col1);}
#if TRACE_ELIM
  cout<<"Step 1, list size = "<<L.num<<endl;  
#endif
  while( (col0 = L.next()) != -1 ) {
    if( (column+col0-1)->num < 1 ) continue;
    int row = (column+col0-1)->next();
    normalize( row, col0 );
    
    /* update column */
    int *pos = col[row];
    int d = *pos++;
    while( d-- ) {
      int c = *pos++ - 1;
      (column + c)->remove(row);
      if((column + c)->num == 1) {col1=c+1; L.put( col1 );
      //    cout<<"List size increases to "<<L.num<<endl;
      }
    }
    eliminate( row, col0 );
  }
}

void smat_elim::step2()
{
  /*  eliminates all rows with 1 or 2 entries  */
  
  list L(nro);
  int row;
  for( row = 0; row < nro; row++ )
    if( *col[row] < 3 && position[row] == -1 ) L.put( row );
  
   while( (row = L.next()) != -1 ) {
    if( position[row] != -1 ) continue;
    int colr = col[row][1];
    normalize( row, colr );
    clear_col ( row, colr, L, 1 );      
    eliminate( row, colr );
    free_space( colr );
  }
}

void smat_elim::step3()
{
  /* eliminates all rows which cut a column which have either one or two
     entries */

  list L(nco);
  int col0,col1;
  //  for( col0 = 0; col0 < nco; col0++ ) {
  for( col0 = nco-1; col0 >=0; col0-- ) {
    int vali = (column+col0)->num;
    if( vali == 2 || vali == 1 ) {col1=col0+1; L.put(col1);}
  }

  while( (col0 = L.next()) != -1 ) {
    if( (column+col0-1)->num < 1 ) continue;
    int row = (column+col0-1)->next();
    normalize( row, col0 );
    clear_col( row, col0, L, 0, 1 );
    eliminate( row, col0 );
    free_space( col0 );
  }
}

void smat_elim::step4 ( )
{
  int* lightness = new int[nco];
  int M, i, wt, r, row;

  // Find maximum column weight
  int maxcolwt=0;
  for( i = 0; i < nco; i++ ) 
    {
      wt = (column+i)->num;
      if( maxcolwt < wt) maxcolwt=wt;
    }
  int M0 = maxcolwt; // max(20, int(maxcolwt/10)); // 20;
  int Mstep = int(maxcolwt/100);
  if (Mstep==0) Mstep=1;
#if TRACE_ELIM
  cout<<"Step 4, max column weight = "<<maxcolwt<<endl;  
#endif

  //float Mscale = 0.9;
  //for( M = M0; M >= 4; M--)
  //for( M = M0; M >= 3; M*=Mscale)
  for( M = M0; M >= 3; M-=Mstep)
    {
#if TRACE_ELIM
  cout<<"Step 4, M = "<<M;
#endif
      /* divides columns in `light' and `heavy' */
      int nlight=0;
      int *l = lightness;
      for( i = 0; i < nco; i++ ) {
	wt = (column+i)->num;
	if( 0 < wt && wt <= M ) {*l++ = 1; nlight++;} //light
	else *l++ = 0;     //heavy; includes columns already eliminated
      }
#if TRACE_ELIM
      cout<<", "<<nlight<<" light columns; ";  
      report();
#endif
  if (nlight==0) break; // from the loop over M
  if (nlight<(nco/2)) break; // from the loop over M
  //if (nlight<=(nco/4)) break; // from the loop over M
      while(1)
	{
	  /* eliminates rows with weight 1 */
	  for( r = 0, row = -1; r < nro; r++ ) {
	    if(has_weight_one(r, lightness) && position[r] == -1) 
	      { row = r; break; }
	  }
	  if( row != -1 ) 
	    {
	      int col0 = 0;       // light col cutting row
	      int d = *col[row]; // weight in the process of eliminating row r.
	      int *pos = col[row] + 1;
	      while( d-- ) {
		int c = *pos++ - 1;
		if(lightness[c] == 1) { col0 = c+1; break; }
	      }
	      if( col0 == 0 ) {cerr << "step4: row doesn't cut light col"<<endl; return;}
	      normalize( row, col0 );
	      list temp(0);
	      clear_col(row,col0, temp, 0, 0, M, lightness);
	      eliminate( row, col0 );
	      free_space( col0 );
	    }
	  else break;
	}       
    }
            
  delete [] lightness;
}

void smat_elim::standard ( ){

  // remaining elimination

  int i, col0, row, wt, mincolwt;
  double density_threshold = 0.2; // 0.25; // 0.5; // 0.75; // 0.9; // 0.7; // 0.5; // 0.2; // 0.04; // 0.15; // 0.2; // 0.1;

 // this threshold can be changed: the code here works fine when the
 // density is low, but when it is higher it's best to switch to using
 // a dense structure for the remaining elimination.  It might also be
 // better not to recompute the density after every single step.


#if TRACE_ELIM
  report();
  cout << "Continuing elimination in sparse structure until density goes above " << density_threshold <<endl;
#endif
  while(active_density() < density_threshold)
    {
  // Find minimum positive column weight
      mincolwt=nro+1; col0=-1;
      for( i = 0; i < nco; i++ ) 
        {
          wt = (column+i)->num;
          if( (wt>0) && (mincolwt > wt) ) {col0=i+1; mincolwt=wt;}
        }
      if (col0==-1) return;  // this is how the while(1) is ended
#if TRACE_ELIM
      //      cout<<"... wt "<<mincolwt<<flush;
#endif
      row = (column+col0-1)->next();
      normalize( row, col0 ); 
      list temp(0);
      clear_col( row, col0, temp );
      eliminate( row, col0 );
      free_space( col0 );
    }
#if TRACE_ELIM
  cout << "Finished elimination in sparse structure, density is now " << active_density() <<endl;
  report();
#endif
  step5dense();
}

void smat_elim::back_sub ( ){

  /* Back substitution */
  for( int n = rank; n; n-- )
    {   
      int row = elim_row[n-1];
      int* pos = col[row] + 1;
      for( int j = 0; j < *col[row]; j++ )
	{
	  int e  = elim_col[*pos++-1];
	  if( e != -1 && e != row )
	    {
	      elim( e, row, -val[row][j] );
	      j = -1;
	      pos = col[row] + 1;
	    }
	}
    }
}

void smat_elim::normalize( int row, int col0)
{
  int d = *col[row];
  int count = find( col0, col[row]+1, d-1 );
  if( col[row][count+1] != col0 ) 
    { cerr << "error in normalize "<<endl; return; }
  if( val[row][count] != 1 ) {
    scalar invValue = invmod( val[row][count], modulus);
    scalar *values = val[row];
    while(d--) { *values = xmm( *values , invValue, modulus ); values++; }
  }
}

void smat_elim::eliminate( const int& row, const int& col0 ) //1<=col0<=nco;
{
  //cout<<"Eliminating (r,c)=("<<row<<","<<col0-1<<")"<<endl;
  elim_col[ col0 - 1 ] = row;
  position[ row ] = col0;
  elim_row[ rank++ ] = row;
}

void 
smat_elim::clear_col( int row,int col0,list& L, int fr, int fc,int M,int* li )
{
  int numRow = (column+col0-1)->num;
  int d = col[row][0];
  int *pos1 = col[row]+1;
  if( numRow == 1 ) {
    for( int s = 0; s < d; s++ ) {
      int c = pos1[s] - 1;
      (column + c)->remove(row);
      if( fc ) check_col( c, L );         // check condition for cols
      if( M ) {
	int l = (column+c)->num;
	if( 0 < l && l <= M ) li[c] = 1;  // col is light
	else li[c] = 0;  // heavy;
      }
    }
    return;
  }
  list::listsize = numRow;

  /* for the d cols in col[row], these lists will contain rows to be taken 
   * in/out of column */
  list *list_row_out = new list [d]; 
  list *list_row_in = new list [d];
  if( !list_row_out ) {cerr << "memory exhausted in elim::clear_col"<<endl; return;};
  if( !list_row_in ) {cerr << "memory exhauted in elim::clear_col"; return;};
  list* lri = list_row_in, *lro = list_row_out;  
  /* eliminate col from other rows cutting col */

  int di = d;
  scalar *veci1 = val[row];
  (column+col0-1)->index = 0;   //reset index for iteration
  for( int l = 0; l < numRow; l++ ) {
    int row2 = (column+col0-1)->next();
    if( row2 == row ) continue;
    int *pos2 = col[row2];
    int d2 = *pos2++;
    int ind = find(col0, pos2, d2-1);
    if( pos2[ind] != col0 ) { cerr << "error in clear_col"<<endl; return; }
    int d2i = d2;
    scalar *oldVal = val[row2]; int *oldMat = col[row2];
    scalar *veci2 = oldVal;
    scalar v2 = mod(modulus-veci2[ind],modulus);
    int *P = col[row2] = new int [ d + d2 + 1 ]; P++;
    scalar *V = val[row2] = new scalar [ d + d2 ];

    /* do row2+= v2*row1 */
    int k = 0;       /*k will be # of non-zero entries of sum*/
    while( d && d2 )
      { 
	if( *pos1 < *pos2 ) { 
	  lri[di-d].put(row2);
	  *P++ = *pos1++; *V++ = xmm( v2,(*veci1++), modulus ); d--; 
	}
	else if(( *P++ = *pos2++ ) < *pos1 ) { *V++ = *veci2++; d2--; }
	else
	  {
	    if( (*V++ = addmod(xmm(v2,(*veci1++),modulus) , (*veci2++), modulus)) == 0)
	      { lro[di-d].put(row2); V--; P--; k--;}
	    pos1++;
	    d--;
	    d2--;
	  }
	k++;  
      }
    if( d == 0 ) while( d2 )
      { *P++ = *pos2++; *V++ = *veci2++; k++; d2--; }
    else while( d ) {
      lri[di-d].put(row2);
      *P++ = *pos1++; *V++ = xmm(v2,(*veci1++), modulus); k++; d--;
    }
    *col[row2] = k;
    if( fr ) check_row(d2i, row2, L);         // check condition for rows
    delete [] oldMat;
    delete [] oldVal;
    d = di;         // reset d, pos1 and veci1
    pos1 -= d;
    veci1 -= d;
  }

  /* update column */
  for( int t = 0; t < di; t++ ) {
    int c = col[row][t+1]-1;
    (column+c)->remove(row);
    column[c].remove(list_row_out[t]);
    (column+c)->put(list_row_in[t]);
    if( fc ) check_col( c, L );            // check condition for cols
    if( M ) {
      int l = (column+c)->num;
      if( 0 < l && l <= M ) li[c] = 1;  // col is light
      else li[c] = 0;  // heavy;
    }
  }
  
  delete [] list_row_out;
  delete [] list_row_in;
}

void smat_elim::free_space( int col0 )
{
  (column+col0-1)->clear();
}

void smat_elim::check_row (int d, int row2, list& L ) 
{
   if( *col[row2] < 3 ) {
      if( *col[row2] == 0 ) position[row2] = 0;
      //if d <= 2 then row2 was already in the list, so
      else if( d > 2 ) L.put(row2);  
  }
}
void smat_elim::check_col( int c, list& L ) 
{
  int vali = (column+c)->num;
  if( vali == 2 || vali == 1 ) {L.put(c+1);}
}

int smat_elim::get_weight( int row, const int* lightness ) 
{
  int wt = 0;
  int *pos = col[row];
  int d = *pos++;
  while( d-- ) wt += lightness[ *pos++ - 1 ];
  return wt;
}

int smat_elim::has_weight_one( int row, const int* lightness )
{
  int wt = 0;
  int *pos = col[row];
  int d = *pos++;
  while( d-- )
    {
      wt += lightness[ *pos++ - 1 ];
      if (wt>1) return 0;
    }
  return (wt==1);
}

int smat_elim::n_active_cols() // number of active columns
{
  // Remaining cols are those with positive column weight
  int j, nrc;
  for(j=nrc=0; j<nco; j++) 
    if (((column+j)->num)>0) 
      nrc++;
  return nrc;
}

int smat_elim::n_active_rows() // number of active rows
{
  // Remaining rows are those with "position" code -1 or those which are empty
  int i, nrr;
  for(i=nrr=0; i<nro; i++)
    if( (*col[i] >0) && (position[i] == -1) )
      nrr++;
  return nrr;
}

long smat_elim::n_active_entries() // number of active entries
{
  // Remaining cols are those with positive column weight
  int j; long n=0;
  for(j=0; j<nco; j++)
    n += ((column+j)->num); 
  return n;
}

double smat_elim::active_density() // density of non-eliminated part
{
  double d = n_active_entries();
  int n = n_active_cols();
  if (!n) return 0;
  d /= n;
  n = n_active_rows();
  if (!n) return 0;
  d /= n;
  return d;
}

void smat_elim::report()
{
  cerr << n_active_entries() << " active entries in ("
       << n_active_rows() << "," << n_active_cols()
       << ") active (rows, cols).  Active density = "
       << active_density() << endl;
  cerr<<"Rank so far = "<<rank<<endl;
}

/* old fashioned elim function for back elimination. Have to change this 
 * later.
 * Do row2+= v2*row1 */
void smat_elim::elim( int row1, int row2, scalar v2 )
{
  int d = *col[row1], d2 = *col[row2];
  scalar *oldVal = val[row2]; int *oldMat = col[row2];
  int *pos1 = col[row1]+1, *pos2 = oldMat + 1;
  scalar *veci1 = val[row1], *veci2 = oldVal;
  int *P = col[row2] = new int [ d + d2 + 1 ]; P++;
  scalar *V = val[row2] = new scalar [ d + d2 ];
  int k = 0;       /*k will be # of non-zero entries of sum*/
  while( d && d2 )
    { 
      if( *pos1 < *pos2 ) 
	{*P++ = *pos1++; *V++ = xmm( v2,(*veci1++),modulus ); d--; }
      else if(( *P++ = *pos2++ ) < *pos1 ) { *V++ = *veci2++; d2--; }
      else
	{
	  if( (*V++ = addmod(xmm(v2,(*veci1++),modulus) , (*veci2++), modulus)) == 0)
	    { V--; P--; k--;}
	  pos1++; // unused, but prevents compiler warning
	  d--;
	  d2--;
	}
      k++;  
    }
  if( d == 0 ) while( d2 )
    { *P++ = *pos2++; *V++ = *veci2++; k++; d2--; }
  else if( d2 == 0 ) while( d )
                       { *P++ = *pos1++; *V++ = xmm(v2,(*veci1++),modulus); k++; d--; }
  *col[row2] = k;
  delete [] oldMat;
  delete [] oldVal;
}

void smat_elim::step5dense()
{
#if TRACE_DENSE
  report();
  cerr<<"switching to dense elimination"<<endl;
#endif
  // (1) Extract the uneliminated "dense" part

  vector<int> remaining_rows, remaining_cols;

  // Remaining rows are those with "position" code -1 or those which are empty
  int i, j;
  for(i=0; i<nro; i++)
    if( (*col[i] >0) && (position[i] == -1) )
      remaining_rows.push_back(i+1);
  int nrr = remaining_rows.size();

  // Remaining cols are those with positive column weight
  for(j=0; j<nco; j++) 
    if (((column+j)->num)>0) 
      remaining_cols.push_back(j+1);
  int nrc = remaining_cols.size();
#if TRACE_DENSE
  cout<<nrr<<" remaining rows, " <<nrc<<" remaining cols"<<endl;
  //    cout<<" remaining rows: " << remaining_rows<<endl;
  //    cout<<" remaining cols: " << remaining_cols<<endl;
#endif
  if(nrr*nrc==0) //(nrr*nrc<10000) // avoid overheads of switching to dense if
                    // there's not a lot left to do
    {
      standard();
      return;
    }
  mat dmat(nrr, nrc);
  for (i=0; i<nrr; i++)
    {
      svec v = row(remaining_rows[i]);
      j=0;
      for( const auto& vi : v.entries)
        {
          while (remaining_cols[j]<(vi.first)) j++;
          dmat.set(i+1,j+1,vi.second);
        }
    }

  // (2) reduce this to echelon form

#if TRACE_DENSE
  cout<<"Constructed dense matrix, starting dense elimination step..." <<endl;
#endif
  vec_i pc,npc; long rk,ny;

  // this will call ref_via_ntl() if FLINT, else ref_via_flint()
  dmat = rref(dmat,pc,npc,rk,ny,modulus);

#if TRACE_DENSE
    cout<<"...finished dense elimination, rank = "<<rk;
    cout<<", nullity = "<<ny<<endl;
    //cout<<"Pivotal columns:    "<<pc<<endl;
    // cout<<"Nonpivotal columns: "<<npc<<endl;
#endif

  // (3) put it back into the sparse structure

  // the (i,j) entry of dmat goes in the remaining_rows[i-1]'th row,
  // remaining_cols[j-1] column.  For simplicity of coding, we create
  // the new rows as svecs and the use setrow().
  int nrd = dmat.nrows(); // may be less than nrr since 0 rows are trimmed
  svec rowi(nco);
  for(i=1; i<=nrd; i++)
    {
      rowi.clear();
      for(j=1; j<=nrc; j++)
        rowi.set(remaining_cols[j-1],dmat(i,j));
      setrow(remaining_rows[i-1],rowi);
    }
  rowi.clear();
  for(i=nrd+1; i<=nrr; i++)
    setrow(remaining_rows[i-1],rowi);

  // (4) Use the known echelon form for these changed rows to eliminate them

#if TRACE_DENSE
    cout<<"remaining elimination within sparse structure"<<endl;
#endif
  for(i=1; i<=nrd; i++) 
    {
      if (is_nonzero(xmod(dmat(i,pc[i])-1,modulus)))
        cout<<"Bad pivot #"<<i<<" ("<<dmat(i,pc[i])<<")"<<endl;
      int r = remaining_rows[i-1]-1;
      int c = remaining_cols[pc[i]-1];
      //      cout<<"Eliminating (r,c)=("<<r<<","<<c<<")"<<endl;
      eliminate(r,c);
      free_space(remaining_cols[pc[i]-1]);
    }
#if TRACE_DENSE
    cout<<"finished dense step"<<endl;
#endif
}

long smat::rank(scalar mod)
{
  smat_elim sme(*this,mod);
  (void) sme.sparse_elimination();
  return sme.get_rank();
}

ssubspace::ssubspace(int n)
  :pivots(vec_i::iota(n)), basis(smat::identity_matrix(n))
{}

ssubspace::ssubspace(const smat& b, const vec_i& p, scalar mod)
  :modulus(mod),pivots(p),basis(b)
{}

ssubspace::ssubspace(const ssubspace& s)
  :modulus(s.modulus),pivots(s.pivots),basis(s.basis)
{}

// assignment
void ssubspace::operator=(const ssubspace& s) 
{
  pivots=s.pivots;
  basis=s.basis;
  modulus=s.modulus;
}

// Definitions of nonmember, nonfriend operators and functions:

ssubspace combine(const ssubspace& s1, const ssubspace& s2)
{
  scalar m = s1.modulus;
  return ssubspace(mult_mod_p(s1.basis,s2.basis,m),s1.pivots[s2.pivots],m);
}
 
smat restrict_mat(const smat& m, const ssubspace& s)
{ 
  return mult_mod_p(m.select_rows(pivots(s)),basis(s),s.modulus);
}
 
ssubspace kernel(const smat& sm, scalar m)
{
  vec_i pivs, npivs;
  smat kern = smat_elim(sm,m).kernel(npivs,pivs);
  return ssubspace(kern,pivs,m);
}
 
ssubspace eigenspace(const smat& sm, scalar lambda, scalar m)
{
  smat m1 = sm; m1.sub_mod_p(lambda, m);
  return kernel(m1, m);
}
 
ssubspace subeigenspace(const smat& sm, scalar l, const ssubspace& s, scalar m)
{
  return combine(s,eigenspace(restrict_mat(sm,s), l, m));
}

ssubspace make1d(const vec& bas, scalar&piv, scalar m)
// make a 1-D ssubspace with basis bas
{
  smat tbasis(1,dim(bas));
  svec sbas(bas);
  tbasis.setrow(1,sbas);
  vec_i pivs(1); // initialised to 0
  pivs[1]=sbas.first_index();
  piv=sbas.elem(pivs[1]);
  return ssubspace(transpose(tbasis),pivs,m);
}
