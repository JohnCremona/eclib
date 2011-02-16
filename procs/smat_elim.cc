// smat_elim.cc: implementation of class smat_elim
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
 
//  implements structured modular elimination
//  original written by Luiz Figueiredo

inline scalar xmm(scalar a, scalar b, scalar m)
{
  //return xmodmul(a,b,m);
  //return (a*b) % m;
  return (a*(int64_t)b) % m;
  //return (scalar)(((long)a*(long)b) % (long)m);
}
inline scalar xmm0(scalar a, scalar b)
{
  //return xmodmul(a,b,m);
  //return (a*b) % m;
  return (a*(int64_t)b) % BIGPRIME;
  //return (scalar)(((long)a*(long)b) % (long)m);
}



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
  if( !new_array ) { cerr << "memory exhausted"; abort(); }      
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
smat_elim::list::find( type& X, int ub, int lb )
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
      abort(); 
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
      abort(); 
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
	abort(); 
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
  if( !column ) { cerr << "memory exhausted"; abort(); }
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

//#define TRACE_ELIM

void smat_elim::sparse_elimination( )
{
#ifdef TRACE_ELIM
  int pop=get_population(*this);
  double density = (double)pop/(nco*nro);
  cout<<"Starting sparse elimination: "<<nro<<" rows, "<<nco<<" columns, "<<pop<<" entries (density = "<<density<<")\n";
  //  cout<<"row weights:\n";
  //  for(int i=0; i<nro; i++) cout<<col[i][0]<<" ";
  //  cout<<endl;
  //  cout<<(*this)<<endl;
  cout<<"Starting step 0..."<<flush;
#endif
  step0();
#ifdef TRACE_ELIM
  cout<<"finished\n"; 
  //  cout<<(*this)<<endl;
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
  cout<<"Starting step 1..."<<flush;
#endif
  step1();
#ifdef TRACE_ELIM
  cout<<"finished\n";
  //  cout<<(*this)<<endl;
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
  cout<<"Starting step 2..."<<flush;
#endif
  step2();
#ifdef TRACE_ELIM
  cout<<"finished\n";
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
  cout<<"Starting step 3..."<<flush;
#endif
  step3();
#ifdef TRACE_ELIM
  cout<<"finished\n";
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
  cout<<"Starting step 4..."<<flush;
#endif
  step4();
#ifdef TRACE_ELIM
  cout<<"finished\n";
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
  cout<<"Starting step 4..."<<flush;
#endif
  step4();
#ifdef TRACE_ELIM
  cout<<"finished\n";
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
#endif
#if(0)  // use dense method for final elimination
#ifdef TRACE_ELIM
  cout << "Switching to dense mode..."<<endl;
#endif
  step5dense();
#ifdef TRACE_ELIM
  cout<<"finished"<<endl;
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
#endif
#else  // use sparse method for final elimination
#ifdef TRACE_ELIM
    cout<<"Starting step 5 (remaining elimination)..."<<flush;
#endif
     standard( );
#ifdef TRACE_ELIM
  cout<<"finished"<<endl;
  pop=get_population(*this);
  density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
#endif
#endif
}

smat smat_elim::kernel( vec& pc, vec& npc)
{
  int i,n,r,denom = 1;
#ifdef TRACE_ELIM
  cout<<"Starting sparse_elimination()..."<<flush;
#endif
  sparse_elimination( );
#ifdef TRACE_ELIM
  cout<<"finished sparse_elimination()"<<endl;
#endif


  int nullity = nco - rank;
  if (nullity>0)
    {
#ifdef TRACE_ELIM
  cout<<"Starting back-substitution..."<<flush;
#endif
      back_sub();
#ifdef TRACE_ELIM
  cout<<"finished back-substitution"<<endl;
#endif

    }
  smat basis( nco, nullity );
  pc.init( rank );
  npc.init( nullity );

  /* set-up vecs pc & npc */
#ifdef TRACE_ELIM
  cout<<"Setting up pc and npc..."<<flush;
#endif
  int ny = 0, k = 0;
  long *new_row = new long [ rank ];
  for( i = 1; i <= nco; i++ )
    {
      if( elim_col[i-1] > -1 ) { pc[++k] = i; new_row[k-1] = elim_col[i-1]; }
      else npc[++ny] = i;
    }

  /* write basis for kernel */
#ifdef TRACE_ELIM
  cout<<"Constructing basis for kernel..."<<flush;
#endif
  for( n = 1; n <= nullity; n++ )
    { 
      int i = npc[n]-1;
      basis.col[i][0] = 1;      //this much storage was granted in the 
      basis.col[i][1] = n;      // in the constructor.
      basis.val[i][0] = denom;
    }

  scalar *aux_val = new scalar [nco];
  int *aux_col = new int [nco];
  for ( r=1; r<=rank; r++)
    { 
      int i = pc[r]-1;
      int count = 0;
      int *axp = aux_col; scalar *axv = aux_val;
      int *posB = col[new_row[r-1]];
      int d = *posB++-1;
      scalar *valB = val[new_row[r-1]];
      for (int j = 1, h = 0; j<=nullity; j++) {
	while( *posB < npc[j] && h < d ) { posB++; h++; }
	if( *posB == npc[j] ) {	*axp++ = j; *axv++ = -valB[h]; count++; }
      }
      delete [] basis.col[i];
      delete [] basis.val[i];
      basis.col[i] = new int [count + 1];
      basis.val[i] = new scalar [count];
      int *pos = basis.col[i];
      scalar *val = basis.val[i];
      axp = aux_col;
      axv = aux_val;
      *pos++ = count;
      for( n = 0; n < count; n++ ) { *pos++ = *axp++; *val++ = *axv++; }
    }
  delete[]new_row;
  delete[]aux_val;
  delete[]aux_col;
#ifdef TRACE_ELIM
  cout<<"Finished constructing basis for kernel"<<endl;
  //  cout<<"Basis = "<<basis<<endl;
#endif
  return basis;
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
	if( *pos != colr ) { cerr << "error in step0!\n"; abort();}
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
#ifdef TRACE_ELIM
  cout<<"Step 1, column weights:"<<endl;  
  //  for( col0 = 0; col0 < nco; col0++ ) cout<<(column+col0)->num<<" ";
  //  cout<<endl;
#endif
  for( col0 = 0; col0 < nco; col0++ )
    if( (column+col0)->num == 1 ) {col1=col0+1; L.put(col1);}
#ifdef TRACE_ELIM
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
    scalar val = (column+col0)->num;
    if( val == 2 || val == 1 ) {col1=col0+1; L.put(col1);}
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
  int Mstep = int(maxcolwt/100);
  if (Mstep==0) Mstep=1;
#ifdef TRACE_ELIM
  cout<<"Step 4, max column weight = "<<maxcolwt<<endl;  
#endif
  
  //  for( M = 20; M >= 4; M--)
  for( M = maxcolwt; M >= 3; M-=Mstep)
    {
#ifdef TRACE_ELIM
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
#ifdef TRACE_ELIM
      cout<<", "<<nlight<<" light columns; ";  
  int pop=get_population(*this);
  double density = (double)pop/(nco*nro);
  cout<<"Rank so far = "<<rank<<"; "
      <<pop<<" entries (density = "<<density<<")\n";
#endif
  if (nlight==0) break; // from the loop over M
  if (nlight<(nco/2)) break; // from the loop over M
      while(1)
	{
	  /* eliminates rows with weight 1 */
	  for( r = 0, row = -1; r < nro; r++ ) {
	    if(get_weight(r, lightness) == 1 && position[r] == -1) 
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
	      if( col0 == 0 ) {cerr << "step4: row doesn't cut light col"; abort();}
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
  while(1)
    {
  // Find minimum positive column weight
      mincolwt=nro+1; col0=-1;
      for( i = 0; i < nco; i++ ) 
        {
          wt = (column+i)->num;
          if( (wt>0) && (mincolwt > wt) ) {col0=i+1; mincolwt=wt;}
        }
      if (col0==-1) break;  // this is how the while(1) is ended
#ifdef TRACE_ELIM
      //      cout<<"... wt "<<mincolwt<<flush;
#endif
      row = (column+col0-1)->next();
      normalize( row, col0 ); 
      list temp(0);
      clear_col( row, col0, temp );
      eliminate( row, col0 );
      free_space( col0 );
    }

  // for( int col0 = 1; col0 <= nco; col0++ )
  //   {
  //     if( (column+col0-1)->num > 0 ) {
  //       int row = (column+col0-1)->next();
  //       normalize( row, col0 ); 
  //       list temp(0);
  //       clear_col( row, col0, temp );
  //       eliminate( row, col0 );
  //       free_space( col0 );
  //     }
  //   }
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
    { cerr << "error in normalize "; abort(); }
  if( val[row][count] != 1 ) {
    scalar invValue = invmod0( val[row][count]);
    scalar *values = val[row];
    while(d--) { *values = xmm0( *values , invValue ); values++; }
  }
}

void smat_elim::eliminate( int& row, int& col0 ) //1<=col0<=nco;
{
  //  cout<<"Eliminating (r,c)=("<<row<<","<<col0<<")"<<endl;
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
  if( !list_row_out ) {cerr << "memory exhausted"; abort();};
  if( !list_row_in ) {cerr << "memory exhauted"; abort();};
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
    if( pos2[ind] != col0 ) { cerr << "error in clear_col"; abort(); }
    int d2i = d2;
    scalar *oldVal = val[row2]; int *oldMat = col[row2];
    scalar *veci2 = oldVal;
    scalar v2 = -1*veci2[ind];
    int *P = col[row2] = new int [ d + d2 + 1 ]; P++;
    scalar *V = val[row2] = new scalar [ d + d2 ];

    /* do row2+= v2*row1 */
    int k = 0;       /*k will be # of non-zero entries of sum*/
    while( d && d2 )
      { 
	if( *pos1 < *pos2 ) { 
	  lri[di-d].put(row2);
	  *P++ = *pos1++; *V++ = xmm0( v2,(*veci1++) ); d--; 
	}
	else if(( *P++ = *pos2++ ) < *pos1 ) { *V++ = *veci2++; d2--; }
	else
	  {
	    if( (*V++ = xmod0(xmm0(v2,(*veci1++)) + (*veci2++))) == 0)
	      { lro[di-d].put(row2); V--; P--; k--;}
	    *pos1++;
	    d--;
	    d2--;
	  }
	k++;  
      }
    if( d == 0 ) while( d2 )
      { *P++ = *pos2++; *V++ = *veci2++; k++; d2--; }
    else while( d ) {
      lri[di-d].put(row2);
      *P++ = *pos1++; *V++ = xmm0(v2,(*veci1++)); k++; d--;
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
  int c1, val = (column+c)->num;
  if( val == 2 || val == 1 ) {c1=c+1; L.put(c1);}
}

int smat_elim::get_weight( int row, int* lightness ) 
{
  int wt = 0;
  int *pos = col[row];
  int d = *pos++;
  while( d-- ) wt += lightness[ *pos++ - 1 ];
  return wt;
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
	{*P++ = *pos1++; *V++ = xmm0( v2,(*veci1++) ); d--; }
      else if(( *P++ = *pos2++ ) < *pos1 ) { *V++ = *veci2++; d2--; }
      else
	{
	  if( (*V++ = xmod0(xmm0(v2,(*veci1++)) + (*veci2++))) == 0)
	    { V--; P--; k--;}
	  *pos1++;
	  d--;
	  d2--;
	}
      k++;  
    }
  if( d == 0 ) while( d2 )
    { *P++ = *pos2++; *V++ = *veci2++; k++; d2--; }
  else if( d2 == 0 ) while( d )
    { *P++ = *pos1++; *V++ = xmm0(v2,(*veci1++)); k++; d--; }
  *col[row2] = k;
  delete [] oldMat;
  delete [] oldVal;
}

void smat_elim::step5dense()
{
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

  cout<<nrr<<" remaining rows, " <<nrc<<" remaining cols"<<endl;
  //  cout<<" remaining rows: " << remaining_rows<<endl;
  //  cout<<" remaining cols: " << remaining_cols<<endl;
  
  mat dmat(nrr, nrc);
  map<int,scalar>::const_iterator vi;
  vector<int>::const_iterator rci;
  for (i=0; i<nrr; i++)
    {
      svec v = row(remaining_rows[i]);
      j=0;
      for(vi=v.entries.begin();
          vi!=v.entries.end(); vi++)
        {
          while (remaining_cols[j]<(vi->first)) j++;
          dmat.set(i+1,j+1,vi->second);
        }
    }

  // (2) reduce this to echeclon form

  vec pc,npc; long rk,ny;
  dmat = echmodp_uptri(dmat,pc,npc,rk,ny,BIGPRIME);
  // cout<<"Rank = "<<rk<<endl;
  // cout<<"Nullity = "<<ny<<endl;
  // cout<<"Pivotal columns:    "<<pc<<endl;
  // cout<<"Nonpivotal columns: "<<npc<<endl;
  
  // (3) put it back into the sparse structure

  // the (i,j) entry of dmat goes in the remaining_rows[i-1]'the rowm
  // remaining_cols[j-1] column.  For simplicity of coding, we create
  // the new rows as svecs and the use setrow().
  int nrd = nrows(dmat); // may be less that nrr since 0 rows are trimmed
  int ncd = ncols(dmat);
  svec rowi(nco);
  for(i=1; i<=nrd; i++)
    {
      rowi.clear();
      for(int j=1; j<=nrc; j++)
        rowi.set(remaining_cols[j-1],dmat(i,j));
      setrow(remaining_rows[i-1],rowi);
    }
  rowi.clear();
  for(i=nrd+1; i<=nrr; i++)
    setrow(remaining_rows[i-1],rowi);

  // (4) Use the known echelon form for these changed rows to eliminate them

  for(i=1; i<=nrd; i++) 
    {
      if (xmod0(dmat(i,pc[i])-1)) cout<<"Bad pivot #"<<i<<" ("<<dmat(i,pc[i])<<")"<<endl;
      int r = remaining_rows[i-1]-1;
      int c = remaining_cols[pc[i]-1];
      //      cout<<"Eliminating (r,c)=("<<r<<","<<c<<")"<<endl;
      eliminate(r,c);
      free_space(remaining_cols[pc[i]-1]);
    }
}


long rank(smat& sm)
{
  smat_elim sme(sm);
  vec pivs, npivs;
  (void) sme.kernel(npivs,pivs);
  return sme.get_rank();
}

ssubspace::ssubspace(int n) 
:pivots(iota((scalar)n)),basis(sidmat((scalar)n))
{}

ssubspace::ssubspace(const smat& b, const vec& p)
:pivots(p),basis(b)
{}

ssubspace::ssubspace(const ssubspace& s) 
:pivots(s.pivots),basis(s.basis) 
{}

// destructor -- no need to do anything as componenets have their own
ssubspace::~ssubspace() 
{}

// assignment
void ssubspace::operator=(const ssubspace& s) 
{
  pivots=s.pivots; 
  basis=s.basis; 
}

// Definitions of nonmember, nonfriend operators and functions:

ssubspace combine(const ssubspace& s1, const ssubspace& s2)
{ 
  return ssubspace(mult_mod_p(s1.basis,s2.basis,BIGPRIME),s1.pivots[s2.pivots]);
}
 
smat restrict(const smat& m, const ssubspace& s)
{ 
  return mult_mod_p(m.select_rows(pivots(s)),basis(s),BIGPRIME);
}
 
ssubspace kernel(const smat& sm)
{
  vec pivs, npivs;
  smat kern = smat_elim(sm).kernel(npivs,pivs);
  return ssubspace(kern,pivs);
}
 
ssubspace eigenspace(const smat& m1, scalar lambda)
{
  smat m = m1; m.sub_mod_p(lambda);
  return kernel(m);
}
 
ssubspace subeigenspace(const smat& m1, scalar l, const ssubspace& s)
{
  return combine(s,eigenspace(restrict(m1,s), l));
}


ssubspace make1d(const vec& bas, long&piv) 
// make a 1-D ssubspace with basis bas
{
  smat tbasis(1,dim(bas));
  svec sbas(bas);
  tbasis.setrow(1,sbas);
  vec pivs(1); // initialised to 0
  pivs[1]=sbas.first_index();
  piv=sbas.elem(pivs[1]);
  return ssubspace(transpose(tbasis),pivs);
}
