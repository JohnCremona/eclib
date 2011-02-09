// smat.cc: implementation of sparse integer matrix class smat
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
 
// ONLY to be included by smatrix.cc

void showrow(int*pos, scalar*val) // for debugging
{
  int d=pos[0]; 
  cout<<"[";
  int* posi=pos; posi++;
  scalar* vali=(scalar*)val;
  while(d--) cout<<"("<<(*posi++)<<","<<(*vali++)<<")";
  /*
  if(sizeof(scalar)==sizeof(int)) 
    {
      int* vali=(int*)val;
      while(d--) cout<<"("<<(*posi++)<<","<<(*vali++)<<")";
    }
  else 
    {
      long* vali=(long*)val;
      while(d--) cout<<"("<<(*posi++)<<","<<(*vali++)<<")";
    }
  */
  cout<<"]";
}

//#define DEBUG_MEM

// Definitions of member operators and functions:

smat::smat(int nr, int nc)
{
  nco = nc;
  nro = nr;
  col = new int * [nr];
  val = new scalar * [nr];
#ifdef DEBUG_MEM
  cout<<"Constructed an smat with (nr,nc)=("<<nr<<","<<nc<<"), with col="<<col<<", val="<<val<<endl;
#endif
  for( int i = 0; i < nr; i++ )
    {
      col[i] = new int [ 2 ];
      val[i] = new scalar [ 1 ];
      col[i][1] = col[i][0] = val[i][0] = 0;
    }
}

smat::smat(const smat& sm)
{
  nco = sm.nco;
  nro = sm.nro;
  col = new int * [nro];
  val = new scalar * [nro];
#ifdef DEBUG_MEM
  cout<<"Constructed an smat (copy constructor) with col="<<col<<", val="<<val<<endl;
#endif
  for( int i = 0; i < nro; i++ )
    { 
      int d = sm.col[i][0];
      col[i] = new int[ d+1 ];  
      val[i] = new scalar[ d ];  
      int *pos = col[i], *pi = sm.col[i];
      scalar *values = val[i], *vi = sm.val[i];
      *pos++ = *pi++;
      while (d--) { *values++ = *vi++; *pos++ = *pi++; }
    }
}

smat::smat(const mat& m)
{
  //  cout<<"Converting mat("<<m.nro<<"x"<<m.nco<<") to smat"<<endl;
  nco = m.nco;
  nro = m.nro;
  col = new int * [nro];
  val = new scalar * [nro];
#ifdef DEBUG_MEM
  cout<<"Constructed an smat (from a mat) with col="<<col<<", val="<<val<<endl;
#endif
  int i, j, k, l, p;
  for( i = 0; i < nro; i++ )
    {
      scalar *veci = m.entries + i*nco;
      for( j = 0, k = 0; j < nco; j++ ) if( *veci++ ) k++;
      col[i] = new int[ k+1 ];  
      val[i] = new scalar[ k ];  
      scalar *values = val[i]; int *pos = col[i];
      veci = m.entries + i*nco;
      *pos++ = k;
      for( l = 0, p = 1;  l < nco; l++, p++,veci++ ) 
	if( *veci ) { *values++ = *veci; *pos++ = p; }
    }
}

smat::~smat()
{
  for( int i = 0; i < nro; i++ ) { delete [] col[i]; delete [] val[i]; }
#ifdef DEBUG_MEM
  cout<<"Destroying an smat with col="<<col<<", val="<<val<<endl;
#endif
  delete [] col;
  delete [] val;
}

// member functions and operators

void smat::set_row( int i, int d, int* pos, scalar* values)
{
  if( col[i][0] != d ) {
    delete [] col[i]; delete [] val[i]; 
    col[i] = new int [d+1];
    val[i] = new scalar [d];
    col[i][0] = d;
  }
  for( int j = 0; j < d; j++ ) {
    col[i][j+1] = *pos++;
    val[i][j] = *values++;
  }
}

void smat::setrow ( int i, const svec& v) // i counts from 1
{
  int j, d=v.entries.size();
  i--;
  if( col[i][0] != d ) {
    delete [] col[i]; delete [] val[i]; 
    col[i] = new int [d+1];
    val[i] = new scalar [d];
    col[i][0] = d;
  }
  map<int,scalar>::const_iterator vi;
  for(vi=v.entries.begin(), j=0;
      vi!=v.entries.end(); vi++, j++)
    {
      col[i][j+1] = vi->first;
      val[i][j] = vi->second;
    }
}

smat smat::select_rows(const vec& rows) const
  {
    int i,r, n=dim(rows);
    smat ans(n,nco);
    for(i=0; i<n; i++)
      {
	r=rows[i+1]-1;
	ans.set_row(i,col[r][0],col[r]+1,val[r]);
      }
    return ans;
  }

mat smat::as_mat( ) const
{
  //  cout<<"Converting smat to mat("<<nro<<"x"<<nco<<")"<<endl;
  mat ans( nro, nco ); 
  scalar *mi = ans.entries;
  for( int i = 0; i < nro; i++ )
    {
      int d = *col[i];
      scalar *values = val[i];
      int *posi = col[i] + 1;
      while( d-- )
	mi[ i*nco + (*posi++) - 1 ] = *values++;
    }
  return ans;
}

svec smat::row(int i) const // extract row i as an svec, i counts from 1
{
  i--;
  svec ans(nco);
  int d = *col[i];
  scalar *values = val[i];
  int *posi = col[i] + 1;
  while( d-- )
    ans.set( (*posi++),  (*values++));
  return ans;
}

scalar smat::elem(int i, int j)  const   /*returns (i,j) entry, 1 <= i <= nro
        				  * can only be used as a rvalue  */
{
 if( (0<i) && (i<=nro) && (0<j) && (j<=nco) )           
   {
     int d = *col[i-1];
     int *posi = col[i-1] + 1;
     scalar *veci = val[i-1]; 
     while( d-- )
       { 
	 if( j == *posi++ ) return *veci;
	 veci++;
       }
     return 0;
   }
 else 
   {
     cerr << "Bad indices in smat::operator ()\n"; 
     return 0;
   }
}


smat& smat::operator=(const smat& sm)
{
 if (this==&sm) return *this;
 nco = sm.nco;
 int i, n = sm.nro; 
 if (nro != n) // delete old space and replace with new.
   {            
     for( i = 0; i < nro; i++ ) { delete [] col[i]; delete [] val[i]; }
#ifdef DEBUG_MEM
  cout<<"in smat.=, destroying old smat with col="<<col<<", val="<<val<<endl;
#endif
     delete [] col;
     delete [] val;
     nro = n;
     col = new int * [nro]; 
     val = new scalar * [nro]; 
#ifdef DEBUG_MEM
  cout<<"in smat.=, creating new smat with col="<<col<<", val="<<val<<endl;
#endif
  for( i = 0; i < nro; i++ )
    {
      col[i] = new int [ 2 ];
      val[i] = new scalar [ 1 ];
      col[i][1] = col[i][0] = val[i][0] = 0;
    }
   }
 for( i = 0; i < nro; i++ )
   { 
     int d = *sm.col[i];
     if(d!=col[i][0])
       {
	 delete[]col[i]; delete[]val[i];
	 col[i] = new int [ d+1 ];
	 val[i] = new scalar [ d ];
	 col[i][0]=d;
       }
     scalar *values = val[i]; int *pos = col[i];
     scalar *vi = sm.val[i]; int *pi = sm.col[i];
     *pos++ = *pi++;
     while (d--) { *values++ = *vi++; *pos++ = *pi++; }
   }
  return *this;
}

smat& smat::operator+=(const smat& mat2)
{
  if ((nro==mat2.nro))
    {
      for(int i = 0; i < nro; i++ )
	{
// 	  cout<<"Adding rows ";
// 	  showrow(col[i],val[i]); cout<<" and "; showrow(mat2.col[i],mat2.val[i]);
// 	  cout<<endl;
	  int d = *col[i], d2 = *mat2.col[i];
	  int *pos1 = col[i] + 1, *pos2 = mat2.col[i]+1;
 	  scalar *val1 = val[i], *val2 = mat2.val[i];
	  int *P = new int [ d + d2 + 1 ]; int* Pi=P+1;
	  scalar *V = new scalar [ d + d2 ]; scalar* Vi=V;
	  int k = 0;       /*k will be # of non-zero entries of sum*/
	  while( d && d2 )
	    { 
	      if( *pos1 < *pos2 ) 
		{  *Pi++ = *pos1++; *Vi++ = *val1++; d--; k++; }
	      else if(*pos2 < *pos1 )
		{  *Pi++ = *pos2++; *Vi++ = *val2++; d2--; k++; }
	      else {
// 		cout<<"two entries at position "<<(*pos1)<<", "<<(*val1)<<" and "<<(*val2)<<endl;
		*Pi = *pos1;
		scalar newval = (*val1++) + (*val2++);
// 		cout<<"sum = "<<newval<<endl;
		if( newval!=0 ) 
		  {
		    *Vi=newval;
// 		    cout<<"nonzero, putting "<<(*V)<<" into sum row in position "<<(*Pi)<<endl; 
		    Vi++; Pi++; k++;
		  }
// 		else cout<<"zero, omitting from sum row"<<endl;
		pos1++; pos2++; d--; d2--;
	      }
	    }
// 	  P[0] = k;
// 	  cout<<"intermediate result is ";
// 	  showrow(P,V);
// 	  cout<<endl;
	  while( d2-- )
	    { *Pi++ = *pos2++; *Vi++ = *val2++; k++; }
	  while( d -- )
	    { *Pi++ = *pos1++; *Vi++ = *val1++; k++; }
	  P[0] = k;
	  delete [] col[i]; col[i]=P;
	  delete [] val[i]; val[i]=V;
// 	  cout<<"Result is ";
// 	  showrow(col[i],val[i]);
// 	  cout<<endl;
	}
    }
  else cerr << "Incompatible smatrices in operator +=\n";
  return *this;
}

smat& smat::operator-=(const smat& mat2)
{
  if ((nro==mat2.nro))
    {
      for(int i = 0; i < nro; i++ ) 
	{
	  int d = *col[i], d2 = *mat2.col[i];
	  int *pos1 = col[i] + 1, *pos2 = mat2.col[i]+1;
 	  scalar *val1 = val[i], *val2 = mat2.val[i];
	  int *P = new int [ d + d2 + 1 ]; int* Pi=P+1;
	  scalar *V = new scalar [ d + d2 ]; scalar* Vi=V;
	  int k = 0;       /*k will be # of non-zero entries of sum*/
	  while( d && d2 )
	    { 
	      if( *pos1 < *pos2 ) 
		{ *Pi++ = *pos1++; *Vi++ = *val1++; d--; k++; }
	      else if(*pos2 < *pos1 )
		{  *Pi++ = *pos2++; *Vi++ = -(*val2++); d2--; k++; }
	      else {
		*Pi = *pos1;
		scalar newval = (*val1++) - (*val2++);
		if( newval!=0 ) { *Vi++=newval; Pi++; k++; }
		pos1++; pos2++; d--; d2--;
	      }
	    }
	  while( d2-- )
	    { *Pi++ = *pos2++; *Vi++ = -(*val2++); k++; }
	  while( d -- )
	    { *Pi++ = *pos1++; *Vi++ = *val1++; k++; }
	  P[0] = k;
	  delete [] col[i]; col[i]=P;
	  delete [] val[i]; val[i]=V;
	}
    }
  else cerr << "Incompatible matrices in operator -=\n";
  return *this;
}

smat& smat::operator+= (const scalar& scal) // adds scalar*identity
{
  if(scal==0) return *this;
  int i, d, k;
  for(i = 0; i < nro; i++ )
    {
      d = *col[i];                      // length of old row
      int *pos1 = col[i] + 1;           // pointer to run along position vector
      scalar *val1 = val[i];            // pointer to run along value vector
      int *P = new int [ d + 2 ];       //  new position vector
      scalar *V = new scalar [ d + 1 ]; //  new value vector
      int* Pi=P+1;
      scalar* Vi=V;
      scalar newval;
      k = 0;           // k will be # of non-zero entries of new row
      while((d)&&(*pos1<(i+1)))  // just copy entries
	{ 
	  *Pi++ = *pos1++; *Vi++ = *val1++; k++; d--;
	}
      if(d&&(*pos1==(i+1)))     // add the scalar, see if it's zero
	{
	  newval = (*val1)+scal;
	  if( newval!=0) { *Vi++ = newval; *Pi++=*pos1; k++; }
	  pos1++; val1++; d--;
	}
      else // insert new entry
	{
	  *Vi++ = scal; *Pi++=(i+1); k++; 
	}
      while(d--) // copy remaining entries if necessary
	{
	  *Pi++ = *pos1++; *Vi++ = *val1++; k++; 
	}
      P[0] = k;
      delete [] col[i]; col[i]=P;
      delete [] val[i]; val[i]=V;
    }
  return *this;
}

void smat::sub_mod_p(const scalar& lambda, const scalar& p) 
// subtracts scalar*identity mod p
{
  this->operator-=(lambda);
  this->reduce_mod_p(p);
}

void smat::reduce_mod_p(const scalar& p)
{
  svec rowi;
  for(int i=1; i<=nro; i++)
    {
      rowi = row(i);
      rowi.reduce_mod_p(p);
      setrow(i,rowi);
    }
}

smat& smat::operator*=(scalar scal)
{
  if(scal==0) cerr<<"Attempt to multiply smat by 0\n"<<endl;
  int i, d; scalar *veci;
  for( i = 0; i < nro; i++)
    {
      d = *col[i];
      veci = val[i];
      while(d--) (*veci++) *= scal;
    }
  return *this;
}

smat& smat::mult_by_scalar_mod_p (scalar scal, const scalar& p)
{
  if(xmod(scal,p)==0) cerr<<"Attempt to multiply smat by 0\n"<<endl;
  int i, d; scalar *veci;
  for( i = 0; i < nro; i++)
    {
      d = *col[i];
      veci = val[i];
      while(d--) {(*veci) = xmodmul(*veci,scal,p); veci++;}
    }
  return *this;
}

smat& smat::operator/=(scalar scal)
{
  if(scal==0) cerr<<"Attempt to divide smat by 0\n"<<endl;
  int i, d; scalar *veci;
   for(  i = 0; i < nro; i++)
    {
      d = *col[i];
      veci = val[i];
      while(d--) (*veci++) /= scal;
    }
  return *this;
 }

mat smat::operator*( const mat& m )
{
  if( nco != nrows(m) )
    {
      cerr << "incompatible smat & mat in operator*\n";
      abort();
    }
  mat product( nro, ncols(m) );
  int i, j, d, t;
  scalar ans;
  for( i = 1; i <= nro; i++ ) 
    {
      d = col[i-1][0];
      for( j = 1; j <= ncols(m); j++ ) 
	{
	  ans = 0;
	  for( t = 0; t < d; t++ ) ans += val[i-1][t]*m(col[i-1][t+1],j);
	  product(i,j) = ans;
	}
    }
  return product;
}


// Definitions of non-member, friend operators and functions

svec operator* ( const smat& A, const svec& v )
{
  if( A.nco != dim(v) ) 
    { 
      cout << "incompatible smat*svec\n"; 
      cout << "Dimensions "<<dim(A)<<" and "<<dim(v)<<endl;
      abort();
    }
  int n = A.nro, j; scalar s;
  svec prod(n);
  for(j = 1; j<=n; j++)  
    {
      s = (A.row(j))*v;
      if(s) prod.entries[j]=s;
    }
  return prod;
}

vec operator*  (smat& m, const vec& v)
{
  int r = nrows(m), c=ncols(m);
  if(c!=dim(v))
    {
      cout<<"Error in smat*vec:  wrong dimensions ("<<r<<"x"<<c<<")*"<<dim(v)<<endl;
      abort();
    }
  vec w(r);
  for(int i=1; i<=r; i++) w.set(i,m.row(i)*v);
  return w;
}

// (col) svec * smat

svec operator* ( const svec& v, const smat& A )
{
  if( v.d != nrows(A) ) 
    { 
      cout << "incompatible sizes in v*A\n"; 
      cout << "Dimensions "<<v.d<<" and "<<dim(A)<<endl;
      abort();
    }
  svec prod(ncols(A));
  map<int,scalar>::const_iterator vi;
  for(vi=v.entries.begin(); vi!=v.entries.end(); vi++)
    prod += (vi->second)*(A.row(vi->first));
  return prod;
}

svec mult_mod_p( const svec& v, const smat& A, const scalar& p  )
{
  if( v.d != nrows(A) ) 
    { 
      cout << "incompatible sizes in v*A\n"; 
      cout << "Dimensions "<<v.d<<" and "<<dim(A)<<endl;
      abort();
    }
  svec prod(ncols(A));
  map<int,scalar>::const_iterator vi;
  for(vi=v.entries.begin(); vi!=v.entries.end(); vi++)
    prod.add_scalar_times_mod_p(A.row(vi->first), vi->second,p);
  return prod;
}

#if(1)
smat operator* ( const smat& A, const smat& B )
{
  if( A.nco != B.nro ) { cerr << "incompatible smats in operator *\n"; abort();}
  int nro = A.nro, nco = B.nco;
  smat prod( nro, nco );

  for (int i=1; i<=nro; i++)
    prod.setrow(i, A.row(i)*B);
  return prod;
}

smat mult_mod_p ( const smat& A, const smat& B, const scalar& p )
{
  if( A.nco != B.nro ) { cerr << "incompatible smats in operator *\n"; abort();}
  int nro = A.nro, nco = B.nco;
  smat prod( nro, nco );

  for (int i=1; i<=nro; i++)
    prod.setrow(i, mult_mod_p(A.row(i),B,p));
  return prod;
}
#endif

#if(0)
smat operator* ( const smat& A, const smat& B )
{
  if( A.nco != B.nro ) { cerr << "incompatible smats in operator *\n"; abort();}
  int nro = A.nro, nco = B.nco;
  smat prod( nro, nco );
  
  /* writes columns of B */
  
  int* ncol = new int [nco];
  int *nc = ncol;
  int l, r, s;
  for( l = 0; l < nco; l++ ) *nc++ = 0;
  
  for( r = 0; r <B.nro; r++ ) // counts # of elements in each col
    {
      int d = *B.col[r];
      int *p = B.col[r] + 1;
      while( d-- ) ncol[*p++ - 1]++;
    }
  scalar **colB_val = new scalar * [nco];
  int **colB_mat = new int * [nco];
  for( l = 0; l < nco; l++ ) {
    colB_val[l] = new scalar [ncol[l]];
    colB_mat[l] = new int [ncol[l]];
  }
  int * aux = new int [nco];
  int *a=aux;
  for( s = 0; s < nco; s++ ) *a++ = 0;

  for( r = 0; r < B.nro; r++ ) {
    int d = *B.col[r];
    scalar *v = B.val[r];
    int *p = B.col[r] + 1;
    while( d-- ) {
      int col = *p++ - 1;
      colB_mat[col][aux[col]] = r+1;
      colB_val[col][aux[col]++] = *v++ ;
    }
  }

  /* multiply A and B */
 
  int * aux_pos = new int [nco];
  
  for( r = 0; r < nro; r++ ) {
    int count = 0;
    int d = *A.col[r];
    int *ax = aux, *axp = aux_pos; //reusing aux. It will now hold values of *
    for( int m = 0; m < nco; m++ ) { *ax++ = 0; *axp++ = 0; }
    int *posA = A.col[r] +1;
    for( int i, c = 0; c < nco; c++ ) { 
      scalar soma = 0;
      int *posB = colB_mat[c];
      for( l = 0, i = 0; l < ncol[c] && i < d; ) {
	if( posA[i] < posB[l] ) i++;
	else if( posA[i] >  posB[l] ) l++;
	else soma= xmod0(soma + xmodmul0(A.val[r][i++] , colB_val[c][l++]));
      }
      if( soma != 0 ) {
	aux_pos[count] = c+1;
	aux[count++] = mod0(soma);
      }
    }
    delete [] prod.col[r];
    delete [] prod.val[r];
    int *pos = prod.col[r] = new int [count+1];
    scalar *val = prod.val[r] = new scalar [count];
    *pos++ = count;
    for( ax = aux, axp = aux_pos, l = 0; l < count; l++ ) {
      *pos++ = *axp++;
      *val++ = *ax++;
    }
  }

  delete [] aux;
  delete [] aux_pos;
  delete [] ncol;
  for( l = 0; l < nco; l++ ) { delete [] colB_mat[l]; delete [] colB_val[l]; }
  delete [] colB_mat; delete [] colB_val;

  return prod;
}

#endif


smat transpose ( const smat& A )
{
  // 1. Count the number of entries in each column (as in operator*() below):
  int *colwts = new int[A.nco];
  int i, r;
  for(i=0; i<A.nco; i++) colwts[i]=0;
  for( r = 0; r <A.nro; r++ ) // counts # of elements in each col
    {
      int d = *A.col[r];
      int *p = A.col[r] + 1;
      while( d-- ) colwts[*p++ - 1]++;
    }
#if(0)
  cout<<"Column weights of A:\n";
  for(i=0; i<A.nco; i++) cout<<colwts[i]<<" ";
  cout<<endl;
#endif
  // 2. Make space for the output matrix:
  smat B(A.nco,A.nro);
  // Remove the default entries in B:
  for( int i = 0; i < B.nro; i++ ) { delete [] B.col[i]; delete [] B.val[i]; }
  // Replace with the correct sizes:
  for( i = 0; i < B.nro; i++ )
    {
      int d = colwts[i];
      B.col[i] = new int[ d+1 ];  
      B.val[i] = new scalar[ d ];  
      B.col[i][0] = d;
    }
  delete[]colwts;
  //3. Copy entries over.  aux[i] holds the number of entries so far
  //   put into row i of the transpose
  int * aux = new int [B.nro];
  int *a=aux;
  for( r = 0; r < B.nro; r++ ) *a++ = 0;
  for( r = 0; r < A.nro; r++ ) {
    int d = *A.col[r];
    //    cout<<"row "<<r<<" of A has "<<d<<" entries\n";
    scalar *v = A.val[r];
    int *p = A.col[r] + 1;
    while( d-- ) {
      int c = *p++ - 1;
      B.col[c][aux[c]+1] = r+1;
      B.val[c][aux[c]] = *v++;
      aux[c]++;
    }
#if(0)
    cout<<"After processing that row, aux = \n";
    for(i=0; i<A.nco; i++) cout<<aux[i]<<" ";
    cout<<endl;
#endif
  }
  delete[]aux;
  return B;
}

int operator==(const smat& sm1, const smat& sm2)
{
   int nr = sm1.nro, i;
   int equal = ( nr == sm2.nro );
   if(!equal) return 0;

   for( i = 0; i < nr && equal; i++ )
     {
       int d1 = *sm1.col[i], d2 = *sm2.col[i];
       if( d1 != d2 ) {return 0;}
       
       scalar *sm1val = sm1.val[i], *sm2val= sm2.val[i];
       int *sm1pos= sm1.col[i] + 1;
       int *sm2pos = sm2.col[i] + 1;
       while (equal && d1--) equal = ((*sm1val++)==(*sm2val++));
       while (equal && d2--) equal = ((*sm1pos++)==(*sm2pos++));
     }
   return equal;
}

int eqmodp(const smat& sm1, const smat& sm2, const scalar& p)
{
   int nr = sm1.nro, i;
   int equal = ( nr == sm2.nro );
   if(!equal) return 0;

   for( i = 0; i < nr && equal; i++ )
     {
       int d1 = *sm1.col[i], d2 = *sm2.col[i];
       if( d1!=d2 ) {return 0;}
       
       scalar *sm1val = sm1.val[i], *sm2val= sm2.val[i];
       int *sm1pos= sm1.col[i] + 1;
       int *sm2pos = sm2.col[i] + 1;
       while (equal && d2--) equal = ((*sm1pos++)==(*sm2pos++));
       while (equal && d1--) equal = (xmod((*sm1val++)-(*sm2val++),p)==0);
     }
   return equal;
}

ostream& operator << (ostream& s, const smat& sm)
{
  for( int i = 0; i < sm.nro; i++ )
    {
      cout << "row[" << i+1 << "] =";
      int d = *sm.col[i]; 
      int *posi = sm.col[i] + 1; scalar *veci = sm.val[i];
      int n = d-1 > 0 ? d-1 : 0;
      s << "{ ";
      s << "values " << "[";
      if( d > 0 ) s << *veci++; 
      while ( n-- ) s << "," << (*veci++); 
      s << "]";
      s << "   positions: " << "[";
      if( d > 0 ) { s << *posi++; d = d-1; }
      while ( d-- ) { s << "," << (*posi++); }
      s << "]    }" << endl;
    }
  return s;
}


istream& operator >> (istream& s, smat& sm)
{
  int *pos = new int [ sm.nco ];
  scalar *values = new scalar [ sm.nco ];
  int r, k, count;
  for( r = 0; r < sm.nro; r++ )
    { 
      cout << "input row " << r+1 << endl;
      int *p = pos; scalar *v = values;
      s >> k;
      for( count = 0; k != 0; s >> k )
	{
	  *v++ = k;
	  s >> k;
	  if( k ) *p++ = k;
	  else { cerr << "enter zero as a value!!!\n"; abort(); }
	  count++;
	}
      
      delete [] sm.col[r];
      delete [] sm.val[r];
      sm.col[r] = new int [ count + 1 ];
      sm.val[r] = new scalar [ count ];
      
      sm.col[r][0] = count;
      p = pos;
      v = values;
      for( k = 0; k < count; k++ ) 
	{ sm.col[r][k+1] = *p++; sm.val[r][k] = *v++; }
    }
  delete [] pos;
  delete [] values;
  return s;
}


    
// Definition of non-friend functions

smat operator+(const smat& sm)
{
        return sm;
}

smat operator-(const smat& sm)
{
        return (-1)*sm;
}

smat operator+(const smat& sm1, const smat& sm2)
{
  smat ans(sm1); ans+=sm2;  return ans;
}

smat operator-(const smat& sm1, const smat& sm2) 
{
  smat ans(sm1); ans-=sm2;  return ans;
}

smat operator*(scalar scal, const smat& sm)
{
  smat ans(sm); ans*=scal;  return ans;
}

smat operator/(const smat& sm, scalar scal)
{
  smat ans(sm); ans/=scal;  return ans;
}

int operator!=(const smat& sm1, const smat& sm2)
{
        return !(sm1==sm2);
}

int get_population(const smat& m )
{
  int r,d,count = 0;
  for( r = 0; r < m.nro; r++ ) 
    {
      d = *(m.col[r]);
      if(d==0) continue;
      int *pos = m.col[r] + 1;
      while( d-- ) { count += ( *pos++ != 0 );}
    }
  return count;
}

smat sidmat(scalar n)  // identity matrix
{
  smat I(n,n); // creates enough space
  for( int i = 0; i < n; i++ )
    {
      I.col[i][0] = 1;   // one entry in this row
      I.col[i][1] = i+1; // ...it's in column i+1
      I.val[i][0] = 1;   // ...its value is 1
    }
  return I;
}

smat liftmat(const smat& mm, scalar pr, scalar& dd, int trace)
{
  scalar modulus=pr,n,d; long nr,nc; dd=1;
  int succ=0,success=1;
  float lim=floor(sqrt(pr/2.0));
  smat m = mm;
  if(trace)
    {
      cout << "Lifting mod-p smat;  smat mod "<<pr<<" is:\n";
      cout << m.as_mat();
      cout << "Now lifting back to Q.\n";
      cout << "lim = " << lim << "\n";
    }
  for(nr=0; nr<m.nro; nr++)
    for(nc=0; nc<m.col[nr][0]; nc++)
      {
	succ = modrat(m.val[nr][nc],modulus,lim,n,d);
	success = success && succ;
	dd=lcm(d,dd);
      }
  if(!success) 
    cout << "Problems encountered with modrat lifting of smat." << endl;
  dd=abs(dd);
  if(trace) cout << "Common denominator = " << dd << "\n";
  for(nr=0; nr<m.nro; nr++)
    for(nc=0; nc<m.col[nr][0]; nc++)
      m.val[nr][nc] = mod(xmodmul(dd,(m.val[nr][nc]),pr),pr);
  if(trace) cout << "Lifted smat = " << m.as_mat() << "\n";
  return m;
}
