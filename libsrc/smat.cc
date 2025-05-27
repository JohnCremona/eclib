// smat.cc: implementation of sparse integer matrix class smat
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
 
// ONLY to be included by smatrix.cc

void showrow(int*pos, scalar*val) // for debugging
{
  int d=pos[0]; 
  cout<<"[";
  int* posi=pos; posi++;
  scalar* vali=(scalar*)val;
  while(d--) cout<<"("<<(*posi++)<<","<<(*vali++)<<")";
  cout<<"]";
}

//#define DEBUG_MEM

// Definitions of member operators and functions:

smat::smat(int nr, int nc)
{
  static const scalar zero(0);
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
      col[i][1] = col[i][0] = 0;
      val[i][0] = zero;
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
      auto veci = m.entries.begin() + i*nco;
      for( j = 0, k = 0; j < nco; j++ ) if( is_nonzero(*veci++) ) k++;
      col[i] = new int[ k+1 ];
      val[i] = new scalar[ k ];
      scalar *values = val[i]; int *pos = col[i];
      veci = m.entries.begin() + i*nco;
      *pos++ = k;
      for( l = 0, p = 1;  l < nco; l++, p++,veci++ )
	if( is_nonzero(*veci) ) { *values++ = *veci; *pos++ = p; }
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

void smat::set_row( int i, int d, int* pos, scalar* values) // i counts from 0
{
  scalar *vali = val[i];
  int *coli = col[i];
  // allocate space for d entries, but only set the non-zero ones
  if( coli[0] != d ) {
    delete [] coli; delete [] vali;
    coli = col[i] = new int [d+1];
    vali = val[i] = new scalar [d];
  }
  coli++;
  int j=d;
  while(j--) {
    scalar v = *values++;
    int c = *pos++;
    if (is_nonzero(v))
      {
        *coli++ = c;
        *vali++ = v;
      }
  }
  col[i][0] = (coli - col[i]) - 1; // difference of pointers
}

void smat::setrow ( int i, const vec& v) // i counts from 1
{
  int j, m, n, dimv=dim(v);
  scalar *vali;
  int *coli;

  // count nonzero entries of v:
  auto vi = v.entries.begin();
  j = dimv;
  n = 0;
  while(j--)
    if (is_nonzero(*vi++))
      n++;

  // (re)allocate position and value arrays
  i--; // so it starts from 0
  coli = col[i];
  vali = val[i];
  if( coli[0] != n ) {
    delete [] coli; delete [] vali;
    coli = col[i] = new int [n+1];
    vali = val[i] = new scalar [n];
    *coli++ = n;
  }

  // copy nonzero entries
  coli++;
  vi = v.entries.begin();
  for(m=1; m<=dimv; m++)
    {
      scalar e = *vi++;
      if (is_nonzero(e))
        {
          *coli++ = m;
          *vali++ = e;
        }
    }
}

void smat::setrow ( int i, const svec& v) // i counts from 1
{
  int d=v.entries.size();
  i--;
  scalar *vali = val[i];
  int *coli = col[i];
  if( coli[0] != d ) {
    delete [] coli;
    delete [] vali;
    coli = col[i] = new int [d+1];
    vali = val[i] = new scalar [d];
    coli[0] = d;
  }

  coli++;
  for( const auto& vi : v.entries)
    {
      *coli++ = vi.first;
      *vali++ = vi.second;
    }

  // check
  // if (row(i+1) != v)
  //   {
  //     cerr << "error in smat::setrow(int, svec):";
  //     cerr << "v = "<<v.as_vec()<<endl;
  //     cerr << "row was set to "<< row(i+1)<<endl;
  //   }
}

smat smat::select_rows(const vec_i& rows) const
{
  int n=dim(rows);
  smat ans(n,nco);
  for(int i=0; i<n; i++)
    {
      int r=rows[i+1]-1;
      ans.set_row(i,col[r][0],col[r]+1,val[r]);
    }
  return ans;
}

mat smat::as_mat( ) const
{
  //  cout<<"Converting smat to mat("<<nro<<"x"<<nco<<")"<<endl;
  mat ans( nro, nco );
  auto mi = ans.entries.begin();
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
  static const scalar zero(0);
  if( (0<i) && (i<=nro) && (0<j) && (j<=nco) )
   {
     int d = *col[i-1];
     if (d==0) return zero;
     int *first = col[i-1]+1;
     int *last = col[i-1]+1+d;
     int *p = std::lower_bound(first, last, j);
     if ((p==last) || (*p!=j))
       return zero;
     return val[i-1][p-first];
   }
 else
   {
     cerr << "Bad indices ("<<i<<","<<j<<") in smat::operator ()! (nro,nco)=("<<nro<<","<<nco<<")\n";
     return zero;
   }
}


smat& smat::operator=(const smat& sm)
{
 static const scalar zero(0);
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
      col[i][1] = col[i][0] = 0;
      val[i][0] = zero;
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
  if (nro==mat2.nro)
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
  if (nro==mat2.nro)
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

  for(int i = 0; i < nro; i++ )
    {
      int d = *col[i];                      // length of old row
      int *pos1 = col[i] + 1;           // pointer to run along position vector
      scalar *val1 = val[i];            // pointer to run along value vector
      int *P = new int [ d + 2 ];       //  new position vector
      scalar *V = new scalar [ d + 1 ]; //  new value vector
      int* Pi=P+1;
      scalar* Vi=V;
      int k = 0;           // k will be # of non-zero entries of new row
      while((d)&&(*pos1<(i+1)))  // just copy entries
	{
	  *Pi++ = *pos1++; *Vi++ = *val1++; k++; d--;
	}
      if(d&&(*pos1==(i+1)))     // add the scalar, see if it's zero
	{
	  scalar newval = (*val1)+scal;
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
  for( int i = 0; i < nro; i++)
    {
      int d = *col[i];
      scalar *veci = val[i];
      while(d--) (*veci++) *= scal;
    }
  return *this;
}

smat& smat::mult_by_scalar_mod_p (scalar scal, const scalar& p)
{
  if(xmod(scal,p)==0) cerr<<"Attempt to multiply smat by 0\n"<<endl;
  for(int i = 0; i < nro; i++)
    {
      int d = *col[i];
      scalar *veci = val[i];
      while(d--) {(*veci) = xmodmul(*veci,scal,p); veci++;}
    }
  return *this;
}

smat& smat::operator/=(scalar scal)
{
  if(scal==0) cerr<<"Attempt to divide smat by 0\n"<<endl;
  for(int i = 0; i < nro; i++)
    {
      int d = *col[i];
      scalar *veci = val[i];
      while(d--) (*veci++) /= scal;
    }
  return *this;
 }

mat smat::operator*( const mat& m )
{
  if( nco != m.nrows() )
    {
      cerr << "incompatible smat & mat in operator*"<<endl;
      return mat();
    }
  mat product( nro, m.ncols() );
  scalar ans;
  for(int i = 1; i <= nro; i++ )
    {
      int d = col[i-1][0];
      for(int j = 1; j <= m.ncols(); j++ ) 
	{
	  ans = 0;
	  for(int t = 0; t < d; t++ ) ans += val[i-1][t]*m(col[i-1][t+1],j);
	  product(i,j) = ans;
	}
    }
  return product;
}

long smat::nullity(const scalar& lambda, scalar m) // nullity of this-lambda*I modulo m
{
  smat sma(*this); sma-=lambda;  return sma.ncols()-sma.rank(m);
}

// Definitions of non-member, friend operators and functions

svec operator* ( const smat& A, const svec& v )
{
  if( A.nco != dim(v) ) 
    { 
      cerr << "incompatible smat*svec\n"; 
      cerr << "Dimensions "<<dim(A)<<" and "<<dim(v)<<endl;
      return svec();
    }
  int n = A.nro;
  svec prod(n);
  for(int j = 1; j<=n; j++)
    {
      scalar s = (A.row(j))*v;
      if(is_nonzero(s)) prod.entries[j]=s;
    }
  return prod;
}

vec operator* (const smat& m, const vec& v)
{
  int r = m.nrows(), c=m.ncols();
  vec w(r);
  if(c!=dim(v))
    cerr<<"Error in smat*vec:  wrong dimensions ("<<r<<"x"<<c<<")*"<<dim(v)<<endl;
  else
    for(int i=1; i<=r; i++) w.set(i,m.row(i)*v);
  return w;
}

// (col) svec * smat

svec operator* ( const svec& v, const smat& A )
{
  svec prod(A.ncols());
  if( v.d != A.nrows() )
    {
      cerr << "incompatible sizes in v*A\n"
           << "Dimensions "<<v.d<<" and "<<dim(A)<<endl;
    }
  else
    {
      std::for_each(v.entries.cbegin(), v.entries.cend(),
                    [&prod, A](auto vi){prod += (vi.second)*(A.row(vi.first));});
    }
  return prod;
}

svec mult_mod_p( const svec& v, const smat& A, const scalar& p  )
{
  vec prod(A.ncols());
  if( v.d != A.nrows() )
    {
      cerr << "incompatible sizes in v*A\n"
           << "Dimensions "<<v.d<<" and "<<dim(A)<<endl;
    }
  else
    {
      for( const auto& vi : v.entries)
        {
          // prod.add_scalar_times_mod_p(A.row(vi->first), vi->second,p);
          int i = (vi.first)-1;         // the row of A to use (from 0)
          scalar c = vi.second;         // the coefficient to multiply it by
          int d = A.col[i][0];          // #nonzero entries in this row
          int *posi = A.col[i] +1;      // pointer to array of columns
          scalar *values = A.val[i];    // pointer to array of values
          while (d--)
            prod.add_modp(*posi++,xmodmul(c,*values++,p),p);
        }
    }
  return svec(prod);
}


svec mult_mod_p( const smat& A, const svec& v, const scalar& p  )
{
  svec w(A.nrows());
  if( v.d != A.ncols() )
    {
      cerr << "incompatible sizes in A*v\n";
      cerr << "Dimensions "<<dim(A)<<" and "<<v.d<<endl;
    }
  else
    {
      for(int i=1; i<=A.nrows(); i++)
        w.set(i,dotmodp(A.row(i),v,p));
    }
  return w;
}

vec mult_mod_p( const smat& A, const vec& v, const scalar& p  )
{
  vec w(A.nrows());
  if( dim(v) != A.ncols() )
    {
      cerr << "incompatible sizes in A*v\n";
      cerr << "Dimensions "<<dim(A)<<" and "<<dim(v)<<endl;
    }
  else
    {
      for(int i=1; i<=A.nrows(); i++)
        w.set(i,dotmodp(A.row(i),v,p));
    }
  return w;
}

smat operator* ( const smat& A, const smat& B )
{
  int nro = A.nro, nco = B.nco;
  smat prod( nro, nco );
  if( A.nco != B.nro )
    {
      cerr << "incompatible smats in operator *"<<endl;
    }
  else
    {
      for (int i=1; i<=nro; i++)
        prod.setrow(i, A.row(i)*B);
    }
  return prod;
}

smat mult_mod_p ( const smat& A, const smat& B, const scalar& p )
{
  int nro = A.nro, nco = B.nco;
  smat prod( nro, nco );
  if( A.nco != B.nro )
    {
      cerr << "incompatible smats in operator *"<<endl;
    }
  else
    {
      for (int i=1; i<=nro; i++)
        prod.setrow(i, mult_mod_p(A.row(i),B,p));
    }
  return prod;
}

smat transpose ( const smat& A )
{
  // 1. Count the number of entries in each column (as in operator*() below):
  int *colwts = new int[A.nco];
  for( int i=0; i<A.nco; i++) colwts[i]=0;
  for( int r=0; r <A.nro; r++ ) // counts # of elements in each col
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
  for( int i = 0; i < B.nro; i++ )
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
  for( int r = 0; r < B.nro; r++ ) *a++ = 0;
  for( int r = 0; r < A.nro; r++ ) {
    int d = *A.col[r];
    //    cout<<"row "<<r<<" of A has "<<d<<" entries\n";
    const scalar *v = A.val[r];
    int *p = A.col[r] + 1;
    while( d-- ) {
      int c = *p++ - 1;
      B.col[c][aux[c]+1] = r+1;
      B.val[c][aux[c]] = *v++;
      aux[c]++;
    }
#if(0)
    cout<<"After processing that row, aux = \n";
    for(int i=0; i<A.nco; i++) cout<<aux[i]<<" ";
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
	  else { cerr << "invalid entry value 0 in smat input"<<endl; }
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
  return scalar(-1)*sm;
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
  int count = 0;
  for(int r = 0; r < m.nro; r++ )
    {
      int d = *(m.col[r]);
      if(d==0) continue;
      int *pos = m.col[r] + 1;
      while( d-- ) { count += ( is_nonzero(*pos++) );}
    }
  return count;
}

scalar maxabs(const smat& m )
{
  scalar a(0);
  for(int r = 0; r < m.nro; r++ )
    {
      int d = *(m.col[r]);
      const scalar *values = m.val[r];
      int *pos = m.col[r] + 1;
      while( d-- ) { a = max(a, abs(*values++)); pos++;}
    }
  return a;
}

smat smat::scalar_matrix(int n, const scalar& a)  // nxn scalar matrix a*I
{
  smat D(n,n); // creates enough space
  for( int i = 0; i < n; i++ )
    {
      D.col[i][0] = 1;   // one entry in this row
      D.col[i][1] = i+1; // ...it's in column i+1
      D.val[i][0] = a;   // ...its value is a
    }
  return D;
}

int liftmat(const smat& mm, scalar pr, smat& m, scalar& dd)
{ int trace=0;
  scalar n,d; long nr,nc; dd=1;
  scalar lim = sqrt(pr>>1);
  m = mm;
  m.reduce_mod_p(pr);
  if(trace)
    {
      cout << "Lifting mod-p smat" << endl;
      if (trace>1) cout << "smat mod "<<pr<<" is:\n" << m.as_mat() <<endl;
      cout << "Now lifting back to Q.\n";
    }
  scalar ma = maxabs(m);
  if (ma < lim)
    {
      if (trace) cout<<"Nothing to do, max entry "<<ma<<" < "<<lim<<endl;
      return 1;
    }
  for(nr=0; nr<m.nro; nr++)
    for(nc=0; nc<m.col[nr][0]; nc++)
      {
        scalar v = m.val[nr][nc];
        if (abs(v) < lim) continue;
	int ok = modrat(v,pr,n,d);
        scalar newdd=lcm(abs(d),dd);
        if (newdd!=dd)
          {
            dd = newdd;
            if(trace)
              cout <<"denom = "<<abs(d)<<", common denom so far is "<<dd<<endl;
          }
        if (!ok)
          {
            if (trace) cerr << "Failed to lift "<<v<<" mod "<<pr<<" to Q"<<endl;
            return 0;
          }
      }
  dd=abs(dd);
  if(trace) cout << "Common denominator = " << dd << "\n";
  for(nr=0; nr<m.nro; nr++)
    for(nc=0; nc<m.col[nr][0]; nc++)
      m.val[nr][nc] = mod(xmodmul(dd,(m.val[nr][nc]),pr),pr);
  if(trace)
    {
      if (trace>1) cout << "Lifted smat = " << m.as_mat() << "\n";
      cout << " Lift has denominator "<<dd<<endl;
    }
  return 1;
}

//#define DEBUG_CHINESE

int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2,
                      smat& m, scalar& dd)
{
  scalar modulus=pr1*pr2;
  scalar n,d,u,v;

  dd = bezout(pr1,pr2,u,v); //==1
  if (dd!=1) return 0;

  // First time through: compute CRTs, common denominator and success flag
  m = m1; // NB We assume that m1 and m2 have nonzero entries in the same places
  for(int nr=0; nr<m1.nro; nr++)
    for(int nc=0; nc<m1.col[nr][0]; nc++)
      {
        scalar mij = xmodmul(v,m1.val[nr][nc],pr1)*pr2 + xmodmul(u,m2.val[nr][nc],pr2)*pr1;
        mij = mod(mij,modulus);
#ifdef DEBUG_CHINESE
        if (((mij-m1.val[nr][nc])%pr1)||((mij-m2.val[nr][nc])%pr2))
          {
            cout<< "bad CRT(["<<m1.val[nr][nc]<<","<<m2.val[nr][nc]<<"],["<<pr1<<","<<pr2<<"]) = "<<mij<<endl;
          }
#endif
        m.val[nr][nc] = mij;
        if (modrat(mij,modulus,n,d))
          dd=lcm(d,dd);
        else
          {
#ifdef DEBUG_CHINESE
            cout<<"CRT("<<m1.val[nr][nc]<<","<<m2.val[nr][nc]<<")="<<mij<<" (mod "<<modulus<<") fails to lift (lim="<<lim<<")\n";
            cout << "Problems encountered in chinese lifting of smat modulo "<<pr1<<" and "<<pr2<< endl;
#endif
            return 0;
          }
      }
  dd=abs(dd);
#ifdef DEBUG_CHINESE
  cout << "Common denominator = " << dd << "\n";
#endif
  // Second time through: rescale
  for(int nr=0; nr<m.nro; nr++)
    for(int nc=0; nc<m.col[nr][0]; nc++)
      {
        m.val[nr][nc] = mod(xmodmul((dd/d),m.val[nr][nc],modulus),modulus);
      }
  return 1;
}

#if FLINT

#include "eclib/flinterface.h"

// FLINT has more than one type for modular matrices: standard in
// FLINT-2.3..2.9 was nmod_mat_t with entries of type mp_limb_t
// (unsigned long) while non-standard was hmod_mat_t, with entries
// hlimb_t (unsigned int).  From FLINT-3 the latter is emulated via a
// wrapper.  We use the former when scalar=long and the latter when
// scalar=int and the FLINT version is at least 3.  The unsigned
// scalar types are #define'd as uscalar.

void mod_mat_from_smat(mod_mat& A, const smat& M, scalar pr)
{
  long nr=M.nrows(), nc=M.ncols();
  long i, j;

  // copy of the modulus for FLINT
  uscalar p = (uscalar)I2long(pr);

  // create flint matrix copy of M:
  mod_mat_init(A, nr, nc, p);
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++)
      mod_mat_entry(A,i,j) = (uscalar)I2long(posmod(M.elem(i+1,j+1),pr));
}

smat smat_from_mod_mat(const mod_mat& A, const scalar& p) //scalar just to fix return type
{
  long nr=mod_mat_nrows(A), nc=mod_mat_ncols(A);

  // create matrix copy of A:
  smat M(nr, nc);
  long i, j;
  for(i=0; i<nr; i++)
    {
      svec rowi(nc);
      for(j=0; j<nc; j++)
        rowi.set(j+1, scalar(mod_mat_entry(A,i,j)));
      M.setrow(i+1,rowi);
    }
  return M;
}

smat mult_mod_p_flint ( const smat& A, const smat& B, const scalar& pr )
{
  if( A.ncols() != B.nrows() )
    {
      cerr << "incompatible smats in operator *"<<endl;
      return smat();
    }
  mod_mat A1, B1, C1;
  mod_mat_from_smat(A1,A,pr);
  mod_mat_from_smat(B1,B,pr);
  mod_mat_init(C1, A.nrows(), B.ncols(), I2long(pr));
  // timer T;
  // T.start();
  mod_mat_mul(C1,A1,B1);
  // T.stop();
  // cout<<"mult_mod_p_flint time (size "<<dim(A)<<"x"<<dim(B)<<"): ";
  // T.show();
  smat C = smat_from_mod_mat(C1, pr);
  mod_mat_clear(A1);
  mod_mat_clear(B1);
  mod_mat_clear(C1);
  return C;
}

#endif
