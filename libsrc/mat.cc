// mat.cc: implementation of integer matrix classes
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
 
// Only to be included by matrix.cc

#define LONGLONG scalar   // long longs caused memory alocation bugs
                        // N.B. if = long long, uncomment abs() below

// Definitions of member operators and functions:

mat::mat(long nr, long nc)
{
 nro=nr;
 nco=nc;
 long n=nr*nc;
 entries=new scalar[n];  
 if (!entries)
   cerr<<"Out of memory in mat constructor!"<<endl;
 else
   memset(entries, 0, n*sizeof(scalar));
}

mat::mat(const mat& m)
{
 nro=m.nro;
 nco=m.nco;
 long n=nro*nco;
 entries=new scalar[n];
 if (!entries)
   cerr<<"Out of memory in mat constructor!"<<endl;
 else
   memcpy(entries, m.entries, n*sizeof(scalar));
}

mat::~mat()
{
  delete[] entries;
}

void mat::init(long nr, long nc) // assigns to zero mat of given size;
{                                 // with defaults (0,0) releases all space.
 long n = nr*nc;
 if (nro*nco!=n) // delete old space; 
   {            // replace with new.
     delete[] entries;
     entries = new scalar[n]; 
   }
 if (!entries)
   cerr<<"Out of memory in mat::init"<<endl;
 else
   {
     nro = nr;
     nco = nc;
     memset(entries, 0, n*sizeof(scalar));
   }
}

scalar& mat::operator()(long i, long j)  const   // returns ref to (i,j) entry
{
   return entries[(i-1)*nco+(j-1)];
}


mat mat::slice(long r1,long r2,long c1,long c2) const
{
  if(c1<0) // abbreviated form with firsts=1
    {
      c2=r2-1; r2=r1-1; r1=c1=0;
    }
  else 
    {
      r1--; c1--; r2--; c2--;
    }
 long n=r2-r1+1,c=c2-c1+1; long cc=c;
 mat ans(n,cc);
 scalar* ap=ans.entries, *mp=entries+r1*nco+c1;
 while(n--)
   {
     memcpy(ap, mp, cc*sizeof(scalar));
     ap += cc;
     mp += nco;
   }
 return ans;
}

mat& mat::operator=(const mat& m)
{
 if (this==&m) return *this;
 long n = m.nro*m.nco;
 if (nro*nco!=n) // delete old space; 
   {            // replace with new.
     delete[] entries;
     entries = new scalar[n]; 
   }
 if (!entries)
   cerr<<"Out of memory in mat assignment!"<<endl;
 else
   {
     nro = m.nro;
     nco = m.nco;
     scalar *m1=entries, *m2=m.entries;
     while(n--) *m1++ = *m2++;
   }
 return *this;
}

scalar mat::sub(long i, long j) const
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) return entries[(i-1)*nco+(j-1)];
  else 
    {
      cerr << "Bad indices ("<<i<<","<<j<<") in mat::sub (nro="<<nro
	   <<", nco="<<nco<<endl;
      return 0;
    }
}

void mat::set(long i, long j, scalar x)
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) entries[(i-1)*nco+(j-1)] = x;
 else 
   {
     cerr << "Bad indices ("<<i<<","<<j<<") in mat::set (nro="<<nro
	  <<", nco="<<nco<<endl;
   }
}

void mat::add(long i, long j, scalar x)
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) entries[(i-1)*nco+(j-1)] += x;
 else 
   {
     cerr << "Bad indices ("<<i<<","<<j<<") in mat::add (nro="<<nro
	  <<", nco="<<nco<<endl;
   }
}

void mat::setrow(long i, const vec& v)
{
 if ((0<i) && (i<=nro) && (dim(v)==nco))
  {
    scalar * rowi = entries + (i-1)*nco;
    scalar * vec = v.entries;
    long c=nco; 
    while(c--) *rowi++ = *vec++;
  }
 else 
   {
     cerr << "Bad indices in mat::setrow (i="<<i<<", nro="<<nro
	  <<", dim(v)="<<dim(v)<<", nco="<<nco<<")"<<endl;
   }
}

void mat::setcol(long j, const vec& v)
{
 if ((0<j) && (j<=nco) && (dim(v)==nro))
  {
   scalar * colj = entries+(j-1);
   scalar * vec = v.entries;
   long n=nro;
   while(n--) {*colj = *vec++; colj+=nco;}
 }
 else
   {
     cerr << "Bad indices in mat::setcol (j="<<j<<", nco="<<nco
	  <<", dim(v)="<<dim(v)<<", nco="<<nco<<")"<<endl;
   }
}

vec mat::row(long i) const
{
 vec mi(nco);
 if ((0<i) && (i<=nro))
   memcpy(mi.entries, entries+(i-1)*nco, nco*sizeof(scalar));
 else
   cerr << "Bad row number "<<i<<" in function mat::row (nro="<<nro<<")"<<endl;
 return mi;
}

vec mat::col(long j) const
{
 vec mj(nro);
 long i=nro; scalar *entriesij=entries+(j-1), *v=mj.entries;
 if ((0<j) && (j<=nco))
   while(i--) {*v++ = *entriesij; entriesij+=nco;}
 else
   cerr << "Bad column number "<<j<<" in function mat::col (nco="<<nco<<")"<<endl;
 return mj;
}

void mat::swaprows(long r1, long r2)
{
  if ((0<r1)&&(0<r2)&&(r1<=nro)&&(r2<=nro))
    {
      scalar *mr1 = entries + (r1-1)*nco;
      scalar *mr2 = entries + (r2-1)*nco;
      long nc=nco;
      while(nc--) {long a = *mr1; *mr1++ = *mr2; *mr2++ = a; }
    }
  else
    {
      cerr << "Bad row numbers " << r1 << "," << r2 << " in swaprow (nro="<<nro<<")"<<endl;
    }
}

void mat::multrow(long r, scalar scal)
{
  if ((0<r)&&(r<=nro))
    {
      long nc=nco; scalar *mij = entries+(r-1)*nco;
      while(nc--) (*mij++) *= scal;
    }
  else
    {
      cerr << "Bad row number " << r << " in multrow (nro="<<nro<<")"<<endl;
    }
}

void mat::divrow(long r, scalar scal)
{
  if ((0<r)&&(r<=nro))
    {
      long nc=nco; scalar *mij = entries+(r-1)*nco;
      while(nc--) (*mij++) /= scal;
    }
  else
    {
      cerr << "Bad row number " << r << " in divrow (nro="<<nro<<")"<<endl;
    }
}

void mat::clearrow(long r)
{
  if ((0<r)&&(r<=nro))
    {
      scalar g=0; long nc=nco; scalar * mij = entries+(r-1)*nco;
      while((nc--)&&(g!=1)) g = gcd(g,(*mij++));
      if (g>1)
	{
	  nc=nco; mij = entries+(r-1)*nco;
	  while(nc--) (*mij++) /= g;
	}
    }
  else 
    {
      cerr << "Bad row number " << r << " in clearrow (nro="<<nro<<")"<<endl;
    }
}

mat& mat::operator+=(const mat& entries2)
{
  if ((nro==entries2.nro) && (nco=entries2.nco)) 
    {
      long n=nro*nco; scalar *m1=entries, *m2=entries2.entries;
      while(n--) (*m1++) += (*m2++);
    }
  else 
    {
      cerr << "Incompatible matrices in operator +="<<endl;
    }
  return *this;
}

mat& mat::operator-=(const mat& entries2)
{
  if ((nro==entries2.nro) && (nco=entries2.nco)) 
    {
      long n=nro*nco; scalar *m1=entries, *m2=entries2.entries;
      while(n--) (*m1++) -= (*m2++);
    }
  else 
    {
      cerr << "Incompatible matrices in operator -="<<endl;
    }
  return *this;
}

mat& mat::operator*=(scalar scal)
{
  scalar* mij = entries; long n=nco*nro;
  while(n--) (*mij++) *= scal;
  return *this;
}

mat& mat::operator/=(scalar scal)
{
  scalar* mij = entries; long n=nco*nro;
  while(n--) (*mij++) /= scal;
  return *this;
}


// Definitions of non-member, friend operators and functions

// add/sub row i of mat to v (implemented in mat.cc)
void add_row_to_vec(const vec& v, const mat& m, long i)
{
  scalar* vi=v.entries, *mij=m.entries+(i-1)*m.nco; long j=v.d;
  while(j--)(*vi++)+=(*mij++);
}

void sub_row_to_vec(const vec& v, const mat& m, long i)
{
  scalar* vi=v.entries, *mij=m.entries+(i-1)*m.nco; long j=v.d;
  while(j--)(*vi++)-=(*mij++);
}

mat operator*(const mat& m1, const mat& m2)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 mat m3(m,p);
 scalar *a=m1.entries, *b=m2.entries, *c=m3.entries;
 if (n==m2.nro)  // algorithm from Dr Dobb's Journal August 1993
   {
     while(m--)
       {
	 scalar *bp=b;
         long k=n;
	 while(k--)
	   {
	     scalar *cp=c;
             long j=p;
	     while(j--)
	       {
		 *cp++ += *a * *bp++;
	       }
	     a++;
	   }
	 c += p;
       }
   }
 else
   {
     cerr << "Incompatible sizes in mat product"<<endl;
   }
 return m3;
}

int operator==(const mat& m1, const mat& m2)
{
   long nr=m1.nro, nc=m1.nco;
   int equal = ((nr==m2.nro) && (nc==m2.nco)); 
   if(!equal) return 0;
   scalar *m1ij=m1.entries, *m2ij=m2.entries; long n=nr*nc;
   while((n--)&&equal) equal=((*m1ij++)==(*m2ij++));
   return equal;
}

void mat::output(ostream& s) const
{
  scalar* mij=entries;
  s << "\n[";
  long nr=nro;
  while(nr--)
    {
      long nc=nco;
      s<<"[";
      while(nc--) {s<<(*mij++); if(nc) s<<",";}
      s<<"]"; if(nr) s<<",\n";
    }
  s << "]\n";
}

void mat::output_pari(ostream& s) const
{
  scalar* mij=entries;
  s << "\n[";
  long nr=nro;
  while(nr--)
    {
      long nc=nco;
      while(nc--) {s<<(*mij++); if(nc) s<<",";}
      if(nr) s<<";";
    }
  s << "]\n";
}

long ndigits(scalar a) 
{
  static double log10 = log((double)10);
  if(a==0) return 1;
  long s=0;
  if(a<0) {a=-a; s=1;}
  return s+(long)floor(log((double)a)/log10)+1;
}

void mat::output_pretty(ostream& s) const
{
  long i,j; scalar* mij;
  long nc=nco,nr=nro;
  int* colwidths = new int[nc];
  for(j=0; j<nco; j++)
    {
      scalar ma=0, mi=0;
      for(i=0, mij=entries+j; i<nro; i++, mij+=nc)
	{
	  if (*mij>ma) ma=*mij;
	  else if (*mij<mi) mi=*mij;
	}
      ma=ndigits(ma);
      mi=ndigits(mi);
      if(mi>ma)ma=mi;
      colwidths[j]=ma;
    }
  mij=entries;
  while(nr--)
    {
      s << "[";
      nc=nco;  j=0;
      while(nc--) 
	{
	  s.width(colwidths[j]); s<<(*mij);
	  //	  s.form("%*d",colwidths[j],*mij); 
	  if(nc) s<<" ";
	  mij++; j++; 
	}
      s<<"]\n";
    }
  delete[]colwidths;
}

void mat::dump_to_file(string filename) const
{
  ofstream fout(filename.c_str(),ofstream::binary);
  fout.write((char*)&nro,sizeof(nro));
  fout.write((char*)&nco,sizeof(nco));
  fout.write((char*)entries,nro*nco*sizeof(scalar));
  fout.close();
}

void mat::read_from_file(string filename)
{
  ifstream fin(filename.c_str());
  fin.read((char*)&nro,sizeof(nro));
  fin.read((char*)&nco,sizeof(nco));
  delete[] entries;
  entries = new scalar[nro*nco];
  fin.read((char*)entries,nro*nco*sizeof(scalar));
  fin.close();
}

istream& operator>>(istream& s, mat& m)
{
 long n=m.nro*m.nco;
 scalar* mij=m.entries;
 while(n--) s >> (*mij++);
 return s;
}

mat colcat(const mat& a, const mat& b)
{
 long nr = a.nro, nca = a.nco, ncb = b.nco;
 mat ans(nr,nca+ncb);
 scalar *ansij=ans.entries, *aij=a.entries, *bij=b.entries;
 if (nr==b.nro)
   while(nr--)
     {
       long nc=nca;
       while(nc--) *ansij++ = *aij++;
       nc=ncb;
       while(nc--) *ansij++ = *bij++;
     }
 else
   {
     cerr << "colcat: matrices have different number of rows!" << endl;
   }
 return ans;
}

mat rowcat(const mat& a, const mat& b)
{
 long nra = a.nro, nc = a.nco, nrb = b.nro;
 mat ans(nra+nrb,nc);
 scalar *ansij=ans.entries, *aij=a.entries, *bij=b.entries;
 if (nc==b.nco)
 {
   long n = nra*nc;
   while(n--) *ansij++ = *aij++;
   n = nrb*nc;
   while(n--) *ansij++ = *bij++;
 }
 else
   {
     cerr << "rowcat: matrices have different number of columns!" << endl;
   }
 return ans;
}

mat directsum(const mat& a, const mat& b)
{
  long n,c, nra=a.nro, nca=a.nco, nrb=b.nro, ncb=b.nco;
  mat ans(nra+nrb,nca+ncb);
  scalar* ansij=ans.entries, *aij=a.entries, *bij=b.entries;
  n=nra; 
  while(n--) 
    {
      c=nca; while(c--) *ansij++ = *aij++;
      c=ncb; while(c--) *ansij++ = 0;
    }
  n=nrb; 
  while(n--) 
    {
      c=nca; while(c--) *ansij++ = 0;
      c=ncb; while(c--) *ansij++ = *bij++;
    }
  return ans;
}

void elimrows(mat& m, long r1, long r2, long pos)   //plain elimination, no clearing
{
  long nc=m.nco;
  scalar *mr1 = m.entries + (r1-1)*nc,
    *mr2 = m.entries + (r2-1)*nc;
  scalar p = mr1[pos-1], q=mr2[pos-1];
  while(nc--)
    {
      (*mr2)= (p*(*mr2))-(q*(*mr1));
      mr1++; mr2++;
    }
}

void elimrows1(mat& m, long r1, long r2, long pos)      //elimination + clearing
{
  elimrows(m,r1,r2,pos);
  m.clearrow(r2);
}

void elimrows2(mat& m, long r1, long r2, long pos, scalar last) //elimination + divide by last pivot
{
  elimrows(m,r1,r2,pos);
  m.divrow(r2,last);
}

// Definition of non-friend functions

mat operator+(const mat& m)
{
        return m;
}

mat operator-(const mat& m)
{
        return (-1)*m;
}

mat operator+(const mat& m1, const mat& m2)
{
  mat ans(m1); ans+=m2;  return ans;
}

mat operator-(const mat& m1, const mat& m2) 
{
  mat ans(m1); ans-=m2;  return ans;
}

mat operator*(scalar scal, const mat& m)
{
  mat ans(m); ans*=scal;  return ans;
}

mat operator/(const mat& m, scalar scal)
{
  mat ans(m); ans/=scal;  return ans;
}

int operator!=(const mat& m1, const mat& m2)
{
  return !(m1==m2);
}

vec operator*(const mat& m, const vec& v)
{
 long r=m.nro, c=m.nco;
 vec w(r);
 if (c==v.d)
   {
     scalar *mp=m.entries, *wp=w.entries;
     while(r--)
       {
	 scalar *vp=v.entries;
         c=m.nco;
	 while(c--)
           *wp += (*mp++)*(*vp++);
	 wp++;
       }
   }
 else
   {
     cerr << "Incompatible sizes in *(mat,vec)"<<endl;
   }
 return w;
}

mat idmat(scalar n)
{
 mat ans(n,n);
 long i;
 for (i=1; i<=n; i++) ans.set(i,i,1);
 return ans;
}

mat transpose(const mat& m)
{
 long i,j,nr,nc;
 nr=m.ncols(); nc=m.nrows();
 mat ans(nr, nc);
 for (i=1; i<=nr; i++)
  for (j=1; j<=nc; j++)
   ans.set(i,j,  m(j,i));
 return ans;
}

mat submat(const mat& m, const vec& iv, const vec& jv)
{long i,j;
 long nr = dim(iv);
 long nc = dim(jv);
 mat ans(nr,nc);
 for (i=1; i<=nr; i++)
  for (j=1; j<=nc; j++)
   ans.set(i,j, m(iv[i],jv[j]));
 return ans;
}

mat echelon(const mat& entries, vec& pcols, vec& npcols,
               long& rk, long& ny, scalar& d, int method)
{
  switch (method)
    {case 0: return echelon0(entries,pcols,npcols,rk,ny,d);
     case 1: return echelonl(entries,pcols,npcols,rk,ny,d);
     case 2: return echelonp(entries,pcols,npcols,rk,ny,d,DEFAULT_MODULUS);
     default: return echelon0(entries,pcols,npcols,rk,ny,d);
    }        
}

//#define DEBUG_ECH_0

//N.B. if(q==0) the following multiplies row r2 by p, which looks
//redundant.  However, it is important to keep this in as in echelon0
//we must guarentee divisibility by "lastpivot".  We do not want to keep 
//computing contents of rows as this is slower.
// Used in forward elimination in echelon0

void conservative_elim(scalar *m, long nc, long r1, long r2, long pos)
{scalar *mr1=m+r1*nc, *mr2=m+r2*nc;
 scalar p = mr1[pos], q = mr2[pos];
 if(p==1)
   if(q==0) {;} // nothing to do
   else
     if(q==1)
       while(nc--) 
	 {
	   (*mr2)-=(*mr1); 
	   mr1++; mr2++; 
	 } 
     else // general q
       while(nc--) 
	 {
	   (*mr2)-=(q*(*mr1)); 
	   mr1++; mr2++; 
	 } 
 else  // p!=1; we cannot assume p>0
   if(q==0) // must still multiply r2 by p
     while(nc--) 
       {
	 (*mr2)*=p; 
	 mr2++; 
       } 
   else
     if(q==1)
       while(nc--) 
	 {
	   (*mr2)=(p*(*mr2))-(*mr1); 
	   mr1++; mr2++; 
	 } 
     else // general q
       while(nc--) 
	 {
	   (*mr2)=(p*(*mr2))-(q*(*mr1)); 
	   mr1++; mr2++; 
	 } 
}

// This version does not multiply row r1 by p unnecessarily 
// (used in back substitution)

void elim(scalar *m, long nc, long r1, long r2, long pos)
{scalar *mr1=m+r1*nc, *mr2=m+r2*nc;
 scalar p = mr1[pos], q = mr2[pos];
#ifdef DEBUG_ECH_0
long n;
cout<<"In elim with p = "<<p<<" and q = " << q << endl;
cout<<"row 1: "; for(n=0; n<nc; n++) cout<<mr1[n]<<",";  cout<<endl;
cout<<"row 2: "; for(n=0; n<nc; n++) cout<<mr2[n]<<",";  cout<<endl;
#endif
 if(p==1)
   if(q==0) {;} // nothing to do
   else
     if(q==1)
       while(nc--) 
	 {
	   (*mr2)-=(*mr1); 
	   mr1++; mr2++; 
	 } 
     else // general q
       while(nc--) 
	 {
	   (*mr2)-=(q*(*mr1)); 
	   mr1++; mr2++; 
	 } 
 else  // p!=1; we cannot assume p>0
   if(q==0)    {;} // nothing to do
   else
     if(q==1)
       while(nc--) 
	 {
	   (*mr2)=(p*(*mr2))-(*mr1); 
	   mr1++; mr2++; 
	 } 
     else // general q
       while(nc--) 
	 {
	   (*mr2)=(p*(*mr2))-(q*(*mr1)); 
	   mr1++; mr2++; 
	 } 
}

void clear(scalar* row, long nc)
{long n=nc; scalar *rowi=row; scalar g = 0;
 while((n--)&&(g!=1)) g=gcd(g,*rowi++);
 if (g<0)g=-g;
 if (g>1) 
   {
     n=nc; rowi=row; while(n--) (*rowi++) /= g;
   }
}

//#ifndef DEBUG_ECH_0
//#define DEBUG_ECH_0
//#endif

#ifdef DEBUG_ECH_0
void show(scalar* m, long nr, long nc)
{
  long i,j; scalar* mij = m;
  for(i=0; i<nr; i++)
    {
      for(j=0; j<nc; j++)
	{cout<<(*mij)<<"\t"; mij++;}
      cout<<"\n";
    }
}
#endif

mat echelon0(const mat& entries, vec& pc, vec& npc,
                long& rk, long& ny, scalar& d)
{
  long nr, nc, r,c,r2,r3,rmin,i;
  scalar min, mr2c, lastpivot=1;
  rk=0; ny=0; r=0;
  nc=entries.nco; nr=entries.nro;
  scalar *m, *mi1, *mi2, *mij; scalar temp;
  m = new scalar[nr*nc];
  long n=nr*nc; mij=m; mi1=entries.entries;
  while(n--) *mij++ = *mi1++;

  scalar *pcols = new scalar[nc];
  scalar *npcols = new scalar[nc];
  for (c=0; (c<nc)&&(r<nr); c++)
    {
      mij=m+r*nc+c;  // points to column c in row r
      min = abs(*mij);  rmin = r;
      for (r2=r+1, mij+=nc; (r2<nr)&&(min!=1); r2++, mij+=nc)
       { mr2c = abs(*mij);
         if ((0<mr2c) && ((mr2c<min) || (min==0))) { min=mr2c; rmin=r2 ;}
       }
     if (min==0) npcols[ny++] = c;
     else
       {pcols[rk++] = c;
#ifdef DEBUG_ECH_0
       cout<<"Using col "<<c<<" as pivotal col"<<endl;
#endif
        if (rmin>r) //swap rows
	  {
#ifdef DEBUG_ECH_0
	    cout<<"Swapping rows "<<r<<" and "<<rmin<<endl;
#endif
	    mi1=m+r*nc; mi2=m+rmin*nc; n=nc;
	    while(n--) {temp = *mi1; *mi1++ = *mi2; *mi2++ = temp;}
	  } 
        for (r3 = r+1 ; r3<nr; r3++)
          {
#ifdef DEBUG_ECH_0
	    cout<<"Eliminating from row "<<r3<<endl;
#endif
            conservative_elim(m,nc,r,r3,c);
	    if(lastpivot>1)
	      {
		mi1 = m+r3*nc; n=nc; 
		while(n--) 
		  {
		    if(*mi1%lastpivot)
		      cout<<"Error in echelon0!  Entry "<<(*mi1)
			  <<" not divisible by lastpivot "<<lastpivot<<endl;
		    *mi1++ /= lastpivot;
		  }
	      }
          }
        lastpivot=min;
#ifdef DEBUG_ECH_0
cout<<"r="<<r<<": pivot = "<<min<<endl;
#endif
        r++;
       }
#ifdef DEBUG_ECH_0
     //     cout<<"Current mat is:\n";show(m,nr,nc);
#endif
    }
  for (c = rk+ny; c<nc; c++) npcols[ny++] = c;
#ifdef DEBUG_ECH_0
cout<<"After forward elimination, rank = "<<rk<<"; pivots are:"<<endl;
for(r3=0; r3<rk; r3++) cout<<(m+r3*nc)[pcols[r3]]<<",";
cout<<endl;
#endif
  d=1;
  if (ny>0)   // Back-substitute and even up pivots
    {for (r=0; r<rk; r++) clear(m+r*nc,nc);
#ifdef DEBUG_ECH_0
cout<<"After clearing, pivots are:"<<endl;
for(r3=0; r3<rk; r3++) cout<<(m+r3*nc)[pcols[r3]]<<",";
cout<<endl;
#endif
     for (r=0; r<rk; r++)
       {
	 mi1=m+r*nc;
#ifdef DEBUG_ECH_0
cout<<"Before back-subst, row "<<r<<" is:"<<endl;
for(r3=0; r3<nc; r3++) cout<<mi1[r3]<<",";
cout<<": pivot = "<<mi1[pcols[r]]<<endl;
#endif
	 for (r2=r+1; r2<rk; r2++) elim(m,nc,r2,r,pcols[r2]);  
#ifdef DEBUG_ECH_0
cout<<"After back-subst, row "<<r<<" is:"<<endl;
for(r3=0; r3<nc; r3++) cout<<mi1[r3]<<",";
cout<<": pivot = "<<mi1[pcols[r]]<<endl;
#endif
	 clear(mi1,nc);
#ifdef DEBUG_ECH_0
cout<<"After clearing, row "<<r<<" is:"<<endl;
for(r3=0; r3<nc; r3++) cout<<mi1[r3]<<",";
x cout<<": pivot = "<<mi1[pcols[r]]<<endl;
#endif
	 d = lcm(d,mi1[pcols[r]]);
       }
     d = abs(d);
     // cout << "d = " << d << "\n";
     for (r=0, mij=m; r<rk; r++)
       {
	 n=nc;
	 scalar fac = d/mij[pcols[r]];  
	 while(n--) *mij++ *= fac;
       }
   }
  else 
    {
      mij=m;
      for (r=0; r<rk; r++)
	for (c=0; c<nc; c++)
	  *mij++ = (c==pcols[r]);  // 0 or 1 !
    }
  // Copy back into mat
  mat ans(rk,nc);
  n=rk*nc; scalar* ansij=ans.entries; mij=m;
  while(n--) *ansij++ = *mij++; 
  
  delete[] m;
  // fix vectors
  pc.init(rk); npc.init(ny);
  for (i=0; i<rk; i++)  pc[i+1]= pcols[i]+1;
  for (i=0; i<ny; i++) npc[i+1]=npcols[i]+1;
  delete[] pcols;
  delete[] npcols;

  return ans;
}

long mat::rank() const
{
 long rk=0, nr,nc,r=1,c,r2,r3;
 long mr2c,lastpivot=1;
 mat m(*this);
 nc=m.ncols(); nr=m.nrows();
 for (c=1; (c<=nc)&&(r<=nr); c++)
 { long min = abs(m(r,c));
   long rmin = r;
   for (r2=r+1; (r2<=nr)&&(min!=1); r2++)
   { mr2c = abs(m(r2,c));
     if ((0<mr2c) && ((mr2c<min) || (min==0))) { min=mr2c; rmin=r2 ;}
   }
   if (min!=0)
     {rk++;
      if (rmin>r) m.swaprows(r,rmin);
      for (r3 = r+1 ; r3<=nr; r3++)
         elimrows2(m,r,r3,c,lastpivot);
      lastpivot=min;
      r++;
     }
 }
 return rk;
}

long mat::nullity() const
{
 return nco-rank();
}

long mat::trace() const
{ long i=0; scalar* aii=entries; long ans=0;
  for (; i<nro; i++, aii+=(nco+1))
    ans += *aii;
  return ans;
}

// FADEEV'S METHOD

vector<long> mat::charpoly() const
{ long n = nrows();
  mat b(*this);
  mat id(idmat((scalar)n));
  vector<long> clist(n+1);
  long t = trace();
  clist[n]   =  1;
  clist[n-1] = -t;
  for (long i=2; i<=n; i++)
    { b=(*this)*(b-t*id);          //     cout << b;   // (for testing only)
        t=b.trace()/i;
        clist[n-i] = -t;
      }
  if (!(b==t*id)) 
    {
      cerr << "Error in charpoly: final b = " << (b-t*id) << endl;
    }
  return clist;
}

long mat::determinant() const
{
 long det = charpoly()[0];
 if (nrows()%2==1)
   return -det;
 else
   return det;
}

void vec::sub_row(const mat& m, int i)
{
  scalar* vi=entries, *wi=m.entries+(i-1)*d; long n=d;
  if (d==m.ncols()) {while(n--)(*vi++)-=(*wi++);}
  else {cerr << "Incompatible vecs in vec::sub_row"<<endl;}
}

void vec::add_row(const mat& m, int i)
{
  scalar* vi=entries, *wi=m.entries+(i-1)*d; long n=d;
  if (d==m.ncols()) {while(n--)(*vi++)+=(*wi++);}
  else {cerr << "Incompatible vecs in vec::add_row(): d="<<d<<" but m has "<<m.ncols()<<"cols"<<endl;}
}

mat addscalar(const mat& m, scalar c)
{
  return m+(c*idmat((scalar)m.nrows()));
}
 
vec apply(const mat& m, const vec& v)    // same as *(mat, vec)
{
 long nr=m.nrows(), nc=m.ncols();
 vec ans(nr);
 if (nc==dim(v))
   for (long i=1; i<=nr; i++) 
     ans[i] = m.row(i)*v;
 else 
   {
     cerr << "Incompatible sizes in *(mat,vec)"<<endl;
   }
 return ans;
}

/*  Need this when LONGLONG = long long!
LONGLONG abs(LONGLONG a)
{return ((a<0) ? -a : a);
}
*/

LONGLONG lgcd(LONGLONG aa, LONGLONG bb)
{LONGLONG a=aa,b=bb;
 while (b!=0) {LONGLONG c=a%b; a=b; b=c;}
 return ((a<0) ? -a : a);
}

LONGLONG llcm(LONGLONG a, LONGLONG b)
{LONGLONG g = lgcd(a,b);
 if (g==0) return 0;
 return a*(b/g);
}

void lelim(LONGLONG *m, long nc, long r1, long r2, long pos)
{LONGLONG *mr1=m+r1*nc, *mr2=m+r2*nc;
 scalar p = mr1[pos], q = mr2[pos];
 while(nc--)
   {
     (*mr2)=(p*(*mr2))-(q*(*mr1)); 
     mr1++; mr2++;
   }
}

void lclear(LONGLONG* row, long nc)
{long n=nc; LONGLONG *rowi=row; LONGLONG g = 0;
 while((n--)&&(g!=1)) g=lgcd(g,*rowi++);
 if (g<0)g=-g;
 if (g>1) 
   {
     n=nc; rowi=row; while(n--) (*rowi++) /= g;
   }
}

// The following version of echelon uses long-long-integers for internal
// calculation and a minimum of function calls; the structures vec and
// mat are used for input/output but not internally.

mat echelonl(const mat& entries, vec& pc, vec& npc,
                long& rk, long& ny, scalar& d)
{
  rk=0; ny=0;
  long nr=entries.nro, nc=entries.nco, r=0, c,r2,r3,i, mr2c,lastpivot=1;
  LONGLONG *m, *mi1, *mi2, *mij; LONGLONG temp;
  m = new LONGLONG[nr*nc];
  long n=nr*nc; mij=m; mi1=entries.entries;
  while(n--) *mij++ = (LONGLONG)(*mi1++);

  scalar *pcols = new scalar[nc];
  scalar *npcols = new scalar[nc];
  for (c=0; (c<nc)&&(r<nr); c++)
    {
      mij=m+r*nc+c;  // points to column c in row r
      long min = abs(*mij);
      long rmin = r;
      for (r2=r+1, mij+=nc; (r2<nr)&&(min!=1); r2++, mij+=nc)
       { mr2c = abs(*mij);
         if ((0<mr2c) && ((mr2c<min) || (min==0))) { min=mr2c; rmin=r2 ;}
       }
     if (min==0) npcols[ny++] = c;
     else
       {pcols[rk++] = c;
        if (rmin>r) //swap rows
	  {
	    mi1=m+r*nc; mi2=m+rmin*nc; n=nc;
	    while(n--) {temp = *mi1; *mi1++ = *mi2; *mi2++ = temp;}
	  } 
        for (r3 = r+1 ; r3<nr; r3++)
          {
            lelim(m,nc,r,r3,c);
	    mi1 = m+r3*nc; n=nc; while(n--) *mi1++ /= lastpivot;
          }
        lastpivot=min;
        r++;
      }
    }
  for (c = rk+ny; c<nc; c++) npcols[ny++] = c;
  d=1;
  if (ny>0)   // Back-substitute and even up pivots
    {for (r=0; r<rk; r++) lclear(m+r*nc,nc);
     for (r=0; r<rk; r++)
       {
	 for (r2=r+1; r2<rk; r2++) lelim(m,nc,r2,r,pcols[r2]);  
	 mi1=m+r*nc;
	 lclear(mi1,nc);
	 d = llcm(d,mi1[pcols[r]]);
       }
     d = abs(d);
     // cout << "d = " << d << "\n";
     for (r=0, mij=m; r<rk; r++)
       {
	 n=nc;
	 scalar fac = d/mij[pcols[r]];  
	 while(n--) *mij++ *= fac;
       }
   }
  else 
    {
      mij=m;
      for (r=0; r<rk; r++)
	for (c=0; c<nc; c++)
	  *mij++ = (c==pcols[r]);  // 0 or 1 !
    }
  // Copy back into mat
  mat ans(rk,nc);
  n=rk*nc; scalar* ansij=ans.entries; mij=m;
  while(n--)
      *ansij++=*mij++;

  delete[] m;
  // fix vectors
  pc.init(rk); npc.init(ny);
  for (i=0; i<rk; i++)  pc[i+1]= pcols[i]+1;
  for (i=0; i<ny; i++) npc[i+1]=npcols[i]+1;
  delete[] pcols;
  delete[] npcols;

  return ans;
}


void elimp(const mat& m, long r1, long r2, long pos, scalar pr)
{
  long nc=m.nco;
  scalar *mr1 = m.entries + (r1-1)*nc + (pos-1),
         *mr2 = m.entries + (r2-1)*nc + (pos-1);
  scalar p = xmod(*mr1,pr), q=xmod(*mr2,pr);
  nc -= (pos-1);
 if(p==1)
   { 
     if(q==0) {return;} // nothing to do
     if(q==1)
       {
         while(nc--)
           {
             (*mr2)= xmod(*mr2-*mr1,pr);
             mr1++; mr2++;
           }
         return;
       }
     if(q==-1)
       {
         while(nc--)
           {
             (*mr2)= xmod(*mr2+*mr1,pr);
             mr1++; mr2++;
           }
         return;
       }
     // general q
     while(nc--)
       {
         (*mr2)= xmod(*mr2-xmodmul(q,*mr1,pr),pr);
         mr1++; mr2++;
       }
     return;
   }
 // general p (p!=1)
 if(q==0) {return;} // nothing to do
 if(q==1)
   { 
     while(nc--)
       {
	 (*mr2)= xmod(xmodmul(p,*mr2,pr)-*mr1,pr);
	 mr1++; mr2++;
       }
     return;
   }
 if(q==-1)
   { 
     while(nc--)
       {
	 (*mr2)= xmod(xmodmul(p,*mr2,pr)+*mr1,pr);
	 mr1++; mr2++;
       }
     return;
   }
 // general q
 while(nc--)
   {
     (*mr2)= xmod(xmodmul(p,*mr2,pr)-xmodmul(q,*mr1,pr),pr);
     mr1++; mr2++;
   }
}

void elimp1(const mat& m, long r1, long r2, long pos, scalar pr)
//same as elimp except assumes pivot is 1
{
 long nc=m.nco;
 scalar *mr1 = m.entries + (r1-1)*nc,
        *mr2 = m.entries + (r2-1)*nc;
 scalar q=xmod(mr2[pos-1],pr);
 if(q==0) return;
 if(q==1)
   { 
     while(nc--)
       {
         (*mr2)= xmod(*mr2-*mr1,pr);
         mr1++; mr2++;
       }
     return;
   }
 if(q==-1)
   { 
     while(nc--)
       {
         (*mr2)= xmod(*mr2+*mr1,pr);
         mr1++; mr2++;
       }
     return;
   }
  // general q
 while(nc--)
   {
     if(*mr1)
       (*mr2)= xmod((*mr2)-xmodmul(q,*mr1,pr),pr);
     mr1++; mr2++;
   }
}

//#define TRACE 1 

// This method uses mod-p arithmetic internally but returns the
// "characteristic zero" echelon form of the mat.  It will only give
// the wrong answer if (a) the rank mod pr is not the actual rank, or (b)
// the actual echelon form has entries which are too big.

mat echelonp(const mat& entries, vec& pcols, vec& npcols,
             long& rk, long& ny, scalar& d, scalar pr)
{
#ifdef TRACE
                cout << "In echelonp\n";
#endif /* TRACE */
 long nc,nr,r,c,r2,r3,rmin;
 scalar min, mr2c;
 nr=entries.nrows(), nc=entries.ncols();
 mat m(nr,nc);
 scalar *mij=m.entries, *entriesij=entries.entries;
 long n=nr*nc;
 while(n--) *mij++ = xmod(*entriesij++,pr);
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0; r=1;
 for (c=1; (c<=nc)&&(r<=nr); c++)
 {
   min = m(r,c);   rmin = r;
   for (r2=r+1; (r2<=nr)&&(min==0); r2++)
   { mr2c = m(r2,c);
     if (0!=mr2c) { min=mr2c; rmin=r2 ;}
   }
   if (min==0) npcols[++ny] = c;
   else       
     {
      pcols[++rk] = c;
      if (rmin>r) m.swaprows(r,rmin);
      for (r3 = r+1 ; r3<=nr; r3++) elimp(m,r,r3,c,pr);
      r++;
     }
 }
 for (c = rk+ny+1; c<=nc; c++) npcols[++ny] = c ;
#ifdef TRACE
       cout << "Finished first stage; rk = " << rk;
       cout << ", ny = " << ny << "\n";
       cout << "Back substitution.\n";
#endif /* TRACE */
 pcols  =  pcols.slice(1,rk);
 npcols =  npcols.slice(1,ny);    // truncate index vectors
 if (ny>0)
 { for (r=1; r<=rk; r++)
      for (r2=r+1; r2<=rk; r2++)  
          elimp(m,r2,r,pcols[r2],pr); 
   for (r=1; r<=rk; r++)
      { scalar fac = xmod(invmod(m(r,pcols[r]),pr),pr);
        for (c=1; c<=nc; c++) m(r,c)=xmodmul(fac,m(r,c),pr);
      }
 }
 else 
       for (r=1; r<=rk; r++)
        {for (c=1; c<=nc; c++) m(r,c)=(c==pcols[r]);    // 0 or 1 !
        }
 scalar modulus=pr;
 float lim=floor(sqrt(pr/2.0));
#ifdef TRACE
     cout << "Finished second stage.\n Echelon mat mod "<<pr<<" is:\n";
     cout << m;
     cout << "Now lifting back to Q.\n";
     cout << "lim = " << lim << "\n";
#endif /* TRACE */
 scalar dd = 1;
 mat nmat(rk,nc);
 mat dmat(rk,nc);

#ifdef TRACE
     cout << "rk = " << rk << "\n";
     cout << "ny = " << ny << "\n";
#endif /* TRACE */
 for (long i=1; i<=rk; i++)
 {
   for (long j=1; j<=rk; j++)
     {
       nmat(i,pcols[j])=(i==j);
       dmat(i,pcols[j])=1;
     }
   for (long j=1; j<=ny; j++)
    {
      scalar n1,d1;
      long jj = npcols[j];
      modrat(m(i,jj),modulus,lim,n1,d1);
      nmat(i,jj)=n1;
      dmat(i,jj)=d1;
      dd=(dd*d1)/gcd(dd,d1);
    }
 }
 dd=abs(dd);
#ifdef TRACE
      cout << "Numerator mat = " << nmat;
      cout << "Denominator mat = " << dmat;
      cout << "Common denominator = " << dd << "\n";
#endif /* TRACE */
      for (long i=1; i<=rk; i++)
        for (long j=1; j<=nc; j++)
          m(i,j)=(dd*nmat(i,j))/dmat(i,j);
 d=dd;
 return m;
}


// The following function computes the echelon form of m modulo the prime pr.

mat echmodp(const mat& entries, vec& pcols, vec& npcols, long& rk, long& ny, scalar pr)
{
// cout << "In echmodp with entries = " << entries;
 long r,c,r2,r3, nr=entries.nro, nc=entries.nco;
 mat m(nr,nc);
 scalar *mij=m.entries, *entriesij=entries.entries;
 long n=nr*nc;
 while(n--) *mij++ = xmod(*entriesij++,pr);
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0; r=1;
 for (c=1; (c<=nc)&&(r<=nr); c++)
   {
     mij=m.entries+(r-1)*nc+c-1;
     scalar min = *mij;
     long rmin = r;
     for (r2=r+1, mij+=nc; (r2<=nr)&&(min==0); r2++, mij+=nc)
       {
	 scalar mr2c = *mij;
	 if (0!=mr2c) { min=mr2c; rmin=r2 ;}
       }
     if (min==0) npcols[++ny] = c;
     else
       {
	 pcols[++rk] = c;
	 if (rmin>r) m.swaprows(r,rmin);
	 entriesij = m.entries+(r-1)*nc;
         if (min!=1)
           {
             if(min==-1)
               {
                 n=nc;
                 while(n--)
                   {
                     *entriesij = - (*entriesij);
                     entriesij++;
                   }
               }
             else
               {
                 scalar fac = xmod(invmod(min,pr),pr);
                 n=nc;
                 while(n--)
                   {
                     *entriesij = xmodmul(fac , *entriesij, pr);
                     entriesij++;
                   }
               }
           }
	 for (r3 = r+1 ; r3<=nr; r3++) elimp1(m,r,r3,c,pr);
	 // for (r3 = r+1 ; r3<=nr; r3++) elimp(m,r,r3,c,pr);
	 r++;
       }
   }
 for (c = rk+ny+1; c<=nc; c++) npcols[++ny] = c ;
 pcols  =  pcols.slice(rk);
 npcols =  npcols.slice(ny);    // truncate index vectors
 // cout << "Rank = " << rk << ".  Nullity = " << ny << ".\n";
 if (ny>0)
   { 
     for (r=1; r<=rk; r++)
       for (r2=r+1; r2<=rk; r2++)  
	 elimp(m,r2,r,pcols[r2],pr); 
     for (r=1; r<=rk; r++)
       { 
	 mij = m.entries+(r-1)*nc;
	 scalar fac = xmod(invmod(mij[pcols[r]-1],pr),pr);
	 n=nc; while(n--) {*mij = xmodmul(fac , *mij, pr); mij++;}
       }
   }
 else 
   {
     mij=m.entries;
     for (r=1; r<=rk; r++)
       for (c=1; c<=nc; c++) 
	 *mij++ = (c==pcols[r]);    // 0 or 1 !
   }
 return m.slice(rk,nc);
}

//#define TRACE_NTL_REF

#include <NTL/mat_lzz_p.h>
#ifdef TRACE_NTL_REF
#include <eclib/timer.h>
#endif

// Construct an NTL mat_lzz_p (matrix mod p) from a mat mod pr

mat_zz_p mat_zz_p_from_mat(const mat& M, scalar pr)
{
  long nr=M.nrows(), nc=M.ncols();
  long i, j;

#ifdef TRACE_NTL_REF
  cout<<"Creating an NTL mat_zz_p from a matrix with " << nr <<" rows and "<<nc<<" columns, mod "<<pr<<endl;
#endif
  // create NTL matrix copy of M:
  zz_pPush push(pr);
  mat_zz_p A(INIT_SIZE, nr, nc);
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++)
      {
        A.put(i,j, conv<zz_p>(M(i+1,j+1)));
      }
#ifdef TRACE_NTL_REF
  cout<<"--done."<<endl;
#endif
  return A;
}

// Construct a mat (scalar type same as pr) from an NTL mat_lzz_p

mat mat_from_mat_zz_p(const mat_zz_p& A, scalar pr) // type of scalar fixes return type
{
 long nr = A.NumRows(), nc = A.NumCols();
 long i, j;

#ifdef TRACE_NTL_REF
  cout<<"Creating a mat from an NTL mat_zz_p with " << nr <<" rows and "<<nc<<" columns, mod "<<pr<<endl;
#endif
 // create matrix copy of A:
 mat M(nr, nc);
 for(i=0; i<nr; i++)
   for(j=0; j<nc; j++)
     M(i+1,j+1) = mod(conv<scalar>(A.get(i,j)), pr);
#ifdef TRACE_NTL_REF
  cout<<"--done."<<endl;
#endif
 return M;
}

// compute ref of M mod pr via NTL, setting rk=rank, ny=nullity,
// pivotal columns pcols, non-pivotal columns npcols

mat ref_via_ntl(const mat& M, vec& pcols, vec& npcols,
                long& rk, long& ny, scalar pr)
{
 long nc=M.ncols();
 long i, j, k;
#ifdef TRACE_NTL_REF
 timer ntl_timer;
 ntl_timer.start();
#endif
 zz_pPush push(pr);
 mat_zz_p A = mat_zz_p_from_mat(M, pr);

#ifdef TRACE_NTL_REF
 cout<<"--calling NTL's gauss()..."<<flush;
#endif
 rk = gauss(A); // reduce to echelon form in place; rk is the rank
#ifdef TRACE_NTL_REF
 cout<<"done." << endl;
#endif
 ny = nc-rk;
#ifdef TRACE_NTL_REF
 cout<<"Rank = " << rk <<", nullity = "<<ny<<endl;
#endif

 // Find pivots, rescale rows so pivots are 1

 pcols.init(rk);
 npcols.init(ny);
 zz_p zero = conv<zz_p>(0);
 zz_p one = conv<zz_p>(1);
 zz_p piv, inv_piv;

 for (i = j = k = 0; i < rk; i++)
   {
     while (A.get(i,j) == zero)
       {
         npcols[k+1] = j+1;
         k++;
         j++;
       }
     piv = A.get(i,j);
     pcols[i+1] = j+1;
     j++;
     if (piv != one)
       {
         inv(inv_piv, piv);
         A[i] = inv_piv*A[i];
       }
   }
 while (k < ny)
   {
     npcols[k+1] = j+1;
     k++;
     j++;
   }

 // copy back to a new matrix for return:
 mat ans = mat_from_mat_zz_p(A, pr).slice(rk,nc);
#ifdef TRACE_NTL_REF
 ntl_timer.start();
 ntl_timer.show();
 cout<<endl;
#endif
 return ans;
}

long rank_via_ntl(const mat& M, scalar pr)
{
#ifdef TRACE_NTL_REF
  cout << "Computing rank mod "<<pr<<" of a matrix of size ("<<M.nrows()<<", "<<M.ncols()<<")..."<<flush;
  timer ntl_timer;
 ntl_timer.start();
#endif
 zz_pPush push(pr);
 mat_zz_p A = mat_zz_p_from_mat(M, pr);
 long rk = gauss(A); // reduce to echelon form in place; rk is the rank
#ifdef TRACE_NTL_REF
 cout << "done: "<<flush;
 ntl_timer.start();
 ntl_timer.show();
 cout<<endl;
#endif
 return rk;
}

long det_via_ntl(const mat& M, scalar pr)
{
#ifdef TRACE_NTL_REF
  cout << "Computing determinant mod "<<pr<<" of a matrix of size ("<<M.nrows()<<", "<<M.ncols()<<")..."<<flush;
  timer ntl_timer;
 ntl_timer.start();
#endif
 zz_pPush push(pr);
 mat_zz_p A = mat_zz_p_from_mat(M, pr);
 zz_p det = determinant(A);
#ifdef TRACE_NTL_REF
 cout << "done: "<<flush;
 ntl_timer.start();
 ntl_timer.show();
 cout<<endl;
#endif
 return mod(conv<scalar>(det), pr);
}

#if FLINT

#include "eclib/flinterface.h"

// FLINT has more than one type for modular matrices: standard in
// FLINT-2.3..2.9 was nmod_mat_t with entries of type mp_limb_t
// (unsigned long) while non-standard was hmod_mat_t, with entries
// hlimb_t (unsigned int).  From FLINT-3 the latter is emulated via a
// wrapper.  We use the former when scalar=long and the latter when
// scalar=int and the FLINT versin is at least 3.  The unsigned
// scalar types are #define'd as uscalar.

void mod_mat_from_mat(mod_mat& A, const mat& M, scalar pr)
{
  long nr=M.nrows(), nc=M.ncols();
  long i, j;

  // copy of the modulus for FLINT
  uscalar mod = (uscalar)pr;

  // create flint matrix copy of M:
  mod_mat_init(A, nr, nc, mod);
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++)
      mod_mat_entry(A,i,j) = (uscalar)posmod(M(i+1,j+1),pr);
}

mat mat_from_mod_mat(const mod_mat& A, scalar a) // scalar just to fix return type
{
 long nr=mod_mat_nrows(A), nc=mod_mat_ncols(A);

  // create matrix copy of A:
  mat M(nr, nc);
  long i, j;
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++)
      M(i+1,j+1) = mod_mat_entry(A,i,j);
  return M;
}

mat ref_via_flint(const mat& M, scalar pr)
{
  // create flint matrix copy of M:
  mod_mat A;
  mod_mat_from_mat(A,M,pr);

  // reduce A to rref:
#ifdef TRACE_FLINT_RREF
  timeit_t t;
  timeit_start(t);
  long nc=M.ncols(), nr=mod_mat_nrows(A);
  cerr<<"(nr,nc)=("<<nr<<","<<nc<<"): "<<flush;
#endif
  mod_mat_rref(A);
#ifdef TRACE_FLINT_RREF
  timeit_stop(t);
  cerr<<" cpu = "<<(t->cpu)<<" ms, wall = "<<(t->wall)<<" ms"<<endl;
#endif

  // copy back to a new matrix for return:
  mat ans = mat_from_mod_mat(A, pr);

  // clear the flint matrix and return:
  mod_mat_clear(A);
  return ans;
}

// The following function computes the reduced echelon form
// of M modulo the prime pr, calling FLINT's nmod_mat_rref function.

mat ref_via_flint(const mat& M, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar pr)
{
  long nc=M.ncols();
  long i, j, k;

#ifdef TRACE_FLINT_RREF
#if (SCALAR_OPTION==1)
  cout << "In ref_via_flint(M) with M having "<<nr<<" rows and "<<nc<<" columns, using hmod_mat and modulus "<<pr<<"."<<endl;
#else
  cout << "In ref_via_flint(M) with M having "<<nr<<" rows and "<<nc<<" columns, using nmod_mat and modulus "<<pr<<"."<<endl;
#endif
  //  cout << "Size of  scalar = "<<8*sizeof(scalar)<<" bits"<<endl;
  //  cout << "Size of uscalar = "<<8*sizeof(uscalar)<<" bits"<<endl;
#endif

  // create flint matrix copy of M:
  mod_mat A;
  mod_mat_from_mat(A,M,pr);

#ifdef TRACE_FLINT_RREF
  timeit_t t;
  timeit_start(t);
  long nr=M.nrows();
  cerr<<"(nr,nc)=("<<nr<<","<<nc<<"): "<<flush;
#endif

  // reduce A to rref:
  rk = mod_mat_rref(A);
#ifdef TRACE_FLINT_RREF
  timeit_stop(t);
  cerr<<"rank = "<<rk<<". cpu = "<<(t->cpu)<<" ms, wall = "<<(t->wall)<<" ms"<<endl;
#endif

  // construct vectors of pivotal and non-pivotal columns
  ny = nc-rk;
  pcols.init(rk);
  npcols.init(ny);
  for (i = j = k = 0; i < rk; i++)
    {
      while (mod_mat_entry(A, i, j) == 0UL)
        {
          npcols[k+1] = j+1;
          k++;
          j++;
        }
      pcols[i+1] = j+1;
      j++;
    }
  while (k < ny)
    {
      npcols[k+1] = j+1;
      k++;
      j++;
    }

  // copy back to a new matrix for return:
  mat ans = mat_from_mod_mat(A,pr).slice(rk,nc);

  // clear the flint matrix and return:
  mod_mat_clear(A);
  return ans;
}
#endif // FLINT

// The following function computes the upper-triangular echelon form
// of m modulo the prime pr.

mat echmodp_uptri(const mat& entries, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar pr)
{
// cout << "In echmodp_uptri with matrix = " << entries; 
long nc,nr,r,c,r2,r3;
 scalar min, mr2c;
 nr=entries.nro, nc=entries.nco;
 mat m(nr,nc);
 scalar *mij=m.entries, *entriesij=entries.entries;
 long n=nr*nc;
 while(n--) *mij++ = xmod(*entriesij++,pr);
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0; r=1;
 for (c=1; (c<=nc)&&(r<=nr); c++)
   {
     mij=m.entries+(r-1)*nc+c-1;
     min = *mij;
     long rmin = r;
     for (r2=r+1, mij+=nc; (r2<=nr)&&(min==0); r2++, mij+=nc)
       {
	 mr2c = *mij;
	 if (0!=mr2c) { min=mr2c; rmin=r2 ;}
       }
     if (min==0) npcols[++ny] = c;
     else
       {
	 pcols[++rk] = c;
	 if (rmin>r) m.swaprows(r,rmin);
	 entriesij = m.entries+(r-1)*nc;
         scalar fac;
         if (min!=1)
           {
             if(min==-1)
               {
                 n=nc; while(n--) {*entriesij = - (*entriesij); entriesij++;}
               }
             else
               {
                 //                 cout<<"pivot = "<<min<<endl;
                 fac = xmod(invmod(min,pr),pr);
                 n=nc; while(n--) {*entriesij = xmodmul(fac , *entriesij, pr); entriesij++;}
               }
           }
	 for (r3 = r+1 ; r3<=nr; r3++) elimp1(m,r,r3,c,pr);
	 r++;
       }
   }
 for (c = rk+ny+1; c<=nc; c++) npcols[++ny] = c ;
 pcols  =  pcols.slice(rk);
 npcols =  npcols.slice(ny);    // truncate index vectors
 // cout << "Rank = " << rk << ".  Nullity = " << ny << ".\n";
 return m.slice(rk,nc);
}

mat matmulmodp(const mat& m1, const mat& m2, scalar pr)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 mat m3(m,p);
 scalar *a=m1.entries, *b=m2.entries, *c=m3.entries;
 if (n==m2.nro)  // algorithm from Dr Dobb's Journal August 1993
   {
     while(m--)
       {
	 scalar *bp=b;
         long k=n;
	 while(k--)
	   {
	     scalar *cp=c;
             long j=p;
	     while(j--)
	       {
		 *cp += xmodmul(*a , *bp++, pr);
		 *cp = xmod(*cp,pr);
		 cp++;
	       }
	     a++;
	   }
	 c += p;
       }
   }
 else
   {
     cerr << "Incompatible sizes in mat product"<<endl;
   }
 return m3;
}

int liftmat(const mat& mm, scalar pr, mat& m, scalar& dd, int trace)
{
  scalar modulus=pr,n,d; long nr,nc,nrc; dd=1;
  int success=1;
  float lim=floor(sqrt(pr/2.0));
  m = mm; scalar *mp;
  if(trace)
    {
      cout << "Lifting mod-p mat;  mat mod "<<pr<<" is:\n";
      cout << m;
      cout << "Now lifting back to Q.\n";
      cout << "lim = " << lim << "\n";
    }
  nr = m.nro; nc = m.nco;  nrc = nr*nc; mp=m.entries;
  while(nrc--)
    {
      int succ = modrat(*mp++,modulus,lim,n,d);
      success = success && succ;
      dd=lcm(d,dd);
    }
  if(!success)
    {
      //cout << "Problems encountered with modrat lifting of mat." << endl;
      return 0;
    }
   dd=abs(dd);
  if(trace) cout << "Common denominator = " << dd << "\n";
  nrc=nr*nc; mp=m.entries;
  while(nrc--)
      {
	*mp=mod(xmodmul(dd,(*mp),pr),pr);
	mp++;
      }
  return 1;
}

double sparsity(const mat& m)
{
    long nr=m.nrows(), nc=m.ncols();
    if(nr==0) return 1;
    if(nc==0) return 1;
    scalar* matij=m.entries;
    double count=0;
    long n=nr*nc;
    long i=n;
    while(i--) {if(*matij++) count+=1;}
    return count/n;
}

#if (FLINT==1)&&(__FLINT_VERSION>2)&&(SCALAR_OPTION==1)

// Implementation of wrapper functions declared in flinterface.h
// written by Fredrik Johansson

#include <flint/gr.h>
#include <flint/gr_mat.h>

void
hmod_mat_init(hmod_mat_t mat, slong rows, slong cols, hlimb_t n)
{
    gr_ctx_t ctx;
    gr_ctx_init_nmod32(ctx, n);
    gr_mat_init((gr_mat_struct *) mat, rows, cols, ctx);
    nmod_init(&(mat->mod), n);
}

void
hmod_mat_clear(hmod_mat_t mat)
{
    if (mat->entries)
    {
        flint_free(mat->entries);
        flint_free(mat->rows);
    }
}

void
hmod_mat_mul(hmod_mat_t C, const hmod_mat_t A, const hmod_mat_t B)
{
    gr_ctx_t ctx;
    gr_ctx_init_nmod32(ctx, C->mod.n);
    GR_MUST_SUCCEED(gr_mat_mul((gr_mat_struct *) C, (gr_mat_struct *) A, (gr_mat_struct *) B, ctx));
}

slong
hmod_mat_rref(hmod_mat_t mat)
{
    slong rank;
    gr_ctx_t ctx;
    gr_ctx_init_nmod32(ctx, mat->mod.n);
    GR_MUST_SUCCEED(gr_mat_rref_lu(&rank, (gr_mat_struct *) mat, (gr_mat_struct *) mat, ctx));
    return rank;
}

#endif
