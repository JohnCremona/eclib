// mat.cc: implementation of integer matrix classes
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
 
// Only to be included by matrix.cc

#define LONGLONG scalar   // long longs caused memory alocation bugs
                        // N.B. if = long long, uncomment abs() below

// Definitions of member operators and functions:

mat::mat(long nr, long nc)
{
 nro=nr;
 nco=nc;
 long n=nr*nc;
 entries=new scalar[n];  if (!entries) {cerr<<"Out of memory!\n"; abort();}
 scalar* m1=entries;
 while(n--) *m1++ = 0;
}

mat::mat(const mat& m)
{
 nro=m.nro;
 nco=m.nco;
 long n = nro*nco;
 entries=new scalar[n]; if (!entries) {cerr<<"Out of memory!\n"; abort();}
 scalar *m1=entries, *m2=m.entries;
 while(n--) *m1++ = *m2++;
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
     if (!entries) {cerr<<"Out of memory!\n"; abort();}
   }
 nro = nr;
 nco = nc;
 scalar *m1=entries;
 while(n--) *m1++ = 0;
}

scalar& mat::operator()(long i, long j)  const   // returns ref to (i,j) entry
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) 
   return entries[(i-1)*nco+(j-1)];
 else 
   {
     cerr << "Bad indices in mat::sub\n"; 
     return entries[0];
   }
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
     c=cc;
     while(c--) *ap++ = *mp++; 
     mp+=(nco-cc);
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
     if (!entries) {cerr<<"Out of memory!\n"; abort();}
   }
 nro = m.nro;
 nco = m.nco;
 scalar *m1=entries, *m2=m.entries;
 while(n--) *m1++ = *m2++;
 return *this;
}

scalar mat::sub(long i, long j) const
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) return entries[(i-1)*nco+(j-1)];
  else {cerr << "Bad indices in mat::sub\n"; return 0;}
}

void mat::set(long i, long j, scalar x)
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) entries[(i-1)*nco+(j-1)] = x;
  else cerr << "Bad indices in mat::set\n";
}

void mat::add(long i, long j, scalar x)
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) entries[(i-1)*nco+(j-1)] += x;
  else cerr << "Bad indices in mat::add\n";
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
 else cerr << "Bad indices in mat::setrow (i="<<i<<", nro="<<nro<<", dim(v)="<<dim(v)<<", nco="<<nco<<")\n";
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
 else cerr << "Bad indices in mat::setcol\n";
}

vec mat::row(long i) const
{
 vec mi(nco);
 long j=nco; scalar *entriesij=entries+(i-1)*nco, *v=mi.entries;
 if ((0<i) && (i<=nro)) 
   while(j--) *v++ = *entriesij++;
 else 
   cerr << "Bad row number in function mat::row\n";
 return mi;
}

vec mat::col(long j) const
{
 vec mj(nro);
 long i=nro; scalar *entriesij=entries+(j-1), *v=mj.entries;
 if ((0<j) && (j<=nco))
   while(i--) {*v++ = *entriesij; entriesij+=nco;}
 else 
   cerr << "Bad column number in function mat::col\n";
 return mj;
}

void mat::swaprows(long r1, long r2)
{
  if ((0<r1)&&(0<r2)&&(r1<=nro)&&(r2<=nro))
    {
      scalar *mr1 = entries + (r1-1)*nco;
      scalar *mr2 = entries + (r2-1)*nco;
      long nc=nco, a;
      while(nc--) {a = *mr1; *mr1++ = *mr2; *mr2++ = a; }
    }
  else cerr << "Bad row numbers " << r1 << "," << r2 << " in swaprow\n";
}

void mat::multrow(long r, scalar scal)
{
  if ((0<r)&&(r<=nro))
    {
      long nc=nco; scalar *mij = entries+(r-1)*nco;
      while(nc--) (*mij++) *= scal;
    }
  else cerr << "Bad row number " << r << " in multrow\n";
}

void mat::divrow(long r, scalar scal)
{
  if ((0<r)&&(r<=nro))
    {
      long nc=nco; scalar *mij = entries+(r-1)*nco;
      while(nc--) (*mij++) /= scal;
    }
  else cerr << "Bad row number " << r << " in divrow\n";
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
  else cerr << "Bad row number " << r << " in clearrow\n";
}

mat& mat::operator+=(const mat& entries2)
{
  if ((nro==entries2.nro) && (nco=entries2.nco)) 
    {
      long n=nro*nco; scalar *m1=entries, *m2=entries2.entries;
      while(n--) (*m1++) += (*m2++);
    }
  else cerr << "Incompatible matrices in operator +=\n";
  return *this;
}

mat& mat::operator-=(const mat& entries2)
{
  if ((nro==entries2.nro) && (nco=entries2.nco)) 
    {
      long n=nro*nco; scalar *m1=entries, *m2=entries2.entries;
      while(n--) (*m1++) -= (*m2++);
    }
  else cerr << "Incompatible matrices in operator -=\n";
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

long nrows(const mat& m) {return m.nro;}
long ncols(const mat& m) {return m.nco;}

// add/sub row i of mat to v (implemented in mat.cc)
void add_row_to_vec(vec& v, const mat& m, long i)
{
  scalar* vi=v.entries, *mij=m.entries+(i-1)*m.nco; long j=v.d;
  while(j--)(*vi++)+=(*mij++);
}

void sub_row_to_vec(vec& v, const mat& m, long i)
{
  scalar* vi=v.entries, *mij=m.entries+(i-1)*m.nco; long j=v.d;
  while(j--)(*vi++)-=(*mij++);
}

mat operator*(const mat& m1, const mat& m2)
{
 long j,k, m=m1.nro, n=m1.nco, p=m2.nco;
 mat m3(m,p); 
 scalar *a=m1.entries, *b=m2.entries, *c=m3.entries, *bp, *cp;
 if (n==m2.nro)  // algorithm from Dr Dobb's Journal August 1993
   {
     while(m--)
       {
	 bp=b; k=n;
	 while(k--)
	   {
	     cp=c; j=p;
	     while(j--)
	       {
		 *cp++ += *a * *bp++;
	       }
	     a++;
	   }
	 c += p;
       }
   }
 else cerr << "Incompatible sizes in mat product\n";
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
  long nc,nr=nro;
  while(nr--)
    {
      nc=nco;
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
  long nc,nr=nro;
  while(nr--)
    {
      nc=nco;
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
  scalar ma, mi;
  int* colwidths = new int[nc];
  for(j=0; j<nco; j++)
    {
      ma=mi=0;
      mij = entries+j;
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

//binary file I/O

// Not supported by libstdc++-v3, it seems

//#ifdef THIS_USED_TO_WORK_WITH_GCC_2_95
#if(1)

void mat::dump_to_file(char* filename) const
{
  ofstream fout(filename,ofstream::binary);
  fout.write((char*)&nro,sizeof(nro));
  fout.write((char*)&nco,sizeof(nco));
  fout.write((char*)entries,nro*nco*sizeof(scalar));
  fout.close();
}

void mat::read_from_file(char* filename)
{
  ifstream fin(filename);
  fin.read((char*)&nro,sizeof(nro));
  fin.read((char*)&nco,sizeof(nco));
  delete[] entries;
  entries = new scalar[nro*nco]; 
  fin.read((char*)entries,nro*nco*sizeof(scalar));
  fin.close();
}

#else

void mat::dump_to_file(char* filename) const
{
  ofstream fout(filename);
  fout<<nro<<" "<<nco<<" ";
  long size=nro*nco; 
  scalar* m=entries;
  while(size--) fout<<*m++<<" ";
  fout.close();
}

void mat::read_from_file(char* filename)
{
  ifstream fin(filename);
  fin>>nro>>nco;
  long size=nro*nco; 
  delete[] entries;
  entries = new scalar[size]; 
  scalar* m=entries;
  while(size--) fin>>*m++;
  fin.close();
}

#endif

istream& operator>>(istream& s, mat& m)
{
 long n=m.nro*m.nco;
 scalar* mij=m.entries;
 while(n--) s >> (*mij++);
 return s;
}

mat colcat(const mat& a, const mat& b)
{
 long nc, nr = a.nro, nca = a.nco, ncb = b.nco;
 mat ans(nr,nca+ncb);
 scalar *ansij=ans.entries, *aij=a.entries, *bij=b.entries;
 if (nr==b.nro)
   while(nr--)
     {
       nc=nca; while(nc--) *ansij++ = *aij++;
       nc=ncb; while(nc--) *ansij++ = *bij++;
     }
 else cerr << "colcat: matrices have different number of rows!" << "\n";
 return ans;
}

mat rowcat(const mat& a, const mat& b)
{
 long n, nra = a.nro, nc = a.nco, nrb = b.nro;
 mat ans(nra+nrb,nc);
 scalar *ansij=ans.entries, *aij=a.entries, *bij=b.entries;
 if (nc==b.nco)
 {
   n = nra*nc; while(n--) *ansij++ = *aij++;
   n = nrb*nc; while(n--) *ansij++ = *bij++;
 }
 else cerr << "rowcat: matrices have different number of columns!" << "\n";
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
 long r=m.nro, c=m.nco; scalar *mp,*vp,*wp;
 vec w(r);
 if (c==v.d)
   {
     mp=m.entries; wp=w.entries;
     while(r--)
       {
	 vp=v.entries; c=m.nco;
	 while(c--) *wp += (*mp++)*(*vp++);
	 wp++;
       }
   }
 else cerr << "Incompatible sizes in *(mat,vec)\n";
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
 nr=ncols(m); nc=nrows(m);
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
     case 2: return echelonp(entries,pcols,npcols,rk,ny,d,BIGPRIME);
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

//  Original version 

// mat echelon0(const mat& entries, vec& pcols, vec& npcols,
//                long& rk, long& ny, long& d)
// {
//  long nc, nr;
//  long r,c,r2,r3,rmin;
//  long min, mr2c,lastpivot;
//  rk=0; ny=0; r=1; lastpivot=1;
// 
//  mat m(entries);
//  nc=ncols(m); nr=nrows(m);
//  pcols.init(nc);
//  npcols.init(nc);
//  for (c=1; (c<=nc)&&(r<=nr); c++)
//   {min = abs(m(r,c));
//    rmin = r;
//    for (r2=r+1; (r2<=nr)&&(min!=1); r2++)
//    { mr2c = abs(m(r2,c));
//      if ((0<mr2c) && ((mr2c<min) || (min==0))) { min=mr2c; rmin=r2 ;}
//    }
//    if (min==0) npcols[++ny] = c;
//    else
//      {pcols[++rk] = c;
//       if (rmin>r) m.swaprows(r,rmin);
//       for (r3 = r+1 ; r3<=nr; r3++)
//          elimrows2(m,r,r3,c,lastpivot);
//       lastpivot=min;
//       r++;
//      }
//  }
//  for (c = rk+ny+1; c<=nc; c++) npcols[++ny] = c;
//  pcols  =  pcols.slice(rk);
//  npcols =  npcols.slice(ny);    // truncate index vectors
// // cout << "In echelon:\n";
// // cout << "pcols = " << pcols << "\n";
// // cout << "npcols = " << npcols << "\n";
//  d=1;
//  lastpivot=1;
//  if (ny>0)
//   {for (r=1; r<=rk; r++)  m.clearrow(r);
//    for (r=1;r<=rk; r++)
//       for (r2=r+1; r2<=rk; r2++)
//         elimrows1(m,r2,r,pcols[r2]);
//    for (r=1; r<=rk; r++)
//       d = (d*m(r,pcols[r]))/gcd(d,m(r,pcols[r]));
//    d = abs(d);
//    for (r=1; r<=rk; r++)
//       {long fac = d/m(r,pcols[r]);
//        m.multrow(r,fac);
//       }
//  }
//  else 
//        for (r=1; r<=rk; r++)
//         for (c=1; c<=nc; c++)
//           m.set(r,c, (c==pcols[r]));  // 0 or 1 !
//  return m.slice(rk,nc);
// }

long rank(const mat& entries)
{
 long rk,nr,nc,r,c,r2,r3,rmin;
 long min, mr2c,lastpivot;
 rk=0; r=1; lastpivot=1;
 mat m(entries);
 nc=ncols(m); nr=nrows(m);
 for (c=1; (c<=nc)&&(r<=nr); c++)
 { min = abs(m(r,c));
   rmin = r;
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

long nullity(const mat& m)
{
 return ncols(m)-rank(m);
}

long trace(const mat& a)
{ long i; long ans=0;
  for (i=1; i<=nrows(a); i++) ans += a(i,i);
  return ans;
}
 
/*    OLD VERSION
vector<long> charpoly(const mat& a)
{ long i,k,r,n = nrows(a);
  vec tlist = vec(n);
  mat apower(a);
  vector<long> clist(n+1);
  tlist[1] = trace(a);
  for (i=2; i<=n; i++)
      { apower*=a;
        tlist[i] = trace(apower);
      }
   clist[n]=1;
   for (k=1; k<=n; k++)
      { long temp = 0;
         for (r=1; r<=k; r++)  temp+= tlist[r]*clist[n+r-k];
         clist[n-k]= -temp/k;
      }
  return clist;
}
*/

// NEW VERSION -- FADEEV'S METHOD

vector<long> charpoly(const mat& a)
{ long n = nrows(a);
  mat b(a);
  mat id(idmat((scalar)n));
  vector<long> clist(n+1);
  long t = trace(a);
  clist[n]   =  1;
  clist[n-1] = -t;
  for (long i=2; i<=n; i++)
      { b=a*(b-t*id);          //     cout << b;   // (for testing only)
        t=trace(b)/i;
        clist[n-i] = -t;
      }
  if (!(b==t*id)) cerr << "Error in charpoly: final b = " << (b-t*id);
  return clist;
}

 
long determinant(const mat& m)
{
 vector<long> cp = charpoly(m);
 long det = cp[0];
 if (nrows(m)%2==1) det=-det;
 return det;
}

mat addscalar(const mat& m, scalar c)
{
  return m+(c*idmat((scalar)nrows(m)));
}
 
vec apply(const mat& m, const vec& v)    // same as *(mat, vec)
{
 long nr=nrows(m), nc=ncols(m);
 vec ans(nr);
 if (nc==dim(v))
   for (long i=1; i<=nr; i++) 
     ans[i] = m.row(i)*v;
 else cerr << "Incompatible sizes in *(mat,vec)\n";
 return ans;
}

/*  Need this when LONGLONG = long long!
LONGLONG abs(LONGLONG a)
{return ((a<0) ? -a : a);
}
*/

LONGLONG lgcd(LONGLONG aa, LONGLONG bb)
{LONGLONG a,b,c;
 a=aa; b=bb;
 while (b!=0) {c=a%b; a=b; b=c;}
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
  long nr, nc, r,c,r2,r3,rmin,i;
  long min, mr2c,lastpivot;
  rk=0; ny=0; r=0; lastpivot=1;
  nc=entries.nco; nr=entries.nro;
  LONGLONG *m, *mi1, *mi2, *mij; LONGLONG temp;
  m = new LONGLONG[nr*nc];
  long n=nr*nc; mij=m; mi1=entries.entries;
  while(n--) *mij++ = (LONGLONG)(*mi1++);

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
    {
      temp = *mij++; 
      if ((INT_MIN<=temp)&&(temp<=INT_MAX))
	{ 
	  *ansij++=temp;
	}
      else 
	{
	  cerr << "Problem in echelonl: entry " << temp << " too big!\n";
	  *ansij++=0;
	}
    }
  
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
 scalar *mr1 = m.entries + (r1-1)*nc,
      *mr2 = m.entries + (r2-1)*nc;
 scalar p = mod(mr1[pos-1],pr), q=mod(mr2[pos-1],pr);
 if(p==1)
   if(q==0) {;} // nothing to do
   else
     if(q==1)
       while(nc--)
	 {
	   (*mr2)= mod(*mr2-*mr1,pr);
	   mr1++; mr2++;
	 }
     else // general q
       while(nc--)
	 {
	   (*mr2)= mod(*mr2-xmodmul(q,*mr1,pr),pr);
	   mr1++; mr2++;
	 }
 else // general p (p!=1)
   if(q==0) {;} // nothing to do
   else 
     if(q==1)
     while(nc--)
       {
	 (*mr2)= mod(xmodmul(p,*mr2,pr)-*mr1,pr);
	 mr1++; mr2++;
       }
     else // general q
       while(nc--)
	 {
	   (*mr2)= mod(xmodmul(p,*mr2,pr)-xmodmul(q,*mr1,pr),pr);
	   mr1++; mr2++;
	 }
}

void elimp1(const mat& m, long r1, long r2, long pos, scalar pr)
//same as elimp except assumes pivot is 1
{
 long nc=m.nco;
 scalar *mr1 = m.entries + (r1-1)*nc,
      *mr2 = m.entries + (r2-1)*nc;
 scalar q=mod(mr2[pos-1],pr);
 if(q==0) {;}
 else
   if(q==1)
   while(nc--)
     {
       (*mr2)= mod(*mr2-*mr1,pr);
       mr1++; mr2++;
     }
   else  // general q
     while(nc--)
       {
	 (*mr2)= mod((*mr2)-xmodmul(q,*mr1,pr),pr);
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
 scalar min, mr2c,lastpivot;
 nr=nrows(entries), nc=ncols(entries);
 mat m(nr,nc);
 scalar *mij=m.entries, *entriesij=entries.entries;
 long n=nr*nc;
 while(n--) *mij++ = xmod(*entriesij++,pr);
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0; r=1; lastpivot=1;
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
 long i,j;
 for (i=1; i<=rk; i++) 
 {  
    for (j=1; j<=rk; j++) 
     {
//    nmat.set(i,pcols[j],(i==j));
      nmat(i,pcols[j])=(i==j);
//    dmat.set(i,pcols[j],1);
      dmat(i,pcols[j])=1;
     }
    for (j=1; j<=ny; j++)
    {scalar n,d;
     long jj = npcols[j];
     modrat(m(i,jj),modulus,lim,n,d);
//   nmat.set(i,jj,n);
     nmat(i,jj)=n;
//   dmat.set(i,jj,d);
     dmat(i,jj)=d;
     dd=(dd*d)/gcd(dd,d);
    }
 }
 dd=abs(dd);
#ifdef TRACE
      cout << "Numerator mat = " << nmat;
      cout << "Denominator mat = " << dmat;
      cout << "Common denominator = " << dd << "\n";
#endif /* TRACE */
 for (i=1; i<=rk; i++)
 {
  for (j=1; j<=nc; j++) m(i,j)=(dd*nmat(i,j))/dmat(i,j);
 }
 d=dd;
 return m; 
}


// The following function computes the echelon form of m modulo the prime pr.

mat echmodp(const mat& entries, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar pr)
{
// cout << "In echmodp with entries = " << entries;
 long nc,nr,r,c,r2,r3,rmin;
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
     min = *mij;   rmin = r;
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
	 for (r3 = r+1 ; r3<=nr; r3++) elimp(m,r,r3,c,pr);
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

// The following function computes the upper-triangular echelon form of m modulo the prime pr.

mat echmodp_uptri(const mat& entries, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar pr)
{
// cout << "In echmodp_uptri with entries = " << entries;
 long nc,nr,r,c,r2,r3,rmin;
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
     min = *mij;   rmin = r;
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
	 scalar fac = xmod(invmod(min,pr),pr);
	 n=nc; while(n--) {*entriesij = xmodmul(fac , *entriesij, pr); entriesij++;}
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
 long j,k, m=m1.nro, n=m1.nco, p=m2.nco;
 mat m3(m,p); 
 scalar *a=m1.entries, *b=m2.entries, *c=m3.entries, *bp, *cp;
 if (n==m2.nro)  // algorithm from Dr Dobb's Journal August 1993
   {
     while(m--)
       {
	 bp=b; k=n;
	 while(k--)
	   {
	     cp=c; j=p;
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
 else cerr << "Incompatible sizes in mat product\n";
 return m3;
}

mat liftmat(const mat& mm, scalar pr, scalar& dd, int trace)
{
  scalar modulus=pr,n,d; long nr,nc,nrc; dd=1;
  int succ,success=1;
  float lim=floor(sqrt(pr/2.0));
  mat m = mm; scalar *mp;
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
      succ = modrat(*mp++,modulus,lim,n,d);
      success = success && succ;
      dd=lcm(d,dd);
    }
  if(!success) 
    cout << "Problems encountered with modrat lifting of mat." << endl;
  dd=abs(dd);
  if(trace) cout << "Common denominator = " << dd << "\n";
  nrc=nr*nc; mp=m.entries;
  while(nrc--)
      {
	*mp=mod(xmodmul(dd,(*mp),pr),pr); 
	mp++;
      }
  return m;
}

double sparsity(const mat& m)
{
    long nr=nrows(m), nc=ncols(m);
    if(nr==0) return 1;
    if(nc==0) return 1;
    scalar* matij=m.entries;
    double count=0;
    long n=nr*nc;
    long i=n;
    while(i--) {if(*matij++) count+=1;}
    return count/n;
}
