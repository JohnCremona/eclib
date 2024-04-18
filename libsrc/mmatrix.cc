// mmatrix.cc: implementation of multiprecision integer matrix class
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
 
#include <eclib/marith.h>
#include <eclib/mmatrix.h>

const bigint MBIGPRIME=atoI(string("6074000003").c_str());
// will convert this string to an bigint
//This is nearly the largest p such that (p/2)^2 < 2^63.

// Definitions of member operators and functions:

mat_m::mat_m(long nr, long nc)
{
 nro=nr;
 nco=nc;
 long n=nr*nc;
 entries=new bigint[n];
 if (!entries)
   cerr<<"Out of memory in mat_m constructor"<<endl;
 else
   {
     bigint* m1=entries;
     while(n--) *m1++ = 0;
   }
}

mat_m::mat_m(const mat_m& m)
{
 nro=m.nro;
 nco=m.nco;
 long n = nro*nco;
 entries=new bigint[n];
 if (!entries)
   cerr<<"Out of memory in mat_m constructor"<<endl;
 else
   {
     bigint *m1=entries, *m2=m.entries;
     while(n--) *m1++ = *m2++;
   }
}

mat_m::mat_m(const mat_i& m)
{
 nro=m.nro;
 nco=m.nco;
 long n = nro*nco;
 entries=new bigint[n];
 if (!entries)
   cerr<<"Out of memory in mat_m constructor"<<endl;
 {
   bigint *m1=entries; int *m2=m.entries;
   while(n--) *m1++ = *m2++;
 }
}

mat_m::mat_m(const mat_l& m)
{
 nro=m.nro;
 nco=m.nco;
 long n = nro*nco;
 entries=new bigint[n];
 if (!entries)
   cerr<<"Out of memory in mat_m constuctor"<<endl;
 {
   bigint *m1=entries; long *m2=m.entries;
   while(n--) *m1++ = *m2++;
 }
}

mat_m::~mat_m()
{
  delete[] entries;
}

void mat_m::init(long nr, long nc) // assigns to zero matrix of given size;
{                                 // with defaults (0,0) releases all space.
 long n = nr*nc;
 if (nro*nco!=n) // delete old space; 
   {            // replace with new.
     delete[] entries;
     entries = new bigint[n]; 
   }
 if (!entries)
   cerr<<"Out of memory in mat_m::init"<<endl;
 else
   {
     nro = nr;
     nco = nc;
     bigint *m1=entries;
     while(n--) *m1++ = 0;
   }
}

bigint& mat_m::operator()(long i, long j)  const   // returns ref to (i,j) entry
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) 
   return entries[(i-1)*nco+(j-1)];
 else
   {
     cerr << "Bad indices in mat_m::sub"<<endl;
     return entries[0];
   }
}


mat_m mat_m::slice(long r1,long r2,long c1,long c2) const
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
 mat_m ans(n,cc);
 bigint* ap=ans.entries, *mp=entries+r1*nco+c1;
 while(n--)
   {
     c=cc;
     while(c--) *ap++ = *mp++; 
     mp+=(nco-cc);
   }
 return ans;
}

mat_m& mat_m::operator=(const mat_m& m)
{
 if (this==&m) return *this;
 long n = m.nro*m.nco;
 if (nro*nco!=n) // delete old space; 
   {            // replace with new.
     delete[] entries;
     entries = new bigint[n]; 
   }
 if (!entries)
   cerr<<"Out of memory in mat_m assignment"<<endl;
 else
   {
     nro = m.nro;
     nco = m.nco;
     bigint *m1=entries, *m2=m.entries;
     while(n--) *m1++ = *m2++;
   }
 return *this;
}

bigint mat_m::sub(long i, long j) const
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) return entries[(i-1)*nco+(j-1)];
 else {cerr << "Bad indices in mat_m::sub"<<endl; bigint ans; return ans;}
}

void mat_m::set(long i, long j, const bigint& x)
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) entries[(i-1)*nco+(j-1)] = x;
 else {cerr << "Bad indices in mat_m::set"<<endl;}
}

void mat_m::add(long i, long j, const bigint& x)
{
 if ((0<i) && (i<=nro) && (0<j) && (j<=nco)) entries[(i-1)*nco+(j-1)] += x;
 else {cerr << "Bad indices in mat_m::add"<<endl;}
}

void mat_m::setrow(long i, const vec_m& v)
{
 if ((0<i) && (i<=nro) && (dim(v)==nco))
  {
    bigint * rowi = entries + (i-1)*nco;
    bigint * vec = v.entries;
    long c=nco;
    while(c--) *rowi++ = *vec++;
  }
 else {cerr << "Bad indices in mat_m::setrow"<<endl;}
}

void mat_m::setcol(long j, const vec_m& v)
{
 if ((0<j) && (j<=nco) && (dim(v)==nro))
  {
   bigint * colj = entries+(j-1);
   bigint * vec = v.entries;
   long n=nro;
   while(n--) {*colj = *vec++; colj+=nco;}
 }
 else {cerr << "Bad indices in mat_m::setcol"<<endl;}
}

vec_m mat_m::row(long i) const
{
 vec_m mi(nco);
 long j=nco; bigint *matij=entries+(i-1)*nco, *v=mi.entries;
 if ((0<i) && (i<=nro)) 
   while(j--) *v++ = *matij++;
 else
   {
     cerr << "Bad row number in function mat_m::row"<<endl;
   }
 return mi;
}

vec_m mat_m::col(long j) const
{
 vec_m mj(nro);
 long i=nro; bigint *matij=entries+(j-1), *v=mj.entries;
 if ((0<j) && (j<=nco))
   while(i--) {*v++ = *matij; matij+=nco;}
 else
   {
     cerr << "Bad column number in function mat_m::col"<<endl;
   }
 return mj;
}

void mat_m::swaprows(long r1, long r2)
{
  if ((0<r1)&&(0<r2)&&(r1<=nro)&&(r2<=nro))
    {
      bigint *mr1 = entries + (r1-1)*nco;
      bigint *mr2 = entries + (r2-1)*nco;
      long nc=nco; bigint a;
      while(nc--) {a = *mr1; *mr1++ = *mr2; *mr2++ = a; }
    }
  else
    {
      cerr << "Bad row numbers " << r1 << "," << r2 << " in swaprow"<<endl;
    }
}

void mat_m::multrow(long r, const bigint& scal)
{
  if ((0<r)&&(r<=nro))
    {
      long nc=nco; bigint *mij = entries+(r-1)*nco;
      while(nc--) (*mij++) *= scal;
    }
  else
    {
      cerr << "Bad row number " << r << " in multrow"<<endl;
    }
}

void mat_m::divrow(long r, const bigint& scal)
{
  if ((0<r)&&(r<=nro))
    {
      long nc=nco; bigint *mij = entries+(r-1)*nco;
      while(nc--) (*mij++) /= scal;
    }
  else
    {
      cerr << "Bad row number " << r << " in divrow"<<endl;
    }
}

void mat_m::clearrow(long r)
{
  if ((0<r)&&(r<=nro))
    {
      bigint g; long nc=nco; bigint * mij = entries+(r-1)*nco;
      while((nc--)&&(!is_one(g))) g = gcd(g,(*mij++));
      if(is_zero(g)||is_one(g)) return;
      nc=nco; mij = entries+(r-1)*nco;
      while(nc--) (*mij++) /= g;
    }
  else
    {
      cerr << "Bad row number " << r << " in clearrow"<<endl;
    }
}

mat_m& mat_m::operator+=(const mat_m& mat2)
{
  if ((nro==mat2.nro) && (nco=mat2.nco)) 
    {
      long n=nro*nco; bigint *m1=entries, *m2=mat2.entries;
      while(n--) (*m1++) += (*m2++);
    }
  else
    {
      cerr << "Incompatible matrices in operator +="<<endl;
    }
  return *this;
}

mat_m& mat_m::operator-=(const mat_m& mat2)
{
  if ((nro==mat2.nro) && (nco=mat2.nco)) 
    {
      long n=nro*nco; bigint *m1=entries, *m2=mat2.entries;
      while(n--) (*m1++) -= (*m2++);
    }
  else
    {
      cerr << "Incompatible matrices in operator -="<<endl;
    }
  return *this;
}

mat_m& mat_m::operator*=(const bigint& scal)
{
  bigint* mij = entries; long n=nco*nro;
  while(n--) (*mij++) *= scal;
  return *this;
}

mat_m& mat_m::operator/=(const bigint& scal)
{
  bigint* mij = entries; long n=nco*nro;
  while(n--) (*mij++) /= scal;
  return *this;
}

mat_i mat_m::shorten(int x) const
{
  mat_i ans(nro,nco);
  bigint *matij=entries; int *ansij=ans.entries; long n=nro*nco;
  bigint minint; minint=MININT;
  bigint maxint; maxint=MAXINT;
  while(n--)
    {
      bigint mij = *matij++;
// NB gmp's test function here is broken!
      if((mij>=minint)&&(mij<=maxint)) 
	{
	  if(is_zero(mij)) *ansij=0;
	  else
	    {
	      int aij = I2int(mij);
	      *ansij =aij;
	      if(BIGINT(*ansij)!=mij)
		{
		  cerr<<"Problem: I2int("<<mij<<") returns "<<(*ansij)<<endl;
		}
	    }
	}
      else
	{
	  cerr << "Problem shortening bigint " << mij << " to an int!" << endl;
	}
      ansij++;
    }
  return ans;
}

mat_l mat_m::shorten(long x) const
{
  mat_l ans(nro,nco);
  bigint *matij=entries; long *ansij=ans.entries; long n=nro*nco;
  bigint minlong; minlong=MINLONG;
  bigint maxlong; maxlong=MAXLONG;
  while(n--)
    {
      bigint& mij = *matij++;
// NB gmp's test function here is broken!
      if((mij>=minlong)&&(mij<=maxlong)) 
	{
	  if(is_zero(mij)) *ansij=0;
	  else
	    {
	      long aij = I2long(mij);
	      *ansij = aij;
	      if(BIGINT(*ansij)!=mij)
		{
		  cerr<<"Problem: I2int("<<mij<<") returns "<<(*ansij)<<endl;
		}
	    }
	}
      else
	{
	  cerr << "Problem shortening bigint " << mij << " to a long!" << endl;
	}
      ansij++;
      }
  return ans;
}

// Definitions of non-member, friend operators and functions

mat_m operator*(const mat_m& m1, const mat_m& m2)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 mat_m m3(m,p);
 bigint *a=m1.entries, *b=m2.entries, *c=m3.entries;
 if (n==m2.nro)  // algorithm from Dr Dobb's Journal August 1993
   {
     while(m--)
       {
         bigint *bp=b;
         long k=n;
         while(k--)
           {
             bigint *cp=c;
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
     cerr << "Incompatible sizes in mat_m product"<<endl;
   }
 return m3;
}

int operator==(const mat_m& m1, const mat_m& m2)
{
   long nr=m1.nro, nc=m1.nco;
   int equal = ((nr==m2.nro) && (nc==m2.nco)); 
   if(!equal) return 0;
   bigint *m1ij=m1.entries, *m2ij=m2.entries; long n=nr*nc;
   while((n--)&&equal) equal=((*m1ij++)==(*m2ij++));
   return equal;
}

ostream& operator<<(ostream& s, const mat_m& m)
{
  bigint* mij=m.entries;
  s << "\n[";
  long nr=m.nro;
  while(nr--)
    {
      long nc=m.nco;
      s<<"[";
      while(nc--) {s<<(*mij++); if(nc) s<<",";}
      s<<"]"; if(nr) s<<",\n";
    }
  s << "]\n";
  return s;
}

istream& operator>>(istream& s, mat_m& m)
{
 long n=m.nro*m.nco;
 bigint* mij=m.entries;
 while(n--) s >> (*mij++);
 return s;
}

mat_m colcat(const mat_m& a, const mat_m& b)
{
 long nr = a.nro, nca = a.nco, ncb = b.nco;
 mat_m ans(nr,nca+ncb);
 bigint *ansij=ans.entries, *aij=a.entries, *bij=b.entries;
 if (nr==b.nro)
   while(nr--)
     {
       long nc=nca;
       while(nc--)
         *ansij++ = *aij++;
       nc=ncb;
       while(nc--)
         *ansij++ = *bij++;
     }
 else
   {
     cerr << "colcat: matrices have different number of rows!" << endl;
   }
 return ans;
}

mat_m rowcat(const mat_m& a, const mat_m& b)
{
 long nra = a.nro, nc = a.nco, nrb = b.nro;
 mat_m ans(nra+nrb,nc);
 bigint *ansij=ans.entries, *aij=a.entries, *bij=b.entries;
 if (nc==b.nco)
 {
   long n = nra*nc;
   while(n--)
     *ansij++ = *aij++;
   n = nrb*nc;
   while(n--)
     *ansij++ = *bij++;
 }
 else
   {
     cerr << "rowcat: matrices have different number of columns!" << endl;
   }
 return ans;
}

mat_m directsum(const mat_m& a, const mat_m& b)
{
  long n,c, nra=a.nro, nca=a.nco, nrb=b.nro, ncb=b.nco;
  mat_m ans(nra+nrb,nca+ncb);
  bigint* ansij=ans.entries, *aij=a.entries, *bij=b.entries;
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

void elimrows(mat_m& m, long r1, long r2, long pos)   
//plain elimination, no clearing
{
 long nc=m.nco;
 bigint *mr1 = m.entries + (r1-1)*nc,
         *mr2 = m.entries + (r2-1)*nc;
 bigint p = mr1[pos-1], q=mr2[pos-1];
 while(nc--)
   {
     (*mr2)= (p*(*mr2))-(q*(*mr1));
     mr1++; mr2++;
   }
}

void elimrows1(mat_m& m, long r1, long r2, long pos) 
//elimination + clearing
{
  elimrows(m,r1,r2,pos);
  m.clearrow(r2);
}

void elimrows2(mat_m& m, long r1, long r2, long pos, const bigint& last) 
//elimination + divide by last pivot
{
  elimrows(m,r1,r2,pos);
  m.divrow(r2,last);
}

// Definition of non-friend functions

mat_m operator+(const mat_m& m)
{return m;}

mat_m operator-(const mat_m& m)
{return BIGINT(-1)*m;}

mat_m operator+(const mat_m& m1, const mat_m& m2)
{mat_m ans(m1); ans+=m2; return ans;}

mat_m operator-(const mat_m& m1, const mat_m& m2) 
{mat_m ans(m1); ans-=m2; return ans;}

mat_m operator*(const bigint& scal, const mat_m& m)
{mat_m ans(m); ans*=scal; return ans;}

mat_m operator*(int scal, const mat_m& m)
{mat_m ans(m); ans*=BIGINT(scal); return ans;}

mat_m operator*(long scal, const mat_m& m)
{mat_m ans(m); ans*=BIGINT(scal); return ans;}

mat_m operator/(const mat_m& m, const bigint& scal)
{mat_m ans(m); ans/=scal; return ans;}

int operator!=(const mat_m& m1, const mat_m& m2)
{return !(m1==m2);}

vec_m operator*(const mat_m& m, const vec_m& v)
{
 long r=m.nro, c=m.nco;
 vec_m w(r);
 if (c==v.d)
   {
     bigint *mp=m.entries, *wp=w.entries;
     while(r--)
       {
         bigint *vp=v.entries;
         c=m.nco;
         while(c--)
           *wp += (*mp++)*(*vp++);
         wp++;
       }
   }
 else
   {
     cerr << "Incompatible sizes in *(mat_m,vec_m)"<<endl;
   }
 return w;
}

mat_m midmat(long n)
{
 mat_m ans(n,n);
 long i; bigint one; one=1;
 for (i=1; i<=n; i++) ans.set(i,i,one);
 return ans;
}

mat_m transpose(const mat_m& m)
{
 long i,j,nr,nc;
 nr=m.ncols(); nc=m.nrows();
 mat_m ans(nr, nc);
 for (i=1; i<=nr; i++)
  for (j=1; j<=nc; j++)
   ans.set(i,j,  m(j,i));
 return ans;
}

mat_m submatrix(const mat_m& m, const vec_i& iv, const vec_i& jv)
{long i,j;
 long nr = dim(iv);
 long nc = dim(jv);
 mat_m ans(nr,nc);
 for (i=1; i<=nr; i++)
  for (j=1; j<=nc; j++)
   ans.set(i,j, m(iv[i],jv[j]));
 return ans;
}

mat_m echelonp(const mat_m& m, vec_i& pcols, vec_i& npcols,
                                  long& rk, long& ny, bigint& d, const bigint& pr);

mat_m echelon(const mat_m& m, vec_l& pcols, vec_l& npcols,
               long& rk, long& ny, bigint& d, int method)
{
  vec_i pc, npc;
  mat_m ans = echelon(m,pc,npc,rk,ny,d,method);
  pcols.init(rk); npcols.init(ny);
  int i;
  for (i=1; i<=rk; i++)  pcols[i]= pc[i];
  for (i=1; i<=ny; i++) npcols[i]=npc[i];
  return ans;
}

mat_m echelon(const mat_m& m, vec_i& pcols, vec_i& npcols,
               long& rk, long& ny, bigint& d, int method)
{
//N.B. case 1 is for consistency with matrix.cc only: redundant.
  switch (method)
    {case 0: return echelon0(m,pcols,npcols,rk,ny,d);
     case 1: return echelon0(m,pcols,npcols,rk,ny,d);
     case 2: return echelonp(m,pcols,npcols,rk,ny,d,MBIGPRIME);
     default: return echelon0(m,pcols,npcols,rk,ny,d);
    }        
}

void elim(bigint *m, long nc, long r1, long r2, long pos)
{bigint *mr1=m+r1*nc, *mr2=m+r2*nc;
 bigint p = mr1[pos], q = mr2[pos];
 while(nc--)
   {
     (*mr2)=(p*(*mr2))-(q*(*mr1)); 
     mr1++; mr2++;
   }
}

void clear(bigint* row, long nc)
{long n=nc; bigint *rowi=row; bigint g;
 while((n--)&&(!is_one(g))) g=gcd(g,*rowi++);
 if (sign(g)<0) g=-g;
 if(is_zero(g)||is_one(g)) return;
 n=nc; rowi=row; while(n--) (*rowi++) /= g;
}

mat_m echelon0(const mat_m& m1, vec_i& pc, vec_i& npc,
                long& rk, long& ny, bigint& d)
{
  long nr, nc, r,c,r2,r3,i;
  bigint mr2c,lastpivot, *m, *mi2, temp;
  rk=0; ny=0; r=0; lastpivot=1;
  nc=m1.nco; nr=m1.nro;
  m = new bigint[nr*nc];
  long n=nr*nc;
  bigint *mij=m;
  bigint *mi1=m1.entries;
  while(n--)
    *mij++ = *mi1++;

  int *pcols = new int[nc];
  int *npcols = new int[nc];
  for (c=0; (c<nc)&&(r<nr); c++)
    {
      mij=m+r*nc+c;  // points to column c in row r
      bigint min = abs(*mij);
      long rmin = r;
      for (r2=r+1, mij+=nc; (r2<nr)&&(!is_one(min)); r2++, mij+=nc)
	{ mr2c = abs(*mij);
	  if ((sign(mr2c)>0) && ((mr2c<min) || (sign(min)==0))) 
	    { 
	      min=mr2c; 
	      rmin=r2 ;
	    }
	}
      if (sign(min)==0) npcols[ny++] = c;
      else
       {pcols[rk++] = c;
        if (rmin>r) //swap rows
          {
            mi1=m+r*nc; mi2=m+rmin*nc; n=nc;
            while(n--) {temp = *mi1; *mi1++ = *mi2; *mi2++ = temp;}
          } 
        for (r3 = r+1 ; r3<nr; r3++)
          {
            elim(m,nc,r,r3,c);
            mi1 = m+r3*nc; n=nc; while(n--) *mi1++ /= lastpivot;
          }
        lastpivot=min;
        r++;
      }
    }
  for (c = rk+ny; c<nc; c++) npcols[ny++] = c;
  d=1;
  if (ny>0)   // Back-substitute and even up pivots
    {for (r=0; r<rk; r++) clear(m+r*nc,nc);
     for (r=0; r<rk; r++)
       {
         for (r2=r+1; r2<rk; r2++) elim(m,nc,r2,r,pcols[r2]);  
         mi1=m+r*nc;
         clear(mi1,nc);
         d = lcm(d,mi1[pcols[r]]);
       }
     d = abs(d);
     // cout << "d = " << d << "\n";
     for (r=0, mij=m; r<rk; r++)
       {
         n=nc;
         bigint fac = d/mij[pcols[r]];  
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
  // Copy back into matrix
  mat_m ans(rk,nc);
  n=rk*nc; bigint* ansij=ans.entries; mij=m;
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

long mat_m::rank() const
{
 long rk,nr,nc,r,c,r2,r3;
 bigint mr2c,lastpivot;
 rk=0; r=1; lastpivot=1;
 mat_m m(*this);
 nc=m.ncols(); nr=m.nrows();
 for (c=1; (c<=nc)&&(r<=nr); c++)
 { bigint min = abs(m(r,c));
   long rmin = r;
   for (r2=r+1; (r2<=nr)&&(!is_one(min)); r2++)
   { mr2c=abs(m(r2,c));
     if ((sign(mr2c)>0) && ((mr2c<min) || (sign(min)==0)))
       {
	 min=mr2c;
	 rmin=r2 ;
       }
   }
   if (sign(min)!=0)
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

long mat_m::nullity() const
{
 return ncols()-rank();
}

bigint mat_m::trace() const
{ long i=0; bigint* aii=entries; bigint ans=BIGINT(0);
  for (; i<nro; i++, aii+=(nco+1))
    ans += *aii;
  return ans;
}

// FADEEV'S METHOD

vector<bigint> mat_m::charpoly() const
{ long n = nrows();
  mat_m b(*this);
  mat_m id(midmat(n)), tid;
  vector<bigint> clist(n+1);
  bigint t = trace(), ii;
  clist[n]   =  1;
  clist[n-1] = -t;
  for (long i=2; i<=n; i++)
      { tid=t*id;
	b-=tid;
	b=b*(*this);          //     cout << b;   // (for testing only)
	ii=i;
        t=b.trace()/ii;
        clist[n-i] = -t;
      }
  tid=t*id;
  if (b!=tid)
    {
      cerr << "Error in charpoly: final b = " << (b-t*id) << endl;
    }
  return clist;
}

bigint mat_m::determinant() const
{
 bigint det = charpoly()[0];
 if (nrows()%2==1)
   return -det;
 else
   return det;
}

mat_m addscalar(const mat_m& m, const bigint& c)
{
  mat_m ans(midmat(m.nrows()));
  ans*=c;
  ans+=m;
  return ans;
}
 
vec_m apply(const mat_m& m, const vec_m& v)    // same as *(mat_m, vec_m)
{
 long nr=m.nrows(), nc=m.ncols();
 vec_m ans(nr);
 if (nc==dim(v))
   for (long i=1; i<=nr; i++) 
     ans[i] = m.row(i)*v;
 else
   {
     cerr << "Incompatible sizes in *(mat_m,vec_m)"<<endl;
   }
 return ans;
}

void elimp(const mat_m& m, long r1, long r2, long pos, const bigint& pr)
{
 long nc=m.nco;
 bigint *mr1 = m.entries + (r1-1)*nc, *mr2 = m.entries + (r2-1)*nc;
 bigint p = mr1[pos-1], q=mr2[pos-1];
 while(nc--)
   {
     (*mr2)= mod(mod(p*(*mr2),pr)-mod(q*(*mr1),pr),pr);
     mr1++; mr2++;
   }
}

//#define TRACE 1

mat_m echelonp(const mat_m& m1, vec_i& pcols, vec_i& npcols,
                 long& rk, long& ny, bigint& d, const bigint& pr)
{
#ifdef TRACE
  cout << "In echelonp\n";
#endif /* TRACE */
 long nc,nr,r,c,r2,r3,rmin;
 bigint min, mr2c,lastpivot;
 nr=m1.nrows(), nc=m1.ncols();
 mat_m m(nr,nc);
 for (c=1; c<=nc; c++)
   for (r=1; r<=nr; r++)
     m(r,c)=mod(m1(r,c),pr);
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0; r=1; lastpivot=1;
 for (c=1; (c<=nc)&&(r<=nr); c++)
 {
   min = m(r,c);   rmin = r;
   for (r2=r+1; (r2<=nr)&&(sign(min)==0); r2++)
   { mr2c = m(r2,c);
     if (0!=sign(mr2c)) { min=mr2c; rmin=r2 ;}
   }
   if (sign(min)==0) npcols[++ny] = c;
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
      { const bigint& temp = m(r,pcols[r]);
	const bigint& fac = invmod(temp,pr);
        for (c=1; c<=nc; c++) 
	  {
	    const bigint& tmp = fac*m(r,c); 
	    m(r,c)=mod(tmp,pr);
	  }
      }
 }
 else 
   for (r=1; r<=rk; r++)
     {
       for (c=1; c<=nc; c++) m(r,c)=(c==pcols[r]);    // 0 or 1 !
     }
  bigint modulus=pr;
  bigint lim=sqrt(pr>>1);
#ifdef TRACE
  cout << "Finished second stage.\n Echelon matrix mod "<<pr<<" is:\n";
  cout << m;
  cout << "Now lifting back to Q.\n";
  cout << "lim = " << lim << "\n";
#endif /* TRACE */
  bigint dd; dd=1;
  mat_m nmat(rk,nc);
  mat_m dmat(rk,nc);

#ifdef TRACE
  cout << "rk = " << rk << "\n";
  cout << "ny = " << ny << "\n";
#endif /* TRACE */
  long i,j;
  for (i=1; i<=rk; i++) 
    {  
      for (j=1; j<=rk; j++) 
        {
          nmat(i,pcols[j])=(i==j);
          dmat(i,pcols[j])=1;
        }
      for (j=1; j<=ny; j++)
        {bigint num, den;
         long jj = npcols[j];
         modrat(m(i,jj),modulus,lim,num,den);
         nmat(i,jj)=num;
         dmat(i,jj)=den;
         dd=(dd*den)/gcd(dd,den);
       }
    }
  dd=abs(dd);
#ifdef TRACE
  cout << "Numerator matrix = " << nmat;
  cout << "Denominator matrix = " << dmat;
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

mat_m echmodp(const mat_m& m1, vec_i& pcols, vec_i& npcols,
                long& rk, long& ny, const bigint& pr)
{
// cout << "In echmodp with mat = " << m1;
 long nc,nr,r,c,r2,r3;
 nr=m1.nro, nc=m1.nco;
 mat_m m(nr,nc);
 bigint *mij=m.entries, *matij=m1.entries;
 long n=nr*nc;
 while(n--) *mij++ = mod(*matij++,pr);
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0; r=1;
 for (c=1; (c<=nc)&&(r<=nr); c++)
   {
     mij=m.entries+(r-1)*nc+c-1;
     bigint min = *mij;
     long rmin = r;
     for (r2=r+1, mij+=nc; (r2<=nr)&&(sign(min)==0); r2++, mij+=nc)
       {
         bigint mr2c = *mij;
         if (!is_zero(mr2c)) { min=mr2c; rmin=r2 ;}
       }
     if (sign(min)==0) npcols[++ny] = c;
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
         bigint fac = invmod(mij[pcols[r]-1],pr);
         n=nc; while(n--) {*mij = mod(fac * *mij, pr); mij++;}
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

mat_m matmulmodp(const mat_m& m1, const mat_m& m2, const bigint& pr)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 mat_m m3(m,p);
 bigint *a=m1.entries, *b=m2.entries, *c=m3.entries;
 if (n==m2.nro)  // algorithm from Dr Dobb's Journal August 1993
   {
     while(m--)
       {
         bigint *bp=b;
         long k=n;
         while(k--)
           {
             bigint *cp=c;
             long j=p;
             while(j--)
               {
                 *cp += mod(*a * *bp++, pr);
                 *cp = mod(*cp,pr);
                 cp++;
               }
             a++;
           }
         c += p;
       }
   }
 else
   {
     cerr << "Incompatible sizes in mat_m product"<<endl;
   }
 return m3;
}



