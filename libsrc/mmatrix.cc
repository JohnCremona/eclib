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

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef svec
#undef smat
#undef smat_elim

#define scalar bigint
#define vec vec_m
#define mat mat_m
#define subspace subspace_m
#define svec svec_m
#define smat smat_m
#define smat_elim smat_m_elim
#undef svec
#undef smat
#undef smat_elim

#include "mat.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

mat_m to_mat_m(const mat_i& m)
{
  const vector<int> & mij = m.get_entries();
  vector<bigint> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), [](const int& x) {return bigint(x);});
  return mat_m(m.nrows(), m.ncols(), n);
}

mat_m to_mat_m(const mat_l& m)
{
  const vector<long> & mij = m.get_entries();
  vector<bigint> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), [](const long& x) {return bigint(x);});
  return mat_m(m.nrows(), m.ncols(), n);
}

mat_i to_mat_i(const mat_m& m)
{
  const vector<bigint> & mij = m.get_entries();
  auto toint = [](const bigint& a) {return is_int(a)? I2int(a) : int(0);};
  vector<int> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), toint);
  return mat_i(m.nrows(), m.ncols(), n);
}

mat_l to_mat_l(const mat_m& m)
{
  const vector<bigint> & mij = m.get_entries();
  auto tolong = [](const bigint& a) {return is_long(a)? I2long(a) : long(0);};
  vector<long> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), tolong);
  return mat_l(m.nrows(), m.ncols(), n);
}

#if(0)

// Definitions of member operators and functions:

mat_m::mat_m(const mat_i& m)
  :nro(m.nro), nco(m.nco)
{
  entries.resize(nro*nco);
  std::transform(m.entries.begin(), m.entries.end(), entries.begin(), [](const int& x) {return bigint(x);});
}

mat_m::mat_m(const mat_l& m)
  :nro(m.nro), nco(m.nco)
{
  entries.resize(nro*nco);
  std::transform(m.entries.begin(), m.entries.end(), entries.begin(), [](const long& x) {return bigint(x);});
}

void mat_m::init(long nr, long nc) // assigns to zero matrix of given size;
{                                 // with defaults (0,0) releases all space.
  nro = nr;
  nco = nc;
  entries.resize(nro*nco);
}

bigint& mat_m::operator()(long i, long j)   // returns ref to (i,j) entry
{
  return entries.at((i-1)*nco+(j-1));
}

bigint mat_m::operator()(long i, long j)  const   // returns (i,j) entry
{
  return entries.at((i-1)*nco+(j-1));
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
 long n=r2-r1+1,c=c2-c1+1;
 mat_m ans(n,c);
 auto ap=ans.entries.begin();
 auto mp=entries.begin()+r1*nco+c1;
 while(n--)
   {
     std::copy(mp, mp+c, ap);
     ap += c;
     mp += nco;
   }
 return ans;
}

mat_m& mat_m::operator=(const mat_m& m)
{
 if (this==&m) return *this;
 nro=m.nro;
 nco=m.nco;
 entries = m.entries;
 return *this;
}

bigint mat_m::sub(long i, long j) const
{
  return entries.at((i-1)*nco+(j-1));
}

void mat_m::set(long i, long j, const bigint& x)
{
  entries.at((i-1)*nco+(j-1)) = x;
}

void mat_m::add(long i, long j, const bigint& x)
{
  if (!is_zero(x)) entries.at((i-1)*nco+(j-1)) += x;
}

void mat_m::setrow(long i, const vec_m& v)
{
  std::copy(v.entries.begin(), v.entries.end(), entries.begin() + (i-1)*nco);
}

void mat_m::setcol(long j, const vec_m& v)
{
  auto colj = entries.begin()+(j-1);
  for ( const auto vi : v.entries)
    {
      *colj = vi;
      colj += nco;
    }
}

vec_m mat_m::row(long i) const
{
 vec_m mi(nco);
 auto e = entries.begin()+(i-1)*nco;
 std::copy(e, e+nco, mi.entries.begin());
 return mi;
}

vec_m mat_m::col(long j) const
{
 vec_m v(nro);
 auto entriesij = entries.begin()+(j-1);
 for ( auto& vi : v.entries)
   {
     vi = *entriesij;
     entriesij+=nco;
   }
 return v;
}

void mat_m::swaprows(long r1, long r2)
{
  auto mr1 = entries.begin() + (r1-1)*nco;
  auto mr2 = entries.begin() + (r2-1)*nco;
  std::swap_ranges(mr1, mr1+nco, mr2);
}

void mat_m::multrow(long r, const bigint& scal)
{
  if (scal==1) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [scal](const bigint& x) {return x * scal;});
}

void mat_m::divrow(long r, const bigint& scal)
{
  if (scal<=1) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [scal](const bigint& x) {return x / scal;});
}

bigint mat_m::content() const
{
  return std::accumulate(entries.begin(), entries.end(), bigint(0),
                         [](const bigint& x, const bigint& y) {return gcd(x,y);});
}

bigint mat_m::row_content(long r) const
{
  auto mij = entries.begin()+(r-1)*nco;
  return std::accumulate(mij, mij+nco, bigint(0),
                         [](const bigint& x, const bigint& y) {return gcd(x,y);});
}

void mat_m::clearrow(long r)
{
  divrow(r, row_content(r));
}

void mat_m::makeprimitive()
{
  bigint g = content();
  if (g<=1) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [g](const bigint& x) {return x / g;});
}

void mat_m::operator+=(const mat_m& n)
{
  std::transform(n.entries.begin(), n.entries.end(), entries.begin(), entries.begin(),
                 [](const bigint& x, const bigint& y) { return x + y;});
}

void mat_m::operator-=(const mat_m& n)
{
  std::transform(n.entries.begin(), n.entries.end(), entries.begin(), entries.begin(),
                 [](const bigint& x, const bigint& y) { return y - x;});
}

void mat_m::operator*=(const bigint& scal)
{
  if (is_one(scal))
    return;
  if (is_zero(scal))
    std::fill(entries.begin(), entries.end(), bigint(0));
  else
    std::transform(entries.begin(), entries.end(), entries.begin(),
                   [scal](const bigint& x) {return x * scal;});
}

void mat_m::operator/=(const bigint& scal)
{
  if (is_one(scal) || is_zero(scal))
    return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const bigint& x) {return x / scal;});
}

mat_i mat_m::shorten(int x) const
{
  mat_i ans(nro,nco);
  std::transform(entries.begin(), entries.end(), ans.entries.begin(), [](const bigint& x) {return I2int(x);});
  return ans;
}

mat_l mat_m::shorten(long x) const
{
  mat_l ans(nro,nco);
  std::transform(entries.begin(), entries.end(), ans.entries.begin(), [](const bigint& x) {return I2long(x);});
  return ans;
}

// Definitions of non-member, friend operators and functions

mat_m operator*(const mat_m& m1, const mat_m& m2)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 mat_m m3(m,p);
 if (n==m2.nro)  // algorithm from Dr Dobb's Journal August 1993
   {
     auto a=m1.entries.begin();                                     // a points to m1(i,k)
     for (auto c=m3.entries.begin(); c!=m3.entries.end(); c+=p)     // c points to m3(i,_) for 0<=i<m
       {
         for (auto b=m2.entries.begin(); b!=m2.entries.end(); b+=p) // b points to m2(k,_) for 0<=k<n
           { // add m1(i,k)*m2(k,j) to m3(i,j) for 0<=j<p
             bigint m1ik = *a++;
             std::transform(b, b+p, c, c,
                            [m1ik] (const bigint& m2kj, const bigint& m3ij) {return m1ik*m2kj+m3ij;});
           }
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
  return (m1.nro==m2.nro) && (m1.nco==m2.nco) && (m1.entries==m2.entries);
}

ostream& operator<<(ostream& s, const mat_m& m)
{
  auto mij=m.entries.begin();
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

istream& operator>>(istream& s, mat_m& m) // m cannot be const
{
 long n=m.nro*m.nco;
 auto mij=m.entries.begin();
 while(n--) s >> (*mij++);
 return s;
}

mat_m colcat(const mat_m& a, const mat_m& b)
{
 long nr = a.nro, nca = a.nco, ncb = b.nco;
 mat_m c(nr,nca+ncb);
 if (nr==b.nro)
   {
     auto aij = a.entries.begin();
     auto bij = b.entries.begin();
     auto cij = c.entries.begin();
     while (cij!=c.entries.end())
       {
         std::copy(aij, aij+nca, cij);
         aij+=nca;
         cij+=nca;
         std::copy(bij, bij+ncb, cij);
         bij+=ncb;
         cij+=ncb;
       }
   }
 else
   {
     cerr << "colcat: matrices have different number of rows!" << endl;
   }
 return c;
}

mat_m rowcat(const mat_m& a, const mat_m& b)
{
 mat_m c(a.nro+b.nro,a.nco);
 if (a.nco==b.nco)
 {
   auto cij = c.entries.begin();
   std::copy(a.entries.begin(), a.entries.end(), cij);
   cij += a.entries.size();
   std::copy(b.entries.begin(), b.entries.end(), cij);
 }
 else
   cerr << "rowcat: matrices have different number of columns!" << endl;
 return c;
}

mat_m directsum(const mat_m& a, const mat_m& b)
{
  return rowcat(colcat(a,mat_m(a.nro,b.nco)),colcat(mat_m(b.nro,a.nco),b));
}

//plain elimination, no clearing
void elimrows(mat_m& m, long r1, long r2, long pos)  // m cannot be const
{
  long nc=m.nco;
  bigint p = m(r1,pos), q=m(r2,pos);
  auto mr1 = m.entries.begin() + (r1-1)*nc;
  auto mr2 = m.entries.begin() + (r2-1)*nc;
  // replace row2 by p*row2-q*row1
  std::transform(mr1, mr1+nc, mr2, mr2,
                 [p,q] (const bigint& x, const bigint& y) {return p*y-q*x;});
}

//elimination + clearing
void elimrows1(mat_m& m, long r1, long r2, long pos)
{
  elimrows(m,r1,r2,pos);
  m.clearrow(r2);
}

//elimination + divide by last pivot
void elimrows2(mat_m& m, long r1, long r2, long pos, const bigint& last)
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
{return BIGINT(scal)*m;}

mat_m operator*(long scal, const mat_m& m)
{return BIGINT(scal)*m;}

mat_m operator/(const mat_m& m, const bigint& scal)
{mat_m ans(m); ans/=scal; return ans;}

int operator!=(const mat_m& m1, const mat_m& m2)
{return !(m1==m2);}

vec_m operator*(const mat_m& m, const vec_m& v)
{
 long c=m.nco;
 vec_m w(m.nro);
 if (c==dim(v))
   {
     auto mi = m.entries.begin();
     for (auto& wi : w.entries)
       {
         wi = std::inner_product(mi, mi+c, v.entries.begin(), bigint(0));
         mi += c;
       }
   }
 else
   cerr << "Incompatible sizes in *(mat_m,vec_m)"<<endl;
 return w;
}

mat_m midmat(long n)
{
 mat_m ans(n,n);
 bigint one(1);
 for (long i=1; i<=n; i++) ans.set(i,i,one);
 return ans;
}

mat_m transpose(const mat_m& m)
{
 long nr=m.ncols(), nc=m.nrows();
 mat_m ans(nr, nc);
 for (long i=1; i<=nr; i++)
  for (long j=1; j<=nc; j++)
   ans.set(i,j,  m(j,i));
 return ans;
}

// submatrix of rows indexed by v, all columns
mat_m rowsubmat(const mat_m& m, const vec_i& v)
{
  long nr = dim(v), nc = m.ncols();
  mat_m ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(v[i],j));
  return ans;
}

// submatrix of rows indexed by iv, columns indexed by jv
mat_m submatrix(const mat_m& m, const vec_i& iv, const vec_i& jv)
{
  long nr = dim(iv), nc = dim(jv);
  mat_m ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
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

void elim(vector<bigint>& m, long nc, long r1, long r2, long pos)
{
  auto mr1=m.begin()+r1*nc;
  auto mr2=m.begin()+r2*nc;
  bigint p = *(mr1+pos), q = *(mr2+pos);
  if (is_one(p)&&is_zero(q))
    return;
  std::transform(mr1, mr1+nc, mr2, mr2,
                 [p,q](const bigint& x, const bigint& y) {return p*y - q*x;});
}

void clear(vector<bigint>& row, long col1, long col2)
{
  auto row1=row.begin()+col1;
  auto row2=row.begin()+col2;
  bigint g = std::accumulate(row1, row2, bigint(0),
                             [](const bigint& x, const bigint& y) {return gcd(x,y);} );
  if (g>1)
    std::for_each(row1, row2, [g](bigint& x) {x/=g;});
}

//#define DEBUG_ECH_0

#ifdef DEBUG_ECH_0
void show(vector<bigint> m, long nr, long nc)
{
  auto mij = m.begin();
  for(long i=0; i<nr; i++)
    {
      for(long j=0; j<nc; j++)
	cout<<(*mij++)<<"\t";
      cout<<"\n";
    }
}
#endif

mat_m echelon0(const mat_m& m1, vec_i& pc, vec_i& npc,
               long& rk, long& ny, bigint& d)
{
#ifdef DEBUG_ECH_0
  cout<<"In echelon0 with matrix:\n"<<m1<<endl;
#endif
  rk=0; ny=0;
  bigint lastpivot(1);
  long r=0, nc=m1.nco, nr=m1.nro;
  vector<bigint> m = m1.entries;
  vector<int> pcols(nc), npcols(nc);
  for (long c=0; (c<nc)&&(r<nr); c++)
    {
      auto mij=m.begin()+r*nc+c;  // points to column c in row r
      bigint piv = abs(*mij);
      long rmin = r;
      mij+=nc;
      for (long r2=r+1; (r2<nr)&&(piv!=1); r2++, mij+=nc)
       {
         bigint mr2c = abs(*mij);
         if ((0<mr2c) && ((mr2c<piv) || (piv==0)))
           {
             piv=mr2c;
             rmin=r2;
           }
       }
      if (is_zero(piv))
        npcols[ny++] = c;
      else
       {
         pcols[rk++] = c;
#ifdef DEBUG_ECH_0
         cout<<"Using col "<<c<<" as pivotal col; pivot="<<piv<<" in row "<<rmin<<endl;
#endif
         if (rmin>r) //swap rows
	  {
#ifdef DEBUG_ECH_0
	    cout<<"Swapping rows "<<r<<" and "<<rmin<<endl;
#endif
            auto mr1 = m.begin() + r*nc;
            auto mr2 = m.begin() + rmin*nc;
            std::swap_ranges(mr1, mr1+nc, mr2);
	  }
         for (long r3 = r+1 ; r3<nr; r3++)
          {
#ifdef DEBUG_ECH_0
	    cout<<"Eliminating from row "<<r3<<endl;
            cout<<"Before, m is\n"; show(m,nr,nc);
#endif
            elim(m,nc,r,r3,c);
#ifdef DEBUG_ECH_0
            cout<<"After, m is\n"; show(m,nr,nc);
#endif
	    if(lastpivot>1)
	      {
		auto mi1 = m.begin()+r3*nc;
                std::transform(mi1, mi1+nc, mi1, [lastpivot]( const bigint& x) {return x/lastpivot;});
              }
          }
         lastpivot=piv;
#ifdef DEBUG_ECH_0
         cout<<"r="<<r<<": pivot = "<<piv<<endl;
#endif
         r++;
       }
#ifdef DEBUG_ECH_0
      cout<<"Current mat is:\n";show(m,nr,nc);
#endif
    }
  for (long c = rk+ny; c<nc; c++) npcols[ny++] = c;
  d=1;
  if (ny>0)   // Back-substitute and even up pivots
    {
      for (long r1=0; r1<rk; r1++)
        clear(m, r1*nc, (r1+1)*nc); // divides row by its content
#ifdef DEBUG_ECH_0
      cout<<"After clearing, pivots are:"<<endl;
      for(long r3=0; r3<rk; r3++)
        cout<<*(m.begin()+r3*nc+pcols[r3])<<",";
      cout<<endl;
#endif
      for (long r1=0; r1<rk; r1++)
        {
          auto mi1 = m.begin()+r1*nc;
#ifdef DEBUG_ECH_0
          cout<<"Before back-subst, row "<<r<<" is:"<<endl;
          for(long r3=0; r3<nc; r3++)
            cout<<*(mi1+r3)<<",";
          cout<<": pivot = "<<*(mi1+pcols[r1])<<endl;
#endif
          for (long r2=r1+1; r2<rk; r2++)
            elim(m,nc,r2,r1,pcols[r2]);
#ifdef DEBUG_ECH_0
          cout<<"After back-subst, row "<<r<<" is:"<<endl;
          for(long r3=0; r3<nc; r3++)
            cout<<*(mi1+r3)<<",";
          cout<<": pivot = "<<*(mi1+pcols[r1])<<endl;
#endif
          clear(m, r1*nc, (r1+1)*nc);
#ifdef DEBUG_ECH_0
          cout<<"After clearing, row "<<r1<<" is:"<<endl;
          for(long r3=0; r3<nc; r3++)
            cout<<*(mi1+r3)<<",";
          cout<<": pivot = "<<*(mi1+pcols[r1])<<endl;
#endif
          d = lcm(d, *(mi1+pcols[r1]));
        }
      d = abs(d);
      // cout << "d = " << d << "\n";
      auto mij = m.begin();
      for (long r1=0; r1<rk; r1++)
        {
          bigint fac = d/mij[pcols[r1]];
          std::transform(mij, mij+nc, mij, [fac](const bigint& x){return fac*x;});
          mij += nc;
        }
    }
  else
    {
      auto mij = m.begin();
      for (long i=0; i<rk; i++)
	for (long j=0; j<nc; j++)
	  *mij++ = (j==pcols[i]);  // 0 or 1 !
    }

  // fix vectors
  pc.init(rk); npc.init(ny);
  for (long i=0; i<rk; i++)  pc[i+1]= pcols[i]+1;
  for (long i=0; i<ny; i++) npc[i+1]=npcols[i]+1;

  // Copy back into mat
  mat_m ans(rk,nc, m);
  return ans;
}

long mat_m::rank() const
{
  long rk=0;
  bigint lastpivot(1);
  mat_m m(*this); // work with a copy, which will be reduced
  for (long c=1, r=1; (c<=nco)&&(r<=nro); c++)
    {
      bigint piv = abs(m(r,c));
      long rmin = r;
      for (long r2=r+1; (r2<=nro)&&(!is_one(piv)); r2++)
        {
          bigint mr2c=abs(m(r2,c));
          if ((!is_zero(mr2c)) && ((mr2c<piv) || is_zero(piv)))
            {
              piv=mr2c;
              rmin=r2 ;
            }
        }
      if (!is_zero(piv))
        {
          rk++;
          if (rmin>r) m.swaprows(r,rmin);
          for (long r3 = r+1 ; r3<=nro; r3++)
            elimrows2(m,r,r3,c,lastpivot);
          lastpivot=piv;
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
{
  bigint tr(0);
  for (long i=0; i<nro; i++)
    tr += entries.at(i*(nco+1));
  return tr;
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
 return (nro%2? -det :det);
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
  return m*v;
}

void elimp(mat_m& m, long r1, long r2, long pos, const bigint& pr)
{
 long nc=m.nco;
 auto mr1 = m.entries.begin() + (r1-1)*nc + (pos-1);
 auto mr2 = m.entries.begin() + (r2-1)*nc + (pos-1);
 bigint p = mod(*mr1,pr), q=mod(*mr2,pr);
 if(q==0) {return;} // nothing to do
 nc -= (pos-1); // first pos-1 entries are assumed 0 already
 std::transform(mr1, mr1+nc, mr2, mr2,
                [pr,p,q](const bigint& x, const bigint& y) {return mod(mod(p*y,pr)-mod(q*x,pr), pr);});
}

//#define TRACE 1

mat_m echelonp(const mat_m& m1, vec_i& pcols, vec_i& npcols,
                 long& rk, long& ny, bigint& d, const bigint& pr)
{
#ifdef TRACE
  cout << "In echelonp\n";
#endif /* TRACE */
 long nr=m1.nrows(), nc=m1.ncols();
 mat_m m(nr,nc);
 std::transform(m1.entries.begin(), m1.entries.end(), m.entries.begin(),
                [pr] (const bigint& x) {return mod(x,pr);});
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0;
 long r=1;
 for (long c=1; (c<=nc)&&(r<=nr); c++)
 {
   bigint piv = m(r,c);
   long rmin = r;
   for (long r2=r+1; (r2<=nr)&&(piv==0); r2++)
   {
     bigint mr2c = m(r2,c);
     if (!is_zero(mr2c))
       {
         piv=mr2c;
         rmin=r2;
       }
   }
   if (is_zero(piv))
     npcols[++ny] = c;
   else
     {
       pcols[++rk] = c;
       if (rmin>r) m.swaprows(r,rmin);
       for (long r3 = r+1 ; r3<=nr; r3++)
         elimp(m,r,r3,c,pr);
       r++;
     }
 }
 for (long c = rk+ny+1; c<=nc; c++)
   npcols[++ny] = c ;
#ifdef TRACE
 cout << "Finished first stage; rk = " << rk;
 cout << ", ny = " << ny << "\n";
 cout << "Back substitution.\n";
#endif /* TRACE */
 pcols  =  pcols.slice(1,rk);
 npcols =  npcols.slice(1,ny);    // truncate index vectors
 if (ny>0)
 {
   for (long r1=1; r1<=rk; r1++)
     for (long r2=r+1; r2<=rk; r2++)
       elimp(m,r2,r1,pcols[r2],pr);
   for (long r1=1; r1<=rk; r1++)
     {
       bigint fac = mod(invmod(m(r1,pcols[r1]),pr),pr);
       for (long c=1; c<=nc; c++)
         m(r1,c)=mod(fac*m(r1,c),pr);
     }
 }
 else
   for (long i=1; i<=rk; i++)
     for (long j=1; j<=nc; j++)
       m(i,j)=(j==pcols[i]);    // 0 or 1 !

 bigint lim=sqrt(pr>>1);
#ifdef TRACE
 cout << "Finished second stage.\n Echelon mat mod "<<pr<<" is:\n";
 cout << m;
 cout << "Now lifting back to Q.\n";
 cout << "lim = " << lim << "\n";
#endif /* TRACE */
 bigint dd(1);
 mat_m nmat(rk,nc);
 mat_m dmat(rk,nc);

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
         bigint n1,d1;
         long jj = npcols[j];
         modrat(m(i,jj),pr,lim,n1,d1);
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

mat_m echmodp(const mat_m& m1, vec_i& pcols, vec_i& npcols,
                long& rk, long& ny, const bigint& pr)
{
  // cout << "In echmodp with p="<<pr<<" and matrix " << m1 << endl;
 long nr=m1.nrows(), nc=m1.ncols();
 mat_m m(nr,nc);
 std::transform(m1.entries.begin(), m1.entries.end(), m.entries.begin(),
                [pr] (const bigint& x) {return mod(x,pr);});
 // cout << " - after reducing modulo p,  matrix is " << m << endl;
 pcols.init(nc);
 npcols.init(nc);
 rk=ny=0;
 long r=1;
 for (long c=1; (c<=nc)&&(r<=nr); c++)
   {
     auto mij=m.entries.begin()+(r-1)*nc+c-1;
     bigint mmin = *mij;
     long rmin = r;
     mij += nc;
     for (long r2=r+1; (r2<=nr)&&(is_zero(mmin)); r2++, mij+=nc)
       {
	 bigint mr2c = *mij;
	 if (!is_zero(mr2c))
           {
             mmin=mr2c;
             rmin=r2;
           }
       }
     if (is_zero(mmin))
       npcols[++ny] = c;
     else
       {
	 pcols[++rk] = c;
	 if (rmin>r)
           m.swaprows(r,rmin);
	 auto entriesij = m.entries.begin()+(r-1)*nc;
         // cout<<"pivot = "<<mmin<<endl;
         bigint fac = mod(invmod(mmin,pr),pr);
         std::transform(entriesij, entriesij+nc, entriesij,
                        [pr,fac] (const bigint& x) {return mod(fac*x, pr);});
         for (long r3 = r+1 ; r3<=nr; r3++)
           elimp(m,r,r3,c,pr);
	 r++;
       }
     // cout << "After c="<<c<<" elimination, matrix is "<<m<<endl;
   }
 for (long c = rk+ny+1; c<=nc; c++)
   npcols[++ny] = c ;
 pcols  =  pcols.slice(rk);
 npcols =  npcols.slice(ny);    // truncate index vectors
 // cout << "After forward elimination, matrix is "<<m<<endl;
 // cout << "Rank = " << rk << ".  Nullity = " << ny << ".\n";
 if (ny>0)
   {
     for (long r1=1; r1<=rk; r1++)
       for (long r2=r1+1; r2<=rk; r2++)
	 elimp(m,r2,r1,pcols[r2],pr);
     for (long r1=1; r1<=rk; r1++)
       {
	 auto mij = m.entries.begin()+(r1-1)*nc;
	 bigint fac = *(mij+pcols[r1]-1);
	 fac = mod(invmod(fac,pr),pr);
         std::transform(mij, mij+nc, mij,
                        [pr,fac] (const bigint& x) {return mod(fac*x, pr);});
       }
   }
 else
   {
     auto mij=m.entries.begin();
     for (long i=1; i<=rk; i++)
       for (long j=1; j<=nc; j++)
	 *mij++ = (j==pcols[i]);    // 0 or 1 !
   }
 return m.slice(rk,nc);
}

mat_m matmulmodp(const mat_m& m1, const mat_m& m2, const bigint& pr)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 mat_m m3(m,p);
 if (n==m2.nro)
   {
     auto a=m1.entries.begin();                                     // a points to m1(i,k)
     for (auto c=m3.entries.begin(); c!=m3.entries.end(); c+=p)     // c points to m3(i,_) for 0<=i<m
       {
         for (auto b=m2.entries.begin(); b!=m2.entries.end(); b+=p) // b points to m2(k,_) for 0<=k<n
           { // add m1(i,k)*m2(k,j) to m3(i,j) for 0<=j<p
             bigint m1ik = *a++;
             std::transform(b, b+p, c, c,
                            [pr,m1ik] (const bigint& m2kj, const bigint& m3ij)
                            {return mod(mod(m1ik*m2kj,pr)+m3ij, pr);});
           }
       }
   }
 else
   {
     cerr << "Incompatible sizes in mat product"<<endl;
   }
 return m3;
}

int liftmat(const mat_m& mm, const bigint& pr, mat_m& m, bigint& dd, int trace)
{
  bigint lim = sqrt(pr>>1);
  if(trace)
    cout << "Lifting mod-p mat;  mat mod "<<pr<<" is:\n"
         << mm
         << "Now lifting back to Q.\n"
         << "lim = " << lim << endl;

  int success=1;
  bigint n,d;
  m = mm;
  dd=1;
  std::for_each(m.entries.begin(), m.entries.end(),
                [&success,&dd,pr,lim,&n,&d] (const bigint& x)
                {success=success&&modrat(x,pr,lim,n,d); d=lcm(d,dd);});
  if(!success)
    return 0;

  dd=abs(dd);
  if(trace)
    cout << "Common denominator = " << dd << "\n";
  std::transform(m.entries.begin(), m.entries.end(), m.entries.begin(),
                 [pr,dd] (const bigint& x) {return mod(mod(dd*x,pr),pr);});
  return 1;
}


#endif
