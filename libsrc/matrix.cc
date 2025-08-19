// matrix.cc: manage implementation of integer matrix classes
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
 
#include <eclib/matrix.h>

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar int
#define vec vec_i
#define mat mat_i
#define subspace subspace_i
#define scalar_is_int

#include "mat.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace
#undef scalar_is_int

#define scalar long
#define vec vec_l
#define mat mat_l
#define subspace subspace_l

#include "mat.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

#define scalar bigint
#define vec vec_m
#define mat mat_m
#define subspace subspace_m

#include "mat.cc"

#undef scalar
#undef vec
#undef mat
#undef subspace

// The following functions are here and not in mat.cc since they are
// not to be created in 3 versions

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

mat_i to_mat_i(const mat_l& m)
{
  const vector<long> & mij = m.get_entries();
  auto toint = [](const long& a) {return int(a);};
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

mat_l to_mat_l(const mat_i& m)
{
  const vector<int> & mij = m.get_entries();
  auto tolong = [](const int& a) {return long(a);};
  vector<long> n(mij.size());
  std::transform(mij.begin(), mij.end(), n.begin(), tolong);
  return mat_l(m.nrows(), m.ncols(), n);
}

///////////////////////////////////////////////////////////////////////////

// Definitions of member operators and functions:

template<class T>
void matT<T>::init(long nr, long nc) // resets to zero mat of given size;
{                                // with defaults (0,0) releases all space.
  nro = nr;
  nco = nc;
  entries.resize(nro*nco, T(0));
}

template<class T>
T& matT<T>::operator()(long i, long j)   // returns ref to (i,j) entry
{
  return entries.at((i-1)*nco+(j-1));
}

template<class T>
T matT<T>::operator()(long i, long j) const   // returns (i,j) entry
{
  return entries.at((i-1)*nco+(j-1));
}

template<class T>
T matT<T>::sub(long i, long j) const
{
  return entries.at((i-1)*nco+(j-1));
}

template<class T>
matT<T> matT<T>::slice(long r1,long r2,long c1,long c2) const
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
 matT<T> ans(n,c);
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

template<class T>
matT<T>& matT<T>::operator=(const matT<T>& m)
{
 if (this==&m) return *this;
 nro=m.nro;
 nco=m.nco;
 entries = m.entries;
 return *this;
}

template<class T>
void matT<T>::set(long i, long j, const T& x)
{
  entries.at((i-1)*nco+(j-1)) = x;
}

template<class T>
void matT<T>::add(long i, long j, const T& x)
{
  if (is_nonzero(x)) entries.at((i-1)*nco+(j-1)) += x;
}

template<class T>
void matT<T>::setrow(long i, const vecT<T>& v)
{
  std::copy(v.entries.begin(), v.entries.end(), entries.begin() + (i-1)*nco);
}

template<class T>
void matT<T>::setcol(long j, const vecT<T>& v)
{
  auto colj = entries.begin()+(j-1);
  for ( const auto vi : v.entries)
    {
      *colj = vi;
      colj += nco;
    }
}

template<class T>
vecT<T> matT<T>::row(long i) const
{
 vecT<T> mi(nco);
 auto e = entries.begin()+(i-1)*nco;
 std::copy(e, e+nco, mi.entries.begin());
 return mi;
}

template<class T>
vecT<T> matT<T>::col(long j) const
{
 vecT<T> v(nro);
 auto entriesij = entries.begin()+(j-1);
 for ( auto& vi : v.entries)
   {
     vi = *entriesij;
     entriesij+=nco;
   }
 return v;
}

template<class T>
void matT<T>::swaprows(long r1, long r2)
{
  auto mr1 = entries.begin() + (r1-1)*nco;
  auto mr2 = entries.begin() + (r2-1)*nco;
  std::swap_ranges(mr1, mr1+nco, mr2);
}

template<class T>
void matT<T>::multrow(long r, const T& scal)
{
  if (is_one(scal)) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [scal](const T& x) {return x * scal;});
}

template<class T>
void matT<T>::divrow(long r, const T& scal)
{
  if (is_zero(scal)||is_one(scal)) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [scal](const T& x) {return x / scal;});
}

template<class T>
T matT<T>::content() const
{
  return std::accumulate(entries.begin(), entries.end(), T(0),
                         [](const T& x, const T& y) {return gcd(x,y);});
}

template<class T>
T matT<T>::row_content(long r) const
{
  auto mij = entries.begin()+(r-1)*nco;
  return std::accumulate(mij, mij+nco, T(0),
                         [](const T& x, const T& y) {return gcd(x,y);});
}

template<class T>
void matT<T>::clearrow(long r)
{
  divrow(r, row_content(r));
}

template<class T>
void matT<T>::make_primitive()
{
  T g = content();
  if (is_zero(g)||is_one(g)) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [g](const T& x) {return x / g;});
}

template<class T>
void matT<T>::operator+=(const matT<T>& n)
{
  std::transform(n.entries.begin(), n.entries.end(), entries.begin(), entries.begin(),
                 [](const T& x, const T& y) { return x + y;});
}

template<class T>
void matT<T>::operator-=(const matT<T>& n)
{
  std::transform(n.entries.begin(), n.entries.end(), entries.begin(), entries.begin(),
                 [](const T& x, const T& y) { return y - x;});
}

template<class T>
void matT<T>::operator*=(const T& scal)
{
  if (is_one(scal))
    return;
  if (is_zero(scal))
    std::fill(entries.begin(), entries.end(), T(0));
  else
    std::transform(entries.begin(), entries.end(), entries.begin(),
                   [scal](const T& x) {return x * scal;});
}

template<class T>
void matT<T>::operator/=(const T& scal)
{
  if (is_zero(scal)||is_one(scal)) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const T& x) {return x / scal;});
}

// Definitions of non-member, friend operators and functions

// add/sub row i of mat to v
template<class T>
void add_row_to_vec(vecT<T>& v, const matT<T>& m, long i)
{
  std::transform(v.entries.begin(), v.entries.end(),
                 m.entries.begin()+(i-1)*m.nco,
                 v.entries.begin(), std::plus<T>());
}

template<class T>
void sub_row_to_vec(vecT<T>& v, const matT<T>& m, long i)
{
  std::transform(v.entries.begin(), v.entries.end(),
                 m.entries.begin()+(i-1)*m.nco,
                 v.entries.begin(), std::minus<T>());
}

template<class T>
matT<T> operator*(const matT<T>& m1, const matT<T>& m2)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 matT<T> m3(m,p);
 if (n==m2.nro)
   {
     auto a=m1.entries.begin();                                     // a points to m1(i,k)
     for (auto c=m3.entries.begin(); c!=m3.entries.end(); c+=p)     // c points to m3(i,_) for 0<=i<m
       {
         for (auto b=m2.entries.begin(); b!=m2.entries.end(); b+=p) // b points to m2(k,_) for 0<=k<n
           { // add m1(i,k)*m2(k,j) to m3(i,j) for 0<=j<p
             T m1ik = *a++;
             std::transform(b, b+p, c, c,
                            [m1ik] (const T& m2kj, const T& m3ij) {return m1ik*m2kj+m3ij;});
           }
       }
   }
 else
   {
     cerr << "Incompatible sizes in mat product"<<endl;
   }
 return m3;
}

template<class T>
int operator==(const matT<T>& m1, const matT<T>& m2)
{
  return (m1.nro==m2.nro) && (m1.nco==m2.nco) && (m1.entries==m2.entries);
}

template<class T>
void matT<T>::output(ostream& s) const
{
  auto mij=entries.begin();
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

template<class T>
void matT<T>::output_pari(ostream& s) const
{
  auto mij=entries.begin();
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

template<class T>
long ndigits(const T& a)
{
  int digits = 0;
  T aa(a);
  if (aa < 0) digits = 1; // for the '-'
  while (is_nonzero(aa)) { aa /= 10; digits++; }
  return digits;
}

template<class T>
void matT<T>::output_pretty(ostream& s) const
{
  // find max ndgits in each column:
  vector<int> colwidths(nco);
  for(long j=0; j<nco; j++)
    {
      auto mij = entries.begin()+j;
      T ma(0), mi(0); // max and min for column j
      for(long i=0; i<nro; i++, mij+=nco)
	{
	  if (*mij>ma) ma=*mij;
	  else if (*mij<mi) mi=*mij;
	}
      long nma=ndigits(ma),
        nmi=ndigits(mi);
      if(nmi>nma)nma=nmi;
      colwidths[j]=nma;
    }
  long nr=nro;
  auto mij=entries.begin();
  while(nr--)
    {
      s << "[";
      for(long j=0; j<nco; j++)
	{
	  if(j) s<<" ";
	  s.width(colwidths[j]);
          s<<(*mij++);
	}
      s<<"]\n";
    }
}

// The binary file input/output only works for T=int or long, not
// bigint, and is not used anyehere else
// void mat::dump_to_file(string filename) const
// {
//   ofstream fout(filename.c_str(),ofstream::binary);
//   fout.write((char*)&nro,sizeof(nro));
//   fout.write((char*)&nco,sizeof(nco));
//   fout.write((char*)entries.data(),nro*nco*sizeof(T));
//   fout.close();
// }

// void mat::read_from_file(string filename)
// {
//   ifstream fin(filename.c_str());
//   fin.read((char*)&nro,sizeof(nro));
//   fin.read((char*)&nco,sizeof(nco));
//   entries.resize(nro*nco);
//   fin.read((char*)entries.data(),nro*nco*sizeof(T));
//   fin.close();
// }

template<class T>
istream& operator>>(istream& s, matT<T>& m) // m cannot be const
{
 long n=m.nro*m.nco;
 auto mij=m.entries.begin();
 while(n--) s >> (*mij++);
 return s;
}

template<class T>
matT<T> colcat(const matT<T>& a, const matT<T>& b)
{
 long nr = a.nro, nca = a.nco, ncb = b.nco;
 matT<T> c(nr,nca+ncb);
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
   cerr << "colcat: matrices have different number of rows!" << endl;
 return c;
}

template<class T>
matT<T> rowcat(const matT<T>& a, const matT<T>& b)
{
 matT<T> c(a.nro+b.nro,a.nco);
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

template<class T>
matT<T> directsum(const matT<T>& a, const matT<T>& b)
{
  return rowcat(colcat(a,matT<T>(a.nro,b.nco)),colcat(matT<T>(b.nro,a.nco),b));
}

//plain elimination, no clearing
template<class T>
void elimrows(matT<T>& m, long r1, long r2, long pos) // m cannot be const
{
  long nc=m.nco;
  T p = m(r1,pos), q=m(r2,pos);
  auto mr1 = m.entries.begin() + (r1-1)*nc;
  auto mr2 = m.entries.begin() + (r2-1)*nc;
  // replace row2 by p*row2-q*row1
  std::transform(mr1, mr1+nc, mr2, mr2,
                 [p,q] (const T& x, const T& y) {return p*y-q*x;});
}

//elimination + clearing (i.e. divide new row by its content)
template<class T>
void elimrows1(matT<T>& m, long r1, long r2, long pos)
{
  elimrows(m,r1,r2,pos);
  m.clearrow(r2);
}

//elimination + divide by last pivot
template<class T>
void elimrows2(matT<T>& m, long r1, long r2, long pos, const T& last)
{
  elimrows(m,r1,r2,pos);
  m.divrow(r2,last);
}

// Definition of non-friend functions

template<class T>
matT<T> operator+(const matT<T>& m)
{
  return m;
}

template<class T>
matT<T> operator-(const matT<T>& m)
{
  return T(-1)*m;
}

template<class T>
matT<T> operator+(const matT<T>& m1, const matT<T>& m2)
{
  matT<T> ans(m1); ans+=m2;  return ans;
}

template<class T>
matT<T> operator-(const matT<T>& m1, const matT<T>& m2) 
{
  matT<T> ans(m1); ans-=m2;  return ans;
}

template<class T>
matT<T> operator*(const T& scal, const matT<T>& m)
{
  matT<T> ans(m); ans*=scal;  return ans;
}

template<class T>
matT<T> operator/(const matT<T>& m, const T& scal)
{
  matT<T> ans(m); ans/=scal;  return ans;
}

template<class T>
int operator!=(const matT<T>& m1, const matT<T>& m2)
{
  return !(m1==m2);
}

template<class T>
vecT<T> operator*(const matT<T>& m, const vecT<T>& v)
{
 long c=m.nco;
 vecT<T> w(m.nro);
 if (c==dim(v))
   {
     auto mi = m.entries.begin();
     for (auto& wi : w.entries)
       {
         wi = std::inner_product(mi, mi+c, v.entries.begin(), T(0));
         mi += c;
       }
   }
 else
   cerr << "Incompatible sizes in *(mat,vec)"<<endl;
 return w;
}

template<class T>
matT<T> matT<T>::scalar_matrix(long n, const T& a)
{
  matT<T> D(n,n);
  for (long i=1; i<=n; i++) D.set(i,i,a);
  return D;
}

template<class T>
matT<T> transpose(const matT<T>& m)
{
  long nr=m.ncols(), nc=m.nrows();
  matT<T> ans(nr, nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j,  m(j,i));
  return ans;
}

// submatrix of rows indexed by v, all columns
template<class T>
matT<T> rowsubmat(const matT<T>& m, const vecT<T>& v)
{
  long nr = dim(v), nc = m.ncols();
  matT<T> ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(v[i],j));
  return ans;
}

template<class T>
matT<T> rowsubmat(const matT<T>& m, const vecT<long>& v)
{
  long nr = dim(v), nc = m.ncols();
  matT<T> ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(v[i],j));
  return ans;
}

// submatrix of rows indexed by iv, columns indexed by jv
template<class T>
matT<T> submat(const matT<T>& m, const vecT<int>& iv, const vecT<int>& jv)
{
  long nr = dim(iv), nc = dim(jv);
  matT<T> ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(iv[i],jv[j]));
  return ans;
}

template<class T>
matT<T> submat(const matT<T>& m, const vecT<long>& iv, const vecT<long>& jv)
{
  long nr = dim(iv), nc = dim(jv);
  matT<T> ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(iv[i],jv[j]));
  return ans;
}

template<class T>
matT<T> echelon(const matT<T>& entries, vecT<int>& pcols, vecT<int>& npcols,
            long& rk, long& ny, T& d, int method)
{
  switch (method)
    {
    case 0: default: return echelon0(entries,pcols,npcols,rk,ny,d);
    case 2: return echelonp(entries,pcols,npcols,rk,ny,d, T(DEFAULT_MODULUS));
    }
}

//#define DEBUG_ECH_0

//N.B. if(q==0) the following multiplies row r2 by p, which looks
//redundant.  However, it is important to keep this in as in echelon0
//we must guarentee divisibility by "lastpivot".  We do not want to keep
//computing contents of rows as this is slower.
// Used in forward elimination in echelon0

template<class T>
void conservative_elim(vector<T>& m, long nc, long r1, long r2, long pos)
{
  auto mr1=m.begin() + r1*nc + pos;
  auto mr2=m.begin() + r2*nc + pos;
  T p = *mr1, q = *mr2;
  nc -= pos;
#ifdef DEBUG_ECH_0
  cout<<"In conservative_elim with p = "<<p<<" and q = " << q << endl;
  cout<<"row 1: "; for(long n=0; n<nc; n++) cout<<*(mr1+n)<<",";  cout<<endl;
  cout<<"row 2: "; for(long n=0; n<nc; n++) cout<<*(mr2+n)<<",";  cout<<endl;
#endif
  if (is_one(p)&&is_zero(q))
    return;
  // generic function to make y (entry in row2) 0
  std::function<T (const T&, const T&)>
    f = [p,q](const T& x, const T& y) {return p*y - q*x;};
  if(is_one(p)) // now q!=0
    {
      if(is_one(q))
        f = [p,q](const T& x, const T& y) {return y - x;};
      else
        {
          if(is_one(-q))
            f = [p,q](const T& x, const T& y) {return y + x;};
          else
            f = [p,q](const T& x, const T& y) {return y - q*x;};
        }
    }
  else  // p!=1
    {
      if(is_zero(q))
        f = [p,q](const T& x, const T& y) {return p*y;};
      if(is_one(q))
        f = [p,q](const T& x, const T& y) {return p*y - x;};
      if(is_one(-q))
        f = [p,q](const T& x, const T& y) {return p*y + x;};
    }
  std::transform(mr1, mr1+nc, mr2, mr2, f);
}

// This version does not multiply row r1 by p unnecessarily.  Used in
// back substitution, it does not assume that the entries in
// columns<pos are 0.

template<class T>
void elim(vector<T>& m, long nc, long r1, long r2, long pos)
{
  auto mr1=m.begin()+r1*nc;
  auto mr2=m.begin()+r2*nc;
  T p = *(mr1+pos), q = *(mr2+pos);
#ifdef DEBUG_ECH_0
  cout<<"In elim with p = "<<p<<" and q = " << q << endl;
  cout<<"row 1: "; for(long n=0; n<nc; n++) cout<<*(mr1+n)<<",";  cout<<endl;
  cout<<"row 2: "; for(long n=0; n<nc; n++) cout<<*(mr2+n)<<",";  cout<<endl;
#endif
  if (is_one(p)&&is_zero(q))
    return;
  // generic function to make y (entry in row2) 0
  std::function<T (const T&, const T&)>
    f = [p,q](const T& x, const T& y) {return p*y - q*x;};
  if(is_one(p)) // now q!=0
    {
      if(is_one(q))
        f = [p,q](const T& x, const T& y) {return y - x;};
      else
        {
          if(is_one(-q))
            f = [p,q](const T& x, const T& y) {return y + x;};
          else
            f = [p,q](const T& x, const T& y) {return y - q*x;};
        }
    }
  else  // p!=1
    {
      if(is_one(q))
        f = [p,q](const T& x, const T& y) {return p*y - x;};
      if(is_one(-q))
        f = [p,q](const T& x, const T& y) {return p*y + x;};
    }
  std::transform(mr1, mr1+nc, mr2, mr2, f);
}

template<class T>
void clear(vector<T>& row, long col1, long col2)
{
  auto row1=row.begin()+col1;
  auto row2=row.begin()+col2;
  T g = std::accumulate(row1, row2, T(0),
                             [](const T& x, const T& y) {return gcd(x,y);});
  if (g>1)
    std::for_each(row1, row2, [g](T& x) {x/=g;});
}

//#ifndef DEBUG_ECH_0
//#define DEBUG_ECH_0
//#endif

#ifdef DEBUG_ECH_0
template<class T>
void show(vector<T> m, long nr, long nc)
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

template<class T>
matT<T> echelon0(const matT<T>& entries, vecT<int>& pc, vecT<int>& npc,
             long& rk, long& ny, T& d)
{
#ifdef DEBUG_ECH_0
  cout<<"In echelon0 with matrix:\n"<<entries<<endl;
#endif
  rk=0; ny=0;
  T lastpivot(1);
  long r=0, nc=entries.nco, nr=entries.nro;
  vector<T> m = entries.entries;
  vector<int> pcols(nc), npcols(nc);
  for (long c=0; (c<nc)&&(r<nr); c++)
    {
      auto mij=m.begin()+r*nc+c;  // points to column c in row r
      T piv = abs(*mij);
      long rmin = r;
      mij+=nc;
      for (long r2=r+1; (r2<nr)&&(piv!=1); r2++, mij+=nc)
       {
         T mr2c = abs(*mij);
         if ((0<mr2c) && ((mr2c<piv) || (piv==0)))
           {
             piv=mr2c;
             rmin=r2;
           }
       }
      if (piv==0)
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
            conservative_elim(m,nc,r,r3,c);
#ifdef DEBUG_ECH_0
            cout<<"After, m is\n"; show(m,nr,nc);
#endif
	    if(lastpivot>1)
	      {
		auto mi1 = m.begin()+r3*nc;
                std::transform(mi1, mi1+nc, mi1, [lastpivot]( const T& x) {return x/lastpivot;});
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
#ifdef DEBUG_ECH_0
  cout<<"After forward elimination, rank = "<<rk<<"; pivots are:"<<endl;
  for(long r3=0; r3<rk; r3++) cout<<*(m.begin()+r3*nc+pcols[r3])<<",";
  cout<<endl;
#endif
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
          T fac = d/mij[pcols[r1]];
          std::transform(mij, mij+nc, mij, [fac](const T& x){return fac*x;});
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
  matT<T> ans(rk,nc, m);
  return ans;
}

template<class T>
long matT<T>::rank() const
{
  long rk=0;
  T lastpivot(1);
  matT<T> m(*this); // work with a copy, which will be reduced
  long nc=m.ncols(), nr=m.nrows();
  for (long c=1, r=1; (c<=nc)&&(r<=nr); c++)
    {
      T mmin = abs(m(r,c));
      long rmin = r;
      for (long r2=r+1; (r2<=nr)&&(!is_one(mmin)); r2++)
        {
          T mr2c = abs(m(r2,c));
          if ((is_nonzero(mr2c)) && ((mr2c<mmin) || (is_zero(mmin))))
            {
              mmin=mr2c;
              rmin=r2;
            }
        }
      if (mmin!=0)
        {
          rk++;
          if (rmin>r) m.swaprows(r,rmin);
          for (long r3 = r+1 ; r3<=nr; r3++)
            elimrows2(m,r,r3,c,lastpivot);
          lastpivot=mmin;
          r++;
        }
    }
  return rk;
}

template<class T>
long matT<T>::nullity() const
{
 return nco-rank();
}

template<class T>
T matT<T>::trace() const
{
  T tr(0);
  for (long i=0; i<nro; i++)
    tr += entries.at(i*(nco+1));
  return tr;
}

// FADEEV'S METHOD

template<class T>
vector<T> matT<T>::charpoly() const
{ long n = nrows();
  matT<T> b(*this);
  matT<T> id(identity_matrix(n));
  vector<T> clist(n+1);
  T t = trace();
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

template<class T>
T matT<T>::determinant() const
{
 T det = charpoly()[0];
 return (nro%2? -det :det);
}

template<class T>
void vecT<T>::sub_row(const matT<T>& m, int i)
{
  long n=entries.size();
  auto wi = m.entries.begin() + (i-1)*n;
  std::transform(entries.begin(), entries.end(), wi, entries.begin(), std::minus<T>());
}

template<class T>
void vecT<T>::add_row(const matT<T>& m, int i)
{
  long n=entries.size();
  auto wi = m.entries.begin() + (i-1)*n;
  std::transform(entries.begin(), entries.end(), wi, entries.begin(), std::plus<T>());
}

template<class T>
matT<T> addscalar(const matT<T>& mm, const T& c)
{
  return mm + matT<T>::scalar_matrix(mm.nrows(), c);
}

template<class T>
vecT<T> apply(const matT<T>& m, const vecT<T>& v)    // same as *(mat, vec)
{
  return m*v;
}

template<class T>
void matT<T>::reduce_mod_p(const T& p)
{
  if (p==0) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [p](const T& mij) {return mod(mij,p);});
}

template<class T>
void elimp(matT<T>& m, long r1, long r2, long pos, const T& pr)
{
  long nc=m.nco;
  auto mr1 = m.entries.begin() + (r1-1)*nc + (pos-1);
  auto mr2 = m.entries.begin() + (r2-1)*nc + (pos-1);
  T p = mod(*mr1,pr), q=mod(*mr2,pr);
  if(q==0) {return;} // nothing to do
  nc -= (pos-1); // first pos-1 entries are assumed 0 already
  // generic function to make y (entry in row2) 0
  std::function<T (const T&, const T&)>
    f = [pr,p,q](const T& x, const T& y) {return mod(xmodmul(p,y,pr)-xmodmul(q,x,pr), pr);};
  // simpler special cases (for same signature they must also capture both p and q)
  if(is_one(p))
   {
     if(is_one(q))
       f = [pr,p,q](const T& x, const T& y) {return mod(y-x, pr);};
     else
       {
         if(is_one(-q))
           f = [pr,p,q](const T& x, const T& y) {return mod(y+x, pr);};
         else
           // general q
           f = [pr,p,q](const T& x, const T& y) {return mod(y-xmodmul(q,x,pr), pr);};
       }
   }
  else // general p!=1
    {
      if(is_one(q))
        f = [pr,p,q](const T& x, const T& y) {return mod(xmodmul(p,y,pr)-x, pr);};
      if(is_one(-q))
        f = [pr,p,q](const T& x, const T& y) {return mod(xmodmul(p,y,pr)+x, pr);};
      // else the generic f will be used
    }
  std::transform(mr1, mr1+nc, mr2, mr2, f);
}

template<class T>
void elimp1(matT<T>& m, long r1, long r2, long pos, const T& pr)
//same as elimp except assumes pivot is 1
{
  long nc=m.nco;
  auto mr1 = m.entries.begin() + (r1-1)*nc + (pos-1);
  auto mr2 = m.entries.begin() + (r2-1)*nc + (pos-1);
  T q=mod(*mr2,pr);
  if(is_zero(q)) return;
  nc -= (pos-1); // first pos-1 entries are assumed 0 already
  // generic function to make y (entry in row2) 0
  std::function<T (const T&, const T&)>
    f = [pr,q](const T& x, const T& y) {return mod(y-xmodmul(q,x,pr), pr);};
  // simpler special cases
  if (is_one(q))
    f = [pr,q](const T& x, const T& y) {return mod(y-x, pr);};
  if (is_one(-q))
    f = [pr,q](const T& x, const T& y) {return mod(y+x, pr);};
  std::transform(mr1, mr1+nc, mr2, mr2, f);
}

//#define TRACE 1

// This method uses mod-p arithmetic internally but returns the
// "characteristic zero" echelon form of the mat.  It will only give
// the wrong answer if (a) the rank mod pr is not the actual rank, or (b)
// the actual echelon form has entries which are too big.

template<class T>
matT<T> echelonp(const matT<T>& entries, vecT<int>& pcols, vecT<int>& npcols,
             long& rk, long& ny, T& d, const T& pr)
{
#ifdef TRACE
  cout << "In echelonp\n";
#endif /* TRACE */
 long nr=entries.nrows(), nc=entries.ncols();
 matT<T> m(nr,nc);
 std::transform(entries.entries.begin(), entries.entries.end(), m.entries.begin(),
                [pr] (const T& x) {return mod(x,pr);});
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0;
 long r=1;
 for (long c=1; (c<=nc)&&(r<=nr); c++)
 {
   T mmin = m(r,c);
   long rmin = r;
   for (long r2=r+1; (r2<=nr)&&(mmin==0); r2++)
   {
     T mr2c = m(r2,c);
     if (0!=mr2c)
       {
         mmin=mr2c;
         rmin=r2;
       }
   }
   if (mmin==0)
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
       T fac = xmod(invmod(m(r1,pcols[r1]),pr),pr);
       for (long c=1; c<=nc; c++)
         m(r1,c)=xmodmul(fac,m(r1,c),pr);
     }
 }
 else
   for (long i=1; i<=rk; i++)
     for (long j=1; j<=nc; j++)
       m(i,j)=(j==pcols[i]);    // 0 or 1 !

#ifdef TRACE
 cout << "Finished second stage.\n Echelon mat mod "<<pr<<" is:\n";
 cout << m;
 cout << "Now lifting back to Q.\n";
#endif /* TRACE */
 T dd(1);
 matT<T> nmat(rk,nc);
 matT<T> dmat(rk,nc);

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
         T n1,d1;
         long jj = npcols[j];
         int ok = modrat(m(i,jj), pr,n1,d1);
         nmat(i,jj)=n1;
         dmat(i,jj)=d1;
         if (ok)
           dd=(dd*d1)/gcd(dd,d1);
         else
           cerr<<"Failed to lift "<<m(i,jj)<<" mod "<<pr<<" to Q"<<endl;
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

template<class T>
matT<T> echmodp(const matT<T>& entries, vecT<int>& pcols, vecT<int>& npcols, long& rk, long& ny, const T& pr)
{
 // cout << "In echmodp with p="<<pr<<" and matrix " << entries << endl;
 long nr=entries.nrows(), nc=entries.ncols();
 matT<T> m(nr,nc);
 std::transform(entries.entries.begin(), entries.entries.end(), m.entries.begin(),
                [pr] (const T& x) {return mod(x,pr);});
 // cout << " - after reducing modulo p,  matrix is " << m << endl;
 pcols.init(nc);
 npcols.init(nc);
 rk=ny=0;
 long r=1;
 for (long c=1; (c<=nc)&&(r<=nr); c++)
   {
     auto mij=m.entries.begin()+(r-1)*nc+c-1;
     T mmin(*mij);
     long rmin = r;
     mij += nc;
     for (long r2=r+1; (r2<=nr)&&(is_zero(mmin)); r2++, mij+=nc)
       {
	 T mr2c(*mij);
	 if (is_nonzero(mr2c))
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
         // cout<<"c = "<<c<<", pivot = "<<mmin<<endl;
         T fac = xmod(invmod(mmin,pr),pr);
         std::transform(entriesij, entriesij+nc, entriesij,
                        [pr,fac] (const T& x) {return mod(xmodmul(fac,x, pr), pr);});
         for (long r3 = r+1 ; r3<=nr; r3++)
           elimp1(m,r,r3,c,pr);
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
	 T fac = *(mij+pcols[r1]-1);
	 fac = mod(invmod(fac,pr),pr);
         std::transform(mij, mij+nc, mij,
                        [pr,fac] (const T& x) {return mod(xmodmul(fac,x, pr), pr);});
       }
   }
 else
   {
     auto mij=m.entries.begin();
     for (long i=1; i<=rk; i++)
       for (long j=1; j<=nc; j++)
	 *mij++ = T(j==pcols[i]);    // 0 or 1 !
   }
 return m.slice(rk,nc);
}

template<class T>
matT<T> echmodp_uptri(const matT<T>& entries, vecT<int>& pcols, vecT<int>& npcols,
                                  long& rk, long& ny, const T& pr)
{
// cout << "In echmodp_uptri with matrix = " << entries;
 long nr=entries.nrows(), nc=entries.ncols();
 matT<T> m(nr,nc);
 std::transform(entries.entries.begin(), entries.entries.end(), m.entries.begin(),
                [pr] (const T& x) {return mod(x,pr);});
 pcols.init(nc);
 npcols.init(nc);
 rk=ny=0;
 long r=1;
 for (long c=1; (c<=nc)&&(r<=nr); c++)
   {
     auto mij=m.entries.begin()+(r-1)*nc+c-1;
     T mmin = *mij;
     long rmin = r;
     mij += nc;
     for (long r2=r+1; (r2<=nr)&&(mmin==0); r2++, mij+=nc)
       {
	 T mr2c = *mij;
	 if (0!=mr2c)
           {
             mmin=mr2c;
             rmin=r2;
           }
       }
     if (mmin==0)
       npcols[++ny] = c;
     else
       {
	 pcols[++rk] = c;
	 if (rmin>r)
           m.swaprows(r,rmin);
	 auto entriesij = m.entries.begin()+(r-1)*nc;
         T fac = mod(invmod(mmin,pr),pr);
         std::transform(entriesij, entriesij+nc, entriesij,
                        [pr,fac] (const T& x) {return mod(fac*x, pr);});
         for (long r3 = r+1 ; r3<=nr; r3++)
           elimp1(m,r,r3,c,pr);
	 r++;
       }
   }
 for (long c = rk+ny+1; c<=nc; c++)
   npcols[++ny] = c ;
 pcols  =  pcols.slice(rk);
 npcols =  npcols.slice(ny);    // truncate index vectors
 // cout << "Rank = " << rk << ".  Nullity = " << ny << ".\n";
 return m.slice(rk,nc);
}

//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Interface with NTL matrices
//
//////////////////////////////////////////////////////////////////////////////////////////////

//#define TRACE_NTL_REF

#include <NTL/mat_lzz_p.h>
#ifdef TRACE_NTL_REF
#include <eclib/timer.h>
#endif

// Construct an NTL mat_lzz_p (matrix mod p) from a mat mod pr

template<class T>
mat_zz_p mat_zz_p_from_mat(const matT<T>& M, const T& pr)
{
  long nr=M.nrows(), nc=M.ncols();
#ifdef TRACE_NTL_REF
  cout<<"Creating an NTL mat_zz_p from a matrix with " << nr <<" rows and "<<nc<<" columns, mod "<<pr<<endl;
#endif
  // create NTL matrix copy of M:
  zz_pPush push(I2long(pr));
  mat_zz_p A(INIT_SIZE, nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      A.put(i,j, conv<zz_p>(M(i+1,j+1)));
#ifdef TRACE_NTL_REF
  cout<<"--done."<<endl;
#endif
  return A;
}

// Construct a mat (T type same as pr) from an NTL mat_lzz_p

template<class T>
matT<T> mat_from_mat_zz_p(const mat_zz_p& A, const T& pr) // type of T fixes return type
{
 long nr = A.NumRows(), nc = A.NumCols();
#ifdef TRACE_NTL_REF
  cout<<"Creating a mat from an NTL mat_zz_p with " << nr <<" rows and "<<nc<<" columns, mod "<<pr<<endl;
#endif
 // create matrix copy of A:
 matT<T> M(nr, nc);
 for(long i=0; i<nr; i++)
   for(long j=0; j<nc; j++)
     M(i+1,j+1) = mod(conv<T>(A.get(i,j)), pr);
#ifdef TRACE_NTL_REF
  cout<<"--done."<<endl;
#endif
 return M;
}

// compute ref of M mod pr via NTL, setting rk=rank, ny=nullity,
// pivotal columns pcols, non-pivotal columns npcols

template<class T>
matT<T> ref_via_ntl(const matT<T>& M, vecT<int>& pcols, vecT<int>& npcols,
                    long& rk, long& ny, const T& pr)
{
 long nc=M.ncols();
 long i, j, k;
#ifdef TRACE_NTL_REF
 timer ntl_timer;
 ntl_timer.start();
#endif
 zz_pPush push(I2long(pr));
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
 matT<T> ans = mat_from_mat_zz_p(A, pr).slice(rk,nc);
#ifdef TRACE_NTL_REF
 ntl_timer.start();
 ntl_timer.show();
 cout<<endl;
#endif
 return ans;
}

template<class T>
long rank_via_ntl(const matT<T>& M, const T& pr)
{
#ifdef TRACE_NTL_REF
  cout << "Computing rank mod "<<pr<<" of a matrix of size ("<<M.nrows()<<", "<<M.ncols()<<")..."<<flush;
  timer ntl_timer;
  ntl_timer.start();
#endif
  zz_pPush push(I2long(pr));
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

template<class T>
T det_via_ntl(const matT<T>& M, const T& pr)
{
#ifdef TRACE_NTL_REF
  cout << "Computing determinant mod "<<pr<<" of a matrix of size ("<<M.nrows()<<", "<<M.ncols()<<")..."<<flush;
  timer ntl_timer;
  ntl_timer.start();
#endif
  zz_pPush push(I2long(pr));
  mat_zz_p A = mat_zz_p_from_mat(M, pr);
  zz_p det = determinant(A);
#ifdef TRACE_NTL_REF
  cout << "done: "<<flush;
  ntl_timer.start();
  ntl_timer.show();
  cout<<endl;
#endif
  return mod(conv<T>(det), pr);
}
