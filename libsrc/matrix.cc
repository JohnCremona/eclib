// matrix.cc: manage implementation of integer matrix classes
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

#include "eclib/convert.h"
#include "eclib/linalg.h"
#include "eclib/polys.h"

#include "eclib/flinterface.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>

// Definitions of member operators and functions:

template<class T>
void Zmat<T>::init(long nr, long nc) // resets to zero mat of given size;
{                                // with defaults (0,0) releases all space.
  nro = nr;
  nco = nc;
  entries.resize(nro*nco, T(0));
}

template<class T>
T& Zmat<T>::operator()(long i, long j)   // returns reference to (i,j) entry
{
  return entries.at((i-1)*nco+(j-1));
}

template<class T>
T Zmat<T>::operator()(long i, long j) const   // returns (i,j) entry
{
  return entries.at((i-1)*nco+(j-1));
}

template<class T>
T Zmat<T>::sub(long i, long j) const
{
  return entries.at((i-1)*nco+(j-1));
}

template<class T>
Zmat<T> Zmat<T>::slice(long r1,long r2,long c1,long c2) const
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
 Zmat<T> ans(n,c);
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
Zmat<T>& Zmat<T>::operator=(const Zmat<T>& m)
{
 if (this==&m) return *this;
 nro=m.nro;
 nco=m.nco;
 entries = m.entries;
 return *this;
}

template<class T>
void Zmat<T>::set(long i, long j, const T& x)
{
  entries.at((i-1)*nco+(j-1)) = x;
}

template<class T>
void Zmat<T>::add(long i, long j, const T& x)
{
  if (is_nonzero(x)) entries.at((i-1)*nco+(j-1)) += x;
}

template<class T>
void Zmat<T>::setrow(long i, const Zvec<T>& v)
{
  std::copy(v.entries.begin(), v.entries.end(), entries.begin() + (i-1)*nco);
}

template<class T>
void Zmat<T>::setcol(long j, const Zvec<T>& v)
{
  auto colj = entries.begin()+(j-1);
  for ( const auto vi : v.entries)
    {
      *colj = vi;
      colj += nco;
    }
}

template<class T>
Zvec<T> Zmat<T>::row(long i) const
{
 Zvec<T> mi(nco);
 auto e = entries.begin()+(i-1)*nco;
 std::copy(e, e+nco, mi.entries.begin());
 return mi;
}

template<class T>
Zvec<T> Zmat<T>::col(long j) const
{
 Zvec<T> v(nro);
 auto entriesij = entries.begin()+(j-1);
 for ( auto& vi : v.entries)
   {
     vi = *entriesij;
     entriesij+=nco;
   }
 return v;
}

template<class T>
void Zmat<T>::swaprows(long r1, long r2)
{
  auto mr1 = entries.begin() + (r1-1)*nco;
  auto mr2 = entries.begin() + (r2-1)*nco;
  std::swap_ranges(mr1, mr1+nco, mr2);
}

template<class T>
void Zmat<T>::multrow(long r, const T& scal)
{
  if (is_one(scal)) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [scal](const T& x) {return x * scal;});
}

template<class T>
void Zmat<T>::divrow(long r, const T& scal)
{
  if (::is_zero(scal)||::is_one(scal)) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [r,scal](const T& x)
  {
    T q,rem;
    if (divrem(x,scal,q,rem)) return q;
    cerr << "Error in dividing row " << r << " by " << scal << endl;
    exit(1);
    return T(0);
  });
}

template<class T>
int Zmat<T>::is_zero() const
{
  return std::all_of(entries.begin(), entries.end(),
                     [](const T& x) {return ::is_zero(x);});
}

template<class T>
T Zmat<T>::content() const
{
  return std::accumulate(entries.begin(), entries.end(), T(0),
                         [](const T& x, const T& y) {return gcd(x,y);});
}

template<class T>
T Zmat<T>::row_content(long r) const
{
  auto mij = entries.begin()+(r-1)*nco;
  return std::accumulate(mij, mij+nco, T(0),
                         [](const T& x, const T& y) {return gcd(x,y);});
}

template<class T>
void Zmat<T>::clearrow(long r)
{
  divrow(r, row_content(r));
}

template<class T>
void Zmat<T>::make_primitive()
{
  T g = content();
  if (::is_zero(g)||is_one(g)) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [g](const T& x) {return x / g;});
}

template<class T>
void Zmat<T>::operator+=(const Zmat<T>& n)
{
  std::transform(n.entries.begin(), n.entries.end(), entries.begin(), entries.begin(),
                 [](const T& x, const T& y) { return x + y;});
}

template<class T>
void Zmat<T>::operator-=(const Zmat<T>& n)
{
  std::transform(n.entries.begin(), n.entries.end(), entries.begin(), entries.begin(),
                 [](const T& x, const T& y) { return y - x;});
}

template<class T>
void Zmat<T>::operator*=(const T& scal)
{
  if (is_one(scal))
    return;
  if (::is_zero(scal))
    std::fill(entries.begin(), entries.end(), T(0));
  else
    std::transform(entries.begin(), entries.end(), entries.begin(),
                   [scal](const T& x) {return x * scal;});
}

template<class T>
void Zmat<T>::operator/=(const T& scal)
{
  if (::is_zero(scal)||is_one(scal)) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const T& x) {return x / scal;});
}

// append n zero rows
template<class T>
void Zmat<T>::append_rows(int n)
{
  entries.insert(entries.end(), n*nco, T(0));
  nro += n;
}

// append given new row
template<class T>
void Zmat<T>::append_row(const Zvec<T>& new_row)
{
  entries.insert(entries.end(), new_row.entries.begin(), new_row.entries.end());
  nro += 1;
}

// append given new rows
template<class T>
void Zmat<T>::append_rows(const vector<Zvec<T>>& new_rows)
{
  for (auto v: new_rows)
    entries.insert(entries.end(), v.entries.begin(), v.entries.end());
  nro += new_rows.size();
}

// delete the last n rows
template<class T>
void Zmat<T>::delete_rows(int n)
{
  entries.erase(entries.end()-n*nco, entries.end());
  nro -= n;
}

template<class T>
void Zmat<T>::output_raw(ostream& s) const
{
  std::copy(entries.begin(), entries.end(), ostream_iterator<T>(s," "));
}

template<class T>
void Zmat<T>::output(ostream& s) const
{
  auto mij=entries.begin();
  s << "[";
  long nr=nro;
  while(nr--)
    {
      long nc=nco;
      s<<"[";
      while(nc--)
        {
          s<<(*mij++);
          if(nc)
            s<<",";
        }
      s<<"]";
      if(nr)
        s<<",\n";
    }
  s << "]";
}

// same as m.output(cout) except no newlines between rows
template<class T>
void Zmat<T>::output_flat(ostream&s) const
{
  auto mij = entries.begin();
  s << "[";
  long nr=nro;
  while(nr--)
    {
      long nc = nco;
      s<<"[";
      while(nc--)
        {
          s<<(*mij++);
          if(nc)
            s<<",";
        }
      s<<"]";
      if(nr)
        s<<",";
    }
  s << "]";
}

template<class T>
void Zmat<T>::output_pari(ostream& s) const
{
  auto mij=entries.begin();
  s << "[";
  long nr=nro;
  while(nr--)
    {
      long nc=nco;
      while(nc--) {s<<(*mij++); if(nc) s<<",";}
      if(nr) s<<";";
    }
  s << "]";
}

template<class T>
long ndigits(const T& a)
{
  int digits = (a<=0); // for the '-' or 0
  T aa(abs(a));
  while (is_nonzero(aa)) { aa /= 10; digits++; }
  return digits;
}

template<class T>
void Zmat<T>::output_pretty(ostream& s) const
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

template<class T>
Zmat<T> Zmat<T>::scalar_matrix(long n, const T& a)
{
  Zmat<T> D(n,n);
  for (long i=1; i<=n; i++) D.set(i,i,a);
  return D;
}

template<class T>
long Zmat<T>::rank() const
{
  fmpz_mat_t A;
  fmpz_mat_init(A, nro, nco);
  flint_mat_from_mat(A,*this);
  long rk = fmpz_mat_rank(A);
  fmpz_mat_clear(A);
  return rk;
}

template<class T>
long Zmat<T>::nullity() const
{
 return nco-rank();
}

template<class T>
T Zmat<T>::trace() const
{
  T tr(0);
  for (long i=0; i<nro; i++)
    tr += entries.at(i*(nco+1));
  return tr;
}

// FADEEV'S METHOD

// CHANGE
template<class T>
vector<T> Zmat<T>::charpoly() const
{ long n = nrows();
  Zmat<T> b(*this);
  Zmat<T> id(identity_matrix(n));
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
      exit(1);
    }
  return clist;
}

// CHANGE
template<class T>
T Zmat<T>::determinant() const
{
 T det = charpoly()[0];
 return (nro%2? -det :det);
}

template<class T>
void Zvec<T>::sub_row(const Zmat<T>& m, int i)
{
  long n=entries.size();
  auto wi = m.entries.begin() + (i-1)*n;
  std::transform(entries.begin(), entries.end(), wi, entries.begin(), std::minus<T>());
}

// The next two need to be instantated hree as they are methods of the
// Zvec class which is mostly instantiated in vector.cc
template<class T>
void Zvec<T>::add_row(const Zmat<T>& m, int i)
{
  long n=entries.size();
  auto wi = m.entries.begin() + (i-1)*n;
  std::transform(entries.begin(), entries.end(), wi, entries.begin(), std::plus<T>());
}
template void Zvec<int>::add_row(const Zmat<int>& m, int i);
template void Zvec<long>::add_row(const Zmat<long>& m, int i);
template void Zvec<ZZ>::add_row(const Zmat<ZZ>& m, int i);
template void Zvec<INT>::add_row(const Zmat<INT>& m, int i);

template<class T>
void Zmat<T>::reduce_mod_p(const T& p)
{
  if (p==0) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [p](const T& mij) {return mod(mij,p);});
}
template void Zvec<int>::sub_row(const Zmat<int>& m, int i);
template void Zvec<long>::sub_row(const Zmat<long>& m, int i);
template void Zvec<ZZ>::sub_row(const Zmat<ZZ>& m, int i);
template void Zvec<INT>::sub_row(const Zmat<INT>& m, int i);

// Definitions of non-member operators and functions

// add/sub row i of mat to v
template<class T>
void add_row_to_vec(Zvec<T>& v, const Zmat<T>& m, long i)
{
  std::transform(v.entries.begin(), v.entries.end(),
                 m.entries.begin()+(i-1)*m.nco,
                 v.entries.begin(), std::plus<T>());
}
template void add_row_to_vec<int>(Zvec<int>& v, const Zmat<int>& m, long i);
template void add_row_to_vec<long>(Zvec<long>& v, const Zmat<long>& m, long i);
template void add_row_to_vec<ZZ>(Zvec<ZZ>& v, const Zmat<ZZ>& m, long i);
template void add_row_to_vec<INT>(Zvec<INT>& v, const Zmat<INT>& m, long i);

template<class T>
void sub_row_to_vec(Zvec<T>& v, const Zmat<T>& m, long i)
{
  std::transform(v.entries.begin(), v.entries.end(),
                 m.entries.begin()+(i-1)*m.nco,
                 v.entries.begin(), std::minus<T>());
}
template void sub_row_to_vec<int>(Zvec<int>& v, const Zmat<int>& m, long i);
template void sub_row_to_vec<long>(Zvec<long>& v, const Zmat<long>& m, long i);
template void sub_row_to_vec<ZZ>(Zvec<ZZ>& v, const Zmat<ZZ>& m, long i);
template void sub_row_to_vec<INT>(Zvec<INT>& v, const Zmat<INT>& m, long i);

template<class T>
Zmat<T> operator*(const Zmat<T>& m1, const Zmat<T>& m2)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 Zmat<T> m3(m,p);
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
     exit(1);
   }
 return m3;
}
template Zmat<int> operator*<int>(const Zmat<int>&, const Zmat<int>&);
template Zmat<long> operator*<long>(const Zmat<long>&, const Zmat<long>&);
template Zmat<ZZ> operator*<ZZ>(const Zmat<ZZ>&, const Zmat<ZZ>&);
template Zmat<INT> operator*<INT>(const Zmat<INT>&, const Zmat<INT>&);

template<class T>
int operator==(const Zmat<T>& m1, const Zmat<T>& m2)
{
  return (m1.nro==m2.nro) && (m1.nco==m2.nco) && (m1.entries==m2.entries);
}
template int operator==<int>(const Zmat<int>&, const Zmat<int>&);
template int operator==<long>(const Zmat<long>&, const Zmat<long>&);
template int operator==<ZZ>(const Zmat<ZZ>&, const Zmat<ZZ>&);
template int operator==<INT>(const Zmat<INT>&, const Zmat<INT>&);


template<class T>
int operator!=(const Zmat<T>& m1, const Zmat<T>& m2)
{
  return !(m1==m2);
}
template int operator!=<int>(const Zmat<int>& m1, const Zmat<int>& m2);
template int operator!=<long>(const Zmat<long>& m1, const Zmat<long>& m2);
template int operator!=<ZZ>(const Zmat<ZZ>& m1, const Zmat<ZZ>& m2);
template int operator!=<INT>(const Zmat<INT>& m1, const Zmat<INT>& m2);

template<class T>
istream& operator>>(istream& s, Zmat<T>& m) // m cannot be const
{
 long n=m.nro*m.nco;
 auto mij=m.entries.begin();
 while(n--) s >> (*mij++);
 return s;
}
template istream& operator>> <int>(istream&s, Zmat<int>&);
template istream& operator>> <long>(istream&s, Zmat<long>&);
template istream& operator>> <ZZ>(istream&s, Zmat<ZZ>&);
template istream& operator>> <INT>(istream&s, Zmat<INT>&);

template<class T>
Zmat<T> colcat(const Zmat<T>& a, const Zmat<T>& b)
{
 long nr = a.nro, nca = a.nco, ncb = b.nco;
 Zmat<T> c(nr,nca+ncb);
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
     exit(1);
   }
 return c;
}
template Zmat<int> colcat<int>(const Zmat<int>& a, const Zmat<int>& b);
template Zmat<long> colcat<long>(const Zmat<long>& a, const Zmat<long>& b);
template Zmat<ZZ> colcat<ZZ>(const Zmat<ZZ>& a, const Zmat<ZZ>& b);
template Zmat<INT> colcat<INT>(const Zmat<INT>& a, const Zmat<INT>& b);

template<class T>
Zmat<T> rowcat(const Zmat<T>& a, const Zmat<T>& b)
{
 Zmat<T> c(a.nro+b.nro,a.nco);
 if (a.nco==b.nco)
 {
   auto cij = c.entries.begin();
   std::copy(a.entries.begin(), a.entries.end(), cij);
   cij += a.entries.size();
   std::copy(b.entries.begin(), b.entries.end(), cij);
 }
 else
   {
     cerr << "rowcat: matrices have different number of columns!" << endl;
     exit(1);
   }
 return c;
}
template Zmat<int> rowcat<int>(const Zmat<int>& a, const Zmat<int>& b);
template Zmat<long> rowcat<long>(const Zmat<long>& a, const Zmat<long>& b);
template Zmat<ZZ> rowcat<ZZ>(const Zmat<ZZ>& a, const Zmat<ZZ>& b);
template Zmat<INT> rowcat<INT>(const Zmat<INT>& a, const Zmat<INT>& b);

template<class T>
Zmat<T> directsum(const Zmat<T>& a, const Zmat<T>& b)
{
  return rowcat(colcat(a,Zmat<T>(a.nro,b.nco)),colcat(Zmat<T>(b.nro,a.nco),b));
}
template Zmat<int> directsum<int>(const Zmat<int>& a, const Zmat<int>& b);
template Zmat<long> directsum<long>(const Zmat<long>& a, const Zmat<long>& b);
template Zmat<ZZ> directsum<ZZ>(const Zmat<ZZ>& a, const Zmat<ZZ>& b);
template Zmat<INT> directsum<INT>(const Zmat<INT>& a, const Zmat<INT>& b);


template<class T>
Zmat<T> operator+(const Zmat<T>& m)
{
  return m;
}
template Zmat<int> operator+<int>(const Zmat<int>&);                   // unary
template Zmat<long> operator+<long>(const Zmat<long>&);                   // unary
template Zmat<ZZ> operator+<ZZ>(const Zmat<ZZ>&);                   // unary
template Zmat<INT> operator+<INT>(const Zmat<INT>&);                   // unary

template<class T>
Zmat<T> operator-(const Zmat<T>& m)
{
  return T(-1)*m;
}
template Zmat<int> operator-<int>(const Zmat<int>&);                   // unary
template Zmat<long> operator-<long>(const Zmat<long>&);                   // unary
template Zmat<ZZ> operator-<ZZ>(const Zmat<ZZ>&);                   // unary
template Zmat<INT> operator-<INT>(const Zmat<INT>&);                   // unary

template<class T>
Zmat<T> operator+(const Zmat<T>& m1, const Zmat<T>& m2)
{
  Zmat<T> ans(m1); ans+=m2;  return ans;
}
template Zmat<int>  operator+<int> (const Zmat<int>&,  const Zmat<int>&);
template Zmat<long> operator+<long>(const Zmat<long>&, const Zmat<long>&);
template Zmat<ZZ>   operator+<ZZ>  (const Zmat<ZZ>&,   const Zmat<ZZ>&);
template Zmat<INT>  operator+<INT> (const Zmat<INT>&,  const Zmat<INT>&);

template<class T>
Zmat<T> operator-(const Zmat<T>& m1, const Zmat<T>& m2)
{
  Zmat<T> ans(m1); ans-=m2;  return ans;
}
template Zmat<int>  operator-<int> (const Zmat<int>&,  const Zmat<int>&);
template Zmat<long> operator-<long>(const Zmat<long>&, const Zmat<long>&);
template Zmat<ZZ>   operator-<ZZ>  (const Zmat<ZZ>&,   const Zmat<ZZ>&);
template Zmat<INT>  operator-<INT> (const Zmat<INT>&,  const Zmat<INT>&);

template<class T>
Zmat<T> operator*(const T& scal, const Zmat<T>& m)
{
  Zmat<T> ans(m); ans*=scal;  return ans;
}
template Zmat<int > operator*<int> (const int&,  const Zmat<int>&);
template Zmat<long> operator*<long>(const long&, const Zmat<long>&);
template Zmat<ZZ>   operator*<ZZ>  (const ZZ&,   const Zmat<ZZ>&);
template Zmat<INT>  operator*<INT> (const INT&,  const Zmat<INT>& m);

template<class T>
Zmat<T> operator/(const Zmat<T>& m, const T& scal)
{
  Zmat<T> ans(m); ans/=scal;  return ans;
}
template Zmat<int>  operator/<int> (const Zmat<int>&,  const int&);
template Zmat<long> operator/<long>(const Zmat<long>&, const long&);
template Zmat<ZZ>   operator/<ZZ>  (const Zmat<ZZ>&,   const ZZ&);
template Zmat<INT>  operator/<INT> (const Zmat<INT>&,  const INT&);

template<class T>
Zvec<T> operator*(const Zmat<T>& m, const Zvec<T>& v)
{
 long c=m.nco;
 Zvec<T> w(m.nro);
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
   {
     cerr << "Incompatible sizes in *(mat,vec)"<<endl;
     exit(1);
   }
 return w;
}
template Zvec<int> operator*<int>(const Zmat<int>&, const Zvec<int>&);
template Zvec<long> operator*<long>(const Zmat<long>&, const Zvec<long>&);
template Zvec<ZZ> operator*<ZZ>(const Zmat<ZZ>&, const Zvec<ZZ>&);
template Zvec<INT> operator*<INT>(const Zmat<INT>&, const Zvec<INT>&);

template<class T>
Zmat<T> transpose(const Zmat<T>& m)
{
  long nr=m.ncols(), nc=m.nrows();
  Zmat<T> ans(nr, nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j,  m(j,i));
  return ans;
}
template Zmat<int> transpose<int>(const Zmat<int>& m);
template Zmat<long> transpose<long>(const Zmat<long>& m);
template Zmat<ZZ> transpose<ZZ>(const Zmat<ZZ>& m);
template Zmat<INT> transpose<INT>(const Zmat<INT>& m);

// submatrix of rows indexed by v, all columns
template<class T>
Zmat<T> rowsubmat(const Zmat<T>& m, const Zvec<int>& v)
{
  long nr = dim(v), nc = m.ncols();
  Zmat<T> ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(v[i],j));
  return ans;
}
template Zmat<int> rowsubmat<int>(const Zmat<int>& m, const Zvec<int>& v);
template Zmat<long> rowsubmat<long>(const Zmat<long>& m, const Zvec<int>& v);
template Zmat<ZZ> rowsubmat<ZZ>(const Zmat<ZZ>& m, const Zvec<int>& v);
template Zmat<INT> rowsubmat<INT>(const Zmat<INT>& m, const Zvec<int>& v);

template<class T>
Zmat<T> rowsubmat(const Zmat<T>& m, const Zvec<long>& v)
{
  long nr = dim(v), nc = m.ncols();
  Zmat<T> ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(v[i],j));
  return ans;
}
template Zmat<int> rowsubmat<int>(const Zmat<int>& m, const Zvec<long>& v);
template Zmat<long> rowsubmat<long>(const Zmat<long>& m, const Zvec<long>& v);
template Zmat<ZZ> rowsubmat<ZZ>(const Zmat<ZZ>& m, const Zvec<long>& v);
template Zmat<INT> rowsubmat<INT>(const Zmat<INT>& m, const Zvec<long>& v);

// submatrix of rows indexed by iv, columns indexed by jv
template<class T>
Zmat<T> submat(const Zmat<T>& m, const Zvec<int>& iv, const Zvec<int>& jv)
{
  long nr = dim(iv), nc = dim(jv);
  Zmat<T> ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(iv[i],jv[j]));
  return ans;
}
template Zmat<int> submat<int>(const Zmat<int>& m, const Zvec<int>& iv, const Zvec<int>& jv);
template Zmat<long> submat<long>(const Zmat<long>& m, const Zvec<int>& iv, const Zvec<int>& jv);
template Zmat<ZZ> submat<ZZ>(const Zmat<ZZ>& m, const Zvec<int>& iv, const Zvec<int>& jv);
template Zmat<INT> submat<INT>(const Zmat<INT>& m, const Zvec<int>& iv, const Zvec<int>& jv);

template<class T>
Zmat<T> submat(const Zmat<T>& m, const Zvec<long>& iv, const Zvec<long>& jv)
{
  long nr = dim(iv), nc = dim(jv);
  Zmat<T> ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(iv[i],jv[j]));
  return ans;
}
template Zmat<int> submat<int>(const Zmat<int>& m, const Zvec<long>& iv, const Zvec<long>& jv);
template Zmat<long> submat<long>(const Zmat<long>& m, const Zvec<long>& iv, const Zvec<long>& jv);
template Zmat<ZZ> submat<ZZ>(const Zmat<ZZ>& m, const Zvec<long>& iv, const Zvec<long>& jv);
template Zmat<INT> submat<INT>(const Zmat<INT>& m, const Zvec<long>& iv, const Zvec<long>& jv);

template<class T>
Zmat<T> ref(const Zmat<T>& M, T& d,
            Zvec<int>& pcols, Zvec<int>& npcols,
            long& rk, long& ny)
{
  Zmat<T> R = ref(M, d);
  pnpcols(R, pcols, npcols, rk, ny);
  return R;
}
template Zmat<int> ref<int>(const Zmat<int>&, int&, Zvec<int>&, Zvec<int>&, long&, long&);
template Zmat<long> ref<long>(const Zmat<long>&, long&, Zvec<int>&, Zvec<int>&, long&, long&);
template Zmat<ZZ> ref<ZZ>(const Zmat<ZZ>&, ZZ&, Zvec<int>&, Zvec<int>&, long&, long&);
template Zmat<INT> ref<INT>(const Zmat<INT>&, INT&, Zvec<int>&, Zvec<int>&, long&, long&);

template<class T>
Zmat<T> addscalar(const Zmat<T>& mm, const T& c)
{
  return mm + Zmat<T>::scalar_matrix(mm.nrows(), c);
}
template Zmat<int> addscalar<int>(const Zmat<int>&, const int&);
template Zmat<long> addscalar<long>(const Zmat<long>&, const long&);
template Zmat<ZZ> addscalar<ZZ>(const Zmat<ZZ>&, const ZZ&);
template Zmat<INT> addscalar<INT>(const Zmat<INT>&, const INT&);

template<class T>
Zvec<T> apply(const Zmat<T>& m, const Zvec<T>& v)    // same as *(mat, vec)
{
  return m*v;
}
template Zvec<int> apply<int>(const Zmat<int>&, const Zvec<int>&);
template Zvec<long> apply<long>(const Zmat<long>&, const Zvec<long>&);
template Zvec<ZZ> apply<ZZ>(const Zmat<ZZ>&, const Zvec<ZZ>&);
template Zvec<INT> apply<INT>(const Zmat<INT>&, const Zvec<INT>&);

//CHANGE
// Assigns d*A^-1 to Ainv and returns d, assuming A square and det(A) nonzero
template<class T>
T inverse(const Zmat<T>& A, Zmat<T>& Ainv, int method)
{
  long n = A.nrows();
  const Zmat<T>& aug=colcat(A, Zmat<T>::identity_matrix(n));
  long rk, ny; Zvec<int> pc,npc; T d;
  const Zmat<T>& R = ref(aug, d, pc, npc, rk, ny);
  Ainv = R.slice(1,n,n+1,2*n);
  return d;
}
template int inverse<int>(const Zmat<int>&, Zmat<int>&, int);
template long inverse<long>(const Zmat<long>&, Zmat<long>&, int);
template ZZ inverse<ZZ>(const Zmat<ZZ>&, Zmat<ZZ>&, int);
template INT inverse<INT>(const Zmat<INT>&, Zmat<INT>&, int);


template<class T>
Zmat<T> matmulmodp(const Zmat<T>& m1, const Zmat<T>& m2, const T& pr)
{
 int m=m1.nro, n=m1.nco, p=m2.nco;
 Zmat<T> m3(m,p);
 if (n==m2.nro)
   {
     auto a=m1.entries.begin();                                     // a points to m1(i,k)
     for (auto c=m3.entries.begin(); c!=m3.entries.end(); c+=p)     // c points to m3(i,_) for 0<=i<m
       {
         for (auto b=m2.entries.begin(); b!=m2.entries.end(); b+=p) // b points to m2(k,_) for 0<=k<n
           { // add m1(i,k)*m2(k,j) to m3(i,j) for 0<=j<p
             T m1ik = *a++;
             std::transform(b, b+p, c, c,
                            [pr,m1ik] (const T& m2kj, const T& m3ij)
                            {return mod(xmodmul(m1ik,m2kj,pr)+m3ij, pr);});
           }
       }
   }
 else
   {
     cerr << "Incompatible sizes in mat product"<<endl;
     exit(1);
   }
 return m3;
}
template Zmat<int> matmulmodp<int>(const Zmat<int>&, const Zmat<int>&, const int& pr);
template Zmat<long> matmulmodp<long>(const Zmat<long>&, const Zmat<long>&, const long& pr);
template Zmat<ZZ> matmulmodp<ZZ>(const Zmat<ZZ>&, const Zmat<ZZ>&, const ZZ& pr);
template Zmat<INT> matmulmodp<INT>(const Zmat<INT>&, const Zmat<INT>&, const INT& pr);

template<class T>
Zvec<T> matvecmulmodp(const Zmat<T>& M, const Zvec<T>& v, const T& pr)
{
  int m=M.nro, n=M.nco;
  if (n!=dim(v))
   {
     cerr << "Incompatible sizes in mat*vec product: "<<m<<"x"<<n<<" and "<<dim(v)<<endl;
     exit(1);
     return Zvec<T>(m);
   }
  vector<T> ve = v.get_entries(), w;
  for (auto a=M.entries.begin(); a!=M.entries.end(); a+=n)
    {
      T wi(0);
      for (int j=0; j<n; j++)
        wi = xmod(wi +xmodmul(a[j],ve[j], pr), pr);
      w.push_back(wi);
    }
  return Zvec<T>(w);
}
template Zvec<int> matvecmulmodp<int>(const Zmat<int>&, const Zvec<int>&, const int& pr);
template Zvec<long> matvecmulmodp<long>(const Zmat<long>&, const Zvec<long>&, const long& pr);
template Zvec<ZZ> matvecmulmodp<ZZ>(const Zmat<ZZ>&, const Zvec<ZZ>&, const ZZ& pr);
template Zvec<INT> matvecmulmodp<INT>(const Zmat<INT>&, const Zvec<INT>&, const INT& pr);

template<class T>
int liftmat(const Zmat<T>& mm, const T& pr, Zmat<T>& m, T& dd)
{
  int trace=0;
  if(trace)
    cout << "Lifting mod-p mat;  mat mod "<<pr<<" is:\n"
         << mm << "\n"
         << "Now lifting back to Q." << endl;

  T n,d, lim = isqrt(pr>>1);
  dd=1;
  m = mm;
  m.reduce_mod_p(pr);
  if (maxabs(m) < lim)
    {
      if(trace)
        cout << "Lifted matrix is " << m << "\n with denominator 1" << endl;
      return 1;
    }
  int success = 1;
  std::for_each(m.entries.begin(), m.entries.end(),
                [&success,lim,&dd,pr,&n,&d] (const T& x)
                {
                  if (abs(x)>lim)
                    {
                      int succ = modrat(x,pr,lim,n,d);
                      if(succ) dd=lcm(d,dd);
                      else
                        {
                          cerr<<"Failed to lift "<<x<<" mod "<<pr<< " in liftmat() " << endl;
                          exit(1);
                          success=0;
                        }
                    }
                });
  dd=abs(dd);
  if(trace)
    cout << "Common denominator = " << dd << "\n";
  std::transform(m.entries.begin(), m.entries.end(), m.entries.begin(),
                 [pr,dd] (const T& x) {return mod(xmodmul(dd,x,pr),pr);});
  if (!success)
    {
      cerr<<"Failed to lift some entries mod "<<pr<< " in liftmat() " << endl;
      return 0;
    }
  if(trace)
    cout << "Lifted matrix is " << m << "\n";
  return 1;
}
template int liftmat<int> (const Zmat<int>&,  const int&,  Zmat<int>&,  int&);
template int liftmat<long>(const Zmat<long>&, const long&, Zmat<long>&, long&);
template int liftmat<ZZ>  (const Zmat<ZZ>&,   const ZZ&,   Zmat<ZZ>&,   ZZ&);
template int liftmat<INT> (const Zmat<INT>&,  const INT&,  Zmat<INT>&,  INT&);

template<class T>
T maxabs(const Zmat<T>& m) // max entry
{
  return m.entries.empty()?
    T(0) :
    std::accumulate(m.entries.begin(), m.entries.end(), T(0),
                    [](const T& x, const T& y) {return max(x,abs(y));});
}
template int  maxabs<int> (const Zmat<int>&);
template long maxabs<long>(const Zmat<long>&);
template ZZ   maxabs<ZZ>  (const Zmat<ZZ>&);
template INT  maxabs<INT> (const Zmat<INT>&);

template<class T>
long population(const Zmat<T>& m) // #nonzero entries
{
  if (m.entries.empty()) return 0;
  return std::count_if(m.entries.begin(), m.entries.end(), [](const T& x) {return is_nonzero(x);});
}
template long population<int>(const Zmat<int>&);
template long population<long>(const Zmat<long>&);
template long population<ZZ>(const Zmat<ZZ>&);
template long population<INT>(const Zmat<INT>&);

template<class T>
double sparsity(const Zmat<T>& m)
{
  if (m.entries.empty()) return 1;
  return double(population(m))/m.entries.size();
}
template double sparsity<int>(const Zmat<int>&);
template double sparsity<long>(const Zmat<long>&);
template double sparsity<ZZ>(const Zmat<ZZ>&);
template double sparsity<INT>(const Zmat<INT>&);

template<class T>
int is_permutation_matrix(const Zmat<T>& m)
{
  if (m.entries.empty())
    return 1;
  int nr = m.nro, nc=m.nco;
  if (nr!=nc)
    return 0;

  auto test =  [](const T& x){return is_zero(x) || is_one(x);};
  if (!all_of(m.entries.begin(), m.entries.end(), test))
    return 0;

  auto test0 = [](const T& x){return is_zero(x);};
  for (int i=1; i<=nr; i++)
    {
      const auto& v = m.row(i).get_entries();
      if (std::count_if(v.begin(), v.end(), test0) != nc-1)
        return 0;
      const auto& w = m.col(i).get_entries();
      if (std::count_if(w.begin(), w.end(), test0) != nr-1)
        return 0;
    }
  return 1;
}
template int is_permutation_matrix<int>(const Zmat<int>& m);
template int is_permutation_matrix<long>(const Zmat<long>& m);
template int is_permutation_matrix<ZZ>(const Zmat<ZZ>& m);
template int is_permutation_matrix<INT>(const Zmat<INT>& m);

template<class T>
int is_signed_permutation_matrix(const Zmat<T>& m)
{
  if (m.entries.empty())
    return 1;
  int nr = m.nro, nc=m.nco;
  if (nr!=nc)
    return 0;

  // The only difference between this and the previous function is the
  // abs() here:
  auto test = [](const T& x){return is_zero(x) || is_one(abs(x));};
  if (!all_of(m.entries.begin(), m.entries.end(), test))
    return 0;

  auto test0 = [](const T& x){return is_zero(x);};
  for (int i=1; i<=nr; i++)
    {
      const auto& v = m.row(i).get_entries();
      if (std::count(v.begin(), v.end(), 0) != nc-1)
        return 0;
      const auto& w = m.col(i).get_entries();
      if (std::count(w.begin(), w.end(), 0) != nr-1)
        return 0;
    }
  return 1;
}
template int is_signed_permutation_matrix<int>(const Zmat<int>& m);
template int is_signed_permutation_matrix<long>(const Zmat<long>& m);
template int is_signed_permutation_matrix<ZZ>(const Zmat<ZZ>& m);
template int is_signed_permutation_matrix<INT>(const Zmat<INT>& m);

template<class T>
int is_identity_matrix(const Zmat<T>& m)
{
  return m == Zmat<T>::identity_matrix(m.nro);
}
template int is_identity_matrix<int>(const Zmat<int>& m);
template int is_identity_matrix<long>(const Zmat<long>& m);
template int is_identity_matrix<ZZ>(const Zmat<ZZ>& m);
template int is_identity_matrix<INT>(const Zmat<INT>& m);

template<class T>
mat_m to_mat_m(const Zmat<T>& M)
{
  const vector<T> & Mij = M.get_entries();
  vector<ZZ> entries(Mij.size());
  std::transform(Mij.begin(), Mij.end(), entries.begin(), [](const T& x) {return to_ZZ(x);});
  return mat_m(M.nrows(), M.ncols(), entries);
}
template mat_m to_mat_m(const Zmat<int>& M);
template mat_m to_mat_m(const Zmat<long>& M);
template mat_m to_mat_m(const Zmat<ZZ>& M);
template mat_m to_mat_m(const Zmat<INT>& M);

template<class T>
mat_i to_mat_i(const Zmat<T>& M)
{
  const vector<T> & Mij = M.get_entries();
  vector<int> entries(Mij.size());
  std::transform(Mij.begin(), Mij.end(), entries.begin(),
                 [](const T& i){return to_int(i);});
  return mat_i(M.nrows(), M.ncols(), entries);
}
template mat_i to_mat_i(const Zmat<int>& M);
template mat_i to_mat_i(const Zmat<long>& M);
template mat_i to_mat_i(const Zmat<ZZ>& M);

template<class T>
mat_l to_mat_l(const Zmat<T>& M)
{
  const vector<T> & Mij = M.get_entries();
  vector<long> entries(Mij.size());
  std::transform(Mij.begin(), Mij.end(), entries.begin(),
                 [](const T& i){return to_long(i);});
  return mat_l(M.nrows(), M.ncols(), entries);
}
template mat_l to_mat_l(const Zmat<int>& M);
template mat_l to_mat_l(const Zmat<long>& M);
template mat_l to_mat_l(const Zmat<ZZ>& M);

template<class T>
mat_I to_mat_I(const Zmat<T>& M)
{
  const vector<T> & Mij = M.get_entries();
  vector<INT> entries(Mij.size());
  std::transform(Mij.begin(), Mij.end(), entries.begin(), [](const T& x) {return to_INT(x);});
  return mat_I(M.nrows(), M.ncols(), entries);
}
template mat_I to_mat_I(const Zmat<int>& M);
template mat_I to_mat_I(const Zmat<long>& M);
template mat_I to_mat_I(const Zmat<ZZ>& M);

///////////////////////////////////////////////////////////////////////////

// FLINT has more than one type for modular matrices: standard in
// FLINT-2.3..2.9 was nmod_mat_t with entries of type mp_limb_t
// (unsigned long) while non-standard was hmod_mat_t, with entries
// hlimb_t (unsigned int).  From FLINT-3 the latter is emulated via a
// wrapper.  We use the former when scalar=long and the latter when
// scalar=int and the FLINT version is at least 3.  The unsigned
// scalar types are #define'd as uscalar.

// Implementation of wrapper functions declared in flinterface.h
// written by Fredrik Johansson

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
#if (__FLINT_VERSION==3)&&(__FLINT_VERSION_MINOR<3)
        flint_free(mat->rows);
#endif
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
  slong rk;
  gr_ctx_t ctx;
  gr_ctx_init_nmod32(ctx, mat->mod.n);
  GR_MUST_SUCCEED(gr_mat_rref_lu(&rk, (gr_mat_struct *) mat, (gr_mat_struct *) mat, ctx));
  return rk;
}

// create flint matrix (type hmod_mat_t) copy of a Zmat<int> modulo pr:
void mod_mat_from_mat(hmod_mat_t& A, const Zmat<int>& M, const int& pr)
{
  long nr=M.nrows(), nc=M.ncols();
  hmod_mat_init(A, nr, nc, (hlimb_t)pr);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      hmod_mat_entry(A,i,j) = (hlimb_t)posmod(M(i+1,j+1),pr);
}

// create flint matrix (type nmod_mat_t) copy of a Zmat<long> modulo pr:
void mod_mat_from_mat(nmod_mat_t& A, const Zmat<long>& M, const long& pr)
{
  long nr=M.nrows(), nc=M.ncols();
  nmod_mat_init(A, nr, nc, (mp_limb_t)pr);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      nmod_mat_entry(A,i,j) = (mp_limb_t)posmod(M(i+1,j+1),pr);
}

// create flint matrix (type fmpz_mat_t) copy of a Zmat<T>:
template<class T>
void flint_mat_from_mat(fmpz_mat_t& A, const Zmat<T>& M)
{
  long nr=M.nrows(), nc=M.ncols();
  fmpz_mat_init(A, nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      fmpz_set(fmpz_mat_entry(A,i,j), *to_FLINT(M(i+1,j+1)));
}
template void flint_mat_from_mat<int>(fmpz_mat_t& A, const Zmat<int>& M);
template void flint_mat_from_mat<long>(fmpz_mat_t& A, const Zmat<long>& M);
template void flint_mat_from_mat<ZZ>(fmpz_mat_t& A, const Zmat<ZZ>& M);
template void flint_mat_from_mat<INT>(fmpz_mat_t& A, const Zmat<INT>& M);

// convert a flint matrix (type fmpz_mat_t) to a Zmat<T>.  The dummy
// variable is needed to determine the return type
template<class T>
Zmat<T> mat_from_flint_mat(fmpz_mat_t& A, const T& dummy)
{
  long nr = fmpz_mat_nrows(A), nc = fmpz_mat_ncols(A);
  Zmat<T> M(nr,nc);
  fmpz_t Mij;
  fmpz_init(Mij);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      {
        fmpz_set(Mij, fmpz_mat_entry(A,i,j));
        set(M(i+1,j+1), Mij);
      }
  fmpz_clear(Mij);
  return M;
}
template Zmat<int> mat_from_flint_mat<int>(fmpz_mat_t& A, const int& dummy);
template Zmat<long> mat_from_flint_mat<long>(fmpz_mat_t& A, const long& dummy);
template Zmat<ZZ> mat_from_flint_mat<ZZ>(fmpz_mat_t& A, const ZZ& dummy);
template Zmat<INT> mat_from_flint_mat<INT>(fmpz_mat_t& A, const INT& dummy);

// create flint matrix (type fmpz_mod_mat_t) copy of a Zmat<ZZ> modulo
// pr. The context flint_modulus must have been preset to match pr and
// fmpz_mod_mat_init(A, flint_modulus) should have been called already.  We pass
// pr as a parameter since the FLINT matrix must have integers in
// therange 0..pr-1.
void mod_mat_from_mat(fmpz_mod_mat_t& A, const fmpz_mod_ctx_t& flint_modulus, const Zmat<ZZ>& M, const ZZ& pr)
{
  long nr=M.nrows(), nc=M.ncols();
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      fmpz_mod_mat_set_entry(A,i,j, *to_FLINT(posmod(M(i+1,j+1),pr)), flint_modulus);
}

// create flint matrix (type fmpz_mod_mat_t) copy of a Zmat<ZZ> modulo
// pr. The context flint_modulus must have been preset to match pr and
// fmpz_mod_mat_init(A, flint_modulus) should have been called already.  We pass
// pr as a parameter since the FLINT matrix must have integers in
// therange 0..pr-1.
void mod_mat_from_mat(fmpz_mod_mat_t& A, const fmpz_mod_ctx_t& flint_modulus, const Zmat<INT>& M, const INT& pr)
{
  fmpz_t tmp;
  long nr=M.nrows(), nc=M.ncols();
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      {
        INT Aij = posmod(M(i+1,j+1),pr);
        Aij.get_fmpz(tmp);
        fmpz_mod_mat_set_entry(A,i,j, tmp, flint_modulus);
      }
  fmpz_clear(tmp);
}

Zmat<int> mat_from_mod_mat(const hmod_mat_t& A)
{
  long nr=hmod_mat_nrows(A), nc=hmod_mat_ncols(A);
  Zmat<int> M(nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      M(i+1,j+1) = (int)hmod_mat_entry(A,i,j);
  return M;
}

Zmat<long> mat_from_mod_mat(const nmod_mat_t& A)
{
  long nr=nmod_mat_nrows(A), nc=nmod_mat_ncols(A);
  Zmat<long> M(nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      M(i+1,j+1) = (long)nmod_mat_entry(A,i,j);
  return M;
}

Zmat<ZZ> mat_from_mod_mat(const fmpz_mod_mat_t& A, const fmpz_mod_ctx_t& flint_modulus, const ZZ& dummy)
{
  long nr=fmpz_mod_mat_nrows(A, flint_modulus), nc=fmpz_mod_mat_ncols(A, flint_modulus);
  Zmat<ZZ> M(nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      {
        fmpz_t Aij = {*fmpz_mod_mat_entry(A,i,j)};
        M(i+1,j+1) = FLINT_to_NTL(Aij);
      }
  return M;
}

Zmat<INT> mat_from_mod_mat(const fmpz_mod_mat_t& A, const fmpz_mod_ctx_t& flint_modulus, const INT& dummy)
{
  long nr=fmpz_mod_mat_nrows(A, flint_modulus), nc=fmpz_mod_mat_ncols(A, flint_modulus);
  Zmat<INT> M(nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      {
        fmpz_t Aij = {*fmpz_mod_mat_entry(A,i,j)};
        M(i+1,j+1) = INT(Aij);
      }
  return M;
}

// The following four versions of ref_mod_p() are not templated.  The
// int and long versions use hmod_mat_rref() and nmod_mat_rref()
// respectively; the ZZ and INT versions are almost identical so could
// be templated, both use fmpz_mod_mat_rref().

template<>
Zmat<int> ref_mod_p(const Zmat<int>& M, const int& pr)
{
  hmod_mat_t A;
  mod_mat_from_mat(A,M,pr);
  long rk = hmod_mat_rref(A);
  Zmat<int> B = mat_from_mod_mat(A).slice(rk, M.ncols());
  hmod_mat_clear(A);
  B.reduce_mod_p(pr);
  return B;
}

template<>
Zmat<long> ref_mod_p(const Zmat<long>& M, const long& pr)
{
  nmod_mat_t A;
  mod_mat_from_mat(A,M,pr);
  long rk = nmod_mat_rref(A);
  Zmat<long> B = mat_from_mod_mat(A).slice(rk, M.ncols());
  nmod_mat_clear(A);
  B.reduce_mod_p(pr);
  return B;
}

template<>
Zmat<ZZ> ref_mod_p(const Zmat<ZZ>& M, const ZZ& pr)
{
  fmpz_mod_ctx_t flint_modulus;
  fmpz_mod_ctx_init(flint_modulus, *to_FLINT(pr));

  long nr=M.nrows(), nc=M.ncols(), rk;
  fmpz_mod_mat_t A, R;
  fmpz_mod_mat_init(A, nr, nc, flint_modulus);
  fmpz_mod_mat_init(R, nr, nc, flint_modulus);

  mod_mat_from_mat(A,flint_modulus,M,pr);
  rk = fmpz_mod_mat_rref(R, A, flint_modulus);
  // the pr here is just to disambiguate the function
  Zmat<ZZ> B = mat_from_mod_mat(R, flint_modulus, pr).slice(rk, nc);
  B.reduce_mod_p(pr);

  fmpz_mod_mat_clear(A,flint_modulus);
  fmpz_mod_mat_clear(R,flint_modulus);
  fmpz_mod_ctx_clear(flint_modulus);

  return B;
}

template<>
Zmat<INT> ref_mod_p(const Zmat<INT>& M, const INT& pr)
{
  fmpz_mod_ctx_t flint_modulus;
  fmpz_t tmp;
  pr.get_fmpz(tmp);
  fmpz_mod_ctx_init(flint_modulus, tmp);
  fmpz_clear(tmp);

  long nr=M.nrows(), nc=M.ncols(), rk;
  fmpz_mod_mat_t A, R;
  fmpz_mod_mat_init(A, nr, nc, flint_modulus);
  fmpz_mod_mat_init(R, nr, nc, flint_modulus);

  mod_mat_from_mat(A,flint_modulus,M,pr);
  rk = fmpz_mod_mat_rref(R, A, flint_modulus);
  // the pr here is just to disambiguate the function
  Zmat<INT> B = mat_from_mod_mat(R, flint_modulus, pr).slice(rk, nc);
  B.reduce_mod_p(pr);

  fmpz_mod_mat_clear(A,flint_modulus);
  fmpz_mod_mat_clear(R,flint_modulus);
  fmpz_mod_ctx_clear(flint_modulus);

  return B;
}

template<class T>
Zmat<T> ref_mod_p(const Zmat<T>& M, const T& pr,
                  Zvec<int>& pcols, Zvec<int>& npcols,
                  long& rk, long& ny)
{
  Zmat<T> R = ref_mod_p(M, pr);
  pnpcols(R, pcols, npcols, rk, ny);
  return R;
}
template Zmat<int> ref_mod_p<int>(const Zmat<int>&, const int&, Zvec<int>&, Zvec<int>&, long&, long&);
template Zmat<long> ref_mod_p<long>(const Zmat<long>&, const long&, Zvec<int>&, Zvec<int>&, long&, long&);
template Zmat<ZZ> ref_mod_p<ZZ>(const Zmat<ZZ>&, const ZZ&, Zvec<int>&, Zvec<int>&, long&, long&);
template Zmat<INT> ref_mod_p<INT>(const Zmat<INT>&, const INT&, Zvec<int>&, Zvec<int>&, long&, long&);

template<class T>
Zmat<T> ref(const Zmat<T>& M, T& dd)
{
  int trace=0;
  long nr=M.nrows(), nc=M.ncols(), rk;
  if (trace) cout << "In ref() with M =\n" << M << endl;

  fmpz_mat_t A, R;
  fmpz_mat_init(A, nr, nc);
  fmpz_mat_init(R, nr, nc);
  flint_mat_from_mat(A,M);
  fmpz_t d;
  fmpz_init(d);
  rk = fmpz_mat_rref(R, d, A);
  set(dd,d);
  Zmat<T> B = mat_from_flint_mat(R, dd).slice(rk, nc);
  if (trace) cout << "ref_via_flint() returns d = " << dd << " and REF =\n" << B << endl;

  fmpz_mat_clear(A);
  fmpz_mat_clear(R);
  fmpz_clear(d);

  // FLINT does not guarantee that content(R) and d are coprime or that d>0
  if (dd<0)
    {
      dd = -dd;
      B = -B;
    }

  if (!is_one(dd))
    {
      T g = gcd(dd, B.content());
      if (!is_one(g))
        {
          dd /= g;
          B /= g;
        }
    }
  return B;
}

// From a matrix in REF, extract the pivotal and non-pivotal columns
template<class T>
void pnpcols(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols, long& rk, long& ny)
{
  vector<int> pc, npc;
  int nc = M.ncols(), nr = M.nrows(), i = 1, j;
  for (j=1; j<=nc; j++)
    {
      if ((i>nr) || (is_zero(M(i,j))))
        {
          npc.push_back(j);
        }
      else
        {
          pc.push_back(j);
          i++;
        }
    }
  pcols  = Zvec<int>(pc);
  npcols = Zvec<int>(npc);
  rk = pc.size();
  ny = npc.size();
}
template void pnpcols<int> (const Zmat<int>&,  Zvec<int>&, Zvec<int>&, long&, long&);
template void pnpcols<long>(const Zmat<long>&, Zvec<int>&, Zvec<int>&, long&, long&);
template void pnpcols<ZZ>  (const Zmat<ZZ>&,   Zvec<int>&, Zvec<int>&, long&, long&);
template void pnpcols<INT> (const Zmat<INT>&,  Zvec<int>&, Zvec<int>&, long&, long&);

// Hermite and Smith Normal Forms (via FLINT)
template<class T>
Zmat<T> HNF(const Zmat<T>& M)
{
  fmpz_mat_t A;//  fmpz_mat_init() is done in the next function
  flint_mat_from_mat(A, M);

  fmpz_mat_hnf(A, A);
  Zmat<T> H = mat_from_flint_mat(A, (T)(0));

  fmpz_mat_clear(A);
  return H;
}
template Zmat<int> HNF<int>(const Zmat<int>& M);
template Zmat<long> HNF<long>(const Zmat<long>& M);
template Zmat<ZZ> HNF<ZZ>(const Zmat<ZZ>& M);
template Zmat<INT> HNF<INT>(const Zmat<INT>& M);

template<class T>
Zmat<T> HNF(const Zmat<T>& M, Zmat<T>& U)
{
  fmpz_mat_t A;//  fmpz_mat_init() is done in the next function
  flint_mat_from_mat(A, M);
  fmpz_mat_t UU;
  fmpz_mat_init(UU, M.nrows(), M.nrows());
  fmpz_mat_one(UU);

  fmpz_mat_hnf_transform(A, UU, A);

  Zmat<T> H = mat_from_flint_mat(A, (T)(0));
  U = mat_from_flint_mat(UU, (T)(0));

  fmpz_mat_clear(A);
  fmpz_mat_clear(UU);

  return H;
}
template Zmat<int> HNF<int>(const Zmat<int>& M, Zmat<int>& U);
template Zmat<long> HNF<long>(const Zmat<long>& M, Zmat<long>& U);
template Zmat<ZZ> HNF<ZZ>(const Zmat<ZZ>& M, Zmat<ZZ>& U);
template Zmat<INT> HNF<INT>(const Zmat<INT>& M, Zmat<INT>& U);

template<class T>
Zmat<T> SNF(const Zmat<T>& M)
{
  fmpz_mat_t A;//  fmpz_mat_init() is done in the next function
  flint_mat_from_mat(A, M);

  fmpz_mat_snf(A, A);
  Zmat<T> S = mat_from_flint_mat(A, (T)(0));

  fmpz_mat_clear(A);

  return S;
}
template Zmat<int> SNF<int>(const Zmat<int>& M);
template Zmat<long> SNF<long>(const Zmat<long>& M);
template Zmat<ZZ> SNF<ZZ>(const Zmat<ZZ>& M);
template Zmat<INT> SNF<INT>(const Zmat<INT>& M);

template<class T>
Zmat<T> SNF(const Zmat<T>& M, Zmat<T>& U, Zmat<T>& V)
{
  fmpz_mat_t A;//  fmpz_mat_init() is done in the next function
  flint_mat_from_mat(A, M);
  fmpz_mat_t UU;
  fmpz_mat_init(UU, M.nrows(), M.nrows());
  fmpz_mat_one(UU);
  fmpz_mat_t VV;
  fmpz_mat_init(VV, M.ncols(), M.ncols());
  fmpz_mat_one(VV);

  fmpz_mat_snf_transform(A, UU, VV, A);

  Zmat<T> S = mat_from_flint_mat(A, (T)(0));
  U = mat_from_flint_mat(UU, (T)(0));
  V = mat_from_flint_mat(VV, (T)(0));

  fmpz_mat_clear(A);
  fmpz_mat_clear(UU);
  fmpz_mat_clear(VV);

  return S;
}
template Zmat<int> SNF<int>(const Zmat<int>& M, Zmat<int>& U, Zmat<int>& V);
template Zmat<long> SNF<long>(const Zmat<long>& M, Zmat<long>& U, Zmat<long>& V);
template Zmat<ZZ> SNF<ZZ>(const Zmat<ZZ>& M, Zmat<ZZ>& U, Zmat<ZZ>& V);
template Zmat<INT> SNF<INT>(const Zmat<INT>& M, Zmat<INT>& U, Zmat<INT>& V);

// LLL row-reduction (via FLINT)
template<class T>
Zmat<T> LLL(const Zmat<T>& M)
{
  Zmat<T> U;
  return LLL(M, U); // U is discarded
}
template Zmat<int> LLL<int>(const Zmat<int>& M);
template Zmat<long> LLL<long>(const Zmat<long>& M);
template Zmat<ZZ> LLL<ZZ>(const Zmat<ZZ>& M);
template Zmat<INT> LLL<INT>(const Zmat<INT>& M);

//#define DEBUG_LLL

template<class T>
Zmat<T> LLL(const Zmat<T>& M, Zmat<T>& U) // U*M=L with U unimodular
{
  fmpz_mat_t A;//  fmpz_mat_init() is done in the next function
  flint_mat_from_mat(A, M);
#ifdef DEBUG_LLL
  flint_printf("fmpz_mat A = \n");
  fmpz_mat_print_pretty(A); cout << endl;
#endif
  auto nr = fmpz_mat_nrows(A);
  fmpz_mat_t UU;
  fmpz_mat_init(UU, nr, nr);
  fmpz_mat_one(UU);
#ifdef DEBUG_LLL
  flint_printf("After fmpz_mat_init, UU = \n");
  fmpz_mat_print_pretty(UU); cout << endl;
#endif
  fmpz_lll_t fl;
  fmpz_lll_context_init_default(fl);

  fmpz_lll(A, UU, fl);
#ifdef DEBUG_LLL
  flint_printf("fmpz_lll sets A = \n");
  fmpz_mat_print_pretty(A); cout << endl;
  flint_printf("fmpz_lll sets UU = \n");
  fmpz_mat_print_pretty(UU); cout << endl;
#endif
  Zmat<T> L = mat_from_flint_mat(A, (T)(0));
  U = mat_from_flint_mat(UU, (T)(0));
#ifdef DEBUG_LLL
  cout << " which converts to U = \n" << U << endl;
#endif
  fmpz_mat_clear(A);
  fmpz_mat_clear(UU);

  return L;
}
template Zmat<int> LLL<int>(const Zmat<int>& M, Zmat<int>& U);
template Zmat<long> LLL<long>(const Zmat<long>& M, Zmat<long>& U);
template Zmat<ZZ> LLL<ZZ>(const Zmat<ZZ>& M, Zmat<ZZ>& U);
template Zmat<INT> LLL<INT>(const Zmat<INT>& M, Zmat<INT>& U);

template<class T>
mat_ZZ mat_to_mat_ZZ(Zmat<T> A)
{
  int d = A.nrows();
  // copy into an NTL matrix:
  mat_ZZ ntl_A;
  ntl_A.SetDims(d,d);
  for(int i=1; i<=d; i++)
    for(int j=1; j<=d; j++)
      ntl_A(i,j) = to_ZZ(A(i,j));
  return ntl_A;
}
template mat_ZZ mat_to_mat_ZZ<int>(Zmat<int> A);
template mat_ZZ mat_to_mat_ZZ<long>(Zmat<long> A);
template mat_ZZ mat_to_mat_ZZ<ZZ>(Zmat<ZZ> A);
template mat_ZZ mat_to_mat_ZZ<INT>(Zmat<INT> A);

template<class T>
mat_ZZ_p mat_to_mat_ZZ_p(Zmat<T> A)
{
  int d = A.nrows();

  // copy into an NTL matrix:
  mat_ZZ_p ntl_A;
  ntl_A.SetDims(d,d);
  for(int i=1; i<=d; i++)
    for(int j=1; j<=d; j++)
      ntl_A(i,j) = NTL::conv<ZZ_p>(to_ZZ(A(i,j)));
  return ntl_A;
}

template mat_ZZ_p mat_to_mat_ZZ_p<int>(Zmat<int> A);
template mat_ZZ_p mat_to_mat_ZZ_p<long>(Zmat<long> A);
template mat_ZZ_p mat_to_mat_ZZ_p<ZZ>(Zmat<ZZ> A);
template mat_ZZ_p mat_to_mat_ZZ_p<INT>(Zmat<INT> A);

// compute char poly of A:
ZZX charpoly(const mat_ZZ& A)
{
  ZZX charpol;
  CharPoly(charpol, A);
  return charpol;
}

// compute char poly of A:
template<class T>
ZZX charpoly(const Zmat<T>& A)
{
  return charpoly(mat_to_mat_ZZ(A));
}

template ZZX charpoly<int>(const Zmat<int>& A);
template ZZX charpoly<long>(const Zmat<long>& A);
template ZZX charpoly<ZZ>(const Zmat<ZZ>& A);
template ZZX charpoly<INT>(const Zmat<INT>& A);

ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den, const ZZ& m)
{
  ZZX charpol;
  CharPoly(charpol, A);
  return reduce_poly(scale_poly_down(charpol, den), m);
}

// return A mod m (or just A if m==0)
mat_ZZ reduce_mat(const mat_ZZ& A, const ZZ& m)
{
  if (m==0) return A;
  int nr= A.NumRows(), nc = A.NumCols();
  mat_ZZ B = A;
  for (int i=1; i<=nr; i++)
    for (int j=1; j<=nc; j++)
      B(i,j) = mod(B(i,j), m);
  return B;
}

// evaluate f(A) (assumes f monic)
mat_ZZ evaluate(const ZZX& f, const mat_ZZ& A)
{
  long d = deg(f);
  mat_ZZ fA = A, I;
  ident(I, A.NumRows());
  for(int i=d-1; i>=0; i--)
    {
      fA += coeff(f,i)*I;
      if(i)
        fA *= A;
    }
  return fA;
}

// evaluate f(A) (assumes f monic)
template<class T>
mat_m evaluate(const ZZX& f, const Zmat<T>& A)
{
  long d = deg(f);
  mat_m mA(to_mat_m(A));
  mat_m fA(mA);
  for(int i=d-1; i>=0; i--)
    {
      fA = addscalar(fA, coeff(f,i));
      if(i)
        fA = fA*mA;
    }
  return fA;
}

template mat_m evaluate<int>(const ZZX& f, const Zmat<int>& A);
template mat_m evaluate<long>(const ZZX& f, const Zmat<long>& A);
template mat_m evaluate<ZZ>(const ZZX& f, const Zmat<ZZ>& A);
template mat_m evaluate<INT>(const ZZX& f, const Zmat<INT>& A);

// evaluate f(A) where the coeffs of f are in co
// (co[1]=const term, co[d+1]=coeff of x^d)
template<class T>
Zmat<T> evaluate(const Zvec<T>& co, const Zmat<T>& A)
{
  long d = dim(co)-1;
  Zmat<T> mA = A;
  Zmat<T> fA = co[d+1]*mA;
  for(int i=d; i>0; i--)
    {
      fA = addscalar(fA, co[i]);
      if(i>1)
        fA = fA * mA;
    }
  return fA;
}

template Zmat<int> evaluate<int>(const Zvec<int>& co, const Zmat<int>& A);
template Zmat<long> evaluate<long>(const Zvec<long>& co, const Zmat<long>& A);
template Zmat<ZZ> evaluate<ZZ>(const Zvec<ZZ>& co, const Zmat<ZZ>& A);
template Zmat<INT> evaluate<INT>(const Zvec<INT>& co, const Zmat<INT>& A);

// p should be monic:
mat_m CompanionMatrix(const ZZX& p)
{
  int d = deg(p);
  mat_m A(d,d);
  ZZ one(1);
  for(int i=1; i<d; i++)
    {
      A(i+1,i) = one;
      A(i,d) = -coeff(p, i-1);
    }
  A(d,d) = -coeff(p, d-1);
  return A;
}

int check_involution(const mat_ZZ& A, const ZZ& den, const ZZ& m, int verbose)
{
  int res = IsDiag(reduce_mat(sqr(A), to_ZZ(m)), A.NumRows(), to_ZZ(den*den));
  if (verbose)
    cout << (res? "Involution!": "NOT an involution....") << "\n";
  return res;
}

int commute(const mat_ZZ& A, const mat_ZZ& B, const ZZ& m)
{
  return IsZero(reduce_mat(A*B-B*A, to_ZZ(m)));
}

// check that a matrix commutes with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist, const ZZ& modulus)
{
  return std::all_of(Blist.begin(), Blist.end(), [A, modulus] (const mat_ZZ& B) {return commute(A,B,modulus);});
}

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vec_m& co, const vector<mat_m>& mats)
{
  int n = mats.size(), d = mats[0].nrows();
  mat_m a(d,d);
  for (int i=0; i<n; i++)
    {
      ZZ c = co[i+1];
      if (c!=0)
        a += c*mats[i];
    }
  return a;
}

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vector<ZZ>& co, const vector<mat_m>& mats)
{
  int n = mats.size(), d = mats[0].nrows();
  mat_m a(d,d);
  for (int i=0; i<n; i++)
    {
      ZZ c = co[i];
      if (c!=0)
        a += c*mats[i];
    }
  return a;
}

// rank of an NTL matrix:
long rank(mat_ZZ A)
{
  ZZ d2;
  return NTL::image(d2, A);
}

// nullity of an NTL matrix:
long nullity(mat_ZZ A)
{
  return A.NumRows()-rank(A);
}

///////////////////////////////////////////////////////////////////////////

// Instantiate Zmat template classes for T=int, long, ZZ, INT

template class Zmat<int>;
template class Zmat<long>;
template class Zmat<ZZ>;
template class Zmat<INT>;

// Instantiate inline Zmat template functions for T=int, long, ZZ, INT
template ostream& operator<< <int>(ostream&s, const Zmat<int>&m);
template ostream& operator<< <long>(ostream&s, const Zmat<long>&m);
template ostream& operator<< <ZZ>(ostream&s, const Zmat<ZZ>&m);
template ostream& operator<< <INT>(ostream&s, const Zmat<INT>&m);

#if(0)
//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Interface with NTL matrices
//
//////////////////////////////////////////////////////////////////////////////////////////////

#define TRACE_NTL_REF

#include <NTL/mat_lzz_p.h>
#ifdef TRACE_NTL_REF
#include "eclib/timer.h"
#endif

// Construct an NTL mat_lzz_p (matrix mod p) from a mat mod pr

template<class T>
mat_zz_p mat_zz_p_from_mat(const Zmat<T>& M, const T& pr)
{
  long nr=M.nrows(), nc=M.ncols();
#ifdef TRACE_NTL_REF
  cout<<"Creating an NTL mat_zz_p from a matrix with " << nr <<" rows and "<<nc<<" columns, mod "<<pr<<endl;
#endif
  // create NTL matrix copy of M:
  NTL::zz_pPush push(I2long(pr));
  mat_zz_p A(NTL::INIT_SIZE, nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      A.put(i,j, NTL::conv<zz_p>(M(i+1,j+1)));
#ifdef TRACE_NTL_REF
  cout<<"--done."<<endl;
#endif
  return A;
}

// Construct a mat (T type same as pr) from an NTL mat_lzz_p

template<class T>
Zmat<T> mat_from_mat_zz_p(const mat_zz_p& A, const T& pr) // type of T fixes return type
{
 long nr = A.NumRows(), nc = A.NumCols();
#ifdef TRACE_NTL_REF
  cout<<"Creating a mat from an NTL mat_zz_p with " << nr <<" rows and "<<nc<<" columns, mod "<<pr<<endl;
#endif
 // create matrix copy of A:
 Zmat<T> M(nr, nc);
 for(long i=0; i<nr; i++)
   for(long j=0; j<nc; j++)
     M(i+1,j+1) = mod(NTL::conv<T>(A.get(i,j)), pr);
#ifdef TRACE_NTL_REF
  cout<<"--done."<<endl;
#endif
 return M;
}

// compute ref of M mod pr via NTL, setting rk=rank, ny=nullity,
// pivotal columns pcols, non-pivotal columns npcols

template<class T>
Zmat<T> ref_via_ntl(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                    long& rk, long& ny, const T& pr)
{
 long nc=M.ncols();
 long i, j, k;
#ifdef TRACE_NTL_REF
 timer ntl_timer;
 ntl_timer.start();
#endif
 NTL::zz_pPush push(I2long(pr));
 mat_zz_p A = mat_zz_p_from_mat(M, pr);

#ifdef TRACE_NTL_REF
 cout<<"--calling NTL's gauss()..."<<flush;
#endif
 rk = NTL::gauss(A); // reduce to echelon form in place; rk is the rank
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
 zz_p zero = NTL::conv<zz_p>(0);
 zz_p one = NTL::conv<zz_p>(1);
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
         NTL::inv(inv_piv, piv);
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
 Zmat<T> ans = mat_from_mat_zz_p(A, pr).slice(rk,nc);
#ifdef TRACE_NTL_REF
 ntl_timer.start();
 ntl_timer.show();
 cout<<endl;
#endif
 return ans;
}

template<class T>
long rank_via_ntl(const Zmat<T>& M, const T& pr)
{
#ifdef TRACE_NTL_REF
  cout << "Computing rank mod "<<pr<<" of a matrix of size ("<<M.nrows()<<", "<<M.ncols()<<")..."<<flush;
  timer ntl_timer;
  ntl_timer.start();
#endif
  NTL::zz_pPush push(I2long(pr));
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
T det_via_ntl(const Zmat<T>& M, const T& pr)
{
#ifdef TRACE_NTL_REF
  cout << "Computing determinant mod "<<pr<<" of a matrix of size ("<<M.nrows()<<", "<<M.ncols()<<")..."<<flush;
  timer ntl_timer;
  ntl_timer.start();
#endif
  NTL::zz_pPush push(I2long(pr));
  mat_zz_p A = mat_zz_p_from_mat(M, pr);
  zz_p det = determinant(A);
#ifdef TRACE_NTL_REF
  cout << "done: "<<flush;
  ntl_timer.start();
  ntl_timer.show();
  cout<<endl;
#endif
  return mod(NTL::conv<T>(det), pr);
}
template Zmat<int> ref_via_ntl<int>(const Zmat<int>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const int& pr);
template long rank_via_ntl<int>(const Zmat<int>& M, const int& pr);
template int det_via_ntl<int>(const Zmat<int>& M, const int& pr);
template Zmat<long> ref_via_ntl<long>(const Zmat<long>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                         long& rk, long& ny, const long& pr);
template long rank_via_ntl<long>(const Zmat<long>& M, const long& pr);
template long det_via_ntl<long>(const Zmat<long>& M, const long& pr);
template Zmat<ZZ> ref_via_ntl<ZZ>(const Zmat<ZZ>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                         long& rk, long& ny, const ZZ& pr);
template long rank_via_ntl<ZZ>(const Zmat<ZZ>& M, const ZZ& pr);
template ZZ det_via_ntl<ZZ>(const Zmat<ZZ>& M, const ZZ& pr);

#endif

