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
 
#include "eclib/convert.h"
#include "eclib/linalg.h"

// Instantiate Zmat template classes for T=int, long, bigint

template class Zmat<int>;
template class Zmat<long>;
template class Zmat<bigint>;

// Definitions of member operators and functions:

template<class T>
void Zmat<T>::init(long nr, long nc) // resets to zero mat of given size;
{                                // with defaults (0,0) releases all space.
  nro = nr;
  nco = nc;
  entries.resize(nro*nco, T(0));
}

template<class T>
T& Zmat<T>::operator()(long i, long j)   // returns ref to (i,j) entry
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
  if (is_zero(scal)||is_one(scal)) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [scal](const T& x) {return x / scal;});
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
  if (is_zero(g)||is_one(g)) return;
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
  if (is_zero(scal))
    std::fill(entries.begin(), entries.end(), T(0));
  else
    std::transform(entries.begin(), entries.end(), entries.begin(),
                   [scal](const T& x) {return x * scal;});
}

template<class T>
void Zmat<T>::operator/=(const T& scal)
{
  if (is_zero(scal)||is_one(scal)) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const T& x) {return x / scal;});
}

// Definitions of non-member, friend operators and functions

// add/sub row i of mat to v
template<class T>
void add_row_to_vec(Zvec<T>& v, const Zmat<T>& m, long i)
{
  std::transform(v.entries.begin(), v.entries.end(),
                 m.entries.begin()+(i-1)*m.nco,
                 v.entries.begin(), std::plus<T>());
}

template<class T>
void sub_row_to_vec(Zvec<T>& v, const Zmat<T>& m, long i)
{
  std::transform(v.entries.begin(), v.entries.end(),
                 m.entries.begin()+(i-1)*m.nco,
                 v.entries.begin(), std::minus<T>());
}

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
   }
 return m3;
}

template<class T>
int operator==(const Zmat<T>& m1, const Zmat<T>& m2)
{
  return (m1.nro==m2.nro) && (m1.nco==m2.nco) && (m1.entries==m2.entries);
}

template<class T>
void Zmat<T>::output(ostream& s) const
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
void Zmat<T>::output_pari(ostream& s) const
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
istream& operator>>(istream& s, Zmat<T>& m) // m cannot be const
{
 long n=m.nro*m.nco;
 auto mij=m.entries.begin();
 while(n--) s >> (*mij++);
 return s;
}

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
   cerr << "colcat: matrices have different number of rows!" << endl;
 return c;
}

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
   cerr << "rowcat: matrices have different number of columns!" << endl;
 return c;
}

template<class T>
Zmat<T> directsum(const Zmat<T>& a, const Zmat<T>& b)
{
  return rowcat(colcat(a,Zmat<T>(a.nro,b.nco)),colcat(Zmat<T>(b.nro,a.nco),b));
}

//plain elimination, no clearing
template<class T>
void elimrows(Zmat<T>& m, long r1, long r2, long pos) // m cannot be const
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
void elimrows1(Zmat<T>& m, long r1, long r2, long pos)
{
  elimrows(m,r1,r2,pos);
  m.clearrow(r2);
}

//elimination + divide by last pivot
template<class T>
void elimrows2(Zmat<T>& m, long r1, long r2, long pos, const T& last)
{
  elimrows(m,r1,r2,pos);
  m.divrow(r2,last);
}

// Definition of non-friend functions

template<class T>
Zmat<T> operator+(const Zmat<T>& m)
{
  return m;
}

template<class T>
Zmat<T> operator-(const Zmat<T>& m)
{
  return T(-1)*m;
}

template<class T>
Zmat<T> operator+(const Zmat<T>& m1, const Zmat<T>& m2)
{
  Zmat<T> ans(m1); ans+=m2;  return ans;
}

template<class T>
Zmat<T> operator-(const Zmat<T>& m1, const Zmat<T>& m2) 
{
  Zmat<T> ans(m1); ans-=m2;  return ans;
}

template<class T>
Zmat<T> operator*(const T& scal, const Zmat<T>& m)
{
  Zmat<T> ans(m); ans*=scal;  return ans;
}

template<class T>
Zmat<T> operator/(const Zmat<T>& m, const T& scal)
{
  Zmat<T> ans(m); ans/=scal;  return ans;
}

template<class T>
int operator!=(const Zmat<T>& m1, const Zmat<T>& m2)
{
  return !(m1==m2);
}

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
   cerr << "Incompatible sizes in *(mat,vec)"<<endl;
 return w;
}

template<class T>
Zmat<T> Zmat<T>::scalar_matrix(long n, const T& a)
{
  Zmat<T> D(n,n);
  for (long i=1; i<=n; i++) D.set(i,i,a);
  return D;
}

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

template<class T>
Zmat<T> echelon(const Zmat<T>& entries, Zvec<int>& pcols, Zvec<int>& npcols,
            long& rk, long& ny, T& d, int method)
{
  switch (method)
    {
    case 0: default: return echelon0(entries,pcols,npcols,rk,ny,d);
    case 2: return echelonp(entries,pcols,npcols,rk,ny,d, default_modulus<T>());
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
Zmat<T> echelon0(const Zmat<T>& entries, Zvec<int>& pc, Zvec<int>& npc,
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
  Zmat<T> ans(rk,nc, m);
  return ans;
}

template<class T>
long Zmat<T>::rank() const
{
  long rk=0;
  T lastpivot(1);
  Zmat<T> m(*this); // work with a copy, which will be reduced
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
    }
  return clist;
}

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

template<class T>
void Zvec<T>::add_row(const Zmat<T>& m, int i)
{
  long n=entries.size();
  auto wi = m.entries.begin() + (i-1)*n;
  std::transform(entries.begin(), entries.end(), wi, entries.begin(), std::plus<T>());
}

template<class T>
Zmat<T> addscalar(const Zmat<T>& mm, const T& c)
{
  return mm + Zmat<T>::scalar_matrix(mm.nrows(), c);
}

template<class T>
Zvec<T> apply(const Zmat<T>& m, const Zvec<T>& v)    // same as *(mat, vec)
{
  return m*v;
}

template<class T>
void Zmat<T>::reduce_mod_p(const T& p)
{
  if (p==0) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [p](const T& mij) {return mod(mij,p);});
}

template<class T>
void elimp(Zmat<T>& m, long r1, long r2, long pos, const T& pr)
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
void elimp1(Zmat<T>& m, long r1, long r2, long pos, const T& pr)
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
Zmat<T> echelonp(const Zmat<T>& entries, Zvec<int>& pcols, Zvec<int>& npcols,
             long& rk, long& ny, T& d, const T& pr)
{
#ifdef TRACE
  cout << "In echelonp\n";
#endif /* TRACE */
 long nr=entries.nrows(), nc=entries.ncols();
 Zmat<T> m(nr,nc);
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
 Zmat<T> nmat(rk,nc);
 Zmat<T> dmat(rk,nc);

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
Zmat<T> echmodp(const Zmat<T>& entries, Zvec<int>& pcols, Zvec<int>& npcols, long& rk, long& ny, const T& pr)
{
 // cout << "In echmodp with p="<<pr<<" and matrix " << entries << endl;
 long nr=entries.nrows(), nc=entries.ncols();
 Zmat<T> m(nr,nc);
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
Zmat<T> echmodp_uptri(const Zmat<T>& entries, Zvec<int>& pcols, Zvec<int>& npcols,
                                  long& rk, long& ny, const T& pr)
{
// cout << "In echmodp_uptri with matrix = " << entries;
 long nr=entries.nrows(), nc=entries.ncols();
 Zmat<T> m(nr,nc);
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
mat_zz_p mat_zz_p_from_mat(const Zmat<T>& M, const T& pr)
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
     M(i+1,j+1) = mod(conv<T>(A.get(i,j)), pr);
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
T det_via_ntl(const Zmat<T>& M, const T& pr)
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
                            {return xmod(xmodmul(m1ik,m2kj,pr)+m3ij, pr);});
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
Zvec<T> matvecmulmodp(const Zmat<T>& M, const Zvec<T>& v, const T& pr)
{
  int m=M.nro, n=M.nco;
  if (n!=dim(v))
   {
     cerr << "Incompatible sizes in mat*vec product: "<<m<<"x"<<n<<" and "<<dim(v)<<endl;
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

template<class T>
int liftmat(const Zmat<T>& mm, const T& pr, Zmat<T>& m, T& dd)
{
  int trace=0;
  if(trace)
    cout << "Lifting mod-p mat;  mat mod "<<pr<<" is:\n"
         << mm
         << "Now lifting back to Q." << endl;

  T n,d;
  T lim = sqrt(pr>>1);
  m = mm;
  m.reduce_mod_p(pr);
  if (maxabs(m) < lim) return 1;
  int success = 1;
  dd=1;
  std::for_each(m.entries.begin(), m.entries.end(),
                [&success,lim,&dd,pr,&n,&d] (const T& x)
                {if (abs(x)>lim) {int succ = modrat(x,pr,n,d); if(succ) d=lcm(d,dd); else success=0;}});
  dd=abs(dd);
  if(trace)
    cout << "Common denominator = " << dd << "\n";
  std::transform(m.entries.begin(), m.entries.end(), m.entries.begin(),
                 [pr,dd] (const T& x) {return mod(xmodmul(dd,x,pr),pr);});
  if (!success)
    {
      cerr<<"liftmat() failed to lift some entries mod "<<pr<<endl;
      return 0;
    }
  if(trace)
    cout << "Lifted matrix is " << m << "\n";
  return 1;
}

template<class T>
T maxabs(const Zmat<T>& m) // max entry
{
  T a(0);
  std::for_each(m.entries.begin(), m.entries.end(), [&a](const T& x) {return max(a,abs(x));});
  return a;
}

template<class T>
long population(const Zmat<T>& m) // #nonzero entries
{
  if (m.entries.empty()) return 0;
  return std::count_if(m.entries.begin(), m.entries.end(), [](const T& x) {return is_nonzero(x);});
}

template<class T>
double sparsity(const Zmat<T>& m)
{
  if (m.entries.empty()) return 1;
  return double(population(m))/m.entries.size();
}


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

#include "eclib/flinterface.h"
#include "flint/gr.h"
#include "flint/gr_mat.h"
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>

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
    slong rank;
    gr_ctx_t ctx;
    gr_ctx_init_nmod32(ctx, mat->mod.n);
    GR_MUST_SUCCEED(gr_mat_rref_lu(&rank, (gr_mat_struct *) mat, (gr_mat_struct *) mat, ctx));
    return rank;
}

// create flint matrix (type hmod_mat_t) copy of a Zmat<int>:
void mod_mat_from_mat(hmod_mat_t& A, const Zmat<int>& M, const int& pr)
{
  long nr=M.nrows(), nc=M.ncols();
  hmod_mat_init(A, nr, nc, (hlimb_t)pr);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      hmod_mat_entry(A,i,j) = (hlimb_t)posmod(M(i+1,j+1),pr);
}

// create flint matrix (type nmod_mat_t) copy of a Zmat<long>:
void mod_mat_from_mat(nmod_mat_t& A, const Zmat<long>& M, const long& pr)
{
  long nr=M.nrows(), nc=M.ncols();
  nmod_mat_init(A, nr, nc, (mp_limb_t)pr);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      nmod_mat_entry(A,i,j) = (mp_limb_t)posmod(M(i+1,j+1),pr);
}

// create flint matrix (type fmpz_mod_mat_t) copy of a Zmat<bigint>:
void mod_mat_from_mat(fmpz_mod_mat_t& A, fmpz_mod_ctx_t& mod, const Zmat<bigint>& M, const bigint& pr)
{
  long nr=M.nrows(), nc=M.ncols();
  fmpz_mod_ctx_init(mod, *NTL_to_FLINT(pr));
  fmpz_mod_mat_init(A, nr, nc, mod);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      fmpz_mod_mat_set_entry(A,i,j, *NTL_to_FLINT(posmod(M(i+1,j+1),pr)), mod);
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

Zmat<bigint> mat_from_mod_mat(const fmpz_mod_mat_t& A, const fmpz_mod_ctx_t& mod)
{
  long nr=fmpz_mod_mat_nrows(A, mod), nc=fmpz_mod_mat_ncols(A, mod);
  Zmat<bigint> M(nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      {
        fmpz_t Aij = {*fmpz_mod_mat_entry(A,i,j)};
        M(i+1,j+1) = FLINT_to_NTL(Aij);
      }
  return M;
}

Zmat<int> ref_via_flint(const Zmat<int>& M, const int& pr)
{
  hmod_mat_t A;
  mod_mat_from_mat(A,M,pr);
  long rk = hmod_mat_rref(A);
  Zmat<int> B = mat_from_mod_mat(A).slice(rk, M.ncols());
  hmod_mat_clear(A);
  return B;
}

Zmat<long> ref_via_flint(const Zmat<long>& M, const long& pr)
{
  nmod_mat_t A;
  mod_mat_from_mat(A,M,pr);
  long rk = nmod_mat_rref(A);
  Zmat<long> B = mat_from_mod_mat(A).slice(rk, M.ncols());
  nmod_mat_clear(A);
  return B;
}

Zmat<bigint> ref_via_flint(const Zmat<bigint>& M, const bigint& pr)
{
  long nr=M.nrows(), nc=M.ncols();

  fmpz_mod_ctx_t mod;
  fmpz_mod_ctx_init(mod, *NTL_to_FLINT(pr));

  fmpz_mod_mat_t A, R;
  fmpz_mod_mat_init(A, nr, nc, mod);
  fmpz_mod_mat_init(R, nr, nc, mod);

  mod_mat_from_mat(A,mod,M,pr);

  long rk = fmpz_mod_mat_rref(R, A, mod);
  Zmat<bigint> B = mat_from_mod_mat(R, mod).slice(rk, nc);
  fmpz_mod_mat_clear(A,mod);
  fmpz_mod_mat_clear(R,mod);
  fmpz_mod_ctx_clear(mod);
  return B;
}

// The following function computes the reduced echelon form of M
// modulo the prime pr, using the appropriate rref function from
// FLINT.

template<class T>
Zmat<T> ref_via_flint(const Zmat<T>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                      long& rk, long& ny, const T& pr)
{
  Zmat<T> R = ref_via_flint(M, pr);

  // construct vectors of pivotal and non-pivotal columns
  rk = R.nrows();
  ny = M.ncols()-rk;
  pcols.init(rk);
  npcols.init(ny);
  long i, j, k;
  T zero(0);
  for (i = j = k = 1; i <= rk; i++)
    {
      while (R(i,j) == zero)
        {
          npcols[k] = j;
          k++;
          j++;
        }
      pcols[i] = j;
      j++;
    }
  while (k <= ny)
    {
      npcols[k] = j;
      k++;
      j++;
    }

  return R;
}

///////////////////////////////////////////////////////////////////////////


// Instantiate Zmat template functions for T=int
template void add_row_to_vec<int>(Zvec<int>& v, const Zmat<int>& m, long i);
template void sub_row_to_vec<int>(Zvec<int>& v, const Zmat<int>& m, long i);
template Zmat<int> operator*<int>(const Zmat<int>&, const Zmat<int>&);
template Zvec<int> operator*<int>(const Zmat<int>&, const Zvec<int>&);
template int operator==<int>(const Zmat<int>&, const Zmat<int>&);
template istream& operator>> <int>(istream&s, Zmat<int>&);
template Zmat<int> colcat<int>(const Zmat<int>& a, const Zmat<int>& b);
template Zmat<int> rowcat<int>(const Zmat<int>& a, const Zmat<int>& b);
template Zmat<int> directsum<int>(const Zmat<int>& a, const Zmat<int>& b);
template void elimrows<int>(Zmat<int>& m, long r1, long r2, long pos);
template void elimrows1<int>(Zmat<int>& m, long r1, long r2, long pos);
template void elimrows2<int>(Zmat<int>& m, long r1, long r2, long pos, const int& last);
template Zmat<int> echelon0<int>(const Zmat<int>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                            long& rk, long& ny, int& d);
template void elimp<int>(Zmat<int>& m, long r1, long r2, long pos, const int& pr);
template void elimp1<int>(Zmat<int>& m, long r1, long r2, long pos, const int& pr);
template Zmat<int> echelonp<int>(const Zmat<int>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                            long& rk, long& ny, int& d, const int& pr);
template Zmat<int> echmodp<int>(const Zmat<int>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                           long& rk, long& ny, const int& pr);
template Zmat<int> echmodp_uptri<int>(const Zmat<int>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                                 long& rk, long& ny, const int& pr);
template Zmat<int> ref_via_ntl<int>(const Zmat<int>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const int& pr);
template Zmat<int> ref_via_flint<int>(const Zmat<int>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const int& pr);
template Zmat<int> rref<int>(const Zmat<int>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                        long& rk, long& ny, const int& pr);
template long rank_via_ntl<int>(const Zmat<int>& M, const int& pr);
template int det_via_ntl<int>(const Zmat<int>& M, const int& pr);
template Zmat<int> transpose<int>(const Zmat<int>& m);
template Zmat<int> matmulmodp<int>(const Zmat<int>&, const Zmat<int>&, const int& pr);
template Zvec<int> matvecmulmodp<int>(const Zmat<int>&, const Zvec<int>&, const int& pr);
template long population<int>(const Zmat<int>& m); // #nonzero entries
template int maxabs<int>(const Zmat<int>& m); // max entry
template double sparsity<int>(const Zmat<int>& m); // #nonzero entries/#entries
template ostream& operator<< <int>(ostream&s, const Zmat<int>&m);
template Zmat<int> operator+<int>(const Zmat<int>&);                   // unary
template Zmat<int> operator-<int>(const Zmat<int>&);                   // unary
template Zmat<int> operator+<int>(const Zmat<int>& m1, const Zmat<int>& m2);
template Zmat<int> operator-<int>(const Zmat<int>& m1, const Zmat<int>& m2);
template Zmat<int> operator*<int>(const int& scal, const Zmat<int>& m);
template Zmat<int> operator/<int>(const Zmat<int>& m, const int& scal);
template int operator!=<int>(const Zmat<int>& m1, const Zmat<int>& m2);
template Zmat<int> rowsubmat<int>(const Zmat<int>& m, const Zvec<int>& v);
template Zmat<int> rowsubmat<int>(const Zmat<int>& m, const Zvec<long>& v);
template Zmat<int> submat<int>(const Zmat<int>& m, const Zvec<int>& iv, const Zvec<int>& jv);
template Zmat<int> submat<int>(const Zmat<int>& m, const Zvec<long>& iv, const Zvec<long>& jv);
template Zmat<int> echelon<int>(const Zmat<int>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                          long& rk, long& ny, int& d, int method=0);
template Zmat<int> addscalar<int>(const Zmat<int>&, const int&);
template Zvec<int> apply<int>(const Zmat<int>&, const Zvec<int>&);
template int liftmat<int>(const Zmat<int>& mm, const int& pr, Zmat<int>& m, int& dd);

template void Zvec<int>::sub_row(const Zmat<int>& m, int i);
template void Zvec<int>::add_row(const Zmat<int>& m, int i);

// Instantiate Zmat template functions for T=long
template void add_row_to_vec<long>(Zvec<long>& v, const Zmat<long>& m, long i);
template void sub_row_to_vec<long>(Zvec<long>& v, const Zmat<long>& m, long i);
template Zmat<long> operator*<long>(const Zmat<long>&, const Zmat<long>&);
template Zvec<long> operator*<long>(const Zmat<long>&, const Zvec<long>&);
template int operator==<long>(const Zmat<long>&, const Zmat<long>&);
template istream& operator>> <long>(istream&s, Zmat<long>&);
template Zmat<long> colcat<long>(const Zmat<long>& a, const Zmat<long>& b);
template Zmat<long> rowcat<long>(const Zmat<long>& a, const Zmat<long>& b);
template Zmat<long> directsum<long>(const Zmat<long>& a, const Zmat<long>& b);
template void elimrows<long>(Zmat<long>& m, long r1, long r2, long pos);
template void elimrows1<long>(Zmat<long>& m, long r1, long r2, long pos);
template void elimrows2<long>(Zmat<long>& m, long r1, long r2, long pos, const long& last);
template Zmat<long> echelon0<long>(const Zmat<long>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                      long& rk, long& ny, long& d);
template void elimp<long>(Zmat<long>& m, long r1, long r2, long pos, const long& pr);
template void elimp1<long>(Zmat<long>& m, long r1, long r2, long pos, const long& pr);
template Zmat<long> echelonp<long>(const Zmat<long>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                      long& rk, long& ny, long& d, const long& pr);
template Zmat<long> echmodp<long>(const Zmat<long>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                     long& rk, long& ny, const long& pr);
template Zmat<long> echmodp_uptri<long>(const Zmat<long>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                     long& rk, long& ny, const long& pr);
template Zmat<long> ref_via_ntl<long>(const Zmat<long>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                         long& rk, long& ny, const long& pr);
template Zmat<long> ref_via_flint<long>(const Zmat<long>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                         long& rk, long& ny, const long& pr);
template Zmat<long> rref<long>(const Zmat<long>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const long& pr);
template long rank_via_ntl<long>(const Zmat<long>& M, const long& pr);
template long det_via_ntl<long>(const Zmat<long>& M, const long& pr);
template Zmat<long> transpose<long>(const Zmat<long>& m);
template Zmat<long> matmulmodp<long>(const Zmat<long>&, const Zmat<long>&, const long& pr);
template Zvec<long> matvecmulmodp<long>(const Zmat<long>&, const Zvec<long>&, const long& pr);
template long population<long>(const Zmat<long>& m); // #nonzero entries
template long maxabs<long>(const Zmat<long>& m); // max entry
template double sparsity<long>(const Zmat<long>& m); // #nonzero entries/#entries
template ostream& operator<< <long>(ostream&s, const Zmat<long>&m);
template Zmat<long> operator+<long>(const Zmat<long>&);                   // unary
template Zmat<long> operator-<long>(const Zmat<long>&);                   // unary
template Zmat<long> operator+<long>(const Zmat<long>& m1, const Zmat<long>& m2);
template Zmat<long> operator-<long>(const Zmat<long>& m1, const Zmat<long>& m2);
template Zmat<long> operator*<long>(const long& scal, const Zmat<long>& m);
template Zmat<long> operator/<long>(const Zmat<long>& m, const long& scal);
template int operator!=<long>(const Zmat<long>& m1, const Zmat<long>& m2);
template Zmat<long> rowsubmat<long>(const Zmat<long>& m, const Zvec<int>& v);
template Zmat<long> rowsubmat<long>(const Zmat<long>& m, const Zvec<long>& v);
template Zmat<long> submat<long>(const Zmat<long>& m, const Zvec<int>& iv, const Zvec<int>& jv);
template Zmat<long> submat<long>(const Zmat<long>& m, const Zvec<long>& iv, const Zvec<long>& jv);
template Zmat<long> echelon<long>(const Zmat<long>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                          long& rk, long& ny, long& d, int method=0);
template Zmat<long> addscalar<long>(const Zmat<long>&, const long&);
template Zvec<long> apply<long>(const Zmat<long>&, const Zvec<long>&);
template int liftmat<long>(const Zmat<long>& mm, const long& pr, Zmat<long>& m, long& dd);

template void Zvec<long>::sub_row(const Zmat<long>& m, int i);
template void Zvec<long>::add_row(const Zmat<long>& m, int i);

// Instantiate Zmat template functions for T=bigint
template void add_row_to_vec<bigint>(Zvec<bigint>& v, const Zmat<bigint>& m, long i);
template void sub_row_to_vec<bigint>(Zvec<bigint>& v, const Zmat<bigint>& m, long i);
template Zmat<bigint> operator*<bigint>(const Zmat<bigint>&, const Zmat<bigint>&);
template Zvec<bigint> operator*<bigint>(const Zmat<bigint>&, const Zvec<bigint>&);
template int operator==<bigint>(const Zmat<bigint>&, const Zmat<bigint>&);
template istream& operator>> <bigint>(istream&s, Zmat<bigint>&);
template Zmat<bigint> colcat<bigint>(const Zmat<bigint>& a, const Zmat<bigint>& b);
template Zmat<bigint> rowcat<bigint>(const Zmat<bigint>& a, const Zmat<bigint>& b);
template Zmat<bigint> directsum<bigint>(const Zmat<bigint>& a, const Zmat<bigint>& b);
template void elimrows<bigint>(Zmat<bigint>& m, long r1, long r2, long pos);
template void elimrows1<bigint>(Zmat<bigint>& m, long r1, long r2, long pos);
template void elimrows2<bigint>(Zmat<bigint>& m, long r1, long r2, long pos, const bigint& last);
template Zmat<bigint> echelon0<bigint>(const Zmat<bigint>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                      long& rk, long& ny, bigint& d);
template void elimp<bigint>(Zmat<bigint>& m, long r1, long r2, long pos, const bigint& pr);
template void elimp1<bigint>(Zmat<bigint>& m, long r1, long r2, long pos, const bigint& pr);
template Zmat<bigint> echelonp<bigint>(const Zmat<bigint>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                      long& rk, long& ny, bigint& d, const bigint& pr);
template Zmat<bigint> echmodp<bigint>(const Zmat<bigint>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                     long& rk, long& ny, const bigint& pr);
template Zmat<bigint> echmodp_uptri<bigint>(const Zmat<bigint>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                     long& rk, long& ny, const bigint& pr);
template Zmat<bigint> ref_via_ntl<bigint>(const Zmat<bigint>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                         long& rk, long& ny, const bigint& pr);
template Zmat<bigint> ref_via_flint<bigint>(const Zmat<bigint>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                         long& rk, long& ny, const bigint& pr);
template Zmat<bigint> rref<bigint>(const Zmat<bigint>& M, Zvec<int>& pcols, Zvec<int>& npcols,
                               long& rk, long& ny, const bigint& pr);
template long rank_via_ntl<bigint>(const Zmat<bigint>& M, const bigint& pr);
template bigint det_via_ntl<bigint>(const Zmat<bigint>& M, const bigint& pr);
template Zmat<bigint> transpose<bigint>(const Zmat<bigint>& m);
template Zmat<bigint> matmulmodp<bigint>(const Zmat<bigint>&, const Zmat<bigint>&, const bigint& pr);
template Zvec<bigint> matvecmulmodp<bigint>(const Zmat<bigint>&, const Zvec<bigint>&, const bigint& pr);
template long population<bigint>(const Zmat<bigint>& m); // #nonzero entries
template bigint maxabs<bigint>(const Zmat<bigint>& m); // max entry
template double sparsity<bigint>(const Zmat<bigint>& m); // #nonzero entries/#entries
template ostream& operator<< <bigint>(ostream&s, const Zmat<bigint>&m);
template Zmat<bigint> operator+<bigint>(const Zmat<bigint>&);                   // unary
template Zmat<bigint> operator-<bigint>(const Zmat<bigint>&);                   // unary
template Zmat<bigint> operator+<bigint>(const Zmat<bigint>& m1, const Zmat<bigint>& m2);
template Zmat<bigint> operator-<bigint>(const Zmat<bigint>& m1, const Zmat<bigint>& m2);
template Zmat<bigint> operator*<bigint>(const bigint& scal, const Zmat<bigint>& m);
template Zmat<bigint> operator/<bigint>(const Zmat<bigint>& m, const bigint& scal);
template int operator!=<bigint>(const Zmat<bigint>& m1, const Zmat<bigint>& m2);
template Zmat<bigint> rowsubmat<bigint>(const Zmat<bigint>& m, const Zvec<int>& v);
template Zmat<bigint> rowsubmat<bigint>(const Zmat<bigint>& m, const Zvec<long>& v);
template Zmat<bigint> submat<bigint>(const Zmat<bigint>& m, const Zvec<int>& iv, const Zvec<int>& jv);
template Zmat<bigint> submat<bigint>(const Zmat<bigint>& m, const Zvec<long>& iv, const Zvec<long>& jv);
template Zmat<bigint> echelon<bigint>(const Zmat<bigint>& m, Zvec<int>& pcols, Zvec<int>& npcols,
                          long& rk, long& ny, bigint& d, int method=0);
template Zmat<bigint> addscalar<bigint>(const Zmat<bigint>&, const bigint&);
template Zvec<bigint> apply<bigint>(const Zmat<bigint>&, const Zvec<bigint>&);
template int liftmat<bigint>(const Zmat<bigint>& mm, const bigint& pr, Zmat<bigint>& m, bigint& dd);

template void Zvec<bigint>::sub_row(const Zmat<bigint>& m, int i);
template void Zvec<bigint>::add_row(const Zmat<bigint>& m, int i);
