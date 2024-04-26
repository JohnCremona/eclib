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

// Definitions of member operators and functions:

void mat::init(long nr, long nc) // resets to zero mat of given size;
{                                // with defaults (0,0) releases all space.
  nro = nr;
  nco = nc;
  entries.resize(nro*nco, 0);
}

scalar& mat::operator()(long i, long j)   // returns ref to (i,j) entry
{
  return entries.at((i-1)*nco+(j-1));
}

scalar mat::operator()(long i, long j) const   // returns (i,j) entry
{
  return entries.at((i-1)*nco+(j-1));
}

scalar mat::sub(long i, long j) const
{
  return entries.at((i-1)*nco+(j-1));
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
 long n=r2-r1+1,c=c2-c1+1;
 mat ans(n,c);
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

mat& mat::operator=(const mat& m)
{
 if (this==&m) return *this;
 nro=m.nro;
 nco=m.nco;
 entries = m.entries;
 return *this;
}

void mat::set(long i, long j, scalar x)
{
  entries.at((i-1)*nco+(j-1)) = x;
}

void mat::add(long i, long j, scalar x)
{
  if (x) entries.at((i-1)*nco+(j-1)) += x;
}

void mat::setrow(long i, const vec& v)
{
  std::copy(v.entries.begin(), v.entries.end(), entries.begin() + (i-1)*nco);
}

void mat::setcol(long j, const vec& v)
{
  auto colj = entries.begin()+(j-1);
  for ( const auto vi : v.entries)
    {
      *colj = vi;
      colj += nco;
    }
}

vec mat::row(long i) const
{
 vec mi(nco);
 auto e = entries.begin()+(i-1)*nco;
 std::copy(e, e+nco, mi.entries.begin());
 return mi;
}

vec mat::col(long j) const
{
 vec v(nro);
 auto entriesij = entries.begin()+(j-1);
 for ( auto& vi : v.entries)
   {
     vi = *entriesij;
     entriesij+=nco;
   }
 return v;
}

void mat::swaprows(long r1, long r2)
{
  auto mr1 = entries.begin() + (r1-1)*nco;
  auto mr2 = entries.begin() + (r2-1)*nco;
  std::swap_ranges(mr1, mr1+nco, mr2);
}

void mat::multrow(long r, scalar scal)
{
  if (scal==1) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [scal](const scalar& x) {return x * scal;});
}

void mat::divrow(long r, scalar scal)
{
  if (scal<=1) return;
  auto mij = entries.begin()+(r-1)*nco;
  std::transform(mij, mij+nco, mij, [scal](const scalar& x) {return x / scal;});
}

scalar mat::content() const
{
  return std::accumulate(entries.begin(), entries.end(), 0,
                         [](const scalar& x, const scalar& y) {return gcd(x,y);});
}

scalar mat::row_content(long r) const
{
  auto mij = entries.begin()+(r-1)*nco;
  return std::accumulate(mij, mij+nco, 0,
                         [](const scalar& x, const scalar& y) {return gcd(x,y);});
}

void mat::clearrow(long r)
{
  divrow(r, row_content(r));
}

void mat::makeprimitive()
{
  scalar g = content();
  if (g<=1) return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [g](const scalar& x) {return x / g;});
}

void mat::operator+=(const mat& n)
{
  std::transform(n.entries.begin(), n.entries.end(), entries.begin(), entries.begin(),
                 [](const scalar& x, const scalar& y) { return x + y;});
}

void mat::operator-=(const mat& n)
{
  std::transform(n.entries.begin(), n.entries.end(), entries.begin(), entries.begin(),
                 [](const scalar& x, const scalar& y) { return y - x;});
}

void mat::operator*=(scalar scal)
{
  if (scal==1)
    return;
  if (scal==0)
    std::fill(entries.begin(), entries.end(), 0);
  else
    std::transform(entries.begin(), entries.end(), entries.begin(),
                   [scal](const scalar& x) {return x * scal;});
}

void mat::operator/=(scalar scal)
{
  if (scal<=1)
    return;
  std::transform(entries.begin(), entries.end(), entries.begin(),
                 [scal](const scalar& x) {return x / scal;});
}

// Definitions of non-member, friend operators and functions

// add/sub row i of mat to v
void add_row_to_vec(vec& v, const mat& m, long i)
{
  std::transform(v.entries.begin(), v.entries.end(),
                 m.entries.begin()+(i-1)*m.nco,
                 v.entries.begin(), std::plus<scalar>());
}

void sub_row_to_vec(vec& v, const mat& m, long i)
{
  std::transform(v.entries.begin(), v.entries.end(),
                 m.entries.begin()+(i-1)*m.nco,
                 v.entries.begin(), std::minus<scalar>());
}

mat operator*(const mat& m1, const mat& m2)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 mat m3(m,p);
 if (n==m2.nro)
   {
     auto a=m1.entries.begin();                                     // a points to m1(i,k)
     for (auto c=m3.entries.begin(); c!=m3.entries.end(); c+=p)     // c points to m3(i,_) for 0<=i<m
       {
         for (auto b=m2.entries.begin(); b!=m2.entries.end(); b+=p) // b points to m2(k,_) for 0<=k<n
           { // add m1(i,k)*m2(k,j) to m3(i,j) for 0<=j<p
             scalar m1ik = *a++;
             std::transform(b, b+p, c, c,
                            [m1ik] (const scalar& m2kj, const scalar& m3ij) {return m1ik*m2kj+m3ij;});
           }
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
  return (m1.nro==m2.nro) && (m1.nco==m2.nco) && (m1.entries==m2.entries);
}

void mat::output(ostream& s) const
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

void mat::output_pari(ostream& s) const
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

long ndigits(scalar a)
{
  int digits = 0;
  if (a < 0) digits = 1; // for the '-'
  while (a) { a /= 10; digits++; }
  return digits;
}

void mat::output_pretty(ostream& s) const
{
  // find max ndgits in each column:
  vector<int> colwidths(nco);
  for(long j=0; j<nco; j++)
    {
      auto mij = entries.begin()+j;
      scalar ma=0, mi=0; // max and min for column j
      for(long i=0; i<nro; i++, mij+=nco)
	{
	  if (*mij>ma) ma=*mij;
	  else if (*mij<mi) mi=*mij;
	}
      ma=ndigits(ma);
      mi=ndigits(mi);
      if(mi>ma)ma=mi;
      colwidths[j]=ma;
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

void mat::dump_to_file(string filename) const
{
  ofstream fout(filename.c_str(),ofstream::binary);
  fout.write((char*)&nro,sizeof(nro));
  fout.write((char*)&nco,sizeof(nco));
  fout.write((char*)entries.data(),nro*nco*sizeof(scalar));
  fout.close();
}

void mat::read_from_file(string filename)
{
  ifstream fin(filename.c_str());
  fin.read((char*)&nro,sizeof(nro));
  fin.read((char*)&nco,sizeof(nco));
  entries.resize(nro*nco);
  fin.read((char*)entries.data(),nro*nco*sizeof(scalar));
  fin.close();
}

istream& operator>>(istream& s, mat& m) // m cannot be const
{
 long n=m.nro*m.nco;
 auto mij=m.entries.begin();
 while(n--) s >> (*mij++);
 return s;
}

mat colcat(const mat& a, const mat& b)
{
 long nr = a.nro, nca = a.nco, ncb = b.nco;
 mat c(nr,nca+ncb);
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

mat rowcat(const mat& a, const mat& b)
{
 mat c(a.nro+b.nro,a.nco);
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

mat directsum(const mat& a, const mat& b)
{
  return rowcat(colcat(a,mat(a.nro,b.nco)),colcat(mat(b.nro,a.nco),b));
}

//plain elimination, no clearing
void elimrows(mat& m, long r1, long r2, long pos) // m cannot be const
{
  long nc=m.nco;
  scalar p = m(r1,pos), q=m(r2,pos);
  auto mr1 = m.entries.begin() + (r1-1)*nc;
  auto mr2 = m.entries.begin() + (r2-1)*nc;
  // replace row2 by p*row2-q*row1
  std::transform(mr1, mr1+nc, mr2, mr2,
                 [p,q] (const scalar& x, const scalar& y) {return p*y-q*x;});
}

//elimination + clearing (i.e. divide new row by its content)
void elimrows1(mat& m, long r1, long r2, long pos)
{
  elimrows(m,r1,r2,pos);
  m.clearrow(r2);
}

//elimination + divide by last pivot
void elimrows2(mat& m, long r1, long r2, long pos, scalar last)
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
 long c=m.nco;
 vec w(m.nro);
 if (c==dim(v))
   {
     auto mi = m.entries.begin();
     for (auto& wi : w.entries)
       {
         wi = std::inner_product(mi, mi+c, v.entries.begin(), 0);
         mi += c;
       }
   }
 else
   cerr << "Incompatible sizes in *(mat,vec)"<<endl;
 return w;
}

mat idmat(scalar n)
{
  mat ans(n,n);
  for (long i=1; i<=n; i++) ans.set(i,i,1);
  return ans;
}

mat transpose(const mat& m)
{
  long nr=m.ncols(), nc=m.nrows();
  mat ans(nr, nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j,  m(j,i));
  return ans;
}

// submatrix of rows indexed by v, all columns
mat rowsubmat(const mat& m, const vec& v)
{
  long nr = dim(v), nc = m.ncols();
  mat ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(v[i],j));
  return ans;
}

// submatrix of rows indexed by iv, columns indexed by jv
mat submat(const mat& m, const vec& iv, const vec& jv)
{
  long nr = dim(iv), nc = dim(jv);
  mat ans(nr,nc);
  for (long i=1; i<=nr; i++)
    for (long j=1; j<=nc; j++)
      ans.set(i,j, m(iv[i],jv[j]));
  return ans;
}

mat echelon(const mat& entries, vec& pcols, vec& npcols,
            long& rk, long& ny, scalar& d, int method)
{
  switch (method)
    {
    case 0: default: return echelon0(entries,pcols,npcols,rk,ny,d);
    case 2: return echelonp(entries,pcols,npcols,rk,ny,d,DEFAULT_MODULUS);
    }
}

//#define DEBUG_ECH_0

//N.B. if(q==0) the following multiplies row r2 by p, which looks
//redundant.  However, it is important to keep this in as in echelon0
//we must guarentee divisibility by "lastpivot".  We do not want to keep
//computing contents of rows as this is slower.
// Used in forward elimination in echelon0

void conservative_elim(vector<scalar>& m, long nc, long r1, long r2, long pos)
{
  auto mr1=m.begin() + r1*nc + pos;
  auto mr2=m.begin() + r2*nc + pos;
  scalar p = *mr1, q = *mr2;
  nc -= pos;
#ifdef DEBUG_ECH_0
  cout<<"In conservative_elim with p = "<<p<<" and q = " << q << endl;
  cout<<"row 1: "; for(long n=0; n<nc; n++) cout<<*(mr1+n)<<",";  cout<<endl;
  cout<<"row 2: "; for(long n=0; n<nc; n++) cout<<*(mr2+n)<<",";  cout<<endl;
#endif
  if ((p==1)&&(q==0))
    return;
  // generic function to make y (entry in row2) 0
  std::function<scalar (scalar, scalar)>
    f = [p,q](const scalar& x, const scalar& y) {return p*y - q*x;};
  if(p==1) // now q!=0
    {
      if(q==1)
        f = [p,q](const scalar& x, const scalar& y) {return y - x;};
      else
        {
          if(q==-1)
            f = [p,q](const scalar& x, const scalar& y) {return y + x;};
          else
            f = [p,q](const scalar& x, const scalar& y) {return y - q*x;};
        }
    }
  else  // p!=1
    {
      if(q==0)
        f = [p,q](const scalar& x, const scalar& y) {return p*y;};
      if(q==1)
        f = [p,q](const scalar& x, const scalar& y) {return p*y - x;};
      if(q==-1)
        f = [p,q](const scalar& x, const scalar& y) {return p*y + x;};
    }
  std::transform(mr1, mr1+nc, mr2, mr2, f);
}

// This version does not multiply row r1 by p unnecessarily.  Used in
// back substitution, it does not assume that the entries in
// columns<pos are 0.

void elim(vector<scalar>& m, long nc, long r1, long r2, long pos)
{
  auto mr1=m.begin()+r1*nc;
  auto mr2=m.begin()+r2*nc;
  scalar p = *(mr1+pos), q = *(mr2+pos);
#ifdef DEBUG_ECH_0
  cout<<"In elim with p = "<<p<<" and q = " << q << endl;
  cout<<"row 1: "; for(long n=0; n<nc; n++) cout<<*(mr1+n)<<",";  cout<<endl;
  cout<<"row 2: "; for(long n=0; n<nc; n++) cout<<*(mr2+n)<<",";  cout<<endl;
#endif
  if ((p==1)&&(q==0))
    return;
  // generic function to make y (entry in row2) 0
  std::function<scalar (scalar, scalar)>
    f = [p,q](const scalar& x, const scalar& y) {return p*y - q*x;};
  if(p==1) // now q!=0
    {
      if(q==1)
        f = [p,q](const scalar& x, const scalar& y) {return y - x;};
      else
        {
          if(q==-1)
            f = [p,q](const scalar& x, const scalar& y) {return y + x;};
          else
            f = [p,q](const scalar& x, const scalar& y) {return y - q*x;};
        }
    }
  else  // p!=1
    {
      if(q==1)
        f = [p,q](const scalar& x, const scalar& y) {return p*y - x;};
      if(q==-1)
        f = [p,q](const scalar& x, const scalar& y) {return p*y + x;};
    }
  std::transform(mr1, mr1+nc, mr2, mr2, f);
}

void clear(vector<scalar>& row, long col1, long col2)
{
  auto row1=row.begin()+col1;
  auto row2=row.begin()+col2;
  scalar g = std::accumulate(row1, row2, 0,
                             [](const scalar& x, const scalar& y) {return gcd(x,y);});
  if (g>1)
    std::for_each(row1, row2, [g](scalar& x) {x/=g;});
}

//#ifndef DEBUG_ECH_0
//#define DEBUG_ECH_0
//#endif

#ifdef DEBUG_ECH_0
void show(vector<scalar> m, long nr, long nc)
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

mat echelon0(const mat& entries, vec& pc, vec& npc,
             long& rk, long& ny, scalar& d)
{
#ifdef DEBUG_ECH_0
  cout<<"In echelon0 with matrix:\n"<<entries<<endl;
#endif
  rk=0; ny=0;
  scalar lastpivot=1;
  long r=0, nc=entries.nco, nr=entries.nro;
  vector<scalar> m = entries.entries;
  vector<scalar> pcols(nc), npcols(nc);
  for (long c=0; (c<nc)&&(r<nr); c++)
    {
      auto mij=m.begin()+r*nc+c;  // points to column c in row r
      scalar piv = abs(*mij);
      long rmin = r;
      mij+=nc;
      for (long r2=r+1; (r2<nr)&&(piv!=1); r2++, mij+=nc)
       {
         scalar mr2c = abs(*mij);
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
                std::transform(mi1, mi1+nc, mi1, [lastpivot]( const scalar& x) {return x/lastpivot;});
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
          scalar fac = d/mij[pcols[r1]];
          std::transform(mij, mij+nc, mij, [fac](const scalar& x){return fac*x;});
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
  mat ans(rk,nc, m);
  return ans;
}

long mat::rank() const
{
  long rk=0, lastpivot=1;
  mat m(*this); // work with a copy, which will be reduced
  long nc=m.ncols(), nr=m.nrows();
  for (long c=1, r=1; (c<=nc)&&(r<=nr); c++)
    {
      long mmin = abs(m(r,c));
      long rmin = r;
      for (long r2=r+1; (r2<=nr)&&(mmin!=1); r2++)
        {
          long mr2c = abs(m(r2,c));
          if ((0<mr2c) && ((mr2c<mmin) || (mmin==0)))
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

long mat::nullity() const
{
 return nco-rank();
}

long mat::trace() const
{
  long tr=0;
  for (long i=0; i<nro; i++)
    tr += entries.at(i*(nco+1));
  return tr;
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
 return (nro%2? -det :det);
}

void vec::sub_row(const mat& m, int i)
{
  long n=entries.size();
  auto wi = m.entries.begin() + (i-1)*n;
  std::transform(entries.begin(), entries.end(), wi, entries.begin(), std::minus<scalar>());
}

void vec::add_row(const mat& m, int i)
{
  long n=entries.size();
  auto wi = m.entries.begin() + (i-1)*n;
  std::transform(entries.begin(), entries.end(), wi, entries.begin(), std::plus<scalar>());
}

mat addscalar(const mat& mm, scalar c)
{
  mat m(mm);
  for (long i=1; i<=m.nrows(); i++)
    m.add(i,i,c);
  return m;
}

vec apply(const mat& m, const vec& v)    // same as *(mat, vec)
{
  return m*v;
}

mat reduce_modp(const mat& m, const scalar& p)
{
  if (p==0) return m;
  mat w(m.nrows(), m.ncols());
  std::transform(m.entries.begin(), m.entries.end(), w.entries.begin(),
                 [p](const scalar& mij) {return mod(mij,p);});
  return w;
}

void elimp(mat& m, long r1, long r2, long pos, scalar pr)
{
  long nc=m.nco;
  auto mr1 = m.entries.begin() + (r1-1)*nc + (pos-1);
  auto mr2 = m.entries.begin() + (r2-1)*nc + (pos-1);
  scalar p = mod(*mr1,pr), q=mod(*mr2,pr);
  if(q==0) {return;} // nothing to do
  nc -= (pos-1); // first pos-1 entries are assumed 0 already
  // generic function to make y (entry in row2) 0
  std::function<scalar (scalar, scalar)>
    f = [pr,p,q](const scalar& x, const scalar& y) {return mod(xmodmul(p,y,pr)-xmodmul(q,x,pr), pr);};
  // simpler special cases (for same signature they must also capture both p and q)
  if(p==1)
   {
     if(q==1)
       f = [pr,p,q](const scalar& x, const scalar& y) {return mod(y-x, pr);};
     else
       {
         if(q==-1)
           f = [pr,p,q](const scalar& x, const scalar& y) {return mod(y+x, pr);};
         else
           // general q
           f = [pr,p,q](const scalar& x, const scalar& y) {return mod(y-xmodmul(q,x,pr), pr);};
       }
   }
  else // general p!=1
    {
      if(q==1)
        f = [pr,p,q](const scalar& x, const scalar& y) {return mod(xmodmul(p,y,pr)-x, pr);};
      if(q==-1)
        f = [pr,p,q](const scalar& x, const scalar& y) {return mod(xmodmul(p,y,pr)+x, pr);};
      // else the generic f will be used
    }
  std::transform(mr1, mr1+nc, mr2, mr2, f);
}

void elimp1(mat& m, long r1, long r2, long pos, scalar pr)
//same as elimp except assumes pivot is 1
{
  long nc=m.nco;
  auto mr1 = m.entries.begin() + (r1-1)*nc + (pos-1);
  auto mr2 = m.entries.begin() + (r2-1)*nc + (pos-1);
  scalar q=mod(*mr2,pr);
  if(q==0) return;
  nc -= (pos-1); // first pos-1 entries are assumed 0 already
  // generic function to make y (entry in row2) 0
  std::function<scalar (scalar, scalar)>
    f = [pr,q](const scalar& x, const scalar& y) {return mod(y-xmodmul(q,x,pr), pr);};
  // simpler special cases
  if (q==1)
    f = [pr,q](const scalar& x, const scalar& y) {return mod(y-x, pr);};
  if (q==-1)
    f = [pr,q](const scalar& x, const scalar& y) {return mod(y+x, pr);};
  std::transform(mr1, mr1+nc, mr2, mr2, f);
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
 long nr=entries.nrows(), nc=entries.ncols();
 mat m(nr,nc);
 std::transform(entries.entries.begin(), entries.entries.end(), m.entries.begin(),
                [pr] (const scalar& x) {return mod(x,pr);});
 pcols.init(nc);
 npcols.init(nc);
 rk=0; ny=0;
 long r=1;
 for (long c=1; (c<=nc)&&(r<=nr); c++)
 {
   scalar mmin = m(r,c);
   long rmin = r;
   for (long r2=r+1; (r2<=nr)&&(mmin==0); r2++)
   {
     scalar mr2c = m(r2,c);
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
       scalar fac = xmod(invmod(m(r1,pcols[r1]),pr),pr);
       for (long c=1; c<=nc; c++)
         m(r1,c)=xmodmul(fac,m(r1,c),pr);
     }
 }
 else
   for (long i=1; i<=rk; i++)
     for (long j=1; j<=nc; j++)
       m(i,j)=(j==pcols[i]);    // 0 or 1 !

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
  // cout << "In echmodp with p="<<pr<<" and matrix " << entries << endl;
 long nr=entries.nrows(), nc=entries.ncols();
 mat m(nr,nc);
 std::transform(entries.entries.begin(), entries.entries.end(), m.entries.begin(),
                [pr] (const scalar& x) {return mod(x,pr);});
 // cout << " - after reducing modulo p,  matrix is " << m << endl;
 pcols.init(nc);
 npcols.init(nc);
 rk=ny=0;
 long r=1;
 for (long c=1; (c<=nc)&&(r<=nr); c++)
   {
     auto mij=m.entries.begin()+(r-1)*nc+c-1;
     scalar mmin = *mij;
     long rmin = r;
     mij += nc;
     for (long r2=r+1; (r2<=nr)&&(mmin==0); r2++, mij+=nc)
       {
	 scalar mr2c = *mij;
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
         // cout<<"pivot = "<<mmin<<endl;
         scalar fac = xmod(invmod(mmin,pr),pr);
         std::transform(entriesij, entriesij+nc, entriesij,
                        [pr,fac] (const scalar& x) {return mod(fac*x, pr);});
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
	 scalar fac = *(mij+pcols[r1]-1);
	 fac = mod(invmod(fac,pr),pr);
         std::transform(mij, mij+nc, mij,
                        [pr,fac] (const scalar& x) {return mod(fac*x, pr);});
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

mat echmodp_uptri(const mat& entries, vec& pcols, vec& npcols,
                                  long& rk, long& ny, scalar pr)
{
// cout << "In echmodp_uptri with matrix = " << entries;
 long nr=entries.nrows(), nc=entries.ncols();
 mat m(nr,nc);
 std::transform(entries.entries.begin(), entries.entries.end(), m.entries.begin(),
                [pr] (const scalar& x) {return mod(x,pr);});
 pcols.init(nc);
 npcols.init(nc);
 rk=ny=0;
 long r=1;
 for (long c=1; (c<=nc)&&(r<=nr); c++)
   {
     auto mij=m.entries.begin()+(r-1)*nc+c-1;
     scalar mmin = *mij;
     long rmin = r;
     mij += nc;
     for (long r2=r+1; (r2<=nr)&&(mmin==0); r2++, mij+=nc)
       {
	 scalar mr2c = *mij;
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
         scalar fac = mod(invmod(mmin,pr),pr);
         std::transform(entriesij, entriesij+nc, entriesij,
                        [pr,fac] (const scalar& x) {return mod(fac*x, pr);});
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

mat_zz_p mat_zz_p_from_mat(const mat& M, scalar pr)
{
  long nr=M.nrows(), nc=M.ncols();
#ifdef TRACE_NTL_REF
  cout<<"Creating an NTL mat_zz_p from a matrix with " << nr <<" rows and "<<nc<<" columns, mod "<<pr<<endl;
#endif
  // create NTL matrix copy of M:
  zz_pPush push(pr);
  mat_zz_p A(INIT_SIZE, nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      A.put(i,j, conv<zz_p>(M(i+1,j+1)));
#ifdef TRACE_NTL_REF
  cout<<"--done."<<endl;
#endif
  return A;
}

// Construct a mat (scalar type same as pr) from an NTL mat_lzz_p

mat mat_from_mat_zz_p(const mat_zz_p& A, scalar pr) // type of scalar fixes return type
{
 long nr = A.NumRows(), nc = A.NumCols();

#ifdef TRACE_NTL_REF
  cout<<"Creating a mat from an NTL mat_zz_p with " << nr <<" rows and "<<nc<<" columns, mod "<<pr<<endl;
#endif
 // create matrix copy of A:
 mat M(nr, nc);
 for(long i=0; i<nr; i++)
   for(long j=0; j<nc; j++)
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

//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Interface with FLINT matrices
//
//////////////////////////////////////////////////////////////////////////////////////////////

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

  // copy of the modulus for FLINT
  uscalar mod = (uscalar)pr;

  // create flint matrix copy of M:
  mod_mat_init(A, nr, nc, mod);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
      mod_mat_entry(A,i,j) = (uscalar)posmod(M(i+1,j+1),pr);
}

mat mat_from_mod_mat(const mod_mat& A, scalar a) // scalar just to fix return type
{
  long nr=mod_mat_nrows(A), nc=mod_mat_ncols(A);

  // create matrix copy of A:
  mat M(nr, nc);
  for(long i=0; i<nr; i++)
    for(long j=0; j<nc; j++)
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

//////////////////////////////////////////////////////////////////////////////////////////////

mat matmulmodp(const mat& m1, const mat& m2, scalar pr)
{
 long m=m1.nro, n=m1.nco, p=m2.nco;
 mat m3(m,p);
 if (n==m2.nro)
   {
     auto a=m1.entries.begin();                                     // a points to m1(i,k)
     for (auto c=m3.entries.begin(); c!=m3.entries.end(); c+=p)     // c points to m3(i,_) for 0<=i<m
       {
         for (auto b=m2.entries.begin(); b!=m2.entries.end(); b+=p) // b points to m2(k,_) for 0<=k<n
           { // add m1(i,k)*m2(k,j) to m3(i,j) for 0<=j<p
             scalar m1ik = *a++;
             std::transform(b, b+p, c, c,
                            [pr,m1ik] (const scalar& m2kj, const scalar& m3ij)
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

int liftmat(const mat& mm, scalar pr, mat& m, scalar& dd, int trace)
{
  float lim=floor(sqrt(pr/2.0));
  if(trace)
    cout << "Lifting mod-p mat;  mat mod "<<pr<<" is:\n"
         << mm
         << "Now lifting back to Q.\n"
         << "lim = " << lim << endl;

  int success=1;
  scalar n,d;
  m = mm;
  dd=1;
  std::for_each(m.entries.begin(), m.entries.end(),
                [&success,&dd,pr,lim,&n,&d] (const scalar& x)
                {success=success&&modrat(x,pr,lim,n,d); d=lcm(d,dd);});
  if(!success)
    return 0;

  dd=abs(dd);
  if(trace)
    cout << "Common denominator = " << dd << "\n";
  std::transform(m.entries.begin(), m.entries.end(), m.entries.begin(),
                 [pr,dd] (const scalar& x) {return mod(xmodmul(dd,x,pr),pr);});
  return 1;
}

long population(const mat& m) // #nonzero entries
{
  if (m.entries.empty()) return 0;
  return std::count_if(m.entries.begin(), m.entries.end(), [](const scalar& x) {return x!=0;});
}

double sparsity(const mat& m)
{
  if (m.entries.empty()) return 1;
  return double(population(m))/m.entries.size();
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
