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
 
// ONLY to be included by smatrix.cc

// Definitions of member operators and functions:

smat::smat(const mat& m)
  :nco(m.nco), nro(m.nro)
{
  rows.resize(nro+1);
  for( int i = 1; i <= nro; i++ )  rows[i]=svec(m.row(i));
}

// member functions and operators

smat& smat::operator= (const smat& m)  // assignment with copy
{
  if (this==&m) return *this;
  nco = m.nco;
  nro = m.nro;
  rows=m.rows;
  return *this;
}

mat smat::as_mat( ) const
{
  mat ans( nro, nco ); 
  for( int i = 1; i <= nro; i++ )
    ans.setrow(i, rows[i].as_vec());
  return ans;
}

smat smat::select_rows(const vec& v) const
  {
    int i, n=dim(v);
    smat ans(n,nco);
    for(i=1; i<=n; i++)
      ans.rows[i]=rows[v[i]];
    return ans;
  }

smat& smat::operator+=(const smat& mat2)
{
  if ((nro!=mat2.nro)||(nco!=mat2.nco))
    {
      cout << "Incompatible smatrices in operator +=\n";
      cout << "Dimensions "<<dim(*this)<<" and "<<dim(mat2)<<endl;
      abort();
    }
  else
    for(int i = 1; i <= nro; i++ ) rows[i]+=mat2.rows[i];
  return *this;
}

smat& smat::operator-=(const smat& mat2)
{
  if ((nro!=mat2.nro)||(nco!=mat2.nco))
    {
      cout << "Incompatible smatrices in operator +=\n";
      cout << "Dimensions "<<dim(*this)<<" and "<<dim(mat2)<<endl;
      abort();
    }
  else
    for(int i = 1; i <= nro; i++ ) rows[i]-=mat2.rows[i];
  return *this;
}

smat& smat::operator+= (const scalar& scal) // adds scalar*identity
{
  if(scal!=0) for(int i = 1; i <= nro; i++ ) rows[i].add(i,scal);
  return *this;
}

smat& smat::operator-= (const scalar& scal) // subtracts scalar*identity
{
  if(scal!=0)  for(int i = 1; i <= nro; i++ ) rows[i].sub(i,scal);
  return *this;
}

void smat::sub_mod_p(const scalar& lambda) // subtracts
					   // scalar*identity mod p
{
  if(lambda!=0)  
    for(int i = 1; i <= nro; i++ ) rows[i].sub_mod_p(i,lambda);
}

smat& smat::operator*=(scalar scal)
{
  //  if(scal==0) cout<<"Attempt to multiply smat by 0\n"<<endl;
  for(int i = 1; i <= nro; i++)  rows[i]*=scal;
  return *this;
}

void smat::reduce_mod_p(const scalar& p)
{
  for(int i = 1; i <= nro; i++)  rows[i].reduce_mod_p(p);
}

smat& smat::mult_by_scalar_mod_p (scalar scal, const scalar& p)
{
  //  if(xmod(scal,p)==0) cout<<"Attempt to multiply smat by 0\n"<<endl;
  for(int i = 1; i <- nro; i++)
    rows[i].mult_by_scalar_mod_p(scal);
  return *this;
}

smat& smat::operator/=(scalar scal)
{
  if(scal==0) {cout<<"Attempt to divide smat by 0\n"<<endl;abort();}
  for(int i = 1; i <= nro; i++)  rows[i]/=scal;
  return *this;
}

smat transpose(const smat& A)
{
  //  cout<<"Transposing an smat of size "<<A.nro<<"x"<<A.nco<<"..."<<flush;
  smat B(A.nco, A.nro);
  for(int j = 1; j <= A.nro; j++)
    for(map<int,scalar>::const_iterator aij = A.rows[j].begin(); 
	aij!=A.rows[j].end(); 
	aij++)
      B.rows[(aij->first)].entries[j] = aij->second;
  //  cout<<"done"<<endl;
  return B;
}

mat smat::operator*( const mat& m )
{
  if( nco != nrows(m) )
    {
      cout << "incompatible smat & mat in operator*\n";
      cout << "Dimensions "<<dim(*this)<<" and ["<<nrows(m)<<" "<<ncols(m)<<"]"<<endl;
      abort();
    }
  int nco = m.nco;
  mat prod( nro, nco );
  for(int i = 1; i <= nro; i++)
    for(int j = 1; j <= nco; j++)
      prod.set(i,j,rows[i]*m.col(j));
  return prod;
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
      s = (A.rows[j])*v;
      if(s) prod.entries[j]=s;
    }
  return prod;
}


smat operator* ( const smat& A, const smat& B )
{
  if( A.nco != B.nro ) 
    { 
      cout << "incompatible smats in operator *\n"; 
      cout << "Dimensions "<<dim(A)<<" and "<<dim(B)<<endl;
      abort();
    }
  int nro = A.nro, nco = B.nco, j;
  smat prod( nro, nco );
  smat tB = transpose(B);
  vector<svec>::const_iterator Ai;
  vector<svec>::iterator Bj, Pi;
  scalar Pij;
  for((Ai = A.rows.begin())++, (Pi = prod.rows.begin())++; Ai != A.rows.end(); Ai++, Pi++)
    {
      if(Ai->size()==0) continue;
      for((Bj = tB.rows.begin())++, j=1; Bj != tB.rows.end(); Bj++, j++)
	if(Bj->size()!=0) 
	  {
	    Pij = (*Ai)*(*Bj);
	    if(Pij) Pi->entries[j]=Pij;
	  }
    }  
  return prod;
}

//#if(0)
smat mult_mod_p ( const smat& A, const smat& B, const scalar& p )
{
  //  cout<<"Multiplying mod-p smats of size "<<A.nro<<"x"<<A.nco<<" and "<<B.nro<<"x"<<B.nco<<"..."<<flush;
  if( A.nco != B.nro ) 
    { 
      cout << "incompatible smats in mult_mod_p(smat,smat,p)\n"; 
      cout << "Dimensions "<<dim(A)<<" and "<<dim(B)<<endl;
      abort();
    }
  int nro = A.nro, nco = B.nco, j;
  smat prod( nro, nco );
  smat tB = transpose(B);
  vector<std::set<int> > Bsups = row_supports(tB);
  vector<svec>::const_iterator Ai;
  vector<svec>::iterator Bj, Pi;
  scalar Pij;
  std::set<int> ABsup;
  std::set<int>::const_iterator ABsupi;
  for((Ai = A.rows.begin())++, (Pi = prod.rows.begin())++; Ai != A.rows.end(); Ai++, Pi++)
    {
      if(Ai->size()==0) continue;
      std::set<int> Asup = Ai->support();
      for((Bj = tB.rows.begin())++, j=1; Bj != tB.rows.end(); Bj++, j++)
	if(Bj->size()!=0) 
	  {
	    ABsup.clear();;
	    set_intersection(Asup.begin(),Asup.end(),Bsups[j].begin(),Bsups[j].end(),inserter(ABsup,ABsup.begin()));
	    //	    cout<<"Intersection of "<<Asup<<" and "<<Bsups[j]<<" is "<<ABsup<<endl;
	    if(ABsup.size()>0)  
	      {
		// Pij = dotmodp(*Ai,*Bj,p);
		Pij=0;
		for(ABsupi=ABsup.begin(); ABsupi!=ABsup.end(); ABsupi++)
		  {
		    Pij=xmod(Pij+xmodmul(Ai->elem(*ABsupi),Bj->elem(*ABsupi),p),p);
		  }
		if(Pij) Pi->entries[j]=Pij;
	      }
	  }
    }  
  //  cout<<"done"<<endl;
  return prod;
}
//#endif
#if(0)
smat mult_mod_p ( const smat& A, const smat& B, const scalar& p )
{
  cout<<"Multiplying mod-p smats of size "<<A.nro<<"x"<<A.nco<<" and "<<B.nro<<"x"<<B.nco<<"..."<<endl;
  //  cout<<"A = "<<A<<endl;
  //  cout<<"B = "<<B<<endl;
  if( A.nco != B.nro ) 
    { 
      cout << "incompatible smats in mult_mod_p(smat,smat,p)\n"; 
      cout << "Dimensions "<<dim(A)<<" and "<<dim(B)<<endl;
      abort();
    }
  int nro = A.nro, nco = B.nco, j;
  smat prod( nro, nco );
  smat tB = transpose(B);
  vector<svec>::const_iterator Ai;
  vector<svec>::iterator Bj, Pi;
  scalar Pij;
  for((Ai = A.rows.begin())++, (Pi = prod.rows.begin())++; Ai != A.rows.end(); Ai++, Pi++)
    {
      if(Ai->size()==0) continue;
      for((Bj = tB.rows.begin())++, j=1; Bj != tB.rows.end(); Bj++, j++)
	if(Bj->size()!=0) 
	  {
	    Pij = dotmodp(*Ai,*Bj,p);
	    if(Pij) Pi->entries[j]=Pij;
	  }
    }  
  cout<<"done"<<endl;
  //  cout<<"product = "<<prod<<endl;
  return prod;
}
#endif

int operator==(const smat& sm1, const smat& sm2)
{
   int nr = sm1.nro;
   if (nr != sm2.nro ) return 0;
   for(int i = 1; i <= nr; i++ )
     if(sm1.rows[i]!=sm2.rows[i]) return 0;
   return 1;
}

int operator==(const smat& sm, const mat& m)
{
   int nr = sm.nro;
   if (nr != nrows(m) ) return 0;
   if (sm.nco != ncols(m) ) return 0;
   for(int i = 1; i <= nr; i++ ) if(sm.rows[i]!=m.row(i)) return 0;
   return 1;
}

int eqmodp(const smat& sm1, const smat& sm2, const scalar& p)
{
   int nr = sm1.nro;
   if (nr != sm2.nro ) return 0;
   for(int i = 1; i <= nr; i++ )
     if(!eqmodp(sm1.rows[i],sm2.rows[i],p)) return 0;
   return 1;
}

ostream& operator << (ostream& s, const smat& sm)
{
  s<<"[";
  for( int i = 1; i <= sm.nro; i++ )
    {
      if(i>1) s<<";";
      s << sm.rows[i];
    }
  s<<"]";
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
  unsigned int count = 0;
  for(int r = 1; r <= m.nro; r++ ) 
    count += m.rows[r].size();
  return count;
}

vector<std::set<int> > row_supports(const smat& A)
{
  vector<std::set<int> > ans(A.nro+1); 
  for(int i=1; i<=A.nro; i++) ans[i]=A.rows[i].support();
  return ans;
}


// fill in the (presumed empty) vector v with at most max random
// numbers in range -10..10
void random_fill_in( svec& v, int max, scalar& seed )
{
  int N = int( (max+1) * ran0(seed) ); //number of entries in vector
  if( N == (max+1) ) N--;  // could occur !
  //  cout<<"random_fill_in(), max="<<max<<", seed="<<seed<<", dim(v)="<<dim(v)<<", N = "<<N<<endl;
  while(N--)
    {
      int i = 1+int( dim(v) * ran0(seed) );
      if(i>dim(v)) i--;
      v.set(i,int( 20 * ran0(seed) ) - 10);	
    }
  //  cout<<"random_fill_in() returns v="<<v<<endl;
}

void random_fill_in( smat& sm, int max, scalar& seed )
{
  for( int r = 1; r <= sm.nro; r++ )
    random_fill_in(sm.rows[r], max, seed);
}

smat sidmat(scalar n)  // identity matrix
{
  smat I(n,n);
  for(int i=1; i<=n; i++) I.rows[i].entries[i]=1;
  return I;
}

smat liftmat(const smat& mm, scalar pr, scalar& dd, int trace)
{
  scalar modulus=pr,n,d; dd=1;
  int succ,success=1;
  float lim=floor(sqrt(pr/2.0));
  smat m = mm;
  if(trace)
    {
      cout << "Lifting mod-p smat;  smat mod "<<pr<<" is:\n";
      cout << mm.as_mat();
      cout << "Now lifting back to Q.\n";
      cout << "lim = " << lim << "\n";
    }
  vector<svec>::iterator ri;
  map<int,scalar>::iterator rij;
  for((ri=m.rows.begin())++; ri!=m.rows.end(); ri++)
    for(rij=ri->begin(); rij!=ri->end(); rij++)
      {  
	succ = modrat(rij->second,modulus,lim,n,d);
	success = success && succ;
	dd=lcm(d,dd);
      }
  if(!success) 
    cout << "Problems encountered with modrat lifting of smat." << endl;
  dd=abs(dd);
  if(trace) cout << "Common denominator = " << dd << "\n";
  for((ri=m.rows.begin())++; ri!=m.rows.end(); ri++)
    for(rij=ri->begin(); rij!=ri->end(); rij++)
      rij->second=mod(xmodmul(dd,(rij->second),pr),pr); 
  if(trace) cout << "liftmat returns " << m.as_mat() << endl;
  return m;
}

