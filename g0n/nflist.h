// FILE NFLIST.H: declaration of class nflist
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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

bool less_ap(const vector<short>v, const vector<short>w);
bool less_nf(const pair<long,vector<short> > f1, 
	     const pair<long,vector<short> > f2);
struct less_newform 
  : public binary_function<const pair<long,vector<short> >,
			   const pair<long,vector<short> >, bool > 
{ 
  bool operator()
    (const pair<long,vector<short> >& f, 
     const pair<long,vector<short> >& g) const
  {
    return less_nf(f,g);
  } 
};

class nflist{
private:
  long n_max;      // top level
  vector<long> ** nf;  // the eigs
  vector<long> nnf;    // the numbers of newforms at each level
  long nap;        // number of eigs for each newform
  vector<long> plist;      // list of nap primes
  long q;                  // Current modulus
  multimap<pair<long, vector<short> >,
	   pair<long,long>, 
	   less_newform > mod_q_table;
  // map from pair(level, eiglist mod q) to pair(level, number)
  int cong1(long n1, const vector<long>& ap1, 
	    long n2, const vector<long>& ap2, long q);

public:
  nflist(long nm, long np);  // constructor, nm=n_max
  ~nflist(void);
  vector<long> form(long level, long number)
    {return nf[level][number];}
  void print1(long level, long number)
    { cout<<level<<"\t"<<(number+1)<<"\t"<<nf[level][number]<<"\n";}
  void printall(long n1, long n2)
    {
      for (long n=n1; n<=n2; n++) 
        for (long i=0; i<nnf[n]; i++)
          cout<<n<<"\t"<<(i+1)<<"\t"<<nf[n][i]<<"\n";
    }
  int prepare(long qq);  // set up mod-q table
  int nflist::find_dups(long qq);  // look for duplicates in mod-q table
  int cong(long n, long iform, long q, long& m, long& jform);
  void cong_bothways(long n, long iform, long q, long& m, long& jform);
  int find_cong(long n, long q);
  int find_cong_all(long q)
  { for(long n=11; n<=n_max; n++)  
    {
      //      cout<<"."<<flush; 
      find_cong(n,q);
      if(n%1000==0)cout<<"["<<n<<"]"<<endl;
    }
  }
};

inline void report(long n, long iform, long m, long jform, long q)
{
  cout << n<<"#"<<(iform+1)<<"\t = "
       << m<<"#"<<(jform+1)<<"\t (mod "<<q<<")\n";

}
