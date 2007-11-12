// FILE NFLIST.CC: implementation of class nflist
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

#include <fstream>
#include <iomanip>
#include "arith.h"
#include "nflist.h"

#define VERBOSE 0
//#define STRICT      // If defined, uses a stricter congruence criterion 
                    // (see below) 

nflist::nflist(long nm, long np) :n_max(nm), nap(np)
{
  plist = primes(nap);
  nf = new vector<long>*[nm+1];  // indexed from 0 to nm but only use from 11
  nnf.resize(nm+1);              // ditto
  long i, j, n, nforms, naq, neigsonfile, dummy;
  short *batch, *batchptr;
  char* name = new char[20];

  for (n=11; n<=nm; n++)
    {  
      //#if VERBOSE>0      
      //      cout<<"."<<flush;
      if(n%1000==0)cout<<"["<<n<<"]"<<endl;
      //#endif
      nforms=0;
      sprintf(name,"newforms/x%d\0",n);
      ifstream datafile(name);
      if(!datafile.is_open())
        {
          cout<<"Unable to open file "<<name<<" for newform input"<<endl;
        }
      else
        {
	  short temp_short;
	  int temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of newforms
	  nforms=temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of bad primes
	  naq=temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of eigs
	  neigsonfile=temp_int;
        }
      nnf[n]=nforms;
#if VERBOSE>1      
      cout<<nforms<<" at level "<<n<<" with "<<neigsonfile<<" eigs on file\n";
#endif
      
      nf[n]=new vector<long>[nforms];  // indexed from 0
      for(i=0; i<nforms; i++)     nf[n][i].resize(nap);
#if VERBOSE>2
      cout<<"created vector<long>s to hold them\n";
      cout<<"Before reading data:\n";
      for(i=0; i<nforms; i++) cout<<nf[n][i]<<endl;
#endif

	  if(nforms==0) continue;
	  else
	    {
	      // skip over extra data for each newform
	      long ntotal = 16*nforms;
	      int* batch_i = new int[ntotal];
	      datafile.read((char*)batch_i,ntotal*sizeof(int));
	      delete batch_i;

	      // skip over aq for each newform
	      ntotal = naq*nforms;
	      batch = new short[ntotal];
	      datafile.read((char*)batch,ntotal*sizeof(short));

	      // read and store ap for each newform
	      ntotal = nap*nforms;
	      delete batch;
	      batch = new short[ntotal];
	      datafile.read((char*)batch,ntotal*sizeof(short));

	    }

      datafile.close();
      batchptr = batch;
      for(j=0; j<nap; j++)
        for(i=0; i<nforms; i++) 
          nf[n][i][j] = *batchptr++;
#if VERBOSE>2
      cout<<"After reading data:\n";
      for(i=0; i<nforms; i++) cout<<nf[n][i]<<endl;
#endif
    } // end of level loop
  delete[] name;

} // end of nflist constructor

nflist::~nflist()
{ 
  for(long n=11; n<=n_max; n++) delete [] nf[n];
  delete [] nf;
}

int nflist::cong1(long n1, const vector<long>& ap1, 
		  long n2, const vector<long>& ap2, long q)
// tests whether newforms ap_i al levels n_i are "almost congruent" mod q
// in the sense of Mazur
{
//cout<<"Comparing "<<ap1<<" and "<<ap2<<"\t";
  long i;  int ok=1;
  for(i=0; ok&&(i<nap); i++)
    {
      long p = plist[i], ap1i=ap1[i], ap2i=ap2[i];
#ifdef STRICT
      long p2=p*p;
      if(::div(p,n1)) {if(::div(p2,n1)) {ap1i=0;} else ap1i=-ap1i;}
      if(::div(p,n2)) {if(::div(p2,n2)) {ap2i=0;} else ap2i=-ap2i;}
      ok = ::div(q,ap1i-ap2i);
#else
      ok = ::div(p,n1) || ::div(p,n2) || ::div(p,q) || ::div(q,ap1i-ap2i);
#endif
    }
//cout<<"ok="<<ok<<endl;
  return ok;
}

int nflist::cong(long n, long iform, long q, long& m, long& jform)
//searches for an earlier form almost congruent to this one mod q,
//returns 1 and level in m, number in jform if found
{
  vector<long> ap=nf[n][iform];
 
  for(m=11; m<n; m++)
    for(jform=0; jform<nnf[m]; jform++)
      {
//	cout<<"Testing n,i,m,j = "<<n<<","<<iform<<","<<m<<","<<jform<<endl;
        if(cong1(n,ap,m,nf[m][jform],q)) return 1;
      }

  m=n;
  for(jform=0; jform<iform; jform++)
    {
//      cout<<"Testing n,i,m,j = "<<n<<","<<iform<<","<<m<<","<<jform<<endl;
      if(cong1(n,ap,m,nf[m][jform],q)) return 1;
    }
  return 0;
}

void nflist::cong_bothways(long n, long iform, long q, long& m, long& jform)
//searches for an earlier or later form almost congruent to this one mod q,
//returns 1 and level in m, number in jform if found
//We check the other forms at the same level first
{
  if((iform<0)||(iform>=nnf[n])) 
    {
      cout<<"Form number must be between 1 and "<<nnf[n]<<"!\n"; 
      return;
    }
  vector<long> ap=nf[n][iform];
  int found=0;

  m=n;
  for(jform=0; jform<nnf[n]; jform++)
    {
      if(jform==iform) continue;
//      cout<<"Testing n,i,m,j = "<<n<<","<<iform<<","<<m<<","<<jform<<endl;
      if(cong1(n,ap,m,nf[m][jform],q)) 
	{
	  found=1;
	  if((q==2)&&cong1(n,ap,m,nf[m][jform],4)) 
	    report(n,iform,m,jform,4);
	  else report(n,iform,m,jform,q);
	  if((q==2)&&(m==14)&&(jform==0)) 
	    {
	      cout<<"Omitting other congruences mod 2\n";
	      return;
	    }
	}
    }
  //if(found) return;  // We only look at other levels if there's nothing
                     // at the same level. 
  int go_on=0;
  cout<<"Search other levels? "; cin>>go_on;
  if(!go_on)return;

  for(m=11; m<=n_max; m++)
    {
      if(m==n) continue;
      for(jform=0; jform<nnf[m]; jform++)
	{
//	  cout<<"Testing n,i,m,j = "<<n<<","<<iform<<","<<m<<","<<jform<<endl;
	  if(cong1(n,ap,m,nf[m][jform],q)) 
	    {
	      if((q==2)&&cong1(n,ap,m,nf[m][jform],4)) 
		report(n,iform,m,jform,4);
	      else report(n,iform,m,jform,q);
	      if((q==2)&&(m==14)&&(jform==0)) 
		{
		  cout<<"Omitting other congruences mod 2\n";
		  return;
		}
	    }
	}
    }
}

int nflist::find_cong(long n, long q)
//for each form at level n, 
//searches for an earlier form almost congruent to this one mod q,
//outputting any found
{
  long m, jform, iform;
  int found, foundany=0;
  for (iform=0; iform<nnf[n]; iform++)
   {
     found = cong(n,iform,q,m,jform);
     if(found) 
      {
       foundany=1;
       //       cout<<endl;
       cout <<setw(6)<< n<<"#"<<(iform+1)<<" = "
	    <<setw(6)<< m<<"#"<<(jform+1)<<"\t (mod "<<q<<")\n";
       /*
       cout << "level "<<n<<", form "<<(iform+1)<<" congruent mod "<<q<<" to "
            << "level "<<m<<", form "<<(jform+1)<<endl;
       */
      }
   }    
//  if(!foundany&&((n%100)==0)) 
//    cout<<"Level "<<n<<": no congruences found."<<endl;
}

int nflist::prepare(long qq)  // set up mod-q table
{
  q = qq;
  long n,i,j;
  vector<long> ap(nap);
  vector<long> pr=primes(nap);
  vector<short> ap_mod_q(nap);
#if VERBOSE>0
  cout<<"Preparing mod "<<q<<" table..."<<flush;
#endif
  for(n=11; n<= n_max; n++)
    {
      for(i=0; i<nnf[n]; i++)
	{
	  ap=nf[n][i];
	  for(j=0; j<nap; j++) 
	    {
	      ap_mod_q[j]=(short)(posmod(ap[j],q));
	      if(((q*n)%pr[j])==0) ap_mod_q[j]=(short)(0);
	    }
	  mod_q_table.insert(pair< pair<long,vector<short> >,pair<long,long> >
			     (
			     pair<long,vector<short> >(n,ap_mod_q),
			     pair<long,long>(n,i+1)
			     ));
	} 
    }
#if VERBOSE>0
  cout<<"done"<<endl;
#endif
}

bool less_ap(const vector<short>v, const vector<short>w)
{
  vector<short>::const_iterator vi=v.begin(), wi=w.begin();
#if VERBOSE>2
  cout<<"Comparing "<<v<<" and "<<w<<"..."<<flush;
#endif
  while(vi!=v.end())
    {
      bool s = (*vi++)<(*wi++); 
      if(s) 
	{
#if VERBOSE>2
	  cout<<"result = "<<s<<endl;
#endif
	  return s; 
	}
    }
#if VERBOSE>2
  cout<<"result = "<<0<<endl;
#endif
  return 0;
}

bool less_nf(const pair<long,vector<short> > f1, 
	     const pair<long,vector<short> > f2)
{
  return less_ap(f1.second,f2.second);  
}

int nflist::find_dups(long qq)  // look for duplicates in mod-q table
{
#if VERBOSE>0
  cout<<"Looking for duplicates mod "<<q<<" ..."<<flush;
#endif
  for (multimap<pair<long,vector<short> >,pair<long,long>, less_newform >::iterator 
	 it = mod_q_table.begin(); it != mod_q_table.end();  ++it)
    {
      pair<long, vector<short> > N_ap = it->first;
      pair<long,long> N_i = it->second;
#if VERBOSE>1
      cout<<"Counting N= "<<N_ap.first<<", ap="<<N_ap.second<<" ..."<<flush;;
#endif
      int count = mod_q_table.count(N_ap);
#if VERBOSE>1
      cout<<"done, number is "<<count<<endl;
#endif
      cout<<N_ap.first<<": "<<N_ap.second
	  <<" ("<<N_i.first<<","<<N_i.second<<")"
	  << ": # = " << count
	  << endl;
      if(count>1)
	{
	  cout << "  [" << N_ap.second << ", ";
	  cout << "]" << endl;
	}
    }
#if VERBOSE>0
  cout<<"done"<<endl;
#endif
 
}
