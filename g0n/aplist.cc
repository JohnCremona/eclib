// FILE APLIST.CC: program for listing ap
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
#define BOOKORDER       // if defined, sorts newforms & curves into order
                        // in the Book (relevant up to 500 only)
#include "curvesort.cc"

int main()
{
  char* code = new char[20];
  char filename[20];
  cout<<"Enter output filename: ";
  cin>>filename;
  ofstream output(filename);
  int n1,n2; short temp;
  cout<<"Enter first and last values of n: ";
  cin>>n1>>n2;
  for (int n=n1; n<=n2; n++)
    { 
      sprintf(filename,"%s/x%d",NF_DIR,n);
      ifstream datafile(filename);
      if (!datafile) cout<<"No file "<<filename<<" exists!"<<endl;
      else
        {
//        cout<<"reading eigs from file "<<filename<<"; ";
          int i,ic,xic,ip,jp,p,iq,nnf,naq,nap,nx; 
	  int nbigprimes = 5;  // max no. of primes over 100

	  short temp_short;
	  int temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of newforms
	  nnf=temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of bad primes
	  naq=temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of eigs
	  nap=temp_int;
	  if(nnf==0) continue;

	      // skip over extra data for each newform
	      long ntotal = 16*nnf;
	      int* batch_i = new int[ntotal];
	      datafile.read((char*)batch_i,ntotal*sizeof(int));
	      delete batch_i;

	      // read and store aq for each newform
	      ntotal = naq*nnf;
	      short* batch = new short[ntotal];
	      datafile.read((char*)batch,ntotal*sizeof(short));
	      short* batchptr = batch;
	      vector< vector<long> > aqtable;
	      for(ic=0; ic<nnf; ic++) aqtable.push_back(vector<long>(naq));
	      for (iq=0; iq<naq; iq++)
		for (ic=0; ic<nnf; ic++) aqtable[ic][iq]=*batchptr++;
	      
	      // read and store ap for each newform
	      ntotal = 25*nnf;
	      delete batch;
	      batch = new short[ntotal];
	      datafile.read((char*)batch,ntotal*sizeof(short));
	      batchptr = batch;
	      vector< vector<long> > aptable;
	      vector<long> bigptable;
	      for(ic=0; ic<nnf; ic++) 
		aptable.push_back(vector<long>(25+nbigprimes));
	      for (ip=0; ip<25; ip++)
		for (ic=0; ic<nnf; ic++) 
		  aptable[ic][ip]=*batchptr++;

	      // deal with big bad primes (>100) if any:
	      int nn=n;
	      ip=25; // index into aptable
	      jp=26; // current prime being considered
	      p=prime_number(jp);
	      for (ip=0; ip<25; ip++)
		{
		  int p = prime_number(ip+1);
		  while (nn%p==0) nn/=p;
		}
	      while (nn>1)     // then there are primes>100 dividing n 	       
		{
		  while (nn%p!=0) // skip over data for p & increment p
		    {
		      datafile.read((char*)batch,nnf*sizeof(short));
		      p=prime_number(++jp);
		    }
		  // now p is a "large" prime divisor of n
		  bigptable.push_back(p);
		  while (nn%p==0) nn/=p;		    
		  datafile.read((char*)batch,nnf*sizeof(short));
		  batchptr=batch;
		  for(ic=0; ic<nnf; ic++) aptable[ic][ip]=*batchptr++;
		  ip++;
		}
	      delete batch;
	      datafile.close();
	      //        cout<<"finished reading eigs, closed data file."<<endl;
          for (xic=0; xic<nnf; xic++)
            { output<<setw(3)<<n;
	      ic=xic;
	      codeletter(xic,code);
	      output<<" "<<code;
#ifdef BOOKORDER
	      ic=booknumber0(n,ic);
#endif
	      iq=0;
              for (int jp=0; jp<25; jp++)
                { int p = prime_number(jp+1);  
		  int ap = aptable[ic][jp];
                  if (n%p==0) // W_q-eig
		    {
		      ap=aqtable[ic][iq++];
		      if(jp>8)  output<<" ";
		      if(ap==1) output<<"  +"; else output<<"  -";
		    }
                  else        // T_p-eig
		    {
		      if (jp>8) output<<setw(4)<<ap; 
		      else output<<setw(3)<<ap;
		    }
                }
	      for(i=0; i<bigptable.size(); i++)
                {
                  if (aqtable[ic][iq++]==1) output<<" +"; else output<<" -";
                  output<<"("<<bigptable[i]<<")";
                }
              output<<endl;
            } // end of xic/ic output loop
        }     // end of if(!data) ... else ...
    }         // end of main n loop
}             // end of main()
