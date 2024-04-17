// FILE APLIST.CC: program for listing ap
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

#include <fstream>
#include <iomanip>
#include <eclib/arith.h>
#include <eclib/moddata.h> // for nf_filename
#define BOOKORDER       // if defined, sorts newforms & curves into order
                        // in the Book (relevant up to 500 only)
#include <eclib/curvesort.h>

int main()
{
  string filename;
  cout<<"Enter output filename: ";
  cin>>filename;
  ofstream output(filename.c_str());
  int n1,n2;
  cout<<"Enter first and last values of n: ";
  cin>>n1>>n2;
  for (int n=n1; n<=n2; n++)
    {
      string data = nf_filename(n,'x');
      ifstream datafile(data.c_str());
      if (!datafile) cout<<"No file "<<data<<" exists!"<<endl;
      else
        {
          int i,ic,xic,ip,jp,p,iq,nnf,naq;
	  int nbigprimes = 5;  // max no. of primes over 100

	  int temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of newforms
	  nnf=temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of bad primes
	  naq=temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of eigs
	  if(nnf==0) continue;

	      // skip over extra data for each newform
	      int ntotal = 16*nnf;
	      int* batch_i = new int[ntotal];
	      datafile.read((char*)batch_i,ntotal*sizeof(int));
	      delete[] batch_i;

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
	      delete[] batch;
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
	      for (int kp=0; kp<25; kp++)
		{
		  int pk = prime_number(kp+1);
		  while (nn%pk==0) nn/=pk;
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
	      delete[] batch;
	      datafile.close();
	      //        cout<<"finished reading eigs, closed data file."<<endl;
          for (xic=0; xic<nnf; xic++)
            { output<<setw(3)<<n;
	      ic=xic;
              output<<" "<<codeletter(xic);
#ifdef BOOKORDER
	      ic=booknumber0(n,ic);
#endif
	      iq=0;
              for (int kp=0; kp<25; kp++)
                { int pk = prime_number(kp+1);
		  int ap = aptable[ic][kp];
                  if (n%pk==0) // W_q-eig
		    {
		      ap=aqtable[ic][iq++];
		      if(kp>8)  output<<" ";
		      if(ap==1) output<<"  +"; else output<<"  -";
		    }
                  else        // T_p-eig
		    {
		      if (kp>8) output<<setw(4)<<ap;
		      else output<<setw(3)<<ap;
		    }
                }
	      for(i=0; i<(int)bigptable.size(); i++)
                {
                  if (aqtable[ic][iq++]==1) output<<" +"; else output<<" -";
                  output<<"("<<bigptable[i]<<")";
                }
              output<<endl;
            } // end of xic/ic output loop
        }     // end of if(!data) ... else ...
    }         // end of main n loop
}             // end of main()
