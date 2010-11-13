// FILE CHECKAP.CC: program for checking ap are in valid range
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

// (Now essentially redundant as we do the check as they are computed)

#include <fstream>
using namespace std;
#include "arith.h"

char* nf_filename(long n, char c)
{
  char* nf_dir = getenv("NF_DIR"); 
  string nf_file;
  if (nf_dir==NULL) 
    nf_file = string("./newforms");
  else
    nf_file = string(nf_dir);
  char* filename=new char[20];
  sprintf(filename,"%s/%c%ld",nf_file.c_str(),c,n);
  return filename;
}

const long MAXNAP = 100000;

int main(void)
{
  char* name = new char[20];
  long firstn, lastn, n, nnf, naq, nap, n2ds, ap, i, j, p;
  int ok, allok=1;
  short temp;
  cout<<"Enter first and last N: "; cin>>firstn>>lastn;
  for(n=firstn; n<=lastn; n++)
    {
      ok=1;
      name = nf_filename(n,'x');
      ifstream datafile(name);
      if(!datafile.is_open())
	{
	  cout<<"\nFile "<<name<<" does not exist!"<<endl;
	}
      else 
	{
	  short temp_short;
	  int temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of newforms
	  nnf=temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of bad primes
	  naq=temp_int;
	  datafile.read((char*)&temp_int,sizeof(int));   // = number of eigs
	  nap=temp_int;
	  if(nnf>0) 
	    {
	      // skip over extra data for each newform
	      long ntotal = 16*nnf;
	      int* batch_i = new int[ntotal];
	      datafile.read((char*)batch_i,ntotal*sizeof(int));

	      // skip over aq for each newform
	      ntotal = naq*nnf;
	      short* batch = new short[ntotal];
	      datafile.read((char*)batch,ntotal*sizeof(short));
	      
	      // read and check ap for each newform
	      ntotal = nap*nnf;
	      delete batch;
	      batch = new short[ntotal];
	      datafile.read((char*)batch,ntotal*sizeof(short));
	      short* batchptr = batch;
	      primevar pr;
	      for(j=0; j<nap; j++, pr++) 
		{
		  p = (long)pr;
		  for(i=0; i<nnf; i++) 
		    {
		      long ap=*batchptr++;
		      if(ap*ap>4*p)
			{
			  cout << "\nError: N = " << n << ", form # " << (j+1) 
			       << ", p = " << p << " has ap = " << ap << flush;
			  ok=0;
			}
		    }
		}
	      delete batch;
	    } // ends if(nnf>0)
	} // ends if(datafile)
	    
      datafile.close();
      if(ok) cout << n << " ";
      else
	{
	  cout << "\nN = " << n << ": errors found!\n";
	  allok=0;
	}
      if(n%10==0) cout<<endl;
    }
  
  if(allok) cout << "\nAll checked ok\n"; else cout << "\nerrors found\n";
  delete[]name;
}
