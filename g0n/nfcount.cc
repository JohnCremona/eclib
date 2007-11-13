// FILE NFCOUNT.CC:  program to count newforms (from data files)
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
  sprintf(filename,"%s/%c%d",nf_file.c_str(),c,n);
  return filename;
}

int main(void)
{
  int limit,n=1,count=0,firstn, verbose; 
  char* name = new char[20];
  int getindex=0;
  cout << "Verbose? (0/1) " << endl;  cin >> verbose;
  cout<<"Enter first and last N: ";cin>>firstn>>limit; 
  n=firstn-1; cout<<endl;
  while (n<limit) { n++;
  if (n>0)
{
  if(verbose) cout << ">>> Level " << n << " <<< ";
  int num=0, nap=0, naq=0;
  name = nf_filename(n,'x');
  ifstream in(name);
  int eig_file_exists = in.is_open();
  if(!eig_file_exists)
    {
      //      cout<<"Unable to open file "<<name<<" for newform input"<<endl;
    }
  else
    {
      in.read((char*)&num,sizeof(int));
      count+=num;
      in.read((char*)&naq,sizeof(int));   // = number of bad primes
      in.read((char*)&nap,sizeof(int));   // = number of eigs
      in.close();
    }
  if(verbose) 
    {
      cout << num << " newform(s), "<<nap<< " eigs on file.";
      if((num>0)&&(nap<168)) cout<<" !!!";
      if(!eig_file_exists) cout<<" [No eig file exists]";
      cout<< endl;
    }
}       // end of if(n)
}       // end of while()
cout << "\n" << count << " newform(s) in range " << firstn << "..." << limit << endl;
}       // end of main()
