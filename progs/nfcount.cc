// FILE NFCOUNT.CC:  program to count newforms (from data files)
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
#include <eclib/arith.h>
#include <eclib/moddata.h> // for nf_filename

int main(void)
{
  int limit,n=1,count=0,firstn, verbose; 
  int getindex=0;
  cout << "Verbose? (0/1) " << endl;  cin >> verbose;
  cout<<"Enter first and last N: ";cin>>firstn>>limit; 
  n=firstn-1; cout<<endl;
  while (n<limit) { n++;
  if (n>1)
{
  if(verbose) cout << ">>> Level " << n << " <<< ";
  int num=0, nap=0, naq=0;
  string name = nf_filename(n,'x');
  ifstream in(name.c_str());
  int eig_file_exists = in.is_open();
  if(!eig_file_exists)
    {
      cout<<"Unable to open file "<<name<<" for newform input"<<endl;
    }
  else
    {
      //      cout<<"Opened file "<<name<<" for newform input"<<endl;
      in.read((char*)&num,sizeof(int));
      count+=num;
      in.read((char*)&naq,sizeof(int));   // = number of bad primes
      in.read((char*)&nap,sizeof(int));   // = number of eigs
      in.close();
    }
  if(verbose) 
    {
      cout << num << " newform(s), "<<nap<< " eigs on file.";
      if((num>0)&&(nap<25)) cout<<" !!!";
      if(!eig_file_exists) cout<<" [No eig file exists]";
      cout<< endl;
    }
}       // end of if(n)
}       // end of while()
cout << "\n" << count << " newform(s) in range " << firstn << "..." << limit << endl;
}       // end of main()
