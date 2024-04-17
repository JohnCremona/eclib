// tilll.cc: test program for illl integer lll reduction
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
 
#include <eclib/mvector.h>
#include <eclib/mmatrix.h>
#include <eclib/illl.h>

int main()
{
  int n;
  cout<<"Enter size n: "; cin>>n;
  vector<vec_m> b(n+1);
  b[0]=vec_m(n);
  for(int i=1; i<=n; i++)
    {
      b[0][i]=1;  // these are the weights
      b[i]=vec_m(n);
      cout<<"Enter vector number "<<(i)<<": ";
      cin>>b[i];
    }
  cout<<"Before reduction, vectors are:\n";
  for(int i=1; i<=n; i++) cout<<b[i]<<endl;
  cout<<endl;

  cout << "FIRST METHOD: JC'S implementation of integer LLL from HC's book\n";

  lll_reduce(n,b);

  cout<<"After reduction, vectors are:\n";
  for(int i=1; i<=n; i++) cout<<b[i]<<endl;
  cout<<endl;


  vec_m shortest=b[1];
  bigint min_length=sdot(b,1,1);
  for(int i=1; i<=n; i++)
    cout<<"Square length of vector "<<i<<" is "<<sdot(b,i,i)<<endl;

  if(n==3) // then we know that any vector shorter than b[1] must be
           // of the form a1*b[1]+a2*b[2]+a3*b[3] with all ai in {-1,0,1}
    {
      cout<<"Candidates for shortest vectors:\n";
      int ok=1, better=0;
      for(int i=1; ok&&(i>-2); i--)
	for(int j=1; ok&&(j>-2); j--)
	  for(int k=1; ok&&(k>-2); k--)
	    {
	      if((i==0)&&(j==0)&&(k==0)) {ok=0;break;}
	      vec_m v=i*b[1]+j*b[2]+k*b[3];
	      bigint norm = sqr(v[1])+sqr(v[2])+sqr(v[3]);
	      cout<<"("<<i<<","<<j<<","<<k<<"): "<<v<<", norm = "<<norm<<endl;
	      if(norm<min_length) {min_length=norm; shortest=v; better=1;}
	    }
      cout<<"The shortest vector is "<<shortest
	  <<" with square length "<<min_length<<endl;
      if(better)
	cout<<"-- shorter than b[1]!"<<endl;
    }
}
