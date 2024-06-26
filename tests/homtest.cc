// FILE HOMTEST.CC: Test program for homspace class
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
//
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/timer.h>

#define AUTOLOOP
//#define SHOW_TIMES

int main(void)
{
  init_time();
 int n=1; 
 int plus=1;
 int verbose=0;
 int cuspidal=0;
 vector<long> dims(3), cdims(3);
 // cout << "Verbose? "; cin >> verbose;
 // cout << "Plus space, minus space or full space (+1,-1,0)? "; cin >> plus;
 int s,limit; 
#ifdef AUTOLOOP
  cout<<"Enter first and last levels: ";cin>>n>>limit; n--;
  while (n<limit) { n++;
#else
  while (n>1) { cout<<"Enter level: "; cin>>n;
#endif
 if (n>1)
{
  {
    cout << "\n>>>Level " << n << "\n";
    for (s=0; s<3; s++)
      {
        plus=(s==0? 1 : (s==1? -1: 0));
        if(verbose)
          {
            cout<<"Computing sign="<<plus<<" space"<<endl;
          }
	start_time();
	homspace hplus(n,plus, cuspidal,verbose);
	stop_time();
	int dim = hplus.h1dim();
	int cdim = hplus.h1cuspdim();
	dims[plus+1]=dim;
	cdims[plus+1]=cdim;
	int d = hplus.h1denom();
	int cd = hplus.h1cdenom();
	cout << "Sign = "; if (plus!=-1) cout<<" ";
	cout << plus << ": ";
	cout << "\tDimension = " << dim;
	cout << "\tCuspidal dimension = " << cdim;
	if(d*cd>1 &&verbose) cout<<" denoms ("<<d<<","<<cd<<")";
#ifdef SHOW_TIMES
	show_time();
#endif
 	cout<<endl;
#ifdef SHOW_TIMES
	cout<<"***************************************************"<<endl;
#endif
      }
#ifdef SHOW_TIMES
    cout<<"***************************************************"<<endl;
#endif
    if (dims[1]==dims[0]+dims[2])
      continue; //cout<<"\tDimensions add up OK at level "<<n;
    else
      cout<<"****************Dimensions inconsistent for level "<<n<<" ***************"<<endl;
    if ((cdims[1]==cdims[0]+cdims[2])&&(cdims[0]==cdims[2]))
      continue; //cout<<"\tCuspidal dimensions add up OK at level "<<n<<endl;
    else
      cout<<"\t****************Cuspidal dimensions inconsistent for level "<<n<<" ***************"<<endl;
  }
  cout<<endl;
  //  stop_time();
  //  show_time();
}       // end of if(n)
}       // end of while()
}       // end of main()
