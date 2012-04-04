// tlegcert.cc: test program for solving legendre equations with certificates
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
 
#include "marith.h"
#include "quadratic.h"
#include "conic.h"
#include "legendre.h"

#ifndef CONIC_METHOD
#define CONIC_METHOD 4
#endif
//#define TEST_PARAM
//#define HOLZER_MEASURES

int main()
{
  initprimes(string("PRIMES").c_str());
  cout<<"Solving ax^2 + by^2 + cz^2 = 0\n";
  cout<<"Using method "<<CONIC_METHOD<<endl<<endl;

  bigint a,b,c,x,y,z;
  bigint u,k1,k2,k3;
  quadratic qx, qy, qz;

  while(1) 
    {
      cout << "Enter coefficients a b c: ";
      cin >> a >> b >> c;
      cout<<a<<" "<<b<<" "<<c<<endl; 
      if(a==0) abort();
      cout << "Enter certificate k1 k2 k3: ";
      cin >> k1 >> k2 >> k3;
      int use_lll=(CONIC_METHOD==5);
      int res=!use_lll;
      if(use_lll) legendre_via_lll(a,b,c,k1,k2,k3,x,y,z);
      else        res = legendre_solve_cert_1(a,b,c,k1,k2,k3,x,y,z,u);
      if(res)
//    if(!solve_conic(a,0,c,-b,x,y,z,CONIC_METHOD))
	{
	  cout<<"No solution!\n";
	}
      else
	{
#ifdef HOLZER_MEASURES
	  cout<<"Before reduction of solution "; 
	  show_xyz(x,y,z); cout<<endl;
	  cout<<"Holzer measure = "<<holzer_measure(a,b,c,x,y,z)<<endl;
	  cancel1(x,y,z);
	  new_legendre_reduce(a,b,c,x,y,z,0);
	  cout<<"After reduction: ";show_xyz(x,y,z);cout<<endl;
	  cout<<"Holzer measure = "<<holzer_measure(a,b,c,x,y,z)<<endl;
#else
	  cout << "Solution: "; show_xyz(x,y,z); cout<<endl;
#endif // HOLZER_MEASURES
	  

      if(check_leg(a,b,c,x,y,z))
	{
	  cout<<" --OK\n";
#ifdef TEST_PARAM
	  legendre_param(a,b,c,x,y,z,qx,qy,qz);
	  //	  cout<<"Parametric solution:\n";
	  cout<<"x = "<<qx<<" * [u^2,uv,v^2]\n";
	  cout<<"y = "<<qy<<" * [u^2,uv,v^2]\n";
	  cout<<"z = "<<qz<<" * [u^2,uv,v^2]\n";
	  if(testparamsol(a,0,c,-b,qx,qy,qz,0))
	    cout<<"Parametric solution is OK\n";
	  else 
	    cout<<"Parametric solution is wrong!\n";
#endif
	}
      else
	cout<<" -- wrong!\n";
    }      
    }
  abort();
}



