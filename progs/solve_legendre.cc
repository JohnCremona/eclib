// solve_legendre.cc:  program for solving legendre equations
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
 
#include <eclib/marith.h>
#include <eclib/quadratic.h>
#include <eclib/conic.h>
#include <eclib/legendre.h>

#ifndef CONIC_METHOD
#define CONIC_METHOD 4
#endif
#define TEST_PARAM

int main()
{
  initprimes("PRIMES");
  cout<<"Solving ax^2 + by^2 + cz^2 = 0\n";
  cout<<"Using method "<<CONIC_METHOD<<endl<<endl;

  bigint a,b,c,x,y,z,zero; zero=0;
  quadratic qx, qy, qz;

  while(1) 
    {
      cout << "Enter coefficients a b c: ";
      cin >> a >> b >> c;
      cout<<a<<" "<<b<<" "<<c<<endl; 
      if(a==0) exit(0);
      int use_lll = (CONIC_METHOD==5);
      if(!legendre_solve(a,b,c,x,y,z,use_lll))
//    if(!solve_conic(a,zero,c,-b,x,y,z,CONIC_METHOD))
	{
	  cout<<"No solution!\n";
	}
      else
	{
      cout << "Solution: "; show_xyz(x,y,z);
      if(check_leg(a,b,c,x,y,z))
	{
	  cout<<" --OK\n";
#ifdef TEST_PARAM
	  legendre_param(a,b,c,x,y,z,qx,qy,qz);
	  //	  cout<<"Parametric solution:\n";
	  cout<<"x = "<<qx<<" * [u^2,uv,v^2]\n";
	  cout<<"y = "<<qy<<" * [u^2,uv,v^2]\n";
	  cout<<"z = "<<qz<<" * [u^2,uv,v^2]\n";
	  cout<<"Disc(qx) = "<<qx.disc()<<endl;
	  cout<<"Disc(qy) = "<<qy.disc()<<endl;
	  cout<<"Disc(qz) = "<<qz.disc()<<endl;
	  if(testparamsol(a,zero,c,-b,qx,qy,qz,0))
	    cout<<"Parametric solution is OK\n";
	  else 
	    cout<<"Parametric solution is wrong!\n";
#endif
	}
      else
	cout<<" -- wrong!\n";
    }      
    }
}



