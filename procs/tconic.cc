// tconic.cc: conic test program
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

#ifndef VERBOSITY
#define VERBOSITY 0
#endif
#ifndef CONIC_METHOD
#define CONIC_METHOD 4
#endif
#define TEST_PARAM

int main()
{
  initprimes("PRIMES",VERBOSITY);
  cout<<"Solving ax^2 + bxz + cz^2 = dy^2\n";
  cout<<"Using method "<<CONIC_METHOD<<endl<<endl;

  bigint a,b,c,d,x0,y0,z0,disc;
  quadratic q, qx, qy, qz;

  while(1) {
  cout << "Enter coefficients a b c d: ";
  cin >> a >> b >> c >> d;
  cout<<a<<" "<<b<<" "<<c<<" "<<d<<endl; 
  if(d==0) {break;}
  q=quadratic(a,b,c);
  int res = solve_conic(q,d,x0,y0,z0,CONIC_METHOD);
  if(res) 
    {
      cout << "Solution: "; show_xyz(x0,y0,z0); cout<<endl;
      if(testsol(a,b,c,d,x0,y0,z0,0))
	cout<<"Solution is correct!\n";
      else
	cout<<"Solution is WRONG\n";
#ifdef TEST_PARAM
      res = solve_conic_param(q,d,qx,qy,qz,CONIC_METHOD,VERBOSITY);
      x0=qx[0]; y0=qy[0]; z0=qz[0]; cancel(x0,y0,z0);
      cout << "Solution: "; show_xyz(x0,y0,z0); cout<<endl;
      cout << "Parametric solution:\n";
      cout << "x = ["<<qx[0]<<","<<qx[1]<<","<<qx[2]<<"]*[u^2,uv,v^2]\n";
      cout << "y = ["<<qy[0]<<","<<qy[1]<<","<<qy[2]<<"]*[u^2,uv,v^2]\n";
      cout << "z = ["<<qz[0]<<","<<qz[1]<<","<<qz[2]<<"]*[u^2,uv,v^2]\n";
      
      bigint dqx = qx.disc();
      cout<<"disc(qx) = "<<dqx;
      if(dqx==4*c*d) cout<<" = 4cd\n";
      else cout<<" --NOT equal to 4cd = " << (4*c*d) <<endl;
      bigint dqz = qz.disc();
      cout<<"disc(qz) = "<<dqz;
      if(dqz==4*a*d) cout<<" = 4ad\n";
      else cout<<" --NOT equal to 4ad = " << (4*a*d) <<endl;
      bigint result = resultant(qx,qz);
      cout<<"resultant(qx,qz)   = "<<result;
      bigint res2 = sqr(d)*q.disc();
      if(result==res2) cout<<" = d^2(b^2-4ac)\n";
      else cout<<" --NOT equal to d^2(b^2-4ac) = "<<res2<<endl;
      if(testparamsol(a,b,c,d,qx,qy,qz,0))
	cout<<"Parametric solution is correct!\n";
      else
	cout<<"Parametric solution is WRONG\n";
#endif
    }
  else cout << "No solution.\n";
  }
  abort();
}



