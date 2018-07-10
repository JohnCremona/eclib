// reduce_cubics.cc: Program for reducing integer binary cubic forms
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>

int main()
{
  initprimes("PRIMES");
  bigint a, b, c, d, disc;
  unimod m;
  while(cout << "Enter cubic coeffs a, b, c, d: ",
	cin >> a >> b >> c >> d,
	!(is_zero(a)&&is_zero(b)&&is_zero(c)&&is_zero(d)))
    {
      cubic g0(a,b,c,d);
      cubic g(g0);
      cout << "Input cubic = "<<g<<endl;;
      disc = g.disc();
      cout << "Discriminant = " << disc << endl;

      if(disc>0)
	{
	  cout << "Using Hessian to reduce...\n";
          g.hess_reduce(m);
	  cout << "Hessian reduced cubic = "<<g<<endl;
	  cout << "after transform by "<<m<<endl;
	  cout << "Root of Hessian = " << g.hess_root() << endl;
	}
      else
	{
	  bigfloat alpha = g.real_root();
	  //cout << "Real root alpha = " << alpha << endl;
	  bigfloat xa=I2bigfloat(a), xb=I2bigfloat(b), xc=I2bigfloat(c), xd=I2bigfloat(d);
	  bigfloat ga = ((xa*alpha+xb)*alpha+xc)*alpha+xd;
	  //cout << "g(alpha) = " << ga << endl;
// First use Mathews reduction
	  cout << "Using Mathews reduction ...\n";
          g.mathews_reduce(m);
	  cout << "Mathews reduced cubic = "<<g<<endl;
	  cout << "after transform by "<<m<<endl;

// Now use JC/Julia reduction
	  cout << "Using JC/Julia reduction ...\n";
	  g=g0;
          g.jc_reduce(m);
	  cout << "JC/Julia reduced cubic = "<<g<<endl;
	  cout << "after transform by "<<m<<endl;
	}
    }
  cout<<endl;
}

