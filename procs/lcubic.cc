// lcubic.cc: Program for listing integer cubics with given discriminant
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
#include "unimod.h"
#include "cubic.h"

int main()
{
  initprimes(string("PRIMES").c_str());
  bigint a, b, c, d, disc;
  bigint absdisc, maxdisc;
  bigint a0,b0,c0,d0;
  bigint alim, blim, a2, b2, b3, cmin, cmax, r;
  bigint P, U, U2;
  bigfloat rdisc, ax, cx, cy;
  int neg;
  bigfloat fac1, fac2;
  unimod m;

  while(cout << "Enter discriminant bound (positive or negative): ",	cin >> maxdisc, !is_zero(maxdisc))
    {
      neg=(maxdisc<0);
      if(neg) 
	{
	  ::negate(maxdisc);
#ifdef LiDIA_ALL
	  fac1 = sqrt(bigfloat(8,27));
	  fac2 = power(bigfloat(2),inverse(bigfloat(3)));
#else
	  fac1 = sqrt((double)8)/sqrt((double)27);
	  fac2 = 1.2599210498948731647672106072782283505;
#endif
	  cout << "Negative discriminants down to " << maxdisc << endl;
	}
      else
	{
	  fac1 = sqrt(to_bigfloat(8))/to_bigfloat(3);
	  fac2 = to_bigfloat(1);
	  cout << "Positive discriminants  up  to " << maxdisc << endl;
	}
      
      for(absdisc=1; absdisc<=maxdisc; absdisc++)
	    {
	      disc=absdisc;
	      if(neg) ::negate(disc);
//	      cout << "Discriminant = " << disc << endl;
	      rdisc = sqrt(I2bigfloat(absdisc));
	      ax = fac1 * sqrt(rdisc);
	      alim=Ifloor(ax);
//	      cout<<"Bound on a = " << alim << endl;
	      for(a=1; a<=alim; a++)
		{
		  a2=a*a;
		  blim=(3*a)/2;
//	          cout<<"a="<<a<<": bound on b = "<<blim<<endl;
		  for(b=-blim; b<=blim; b++)
		    {
		      b2=b*b; b3=b*b2;
		      bigfloat i3a = to_bigfloat(1)/I2bigfloat(3*a);
		      cy=I2bigfloat(b2)*i3a;
		      cx=(I2bigfloat(b2)-fac2*rdisc)*i3a;
		      cmin=Iceil(cx);
		      cmax=Ifloor(cy);
//		      cout<<"a="<<a<<", b="<<b<<": bounds on c: "<<cmin<<","<<cmax<<endl;
		      for(c=cmin; c<=cmax; c++)
			{
			  P = b2-3*a*c;
			  U2 = 4*P*P*P-27*disc*a2;
			  if(isqrt(U2,U))
			    {
			      if(::divides(U-2*b3+9*a*b*c,27*a2,d,r))
				{
				  cout<<disc<<"\t";
				  cubic g(a,b,c,d);
				  cout<<g<<"\t----------->\t";
				  if(neg)
				    g.jc_reduce(m);
				  else
				    g.hess_reduce(m);
				  cout<<g<<endl;
				}
			    }
			}
		    }
		}
	    }
    }
  cout<<endl;
}
