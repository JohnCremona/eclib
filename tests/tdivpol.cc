// tdivpol.cc -- test for division poly functions in divpol.h/cc 
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
//

#include <eclib/curve.h>

#include <eclib/polys.h>
#include <eclib/divpol.h>

int main()
{
  //  set_precision("Enter precision in bits");
  initprimes("PRIMES",0);

  Curve E(BIGINT(0),BIGINT(0),BIGINT(1),BIGINT(-7),BIGINT(6));

  Curvedata C(E);
  cout << "Curve " << E << endl;

  cout<<"Division Poly (2) = \t"<<makepdivpol(&C,2)<<endl;
  int i;
  for(i=3; i<12; i+=2)
    cout<<"Division Poly ("<<i<<") = \t"<<makepdivpol(&C,i)<<endl;


}


//end of file tdivpol.cc





