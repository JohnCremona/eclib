// tdivpol.cc -- test for division poly functions in divpol.h/cc 
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
//

#include "curve.h"

#include "polys.h"
#include "divpol.h"

int main()
{
  //  set_precision("Enter number of decimal places");
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





