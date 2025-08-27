// FILE CUSP.H
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

#ifndef _ECLIB_CUSP_H
#define _ECLIB_CUSP_H      1
                           //flags that this file has been included

#include "moddata.h"
#include "rat.h"

class cusplist {
 private:
    const moddata* N;
    vector<rational> list;
  int cuspeq(const rational& c1, const rational& c2, int plusflag=0) const;
 public:
  cusplist(long n=0, const moddata* iN=0) :N(iN) {list.reserve(n);}
  long index(const rational& a);
  long index_1(const rational& a);
  long index_2(const rational& a);
  rational item(long n) const {
    if (n<list.size())
      return list[n];
    else
      {
        cerr<< "invalid cusp index "<<n<<": #cusps is "<<list.size()<<endl;
        return rational(0);
      }
  }  //should check n really
  void display() const {
    for(unsigned int i=0; i<list.size(); i++) cout<<i<<"\t"<<list[i]<<endl;
  }
  long count() const {return list.size();}
};

#endif
