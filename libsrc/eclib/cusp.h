// FILE CUSP.H
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

#ifndef _CUSP_H
#define _CUSP_H      1
                           //flags that this file has been included

#include <eclib/moddata.h>
#include <eclib/rat.h>

class cusplist {
 private:
    const moddata* N;
    rational *list;
    long number,maxnumber;
  int cuspeq(const rational& c1, const rational& c2, int plusflag=0) const;
 public:
    cusplist(long n=0, const moddata* iN=0) :N(iN), number(0), maxnumber(n)
      { list=new rational[n];}
    ~cusplist() {delete[] list;}
    long index(const rational& a);
    long index_1(const rational& a);
    long index_2(const rational& a);
    rational item(long n) const {return list[n];}  //should check n really
    void display() const 
      {for(long i=0; i<number; i++) cout<<i<<"\t"<<list[i]<<endl;}
    long count() const {return number;}
};   

#endif
