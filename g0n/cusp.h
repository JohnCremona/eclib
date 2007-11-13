// FILE CUSP.H
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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

class cusplist {
 private:
    const moddata* N;
    rational *list;
    long number,maxnumber;
    int cuspeq(const rational& c1, const rational& c2) const;
 public:
    cusplist(long n=0, const moddata* iN=0) :N(iN), number(0), maxnumber(n)
      { list=new rational[n];}
    ~cusplist() {delete[] list;}
    long index(const rational& a);
    rational item(long n) const {return list[n];}  //should check n really
    void display() const 
      {for(long i=0; i<number; i++) cout<<i<<"\t"<<list[i]<<endl;}
    long count() const {return number;}
};   
