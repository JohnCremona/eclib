// FILE MODDATA.H: Declaration of class moddata
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

#include <fstream>
#include "rat.h"
#include "xsplit.h"

class level {
  friend class symb;
  friend class cusplist;
  friend class oldforms;
  friend class eigdata;  
  friend class newform;  
  friend class h1newform;  
  friend class summer;  
//protected:
public:     // got tired of making everything friends...
 long modulus;
 int plusflag; int squarelevel;
 vector<long> plist,dlist,primelist;
 long p0;  // first good prime
 long npdivs,ndivs,sqfac,nap;
 long reduce(long res)const {res%=modulus; return (res<0)?modulus+res:res;}
public:
 level(long n, long neigs=20);
};

class moddata :public level {
protected:
 long phi,psi,nsymb1,nsymb2;
 vector<long> invlist,noninvlist,noninvdlist,dstarts,gcdtable,unitdivlist;
 long code(long res) const {return invlist[reduce(res)];}
public:
  long nsymb;
 moddata(long n);                                //constructor
 void display() const;
 void abort(const string mess) const 
  {
    cout<<"Out of memory ("<<mess<<")\n";  
    ::abort();
  }
 long gcd(long res) const {return gcdtable[reduce(res)];}
 long unitdiv(long res) const {return unitdivlist[reduce(res)];}
};

// Returns oldform filename. Default is "newforms/", but can
// be changed via OF_DIR environment variable.
string of_filename(long n, char c);

// Function to construct a string representing newform filename,
// including a prefix (default "newforms/" but taken from environment
// variable NF_DIR if defined), followed by a single character
// (usually 'x') and the decimal digits of n (positive integer).
string nf_filename(long n, char c);
