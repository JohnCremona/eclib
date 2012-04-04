// FILE SYMB.H: Declarations for M-symbols, modular symbols
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

class symb {
 private:
   long c,d;
   const moddata * N;  // needs to be creatable from a const thing
 public:
   symb() {c=d=0; N=0;}
   symb(long ic, long id, const moddata * iN)   {  c=ic; d=id; N=iN;}
   long cee() const        {return c;}
   long dee() const        {return d;}
   long ceered() const        {return N->reduce(c);}
   long deered() const        {return N->reduce(d);}
   long modulus() const       {return N->modulus;}
   int operator==(const symb& s) const 
  {return 0==((xmodmul(c,s.d,N->modulus)-xmodmul(s.c,d,N->modulus))%(N->modulus));}
   int eq(const symb& s) const 
       {return ((c==s.c)&&(d==s.d))||((c==-s.c)&&(d==-s.d));}
   symb normalize() const;
   friend ostream& operator<< (ostream&s, const symb&);
   long orbitlength() const
     {long n=N->modulus, cr=N->reduce(c); cr=cr*cr; 
     return n/N->gcd(cr);}
};

class modsym {
 private:
    rational a,b;
 public:
    modsym()                                {a=rational(0); b=rational(0);}
    modsym(const rational& ra, const rational& rb) {a=ra; b=rb;}
    modsym(const symb&);                        //conversion
    rational alpha() const {return a;}
    rational beta() const {return b;}
    friend ostream& operator<< (ostream& s, const modsym& m);
};

#include <map>

class symblist {
 private:
    symb *list;
    map<pair<long,long>,long> hashtable;
    long num,maxnum;
 public:
  symblist(long n=0);
  ~symblist();
    void add(const symb& s, long start=0);
    long index(const symb& s, long start=0) const;
    symb operator[](long n) const {return list[n];}
    symb item(long n) const;
    void display() const {for(long i=0; i<num; i++) cout<<i<<"\t"<<list[i]<<"\n";}
    long count() const {return num;}
};

class symbdata :public moddata {
 private:
    symblist specials;         // The list of "special" symbols
 public:
    symbdata(long);             // The constructor
    long index2(long c, long d) const;
    long index(const symb& s) const {return index2(s.cee(),s.dee());}
    symb symbol(long i) const;
    void display() const;
    void check() const;
    long rof(long i) const {symb s=symbol(i); return index2(s.dee(), s.cee());}
    long rsof(long i) const {symb s=symbol(i); return index2(-s.cee(),s.dee());}
    long sof(long i) const {symb s=symbol(i); return index2(-s.dee(), s.cee());}
    long tof(long i) const {symb s=symbol(i); long c=s.cee(), d=s.dee(); return index2(c-d, c);} 
};

modsym jumpsymb(symb s1, symb s2);

