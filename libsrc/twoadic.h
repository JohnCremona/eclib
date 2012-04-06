// twoadic.h: declarations of functions for existence of 2-adic points
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
 

// Case 1 is A=0, B=0, x=2 (mod 4)
// Case 2 is A=1, B=2, x=1 (mod 4)

// Here A = -27*I, B = -27*J

// The following macros are due to Michael Stoll:

#define HIGH 0x7FFF
/* 2-adic valuation function */
#define val1(a,s) (((a)&(0x1<<(s))) ? (s) : (s)+1)
#define val2(a,s) (((a)&(0x3<<(s))) ? val1(a,s) : val1(a,(s)+2))
#define val4(a,s) (((a)&(0xF<<(s))) ? val2(a,s) : val2(a,(s)+4))
#define val8(a,s) (((a)&(0xFF<<(s))) ? val4(a,s) : val4(a,(s)+8))
#define val16(a) (((a)&0xFFFF) ? val8(a,0) : val8(a,16))
inline long val(long a){long r = (a) ? val16(a) : HIGH;  return r;}
inline long val(bigint a){long r = (a==0) ? HIGH: val(2,a);  return r;}

// try1(poly), with poly a deg 3 polynomial in x, determines if 
// there is a 2-adic integer a such that poly(a) is a square    
// in Q_2. Returns 1 if successful, 0 otherwise.                

// These were originally used to determine the index in the ambiguous
// cases, but were much slower than the non-recursive functions,
// called case1() and case2()

// In all cases the return value is 0 if the index is 1 
//                              and 1 if the index is 2
// so add 1 to get the index, which is also the "number of I,J pairs"


long try1(long poly[4]);
long try1(bigint poly[4]);
long case1(long a, long b); // A=4a, B=4b
long case2(long a, long b); // A=4a+1, B=4b+2

long case1(bigint a, bigint b); // A=4a, B=4b
long case2(bigint a, bigint b); // A=4a+1, B=4b+2

