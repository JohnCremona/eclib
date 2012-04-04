// file fixc6.h: declaration of class for "fixing" large c6 values
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

// Background: for curves with large c6, the value cannot be computed
// with sufficient precision without using multiprecision floating
// point arithmetic, which slows down all the other cases
// unnecessarily.  To avoid this, after first computing the "correct"
// c6 value once in multiprecision, we add that value to a table
// indexed by (level, form#), and look up in the table when we need
// the value.

// The fixc6 class has two static data members, of type map<
// pair<long,int>, bigint> such that an entry ((N,i),c6) or ((N,i),c4)
// says that the c6 or c4 value for form i at level N is c6 or c4.  Of
// course, for most (N,i) pairs this is blank -- and we must avoid
// inserting wrong dummy entries of the form ((N,i),0).

// April 2005: added facility for fixing c4 as well as c6, but the
// class name is unchanged

class fixc6 {

  static map< pair<long,int>, bigint > fixc4table;
  static map< pair<long,int>, bigint > fixc6table;

public:

  fixc6();  // global initializer, see fixc6.cc
  void operator()(long N, int i, bigint& c4, bigint& c6);  
// look up value, changes c4 and/or c6 if there's an entry in the table,
// otherwise leaves unchanged

  };

extern fixc6 c4c6fixer;  // the one and only instance of the class

