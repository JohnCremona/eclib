// desc2.h:  declaration of second descent (via 2-isogeny) procedure
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
 
int desc2(const bigint& c, const bigint& d1, const bigint& d2,
	  const vector<bigint>& plist, const vector<bigint>& supp, const vector<bigint>& bgens,
	  long mask,  double hlim,
	  bigint& x, bigint& y, bigint& z, int verb, int selmer_only=0, int alldesc=0);
// Works on homogeneous space (d1,0,c,0,d2)
// Returns 
//   -1 if it certainly has no points (if no ELS descendents)
//   +1 if it has a point (coordinates returned in x, y, z)
//    0 if undecided (ELS descendents exist but no rational points were found)
// if alldesc==1 it does not stop when it finds one descendent with a point on it,
// but goes on to look at all the others.
