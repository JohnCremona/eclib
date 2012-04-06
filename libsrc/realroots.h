// realroots.h: declarations of funtions for real roots of polynomials
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
 
bigfloat safe_sqrt(const bigfloat& x);
bigfloat cube_root(const bigfloat& x);

// coeff contains deg+1 reals starting with the leading coefficient
// which must be nonzero
// 
// we assume the roots are distinct

vector<bigfloat> realroots( const vector<bigfloat>& coeff );

// As above but only root in the interval [-1,1]

vector<bigfloat> realroots11( const vector<bigfloat>& coeff );



