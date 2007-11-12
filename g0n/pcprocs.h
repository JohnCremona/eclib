// FILE PCPROCS.H: declaration of functions used to compute periods
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

void ratapprox(bigfloat x, long& a, long& b);

int get_curve(long n, long fac, long maxnx, long maxny,
	      const bigfloat& x0, const bigfloat& y0, 
	      long& nx, long& ny, int& type, int detail=0);
