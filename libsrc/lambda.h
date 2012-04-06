// lambda.h   Declarations of functions which compute Silverman's
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
 //            finite set Lambda_bad for a curve


// N.B. (1) Uses my height normalization, double S's.
// (3) Uses the local height normalization WITHOUT the log|Delta|
// (2) Intended for use in computing Heegner points (not yet implemented)

vector<bigfloat> lambda_bad_1(const bigint& p, long kcode, long npd, long& nlambda);
vector<bigfloat> lambda_bad(const CurveRed& C, long& nlambda, int verbose=0);

// kcode is Kodaira Code with usual coding 10*n for I_n, 10*m+1 for I*_m, 
// 1,2,3,4,5,6,7 for I,II,III,IV,IV*.III*,II* respectively

// caller must delete returned array; first nlambda entries are filled.

int make_point_from_x(Curvedata* CD, const bigint& a, const bigint& d, Point* P);
int make_point_from_x(Curvedata* CD, const bigfloat& x, long maxdd, Point* P);
int make_point_from_x_and_ht(Curvedata* CD, vector<bigfloat> lambdas, const bigfloat& xp, const bigfloat& ht, Point* P);

