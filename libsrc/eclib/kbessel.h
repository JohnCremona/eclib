// kbessel.h: implementation of K-Bessel function for arbitrary real
//            parameter nu. Adapted from the PARI function kbessel.
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
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

#ifndef _ECLIB_KBESSEL_H
#define _ECLIB_KBESSEL_H      1
                           //flags that this file has been included


double kbessel(double nu, double gx, int debug=0);

inline double K0(double x) {return kbessel(0,x);}
inline double K1(double x) {return kbessel(1,x);}

#endif
