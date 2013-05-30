// file unimod.cc: implementation of functions in unimod.h
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
 
#include <eclib/interface.h>
#include <eclib/marith.h>
#include <eclib/unimod.h>

unimod operator*(const unimod& a, const unimod& b)
{
  return unimod(
		a.m11*b.m11 + a.m12*b.m21,
		a.m11*b.m12 + a.m12*b.m22,
		a.m21*b.m11 + a.m22*b.m21,
		a.m21*b.m12 + a.m22*b.m22
		);
}

scaled_unimod operator*(const scaled_unimod& a, const scaled_unimod& b)
{
  return scaled_unimod(
		       a.m00*b.m00,
		       a.m11*b.m11 + a.m12*b.m21,
		       a.m11*b.m12 + a.m12*b.m22,
		       a.m21*b.m11 + a.m22*b.m21,
		       a.m21*b.m12 + a.m22*b.m22
		       );
}
