// random.cc: implementation of random functions used in smat testing
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
 
/* codes were copied from the book 
 * Numerical recipes in C,  
 * under "immediate license" agreement,
 * ( see page xvi of the second edition ).
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0( long& idum ) 
{
  long k;
  float ans;

  idum ^= MASK;
  k = (idum)/IQ;
  idum = IA*(idum - k*IQ) - IR*k;
  if( idum < 0 ) idum += IM;
  ans = AM*(idum);
  idum ^= MASK;
  return ans;
}

float ran0( int& idum ) 
{
  int k;
  float ans;

  idum ^= MASK;
  k = (idum)/IQ;
  idum = IA*(idum - k*IQ) - IR*k;
  if( idum < 0 ) idum += IM;
  ans = AM*(idum);
  idum ^= MASK;
  return ans;
}
