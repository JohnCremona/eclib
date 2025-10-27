// modulus.h: declarations of modulus handling class and functions
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
 
#if     !defined(_ECLIB_MODULUS_H)
#define _ECLIB_MODULUS_H      1       //flags that this file has been included

#include "marith.h"

template<class T> class modulus_factory {
public:
  modulus_factory();
  T set_modulus(const T& new_modulus); // returns old modulus
  T get_modulus() const {return modulus;}
private:
  T modulus;
};

extern modulus_factory<int> modulus_factory_int;
extern modulus_factory<long> modulus_factory_long;
extern modulus_factory<ZZ> modulus_factory_ZZ;

template<class T> T default_modulus();
template<class T> T set_default_modulus(const T& new_modulus); // returns old modulus

#endif
