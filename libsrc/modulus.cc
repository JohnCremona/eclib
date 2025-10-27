// modulus.cc: implementations of modulus handling class and functions
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
 
#include "eclib/modulus.h"

// Instantiate/specialise modulus_factory template classes for T=int, long, ZZ

// PRIME30 is (as defined in xmod.h) 1073741789,  the largest p such that p < 2^30
template<> int default_modulus<int>()       { return modulus_factory_int.get_modulus();}
template<> long default_modulus<long>()     { return modulus_factory_long.get_modulus();}
template<> ZZ default_modulus<ZZ>() { return modulus_factory_ZZ.get_modulus();}

template<> modulus_factory<int>::modulus_factory()    :modulus(1073741789) { ; }
template<> modulus_factory<long>::modulus_factory()   :modulus(1073741789) { ; }
template<> modulus_factory<ZZ>::modulus_factory()
  //  :modulus(to_ZZ("6074000003"))
   :modulus(to_ZZ("1000000000000000000000000000057"))
{ ; }

template<class T> T modulus_factory<T>::set_modulus(const T& new_modulus)
{
  T old_modulus = modulus;
  modulus = new_modulus;
  return old_modulus;
}

modulus_factory<int> modulus_factory_int;
modulus_factory<long> modulus_factory_long;
modulus_factory<ZZ> modulus_factory_ZZ;

template class modulus_factory<int>;
template class modulus_factory<long>;
template class modulus_factory<ZZ>;

template<>
int set_default_modulus<int>(const int& new_modulus)
{
  return modulus_factory_int.set_modulus(new_modulus);
}

template<>
long set_default_modulus<long>(const long& new_modulus)
{
  return modulus_factory_long.set_modulus(new_modulus);
}

template<>
ZZ set_default_modulus<ZZ>(const ZZ& new_modulus)
{
  return modulus_factory_ZZ.set_modulus(new_modulus);
}
