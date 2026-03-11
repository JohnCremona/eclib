// modulus.cc: implementations of modulus handling class and functions
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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
template<> INT default_modulus<INT>() { return modulus_factory_INT.get_modulus();}

const int default_modulus_int(1073741789);
const long default_modulus_long(1073741789);
const string bigprime("1000000000000000000000000000057");
const ZZ default_modulus_ZZ(to_ZZ(bigprime.c_str()));
const INT default_modulus_INT(bigprime);

template<> modulus_factory<int>::modulus_factory()    :modulus(default_modulus_int) { ; }
template<> modulus_factory<long>::modulus_factory()   :modulus(default_modulus_long) { ; }
template<> modulus_factory<ZZ>::modulus_factory()     :modulus(default_modulus_ZZ) { ; }
template<> modulus_factory<INT>::modulus_factory()    :modulus(default_modulus_INT) { ; }

template<class T> T modulus_factory<T>::set_modulus(const T& new_modulus)
{
  T old_modulus = modulus;
  modulus = new_modulus;
  return old_modulus;
}

modulus_factory<int> modulus_factory_int;
modulus_factory<long> modulus_factory_long;
modulus_factory<ZZ> modulus_factory_ZZ;
modulus_factory<INT> modulus_factory_INT;

template class modulus_factory<int>;
template class modulus_factory<long>;
template class modulus_factory<ZZ>;
template class modulus_factory<INT>;

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

template<>
INT set_default_modulus<INT>(const INT& new_modulus)
{
  return modulus_factory_INT.set_modulus(new_modulus);
}
