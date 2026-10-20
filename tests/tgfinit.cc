// Check default finite-field construction before the context cache exists.
// Copyright 2026 Dima Pasechnik
// SPDX-License-Identifier: GPL-2.0-or-later

#include "eclib/gf.h"

// Reproduce the initialization order possible for ffmodq::Fq in a shared
// library. Including the implementation below places the cache after this
// object in one translation unit, so the test does not depend on link order.
galois_field field_before_cache;

#include "../libsrc/gf.cc"

int main()
{
  if (field_before_cache.characteristic() != to_ZZ(2)
      || ZZ_p::modulus() != to_ZZ(2))
    return 1;

  // The explicit constructor must still initialize and restore cached fields.
  galois_field five(to_ZZ(5));
  if (five.characteristic() != to_ZZ(5) || ZZ_p::modulus() != to_ZZ(5))
    return 1;
  galois_field seven(to_ZZ(7));
  if (seven.characteristic() != to_ZZ(7) || ZZ_p::modulus() != to_ZZ(7))
    return 1;
  galois_field five_again(to_ZZ(5));
  return five_again.characteristic() != to_ZZ(5)
         || ZZ_p::modulus() != to_ZZ(5);
}
