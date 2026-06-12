// ssubspace.cc: implementation of template class ssubZspace
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2025 John Cremona
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

#include "eclib/linalg.h"

template<class T>
ssubZspace<T>::ssubZspace(int n, T mod)
  :modulus(mod), pivots(Zvec<int>::iota(n)), basis(sZmat<T>::identity_matrix(n))
{}

template<class T>
ssubZspace<T>::ssubZspace(const sZmat<T>& b, const Zvec<int>& p, T mod)
  :modulus(mod), pivots(p), basis(b)
{}

template<class T>
ssubZspace<T>::ssubZspace(const ssubZspace<T>& s)
  :modulus(s.modulus), pivots(s.pivots), basis(s.basis)
{}

// assignment
template<class T>
void ssubZspace<T>::operator=(const ssubZspace<T>& s)
{
  modulus=s.modulus;
  pivots=s.pivots;
  basis=s.basis;
}

// Definitions of nonmember, nonfriend operators and functions:

template<class T>
ssubZspace<T> combine(const ssubZspace<T>& s1, const ssubZspace<T>& s2)
{
  T m = s1.modulus;
  return ssubZspace<T>(mult_mod_p(s1.basis,s2.basis,m),s1.pivots[s2.pivots],m);
}
template ssubZspace<int> combine<int>(const ssubZspace<int>&, const ssubZspace<int>&);
template ssubZspace<long> combine<long>(const ssubZspace<long>&, const ssubZspace<long>&);
template ssubZspace<ZZ> combine<ZZ>(const ssubZspace<ZZ>&, const ssubZspace<ZZ>&);
template ssubZspace<INT> combine<INT>(const ssubZspace<INT>&, const ssubZspace<INT>&);

template<class T>
sZmat<T> restrict_mat(const sZmat<T>& m, const ssubZspace<T>& s)
{
  return mult_mod_p(m.select_rows(s.pivots), s.basis, s.modulus);
}
template sZmat<int> restrict_mat<int>(const sZmat<int>&, const ssubZspace<int>&);
template sZmat<long> restrict_mat<long>(const sZmat<long>&, const ssubZspace<long>&);
template sZmat<ZZ> restrict_mat<ZZ>(const sZmat<ZZ>&, const ssubZspace<ZZ>&);
template sZmat<INT> restrict_mat<INT>(const sZmat<INT>&, const ssubZspace<INT>&);

template<class T>
ssubZspace<T> eigenspace(const sZmat<T>& sm, T lambda, T m)
{
  sZmat<T> m1 = sm; m1.sub_mod_p(lambda, m);
  return kernel(m1, m);
}
template ssubZspace<int> eigenspace<int>(const sZmat<int>&, int, int);
template ssubZspace<long> eigenspace<long>(const sZmat<long>&, long, long);
template ssubZspace<ZZ> eigenspace<ZZ>(const sZmat<ZZ>&, ZZ, ZZ);
template ssubZspace<INT> eigenspace<INT>(const sZmat<INT>&, INT, INT);

template<class T>
ssubZspace<T> subeigenspace(const sZmat<T>& sm, T l, const ssubZspace<T>& s, T m)
{
  return combine(s,eigenspace(restrict_mat(sm,s), l, m));
}
template ssubZspace<int> subeigenspace<int>(const sZmat<int>&, int, const ssubZspace<int>&, int);
template ssubZspace<long> subeigenspace<long>(const sZmat<long>&, long, const ssubZspace<long>&, long);
template ssubZspace<ZZ> subeigenspace<ZZ>(const sZmat<ZZ>&, ZZ, const ssubZspace<ZZ>&, ZZ);
template ssubZspace<INT> subeigenspace<INT>(const sZmat<INT>&, INT, const ssubZspace<INT>&, INT);

template<class T>
ssubZspace<T> make1d(const Zvec<T>& bas, T&piv, T m)
// make a 1-D ssubspace with basis bas
{
  sZmat<T> tbasis(1,dim(bas));
  sZvec<T> sbas(bas);
  tbasis.setrow(1,sbas);
  Zvec<int> pivs(1); // initialised to 0
  pivs[1]=sbas.first_index();
  piv=sbas.elem(pivs[1]);
  return ssubZspace<T>(transpose(tbasis),pivs,m);
}
template ssubZspace<int> make1d<int>(const Zvec<int>&, int&, int);
template ssubZspace<long> make1d<long>(const Zvec<long>&, long&, long);
template ssubZspace<ZZ> make1d<ZZ>(const Zvec<ZZ>&, ZZ&, ZZ);
template ssubZspace<INT> make1d<INT>(const Zvec<INT>&, INT&, INT);

// Instantiate ssubspace template classes for T=int, long, ZZ, INT

template class ssubZspace<int>;
template class ssubZspace<long>;
template class ssubZspace<ZZ>;
template class ssubZspace<INT>;
