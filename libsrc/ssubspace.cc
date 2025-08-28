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

// Instantiate ssubspace template classes for T=int, long, bigint

template class ssubZspace<int>;
template class ssubZspace<long>;
template class ssubZspace<bigint>;

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
 
template<class T>
sZmat<T> restrict_mat(const sZmat<T>& m, const ssubZspace<T>& s)
{ //cout<<"In restrict_mat(), s.modulus = "<<s.modulus<<endl;
  return mult_mod_p(m.select_rows(pivots(s)),basis(s),s.modulus);
}

template<class T>
ssubZspace<T> eigenspace(const sZmat<T>& sm, T lambda, T m)
{
  sZmat<T> m1 = sm; m1.sub_mod_p(lambda, m);
  return kernel(m1, m);
}
 
template<class T>
ssubZspace<T> subeigenspace(const sZmat<T>& sm, T l, const ssubZspace<T>& s, T m)
{
  return combine(s,eigenspace(restrict_mat(sm,s), l, m));
}

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

// Instantiate template functions for T=int
template ssubZspace<int> make1d<int>(const Zvec<int>& bas, int&piv, int m);
template sZmat<int> restrict_mat<int>(const sZmat<int>& m, const ssubZspace<int>& s);
template int dim<int>(const ssubZspace<int>& s);
template Zvec<int> pivots<int>(const ssubZspace<int>& s);
template sZmat<int> basis<int>(const ssubZspace<int>& s);
template ssubZspace<int> combine<int>(const ssubZspace<int>& s1, const ssubZspace<int>& s2);
template ssubZspace<int> eigenspace<int>(const sZmat<int>& sm, int lambda, int m);
template ssubZspace<int> subeigenspace<int>(const sZmat<int>& sm, int l, const ssubZspace<int>& s, int m);

// Instantiate template functions for T=long
template ssubZspace<long> make1d<long>(const Zvec<long>& bas, long&piv, long m);
template sZmat<long> restrict_mat<long>(const sZmat<long>& m, const ssubZspace<long>& s);
template int dim<long>(const ssubZspace<long>& s);
template Zvec<int> pivots<long>(const ssubZspace<long>& s);
template sZmat<long> basis<long>(const ssubZspace<long>& s);
template ssubZspace<long> combine<long>(const ssubZspace<long>& s1, const ssubZspace<long>& s2);
template ssubZspace<long> eigenspace<long>(const sZmat<long>& sm, long lambda, long m);
template ssubZspace<long> subeigenspace<long>(const sZmat<long>& sm, long l, const ssubZspace<long>& s, long m);

// Instantiate template functions for T=bigint
template ssubZspace<bigint> make1d<bigint>(const Zvec<bigint>& bas, bigint&piv, bigint m);
template sZmat<bigint> restrict_mat<bigint>(const sZmat<bigint>& m, const ssubZspace<bigint>& s);
template int dim<bigint>(const ssubZspace<bigint>& s);
template Zvec<int> pivots<bigint>(const ssubZspace<bigint>& s);
template sZmat<bigint> basis<bigint>(const ssubZspace<bigint>& s);
template ssubZspace<bigint> combine<bigint>(const ssubZspace<bigint>& s1, const ssubZspace<bigint>& s2);
template ssubZspace<bigint> eigenspace<bigint>(const sZmat<bigint>& sm, bigint lambda, bigint m);
template ssubZspace<bigint> subeigenspace<bigint>(const sZmat<bigint>& sm, bigint l, const ssubZspace<bigint>& s, bigint m);
