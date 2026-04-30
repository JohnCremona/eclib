// File RREF.CC: functions for working with number fields for Hecke eigenvalues
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


#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#define uscalar mp_limb_t // unsigned long

// Convert a mat_m to a flint fmpz_mat

void flint_mat_from_mmat(fmpz_mat_t& A, const mat_m& M)
{
  long nr=M.nrows(), nc=M.ncols();
  long i, j;

  // create flint matrix copy of M:
  fmpz_mat_init(A, nr, nc);
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++)
      fmpz_set_si(fmpz_mat_entry(A,i,j), I2long(M(i+1,j+1)));
}

mat_m mmat_from_flint_mat(const fmpz_mat_t& A)
{
  long nr=fmpz_mat_nrows(A), nc=fmpz_mat_ncols(A);

  // create matrix copy of A:
  mat_m M(nr, nc);
  long i, j;
  for(i=0; i<nr; i++)
    for(j=0; j<nc; j++)
      M(i+1,j+1) = fmpz_get_si(fmpz_mat_entry(A,i,j));
  return M;
}

