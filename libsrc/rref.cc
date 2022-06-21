#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
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

// mat_m echelon_via_flint(const mat_m& m1, vec_i& pc, vec_i& npc,
//                 long& rk, long& ny, bigint& d)
// {

// }

