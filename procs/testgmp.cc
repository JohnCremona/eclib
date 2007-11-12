// testgmp.cc: test program for gmp library

#include	<gmp.h>

int main ()
{
	mpz_t	x, y, z;

	mpz_init(x);
	mpz_init(y);
	mpz_init(z);

	mpz_gcd(z, x, y);

	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(z);
	return 0;
}
