/* parifact.cc: interface to libpari functions */
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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

#include <eclib/parifact.h>
#include <eclib/interface.h>
#include <pari/pari.h>
#include <iostream>
#include <sstream>

//#define DEBUG_GPFACT

namespace PARI
{
	using std::string;
	using std::cout;
	using std::stringstream;

	string Z_to_string(const bigint& n)
	{
		// convert NTL/ZZ to string;
		stringstream ssa;
		ssa << n;
		return ssa.str().c_str();
	}

	GEN Z_to_pari(const bigint& n)
	{
		// convert NTL/ZZ to t_INT libpari;
		GEN s = strtoi(Z_to_string(n*sign(n)).c_str());
		return gmulgs(s, sign(n));
	}

	vector<bigint> factor(const bigint& n)
	{
		//factor bigint via libpari(eg pseudoprimes);
		if (!avma)
		{
    		long pari_size = strtol(getenv_with_default("PARI_SIZE", "1000000000").c_str(), NULL, 0);
    		if (pari_size==0) // e.g. syntax error in the environment variable PARI_SIZE
      		pari_size = 1000000000;
#ifdef DEBUG_GPFACT
    std::cout<<"calling pari_init with pari_size = "<<pari_size<<endl;
#endif
    		// the first parameter is the maximum stack size in bytes
    		// the second parameter is the maximum precomputed prime
    		pari_init(pari_size, 1000000);
  		}
#ifdef DEBUG_GPFACT
  std::cout<<"factor called with "<<n<<endl;
#endif
		pari_sp av=avma;  // store pari stack pointer

		GEN x = Z_to_pari(n);
		setsigne(x,1);
		x = gel(Z_factor(x),1);
		settyp(x,t_VEC);

		int i;
		vector<bigint> res;
		for (i=1; i<lg(x); i++)
		{
			res.push_back(BIGINT(GENtostr_raw(gel(x,i))));
		}
#ifdef DEBUG_GPFACT
  std::cout<<"factor returns "<<res<<endl;
#endif
		avma=av;         // restore pari stackpointer

		return res;
	}


	int is_prime(const bigint& n)
	{
		//check if is prime bigint via libpari;
		if (!avma)
		{
			long pari_size = strtol(getenv_with_default("PARI_SIZE", "1000000000").c_str(), NULL, 0);
#ifdef DEBUG_GPFACT
    std::cout<<"calling pari_init with pari_size = "<<pari_size<<endl;
#endif
    		if (pari_size==0) // e.g. syntax error in the environment variable PARI_SIZE
      		pari_size = 1000000000;
    		// the first parameter is the maximum stack size in bytes
    		// the second parameter is the maximum precomputed prime
    		pari_init(pari_size, 1000000);
  		}
  		pari_sp av=avma;  // store pari stack pointer

#ifdef DEBUG_GPFACT
  std::cout<<"is_prime called with "<<n<<"..."<<flush;
#endif

		int ans;
		ans = (gisprime(Z_to_pari(n),1),1);

#ifdef DEBUG_GPFACT
  std::cout<<"and returns "<<ans<<std::endl;
#endif

		avma=av;         // restore pari stackpointer

		return ans;
	}

	vector<bigint> solve_conic(const vector<bigint>& L)
	{
		// Input: Vector of bigint of lenght 6, that represent the coeff of x^2, x*y, x*z, y^2, y*z, z^2, respectivament;
		// Output: integer vector of lenght 3 that is solution of quadric;

		if (!avma)
		{
    		long pari_size = strtol(getenv_with_default("PARI_SIZE", "1000000000").c_str(), NULL, 0);
    		if (pari_size==0) // e.g. syntax error in the environment variable PARI_SIZE
      		pari_size = 1000000000;
#ifdef DEBUG_GPFACT
    std::cout<<"calling pari_init with pari_size = "<<pari_size<<endl;
#endif
    		// the first parameter is the maximum stack size in bytes
    		// the second parameter is the maximum precomputed prime
    		pari_init(pari_size, 1000000);
  		}
#ifdef DEBUG_GPFACT
  std::cout<<"solve_conic called with "<<L<<"..."<<flush;
#endif
		pari_sp av=avma;  // store pari stack pointer

		GEN M = zeromatcopy(3, 3);
		gcoeff(M,1,1) = Z_to_pari(L[0]);
		gcoeff(M,1,2) = gdiv(Z_to_pari(L[1]), gen_2);
		gcoeff(M,1,3) = gdiv(Z_to_pari(L[2]), gen_2);
		gcoeff(M,2,1) = gdiv(Z_to_pari(L[1]), gen_2);
		gcoeff(M,2,2) = Z_to_pari(L[3]);
		gcoeff(M,2,3) = gdiv(Z_to_pari(L[4]), gen_2);
		gcoeff(M,3,1) = gdiv(Z_to_pari(L[2]), gen_2);
		gcoeff(M,3,2) = gdiv(Z_to_pari(L[4]), gen_2);
		gcoeff(M,3,3) = Z_to_pari(L[5]);
		GEN sol = gtomat(qfsolve(M));

		vector<bigint> res(3);
		int i;
		GEN den;
		den = gen_1;

		for (i=1; i<4; i++)
		{
			if (!(is_intreal_t(typ(gcoeff(sol,i,1)))))
			{
				den = glcm0(den, Q_denom(gcoeff(sol,i,1)));
			}
		}
		for (i=1; i<4; i++)
		{
			res[i-1] = BIGINT(GENtostr_raw(gmul(gcoeff(sol,i,1), den)));
		}
#ifdef DEBUG_GPFACT
  std::cout<<"and returns "<<res<<std::endl;
#endif
		avma=av;         // restore pari stackpointer

		return res;
	}

	vector<vector<bigint>> param_conic(const vector<bigint>& L, const vector<bigint>& sol)
	{
		// Input: Vector of bigint of lenght 6, that represent the coeff of x^2, x*y, x*z, y^2, y*z, z^2, respectivament;
		//		   and sol for conic associated with L;
		// Output: Vector representing the parametrization for solutions of L;

		if (!avma)
		{
    		long pari_size = strtol(getenv_with_default("PARI_SIZE", "1000000000").c_str(), NULL, 0);
    		if (pari_size==0) // e.g. syntax error in the environment variable PARI_SIZE
      		pari_size = 1000000000;
#ifdef DEBUG_GPFACT
    std::cout<<"calling pari_init with pari_size = "<<pari_size<<endl;
#endif
    		// the first parameter is the maximum stack size in bytes
    		// the second parameter is the maximum precomputed prime
    		pari_init(pari_size, 1000000);
  		}
#ifdef DEBUG_GPFACT
  std::cout<<"solve_conic called with "<<L<<"..."<<flush;
#endif
		pari_sp av=avma;  // store pari stack pointer

		GEN M = zeromatcopy(3, 3);
		gcoeff(M,1,1) = Z_to_pari(L[0]);
		gcoeff(M,1,2) = gdiv(Z_to_pari(L[1]), gen_2);
		gcoeff(M,1,3) = gdiv(Z_to_pari(L[2]), gen_2);
		gcoeff(M,2,1) = gdiv(Z_to_pari(L[1]), gen_2);
		gcoeff(M,2,2) = Z_to_pari(L[3]);
		gcoeff(M,2,3) = gdiv(Z_to_pari(L[4]), gen_2);
		gcoeff(M,3,1) = gdiv(Z_to_pari(L[2]), gen_2);
		gcoeff(M,3,2) = gdiv(Z_to_pari(L[4]), gen_2);
		gcoeff(M,3,3) = Z_to_pari(L[5]);

		GEN g_sol = zerocol(3);
		gel(g_sol,1) = Z_to_pari(sol[0]);
		gel(g_sol,2) = Z_to_pari(sol[1]);
		gel(g_sol,3) = Z_to_pari(sol[2]);

		GEN param = qfparam(M,g_sol,1);

		vector<vector<bigint>> res(3,vector<bigint> (3,BIGINT(0)));
		res[0][0] = BIGINT(GENtostr_raw(gcoeff(param,1,1)));
		res[0][1] = BIGINT(GENtostr_raw(gcoeff(param,1,2)));
		res[0][2] = BIGINT(GENtostr_raw(gcoeff(param,1,3)));
		res[1][0] = BIGINT(GENtostr_raw(gcoeff(param,2,1)));
		res[1][1] = BIGINT(GENtostr_raw(gcoeff(param,2,2)));
		res[1][2] = BIGINT(GENtostr_raw(gcoeff(param,2,3)));
		res[2][0] = BIGINT(GENtostr_raw(gcoeff(param,3,1)));
		res[2][1] = BIGINT(GENtostr_raw(gcoeff(param,3,2)));
		res[2][2] = BIGINT(GENtostr_raw(gcoeff(param,3,3)));

#ifdef DEBUG_GPFACT
  std::cout<<"and returns "<<res<<std::endl;
#endif
		avma=av;         // restore pari stackpointer

		return res;
	}
}
