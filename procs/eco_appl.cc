//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id: eco_appl.cc,v 1.1 2002/08/23 07:35:19 jec Exp $
//
//	Author	: Markus Maurer, Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#include	<LiDIA/eco_prime.h>
#include	<LiDIA/trace_list.h>
#include	<LiDIA/elliptic_curve.h>
#include	<LiDIA/point.h>
#include	<LiDIA/galois_field.h>
#include	<LiDIA/gf_element.h>
#include	<LiDIA/timer.h>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main(int, char**)
{
	bigint p;
	int info;

	//
	// Read characteristic and curve parameters
	//
	info = 1000;

	std::cout << "\n\nECO Point Counting Program (LiDIA Version)";
	std::cout << "\n==========================================\n\n";

	std::cout << "\nCharacteristic : ";
	std::cin >> p;
	p = next_prime(p-1);
	std::cout << "\nChoose next prime, p = "<< p;
	std::cout << " (" << decimal_length(p) << ")" << std::flush;

	bigmod::set_modulus(p);
	galois_field theField(p);

	gf_element a4(theField), a6(theField);

	std::cout << "\n\nCoefficient a4 : ";
	std::cin >> a4;

	std::cout << "\nCoefficient a6 : ";
	std::cin >> a6;

	eco_prime ep(a4, a6);
	int strat;

	std::cout<<"\n\nGeneral strategies: \n";
	std::cout<<"\n(1) Splitting degree computation always";
	std::cout<<"\n(2) Splitting degree computation only for Elkies primes";
	std::cout<<"\n(3) Splitting degree computation only for Atkin primes";
	std::cout<<"\n(4) Splitting degree computation never\n";
	std::cout<<"\nPlease input index of strategy : "; std::cin >> strat;

	switch (strat) {
	case 1:
		ep.set_strategy(eco_prime::COMPUTE_SP_DEGREE);
		break;
	case 2:
		ep.set_strategy(eco_prime::COMPUTE_SP_DEGREE_IF_ELKIES);
		break;
	case 3:
		ep.set_strategy(eco_prime::COMPUTE_SP_DEGREE_IF_ATKIN);
		break;
	case 4:
		ep.set_strategy(eco_prime::DONT_COMPUTE_SP_DEGREE);
		break;
	default:
		std::cerr <<"\n\n No implemented strategy, aborting ...\n\n";
		return 1;
	}

	std::cout<<"\n\nStrategies for finding the eigenvalue: \n";
	std::cout<<"\n(1) Division polynomials with test of X-coordinate";
	std::cout<<"\n(2) Division polynomials with test of Y-coordinate";
	std::cout<<"\n(3) Rational functions with test of X-coordinate";
	std::cout<<"\n(4) Rational functions with test of Y-coordinate";
	std::cout<<"\n(5) Rational functions with test of X-coordinate and table";
	std::cout<<"\n(6) Rational functions with test of Y-coordinate and table";
	std::cout<<"\n(7) Babystep Giantstep for l > 20, Division Polynomials else";
	std::cout<<"\n(8) Funny Babystep Giantstep for l > 20, Division Polynomials else";
	std::cout<<"\n(9) Built-in strategy (Combination Funny BG and Division Polynomials)\n";
	std::cout<<"\nPlease input index of eigenvalue search strategy : "; std::cin >> strat;

	switch (strat) {
	case 1:
		ep.set_ev_search_strategy(eco_prime::EV_DIVISION_POLYNOMIAL |
					  eco_prime::TEST_X_COORDINATE);
		break;
	case 2:
		ep.set_ev_search_strategy(eco_prime::EV_DIVISION_POLYNOMIAL |
					  eco_prime::TEST_Y_COORDINATE);
		break;
	case 3:
		ep.set_ev_search_strategy(eco_prime::EV_RATIONAL_FUNCTION |
					  eco_prime::TEST_X_COORDINATE);
		break;
	case 4:
		ep.set_ev_search_strategy(eco_prime::EV_RATIONAL_FUNCTION |
					  eco_prime::TEST_Y_COORDINATE);
		break;
	case 5:
		ep.set_ev_search_strategy(eco_prime::EV_RATIONAL_FUNCTION_TABLE |
					  eco_prime::TEST_X_COORDINATE);
		break;
	case 6:
		ep.set_ev_search_strategy(eco_prime::EV_RATIONAL_FUNCTION_TABLE |
					  eco_prime::TEST_Y_COORDINATE);
		break;
	case 7:
		ep.set_ev_search_strategy(eco_prime::EV_BABYSTEP_GIANTSTEP);
		break;
	case 8:
		ep.set_ev_search_strategy(eco_prime::EV_FUNNY_BABYSTEP_GIANTSTEP);
		break;
	case 9:
		ep.set_ev_search_strategy(eco_prime::EV_OPTIMAL);
		break;
	default:  std::cerr <<"\n\n No implemented strategy, aborting ...\n\n";
		return 1;
	}

	ep.compute_group_order(info);
	return 0;
}
