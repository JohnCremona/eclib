// ct.h: Functions for Cassels-Tate paring between two quartics equivariants.
// Implementation based On Binary Quartics and the Cassels-Tate Pairing (Tom Fisher) //
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _CT_H_
#define _CT_H_
#include <eclib/interface.h>

void hom_inv(const bigint& ii, const bigint& jj, vector<vector<bigint>>& qelsgens);
	// ii, jj are invariants of minimal elliptic curve associate to quartics on qelsgens;
	// if necessary hom_inv doubleup for that all quartics has the same invariant

void hom_CT(const bigint& ii, const bigint& jj, vector<vector<bigint>>& qgens);
	// retires of qelsgens quartics with non-trivial pairing CT;

#endif
