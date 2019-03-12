// mlocsol.h: declaration of functions for local solubility of quartics
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

#if     !defined(_ECLIB_MLOCSOL_H)
#define _ECLIB_MLOCSOL_H      1       //flags that this file has been included

#include <eclib/mquartic.h>
#include <eclib/compproc.h>
#include <eclib/quadratic.h>


static vector<bigint> v0_;
static vector< vector<bigint> > v00_;
static quadratic q0_=quadratic(BIGINT(0), BIGINT(0), BIGINT(0));

// Checks for solublility in Qp:
// If return 1 and s_!=0, then add at xplist one solution (xp, zp) for 
//                             y^2=g(x,z) with val(p, q(x,z))<=precision of (xp,zp)
int qpsoluble(const quartic& g, const bigint& p, int s_=0, vector< vector<bigint> >& xplist=v00_,
				const quadratic& q=q0_);
int qpsoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
	      		const bigint& e, const bigint& p, int s_=0, vector< vector<bigint> >& xplist=v00_,
				const quadratic& q=q0_);
int qpsoluble(const bigint& a, const bigint& c, const bigint& e, const bigint& p, int s_=0, 
				vector< vector<bigint> >& xplist=v00_, const quadratic& q=q0_);
// latter assumes b=d=0


// Check R-solubility, if it's and s_!=0, then add one solution at xplist with q(x,z)!=0;
int Rsoluble(const quartic& g,int s_=0, vector< vector<bigint> >& xplist=v00_, 
				const quadratic& q=q0_);
int Rsoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
	      		const bigint& e, int s_=0, vector< vector<bigint> >& xplist=v00_, 
				const quadratic& q=q0_);


// Checks for local solubility in Qp for all p in plist;
// If it's soluble and s_!=0 then xplist[0] is solution with q(x,z)!=0 and xplist[i] is solution 
//     in Qp, where p=plist[i-1], with val(p, q(x,z))<precision of (xp,zp);
// If not, badp will hold the first p for which NOT soluble in Qp:
int locallysoluble(const quartic& g, const vector<bigint>& plist, bigint& badp, int s_=0,
                   vector< vector<bigint> >& xplist=v00_, const quadratic& q=q0_);
int locallysoluble(const bigint& a, const bigint& b, const bigint& c, const bigint& d,
					 const bigint& e, const vector<bigint>& plist, bigint& badp, int s_=0, 
					 vector< vector<bigint> >& xplist=v00_, const quadratic& q=q0_);
int locallysoluble(const bigint& a, const bigint& c, const bigint& e,
		           const vector<bigint>& plist, bigint& badp, int s_=0,
                   vector< vector<bigint> >& xplist=v00_, const quadratic& q=q0_);
// latter assumes b=d=0


/* Samir Siksek's Local Solubility Test for odd p */

int local_sol(const bigint& p,bigint *c, int verbose=0);

 // Checks for solublility in Qp
int new_qpsoluble(const quartic& g, const bigint& p, int verbose=0);
int new_qpsoluble(const bigint& a, const bigint& b, const bigint& c,
		  const bigint& d, const bigint& e,
		  const bigint& p, int verbose=0);
int new_qpsoluble(const bigint& a, const bigint& c, const bigint& e,
		  const bigint& p, int verbose=0);

int new_zpsol(const bigint& a,const bigint& b,const bigint& c,const bigint& d,
	      const bigint& e, const bigint& p, int verbose=0);

#endif
