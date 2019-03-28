// ct.cc: Functions for Cassels-Tate paring between two quartics equivariants.
// Implementation based On Binary Quartics and the Cassels-Tate Pairing (Tom Fisher) //
///////////////////////////////////////////////////////////////////////////////////////

#include <eclib/ct.h>
#include <eclib/mlocsol.h>
#include <eclib/hilbert.h>
#include <eclib/transform.h>
#include <eclib/parifact.h>
#include <eclib/quadratic.h>


void doubleup(vector<bigint>& v)
{
	v[1]*=2;
	v[2]*=4;
	v[3]*=8;
	v[4]*=16;
	//b*=2; c*=4; d*=8; e*=16; ii*=16; jj*=64; disc*=4096;
}

void hom_inv(const bigint& ii, const bigint& jj, vector<vector<bigint>>& qelsgens)
{
	// Set ii, jj invariants of elliptic curve;
	// doubleup for that all quartics has the same invariant
	long i, j, all_big=0;
	vector<bigint> v;
	for (i=0; i<qelsgens.size(); i++)
	{
		if (all_big)
		{
			v = qelsgens[i];
			if (!((II(v)==ii)&&(JJ(v)==2*jj)))
			{
				doubleup(v);
				qelsgens[i] = v;
			}
		}
		else
		{
			v = qelsgens[i];
			if ((II(v)==ii)&&(JJ(v)==2*jj))
			{
				all_big = 1;
				for (j=0; j<i; j++)
				{
					v = qelsgens[j];
					doubleup(v);
					qelsgens[j] = v;
				}
			}
		}
	}

}

bigint norm_z(const bigint& ii, const bigint& jj, const vector<bigint>& z0)
{
	// norm of z0 on L=Q[x]/(x^3-3*ii*x+jj)
	vector<bigint> L;
	if (z0.size()==2) { L={BIGINT(0),z0[0],z0[1]};}
	else { L=z0;}

	return (pow(jj,2)*pow(L[0], 3) + 3*ii*jj*pow(L[0], 2)*L[1] + 9*pow(ii,2)*pow(L[0],2)*L[2] 
			- jj*pow(L[1], 3) + 3*jj*L[0]*L[1]*L[2] - 3*ii*pow(L[1],2)*L[2] 
			+ 6*ii*L[0]*pow(L[2],2) + pow(L[2],3));
}

vector<bigint> z_inv_mult27(const bigint& ii, const bigint& jj, vector<bigint>& v)
{
	//return 27*z(g), for that z(g) and z(g2)*z_*H(g) has integer coeff;
	vector<bigint> z;
	z = {4*v[0]*27, (3*v[1]*v[1]-8*v[0]*v[2])*27};//(phi^n,...,phi,const)
	if (norm_z(ii, jj, z)!=BIGINT(0)) return z;

	//here norm(z)=0
	unimod m1{BIGINT(0), BIGINT(1), BIGINT(1), BIGINT(2)},
			m2{BIGINT(1), BIGINT(1), BIGINT(-1), BIGINT(0)},
			m3{BIGINT(1), BIGINT(1), BIGINT(1), BIGINT(2)},
			m4{BIGINT(1), BIGINT(1), BIGINT(2), BIGINT(1)},
			m5{BIGINT(-1), BIGINT(-1), BIGINT(2), BIGINT(1)},
			m6{BIGINT(2), BIGINT(1), BIGINT(-1), BIGINT(0)};
	vector<unimod> M;
	M = {m1,m2,m3,m4,m5,m6};
	int i;
	vector<bigint> v0;
	for (i=0; i<6; i++)
	{
		v0 = v;
		apply_transform(v0[0], v0[1], v0[2], v0[3], v0[4], M[i]);
		z = {4*v0[0]*27, (3*v0[1]*v0[1]-8*v0[0]*v0[2])*27};
		if (norm_z(ii, jj, z)!=BIGINT(0))
		{
			v = v0;
			return z;
		}
	}
}

vector<bigint> z_(const bigint& ii, const bigint& jj, const vector<bigint>& z1_, const vector<bigint>& z2_)
{
	//return <a, b, c> that z_ = a*phi^2 + b*phi + c and z3=z1*z2*z_^2 is linear on L*/L*^2;
	bigint z0, z1, z2;
	z2 = z2_[0]*z1_[0];
	z1 = z1_[0]*z2_[1]+z1_[1]*z2_[0];
	z0 = z1_[1]*z2_[1];

	vector<bigint> L, res;
	L = {(3*ii*(z0 + 3*z2*ii) - z1*jj), 2*(3*z1*ii - z2*jj), 2*(z0 + 3*z2*ii), (z0 + 3*z2*ii), 2*z1, z2};
	res = PARI::solve_conic(L);
	// need check if non-null norm
	return res;
}

quadratic gamma_(const bigint& ii, const bigint& jj, const vector<bigint>& qv1, const vector<bigint>& z2, const vector<bigint>& z_)
{
	// set z2 the z invariant of g2, z_ as above and qv a quartic q1 with invariants ii, jj;
	// return the gamma associated;
	const bigint &a=qv1[0], &b=qv1[1], &c=qv1[2], &d=qv1[3], &e=qv1[4];
	const bigint &z21=z2[0], &z20=z2[1];
	const bigint &z_2=z_[0], &z_1=z_[1], &z_0=z_[2];
	quadratic qq;
	qq.set(4*a*ii*z21*z_2 + b*b*z20*z_2 - (8*a*c*z20*z_2)/3 + b*b*z21*z_1 - (8*a*c*z21*z_1)/3
			+ (4*a*z20*z_1)/3 + (4*a*z21*z_0)/3,
			2*b*ii*z21*z_2 + (2*b*c*z20*z_2)/3 - 4*a*d*z20*z_2 + (2*b*c*z21*z_1)/3 - 4*a*d*z21*z_1
			+ (2*b*z20*z_1)/3 + (2*b*z21*z_0)/3,
			(2*c*ii*z21*z_2)/3 + (2*c*c*z20*z_2)/9 - (b*d*z20*z_2)/3 - (8*a*e*z20*z_2)/3 + (2*c*c*z21*z_1)/9
			- (b*d*z21*z_1)/3 - (8*a*e*z21*z_1)/3 + (2*jj*z21*z_2)/9 - (4*ii*z20*z_2)/9 - (4*ii*z21*z_1)/9
			+ (2*c*z20*z_1)/9 + (2*c*z21*z_0)/9 - (2*z20*z_0)/9);
	return qq;

}

vector<bigint> hil_bad_p(const quadratic& gmm, const bigint& a_g2)
{
	//return vector with primes where hilbert symbol of pairing may be non-trivial;
	vector<bigint> hil_bad_p;
	if (!(gmm.disc()==0)) {	hil_bad_p = factor(2*a_g2*gmm.disc());}
	else { hil_bad_p = factor(a_g2*2*(gmm[0]+gmm[2]));}

	return hil_bad_p;
}

long loc_hil(const bigint& a, const bigint& b, const bigint& p)
{
	// standart hilbert symbol;
	if (a*b!=BIGINT(0))
	{
		if (local_hilbert(a,b,p)) {	return -1;}
		else { return 1;}
	}
	return 0;
}

void hom_CT(const bigint& ii, const bigint& jj, vector<vector<bigint>>& qgens)
{
	// retires of qelsgens quartics with non-trivial pairing CT;
	int i=qgens.size()-1, j, rr;
	if (i>0)
	{
		vector<vector<bigint>> qelsgensl=qgens;
		quadratic gmm;
		vector<bigint> z10, z20, z_0, hil_bad_places;
		bigint non;
		vector< vector<bigint> > xplist;
		vector<int> mask(i+1);
		int pairing;

		for (; i+1>0; i--)
		{
			z10 = z_inv_mult27(ii, jj, qelsgensl[i]);
			if (qelsgensl[i][0]==BIGINT(0)) continue;
			for (j=i-1; j+1>0; j--)
			{
				pairing = 1;
				z20 = z_inv_mult27(ii, jj, qelsgensl[j]);
				if (qelsgensl[j][0]==BIGINT(0)) continue;
				z_0 = z_(ii, jj, z10, z20);
				gmm = gamma_(ii, jj, qelsgensl[i], z20, z_0);

				hil_bad_places = hil_bad_p(gmm, qelsgensl[j][0]);
				xplist = v00_;

				locallysoluble(qelsgensl[i][0], qelsgensl[i][1], qelsgensl[i][2], qelsgensl[i][3], qelsgensl[i][4], hil_bad_places, non, 1, xplist, gmm);

				pairing *= loc_hil(gmm(xplist[0][0], xplist[0][1]), qelsgensl[j][0], BIGINT(0));
				for (rr=1; rr<xplist.size(); rr++)
				{
					pairing *= loc_hil(gmm(xplist[rr][0], xplist[rr][1]), qelsgensl[j][0], hil_bad_places[rr-1]);
				}
				if (pairing==-1)
				{
					mask[i]=1;
					mask[j]=1;
				}
			}
		}
		qelsgensl = v00_;
		for (i=0; i<mask.size(); i++)
		{
			if (!mask[i])
			{
				qelsgensl.push_back(qgens[i]);
			}
		}
		qgens = qelsgensl;
	}
}
