// ct.cc: Functions for Cassels-Tate paring between two quartics with same invariants.
// Implementation based On Binary Quartics and the Cassels-Tate Pairing (Tom Fisher) //
///////////////////////////////////////////////////////////////////////////////////////

#include <eclib/ct.h>
#include <eclib/mlocsol.h>
#include <eclib/hilbert.h>
#include <eclib/transform.h>
#include <eclib/parifact.h>
#include <eclib/minim.h>


void doubleup(vector<bigint>& v)
{
	v[1]*=2;
	v[2]*=4;
	v[3]*=8;
	v[4]*=16;
	//b*=2; c*=4; d*=8; e*=16; ii*=16; jj*=64; disc*=4096;
}

void hom_inv(vector<vector<bigint>>& qelsgens)
{
	// doubleup if necessary for that all quartics has the same invariant
	if (1<qelsgens.size())
	{
		long all_big=0;
		bigint ii0=II(qelsgens[0]), jj0=JJ(qelsgens[0]), ii, jj;
		for (long i=1; i<qelsgens.size(); i++)
		{
			ii = II(qelsgens[i]);
			jj = JJ(qelsgens[i]);
			if (!((ii==ii0)&&(jj==jj0)))
			{
				if (all_big) {doubleup(qelsgens[i]);}
				else
				{
					ii0 = max(abs(ii), abs(ii0))*sign(ii);
					jj0 = max(abs(jj), abs(jj0))*sign(jj);
					all_big = 1;
					for (long j=0; j<=i; j++)
					{
						ii = II(qelsgens[j]);
						jj = JJ(qelsgens[j]);
						if (!((ii==ii0)&&(jj==jj0))) {doubleup(qelsgens[j]);}
					}
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

vector<bigint> oper_L(const bigint& ii, const bigint& jj, const vector<bigint>& z1_, const vector<bigint>& z2_)
{
	// product (z1_*z2_) on quotient L;
	vector<bigint> z1, z2;

	if (z1_.size()==2) {z1 = {BIGINT(0), z1_[0], z1_[1]};} else {z1 = z1_;}
	if (z2_.size()==2) {z2 = {BIGINT(0), z2_[0], z2_[1]};} else {z2 = z2_;}

	return vector<bigint> { (3*ii*z1[0]*z2[0] + z1[2]*z2[0] + z1[1]*z2[1] + z1[0]*z2[2]),
							(-jj*z1[0]*z2[0] + 3*ii*z1[1]*z2[0] + 3*ii*z1[0]*z2[1] + z1[2]*z2[1] + z1[1]*z2[2]),
							- jj*z1[1]*z2[0] - jj*z1[0]*z2[1] + z1[2]*z2[2]};
}

vector<bigint> linearizer_z(const bigint& ii, const bigint& jj, const vector<bigint>& z)
{
	// return z_ that (z*z_^2) is a linear representation for z in L/L^2;
	vector<bigint> L;
	L = {(3*ii*(z[2] + 3*z[0]*ii) - z[1]*jj), 2*(3*z[1]*ii - z[0]*jj), 2*(z[2] + 3*z[0]*ii), (z[2] + 3*z[0]*ii), 2*z[1], z[0]};
	return PARI::solve_conic(L);
}

vector<bigint> quartic_of_z(const bigint& ii, const bigint& jj, const vector<bigint>& z)
{
	bigint r;
	if (!isqrt(norm_z(ii, jj, z), r)) {cout << "\nError: norm of z invariant not square. \n"; abort();}

	vector<bigint> g={3*z[0], BIGINT(0), -18*z[1]*z[0], 24*z[0]*r, 36*ii*z[0]*z[0]*z[0]-9*z[1]*z[1]*z[0]};
	scaled_unimod m;
	bigint I=II(g), J=JJ(g);
	vector<bigint> badprimes0 = factor(6*z[0]);

	minim_all(g[0], g[1], g[2], g[3], g[4], I, J, badprimes0, m, 1);

	return g;
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
		if (local_hilbert(a,b,p)) { return -1;}
		else { return 1;}
	}
	return 0;
}

void hom_CT(const bigint& ii, const bigint& jj, vector<vector<bigint>>& qgens)
{
	// retires of qelsgens quartics with non-trivial pairing CT;
	int lg=qgens.size();
	if (lg>1)
	{
		long i, j, l;
		vector<vector<bigint>> qelsgensl=qgens;
		quadratic gmm;
		vector<bigint> z10, z20, z_0, hil_bad_places;
		bigint non=BIGINT(0);
		vector< vector<bigint> > xplist;
		vector<vector<int>> M(lg,vector<int> (lg,0));
		int pairing;

		for (i=lg-1; i+1>0; i--)
		{
			z10 = z_inv_mult27(ii, jj, qelsgensl[i]);
			if (qelsgensl[i][0]==BIGINT(0)) continue;
			for (j=i-1; j+1>0; j--)
			{
				pairing = 1;
				z20 = z_inv_mult27(ii, jj, qelsgensl[j]);
				if (qelsgensl[j][0]==BIGINT(0)) continue;
				z_0 = linearizer_z(ii, jj, oper_L(ii, jj, z10, z20));

				gmm = gamma_(ii, jj, qelsgensl[i], z20, z_0);

				hil_bad_places = hil_bad_p(gmm, qelsgensl[j][0]);
				xplist = v00_;

				locallysoluble(qelsgensl[i][0], qelsgensl[i][1], qelsgensl[i][2], qelsgensl[i][3], qelsgensl[i][4], hil_bad_places, non, 1, xplist, gmm);

				pairing *= loc_hil(gmm(xplist[0][0], xplist[0][1]), qelsgensl[j][0], BIGINT(0));

				for (l=1; l<xplist.size(); l++)
				{
					pairing *= loc_hil(gmm(xplist[l][0], xplist[l][1]), qelsgensl[j][0], hil_bad_places[l-1]);
				}
				if (pairing==-1) { M[i][j]=M[j][i]=1;}
			}
		}

		qgens=qelsgensl;
		qelsgensl={};
		long ngens=1<<lg, mask_S4=0;

		for (i=1; i<ngens; i++)
		{
			if (i&mask_S4) continue;
			vector<int> vp(lg,0);

			l=-1; j=i;
			do {
				l++;
				if(j&1) {for (long r=0; r<lg; r++) {vp[r]^=M[l][r];}}
			} while (j>>=1);

			if (vp==vector<int> (lg,0))
			{
				if (i==(1<<l)) {mask_S4|=(1<<l); qelsgensl.push_back(qgens[l]);}
				else
				{
					l=-1; j=i;
					vector<bigint> product_z={BIGINT(0), BIGINT(0), BIGINT(1)};
					int need_complete_square=0;

					do {
						l++;
						if(j&1)
						{
							vector<bigint> z_j = z_inv_mult27(ii, jj, qgens[l]);
							need_complete_square^=1;
							product_z = oper_L(ii, jj, product_z, z_j);
						}
					} while (j>>=1);

					if (need_complete_square) {for (long r=0; r<3; r++) {product_z[r]*=BIGINT(3);}}

					vector<bigint> lin_prod=linearizer_z(ii, jj, product_z);
					product_z = oper_L(ii, jj, product_z, oper_L(ii, jj, lin_prod, lin_prod));
					product_z = {product_z[1], product_z[2]}; // remove coeff of phi^2 that here is 0

					qelsgensl.push_back(quartic_of_z(ii, jj, product_z));

					mask_S4|=(1<<l);

					for (long r=0; r<lg; r++) { M[l][r]=M[r][l]=0;}
				}
			}
		}

		//cout << "\n";
		//for (i=0; i<lg; i++) { cout  << M[i] << "\n";}

		hom_inv(qelsgensl);
		qgens = qelsgensl;
	}
}
