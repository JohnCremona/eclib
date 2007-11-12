/* Samir's Local Solubility Test for odd p */

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

