// poltest.cc: test program for LiDIA polynomials

#include	<LiDIA/bigint.h>
#include	<LiDIA/polynomial.h>
#include	<LiDIA/gf_polynomial.h>
using namespace LiDIA;
using namespace std;

#define MODULAR

#ifdef MODULAR

int main()
{
  cout<<"Test of multiplication of two polynomials modulo p"<<endl;

	bigint p = 1039;

	Fp_polynomial f,g,h;

	f.set_modulus(p);
	f.set_coefficient(-376,0);
	f.set_coefficient(1,1);
	cout << "f = " << f << endl;
	cout << "deg(f) = " << f.degree() << endl;

	g.set_modulus(p);
	g.set_coefficient(-663,0);
	g.set_coefficient(1,1);
	cout << "g = " << g << endl;
	cout << "deg(g) = " << g.degree() << endl;

	h = f*g;

	cout << "h = f*g = " << h << endl;
	cout << "deg(h) = " << h.degree() << endl;

	return 0;
}
#else
int main()
{
  cout<<"Test of multiplication of two integer polynomials"<<endl;

	polynomial<bigint> f,g,h;

	f.assign_x();
	f-=polynomial<bigint>(376);
	cout << "f = " << f << endl;
	cout << "deg(f) = " << f.degree() << endl;

	g.assign_x();
	g-=polynomial<bigint>(663);
	cout << "g = " << g << endl;
	cout << "deg(g) = " << g.degree() << endl;

	h = f*g;

	cout << "h = f*g = " << h << endl;
	cout << "deg(h) = " << h.degree() << endl;

	return 0;
}
#endif
