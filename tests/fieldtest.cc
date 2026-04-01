// FILE FIELDTEST.CC: test program for classes Field, FieldElement, FieldIso

#include <cassert>
#include "eclib/polys.h"
#include "eclib/field.h"
#include <eclib/pari_init.h>

void test_field(const Field& F);

int main()
{
  cout << "Program fieldtest: test program for classes Field, FieldElement, FieldIso" << endl
       << "\nInput format for a field: either Q or a variable name followed by"
       << "\n a list of coefficients of a monic irreducible defining polynomial (constant term first)"
       << "\n e.g. i [1 0 1]  for Q(i)" << endl;

  eclib_pari_init();
  Field F;
  cout << "Default constructor:\n";
  test_field(F);
  while(1)
    {
      cerr << "Enter a field F (Q to end): ";
      cin >> F;
      if (F.isQ())
        break;
      test_field(F);
    }
  exit(0);
}

void test_field(const Field& F)
{
  cout << "F = " << F << endl;
  F.display();
  cout << "F is Q? " << F.isQ() << endl;
  int d = F.degree();
  cout << "min poly = " << str(F.poly()) << " of degree " << d << endl;
  string v = F.get_var();
  if (d>1) cout << "generator name = " << v << endl;
  FieldElement a = F.gen();
  cout << "generator = " << a << endl;
  cout << "0 in F is " << F.zero() << endl;
  cout << "1 in F is " << F.one() << endl;
  cout << "-1 in F is " << F.minus_one() << endl;
  cout << "2 in F is " << F.two() << endl;
  cout << "-2 in F is " << F.minus_two() << endl;

  cout << "Applying polred..." << endl;
  FieldIso iso = F.reduction_isomorphism(F.get_var()+"0", 1);
  cout << iso << endl;
  cout << "Reduced field is ";
  iso.codom()->display();
  cout << endl;

  FieldElement b = a*a;
  ZZX cp = b.charpoly();
  ZZX mp = b.minpoly();
  int db = deg(mp);
  int e = 2;
  cout << "b = " << a << "^" << e << " = " << b
       << ", with degree " << db << ", char poly = " << str(cp)
       << " and min poly " << str(mp) << endl;
  while (db<d)
    {
      b *= a;
      e += 1;
      cp = b.charpoly();
      mp = b.minpoly();
      db = deg(mp);
      cout << "b = " << a << "^" << e << " = " << b
           << ", with degree " << db << ", char poly = " << str(cp)
           << " and min poly " << str(mp) << endl;
    }
  cout << "Automorphism mapping " << v << " to " << b << "..."<<endl;
  FieldIso aut(F.change_generator(b));
  cout << iso << endl;

  b = a;
  FieldElement r(&F);
  int sq = b.is_square(r);
  cout << "b = " << b << (sq? " is" : " is not") << " a square in F" << endl;
  while (sq)
    {
      cout << "sqrt(b) = " << r << endl;
      assert (r*r==b);
      b += 1;
      sq = b.is_square(r);
      cout << "b = " << b << (sq? " is" : " is not") << " a square in F" << endl;
    }
  cout << "adjoining sqrt(b) to F..."<<endl;
  FieldIso emb = F.sqrt_embedding(b, "r", r, 1); // reduce=1 for polredbest
  cout << emb << endl;
  Field Frootb = *(iso.codom());
  FieldElement eb = emb(b);
  cout << "The image of b is " << eb << " with sqrt("<<eb<<") = " << r << endl;
  assert (r*r==eb);
  cout << endl;
}
