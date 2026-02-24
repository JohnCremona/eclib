// FILE FIELDTEST.CC: test program for classes Field, FieldElement, FieldIso

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
      cin >> &F;
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
  cout << "min poly = " << str(F.poly()) << endl;
  cout << "generator name = " << F.get_var() << endl;
  cout << "generator = " << F.gen() << endl;
  cout << "0 in F is " << F.zero() << endl;
  cout << "1 in F is " << F.one() << endl;
  cout << "-1 in F is " << F.minus_one() << endl;
  cout << "2 in F is " << F.two() << endl;
  cout << "-2 in F is " << F.minus_two() << endl;
  cout << "Applying polredabs..." << endl;
  FieldIso iso = F.reduction_isomorphism(F.get_var()+"0");
  cout << iso << endl;
  cout << "Reduced field is ";
  iso.codom()->display();
  cout << endl;
}
