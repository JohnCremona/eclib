// FILE FIELDTEST.CC: test program for classes Field, FieldElement, FieldIso
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

#include <cassert>
#include "eclib/polys.h"
#include "eclib/field.h"
#include <eclib/pari_init.h>
#include "eclib/polred.h"

void test_field(const Field& F); // general test for any field
void test_mqfield(const Field& F); // test of multiquadratic extensions

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
      if (F.degree() < 10)
        test_mqfield(F);
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
  for (auto n: {0,1,-1,2,-2})
    cout << n << " in F is " << F(n) << endl;

  cout << "Applying polred..." << endl;
  Field Fred;
  FieldIso iso = F.reduction_isomorphism(F.get_var()+"0", Fred, 1);
  cout << iso << endl;
  cout << "Reduced field is ";
  Fred.display();
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
  cout << "Automorphism mapping generator " << v << " to " << b << "..."<<endl;
  Field Qb;
  FieldIso aut(F.change_generator(b, Qb));
  cout << aut << endl;

  b = a;
  FieldElement r;
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
  Field Frootb;
  FieldIso emb = F.sqrt_embedding(b, "sqrt(r)", Frootb, r, 1); // reduce=1 for polredabs
  cout << emb << endl;
  FieldElement eb = emb(b);
  cout << "The image of b is " << eb << " with sqrt("<<eb<<") = " << r << endl;
  assert (r*r==eb);
  cout << endl;

  cout << "Testing nfinit interface for reduced polynomial " << str(Fred.poly()) << ":\n";
  ZZ ind;
  vector<Qvec> zbasis_coords;
  mat_m bcm;
  nfinit(Fred.poly(), ind, zbasis_coords, bcm);
  cout << "Equation order has index " << ind << " in maximal order " << endl;
  vector<FieldElement> Zbasis;
  for (int i=0; i<d; i++)
    Zbasis.push_back(FieldElement(Fred, zbasis_coords[i]));
  cout << "Integral basis: " << Zbasis << " // random" << endl;
  cout << "Matrix of coefficients of powers w.r.t. integral basis:\n";
  output_flat_matrix(bcm);
  cout << " // random" << endl;
  cout << endl;
  cout << endl;
}

void test_mqfield(const Field& F) // test of multiquadratic extensions
{
  cout << "\nTest of multiquadratic extensions" <<endl;
  vector<int> ri = {-1,2,3};
  size_t nr = ri.size();
  vector<FieldElement> r_list(nr);
  std::transform(ri.begin(), ri.end(), r_list.begin(), [&F](const int& r){return F(r);});
  cout << "Adjoining sqrt(r) to " << F << " for all r in " << ri << endl;
  vector<FieldElement> sqrt_r_list;
  Field K;
  FieldIso emb = F.sqrt_embedding(r_list, "a", K, sqrt_r_list, 1);
  cout << "Extension field is " << K << endl;
  for (int i=0; i<nr; i++)
    {
      cout << "sqrt(" << ri[i] << ") = " << sqrt_r_list[i] << endl;
      assert (sqrt_r_list[i]*sqrt_r_list[i] == emb(r_list[i]));
    }
}
