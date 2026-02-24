// FILE POLREDTEST.CC: test program for functions for reducing ZZX polynomials (polredabs) via libpari

#include "eclib/polred.h"

int main()
{
  cout << "Program polredtest: conversions between ZZX and t_POL and reduction of polynomials." << endl;

  eclib_pari_init();
  int d;
  ZZ den;
  while (1)
    {
      cerr << "Enter degree d: ";
      cin >> d;
      if (d<1) exit(0);
      cerr << "Enter " << d+1 << " integers, starting with the leading coefficient: ";
      ZZX f;
      ZZ c;
      for (int i=0; i<=d; i++)
        {
          cin >> c;
          SetCoeff(f, d-i, c);
        }
      cout << "f = " << str(f) << endl;
      GEN P = ZZX_to_t_POL(f);
      pari_printf(" -as a t_POL: %Ps\n", P);
      ZZX g = t_POL_to_ZZX(P, den);
      assert (den==1);
      cout << " -back to ZZX: " << str(g) << endl;
      if (f==g)
        cout << " OK " << endl;
      else
        cout << " *** WRONG *** " << endl;

      if (IsIrreducible(f))
        {
          cout << "Applying polredabs..." << endl;
          ZZX h;
          g = polredabs(f, h, den);
          cout << "... reduced polynomial is g = " << str(g);
          if (f==g)
            cout << " -- no change"<< endl;
          else
            cout << " -- polynomial has been reduced"<< endl;
          if (den==1)
            cout << "A root of f is a = " << str(h, "b") << " where g(b)=0" << endl;
          else
            cout << "A root of f is a = (" << str(h, "b") << ") / " << den << " where g(b)=0" << endl;
        }
    }
  exit(0);
}
