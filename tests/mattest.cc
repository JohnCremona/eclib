// mattest.cc: Matrix package test program
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

//#define TIMER  // not defined for automatic tests

#include <cassert>
#include <iostream>
#ifdef TIMER
#include <eclib/timer.h>
#endif
#include <eclib/linalg.h>

const scalar modulus(default_modulus<scalar>());

// test with a size nr x nc matrix of type scalar whose entries are in Mij
void testHNFetc(int nr, int nc, const vector<int>& Mij);

int main(void)
{
#ifdef TIMER
  init_time();
  start_time();
#endif
  cout << "Matrix test program with scalar type " << scalar_type << ".\n\n";

  long i; int r;
  vec_i pc(1),npc(1);
  vec poly(1);
  scalar two(2), three(3);
  cout << "Enter size of a square matrix A: "; cin >> r;
  mat a(r,r);
  cout << "Enter entries of A: "; cin >> a;
  cout << "A = \n" << a << endl;
  cout << "Using A.output(cout): \n";  a.output(cout); cout<<endl;
  cout << "Using A.output_pari(cout): \n";  a.output_pari(cout); cout<<endl;
  cout << "Using A.output_pretty(cout): \n";  a.output_pretty(cout);// cout<<endl;

  cout << "Creating an array of 3 matrices\n";
  vector<mat> matlist(3);
  matlist[0] = a;
  matlist[1] = two*a;
  matlist[2] = three*a;
  cout << " A=\n" << matlist[0] << endl;
  cout << "2A=\n" << matlist[1] << endl;
  cout << "3A=\n" << matlist[2] << endl;

  for (i=1; i<=r; i++)
    cout << "row(A,"<<i<<") = " << a.row(i) << endl;
  cout << "A = \n" << a << endl;
  for (int j=1; j<=r; j++)
    cout << "col(A,"<<j<<") = " << a.col(j) << endl;
  cout << "A = \n" << a << endl;
  cout << "directsum(A,A) = \n" << directsum(a,a) << endl;
  cout << "Is A zero? " << (a.is_zero()?"yes":"no") << endl;

  mat z = a-a;
  cout << "Z=A-A=\n"<<z<<endl;
  cout << "Is Z zero? " << (z.is_zero()?"yes":"no") << endl;
  mat b = a;
  cout << "B = A = \n" << b << endl;
  cout << "B==A?" << (b==a) << endl;
  cout << "B!=A?" << (b!=a) << endl;
  b+=a;
  cout << "after B+:=A, A = \n" << a << "\nand B = \n" << b << endl;
  b-=a;
  cout << "after B-:=A, A = \n" << a << "\nand B = \n" << b << endl;
  b*=two;
  cout << "after B*:=2, A = \n" << a << "\nand B = \n" << b << endl;
  b/=two;
  cout << "after B/:=2, A = \n" << a << "\nand B = \n" << b << endl;
  cout << "A+B=\n" << (a+b) << endl;
  cout << "Now A = \n" << a << "\nand B = \n" << b << endl;
  cout << "A-B=\n" << (a-b) << endl;
  cout << "Now A = \n" << a << "\nand B = \n" << b << endl;
  cout << "A*B=\n" << (a*b) << endl;
  cout << "Now A = \n" << a << "\nand B = \n" << b << endl;
  cout << "-A=\n" << (-a) << endl;
  cout << "Now A = \n" << a << endl;
  cout << "-A=\n" << (-a) << endl;
  cout << "Now A = \n" << a << endl;
  vector<scalar> cp = a.charpoly();
  cout << "char. poly. of A has coefficients " << cp << endl;
  cout << "det(A) = " << a.determinant() << endl;
  mat aug = colcat(a,mat::identity_matrix(r));
  cout << "Augmented matrix = \n" << aug << endl << endl;

  long rk, ny;
  scalar denom;
  // int method=0;
  // cout << "Which echelon method? (0=standard,1=longlong,2=modular) ";
  // cin>>method;
  for (auto method: {0, 2, 3})
    {
  cout << "\nUsing method " << method;
  if(method>0) cout << " (modulus = " << modulus << ")";
  cout << endl;
  mat ref = echelon(aug, pc, npc, rk, ny, denom, method);
  cout << "Echelon matrix = \n" << ref << endl;
  cout << "pivotal columns: " << pc << endl;
  cout << "nonpivotal columns: " << npc << endl;
  cout << "Denom = " << denom << endl;

  for (i=1,rk=0; (i<=r)&&(pc[i]<=r); i++,rk++) ;
  ny = r-rk;
  cout << "Rank = " << rk << endl;
  cout << "Nullity = " << ny << endl;
  if (rk<r) // non-invertible!
    {
      cout << "A is not invertible; rk = " << rk << endl;
    }
  else
    {
      mat ainv = ref.slice(1,r,r+1,r+r);
      cout << "A has inverse ";
      if (denom>1)
        cout << "(1/" << denom << ")*";
      cout << endl << ainv << endl;
      cout << "Check: A.A^(-1) = I ?";
      int ok = (a*ainv == mat::scalar_matrix(r,denom));
      cout << (ok? " True!": " False!") << endl;
      //      assert (ok);
    }
    } // loop over methods

  cout << "Calling ref_via_flint()..." << endl;
  mat R = ref_via_flint(a, modulus);
  cout << " ref_via_flint() returns\n" << R << endl;

  // test conversion to/from FLINT matrices:

  cout << "Testing conversion to and from FLINT matrices: " << flush;
  fmpz_mat_t A;
  flint_mat_from_mat(A,a);
  b = mat_from_flint_mat(A, (scalar)(0));
  if (a == b)
    cout << "OK" << endl;
  else
    {
      cout << "***** not OK: from "; a.output_flat(cout);
      cout << " to "; b.output_flat(cout); cout << endl;
    }

  // Test HNF, SNF and LLL via FLINT
  cout << "Testing HNF, SNF and LLL via FLINT: " << endl;

  // 3x3 example:
  testHNFetc(3, 3, {2,4,4,  -6,6,12,  10,-4,-16});

  // larger examples: entries were generated randomly, then hard-coded
  // in (as otherwise output is random).
#if (SCALAR_OPTION==1) // int
  int d = 5, m=5;
  vector<int> entries = {2, -1, 2, 3, -5, 1, 5, -1, -1, 2, 0, 1, 5, 1, -1, 3, -3, 5, -4, -3,
                         -3, -2, -2, -1, -4};
#endif
#if (SCALAR_OPTION==2) // long
  int d = 10, m=5;
  vector<int> entries = {-3, -3, 2, -1, 3, 3, 0, -1, 4, -5, 4, -3, -3, 3, 2, 3, 2, -5, 2, 3, 2, 0, -3,
                         -2, 5, -2, -1, 0, 2, 1, 1, 1, 0, -4, 3, 1, -5, 5, -3, -1, 0, 5, -3, -1, 0, 1,
                         1, -5, -5, 5, 4, -1, -2, -3, 0, 3, 2, -4, 2, 0, 5, 5, -4, -5, 4, 3, -1, -2, -1,
                         5, -4, 2, 3, 4, 5, 0, -2, 3, -2, 5, 0, -3, 2, -2, -4, 5, 5, -1, -3, 5, -3, 1, 2,
                         5, 2, -2, -2, 4, 0, 5};
#endif
#if (SCALAR_OPTION==3) // ZZ (NTL)
  int d = 20, m=1000;
  vector<int> entries = {60, -862, 248, 707, -8, 991, -265, -423, 299, -84, -298, -268, 253, -96, 422,
                         -590, -948, 106, -831, -910, 23, -134, -664, 534, -224, 235, 53, -151, 522, -84,
                         -599, -287, 54, 271, -42, 464, 193, -332, 289, 847, -124, 839, -135, -38, 445,
                         -83, -566, 45, -828, 588, 653, 495, -254, 988, 649, -766, 691, -551, 567, -966,
                         -177, 859, 469, 137, 363, 852, 801, -905, 645, -467, 928, -418, 823, 889, -154,
                         889, 852, 17, 203, -271, 931, 282, 30, -350, -993, -350, 338, -818, 251, -910,
                         798, -104, 952, -397, 228, -189, 442, -478, 883, 227, -992, 894, 971, 43, 37,
                         446, 745, 756, 128, -894, -24, 172, 83, -927, 464, 361, -121, -993, 292, 338,
                         -799, -692, -200, -327, 933, 598, 121, 674, 398, -306, 999, -499, -364, -762,
                         854, -756, -944, 297, 589, 486, -514, 437, 680, -18, -177, 693, -747, 342, -199,
                         -929, -539, 902, 860, 808, -163, -874, -155, 622, 49, 12, -557, 41, 415, -698,
                         157, -962, -279, -825, -55, -38, 892, -223, 786, 887, -500, -469, 222, 91, -596,
                         688, -13, -347, -324, -696, 875, 327, -935, -23, 661, -816, 89, -457, 142, -831,
                         -351, -643, 605, -382, -17, 355, -900, -602, 650, -338, 432, 487, 229, 268, 388,
                         535, 243, -947, -230, 832, 503, -299, -493, -527, 912, -469, 816, -747, -749, -63,
                         586, -587, -276, -308, 14, 858, 306, -283, -414, 991, -367, -878, -393, 51, 443,
                         660, -410, 851, 368, -104, 778, 768, 9, -600, -175, 658, -205, 997, -297, -412,
                         519, -845, 680, -59, -159, -243, -898, -341, 543, 577, -666, -154, 947, 869,
                         -707, 424, -159, 224, 553, 544, -694, -187, 927, -572, 866, -156, 247, 915, -687,
                         -657, -138, 150, 158, 246, 59, 100, -42, -209, -886, -553, 184, 262, -133, -181,
                         -691, -791, -151, 317, -444, 522, -347, 299, 233, -613, -411, -570, 331, -121,
                         -126, -377, 792, 159, 339, -364, -387, 7, -123, 231, -379, -971, 917, -401, 776,
                         -619, -565, -992, 524, -784, -407, 894, 169, -861, -310, -129, 348, 997, -815,
                         278, -36, 688, -41, -935, 784, -452, 85, -488, 256, 639, -974, 755, 871, -997,
                         -29, -881, -666, -409, 637, -614, -241, -921, -373, -752, 151, -75, 410, -864,
                         -814, -335, 136, -543, 114, -270, -37, 564, 366, 456, -765, 496, -753, 560, 888,
                         467, 147, -484, 541, -645, 103, 906, -504, 169, -345, 555, -235, 297, -57, -981};
#endif
#if (SCALAR_OPTION==4) // INT (FLINT fmpz)
  int d = 10, m=10000;
  vector<int> entries = {-8660, 1243, 8012, -3806, -7123, -7463, 2229, 869, -3386, -2202, 2445, 8371,
                         -3674, -6146, 2364, 7997, 6869, -988, 3247, -6272, -5793, 6501, -6476, 5132,
                         53, 4076, 6564, -5381, -6155, -8797, -3935, -7048, 5337, -4946, -5946, -1346,
                         19, 2033, -3616, 7212, -1112, -2249, -1138, -6880, 755, -2070, -2621, 5965,
                         8969, -3303, 1017, 9459, 605, 8619, 8870, 3551, -1408, 9109, 7957, -7951, -909,
                         988, 8053, -7808, 1200, 2807, -8719, -1012, 3842, 3042, 3962, -5776, 2477, -3680,
                         2942, -8977, 8970, 5122, -4710, 2405, 5578, 6035, 6029, -1487, 1940, -5480, 8383,
                         1604, 1588, -1378, 2175, -8601, 473, 3407, -9877, -1266, 5230, -7472, 8010, 9255};
#endif
  testHNFetc(d, d, entries);  // random_vector(d*d, -m,m, 0)); // 0: not primitive

#ifdef TIMER
  stop_time();
  cout << "cpu time = "; show_time(); cout << endl;
#endif
}

// test with a size nr x nc matrix of type scalar whose entries are in Mij
void testHNFetc(int nr, int nc, const vector<int>& Mij)
{
  vector<scalar> entries(Mij.size());
  std::transform(Mij.begin(), Mij.end(), entries.begin(), [](const int& aij){return scalar(aij);});
  mat M = mat(nr,nc,entries);
  cout << "M = \n" << M << endl;
  cout << "max entry of M: " << maxabs(M) << endl;
  mat H = HNF(M);
  cout << "HNF(M) =\n" << H << endl;
  mat S = SNF(M);
  cout << "SNF(M) =\n" << S << endl;
  mat L = LLL(M);
  cout << "LLL(M) = L =\n" << L << endl;
  cout << "max entry of L: " << maxabs(L) << endl;

#if (SCALAR_OPTION==3) // ZZ
  // Some tests of Qmat class (rational matrices, only implemented for scalar ZZ)
  cout << "Testing Qmat class: " << endl;

  Qmat QM(M, scalar(6));
  cout << "M/6 = \n" << QM << endl;
  cout << "with HNF\n" << HNF(QM) << endl;
  cout << "and  SNF\n" << SNF(QM) << endl;
#endif
}
