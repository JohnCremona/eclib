// TEGROS.CC, test program for egros functions
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2025 John Cremona
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
//

#include <eclib/curve.h>
#include <eclib/egros.h>

int main(void)
{
  set_precision(30);
  initprimes("PRIMES",0);
  //the_primes.init(25000000);
  Curve E;

  const vector<bigint> S5 = {bigint(2), bigint(3), bigint(5)};
  const vector<bigint> S23 = {bigint(2), bigint(3)};
  const vector<bigint> S35 = {bigint(3), bigint(5)};
  const vector<bigint> S11 = {bigint(11)};

  vector<bigint> Q_S5_2 = twist_factors(S5,2);
  cout << "Q("<<S5<<", 2) = "<<Q_S5_2<<" has "<<Q_S5_2.size() << " elements (should be "<<2*pow(2,S5.size())<<")\n";
  vector<bigint> Q_S5_4 = twist_factors(S5,4);
  cout << "Q("<<S5<<", 4) has "<<Q_S5_4.size() << " elements (should be "<<2*pow(4,S5.size())<<")\n";
  vector<bigint> Q_S5_6 = twist_factors(S5,6);
  cout << "Q("<<S5<<", 6) has "<<Q_S5_6.size() << " elements (should be "<<2*pow(6,S5.size())<<")\n";
  cout<<endl;

  bigrational j(bigint(272223782641), bigint(164025));
  cout<<"Curves with good reduction outside "<<S5<<" and j = "<<j;
  if (is_j_possible(j,S5))
    cout <<" do exist"<<endl;
  else
    cout<<" do NOT exist"<<endl;
  cout<<endl;

  vector<CurveRed> egr_S23_0 = egros_from_j_0(S23);
  std::sort(egr_S23_0.begin(), egr_S23_0.end());
  cout << egr_S23_0.size()<< " curves with j=0 and good reduction outside "<<S23<<" (should be 72):\n";
  for (auto E: egr_S23_0) cout<<(Curve)E<<" conductor "<<E.conductor()<<" sort key "<<E.sort_key()<<endl;
  cout<<endl;

  vector<CurveRed> egr_S23_1728 = egros_from_j_1728(S23);
  std::sort(egr_S23_1728.begin(), egr_S23_1728.end());
  cout << egr_S23_1728.size()<< " curves with j=1728 and good reduction outside "<<S23<<" (should be 32):\n";
  for (auto E: egr_S23_1728) cout<<(Curve)E<<" conductor "<<E.conductor()<<" sort key "<<E.sort_key()<<endl;
  cout<<endl;

  cout << "Elliptic curves with conductor a power of 11, from their known j-invariants" << endl;
  vector<bigrational> j11 = {
    bigrational(bigint(-122023936), bigint(161051)),
    bigrational(bigint(-52893159101157376), bigint(11)),
    bigrational(bigint(-4096), bigint(11)),
    bigrational(bigint(-121)),
    bigrational(bigint(-32768)),
    bigrational(bigint(-24729001))
  };
  vector<CurveRed> egr_11;
  for (auto j: j11)
    {
      auto EE = egros_from_j(j, S11);
      cout << EE.size() << " curves with j = " << j << ":";
      for ( auto E: EE) cout << " " << (Curve)E;
      cout << endl;
      egr_11.insert(egr_11.end(), EE.cbegin(), EE.cend());
    }
  std::sort(egr_11.begin(), egr_11.end());
  cout<<"Sorted list:\n";
  for (auto E: egr_11)
    cout << "conductor " << E.conductor() << "\t" << (Curve)E << "\tj = " << j_invariant(E) << endl;
  cout << endl;

  vector<bigint> N_j_0, N_j_1728;
  vector<CurveRed> E_j_0, E_j_1728;
  for (int n=1; n<=100; n++)
    {
      if (!is_valid_conductor(n)) continue;
      bigint N(n);
      vector<bigint> S = pdivs(N);
      if (is_N_possible_j_0(N, S))
        {
          N_j_0.push_back(N);
          auto EE = egros_from_j_0(S);
          for (auto E: EE)
            if (E.conductor()==N)
              E_j_0.push_back(E);
        }
      if (is_N_possible_j_1728(N, S))
        {
          N_j_1728.push_back(N);
          auto EE = egros_from_j_1728(S);
          for (auto E: EE)
            if (E.conductor()==N)
              E_j_1728.push_back(E);
        }
    }
  cout << "Possible conductors <= 100 of curves with j=0:    " << N_j_0 << endl;
  cout << "Actual conductors and curves:\n";
  for (auto E: E_j_0)
    cout << E.conductor() << "\t" << (Curve)E << endl;
  cout << "Possible conductors < 100 of curves with j=1728: " << N_j_1728 << endl;
  cout << "Actual conductors and curves:\n";
  for (auto E: E_j_1728)
    cout << E.conductor() << "\t" << (Curve)E << endl;

  return 0;
}
