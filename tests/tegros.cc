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

  const vector<ZZ> S5 = {ZZ(2), ZZ(3), ZZ(5)};
  const vector<ZZ> S23 = {ZZ(2), ZZ(3)};
  const vector<ZZ> S35 = {ZZ(3), ZZ(5)};
  const vector<ZZ> S11 = {ZZ(11)};

  vector<ZZ> Q_S5_2 = twist_factors(S5,2);
  cout << "Q("<<S5<<", 2) = "<<Q_S5_2<<" has "<<Q_S5_2.size() << " elements (should be "<<2*pow(2,S5.size())<<")\n";
  vector<ZZ> Q_S5_4 = twist_factors(S5,4);
  cout << "Q("<<S5<<", 4) has "<<Q_S5_4.size() << " elements (should be "<<2*pow(4,S5.size())<<")\n";
  vector<ZZ> Q_S5_6 = twist_factors(S5,6);
  cout << "Q("<<S5<<", 6) has "<<Q_S5_6.size() << " elements (should be "<<2*pow(6,S5.size())<<")\n";
  cout<<endl;

  bigrational j(to_ZZ("272223782641"), ZZ(164025));
  cout<<"Curves with good reduction outside "<<S5<<" and j = "<<j;
  if (is_j_possible(j,S5))
    cout <<" do exist"<<endl;
  else
    cout<<" do NOT exist"<<endl;
  cout<<endl;

  vector<CurveRed> egr_S23_0 = egros_from_j_0(S23);
  cout << egr_S23_0.size()<< " curves with j=0 and good reduction outside "<<S23<<" (should be 72):\n";
  for (auto E1: egr_S23_0) cout<<(Curve)E1<<" conductor "<<E1.conductor()<<" sort key "<<E1.sort_key()<<endl;
  cout<<endl;

  vector<CurveRed> egr_S23_1728 = egros_from_j_1728(S23);
  cout << egr_S23_1728.size()<< " curves with j=1728 and good reduction outside "<<S23<<" (should be 32):\n";
  for (auto E1: egr_S23_1728) cout<<(Curve)E1<<" conductor "<<E1.conductor()<<" sort key "<<E1.sort_key()<<endl;
  cout<<endl;

  cout << "Elliptic curves with conductor a power of 11, from their known j-invariants" << endl;
  vector<bigrational> j11 = {
    bigrational(ZZ(-122023936), ZZ(161051)),
    bigrational(to_ZZ("-52893159101157376"), ZZ(11)),
    bigrational(ZZ(-4096), ZZ(11)),
    bigrational(ZZ(-121)),
    bigrational(ZZ(-32768)),
    bigrational(ZZ(-24729001))
  };
  vector<CurveRed> egr_11;
  for (auto ji: j11)
    {
      auto EE = egros_from_j(ji, S11);
      cout << EE.size() << " curves with j = " << ji << ":";
      for ( auto E1: EE) cout << " " << (Curve)E1;
      cout << endl;
      egr_11.insert(egr_11.end(), EE.cbegin(), EE.cend());
    }
  cout<<endl;
  std::sort(egr_11.begin(), egr_11.end());
  cout<<"Full sorted list:\n";
  for (auto E1: egr_11)
    cout << "conductor " << E1.conductor() << "\t" << (Curve)E1 << "\tj = " << j_invariant(E1) << endl;
  cout << endl;

  vector<ZZ> N_j_0, N_j_1728;
  vector<CurveRed> E_j_0, E_j_1728;
  for (int n=1; n<=100; n++)
    {
      if (!is_valid_conductor(n)) continue;
      ZZ N(n);
      vector<ZZ> S = pdivs(N);
      if (is_N_possible_j_0(N, S))
        {
          N_j_0.push_back(N);
          auto EE = egros_from_j_0(S);
          for (auto E1: EE)
            if (E1.conductor()==N)
              E_j_0.push_back(E1);
        }
      if (is_N_possible_j_1728(N, S))
        {
          N_j_1728.push_back(N);
          auto EE = egros_from_j_1728(S);
          for (auto E1: EE)
            if (E1.conductor()==N)
              E_j_1728.push_back(E1);
        }
    }
  cout << "Possible conductors <= 100 of curves with j=0:    " << N_j_0 << endl;
  cout << "Actual conductors and curves:\n";
  for (auto E1: E_j_0)
    cout << E1.conductor() << "\t" << (Curve)E1 << endl;
  cout << endl;
  cout << "Possible conductors <= 100 of curves with j=1728: " << N_j_1728 << endl;
  cout << "Actual conductors and curves:\n";
  for (auto E1: E_j_1728)
    cout << E1.conductor() << "\t" << (Curve)E1 << endl;

  return 0;
}
