// FILE tversion.cc : Simple test of eclib versioning utilities
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 Marcus Mo
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

#include <iostream>
#include <eclib/version.h>

int main( int argc, char **argv  )
{
  cout << "This eclib version is " << eclib_version() << endl;
  vector<int> date = eclib_date();
  cout << "[ year month day ] = " << date << endl;

  cout << "\nConfigure options and compilation date:\n" << endl;
  show_version();

  cout << "\nTesting comparison of version date and various (year, month, day) date triples:\n" << endl;
  // test version date comparison
  int y, m, d, s;
  y=2020; m=1; d=1;
  s = compare_eclib_version(y,m,d);
  cout << "Date ";
  cout << setw(4)<<std::setfill('0')<<y;
  cout << setw(2)<<std::setfill('0')<<m;
  cout << setw(2)<<std::setfill('0')<<d;
  cout << " is ";
  if (s>0) cout << "before"; else
    if (s<0) cout << "after"; else
      cout << "the same as";
  cout << " the version date ";
  cout << setw(4)<<std::setfill('0')<<date[0];
  cout << setw(2)<<std::setfill('0')<<date[1];
  cout << setw(2)<<std::setfill('0')<<date[2];
  cout << endl;

  y=2030; m=6; d=11;
  s = compare_eclib_version(y,m,d);
  cout << "Date ";
  cout << setw(4)<<std::setfill('0')<<y;
  cout << setw(2)<<std::setfill('0')<<m;
  cout << setw(2)<<std::setfill('0')<<d;
  cout << " is ";
  if (s>0) cout << "before"; else
    if (s<0) cout << "after"; else
      cout << "the same as";
  cout << " the version date ";
  cout << setw(4)<<std::setfill('0')<<date[0];
  cout << setw(2)<<std::setfill('0')<<date[1];
  cout << setw(2)<<std::setfill('0')<<date[2];
  cout << endl;

  y=2021; m=3; d=17;
  s = compare_eclib_version(y,m,d);
  cout << "Date ";
  cout << setw(4)<<std::setfill('0')<<y;
  cout << setw(2)<<std::setfill('0')<<m;
  cout << setw(2)<<std::setfill('0')<<d;
  cout << " is ";
  if (s>0) cout << "before"; else
    if (s<0) cout << "after"; else
      cout << "the same as";
  cout << " the version date ";
  cout << setw(4)<<std::setfill('0')<<date[0];
  cout << setw(2)<<std::setfill('0')<<date[1];
  cout << setw(2)<<std::setfill('0')<<date[2];
  cout << endl;

  exit(0);
}
