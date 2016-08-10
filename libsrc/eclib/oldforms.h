// FILE OLDFORMS.H: declaration of class oldforms
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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

class oldforms {
 public:
  long noldclasses, nap, ntp;
  long totalolddim;
 private:
  const level* N;
  int plusflag;
  vector< vector<long> > oldformap;
  vector<long> oldclassdims, oldlevels;
  void getoldclasses(long d, int verbose);
 public:
  oldforms(long intp, const level* iN, int verbose=0, int plus=1);   
             //intp = input value of ntp = max. number of Tp to use
  ~oldforms(){;}
  long dimoldpart(vector<long> aplist) const;
  void display(void) const;
};

// Utilities for creating filenames for curve files
bool file_exists(string filename);
string ecdb_filename(long N);
string single_curve_filename(long N);
string curve_filename(long N);


