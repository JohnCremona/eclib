// FILE OLDFORMS.H: declaration of class oldforms
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////

#define BINARY
#ifdef BINARY
#define EIG_FILE_PREFIX "eigs/x"
#else
#define EIG_FILE_PREFIX "eigs/e"
#endif

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

char* eigfile(long d);    //returns filename for eigs at level d
