// FILE MANIN.H: declaration of class manin
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

class manin :public newforms {
private:
  int easy; long nq, dq;
  vector<long> dp0, pdotlist, qdotlist; // , ipdotlist, iqdotlist;
  vec initvec;
  void findq();    //Computes nq, dq, qdotlist
  void getoneap(long p, int output, ofstream& out, int verbose=0);
public:
//Constructor:
  manin(long n, int useolddata=0, long ntp=5, int cuspidalflag=0, int disp=0);
//Destructor:
  ~manin() {;}
//Main function to produce Hecke eigenvals
  void getap(long first, long last, int output, char* eigfile, int verbose=0);
};
