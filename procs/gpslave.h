// parislave.h: class for starting up a "slave" background gp process 
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
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
 
#include "interface.h"

// The slave  can be sent input via a fifo file and the output
// recovered.  Intended for integer factorization only at present.

class parislave {
  char gpinfilename[18];
  char gpoutput[1000];
  FILE* gp;  // pointer to file from which gp takes its input
  ofstream gpin; // ofstream attached to latter
  int dummy;  // =1 iff we have no gp to use
 public:
  parislave();
  ~parislave();
  vector<bigint> factor(const bigint& n, int proof=1);
  int is_prime(const bigint& n);
  int is_running() {return !dummy;}
};

extern parislave the_pari_slave;  // The one and only instance
