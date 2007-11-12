// FILE GETOLDEIGS.CC
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
//
#include <iostream.h>
#include <streambuf.h>
#include "arith.h"

char command[40];
char eigfilename[20];

int main(void)
{
  long level;
  cout<<"Enter level: "; cin>>level;
  longlist oldlevels = posdivs(level);
  for (int i=1; i<oldlevels.length-1; i++)
    {
      long m = oldlevels[i];
//      cout<<"oldlevel "<<m<<"\n";
      if(m<11) continue;
      sprintf(eigfilename,"eigs/ex%d\0",m);
      sprintf(command,"cd eigs; extract %d\0",m);
//      cout << "Extracting " << eigfilename << " from tar file.\n";
      system(command);
    }
}
