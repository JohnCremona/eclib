// timer.cc:  implementations of timer functions
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
 
#include <iostream>
#include <iomanip>
#include "interface.h"
#include "timer.h"

#ifndef CLK_TCK
#define CLK_TCK	CLOCKS_PER_SEC
#endif

double starttime,stoptime;
void init_time() { ;}
void start_time() { starttime = GetTime();}
void stop_time()  { stoptime  = GetTime();}
void show_time(ostream& s) 
{
  s<<" ("<<(stoptime-starttime)<<" seconds)" << flush; 
}
