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

#ifdef LiDIA_INTS
timer T;
void init_time() { T.set_print_mode(HMS_MODE);}
void start_time() { T.start_timer();}
void stop_time()  { T.stop_timer();}
void show_time(ostream& s) 
{ 
  s<<" ("<<T.user_time()/1000.0<<" seconds)"<<flush; 
}

timer conic_timer;
double conic_time;
void init_conic_timer() { conic_timer.set_print_mode(HMS_MODE); conic_time=0;}
void start_conic_timer() { conic_timer.start_timer();}
void stop_conic_timer()  { conic_timer.stop_timer(); conic_time+=conic_timer.user_time()/100.0;}
void show_conic_timer(ostream& s) 
{
  s<<"Time spent so far in solving conics was "<<conic_time<<" seconds\n"; 
}

#else
#ifdef NTL_INTS
double starttime,stoptime;
void init_time() { ;}
void start_time() { starttime = GetTime();}
void stop_time()  { stoptime  = GetTime();}
void show_time(ostream& s) 
{
  s<<" ("<<(stoptime-starttime)<<" seconds)" << flush; 
}

#else // not LiDIA

struct tms buffer;
long starttime,stoptime;
void init_time() { ;}
void start_time() { times(&buffer); starttime = buffer.tms_utime;}
void stop_time()  { times(&buffer); stoptime = buffer.tms_utime;}
void show_time(ostream& s) { s<<setprecision(3)<<" ("<<(stoptime-starttime)*1./CLK_TCK<<" seconds)"<<flush; }

void init_conic_timer() {;}
void start_conic_timer() {;}
void stop_conic_timer()  {;}
void show_conic_timer() {;}


#endif

#endif
