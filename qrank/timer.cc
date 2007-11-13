// timer.cc: implementation of timer functions
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
#include "timer.h"

#ifdef LiDIA_INTS
timer T;
void init_time() { T.set_print_mode(HMS_MODE);}
void start_time() { T.start_timer();}
void stop_time()  { T.stop_timer();}
void show_time() { cerr.form(" (%.2f seconds)\n", T.user_time()/100.0); }
//void show_time()   { cerr << " (" << T << ")\n";}

timer conic_timer;
double conic_time;
void init_conic_timer() { conic_timer.set_print_mode(HMS_MODE); conic_time=0;}
void start_conic_timer() { conic_timer.start_timer();}
void stop_conic_timer()  { conic_timer.stop_timer(); conic_time+=conic_timer.user_time()/100.0;}
void show_conic_timer() { cerr.form("Time spent so far in solving conics was %.2f seconds\n", conic_time); }

#else
struct tms buffer;
#ifndef HZ
#define HZ 60  // fix for alpha
#endif
long starttime,stoptime;
void init_time() { ;}
void start_time() { times(&buffer); starttime = buffer.tms_utime;}
void stop_time()  { times(&buffer); stoptime = buffer.tms_utime;}
void show_time()   { cerr.form(" (%.2f seconds)\n",(stoptime-starttime)*1./HZ);}

void init_conic_timer() {;}
void start_conic_timer() {;}
void stop_conic_timer()  {;}
void show_conic_timer() {;}


#endif
