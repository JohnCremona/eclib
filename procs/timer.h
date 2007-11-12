// timer.h:  declarations of timer functions
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
 
#ifndef _TIMER_H
#define _TIMER_H      1
                           //flags that this file has been included

#ifdef LiDIA_INTS
#include "LiDIA/timer.h"
#else
#ifdef NTL_INTS
#include "NTL/tools.h"
#else
#include <time.h>
#include <sys/times.h>
#endif
#endif

void init_time();
void start_time();
void stop_time();
void show_time();

#ifndef TIME_CONICS
#define TIME_CONICS 0
#endif

void init_conic_timer();
void start_conic_timer();
void stop_conic_timer();
void show_conic_timer();

#endif
