// timer.h:  declarations of timer functions
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
//                     Marcus Mo    (timer class)
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
 
#ifndef _ECLIB_TIMER_H
#define _ECLIB_TIMER_H      1

// Timer header files
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <unordered_map>

#include "interface.h"

void init_time();
void start_time();
void stop_time();
void show_time(ostream& s = cout);

// New class object for modularity (same code)
class timer {
  public:
    timer();
    explicit timer( string filename );
    ~timer();                   

    void    stream( string filename = "" );
    void    start( string name = "default" );
    void    split( string name = "default" );
    void    stop( string name = "default" );
    void    stopAll();
    void    write( string message );
    void    show( int nline = 0, string name = "default", 
                  int idx1 = 0, int idx2 = -1 );
    void    showAll( int nline = 0 );
    void    clear( string name = "default" );
    void    clearAll();
    void    add( string name );
    void    list();

    int     count( string name = "default" );
    double  total( string name = "default" );
    double  average( string name = "default" );

    typedef unordered_map< string, vector<double> > timers;

  private:
    ostream* s_;
    ofstream file_;
    timers   times_;
    
    double   getWallTime();

    template< typename T > 
    string toString( T el );
};

#ifndef TIME_CONICS
#define TIME_CONICS 0
#endif

void init_conic_timer();
void start_conic_timer();
void stop_conic_timer();
void show_conic_timer(ostream& s = cout);

#endif
