// FILE logger.h : Declaration of class logger
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 Marcus Mo
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

/**
 * logger.h
 *
 * Declarations for logger class. 
 * Verbosity levels are integers:
 *    -1 : Error. Always print
 *     0 : Print if verbose > 0
 *     1 : Print if verbose > 1
 *    ... etc ...
 */

#ifndef LOGGER_H
#define LOGGER_H

// Include headers
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>

// Disable multithreading
// #undef ECLIB_MULTITHREAD

#ifdef ECLIB_MULTITHREAD
#include <boost/thread/thread.hpp>

// Uncomment if multithreading debug messages required
//#define ECLIB_MULTITHREAD_DEBUG
#endif

// Logging configurations
#ifndef ECLOG1
#define ECLOG1(v) if(eclogger::level() <= (v)) ; else eclogger().stream()
#endif
#ifndef ECLOG2
#define ECLOG2(v) if(eclogger::level() <= (v)) ; else eclogger().stream(__FILE__,__LINE__)
#endif

// Expose logging macro
#ifndef ECLOG
#define ECLOG ECLOG1
#endif

class eclogger {
  public:
    eclogger();
    ~eclogger();

    std::ostringstream& stream();
    std::ostringstream& stream( const char *file, const unsigned long line );

    static int  level();
    static void setLevel( const int verbose );

  private:
    static int level_;
    std::ostringstream s;
};

#endif
