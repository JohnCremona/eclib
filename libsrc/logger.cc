// FILE logger.cc : Implementation of member functions for class logger
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 Marcus Mo
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
 * logger.cc
 *
 * Simple thread-safe logger class. Prints to stdout.
 */

#include "eclib/logger.h"

int eclogger::level_ = 0;

/**
 * eclogger()
 *
 * Main constructor. Defaults verbosity to 0.
 */
eclogger::eclogger() {}

/**
 * ~eclogger()
 *
 * Destructor. Appends newline character.
 * Uses thread-safe fprintf, and flushes.
 */
eclogger::~eclogger() {
  //s "< std::endl"
  fprintf( stdout, "%s", s.str().c_str() );  
  fflush( stdout );
}

/**
 * level()
 *
 * Returns logger verbosity level.
 */
int eclogger::level() {
  return level_;
}

/**
 * setLevel()
 *
 * Allows for global change of verbosity level.
 */
void eclogger::setLevel( const int verbose ) {
  level_ = verbose;
}

/**
 * stream()
 *
 * Returns stream object.
 */
std::ostringstream& eclogger::stream() {
#if defined(ECLIB_MULTITHREAD) && defined(ECLIB_MULTITHREAD_DEBUG)
  s << "Thread " << boost::this_thread::get_id() << "\t"; 
#endif

  return s;
}

/**
 * stream()
 *
 * Returns stream object with FILE and LINE details.
 */
std::ostringstream& eclogger::stream( const char *file, const unsigned long line ) {
  std::string filename = file;
  s << std::setw(20) << filename << std::setw(5) << line << " ";
  return stream();  
}
