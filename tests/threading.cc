// FILE threading.cc : Simple test program for threadpool class
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

// Enable/disable multithreading
// #undef MULTITHREAD

// Include headers
#include <iostream>
#include <iomanip>
#include <cmath>
#include <eclib/interface.h>
#include <eclib/timer.h>
#include <eclib/threadpool.h>

// Primality testing task class
class isprime {
  public:
    explicit isprime( unsigned int n )
      : n_( n ), prime_( 1 ) {}
    isprime( unsigned int min, unsigned int max )
      : min_( min ), max_( max ) {}
    ~isprime() {;}

    void operator()() {
      // Find square root
      unsigned int sr = sqrt( (double) n_ );

      // Loop
      for( unsigned int i = 2; i <= sr; i++ ) {
        if( ( n_ % i ) == 0 ) {
          prime_ = 0;
          break;
        }
      }
    }

    unsigned int n() {
      return n_;
    }

    bool isPrime() {
      return prime_;
    }

  private:
    unsigned int n_;
    unsigned int min_;
    unsigned int max_;
    bool prime_;
};

// Main function.
// Creates tasks and posts to job queue
int main( int argc, char **argv ) {

  // Variables (with default values)
  unsigned int N = 10000;                 // Number of tasks
  unsigned int nt = 1;                    // Number of threads
  unsigned int v = 0;                     // Verbosity
  std::vector< isprime* > tasks;          // Array to hold tasks
  int count = 0;                          // Counter

  // Read in command line arguments
  if( argc > 1 )  N = atoi( argv[1] );     // Read in number of tasks
  if( argc > 2 ) nt = atoi( argv[2] );     // Read in number of threads

  // Initiate timer
  timer profile;

  // Start default timer
  profile.start();

#ifdef MULTITHREAD
  // Initiate threadpool/job queue with verbose output
  threadpool pool( nt, v );
#endif

  // Create tasks and add to array for later reference
  for( unsigned int i = 2; i < N; i++ ) {
    // Create a new task object
    isprime* task = new isprime( i );

    // Add task to container so it can be accessed later
    tasks.push_back( task );

#ifdef MULTITHREAD
    // Add task to job queue
    pool.post< isprime >( *tasks.back() );
#else
    // Serial version
    task -> operator()();
#endif
  }

#ifdef MULTITHREAD
  // Wait for all jobs to be complete
  pool.close();
#endif

  // Stop timer
  profile.stop();

  // Print out results
  std::cout << "Primes up to " << N << std::endl;
  for (const auto& t : tasks) {
    if( t -> isPrime() ) {
      std::cout << std::setw(10) << t -> n() << " ";
      count++;
      if( (count % 10) == 0 ) std::cout << std::endl;
    }
  }
  std::cout << std::endl;

  // Print out run time
  std::cout << "Running with " << nt << " threads" << std::endl;
  std::cout << "There are " << count << " primes up to " << N << std::endl;
  if( v ) profile.show( 1 );

  // Delete tasks
  for( auto it : tasks)
     delete it;

  exit( EXIT_SUCCESS );

}
