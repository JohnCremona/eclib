/**
 * timing.cc
 *
 * Testing timing class (found in libsrc/timer.cc).
 * Features are a simple, semantic interface;
 * subtimers; basic statistics.
 */

// Include headers
#include <iostream>
#include <sstream>
#include <eclib/interface.h>
#include <eclib/timer.h>

// Basic platform detection
#ifdef _WIN32
#include <windows.h>
#define powernap(x) Sleep(x)
#else
#include <unistd.h>
#define powernap(x) usleep((x)*1000)
#endif

int main( int argc, char **argv  ) {

  // Variables
  int         numTimers = 1;
  std::string prefix    = "timer";
  std::vector< std::string > names;
  int         verbose   = 0;

  // Verbosity
  std::cout << "Verbose (0/1): "; 
  std::cin  >> verbose;

  // Introduction
  std::cout << "Program timing - demonstration and testing." << std::endl;
  // Read in number of subtimers
  std::cout << "How many subtimers? ";
  std::cin  >> numTimers;
  
  // Initiate single timer object instance
  // Timings will be written to std::cout
  timer logbook;
  
  // OR: Write timings to a file. Pass a filename (string)
  // timer logbook( "times.dat" );

  // A default timer is automatically setup
  std::cout << "Start default timer" << std::endl;
  logbook.start();

  // Initiate the subtimers
  for( int t = 0; t < numTimers; t++ ) {
    // Create a name for each subtimer
    // Usually hard-coded into program but we just 
    // increment a counter and attach to a prefix.
    stringstream ss;
    ss << prefix << t;
    names.push_back( ss.str() );

    // add() takes in a string object
    std::cout << "Adding " << names[t] << std::endl;
    logbook.add( names[t] );
  }

  // Slumber for one second
  powernap(1000);

  // Start all subtimers
  for( int t = 0; t < numTimers; t++ ) {
    std::cout << "Starting " << names[t] << std::endl;
    logbook.start( names[t] );

    // Slumber for one second
    powernap(1000);
  }
 
  // Stop all timers, including default timer
  std::cout << "Stopping all timers" << std::endl;
  logbook.stopAll();

  // Write statistics to terminal
  // Commented out for `make check` to pass
  if( verbose ) logbook.showAll();

  exit( EXIT_SUCCESS );
}
