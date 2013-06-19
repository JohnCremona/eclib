/**
 * threadpool.h
 *
 * Declarations for threadpool class
 * for C++11 threads or Boost threads
 */

#ifndef THREADPOOL_H
#define THREADPOOL_H

// Include headers based on
// available packages on system
#include <iostream>
#include <stdlib.h>
#include <boost/thread/thread.hpp>
#include <boost/thread/future.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

template< class Task >
class Threadpool {
  public:
    Threadpool( unsigned int numThreads );  
    ~Threadpool();

    void post( Task task );

    unsigned int getThreadCount();
    unsigned int getMaxThreads();

  private:
    unsigned int maxThreads;
    unsigned int threadCount;
             int verbose;

    boost::asio::io_service       io_service;
    boost::asio::io_service::work work;
    boost::thread_group           threads;
    
}

#endif
