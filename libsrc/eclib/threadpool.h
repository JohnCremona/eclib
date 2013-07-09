/**
 * threadpool.h
 *
 * Declarations for threadpool class
 * for C++11 threads or Boost threads
 */

// Include only if Boost installed
#ifdef MULTITHREAD

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

class threadpool {
  public:
    threadpool( unsigned int numThreads = 0, int verbose = 0 );
    ~threadpool();

    /**
     * post()
     *
     * Add a job to the job queue so that the next
     * idle thread in the threadpool may execute it.
     * Templated function must reside in header file.
     */
    template< class Task >
    void post( Task &task ) {
      // Add reference to new task to job queue
      io_service_.post( boost::bind< void >( boost::ref( task ) ) );
    }

    void close();

    unsigned int getThreadCount();
    unsigned int getMaxThreads();

  private:
    unsigned int maxThreads;
    unsigned int threadCount;
             int verbose;

    boost::asio::io_service io_service_;
    boost::shared_ptr< boost::asio::io_service::work > work_;
    boost::thread_group     threads_;
    
};

#endif // THREADPOOL_H

#endif // MULTITHREAD
