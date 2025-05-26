// FILE xsplit_data.h : Declaration of class ff_data
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

#ifndef _ECLIB_XSPLIT_DATA_H
#define _ECLIB_XSPLIT_DATA_H

// Include headers
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>

// Disable multithreading
// #undef ECLIB_MULTITHREAD

#ifdef ECLIB_MULTITHREAD 
#include <boost/thread/mutex.hpp>
#endif

// Header files for macros and custom data structures
#include <eclib/method.h>

// Forward declaration of classes ... prevents circular includes.
// form_finder only required as pointer
class form_finder;

// Ennumerate node and child flags, exposed to global scope for outside access
enum nodestatus  { INTERNAL, ALL_OLD, FOUND_NEW, MAX_DEPTH };
enum childstatus { NOT_COMPLETE, COMPLETE, DESTROYED };

/**
 * ff_data (form_finder_data)
 *
 * Class name follows convention of that contained in xsplit.cc for ease of use.
 */
class ff_data {
  public:
    // Constructor, destructor.
    explicit ff_data( form_finder* ff );
    ~ff_data();

#ifdef ECLIB_MULTITHREAD
    void operator()();                              // Executed upon submission to job queue
#endif

    // Getters (to maintain consistency)
    nodestatus     status();                        // Return status of current node
    ssubspace*     abs_space();                     // Return parent absolute subspace
    ssubspace*     rel_space();                     // Return parent relative subspace
    long           depth();                         // Return current depth
    long           subdim();                        // Return subdimension
    long           eig();                           // Return associated eigenvalue
    vector< long > eiglist();                       // Return sequence of eigenvalues
    ff_data*       child( long eig );               // Return pointer to child
    int            numCompleteChildren();           // Return number of complete children
    bool           complete();                      // Return true if all children complete
  
    // Modifiers (available if needed i.e. not a friend class)
    void setStatus( nodestatus s );                 // Store status of current node

    void increaseDepth( long delta = 1 );           // Increase current depth
    void decreaseDepth( long delta = 1 );           // Decrease current depth
    
    void increaseSubmatUsage();                     // Locked counter incrementer

    void storeBplus( vec bp );                      // Store vector bplus
    void storeBminus( vec bm );                     // Store vector bminus
    
    void addChild( long eig, ff_data &child );      // Store new child of current node 
    void eraseChild( long eig );                    // Destroy child given eigenvalue
    void eraseChild( int idx );                     // Destroy child given index
    
    void setParent( ff_data *parent );              // Store parent node pointer
    void setEigenvalue( long eig );                 // Store eigenvalue

    void setChildren( vector<long> eigs );          // Stores number of children
                                                    //  and eigrange
    void childStatus( long eig, childstatus flag ); // Sets completed children
    
    void eraseChildren();                           // Destroys all children recursively

    // Make form_finder class a friend to gain access to protected/private methods
    friend class form_finder;

  private:
    form_finder*        ff_;                 // Back-pointer to form_finder class
                                             // Allows access to form_finder methods      

    nodestatus          status_;             // Status of current node
    long                depth_;              // Indicator of current depth
    long                subdim_;             // Dimension of current subspace
    long                eigenvalue_;         // Corresponding eigenvalue
    vector< long >      eigrange_;           // List of all of all eigenvalues at this depth (used by map)
    vector< long >      eiglist_;            // Sequence of eigenvalues leading to current
    vec                 bplus_, bminus_;
    ssubspace*          abs_space_;          // Current absolute subspace (dynamically created)
    ssubspace*          rel_space_;          // Current relative subspace (dynamically created)
    smat                conjmat_;            // Used only when plus==0 and bigmats==1
    smat                the_opmat_;
    smat                submat_;

    ff_data*            parent_;             // Pointer to parent data node
    vector<ff_data*>    children_;           // Pointers to corresponding data nodes
    ff_data*            child_;              // Pointer to favoured child
    int                 numChildren_;        // Number of children
    vector<childstatus> completedChildren_;  // Flags for child completion
    int                 submatUsage_;        // Counter for submat

#ifdef ECLIB_MULTITHREAD
    boost::mutex        childComplete_lock_; // Lock for completed children
    boost::mutex        submat_lock_;        // Lock for submat usage
    boost::mutex        go_up_lock_;         // Lock for go_up() method
#endif

    // Helper methods
    int map( long eig );                    // Map eigenvalue to index
};

#endif

// end of XSPLIT_DATA.H
