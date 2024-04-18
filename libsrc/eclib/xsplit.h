// xsplit.h: Declaration of class form_finder
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
//                     Marcus Mo     (parallel code)
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
 
#if     !defined(_ECLIB_XSPLIT_H)
#define _ECLIB_XSPLIT_H      1       //flags that this file has been included

// Disable multithreading
// #undef ECLIB_MULTITHREAD

#ifdef ECLIB_MULTITHREAD 
#include <boost/thread/mutex.hpp>
#endif

#include <eclib/xsplit_data.h>
#include <eclib/splitbase.h>
#ifdef ECLIB_MULTITHREAD
#include <eclib/threadpool.h>
#endif

// flags set on construction:

// dual:  use dual bases throughout (the usual case)
// plusflag: if 1 then targetdim=1 and conjmat not used
//           if 0 then targetdim=2 and conjmat is used

// bigmats: if 1 then uses opmat() to get operators on ambient space,
// does its own restriction to the current subspace; if 0 then uses
// opmat_restricted() to get operators restricted to current space
// directly -- ONLY possible when dual=1

class form_finder {
  public:
    form_finder(splitter_base* hh, int plus, int maxd, int mind=0, 
                int dualflag=1, int bigmatsflag=0, int v=0);
    ~form_finder(void); 
    
    void find();
    void find(ff_data &data);
    void recover(vector< vector<long> > eigs);
    void splitoff(const vector<long>& eigs);
    void store(vec bp, vec bm, vector<long> eigs);
    
    vec  getbasis( const ff_data &data ) const {return data.bplus_;}
    vec  getbasisplus( const ff_data &data ) const {return data.bplus_;}
    vec  getbasisminus( const ff_data &data ) const {return data.bminus_;}

    friend class ff_data; 
  
  protected:
    splitter_base* h;
    
    int            plusflag, dual, bigmats, verbose, targetdim;
    int            gnfcount;                  // Global newform counter
    long           maxdepth, mindepth, dimen;
    SCALAR         denom1;
    vector< vector<long> > gaplist;           // Vector to hold all (sub)eiglists
    vector<vec>    gbplus, gbminus;           // Vector to hold all bplus/bminus

    ff_data* root;                            // Always points to root data node

    void make_opmat(long i, ff_data &data);   // Puts it in the_opmat
    void make_submat(ff_data &data);
    smat make_nested_submat(long ip, ff_data &data);
    void go_down(ff_data &data, long eig, int last=0);
    void go_up( ff_data &data );
    void make_basis(ff_data &data);
    vec make_basis1(ff_data &data);  // current space has dim. 1
    vec make_basis2(ff_data &data, const vec& v);  // current space has dim. >1 and
                                                   // relative basis vector v

#ifdef ECLIB_MULTITHREAD
    threadpool   pool;                        // Job queue
    boost::mutex store_lock;                  // Lock for store() function
#endif
};

vec  getbasis1(const ssubspace* s);       // Assuming dim(s)=1, get basis vector
vec  lift(const vec& v);                  // Lift basis vector

#endif
