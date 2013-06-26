// xsplit.h: Declaration of class form_finder
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
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
 
#if     !defined(_XSPLIT_H)
#define _XSPLIT_H      1       //flags that this file has been included

//#include <boost/threads/mutex.hpp>

#include "method.h"  // #defines form_finder=form_finder0/1/2/3/4
#include "splitbase.h"
#include "xsplit_data.h"

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
    
    vec  getbasis() const {return bplus;}
    vec  getbasisplus() const {return bplus;}
    vec  getbasisminus() const {return bminus;}
  
  protected:
    splitter_base* h;
    
    int            plusflag, dual, bigmats, verbose, targetdim;
    int            gnfcount;              // Global newform counter
    long           maxdepth, mindepth, depth, subdim, dimen;
    SCALAR         denom1;
    vector<long>   eiglist;
    vec            bplus, bminus;
    vector< vector<long> > gaplist;       // Vector to hold all (sub)eiglists
    vector<vec>    gbplus, gbminus;       // Vector to hold all bplus/bminus
    ssubspace**    nest;                  // Array of pointers to subspaces.
                                          // "Current" subspace is *nest[depth] 
                                          // of dimension subdim

    int*           havemat;
    vector<string> opfilenames;           // Temp filenames
    smat           conjmat;               // Only used if plus==0 and bigmats==1
    smat           the_opmat;
    smat*          submats;               // Holds current restriction for i>0
   
    // NEW THINGS REQUIRED
    //threadpool/job queue
    ff_data* root; 
    ff_data* current;
    //boost::mutex store_lock;              // Lock for store() function 
    
    void make_opmat(long i, ff_data *data); // Puts it in the_opmat
    void make_submat(ff_data *data);
    void go_down(ff_data &data, long eig, int last=0);
    void go_up();
    void make_basis(ff_data &data);
    vec  getbasis1(const ssubspace* s);   // Assuming dim(s)=1, get basis vector
};

#endif
