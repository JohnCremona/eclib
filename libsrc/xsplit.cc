// xsplit.cc: implementation of class form_finder
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
 
#define ECLIB_INT_NUM_THREADS 8
#define ECLIB_RECURSION_DIM_LIMIT 5821
//#define ECLIB_MULTITHREAD_DEBUG

#include <eclib/logger.h>
#include <eclib/method.h>
#include <eclib/xsplit.h>

template class form_finderT<int>;
template class form_finderT<long>;
template class form_finderT<bigint>;

template<class T>
Zvec<T> lift(const Zvec<T>& v, T mod)
{
  Zvec<T> w;
  if ( lift(v, mod, w) )
    return w;
  else
    {
      cout << "Unable to lift eigenvector " << v << " from Z/" << mod << " to Z" << endl;
      return v;
    }
  return w;
}

// CLASS FORM_FINDER (was called splitter)

template<class T>
form_finderT<T>::form_finderT(splitter_base<T>* hh, T mod,
                              int plus, int maxd, int mind,
                              int dualflag, int bigmatsflag, int v)
  :h(hh), modulus(mod), plusflag(plus), dual(dualflag), bigmats(bigmatsflag), verbose(v),
   gnfcount(0), maxdepth(maxd), mindepth(mind)
{
  //cout<<"In form_finder constructor, modulus="<<modulus<<", plusflag="<<plus<<", maxd="<<maxd<<", mind="<<mind<<", dualflag="<<dualflag<<", bigmatsflag="<<bigmatsflag<<endl;
  eclogger::setLevel( verbose );
  denom1 = h->matden();
  dimen  = h->matdim();

  // Create and initialise new data object as root node
  // passing a constant pointer of current form_finder object
  // to data class constructor
  root = new ff_data<T>( this );

  // Set initial values
  // form_finder class is a friend of ff_data class
  // so may access private members
  root -> subdim_ = dimen;
 
  targetdim = 1;
  if( !plusflag ) {           // full conjmat not needed when plusflag is true
    targetdim=2;
    if( bigmats ) root -> conjmat_ = h -> s_opmat(-1,dual); 
  }
}

template<class T>
form_finderT<T>::~form_finderT(void) {
  // Decendants of root node will be recursively deleted
  // if they have not already been deleted during find()
  // All dynamically created objects (subspaces) held 
  // in each data node will also be deleted.
  delete root;
}

// This is only used when bigmats==1 and we compute opmats on the entire ambient space:
template<class T>
void form_finderT<T>::make_opmat(long i, ff_data<T> &data) { 
  data.the_opmat_ = h -> s_opmat(i,dual,verbose); 
}

template<class T>
void form_finderT<T>::make_submat( ff_data<T> &data ) {
  // Cache current data node depth
  long depth = data.depth_;

  if( bigmats ) { 
    // fetch the_opmat from file, or compute
    make_opmat(depth,data);

    if( depth == 0 ) data.submat_ = data.the_opmat_;
    else {
	    ECLOG(1) << "restricting the_opmat to subspace...";
	    data.submat_ = restrict_mat(data.the_opmat_,*(data.abs_space_));
	    ECLOG(1) << "done." << endl;
	  }

    data.the_opmat_ = sZmat<T>(0,0); // releases its space
  }
  else {
    if( data.submat_.nrows() == 0 ) // else we have it already
      {
        if( depth == 0 )
          data.submat_ = h -> s_opmat(depth,1,verbose);
        else
          {
            //data.submat_ = h -> s_opmat_restricted(depth,*(data.abs_space_),1,verbose);
            data.submat_ = make_nested_submat(depth,data);
          }
      }
  }
}

/**
 * make_nested_submat()
 *
 * Computes and returns the submat -- new nested version
 */

template<class T>
sZmat<T> form_finderT<T>::make_nested_submat(long ip, ff_data<T> &data)
{
  long depth = data.depth_;  // current depth
  long subdim = data.subdim_;  // current dimension
  ff_data<T> *d = &data; // Pointer to nodes

  ECLOG(1) << "Computing operator of size " << subdim
           << " at depth " << depth << "..." << flush;

  // first we go up the chain, composing pivotal indices

  Zvec<int> jlist = Zvec<int>::iota(subdim);
  sZmat<T> b = d->rel_space_->bas();
  int level = depth;
  while (level--)
    {
      ECLOG(2) << "["<<level<<"]" << flush;
      jlist = d->rel_space_->pivs()[jlist];
      d->parent_->child_ = d;
      d = d->parent_;
      if(level) b = mult_mod_p(d->rel_space_->bas(), b, modulus);
    }

  // now compute the matrix of images of the j'th generator for j in jlist
  ECLOG(2) << " basis done..." << flush;
  sZmat<T> m = h -> s_opmat_cols(ip, jlist, 0);
  ECLOG(2) << " sub-opmat done..." << flush;
  m = mult_mod_p(m,b,modulus);
  ECLOG(1) <<" opmat done."<<endl;
  return m;
}


/**
 * go_down()
 *
 * Initiates creation of new subspace; data stored in new
 * data node. Data node passed as parameter will become the
 * _parent_ of this new data node.
 */
template<class T>
void form_finderT<T>::go_down(ff_data<T> &data, long eig, int last) {
  // Cache current depth
  long depth = data.depth_;

  // Locate required child w.r.t test eigenvalue
  ff_data<T> *child = data.child( eig ); 

  // Set new depth
  child -> depth_ = depth + 1;
   
  T eig2 = eig*denom1;

  ECLOG(1) << "Increasing depth to " << depth+1 << ", "
           << "trying eig = " << eig << "..."
           << "after scaling, eig =  " << eig2 << "..." << endl;
  ssubZspace<T> s(0, modulus);

  vector<int> submat_dim = dim(data.submat_);
  stringstream submat_dim_ss;
  std::copy(submat_dim.begin(),submat_dim.end(),ostream_iterator<int>(submat_dim_ss," "));

  ECLOG(1) << "Using sparse elimination (size = [ "
                     << submat_dim_ss.str() << "], density ="
		                 << density(data.submat_) << ")..." << flush;
  ECLOG(3) << "submat = " << data.submat_ << flush;

  s = eigenspace(data.submat_,eig2, modulus); // the relative eigenspace

  // Increment data usage counter for parent node
  data.increaseSubmatUsage();

  // Reset current submat if all children have used it
  // Save space (will recompute when needed)
  //if(    ( depth == 0 )
  //    && ( dim(s) > 0 )
  //    && ( data.submat_.nrows() > 1000 )
  //    && ( data.submatUsage_ == data.numChildren_ ) ) {
  //  data.submat_ = smat(0,0); 
  //}

  ECLOG(1) << "done (dim = " << dim(s) << ")"<<endl;
  // ECLOG(1) << ", combining subspaces..." << flush;

  child -> rel_space_ = new ssubZspace<T>(s);
  // if( depth == 0 )
  //   child -> abs_space_ = new ssubZspace<T>(s);
  // else
  //   child -> abs_space_ = new ssubZspace<T>(combine( *(data.abs_space_),s ));
  // ECLOG(1) << "done." << endl;
  
  depth++; // Local depth increment (does not effect data nodes)

  child -> subdim_ = dim( *(child -> rel_space_) );

  ECLOG(1) << "Eigenvalue " << eig 
           << " has multiplicity " << child -> subdim_ << endl;
  if(child -> subdim_>0) {
    ECLOG(0) << " eig " << eig 
             << " gives new subspace at depth " << depth
             << " of dimension " << child -> subdim_ << endl;
  }
}

template<class T>
void form_finderT<T>::go_up( ff_data<T> &data ) {
  // Cache pointer to parent data node for access after current node is deleted
  ff_data<T> *parent = data.parent_;

#ifdef ECLIB_MULTITHREAD
  // Lock parent node with scoped lock 
  boost::mutex::scoped_lock lock( parent -> go_up_lock_ );
#ifdef ECLIB_MULTITHREAD_DEBUG
  ECLOG(1) << "in go_up for eig=" << data.eigenvalue_ 
           << " depth=" << data.depth_
           <<" status=" << data.status_ << std::endl;
#endif
#endif

  // Erasing node via children array of parent which calls destructor
  // of object (ff_data<T>), using eigenvalue as key
  parent -> childStatus( data.eigenvalue_, COMPLETE );
  parent -> eraseChild( data.eigenvalue_ );

#ifdef ECLIB_MULTITHREAD
  lock.unlock();
  
  // Only last child to complete will execute the following (Detects if parent is root node)
  if( parent -> complete() && parent -> parent_ != NULL ) go_up( *parent );
#endif
}

template<class T>
void form_finderT<T>::make_basis( ff_data<T> &data ) {
  // Cache data values
  long depth  = data.depth();
  long subdim = data.subdim();
  vector< long > eiglist = data.eiglist();

  if( subdim != targetdim ) {
    cout << "error in form_finder::make_basis with eiglist = ";
    for(int i=0; i<depth; i++) 
      cout << eiglist[i] << ",";
    cout << "\nfinal subspace has dimension " << subdim << endl;
    cout << "aborting this branch!" << endl;
    return;
  }

  if(plusflag) {
    // must treat separately since we did not
    // define abs_space[0] in order to save space
    if(depth==0) {
      data.bplus_    = Zvec<T>(dimen);
      data.bplus_[1] = 1;
	  }
    else {
      data.bplus_ = make_basis1(data);
    }

    return;
  }

  ssubZspace<T> *spm_rel; //, *spm_abs;
  T eig = denom1;
  sZmat<T> subconjmat;            // only used when depth>0
  if( bigmats ) {
    ssubZspace<T>* s;
    s = data.abs_space_;  // only used when depth>0
    cout<<"data.abs_space_ = (pointer) "<<data.abs_space_<<endl;
    cout<<"s->modulus = "<<s->mod()<<endl;
    subconjmat = (depth) ? restrict_mat(data.conjmat_, *s) : data.conjmat_;
    // will only be a 2x2 in this case (genus 1 only!)
  }
  else {
    subconjmat = make_nested_submat(-1,data);
  }

  // C++11 loop over two variables (similar to python)
  // for( int b : { -1,+1 } ) { /* use b as -1 or +1 */ }
  for(long signeig=+1; signeig>-2; signeig-=2) {
    T seig;
    seig = eig;

    if(signeig<0) seig = -eig;

    if(depth) {
      spm_rel = new ssubZspace<T>(eigenspace(subconjmat,seig, modulus));
	    //spm_abs  = new ssubZspace<T>(combine(*s,*spm_rel));
    }
    else {
      spm_rel = new ssubZspace<T>(eigenspace(subconjmat,seig, modulus));
      //spm_abs = spm_rel;
    }

    if(dim(*spm_rel)!=1) {
      cout << "error in form_finder::makebasis; ";
      cout << "\nfinal (";

      if(signeig>0) cout << "+";
      else cout << "-";

      cout << ") subspace has dimension " << dim(*spm_rel) << endl;
      cout << "aborting this branch!" << endl;
      //delete spm_abs;
      delete spm_rel;
      return;
    }

    Zvec<T> w = make_basis2(data, spm_rel->bas().as_mat().col(1));

    if(signeig>0) data.bplus_  = w;
    else          data.bminus_ = w;

    //delete spm_abs;
    delete spm_rel;
  }
}

template<class T>
Zvec<T> form_finderT<T>::make_basis2(ff_data<T> &data, const Zvec<T>& v)
{
  ff_data<T> *d = &data;
  int level = data.depth_;
  Zvec<T> w = v;
  while (level--)
    {
      w = mult_mod_p(d->rel_space_->bas(), w, modulus);
      d = d->parent_;
    }
  return lift(w, modulus);
}

template<class T>
Zvec<T> form_finderT<T>::make_basis1(ff_data<T> &data)
{
  Zvec<T> v(1);  v.set(1,T(1));
  return make_basis2(data, v);
}

template<class T>
void form_finderT<T>::recover(vector< vector<long> > eigs) {
  for(unsigned int iform=0; iform<eigs.size(); iform++) {
    if(verbose) {
	    cout << "Form number " << iform+1 << " with eigs ";

      int n = eigs[iform].size(); 
      if(n>10) n = 10;

	    copy(eigs[iform].begin(), eigs[iform].begin() + n, 
	       ostream_iterator<long>(cout, " ")); 
	    cout << "..." << endl;
	  }
 
    splitoff(eigs[iform]);
  }  
  // Clears all nodes.  This cannot be done automatically since we
  // don't know how many eigs lists were to be splt off.
  root -> eraseChildren();
}

template<class T>
void form_finderT<T>::splitoff(const vector<long>& eigs) {

  // Always start at root node
  ff_data<T> *current = root;

  // Temporary variables
  long depth  = current -> depth_;
  long subdim = current -> subdim_; 

  if( verbose ) {
    cout << "Entering form_finder, depth = " << depth 
         << ", dimension " << subdim << endl;
  }
 
  // Walk down nodes (if any already created) for common branches
  while( current -> numChildren_ > 0
         && current -> child(eigs[depth]) != NULL ) {
 
    // Update current node pointer
    current = current -> child(eigs[depth]);

    // Update data
    depth  = current -> depth_;
    subdim = current -> subdim_;
    if (verbose) {
      cout << "...increasing depth to " << depth 
           << ", dimension " << subdim << endl;
    }
  }
  
  // Current node is new branch point
  // We trim all sub-branches ...
  current -> eraseChildren();
  
  if( verbose ) {
    cout << "restarting at depth = " << depth << ", "
         << "dimension " << subdim << endl;
  }
  
  // ... and grow a new branch down to required depth.
  while( (subdim > targetdim) && (depth < maxdepth) ) {
    // Get number of possible eigenvalues
    if( current -> numChildren_ <= 0 ) {
      vector<long> t_eigs = h->eigrange(depth);
      current -> setChildren( t_eigs );
    }

    // Create new child node
    ff_data<T> *child = new ff_data<T>( this );
    
    // Configure data node ancestry
    current -> addChild( eigs[depth], *child );
 
    // Create submat for current node
    make_submat( *current );
    
    // Proceed to go down
    go_down(*current,eigs[depth],1);

    // Update to new values
    current = child;
    depth   = current -> depth_;
    subdim  = current -> subdim_;
  }

  // Creating newforms
  make_basis(*current);
  h->use(current->bplus_,current->bminus_,eigs); 

  return;
}

template<class T>
void form_finderT<T>::find() {
#ifdef ECLIB_MULTITHREAD
  // Set number of threads to use either through default
  // ECLIB_INT_NUM_THREADS macro defined above, or
  // ECLIB_EXT_NUM_THREADS environment variable.
  unsigned int eclib_num_threads = ECLIB_INT_NUM_THREADS;

  stringstream s;
  s << getenv("ECLIB_EXT_NUM_THREADS");
  if( !s.str().empty() ) eclib_num_threads = atoi(s.str().c_str());

  // Start job queue. We keep job queue local to ensure threads are 
  // not kept busy for longer than necessary.
  pool.start( eclib_num_threads, verbose );
#endif

  // Proceed in recursive find, passing a node through
  find( *root );
  
#ifdef ECLIB_MULTITHREAD
  // Join all threads in threadpool to wait for all jobs to finish
  // Or detect when all branches of the tree has been traversed
  pool.close();
#endif

  // Clear all nodes.  This should have been be done automatically but not all nodes are deleted when running in multithreaded mode.
  root -> eraseChildren();

  // Now compute all newforms only if recursion has finished
  if(verbose>1) cout << "Now performing use() on all lists at once" << endl;
  for( int nf = 0; nf < gnfcount; nf++ ) {
    h-> use(gbplus[nf],gbminus[nf],gaplist[nf]);
  }
}

template<class T>
void form_finderT<T>::find( ff_data<T> &data ) {
  // Cache values of current data
  long depth  = data.depth();
  long subdim = data.subdim();
  vector< long > eiglist = data.eiglist();

  vector<long> subeiglist(eiglist.begin(),eiglist.begin()+depth);

  int dimold = h->dimoldpart(subeiglist);

  stringstream subeiglist_ss;
  std::copy(subeiglist.begin(),subeiglist.end(),ostream_iterator<long>(subeiglist_ss," "));

  ECLOG(0) << "In form_finder, depth = " << depth 
           << ", aplist = [ " << subeiglist_ss.str() << "];\t"
           << "dimsofar=" << subdim
           << ", dimold=" << dimold
           << ", dimnew=" << subdim-dimold << "\n";
 
  if( dimold == subdim ) {
    data.setStatus( ALL_OLD );     // Set status of current node
    ECLOG(0) << "Abandoning a common eigenspace of dimension " << subdim
             << " which is a sum of oldclasses." << endl;
    return;   // This branch of the recursion ends: all is old
  }

  if( ( subdim == targetdim ) && ( depth > mindepth ) ) { 
    data.setStatus( FOUND_NEW );   // Set status of current node
    make_basis( data );
    store(data.bplus_,data.bminus_,subeiglist);
    return;
  }

  if( depth == maxdepth ) { 
    data.setStatus( MAX_DEPTH );
    if(1) {       // we want to see THIS message whatever the verbosity level! 
      cout << "\nFound a " << subdim << "D common eigenspace\n";
      cout << "Abandoning, even though oldforms only make up ";
      cout << dimold << "D of this." << endl;
    }
    return;
  }

  // Pass data node through to make_submat()
  // NOTE originally called in go_down(), but relocated here since
  // it only needs to be called once per node. 
  make_submat(data);

  // The recursive part:
  vector<long> t_eigs = h->eigrange(depth);
  auto apvar = t_eigs.begin();

  stringstream t_eigs_ss;
  std::copy(t_eigs.begin(),t_eigs.end(),ostream_iterator<long>(t_eigs_ss," "));

  ECLOG(0) << "Testing eigenvalues [ " << t_eigs_ss.str() 
           << "] at level " << (depth+1) << endl;

  // Set children counter
  data.setChildren( t_eigs );

  while( apvar != t_eigs.end() ) { 
    ECLOG(1) << "Going down with ap = " << (*apvar) <<endl;
    long eig = *apvar++;

    // Initiate new data node, passing a constant reference of the current 
    // form_finder object to the data class constructor
    ff_data<T> *child = new ff_data<T>( this );

    // Configure data node ancestry
    data.addChild( eig, *child );

#ifdef ECLIB_MULTITHREAD
    if( data.subdim_ > ECLIB_RECURSION_DIM_LIMIT ) {
      // Post newly created child node to threadpool
      pool.post< ff_data<T> >( *child );
    }
    else {
      // Parallel granularity control. Continue in serial.
      go_down( data, eig, apvar==t_eigs.end() );
      if( child -> subdim_ > 0 ) find( *child );
      //if( child -> status_ != INTERNAL || child -> subdim_ == 0 ) go_up( *child );
    }
#else   
    // Pass through current data node and new test eigenvalue to go_down()
    go_down( data, eig, apvar==t_eigs.end() );
    
    // We pass new child node to find() 
    if( child -> subdim_ > 0 ) find( *child );

    go_up( *child );
#endif
  }  

#ifndef ECLIB_MULTITHREAD
  ECLOG(0) << "Finished at level " << (depth+1) << endl;
#endif
}

template<class T>
void form_finderT<T>::store(Zvec<T> bp, Zvec<T> bm, vector<long> eigs) {
#ifdef ECLIB_MULTITHREAD
  // Lock function
  boost::mutex::scoped_lock lock( store_lock );
#endif

  // Store sub-bplus,bminus,eiglists in class level containers
  gbplus.push_back(bp);
  gbminus.push_back(bm);
  gaplist.push_back(eigs);

  // Increment global counter
  gnfcount++;

  // Inform about newform count
  ECLOG(1) << "Current newform subtotal count at " << gnfcount << endl;
}

// end of XSPLIT.CC
