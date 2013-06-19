// xsplit.cc: implementation of class form_finder
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
 
#include <unistd.h>  // for unlink() (not needed on linux)

#define USE_SPARSE 1
#include <eclib/xsplit.h>
#include <eclib/xsplit_data.h>

#include <eclib/smatrix_elim.h>
subspace sparse_combine(const subspace& s1, const subspace& s2);
mat sparse_restrict(const mat& m, const subspace& s);
smat restrict_mat(const smat& m, const subspace& s);

// CLASS FORM_FINDER (was called splitter)

//#define DEBUG_NEST

form_finder::form_finder(splitter_base* hh, int plus, int maxd, int mind, int dualflag, int bigmatsflag, int v)
:h(hh), plusflag(plus), dual(dualflag), bigmats(bigmatsflag), verbose(v), 
 maxdepth(maxd), mindepth(mind), gnfcount(0)
{
  denom1 = h->matden();
  dimen  = h->matdim();
 
  // Create and initialise new data object as root node
  root = new ff_data();

  // Set initial values
  // form_finder class is a friend of ff_data class
  // so may access private members
  r.subdim_ = dimen;
  r.eiglist_ = vector< long >(maxd);
  
  targetdim = 1;
  if( !plusflag ) {           // full conjmat not needed when plusflag is true
    targetdim=2;
    if( bigmats ) r.conjmat = h->s_opmat(-1,dual); 
  }
}

form_finder::~form_finder(void) {
  // Decendants of root node will be recursively deleted
  // if they have not already been deleted during find()
  // All dynamically created objects (subspaces) held 
  // in each data node will also be deleted.
  delete root;
}

void form_finder::make_opmat(long i, ff_data &data) { 
  data.the_opmat_ = h->s_opmat(i,dual,verbose); 
}

void form_finder::make_submat( ff_data &data ) {
  // Cache current data node depth
  long depth = data.depth_;

  if( bigmats ) { 
    // fetch the_opmat from file, or compute
    make_opmat(depth,data);
    
    if( depth == 0 ) {
	    data.submat_ = the_opmat;
	  }
    else {
	    if( verbose > 1 ) cout << "restricting the_opmat to subspace..." << flush;
#ifdef DEBUG_NEST
      cerr << "Accessing nest[" << depth << "]" << endl;
#endif
	    data.submat_ = restrict_mat(data.the_opmat_,data.nest_);
	    if( verbose > 1 ) cout << "done." << endl;
	  }
      
    data.the_opmat_ = smat(0,0); // releases its space
  }
  else {
    if( nrows(data.submat_) == 0 ) {
	    if( depth == 0 ) {  
        data.submats_ = h->s_opmat(depth,1,verbose);
	    }
      else {
#ifdef DEBUG_NEST
        cerr << "accessing nest[" << depth << "]" << endl;
#endif
	      data.submat_ = h->s_opmat_restricted(depth,data.nest_,1,verbose);
	    }
	  }
  }
}

/**
 * go_down()
 *
 * Initiates creation of new subspace; data stored in new
 * data node. Data node passed as parameter will become the
 * _parent_ of this new data node.
 */
void form_finder::go_down(long eig, int last, ff_data &data) {
  // Cache current depth
  long depth = data.depth_;

  if(verbose>1)
    cout << "Increasing depth to " << depth+1 << ", trying eig = " << eig << "..." << flush;
 
  // Initiate new data node
  ff_data *child = new ff_data();

  // Configure data node ancestry 
  data.addChild( child );
  child.parent_ = *data;

  // Set new depth
  child.depth_ = depth + 1;
  
  // Store new test eigenvalue in new node
  child.eigenvalue_ = eig;
  
  SCALAR eig2 = eig*denom1;
  if(verbose>1)
    cout << "after scaling, eig =  " << eig2 << "..." << flush;
  // if(depth) eig2*= denom(*nest[depth]); // else latter is 1 anyway
  ssubspace s(0);

  // Pass new data node through to make_submat()
  make_submat(child);

  if(verbose>1) cout << "Using sparse elimination (size = "
                     << dim(child.submat_) << ", density ="
		                 << density(child.submat_) << ")..." << flush;
  if(verbose>3) cout << "submat = " << child.submat_ << flush;

  s = eigenspace(child.submat_,eig2);

  // Save space (will recompute when needed)
  if(((depth==0)&&(dim(s)>0)&&(nrows(child.submat_)>1000))||last)
    child.submat_ = smat(0,0); 
     
  if(verbose>1) cout << "done (dim = " << dim(s) << "), combining subspaces..." << flush;
  
  if( depth == 0 ) {
#ifdef DEBUG_NEST
    cerr << "creating nest[" << (depth+1) << "]" << endl;
#endif
    child.nest_ = new ssubspace(s);
  }
  else {
#ifdef DEBUG_NEST
    cerr << "accessing nest[" << depth     << "]" << endl;
    cerr << "creating nest["  << (depth+1) << "]" << endl;
#endif
    child.nest_ = new ssubspace(combine( *(data.nest_),s ));
  }
  
  if(verbose>1) cout << "done." << endl;
  
  // Local depth increment (does not effect data nodes)
  depth++;

#ifdef DEBUG_NEST
  cerr << "accessing nest[" << depth << "]" << endl;
#endif
  child.subdim_ = dim( *(child.nest_) );
  
  if(verbose>1) {
    cout << "Eigenvalue " << eig 
         << " has multiplicity " << child.subdim_ << endl;
  }
  if(verbose && (child.subdim_>0)) {
    cout << " eig " << eig 
         << " gives new subspace at depth " << depth
         << " of dimension " << child.subdim_ << endl;
  }
}

void form_finder::go_up() {
  // May not be required anymore
  // Potentially change to return parent data node
  // of given node, and perform memory management.
}

void form_finder::make_basis( ff_data &data ) {
  // Cache data values
  long depth  = data.depth();
  long subdim = data.subdim();
  vector< long > eiglist = data.eiglist();

  if(subdim!=targetdim)
  {
    cout << "error in form_finder::make_basis with eiglist = ";
    for(int i=0; i<depth; i++) 
      cout << eiglist[i] << ",";
    cout << "\nfinal subspace has dimension " << subdim << endl;
    cout << "aborting this branch!" << endl;
    return;
  }

  if(plusflag) {
    // must treat separately since we did not
    // define nest[0] in order to save space
    if(depth==0) {
      data.bplus_    = vec(dimen); 
      data.bplus_[1] = 1;
	  }
    else {
#ifdef DEBUG_NEST
      cerr << "accessing nest[" << depth << "]" << endl;
#endif
      bplus = getbasis1(data.nest_);
    }
     
    return;
  }

#ifdef DEBUG_NEST
  cerr << "accessing nest[" << depth << "]" << endl;
#endif

  ssubspace* s = data.nest_;  // only used when depth>0
  ssubspace *spm0, *spm;
  SCALAR eig = denom1;
  //  if(depth) eig*=denom(*s);
  smat subconjmat;          // only used when depth>0
  if(bigmats)
  {
    if(depth) subconjmat = restrict_mat(data.conjmat_,*s);
    else      subconjmat = data.conjmat_;  
    // will only be a 2x2 in this case (genus 1 only!)
  }
  else
  {
    subconjmat=h->s_opmat_restricted(-1,*s,1,verbose);
  }

  for(long signeig=+1; signeig>-2; signeig-=2)
  {
    SCALAR seig; 
           seig = data.eigenvalue_;
    
    if(signeig<0) seig =- data.eigenvalue_;
    
    if(depth)
    {
	    spm0 = new ssubspace(eigenspace(subconjmat,seig));
	    spm = new ssubspace(combine(*s,*spm0));
	    delete spm0;
    }
    else 
        spm = new ssubspace(eigenspace(subconjmat,seig));
    
    if(dim(*spm)!=1)
    {
      cout << "error in form_finder::makebasis; ";
      cout << "\nfinal (";
      
      if(signeig>0) cout << "+"; 
      else cout << "-";
        
      cout << ") subspace has dimension " << dim(*spm) << endl;
      cout << "aborting this branch!" << endl;
	    delete spm;
      return;
    }
    
    if(signeig>0) data.bplus_  = getbasis1(spm);
    else          data.bminus_ = getbasis1(spm);

    delete spm;
  }
}

vec form_finder::getbasis1(const ssubspace* s)
{
  VEC b = basis(*s).as_mat().col(1);
#ifdef MODULAR
  if(!liftok(b,MODULUS)) 
     cout << "Unable to lift eigenvector from mod " << MODULUS << endl;
#else
  makeprimitive(b);
#endif
#ifdef MULTI
  scalar n=0; // dummy variable to gt the right type in next line
  return b.shorten(n);
#else      
  return b;
#endif
}

void form_finder::recover(vector< vector<long> > eigs)
{
  for(unsigned int iform=0; iform<eigs.size(); iform++)
  {
    if(verbose)
  	{
	    cout << "Form number " << iform+1 << " with eigs ";
	    
      int n = eigs[iform].size(); 
      if(n>10) n = 10;

	    copy(eigs[iform].begin(), eigs[iform].begin() + n, 
	       ostream_iterator<long>(cout, " ")); 
	    cout << "..." << endl;
	  }
 
    splitoff(eigs[iform]);
  }  
}

void form_finder::splitoff(const vector<long>& eigs)
{
  if(verbose)
    cout<<"Entering form_finder, depth = "<<depth<<", dimension "<<subdim<<endl;

// BACKTRACK:

  while(!startswith(eiglist,eigs,depth,0)) go_up();

// GO DEEPER, UNTIL DIMENSION = TARGET DIMENSION:

  if(verbose)
    cout<<"restarting at depth = "<<depth<<", dimension "<<subdim<<endl;

  while((subdim>targetdim) && (depth<maxdepth))
  {
    go_down(eigs[depth],1);
  }

  make_basis();
  h->use(bplus,bminus,eigs); 
  
  return;
}

void form_finder::find() {
  find( data );
  
  // Now compute all newforms only if recursion has finished
  cout << "Now performing use() on all lists at once" << endl;
  for( int nf = 0; nf < gnfcount; nf++ ) {
    h-> use(gbplus[nf],gbminus[nf],gaplist[nf]);
  }
}

void form_finder::find( ff_data &data )
{
  // Cache values of current data
  long depth  = data.depth();
  long subdim = data.subdim();
  vector< long > eiglist = data.eiglist();

  vector<long> subeiglist(eiglist.begin(),eiglist.begin()+depth);
 
  if(verbose) cout << "In formfinder, depth = " << depth 
                   << ", aplist = " << subeiglist << ";\t";

  int dimold = h->dimoldpart(subeiglist);
 
  if(verbose) cout << "dimsofar=" << subdim
                   << ", dimold=" << dimold
                   << ", dimnew=" << subdim-dimold << "\n";
 
  if( dimold == subdim ) {
    if(verbose) {
      cout << "Abandoning a common eigenspace of dimension " << subdim;
      cout << " which is a sum of oldclasses." << endl;
    }
    return;   // This branch of the recursion ends: all is old
  }

  if( ( subdim == targetdim ) && ( depth > mindepth ) ) { 
    make_basis( data );
    store(data.bplus_,data.bminus_,subeiglist);
    return;
  }

  if( depth == maxdepth ) { 
    if(1) {       // we want to see THIS message whatever the verbosity level! 
      cout << "\nFound a " << subdim << "D common eigenspace\n";
      cout << "Abandoning, even though oldforms only make up ";
      cout << dimold << "D of this." << endl;
    }
    return;
  }
 
  // The recursive part:
  vector<long> t_eigs = h->eigrange(depth);
  vector<long>::const_iterator apvar = t_eigs.begin();
 
  if(verbose) cout << "Testing eigenvalues " << t_eigs 
                   << " at level " << (depth+1) << endl;

  while( apvar != t_eigs.end() ) { 
    // cout << "Going down with ap = " << (*apvar) <<endl;
    long eig = *apvar++;

    go_down(eig,apvar==t_eigs.end(),data);
    
    // Only call following statements in serial version
    if( data.subdim_ > 0 ) find( data );
    go_up();
  }  
  
  // NOTE - Remove this statement with multithreading ...
  if(verbose) cout << "Finished at level " << (depth+1) << endl;
}

void form_finder::store(vec bp, vec bm, vector<long> eigs) {
  // Lock function
  boost::mutex::scoped_lock lock( store_lockhere );

  // Store sub-bplus,bminus,eiglists in class level containers
  gbplus.push_back(bp);
  gbminus.push_back(bm);
  gaplist.push_back(eigs);

  // Increment global counter
  gnfcount++;

  // Inform about newform count
  if(verbose) 
    cout << "Current newform subtotal count at " << gnfcount << endl;
}

#if (METHOD==2)
subspace sparse_combine(const subspace& s1, const subspace& s2)
{
  // we assume s1, s2 are subspace mod BIGPRIME!
   scalar d = denom(s1)*denom(s2);
   const smat& sm1(basis(s1)), sm2(basis(s2));
   const mat&  b = (sm1*sm2).as_mat(); 
   const vec&  p = pivots(s1)[pivots(s2)];
   return subspace(b,p,d);
   //  return COMBINE(s1,s2);
}

mat sparse_restrict(const mat& m, const subspace& s)
{
  if(dim(s)==nrows(m)) return m; // trivial special case, s is whole space
  scalar dd = denom(s);  // will be 1 if s is a mod-p subspace
  mat b(basis(s));
  smat sm(m), sb(b);
  vec piv=pivots(s);
  smat smr = sm.select_rows(piv);
  smat ans = smr*sb;
  int check=0;
  if(check) {
    smat left = sm*sb; 
    if(dd!=1) {cout<<"(dd="<<dd<<")"; left.mult_by_scalar_mod_p(dd);}
    smat right = sb*ans;
    int ok = eqmodp(left,right);
    if (!ok) 
    {
      cout<<"Warning from sparse_restrict: subspace not invariant!\n";
      cout<<"Difference = \n"<<left-right<<endl;
      abort();
    }
  }
  check=0;
  if(check) {
    int ok = (ans.as_mat()==RESTRICT(m,s));
    if (!ok) 
    {
      cout<<"Error in sparse_restrict: sparse result differs fromnormal!\n";
      abort();
    }
  }
  return ans.as_mat();
}

smat restrict_mat(const smat& m, const subspace& s)
{
  if(dim(s)==nrows(m)) return m; // trivial special case, s is whole space
  return mult_mod_p(m.select_rows(pivots(s)),smat(basis(s)),BIGPRIME);
}

#endif

// end of XSPLIT.CC
