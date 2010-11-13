// xsplit.cc: implementation of class form_finder
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2005 John Cremona
// 
// This file is part of the mwrank package.
// 
// mwrank is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// mwrank is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with mwrank; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include <unistd.h>  // for unlink() (not needed on linux)

#define USE_SPARSE 1
#include "xsplit.h"


// STORE_OPMATS flag: only relevant when bigmats=1.  Each opmat (on
// ambient space) is stored in a temporary disk file when first
// computed, then read in when needed later.  Otherwise recomputed
// each time -- saves on memory at the expense of speed

//#define STORE_OPMATS

#include "smatrix_elim.h"
subspace sparse_combine(const subspace& s1, const subspace& s2);
mat sparse_restrict(const mat& m, const subspace& s);
smat restrict(const smat& m, const subspace& s);

// CLASS FORM_FINDER (was called splitter)

#ifdef STORE_OPMATS
static char* opmatfilename()
{
  char* filename = new char[40];
  char* tmpmatdir = getenv("TMPMATDIR");
  if (tmpmatdir==NULL) 
    {
      sprintf(filename,"%s","/tmp/opmatXXXXXX");
    }
  else
    {
      sprintf(filename,"%s%s",tmpmatdir,"/opmatXXXXXX");
    }
  mkstemp(filename);
  delete tmpmatdir;
  //  cout << "opmatfilename() returns "<<filename<<endl;
  return filename;
}
#endif

form_finder::form_finder(splitter_base* hh, int plus, int maxd, int mind, int dualflag, int bigmatsflag, int v)
:h(hh), plusflag(plus), dual(dualflag), bigmats(bigmatsflag), verbose(v), 
 maxdepth(maxd), mindepth(mind)
{
  denom1 = h->matden();
  dimen  = h->matdim();
  depth=0;
  nest = new ssubspace*[maxd];
  //  nest[0] = new ssubspace(dimen); // this just wastes space, but we must now
                                  // treat depth 0 differently in go_down() 
  subdim=dimen;
  eiglist = vector<long>(maxd);
  submats = new smat[maxd];
#ifdef STORE_OPMATS
  havemat = new int[maxd]; int i=maxd; while(i--) havemat[i]=0;
  opfilenames=new char*[maxd+1];
  // filenames and files will be created as needed
#endif

  targetdim = 1;
  if(!plusflag) // full conjmat not needed when plusflag is true
    {
      targetdim=2;
      if(bigmats)
	{
	  conjmat=h->s_opmat(-1,dual); 
#ifdef STORE_OPMATS
	  opfilenames[maxd]=opmatfilename();
	  // cout<<"Writing conjmat to file " << opfilenames[maxd] << "\n";
	  conjmat.dump_to_file(opfilenames[maxd]);
	  conjmat.init(0,0); // release space until needed 
#endif
	}
    }
}

form_finder::~form_finder(void) 
{
  //  while(depth+1) delete nest[depth--]; 
#ifdef STORE_OPMATS
  long i;
  for(i=0; i<maxdepth; i++) 
    {
      if(havemat[i]) 
	{
	  unlink(opfilenames[i]);
	  delete[] opfilenames[i];
	}
    }
  if(!plusflag)
    {
      unlink(opfilenames[maxdepth]);
      delete[] opfilenames[maxdepth];
    }
  delete[] opfilenames;
  delete[] havemat; 
#endif
  while(depth) delete nest[depth--];   // nest[0] was not created!
  delete[] nest; 
  delete[] submats;
}


void form_finder::make_opmat(long i) 
{ 
#ifdef STORE_OPMATS
  if(!havemat[i]) 
    {
      the_opmat=h->s_opmat(i,dual,verbose); 
      opfilenames[i]=opmatfilename();
//    cout<<"Writing matrix "<<i<<" to file " << opfilenames[i] << "\n";
      the_opmat.dump_to_file(opfilenames[i]);
      havemat[i]=1;
    }
  else
    {
//    cout<<"Reading matrix "<<i<<" from file " << opfilenames[i] << "\n";
      the_opmat.read_from_file(opfilenames[i]);
    }
#else
  the_opmat=h->s_opmat(i,dual,verbose); 
#endif
}

void form_finder::make_submat()
{
  if(bigmats)
    {
      make_opmat(depth);  // fetch the_opmat from file, or compute
      if(depth==0)
	{
	  submats[depth]=the_opmat;
	}
      else
	{
	  if(verbose>1) cout<<"restricting the_opmat to subspace..."<<flush;
	  submats[depth] = restrict(the_opmat,*nest[depth]);
	  if(verbose>1) cout<<"done."<<endl;
	}
      the_opmat=smat(0,0); // releases its space
    }
  else
    {
      if(nrows(submats[depth])==0)
	{
	  if(depth==0)
	    submats[depth]= h->s_opmat(depth,1,verbose);
	  else
	    {
	      submats[depth] = h->s_opmat_restricted(depth,*nest[depth],1,verbose);
	    }
	}
    }
}

void form_finder::go_down(long eig, int last) 
{
  if(verbose>1)
    cout<<"Increasing depth to "<<depth+1<<", trying eig = "<<eig<<"..."<<flush;
  eiglist[depth]=eig;
  SCALAR eig2 = eig*denom1;
  if(verbose>1)
    cout<<"after scaling, eig =  "<<eig2<<"..."<<flush;
  //  if(depth) eig2*= denom(*nest[depth]); // else latter is 1 anyway
  ssubspace s(0);
  make_submat();

  if(verbose>1) cout<<"Using sparse elimination (size = "
                    <<dim(submats[depth])<<", density ="
		    <<density(submats[depth])<<")..."<<flush;
  if(verbose>3) cout<<"submat = "<<submats[depth]<<flush;

  s = eigenspace(submats[depth],eig2);
  if(((depth==0)&&(dim(s)>0)&&(nrows(submats[depth])>1000))||last)
    submats[depth]=smat(0,0); // to save space (will recompute when needed)
  if(verbose>1) cout<<"done (dim = "<<dim(s)<<"), combining subspaces..."<<flush;
  if(depth==0)
    nest[depth+1] = new ssubspace(s);
  else
    nest[depth+1] = new ssubspace(combine(*nest[depth],s));
  if(verbose>1) cout<<"done."<<endl;
  depth++;
  subdim = dim(*(nest[depth]));
  if(verbose>1) cout<<"Eigenvalue "<<eig<<" has multiplicity "<<subdim<<endl;
  if(verbose&&(subdim>0))
    cout<<" eig " << eig << " gives new subspace at depth "<<depth
        <<" of dimension "<<subdim<<endl;
}

void form_finder::go_up()
{
  if(depth>0) 
    {
      delete nest[depth]; 
      submats[depth]=smat(0,0); 
      depth--;
    }
  if(depth) subdim = dim(*nest[depth]);
  else      subdim = dimen;
}

void form_finder::make_basis()
{
  if(subdim!=targetdim)
    {
      cout<<"error in form_finder::make_basis with eiglist = ";
      for(int i=0; i<depth; i++) cout<<eiglist[i]<<",";
      cout<<"\nfinal subspace has dimension "<<subdim<<endl;
      cout<<"aborting this branch!"<<endl;
      return;
    }

  if(plusflag) 
    {
      if(depth==0) // must treat separately since we did not define nest[0]
                   // in order to save space
	{
	  bplus = vec(dimen); 
	  bplus[1]=1;
	}
      else 
	bplus = getbasis1(nest[depth]);
      return;
    }

  ssubspace* s = nest[depth];  // only used when depth>0
  ssubspace *spm0, *spm;
  SCALAR eig = denom1;
  //  if(depth) eig*=denom(*s);
  smat subconjmat;          // only used when depth>0
  if(bigmats)
    {
#ifdef STORE_OPMATS
      conjmat.read_from_file(opfilenames[maxdepth]);
#endif
      if(depth) subconjmat=restrict(conjmat,*s);
      else      subconjmat=conjmat;  
      // will only be a 2x2 in this case (genus 1 only!)
#ifdef STORE_OPMATS
      conjmat=smat(0,0);
#endif
    }
  else
    {
      subconjmat=h->s_opmat_restricted(-1,*s,1,verbose);
    }

  for(long signeig=+1; signeig>-2; signeig-=2)
  {
    SCALAR seig; seig=eig;
    if(signeig<0) {seig=-eig;}
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
        cout<<"error in form_finder::makebasis; ";
        cout<<"\nfinal (";
        if(signeig>0)cout<<"+"; else cout<<"-";
        cout<<") subspace has dimension "<<dim(*spm)<<endl;
        cout<<"aborting this branch!"<<endl;
	delete spm;
        return;
      }
    if(signeig>0)  bplus  = getbasis1(spm);
    else           bminus = getbasis1(spm);

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
	  cout<<"Form number "<<iform+1<<" with eigs ";
	  int n=eigs[iform].size(); if(n>10) n=10;
	  copy(eigs[iform].begin(), eigs[iform].begin() + n, 
	       ostream_iterator<long>(cout, " ")); 
	  cout<<"..."<<endl;
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

void form_finder::find()
{
   vector<long> subeiglist(eiglist.begin(),eiglist.begin()+depth);
   if (verbose) 
      cout << "In formfinder, depth = " << depth 
           << ", aplist = " << subeiglist << ";\t";
   int dimold = h->dimoldpart(subeiglist);
   if (verbose) cout<<"dimsofar="<<subdim
                    <<", dimold="<<dimold
                    <<", dimnew="<<subdim-dimold<<"\n";
   if (dimold==subdim)
     {
       if (verbose)
	 {
	   cout<<"Abandoning a common eigenspace of dimension "<<subdim;
	   cout<<" which is a sum of oldclasses."<<endl;
	 }
       return;   // This branch of the recursion ends: all is old
     }
   if ((subdim==targetdim) && (depth>mindepth))
     { 
       make_basis();
       h->use(bplus,bminus,subeiglist); 
       return;
     }
   if (depth==maxdepth)
     { 
       if (1) // we want to see THIS message whatever the verbosity level!
	 { 
	   cout << "\nFound a "<<subdim<<"D common eigenspace\n";
	   cout << "Abandoning, even though oldforms only make up ";
	   cout << dimold << "D of this." << endl;
	 }
       return;
     }
// The recursive part:
   vector<long> t_eigs = h->eigrange(depth);
   vector<long>::const_iterator apvar=t_eigs.begin();
   if(verbose) 
     cout<<"Testing eigenvalues "<<t_eigs<<" at level "<<(depth+1)<<endl;
   while(apvar!=t_eigs.end())
     { 
       //       cout<<"Going down with ap = "<<(*apvar)<<endl;
       long eig=*apvar++;
       go_down(eig,apvar==t_eigs.end());
       if(subdim>0) find();
       go_up();
     }
   if(verbose) cout<<"Finished at level "<<(depth+1)<<endl;
}

#if (METHOD==2)
subspace sparse_combine(const subspace& s1, const subspace& s2)
{
  // we assume s1, s2 are subspace mod BIGPRIME!
   scalar d=denom(s1)*denom(s2);
   const smat& sm1(basis(s1)), sm2(basis(s2));
   const mat& b = (sm1*sm2).as_mat(); 
   const vec& p = pivots(s1)[pivots(s2)];
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

smat restrict(const smat& m, const subspace& s)
{
  if(dim(s)==nrows(m)) return m; // trivial special case, s is whole space
  return mult_mod_p(m.select_rows(pivots(s)),smat(basis(s)),BIGPRIME);
}

#endif

// end of XSPLIT.CC
