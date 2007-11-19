// mrank2.cc: implementation of class rank2 for descent via 2-isogeny
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
 
#include <iomanip> // for setbase(), used for hex output of codes
#include "bitspace.h"
#include "compproc.h"
#include "points.h"
#include "mwprocs.h"
#include "mquartic.h"
#include "msoluble.h"
#include "descent.h"
#include "mrank2.h"
#include "sqfdiv.h"
#include "desc2.h"

#ifndef QSIEVE_OPT
#define QSIEVE_OPT 0 // uses Stoll's sieve
#endif

//#define DEBUG_ROOTS
//#define DEBUG_SQF
//#define DEBUG_ELS
//#define DEBUG_GLS

// bound on first descent search when a second descent is going to be done
#define FIRST_DESCENT_BOUND 8

#define MAX_R 5   // will not attempt to list all coset reps
                  // for 2E(Q) in E(Q) if rank is more than this

void rank2::makepoint(const bigint& c,const bigint& d1,const bigint& d2,
	       const bigint& x, const bigint& y, const bigint& z, 
	       int which)
{
  Point P(ee);
  if (verbose) cout<<" (x:y:z) = ("<<x<<":"<<y<<":"<<z<<")\n";
  if(which==0)
    {
      if(sign(z)!=0)	P.init(ee, d1*x*x*z, d1*x*y, pow(z,3));
      if(verbose) 
	{
	  cout<<"\tCurve E " /* <<(Curve)ee */ <<"\tPoint "<< P;
	  bigfloat h = height(P);
	  cout << ", height = " << h;
	  if(!P.isvalid()) {cout << " --warning: NOT on curve! " ;}
	  cout << "\n";
	}
    }
  else
    {
      if(verbose) // display point on isogenous curve too
	{
	  Point Q(eedash);
	  if(sign(z)!=0) Q.init(eedash, d1*x*x*z, d1*x*y, pow(z,3));
	  cout<<"\tCurve E' " /* <<(Curve)eedash */ <<"\tPoint "<< Q;
	  bigfloat h = height(Q);
	  cout << ", height = " << h;
	  if(!Q.isvalid()) {cout << " --warning: NOT on curve! " ;}
	  cout << "\n";
	}
      const bigint& xz=x*z, x2=x*x, z2=z*z;
      if(sign(xz)!=0) P.init(ee,2*y*y*xz,y*(d1*x2*x2-d2*z2*z2),pow(2*xz,3));
      if(verbose) 
	{
	  cout<<"\tCurve E " /* <<(Curve)ee */ <<"\tPoint "<< P;
	  bigfloat h = height(P);
	  cout << ", height = " << h;
	  if(!P.isvalid()) {cout << " --warning: NOT on curve! " ;}
	  cout << "\n";
	  
	}
    }
  if(order(P)<0)  {pointlist.push_back(P); npoints++;} // else torsion so ignore
}

int rank2::testquartic(const bigint& c,const bigint& d1,const bigint& d2,int which)
{
  // creates quartic (d1,0,c,0,d2), assumed els, and tries to find a rational point
  // returns +1 if rational point found (handled my makepoint())
  //          0 if undecided
  static const bigint zero = BIGINT(0);
  static const bigint one  = BIGINT(1);
  quartic q(d1,  zero, c,  zero, d2);
  if (verbose) cout<<q<<": ";

  bigint x,y,z;

  // First a quick search for a small point:
  if (ratpoint(q,one,BIGINT(lim1),x,y,z))
    { 
      makepoint(c,d1,d2,x,y,z,which);
      return 1;
    }

  // Next a more serious search using a sieve,
  // but not too far if we are going to use second descent
  quartic_sieve qs(&q,QSIEVE_OPT,0);
  long lim3=lim2; 
  if(do_second_descent)    // save the long search till later...
    if(lim3>FIRST_DESCENT_BOUND) lim3=FIRST_DESCENT_BOUND;
  if(qs.search(lim3, 1, 1)) // maxnpoints, pos-x-only
    { 
      qs.getpoint(x,y,z);
      makepoint(c,d1,d2,x,y,z,which);
      return 1;     
    }
  if (verbose) 
    cout<<" no rational point found (hlim="<<lim3<<")\n";
  return 0;
}  /* end of testquartic*/
 
int rank2::second_descent(const bigint& c, const bigint& d1, const bigint& d2, int which)
  // creates quartic (d1,0,c,0,d2), assumed els, and tries to find a rational point
  // by second descent (after testquartic() has failed using direct search)
  // returns +1 if rational point found (handled by makepoint())
  //          0 if undecided
  //         -1 if no rational point exists (proved by second descent)
{
  int res, verb=verbose; if(verb) verb--;  // reduced verbosity level within desc2()
  bigint x,y,z;
  if(verbose) cout<<"d1="<<d1<<": "<<flush;
  if(which)
    res = desc2(c,d1,d2,badprimes, supp0, elsgens0, mask0,  lim2,x,y,z,verb);
  else
    res = desc2(c,d1,d2,badprimes, supp1, elsgens1, mask1,  lim2,x,y,z,verb);

  if(verbose==1) cout<<endl;

  switch(res)
    {
    case 1: // Positive success, point found
      makepoint(c,d1,d2,x,y,z,which);
      if(verbose)
	cout<<"Second descent successfully found rational point for d1="<<d1<<"\n\n";
      break;
    case -1: // Negative success, no point exists
      if(verbose)
	cout<<"Second descent shows that no rational point exists for d1="<<d1<<"\n\n";
      break;
    case 0: // Undecided
      if (verbose) 
	cout<<"Second descent inconclusive for d1="<<d1<<": ELS descendents exist but no rational point found\n\n";
    }
  return res;
}

// First local descent: determine 
// S^(phi)(E')   if which==0
// S^(phi')(E)   if which==1

void rank2::find_elsgens(int which, const bigint& c, const bigint& d)
{
  static const bigint zero = BIGINT(0);
  if (verbose>1) 
    {
      if(which) cout<<"\n";
      cout<<"Finding els gens for E"; if(which) cout<<"'";
      cout<<" (c"; if(which) cout<<"'";
      cout<<"= "<<c<<", d"; if(which) cout<<"'";
      cout<<"= "<<d<<")"<<endl;
    }
  bigint ddash = c*c-4*d;
  int posd1, istep=1;
  if (sign(ddash)< 0) 
    posd1=1; 
  else 
    posd1= (c+sqrt(ddash)< 0);
  if(posd1) istep=2;  // this will skip negative d1
  
  int extra2torsion=0;
  bigint ee2, ee3;
  if(which)
    {
      if(d_is_sq)
	{
	  extra2torsion=1;
	  ee2 = e2dash;
	  ee3 = e3dash;
	}
    }
  else
    {
      if(ddash_is_sq)
	{
	  extra2torsion=1;
	  ee2 = e2;
	  ee3 = e3;
	}      
    }

  bigint d1, d2, badp;
  const vector<bigint>& supp = (which? supp1 : supp0);
  long ns = supp.size();
  if(verbose>1)
    {
      cout<<"Support (length "<<ns<<"): "<<supp<<endl;
    }
  if(ns>31) {cout<<"Too many primes dividing d!\n"; abort();}

  long mask=0, maxn=1<<ns, index, j, nelsgens=0;
  vector<bigint> elsgens;

// use all torsion: added 24/6/02
// Find and process torsion

  vector<Point> torsion;
  if(which)
    torsion=torsion_points(eedash); 
  else 
    torsion=torsion_points(ee); 
  long it,ntorsion=torsion.size();
#ifdef DEBUG_ELS
  cout<<"Number of torsion points = "<<ntorsion<<endl;
  cout<<"Torsion subgroup is ";
  if(extra2torsion) cout<<"not "; 
  cout<<"cyclic"<<endl;
#endif

  // find generators of torsion mod 2*torsion

  Point T1,T2; long orderT1=1;
  for(it=0; it<ntorsion; it++)
    {  
      Point T = torsion[it];
      long orderT = order(T);
      if(orderT>orderT1) {T1=T; orderT1=orderT;}
    }
#ifdef DEBUG_ELS
  cout<<"Generator of maximal order = "<<T1<<" of order "<<orderT1<<endl;
#endif
  if(extra2torsion) 
// then we need a second generator, for which we can take either of
// the 2-torsion points which is not the multiple of T1
    {
      T2=(ntorsion/4)*T1;
      bigint xT2=getX(T2);
      if(xT2==zero) 
	T2.init(ee,ee2,zero); 
      else 
	T2.init(ee,zero,zero);
#ifdef DEBUG_ELS
  cout<<"Second generator = "<<T2<<" of order 2"<<endl;
#endif
      ntorsion=2;
      torsion[0]=T1;
      torsion[1]=T2;
    }
  else 
    {
      ntorsion=1;
      torsion[0]=T1;
    }
#ifdef DEBUG_ELS
  cout<<"Processing generating torsion points"<<endl;  
#endif
  
  bigint d1x;
  for(it=0; it<ntorsion; it++)
    {
      Point T = torsion[it];
      if(T.iszero()) continue;
      d1=getX(T);
      if(d1==zero) d1=d;
#ifdef DEBUG_ELS
      if(verbose>1) cout<<"Processing torsion d1 = " << d1 << endl;
#endif
      index=makeindex(supp, d1, d1x);
      d1=d1x;  // the square-free part, = d1/square
#ifdef DEBUG_ELS
      if(verbose>1) cout<<"Squarefree part of d1 = " << d1 << endl;
#endif
      if(mask&index) continue;
      if(index)
	{
	  elsgens.push_back(d1); nelsgens++;
	  if(verbose>1) 
	    cout<<"Adding (torsion) els generator #"<<nelsgens<<": d1 = " << d1;
	  for(j=ns-1; j>=0; j--)
//	    if(testbit(index,j)) {setbit(mask,j); break;}
	    if(index&(1<<j)) {mask|=(1<<j); break;}
#ifdef DEBUG_ELS
	  if(verbose>1) cout<<" (pivot = "<<j<<": "<<supp[j]<<")";
#endif
	  if(verbose>1) cout<<endl;
	}
    }

// Record that the first 0, 1 or 2 elsgens come from known torsion points, 
// so can be ingored in the gls step later:
  if(which)
    {
      nt2gens1=nelsgens;
    }
  else
    {
      nt2gens0=nelsgens;
    }

// Now systematically go through square-free divisors d1 of d, always modulo 
// the subgroup of those which we know are els

  int res;
  for(index=istep; index<maxn; index+=istep)
    {
      if(mask&index) continue;
      d1 = makenum(supp,index);
      d2 = d/d1;
#ifdef DEBUG_ELS
      if(verbose>1) cout<<"Testing d1 = "<<d1<<":\t";
#endif
      res = locallysoluble(d1,c,d2,badprimes,badp);
      if(res)
	{
	  elsgens.push_back(d1); nelsgens++;
	  if(verbose>1) 
	    cout<<"Adding els generator #"<<nelsgens<<": d1 = " << d1;
	  for(j=ns-1; j>=0; j--)
//	    if(testbit(index,j)) {setbit(mask,j); break;}
	    if(index&(1<<j)) {mask|=(1<<j); break;}
#ifdef DEBUG_ELS
	  if(verbose>1) cout<<" (pivot = "<<j<<": "<<supp[j]<<")";
#endif
	  if(verbose>1) cout<<endl;
	}
      else 
	{
#ifdef DEBUG_ELS
	  if(verbose>1) cout<<"not locally soluble at p = "<<badp<<endl;
#endif
	}
    }
  if(verbose>1)
    {
      cout<<"After els sieving, nelsgens = " << nelsgens;
      cout << endl;

      cout<<"2-rank of S^{phi";      if(which) cout<<"'";
      cout<<"}(E";                   if(!which) cout<<"'";
      cout<<") = "<<nelsgens<<endl;

      if(nelsgens>0) cout<<"(els)gens: "<<elsgens<<endl;
      //      cout<<"mask = "<<mask<<endl;
    }
  if(which)
    {
      elsgens1=elsgens; mask1=mask; els1=nelsgens;
    }
  else
    {
      elsgens0=elsgens; mask0=mask; els0=nelsgens;
    }
}

// Second local descent: determine 
// phi'(S^2(E)) in S^(phi)(E')   if which==0
// phi(S^2(E')) in S^(phi')(E)   if which==1

void rank2::find_els2gens(int which, const bigint& c, const bigint& d)
{
  if (verbose>1) 
    {
      if(which) cout<<"\n";
      cout<<"Finding els2 gens for E"; if(which) cout<<"'";
      cout<<" (c"; if(which) cout<<"'";
      cout<<"= "<<c<<", d"; if(which) cout<<"'";
      cout<<"= "<<d<<") which lift to S^2(E"; 
      if(which) cout<<"'"; cout<<")"<<endl;
    }

  vector<bigint>& elsgens = (which? elsgens1: elsgens0);
  long nelsgens        = (which? els1: els0);
  long nt2gens         = (which? nt2gens1: nt2gens0);

  bigint d1, d2, badp, x,y,z;
  unsigned long els2mask; long index;
  long maxn=1<<nelsgens, nels2gens=0, els2piv;
  vector<bigint> els2gens;
  bitspace els2_space(nelsgens);

// first record the torsion contribution:

  for(index=0; index<nt2gens; index++)
    {
#ifdef DEBUG_ELS
      d1=elsgens[index];
      if(verbose>1) cout<<"Processing torsion d1 = " << d1 << ":"<<endl;
#endif
      els2mask=(1<<index);
      if(els2_space.mask(els2mask)) continue;  // we work mod the els2 subgp
      els2piv=els2_space.reduce(els2mask);
      if(els2piv<0) continue; // means we're in the subgp; won't happen
#ifndef DEBUG_ELS // else done earlier
      d1=elsgens[index];
#endif
      els2gens.push_back(d1); nels2gens++;
      els2_space.augment(els2mask,els2piv);
      if(verbose>1) cout<<"Adding (torsion) els2 generator #"<<(nels2gens)
			<<": d1 = " << d1 <<endl;
#ifdef DEBUG_ELS
      else
	cout<<"Just incremented nels2gens  to "<<nels2gens<<endl;
      cout<<"now bitmask = "<<els2_space.getbitmask()<<endl;
#endif	      
    } // end of torsion loop

  int res, verb=0; if(verbose>2)verb=verbose-2;
  // reduced verbosity level within desc2()

  for(index=1; (index<maxn)&&(nels2gens<nelsgens); index++)
    {
// First mask against known els2 subgroup:
      if(els2_space.mask(index)) continue;  // we work mod the els2 subgp  
      els2mask=index;
      els2piv=els2_space.reduce(els2mask);
#ifdef DEBUG_ELS
      cout<<"index = "<<index<<endl;
      cout<<"els2mask = "<<els2mask<<endl;
      cout<<"els2piv  = "<<els2piv <<endl;
#endif
      if(els2piv<0) continue; // means w're in the els2 subgp; won't happen   
      d1 = makenum(elsgens,index);
      d2 = d/d1;
#ifdef DEBUG_ELS
      if(verbose>1) cout<<"Processing d1 = " << d1 << ":\t";
#endif
      if(which)
	res = 1+desc2(c,d1,d2,badprimes,supp0,elsgens0,mask0,lim2,x,y,z,verb,1);
      else
	res = 1+desc2(c,d1,d2,badprimes,supp1,elsgens1,mask1,lim2,x,y,z,verb,1);
#ifdef DEBUG_ELS
      if(verbose>1) cout<<"res = " << res << endl;
#endif
      if(res)
	{
	  els2gens.push_back(d1); nels2gens++;
	  els2_space.augment(els2mask,els2piv);
#ifdef DEBUG_ELS
	  cout<<"now bitmask = "<<els2_space.getbitmask()<<endl;
#endif
	  if(verbose>1) 
	    cout<<"Adding els2 generator #"<<nels2gens<<": d1 = " << d1<<endl;
	}
    }

  if(verbose>1)
    {
      cout<<"After els2 sieving, nels2gens = " << nels2gens;
      cout << endl;

      cout<<"2-rank of phi";if(which) cout<<"'";
      cout<<"(S^2(E";       if(!which) cout<<"'";
      cout<<")) = "<<nels2gens<<endl;

      if(nels2gens>0) cout<<"(els2)gens: "<<els2gens<<endl;
    }
  if(which)
    {
      els2gens1=els2gens; els21=nels2gens;
    }
  else
    {
      els2gens0=els2gens; els20=nels2gens;
    }
}

void rank2::find_glsgens(int which, const bigint& c, const bigint& d)
{
  vector<bigint>& elsgens = (which? els2gens1: els2gens0);
  long nelsgens        = (which? els21: els20);
  long nt2gens         = (which? nt2gens1: nt2gens0);
  vector<bigint> gls_gens;
  bitspace gls_space(nelsgens);
  long glspiv, maxn = 1<<nelsgens;
  long nglsgens=0;
  unsigned long glsmask;
  long index, stage, nstages=1; if(do_second_descent) nstages=2;
  long shortfall1, shortfall2;
  bigint d1, d2;
  int res;
  
// first record the torsion contribution:

  for(index=0; index<nt2gens; index++)
    {
      glsmask=(1<<index);
      d1=elsgens[index];
      glspiv=gls_space.reduce(glsmask);
      if(glspiv<0) 
	{
#ifdef DEBUG_GLS
	  cout<<"d1="<<d1<<": known gls (torsion)"<<endl;
#endif	      
	  continue; // as we are certainly in the gls subgroup
	}
      if(verbose>1) 
	{
	  cout<<"Adding (torsion) gls generator #"<<(nglsgens+1)
	      <<": d1 = " << d1 <<endl;
	}
      gls_gens.push_back(d1); nglsgens++;
      gls_space.augment(glsmask,glspiv);
#ifdef DEBUG_GLS
      cout<<"Just incremented nglsgens to "<<nglsgens<<endl;
#endif	      
    } // end of torsion loop

// Next test all first descent curves for global solubility (stage 1)
// and (optionally) do a second descent (stage 2) on unsuccessful ones.

// The next lines are in case torsion accounts for everything, 
// as then the following loop is not executed at all: 
  shortfall1 = 0;  shortfall2 = 0;
  if(which)
    {gls1=gls21=nglsgens;}
  else
    {gls0=gls20=nglsgens;}

  if((nglsgens==nelsgens)&&verbose)
    cout<<"This component of the rank is 0\n";

  for(stage=1; (stage<=nstages)&&(nglsgens<nelsgens); stage++)
    {
      if(verbose&&do_second_descent) 
	{
	  if(stage==1) 
	    cout << "First stage (no second descent yet)...\n";
	  else        
	    cout << "Second stage (using second descent)...\n";
	}
      for(index=1; (index<maxn)&&(nglsgens<nelsgens); index++)
	{
#ifdef DEBUG_GLS
	  d1 = makenum(elsgens,index);  // do it now for output purposes
	  cout<<"index="<<setbase(16)<<index<<setbase(10);
	  cout<<"\td1 = "<<d1<<":\t";
#endif
// First mask against gls subgroup:
	  glsmask=index;
	  glspiv=gls_space.reduce(glsmask);
	  if(glspiv<0) 
	    {
#ifdef DEBUG_GLS
	      cout<<"known to be gls"<<endl;
#endif	      
	      continue; // as we are certainly in the gls subgroup
	    }
#ifdef DEBUG_GLS
	  cout<<"not known gls"
	      <<", mask="<<setbase(16)<<glsmask<<setbase(10)
	      <<endl;
	  cout<<"Keeping j="<<glspiv
	      <<" with mask "<<setbase(16)<<glsmask<<setbase(10)
	      <<" for possible use as pivot"<<endl;
#else
	  d1 = makenum(elsgens,index);  // else we did it earlier
#endif
	  d2 = d/d1;
	  if(stage==1) res = testquartic(c,d1,d2,which);
	  else         res = second_descent(c,d1,d2,which);
#ifdef DEBUG_GLS
	  cout<<"result = "<<res<<endl;
#endif	      
	  switch(res){
	  case -1: // should not happen
	    {
	      cout<<"Problem in 2nd descent!" <<endl;
	    }
	  case 0: // no point found, nothing to do
	    {
	      break;
	    }
	  case 1:
	    {
	      gls_space.augment(glsmask, glspiv);
	      gls_gens.push_back(d1);
	      nglsgens++; 
#ifdef DEBUG_GLS
	      cout<<"Just incremented nglsgens to "<<nglsgens<<endl;
#endif	      
	      if(verbose>1) 
		{
		  cout<<"Adding gls generator #"<<(nglsgens)
		      <<": d1 = " << d1;
#ifdef DEBUG_GLS
		  cout<<" (g_pos = "<<glspiv<<", g_gen = "
		      <<setbase(16)<<glsmask<<setbase(10)<<")";
#endif
		  cout<<endl;
		}
	      break;
	    }
	  } // end of switch(res)
	}   // end of index loop
#ifdef DEBUG_GLS
      cout<<"At bottom of index loop\n";
      cout<<"nglsgens = "<<nglsgens<<"\tgls_gens = ";
      cout << gls_gens << endl;
#endif   
      if(stage==1) // finished stage 1, i.e. first descent
	{
	  shortfall1 = nelsgens-nglsgens;
	  if(verbose)
	    {
	      cout << "After first global descent, this component of the rank";
	      if(shortfall1)
		{
		  cout << "\n\thas lower bound " << nglsgens-nt2gens 
		       << "\n\tand upper bound " << nelsgens-nt2gens 
		       << "\n\t(difference =   " << shortfall1 << ")\n";
		  if(nstages==2)
		    cout<<"Second descent will attempt to reduce this\n";
		}
	      else
		cout << " = " << nelsgens <<endl;
	      if(verbose>1) cout<<"\n";
	    }
	  if(which)
	      gls1=gls21=nglsgens;
	  else
	      gls0=gls20=nglsgens;
	}
      else // finished stage==2, i.e. 2nd descent
	{
	  shortfall2 = nelsgens-nglsgens;
	  if(verbose)
	    {
	      cout << "After second global descent, this component of the rank";
	      if(shortfall2)
		{
		  cout << "\n\thas lower bound " << nglsgens-nt2gens 
		       << "\n\tand upper bound " << nelsgens-nt2gens 
		       << "\n\t(difference =   " << shortfall2 << ")\n";
		}
	      else
		cout << " = " << nelsgens <<endl;
	      if(verbose>1) cout<<"\n";
	    }
	  if(which)
	    gls21=nglsgens;
	  else
	    gls20=nglsgens;
	}

    }  // end of stage loop
   
  if(which)
      glsgens1=gls_gens; 
  else
      glsgens0=gls_gens; 

  if(verbose>1)
    {
      cout<<"After gls sieving, nglsgens = " << nglsgens << endl;
      cout<< "shortfall in rank from first  descent = "<<shortfall1<<endl;
      if(do_second_descent)
	cout<< "shortfall in rank from second descent = "<<shortfall2<<endl;
      if(nglsgens>0) cout<<"(gls)gens: "<<gls_gens<<endl;
    }
}
 
void rank2::local_descent(const bigint& x0)
{
  bigint c,d,cdash,ddash,disc,rootd;
  const bigint zero = BIGINT(0), two=BIGINT(2), minusone=BIGINT(-1);

  c =  3 * x0 + s2;
  d = x0*(c+s2) + s4;
  cdash = - 2 * c;
  ddash = c*c -  4*d;
  disc = 2*d*ddash;
  if (is_zero(disc)) // this should have been caught by the calling program!
    {
      cout << "Curve is singular\n"; 
      success = 0;
      return;
    }

  if(verbose>1)
    {
      cout<<"(c,d)  =("<<c<<","<<d<<")"<<endl;
      cout<<"(c',d')=("<<cdash<<","<<ddash<<")"<<endl;
    }

  ee     = Curvedata(zero,  c,  zero,  d,  zero);
  eedash = Curvedata(zero,cdash,zero,ddash,zero);
  if(verbose) 
    cout<<"Using 2-isogenous curve "<<(Curve)(Curvedata(eedash,0))
	<<endl;

  supp0 = pdivs(d);
  supp1 = pdivs(ddash);
  badprimes.clear();
  set_union(supp0.begin(),supp0.end(),supp1.begin(),supp1.end(),back_inserter(badprimes));
  if(find(badprimes.begin(),badprimes.end(),two)==badprimes.end()) 
    badprimes.insert(badprimes.begin(),two);
  //  cout<<"supp0="<<supp0<<", supp1="<<supp1<<endl;
  //   cout<<"badprimes="<<badprimes<<endl;
  // NB -1 MUST be the first entry in the supports!
  supp0.insert(supp0.begin(),minusone);
  supp1.insert(supp1.begin(),minusone);
  //   cout<<"supp0="<<supp0<<", supp1="<<supp1<<endl;
  
  d_is_sq = ddash_is_sq = 0;
  if(isqrt(ddash,rootd))
  {
    ddash_is_sq=1;
    e2 = (-c-rootd)/2;
    e3 = (-c+rootd)/2;
  }
  if(isqrt(d,rootd))
  {
    d_is_sq=1;
    e2dash = c-2*rootd;
    e3dash = c+2*rootd;
  }

  if(verbose) 
    {
      cout<<"-------------------------------------------------------\n";
      cout<<"First step, determining 1st descent Selmer groups\n";
      cout<<"-------------------------------------------------------\n";
    }
  //  cout<<"Calling find_elsgens(0,...)"<<endl;
  find_elsgens(0,c,d);
  //  cout<<"Calling find_elsgens(1,...)"<<endl;
  find_elsgens(1,cdash,ddash);
  rank_bound = els0+els1-2;

  if(verbose) 
    {
      cout<<"After first local descent, rank bound = "<<rank_bound<<"\n";
      cout<<"rk(S^{phi}(E'))=\t"<<els0<<endl;
      cout<<"rk(S^{phi'}(E))=\t"<<els1<<endl;
      cout<<endl;
      cout<<"-------------------------------------------------------\n";
      cout<<"Second step, determining 2nd descent Selmer groups\n";
      cout<<"-------------------------------------------------------\n";
    }

  if((rank_bound==0)||(!do_second_descent))
    {
      els20=els0; els21=els1;
      els2gens0=elsgens0; els2gens1=elsgens1;
      if(verbose)
	{
	  if(do_second_descent)
	    cout<<"...skipping since we already know rank=0"<<endl;
	  else
	    cout<<"...skipping -- no second descent requested"<<endl;
	}      
    }
  else
    {
      find_els2gens(0,c,d);
      find_els2gens(1,cdash,ddash);
      rank_bound = els20+els21-2;
    }
  if(verbose)
    {
      if(do_second_descent) 
	{
	  cout<<"After second local descent, rank bound = "<<rank_bound<<"\n";
	  cout<<"rk(phi'(S^{2}(E)))=\t"<<els20<<endl;
	  cout<<"rk(phi(S^{2}(E')))=\t"<<els21<<endl;
	}
      cout<<"rk(S^{2}(E))=\t"<<els1+els20-1+ddash_is_sq<<endl;
      cout<<"rk(S^{2}(E'))=\t"<<els0+els21-1+d_is_sq<<endl;
      cout<<endl;
    }
}

rank2::rank2(Curvedata* ec, int verb, int sel, long l1, long l2, int second)
  : rank12(ec,verb,sel,l1,l2,0,second)
{
  bigint a1, a2, a3, a4, a6;
  ec->getai(a1,a2,a3,a4,a6);
  fullnpoints = npoints = 0;
  rank = 0;            // default value if failure occurs
  int best_isogeny=0, best_rank_bound=999999;

  success = 1;
  int n, scaled=0;
  if (odd(a1) || odd(a3)) 
    { s2= a1*a1+ 4*a2;
      s4=  8*(a1*a3+ 2*a4);
      s6= 16*(a3*a3+ 4*a6);
      scaled=1;
    }
  else 
    { s2=a2; s4=a4; s6=a6;
    }

  vector<bigint> xlist = Introotscubic(s2,s4,s6);
  ntwo_torsion = xlist.size();
  if (ntwo_torsion==0)
    { 
      success=0;
      if (verbose) cout << "No points of order 2\n";
      return;
    }

  if(verbose) cout << "\n" << ntwo_torsion << " points of order 2:\n";

  // If there are 3 points of order 2, we order them (for consistency:
  // otherwise the order can be machine dependent)
  if(ntwo_torsion==3) sort(xlist.begin(),xlist.end());

  two_torsion.resize(ntwo_torsion);
  for(n=0; n<ntwo_torsion; n++)
    {
      bigint ei = xlist[n];
      if(scaled) two_torsion[n].init(the_curve,2*ei,-a1*ei-4*a3,BIGINT(8));
      else two_torsion[n].init(the_curve,ei,BIGINT(0),BIGINT(1));
      if(verbose)
	{ if(n>0) cout<<", "; cout<<two_torsion[n];}
    }
  if(verbose)cout<<endl<<endl;

  long mw0,mw1,sel0,sel1,sha0,sha1,sel20,sel21,sha20,sha21,lindex;
  
 // we loop over the different 2-isogenies (if there are 3), doing the
 // local descent 3 times to get the best rank bound (or stopping if
 // that bound is 0).  

  // Fudge: we have to repeat the best one if it was not the last in
  // order to get the class-global variables right for the global
  // descent.  This should obviously be fixed.  

  // The relevant variables are:
  // ee,eedash,supp0,supp1,badprimes,d_is_sq,ddash_is_sq,
  // e2,e3,e2dash,e3dash,els20,els21,elsgens0,elsgens1

  for(n=0; n<ntwo_torsion; n++) 
    {
      if(verbose&&(ntwo_torsion>1)) 
	{
	  cout<<"****************************"<<endl;
	  cout<<"* Using 2-isogeny number "<<(n+1)<<" *"<<endl;
	  cout<<"****************************\n"<<endl;
	}
      local_descent(xlist[n]);
      if((n==0)||(rank_bound<=best_rank_bound)) 
	{
	  best_rank_bound=rank_bound;
	  best_isogeny=n;
	}
      if(best_rank_bound==0) break;
    }
  
  int rerun_needed=(rank_bound>best_rank_bound);
  rank_bound=best_rank_bound;
  selmer_rank=rank_bound; // for compatibility with full 2-descent

  if(verbose&&best_isogeny>0)
    {      
      cout<<"After second local descent, combined upper bound on rank = "
	  <<rank_bound<<"\n";
    }

  if(selmer_only) // Nothing more to do in this case
    {
      rank=0;
      certain=(rank_bound==0);
    }
  else // do the global stages...
    {

      // Redo the best local descent if necessary
      bigint x0=xlist[best_isogeny];
      if(rerun_needed) 
	{
	  if(verbose) cout<<"Rerunning the local descent for isogeny number "
			  <<(best_isogeny+1) 
			  << ", which gave the best upper bound on the rank"
			  <<endl;
	  local_descent(x0);
	}

  bigint c =  3 * x0 + s2;
  bigint d = x0*(c+s2) + s4;
  bigint cdash = - 2 * c;
  bigint ddash = c*c -  4*d;
  
  if(verbose)
    {
      cout<<"Third step, determining E(Q)/phi(E'(Q)) and E'(Q)/phi'(E(Q))\n";
      cout<<"-------------------------------------------------------\n";
      cout<<"1. E(Q)/phi(E'(Q))\n";
      cout<<"-------------------------------------------------------\n";
    }
  
  if(verbose)
    {
      cout<<"(c,d)  =("<<c<<","<<d<<")"<<endl;
      cout<<"(c',d')=("<<cdash<<","<<ddash<<")"<<endl;
    }

  find_glsgens(0,c,d);
  
  npoints1=npoints;
  if(verbose)
    {
      cout<<"-------------------------------------------------------\n";
      cout<<"2. E'(Q)/phi'(E(Q))\n";
      cout<<"-------------------------------------------------------\n";
    }
  find_glsgens(1,cdash,ddash);

//Debug only:
/*
cout<<"gls20 = "<<gls20<<endl;
cout<<"gls21 = "<<gls21<<endl;
cout<<"els0 = "<<els0<<endl;
cout<<"els1 = "<<els1<<endl;
cout<<"els20 = "<<els20<<endl;
cout<<"els21 = "<<els21<<endl;
*/

  rank        = gls20+gls21-2; // lower bound for rank
//rank_bound  = els20+els21-2; // upper bound for rank, from above

  mw0 = 1<<gls20;  
  sel0 = 1<<els0;    sha0=sel0/mw0;      
  sel20 = 1<<els20;  sha20=sel20/mw0;

// gls20  is the 2-rank of E/phi(E')                    (lower bound)
// els0   is the 2-rank of S^(phi)(E')                  (exact)
// els20  is the 2-rank of its subgroup phi'(S^(2)(E))  (exact if 2nd descent)
// els0-gls20  is the 2-rank of III(E')[phi)                 (upper bound)
// els20-gls20 is the 2-rank of its subgroup phi'(III(E)[2]) (upper bound)

  mw1 = 1<<gls21;  
  sel1 = 1<<els1;  sha1=sel1/mw1;
  sel21 = 1<<els21;  sha21=sel21/mw1;

// gls21   is the 2-rank of E'/phi'(E)                   (lower bound)
// els1    is the 2-rank of S^(phi')(E)                  (exact)
// els21   is the 2-rank of its subgroup phi(S^(2)(E'))  (exact if 2nd descent)
// els1-gls21  is the 2-rank of III(E)[phi']                 (upper bound)
// els21-gls21 is the 2-rank of its subgroup phi(III(E')[2]) (upper bound)

  index2=mw0*mw1;
  lindex=gls20+gls21;
  if(!ddash_is_sq) {index2/=2; lindex-=1;}
  
  int certain0 = (els20==gls20); // i.e. descent was conclusive
  int certain1 = (els21==gls21); //
  certain = certain0&&certain1;
  
  if (verbose) 
    {
      cout<<"\n";
      cout<<"-------------------------------------------------------\n";
      cout<<"Summary of results:\n";
      cout<<"-------------------------------------------------------\n";
      
      if(certain) cout << "\trank(E) " << "= " << rank;
      else        cout << "\t"<<rank<< " <= rank(E) <= " << rank_bound;
      cout<<"\n";
      cout << "\t#E(Q)/2E(Q) ";
      if(!certain) cout<<">";
      cout<<"= " << index2 << "\n\n";
      
      long lb, ub; int lb1, eq;

      if(do_second_descent)
	{
      cout<<"Information on III(E/Q):\n";

      lb = 1; lb1=1; ub = sha1; eq=(lb==ub);
      cout<<"\t#III(E/Q)[phi']    ";
      if(eq) cout<<"= "<<ub<<"\n";
      else
	{
	  if(lb1) cout<<"<= "<<ub<<"\n";
	  else    cout<<"is between "<<lb<<" and "<<ub<<"\n";
	}

      lb = sel1/sel21; lb1=(lb==1); ub = sha1*sha20; eq=(lb==ub);
      cout<<"\t#III(E/Q)[2]       ";
      if(eq) cout<<"= "<<ub<<"\n";
      else
	{
	  if(lb1) cout<<"<= "<<ub<<"\n";
	  else    cout<<"is between "<<lb<<" and "<<ub<<"\n";
	}
      cout<<endl;

      cout<<"Information on III(E'/Q):\n";

      lb = 1; lb1=1; ub = sha20; eq = (lb==ub);
      cout<<"\t#phi'(III(E/Q)[2]) ";
      if(!eq) cout<<"<";
      cout<<"= " << ub << "\n";

      lb = sel0/sel20; lb1=(lb==1); ub = sha0; eq=(lb==ub);
      cout<<"\t#III(E'/Q)[phi]    ";
      if(eq) cout<<"= "<<ub<<"\n";
      else
	{
	  if(lb1) cout<<"<= "<<ub<<"\n";
	  else    cout<<"is between "<<lb<<" and "<<ub<<"\n";
	}

      lb = sel0/sel20; lb1=(lb==1); ub = sha0*sha21; eq=(lb==ub);
      cout<<"\t#III(E'/Q)[2]      ";
      if(eq) cout<<"= "<<ub<<"\n";
      else
	{
	  if(lb1) cout<<"<= "<<ub<<"\n";
	  else    cout<<"is between "<<lb<<" and "<<ub<<"\n";
	}
      cout<<endl;
	}  // end of do_second_descent branch
      else
	{
      cout<<"Information on III(E/Q):\n";

      ub = sha1;
      if(ub==1) cout<<"\t III(E/Q)[phi'] is trivial\n";
      else cout<<"\t#III(E/Q)[phi'] <= " << ub << "\n";

      ub = sha0*sha1;
      if(ub==1) cout<<"\t III(E/Q)[2]    is trivial\n";
      else cout<<"\t#III(E/Q)[2]    <= " << ub << "\n";

      cout<<"Information on III(E'/Q):\n";

      ub = sha0;
      if(ub==1) cout<<"\t III(E'/Q)[phi] is trivial\n";
      else cout<<"\t#III(E'/Q)[phi] <= " << ub << "\n";

      ub = sha0*sha1;
      if(ub==1) cout<<"\t III(E'/Q)[2]   is trivial\n";
      else cout<<"\t#III(E'/Q)[2]   <= " << ub << "\n";

      cout<<endl;
	}  // end of no_second_descent branch

    }
  
  if(npoints>0) makegens();
  if(rank<=MAX_R)makepoints();
    }  // end of "else" clause after if(selmer_only)
}

void rank2::makegens()
{
  Curvedata ee_min;
  bigint u, r, s, t, x, y, z; int i;
  ee_min=ee.minimalize(u,r,s,t);
  if(verbose)
    {
      cout<<"-------------------------------------------------------\n";
      cout << "\nList of points on E = " << (Curve)ee_min << ":\n";
      cout<<"\nI.  Points on E mod phi(E')\n";
      if(npoints1==0) cout << "--none (modulo torsion).\n\n";
    }
  for(i=0; i<npoints; i++)
    {
      if(verbose&&(i==npoints1)) {cout<<"\nII. Points on phi(E') mod 2E\n";}
      Point q = shift(pointlist[i],the_curve,u,r,s,t);
      bigfloat h = height(q);
      int valid = q.isvalid();
      if(verbose||!valid) cout << "Point " << q << ", height = " << h;
      if(!valid) cout << " --warning: NOT on curve!";
      if(verbose||!valid) cout << "\n";
      pointlist[i]=q;
    }
  if(verbose&&(npoints1==npoints)) 
    {
      cout<<"\nII.  Points on phi(E') mod 2E\n";
      cout << "--none (modulo torsion).\n\n";
    }
}

void rank2::listgens()
{
  long i;
  cout << "List of generators of E(Q)/2E(Q) for E = " 
       << (Curve)(*the_curve) << ": \n";
  for(i=0; i<npoints; i++)
    {
      Point p = pointlist[i];
      cout << "Point " <<
	//	    "on " <<  (Curve)(*CD_orig) << ": " << 
	p;
      bigfloat h = height(p);
      cout << ", height = " << h;
      if(!p.isvalid()) cout << " --warning: NOT on curve!";
      cout << "\n";
    }
}

void rank2::listgens(Curvedata* CD_orig, const bigint& u, const bigint& r, 
		     const bigint& s, const bigint& t)
{
  long i;
  cout << "List of generators of E(Q)/2E(Q) (mod torsion) for E = " 
       << (Curve)(*CD_orig) << ": \n";
  for(i=0; i<npoints; i++)
    {
      Point p = shift(pointlist[i],CD_orig,u,r,s,t,1);
      cout << "Point " << (i+1) <<
	//	    "on " <<  (Curve)(*CD_orig) << 
	": " << p;
      bigfloat h = height(p);
      cout << ", height = " << h;
      if(!p.isvalid()) cout << " --warning: NOT on curve!";
      cout << "\n";
    }
}

void rank2::makepoints()
{
  if(fullnpoints>0) return;   // avoids calling this twice
  int i, j;
  long smallindex = index2/(1+ntwo_torsion);
  fullnpoints=1;         // will be smallindex
  fullpointlist.resize(smallindex);
  fullpointlist[0]=Point(the_curve);

  if(verbose&&(rank>0))
    {
      cout<<"-------------------------------------------------------\n";
      cout << "Computing full set of "<<smallindex<<" coset representatives for\n";
      cout << "2E(Q) in E(Q) (modulo torsion), and sorting into height order...." << flush;
    }
  for(i=0; i<rank; i++)
    {
      for(j=0; j<fullnpoints; j++)
	fullpointlist[j+fullnpoints] = fullpointlist[j]+pointlist[i];
      fullnpoints*=2;
    }
  if(fullnpoints!=smallindex)
    {
      cout << "Problem: index = " << index2 << " but " <<fullnpoints<<" cosets\n";
    }
//
// Now reorder points into increasing height order:
//
  for(i=0; i<fullnpoints; i++)
    for(j=i+1; j<fullnpoints; j++)
      if(height(fullpointlist[j])<height(fullpointlist[i]))
	{
	  Point temp = fullpointlist[i]; 
	  fullpointlist[i]=fullpointlist[j]; 
	  fullpointlist[j]=temp;
	}
  if(verbose&&(rank>0))  cout << "done.\n" << endl;
}

void rank2::listpoints()
{
  makepoints();
  cout << "Points on curve E = " << (Curve)(*the_curve) 
       << " covering E(Q)/2E(Q), modulo torsion:";
  if(rank==0) cout<<" none.";
  else
    if(rank>MAX_R)
      cout << "Too many to list ("<<fullnpoints<<" mod torsion)\n";
    else
    {
      cout<<"\n"<<fullnpoints<<" points, [0:1:0] and:\n";
      for (long i=1; i<fullnpoints; i++)
	{
	  Point p = fullpointlist[i];
	  cout << "Point " << p;
	  bigfloat h = height(p);
	  cout << ", height = " << h;
	  if(!p.isvalid()) {cout << " --warning: NOT on curve! " ;}
	  cout << "\n";
	}
    }
  cout<<"\n\n";
}

void rank2::listpoints(Curvedata* CD_orig, const bigint& u, const bigint& r, 
		                           const bigint& s, const bigint& t)
{
  makepoints();
  cout << "Points on original curve E = " << (Curve)(*CD_orig) 
       << " covering E(Q)/2E(Q), modulo torsion:";
  if(rank==0) cout<<" none.";
  else if(rank>MAX_R) 
    cout << "Too many to list ("<<fullnpoints<<" mod torsion)\n";
  else
    {
      cout<<"\n"<<fullnpoints<<" points:"<<endl;
      for (long i=1; i<fullnpoints; i++)
	{
	  Point p0 = fullpointlist[i];
//	  cout << "Point on " <<  (Curve)(*the_curve) << ": " << p0;
//	  if(!p0.isvalid()) cout << " --warning: NOT on curve!";
//	  cout << "\n";
	  Point p = shift(p0,CD_orig,u,r,s,t,1);
	  cout << "Point " <<
	    //	    "on " <<  (Curve)(*CD_orig) << ": " << 
	    p;
	  bigfloat h = height(p);
	  cout << ", height = " << h;
	  if(!p.isvalid()) cout << " --warning: NOT on curve!";
	  cout << "\n";
	}
    }
  cout<<"\n\n";
}
