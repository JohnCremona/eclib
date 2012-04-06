// mrank1.h -- implementation of class rank1 for general 2-descent
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
 
#include "points.h"   // from qcurves library
#include "mwprocs.h"   // from qcurves library
#include "mquartic.h"
#include "mequiv.h"
#include "msoluble.h"
#include "qc.h"
#define USE_BIGINTS
#include "descent.h"
#include "mrank1.h"
#include "twoadic.h"

#ifdef USE_BIGINTS
//#define DEFAULT_NAUX 12
#define DEFAULT_NAUX 8
#else
#define DEFAULT_NAUX 5
#endif

// For testing: defines the strategy for dealing with large quartics
// 1 for original BSD criteria
// 2 for simple version of JC+MS criteria (Lemma 5.1 only)
// 3 for optimal JC+MS criteria, with 2-adic refinement (uses twoadic.h)

// So far, the options only affect the criteria for using small
// quartics only, not the handling of large quartics using the exact
// 2-adic index.

// 4 for intelligent handling of large quartics when 2-adic index is 2
// or 4.  When index=2 we can abort large quartic search as soon as
// one is found; when index=2 we can abort after two are found
// provided that they are independent modulo small quartics.

#ifndef LARGE_Q
#define LARGE_Q 4
#endif

#ifndef QSIEVE_OPT
#define QSIEVE_OPT 0 // uses Stoll's sieve
#endif


#define SQUARE_A_FIRST
//#define NO_PADIC_FILTERING

//#define SHOW_ABC_RANGES
//#define DEBUG_AH

#define maxnquartics 2048 // OK for curves of rank 12 or less.
#define abceps 0.001      // used in abc-test

#define ROUNDADJUST 0.001 // Amount add/subtracted before rounding for safety
long roundtemp;
inline int FITS_IN_LONG(const bigfloat& x) {return (x<=MAXLONG)&&(x>=MINLONG);}
inline void ROUNDUP(long& a,const bigfloat& x) 
{
  if(FITS_IN_LONG(x)) return Iasb(a,ceil((x)-ROUNDADJUST));
  cout<<"Attempt to round "<<x<<" to a long int fails, aborting!\n";
  abort();
}
inline void ROUNDDOWN(long& a,const bigfloat& x) 
{
  if(FITS_IN_LONG(x)) return Iasb(a,floor((x)+ROUNDADJUST));
  cout<<"Attempt to round "<<x<<" to a long int fails, aborting!\n";
  abort();
}

#define NEQPLIST 5        // Number of primes for equiv-test sieving

#ifndef USE_BIGINTS
int xsqrt(bigfloat a, bigfloat &b)
{
  if(a<0) return 0;
  b=floor(sqrt(a)+0.1);
  return is_zero(b*b-a);
}
#endif

long * rank1::qeps(const quartic& q, int x2)
{
  long* vec = new long[num_aux];  // position 0 is not used
  long i; vec[0]=0;
  for(i=1; i<num_aux; i++)
    {
      long a = posmod(q.geta(),auxs[i]);
      long H = posmod(q.getH(),auxs[i]);
      if(x2) H=posmod(hscalemod[i]*H,auxs[i]);
      vec[i] = flags[i][a][H];
    }
  return vec;
}

void rank1::show_eps_vec(long * vec)
{
  long i;
  cout<<"(";
  for(i=1; i<num_aux; i++) {
    if(i>1) cout<<":";
    if(aux_types[i]==1)
      switch(vec[i]) {
      case 15: cout<<"0"; break;
      case  5: cout<<"1"; break;
      default: cout<<"?";
    }
    else
      switch(vec[i]) {
      case 15: cout<<"00"; break;
      case  5: cout<<"01"; break;
      case  3: cout<<"10"; break;
      case  1: cout<<"11"; break;
      default: cout<<"??";
    }
    //    cout<<"("<<vec[i]<<")";
  }
  cout<<")";
}

 // process latest quartic found
void rank1::addquartic(const bigint& a, const bigint& b, const bigint& c, 
		       const bigint& d, const bigint& e)
{
  long firsti, i, oldnumber=0, thisnumber, nfl; 
  char ab; unsigned long code;  
  int trivial=0, newone=1, gls=0, els=0;
  quartic * qlist, *thisq;;
  bigint x,y,z, badp; 	  Point Ptemp;
  int btype = 0;    
  int pivtype=-1; // set to 0 for \infty, 1 for odd prime, 2 for 2
  if (type==1) // then we have an egg point, i.e. \infty is pivotal
    {btype=1; pivtype=0;}
#if LARGE_Q>3
  if(!btype)
    {
      if (extra2&&(twoadic_index==2)) // then 2 is pivotal
	{btype=1; pivtype=2;}
    }
#endif
  if(!btype)
    {
      if (ipivot>=0)  // then we have a pivotal odd prime
	{btype=1; pivtype=1;}
    }
  int atype = !btype;

  // NB We do not use 2 as a pivotal prime when the 2-adic index is 4,
  // since we have not implemented the corresponding local map from
  // quartics to (Z/2Z)^2 which maps large quartics to one of the
  // three non-trivial elements.  So when the 2-adic index is 4, we
  // only use large quartics for pivoting when there is an odd pivotal
  // prime suitable.

  if(atype) 
    {qlist=qlista; thisnumber=nquarticsa; nfl=nfirstlota; ab='A';}
  else 
    {qlist=qlistb; thisnumber=nquarticsb; nfl=nfirstlotb; ab='B';}

  qlist[thisnumber].assign(a,b,c,d,e,croots,type,ii,jj,disc);
  thisq=qlist+thisnumber;

  if (verbose) cout << (*thisq) << "\t";
  if (verbose>1) 
    {
      cout << "(ipivot = "<<ipivot<<", type = "<<ab<<") \t";  
      long * vec = qeps(*thisq,extra2);
      show_eps_vec(vec);
      cout<<"\t";
      delete[] vec;
    }

  // Check triviality

  if(atype) trivial = thisq->trivial(); // else certainly nontrivial
  if (trivial)
    {
      if (verbose) cout << "--trivial"<<endl;
      return;
    }
  if (verbose) cout<<"--nontrivial..."<<flush;

  // Check current is inequivalent to previous

  code = thisq->set_equiv_code(eqplist);
  firsti = (extra2==1 ? nfl : 0);
  for (i=firsti; newone && (i<thisnumber); i++)
    { 
      if((!atype)&&(!qlistbflag[i])) continue;
      if(traceequiv) 
	cout << "\nTesting equiv with number " << ab<< i+1 << endl;
#ifdef NEW_EQUIV
      newone = ! new_equiv(thisq,qlist+i,traceequiv);
#else
      newone = !     equiv(thisq,qlist+i,dlist,traceequiv);
#endif
      if (!newone) oldnumber=i+1;
    }

  if(newone)      // Check local and global solubility:
    {
      if(atype) 
	{
	  nquarticsa++;  
	  if(verbose) cout<<"--new (A) #"<<nquarticsa<<"\t"<<flush;
	}
      // but we do not increment nquarticsb unless the quartic is
      // els with no point, as otherwise we do not need to keep it

      if(selmer_only)
	{
	  gls = ratpoint(*thisq,BIGINT(1),BIGINT(lim1),x,y,z);
	  if(gls) els=1; 
	  else els=locallysoluble(*thisq,plist,badp);
          if(verbose) 
	    {
	      if(els) cout<<"locally soluble\n";
	      else cout<<"not locally soluble (p = "<<badp<<")\n";
	    }
	}
      else

      if (ratpoint(*thisq,BIGINT(1),BIGINT(lim1),x,y,z))
	{ 
	  gls=els=1;
	  if (verbose) cout<<"(x:y:z) = ("<<x<<" : "<<y<<" : "<<z<<")\n";
	}
      else
	{
// 	  cout<<"\nChecking "<<(*thisq)<<" for local solubility at "<<plist<<endl;
	  if (locallysoluble(*thisq,plist,badp))
	    { 
	      els=1;
	      if (verbose) cout<<"locally soluble..."<<flush;
	      quartic_sieve qs(thisq,QSIEVE_OPT,0);
	      if(qs.search(lim2))
		{ 
		  qs.getpoint(x,y,z); gls=1;
		  if (verbose) 
		    cout<<"(x:y:z) = ("<<x<<" : "<<y<<" : "<<z<<")\n";
		}
	      else 
		if (verbose) 
		  cout<<"no rational point found (limit "<<lim2<<")"<<flush;
	    }
	  else 
	    if (verbose) 
	      cout<<"not locally soluble (p = "<<badp<<")\n";
	}

      if(gls||(els&&selmer_only))
	{
	  if(!selmer_only)
	    {
  //cout<<"Calling qc() with (x:y:z) = ("<<x<<" : "<<y<<" : "<<z<<")\n";  
	      qc(*thisq,x,y,z,the_curve,&IJ_curve,tr_u,tr_r,tr_s,tr_t, Ptemp,verbose);
  //cout<<"qc() returns giving point " << Ptemp << "\n";
	    }
	  if(atype) // we have a quartic in A = ker(eps)
	    {
	      if(!selmer_only)
		{
		  pointlist1.push_back(Ptemp);
		  npoints1++; 
		  n1++; 
		}
	      n2++;
	      if(verbose) 
		{
		  if(selmer_only) cout<<"Selmer rank increases to "<<n2<<endl;
		  else cout<<"Size of A=ker(eps) increases to "<<n1<<endl;
		}
	    }
	  else // we have a quartic in B = im(eps)
	       // pivtype tells which type of prime is pivotal
	    {
	      if(!selmer_only)
		{
		  pointlist2.push_back(Ptemp);
		  npoints2++; 
		}
	      switch(pivtype) {
	      case 2:  // this is a large quartic
		{
		  global_index*=2;
		  if(verbose)
		    {
		      cout<<"Doubling global 2-adic index to "<<global_index<<endl;
		    }
		  if(global_index==twoadic_index) 
		    {
		      if(verbose)
			{
			  cout<<"global 2-adic index is equal to local index\n";
			  cout<<"so we abort the search for large quartics"<<endl;
			}
		      have_large_quartics=1;
// now go through qlistb to see if any of them were large 
// and now redundant, by comparing discriminants
		      
		      for(i=0; i<nquarticsb; i++)
			{
			  if(!qlistbflag[i]) continue;  
			  // i'th is already redundant
			  qlistbflag[i]=(qlistb[i].getdisc()!=disc);
			  if((!qlistbflag[i])&&verbose) 
			    cout<<"Quartic B #"<<(i+1)<<" is now redundant\n";
			}
		    }
		  break;
		} // end of case 2
	      case 0: //  on the egg
		{
		  have_eggpoint=1;
// now go through qlistb to see if any of them were Type 1
// and now redundant
		  for(i=0; i<nquarticsb; i++)
		    {
		      if(!qlistbflag[i]) continue;  
		                   // i'th is already redundant
		      qlistbflag[i]=(qlistb[i].gettype()!=1);
		      if(!qlistbflag[i]) 
			if(verbose)
			  cout<<"Quartic B #"<<(i+1)<<" is now redundant\n";
		    }
		  break;
		} // end of case 0:
	      case 1: //odd pivotal prime # ipivot
		{
		  int oldflag=aux_flags[ipivot], newflag=8;
		  if((oldflag==1)&&(pivflag==5)) newflag=2;
		  if((oldflag==1)&&(pivflag==3)) newflag=4;
		  if((oldflag==1)&&(pivflag==1)) newflag=4;
		  aux_flags[ipivot]=newflag;
		  if(verbose>1) cout<<"\nipivot = "<<ipivot
				    <<", changing mask from "<<oldflag
				    <<" to "<<newflag<<endl;
 // go back through qlistb to see if any of them would now be sieved out, 
 // as they are now redundant
		  long auxpiv = auxs[ipivot];
		  long hscale = hscalemod[ipivot];
		  for(i=0; i<nquarticsb; i++)
		    {
		      if(!qlistbflag[i]) continue;  // i'th already redundant
		      // compute (a,h) of i'th quartic:
		      long a = posmod(qlist[i].geta(),auxpiv);	      
		      long H = posmod(qlist[i].getH(),auxpiv);
		      if(extra2) if(i>=nfl) H=posmod(hscale*H,auxpiv);
		      
		      // Check if it would now be sieved out:
		      long fl = flags[ipivot][a][H];
		      if(verbose>1) 
			cout<<(i+1)<<"-th quartic in list has flag = "<<fl<<endl;
		      qlistbflag[i] = (fl & newflag)!=0;
		      if(!qlistbflag[i]) 
			if(verbose)
			  cout<<"Quartic B #"<<(i+1)<<" is now redundant\n";
		    }
		} // end of case 1
	      } // end of switch between cases
	      rank_B++; 
	      if(verbose)
		{
		  if(selmer_only)
		    cout<<"Selmer rank increases to "<<rank_B<<endl;
		  else
		    {
		      cout<<"Rank of B=im(eps) increases to "<<rank_B;
		      if(type==1) 
			cout<<" (The previous point is on the egg)";
		      else {if(verbose>1) 
			cout<<" (pivotal prime =" << auxs[ipivot] << ")";}
		      cout<<endl;
		    }
		}
 // in this case we do not need to keep the quartic,  so do not increment nquarticsb
	    } // end of B=im(eps) case
	}  // end of gls case
      else // no rational point was found (& we are not doing selmer_only)
	if(els) // have a possible "Selmer point"
	  {
	    if(atype) 
	      {
		n2++;
		if(verbose) cout<<endl;
	      }
	    else
	      {
		qlistbflag[nquarticsb]=1;
		nquarticsb++;
		if(verbose) cout<<" --new (B) #"<<nquarticsb<<endl;
	      }
	  }
      return;
    }
  else
    {
      if (verbose) cout << "--equivalent to ("<<ab<<") #"<<oldnumber<<endl;
    }
} // end of addquartic()


void rank1::getquartics()
{
  nquarticsa = 0; nfirstlota = 0; 
  nquarticsb = 0; nfirstlotb = 0; 
  have_eggpoint = 0;
  have_large_quartics = 0;
  rank_B = 0;
  ah_count = ah_sieve_0 = ah_sieve_1 = ah_sieve_2 = 0;
  ah_rfail = ah_dfail   = ah_efail   = ah_extra2fail = ah_pass = 0;
  ii=c4; jj=2*c6;   disc = 4*d1728;

  if(div(16,ii)&&div(64,jj)) 
    {
      ii/=16; jj/=64; disc/=4096;
      tr_u/=2; tr_r/=4; tr_s/=2; tr_t/=8;
    }
  if(verbose)
    {
      cout<<"Basic pair: I="<<ii<<", J="<<jj<<endl;
      cout<<"disc="<<disc<<endl;
    }

  xii=I2bigfloat(ii); xjj=I2bigfloat(jj);

  long Imod4 = posmod(ii,4), Jmod4 = posmod(jj,4);
  npairs=2; // default
  global_index=1;  // only gets increased when large quartics are found
  twoadic_index = 2;
  if(Jmod4==0) if((Imod4==2)||(Imod4==3)) twoadic_index=4;
  if(verbose) 
    cout<<"2-adic index bound = "<<twoadic_index<<endl;

  bsd_npairs = 2;
  if (div(4,ii)&&div(8,jj)&&div(16,2*ii+jj))
    {
      bsd_npairs=1;
    }
#if LARGE_Q==1
  npairs=twoadic_index=bsd_npairs;
#else // LARGE_Q>1
  if( (Imod4==0)&&(Jmod4==0) )
    { 
      if  (div(16,2*ii+jj)||div(16,2*ii+jj-4))
// Case covered by Lemma 5.1(a)!	
	{
	  npairs=twoadic_index=1;
	  if(verbose)
	    {
	      cout<<"By Lemma 5.1(a), 2-adic index = ";
	      if(npairs==1) cout<<"1\n"; else cout<<"2\n";
	    }
	}
#if LARGE_Q>2
      else // use 2-adic refinement to determine index (case 1)
	{
	  bigint a = -27*ii/4;
	  bigint b = -27*jj/4;
	  if(verbose>1) 
	    cout<<"Case 1 with a = I/4 = "<<a<<", b = J/4 = "<<b<<endl;
	  twoadic_index = npairs = 1 + case1(a,b);
	  if(verbose)
	    {
	      cout<<"After 2-adic refinement (case 1); 2-adic index = ";
	      if(npairs==1) cout<<"1\n"; else cout<<"2\n";
	    }
	}
#endif
    }
  if( (Imod4==1)&&(Jmod4==2) )
    {
      if  (div(16,ii+jj+5)||div(16,ii+jj+1))
// Case covered by Lemma 5.1(b)!	
	{
	  npairs=twoadic_index=1;
	  if(verbose)
	    {
	      cout<<"By Lemma 5.1(b), 2-adic index = ";
	      if(npairs==1) cout<<"1\n"; else cout<<"2\n";
	    }
	}
#if LARGE_Q>2
      else  // use 2-adic refinement to determine index (case 2)
	{
	  bigint a = -7-27*(ii-1)/4;
	  bigint b = -14-27*(jj-2)/4;
	  if(verbose>1)
	    cout<<"Case 2 with a = (I-1)/4 = "<<a<<", b = (J-2)/4 = "<<b<<endl;
	  twoadic_index = npairs = 1 + case2(a,b);
	  if(verbose)
	    {
	      cout<<"After 2-adic refinement (case 2); 2-adic index = ";
	      if(npairs==1) cout<<"1\n"; else cout<<"2\n";
	    }
	}
#endif
    }
  if (verbose) 
    {
#if LARGE_Q>2
      cout<<"2-adic index = "<<twoadic_index<<endl;
#endif
      if(npairs==2) cout<<"Two (I,J) pairs";
      else cout<<"One (I,J) pair";
      cout<<endl;
    }
  
  if (verbose) 
    {
      if (div(4,ii)&&div(8,jj)&&div(16,2*ii+jj)) // BSD say 1 pair
	{
	  if(npairs==2) // then the new result is worse -- should not happen
	    {
	      cout<<"!!! BSD give one (I,J) pair "
		  <<"-- this should NOT happen"<<endl;
	    }
	}
      else // BSD say 2 pairs
	{
	  if(npairs==1) // then the new result is an improvement
	    {
	      cout<<"*** BSD give two (I,J) pairs"<<endl;
	    }
	}
    }
#endif // LARGE_Q==1 or >1

  vector<bigint> plist0 = getbad_primes(*the_curve); // sorted by construction
  // now make sure 2 and 3 are in the list of primes
  vector<bigint> p23; p23.push_back(BIGINT(2)); p23.push_back(BIGINT(3));
  set_union(plist0.begin(),plist0.end(),p23.begin(),p23.end(),back_inserter(plist));
//   cout<<"\nplist0 = "<<plist0<<", p23="<<p23<<endl;
//   cout<<"\tplist = "<<plist<<endl;

#ifndef NEW_EQUIV
  dlist = sqdivs(disc,plist);
#endif
  threediv = div(3,ii);

  aux_init();
  extra2=0;

  //  bigcomplex c1(0), c2(-3*xii), c3(xjj);
  //  cphi = solvecubic( c1, c2, c3);
  bigcomplex w =  bigcomplex(to_bigfloat(-1), sqrt(to_bigfloat(3)))/to_bigfloat(2);
  bigfloat one_third = to_bigfloat(1)/to_bigfloat(3);
  cphi = new bigcomplex[3];
  if(is_zero(ii))
    {
      if(xjj>0) 
	cphi[2]=-exp(one_third*log(xjj));
      else 
	cphi[2]=exp(one_third*log(-xjj));
      cphi[1] = w*cphi[2];
      cphi[0] = conj(cphi[1]);
    }
  else
    {
      bigfloat xmdisc = I2bigfloat(-disc);
      bigcomplex t;
      if(posdisc)
	{
	  bigcomplex t3 = sqrt(bigcomplex(xmdisc))-xjj;
	  t3 /= to_bigfloat(2);
	  t  = to_bigfloat(3)*exp(one_third*log(t3));
	}
      else
	{
	  bigfloat t3 = sqrt(xmdisc)-xjj;
	  t3 /= to_bigfloat(2);
	  if(t3>0)
	    t  = to_bigfloat(3)*exp(one_third*log(t3));
	  else
	    t  = -to_bigfloat(3)*exp(one_third*log(-t3));
	}
      cphi[2] = (t+9*xii/t)*one_third; // real when disc<0
      t*=w;
      cphi[1] = (t+9*xii/t)*one_third;
      t*=w;
      cphi[0] = (t+9*xii/t)*one_third;
    }

  Imod2=odd(ii);
  Jmod2=odd(jj);
  if(verbose>1)
    cout<<"Before sorting, phi = "<<cphi[0]<<","<<cphi[1]<<","<<cphi[2]<<endl;
  flag_init();
  getquartics1();
  if (npairs==2) 
    {
      nfirstlota = nquarticsa;
      nfirstlotb = nquarticsb;
      extra2=1;  ii*= 16;  jj*= 64;  disc*= 4096; xii*=16; xjj*=64;
      tr_u*=2; tr_r*=4; tr_s*=2; tr_t*=8;
      Imod2=Jmod2=0;
#ifndef NEW_EQUIV
      dlist = sqdivs(disc,plist);
#endif
      bigfloat four = to_bigfloat(4);
      cphi[0]*=four; cphi[1]*=four; cphi[2]*=four;
      getquartics1();
    }
  delete[] cphi;
  if(verbose>1) 
    {
      cout << ah_count      << "\t (a,b,c) triples in search region\n";
      cout << ah_sieve_1    << "\t failed c-divisiblity,\n";
      cout << ah_sieve_2    << "\t failed syzygy sieve,\n";
      cout << ah_sieve_0    << "\t passed sieve.\n";
      cout << ah_rfail      << "\t failed syzygy after sieving,\n";
      cout << ah_dfail      << "\t failed d-integrality,\n";
      cout << ah_efail      << "\t failed e-integrality,\n";
      cout << ah_extra2fail << "\t failed extra-2 divisibility conditions,\n";
      cout << ah_pass       << "\t passed all and produced quartics.\n";
    }
  clear_sieve(); // clears memory allocated by aux_init()
  return;
} // End of getquartics()

void rank1::getquartics1()
{
  if (verbose) 
    cout<<"Looking for quartics with I = "<< ii << ", J = " << jj << endl;
  
  static bigint zero = BIGINT(0);
  IJ_curve = Curvedata(zero,zero,zero,-27*ii,-27*jj,0);  // don't minimise

  if (posdisc) 
    {
      gettype(2);  // get type 2s first as they are a subgroup of index 1 or 2
      if(!have_eggpoint) 
	gettype(1);
    }  
  else 
    {
      gettype(3);
    }
}  //  of getquartics1()


void rank1::gettype(int t) // new hybrid version 13/2/96
{
  type=t;
  long a,astep, amin=0, amax=0, firsta, lasta;
  long b,bstep;
  long c,cstep,cmod3;  
  int a_is_odd, b_is_odd, a_div_by_4;
  int a_positive;
  bigint I48=48*ii, J64=64*jj;
  static const bigint m27=BIGINT(-27);
  static const bigfloat root27=sqrt(to_bigfloat(27));
  static const bigfloat zero=to_bigfloat(0);
  bigint rsq, r, rem, h, d, e, ee;
  bigfloat r1, r2, r3, xr;
  long efactor, cfactor;
  bigcomplex c1;
  // Unnecessary initializations to keep -Wall happy:
  bigfloat phi=zero,phi1,phi2,phi3; 
  bigfloat amax0=zero, amin0, amax2, amin2, amax3, amin3;
  bigfloat hmin=zero, hmax=zero, hmin0=zero, hmax0=zero, hmin2, hmax2, htemp;
  bigfloat const6=zero,const5=zero,const3=zero,const2=zero;
  

  int extraextra2 = div(64,ii)&&div(128,jj);  // Pascale's extra condition
  cfactor = 1;
  cmod3 = mod(jj,3);

  if (verbose) cout << "Looking for Type " << t << " quartics:\n";

// Set phi to be the real root in type 3, 
// else set the phi_i to be the three real roots 
//   in descending order phi1 > phi2 > phi3:

  switch(type) {
  case 1:
    phi1 = real(cphi[0]);
    phi2 = real(cphi[1]);
    phi3 = real(cphi[2]);
    orderreal(phi1,phi2,phi3);  // decreasing order
    hmax0 = 4*(phi3*phi3-xii)/3;
    if(verbose>1)
      cout<<"phi1 = "<<phi1<<"\nphi2 = "<<phi2<<"\nphi3 = "<<phi3<<"\nhmax0 = "<<hmax0<<endl;
    break;
  case 2:  // N.B. repetition here to allow types 1 and 2 in either order
    phi1 = real(cphi[0]);
    phi2 = real(cphi[1]);
    phi3 = real(cphi[2]);
    orderreal(phi1,phi2,phi3);  // decreasing order
    hmin0 = 4*(phi2*phi2-xii)/3;
    break;
  case 3:
    // find the real phi and make it phi[2]
    if (is_real(cphi[1])) 
      {
	phi=real(cphi[1]);
	cphi[1]=cphi[2];
	cphi[2]=phi;
      }
    else 
      {
	if (is_real(cphi[2])) 
	  {
	    phi=real(cphi[2]);
	  }
	else 
	  {
	    if (is_real(cphi[0])) 
	      {
		phi=real(cphi[0]);
		cphi[0]=cphi[2];
		cphi[2]=phi;
	      }
	    else // error, non are detected to be real!
	      {
		cout<<"ERROR: none are real, quitting"<<endl;
		abort();
	      }
	  }
      }
    const2 = phi*phi-4*xii;
    const5 = phi*phi-xii;
    const6 = 0;
    if(const2>0) const6 = sqrt(const2);  // const2 must be >0 but avoid rounding problems
    hmin0 = 4*const5/3;
    amax0 = (abs(cphi[0]-cphi[1]) + 2*abs(cphi[0]-cphi[2]))/18;
    // (2*sqrt(phi^2-I)+sqrt(phi^2-4I))/6*sqrt(3)
    if(verbose>1)
      {
	cout<<"After  sorting, phi = "<<cphi[0]<<","
	    << cphi[1]<<","<<cphi[2]
	    << "\nBasic a bound = " << amax0 <<"\n";
      }
  } // end of switch(type)


// Set bounds on a loop:

  switch(type) {
  case 1:
    amin = 1;
    ROUNDDOWN(amax,(phi1-phi3) / 9);
    break;
  case 2:
    ROUNDUP(amin,(phi3-phi2)/9);    // negative
    ROUNDDOWN(amax,(phi1-phi2)/9);  // positive
    break;
  case 3:
    if(phi>0)
      {
	amax2 = (phi+const6)/6;    // const6 = sqrt(phi^2-4I)
	amax3 = 2*const5/(9*phi);  // const5 = phi^2-I
	if(verbose>1)
	  {
	if((amax2<amax0)&&(amax0<amax3))  // then we use amax0 but old version used amax2
	  {
	    cout<<"a upper bound = "<<amax0<<" while old version wrongly had "<<amax2<<endl;
	  }
	if((amax2<amax3)&&(amax3<amax0))  // then we use amax3 but old version used amax2
	  {
	    cout<<"a upper bound = "<<amax3<<" while old version wrongly had "<<amax2<<endl;
	  }
	  }
	if(amax3>amax2) amax2=amax3;  // so amax2=max of two previous
	if(amax2>amax0) 
	  {
	    amax2=amax0;              // so amax2=min of two previous
	  }
	ROUNDDOWN(amax,amax2);

	amin2= -const6/root27;
	if(amin2<-amax0)
	  {
	    amin2=-amax0;
	    cout<<"New a lower bound worse, not using -- should NOT happen!\n";
	  }
	ROUNDUP(amin,amin2);
      }
    else // phi<0
      {
	amax2= const6/root27;
	if(amax2>amax0)
	  {
	    amax2=amax0;
	    cout<<"New a upper bound worse, not using -- should NOT happen!\n";
	  }
	ROUNDDOWN(amax,amax2);

	amin0 = -amax0;         // this and the two following are all negative
	amin2 = (phi-const6)/6;
	amin3 = 2*const5/(9*phi);  // const5 = phi^2-I
	if(verbose>1)
	  {
	if((amin0<amin3)&&(amin3<amin2))  // then we use amin3 but old version used amin2
	  {
	    cout<<"a lower bound = "<<amin3<<" while old version wrongly had "<<amin2<<endl;
	  }
	if((amin3<amin0)&&(amin0<amin2))  // then we use amin0 but old version used amin2
	  {
	    cout<<"a lower bound = "<<amin0<<" while old version wrongly had "<<amin2<<endl;
	  }
	  }
	if(amin2>amin3) amin2=amin3;  // so amin2=min of two previous
	if(amin2<amin0)
	  {
	    amin2=amin0;              // so amin2=max of two previous
	  }
	ROUNDUP(amin,amin2);
      }
    if(verbose>1) cout<<"Search range for a: ("<<amin<<","<<amax<<")\n";
  } // end of switch(type)

  // DEBUG:    amin=-19; amax=-19;
  if(verbose>1) cout<<"Search range for a: ("<<amin<<","<<amax<<")\n";

  //negative first:
  //  for(a_positive=0; a_positive<=1; a_positive++)
  //positive first:
  //    for(a_positive=1; a_positive>=0; a_positive--)
#ifdef SQUARE_A_FIRST
  for(a_positive=2; a_positive>=0; a_positive--)
#else
  for(a_positive=1; a_positive>=0; a_positive--)
#endif
    {
#ifdef SQUARE_A_FIRST
      int square_a_only = (a_positive==2);
#endif
      if (a_positive) 
	{
	  firsta=amin; if(firsta<1) firsta=1;
	  lasta =amax;
	  if(firsta>lasta) continue;
	  a=firsta-1;
	  astep=1;
	  if(verbose) 
	    {
	      cout << "Trying positive a from " << firsta << " up to " << lasta;
#ifdef SQUARE_A_FIRST
	      if(square_a_only)
		cout << " (square a first...)";
	      else
		cout << " (...then non-square a)";
#endif
	      cout<<  endl;
	    }
	}
      else  // negative range will be traversed downwards!
	{
	  firsta=amax; if(firsta>-1) firsta=-1;
	  lasta =amin;
	  if(firsta<lasta) continue;
	  a=firsta+1;
	  astep=-1;
	  if(verbose) 
	    cout << "Trying negative a from " << firsta << " down to " << lasta << endl;
	}

#ifdef ABORT_LARGE
      // TEMPORARY CODE INCLUDED ONLY FOR RUNNING MARK WATKINS'S 2
      // MILLION RANK 2 CURVES!

      if(abs(lasta)>10000)
	{
	  if(verbose)
	    cout<<"**************\na range too big, quitting\n*************\n";
	  success=0;
	  return;
	}
#endif

  long iaux; 
  long *amodi, *hmodi, *auxi, *hstepmodi, *hscalemodi, *astepmodi;

  int ***flagsi; int **flagai;
  iaux=num_aux; amodi=amod; auxi=auxs;  astepmodi=astepmod;
  while(iaux--) 
    {
      *amodi++     = posmod(a,    *auxi); 
      *astepmodi++ = posmod(astep,*auxi);
      auxi++;
    }

  while(a!=lasta)
    {
      a+=astep;
      for(iaux=0, amodi=amod, auxi=auxs, flagai=flaga, flagsi=flags, 
	  astepmodi=astepmod;
	  iaux<num_aux; 
	  iaux++, amodi++, auxi++, flagai++, flagsi++, astepmodi++) 
	{
	  (*amodi)+=(*astepmodi); 
	  if((*amodi)>=(*auxi)) (*amodi)-=(*auxi);
	  *flagai = (*flagsi)[*amodi];
	}

#ifdef SQUARE_A_FIRST
      // First time through positive a, only look at square a:
      if(a_positive)
	{
	  long roota=(long)(sqrt((double)a)+0.1);
	  if((square_a_only)!=(a==roota*roota))  continue;
	}
#endif

//
// Tests: not(4|a) if extra2, not(2|a) if extraextra2:
//
      a_is_odd = (a&1);
      a_div_by_4 = !(a&3);
      if ((! (extra2 && (a_div_by_4))) && (! (extraextra2 && !(a_is_odd)))   )
	{
#ifdef SHOW_ABC_RANGES
	  if(verbose)cout<<"a = "<<a<<endl;
#endif
	  efactor = a_is_odd ? 16 : 8;  // only relevant in extra2 case
	  long cfac2 = (extra2?(a_is_odd ? 2 : 4):1);
	  long cfac3 = 3;
	  cstep   = cfac3*cfac2;
	  int b_must_be_odd =  (!extra2) && ((Jmod2 && !a_is_odd) || a_div_by_4);
// latter because else (a,b,c,d,e)==(a/4,b/2,c,2d,4e)
	  bstep = 1;
	  if(extra2) bstep=4; else
	    if(b_must_be_odd) bstep=2;

	  long absa=abs(a); 
	  long absa2=absa<<1;
	  bigfloat xa=to_bigfloat(a);
	  bigfloat xa4 = 4*xa, xa8=8*xa;
	  bigfloat oneover4a = 1/xa4;
	  
// hstep does not depend on b so can be set here:
	  long hstep = 8*a*cstep;
	  iaux=num_aux; hstepmodi=hstepmod; auxi=auxs; hscalemodi=hscalemod;
	  if(extra2)
	    while(iaux--) 
	      {
		*hstepmodi++ = posmod(hstep*(*hscalemodi),*auxi); 
		auxi++; hscalemodi++;
	      }
	  else
	    while(iaux--) 
	      {
		*hstepmodi++ = posmod(hstep,*auxi); 
		auxi++;
	      }

// Set up bounds on H

	  switch(type) {
	  case 1:
	    hmin = xa4*phi2;
	    //if(verbose>1) cout<<"hmin = "<<hmin<<endl;
	    hmax = xa4*phi3 + hmax0;
	    //if(verbose>1) cout<<"hmax1 = "<<hmax<<endl;
	    hmax2 = xa4*phi1;
	    //if(verbose>1) cout<<"hmax2 = "<<hmax2<<endl;
	    if(hmax>hmax2) hmax = hmax2;
	    //if(verbose>1) cout<<"hmax = "<<hmax<<endl;
	    break;
	  case 2:  //  NB these are really the min/max of sgn(a)*H
	    if(a_positive)
	      {
		hmin = xa4*phi2 + hmin0;
		hmax = xa4*phi3;
	      }
	    else
	      {
		hmin = xa4*phi1;
		hmax = xa4*phi2 + hmin0;
	      }
	    break;
	  case 3:
	    {
	      hmax = xa4*phi;
	      const3 = const2*(4*const5-27*xa*xa);
	      if(const3<0) 
		{
		  //		  cout<<"const3 = "<<const3<<endl;
		  const3=0; 
		}
	      else const3=2*sqrt(const3)/3;
	      //	      cout<<"const3 =  "<<const3 << endl;
	      hmax2 = -hmax/2 + const3;
	      //	      cout<<"hmax2 =  "<<hmax2 << endl;

	      hmin = hmax - hmin0;
	      hmin2 = -hmax/2 - const3;

	      if(hmin2>hmin) hmin=hmin2;
	      if(hmax2<hmax) hmax=hmax2;  // don't move this line up!

	      if(!a_positive) // swap the bounds over:
	      {
		htemp=hmin; hmin=hmax; hmax=htemp;
	      }
	    }
	  }

#ifdef SHOW_ABC_RANGES
	  if(verbose>1)
	    cout<<"(a="<<a<<")\tSearch range for H: ("<<hmin<<","<<hmax<<")\n";
#endif
	  if(a_positive)
	    {
	      if(hmin>hmax) 
		{
//	  cout<<"Empty H-range!  hmin = "<<hmin<<", hmax = "<<hmax<<"\n";
		  continue;
		} // skip to next a
	    }
	  else
	    {
	      if(hmax>hmin) 
		{
//	  cout<<"Empty H-range!  hmin = "<<hmax<<", hmax = "<<hmin<<"\n";
		  continue;
		} // skip to next a
	    }

// Set up b loop:  first value used is 0 or 1 after incrementing

	  b=b_must_be_odd-bstep;

	  while(b<=(absa2-bstep))
	    {
	      b+=bstep;
#ifdef SHOW_ABC_RANGES
	      if(verbose)cout<<"\tb = "<<b;
#endif
	      b_is_odd = (b&1);
// bstep handles these conditions
/*	      if (
		  ((!Jmod2) || a_is_odd || b_is_odd)
		  && ((!extra2) || !(b&3))
		  && (b_is_odd||(!a_div_by_4))
		  )
*/
		{
		  int doboth = ((0<b) && (b<absa2));
		  bigfloat xb=to_bigfloat(b);
		  bigfloat xbb3=3*xb*xb;

// set up bounds on c loop: (avoiding rounding error)

#ifdef SHOW_ABC_RANGES
		  cout<<"Before rounding, cmin = "<<(((hmin+xbb3)/xa8))<<endl;
#endif
		  long cmin; ROUNDUP(cmin,((hmin + xbb3)/xa8));
#ifdef SHOW_ABC_RANGES
		  cout<<"Before rounding, cmax = "<<(((hmax+xbb3)/xa8))<<endl;
#endif
		  long cmax; ROUNDDOWN(cmax,((hmax + xbb3)/xa8));

		  if(cmin>cmax+1)
		    cout<<"Empty c-range!  cmin = "<<cmin<<", cmax = "<<cmax<<"\n";
		  //DEBUG: cmin = 804861;

		  while(mod(cmin,3)!=cmod3) 
		    cmin++; // Skip to correct residue mod 3
		  while(cmin%cfac2)       
		    cmin+=3; // Skip to next multiple of cfac2

		  if(cmin>cmax) continue;   // Skip to next b
#ifdef SHOW_ABC_RANGES
		  if(verbose)cout<<":\tcmin = "<<cmin<<", cmax = "<<cmax<<" (cstep = " << cstep << ")\n";
#endif

		  c=cmin-cstep; // So first value used is cmin

		  iaux=num_aux; hmodi=hmod; auxi=auxs; hscalemodi=hscalemod;
		  while(iaux--) 
		    {
		      long aux=(*auxi);
		      long  cmod  = c%aux, bmod  = b%aux;
		      long bb3mod = (3*bmod*bmod)%aux;
		      long h0 = (((8*a*cmod)%aux)-bb3mod)%aux; 
		      if(extra2) h0=posmod((*hscalemodi)*h0,aux);
		      else       h0=posmod(h0,aux);
		      *hmodi++ = h0;
		      hscalemodi++;
		      auxi++;
		    }
		  
		  while(c<=cmax-cstep)
		    {
		      ah_count++;
		      c+=cstep;
		      for(iaux=0, hmodi=hmod, hstepmodi=hstepmod, auxi=auxs; 
			  iaux<num_aux; 
			  iaux++, hmodi++, hstepmodi++, auxi++) 
			{
			  (*hmodi) += (*hstepmodi);
			  if((*hmodi)>=(*auxi)) (*hmodi)-=(*auxi);
			}

		      int flagok = (b_is_odd || even(c-Imod2));
		      if(!flagok) {ah_sieve_1++; continue;}
		      
		      ipivot=-1;

		      for(iaux=0, hmodi=hmod, flagai=flaga; 
			  flagok&&(iaux<num_aux); 
			  iaux++, hmodi++, flagai++)
			{
			  int thisflag = (*flagai)[*hmodi];
			  flagok = aux_flags[iaux] & thisflag;
#ifndef NO_PADIC_FILTERING
			  if(flagok&&(iaux>0)&&(!(thisflag&8))&&(ipivot==-1))
			    {
			      ipivot=iaux; pivflag=thisflag;
			    }
#endif
			}

		      if(!flagok) {ah_sieve_2++; continue;}

// We have an (a,b,c)-triple which passes the sieve test

		      ah_sieve_0++;
		      bigfloat xh=8*xa*c-xbb3;
		      int ok=1;

// Check that rounding has not put us outside the range:
		      if((a_positive&&((xh>hmax)||(xh<hmin)))
			 ||
			 ((!a_positive)&&((xh>hmin)||(xh<hmax))))
			{
			  if(verbose>1)
			    cout<<"(a,b,c)=("<<a<<","<<b<<","<<c<<"): "
			      <<"H = "<<xh<<" is outside range "
				<<hmin<<"..."<<hmax<<"\n";
			  ok=0;
			}
		      if(!ok) continue;

		      bigint biga=BIGINT(a);
		      bigint biga8=biga<<3;
		      bigint bigb=BIGINT(b);
		      bigint bigbsq=sqr(bigb);
		      bigint bigbb3=3*bigbsq;
		      h = biga8*c-bigbb3;

#ifdef USE_BIGINTS
// use bigints from now on
		      bigint asq=sqr(biga);
		      bigint cub = h*(sqr(h)-asq*I48)+biga*asq*J64;
		      ok = ::divides(cub,m27,rsq,rem);
		      if(!ok) 
			{
			  cout<<"cub not divisible by 27 for (a,b,c,h)=("
			    <<a<<","<<b<<","<<c<<","<<h<<")\n";
			  cout<<"cub       = "<<cub<<endl;
			  ah_rfail++;
			  continue;
			}
		      //cout<<"isqrt "<<rsq<<"\n";
		      if(!isqrt(rsq,r)) {ah_rfail++; continue;}
		      xr=I2bigfloat(r);
#else
		      bigfloat asq=xa*xa;
		      bigfloat xrsq=(xh*(xh*xh-48*asq*xii)+64*xa*asq*xjj)/(-27);
		      if(!xsqrt(xrsq,xr)) {ah_rfail++; continue;}
#endif
#ifdef DEBUG_AH
		      cout<<"; r = "<<r<<" "<<flush;
#endif	      
		       bigint bigc = BIGINT(c); 
		       bigint bigcsq = bigc*bigc; 
		       bigint ii_cc = ii-bigcsq;  
#ifdef USE_BIGINTS
		       bigint temp = bigb*(bigbsq-4*biga*bigc);
		       // must compute as bigints
#else
		       bigfloat xc=c, xb=b;
		       bigfloat xcc=xc*xc, xbb=xb*xb;
		       bigfloat temp = xb*(xbb-xa4*xc);  
		       // must compute as bigfloats
#endif
// Loop on  sign of b:
		      long ib, sb=b, sb3=3*b;
		      for(ib=0; ib<1+doboth; ib++)
			{if(ib) {sb=-sb; sb3=-sb3; temp=-temp;}
#ifdef DEBUG_AH
			 cout<<"\na,b,b3,c,temp = "<<a<<","<<sb<<","
			   <<sb3<<","<<c<<","<<temp<<flush;
#endif
#ifdef USE_BIGINTS
			 bigint aa8 = 8*asq;
			 ok = ::divides(r-temp,aa8,d,rem);
#ifdef DEBUG_AH
			 cout<<"; aa8,d,rem = "<<aa8<<","<<d<<","<<rem<<flush;
#endif
			 if(!ok) {ah_dfail++; continue;}
#else
			 bigfloat xd = (xr-temp)/(8*xa*xa);
			 bigfloat xxd = abs(xd-floor(xd+0.5));
			 if (xxd>abceps) {ah_dfail++; continue;}
#endif
#ifdef USE_BIGINTS
			 ee=ii_cc+sb3*d;
#ifdef DEBUG_AH
			 cout<<"\n ii,b3,d,ee = "<<ii<<","<<sb3<<","
			   <<d<<","<<ee;
#endif
			 bigint a12=12*biga;
			 ok = ::divides(ee,a12,e,rem);
#ifdef DEBUG_AH
			 cout<<"\n ee,a12,e,rem = "<<ee<<","<<a12<<","<<
			   e<<","<<rem<<flush;
#endif
			 if(!ok) {ah_efail++; continue;}
#else
			 bigfloat xe = (xii-xcc+sb3*xd)/xa;
			 bigfloat xxe = abs(xe-floor(xe+0.5));
			 if (xxe>abceps) {ah_efail++; continue;}
			 d = Iround(xd);
			 e = Iround(xe);
			 
#endif
#ifdef DEBUG_AH
			 cout << ":\n [" << a<<","<<sb << "," << c 
			   << "," << d << "," << e << "]"<<endl; 
#endif
// 
// Now test divisibility conditions in extra2 case:
// (already know 4|b   since 16|h, and not(4|a))
//
			 int skip=0;
			 if(extra2)
			   {
			     skip = !(ndiv(efactor,e)&&ndiv(efactor,a+sb+c+d+e));
			   }
			 if(skip) {ah_extra2fail++; continue;}
			 ah_pass++;
// Now we have a quartic
// Check the invariants are right (for debugging only):
			 bigb=sb;
			 bigint iiabcde = II(biga,bigb,bigc,d,e);

			 if ( ii != iiabcde )
			   {
			     cout<<"Error: constructed quartic ";
			     
			     cout << "[" << a<<","<<sb << "," << c << "," << d << "," << e << "]"; 
			     cout << " has wrong I-invariant "<<iiabcde<<", not "<<ii<<endl;
			     continue;
			   }
			 bigint jjabcde = JJ(biga,bigb,bigc,d,e);
			 if (jj != jjabcde)
			   {
			     cout<<"Error: constructed quartic ";
			     cout << "[" << a<<","<<sb << "," << c << "," << d << "," << e << "]"; 
			     cout << " has wrong J-invariant "<<jjabcde<<", not "<<jj<<endl;
			     continue;
			   }
// And finally when disc>0 check that the type is correct
			 if(posdisc)
			   {
			     bigint habcde = biga8*bigc-bigbb3;
			     bigint qabcde = habcde*habcde-16*biga*biga*ii;
			                // =3*Q
			     if(type==1)
			       {
			       if((habcde<0)&&(qabcde>0))
				 {
				   cout<<"Error: constructed quartic ";
				   cout<<"[" << a<<","<<sb << "," << c << "," << d << "," << e << "]"; 
				   cout<<" has type 2, not 1!\n";
				   cout<<"Please report"<<endl;
				   continue;				 
				 }
			     }
			     else // type==2
			       if((habcde>=0)||(qabcde<=0))
				 {
				   cout<<"Error: constructed quartic ";
				   cout<<"[" << a<<","<<sb << "," << c << "," << d << "," << e << "]"; 
				   cout<<" has type 1, not 2!\n";
				   cout<<"Please report"<<endl;
				   continue;				 
				 }
			   }
// Now construct the roots of the quartic and add it to the list.

			 bigfloat three=to_bigfloat(3);
			 switch(type) {
			 case 1: // no real roots
			   r1 = sqrt((xa4*phi1-xh)/three);
			   r2 = sqrt(-(xa4*phi2-xh)/three);
			   r3 = sqrt(-(xa4*phi3-xh)/three);
			   croots[0]=bigcomplex( r1-sb,  r2-r3) * oneover4a;
			   croots[1]=bigcomplex( r1-sb, -r2+r3) * oneover4a;
			   croots[2]=bigcomplex(-r1-sb,  r2+r3) * oneover4a;
			   croots[3]=bigcomplex(-r1-sb, -r2-r3) * oneover4a;
			   break;
			 case 2: // all real roots
			   r1 = sqrt((xa4*phi1-xh)/three);
			   r2 = sqrt((xa4*phi2-xh)/three);
			   r3 = sqrt((xa4*phi3-xh)/three);
			   croots[0]=bigcomplex( r1+r2-r3-sb) * oneover4a;
			   croots[1]=bigcomplex( r1-r2+r3-sb) * oneover4a;
			   croots[2]=bigcomplex(-r1+r2+r3-sb) * oneover4a;
			   croots[3]=bigcomplex(-r1-r2-r3-sb) * oneover4a;
			   break;
			 case 3: // roots 2,3 are real
			   c1 = sqrt((xa4*cphi[0]-xh)/three);
			   r3 = sqrt((xa4*phi-xh)/three);
			   if(xr<0) r3=-r3;
			   croots[0]=bigcomplex( r3-sb,  2*imag(c1)) * oneover4a;
			   croots[1]=bigcomplex( r3-sb, -2*imag(c1)) * oneover4a;
			   croots[2]=bigcomplex(-r3-sb  +2*real(c1)) * oneover4a;
			   croots[3]=bigcomplex(-r3-sb  -2*real(c1)) * oneover4a;
			   break;
			 }
			 addquartic(biga,bigb,bigc,d,e);
			 if((type==1)&&(have_eggpoint))
			   {
			     if (verbose) 
			       {
				 cout << "Exiting search for Type 1 quartics after ";
				 cout << "finding one which is globally soluble.\n";
			       }
			     return;
			   }
#if LARGE_Q>3			 
			 if((extra2)&&(have_large_quartics))
			   {
			     if (verbose) 
			       {
				 cout << "Exiting search for large quartics after ";
				 cout << "finding enough globally soluble ones.\n";
			       }
			     return;
			   }
#endif // LARGE_Q>3
		       }  // end of b-sign-loop
		    } // end of c loop
		} // end of b conditions
	    } // end of b loop
	} // end of a conditions
    } // end of main a loop
    } // end of loop on sign of a
   if (verbose) cout << "Finished looking for Type " << t << " quartics.\n";
}  // end of gettype()


rank1::rank1(Curvedata* ec, int verb, int sel, long lim1, long lim2,long n_aux)
  : rank12(ec,verb,sel,lim1,lim2,n_aux,1)
{
  traceequiv=0;
  success=1; // the default!
  if(num_aux==-1) num_aux=DEFAULT_NAUX;
  if(verbose>1) 
    {
      cout << "Using (a,b,c) search with (a,h) sieve and algebraic method\n";
#ifdef USE_BIGINTS
      cout << "(with bigints to solve the syzygy)\n";
#else
      cout << "(with bigfloats to solve the syzygy)\n";
#endif
    }
  
  qlista = new quartic[maxnquartics];
  qlistb = new quartic[maxnquartics];
  qlistbflag = new int[maxnquartics];
  croots = new bigcomplex[4];

  the_curve->getci(c4,c6);
  d1728 = c4*c4*c4-c6*c6;
  if (is_zero(d1728)) {cout<<"Curve is singular\n"; success=0; return;}

  // Set up the transformation [u,r,s,t] from the minimal model to the model
  // [0,0,0,-27*c4,-54*c6]; from these we will later obtain (by simple scaling) 
  // the transformations to the IJ-curve for various I,J
  bigint a1,a2,a3,a4,a6,b2=getb2(*the_curve); 
  the_curve->getai(a1,a2,a3,a4,a6);
  tr_u=6; tr_r=3*b2; tr_s=3*a1; tr_t=108*a3;

  vector<bigint> ir = Introotscubic( BIGINT(0),-27*c4,-54*c6);  
  n0=ir.size()+1;

  long e0,e1,e2;
  if(!intlog2(n0,e0,0))
    {
      success=0; 
      cout<<"!!! Fatal error: n0=#E[2]="<<n0<<" is not a power of 2\n";
      cout<<"Please inform author by email!\n";
      return;
    }
  
  posdisc = is_positive(d1728);
  npoints1=npoints2=0;
  n1=n2=1;

#ifdef TEST_EQCODE
  cout << "Setting eqplist, length = " << NEQPLIST << "\n";
#endif
  eqplist.reserve(NEQPLIST);
  for(primevar pp; eqplist.size()<NEQPLIST; pp++)
    {
      long p = pp;
      if(ndiv(p,d1728)) eqplist.push_back(p);
    }
#ifdef TEST_EQCODE
  cout << "eqplist = " << eqplist << "\n";
#endif

  getquartics();

  //  cout<<"LOCAL INDEX "<<twoadic_index<<" GLOBAL INDEX "<<global_index<<" BSD "<<bsd_npairs<<endl;


  // Compute rank/selmer rank from B=im(eps) first:
  // Must count the quartics in qlistb which are still needed

  long jp, n3 = 0;
  for(jp=0; jp<nquarticsb; jp++) if(qlistbflag[jp]) n3++;
  if(verbose>1)
    {
      cout<<"After getquartics(): \n";
      cout<<"n1 = "<<n1<<endl;
      cout<<"n2 = "<<n2<<endl;
      cout<<"n3 = "<<n3<<endl;
      cout<<"B-rank = "<<rank_B<<endl;
    }
  delete [] qlista;
  delete [] qlistb;
  delete [] croots;
  delete [] qlistbflag;
  
  if(n3>0){
    if(n2>1){
      if(n3%n2==0) {
	n3/=n2;
      }
      else {
	cout<<"\n!!! n3 = "<<n3<<" not a multiple of n2 = "<<n2<<endl;
      }
    }
  }
  n3++;
  long keep_n3 = n3, e3;
  if (!intlog2(n3,e3,1))
    {
      cout<<"\n!!! n3 = "<<keep_n3<<" not a power of 2, rounding up to "<<n3<<"\n";
    }
  long selmer_rank_B = rank_B + e3;

  if(verbose)
    {
      if(!selmer_only) 
	cout << "Mordell rank contribution from B=im(eps) = " << rank_B << endl;
      cout << "Selmer  rank contribution from B=im(eps) = " << selmer_rank_B << endl;
      if(!selmer_only) 
	cout << "Sha     rank contribution from B=im(eps) = " << e3 << endl;
    }
  
  long keep_n1=n1, keep_n2=n2;
  if (!intlog2(n1,e1,1))
    {
      cout<<"\n!!! n1 = "<<keep_n1<<" not a power of 2, rounding up to "<<n1<<"\n";
      cout<<"(Probably due to too small a bound on quartic point search\n";
      cout<<" leading to rational points not being found)\n";
      cout<<" The points listed will be incomplete, but may still\n";
      cout<<" generate a subgroup of finite index.\n";
    }
  if (!intlog2(n2,e2,0))
    {
      success=0; 
      cout<<"\n\n!!! Fatal error: n2 = "<<keep_n2<<" not a power of 2\n";
      cout<<"Please inform author by email!\n";
    }
  
  long rank_A = e1-e0,  selmer_rank_A = e2-e0;
    
  if(verbose)
    {
      if(!selmer_only) 
	cout << "Mordell rank contribution from A=ker(eps) = " << rank_A << endl;
      cout << "Selmer  rank contribution from A=ker(eps) = " << selmer_rank_A << endl;
      if(!selmer_only) 
	cout << "Sha     rank contribution from A=ker(eps) = " << (e2-e1) << endl;
    }

  n1 <<= rank_B;
  n2 <<= selmer_rank_B;
  rank = rank_A + rank_B;
  selmer_rank = selmer_rank_A + selmer_rank_B;
  rank_bound = selmer_rank;
  sha_rank = selmer_rank - rank;
  sha2 = (n2/n1);
  certain = (sha_rank==0)||(selmer_only);
  int strange = odd(sha_rank);

  if(verbose&&strange&&(!selmer_only))
    {
      cout<<"\nWarning: Selmer rank = "<<selmer_rank<<" and program finds \n";
      cout<<"lower bound for rank = "<<rank<<" which differs by an odd\n";
      cout<<"integer from the Selmer rank.   Hence the rank must be 1 more\n";
      cout<<"than reported here.  Try rerunning with a higher bound for\n";
      cout<<"quartic point search.\n";
    }
  
  if (verbose&&(!certain))
    {
      cout << "\nSummary of results (all should be powers of 2):\n\n";
      cout<<"n0 = #E(Q)[2]    = "<<n0<<"\n";
      cout<<"n1 = #E(Q)/2E(Q) ";
      if(!certain) cout<<">";
      cout<<"= "<<n1<<"\n";
      cout<<"n2 = #S^(2)(E/Q) = "<<n2<<"\n";
      cout<<"#III(E/Q)[2]     ";
      if(!certain) cout << "<"; 
      cout<<"= "<<sha2<<"\n\n";

      if(certain) 
	cout << "rank " << "= " << rank;
      else 
	cout << rank << " <= rank <= selmer-rank = " << selmer_rank;
      cout << endl << endl; ;

    }  

  if(verbose&&selmer_only)
    cout << "Selmer rank = " << selmer_rank << endl;

  if(rank<6) sortpoints();

}  // end of rank1 constructor

void rank1::sortpoints()  // reorder points into increasing height order
{
  long i,j;
  for(i=0; i<npoints1; i++)
    for(j=i+1; j<npoints1; j++)
      if(height(pointlist1[j])<height(pointlist1[i]))
	{
	  Point temp = pointlist1[i]; 
	  pointlist1[i]=pointlist1[j]; 
	  pointlist1[j]=temp;
	}
  for(i=0; i<npoints2; i++)
    for(j=i+1; j<npoints2; j++)
      if(height(pointlist2[j])<height(pointlist2[i]))
	{
	  Point temp = pointlist2[i]; 
	  pointlist2[i]=pointlist2[j]; 
	  pointlist2[j]=temp;
	}
}

void showpoint(Point P)
{
  bigfloat h = height(P);
  cout << P << ", height = " << h;
  if(!P.isvalid()) {cout << " --warning: NOT on curve!\n"; abort();}
  cout << "\n";
}

void showpoint(Point P, Curvedata* CD, const bigint& u, const bigint& r, 
	                               const bigint& s, const bigint& t)
{
  showpoint(transform(P,CD,u,r,s,t,1));
}

void rank1::listpoints(Curvedata* CD_orig, const bigint& u, const bigint& r,
		       const bigint& s, const bigint& t)
{
  int explanation_needed = (npoints1>0)&&(npoints2>0);
  if(explanation_needed)
    {
      cout<<"p-adic filtration expresses E(Q)/2E(Q) as a direct sum A+B\n";
      cout<<"where A = E(Q)\\cap\\sum 2E(Q_p) for certain primes p.\n";
      cout<<"We list all nonzero points of A, and generators of B\n";
    }

  if(npoints1>0)
    {
      if(explanation_needed)
	{
	  cout << "Points in A:\n";
	}
      else
	{
	  cout << "Points covering E(Q)/2E(Q):\n";
	}
      for (long i=0; i<npoints1; i++)
	{
	  Point p = pointlist1[i];
	  //	  cout<<p<<" on "<<(Curve)((p.getcurve()))<<endl;
	  cout << "Point "; showpoint(p,CD_orig,u,r,s,t);
	}
    }

  if(npoints2>0)
    {
      if(explanation_needed)
	{
	  cout << "Points generating B:\n";
	}
      else
	{
	  cout << "Points generating E(Q)/2E(Q):\n";
	}
      for (long i=0; i<npoints2; i++)
	{
	  Point p = pointlist2[i];
	  //	  cout<<p<<" on "<<(Curve)((p.getcurve()))<<endl;
	  cout << "Point "; showpoint(p,CD_orig,u,r,s,t);
	}
    }
}

void rank1::listpoints()
{
  int explanation_needed = (npoints1>0)&&(npoints2>0);
  if(explanation_needed)
    {
      cout<<"p-adic filtration expresses E(Q)/2E(Q) as a direct sum A+B\n";
      cout<<"where A = E(Q)\\cap\\sum 2E(Q_p) for certain primes p.\n";
      cout<<"We list all nonzero points of A, and generators of B\n";
    }

  if(npoints1>0)
    {
      if(explanation_needed)
	{
	  cout << "Points in A:\n";
	}
      else
	{
	  cout << "Points covering E(Q)/2E(Q):\n";
	}
      for (long i=0; i<npoints1; i++)
	{
	  Point p = pointlist1[i];
	  //	  cout<<p<<" on "<<(Curve)((p.getcurve()))<<endl;
	  cout << "Point "; showpoint(p);
	}
    }

  if(npoints2>0)
    {
      if(explanation_needed)
	{
	  cout << "Points generating B:\n";
	}
      else
	{
	  cout << "Points generating E(Q)/2E(Q):\n";
	}
      for (long i=0; i<npoints2; i++)
	{
	  Point p = pointlist2[i];
	  //	  cout<<p<<" on "<<(Curve)((p.getcurve()))<<endl;
	  cout << "Point "; showpoint(p);
	}
    }
}

vector<Point> rank1::getpoints()
// We construct a set of coset reps for 2E(Q) in E(Q) given
// reps for the subgroup A in pointlist1 and 
// gens for the complementary subgroup B in pointlist2
{
  long np = (1+npoints1) << npoints2;
  vector<Point> ans;
  long j, k, ip=1+npoints1;
  ans.push_back(Point(the_curve,BIGINT(0),BIGINT(1),BIGINT(0)));
  ans.insert(ans.end(),pointlist1.begin(),pointlist1.end());
  ans.resize(np);
  for(j=0; j<npoints2; j++, ip*=2) 
    for(k=0; k<ip; k++)
      {
	ans[ip+k]=ans[k]+pointlist2[j];
      }
  return ans;
}

vector<Point> rank1::getgens() const
// Returns a set of generators for E(Q) mod 2E(Q)
// (but not necessarily independent)
{
  vector<Point> ans; ans.reserve(pointlist1.size()+pointlist2.size());
  copy(pointlist1.begin(),pointlist1.end(),back_inserter(ans));
  copy(pointlist2.begin(),pointlist2.end(),back_inserter(ans));
  return ans;
}

void rank1::aux_init()  // define  auxiliary moduli and squares
{
  long i, j, a;

  auxs = new long[num_aux];
  phimod = new long*[num_aux];
  aux_flags = new int[num_aux];
  aux_types = new int[num_aux];
  squares = new int*[num_aux];
  flags = new int**[num_aux];
  flaga = new int*[num_aux];
  amod = new long[num_aux];
  hmod = new long[num_aux];
  hstepmod = new long[num_aux];
  astepmod = new long[num_aux];
  hscalemod = new long[num_aux];

  auxs[0]=9;  // treated specially
  aux_flags[0]=1;  aux_types[0]=0;
  for(i=0; i<num_aux; i++) phimod[i] = new long[3];
  i=1;

  // the rest of the auxs must be chosen carefully:  if possible they should
  // be good odd primes p, such that the resolvent cubic is not irreducible mod p.
  // If it has a unique root phi mod p, E(Qp)/2E(Qp) has order 2 and the coset in 
  // which the image of a quartic lies depends on whether it has 0 or 2 roots mod p.
  // If it has 3 roots  mod p, E(Qp)/2E(Qp) has order 4 and the coset in 
  // which the image of a quartic lies depends on whether it has 0 or 4 roots mod p; 
  // if 0, then a further condition determines which non-trivial coset it belongs to

  primevar pr; pr++; pr++;  // skip past 2 and 3

  for(;pr.ok()&&i<num_aux; pr++)
    {
      long p = pr;
      if(div(p,disc)) continue;
      long minus3imodp = mod(-3*ii,p);
      long jmodp = mod(jj,p);
      long nr = nrootscubic(0,minus3imodp,jmodp,p,phimod[i]);
      if(nr>0)
	{
	  auxs[i]=p;
	  aux_flags[i] = 1;
	  aux_types[i] = 1;
	  if(nr>1)  aux_types[i] = 2;
	  i++;
	}
    }

  // report on which primes will be used:

  if((verbose>1)&&(num_aux>0)) 
    {
      cout<<"(a,h) sieving using " <<num_aux<<" moduli: \n";
      cout<<"p:\t";
      for(j=0; j<num_aux; j++) cout<<auxs[j]<<"\t";
      cout<<"\n";
      cout<<"k_p:\t\t";
      for(j=1; j<num_aux; j++) cout<<aux_types[j]<<"\t";
      cout<<"\n";
      cout<<"phi1:\t\t";
      for(j=1; j<num_aux; j++) cout<<phimod[j][0]<<"\t";
      cout<<"\n";
      cout<<"phi2:\t\t";
      for(j=1; j<num_aux; j++) 
	if(aux_types[j]==1) cout<<"*\t";
	else 
	  cout<<phimod[j][1]<<"\t";
      cout<<"\n";
      cout<<"phi3:\t\t";
      for(j=1; j<num_aux; j++) 
	if(aux_types[j]==1) cout<<"*\t";
	else 
	  cout<<phimod[j][2]<<"\t";
      cout<<"\n";
    }

  // initialize flag arrays for squares:

  for (i = 0; i < num_aux; i++)
    {
      long aux = auxs[i];
      long half_aux = ((aux + 1) / 2);
      squares[i] = new int[aux];
      for (j = 0; j < aux; j++)      squares[i][j]=0;
      for (j = 0; j < half_aux; j++) squares[i][posmod( j*j, aux )]=1;

      flags[i] = new int*[aux];
      for(a=0; a<aux; a++)
	flags[i][a] = new int[aux];
    } // end of aux loop

  // initialize scaling factors for use with large I,J pair:
  // NB we use the same sieve for both; (a,h) passes for the larger I,J 
  // iff (a,h/4) passes for the standard I,J.

  for(i=0; i<num_aux; i++)
    {
      hscalemod[i] = invmod(4,auxs[i]);
    }

  if((verbose>1)&&(num_aux>0)) cout<<"finished aux_init()"<<endl;
}

//#define COUNT_CODES

void rank1::flag_init() // set up flag array
{
  if((verbose>1)&&(num_aux>0))  cout<<"starting flag_init()"<<endl;

  int thisflag;
  int ***flagsi=flags;
  int **squaresi=squares;
  long * a4phi= new long[3];
  long * eps  = new long[3];

#ifdef COUNT_CODES
  long * code_count = new long[5]; long icc;
#endif

  for(long i=0; i<num_aux; i++, squaresi++, flagsi++) 
    {
#ifdef COUNT_CODES
      for(icc=0; icc<5; icc++) code_count[icc]=0;
#endif
      
      int case1 = (aux_types[i]==1);  // phi cubic has 1 root  mod p
      int case2 = !case1;             // phi cubic has 3 roots mod p
      long a, h;
      long aux = auxs[i];
      long aux2 = (i==0 ? 27 : aux);
      long I = mod(ii,aux2), J = mod(jj,aux2);
      long I16 = (16*I)%aux2, I48 = (3*I16)%aux2, J64 = (64*J)%aux2;

      int **flagsia = *flagsi;

      for(a=0; a<aux; a++, flagsia++)
	{
	  long a2=(a*a)%aux2;  long a2I48 = (a2*I48)%aux2;  long a2I16 = (a2*I16)%aux2;
	  long a3=(a*a2)%aux2; long a3J64 = (a3*J64)%aux2;
	  if(i>0)
	    {
	      a4phi[0] = (4*a*phimod[i][0])%aux2;
	      if(case2)
		{
		  a4phi[1] = (4*a*phimod[i][1])%aux2;
		  a4phi[2] = (4*a*phimod[i][2])%aux2;
		}
	    }

	  int *flagsiah = *flagsia;

	  for(h=0; h<aux; h++, flagsiah++)
	    {
	      long h2 = (h*h)%aux2;
	      long disc = (((h*((h2-a2I48)%aux2))%aux2)+a3J64)%aux2;
	      if(i==0) // special mod-9 condition to force 27-divisibility
		{
		  *flagsiah = (disc==0);
		}
	      else
		{
		  disc = posmod(-3*disc,aux2);
		  thisflag = (*squaresi)[disc];
		  if(thisflag) 
		    {
	      // look further to see how many roots mod p an (a,h) quartic
	      // could have.  
		      if(case1)
			{
		  //By choice of auxs there must be 0 or 2 roots, and
		  // the flag is set to 15 if there are 2 roots, else 5
			  thisflag=5;
			  if(disc==0) //must look at Q-seminvariant
			    {
			      long q3 = posmod(3*(h2-a2I16),aux2);
			      if((*squaresi)[q3]) thisflag=15;
			    }
			  else
			    {
			      long z3 = posmod(3*(a4phi[0]-h),aux2);
			      if((*squaresi)[z3]) thisflag=15;
			    }
			}
		      else // case 2
			{
		  // By choice of auxs there must be 0 or 4 roots;
		  // the flag is set to 15 if there are 4 roots, else to
		  // one of 5, 3, 1 depending on which element of
		  // E(Qp)/2E(Qp) the quartic maps to:
			  long iz, z;
			  for(iz=0; iz<3; iz++)
			    {
			      z = posmod(3*(a4phi[iz]-h),aux2);
			      eps[iz] = 2*((*squaresi)[z])-(z==0)-1;
			      // = -1, 0, +1
			    }
// At most one eps is 0 (since E has good reduction at p).
// In this case replace 0 by the value which makes the product +1

			  if(eps[0]==0) {eps[0]=eps[1]*eps[2];}
			  else {
			    if(eps[1]==0) {eps[1]=eps[0]*eps[2];}
			    else {
			      if(eps[2]==0) {eps[2]=eps[0]*eps[1];}
			    }
			  }

//Now each eps[i] = +1 or -1 and the product is +1

//(+,+,+) maps to flag 15=8+4+2+1     (+,-,-) maps to flag 5=  4 + 1
//(-,+,-) maps to flag  3=    2+1     (-,-,+) maps to flag 1=      1
			  
			  if(eps[0]==1) thisflag=(eps[1]==1? 15: 5);
			  else          thisflag=(eps[1]==1?  3: 1);
			}
		    }
#ifdef COUNT_CODES
		  if(thisflag==0) code_count[0]++;
		  if(thisflag==15) code_count[1]++;
		  if(thisflag== 5) code_count[2]++;
		  if(thisflag== 3) code_count[3]++;
		  if(thisflag== 1) code_count[4]++;
#endif
		  *flagsiah = thisflag;
		}
	    }
	}
#ifdef COUNT_CODES
      if(i>0)
	{
	  cout << "Code count for p = " << aux << ":\n";
	  cout << 0 << "\t"<< 15 << "\t"<< 5 << "\t"<< 3 << "\t"<< 1 << "\n";
	  for(icc=0; icc<5; icc++) cout<<code_count[icc]<<"\t";
	  cout<<endl;
	  double ratio = ((double)(code_count[0]))/(aux*aux);
	  cout<<"Percentage of (a,H) pairs failing sieve = "<<100*ratio<<endl;
	}
#endif
    }
  delete [] a4phi; delete [] eps;
#ifdef COUNT_CODES
  delete [] code_count;
#endif

  if((verbose>1)&&(num_aux>0)) 
    cout<<"finished flag_init()"<<endl;

}

void rank1::clear_sieve()  // free memory related to sieve;
{
  for(long i=0; i<num_aux; i++) 
    {
      long aux = auxs[i];
      delete[] squares[i];
      for(long a=0; a<aux; a++)
	{
	  delete[] flags[i][a];
	}
      delete[] flags[i];
      delete[] phimod[i];
    }
  delete[] auxs; 
  delete[] phimod; 
  delete[] squares;
  delete[] aux_flags;  delete[] aux_types;
  delete[] flags;  delete[] flaga;
  delete[] amod; delete[] hstepmod; delete[] hscalemod;
  delete[] hmod; delete[] astepmod;
}
