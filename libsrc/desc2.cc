// desc2.cc:  implementation of second descent (via 2-isogeny) procedure
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
 
#include "mquartic.h"
#include "unimod.h"
#include "quadratic.h"
#include "conic.h"
#include "mlocsol.h"
#include "mglobsol.h"
#include "transform.h"
#include "minim.h"
#include "reduce.h"
#include "sqfdiv.h"
#include "desc2.h"
#include "timer.h"  // for timing of conic solving
#include "hilbert.h"

#ifndef QSIEVE_OPT
#define QSIEVE_OPT 0 // uses Stoll's sieve
#endif

int process_d3(const quadratic& q0, const bigint& d3, 
	       const vector<bigint>& plist, const vector<bigint>& factorbase, double hlim,
	       const quadratic& q1, const quadratic& q3,
	       bigint& x, bigint& y, bigint& z, int verb, int selmer_only=0);
// Processes an individual d3
// Returns -1 if not els
//          0 if els but no point found
//         +1 if els and point found
//
// if selmer_only==1, returns 0 or -1, with no point searching
//
// We only pass d3 to process_d3 if the Hilbert symbols tests show
// that both q1=d3 and q3=d3 are soluble (but this does not yet imply
// that they are simultaneously soluble).

int desc2(const bigint& c, const bigint& d1, const bigint& d2,
	  const vector<bigint>& plist, const vector<bigint>& supp, const vector<bigint>& bgens,
	  long mask,  double hlim,
	  bigint& x, bigint& y, bigint& z, int verb, int selmer_only, int alldesc)
// Works on homogeneous space (d1,0,c,0,d2) (assumed ELS)
// Returns 
//   -1 if it certainly has no points (if no ELS descendents)
//   +1 if it has a point (coordinates returned in x, y, z)
//    0 if undecided (ELS descendents exist but no rational points were found)
//
//  plist is a complete list of bad primes (including 2)
//  supp is a list of primes (and -1) dividing d'
//  bgens is a list of generators of the "opposite" Selmer group B', a subgroup
//        of group A' = the span of supp
//  mask enables one to loop over A' modulo B'
//
// if alldesc==1 it does not stop when it finds one descendent with a point on it,
// but goes on to look at all the others.
//
// if selmer_only==1 it only checks whether els descendents exist,
// returning -1 or 0, but does no global point searching.
//
{
  int xverb = (verb>1), res;
  if(verb) 
    {
      cout<<"Using desc2("<<d1<<","<<c<<","<<d2<<")\n";
      cout<<"supp="<<supp<<"; mask="<<mask<<"; bgens="<<bgens<<endl;
    }
  bigint x0,z0,d0, p, result;
  bigint d = d1*d2, cdash = -2*c, ddash = sqr(c)-4*d;
  static const bigint one = BIGINT(1);
  static const bigint two = BIGINT(2);
  int add2tosupp = (val(2,ddash)==4);  
  if(add2tosupp) add2tosupp = (find(supp.begin(),supp.end(),two)==supp.end());  
// For we are on E' and the original d was odd

  // Step 1: solve the conic d1*x0^2 + c*x0*z0 + d2*z0^2 = y0^2 
  // Step 2: with parametrization

  int ddash_is_square = isqrt(ddash,d0);

  quadratic q0(d1,c,d2), q1, q2, q3;
  
  res = solve_conic_param(q0,one,plist,q1,q2,q3);

  if(!res) {cout<<"solve_conic failed -- should not happen!\n"; return -1;}
  x0=q1[0]; z0=q3[0];
  vector<bigint> factorbase = plist;
  int factorbase_enlarged=0;

  // We only need to factorize x0, z0 if some d3 passes the Hilbert
  // symbol test, since we use their factorizations to solve q1=d3 but
  // NOT to test whether q1=d3 and q3=d3 are soluble.  We'll carry out
  // the following only when necessary:

  //  factorbase=vector_union(factorbase,pdivs(x0));
  //  factorbase=vector_union(factorbase,pdivs(z0));

  // The quadratics q1, q3 have discriminants 4*d2, 4*d1, and resultant d'.

  if(xverb) 
    {
      cout<<"q1-coeffs: "<<q1<<"\n";
      cout<<"q3-coeffs: "<<q3<<"\n";
      result = resultant(q1,q3);
      cout<<"resultant(q1,q3)   = "<<result;
      if(result==ddash) cout<<" = d'\n";
      else cout<<" = "<<(result/ddash)<<" * d'\n";
      cout<<"support of d': "<<supp<<endl;
    }

  // Step 3
  long nsupp = supp.size();
  long nd3=1<<nsupp, id3;

  // Step 4: For all square-free divisors of the resultant, attempt to
  //         form a descendent quartic

  bigint d3, keepd3, s0, t0, u0, s1, t1, u1;
  int looking=1, found=0, hres;

  for(id3=0; (id3<nd3)&&(looking||alldesc); id3++) // must start at 0
    {
      d3=makenum(supp,id3);
      if(id3&mask) continue;

// We first test whether q1=d3 and q3=d3 are soluble using Hilbert symbols:
      hres = global_hilbert(q1,d3,factorbase,p) || 
	     global_hilbert(q3,d3,factorbase,p);
      if(!hres)
	{
	  if(verb) cout<<"d3= "<<d3<<" passes Hilbert symbol tests\n";
	  if(xverb) cout<<"About to factorize "<<x0<<" and "<<z0<<endl;
	  if(!factorbase_enlarged)
	    {
   //	      cout<<"Before enlarging, factorbase = "<<factorbase<<endl;
	      factorbase = vector_union(factorbase,pdivs(x0));
	      factorbase = vector_union(factorbase,pdivs(z0));
	      factorbase_enlarged=1;
   //	      cout<<"After enlarging, factorbase = "<<factorbase<<endl;
	    }
	  res = process_d3(q0,d3,plist,factorbase,hlim,q1,q3,x,y,z,verb,selmer_only);
	  if(xverb) cout<<"process_d3("<<d3<<") returns "<<res<<endl;
	}
      else 
	{
	  //	  if(xverb) 
	  //	    cout<<"d3= "<<d3<<" fails Hilbert symbol tests (p="<<p<<")\n";
	  res=-1;
	}
      if(res!=-1)  // descendent is els, lift to S^(2) exists
	{
	  looking=0; 
	  keepd3=d3;
	}
      if(res==1)   // descendent is gls, lift to E/2E exists
	{
	  found=1;
	}

      if((looking||alldesc)&&add2tosupp)
	{
	  d3=2*d3;
// We first test whether q1=d3 and q3=d3 are soluble using Hilbert symbols:
  	  hres = global_hilbert(q1,d3,factorbase,p) || 
  	         global_hilbert(q3,d3,factorbase,p);
	  if(!hres)
	    {
	      if(verb) cout<<"d3= "<<d3<<" passes Hilbert symbol tests\n";
	      if(!factorbase_enlarged)
		{
		  if(xverb)
		    cout<<"About to factorize "<<x0<<" and "<<z0<<endl;
  //	          cout<<"Before enlarging, factorbase = "<<factorbase<<endl;
		  factorbase = vector_union(factorbase,pdivs(x0));
		  factorbase = vector_union(factorbase,pdivs(z0));
		  factorbase_enlarged=1;
  //	          cout<<"After enlarging, factorbase = "<<factorbase<<endl;
		}
	      res = process_d3(q0,d3,plist,factorbase,hlim,q1,q3,x,y,z,verb,selmer_only);
	      if(xverb) cout<<"process_d3("<<d3<<") returns "<<res<<endl;
	    }
	  else 
	    {
	      //	      if(xverb) 
	      //		cout<<"d3= "<<d3<<" fails Hilbert symbol tests (p="<<p<<")\n";
	      res=-1;
	    }
	  if(res!=-1)  // descendent is els, lift to S^(2) exists
	    {
	      looking=0; 
	      keepd3=d3;
	    }
	  if(res==1)   // descendent is gls, lift to E/2E exists 
	    {
	      found=1;
	    }
	}
    }
  if(found) 
    {
      if(verb) 
	cout<<"Found a descendent with a rational point so terminating second descent step.\n";
      return +1;   // Found an els descendent with a point on it
    }
  if(looking) return -1; // No els descendents exist
  
  if(selmer_only) return 0;

  //
  // Go on searching on other descendents...
  //
  if(verb) cout<<"Found an els descendent but no rational point on it...\n";

  nd3=1<<bgens.size();
  long d3step=(2-ddash_is_square);
  long ndesc=nd3/d3step;
  if(nd3==d3step)
    {
      if(verb) cout<<"No further descendents.\n";
      return 0;
    }
  if(verb) 
    cout<<"We now construct and search the "<<(ndesc-1)<<" other descendents...\n";
  for(id3=d3step; (id3<nd3); id3+=d3step)
    {
      d3 = sqfmul(keepd3,makenum(bgens,id3));
      int res = process_d3(q0,d3,plist,factorbase,hlim,q1,q3,x,y,z,verb);
      if(xverb) cout<<"process_d3("<<d3<<") returns "<<res<<endl;
      if(res==1)   // descendent is gls, lift to E/2E exists
	{
	  found=1;
	  if(verb) 
	    cout<<"Found a descendent with a rational point so terminating second descent step.\n";
	  if(!alldesc) return +1;
	}
    }
  if(verb&&(!found)) cout<<"No rational point found on any of the other descendents...\n";
  return found;
  //
  //
  //
}  // end of desc2()


int process_d3(const quadratic& q0, const bigint& d3, 
	       const vector<bigint>& plist, const vector<bigint>& factorbase, double hlim,
	       const quadratic& q1, const quadratic& q3,
	       bigint& x, bigint& y, bigint& z, int verb, int selmer_only)
// Processes an individual d3
// Returns -1 if not els
//          0 if els but no point found
//         +1 if els and point found
//
// if selmer_only==1, returns 0 or -1, with no point searching
//
// We only pass d3 to process_d3 if the Hilbert symbols tests show
// that both q1=d3 and q3=d3 are soluble (but this does not yet imply
// that they are simultaneously soluble).
{
  int xverb=verb>1;
  if(verb) cout<<"Processing d3 = "<<d3<<": \t";
  //  if(xverb) cout<<"\nplist = "<<plist<<"\n";
  bigint s1,u1,t1,ga,gb,gc,gd,ge,cont;
  bigint p,ggI,ggJ,ggD;
  static const bigint one = BIGINT(1), zero = BIGINT(0);
  vector<bigint> ggbadp, ggextrap;
  quadratic Q1, Q2, Q3;
  int resd3 = solve_conic_param(q1,d3,factorbase,Q1,Q2,Q3);

  if(!resd3)
    {
      cout<<"Problem solving q1=d3!\n"<<endl;
      cout<<"q1 = "<<q1<<endl;
      cout<<"d3 = "<<d3<<endl;
      cout<<"factorbase= "<<factorbase<<endl;
      if(verb) cout<<"q1=d3 not soluble.\n";
      return -1;
    } 
  if(verb) cout<<"q1=d3 is soluble.\t";

  resd3 = solve_conic(q3,d3,factorbase,s1,u1,t1);

  if(!resd3)
    {
      cout<<"Problem solving q3=d3!\n"<<endl;
      cout<<"q3 = "<<q3<<endl;
      cout<<"d3 = "<<d3<<endl;
      cout<<"factorbase= "<<factorbase<<endl;
      if(verb) cout<<"q3=d3 not soluble.\n";
      return -1;
    }
  if(verb) cout<<"q1=d3 and q3=d3 both soluble, forming quartic.\n";

  // g(x,y)=q3(Q1(x,y),Q3(x,y)):
	      
  ga = q3(Q1[0],Q3[0]);
	  
  gb = q3[0]*(2*Q1[0]*Q1[1]) + q3[1]*(Q1[0]*Q3[1]+Q1[1]*Q3[0]) + q3[2]*(2*Q3[0]*Q3[1]);
	  
  gc = q3[0]*(sqr(Q1[1]) + 2*Q1[0]*Q1[2]) 
    + q3[1]*(Q1[0]*Q3[2] + Q1[1]*Q3[1] + Q1[2]*Q3[0]) 
    + q3[2]*(sqr(Q3[1]) + 2*Q3[0]*Q3[2]);
	  
  gd = q3[0]*(2*Q1[1]*Q1[2]) + q3[1]*(Q1[1]*Q3[2] + Q1[2]*Q3[1]) + q3[2]*(2*Q3[1]*Q3[2]);
	  
  ge = q3(Q1[2],Q3[2]);

  ga*=d3; gb*=d3; gc*=d3; gd*=d3; ge*=d3;
	  
  // Simplify quartic by dividing by square-part of content
  cont = g_content(ga,gb,gc,gd,ge);
	  
  if(cont>1)
    {
      if(xverb) cout<<"Dividing quartic by "<<cont<<" squared\n";
      cont=sqr(cont);
      ga/=cont; gb/=cont; gc/=cont; gd/=cont; ge/=cont;
    }
	  
  bigint ga0=ga, gb0=gb, gc0=gc, gd0=gd, ge0=ge;
  quartic gg(ga,gb,gc,gd,ge);
  ggI = gg.getI();
  ggJ = gg.getJ();
  ggD = abs(gg.getdisc());
  ggbadp=plist;
  for(unsigned long ip=0; ip<plist.size(); ip++) divide_out(ggD,plist[ip]);
  int extras=(ggD>1);
  if(extras) // then we have introduced some extra bad primes
    {
      if(xverb) cout<<"Having to factorize ggD = "<<ggD<<endl;
      ggextrap = pdivs(ggD);
      ggbadp = vector_union(ggbadp,ggextrap);
    }
  if(xverb)
    {
      cout<<"Quartic gg is " << gg << endl;
      cout<<"I(gg) = "<<ggI<<endl;
      cout<<"J(gg) = "<<ggJ<<endl;
      if(extras) cout<<"extra bad primes = "<<ggextrap<<endl;
      cout<<"bad primes = "<<ggbadp<<endl;
    }
	  
// NB the full minimalization procedure ONLY works for locally soluble quartics
//    so we have to check solubility before minimalizing.
// But we can partially minimise anyway, which helps local solubility test
		  
  scaled_unimod m;
  
  if(xverb) cout << "preliminary minimalization of gg...\n";
  minim_all(ga,gb,gc,gd,ge,ggI,ggJ,ggbadp,m,0,xverb);
  gg.assign(ga,gb,gc,gd,ge);  // must reset roots before searching
  if(xverb) 
    {
      cout<<"transform "<<m<<"\n";
      cout<<"After preliminary minimalizing, gg = "<<gg<<endl;
      cout<<"I(gg) = "<<ggI<<endl;
      cout<<"J(gg) = "<<ggJ<<endl;
      if(check_transform(ga0,gb0,gc0,gd0,ge0,m,ga,gb,gc,gd,ge))
	{cout<<"transform check OK\n";}
      else
	{cout<<"transform check fails!\n";}
    }
  
  if(!locallysoluble(gg,plist,p))
    {
      if(verb) cout<<"Not locally soluble (p="<<p<<")\n";
      return -1;
    }

  if(selmer_only) return 0;

  if(xverb) 
    {
      cout << "Everywhere locally soluble, ";
      cout << "minimalizing gg...\n";
    }
  minim_all(ga,gb,gc,gd,ge,ggI,ggJ,ggbadp,m,1,xverb);
  gg.assign(ga,gb,gc,gd,ge);  // must reset roots before searching
  if(xverb) 
    {
      cout<<"transform "<<m<<"\n";
      cout<<"After minimalizing, gg = "<<gg<<endl;
      cout<<"I(gg) = "<<ggI<<endl;
      cout<<"J(gg) = "<<ggJ<<endl;
      if(check_transform(ga0,gb0,gc0,gd0,ge0,m,ga,gb,gc,gd,ge))
	{cout<<"transform check OK\n";}
      else
	{cout<<"transform check fails!\n";}
      cout<<"Descendent quartic (before reduction) = "
	  <<gg<<endl;
      cout<<"Now reducing gg...\n";
    }
  int better=1, ired;
  bigint oldga(ga), oldgb(gb), oldgc(gc), oldgd(gd), oldge(ge);
  unimod n;
  for(ired=0; (ired<5)&&better; ired++)
    {
      reduce(ga,gb,gc,gd,ge,n);
      better = (abs(ga)<=abs(oldga))&&!((gb==oldgb)&&(gc==oldgc)&&(gd==oldgd)&&(ge==oldge));
      if(better||(ired==0))
	{
	  m*=n;
	  n.reset();
	  oldga=ga; oldgb=gb; oldgc=gc; oldgd=gd; oldge=ge;
	  gg.assign(ga,gb,gc,gd,ge);
	  if(xverb) 
	    {
	      cout<<"Descendent quartic (after "<<(ired+1)
		  <<" reductions) = "<<gg<<endl;
	      if(check_transform(ga0,gb0,gc0,gd0,ge0,m,ga,gb,gc,gd,ge))
		{cout<<"transform check OK\n";}
	      else
		{cout<<"transform check fails!\n";}
	    }
	}
    }
  if(xverb) 
    {
      cout<<"transform "<<m<<"\n";
      if(check_transform(ga0,gb0,gc0,gd0,ge0,m,
			 gg.geta(),gg.getb(),gg.getcc(),gg.getd(),gg.gete()))
	{cout<<"transform check OK\n";}
      else
	{cout<<"transform check fails!\n";}
    }
  if(verb) 
    cout<<"Descendent quartic (after reduction) = "<<gg<<endl;
  
// Step 5: We have an ELS descendent quartic;
//         now attempt to find a rational point on it
  
  quartic_sieve qs(&gg,QSIEVE_OPT,0);
  if(verb) 
    cout << "Searching for points on gg up to height "<<hlim<<endl;
  
  //	      if(qs.search(hlim,1000000))
  if(!qs.search(hlim))
    {
      if(verb) cout << "No point found.\n";
      return 0;
    }
  qs.getpoint(x,y,z);
  if(verb)
    {
      show_xyz(x,y,z);
      cout<<"...mapping this point back to original quartic...\n";
    }
  bigint x3 = m(1,1)*x+m(1,2)*z;
  bigint z3 = m(2,1)*x+m(2,2)*z;
  bigint fac = gcd(x3,z3); if(fac>1) {x3/=fac; z3/=fac;}
  //	  cout<<"x3="<<x3<<"\n";
  //	  cout<<"z3="<<z3<<"\t should lie on original quartic\n";
  bigint q1xz = Q1(x3,z3);
  bigint q3xz = Q3(x3,z3);
  fac=gcd(q1xz,q3xz); if(fac>1) {q1xz/=fac; q3xz/=fac;}
  
  bigint x2 = abs(q1(q1xz,q3xz));
  bigint z2 = abs(q3(q1xz,q3xz));
  //NB These abs() are OK because x2,z2 do have the same sign
  fac=gcd(x2,z2);  if(fac>1) {x2/=fac; z2/=fac;}
  bigint y2 = q0(x2,z2);
  
  if(isqrt(x2,x)&&isqrt(z2,z)&&isqrt(y2,y))
    {
      if(verb) 
	{
	  cout<<"Point on original quartic is "; 
	  show_xyz(x,y,z);
	  cout<<endl;
	}
      return +1;
    }
  cout<<"process_d3 failed! x2="<<x2<<", z2="<<z2<<", y2="<<y2<<".\n";
  return 0;
}
