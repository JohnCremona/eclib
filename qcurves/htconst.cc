// htconst.cc:  implementations of functions for height bounds
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
 
// Here we implement (1) Silverman, (2) CPS (Cremona/Prickett/Siksek)
// bounds on the difference between naive and canonical height.

#include "marith.h"
#include "curve.h"
#include "compproc.h"
#include "matrix.h"
#include "subspace.h"
#include "cperiods.h"
#include "points.h"
#include "sieve_search.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "saturate.h"
#include "mwprocs.h"
#include "elog.h"
#include "egr.h"
#include "htconst.h"

#include "realroots.h"

// Code for Silverman bound

double logplus(double x)
{
  double ax = fabs(x);
  if(ax<1) return 0;
  return log(ax);
}

double hj(const Curvedata& CD, double& realjay)
{
  bigint c4, c6, njay, djay;
  c4=getc4(CD);
  c6=getc6(CD);
  njay = pow(c4,3);
  djay = getdiscr(CD);
  if((djay==0)||(njay==0)) {realjay=0; return 0;}

  double g = I2double(gcd(njay,djay));
  double xnjay = I2double(njay)/g;
  double xdjay = I2double(djay)/g;

  realjay = xnjay/xdjay;

  double x = log(fabs(xnjay));
  double y = log(fabs(xdjay));

  if(x<y) return y;
  else    return x;
}

double silverman_bound(const Curvedata& CD)
{
  static double log2 = log(2.0);
  bigint b2 = getb2(CD);
  bigint delta = getdiscr(CD);
  double realjay;
  double hjay = hj(CD,realjay);

// NB the constant 1.922 = 2*0.961 below is from Bremner's correction 
// to Silverman's paper; Silverman has 0.973 givin 2*0.973 = 1.946.

  double mu = 1.922  + hjay/12
                     + log(fabs(I2double(delta)))/6
                     + logplus(realjay)/6
	             + logplus(I2double(b2)/12);
  
  if(b2!=0) mu += log2;

  return mu;
}


// Cremona-Prickett-Siksek height bound, August 2002
// NB: We assume a minimal model here!

//#define DEBUG_CPS
//#define TEST_CPS

double cps_real(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8);

double egr_height_constant(const Curvedata& CD)
{ 
  return cps_real(I2bigfloat(getb2(CD)),I2bigfloat(getb4(CD)),
		     I2bigfloat(getb6(CD)),I2bigfloat(getb8(CD)));
}

double cps_bound(const Curvedata& CD)
{
  double bd = cps_real(I2bigfloat(getb2(CD)),I2bigfloat(getb4(CD)),
			  I2bigfloat(getb6(CD)),I2bigfloat(getb8(CD)));
#ifdef DEBUG_CPS
  cout<<"In cps_bound() for "<<(Curve)CD<<endl;
  cout<<"cps_real = "<<bd<<endl;
#endif
  CurveRed CR(CD);
  vector<bigint> plist = getbad_primes((Curvedata&)CD);
  for(unsigned int i=0; i<plist.size(); i++)
    {
      bigint q=plist[i];
      if(getc_p(CR,q)==1) 
	{
#ifdef DEBUG_CPS
      cout<<"q = "<<q<<", alpha = 0 since c_q=0"<<endl;
      cout<<"sum so far = "<<bd<<endl;
#endif
	  continue;      
	}
      double alpha =0;
      int m, Kc = getKodaira_code(CR,q).code;
      switch (Kc%10){
      case 0: // Im
	m = Kc/10; 
	alpha = (m%2 ? double(m*m-1)/double(4*m) : double(m)/4);
	break;
      case 1: // I*m
	m = (Kc - 1)/10;  
	alpha = (m==0? 1 : (getc_p(CR,q)==2? 1 : double(m+4)/4));
	break;
      case 3: // III
	alpha = 0.5; 
	break;
      case 4: // IV
	alpha = double(2)/3; 
	break;
      case 5: // IV*
	alpha = double(4)/3; 
	break;
      case 6: // III*
	alpha = 1.5; 
	break;
      default: // II, II*: c_p=1
	break;
      };
      bd += alpha*log(double(I2long(q)));
#ifdef DEBUG_CPS
      cout<<"q = "<<q<<", alpha = "<<alpha<<", q-term = "<< alpha*log(double(I2long(q))) <<endl;
      cout<<"sum so far = "<<bd<<endl;
#endif
    }
  return bd;
}

// Implementation for the real place originally by Nigel Smart, here
// rewritten by JC, 22/8/02

// NB the quantities which here are called dv and dvd are those which
// in Cremona, Prickett and Siksek (JNT 2006, lemma 9) are denoted e,
// e' and *not* those denoted d,d'.  So the value of egr_real() and
// egr_height_const() is -log(eps_infty)/3 (=-log(alpha) in ANTS7 paper).

bigfloat calc_dvd_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8);
bigfloat calc_dv_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8);
bigfloat old_calc_dvd_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8, const bigfloat& del);
bigfloat old_calc_dv_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8, const bigfloat& del);

double cps_real(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8)
{
  bigfloat dv=calc_dv_inf(b2,b4,b6,b8);     
  bigfloat dvd=calc_dvd_inf(b2,b4,b6,b8); 
#ifdef TEST_CPS
  bigfloat del = -b2*b2*b8-8*b4*b4*b4-27*b6*b6+9*b2*b4*b6;
  bigfloat dv2=old_calc_dv_inf(b2,b4,b6,b8,del);
  bigfloat dvd2=old_calc_dvd_inf(b2,b4,b6,b8,del); 
  bigfloat diff = dv-dv2;
  if(!is_zero(diff)) 
    cout<<"old_calc_dv gives "<<dv2<<"\nwhile new gives "<<dv<<endl<<"diff = "<<diff<<endl;
  diff = dvd-dvd2;
  if(!is_zero(diff)) 
    cout<<"old_calc_dvd gives "<<dvd2<<"\nwhile new gives "<<dvd<<endl<<"diff = "<<diff<<endl;
#endif
#ifdef DEBUG_CPS
  cout << "dv=" << dv << endl;
  cout << "dvd=" << dvd << endl;
#endif
  bigfloat htc=to_bigfloat(0);
  if(dv==-1)
    {
      if(dvd==-1) htc = to_bigfloat(0);
      else 
	{
	  if(dvd>0) htc = -log(dvd)/3;
	  else
	    {
	      cout<<"Precision problem in cps_real(): dvd = "<<dvd<<" but should be >0"<<endl;
	      cout<<"Height constant will not be correct"<<endl;
	      abort();
	      htc=0;
	    }
	}
    }
  else
    if(dvd==-1) 
      {
	if(dv>0) htc = -log(dv)/3;
	else
	  {
	    cout<<"Precision problem in cps_real(): dv = "<<dv<<" but should be >0"<<endl;
	    cout<<"Height constant will not be correct"<<endl;
	    abort();
	    htc=0;
	  }
      }
    else 
      {
	bigfloat mindv=min(dv,dvd);
	if(mindv>0) htc = -log(mindv)/3;
	else
	  {
	    cout<<"Precision problem in cps_real(): min(dv,dvd) = "<<mindv<<" but should be >0"<<endl;
	    cout<<"Height constant will not be correct"<<endl;
	    abort();
	    htc=0;
	  }
      }
  
#ifdef DEBUG_CPS
  cout<<"cps_real() returns -log(min(dv,dvd))/3 = "<<htc<<endl;
#endif

#ifdef MPFP
  double ans;
  doublify(htc,ans);
  return ans;
#else
  return htc;
#endif
}

// coeff has length 5 but may start with leading zeros
// returns real roots in [-1,1]
vector<bigfloat> roots11( const vector<bigfloat>& coeff );

inline vector<bigfloat> set_coeff(const bigfloat& c0, const bigfloat& c1, const bigfloat& c2, const bigfloat& c3, const bigfloat& c4)
{
  vector<bigfloat> coeff(5);
  coeff[0]=c0;  coeff[1]=c1;  coeff[2]=c2;  coeff[3]=c3;  coeff[4]=c4;
  return coeff;
}

inline vector<bigfloat> set_coeff(const bigfloat& c0, const bigfloat& c1, const bigfloat& c2, const bigfloat& c3)
{
  vector<bigfloat> coeff(4);
  coeff[0]=c0;  coeff[1]=c1;  coeff[2]=c2;  coeff[3]=c3;
  return coeff;
}

inline vector<bigfloat> set_coeff(const bigfloat& c0, const bigfloat& c1, const bigfloat& c2)
{
  vector<bigfloat> coeff(3);
  coeff[0]=c0;  coeff[1]=c1;  coeff[2]=c2;
  return coeff;
}

inline vector<bigfloat> set_reverse_coeff(const bigfloat& c0, const bigfloat& c1, const bigfloat& c2, const bigfloat& c3, const bigfloat& c4)
{
  return set_coeff(c4,c3,c2,c1,c0);
}

vector<bigfloat> reals_in ( vector<bigcomplex>& v);
vector<bigfloat> reals_in_11 ( vector<bigcomplex>& v);

//  Procedure to calculate dv' for the infinite prime

bigfloat calc_dvd_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8)
{ 
  bigfloat rx;
  bigfloat dvd=to_bigfloat(0),x,x2,Fx,Gx;

#ifdef DEBUG_CPS
  cout<<"\nIn calc_dvd_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif

  vector<bigfloat> rts; // holds real roots of F, G, etc
  std::set<bigfloat> crit_pts; // use a set, so we don't get any duplicates
  std::set<bigfloat> F_roots;

  crit_pts.insert(to_bigfloat(1));
  crit_pts.insert(to_bigfloat(-1));
  
#ifdef DEBUG_CPS
  cout<<"crit_pts = "<<crit_pts<<endl;
#endif

  // Put the roots of G,F',G',F+G,F-G into crit_pts:
  // Keep the roots of F separate (see comment below)

  // Roots of F
  //  rts=roots11(set_reverse_coeff(to_bigfloat(0),to_bigfloat(1),b2/4.0,b4/2.0,b6/4.0));
  rts=realroots11(set_coeff(b6,2*b4,b2,to_bigfloat(4),to_bigfloat(0)));
  F_roots.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of F: "<<rts<<endl;
  cout<<"After adding roots of F, F_roots = "<<F_roots<<endl;
#endif

  // Roots of G
  //  rts=roots11(set_reverse_coeff(to_bigfloat(1),to_bigfloat(0),-b4,-2*b6,-b8));
  rts=realroots11(set_coeff(-b8,-2*b6,-b4,to_bigfloat(0),to_bigfloat(1)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of G: "<<rts<<endl;
  cout<<"After adding roots of G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of F+G
  //  rts=roots11(set_reverse_coeff(to_bigfloat(1),to_bigfloat(4),b2-b4,2.0*(b4-b6),b6-b8));
  rts=realroots11(set_coeff(b6-b8,2.0*(b4-b6),b2-b4,to_bigfloat(4),to_bigfloat(1)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of F+G: "<<rts<<endl;
  cout<<"After adding roots of F+G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of G-F // Change from NPS
  //  rts=roots11(set_reverse_coeff(to_bigfloat(1),to_bigfloat(-4),-(b2+b4),-2.0*(b4+b6),-(b6+b8)));
  rts=realroots11(set_coeff(-(b6+b8),-2.0*(b4+b6),-(b2+b4),to_bigfloat(-4),to_bigfloat(1)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of F-G: "<<rts<<endl;
  cout<<"After adding roots of F-G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of G'
  //  rts=roots11(set_coeff(to_bigfloat(0),-2*b8,-3*b6,-b4,to_bigfloat(0)));
  rts=realroots11(set_coeff(-2*b8,-3*b6,-b4,to_bigfloat(0)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of G': "<<rts<<endl;
  cout<<"After adding roots of G', crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of F'
  //  rts=roots11(set_coeff(to_bigfloat(0),2*b6,3*b4,b2,to_bigfloat(2)));
  rts=realroots11(set_coeff(2*b6,3*b4,b2,to_bigfloat(2)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of F': "<<rts<<endl;
  cout<<"After adding roots of F', crit_pts = "<<crit_pts<<endl;
#endif

  // Evaluate max(|F(x)|,|G(x)|) at each of the x-values in array
  // crit_pts for which F(x)>=0: First take max|G(x)| over the roots
  // of F in [-1,1] -- this avoids the possibility that evaluating f
  // at one of its roots gives a slightly negative value causing a
  // possibly relevant value of G(x) to be ignored: observation of
  // Samir 25/07/04

  int first=1;
  std::set<bigfloat>::const_iterator xi = F_roots.begin();
  while(xi!=F_roots.end())
    {
      x=*xi++;
      x2=x*x;
      Gx=abs(1-b4*x2-2.0*b6*x*x2-b8*x2*x2);
#ifdef DEBUG_CPS
      cout<<"x="<<x<<", Gx="<<Gx<<endl;
#endif
      if (first)  { dvd=Gx; first=0;}
      else if (dvd>Gx) { dvd=Gx; }
#ifdef DEBUG_CPS
      cout<<"dvd so far ="<<dvd<<endl;
#endif
    }
  xi = crit_pts.begin();
  while(xi!=crit_pts.end())
    { 
      x=*xi++;
      x2=x*x;
      Fx=(4.0*x+b2*x2+2.0*b4*x*x2+b6*x2*x2);
#ifdef DEBUG_CPS
      cout<<"x="<<x<<", Fx="<<Fx<<endl;
#endif
      if(Fx<0) continue;
      Gx=abs(1-b4*x2-2.0*b6*x*x2-b8*x2*x2);
#ifdef DEBUG_CPS
      cout<<"x="<<x<<", Gx="<<Gx<<endl;
#endif
      rx=max(Fx,Gx);
      if (first)  { dvd=rx; first=0;}
      else if (dvd>rx) { dvd=rx; }
#ifdef DEBUG_CPS
      cout<<"dvd so far ="<<dvd<<endl;
#endif
    }
  if(first) return to_bigfloat(-1);

  return dvd;

}

// Procedure to calculate dv for the infinite prime

bigfloat calc_dv_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8)
{
  bigfloat rx;
  bigfloat dv=to_bigfloat(0),x,x2,fx,gx;

#ifdef DEBUG_CPS
  cout<<"\nIn calc_dv_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif

  vector<bigfloat> rts; // holds real roots of f, g, etc
  std::set<bigfloat> crit_pts; // use a set, so we don't get any duplicates
  std::set<bigfloat> f_roots;

  crit_pts.insert(to_bigfloat(1));
  crit_pts.insert(to_bigfloat(-1));
  
#ifdef DEBUG_CPS
  cout<<"crit_pts = "<<crit_pts<<endl;
#endif

  // Put the roots of g,f',g',f+g,f-g into crit_pts:
  // Keep the roots of f separate (see comment below)

  // Roots of f
  //  rts=roots11(set_coeff(to_bigfloat(0),to_bigfloat(1),b2/4.0,b4/2.0,b6/4.0));
  rts=realroots11(set_coeff(to_bigfloat(4),b2,2*b4,b6));
  f_roots.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of f: "<<rts<<endl;
  cout<<"After adding roots of f, f_roots = "<<f_roots<<endl;
#endif

  // Roots of g
  //  rts=roots11(set_coeff(to_bigfloat(1),to_bigfloat(0),-b4,-2*b6,-b8));
  rts=realroots11(set_coeff(to_bigfloat(1),to_bigfloat(0),-b4,-2*b6,-b8));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of g: "<<rts<<endl;
  cout<<"After adding roots of g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of f+g
  //  rts=roots11(set_coeff(to_bigfloat(1),to_bigfloat(4),b2-b4,2.0*(b4-b6),b6-b8));
  rts=realroots11(set_coeff(to_bigfloat(1),to_bigfloat(4),b2-b4,2.0*(b4-b6),b6-b8));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of f+g: "<<rts<<endl;
  cout<<"After adding roots of f+g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of g-f // Change from NPS
  //  rts=roots11(set_coeff(to_bigfloat(1),to_bigfloat(-4),-(b2+b4),-2.0*(b4+b6),-(b6+b8)));
  rts=realroots11(set_coeff(to_bigfloat(1),to_bigfloat(-4),-(b2+b4),-2.0*(b4+b6),-(b6+b8)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of f-g: "<<rts<<endl;
  cout<<"After adding roots of f-g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of g'
  //  rts=roots11(set_coeff(to_bigfloat(0),to_bigfloat(1),to_bigfloat(0),-b4/2.0,-b6/2.0));
  rts=realroots11(set_coeff(to_bigfloat(2),to_bigfloat(0),-b4,-b6));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of g': "<<rts<<endl;
  cout<<"After adding roots of g', crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of f'
  //  rts=roots11(set_coeff(to_bigfloat(0),to_bigfloat(0),to_bigfloat(1),b2/6,b4/6));
  rts=realroots11(set_coeff(to_bigfloat(6),b2,b4));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of f': "<<rts<<endl;
  cout<<"After adding roots of f', crit_pts = "<<crit_pts<<endl;
#endif

  // Evaluate max(|f(x)|,|g(x)|) at each of the x-values in array
  // crit_pts for which f(x)>=0: First take max|g(x)| over the roots
  // of f in [-1,1] -- this avoids the possibility that evaluating f
  // at one of its roots gives a slightly negative value causing a
  // possibly relevant value of g(x) to be ignored: observation of
  // Samir 25/07/04

  int first=1;
  std::set<bigfloat>::const_iterator xi = f_roots.begin();
  while(xi!=f_roots.end())
    {
      x=*xi++;
      x2=x*x;
      gx=abs(x2*x2-b4*x2-2.0*b6*x-b8);
#ifdef DEBUG_CPS
      cout<<"x="<<x<<", gx="<<gx<<endl;
#endif
      if (first)  { dv=gx; first=0;}
      else if (dv>gx) { dv=gx; }
#ifdef DEBUG_CPS
      cout<<"dv so far ="<<dv<<endl;
#endif
    }
  xi = crit_pts.begin();
  while(xi!=crit_pts.end())
    { 
      x=*xi++;
      x2=x*x;
      fx=(4.0*x*x2+b2*x2+2.0*b4*x+b6);
#ifdef DEBUG_CPS
      cout<<"x="<<x<<", fx="<<fx<<endl;
#endif
      if(fx<0) continue;
      gx=abs(x2*x2-b4*x2-2.0*b6*x-b8);
#ifdef DEBUG_CPS
      cout<<"x="<<x<<", gx="<<gx<<endl;
#endif
      rx=max(fx,gx);
      if (first)  { dv=rx; first=0;}
      else if (dv>rx) { dv=rx; }
#ifdef DEBUG_CPS
      cout<<"dv so far ="<<dv<<endl;
#endif
    }
  if(first) return to_bigfloat(-1);
  return dv;
}


// PAST THIS POINT the Smart method is implemented

int interval_test(const bigfloat& x, const vector<bigfloat> rts, int debug=0);

vector<bigfloat> reals_in_interval ( vector<bigcomplex>& v, const vector<bigfloat> rts);

// Inserts from C in to S the elements which are real and between -1
// and +1
void include_real_11(std::set<bigfloat>& S, const vector<bigcomplex>& C);


// test if rr lies in one of the first numint intervals out of [i1l.i1r],[i2l.i2r]

int is_in_int(const bigfloat rr,const bigfloat i1l,const bigfloat i1r,
              const bigfloat i2l,const bigfloat i2r,const long numint)
{
if (numint>0)
  { if (rr<=i1r && rr>=i1l) { return 1; }
    if (numint==2 && rr<=i2r && rr>=i2l) { return 1; }
  }
return 0;
}

// test if rr lies in one of the numint intervals  [intervals[i][0],intervals[i][1]]

int is_in_int2(const bigfloat rr,bigfloat** intervals,const long numint)
{
  long i;
  for (i=0; i<numint; i++)
    { if ((rr>=intervals[i][0]) && (rr<=intervals[i][1]))
	{ return 1; }
    }
  return 0;
}

bigfloat old_calc_dv_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8, const bigfloat& del)
{
  bigfloat rr=to_bigfloat(0),rx,rn;

#ifdef DEBUG_CPS
  cout<<"\nIn old_calc_dv_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif
  vector<bigcomplex> rt = solvecubic(b2/4,b4/2,b6/4);

  bigfloat i1l,i1u,i2l=to_bigfloat(0),i2u=to_bigfloat(0);  // Bounds on Intervals
  long i,numint;
  if (del<0)
    { for (i=0; i<3; i++)
	{ if (is_approx_zero(rt[i].imag()))
	     { rr=rt[i].real(); i=3; }
	}
      numint=0; 
      if (rr<=1.0);
	{ numint=1;  
          i1l=rr;  i1u=1.0;
	  if (rr<-1.0) { i1l=-1.0; }
	}
    }
  else
    { rn=rt[0].real();
      rx=rt[1].real();
	// Want rn<rr<rx
      if (rn>rx) { rr=rx; rx=rn; rn=rr; }
      rr=rt[2].real();
      if (rr<rn)
        { rr=rn; rn=rt[2].real(); }
      else if (rx<rr)
        { rr=rx; rx=rt[2].real(); }

      numint=2;
      i1l=rn; i1u=rr; i2l=rx; i2u=1.0;
	// Deal With Right Most Interval
      if (i2l>i2u) { numint=1; }
      if (i2l<-1.0) 
	{ numint=1;
	  i1l=-1.0; i1u=i2u;
        }
	// Now Deal With Left Most Interval
      if (i1u>1.0)  { i1u=1.0; }
      if (i1l<-1.0) { i1l=-1.0; }
      if (i1l>i1u) 
	{ numint=numint-1;
	  i1l=-1.0; i1u=i2u;
	}
    }

  if(numint==0) return to_bigfloat(-1);  // code for "infinity"

  bigfloat* te = new bigfloat[36];
  long tec=2*numint;
  if (numint>0)
    { te[0]=i1l; te[1]=i1u;
      if (numint==2)
        { te[2]=i2l; te[3]=i2u; }
    }

  if (numint!=0)
    { // Roots of g
      rt=solvequartic(to_bigfloat(0),-b4,-2.0*b6,-b8);
      for (i=0; i<4; i++)
        { if (is_approx_zero(rt[i].imag()))
	    { rr=rt[i].real();
              if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
		{ te[tec]=rr; tec=tec+1; }
            }
        }

      // Roots of f (again!)
      rt=solvecubic(b2/4.0,b4/2.0,b6/4.0);
      for (i=0; i<3; i++)
        { if (is_approx_zero(rt[i].imag()))
            { rr=rt[i].real();
              if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
                { te[tec]=rr; tec=tec+1; }
            }
        }

      // Roots of f+g
      rt=solvequartic(to_bigfloat(4),b2-b4,2.0*(b4-b6),b6-b8);
      for (i=0; i<4; i++)
        { if (is_approx_zero(rt[i].imag()))
            { rr=rt[i].real();
              if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
                { te[tec]=rr; tec=tec+1; }
            }
        }

      // Roots of g-f // Change from NPS
      rt=solvequartic(to_bigfloat(-4),-(b2+b4),-2.0*(b4+b6),-(b6+b8));
      for (i=0; i<4; i++)
        { if (is_approx_zero(rt[i].imag()))
            { rr=rt[i].real();
              if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
                { te[tec]=rr; tec=tec+1; }
            }
        }

      // Roots of f'
      if (-24.0*b4+b2*b2>=0.0)
        { rn=sqrt(-24.0*b4+b2*b2);
	  rr=(-b2+rn)/12.0;
	  if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
                { te[tec]=rr; tec=tec+1; }
	  rr=(-b2-rn)/12.0;
	  if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
                { te[tec]=rr; tec=tec+1; }
	}

      // Roots of g'
      rt=solvecubic(to_bigfloat(0),-b4/2.0,-b6/2.0);
      for (i=0; i<3; i++)
        { if (is_approx_zero(rt[i].imag()))
            { rr=rt[i].real();
              if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
                { te[tec]=rr; tec=tec+1; }
            }
        }

    }

  bigfloat dv=to_bigfloat(0);
  // Evaluate |f|, |g| at each of the tec(>0) x-values in array te:
  for (i=0; i<tec; i++)
    { rr=te[i]*te[i];
      rn=4.0*rr*te[i]+b2*rr+2.0*b4*te[i]+b6;
#ifdef DEBUG_CPS
      cout<<"x="<<te[i]<<", fx="<<rn<<endl;
#endif
      rn=abs(rn);
      rx=rr*rr-b4*rr-2.0*b6*te[i]-b8;
      rx=abs(rx);
#ifdef DEBUG_CPS
      cout<<"x="<<te[i]<<", gx="<<rx<<endl;
#endif
      if (rn>rx) { rx=rn; }
      if (i==0)  { dv=rx; }
      else if (dv>rx) { dv=rx; }
#ifdef DEBUG_CPS
      cout<<"dv so far ="<<dv<<endl;
#endif
    }
  delete[] te;

  return dv;
}

bigfloat old_calc_dvd_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8, const bigfloat& del)
{ 
  bigfloat rr=to_bigfloat(0),rn;

#ifdef DEBUG_CPS
  cout<<"\nIn old_calc_dvd_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif
  vector<bigcomplex> rt = solvecubic(b2/4,b4/2,b6/4);
  long i,j,numrt;
  bigfloat rrt[4];
  rrt[0]=to_bigfloat(0);
  if (del<0)
    { for (i=0; i<3; i++)
        { if (is_approx_zero(rt[i].imag()))
             { rr=rt[i].real(); i=3; }
        }
      if (is_approx_zero(rr)) { numrt=1; }
      else                     { rrt[1]=1.0/rr; numrt=2; }
    }
  else
    { numrt=1;
      for (i=0; i<3; i++)
	{ if (!is_approx_zero(rt[i].real()))
	    { rrt[numrt]=1.0/rt[i].real(); 
	      numrt++;
            }
	}
    }
 
  // Sort the real roots 
  long fl=0;
  while (fl==0)
    { fl=1;
      for (i=0; i<numrt-1; i++)
        { if (rrt[i]>rrt[i+1])
	    { rr=rrt[i+1];
	      rrt[i+1]=rrt[i];
	      rrt[i]=rr;
	      fl=0;
	    }
        }
    }

  // Make Intervals
  bigfloat **intervals=new bigfloat*[3];
  for (i=0; i<3; i++)
    { intervals[i]=new bigfloat[2]; }

  long numint=0;
  if (b6<0.0)
    { for (i=0; i<numrt; i=i+2)
	{ intervals[numint][1]=rrt[numrt-i-1];
	  intervals[numint][0]=rrt[numrt-i-2];
	  numint=numint+1;
	}
    }
  else
    { intervals[numint][1]=rrt[numrt-1]+2;
      intervals[numint][0]=rrt[numrt-1];
      numint=numint+1;
      for (i=1; i<numrt; i=i+2)
        { intervals[numint][1]=rrt[numrt-i-1];
	  if (numrt-i-2>0) 
	    { intervals[numint][0]=rrt[numrt-i-2]; }
	  else
	    { intervals[numint][0]=-1.0; }
          numint=numint+1;
        }
    }
  for (i=0; i<numint; i++)
    { if (intervals[i][1]>1.0)  { intervals[i][1]=1.0; }
      if (intervals[i][0]<-1.0) { intervals[i][0]=-1.0; }
    }
  i=0;
  while (i<numint)
    { fl=0;
      while (fl==0 && i<numint)
	{ fl=1;
          if (intervals[i][0]>intervals[i][1])
             { fl=0; 
	       for (j=i; j<numint-1; j++)
		  { intervals[j][0]=intervals[j+1][0];
		    intervals[j][1]=intervals[j+1][1];
		  }
	       numint=numint-1;
	     }
	}
      i=i+1;
    }

  bigfloat *te=new bigfloat[36];
  long tec=2*numint;
  for (i=0; i<numint; i++)
     { te[2*i]=intervals[i][0];
       te[2*i+1]=intervals[i][1];
     }

  if (numint!=0)
    { // Roots of g
      rt=solvequartic(to_bigfloat(0),-b4,-2.0*b6,-b8);
      for (i=0; i<4; i++)
	{ if (!is_approx_zero(rt[i]))
	    { rt[i]=to_bigfloat(1)/rt[i];
              if (is_approx_zero(rt[i].imag()))
                 { rr=rt[i].real();
                   if (is_in_int2(rr,intervals,numint))
                     { te[tec]=rr; tec=tec+1; }
		 }
            }
        }

      // Roots of f
      rr=0.0;
      if (is_in_int2(rr,intervals,numint))
                     { te[tec]=rr; tec=tec+1; }
      rt=solvecubic(b2/4.0,b4/2.0,b6/4.0);
      for (i=0; i<3; i++)
        { if (!is_approx_zero(rt[i]))
            { rt[i]=to_bigfloat(1)/rt[i];
              if (is_approx_zero(rt[i].imag()))
                 { rr=rt[i].real();
                   if (is_in_int2(rr,intervals,numint))
                     { te[tec]=rr; tec=tec+1; }
                 }
            }
        }

      // Roots of f+g
      rt=solvequartic(to_bigfloat(4),b2-b4,2.0*(b4-b6),b6-b8);
      for (i=0; i<4; i++)
        { if (!is_approx_zero(rt[i]))
            { rt[i]=to_bigfloat(1)/rt[i];
              if (is_approx_zero(rt[i].imag()))
                 { rr=rt[i].real();
                   if (is_in_int2(rr,intervals,numint))
                     { te[tec]=rr; tec=tec+1; }
		 }
            }
        }

      // Roots of f-g
      rt=solvequartic(to_bigfloat(-4),-b2-b4,2.0*(-b4-b6),-b6-b8);
      for (i=0; i<4; i++)
        { if (!is_approx_zero(rt[i]))
            { rt[i]=to_bigfloat(1)/rt[i];
	      if (is_approx_zero(rt[i].imag()))
                 { rr=rt[i].real();
                   if (is_in_int2(rr,intervals,numint))
                     { te[tec]=rr; tec=tec+1; }
		 }
            }
        }

      // Roots of f'
      rt=solvecubic(b2/2.0,3.0*b4/2.0,b6);
      for (i=0; i<3; i++)
        { if (!is_approx_zero(rt[i]))
            { if (is_approx_zero(rt[i].imag()))
                { rr=rt[i].real();
                  if (is_in_int2(rr,intervals,numint))
                    { te[tec]=rr; tec=tec+1; }
		}
            }
        }

      // Roots of g'
      rr=0.0;
      if (is_in_int2(rr,intervals,numint))
          { te[tec]=rr; tec=tec+1; }
      bigfloat dr=9.0*b6*b6-8.0*b8*b4;
      if (dr>=0.0)
	{ dr=sqrt(dr);
	  rr=(3.0*b6+dr)/(-4.0*b8);
	  if (is_in_int2(rr,intervals,numint))
                { te[tec]=rr; tec=tec+1; }
	  rr=(3.0*b6-dr)/(-4.0*b8);
          if (is_in_int2(rr,intervals,numint))
                { te[tec]=rr; tec=tec+1; }
	}
    }

  bigfloat dvd=to_bigfloat(0);
  bigfloat f,g;
  for (i=0; i<tec; i++)
    { rr=te[i];
#ifdef DEBUG_CPS
      cout<<"x="<<te[i]<<", fx="<<rr<<endl;
#endif
      f=(((b6*rr+2.0*b4)*rr+b2)*rr+4.0)*rr;
      g=(((-b8*rr-2.0*b6)*rr-b4)*rr)*rr+1.0;
      f=abs(f); g=abs(g);
      rn=f;
#ifdef DEBUG_CPS
      cout<<"x="<<te[i]<<", gx="<<rn<<endl;
#endif
      if (g>f) { rn=g; }
      if (i==0)  { dvd=rn; }
      else if (dvd>rn) { dvd=rn; }
#ifdef DEBUG_CPS
      cout<<"dvd so far ="<<dvd<<endl;
#endif
    }
  delete[] te;
  for (i=0; i<3; i++)
    { delete[] intervals[i]; }
  delete[] intervals;

  return dvd;

}

// coeff has length 5 but may start with leading zeros
// returns real roots in [-1,1]
vector<bigfloat> roots11( const vector<bigfloat>& coeff )
{
#ifdef DEBUG_CPS
  cout<<"In roots11() with coeff =  "<<coeff<<endl;
#endif
  bigfloat a=coeff[0], a1;
#ifdef DEBUG_CPS
  cout<<"coeff[0] = a = "<<a<<endl;
#endif
  if(a!=0) 
    {
#ifdef DEBUG_CPS
      cout<<"inverting nonzero a ="<<a<<endl;
#endif
      a1=1/a;
      vector<bigcomplex> cr = solvequartic(a1*coeff[1],a1*coeff[2],a1*coeff[3],a1*coeff[4]);
      return reals_in_11(cr);
    }
  a=coeff[1];
#ifdef DEBUG_CPS
  cout<<"coeff[1] = a = "<<a<<endl;
#endif
  if(a!=0) 
    {
#ifdef DEBUG_CPS
      cout<<"inverting nonzero a ="<<a<<endl;
#endif
      a1=1/a;
      vector<bigcomplex> cr = solvecubic(a1*coeff[2],a1*coeff[3],a1*coeff[4]);
      return reals_in_11(cr);
    }
  vector<bigfloat> ans; // zero length
  a=coeff[2];
#ifdef DEBUG_CPS
  cout<<"coeff[2] = a = "<<a<<endl;
#endif
  if(a==0) {cout<<"Error in roots: degree<2"<<endl; abort();}
  bigfloat b=coeff[3], c=coeff[4], x;
  bigfloat d=b*b-4*a*c;
  if(d>=0) 
    {
      d=sqrt(d);  a1=1/(2*a);
      x=a1*(-b+d);
      if((x<=1)&&(x>=-1)) ans.push_back(x);
      x=a1*(-b-d);
      if((x<=1)&&(x>=-1)) ans.push_back(x);
    }
  return ans;
}

vector<bigfloat> reals_in ( vector<bigcomplex>& v)
{
  vector<bigfloat> vr;
  vector<bigcomplex>::iterator vi = v.begin();
  while(vi!=v.end()) 
    {
      if(is_real(*vi)) vr.push_back((*vi).real());
      vi++;
    }
  return vr;
}

vector<bigfloat> reals_in_11 ( vector<bigcomplex>& v)
{
  vector<bigfloat> vr;
  vector<bigcomplex>::iterator vi = v.begin();
  while(vi!=v.end()) 
    {
      if(is_real(*vi)) 
	{
	  bigfloat x = (*vi).real();
	  if((x<=1)&&(x>=-1)) vr.push_back(x);
	}
      vi++;
    }
  return vr;
}

int interval_test(const bigfloat& x, const vector<bigfloat> rts, int debug)
{
// rts will have size 1 or 3, and be ordered
  if(debug) cout<<"Interval test("<<x<<"), rts="<<rts<<endl;
  if(x>1) {if(debug) cout<<"\t returns 0\n"; return 0;}
  if(x<-1) {if(debug) cout<<"\t returns 0\n";  return 0;}
  int ans;
  if(rts.size()==1) 
    ans = (x>=rts[0]);
  else 
    ans =  ((x>=rts[0]) && (x<=rts[1])) || (x>=rts[2]); 
  if(debug) cout<<"\t returns "<<ans<<"\n";
  return ans;
}

vector<bigfloat> reals_in_interval ( vector<bigcomplex>& v, const vector<bigfloat> rts)
{
  vector<bigfloat> vr;
  vector<bigcomplex>::iterator vi = v.begin();
  bigfloat x;
  while(vi!=v.end()) 
    {
      if(is_real(*vi)) 
	{
	  x=(*vi).real();
	  if(interval_test(x,rts,1)) vr.push_back(x);
	}
      vi++;
    }
  return vr;
}

// Inserts from C in to S the elements which are real and between -1
// and +1
void include_real_11(std::set<bigfloat>& S, const vector<bigcomplex>& C)
{
  bigfloat x;
  vector<bigcomplex>::const_iterator Ci = C.begin();
  while(Ci!=C.end())
    {
      if(is_real(*Ci))
	{
	  x = (*Ci).real();
	  if((x<=1)&&(x>=-1)) S.insert(x);
	}
      Ci++;
    }
}

inline int in_fund_region(const bigcomplex& z)
{
  return (imag(z)>0)&&(abs(z)>0.999)&&(abs(2*real(z))<1.001);
}

//#define HTB_DEBUG

// This gives a lower bound on non-torsion points; experimental, not
// very good (and not always positive?!!

bigfloat lower_height_bound_alt(const Curvedata& CD)
{
  static bigfloat log2 = log(to_bigfloat(2));
  bigcomplex w1, w2;
  Cperiods pers(CD);
  bigcomplex tau = gettau(pers);
#ifdef HTB_DEBUG
  cout<<"\nAfter normalizing, periods are: "<<pers<<endl;
  if(!in_fund_region(tau)) 
    cout<<"tau is not in the fundamental region!"<<endl;
#endif  
  bigfloat aq = abs(q(tau));
  pers.getwi(w1,w2);  //NB w1 is the minimal period here, not necessarily the real period!
  int realcase = abs(imag(w1))<0.0001;
#ifdef HTB_DEBUG
  cout<<"q  = "<<q(tau)<<endl;
  cout<<"aq = "<<aq<<endl;
  cout<<"w1 = "<< w1 << endl;
  cout<<"w2 = "<< w2 << endl;
  cout<<"real case?  ";if(realcase)cout<<"yes";else cout<<"no"; cout<<endl;
#endif  

  // Either (when w1 is real):
  // (1/6)log|Delta| - 2log2 - (1/6)log|q| - 4(|q|)/(1-|q|)

  // Or (when w1 is not real):
  // (1/6)log|Delta| - 2log2 + (1/12)log|q| - 2(1+|q|^2)/(1-|q|)

  bigfloat term1 = log(abs(I2bigfloat(getdiscr(CD))))/6;
  bigfloat term2 = - 2*log2;
  bigfloat term3, term4;

  if(realcase) 
    {
      term3 = -log(aq)/6;
      term4 = -4*(aq)/(1-aq);
    }
  else
    {
      term3 = log(aq)/12;
      term4 = -2*(1+aq*aq)/(1-aq);
    }
  bigfloat bound = term1 + term2 + term3 + term4;
#ifdef HTB_DEBUG
  cout<<"term 1 = "<<term1<<endl;
  cout<<"term 2 = "<<term2<<endl;
  cout<<"term 3 = "<<term3<<endl;
  cout<<"term 4 = "<<term4<<endl;
#endif
  return bound;
}

// This gives a lower bound on non-torsion points, by searching, given
// the regulator of a known subgroup

// If point search bound is greater than this, output a warning
// message and reduce to this value:
const int max_search_bound = 18;

bigfloat lower_height_bound_search(const Curvedata& CD, const bigfloat& reg)
{
  const int verbose=0;
  // Find optimally x-shifted curve for better point searching...
  bigint x_shift;
  Curvedata C_opt = opt_x_shift(CD,x_shift);
  int shift_flag = !is_zero(x_shift);
  if(shift_flag&&verbose) 
    cout<<"Using shifted model "<<(Curve)C_opt<<" for searching"<<endl;

  double hc = height_constant(C_opt);
  if(verbose) 
    cout<<"height bound constant for shifted curve  = "<<hc<<endl;
  double hc1; 
  doublify(reg,hc1); 
  hc1 = (hc+hc1/(3.9));
  if(hc1>12) hc1=12;    // so hc1 = min(12,R/4+ht.const.)
  double hcx = hc1-hc; // = min(12-ht.const., R/4)
  if(hcx<0) {hcx=0.1; hc1=hcx+hc;}
  if(verbose)
    {
      cout<<"Searching for all points to naive height "<<hc1<<endl;
    }
  if(hc1>max_search_bound) 
    {
      cout<<"\n***Warning: search bound of "<<hc1
	  <<" reduced to "<<max_search_bound
	  <<" -- points may not be saturated***"<<endl;     
      hc1=max_search_bound;
    }
  point_min_height_finder pmh(&C_opt,0,verbose);
  pmh.search(to_bigfloat(hc1));
  bigfloat lambda=pmh.get_min_ht();
  Point Pmin = pmh.get_min_ht_point();
  if(lambda==0)
    {	  
      lambda=hcx;
      if(verbose) 
	cout<<"No points found, lambda = "<<lambda<<endl;
    }
  else
    {
      if(verbose) 
	cout<<"Min height of points found = "<<lambda<<" (point "<<Pmin<<")"<<endl;
      if(lambda>hcx) lambda=hcx;
      if(verbose) 
	cout<<"Using lambda = "<<lambda<<endl;
    }
  return lambda;
}

static long fact_tab[13] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};

long factorial(long n) // OK for n up to 12
{
  if(n<2) return 1;
  if(n>12) 
    {
      cout<<"factorial(long) called with n = "<<n<<" -- too big!"<<endl;
      abort();
    }
  //  return n*factorial(n-1);
  return fact_tab[n];
}

bigfloat factorial(bigfloat n)
{
  if(n<2) return to_bigfloat(1);
  if(n<13) return to_bigfloat(fact_tab[longify(n)]);
  return n*factorial(n-1);
}

bigfloat lattice_const(int r)
{
  static bigfloat one   = to_bigfloat(1);
  static bigfloat two   = to_bigfloat(2);
  static bigfloat three = to_bigfloat(3);
  static bigfloat four  = to_bigfloat(4);
  static bigfloat five  = to_bigfloat(5);
  static bigfloat six   = to_bigfloat(6);
  static bigfloat seven = to_bigfloat(7);
  static bigfloat c1 = one;
  static bigfloat c2 = sqrt(four/three);
  static bigfloat c3 = pow(two,one/three);
  static bigfloat c4 = sqrt(two);
  static bigfloat c5 = pow(two,three/five);
  static bigfloat c6 = two/pow(three,one/six);
  static bigfloat c7 = pow(two,six/seven);
  static bigfloat c8 = two;
  switch(r) {
  case 0: return c1;
  case 1: return c1;
  case 2: return c2;
  case 3: return c3;
  case 4: return c4;
  case 5: return c5;
  case 6: return c6;
  case 7: return c7;
  case 8: return c8;
  default: 
    {
      bigfloat xr = to_bigfloat(r);
      if(r%2)
	{
	  return factorial(xr)*sqrt(Pi())/(factorial((xr-1)/2)*(1<<r));
	}
      else // r is even
	{
	  return factorial(xr/2);
	}
    }
  }
}

point_min_height_finder::point_min_height_finder(Curvedata* EE, int egr, int verb)
  :E(EE), egr_flag(egr), verbose(verb) 
{
  min_ht=to_bigfloat(0);
  Pmin=Point(E);
  E -> getai(a1,a2,a3,a4,a6);
  if(egr_flag) CG=ComponentGroups(*E);
  iso = !((a1==0)&&(a3==0));
  c.resize(4);
  if(iso)
    {
      c[0]=16*getb6(*E);
      c[1]= 8*getb4(*E);
      c[2]=   getb2(*E);
      c[3]=1;
    }
  else
    {
      c[0]=a6;
      c[1]=a4;
      c[2]=a2;
      c[3]=1;
    }
}

int point_min_height_finder::process(const bigint& x, const bigint& y, const bigint& z) 
{
  bigint rz; isqrt(z,rz);
  bigint x1=x*rz, y1=y, z1=z*rz;
  if(iso)
    {
      y1 -= (a1*x1+4*a3*z1);
      x1 *= 2;
      z1 *= 8;
    }
  Point P(E, x1,y1,z1);
  if(P.isvalid())
    {
      if(order(P)<0) {
      int egr=1; bigint p0;
      if(egr_flag) egr=CG.HasGoodReduction(P,p0);
      if(egr)
	{
	  bigfloat hP=height(P);
	  if(is_zero(hP)) return 0;
	  if(verbose) 
	    cout<<"Found point "<<P<<" with height "<<hP<<endl;
	  all_points.push_back(P);
	  
	  if((min_ht==0)||(hP<min_ht)) 
	    {
	      if(verbose) 
		cout<<"New minimum height = "<<hP<<endl;
	      min_ht=hP;
	      Pmin=P;
	    }
	}
      else
	{
	  if(verbose)
	    cout<<"Found point "<<P
		<<" but ignoring as not egr (bad reduction at "<<p0<<")"<<endl;
	}
      }
    }
  else
    {
      cout<<"Raw point       x,y,z = "<<x<<", "<<y<<", "<<z<<endl;
      cout<<"converted point P = "<<P<<" --not on curve!"<<endl;
    }
  return 0;
}

void point_min_height_finder::search(bigfloat h_lim)
{
    if(iso) h_lim+=2.08;
//  if(iso) cout<<"Adding log(8) to h_lim, increasing it to "<<h_lim<<endl;
    qsieve s(this, 3, c, h_lim, (verbose>1)); 
    bigcomplex c1(I2bigfloat(c[2])),
	c2(I2bigfloat(c[1])),
	c3(I2bigfloat(c[0]));
    vector<bigcomplex> roots=solvecubic(c1,c2,c3);
    //  cout<<"solvecubic("<<c1<<","<<c2<<","<<c3<<") returns "<<roots<<endl;
    vector<double> bnd(3);
    int nrr=order_real_roots(bnd,roots);
    s.set_intervals(bnd,nrr,1);
    s.search();
}

//#define debugLB

///////////////////////////////////////////////////////////////////////
//
// class Interval represents a closed interval [lh,rh] where either
// empty=1; or empty=0 and lh <= rh; flags rhinf, lhinf denote
// rh=infty and lh=-infty resp.
//
///////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& os, const Interval& I)
{
  if(I.empty) os<<"[]"; 
  else 
    {
      os<<"[";
      if(I.lhinf) os << "-infty"; else os << I.lh;
      os << ",";
      if(I.rhinf) os << "+infty"; else os << I.rh;
      os << "]";
    }
  return os;
}

void Interval::intersect(const Interval& I)
{
  if(empty) return;
  if(I.empty) {empty=1; return;}
  if(lhinf) 
    {
      lhinf=I.lhinf; lh=I.lh;
    }
  else if(!I.lhinf) lh=max(lh,I.lh); 
  if(rhinf) 
    {
      rhinf=I.rhinf; rh=I.rh;
    }
  else if(!I.rhinf) rh=min(rh,I.rh); 
  if((!lhinf)&&(!rhinf)&&(lh>rh)) empty=1;
}

vector<Interval> intersect(const vector<Interval>& L1, const vector<Interval>& L2)
{
  vector<Interval> ans;
  vector<Interval>::const_iterator I, J;
  for(I=L1.begin(); I!=L1.end(); I++)
    for(J=L2.begin(); J!=L2.end(); J++)
      {
	Interval K = intersect(*I,*J);
	if(!K.is_empty()) ans.push_back(K);
      }
  return ans;
}


///////////////////////////////////////////////////////////////////////
//
// class Interval01 represents a closed subinterval [lh,rh] of [0,1],
// where either empty=1; or empty=0 and lh <= rh.
//
///////////////////////////////////////////////////////////////////////

ostream& operator<< (ostream& os, const Interval01& I)
{
  if(I.empty) 
    os<<"[]"; 
  else 
    os<<"[" << I.lh << "," << I.rh << "]";
  return os;
}

void Interval01::intersect(const Interval01& I)
{
  if(empty) return;
  if(I.empty) {empty=1; return;}
  lh=max(lh,I.lh); 
  rh=min(rh,I.rh); 
  empty=(lh>rh);
}

vector<Interval01> intersect(const vector<Interval01>& L1, 
			     const vector<Interval01>& L2)
{
  vector<Interval01> ans;
  vector<Interval01>::const_iterator I, J;
  for(I=L1.begin(); I!=L1.end(); I++)
    for(J=L2.begin(); J!=L2.end(); J++)
      {
	Interval01 K = intersect(*I,*J);
	if(!K.is_empty()) ans.push_back(K);
      }
  return ans;
}

Interval01 operator/(const Interval01& I, const long n)
{
  if(I.empty) return I;
  return Interval01(I.lh/to_bigfloat(n),I.rh/to_bigfloat(n));
}

Interval01 operator+(const Interval01& I, const bigfloat& shift)
{
  if(I.empty) return I;
  return Interval01(I.lh+shift,I.rh+shift);
}

///////////////////////////////////////////////////////////////////////
//
// Implementation of class CurveHeightConst
//
///////////////////////////////////////////////////////////////////////

CurveHeightConst::CurveHeightConst(const Curvedata& CD)
  : Curvedata(CD), Cperiods(CD)
{
  c = to_bigfloat(egr_height_constant(*this)); // =-log(alpha) in ANTS7
  e3 = get_e3();
  n_max=10;
  n_ann=25; // i.e. all p<100
  ann = annihilators(*this,n_ann);
#ifdef debugLB
  cout<<"e3 = "<<e3<<endl;
  cout<<"archContrib = log(epsilon)/3 = "<<c<<endl;
  cout<<"annihilators = "<<ann<<endl;
  cout<<"n_max = "<<n_max<<endl;
#endif
}

bigfloat CurveHeightConst::D(const long n) const // "denomContrib"
{
  bigfloat denomContrib = to_bigfloat(0);  // = D_E(n) in the paper
  primevar pr; 
  long pmax=(n+1)*(n+1), p=pr;
  for(int i=0; (i<n_ann)&&(p<pmax); i++, pr++)       
    if(n%ann[i]==0) 
      {
#if(0)
//#ifdef debugLB
	cout<<"i="<<i<<", p="<<p
	    <<", ann[i]="<<ann[i]<<"; adding "<<(2*log(p))
	    <<" to denomContrib"<<endl;
#endif
	denomContrib+=2*log((double)p);
	pr++; p=(long)pr;
      }  
  return denomContrib;
}

void CurveHeightConst::compute_phase1()
{
  int success;
  bigfloat target=to_bigfloat(1), fac=to_bigfloat(2);

  // Step 1: find a value of target which succeeds but fac*target fails.
  // The former is called "lower" and the latter "upper":

  success=test_target(target,n_max);

  if(success) 
    {
      lower=target;
      while(success)
	{
	  target*=fac;
	  success=test_target(target,n_max);
	}
      upper=target; lower=target/fac;
    }
  else
    {
      upper=target;		 
      while(!success)
	{
	  target/=fac;
	  n_max+=5;
	  success=test_target(target,n_max);
	}
      lower=target; upper=target*fac;
    }
#ifdef debugLB
  cout<<"Initial interval for LB = ["<<lower<<","<<upper<<"]"<<endl;
#endif
}

void CurveHeightConst::compute_phase2()
{
  int success;

  // Step 2: repeatedly test lower*sqrt(fac) and replace either lower or
  // upper by it and replace fac by sqrt(fac).
  
  int i,nsteps=1000;  // just an upper bound
  bigfloat tolerance=to_bigfloat(0.001); // will stop when |upper-lower|<tolerance
  bigfloat target=to_bigfloat(1), fac=to_bigfloat(2);
  for(i=0; i<nsteps; i++)
    {
      fac=sqrt(fac);
      target=lower*fac;
      success=test_target(target,n_max);
      if(success) {lower=target;} else {upper=target;}
#ifdef debugLB
      cout<<"After "<<(i+1)<<" refinements, interval for LB = ["
	  <<lower<<","<<upper<<"]"<<endl;
#endif
      if(upper-lower<tolerance) return;
    }
}

int CurveHeightConst::test_target(const bigfloat& target, long k) 
{ 
  for(int n=1; n<k; n++)
    if(Bnmu(n,target) < to_bigfloat(1)) 
      return 1;
  return canonicalHeightInterval01(target,k).size()==0;
}


vector<Interval01> CurveHeightConst::canonicalHeightInterval01(const bigfloat& target, long k)
{
  vector<Interval01> solution;
  solution.push_back(Interval01(to_bigfloat(0.5),to_bigfloat(1))); 
  // i.e. image of [e3,+infty] under psi
#ifdef debugLB
  cout<<"testing target bound "<<target<<" with k = "<<k<<endl;
  cout<<"Starting interval: "<<solution<<endl;
#endif
  for(int n=1; n<=k; n++)
    {
#ifdef debugLB
     cout<<"n = "<<n<<endl;
#endif
     bigfloat B=Bnmu(n,target);     // = B_n(target) in the paper
#ifdef debugLB
     cout<<"B_"<<n<<"("<<target<<") = "<<B<<endl;
#endif
     if(B > 1.0e100) break;
     vector<Interval01> t=solveLEQ01(n,B);
     solution=intersect(solution,t);
#ifdef debugLB
     cout<<"intervals from solveLEQ01: "<<t<<endl;
     cout<<"intervals now: "<<solution<<endl;
#endif
     if(solution.size()==0) return solution;
     t=solveGEQ01(n,-B);
     solution=intersect(solution,t);
#ifdef debugLB
     cout<<"intervals from solveGEQ01: "<<t<<endl;
     cout<<"intervals now: "<<solution<<endl;
#endif
     if(solution.size()==0) return solution;
   }
 return  solution;
}

// Returns a list of subintervals of [0,1] containing the elliptic
// logs of P for which x(nP) <= B

vector<Interval01> CurveHeightConst::solveLEQ01(long n, const bigfloat& B)
{
#ifdef debugLB
  cout<<"solveLEQ01("<<n<<","<<B<<")"<<endl;
#endif
  vector<Interval01> ans;
  if(B < e3) return ans;
#ifdef debugLB
  cout<<"B = "<<B<<endl;
#endif
  bigfloat x0=psi(B);
#ifdef debugLB
     cout<<"x0 = psi(B) = "<<x0<<endl;
#endif
  bigfloat oneovern = to_bigfloat(1)/to_bigfloat(n);
  Interval01 I(1-x0,x0);  I=I/n;
#ifdef debugLB
     cout<<"interval: "<<I<<endl;
#endif
  for(int i=0; i<n; i++, I=I+oneovern) {ans.push_back(I);}
  return ans;
}

// Returns a list of subintervals of [0,1] containing the elliptic
// logs of P for which x(nP) >= B

vector<Interval01> CurveHeightConst::solveGEQ01(long n, const bigfloat& B)
{
  if(B <= e3) 
    {
      vector<Interval01> ans;
      ans.push_back(Interval01()); // i.e.[0,1]
      return ans;
    }
  vector<Interval01> ans;
  bigfloat x0=psi(B);
  bigfloat oneovern = to_bigfloat(1)/to_bigfloat(n);
  Interval01 I(to_bigfloat(0),1-x0);  I=I/n;
  for(int i=0; i<n; i++, I=I+oneovern) {ans.push_back(I);}
  I=Interval01(x0,to_bigfloat(1));  I=I/n;
  for(int i=0; i<n; i++, I=I+oneovern) {ans.push_back(I);}
  return ans;
}

// mimic gp's ellordinate(): given a real x, returns a vector of
// length 0,1 or 2 containing those y for which [x,y] is on the curve;
// if there are two such values, the largest is first.

vector<bigfloat> CurveHeightConst::ordinates(const bigfloat& x)
{
  vector<bigfloat> ans;
  static const bigfloat four=to_bigfloat(4), two=to_bigfloat(2);
  bigfloat d = ((four*x+I2bigfloat(b2))*x+(two*I2bigfloat(b4)))*x+I2bigfloat(b6);
  if(d<0) return ans;
  bigfloat y = -(I2bigfloat(a1)*x+I2bigfloat(a3))/two;
  if(d==0) {ans.push_back(y); return ans;}
  d=sqrt(d)/two;  // positive
  ans.push_back(d+y);  // the larger value
  ans.push_back(-d+y); // the smaller value
  return ans;
}

// elliptic log function (called psi in the paper) with domain
// [e3,infty], codomain [0.5,1]

bigfloat CurveHeightConst::psi(const bigfloat& x)
{
  if(x<e3) 
    {
      cout<<"Error in CurveHeightConst::psi(): x="<<x<<" < e3 = "<<e3<<endl;  
      abort();
    }
  //  cout<<"computing psi(x) with x = "<<x<<endl;
  //  cout<<"ordinates: "<<ordinates(x)<<endl;
  bigfloat y = ordinates(x)[0];
  //  cout<<"y = "<<y<<endl;
  bigcomplex z = pointtoz(x,y);
  //  cout<<"z = "<<z<<endl;
  return real(z/get_real_period()); // in [0.5,1]
}

// returns the order of the reduction of CD mod p (i.e. of E^0(Qp)/E^1(Qp))
long annihilator(CurveRed& CR, long p)
{
  bigint pp = BIGINT(p);
  bigint tr = Trace_Frob(CR,pp);
  long ans = p-I2long(tr);
  if(getord_p_N(CR,pp)==0) ans++; 
  return ans;
}

vector<long> annihilators(const Curvedata& CD, long n)
{
  vector<long> ans;
  primevar pr;
  CurveRed CR(CD);
  while(n--) {ans.push_back(annihilator(CR,(long)pr)); pr++;}
  return ans;
}

#if(0) // old code using x-intervals instead of [0,1]-intervals

vector<Interval> CurveHeightConst::canonicalHeightInterval(const bigfloat& target, long k)
{
  vector<Interval> solution;
  solution.push_back(Interval(e3)); // i.e. [e3,+infty]
#ifdef debugLB
  cout<<"testing target bound "<<target<<" with k = "<<k<<endl;
  cout<<"Starting interval: "<<solution<<endl;
#endif
  for(int n=1; n<=k; n++)
    {
     bigfloat B=Bnmu(n,target);      // = B_n(target) in the paper
#ifdef debugLB
     cout<<"n = "<<n<<endl;
     cout<<"B_"<<n<<"("<<target<<") = "<<B<<endl;
#endif
     if(B > 1.0e100) break;
     vector<Interval> t=solveLEQ(n,B);
     solution=intersect(solution,t);
#ifdef debugLB
     cout<<"intervals from solveLEQ: "<<t<<endl;
     cout<<"intervals now: "<<solution<<endl;
#endif
     if(solution.size()==0) return solution;
     t=solveGEQ(n,-B);
     solution=intersect(solution,t);
#ifdef debugLB
     cout<<"intervals from solveGEQ: "<<t<<endl;
     cout<<"intervals now: "<<solution<<endl;
#endif
     if(solution.size()==0) return solution;
   }
 return  solution;
}

// Returns a list of intervals containing those x(P) for which x(nP) <= B

vector<Interval> CurveHeightConst::solveLEQ(long n, const bigfloat& B)
{
  vector<Interval> ans;
  if(B < e3) return ans;
  vector<bigfloat> xlist=solveEllNPower(n,B);
  if(n%2) // n is odd
    {
      ans.push_back(Interval(e3,xlist[0]));
      for(int i=1; i<n; i+=2)
	ans.push_back(Interval(xlist[i],xlist[i+1]));
    }	
  else // n is even
    {
      for(int i=0; i<n; i+=2)
	ans.push_back(Interval(xlist[i],xlist[i+1]));      
    }  
  return ans;
}

// Returns a list of intervals containing those x(P) for which x(nP) >= B

vector<Interval> CurveHeightConst::solveGEQ(long n, const bigfloat& B)
{
  vector<Interval> ans;
  if(B <= e3) 
    {
      ans.push_back(Interval(e3)); // i.e.[e3,infty]
      return ans;
    }
  vector<bigfloat> xlist=solveEllNPower(n,B);
  int i;
  if(n%2) // n is odd
    {
      for(i=1; i<n; i+=2)
	ans.push_back(Interval(xlist[i-1],xlist[i]));
      ans.push_back(Interval(xlist[n-1]));
    }	
  else // n is even
    {
      ans.push_back(Interval(e3,xlist[0]));
      for(i=2; i<n; i+=2)
	ans.push_back(Interval(xlist[i-1],xlist[i]));      
      ans.push_back(Interval(xlist[n-1]));
    }  
  return ans;
}

// solves x(nP)=x1 on $E_0$; returns vector of n reals in [0,1)
// (elliptic logs of solutions)

vector<bigfloat> CurveHeightConst::solveEllNPower01(long n, const bigfloat& x1)
{
  vector<bigfloat> ans;
  if(x1<e3)  return ans; // empty
  bigfloat z=psi(x1);
#ifdef debugLB
  cout<<"psi("<<x1<<") = "<<z<<endl;
#endif
  bigfloat oneovern = to_bigfloat(1)/to_bigfloat(n);
  z*=oneovern;
  for(int i=0; i<n; i++, z+=oneovern) ans.push_back(z);
#ifdef debugLB
  cout<<"solveEllNPower01("<<n<<","<<x1<<") returns "<<ans<<endl;
#endif
  return ans;
}

// solves x(nP)=x1 on $E_0$; returns vector of n reals

vector<bigfloat> CurveHeightConst::solveEllNPower(long n, const bigfloat& x1)
{
  vector<bigfloat> ans=solveEllNPower01(n,x1);
  bigfloat om=real(get_real_period());
  for(int i=0; i<n; i++) ans[i]=real(ztopoint(om*ans[i])[0]);
  sort(ans.begin(), ans.end());
#ifdef debugLB
  cout<<"solveEllNPower("<<n<<","<<x1<<") returns "<<ans<<endl;
#endif
  return ans;
}

#endif
