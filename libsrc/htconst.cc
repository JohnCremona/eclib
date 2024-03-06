// htconst.cc:  implementations of functions for height bounds
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
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
 
// Here we implement (1) Silverman, (2) CPS (Cremona/Prickett/Siksek)
// bounds on the difference between naive and canonical height.

#include <eclib/mwprocs.h> // only needed for order_real_roots
#include <eclib/htconst.h>

#include <eclib/realroots.h>

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
  double bd = cps_real(I2bigfloat(getb2(CD)),I2bigfloat(getb4(CD)),
		     I2bigfloat(getb6(CD)),I2bigfloat(getb8(CD)));
  if (abs(bd) < 1e-30) bd=0; // otherwise the output sometimes prints as "-0"
  return bd;
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
  for( const auto& q : plist)
    {
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
  if (abs(bd) < 1e-30) bd=0; // otherwise the output sometimes prints as "-0"
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
  bigfloat zero=to_bigfloat(0);
  bigfloat htc=zero, dv=zero, dvd=zero;
#ifdef MPFP
  long original_prec, prec;
  prec = original_prec = bit_precision();
  while (dv==0 || dvd==0)
    {
      dv=calc_dv_inf(b2,b4,b6,b8);
      dvd=calc_dvd_inf(b2,b4,b6,b8);
      if (dv==0 || dvd==0)
        {
#ifdef DEBUG_CPS
          cout << "In cps_real(), with bit precision " << prec << " we have dv="<<dv<<", dvd="<<dvd << endl;
#endif
          prec *=2;
          set_bit_precision(prec); // does not change output precision
#ifdef DEBUG_CPS
          cout << " -- doubling bit precision to " << prec << endl;
#endif
        }
    }
  if (prec!=original_prec)
    {
      set_bit_precision(original_prec); // does not change output precision
#ifdef DEBUG_CPS
      cout << " --resetting bit precision back to " << original_prec << endl;
#endif
    }

#else // using C double precision it may be impossible to compute dv, dvd
  dv=calc_dv_inf(b2,b4,b6,b8);
  dvd=calc_dvd_inf(b2,b4,b6,b8);
  if (dv==0 || dvd==0)
    {
      cout << "In cps_real(), using C doubles we have dv="<<dv<<", dvd="<<dvd << endl;
      cout << "Unable to compute height constant." << endl;
      return htc;
    }
#endif

#ifdef DEBUG_CPS
  cout << "dv=" << dv << endl;
  cout << "dvd=" << dvd << endl;
#endif
  if(dv==-1)
    {
      if(dvd==-1) htc = to_bigfloat(0);
      else 
	{
	  if(dvd>0) htc = -log(dvd)/3;
	  else
	    {
	      cerr<<"Precision problem in cps_real(): dvd = "<<dvd<<" but should be >0"<<endl;
	      cerr<<"Height constant will not be correct"<<endl;
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
	    cerr<<"Precision problem in cps_real(): dv = "<<dv<<" but should be >0"<<endl;
	    cerr<<"Height constant will not be correct"<<endl;
	    htc=0;
	  }
      }
    else 
      {
	bigfloat mindv=min(dv,dvd);
	if(mindv>0) htc = -log(mindv)/3;
	else
	  {
	    cerr<<"Precision problem in cps_real(): min(dv,dvd) = "<<mindv<<" but should be >0"<<endl;
	    cerr<<"Height constant will not be correct"<<endl;
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
  vector<bigfloat> coeff = {c0,c1,c2,c3,c4};
  return coeff;
}

inline vector<bigfloat> set_coeff(const bigfloat& c0, const bigfloat& c1, const bigfloat& c2, const bigfloat& c3)
{
  vector<bigfloat> coeff = {c0,c1,c2,c3};
  return coeff;
}

inline vector<bigfloat> set_coeff(const bigfloat& c0, const bigfloat& c1, const bigfloat& c2)
{
  vector<bigfloat> coeff = {c0,c1,c2};
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
  bigfloat zero=to_bigfloat(0), one=to_bigfloat(1), two=to_bigfloat(2), three=to_bigfloat(3), four=to_bigfloat(4);
  bigfloat rx;
  bigfloat dvd=zero,x2,Fx,Gx;

#ifdef DEBUG_CPS
  cout<<"\nIn calc_dvd_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif

  vector<bigfloat> rts; // holds real roots of F, G, etc
  std::set<bigfloat> crit_pts; // use a set, so we don't get any duplicates
  std::set<bigfloat> F_roots;

  crit_pts.insert(one);
  crit_pts.insert(-one);

#ifdef DEBUG_CPS
  cout<<"crit_pts = "<<crit_pts<<endl;
#endif

  // Put the roots of G,F',G',F+G,F-G into crit_pts:
  // Keep the roots of F separate (see comment below)

  // Roots of F
  rts=realroots11(set_coeff(b6,2*b4,b2, four,zero));
  F_roots.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of F: "<<rts<<endl;
  cout<<"After adding roots of F, F_roots = "<<F_roots<<endl;
#endif

  // Roots of G
  rts=realroots11(set_coeff(-b8,-two*b6,-b4,zero,one));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of G: "<<rts<<endl;
  cout<<"After adding roots of G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of F+G
  rts=realroots11(set_coeff(b6-b8,two*(b4-b6),b2-b4, four,one));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of F+G: "<<rts<<endl;
  cout<<"After adding roots of F+G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of G-F // Change from NPS
  rts=realroots11(set_coeff(-(b6+b8),-two*(b4+b6),-(b2+b4), -four,one));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of F-G: "<<rts<<endl;
  cout<<"After adding roots of F-G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of G'
  rts=realroots11(set_coeff(-two*b8,-three*b6,-b4, zero));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of G': "<<rts<<endl;
  cout<<"After adding roots of G', crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of F'
  rts=realroots11(set_coeff(two*b6,three*b4,b2,two));
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
  for (const auto& x : F_roots)
    {
      x2=x*x;
      Gx=abs(1-b4*x2-two*b6*x*x2-b8*x2*x2);
#ifdef DEBUG_CPS
      cout<<"x="<<x<<", Gx="<<Gx<<endl;
#endif
      if (first)  { dvd=Gx; first=0;}
      else if (dvd>Gx) { dvd=Gx; }
#ifdef DEBUG_CPS
      cout<<"dvd so far ="<<dvd<<endl;
#endif
    }
  for (const auto& x : crit_pts)
    {
      x2=x*x;
      Fx=(four*x+b2*x2+two*b4*x*x2+b6*x2*x2);
#ifdef DEBUG_CPS
      cout<<"x="<<x<<", Fx="<<Fx<<endl;
#endif
      if(Fx<0) continue;
      Gx=abs(one-b4*x2-two*b6*x*x2-b8*x2*x2);
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
  if(first)
    return -one;
  else
    return dvd;
}

// Procedure to calculate dv for the infinite prime

bigfloat calc_dv_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8)
{
  bigfloat zero=to_bigfloat(0), one=to_bigfloat(1), two=to_bigfloat(2), six=to_bigfloat(6), four=to_bigfloat(4);
  bigfloat rx;
  bigfloat dv=zero,x2,fx,gx;

#ifdef DEBUG_CPS
  cout<<"\nIn calc_dv_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif

  vector<bigfloat> rts; // holds real roots of f, g, etc
  std::set<bigfloat> crit_pts; // use a set, so we don't get any duplicates
  std::set<bigfloat> f_roots;

  crit_pts.insert(one);
  crit_pts.insert(-one);
  
#ifdef DEBUG_CPS
  cout<<"crit_pts = "<<crit_pts<<endl;
#endif

  // Put the roots of g,f',g',f+g,f-g into crit_pts:
  // Keep the roots of f separate (see comment below)

  // Roots of f
  rts=realroots11(set_coeff(four,b2,two*b4,b6));
  f_roots.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of f: "<<rts<<endl;
  cout<<"After adding roots of f, f_roots = "<<f_roots<<endl;
#endif

  // Roots of g
  rts=realroots11(set_coeff(one, zero,-b4,-two*b6,-b8));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of g: "<<rts<<endl;
  cout<<"After adding roots of g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of f+g
  rts=realroots11(set_coeff(one, four, b2-b4,two*(b4-b6),b6-b8));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of f+g: "<<rts<<endl;
  cout<<"After adding roots of f+g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of g-f // Change from NPS
  rts=realroots11(set_coeff(one,-four,-(b2+b4),-two*(b4+b6),-(b6+b8)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of f-g: "<<rts<<endl;
  cout<<"After adding roots of f-g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of g'
  rts=realroots11(set_coeff(two,zero,-b4,-b6));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_CPS
  cout<<"relevant roots of g': "<<rts<<endl;
  cout<<"After adding roots of g', crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of f'
  rts=realroots11(set_coeff(six,b2,b4));
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
  for( const auto& x : f_roots)
    {
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
  for( const auto& x : crit_pts)
    {
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
  if(first)
    return -one;
  else
    return dv;
}

//#define HTB_DEBUG

// returns lower bound for height of non-torsion points, following
// Cremona & Siksek in ANTS7.  If egr==1, return a lower bound for the
// height of nontorsion pints with everywhere good reduction.

bigfloat lower_height_bound(const Curvedata& CD, int egr)
{
  CurveRed CR(CD);
  CurveHeightConst CHC(CR);
  CHC.compute();
  bigfloat lambda=CHC.get_value(); // assumes egr
  if (!egr)
    {
      long c = I2long(global_Tamagawa_exponent(CR, 1)); // 1 means include R
      lambda /= (c*c);
    }
  return lambda;
}

// This gives a lower bound on non-torsion points, by searching, given
// the regulator of a known subgroup

// If point search bound is greater than this, output a warning
// message and reduce to this value:
const int max_search_bound = 18;

bigfloat lower_height_bound_search(const Curvedata& CD, const bigfloat& reg)
{
  int verbose=0;
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

bigfloat factorial(long n)
{
  if(n<2)
    return to_bigfloat(1);
  if(n<13)
    return to_bigfloat(fact_tab[n]);
  return n*factorial(n-1);
}

bigfloat Gamma_n(long n) // Gamma(n) = (n-1)!
{
  return factorial(n-1);
}

bigfloat Gamma_n_plus_half(long n) // Gamma(n+1/2) = (2n)!sqrt(pi) / (4^n*n!)
{
  // cout<<"n = "<<n<<endl;
  // cout<<"factorial(n) = "<< factorial(n) <<endl;
  // cout<<"factorial(2*n) = "<< factorial(2*n) <<endl;
  // cout<<"2^(2*n) = "<< power2_RR(2*n) <<endl;
  // cout<<"sqrt(pi) = "<< sqrt(Pi()) <<endl;
  return sqrt(Pi()) * factorial(2*n) / (power2_RR(2*n) * factorial(n));
}

static long gam_tab[9] = {1, 1, 4, 2, 4, 8, 64, 64, 256};

bigfloat lattice_const(int r)
// Return gamma_r such that for a lattice of rank r and determinant D,
// the shortest nonzero vector has length at most gamma_r*D^(1/r).
//
// for   r    = 1  2  3 4 5 6     7   8
//  we have
//  gamma_r^r = 1 4/3 2 4 8 64/3 64 256
//
// while for larger r we use Blichfeld's bound
// (2/pi)*Gamma(2+r/2)^(2/r) (see Wikipedia
// https://en.wikipedia.org/wiki/Hermite_constant)
//
//
{
  if (r<=8)
    {
      bigfloat gam = to_bigfloat(gam_tab[r]);
      if (r%4==2)
        gam /= to_bigfloat(3);
      return pow(gam, inv(to_bigfloat(r)));
    }
  else
    {
      bigfloat gam = (r%2? Gamma_n_plus_half((r+3)/2): Gamma_n(2+r/2));
      return 2 * inv(Pi()) * pow(sqr(gam), inv(to_bigfloat(r)));
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
  for( const auto& I : L1)
    for( const auto& J : L2)
      {
	Interval K = intersect(I,J);
	if(!K.is_empty())
          ans.push_back(K);
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
  for( const auto& I : L1)
    for( const auto& J : L2)
      {
	Interval01 K = intersect(I,J);
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

////////////////////////////
//
// Local exponent function
//
////////////////////////////

// returns the exponent of the reduction of CD mod p (i.e. of
// E^0(Qp)/E^1(Qp), or of E^ns(F_p))
//
// NB for good reduction and p>3 we can use the curvemodqbasis class,
// but that is not implemented for p=2, 3.  We also need special code
// for bad reduction.

long exponent(CurveRed& CR, long p)
{
  bigint pp = BIGINT(p);
  int ord_p_N = getord_p_N(CR, pp);

  if (ord_p_N>1)
    // additive reduction, cyclic of order p
    {
      return p;
    }

  if (ord_p_N==1)
    // multiplicative reduction, cyclic:
    // order p-1 if split, i.e. root number -1
    // order p+1 if split, i.e. root number +1
    {
      return p + LocalRootNumber(CR, pp);
    }

  // good reduction
  if (p>3)
    {
      curvemodqbasis Emodq(CR,pp);
      return I2long(Emodq.get_exponent());
    }
  // now p=2 or 3
  int np = 1 + p - I2long(Trace_Frob(CR,pp));
  if (p==2 || np!=4)
    return np; // exponent=order
  // Now p==3, and order=4, test whether we have full 2-torsion
  // The b-invariants are (0, 2, _, 2) for C4 and (0,1,0,2) for
  // C2xC2; so looking at b4 suffices:
  return ((posmod(getb4(CR),3)==1)? 2 : 4);
}

///////////////////////////////////////////////////////////////////////
//
// Implementation of class CurveHeightConst
//
///////////////////////////////////////////////////////////////////////

CurveHeightConst::CurveHeightConst(CurveRed& CR)
  : CurveRed(CR), Cperiods(CR)
{
  c = to_bigfloat(egr_height_constant(*this)); // =-log(alpha) in ANTS7
  e3 = get_e3();
  n_max=10;
#ifdef debugLB
  cout<<"e3 = "<<e3<<endl;
  cout<<"archContrib = log(epsilon)/3 = "<<c<<endl;
  cout<<"n_max = "<<n_max<<endl;
#endif
}

long CurveHeightConst::e_p(long p)
{
  auto pe = ann.find(p);
  if (pe!=ann.end())
    return pe->second;
  long e = exponent(*this,p);
  ann[p] = e;
  return e;
}

bigfloat CurveHeightConst::D(long n) // D_E(n) in the paper
{
  auto DEn = DE.find(n);
  if (DEn!=DE.end())
    {
#ifdef debugLB
      cout << "stored D("<<n<<") = "<< DEn->second <<endl;
#endif
      return DEn->second;
    }
  // else compute and store it:
  bigfloat ans = to_bigfloat(0);
  primevar pr;
  long p, e, pmax = (n+1)*(n+1);
  for (p=pr.value(); p<pmax; pr++, p=pr.value())
    {
      e = e_p(p);
#ifdef debugLB
      //      cout << " p="<<p<<", e_p="<<e<<endl;
#endif
      if(divides(e, n))
        ans+=2*(1+val(p,n/e))*log(p);
    }
#ifdef debugLB
  cout << "D("<<n<<") = "<< ans<<endl;
#endif
  DE[n] = ans;
  return ans;
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
      cerr<<"Error in CurveHeightConst::psi(): x="<<x<<" < e3 = "<<e3<<endl;  
      return to_bigfloat(0);
    }
  //  cout<<"computing psi(x) with x = "<<x<<endl;
  //  cout<<"ordinates: "<<ordinates(x)<<endl;
  bigfloat y = ordinates(x)[0];
  //  cout<<"y = "<<y<<endl;
  bigcomplex z = pointtoz(x,y);
  //  cout<<"z = "<<z<<endl;
  return real(z/get_real_period()); // in [0.5,1]
}
