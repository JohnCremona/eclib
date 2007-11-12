// silver.cc: implementations of (1) Silverman, (2) Siksek bounds on
// the difference between naive and canonical height.  Curently (2) is
// under development

#include "curve.h"
#include "compproc.h"
#include "silver.h"

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
  
  if(b2!=0) mu += log(2.0);

  return mu;
}


// Siksek's height bound, August 2002
// NB: We assume a minimal model here!

//#define DEBUG_SIKSEK
//#define TEST_SIKSEK

double siksek_real(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8);

double siksek_bound(const Curvedata& CD)
{
  bigfloat b2=I2bigfloat(getb2(CD));
  bigfloat b4=I2bigfloat(getb4(CD));
  bigfloat b6=I2bigfloat(getb6(CD));
  bigfloat b8=I2bigfloat(getb8(CD));
  double bd = siksek_real(I2bigfloat(getb2(CD)),I2bigfloat(getb4(CD)),
			  I2bigfloat(getb6(CD)),I2bigfloat(getb8(CD)));
#ifdef DEBUG_SIKSEK
  cout<<"In siksek_bound() for "<<(Curve)CD<<endl;
  cout<<"siksek_real = "<<bd<<endl;
#endif
  CurveRed CR(CD);
  vector<bigint> plist = getbad_primes(CD);
  for(int i=0; i<plist.size(); i++)
    {
      bigint q=plist[i];
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
#ifdef DEBUG_SIKSEK
      cout<<"q = "<<q<<", alpha = "<<alpha<<", q-term = "<< alpha*log(double(I2long(q))) <<endl;
      cout<<"sum so far = "<<bd<<endl;
#endif
    }
  return bd;
}

// Implementation for the real place originally by Nigel Smart, here
// rewritten by JC, 22/8/02

bigfloat calc_dvd_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8);
bigfloat calc_dv_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8);
bigfloat old_calc_dvd_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8, const bigfloat& del);
bigfloat old_calc_dv_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8, const bigfloat& del);

double siksek_real(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8)
{
  bigfloat dv=calc_dv_inf(b2,b4,b6,b8);     
  bigfloat dvd=calc_dvd_inf(b2,b4,b6,b8); 
#ifdef TEST_SIKSEK
  bigfloat del = -b2*b2*b8-8*b4*b4*b4-27*b6*b6+9*b2*b4*b6;
  bigfloat dv2=old_calc_dv_inf(b2,b4,b6,b8,del);
  bigfloat dvd2=old_calc_dvd_inf(b2,b4,b6,b8,del); 
  if(!is_zero(dv-dv2)) cout<<"old_calc_dv gives "<<dv2<<"\nwhile new gives "<<dv<<endl;
  if(!is_zero(dvd-dvd2)) cout<<"old_calc_dvd gives "<<dvd2<<"\nwhile new gives "<<dvd<<endl;
#endif
#ifdef DEBUG_SIKSEK
  cout << "dv=" << dv << endl;
  cout << "dvd=" << dvd << endl;
#endif
  bigfloat htc;
  if(dv==-1)
    {
      if(dvd==-1) htc = 0;
      else htc = -log(dvd)/3;
    }
  else
    if(dvd==-1) -log(dv)/3;
    else htc = -log(min(dv,dvd))/3;

#ifdef DEBUG_SIKSEK
  cout<<"siksek_real() returns -log(min(dv,dvd))/3 = "<<htc<<endl;
#endif

#ifdef LiDIA_ALL
  double ans;
  htc.doublify(ans);
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

inline vector<bigfloat> set_reverse_coeff(const bigfloat& c0, const bigfloat& c1, const bigfloat& c2, const bigfloat& c3, const bigfloat& c4)
{
  return set_coeff(c4,c3,c2,c1,c0);
}

template <class T >
inline ostream& operator<<(ostream& os, const std::set<T>& v);

vector<bigfloat> reals_in ( vector<bigcomplex>& v);
vector<bigfloat> reals_in_11 ( vector<bigcomplex>& v);

//  Procedure to calculate dv' for the infinite prime

bigfloat calc_dvd_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8)
{ 
  bigfloat rr,rx,rn;
  bigfloat d,dvd,x,x2,Fx,Gx;

#ifdef DEBUG_SIKSEK
  cout<<"\nIn calc_dvd_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif

  vector<bigfloat> rts; // holds real roots of F, G, etc
  std::set<bigfloat> crit_pts; // use a set, so we don't get any duplicates

  crit_pts.insert(1);
  crit_pts.insert(-1);
  
#ifdef DEBUG_SIKSEK
  cout<<"crit_pts = "<<crit_pts<<endl;
#endif

  // Put the roots of F,G,F',G',F+G,F-G into crit_pts:

  // Roots of G
  rts=roots11(set_reverse_coeff(1,0,-b4,-2*b6,-b8));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of G: "<<rts<<endl;
  cout<<"After adding roots of G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of F
  rts=roots11(set_reverse_coeff(0,1,b2/4.0,b4/2.0,b6/4.0));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of F: "<<rts<<endl;
  cout<<"After adding roots of F, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of F+G
  rts=roots11(set_reverse_coeff(1,4,b2-b4,2.0*(b4-b6),b6-b8));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of F+G: "<<rts<<endl;
  cout<<"After adding roots of F+G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of G-F // Change from NPS
  rts=roots11(set_reverse_coeff(1,-4,-(b2+b4),-2.0*(b4+b6),-(b6+b8)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of F-G: "<<rts<<endl;
  cout<<"After adding roots of F-G, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of G'
  rts=roots11(set_coeff(0,-2*b8,-3*b6,-b4,0));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of G': "<<rts<<endl;
  cout<<"After adding roots of G', crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of F'
  rts=roots11(set_coeff(0,2*b6,3*b4,b2,2));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of F': "<<rts<<endl;
  cout<<"After adding roots of F', crit_pts = "<<crit_pts<<endl;
#endif

  // Evaluate max(|F(x)|,|G(x)|) at each of the x-values in array crit_pts for which F(x)>=0:
  std::set<bigfloat>::const_iterator xi = crit_pts.begin();
  int first=1;
  while(xi!=crit_pts.end())
    { 
      x=*xi++;
      x2=x*x;
      Fx=abs(4.0*x+b2*x2+2.0*b4*x*x2+b6*x2*x2);
#ifdef DEBUG_SIKSEK
      cout<<"x="<<x<<", Fx="<<Fx<<endl;
#endif
      if(Fx<0) continue;
      Gx=abs(1-b4*x2-2.0*b6*x*x2-b8*x2*x2);
#ifdef DEBUG_SIKSEK
      cout<<"x="<<x<<", Gx="<<Gx<<endl;
#endif
      rx=max(Fx,Gx);
      if (first)  { dvd=rx; first=0;}
      else if (dvd>rx) { dvd=rx; }
#ifdef DEBUG_SIKSEK
      cout<<"dvd so far ="<<dvd<<endl;
#endif
    }
  if(first) return -1;

  return dvd;

}

// Procedure to calculate dv for the infinite prime

bigfloat calc_dv_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8)
{
  bigfloat rr,rx,rn;
  bigfloat d,dv,x,x2,fx,gx;

#ifdef DEBUG_SIKSEK
  cout<<"\nIn calc_dv_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif

  vector<bigfloat> rts; // holds real roots of f, g, etc
  std::set<bigfloat> crit_pts; // use a set, so we don't get any duplicates

  crit_pts.insert(1);
  crit_pts.insert(-1);
  
#ifdef DEBUG_SIKSEK
  cout<<"crit_pts = "<<crit_pts<<endl;
#endif

  // Put the roots of f,g,f',g',f+g,f-g into crit_pts:

  // Roots of g
  rts=roots11(set_coeff(1,0,-b4,-2*b6,-b8));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of g: "<<rts<<endl;
  cout<<"After adding roots of g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of f
  rts=roots11(set_coeff(0,1,b2/4.0,b4/2.0,b6/4.0));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of f: "<<rts<<endl;
  cout<<"After adding roots of f, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of f+g
  rts=roots11(set_coeff(1,4,b2-b4,2.0*(b4-b6),b6-b8));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of f+g: "<<rts<<endl;
  cout<<"After adding roots of f+g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of g-f // Change from NPS
  rts=roots11(set_coeff(1,-4,-(b2+b4),-2.0*(b4+b6),-(b6+b8)));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of f-g: "<<rts<<endl;
  cout<<"After adding roots of f-g, crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of g'
  rts=roots11(set_coeff(0,1,0,-b4/2.0,-b6/2.0));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of g': "<<rts<<endl;
  cout<<"After adding roots of g', crit_pts = "<<crit_pts<<endl;
#endif

  // Roots of f'
  rts=roots11(set_coeff(0,0,1,b2/6,b4/6));
  crit_pts.insert(rts.begin(),rts.end());
#ifdef DEBUG_SIKSEK
  cout<<"relevant roots of f': "<<rts<<endl;
  cout<<"After adding roots of f', crit_pts = "<<crit_pts<<endl;
#endif

  // Evaluate max(|f(x)|,|g(x)|) at each of the x-values in array crit_pts for which f(x)>=0:
  std::set<bigfloat>::const_iterator xi = crit_pts.begin();
  int first=1;
  while(xi!=crit_pts.end())
    { 
      x=*xi++;
      x2=x*x;
      fx=abs(4.0*x*x2+b2*x2+2.0*b4*x+b6);
#ifdef DEBUG_SIKSEK
      cout<<"x="<<x<<", fx="<<fx<<endl;
#endif
      if(fx<0) continue;
      gx=abs(x2*x2-b4*x2-2.0*b6*x-b8);
#ifdef DEBUG_SIKSEK
      cout<<"x="<<x<<", gx="<<gx<<endl;
#endif
      rx=max(fx,gx);
      if (first)  { dv=rx; first=0;}
      else if (dv>rx) { dv=rx; }
#ifdef DEBUG_SIKSEK
      cout<<"dv so far ="<<dv<<endl;
#endif
    }
  if(first) return -1;
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
  bigfloat rr,rx,rn;

#ifdef DEBUG_SIKSEK
  cout<<"\nIn old_calc_dv_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif
  vector<bigcomplex> rt = solvecubic(b2/4,b4/2,b6/4);

  bigfloat i1l,i1u,i2l,i2u;  // Bounds on Intervals
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

  if(numint==0) return -1;  // code for "infinity"

  bigfloat* te = new bigfloat[36];
  long tec=2*numint;
  if (numint>0)
    { te[0]=i1l; te[1]=i1u;
      if (numint==2)
        { te[2]=i2l; te[3]=i2u; }
    }

  if (numint!=0)
    { // Roots of g
      rt=solvequartic(bigfloat(0),-b4,-2.0*b6,-b8);
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
      rt=solvequartic(bigfloat(4),b2-b4,2.0*(b4-b6),b6-b8);
      for (i=0; i<4; i++)
        { if (is_approx_zero(rt[i].imag()))
            { rr=rt[i].real();
              if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
                { te[tec]=rr; tec=tec+1; }
            }
        }

      // Roots of g-f // Change from NPS
      rt=solvequartic(bigfloat(-4),-(b2+b4),-2.0*(b4+b6),-(b6+b8));
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
      rt=solvecubic(bigfloat(0),-b4/2.0,-b6/2.0);
      for (i=0; i<3; i++)
        { if (is_approx_zero(rt[i].imag()))
            { rr=rt[i].real();
              if (is_in_int(rr,i1l,i1u,i2l,i2u,numint))
                { te[tec]=rr; tec=tec+1; }
            }
        }

    }

  bigfloat dv;
  // Evaluate |f|, |g| at each of the tec(>0) x-values in array te:
  for (i=0; i<tec; i++)
    { rr=te[i]*te[i];
      rn=4.0*rr*te[i]+b2*rr+2.0*b4*te[i]+b6;
#ifdef DEBUG_SIKSEK
      cout<<"x="<<te[i]<<", fx="<<rn<<endl;
#endif
      rn=abs(rn);
      rx=rr*rr-b4*rr-2.0*b6*te[i]-b8;
      rx=abs(rx);
#ifdef DEBUG_SIKSEK
      cout<<"x="<<te[i]<<", gx="<<rx<<endl;
#endif
      if (rn>rx) { rx=rn; }
      if (i==0)  { dv=rx; }
      else if (dv>rx) { dv=rx; }
#ifdef DEBUG_SIKSEK
      cout<<"dv so far ="<<dv<<endl;
#endif
    }
  delete[] te;

  return dv;
}

bigfloat old_calc_dvd_inf(const bigfloat& b2, const bigfloat& b4, const bigfloat& b6, const bigfloat& b8, const bigfloat& del)
{ 
  bigfloat rr,rx,rn;

#ifdef DEBUG_SIKSEK
  cout<<"\nIn old_calc_dvd_inf"<<endl;
  cout<<"b2="<<b2<<"\nb4="<<b4<<"\nb6="<<b6<<"\nb8="<<b8<<endl;
#endif
  vector<bigcomplex> rt = solvecubic(b2/4,b4/2,b6/4);
  long i,j,numrt;
  bigfloat rrt[4];
  rrt[0]=0.0;
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

  bigfloat rhs,lhs;
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
      rt=solvequartic(bigfloat(0),-b4,-2.0*b6,-b8);
      for (i=0; i<4; i++)
	{ if (!is_approx_zero(rt[i]))
	    { rt[i]=1.0/rt[i];
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
            { rt[i]=1.0/rt[i];
              if (is_approx_zero(rt[i].imag()))
                 { rr=rt[i].real();
                   if (is_in_int2(rr,intervals,numint))
                     { te[tec]=rr; tec=tec+1; }
                 }
            }
        }

      // Roots of f+g
      rt=solvequartic(bigfloat(4),b2-b4,2.0*(b4-b6),b6-b8);
      for (i=0; i<4; i++)
        { if (!is_approx_zero(rt[i]))
            { rt[i]=1.0/rt[i];
              if (is_approx_zero(rt[i].imag()))
                 { rr=rt[i].real();
                   if (is_in_int2(rr,intervals,numint))
                     { te[tec]=rr; tec=tec+1; }
		 }
            }
        }

      // Roots of f-g
      rt=solvequartic(bigfloat(-4),-b2-b4,2.0*(-b4-b6),-b6-b8);
      for (i=0; i<4; i++)
        { if (!is_approx_zero(rt[i]))
            { rt[i]=1.0/rt[i];
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

  bigfloat dvd;
  bigfloat f,g;
  for (i=0; i<tec; i++)
    { rr=te[i];
      f=(((b6*rr+2.0*b4)*rr+b2)*rr+4.0)*rr;
      g=(((-b8*rr-2.0*b6)*rr-b4)*rr)*rr+1.0;
      f=abs(f); g=abs(g);
      rn=f;
      if (g>f) { rn=g; }
      if (i==0)  { dvd=rn; }
      else if (dvd>rn) { dvd=rn; }
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
  bigfloat a=coeff[0], a1;
  if(a!=0) 
    {
      a1=1/a;
      vector<bigcomplex> cr = solvequartic(a1*coeff[1],a1*coeff[2],a1*coeff[3],a1*coeff[4]);
      return reals_in_11(cr);
    }
  a=coeff[1];
  if(a!=0) 
    {
      a1=1/a;
      vector<bigcomplex> cr = solvecubic(a1*coeff[2],a1*coeff[3],a1*coeff[4]);
      return reals_in_11(cr);
    }
  vector<bigfloat> ans; // zero length
  a=coeff[2];
  if(a==0) cerr<<"Error in roots: degree<2"<<endl;
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

template <class T >
inline ostream& operator<<(ostream& os, const std::set<T>& v)
{
  os <<"[";
  copy(v.begin(),v.end(), ostream_iterator<T>(cout, " "));
  os << "]";
  return os;
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
