//  file esign.cc
//
// Sign of functional equation for an elliptic curve over Q
//
// Taken partly from Odile Lecacheux's GP code, partly translated by
// JEC from the pari C code in elliptic.c from pari-2.1.3
//
// For the theory, see
//
// Halberstadt, Emmanuel. Signes locaux des courbes elliptiques en 2 et 3.
// (French) [Local root numbers of elliptic curves for $p=2$ or $3$]
// C. R. Acad. Sci. Paris Sér. I Math. 326 (1998), no. 9, 1047--1052.
//

#include "curve.h"

//#define DEBUG_ESIGN 1

int kro(const bigint& d, const bigint& n)
{
#ifdef DEBUG_ESIGN
  cout<<"kro("<<d<<","<<n<<") returns "<<flush;
#endif
  int ans = kronecker(d,n);
#ifdef DEBUG_ESIGN
  cout<<ans<<endl;
#endif
  return ans;
}

int kro(const bigint& d, long n)
{
#ifdef DEBUG_ESIGN
  cout<<"kro("<<d<<","<<n<<") returns "<<flush;
#endif
  int ans = kronecker(d,n);
#ifdef DEBUG_ESIGN
  cout<<ans<<endl;
#endif
  return ans;
}

int kro(long d, long n)
{
#ifdef DEBUG_ESIGN
  cout<<"kro("<<d<<","<<n<<") returns "<<flush;
#endif
  int ans = kronecker(d,n);
#ifdef DEBUG_ESIGN
  cout<<ans<<endl;
#endif
  return ans;
}

int kro_m1(long x) // kronecker(-1,x) with x>0 odd
{
  static int kro_m1_tab[4] = {0,1,0,-1};
#ifdef DEBUG_ESIGN
  if (!((x>0)&&(odd(x))))  
    { 
      cout<<"kro_m1() called with x="<<x<<endl;
      return 0;
    }
#endif
  return kro_m1_tab[x&3];
}

int kro_p2(long x) // kronecker(2,x) with x>0 odd
{
  static int kro_p2_tab[8] = {0,1,0,-1,0,-1,0,1};
#ifdef DEBUG_ESIGN
  if (!((x>0)&&(odd(x))))  
    {
      cout<<"kro_p2() called with x="<<x<<endl;
      return 0;
    }
#endif
  return kro_p2_tab[x&7];
}

int kro_m2(long x) // kronecker(-2,x) with x>0 odd
{
  static int kro_m2_tab[8] = {0,1,0,1,0,-1,0,-1};
#ifdef DEBUG_ESIGN
  if (!((x>0)&&(odd(x))))  
    {
      cout<<"kro_m2() called with x="<<x<<endl;
      return 0;
    }
#endif
  return kro_m2_tab[x&7];
}

int kro_3(long x) // kronecker(x,3)
{
  static int kro_3_tab[3] = {0,1,-1};
#ifdef DEBUG_ESIGN
  if (!(x>0))  
    {
      cout<<"kro_3() called with x="<<x<<endl;
      return 0;
    }
#endif
  return kro_3_tab[x%3];
}

// utility function for converting Kodaira codes to the Pari coding

// Kodaira Symbol        My coding    Pari Coding

// I0                    0              1
// I*0                   1             -1
// Im  (m>0)             10*m           m+4
// I*m (m>0)             10*m+1        -(m+4)
// II, III, IV           2, 3, 4        m
// II*. III*, IV*        7, 6, 5       -m


int PariKodairaCode(Kodaira_code Kod)
{
  int ans, code=Kod.code;
  if (code==0) ans = 1;
  else if (code==1) ans = -1;
  else if (code%10 ==0) ans = (code/10)+4;
  else if (code%10 ==1) ans = ((1-code)/10)-4;
  else if (code>4) ans = code-9;
  else ans=code;
#ifdef DEBUG_ESIGN
  cout<<"PariKodairaCode("<<Kod<<") returns "<<ans<<endl;
#endif
  return ans;
}

// LOCAL ROOT NUMBERS, D'APRES HALBERSTADT halberst@math.jussieu.fr 

// translated by JEC from the pari C code in elliptic.c from pari-2.1.3

//  p = 2 or 3 for the neron function

int neron(const CurveRed& E, long p, int kod)
{
  bigint c4,c6; E.getci(c4,c6);
  bigint d=getdiscr(E);
  int v4=val(p,c4);
  int v6=val(p,c6);
  int vd=val(p,d);
#ifdef DEBUG_ESIGN
  cout<<"In neron with p="<<p<<", v4="<<v4<<", v6="<<v6<<", vd="<<vd<<", kod="<<kod<<endl;
#endif
  if (p==3)
  {
    if (abs(kod)>4) return 1;
    else
      {
	switch(kod)
	  {
	  case -1: case 1: return v4&1 ? 2 : 1;
	  case -3: case 3: return (2*v6>vd+3) ? 2 : 1;
	  case -4: case 2:
	    switch (vd%6)
	      {
	      case 4: return 3;
	      case 5: return 4;
	      default: return v6%3==1 ? 2 : 1;
	      }
	  default: /* kod = -2 et 4 */
	    switch (vd%6)
	      {
	      case 0: return 2;
	      case 1: return 3;
	      default: return 1;
	      }
	  }
      }
  }
  if(p==2)
    {
      if (kod>4) return 1;
      else
	{
	  switch(kod)
	    {
	    case 1: return (v6>0) ? 2 : 1;
	    case 2:
	      if (vd==4) return 1;
	      else
		{
		  if (vd==7) return 3;
		  else return v4==4 ? 2 : 4;
		}
	    case 3:
	      switch(vd)
		{
		case 6: return 3;
		case 8: return 4;
		case 9: return 5;
		default: return v4==5 ? 2 : 1;
		}
	    case 4: return v4>4 ? 2 : 1;
	    case -1:
	      switch(vd)
		{
		case 9: return 2;
		case 10: return 4;
		default: return v4>4 ? 3 : 1;
		}
	    case -2:
	      switch(vd)
		{
		case 12: return 2;
		case 14: return 3;
		default: return 1;
		}
	    case -3:
	      switch(vd)
		{
		case 12: return 2;
		case 14: return 3;
		case 15: return 4;
		default: return 1;
		}
	    case -4: return v6==7 ? 2 : 1;
	    case -5: return (v6==7 || v4==6) ? 2 : 1;
	    case -6:
	      switch(vd)
		{
		case 12: return 2;
		case 13: return 3;
		default: return v4==6 ? 2 : 1;
		}
	    case -7: return (vd==12 || v4==6) ? 2 : 1;
	    default:
	      return v4==6 ? 2 : 1;
	    }
	}
    }
  cerr<<"neron() returns 0 -- should not happen!"<<endl;
  return 0; /* should not occur */
}

// Given E, an elliptic curve over Q,
// returns +1 or -1, the "local root number" or local factor at 2
// in the sign of the functional  equation of L(E,s).

int LocalRootNumber2(const CurveRed& E)
{
  static const bigint two = BIGINT(2);
  int kod=PariKodairaCode(getKodaira_code(E,two));
  int n2=neron(E,2,kod); 
  bigint a1,a2,a3,a4,a6,c4,c6;
  E.getai(a1,a2,a3,a4,a6);
  E.getci(c4,c6);

#ifdef DEBUG_ESIGN
  cout<<"\nIn LocalRootNumber2(), n2 = "<<n2<<"..."<<endl;
#endif

  bigint mu,mv; long u,v; int v4,v6;

  if (is_zero(c4)) {v4=12; u=0;}
  else {mu=c4; v4=divide_out(mu,two); u = posmod(mu,64);}
#ifdef DEBUG_ESIGN
  cout<<"c4="<<c4<<", v4="<<v4<<", u="<<u<<endl;
#endif

  if (is_zero(c6)) {v6=12; v=0;}
  else {mv=c6; v6=divide_out(mv,two); v = posmod(mv,64);}
#ifdef DEBUG_ESIGN
  cout<<"c6="<<c6<<", v6="<<v6<<", v="<<v<<endl;
#endif

  if (kod > 4)  return div(2,a2+a3)? -1: 1;
  if (kod < -9) return (n2==2)? -kro_m1(v) : -1;

  bigint tmp = getdiscr(E);
  divide_out(tmp,two);
  long d1=posmod(tmp,64);

  long x1=u+v+v, y1;

  switch(kod)
  {
    case 1: return 1;
    case 2:
      switch(n2)
      {
	case 1:
	  switch(v4)
	  {
	    case 4: return kro_m1(u);
	    case 5: return 1;
	    default: return -1;
	  }
	case 2: return (v6==7) ? 1 : -1;
	case 3: return (v%8==5 || (u*v)%8==5) ? 1 : -1;
	case 4: if (v4>5) return kro_m1(v);
	  return (v4==5) ? -kro_m1(u) : -1;
      }
    case 3:
      switch(n2)
      {
	case 1: return -kro_p2(u*v);
	case 2: return -kro_p2(v);
        case 3: y1=posmod((u-(c6 >> 5)) , 16);
	  return (y1==7 || y1==11) ? 1 : -1;
	case 4: return (v%8==3 || (2*u+v)%8==7) ? 1 : -1;
	case 5: return v6==8 ? kro_p2(x1) : kro_m2(u);
      }
    case -1:
      switch(n2)
      {
	case 1: return -kro_p2(x1);
	case 2: return (v%8==7) || (x1%32==11) ? 1 : -1;
	case 3: return v4==6 ? 1 : -1;
	case 4: if (v4>6) return kro_m1(v);
	  return v4==6 ? -kro_m1(u*v) : -1;
      }
    case -2: return n2==1 ? kro_m2(v) : kro_m1(v);
    case -3:
      switch(n2)
      {
	case 1: y1=posmod((u-2*v),64);
	  return (y1==3) || (y1==19) ? 1 : -1;
      case 2: //return kro(2*kro_m1(u),v);
	if(kro_m1(u)==1) return kro_p2(v); else return kro_m2(v);
      case 3: // return -kro_m1(u)*kro(-2*kro_m1(u),u*v);
	if(kro_m1(u)==1) return -kro_m2(u*v); else return kro_p2(u*v);
	case 4: return v6==11 ? kro_m2(x1) : -kro_m2(u);
      }
    case -5:
      if (n2==1) return x1%32==23 ? 1 : -1;
      else return -kro_p2(2*u+v);
    case -6:
      switch(n2)
      {
	case 1: return 1;
	case 2: return v6==10 ? 1 : -1;
	case 3: return (u%16==11) || ((u+4*v)%16==3) ? 1 : -1;
      }
    case -7:
      if (n2==1) return 1;
      else
      {
        y1= posmod((u+(c6 >> 8)) , 16);
	if (v6==10) return (y1==9) || (y1==13) ? 1 : -1;
	else return (y1==9) || (y1==5) ? 1 : -1;
      }
    case -8: return n2==2 ? kro_m1(v*d1) : -1;
    case -9: return n2==2 ? -kro_m1(d1) : -1;
    default: return -1;
  }
}

// Given E, an elliptic curve over Q,
// returns +1 or -1, the "local root number" or local factor at 3 
// in the sign of the functional  equation of L(E,s).

int LocalRootNumber3(const CurveRed& E)
{
  static const bigint three = BIGINT(3);
  int kod=PariKodairaCode(getKodaira_code(E,three));
  int n2=neron(E,3,kod); 
  bigint a1,a2,a3,a4,a6,c4,c6;
  E.getai(a1,a2,a3,a4,a6);
  E.getci(c4,c6);

#ifdef DEBUG_ESIGN
  cout<<"\nIn LocalRootNumber3()..."<<endl;
#endif
  bigint mu,mv; long u,v; int v4,v6;

  if (is_zero(c4)) { v4=12; u=0;}
  else {mu=c4; v4=divide_out(mu,three); u = posmod(mu,81);}

#ifdef DEBUG_ESIGN
  cout<<"c4="<<c4<<", v4="<<v4<<", u="<<u<<endl;
#endif

  if (is_zero(c6)) { v=0;}
  else {mv=c6; v6=divide_out(mv,three); v = posmod(mv,81);}

#ifdef DEBUG_ESIGN
  cout<<"c6="<<c6<<", v6="<<v6<<", v="<<v<<endl;
#endif

  bigint tmp = getdiscr(E);
  divide_out(tmp,three);
  long d1=posmod(tmp,81);

#ifdef DEBUG_ESIGN
  cout<<"d1="<<d1<<endl;
#endif

  long r6 = posmod(v,9);
  long K4=kro_3(u), K6=kro_3(v);
#ifdef DEBUG_ESIGN
  cout<<"r6="<<r6<<endl;
  cout<<"K4="<<K4<<endl;
  cout<<"K6="<<K6<<endl;
#endif

  if (kod > 4) { return K6;}
  
  switch(kod)
    {
    case 1: case 3: case -3: return 1;
    case 2:
      switch(n2)
	{
	case 1: return (r6==4 || r6>6) ? 1 : -1;
	case 2: return -K4*K6;
	case 3: return 1;
	case 4: return -K6;
	}
    case 4:
      switch(n2)
	{
	case 1: return K6*kro_3(d1);
	case 2: return -K4;
	case 3: return -K6;
	}
    case -2: return n2==2 ? 1 : K6;
    case -4:
      switch(n2)
	{
	case 1:
	  if (v4==4) return (r6==4 || r6==8) ? 1 : -1;
	  else return (r6==1 || r6==2) ? 1 : -1;
	case 2: return -K6;
	case 3: return (r6==2 || r6==7) ? 1 : -1;
	case 4: return K6;
	}
    default: return -1;
    }
}  

// Given E, an elliptic curve over Q, and a prime p not 2 or 3, 
// returns +1 or -1, the "local root number" or local factor in the 
// sign of the functional  equation of L(E,s).

int LocalRootNumber_not_2_or_3(const CurveRed& E, const bigint& p)
{
  if (getord_p_N(E,p) == 1) 
    {
      bigint c4,c6;
      E.getci(c4,c6);
      return -kro(-c6,p);
    }
  
  long sp=posmod(p,24);
  if (getord_p_j_denom(E,p) >0)  return kro_m1(sp); 
  
  long ep=12 / gcd(12,getord_p_discr(E,p));
  if(ep==4) return kro_m2(sp);
  if(odd(ep)) return kro_3(sp); 
  return kro_m1(sp); 
  //   long z= (ep==4)? 2 : ((ep%2)? 3 : 1);
  //  return kro(BIGINT(-z),p);
}

// Given E, an elliptic curve over Q, and a prime p, returns +1 or -1, 
// the "local root number" or local factor in the sign of the functional 
// equation of L(E,s).
//
// For p=0 it gives the factor at the infinte prime, which is always -1
//
// This function just delegates to subsidiary ones for the cases 
// p=2, p=3, and p>=5.
//
int LocalRootNumber(const CurveRed& E, const bigint& p)
{
  if (is_zero(p)) return -1;  // local factor at infinite prime
  if (p==2)       return LocalRootNumber2(E);
  if (p==3)       return LocalRootNumber3(E);
  return LocalRootNumber_not_2_or_3(E,p);  
}

// Given E, an elliptic curve over Q, returns +1 or -1, being the
// "global root number" or sign of the functional equation of L(E,s).

int GlobalRootNumber(const Curvedata e)
{
  CurveRed E(e);
  int ans=-1;
  vector<bigint> bad_primes = getbad_primes(E);
  vector<bigint>::const_iterator pi= bad_primes.begin();
  while(pi!=bad_primes.end())  
    {
      bigint p=*pi++;
      int s = LocalRootNumber(E,p); 
#ifdef DEBUG_ESIGN
      cout<<"Sign at p="<<p<<" is "<<s<<endl;
#endif
      ans *= s;
    }
  return ans;
}
