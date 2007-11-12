// mglobsol.cc: implementation of class quartic_sieve and functions for quartic solubility testing
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
 
// (1) simple search with function ratpoint()
// (2) sieve-assisted search using class quartic_sieve

#include "mglobsol.h"

int test(const bigint& x, const bigint& z, const bigint& y2, bigint&xx, bigint&yy, bigint&zz)
{
  bigint y;
  int ans = isqrt(y2,y);
  if (ans) { xx=x; yy=y; zz=z; }
  return ans;
}

int ratpoint(const quartic& g, const bigint& min, const bigint& max, bigint& xx, bigint&yy, bigint& zz)
{ bigint a=g.geta(), b=g.getb(), c=g.getcc(), d=g.getd(), e=g.gete();
  bigint n,sx,gg; int found = 0;
  // Check for rational points at 0 or infinity:
  if(isqrt(a,n)) {xx=1; yy=n; zz=0; return 1;}
  if(isqrt(e,n)) {xx=0; yy=n; zz=1; return 1;}
  bigint x,x2,x3,x4,z,z2,z3,z4,ax4,bx3z,cx2z2,dxz3,ez4;
  static const bigint zero = BIGINT(0);
  static const bigint one  = BIGINT(1);
  for (n=min; (n<=max) && (!found); ++n)
    { 
      if (n==1) 
	{ 
	  found=test( one, zero, a, xx, yy, zz);  // pt at infty iff a square
	  if (!found) 
	    found=test( zero, one, e, xx, yy, zz); // 0 iff e square
	}
      else
	{
	  for (sx=1; (sx<n) && (! found); ++sx)
	    { 
	      gg=gcd(sx,n);
	      if (gg==1)
		{
		  x=sx;
		  x2=x*x; x3=x*x2; x4=x*x3; z=n-x; z2=z*z; z3=z*z2; z4=z*z3;
		  ax4=a*x4; bx3z=b*x3*z; cx2z2=c*x2*z2; dxz3=d*x*z3; ez4=e*z4;
		  found = test(x,z,ax4+bx3z+cx2z2+dxz3+ez4,xx,yy,zz);
		  if (!found) 
		    found = test(-x,z,ax4-bx3z+cx2z2-dxz3+ez4,xx,yy,zz);
		}
	    }
	}
    }
  return found;
}           /* end of ratpoint */
 

quartic_sieve::quartic_sieve(quartic * gg, int moduli_option, int verb)
: g(gg), verbose(verb) 
{
  a=g->geta();
  b=g->getb();
  c=g->getcc();
  d=g->getd();
  e=g->gete();
  nwprimes= 0;

  easy=0;
  if(isqrt(a,roota)) easy+=1;
  if(isqrt(e,roote)) easy+=2;
  
  //  if(easy&1) cout<<"Easy case: "<<a<<" is square, root = "<<roota<<endl;
  //  if(easy&2) cout<<"Easy case: "<<e<<" is square, root = "<<roote<<endl;
  
  use_stoll = !moduli_option;
  if(use_stoll)  // using Stoll search so rest not needed.
    return; 
  
// set up list of primes p which cannot divide w (=denom(x)) since (a/p)=-1
// (idea of J. Gebel)
// similar list of primes p which cannot divide u (=numer(x)) since (e/p)=-1

  if(!easy) 
    {
      nwprimes=25;
      wprimes = new long[nwprimes];
      long nwp=0;
      long a8=mod(a,8); long b8=mod(2*b,8); long c8=mod(4*c,8);
      long e8=mod(e,8); long d8=mod(2*d,8);
      long z0, t, tt, x, p;
      int two_is_ok = 0;
      for(z0=0; (z0<4) && !two_is_ok; z0++)
	{
	  t = (a8+c8*z0*z0)%8;
	  for(x=1; (x<8) && !two_is_ok; x+=2)
	    {
	      tt = mod((t + b8*x*z0),8);
	      if((tt==0)||(tt==1)||(tt==4)) two_is_ok=1; 
	    }
	}
      if(!two_is_ok) wprimes[nwp++]=2;
      primevar pr; pr++; // to start at 3
      for(;nwp<nwprimes; pr++)
	{
	  p=pr;
	  if(legendre(a,p)==-1) wprimes[nwp++]=p;
	}
      if(verbose) 
	{
	  cout<<"w-primes: ";
	  for(nwp=0; nwp<nwprimes; nwp++) cout<<wprimes[nwp]<<" ";
	  cout<<endl;
	}
// repeat for u-primes:
      uprimes = new long[nwprimes];
      nwp=0;
      two_is_ok = 0;
      for(z0=0; (z0<4) && !two_is_ok; z0++)
	{
	  t = (e8+c8*z0*z0)%8;
	  for(x=1; (x<8) && !two_is_ok; x+=2)
	    {
	      tt = (t + d8*x*z0)%8;
	      if((tt==0)||(tt==1)||(tt==4)) two_is_ok=1; 
	    }
	}
      if(!two_is_ok) uprimes[nwp++]=2;
      pr.init(); pr++; // to start at 3
      for(;nwp<nwprimes; pr++)
	{
	  p=pr;
	  if(legendre(e,p)==-1) uprimes[nwp++]=p;
	}
      if(verbose) 
	{
	  cout<<"u-primes: ";
	  for(nwp=0; nwp<nwprimes; nwp++) cout<<uprimes[nwp]<<" ";
	  cout<<endl;
	}
    }

// set up list of auxiliary moduli
// and a list of which residues modulo each of these are squares

  switch(moduli_option) {
  case 1:
    num_aux = 10; 
    auxs = new long[num_aux];
    auxs[0]=3;
    auxs[1]=5;
    auxs[2]=7;
    auxs[3]=11;
    auxs[4]=13;
    auxs[5]=17;
    auxs[6]=19;
    auxs[7]=23;
    auxs[8]=29;
    auxs[9]=31;
    break;

  case 2:// the following taken from Gebel's scheme
    num_aux = 3; 
    auxs = new long[num_aux];
    auxs[0]=5184;  // = (2^6)*(3^4)   // old: 6624; //  = (2^5)*(3^2)*23
    auxs[1]=5929;  // = (7^2)*(11^2)  // old: 8075; //  = (5^2)*17*19
    auxs[2]=4225;  // = (5^2)*(13^2)  // old: 7007; //  = (7^2)*11*13
    break;

  case 3:
  default:
    num_aux = 9;
    auxs = new long[num_aux];
    auxs[0]=32;
    auxs[1]= 9;
    auxs[2]=25;
    auxs[3]=49;
    auxs[4]=11;
    auxs[5]=13;
    auxs[6]=17;
    auxs[7]=19;
    auxs[8]=23;
    break;
  }
  
  xgood_mod_aux = new int*[num_aux];
//  x1good_mod_aux = new int*[num_aux];
  squares = new int*[num_aux];
  umod = new long[num_aux];

  long i,j;
  for (i = 0; i < num_aux; i++)
    {
      long aux = auxs[i];
      long half_aux = ((aux + 1) / 2);
      squares[i] = new int[aux];
      for (j = 0; j < aux; j++)      squares[i][j]=0;
      for (j = 0; j < half_aux; j++) squares[i][posmod( j*j, aux )]=1;
      xgood_mod_aux[i] = new int[aux];
    }  // end of aux loop
  
  if(verbose>1) 
    {
      cout << "Finished constructing quartic_sieve, using ";
      switch(moduli_option)
	{
	case 1: cout << "ten primes 3..31"; break;
	      case 2: cout << "three composite moduli"; break;
	      case 3: cout << "prime powers"; break;
	      }
      cout << endl;
    }
}

quartic_sieve::~quartic_sieve()
{
  if(nwprimes) { delete[] wprimes; delete[] uprimes;}

  if(use_stoll)  // using Stoll search so rest not needed.
    return; 
  
  delete[] auxs;
  for(long i=0; i<num_aux; i++) 
    {
      delete[] xgood_mod_aux[i];
      delete[] squares[i];
    }
  delete[] xgood_mod_aux;
  delete[] squares;
  delete[] umod;
}

//#define DEBUG_RANGES

void sort(bigfloat& x1, bigfloat& x2, bigfloat& x3, bigfloat& x4)
// sorts into increasing order
{
#ifdef DEBUG_RANGES
  cout << "sort called with roots "<<x1<<", "<<x2<<", "<<x3<<", "<<x4<<"\n";
#endif
  bigfloat t;
  if(x1>x2) {t=x1; x1=x2; x2=t;}
  if(x2>x3) {t=x2; x2=x3; x3=t;}
  if(x3>x4) {t=x3; x3=x4; x4=t;} // now x4 is biggest
  if(x1>x2) {t=x1; x1=x2; x2=t;}
  if(x2>x3) {t=x2; x2=x3; x3=t;} // now x3 is second biggest
  if(x1>x2) {t=x1; x1=x2; x2=t;} 
#ifdef DEBUG_RANGES
  cout << "sort returns roots "<<x1<<", "<<x2<<", "<<x3<<", "<<x4<<"\n";
#endif
} 

long quartic_sieve::search(double h_lim, long maxnpts, int posxonly)
{
  npoints = 0;
  maxnpoints=maxnpts;
  
  if(easy&1) // a is a square with root roota
    {
      pu=1; pv=roota; pw=0; npoints++;
      //      cout<<"Easy case: "<<a<<" is square, root = "<<roota<<endl;
      //      cout<<"Rational point is ("<<pu<<":"<<pv<<":"<<pw<<")"<<endl;
    }
  if(npoints>=maxnpoints) return npoints;
  if(easy&2) // e is a square with root roote
    {
      pu=0; pv=roote; pw=1; npoints++;
      //      cout<<"Easy case: "<<e<<" is square, root = "<<roote<<endl;
      //      cout<<"Rational point is ("<<pu<<":"<<pv<<":"<<pw<<")"<<endl;
    }
  if(npoints>=maxnpoints) return npoints;

  // Now do some harder work:

  if(use_stoll) return stoll_search(h_lim,posxonly);

  // set initial bounds for point coefficients
  ulim = (long)floor(exp(h_lim));

  if (verbose)
    cout << "quartic_sieve::search: trying u,w up to "<<ulim<<endl;

  int type = g->gettype();
  bigcomplex * roots = g->getroots();
  bigfloat x1, x2, x3, x4,t, zero=to_bigfloat(0);
  switch (type) {
  case 0:  default: // no roots info known
#ifdef DEBUG_RANGES
      cout << "sieve::search: no info about real roots.\n";
#endif
      search_range(0,zero,0,zero,posxonly);
      return npoints;
      break;
    case 1: // no real roots
#ifdef DEBUG_RANGES
      cout << "sieve::search: no real roots.\n";
#endif
      search_range(0,zero,0,zero,posxonly);
      return npoints;
      break;
    case 3:  // 2 real roots, one or two ranges
      x1 = real(roots[2]); x2 = real(roots[3]);
      if(x1>x2) {t=x1; x1=x2; x2=t;}
#ifdef DEBUG_RANGES
      cout << "sieve::search: type 3, real roots "<<x1<<", "<<x2<<"\n";
#endif
      if(a>0)
	{
	  search_range(1,x2,0,zero,posxonly);  //  x2<x
	  if(npoints<maxnpoints) 
	    search_range(0,zero,1,x1,posxonly);  //     x<x1
	}
      else
	{
	  search_range(1,x1,1,x2,posxonly);  //  x1<x<x2
	}
      return npoints;
      break;
    case 2:  // 4 real roots, 2 or 3 ranges
      x1 = real(roots[0]); x2 = real(roots[1]);
      x3 = real(roots[2]); x4 = real(roots[3]);
      sort(x1,x2,x3,x4);  // put in increasing order
#ifdef DEBUG_RANGES
      cout << "sieve::search: type 2, real roots "<<x1<<", "<<x2<<", "<<x3<<", "<<x4<<"\n";
#endif
      if(a>0)
	{
	  search_range(1,x2,1,x3,posxonly);  //  x2<x<x3
	  if(npoints<maxnpoints) 
	    search_range(0,zero,1,x1,posxonly);  //      x<x1
	  if(npoints<maxnpoints) 
	    search_range(1,x4,0,zero,posxonly);  //   x4<x
	}
      else
	{
	  search_range(1,x1,1,x2,posxonly);  //  x1<x<x2
	  if(npoints<maxnpoints) 
	    search_range(1,x3,1,x4,posxonly);  //  x3<x<x4
	}
      return npoints;
    }
}

long quartic_sieve::stoll_search(double h_lim, int posxonly) 
{
  //  cout<<"In quartic_sieve::stoll_search()"<<endl;
  vector<bigint> coef(5);
  coef[0]=e; coef[1]=d; coef[2]=c; coef[3]=b; coef[4]=a;
  qsieve s(this, 4, coef, to_bigfloat(h_lim), verbose);
  // Find and set search intervals

  int type = g->gettype(), lb_flag;
  bigcomplex * roots = g->getroots();
  vector<double> real_roots(4); int nrr=0;
  bigfloat x1, x2, x3, x4, t;
  switch (type) {
  case 0:  case 1: default: // no roots info known or no real roots
#ifdef DEBUG_RANGES
    if(type==0)
      cout << "sieve::stoll_search: no info about real roots.\n";
    else
      cout << "sieve::search: no real roots.\n";
#endif
    s.set_intervals(real_roots,0,0,posxonly);
    break;
  case 3:  // 2 real roots, one or two ranges
    nrr=2;
    x1 = real(roots[2]); x2 = real(roots[3]);
    if(x1>x2) {t=x1; x1=x2; x2=t;}
#ifdef DEBUG_RANGES
    cout << "sieve::search: type 3, real roots "<<x1<<", "<<x2<<"\n";
#endif
    if(doublify(x1,real_roots[0])||doublify(x2,real_roots[1]))
      {
	cout<<"##WARNING##: real roots "<<x1<<", "<<x2<<"are not doubles.\n";
	cout<<"Search will be made over [-height,height]."<<endl;
	nrr=0;
      }
    
    lb_flag = (a<0);
    s.set_intervals(real_roots,nrr,lb_flag,posxonly);
    break;
  
  case 2:  // 4 real roots, 2 or 3 ranges
    nrr=4;
    x1 = real(roots[0]); x2 = real(roots[1]);
    x3 = real(roots[2]); x4 = real(roots[3]);
    sort(x1,x2,x3,x4);  // put in increasing order
    if(doublify(x1,real_roots[0])||doublify(x2,real_roots[1])
       ||doublify(x3,real_roots[2])||doublify(x4,real_roots[3]))
      {
	cout<<"##WARNING##: real roots "<<x1<<", "<<x2<<", "<<x3<<", "<<x4
	    <<"are not all doubles.\n";
	cout<<"Search will be made over [-height,height]."<<endl;
	nrr=0;
      }
    lb_flag = (a<0);
    s.set_intervals(real_roots,nrr,lb_flag,posxonly);
  }
  
  //  cout<<"About to start search..."<<flush;
  npoints = s.search();
  //  cout<<"...finished, "<<npoints<<" points found"<<endl;
  return npoints;

}

long quartic_sieve::search_range(int lower, bigfloat lower_bound, 
				int upper, bigfloat upper_bound, int posxonly)
{
// Adjust bounds so exact roots are not missed by rounding error
#if !defined(LiDIA_ALL) && !defined(NTL_ALL)
  bigfloat eps_adjust(to_bigfloat(0.000001));
  if(lower) lower_bound-=eps_adjust;
  if(upper) upper_bound+=eps_adjust;
#endif
  if(posxonly)  // make sure only positive x are used
    {
      if(upper&&upper_bound<0) return npoints;
      if(lower) // revise or create a lower bound
	{if(lower_bound<0) lower_bound=0;}
      else 
	{lower=1; lower_bound=0;}
    }
#ifdef DEBUG_RANGES
  cout<<"sieve::search_range: ";
  if(lower) cout<<"lower bound = "<<lower_bound<<"; ";
  else cout<<"no lower bound; ";
  if(upper) cout<<"upper bound = "<<upper_bound<<".\n";
  else cout<<"no upper bound.\n";
//  cout<<"a="<<a<<"\n";
//  cout<<"b="<<b<<"\n";
//  cout<<"c="<<c<<"\n";
//  cout<<"d="<<d<<"\n";
//  cout<<"e="<<e<<"\n";
#endif
// declare other loop variables
  bigint w2,w3,w4, aw,bw,cw,dw,ew;
  long i, paw,pbw,pcw,pdw,pew, u, w, aux;
  bigint vsq, v, f;

//
// MAIN LOOP on w (denominator)
//
//#define W_START 259000000          // For debugging purposes!
#define W_START 1
  long wstart=W_START;

  int odd_w_only=0;
  if(nwprimes>0) if(wprimes[0]==2) odd_w_only=1;
  long wstep = 1+odd_w_only;
  if(odd_w_only) 
    {
      if(!odd(wstart)) wstart++;
    }
  int* wflag = new int[10000];

  int odd_u_only=0;
  if(nwprimes>0) if(uprimes[0]==2) odd_u_only=1;
  long ustep = 1+odd_u_only;

  for (w = wstart; (w <= ulim) && (npoints<maxnpoints); w+=wstep)
    {
// set up limits for u-loop
      long first_u = -ulim;
      long last_u = ulim;
      long min_u = I2long(Iceil(w*lower_bound));
      long max_u = I2long(Ifloor(w*upper_bound));
#ifdef DEBUG_RANGES
  cout<<"sieve::search_range: (w="<<w<<"), min_u="<<min_u<<", max_u="<<max_u<<endl;;
#endif
      if(lower)	if(first_u < min_u) first_u = min_u;
      if(upper)	if(last_u  > max_u) last_u  = max_u;
#ifdef DEBUG_RANGES
  cout<<"sieve::search_range: first_u="<<first_u<<", last_u="<<last_u<<endl;;
#endif
      if(odd_u_only)
	{
	  if(!odd(first_u)) first_u++;
	  if(!odd(last_u))  last_u--;
	}

      if(first_u>last_u) {continue;}  // to next w value

// Check that w has no impossible prime factors:
      int w_is_ok = 1; long nwp;
      for(nwp=0; (nwp<nwprimes) && w_is_ok; nwp++)
	w_is_ok = ndiv(wprimes[nwp],w);
      if(!w_is_ok) 
	{
//	  if(verbose) cout << " -- skipping w="<<w<<" (bad prime factor) \n";
	  continue;  // skip to next w
	}

      if (verbose) 
	{
	  cout<<"quartic_sieve::search: trying w = "<<w;
	  cout<<"  ("<<first_u<<" <= u <= "<<last_u<<")\n";
	}

      int use_w_sieve = ((last_u-first_u)>10);
      int use_gcd_table = use_w_sieve&&(w<10000)&&((last_u-first_u)>(w/2));
      long umodw=0; 
      int w_vars_set = 0;

      if(use_w_sieve)
	{

// some preliminary calculations of multiples of w etc.
	  w2 = sqr(BIGINT(w));  w3 = w*w2; w4 = w2*w2;
	  aw = a; bw = b*w; cw = c*w2; dw = d*w3; ew = e*w4;

	  for ( i=0; i < num_aux; i++)  
	    umod[i] = posmod(first_u-ustep, auxs[i]);

// set up flag array of residues coprime to w
	  if(use_gcd_table)
	    {
	      umodw = posmod(first_u-ustep,w);
	      wflag[0]=(w==1);
	      for(i=1; i<=w-i; i++) wflag[i] = wflag[w-i] = (gcd(i,w)==1);
	    }

// set the main flag matrix
	  for (long index = 0; index < num_aux; index++)
	    {
	      aux = auxs[index];
	      paw = posmod(aw, aux);
	      pbw = posmod(bw, aux);
	      pcw = posmod(cw, aux);
	      pdw = posmod(dw, aux);
	      pew = posmod(ew, aux);
	  
	      long ddddf= posmod(24*paw , aux);
	      long dddf = posmod(36*paw + 6*pbw , aux);
	      long ddf  = posmod(14*paw + 6*pbw + 2*pcw , aux);
	      long df   = posmod(paw+pbw+pcw+pdw, aux);
	      long f    = posmod(pew , aux);
	  
	      int* flag = xgood_mod_aux[index];
	      int* sqs = squares[index];
	      long x=aux;
	      while(x--)
		{
		  *flag++ = sqs[f];
		  f    +=    df; if(f    >= aux) f    -= aux;
		  df   +=   ddf; if(df   >= aux) df   -= aux;
		  ddf  +=  dddf; if(ddf  >= aux) ddf  -= aux;
		  dddf += ddddf; if(dddf >= aux) dddf -= aux;
		}
	    }  // end of aux loop
	}      // end of if(use_w_sieve)


     for (u=first_u; (u <= last_u) && (npoints<maxnpoints); u+=ustep)
	{
	  int u_is_ok=1;

// check that u is good for all the auxiliaries
	  if(use_w_sieve)
	    {
	      for ( i=0; (i<num_aux); i++)
		{ long& umodi = umod[i];
		  umodi+=ustep;
		  while (umodi >= auxs[i]) umodi -= auxs[i];
		  if(u_is_ok) 
		    {
		      u_is_ok = xgood_mod_aux[i][umodi];
//		      if(!u_is_ok) 
//		      cout<<"(u,w)=("<<u<<","<<w<<") failed sieve mod "<<auxs[i]<<endl;
		    }
		}
//	      if(u_is_ok) cout<<"(u,w)=("<<u<<","<<w<<") passed sieve "<<endl;
//	      else cout<<"(u,w)=("<<u<<","<<w<<") failed sieve "<<endl;
	    }

// check that gcd(u,w)==1
	  if(use_gcd_table)
	    {
	      umodw+=ustep; while(umodw>=w) umodw-=w;
	      u_is_ok = wflag[umodw];  // true if gcd(u,w)=1
	    }
	  if(!u_is_ok) continue;	  

// Check that u has no impossible prime factors:
	  for(nwp=0; (nwp<nwprimes) && u_is_ok; nwp++)
	    u_is_ok = ndiv(uprimes[nwp],u);

	  if(!u_is_ok) continue;	  

	  if(!w_vars_set)
	    {
	      w2 = sqr(BIGINT(w));  w3 = w*w2; w4 = w2*w2;
	      aw = a; bw = b*w; cw = c*w2; dw = d*w3; ew = e*w4;
	      w_vars_set=1;
	    }
	  f=aw; f*=u; f+=bw; f*=u; f+=cw; f*=u; f+=dw; f*=u; f+=ew;
//        f = ew+u*(dw+u*(cw+u*(bw+u*aw))); 
	  if(isqrt(f,v))
	    {
#ifdef DEBUG_RANGES
	      cout<<"u="<<u<<"\n";
	      cout<<"f="<<f<<"\n";
	      cout<<"v="<<v<<", v^2 = " << v*v << "\n";
#endif
	      npoints++;
	      pu=u; pv=v; pw=w;
	    }
	} // ends u-loop
    } // ends w- loop

  delete[] wflag;
  return npoints;

} // end of quartic_sieve::search_range()


/* END OF FILE MGLOBSOL.CC */
