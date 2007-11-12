// tlss.cc: implementation of class TLSS for sieving E(Q)/pE(Q) at one prime q
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
 

// NB: TLSS = Tate--Lichtenbaum--Samir-Siksek: we use a simple
// discrete log a la Siksek when the p-torsion in E(F_q) is cyclic,
// else use the Tate-Lichtenbaum pairing

#include "matrix.h"
#include "subspace.h"

#include "points.h"
#include "polys.h"
#include "curvemod.h"
#include "pointsmod.h"
#include "ffmod.h"
#include "tlss.h"

void TLSS::init(int  pp, int verb)
{
  verbose=verb;
  p=pp;

  Pi=Emodq.get_pbasis(p);
  rank=Pi.size();

  if((verbose>1)&&(rank>0))
    {
      cout<<"Generators of "<<p<<"-torsion mod "<<q<<": \n";
      cout<<"P1 = "<<Pi[0]<<endl;
      if(rank>1) cout<<"P2 = "<<Pi[1]<<endl;
    }  
  if(rank==2) init_tlpolys();
}

void TLSS::init(int  pp, const vector<bigint>& pdivpol, int verb)
{
  verbose=verb;
  p=pp;

  Pi=Emodq.get_pbasis_via_divpol(p, pdivpol);
  rank=Pi.size();

  if((verbose>1)&&(rank>0))
    {
      cout<<"Generators of "<<p<<"-torsion mod "<<q<<": \n";
      cout<<"P1 = "<<Pi[0]<<endl;
      if(rank>1) cout<<"P2 = "<<Pi[1]<<endl;
    }  
  if(rank==2) init_tlpolys();
}

void TLSS::init_tlpolys()
{
 
  if (rank<2) return;

  // case rank=2: prepare to use TL map
  // first find p'th root of 1 mod q 

  // initialize the array of all p'th roots mod q
  q1p=(q-1)/p;
  mu_p = roots_of_unity(Fq,p);

  if(verbose>1) 
    {
      cout<<"q="<<q<<endl;
      cout<<"p'th roots of unity mod q = "<<mu_p<<endl;
      cout<<"Rank of p-torsion mod q = "<<rank<<endl;
    }

  // initialize the ffmodq class
  ffmodq dummy((const curvemodq)Emodq);
  // initialize the TL-functions
  TLpolys.resize(0);
  int i;
  for(i=0; i<rank; i++) TLpolys.push_back(weil_pol(Pi[i],p));
  if(verbose>1) for(i=0; i<rank; i++) cout<<"TL poly: "<<TLpolys[i]<<endl;
}

//#define debugTL

// apply map to P, result in (rank*)[0..p-1]:
vector<int> TLSS::map1point(const Point& P) const
{
#ifdef debugTL
  cout<<"Applying TLSS::map1point (q="<<q<<") to P="<<P<<endl;
#endif
  pointmodq Pmodq = reduce_point(P,(const curvemodq&)Emodq);
#ifdef debugTL
  cout<<"P mod q ="<<Pmodq<<endl;
#endif
  int i, sw;
  vector<int> ans;
  ans.resize(rank);
  for(i=0; i<rank; i++) ans[i]=0;
  if(Pmodq.is_zero()) return ans;
  
  if(rank==1)
    {
      pointmodq P1=Pi[0];
      bigint n1 = Emodq.n1; // order of relevant cyclic factor if known
      if(n1==0) n1 = Emodq.n; // else full group order if not
      Pmodq = I2long(n1/p) * Pmodq;
      if(Pmodq.is_zero()) return ans;

#ifdef debugTL
      cout<<"after multiplying by "<<I2long(n1/p)<<", get "<<Pmodq<<endl;
      cout<<"finding discrete log of "<<Pmodq<<" w.r.t. "<<P1<<endl;
#endif

      // Method 1: use LiDIA's discrete log (was buggy)
      /*
      ans[0]=I2long(bg_algorithm(P1,Pmodq,0,p-1));
      */

      // Method 2: elementary method (OK for small p!)

#ifdef debugTL
      cout<<"discrete log: finding which multiple of "<<P1<<" is "<<Pmodq<<"..."<<endl;
#endif
      pointmodq Q(Emodq);
      pointmodq minusPmodq=-Pmodq;
      for(i=0; i<p; i++) 
	{
	  if(Q==Pmodq) {ans[0]=i; return ans;}
	  if(i) if(Q==minusPmodq) {ans[0]=p-i; return ans;}
	  Q+=P1;
	}

      // Optional check:      

#ifdef debugTL
      int check = ans[0]*P1 == Pmodq;
      if(!check)
	{
	  cout<<"Error: discrete log of "<<Pmodq<<" w.r.t. "<<P1<<" returns "<<ans[0]<<endl;
	  abort();
	}
#endif

      return ans;
    }

  // else apply TL maps

  gf_element xP=Pmodq.get_x();
  gf_element yP=Pmodq.get_y();
  gf_element lambda, mu, t;
  gf_element a1,a2,a3,a4,a6;
  Emodq.get_ai(a1,a2,a3,a4,a6);
  gf_element b2 = a1*a1 + 4*a2; 
  gf_element b4 = 2*a4 + a1*a3;

  for(i=0; i<rank; i++)
    {  
      const pointmodq& T = Pi[i];
      gf_element xT=T.get_x(), yT=T.get_y();
#ifdef debugTL
      cout<<"(xT,yT)=("<<xT<<","<<yT<<")"<<endl;
#endif
      switch(p) {
      case 2:{
	t=xP-xT;
	if(t==0)
	  {
	    lambda=4*xT+b2;
	    mu=xT*lambda+b4+b4;
	    t=4*xP*xP+lambda*xP+mu;
	  }
#ifdef debugTL
	cout<<"calling power_mod with (t,q1p,q) = ("<<mu<<","<<t<<","<<q1p<<","<<q<<")"<<endl;
#endif
	power(mu,t,q1p);
#ifdef debugTL
	cout<<"mu = "<<mu<<endl;
#endif
	ans[i] = ((mu==1)? 0 : 1);	
	break;
      }
      default:{
#ifdef debugTL
	cout<<"applying TL poly = "<<TLpolys[i]<<" to "<<Pmodq<<endl;
#endif
	t=TLpolys[i](Pmodq);
	sw=0;
	if(t==0) 
	  {
	    sw=1; 
	    t=TLpolys[i](-Pmodq);
	  }
	if(t==0) 
	  {
	    cout<<"Error: both P and -P map to 0"<<endl;
	    abort();
	  }
	power(mu,t,q1p);

#ifdef debugTL
	cout<<"t="<<t<<endl;
	cout<<"mu = "<<mu<<endl;
	gf_element t2, mu2;
	if(!(p*Pmodq).is_zero()) t2 = evaluate_weil_pol(T,p,Pmodq);
	else
	  {
	    cout<<"(new method having to use shifting trick)"<<endl;
	    pointmodq R=Pmodq.get_curve().random_point();
	    while((p*R).is_zero())
	      R=Pmodq.get_curve().random_point();
	    t2 = evaluate_weil_pol(T,p,R+Pmodq)/evaluate_weil_pol(T,p,R);
	  }
	power(mu2,t2,q1p);
	if(mu==mu2) cout<<"New method agrees"<<endl;
	else cout<<"New method gives mu value "<<mu2<<" instead!"<<endl;
#endif

	if(mu==1) ans[i]=0; else // discrete log: mu is some power of mu_p
	  {
	    int dl=find(mu_p.begin(),mu_p.end(),mu)-mu_p.begin();
#ifdef debugTL
	    cout<<"discrete log of mu = "<<dl<<endl;
#endif
	    if(dl==p) 
	      {
		cout<<"Error:  mu="<<mu<<" (mod "<<q<<") is not a "<<p<<"'th root of unity!"<<endl;
		abort();
	      }
	    ans[i] = (sw? p-dl: dl);
	  }
	break;
      }  //  end of default case
      } // switch(p)
    }
#ifdef debugTL
  cout<<"ans = "<<ans<<endl;
#endif
  return ans;
}

// apply map to all P in Plist, result is a (rank*#Plist) matrix:
mat TLSS::map_points(const vector<Point>& Plist) const
{
  int npts = Plist.size();
  mat TLim(rank,npts);
  int i,j;
  for(i=0; i<npts; i++)
    {
      Point P=Plist[i];
      vector<int> tlP = map1point(P);
      if(verbose>1) cout<<"P="<<P<<" -> "<<tlP<<endl;
      for(j=0; j<rank; j++)
	{
	  TLim(j+1,i+1)=tlP[j];
	}
    }
  return TLim; 
}
