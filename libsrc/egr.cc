// egr.cc: implementation of functions for reduction of points & component groups
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
 
#include <eclib/matrix.h>
typedef mat_l mat;
typedef long scalar;
#include <eclib/points.h>
#include <eclib/egr.h>

//#define DEBUG_EGR
//#define DEBUG_EGR_EXTRA

// Returns [m] for cyclic of order m, [2,2] for 2*2 (type I*m, m even)

vector<int> ComponentGroups::ComponentGroup(const bigint& p) const
{
  vector<int> ans(1);
  if(p==0) ans[0]=conncomp; // 1 or 2
  else
    {
      ans[0]=1;
      auto ri = reduct_array.find(p);
      if(ri==reduct_array.end()) return ans; // p has good reduction    
      ans[0] = (ri->second).c_p;      // usual case: cyclic of order cp
      int code=(ri->second).Kcode.code;
      if((code%10==1)&&even((code-1)/10)) // Type I*m, m even: [2,2]
	{ans[0]=2; ans.push_back(2);}
    }
  return ans;
}

// Returns 1 iff P and Q have same image in the component group at p:
//
int ComponentGroups::InSameComponent(const Point& P, const Point& Q, const bigint& p) const
{
  if(P==Q) return 1;
  return (P-Q).has_good_reduction(p);
}

// For reduction type Im, multiplicative reduction where component
// group is cyclic of order m, using Tate curve formula from Silverman
// Returns a such that P mod pr maps to +a or -a mod m in the
// component group

long ComponentGroups::ImageInComponentGroup_Im_pm(const Point&P, const bigint& p, int m) const
{
#ifdef DEBUG_EGR
  cout<<"In ImageInComponentGroup_Im_pm() with point "
          <<P
      <<", m="<<m<<"..."<<flush;
#endif
  if(P.has_good_reduction(p)) return 0;
  bigint x, y, z; P.getcoordinates(x,y,z);
  bigint zroot = gcd(x,z); // = cube root of z
  long ans = val(p, 2*y + a1*x + a3*z) - 3*val(p,zroot);
#ifdef DEBUG_EGR
  cout<<"ImageInComponentGroup_Im_pm() first gives ans = "<<ans<<endl;
#endif 
  if(even(m)&&(ans>(m/2))) ans=m/2;
#ifdef DEBUG_EGR
  cout<<"ImageInComponentGroup_Im_pm() returns "<<ans<<endl;
#endif 
  return ans;
}

// For reduction type Im, multiplicative reduction where component
// group is cyclic of order m, using full Tate curve formula.
// Returns a such that P mod pr maps to a mod m in the component group

long ComponentGroups::ImageInComponentGroup_Im(const Point&P, const bigint& p, int m) const
{
#ifdef DEBUG_EGR
  cout<<"In ImageInComponentGroup_Im() with point "
          <<P
      <<", m="<<m<<"..."<<flush;
#endif
  if(P.has_good_reduction(p)) return 0;
  // The following is independent of P
  long N = m;  // to match the write-up
  //  long N2=(N+1)/2;  // =N/2 rounded up
  long N2=N;
  bigint pN = pow(p,N2);
  bigint c4inv = invmod(c4,pN);
  bigint x0 = c4inv*(18*b6-b2*b4);
  bigint y0 = c4inv*(a1*a1*a1*a4 - 2*a1*a1*a2*a3 + 4*a1*a2*a4 + 3*a1*a3*a3 - 36*a1*a6 - 8*a2*a2*a3 + 24*a3*a4);
  x0 = x0 % pN;
  y0 = y0 % pN;
#ifdef DEBUG_EGR_EXTRA
  cout<<"c4inv = "<<c4inv<<endl;
  cout<<"x0 = "<<x0<<", y0="<<y0<<endl;
#endif
  bigint d2 = b2+12*x0;
  bigint d;  sqrt_mod_p_power(d,d2,p,N2);
  bigint alpha1 = (d-a1);
  if(odd(alpha1))
    {
      if(p==2)
        cout<<"Problem in ComponentGroups::ImageInComponentGroup_Im(): "
            <<"quadratic has no roots\n";
      else alpha1+=pN; // saves inverting 2 mod pN
    }
  alpha1/=2;
  bigint alpha2 = (-a1-alpha1);
#ifdef DEBUG_EGR_EXTRA
  cout<<"alpha1 = "<<alpha1<<", alpha2="<<alpha2<<endl;
#endif

  bigint x, y, z, w1, w2;
  P.getcoordinates(x,y,z);
  bigint zinv = invmod(z,pN);
  x=(x*zinv-x0) %pN;
  y=(y*zinv-y0) %pN;
  w1 = (y - alpha1*x) %pN;
  w2 = (y - alpha2*x) %pN;
  long e1 = val(p, w1); if(e1>N2) e1=N2;
  long e2 = val(p, w2); if(e2>N2) e2=N2;
#ifdef DEBUG_EGR_EXTRA
  cout<<"x = "<<x<<", y="<<y<<endl;
  cout<<"(y-y0) - alpha1*(x-x0) = "<<w1<<endl;
  cout<<"(y-y0) - alpha2*(x-x0) = "<<w2<<endl;
  cout<<"e1="<<e1<<", e2="<<e2<<endl;
#endif 
  long ans=0;
  if(e1<e2) 
    {
      if((0<e1)&&(2*e1<N))      ans=-e1; 
      else
	cout<<"Problem! e1 not between 1 and "<<((N+1)/2)-1<<endl;

    }
  else
    if(e2<e1) 
      {
	if((0<e2)&&(2*e2<N)) 	ans=+e2; 
	else
	  cout<<"Problem! e2 not between 1 and "<<((N+1)/2)-1<<endl;
      }
    else
      {
	if(even(N))   ans=N/2;
	else
	  cout<<"Problem! e1=e2="<<e1<<" but N="<<N<<"is not even!"<<endl;
      }
#ifdef DEBUG_EGR
  cout<<"ImageInComponentGroup_Im() gives ans = "<<ans<<endl;
#endif 
  return ans;
}

long ComponentGroups::ImageInComponentGroup(const Point&P, const bigint& p, vector<int> grp) const
{
  int ans=0;      // the default
#ifdef DEBUG_EGR
  cout<<"In ImageInComponentGroup() with point "
          <<P
      <<", p = "<<p<<", group = "<<grp<<"..."<<flush;
#endif
  if(grp.size()==2) // C2xC2, cannot handle
    {
      cerr<<"Error in ComponentGroups::ImageInComponentGroup(): noncyclic case"<<endl;
      return ans;
    }

  long n=grp[0];  // the group is cyclic of order n
  switch(n) {
  case 1: {break;}
  case 2: case 3: {
    ans = (P.has_good_reduction(p)? 0 : 1);
    break;}
  case 4: {
    ans = (P.has_good_reduction(p)? 0 : ((2*P).has_good_reduction(p)? 2: 1));
    break;}
  default:
    // Now we are in the case of I_n
    {
      ans= ImageInComponentGroup_Im(P,p,n);
#if(0) // testing only
      for(int i=0; i<n; i++) 
	{
	  int t = (i*ans-ImageInComponentGroup_Im(i*P,p,n)) % n;
	  if(t) cout<<"Problem!"<<endl;
	}
#endif
    }
  }
#ifdef DEBUG_EGR
  cout<<"ImageInComponentGroup() returns "<<ans<<endl;
#endif 
  return ans;
 
}


//#undef DEBUG_EGR
//#define DEBUG_EGR

// Return least j>0 such that j*P has good reduction at p; the
// component group order is given so we only test j dividing this
//
// Since ComponentGroups are small we use nothing fancy here...

int ComponentGroups::OrderInComponentGroup(const Point& P, const bigint& p, const vector<int>& grp)
{
#ifdef DEBUG_EGR
  cout<<"In OrderInComponentGroup() with point "
          <<P
      <<", p = "<<p<<", group = "<<grp<<"..."<<flush;
#endif
  int ans=1;
  long n=grp[0];  // the group order

  if(grp.size()==2) // C2xC2
    {
      if(!P.has_good_reduction(p)) ans=2;
    }
  else // cyclic of order n
    {
      ans = n/gcd(n,ImageInComponentGroup(P,p,grp));
    }
#ifdef DEBUG_EGR
  cout<<"OrderInComponentGroup() returns "<<ans<<endl;
#endif
  return ans;
}

//#define DEBUG_EGR

// replace (independent) points in Plist with a new set which spans
// the subgroup of the original with good reduction at p, returning
// the index

int ComponentGroups::gr1prime(vector<Point>& Plist, const bigint& p)
{
  int j,k,m,n,n0,n1,npts=Plist.size();
#ifdef DEBUG_EGR
  cout<<"in gr1prime with p="<<p<<endl;
  //  bigfloat reg0 = regulator(Plist);
  //  cout<<"regulator = "<<reg0<<endl;
  //  cout<<n<<" points"<<endl;
#endif
  if(Plist.empty()) return 1;
  Point P0 = Plist[0], P1;
  vector<int> CG = ComponentGroup(p);
  long CGexpo=CG[0];
  long CGOrder=CGexpo;
  if(CG.size()>1) CGOrder*=CG[1]; // =4
#ifdef DEBUG_EGR
  cout<<"Component group structure = "<<CG<<endl;
#endif
  if(CG.size()>1)
    {
#ifdef DEBUG_EGR
      cout<<"Non-cyclic component group "<<endl;
#endif
      int logm=0;  // = no. of gens so far
      m = 1;       // = order so far = 2^logm
      Point Third; // will hold sum of first two gens

      for (k=0; k<npts; k++)
	{
	  Point Pk=Plist[k];
#ifdef DEBUG_EGR
	  cout<< "processing point#"<<k<<": "<<Pk<<endl;
#endif
	  j=OrderInComponentGroup(Pk,p,CG);
#ifdef DEBUG_EGR
	  cout<< "image has order "<<j<<" in component group"<<endl;
#endif
	  if (2%j)
	    {
	      cout<< "Error:  order of "<<Pk
		  <<" in the component group is not 1 or 2!"<<endl;
	      return 0;
	    }

	  if (j==1) continue; // good reduction, nothing to do

	  switch(logm)
	    {
	    case 0:  // no generators yet, we have the first generator
	      {
		if(k>0)	// Make this the first point in list
		  {
		    Plist[k]=Plist[0];
		    Plist[0]=Pk;
		    Pk=Plist[k];
		  }
		logm+=1; m*=2;
		break;
	      } // end of case 0

	    case 1: // one generator so far: do we have a second generator?
	      {
		if (InSameComponent(Pk,Plist[0],p)) // no
		  {
		    Plist[k]=Pk-Plist[0];
		  }
		else                                  // yes
		  {
		    if(k>1)	// Make this the second point in list
		      {
			Plist[k]=Plist[1];
			Plist[1]=Pk;
			Pk=Plist[k];
		      }
		    Third = Plist[0]+Plist[1];
		    logm+=1; m*=2;
		  }
		break;
	      } // end of case 1

	    case 2: // we already have 2 generators
	      {
		if (InSameComponent(Pk,Plist[0],p))
		  {
		    Plist[k]=Pk-Plist[0];
		  }
		else
		  if (InSameComponent(Pk,Plist[1],p))
		    {
		      Plist[k]=Pk-Plist[1];
		    }
		  else
		    {
		      if (InSameComponent(Pk,Third,p))
			{
			  Plist[k]=Pk-Third;
			}
		      else
			{
			  cerr<<"Problem in non-cyclic component group case!"<<endl;
			  return 0;
			}
		    }
	      } // end of case 2
	    } // end of switch on logm
	}  // end of loop on points
      if(logm>0) Plist[0]=2*Plist[0];
      if(logm>1) Plist[1]=2*Plist[1];
      return m;
    }  // end of cyclic switch

#ifdef DEBUG_EGR
  cout<<"Cyclic component group of order "<<CGOrder<<endl;
#endif

  if(CGOrder==1) return 1;
#ifdef DEBUG_EGR
  cout<< "processing point #1: "<<flush;
#endif
  m=OrderInComponentGroup(P0,p,CG); // will hold current index
  for(j=1; j<npts; j++)
    {
      n0 = m;
      P0=Plist[0];  // NB this will change
      P1=Plist[j];
#ifdef DEBUG_EGR
      cout<< "processing point #"<<(j+1)<<": "<<flush;
#endif
      n1 = OrderInComponentGroup(P1,p,CG);
      if(n1==1)
	{
#ifdef DEBUG_EGR
	  cout<< "point has good reduction, no action needed "<<endl;
#endif     
	  continue;
	}
#ifdef DEBUG_EGR
      cout<< "image has order "<<n1<<" in component group"<<endl;
#endif
      // swap points so n0>=n1 for convenience
      if(n1>n0) 
	{
#ifdef DEBUG_EGR
	  cout<<" swapping "<<P0<<" and "<<P1<<endl;
#endif
	  Point Q=P0; P0=P1; P1=Q; 
	  n=n0; n0=n1; n1=n; m=n0;
	}

      if((n0%n1)==0) // P1 is a multiple of P0 in CG...
	{
	  while((!P1.has_good_reduction(p))) {P1=P1-P0;}
#ifdef DEBUG_EGR
	  cout<<"P1 replaced by  "<<P1<<" with good reduction"<<endl;
	  cout<<"index is now "<<m<<endl;
#endif
	  Plist[0]=P0;
	  Plist[j]=P1;
	}
      else // lcm(n0,n1)>n0: we gain something
	{
	  long a,b,g=gcd(n0,n1);
	  // Now find u (coprime to g) s.t. (n1/g)P1 == u* (n0/g)P0
	  Point Q0=(n0/g)*P0; Point Q=Q0; // (P0.getcurve());
	  Point Q1=(n1/g)*P1;
	  long u=1; //0;
#ifdef DEBUG_EGR
	  cout<<"Looking for u"<<endl;
#endif
	  while(gcd(u,CGOrder)>1 || (!InSameComponent(Q,Q1,p))) {u++; Q=Q+Q0;}
#ifdef DEBUG_EGR
	  cout<<" u="<<u<<endl;
#endif
	  g=bezout(u*n0,n1,a,b);
#ifdef DEBUG_EGR
	  cout<<" a="<<a<<", b="<<b<<", g="<<g<<endl;
#endif
	  Point newP0 = b*P0+a*P1;
	  P1=(-u*n0/g)*P0+(n1/g)*P1;
	  P0=newP0;
#ifdef DEBUG_EGR
	  cout<<"new index = lcm="<<((n0*n1)/g)<<endl;
#endif
	  m= ((n0*n1)/g);
	  Plist[0]=P0; // this now has index m
	  Plist[j]=P1; // this has good reduction
	}
    }
  // At this point all points after Plist[0] have good reduction, and
  // Plist[0] has index m
  Plist[0]=m*Plist[0];
#ifdef DEBUG_EGR
  cout<<" gr1prime returns index "<<m<<", points "<<Plist<<endl;
  //  bigfloat reg1 = regulator(Plist);
  //  cout<<" new regulator = "<<reg1<<endl;
  //  cout<<" ratio = "<<reg1/reg0<<endl;
#endif
  return m;
}

// replaces the (independent) points with a new set which spans the
// subgroup of the original with good reduction at all p in plist,
// returning the overall index
int ComponentGroups::grprimes(vector<Point>& Plist, const vector<bigint>& plist)
{
#ifdef DEBUG_EGR
  cout<<"in grprimes with plist="<<plist<<endl;
#endif
  int m = std::accumulate(plist.begin(), plist.end(), 1,
                          [this, &Plist] (int m, const bigint& p) {return m*gr1prime(Plist,p);});
#ifdef DEBUG_EGR
  cout<<" grprimes returns index "<<m<<endl;
#endif
  return m;
}

// replaces the (independent) points with a new set which spans the
// subgroup of the original with good reduction at all p,
// returning the overall index
int ComponentGroups::egr_subgroup(vector<Point>& Plist, int real_too)
{
  static const bigint zero(0);
  if(Plist.size()==0) return 1;
  vector<bigint> plist = the_bad_primes;
  if(real_too && (conncomp==2)) plist.push_back(zero);
#ifdef DEBUG_EGR
  cout<<"Using primes "<<plist<<endl;
#endif
  return grprimes(Plist,plist);
}

bigint comp_map_image(const vector<int> moduli, const mat& image);

//#undef DEBUG_EGR
//#define DEBUG_EGR

bigint egr_index(const vector<Point>& Plist, int real_too)
{
  static const bigint zero(0), one(1);
  if(Plist.size()==0) return one;;

  // Compute minimal model and map points to it if necessary

#ifdef DEBUG_EGR
  cout<<"Plist = "<<Plist<<" (size "<<Plist.size()<<")\n";
#endif
  Point P0 = Plist[0];
  Curve E = P0.getcurve();
  Curvedata E_orig = Curvedata(E, 0); // don't minimise yet
  bigint u, r, s, t;
  Curvedata Emin = E_orig.minimalize(u,r,s,t);
  int is_min = is_one(u) && is_zero(r) && is_zero(s) && is_zero(t);
  // then E was already reduced and minimal, no adjustments needed

#ifdef DEBUG_EGR
  cout<<"In egr_index("<<E<<"): ";
  if (is_min)
    cout << "already minimal and reduced\n";
  else
    cout << " minimal reduced curve is " << (Curve)Emin <<"\n";
#endif

  vector<Point> Plist_min = Plist;
  if (!is_min)
    std::transform(Plist.begin(), Plist.end(), Plist_min.begin(),
                   [&Emin, u, r, s, t] ( const Point& P ) {return transform(P, &Emin, u, r, s, t);});

  // Construct the ComponentGroups class from the minimal model

  ComponentGroups CGS(Emin);
  vector<bigint> plist = getbad_primes(CGS);
  if(real_too && (getconncomp(CGS)==2)) plist.push_back(zero);
#ifdef DEBUG_EGR
  cout<<"Using primes "<<plist<<endl;
#endif

  vector<vector<vector<int> > > imagematrix;
  vector<int> moduli;
  int n=0;

  // Loop through bad primes (possibly including infinity), computing
  // the image in the local component group

  for( const auto& p : plist)
    {
#ifdef DEBUG_EGR
      cout<<"p = "<<p<<endl;
#endif
      vector<vector<int> > im=MapPointsToComponentGroup(CGS,Plist_min,p);
#ifdef DEBUG_EGR
      cout<<"image = ";
      for(unsigned int j=0; j<im.size(); j++) cout << im[j] << " ";
      cout<<endl;
#endif
      imagematrix.push_back(im);
      vector<int> CG=CGS.ComponentGroup(p);
      for(unsigned int ni=0; ni<CG.size(); ni++, n++)
	moduli.push_back(CG[ni]);
    }

  // Compute the kernel

  mat m(Plist.size(),n);
  unsigned int j=0, j1, j2, i;
  bigint imageorder=one;
  for(i=0; i<moduli.size(); i++) imageorder*=moduli[i];
  for(j1=0; j1<plist.size(); j1++)
    for(j2=0; j2<imagematrix[j1][0].size(); j2++)
      {
	for(i=0; i<Plist.size(); i++)
	  m(i+1,j+1)=imagematrix[j1][i][j2];
	j++;
      }
#ifdef DEBUG_EGR
  cout<<"Moduli:      = "<<moduli<<endl;
  cout<<"Image matrix = "<<m<<endl;
  cout<<"Maximum image order = "<<imageorder<<endl;
#endif
  imageorder = comp_map_image(moduli,m);
#ifdef DEBUG_EGR
  cout<<"Actual image order = "<<imageorder<<endl;
#endif
  return imageorder;
}

// Given a list of points P1,...,Pn and a prime p, this returns a
// vector [c1,c2,...,cn] where the image of Pi in the component group
// is ci mod m, where m is the exponent of the component group at p.
// 
// Each ci is a vector of length 1 or 2 (the latter for when the
// component group is C2xC2), not just an integer.
//
// If p=0 then m=1 or 2 (m=2 iff there are two real components and at
// least one P_i is not in the connected component)
//

vector<vector<int> >  MapPointsToComponentGroup(const ComponentGroups& CG, const vector<Point>& Plist,  const bigint& p)
{
  int n=Plist.size();
  vector<vector<int> > images;
  if (n==0) return images;

  // Construct the local component group and find its structure:

  vector<int> G=CG.ComponentGroup(p);
#ifdef DEBUG_EGR
  cout<<"Component group = "<<G<<endl;
#endif
  int m = G.size();
  int cyclic = (m==1);
  int orderG = (cyclic? G[0]: 4);

  // Initialize the images to 0:

  images.resize(n, vector<int>(m,0));

  if (orderG==1) return images;

  if (cyclic) // Now G is cyclic and nontrivial
    {
#ifdef DEBUG_EGR
      cout<< "cyclic case (order "<<orderG<<")"<<endl;
#endif
      for(int i=0; i<n; i++)
        images[i][0]=CG.ImageInComponentGroup(Plist[i],p,G);

      // if order =3 or =4, check for compatibility since our map is
      // only then defined up to sign....
      if((m==3)||(m==4))
	{
	  // Find a point with image +1, if any:
	  int i0=-1;
	  for(int i=0; i<n; i++)
            {
              if(images[i][0]==1)
                {
                  i0=i;
                  break;
                }
            }
	  if(i0!=-1) // else nothing to do
	    {
	      Point P0=Plist[i0];
	      for(int i=i0+1; i<n; i++)
		if(images[i][0]==1)
		  if(!CG.InSameComponent(P0,Plist[i],p))
		    images[i][0]=-1;
	    }
	}  // end of special treatment for m=3,4
    } // end of cyclic case
  else
    {
#ifdef DEBUG_EGR
      cout<< "non-cyclic case"<<endl;
#endif
      // points representing up to 3 nontrivial components:
      vector<Point> PointReps;  
      // the three nonzero images:
      vector<vector<int> > ims(3);
      ims[0]=vector<int>(2); ims[0][0]=1; ims[0][1]=0;
      ims[1]=vector<int>(2); ims[1][0]=0; ims[1][1]=1;
      ims[2]=vector<int>(2); ims[2][0]=1; ims[2][1]=1;
      for (int k=0; k<n; k++)
	{
	  Point Pk=Plist[k];
#ifdef DEBUG_EGR
	  cout<< "processing point#"<<k<<": "<<Pk<<endl;
#endif
	  if(Pk.has_good_reduction(p)) continue;
	  int coset=-1;
	  for(unsigned int j=0; (j<PointReps.size())&&(coset==-1); j++)
	    if(CG.InSameComponent(Pk,PointReps[j],p))
	      coset=j;
	  if(coset==-1) // Pk is in a new coset...
	    {
#ifdef DEBUG_EGR
	      cout<<"Pk is in a new coset"<<endl;
#endif
	      coset=PointReps.size();
	      PointReps.push_back(Pk);
	    }
	  images[k]=ims[coset];
#ifdef DEBUG_EGR
	  cout<<"Pk is in coset #"<<(coset+1)<<", image = "<<images[k]<<endl;
#endif
	}  // loop on points
    } // else ... noncyclic case
  return images;
}

// class Kodaira_code just holds an int which "codes" the type as follows:
// (this coding originally from R.G.E.Pinch)
//
// Im                 -> 10*m
// I*m                -> 10*m+1
// I, II, III, IV     -> 1, 2, 3, 4
// I*, II*. III*, IV* -> 5, 6, 7, 8
//

//#define DEBUG_INDEX

bigint comp_map_image(const vector<int> moduli, const mat& image)
{
  bigint ans; ans=1;
  mat m=image;
#ifdef DEBUG_INDEX
  cout<<"In comp_map_image, m="<<m;
  cout<<"moduli = "<<moduli<<endl;
#endif
  int npts=m.nrows(), np=m.ncols();
  int i, j, jj;
  if(np==0) return ans;
  for(j=1; j<=np; j++)
    {
      scalar modulus=moduli[j-1];
#ifdef DEBUG_INDEX
      cout<<"Working on column "<<j<<", modulus "<<modulus<<endl;
#endif
      if(modulus==1) continue;
#ifdef DEBUG_INDEX
      cout<<"Column "<<j<<" = "<<m.col(j)<<endl;
#endif
      for(i=1; (i<=npts); i++) m(i,j)=m(i,j)%modulus;
      scalar g=0,gm;
      for(i=1; (i<=npts)&&(g!=1); i++) g=gcd(g,m(i,j));
#ifdef DEBUG_INDEX
      cout<<"Column gcd = "<<g<<endl;
#endif
      if(g==0) continue;
      if(g>1)
	{
	  gm=gcd(g,modulus);
	  if(gm>1) 
	    {
	      modulus/=gm; g/=gm;
	      for(i=1; i<=npts; i++) m(i,j)=(m(i,j)/gm)%modulus;
	    }
	  if(g>1) for(i=1; i<=npts; i++) m(i,j)=(m(i,j)/g)%modulus;
	}
#ifdef DEBUG_INDEX
      cout<<"After scaling,  modulus = "<<modulus<<endl;
#endif
      if(modulus==1) continue;
#ifdef DEBUG_INDEX
      cout<<"Column "<<j<<" = "<<m.col(j)<<endl;
#endif
      long colmin=modulus, imin=0, r, q;
      while(abs(colmin)>1)
	{
#ifdef DEBUG_INDEX
	  cout<<"colmin = "<<colmin<<endl;
#endif
	  for(i=1; i<=npts; i++) 
	    if(m(i,j)!=0) 
	      if(abs(m(i,j))<abs(colmin)) colmin=m(i,j);
	  for(i=1; i<=npts; i++) if(m(i,j)==colmin) {imin=i; break;}
#ifdef DEBUG_INDEX
	  cout<<"colmin = "<<colmin<<" at row imin="<<imin<<endl;
#endif
	  for(i=1; i<=npts; i++)
	    if(ndivides(colmin,m(i,j)))
	      {
		r = m(i,j) % colmin;
		q = (m(i,j) - r) / colmin;
#ifdef DEBUG_INDEX
		cout<<"subtracting "<<q<<" times row "<<imin<<" from row "<<i<<endl;
#endif
		for(jj=1; jj<=np; jj++) m(i,jj)-=q*m(imin,jj);
	      }
#ifdef DEBUG_INDEX
	  cout<<"Column "<<j<<" = "<<m.col(j)<<endl;
#endif
	}
      // now |colmin|=1 and we can clear
#ifdef DEBUG_INDEX
      cout<<"colmin = "<<colmin<<", clearing..."<<endl;
      cout<<"Column "<<j<<" = "<<m.col(j)<<endl;
#endif
      for(i=1; i<=npts; i++)
        if(m(i,j)==colmin)
          {
            imin=i; break;
          }
      for(i=1; i<=npts; i++)
	if(i!=imin)
          {
            long piv = m(i,j)*colmin;
            for(jj=1; jj<=np; jj++)
              {
		m(i,jj) = (m(i,jj)-piv*m(imin,jj)) % moduli[jj-1];
              }
          }
      // now update matrix and ans:
#ifdef DEBUG_INDEX
      cout<<"Now matrix = "<<m<<endl;
      cout<<"Column = "<<m.col(j)<<endl;
      cout<<"multiplying ans and row "<<imin<<" by "<<modulus<<endl;
#endif
      ans*=modulus;
      for(jj=1; jj<=np; jj++)
	m(imin,jj)=(modulus*m(imin,jj))%moduli[jj-1];
#ifdef DEBUG_INDEX
      cout<<"Now matrix = "<<m<<endl;
      cout<<"and ans = "<<ans<<endl;
#endif
    }
  return ans;
}
