// mequiv.cc: implementation of quartic equivalence functions
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
 
#include "mequiv.h"

//#ifdef NEW_EQUIV

int new_equiv(quartic* q1, quartic* q2, int info)
{
  if(info)
    {
      cout<<"Checking equivalence of " << *q1 << " and " << *q2 << "\n";
    }
  bigint& ii=q1->ii;
  bigint& jj=q1->jj;
  if (!(ii==q2->ii && jj==q2->jj && q1->disc==q2->disc && q1->type== q2->type))
    { 
      if (info) 
	{
	  cout << "equiv failed on first test!\n";
	  cout << "First  has I="<<q1->ii<<", J="<<q1->jj<<",";
	  cout << " disc="<<q1->disc<<", type="<<q1->type<<endl;
	  cout << "Second has I="<<q2->ii<<", J="<<q2->jj<<",";
	  cout << " disc="<<q2->disc<<", type="<<q2->type<<endl;
	}
      return 0;
    }

  if(q1->equiv_code!=q2->equiv_code)
    {
      if(info) cout << "--equiv_codes not equal\n";
      return 0;
    }

  q1->make_zpol();
  q2->make_zpol();
  bigint& p1=q1->p;      bigint& p2=q2->p;
  bigint& r1=q1->r;      bigint& r2=q2->r;
  bigint& p1sq=q1->psq;  bigint& p2sq=q2->psq;
  bigint& a1=q1->a;      bigint& a2=q2->a;
  bigint& a1sq=q1->asq;  bigint& a2sq=q2->asq;
  const bigint& a1a2=a1*a2;
  const bigint& p1p2=p1*p2;
  const bigint& p = (32*a1a2*ii + p1p2)/3;
  const bigint& s = (-256*jj*a1a2*(a1*p2+a2*p1) 
	       + 64*ii*(p1sq*a2sq+p2sq*a1sq+p1p2*a1a2)
	       - p1sq*p2sq) / 27;
  const bigint& r = r1*r2;
  if(info) cout<<"u-poly = [1,0, " << -2*p<< ", "<<-8*r<<", "<<s<<"]\n";

  // Compute the roots of "u-poly",  see if any are integral, and check:
  // The roots are sqrt(z1*w1)+sqrt(z2*w2)+sqrt(z3*w3) with some choice of signs,
  // where the wi are the z-values of the second quartic.  We compute the zi
  // from the roots of the first, then to ensure that conjugates are matched
  // we compute the wi from the zi (via, invisibly, the phi_i).

  bigcomplex* roots1 = q1->getroots();
  bigfloat xa1=I2bigfloat(a1), xa2=3*I2bigfloat(a2);
  bigcomplex rz1=xa1*(roots1[0]+roots1[1]-roots1[2]-roots1[3]);
  bigcomplex rz2=xa1*(roots1[0]-roots1[1]+roots1[2]-roots1[3]);
  bigcomplex rz3=xa1*(roots1[0]-roots1[1]-roots1[2]+roots1[3]);
  bigcomplex t = I2bigfloat(a1*p2-a2*p1), tt=1/I2bigfloat(3*a1);
  bigcomplex rzw1 = rz1*sqrt((xa2*rz1*rz1+t)*tt);
  bigcomplex rzw2 = rz2*sqrt((xa2*rz2*rz2+t)*tt);
  bigcomplex rzw3 = rz3*sqrt((xa2*rz3*rz3+t)*tt);

  bigcomplex u;
  
  for(long k=0; k<4; k++)
    {
      switch(k) {
      case 0:  u=rzw1+rzw2+rzw3; break;
      case 1:  u=rzw1+rzw2-rzw3; break;
      case 2:  u=rzw1-rzw2+rzw3; break;
      case 3:  u=rzw1-rzw2-rzw3; break;
      }
      if(is_approx_zero(imag(u)))
	{
	  bigint uu = Iround(real(u));
	  bigint uu2=uu*uu;
	  bigint fuu1 = (uu2-2*p)*uu2+s;
	  bigint fuu2 = 8*r*uu;
	  if(fuu1==fuu2) 
	    {
	      if(info) cout<<"Root u = "<<uu<<endl;
	      return 1;
	    }
	  if(fuu1==-fuu2) 
	    {
	      if(info) cout<<"Root u = "<<-uu<<endl;
	      return 1;
	    }
	}
    }

  if(info) cout<<"No integral roots"<<endl;
  return 0;
}

//#else // not using new_equiv so no need to compile this stuff

int testd(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& d, const bigint& e, const bigint& as, 
	  const bigint& bs, const bigint& cs, const bigint& ds, 
	  const bigint& es, const bigint& dd, const bigint& al, 
	  const bigint& be, const bigint& ga, const bigint& de, 
	  int info);

bigcomplex crossratio(const bigcomplex& x1,const bigcomplex& x2,const bigcomplex& x3,const bigcomplex& x4);

int rootsequiv(const quartic* q1, const quartic* q2, int i, const vector<bigint>& dlist, int info);

int allperms[24][4] =
                 {{0,1,2,3},{1,0,2,3},{0,1,3,2},{1,0,3,2},  // Up to here for Type III
		  {2,3,0,1},{2,3,1,0},{3,2,0,1},{3,2,1,0},  // Up to here for Type I
		  {0,2,1,3},{0,2,3,1},{0,3,1,2},{0,3,2,1},
                  {1,2,0,3},{1,2,3,0},{1,3,0,2},{1,3,2,0},
                  {2,0,1,3},{2,0,3,1},{2,1,0,3},{2,1,3,0},
                  {3,0,1,2},{3,0,2,1},{3,1,0,2},{3,1,2,0},}; // All for Type II


int testd(const bigint& a, const bigint& b, const bigint& c, 
	  const bigint& d, const bigint& e, const bigint& as, 
	  const bigint& bs, const bigint& cs, const bigint& ds, 
	  const bigint& es, const bigint& dd, const bigint& al, 
	  const bigint& be, const bigint& ga, const bigint& de, 
	  int info)
{
 bigint d2 = dd*dd;
 bigint al2=al*al, al3=al2*al, al4=al3*al;
 bigint ga2=ga*ga, ga3=ga2*ga, ga4=ga3*ga;
 bigint temp = ga4*es + al*ga3*ds + al2*ga2*cs + al3*ga*bs + al4*as - d2*a;
 if (!is_zero(temp)) return 0;
 bigint de2=de*de, de3=de2*de, de4=de3*de;
 bigint be2=be*be, be3=be2*be, be4=be3*be;
 temp = de4*es + be*de3*ds + be2*de2*cs + be3*de*bs + be4*as - d2*e;
 if (!is_zero(temp)) return 0;
 temp = 4*ga3*de*es + (3*al*ga2*de+be*ga3)*ds
                    + 2*(al2*ga*de+al*be*ga2) * cs
                    + (3*al2*be*ga+al3*de)*bs + 4*al3*be*as - d2*b;
 if (!is_zero(temp)) return 0;
 temp = 4*ga*de3*es + (3*be*ga*de2+al*de3)*ds
                   + 2*(be2*ga*de+al*be*de2)*cs
                   + (be3*ga+ 3*al*be2*de)*bs + 4*al*be3*as - d2*d;
 if (!is_zero(temp)) return 0;
 temp = 6*ga2*de2*es + 3*(be*ga2*de+al*ga*de2) * ds
                + (be2*ga2+ 4*al*be*ga*de+al2*de2) * cs
                + 3*(al*be2*ga+al2*be*de) * bs + 6*al2*be2*as - d2*c;
 if (!is_zero(temp)) return 0;
 return 1;
}  /* end of testd() */


bigcomplex crossratio(const bigcomplex& x1,const bigcomplex& x2,const bigcomplex& x3,const bigcomplex& x4)
{
 return ((x1-x3)*(x2-x4))/((x1-x4)*(x2-x3));  
}

int rootsequiv(const quartic* q1, const quartic* q2, int i, const vector<bigint>& dlist, int info)
{  int ans=0;
   bigcomplex *x = q1->getroots();
   bigcomplex *y = q2->getroots();
   bigcomplex x1=x[0], x2=x[1], x3=x[2], x4=x[3];
   bigcomplex y1=y[allperms[i][0]],y2=y[allperms[i][1]],
           y3=y[allperms[i][2]],y4=y[allperms[i][3]];
//cout<<"X-roots ("<<x<<"): "<<x1<<x2<<x3<<x4<<endl;
//cout<<"Y-roots ("<<y<<"): "<<y1<<y2<<y3<<y4<<endl;
   bigcomplex xtheta = crossratio(x1,x2,x3,x4);   
   bigcomplex ytheta = crossratio(y1,y2,y3,y4);
   if (abs(xtheta-ytheta)>0.1)
     {if (info)
        {cout << i+1 << ": Cross-ratios unequal: "<<xtheta<<" and "<<ytheta<<".\n";
       }
      return 0;
    }
   if(info) cout << i+1 << ": Cross-ratios equal: "<<xtheta<<" and "<<ytheta<<".\n";
   bigcomplex xy1=x1*y1, xy2=x2*y2,
   xy3=x3*y3, // xy4=x4*y4,
   xy23=(x2*y3)-(x3*y2),
   xy13=(y1*x3)-(x1*y3),
   xy12=(x1*y2)-(y1*x2);
   bigcomplex calpha = (xy1*(y2-y3)) + (xy2*(y3-y1)) + (xy3*(y1-y2)),
           cbeta  = (xy1*xy23) + (xy2*xy13) + (xy3*xy12),
           cgamma = xy13 + xy12 + xy23,
           cdelta = (xy1*(x2-x3)) + (xy2*(x3-x1)) + (xy3*(x1-x2));
   bigcomplex scale = calpha;
   if (abs(scale)<abs(cbeta)) scale=cbeta;
   if (abs(scale)<abs(cgamma)) scale=cgamma;
   if (abs(scale)<abs(cdelta)) scale=cdelta;
   if(is_small(scale)) 
     {
       cout << "Warning from rootsequiv(): scale = " << scale << endl;
       cout << "alpha, beta, gamma, delta = " << calpha << cbeta << cgamma << cdelta << endl;
     }
   calpha /= scale;
   cbeta  /= scale;
   cgamma /= scale;
   cdelta /= scale;
   if (!(is_real(calpha) && is_real(cbeta) 
         && is_real(cgamma) 
         && is_real(cdelta)))
     {if (info) 
	{cout << "Transformation not real.\n";
	 cout << "alpha, beta, gamma, delta = " << calpha << cbeta << cgamma << cdelta << endl;
       }
      return 0;
    }
   bigfloat alpha = real(calpha),
          beta  = real(cbeta),
          gamma = real(cgamma),
          delta = real(cdelta);
   bigfloat det = alpha*delta-beta*gamma;

   if (info) 
     {
       cout << "Real transformation has alpha, beta, gamam, delta = ";
       cout <<alpha<<" "<<beta<<" "<<gamma<<" "<<delta<<endl;
       cout << "Testing divisors of "<<q1->getdisc()<<":\n";
//       cout << dlist << endl;
     }
   bigint d;
   vector<bigint>::const_iterator dvar;
   for (dvar=dlist.begin(); dvar!=dlist.end() && (!ans); dvar++)
     {
       d = *dvar;
       bigfloat rscale = sqrt(abs(I2bigfloat(d)/det));
       bigfloat rscaler = floor(rscale+0.5);
       if(abs(rscale-rscaler)<0.001)
	 {
	   bigint al = Iround(rscale*alpha),
	   be = Iround(rscale*beta ),
	   ga = Iround(rscale*gamma),
	   de = Iround(rscale*delta);
	   bigint ddet = abs(al*de-be*ga);
	   if(d==ddet)
	     {
	       if (info)
		 {cout << "d = " << d << endl;
		  cout<<"rscale = "<<rscale<<endl;
		  cout << "al,be,ga,de = "<<al<<" "<<be<<" "<<ga<<" "<<de<<endl;
		}
	       ans=testd(q1->geta(),q1->getb(),q1->getcc(),q1->getd(),q1->gete(),
			 q2->geta(),q2->getb(),q2->getcc(),q2->getd(),q2->gete(),
			 d,al,be,ga,de,info);
	     }
	 }
     }
   return ans ;
 }    // of rootsequiv()
 
 
int equiv(const quartic* q1, const quartic* q2, const vector<bigint>& dlist, int info)
{
   bigint iiq1 = q1->getI(), jjq1 = q1->getJ(), discq1 = q1->getdisc();
   bigint iiq2 = q2->getI(), jjq2 = q2->getJ(), discq2 = q2->getdisc();
   int typeq1 = q1->gettype(), typeq2 = q2->gettype();
   if(info)
     {
       cout<<"Checking equivalence of \n"; q1->dump(cout);
       cout<<"and\n"; q2->dump(cout);
     }
   if (iiq1==iiq2 && jjq1==jjq2 && discq1==discq2 && typeq1== typeq2)
     {
       int nperms = (typeq1==1)? 8 : (typeq1==2)? 24 : 4;
       if (info) cout << "Params agree; calling rootsequiv "<<nperms<<" times.\n";
       int i,ans = 0;
       for (i=0; i<nperms && (!ans); i++)
	 {
	   ans=rootsequiv(q1,q2,i,dlist,info);
	 }
       if(info) {if(!ans) cout << "Not "; cout<<"equiv\n";}
       return ans;
     }
   else 
     { 
       if (info) {cout << "equiv failed on first test!\n";
                  cout << "First  has I="<<iiq1<<", J="<<jjq1<<",";
                  cout << " disc="<<discq1<<", type="<<typeq1<<endl;
                  cout << "Second has I="<<iiq2<<", J="<<jjq2<<",";
                  cout << " disc="<<discq2<<", type="<<typeq2<<endl;
                }
       return 0;
     }
 }  // of equiv()


//#endif
