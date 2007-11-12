// compproc.cc: declarations of functions using complex numbers
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
 
#include "compproc.h" 
 
int is_real(const bigcomplex& z) {return(is_approx_zero(imag(z)));}
 
int is_small(bigfloat x) {return is_approx_zero(x);}
int is_small(const bigcomplex& z){return is_approx_zero(z);}

void orderreal(bigfloat& e1, bigfloat& e2, bigfloat& e3)  // puts in decreasing order
{ 
  bigfloat t;
  if      (e1 < e3) {t=e1; e1=e3; e3=t;}      //swap(e1,e3);
  if      (e1 < e2) {t=e1; e1=e2; e2=t;}      //swap(e1,e2);
  else if (e2 < e3) {t=e2; e2=e3; e3=t;}      //swap(e2,e3);
}
 
//#define DEBUG_CAGM

bigcomplex cagm(const bigcomplex& a, const bigcomplex& b)
{
  bigcomplex x=a, y=b, oldx;
#ifdef DEBUG_CAGM
  cout<<"cagm("<<x<<","<<y<<"):"<<endl;
#endif
  static bigfloat two=to_bigfloat(2);
  bigfloat theta, piby2=Pi()/two;
  while (1)
    { 
      oldx=x;
      x=(x+y)/two;
      y= sqrt(oldx*y);
      theta = arg(y/x);
      if ((theta>piby2) || (theta<=-piby2)) y=-y;   
#ifdef DEBUG_CAGM
      cout<<"Relative error = "<<abs((x-y)/x)<<endl;
#endif
      if(is_approx_zero(abs((x-y)/x))) return x;
    }
  return x;
}

bigcomplex normalize(bigcomplex& w1, bigcomplex& w2)
{
   bigcomplex tau = w1/w2, w3;
   if (tau.imag() < 0) { w1=-w1 ; tau=-tau ; }
   w1=w1-w2*round(tau.real());
   tau=w1/w2;
 
   for (int i=1;i<50 && (abs(tau)<1);i++)  
// {Just to stop infinite loop due to rounding}
   { w3=-w1; w1=w2; w2=w3; tau=w1/w2;
    w1=w1-w2*round(tau.real());
    tau=w1/w2;
   }
   return tau;
}
 
void getc4c6(const bigcomplex& w1, const bigcomplex& w2, 
	     bigcomplex& c4, bigcomplex &c6)
{
  bigcomplex tau= w1/w2;
  static bigfloat zero(to_bigfloat(0)), one(to_bigfloat(1)), two(to_bigfloat(2));
  bigfloat pi(Pi());
  bigfloat x = two*pi*tau.real();
  bigfloat y = two*pi*tau.imag();
  bigcomplex q = exp(-y) *  bigcomplex(cos(x),sin(x));
  bigcomplex f = two*pi/w2;
  bigcomplex f2 = f*f;  
  bigcomplex f4=f2*f2;
  
  bigcomplex term   = bigcomplex(one); 
  bigcomplex qpower = bigcomplex(one);
  bigcomplex sum4   = bigcomplex(zero);
  bigcomplex sum6   = bigcomplex(zero);
  
  bigfloat n, n2;

  for (n=1; 
#ifdef MPFP
       !is_approx_zero(term); 
#else
       !is_zero(term); 
#endif
       n+=1)
    {  n2      = n*n;
       qpower *= q;
       term    = n*n2*qpower/(one-qpower);
       sum4   += term;
       term   *= n2;
       sum6   += term;
     }
  c4= (one +  to_bigfloat(240)*sum4)*f4;
  c6= (one -  to_bigfloat(504)*sum6)*f4*f2;
}
 
bigcomplex discriminant(const bigcomplex& b, const bigcomplex& c, const bigcomplex& d)
{
    bigcomplex bb = b*b, cc = c*c, bc = b*c;
    return to_bigfloat(27)*d*d - bc*bc +  to_bigfloat(4)*bb*b*d 
      -  to_bigfloat(18)*bc*d +  to_bigfloat(4)*c*cc;
}
 
#ifdef LIDIA_ALL

vector<bigcomplex> solvecubic(const bigcomplex& c1, const bigcomplex& c2, const bigcomplex& c3)
{
  polynomial<bigcomplex> f; 
  bigcomplex* croots =  new bigcomplex[3];

  f.set_degree(3);  f[3]=1;  f[2]=c1;  f[1]=c2;  f[0]=c3;
  roots(f,croots);  // LiDIA's roots function wants 2nd parameter to be plain array
  vector<bigcomplex> roots(3);
  for(int i=0; i<3; i++) roots[i]=croots[i];
  delete[]  croots;
  return roots;
}

#else

//#define DEBUG_CUBIC

bigcomplex cube_root(const bigcomplex& z)
{
#ifdef DEBUG_CUBIC
  cout << "Taking complex cube root of z = "<<z<<endl;
#endif
  if(is_zero(z)) return z;
  return exp(log(z)/to_bigfloat(3));
}

vector<bigcomplex> solvecubic(const bigcomplex& c1, const bigcomplex& c2, const bigcomplex& c3)
{
#ifdef DEBUG_CUBIC
cout << "In solvecubic with c1 = "<<c1<<", c2 = "<<c2<<", c3 = "<<c3<<"\n";
#endif
  long i, iter, niter=2; // number of iterations in Newton refinement of roots
  bigfloat three(to_bigfloat(3)), two(to_bigfloat(2)), one(to_bigfloat(1));
  bigfloat third = one/three;
   bigcomplex w =  bigcomplex(to_bigfloat(-1), sqrt(three))/two;
   bigcomplex disc = discriminant(c1,c2,c3);
   bigcomplex p3=  three*c2 - c1*c1;
   bigcomplex mc1over3 = -c1*third;
#ifdef DEBUG_CUBIC
cout << "p3 = "<<p3<<", -c1/3 = "<<mc1over3<<"\n";
#endif
   vector<bigcomplex> roots(3);
 
   if (is_zero(abs(disc))) 
     {if (is_zero(abs(p3)))     // triple root
        {roots[0]=roots[1]=roots[2]= mc1over3;
        }
      else           // double root
        {roots[0]=roots[1]= (c1*c2 -  to_bigfloat(9)*c3)/( p3+p3);
         roots[2]=-( roots[0] + roots[0]+c1);
        }
     }
   else              // distinct roots
     {
      bigcomplex q = (((mc1over3+c1)*mc1over3 +c2)*mc1over3 +c3); 
                    // = F(mc1over3);
      if (is_approx_zero(abs(p3)))       // pure cubic
        { 
	  roots[0]=-cube_root(q);
          roots[1]=w*roots[0];
          roots[2]=w*roots[1];
          roots[0]+=mc1over3;
          roots[1]+=mc1over3;
          roots[2]+=mc1over3;
        }
      else
        {bigcomplex d =  to_bigfloat(729)*q*q+ to_bigfloat(4)*p3*p3*p3;
#ifdef DEBUG_CUBIC
cout << "q = " << q << ", p3 = " << p3 << ", d = " << d << "\n";
#endif
         bigcomplex t1cubed = to_bigfloat(0.5)*(sqrt(d)- to_bigfloat(27)*q);
#ifdef DEBUG_CUBIC
	 cout << "t1cubed = " << t1cubed << endl;
#endif
      if (is_approx_zero(abs(t1cubed)))       // approximately pure cubic
	{
#ifdef DEBUG_CUBIC
	  cout << "t1cubed approx 0 so treating as pure cubic " << endl;
#endif
	  roots[0]=-cube_root(q);
          roots[1]=w*roots[0];
          roots[2]=w*roots[1];
          roots[0]+=mc1over3;
          roots[1]+=mc1over3;
          roots[2]+=mc1over3;
	  for(i=0; i<3; i++)
	    {
	      bigcomplex z = roots[i], fz, fdashz;
	      for(iter=0; iter<niter; iter++)
		{
		  fz = ((z+c1)*z+c2)*z+c3;
		  fdashz = (three*z+two*c1)*z+c2;
		  if(!is_zero(fdashz)) z -= fz/fdashz;
		}
	      roots[i] = z;
	    }
	}
      else
	{
	  bigcomplex t1 = cube_root(t1cubed);
	  bigcomplex t2 = t1*w;
	  bigcomplex t3 = t2*w;
#ifdef DEBUG_CUBIC
	  cout<<"t1cubed = "<<t1cubed<<"\n";
	  cout<<"t1 = "<<t1<<", t2 = "<<t2<<", t3 = "<<t3<<"\n";
#endif
	  roots[0] = (-c1+t1-p3/t1)* third;
	  roots[1] = (-c1+t2-p3/t2)* third;
	  roots[2] = (-c1+t3-p3/t3)* third;
	}
        }
    }
#ifdef DEBUG_CUBIC
cout << "refining roots using Newton with " << niter << "iterations\n";
cout << "unrefined roots: ";
for(i=0; i<3;i++) cout << roots[i] << "\n";
#endif
for(i=0; i<3; i++)
  {
    bigcomplex z = roots[i], fz, fdashz;
    for(iter=0; iter<niter; iter++)
      {
	fz = ((z+c1)*z+c2)*z+c3;
	fdashz = (three*z+two*c1)*z+c2;
	if(!is_zero(fdashz)) z -= fz/fdashz;
      }
    roots[i] = z;
  }
#ifdef DEBUG_CUBIC
cout << "refined roots: ";
for(i=0; i<3;i++) cout << roots[i] << "\n";
#endif

  return roots;
}
 
#endif


vector<bigcomplex> solverealquartic(const bigfloat& a, const bigfloat& b, const bigfloat& c, const bigfloat& d, const bigfloat& e)
{
#ifdef DEBUG
  cout<<"In solverealquartic with (a,b,c,d,e)=("<<a<<","<<b<<","<<c<<","<<d<<","<<e<<")\n";
#endif
  bigfloat three = to_bigfloat(3);
  bigfloat ii = 12*a*e - 3*b*d + c*c;
  bigfloat jj =  (72*a*e + 9*b*d - 2*c*c) * c - 27*(a*d*d + b*b*e);
#ifdef DEBUG
  cout<<"ii="<<ii<<"\njj="<<jj<<"\n";
#endif
  bigfloat disc = 4*ii*ii*ii-jj*jj;
  bigfloat  H = 8*a*c - 3*b*b, R = b*b*b + 8*d*a*a - 4*a*b*c;
  bigfloat  Q = H*H-16*a*a*ii;  // = 3*Q really
  int type, nrr;
  if(disc<0) 
    {type=3; nrr=2;}       // 2 real roots
  else 
    {
      if((H<0)&&(Q>0)) 
	{type=2; nrr=4;}   // 4 real roots
      else 
	{type=1; nrr=0;}   // 0 real roots
    }
  bigcomplex c1(to_bigfloat(0)), c2(-3*ii), c3(jj);
#ifdef DEBUG
  cout<<"Type = " << type << " ("<<nrr<<" real roots)\n";
  cout<<"Coeffs of resolvent cubic are:\n"<<c1<<"\n"<<c2<<"\n"<<c3<<endl;
#endif
  vector<bigcomplex> cphi = solvecubic( c1, c2, c3);
  vector<bigcomplex> roots(4);
  bigfloat a4=4*a;
  bigfloat oneover4a = to_bigfloat(1)/a4;
#ifdef DEBUG
  cout<<"Roots of cubic are:\n"<<cphi<<endl;
#endif
  
  if(type<3) 
    {
#ifdef DEBUG
      cout<<"Positive discriminant\n";
#endif
      // all the phi are real;  order them so that a*phi[i] decreases
      bigfloat phi1 = real(cphi[0]);
      bigfloat phi2 = real(cphi[1]);
      bigfloat phi3 = real(cphi[2]);
      if(a>0)      orderreal(phi1,phi2,phi3); 
      else         orderreal(phi3,phi2,phi1); 
#ifdef DEBUG
      cout<<"phi = "<<phi1<<", "<<phi2<<", "<<phi3<<"\n";
#endif

      if(type==2) // all roots are real
	{
#ifdef DEBUG
	  cout<<"Type 2\n";
#endif
	  bigfloat r1 = sqrt((a4*phi1-H)/three);
	  bigfloat r2 = sqrt((a4*phi2-H)/three);
	  bigfloat r3 = sqrt((a4*phi3-H)/three);
	  if(R<0)  r3 = -r3;
#ifdef DEBUG
	  cout<<"r_i = "<<r1<<", "<<r2<<", "<<r3<<"\n";
	  cout<<"product = "<<r1*r2*r3<<", R = "<<R<<endl;
#endif
	  roots[0] = bigcomplex(( r1 + r2 - r3 - b) * oneover4a);
	  roots[1] = bigcomplex(( r1 - r2 + r3 - b) * oneover4a);
	  roots[2] = bigcomplex((-r1 + r2 + r3 - b) * oneover4a);
	  roots[3] = bigcomplex((-r1 - r2 - r3 - b) * oneover4a);
	  // Those are all real and in descending order of size
	}
      else // no roots are real
	{
#ifdef DEBUG
	  cout<<"Type 1\n";
#endif
	  bigfloat r1 = sqrt((a4*phi1-H)/3);
	  bigfloat ir2 = sqrt(-((a4*phi2-H)/3));
	  bigfloat ir3 = sqrt(-((a4*phi3-H)/3));
	  if(R>0)  r1 = -r1;
#ifdef DEBUG
	  cout<<"r_i = "<<r1<<", "<<ir2<<"i, "<<ir3<<"i\n";
	  cout<<"product = "<<-r1*ir2*ir3<<", R = "<<R<<endl;
#endif
	  roots[0] = bigcomplex( r1-b,  ir2 - ir3) * oneover4a;
	  roots[1] = conj(roots[0]); // bigcomplex( r1-b,  ir2 - ir2) * oneover4a;
	  roots[2] = bigcomplex(-r1-b,  ir2 + ir3 ) * oneover4a;
	  roots[3] = conj(roots[2]); // bigcomplex(-r1-b, -ir2 - ir3 ) * oneover4a;
	}
    }
  else // disc < 0
    {
#ifdef DEBUG
      cout<<"Negative discriminant\nType 3\n";
#endif
      bigfloat realphi;       // will hold the real root, which will be cphi[2]
      if (is_real(cphi[1])) 
	{
	  realphi=real(cphi[1]);
	  cphi[1]=cphi[2];
	  cphi[2]=realphi;
	}
      else 
	if (is_real(cphi[2])) 
	  {
	    realphi=real(cphi[2]);
	  }
	else 
	  {
	    realphi=real(cphi[0]);
	    cphi[0]=cphi[2];
	    cphi[2]=realphi;
	  }
#ifdef DEBUG
      cout<<"Sorted roots of cubic (real one last) are \n";
      cout<<cphi[0]<<"\n"<<cphi[1]<<"\n"<<cphi[2]<<endl;
#endif
      bigcomplex r1 = sqrt((a4*cphi[0]-H)/three);
      bigfloat r3   = sqrt((a4*realphi-H)/three);
      if(R<0)  r3 = -r3;
#ifdef DEBUG
      cout<<"r_i = "<<r1<<", "<<conj(r1)<<", "<<r3<<"\n";
      cout<<"product = "<<r1*conj(r1)*r3<<", R = "<<R<<endl;
#endif
      roots[0] = bigcomplex( r3 - b, 2*imag(r1) ) * oneover4a;
      roots[1] = conj(roots[0]);
      roots[2] = bigcomplex(( 2*real(r1) - r3 - b)) * oneover4a;
      roots[3] = bigcomplex((-2*real(r1) - r3 - b)) * oneover4a;
      // roots[2] and roots[3] are real
    }
  return roots;
}

#ifdef LiDIA_ALL

vector<bigcomplex> solvequartic(const bigcomplex& a, const bigcomplex& b, const bigcomplex& c, const bigcomplex& d)
{
  polynomial<bigcomplex> f; 
  bigcomplex* croots =  new bigcomplex[4];

  f.set_degree(4);  f[4]=1;  f[3]=a;  f[2]=b;  f[1]=c;  f[0]=d;
  roots(f,croots);  // LiDIA's roots function wants 2nd parameter to be plain array
  vector<bigcomplex> roots(4);
  for(int i=0; i<4; i++) roots[i]=croots[i];
  delete[]  croots;
  return roots;
}

#else

vector<bigcomplex> solvequartic(const bigcomplex& a, const bigcomplex& b, const bigcomplex& c, const bigcomplex& d)
{ bigcomplex p,q,r,aa,e,f1,f2;
  static bigfloat zero = to_bigfloat(0);
  static bigfloat two = to_bigfloat(2);
  static bigfloat three = to_bigfloat(3);
  static bigfloat four = to_bigfloat(4);
  static bigfloat eight = to_bigfloat(8);
  static bigfloat x16 = to_bigfloat(16);
  static bigfloat x64 = to_bigfloat(64);
  static bigfloat x256 = to_bigfloat(256);
  bigcomplex a4=a/four;
  vector<bigcomplex> roots(4);
  if(is_zero(d))
    { 
      roots[0]= zero;
      vector<bigcomplex> cuberoots=solvecubic(a,b,c);
      roots[1] = cuberoots[0];
      roots[2] = cuberoots[1];
      roots[3] = cuberoots[2];
    }
    else
      {
	p = b -  three*a*a / eight;
	q = ((a / two) * (a*a4 - b)) + c;
	r = ((x256*d) - (x64*a*c) + (x16*a*a*b) - (three*a*a*a*a)) / x256;
	if( is_approx_zero(q) ) 
	  {
	    bigcomplex s = sqrt(p*p - four*r);
	    roots[0] =  sqrt((-p + s) / two) - a4;
	    roots[1] = -sqrt((-p + s) / two) - a4;
	    roots[2] =  sqrt((-p - s) / two) - a4;
	    roots[3] = -sqrt((-p - s) / two) - a4;
	  }
      else
	{
	  vector<bigcomplex> aaroots= 
	    solvecubic(-p / two,-r,((p * r) / two - (q * q) / eight));
	  aa = aaroots[0];
	  if( is_approx_zero(aa) ) aa=zero;
	  e = sqrt(-p + two*aa);
	  f1 = (aa + q / (two*e));
	  f2 = (aa - q / (two*e));
	  bigcomplex s1 = sqrt(e*e -  four*f1);
	  bigcomplex s2 = sqrt(e*e -  four*f2);
	  roots[0] = (( e + s1) /  two) - a4;
	  roots[1] = (( e - s1) /  two) - a4;
	  roots[2] = ((-e + s2) /  two) - a4;
	  roots[3] = ((-e - s2) /  two) - a4;
	}
      }

  long i, iter, niter=2; // number of iterations in Newton refinement of roots
#ifdef DEBUG_QUARTIC
  cout << "refining roots using Newton with " << niter << "iterations\n";
  cout << "unrefined roots: ";
  for(i=0; i<4;i++) cout << roots[i] << "\n";
#endif
  for(i=0; i<4; i++)
    {
      bigcomplex z = roots[i], fz, fdashz;
      for(iter=0; iter<niter; iter++)
	{
	  fz = (((z+a)*z+b)*z+c)*z+d;
	  fdashz = ((four*z+three*a)*z+two*b)*z+c;
	  if(!is_zero(fdashz)) z -= fz/fdashz;
	}
      roots[i] = z;
    }
#ifdef DEBUG_QUARTIC
  cout << "refined roots: ";
  for(i=0; i<4;i++) cout << roots[i] << "\n";
#endif
  return roots;
}

#endif

void quadsolve(const bigfloat& p, const bigfloat& q, 
	       bigcomplex& root1, bigcomplex& root2)
{
  static bigfloat two = to_bigfloat(2);
  static bigfloat four = to_bigfloat(4);
  bigcomplex disc(p*p- four*q);
  bigcomplex rootdisc = sqrt(disc);
  root1 = ( rootdisc-p)/ two;
  root2 = (-rootdisc-p)/ two;
}
 
vector<long> introotscubic(long a, long b, long c, int& nr)
{ bigcomplex za(to_bigfloat(a)), zb(to_bigfloat(b)), zc(to_bigfloat(c));
  vector<bigcomplex> croots =  solvecubic(za,zb,zc);
  vector<long> iroots;
  int i; long x,cx;
  for (i=0; i<3; i++)
    {
      cout << "Complex root = " << croots[i] << endl;
      bigfloat xx = croots[i].real();
      Iasb(x,xx);
      cout << "Rounds to " << x << endl;
      if (x==0) {if (c==0) iroots.push_back(x);}
      else
        {
	  cx = c/x;
	  if (x*cx==c)
            if (((x+a)*x+b+cx) ==  0)
	      iroots.push_back(x);
        }
     }
  return iroots;  
}

