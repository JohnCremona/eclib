// isogs.cc:  implementation of class IsogenyClass and related functions
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
 
#include "matrix.h"
#include "isogs.h"
#include "points.h"

#define DEBUG

#ifdef MPFP
#define EPS1 0.0001
#else
#define EPS1 0.1 // slower than 0.001 but doesn't miss so many isogenies 
                 // when in ordinary double precision...
#endif

#define close(x,y,e) (abs((x)-(y))<(e))

vector<CurveRed> lisog(const CurveRed& CR, Cperiods& cp, long ell, int verbose)
//INPUT: a curve of type CurveRed, so we know its invariants & conductor
//       the periods of the curve in a Cperiods type
//       a prime ell (primality not checked)
//OUTPUT: an array of the curves isogenous to the given curve (possibly empty)
{
  if(ell==2) return twoisog(CR,verbose);
  if(ell==3) return threeisog(CR,verbose);
  bigcomplex x,X,y,Y,z;
  bigint a1,a2,a3,a4,a6, b2,b4,b6,b8;
  CR.getai(a1,a2,a3,a4,a6);
  CR.getbi(b2,b4,b6,b8);
  bigfloat ra1=I2bigfloat(a1),ra2=I2bigfloat(a2),ra3=I2bigfloat(a3),
         ra4=I2bigfloat(a4),ra6=I2bigfloat(a6),
         rb2=I2bigfloat(b2),rb4=I2bigfloat(b4),rb6=I2bigfloat(b6); 
  // b8 isn't used
  if (verbose>1)
    cout<<"\nra1: "<<ra1<<"\tra2: "<<ra2<<"\tra3: "<<ra3<<"\tra4: "<<ra4<<"\tra6: "<<ra6<<"\nrb2: "<<rb2<<"\trb4: "<<rb4<<"\trb6: "<<rb6<<endl;
  bigint conductor = getconductor(CR);
  bigcomplex w1, w2;  cp.getwi(w1, w2);
  int type = cp.getwRI(w1,w2);
  long subgroup, nsubgroups;
  if (ell==2)
    { if (type==1) nsubgroups = 1; else nsubgroups = 3; }
  else nsubgroups = 2;
  // Now nsubgroups is the number of subgroups of C/L of index ell 
  // defined over R
  vector<CurveRed> ans;

  // other loop variables:
  bigcomplex t, w, ti, ui;
  bigint ell1, ad4, ad6, lad4, lad6;
  bigfloat rad4, rad6;
  bigfloat xell(to_bigfloat(ell));
  static bigfloat two(to_bigfloat(2));
  static bigfloat three(to_bigfloat(3));
  static bigfloat four(to_bigfloat(4));
  static bigfloat five(to_bigfloat(5));
  static bigfloat six(to_bigfloat(6));
  static bigfloat seven(to_bigfloat(7));
  ell1 = ell;
  bigint ell2 = sqr(ell1);
  bigint ell3 = ell1*ell2;
  bigint ell4 = ell2*ell2;
  bigint ell6 = ell2*ell4;

  for (subgroup = 1; subgroup <= nsubgroups; subgroup++)
  {
    if (ell==2)
      { if (subgroup==1) z = w1/two;
	else if (subgroup==2) z = w2/two;
	else z = (w1 + w2)/two; }
    else
      { if (subgroup==1) z = w1/xell;
	else if (type==1) z = (w1 - w2 - w2)/xell;
	else z = w2/xell; }
    t =  w = to_bigfloat(0);
    if (verbose>1) cout<<"Subgroup: "<<subgroup<<": z = "<<z<<endl;

    long iz, ilim;
    if (ell==2) ilim = 1;
    else ilim = (ell - 1)/2;
    for (iz=1; iz <= ilim; iz++)
      { cp.XY_coords(X, Y, to_bigfloat(iz)*z);
	if(ell==2) Y=to_bigfloat(0);  // fix to avoid q(z)=-q(tau) awkward case
//Must convert from Y^2=4X^3 + ... model to our minimal model:	
	x = X - (ra1*ra1 + 4*ra2)/12;
	y = (Y - ra1*x - ra3)/two;
	if (verbose>2)
	  cout<<"i = "<<iz<<":  i*z = " << to_bigfloat(iz)*z 
	      << "\n x = "<<x<<"\n  y = "<<y<<endl ;
	if (ell==2) ti = three*x*x + two*ra2*x + ra4 - ra1*y;
	else ti = six*x*x + rb2*x + rb4;
	ui = four*x*x*x + rb2*x*x + two*rb4*x + rb6;
	t += ti;
	w += ui + x*ti;
      }
    if (verbose>1)
      cout<<"t = "<<t<<";\tw = "<<w<<endl;

    rad4 = ra4 - five*real(t);
    rad6 = ra6 - rb2*real(t) - seven*real(w);
    if (verbose>1)
      cout<<"new a4: "<<rad4<<"\tnew a6: "<<rad6<<endl;

    ad4 = Iround(rad4);
    ad6 = Iround(rad6);
    if (verbose>1)
      cout<<"bigint values are "<<ad4<<" and "<<ad6<<endl;

    if ( close(I2bigfloat(ad4) , rad4 , EPS1) &&  
	 close(I2bigfloat(ad6) , rad6 , EPS1) )
      {
	Curve newcurve(a1,a2,a3,ad4,ad6);
	if(verbose>0) cout << "Testing curve " << newcurve << endl;
	Curvedata newCD(newcurve,1);
	if ((Curve)CR == (Curve)newCD)
	  continue; 
	// cout<<" ## Warning! curve possibly "<<ell<<"-isogenous to itself"<<endl;
	else
	  {
	    CurveRed newCR(newCD);
	    if (conductor==getconductor(newCR))
	      { 
		ans.push_back(newCR);
		if (verbose>1) 
		  cout<<"new curve:\n"<<(Curve)newCD<<endl;
	      }
	  }
      }
    else
      {
	if(verbose>1) cout << "Not close enough\n";
	rad4 = rad4*I2bigfloat(ell4);
	rad6 = rad6*I2bigfloat(ell6);
	if (verbose>1)
	  cout<<"new a4: "<<rad4<<"\tnew a6: "<<rad6<<endl;
	lad4 = Iround(rad4);
	lad6 = Iround(rad6);
	if (verbose>1)
	  cout<<"bigint values are "<<lad4<<" and "<<lad6<<endl;
	if ( close(I2bigfloat(lad4) , rad4 , EPS1) &&
	     close(I2bigfloat(lad6) , rad6 , EPS1) )
	  { 
	    Curve newcurve(ell*a1,ell2*a2,ell3*a3,lad4,lad6);
	    if(verbose>0) cout << "Testing curve " << newcurve << endl;
	    Curvedata newCD(newcurve,1);
	    if ((Curve)CR == (Curve)newCD)
	      cout<<" ## Warning! curve possibly "<<ell<<"-isogenous to itself"<<endl;
	    else 
	      {
		CurveRed newCR(newCD);
		if (conductor==getconductor(newCR))
		  { 
		    ans.push_back(newCR);
		    if (verbose>1) 
		      cout<<"new curve:\n"<<(Curve)newCD<<endl;
		  }
	      }
	  }
	else if(verbose>1) cout << "Not close enough\n";
      }
  } // ends for-loop

  return ans;
}

int semistable(const CurveRed& CR)
{
  int ans=1;
  vector<bigint> plist = getbad_primes(CR);
  vector<bigint>::iterator pvar = plist.begin();
  while(pvar!=plist.end())
    if(getord_p_N(CR,*pvar++)>1) return 0;
  return ans;
}

int comprat(const bigint& n1, const bigint& d1,
	    const bigint& n2, const bigint& d2)
{
  return n1*d2==n2*d1;
}

vector<long> getelllist(const CurveRed& CR)
{
  static const bigint j11a = BIGINT(-32768);
  static const bigint j11b = BIGINT(-121);
  static const bigint j11c = BIGINT(-24729001);

  static const bigint nj17a = BIGINT(-297756989);
  static const bigint dj17a = BIGINT(2);
  static const bigint nj17b = BIGINT(-882216989);
  static const bigint dj17b = BIGINT(131072);

  static const bigint j37a = BIGINT(-9317);
  static const bigint j37b = atoI("-162677523113838677");

  static const bigint j19 = BIGINT(-884736);
  static const bigint j43 = BIGINT(-884736000);
  static const bigint j67 = atoI("-147197952000");
  static const bigint j163 = atoI("-262537412640768000");

  static const bigint one = BIGINT(1);

  vector<long> ans; ans.reserve(4);
  ans.push_back(2);
  ans.push_back(3);
  ans.push_back(5);
  ans.push_back(7);
  bigint N = getconductor(CR);
  if(!semistable(CR))
    {
      ans.push_back(13);
      bigint njay=pow(getc4(CR),3);
      bigint djay=getdiscr(CR);
      bigint g=gcd(njay,djay);
      if(!is_one(g)) {njay/=g; djay/=g;}
      if(djay<0) {djay=-djay; njay=-njay;} // Thanks to Mark Watkins
      if(is_one(djay))
	{
	  if((njay==j11a)||(njay==j11b)||(njay==j11c)) ans.push_back(11); 
	  else {if((njay==j37a)||(njay==j37b))         ans.push_back(37); 
	  else {if(njay==j19)                          ans.push_back(19); 
	  else {if(njay==j43)                          ans.push_back(43); 
	  else {if(njay==j67)                          ans.push_back(67); 
	  else {if(njay==j163)                         ans.push_back(163); 
	  }}}}}
	}  // end if integral cases
      else // check with j17a, j17b
	{
	  if(div(17,N))
	    if(comprat(njay,djay,nj17a,dj17a)||
	       comprat(njay,djay,nj17b,dj17b)) ans.push_back(17);
	}
    }
  return ans;
}

IsogenyClass::IsogenyClass(const CurveRed& C, int verbose)
{
  verb=verbose;
  cp = Cperiods(C);
  if(verb)
    {
      cout << endl;
      C.output(cout);
      cout<<"\nPeriod lattice:\n" << cp << endl;
    }
  llist = getelllist(C);
  ss = semistable(C);
  if(verb)
    {
      cout << "Curve is ";if(!ss)cout<<"NOT ";cout<<"semistable."<<endl;
    }
  nell = llist.size();
  curves.push_back(C); 
  fromlist.push_back(0); 
  isoglist.push_back(0);
  matij = vector<long>(MAXNCURVES*MAXNCURVES,0); // initialized to 0
}

void IsogenyClass::process(long i)  // process i'th curve
{
  vector<long> lworks(nell); // only used when i=0
  CurveRed thisc = curves[i];
  if (verb) cout << "Working on curve " << i+1 << ": " << (Curve)thisc << endl;
  Cperiods pers(thisc);
  vector<long>::iterator lvar=llist.begin();
  long il=0, ell, n;
  while(lvar!=llist.end())
    {
      ell = *lvar++;
      if (verb) cout << "trying l = " << ell << "..." << flush;

      vector<CurveRed> lisogcurves = lisog(thisc,pers,ell,verb);

      if (verb) cout << lisogcurves.size() << " isogenous curves found." << endl;
      if(i==0) 
	{
// 	  cout<<"setting lworks["<<il<<"] to "<< !lisogcurves.empty() << endl;
// 	  cout<<"where llist["<<il<<"] = "<<llist[il]<<endl;
	  lworks[il++] = !lisogcurves.empty();
	}
      vector<CurveRed>::iterator Ci=lisogcurves.begin();
      n=0;
      while(Ci!=lisogcurves.end())
	{ 
	  CurveRed newc = *Ci++; n++;
	  if (verb) cout << "\t"<<n<<": "<<(Curve)newc<<"\t: ";
	  int j=0, isnew = 1;
	  vector<CurveRed>::iterator oldCi=curves.begin();
	  while(oldCi!=curves.end())
	    {	      
	      if ((Curve)newc==(Curve)(*oldCi++))
		{
		  isnew=0;
		  matset(i,j,ell);
		  matset(j,i,ell);
		}
	      j++;
	    }
	  if (isnew) 
	    {
	      curves.push_back(newc);
	      fromlist.push_back(i);
	      isoglist.push_back(ell);
	      matset(i,ncurves,ell);
	      matset(ncurves,i,ell);
	      ncurves++;
	      if (verb) cout << "new # " << ncurves << endl;
	    }
	  else if (verb) cout << "repeat" << endl;
	}
    }  // end of ell loop
  if(i==0) // reset llist to good l only;
    {
      vector<long> goodllist;
      for(long i=0; i<nell; i++) 
	if(lworks[i]) 
	  goodllist.push_back(llist[i]);
      nell=goodllist.size();
      llist=goodllist;
      if(verb)
	{
	  cout << "Number of useful l is " << nell;
	  if(nell) cout << ": " << llist;
	  cout<<endl;
	}
    }
}

void IsogenyClass::grow(void)       // does the work
{
  if(verb) cout << "Trying l values: " << llist << endl;
// N.B. ncurves will increase as more curves are found
  for (ndone=0, ncurves=1; ndone<ncurves; ndone++) 
    {
//       cout << "After processing "<<ndone<<" curve(s), :";
//       display(cout);
      process(ndone);
    }
}

void IsogenyClass::displaycurves(ostream& os)const
{
  os << endl << ncurves << " curve(s) in the isogeny class"<<endl<<endl;
  if(ncurves==0) return;
  long i;
  for (i=0; i<ncurves; i++)
    { Curve ci = (Curve)curves[i];
      os << (i+1) << ": " << ci;
      if (i>0) os << "  is "<< isoglist[i]<<"-isogenous to curve "<<fromlist[i]+1;
      os<<endl;
    }
  os<<endl;
}

void IsogenyClass::displaymat(ostream& os)const
{
  if(ncurves==0) return;
  long i,j;
  os << "Isogeny matrix:\n";
  os << "\t"; for(j=0; j<ncurves; j++) os<<(j+1)<<"\t";  os<<"\n";
  for(i=0; i<ncurves; i++)
    {
      os<<(i+1)<<"\t"; for(j=0; j<ncurves; j++) os<<mat_entry(i,j)<<"\t";  os<<"\n";
    }
  os<<endl;
}

void IsogenyClass::dumpdata(ostream& os, long rank)  
// output for textab to input
{
  os << ncurves << "\n";
  long ic, jc, np, il, nj;
  char sep = ' ';
  for (ic=0; ic<ncurves; ic++)
    {
      CurveRed& C = curves[ic];
      os << C.a1 << sep << C.a2 << sep << C.a3 << sep << C.a4 << sep 
         << C.a6 << sep ;
      os << rank << sep;
//cout << "C.ntorsion = " << C.ntorsion << endl;
      os << C.get_ntorsion() << sep;
//cout << "After C.get_ntorsion(), C.ntorsion = " << C.ntorsion << endl;
      np = C.the_bad_primes.size();
      if(sign(C.discr)>0) os << "+1" <<sep;
      else  os << "-1" <<sep;
      vector<bigint>::const_iterator pi;
      pi=C.the_bad_primes.begin();
      while(pi!=C.the_bad_primes.end()) os << C.reduct_array[*pi++].ord_p_discr << sep;
      pi=C.the_bad_primes.begin();
      while(pi!=C.the_bad_primes.end()) os << C.reduct_array[*pi++].ord_p_j_denom << sep;
      pi=C.the_bad_primes.begin();
      while(pi!=C.the_bad_primes.end()) os << C.reduct_array[*pi++].c_p << sep;
      pi=C.the_bad_primes.begin();
      while(pi!=C.the_bad_primes.end()) os << C.reduct_array[*pi++].Kcode.code << sep;

      os << nell << sep;
      for(il=0; il<nell; il++)
	{
	  long ell = llist[il];
	  for(jc=0, nj=0; jc<ncurves; jc++) 
	    if(mat_entry(ic,jc)==ell) nj++;
	  os << nj << sep;   // # of ell-isogenous curves
	  os << llist[il] << sep;
	  for(jc=0; jc<ncurves; jc++) 
	    if(mat_entry(ic,jc)==ell) os << (jc+1) << sep;
                    	    // id #s of l-isogenous curves
	}
      os << endl;
    }
}

vector<long> IsogenyClass::getmat() const
{
  vector<long> ans(ncurves*ncurves);
  long i,j;
  for(i=0; i<ncurves; i++)
    for(j=0; j<ncurves; j++)
      ans[i*ncurves+j] = mat_entry(i,j);
  return ans;
}

mat IsogenyClass::getmatrix() const
{
  mat ans(ncurves,ncurves);
  long i,j;
  for(i=0; i<ncurves; i++)
    for(j=0; j<ncurves; j++)
      ans.set(i+1,j+1,mat_entry(i,j));
  return ans;
}


vector<CurveRed> twoisog(const CurveRed& CR, int verbose)
//INPUT: a curve of type CurveRed
//OUTPUT: an array of the curves 2-isogenous to the given curve
//(possibly empty)
{
  if(verbose>1) cout<<"In twoisog with CR = "<<CR<<endl;
  Curvedata CD((Curvedata)CR);
  if(verbose>1) cout<<"In twoisog with CD = "<<CD<<endl;
  vector<Point> tt = two_torsion(CD);  // include [0:1:0]
  vector<CurveRed> ans;
  if(tt.size()==1) return ans;
  bigint a1,a2,a3,a4,a6,b2,b4,b6,b8;
  CR.getai(a1,a2,a3,a4,a6);
  CR.getbi(b2,b4,b6,b8);
  unsigned int i;
  for(i=1; i<tt.size(); i++)
    {
      Point T = tt[i];
      bigint x = (4*getX(T))/getZ(T);  // =4x(T)
      bigint t = 3*x*x+2*b2*x+8*b4;
      if(verbose) cout<<"t = "<<t<<endl;
      bigint w = x*t;
      if(verbose) cout<<"w = "<<w<<endl;
      Curve E(2*a1, 4*a2, 8*a3, 16*a4-5*t, 64*a6-4*b2*t-7*w);
      if(verbose) 
	cout<<"raw 2-isogenous curve = "<<E<<endl;
      Curvedata EE(E,1);
      if(verbose) 
	cout<<"after minimising,  2-isogenous curve = "<<(Curve)EE<<endl;
      ans.push_back(CurveRed(EE));
    }
  return ans;
}

vector<CurveRed> threeisog(const CurveRed& CR, int verbose)
//INPUT: a curve of type CurveRed
//OUTPUT: an array of the curves 3-isogenous to the given curve
//(possibly empty)
{
  if(verbose>1) cout<<"In threeisog with CR = "<<CR<<endl;
  Curvedata CD((Curvedata)CR);
  if(verbose>1) cout<<"In threeisog with CD = "<<CD<<endl;
  vector<bigint> xt3 = three_torsion_x(CD);
  if(verbose>1) cout<<"xt3 = "<<xt3<<endl;
  vector<CurveRed> ans;
  if(xt3.size()==0) return ans;
  bigint a1,a2,a3,a4,a6,b2,b4,b6,b8;
  CR.getai(a1,a2,a3,a4,a6);
  CR.getbi(b2,b4,b6,b8);
  for(unsigned int i=0; i<xt3.size(); i++)
    {
      bigint x = xt3[i];  // = 3*x-coord
      bigint t = (2*x+b2)*x+3*b4;
      if(verbose) cout<<"t = "<<t<<endl;
      bigint w = ((10*x+6*b2)*x+27*b4)*x+27*b6;
      if(verbose) cout<<"w = "<<w<<endl;
      Curve E(3*a1, 9*a2, 27*a3, 81*a4-135*t, 729*a6-243*b2*t-189*w);
      if(verbose) 
	cout<<"raw 3-isogenous curve = "<<E<<endl;
      Curvedata EE(E,1);
      if(verbose) 
	cout<<"after minimising,  3-isogenous curve = "<<(Curve)EE<<endl;
      ans.push_back(CurveRed(EE));
    }
  return ans;
}
