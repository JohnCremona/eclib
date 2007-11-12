void rank1::gettype(int t) // new version 31/1/96
{
  type=t;
  long a,absa,b,bb,xc,c,cc;  int a_is_odd;
  bigint a12, aa8, aa, asq, asqI48, acubed, acubedJ64, I48=48*ii, J64=64*jj;
  bigint cub, rsq, r, rem, temp, g, h, x, d, e;
  int extraextra2 = div(64,ii)&&div(128,jj);  // Pascale's extra condition

  if (verbose) cout << "Looking for Type " << t << " quartics:\n";

// Set phi to be the appropriate real root:

  switch(type) {
  case 1:
    phi1 = real(cphi[0]);
    phi2 = real(cphi[1]);
    phi3 = real(cphi[2]);
    orderreal(phi1,phi2,phi3);  // decreasing order
    phi = phi1;
    break;
  case 2:  // N.B. repetition here to allow types 1 and 2 in either order
    phi1 = real(cphi[0]);
    phi2 = real(cphi[1]);
    phi3 = real(cphi[2]);
    orderreal(phi1,phi2,phi3);  // decreasing order
    phi=phi2;
    break;
  case 3:
    if (is_real(cphi[1])) phi=real(cphi[1]);
    else if (is_real(cphi[2])) phi=real(cphi[2]);
    else phi=real(cphi[0]);                  // The real root 
  } // end of switch(type)
#ifdef SHOW_PHI
  if (verbose) cout<<"phi = "<<phi<<endl;
#endif
  const1 = 4*xii-phi*phi;
  
// Set bounds on a loop:

  switch(type) {
  case 1:
    k = const1/ 3.0;  // if(verbose) cout<<"k="<<k<<endl;
    rootk = sqrt(k);
    const4 =  2.0*(k + rootk*phi1);
    const5 =  6.0*rootk +  2.0*phi1;
    amin = 1;
    amax =  ROUNDDOWN((const4 / (const5+ 4.0*phi2)));
    break;
  case 2:
    const6 =  4.0*(xii-phi*phi)/ 3.0;
    amin = ROUNDUP ((const6/( 4*(phi2-phi1))));  // negative
    amax = ROUNDDOWN((const6/( 4*(phi2-phi3))));  // positive
    break;
  case 3:
    double xr =  2.0*sqrt(abs(phi*phi-xii)/ 3.0);
    // The abs is only to avoid rounding errors as phi*phi-ii >= 0 
    amin =   ROUNDUP(((phi-xr)/ 3));
    amax =   ROUNDDOWN(((phi+xr)/ 3));
  } // end of switch(type)

  if (verbose) cout << "Trying a from " << amin << " to " << amax << endl;
  a=amin-1;  // So first value used is amin after incrementing

  long iaux; long *amodi, *hmodi, *auxi;
  int ***flagsi; int **flagai;
  iaux=num_aux; amodi=amod; auxi=auxs;
  while(iaux--) {*amodi++ = posmod(a,*auxi); auxi++;}

  while(a<amax)
    {
      a++;
      for(iaux=0, amodi=amod, auxi=auxs, flagai=flaga, flagsi=flags; 
	  iaux<num_aux; 
	  iaux++, amodi++, auxi++, flagai++, flagsi++) 
	{
	  (*amodi)++; 
	  if((*amodi)==(*auxi)) (*amodi)=0;
	  *flagai = (*flagsi)[*amodi];
	}
//
// Tests: a!=0, not(4|a) if extra2, not(2|a) if extraextra2:
//
      if (a && (! (extra2 && !(a&3))) && (! (extraextra2 && !(a&1)))   )
	{
#ifdef SHOW_ABC_RANGES
	  if(verbose)cout<<"a = "<<a<<endl;
#endif
	  a_is_odd = (a&1);
	  efactor = a_is_odd ? 16 : 8;  // only relevant in extra2 case
	  cfactor = (threediv?3:1)*(extra2?(a_is_odd ? 2 : 4):1);
	  int b_must_be_odd = Jmod2 && !a_is_odd;

	  absa=abs(a);
	  a2=a<<1; a4=a2<<1; a8=a4<<1; a3=a+a2; a6=a3<<1;
	  a12=a8+a4;
	  aa6=a*a6; aa8=a*a8; aa12=2*aa6; aa9=a3*a3;
	  a12phi=12*a*phi; a2phi= a2*phi; a4phi=a4*phi;
#ifdef USE_BIGINTS
	  bigint aa=a;  bigint asq=aa*aa; bigint acubed=asq*aa;  
	  bigint asqI48=asq*I48, acubedJ64=acubed*J64;
#else
	  double asq=a*a; double acubed=a*asq;
	  double asqI48=48*asq*xii, acubedJ64=64*acubed*xjj;
#endif
// Set up extra sieve to ensure that H = -3*square mod 8*a
	  long abs2a = absa<<1, abs8a = absa<<3;;
//	  abs2a*=cfactor;
	  abs8a*=cfactor;
	  int * xflag = new int[abs8a];
	  longlist hlist(abs8a);
	  longlist* blists = new longlist[abs8a];
	  long nhmod8a=0;
	  long ib0=0,ib,ih,bstep=1;
	  if(extra2) bstep=4;
	  if(b_must_be_odd) ib0=1;
	  for(ih=0; ih<abs8a; ih++) {xflag[ih]=0;}
	  for(ib=ib0; ib<=abs2a; ib+=bstep) 
	    {
	      ih = posmod(-3*ib*ib,abs8a);
	      if(!xflag[ih]) 
		{
		  hlist[nhmod8a++]=ih; 
		  blists[ih]=longlist(128);
// 64 is enough for most cases BUT NOT ALL!
		}
	      blists[ih][xflag[ih]++] = ib;
	    }
	  hlist.truncate(nhmod8a);
#ifdef DEBUG_AH
	  for(ih=0; ih<abs8a; ih++) {blists[ih].truncate(xflag[ih]);}
	  cout<<"a = "<<a<<", number of h mod 8a = " << nhmod8a 
	      << ", proportion = " << double(nhmod8a)/abs8a << "\n";
	  cout << "hlist = " << hlist << endl;
	  for(ih=0; ih<abs8a; ih++)
	    if(xflag[ih]>0) 
	      cout << ih<<": "<<xflag[ih]<<": ["<<blists[ih]<<"]\n";
#endif
// Set up H-ranges
	  long long_h, hmin, hmax, hstep, hstep0=1;
	  if(threediv||(a%3==0)) hstep0=3;  // 3|a => 3|H = 8ac-3bb
	  if(extra2) hstep0*=16;    
	  long x,y;
	  long gg = bezout(hstep0,abs8a,x,y);
	  long h0 = hstep0/gg;
	  hstep = h0*abs8a;

	  for(ih=0; ih<nhmod8a; ih++)
	    {
	      if(div(gg,hlist[ih])) 
		hlist[ih] = posmod(h0*hlist[ih]*x,hstep);
	      else 
		{
		  xflag[hlist[ih]]=0; 
		  hlist[ih]=-1;       // code to skip this value;
		}
	    }
#ifdef DEBUG_AH
	  cout << "After adjustment for divisibility by hstep = "<<hstep0
	       <<", hlist = " << hlist << endl;
#endif	  
// Now there are nhmod8a separate A.P.s for h to run through, namely
// h = hlist[i] (mod hstep) (skip if hlist[i]=-1)

	  switch(type) {
	  case 1:
	    hmin = ROUNDUP((a4*phi2));
	    hmax = ROUNDDOWN((a4*phi1));
	    break;
	  case 2:
	    hmin = ROUNDUP((a4*phi2-const6));
	    if(a>0)
	      hmax = ROUNDDOWN((a4*phi3));
	    else
	      hmax = ROUNDDOWN((a4*phi1));
	    break;
	  case 3:
	    hmin = ROUNDUP((12*xii+aa9-a2phi-3*phi*phi));
	    hmax = ROUNDDOWN((a4*phi));
	  } // end of switch(type) on H-ranges
#ifdef SHOW_ABC_RANGES
	  if(verbose)cout<<"hmin = "<<hmin<<", hmax = " << hmax<<endl;
#endif
	  ah_count += (((hmax-hmin)*nhmod8a)/hstep);

	  long nhseq;   // which h A.P. we are in
	  for(nhseq=0; nhseq<nhmod8a; nhseq++)
	{  
	  long h0 = hlist[nhseq];
	  if(h0<0) continue;  // h0=-1 flags to skip this one
	  long hmod8a = posmod(h0,abs8a);
	  longlist blist = blists[hmod8a];
	  long nblist = xflag[hmod8a];
#ifdef DEBUG_AH
	  cout<<"nhseq="<<nhseq<<", h0="<<h0<<", hmod8a="<<hmod8a;
	  cout<<"; nblist="<<nblist<<", blist = ["<<blist<<"]\n";
#endif
	  long_h = h0 + hstep*((hmin-h0)/hstep);
	  long_h-=hstep;         // so first value is correct after increment
	  iaux=num_aux; hmodi=hmod; auxi=auxs;
	  while(iaux--) {*hmodi++ = posmod(long_h,*auxi); auxi++;}

	  while(long_h<=hmax-hstep)
	    {
	      long_h+=hstep;

	      for(iaux=0, hmodi=hmod, auxi=auxs; 
		  iaux<num_aux; 
		  iaux++, hmodi++, auxi++) 
		{
		  (*hmodi) = posmod((*hmodi)+hstep,*auxi);
		}
	      int flagok = nblist;  // Must be >0
	      if(!flagok) {ah_sieve_1++; continue;}
	      for(iaux=0, hmodi=hmod, flagai=flaga; 
		  flagok&&(iaux<num_aux); 
		  iaux++, hmodi++, flagai++)
		flagok = (*flagai)[*hmodi];
	      if(!flagok) {ah_sieve_2++; continue;}

#ifdef DEBUG_AH
	      cout<<"h = "<<long_h<<" ";
#endif
	      ah_sieve_0++;

// Now do the real work: cub should be -27*r^2

#ifdef USE_BIGINTS
	      h = long_h;  // use bigints from now on
	      cub = h*(h*h-asqI48)+acubedJ64;
	      divide(-cub,27,rsq,rem);
	      if(!is_zero(rem)) 
		{
		  cout<<"cub not divisible by 27 for (a,h)=("
		      <<a<<","<<h<<")\n";
		  ah_rfail++;
		  continue;
		}
	      if(!isqrt(rsq,r)) {ah_rfail++; continue;}
#else
	      double xh=long_h;
	      double xcub=xh*(xh*xh-asqI48)+acubedJ64;
	      double xrsq=xcub/(-27);
	      if(xrsq<0) {ah_rfail++; continue;}
	      double xr=sqrt(xrsq);
	      double xxr = abs(xr-floor(xr+0.5));
	      if (xxr>abceps) {ah_rfail++; continue;}
	      r = Iround(xr);
#endif
#ifdef DEBUG_AH
	      cout<<"r = "<<r<<" ";
#endif
// 	     cout<<"a,h,cub = "<<a<<","<<h<<","<<cub;
// 	     bigintArray plistr=pdivs(r);
// 	     cout << "; r = " << r << " has " << plistr.length 
// 		  << " prime divisors: " << plistr << "\t";
// 	     //	     cout << "with exponents: "; 
// 	     //	     for (bigintvar pr(plistr); pr.ok(); pr++) 
// 	     //	       cout << pr.value() <<":"<<val(pr.value(),r) << "\t";
// 	     cout<<endl;
	      
	      for(long k=0; k<nblist; k++)
		{
		  b = blist[k];
		  bb = b*b;
		  b3 = 3*b;
		  c = long_h+3*bb;  // this is 8*a*c so far
		  if(ndiv(a8,c))    // should not happen
		    {
		      cout<<"c not integral for (a,h,b)=("
			  <<a<<","<<h<<","<<b<<")\n";
		      continue;
		    }
		  xc=c>>1;  // = 4*a*c, used below
		  c/=a8;    // = the actual value of c
		  if(!div(cfactor,c)) // shouldn't happen
		    {
		      cout<<"c not divisible by cfactor = "<<cfactor
			<<" for (a,h,b)=("<<a<<","<<h<<","<<b<<")\n";
		      continue;
		    }
		  cc = c*c;
		  temp = b*(bb-xc);
// Loop on  sign of b:
                  int bothb = 1;
                  if((b==0)||(b==abs2a)) bothb=0;
                  for(ib=0; ib<1+bothb; ib++)
		    {if(ib) {b=-b; b3=-b3; temp=-temp;}
//cout<<"\na,b,b3,c,cc,temp = "<<a<<","<<b<<","<<b3<<","<<c<<","<<cc<<","<<temp;
		     divide(r-temp,aa8,d,rem);
//cout<<"; d,rem = "<<d<<","<<rem;
		     if(!is_zero(rem)) {ah_dfail++; continue;}
		     bigint ee=ii+b3*d-cc;
//cout<<"\n ii,b3,d,cc,ee = "<<ii<<","<<b3<<","<<d<<","<<cc<<","<<ee;
		     divide(ee,a12,e,rem);
//cout<<"\n ee,a12,e,rem = "<<ee<<","<<a12<<","<<e<<","<<rem;
		     if(!is_zero(rem)) {ah_rfail++; continue;}
//cout << ":\n [" << a<<","<<b << "," << c << "," << d << "," << e << "]\n"; 
// 
// Now test dividibility conditions in extra2 case:
// (already know 4|b   since 16|h, and not(4|a))
//
		     int skip=0;
		     if(extra2)
		       {
			 skip = !(ndiv(efactor,e)&&ndiv(efactor,a+b+c+d+e));
		       }
		     if(skip) {ah_extra2fail++; continue;}

// Now we have a quartic
// Check the invariants are right (for debugging only):
		     bigint ma(a), mb(b), mc(c);
		     iiabcde = II(ma,mb,mc,d,e);
		     bigint iiabcde = a12*e - b3*d + cc;
		     if ( ii != iiabcde )
		       {
			 cout<<"Error: constructed quartic ";

			 cout << "[" << a<<","<<b << "," << c << "," << d << "," << e << "]"; 
			 cout << " has wrong I-invariant "<<iiabcde<<", not "<<ii<<endl;
			 continue;
		       }
		     jjabcde = JJ(ma,mb,mc,d,e);
		     if (jj != jjabcde)
		       {
			 cout<<"Error: constructed quartic ";
			 
			 cout << "[" << a<<","<<b << "," << c << "," << d << "," << e << "]"; 
			 cout << " has wrong J-invariant "<<jjabcde<<", not "<<jj<<endl;
			 continue;
		       }

// Now construct the root of the quartic and add to the list.
		     bigcomplex za(a), zb((b)), zc((c)), 
		                zd(I2double(d)), ze(I2double(e));
		     bigcomplex* croots =  solvequartic(zb/za,zc/za,zd/za,ze/za);
		     int nrealroots=0, iroot, actual_type;
		     for(iroot=0; iroot<4; iroot++)
		       if(is_real(croots[iroot])) nrealroots++;
		     switch(nrealroots) {
		     case 0: actual_type=1; break;
		     case 2: actual_type=3; break;
		     case 4: actual_type=2; break;
		     default: 
		       cout<<"Quartic [" << a<<","<<b << "," << c << "," << d << "," << e << "] has " << nrealroots << "real roots!\n";
		       for(iroot=0; iroot<4; iroot++) cout<<croots[iroot];
		       cout<<endl;
		     }
		     if (type!=actual_type) // Problem -- wrong type!
		       {
			 cout << "Problem with types: looking for type " 
			   << type << ", but quartic has type " 
			     << actual_type << endl;
			 cout << "Will record it with its actual type.\n";
		       }
		     if(actual_type==3) // make sure two real roots are LAST
		       {
			 bigcomplex tc;
			 for(iroot=0; iroot<2; iroot++)
			   {
			     if(is_real(croots[iroot]))
			       {
				 if(!is_real(croots[2]))
				   {
				     tc=croots[iroot]; 
				     croots[iroot]=croots[2]; 
				     croots[2]=tc;
				   }
				 else
				   if(!is_real(croots[3]))
				     {
				       tc=croots[iroot]; 
				       croots[iroot]=croots[3]; 
				       croots[3]=tc;
				     }
			       }
			   }
		       } // ends if(type==3){}
		     qlist[nquartics].assign(a,b,c,d,e,croots,actual_type,ii,jj,disc);
		     addquartic();
		   }  // end of ib-loop
		} // end of b-loop
	    }  // end of h-loop
	} // end of h-sequence-loop
	  delete[] xflag; delete[] blists;
	} // end of a conditions
    } // end of a loop
  if (verbose) cout << "Finished looking for Type " << t << " quartics.\n";
}  // end of gettype()

