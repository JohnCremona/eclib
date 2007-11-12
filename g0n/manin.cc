// FILE MANIN.CC: implementation of manin class
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2007 John Cremona
// 
// This file is part of the mwrank/g0n package.
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

#include <iomanip>
#include "moddata.h"
#include "symb.h"
#include "cusp.h"
#include "homspace.h"
#include "oldforms.h"
#include "newforms.h"
#include "manin.h"

manin::manin(long n, int useolddata, long ntp, int cuspidalflag, int disp) 
:newforms(n,useolddata, ntp,cuspidalflag,disp)
{
  makeh1(); // In case it wasn't done in the newforms constructor
//  cout << "In manin constructor, after constructing the newforms." << endl;
//  if(disp)newforms::display();
  allproj();  // Compute homspace::projcoord, so projcycle can be used
  if(n1ds)
    {
      //Compute p0, pdotlist, dp0;
      primevar pr; while (::div(pr,modulus)) pr++; 
      p0 = pr;  // First "good" prime
      dp0.resize(n1ds);
      pdotlist.resize(n1ds);
      mvp=h1->projmaninvector(p0);          //length=n1ds but starts at 1!
      easy = 1;
      for (long i=0; i<n1ds; i++)
	{  dp0[i] = 1 + p0 - nflist[i].ap0;
	   pdotlist[i]=mvp[i+1];
	   if (pdotlist[i]==0) easy=0;
	 }

      if(disp){
	cout<<"pdotlist: "<<pdotlist<<endl;
	cout<<"dp0:      "<<dp0<<endl;
      }

      findq(); // Initialize nq, dq, qdotlist, initvec
      if(disp){
	cout<<"q = "<<nq<<"/"<<dq<<endl;
	cout<<"qdotlist: "<<qdotlist<<endl;
      }
    }
  else if(disp) cout<<"manin: No newforms!" << endl;
}

void manin::findq()
{
  //  Look for a rational q for which {q,infinity} is nontrivial
  long i,n,d; int stillok,foundq;
  qdotlist.resize(n1ds);
//  iqdotlist.resize(n1ds);
  for (d = 2, foundq=0; !foundq; d++)
    { if (gcd(d,modulus)==1)
	{ for (n = 1; (n<=d/2) && !foundq; n++)
	    { if (gcd(n,d)==1)
		{   // found a candidate n/d
		  initvec=h1->projcycle(n,d);  //starts at 1
		  for (i=0, stillok=1; (i<n1ds) && stillok; i++)
		    { qdotlist[i]= pdotlist[i]-dp0[i]*initvec[i+1];
		      stillok = (qdotlist[i]!=0);
		    }
		  foundq = stillok;
		  if (foundq) 
		    { 
		      nq=n; dq=d; 
//		      for(i=0; i<n1ds; i++) 
//			iqdotlist[i]=invmod(qdotlist[i],MODULUS);
		    }
		}
	    }
	}
    }
}

void show(long ap) {cout<<setw(4)<<ap<<" ";}

#ifdef BINARY
void putout(ofstream& of, short a)
{
  of.write((char*)&a,sizeof(short));
}
void putout(ofstream& of, int a)
{
  short temp = a;
  of.write((char*)&temp,sizeof(short));
}
void putout(ofstream& of, long a)
{
  short temp = a;
  of.write((char*)&temp,sizeof(short));
}
void nl(ofstream& of) {;}
#else
void putout(ofstream& of, long a)
{
  of<<setw(5)<<a;
}
void nl(ofstream& of) {of<<"\n";}
#endif

void manin::getoneap(long p, int output, ofstream& out, int verbose)
{
  vec v;
  long i,ap; int good = (modulus%p!=0);
  if(verbose) {cout << "p = " << setw(5) << p << ": ";}
  if ( easy)  //   Symbol {0,infinity} non-trivial in all cases
    { 
      if (!good)
	{ 
	  for (i=0; i<n1ds; i++)
	    { 
	      int ip = find(plist.begin(),plist.end(),p)-plist.begin();
	      ap = nflist[i].aplist[ip];
	      if(verbose) show(ap);
	      if(output)  putout(out,ap);
	    }
	  if(verbose) cout << "   *** bad prime\n";
	  if(output) nl(out);
	}
      else
	{  // good prime
	  long maxap=(long)(2*sqrt((double)p));
	  v = h1->projmaninvector(p);   //starts at 1
	  for (i=0; i<n1ds; i++)
	    { 
//	      ap = 1+p-mod(mod(v[i+1]*dp0[i],MODULUS)*ipdotlist[i],MODULUS);
	      ap = 1+p-((v[i+1]*dp0[i])/pdotlist[i]);
	      if((ap>maxap)||(-ap>maxap))
		{
		  cout<<"Error:  eigenvalue "<<ap<<" for p="<<p<<" for form # "<<(i+1)<<" is outside valid range "<<-maxap<<"..."<<maxap<<endl;
		  abort();
		}
	      if(verbose) show(ap);
	      if(output)  putout(out,ap);
	    }
	  if(verbose) cout << endl;
	  if(output) nl(out);
	}
    }        // end of "if easy then..."
  else 
    { 
      if (!good)
	{ 
	  for (i=0; i<n1ds; i++)
	    { 
	      int ip = find(plist.begin(),plist.end(),p)-plist.begin();
	      ap = nflist[i].aqlist[ip];
	      if(verbose) show(ap);
	      if(output) putout(out,ap);
	    }
	  if(verbose)cout << "   *** bad prime\n";
	  if(output) nl(out);
	}
      else
	{
	  long maxap=(long)(2*sqrt((double)p));
	  v = h1->newhecke(p,nq,dq) - (p+1)*initvec;  //starts at 1
	  for (i=0; i<n1ds; i++)
	    { 
//	      ap = 1+p- mod(mod(v[i+1] * dp0[i],MODULUS) * iqdotlist[i] , MODULUS );
	      ap = 1+p- ( (v[i+1]*dp0[i]) / qdotlist[i] );
	      if((ap>maxap)||(-ap>maxap))
		{
		  cout<<"Error:  eigenvalue "<<ap<<" for p="<<p<<" for form # "<<(i+1)<<" is outside valid range "<<-maxap<<"..."<<maxap<<endl;
		  abort();
		}
	      if(verbose) show(ap);
	      if(output) putout(out,ap);
	    }
	  if(verbose)cout << endl;
	  if(output) nl(out);
	}
    }      // end of "if easy then ... else ..."
}    

void manin::getap(long first, long last, int output, char* eigfile, int verbose)
{
  ifstream in;  ofstream out;
  if(output) 
    if(first>1)
      {
// open old eigs file
	short file_nforms, file_nforms2, file_nap, temp;
	in.open(eigfile);
// read in numbers of newforms and ap
#ifdef BINARY
	in.read((char*)&temp,sizeof(short)); file_nforms=temp;
	in.read((char*)&temp,sizeof(short)); file_nforms2=temp;
	in.read((char*)&temp,sizeof(short)); file_nap=temp;
#else
	in >> file_nforms >> file_nforms2 >> file_nap;
#endif
	if(first-file_nap!=1)
	  {
	    cout<<"Problem in getap: existing eigs file has "<<file_nap;
	    cout<<" eigs but we are starting at "<<first<<"!"<<endl;
	    cout<<"Aborting"<<endl;
	    return;
	  }
// read in all the old ap
	long ntotal = file_nap*file_nforms;
	short* batch = new short[ntotal];
	short *batchptr=batch;
#ifdef BINARY
	in.read((char*)batch,ntotal*sizeof(short));
#else
	while(ntotal--) in >> (*batchptr++);
#endif
	in.close();
// write them out again with new nap
	out.open(eigfile);  // will delete existing data
	temp = (n1ds?last:0);
#ifdef BINARY
	putout(out,n1ds);
	putout(out,n2ds);
	putout(out,temp);
#else
	out<<n1ds<<" "<<n2ds<<"\n"<<temp<<"\n\n";
#endif
	long ip,jform; 
	batchptr=batch;
	for(ip=0; ip<file_nap; ip++)
	  {
	    for(jform=0; jform<file_nforms; jform++)
	      {
		putout(out,*batchptr++);
	      }
	    nl(out);
	  }
	delete [] batch;
      }
    else 
      {
	out.open(eigfile);
	short temp = (n1ds?last:0);
#ifdef BINARY
	putout(out,n1ds);
	putout(out,n2ds);
	putout(out,temp);
#else
	out<<n1ds<<" "<<n2ds<<"\n"<<temp<<"\n\n";
#endif
      }
  if(n1ds>0)
    {
      if (easy&&verbose)     //   Symbol {0,infinity} non-trivial in all cases
	cout << "L(f,1) nonzero for each form"<<endl; 
      for (primevar pr(last,first); pr.ok(); pr++) 
	getoneap(pr,output,out,verbose);
    }
  else if(verbose>1) cout<<"getap: No newforms!"<<endl;
  if(output) out.close();
}
