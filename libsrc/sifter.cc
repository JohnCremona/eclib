// sifter.cc: implementation of class for sifting E(Q)/2E(Q)
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
 
// NB This is used for proving that points are independent; now
// largely obsolete, being superceded by general saturation algorithms

#include "points.h"
#include "sifter.h"

void sifter::init()
{
  long iaux, i, j, nr;

  auxs = new long[num_aux];
  nroots = new int[num_aux];
  thetamod = new long*[num_aux];
  squares = new int*[num_aux];
  all_p = new long[2*num_aux];

  for(i=0; i<num_aux; i++) thetamod[i] = new long[3];
  iaux=0;
  max_dim_im=0;

  // the rest of the auxs must be chosen as follows: they should
  // be good primes p>5, such that the resolvent cubic has root(s) mod p.

  primevar pr; pr++; pr++;  // skip past 2 and 3

  for(;pr.ok()&&iaux<num_aux; pr++)
    {
      long p = pr;
      if(div(p,disc)) continue;
      if(verbose>1) cout<<"Trying p = " << p << endl;
      long c1 = mod(-27*I,p);
      long c2 = mod(-27*J,p);
      nr = nrootscubic(0,c1,c2,p,thetamod[iaux]);
      if(verbose>1) cout<<"nr = " << nr << endl;
      if(nr>0) 
	{
	  auxs[iaux]=p;
	  nroots[iaux]=nr;
	  iaux++;
	  all_p[max_dim_im++]=p; 
	  if(nr>1) all_p[max_dim_im++]=p; //again, since p gives 2 bits
	}
    }

  pivcols = new int[max_dim_im];
  eps_mat  = new int*[max_dim_im];
  for(i=0; i<max_dim_im; i++)
    {
      eps_mat[i] = new int[max_dim_im];
    }

  // report on which primes will be used:

  if((verbose>1)&&(num_aux>0)) 
    {
      cout<<"sifting using " <<num_aux<<" moduli: \n";
      cout<<"p:\t";
      for(j=0; j<num_aux; j++) cout<<auxs[j]<<"\t";
      cout<<"\n";
      cout<<"nroots:\t";
      for(j=0; j<num_aux; j++) cout<<nroots[j]<<"\t";
      cout<<"\n";
      cout<<"theta1:\t";
      for(j=0; j<num_aux; j++) cout<<thetamod[j][0]<<"\t";
      cout<<"\n";
      cout<<"theta2:\t";
      for(j=0; j<num_aux; j++) 
	if(nroots[j]==1) cout<<"*\t";
	else 
	  cout<<thetamod[j][1]<<"\t";
      cout<<"\n";
      cout<<"theta3:\t";
      for(j=0; j<num_aux; j++) 
	if(nroots[j]==1) cout<<"*\t";
	else 
	  cout<<thetamod[j][2]<<"\t";
      cout<<"\n";
    }

  // initialize flag arrays for squares:

  for (i = 0; i < num_aux; i++)
    {
      long aux = auxs[i];
      long half_aux = ((aux + 1) / 2);
      squares[i] = new int[aux];
      for (j = 0; j < aux; j++)      squares[i][j]=0;
      for (j = 0; j < half_aux; j++) squares[i][posmod( j*j, aux )]=1;

    } // end of aux loop

  if((verbose>1)&&(num_aux>0)) 
    cout<<"finished sifter::init()"<<endl;
}

void sifter::clear()
{
  long i;
  for(i=0; i<max_dim_im; i++) delete[]eps_mat[i];
  delete[]eps_mat;
  for(i=0; i<num_aux; i++) delete[]thetamod[i];
  delete[]thetamod;
  for(i=0; i<num_aux; i++) delete[]squares[i];
  delete[]squares;
  delete[]auxs;
  delete[]nroots;
  delete[]all_p;
  delete[]pivcols;
}

void sifter::vecout(int* v)
{
  int i,j=0, first=1;
  for(i=0; i<max_dim_im; i++) 
    {
      cout << v[i];
      if(nroots[j]==1) {j++; cout<<" ";}
      else {if(!first) {j++; cout<<" ";} first=!first;}
    }
  cout<<endl;
}

int* sifter::eps(const bigint& x, const bigint& z2)
{
  int *ans = new int[max_dim_im];  // caller is reposnsible for deletion
  int i, c, *ansi=ans;

  for(i=0; i<num_aux; i++)
    {
      c = code(x,z2,i);
      if(nroots[i]==1) 
	*ansi++ = (c&1);
      else 
	{
	  *ansi++ = (c&1); 
	  c>>=1;
	  *ansi++ = (c&1); 
	}
    }
  return ans;
}

void sifter::process(const vector<Point>& Plist) 
{
  vector<Point>::const_iterator P=Plist.begin();
  while(P!=Plist.end()) 
    {
      if(verbose) cout<<"Processing point "<<*P<<endl;
      process(*P++);
    }
}

//#define DEBUG

void sifter::process(const Point& P) 
{
  bigint x0,y0,z0;  P.getcoordinates(x0,y0,z0); //P=[x0:y0:z0] on E
  bigint z=gcd(x0,z0); x0/=z;                // =[x0/z2,y0/z0] now, z0=z^3
  bigint z2=z*z;
  bigint x = (36)*x0+r*z2;
#ifdef DEBUG
  bigint y = (216)*y0+(36)*s*x0*z+t*z0;   // (x/z^2,y/z^3) is on E_{I,J}
  bigint check = y*y-(x*x*x-27*I*x*z2*z2-27*J*z0*z0);
  if(is_zero(check))
    {
      if(verbose) 
	cout<<"Transformed P (on I,J curve) has (x,y,z) = ("<<x<<","<<y<<","<<z<<")\n";
    }
  else
    {
      cout<<"Error in transforming P to I,J curve!\n";
      cout<<"Transformed P has (x,y,z) = ("<<x<<","<<y<<","<<z<<")\n";
    }
#endif

  int i,j;
  int * image = eps(x,z2);
  int *pivrow;

  if(verbose)
    {
      cout << "Image =           \t";
      vecout(image);
    }

  for(i=0; i<rank; i++)
    {
      if(image[pivcols[i]])
	{
	  pivrow = eps_mat[i];
	  for(j=0; j<max_dim_im; j++)
	    {
	      image[j] = image[j] ^ pivrow[j];
	    }
	}
    }

  if(verbose)
    {
      cout << "After elimination:\t";
      vecout(image);
    }

  int newpiv=-1;
  for(j=0; (j<max_dim_im)&&(newpiv<0); j++)
    {
      if(image[j]) newpiv=j;
    }
  if(newpiv<0) // this point is dependent
    {
      if(verbose)
	cout << "eps(P) dependent on previous points!\n"; 
    }
  else
    {
      for(j=0; j<max_dim_im; j++) 
	eps_mat[rank][j]=image[j];
      pivcols[rank++] = newpiv;
      if(verbose)
	{
	  cout << "P independent of previous points (using prime "
	       <<all_p[newpiv]<<")\n";
	  cout << "rank increases to "<<rank<<endl;
	}
    }
  delete[]image;
}

int sifter::code(const bigint& x, const bigint& z2, int i)
{
  long p = auxs[i], theta, alpha, beta;
  int j, ans; int eps[3];
  switch(nroots[i]) {
  case 1: 
    theta=thetamod[i][0];
    alpha = posmod(x - theta*z2 , p);
    if(alpha) return !(squares[i][alpha]);  // 1 or 0
    beta = posmod((3*x*x-27*I*z2*z2) , p);
    return !(squares[i][beta]);
    break;
  case 3:
    for(j=0; j<3; j++)
      {
	theta = thetamod[i][j];
	alpha = posmod(x - theta*z2 , p);
	eps[j] = 2*(squares[i][alpha])-(alpha==0)-1; // =0,-1,+1
      }
    if(eps[0]==0) {eps[0]=eps[1]*eps[2];}
    else {
      if(eps[1]==0) {eps[1]=eps[0]*eps[2];}
      else {
	if(eps[2]==0) {eps[2]=eps[0]*eps[1];}
      }
    }
    if(eps[0]==1) ans=(eps[1]==1? 0: 1);
    else          ans=(eps[1]==1? 2: 3);  // binary coding for 00,01,10,11
    return ans;
    break;
  case 0: default: return 0; break;
  }
}

//end of file sifter.cc





