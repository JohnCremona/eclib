// testing existence of 2-adic points on cubics/quartics

#include "iostream.h"
#include "mlocsol.h"
#include "twoadic.h"

#define  LOOP
//#define COMPARE_LEMMA7
//#define USE_BIGINTS
//#define SHOW_POLY

//
// Case 1 is A=0, B=0, x=2 (mod 4)
// Case 2 is A=1, B=2, x=1 (mod 4)

int main()
{
  int CASE;
  cout<<"Case? ";
  cin>>CASE;
  cout<<"Case "<<CASE<<endl;

  long mod=4;
  long mod1=mod-1;

  int res1, res2;
  long A, B, max;
  long poly[4];
  bigint mpoly[4];
      // Forcing A=B=0(mod 4) and x=2(mod 4)
  poly[3]=64;  mpoly[3]=64;
  switch(CASE) {
  case 1:  poly[2]=96;  mpoly[2]=96; break;
  case 2:  poly[2]=48;  mpoly[2]=48; break;
  }
  long counter[2];  counter[0]=0;  counter[1]=0;
  long countab[32][4][2];
  int i, j;
  for(i=0; i<32; i++) for(j=0; j<4; j++) 
    {countab[i][j][0]=0;countab[i][j][1]=0;}
  bigint a,b,c,d,e,zero,one,two;  zero=0; one=1; two=2; a=zero;

#ifdef LOOP
  cout<<"Enter max: "; cin>>max;
  for(A=0; A<max; A++)
    for(B=0; B<max; B++)
#else
    while(1)
#endif
    {
#ifndef LOOP
      cout<<"Enter A B :";
      cin >> A >> B;
      cout << A << " "<< B << "\t"<<flush;
#endif
      if((A==0)&&(B==0)&&(CASE==1)) 
#ifdef LOOP	
	continue; 
#else
      abort();
#endif
      switch(CASE)
	{
	case 1:
      // Forcing "A"="B"=0(mod 4) and x=2(mod 4)
      poly[1]=16*(A+3);
      poly[0]=4*(B+2*A+2);
      break;
	case 2:
      // Forcing "A"=1, "B"=2(mod 4) and x=1(mod 4)
      poly[1]=16*(A+1);
      poly[0]=4*(B+A+1);
      break;
	}
#ifdef SHOW_POLY      
      cout<<"("<<poly[3]<<","<<poly[2]<<","<<poly[1]<<","<<poly[0]<<")"<<endl;
#endif
      mpoly[0]=poly[0];
      mpoly[1]=poly[1];
#ifdef USE_BIGINTS
      res1 = try1(mpoly);
#else
      if(CASE==1)
	{      
	  //	  res1 = case1(A,B);
	  a=A; b=B; res1 = case1(a,b);
	  if(res1!=(res2=case1(A,B))) 
	    cout<<"bigint case1 returns "<<res1<<", not "<<res2<<endl;
	  /*
	  long C = (2*A+B)%4;
	  int old_res1 = ((C==2)||(C==3));
	  if(res1!=old_res1)
	    {
	      cout<<"(a,b)=("<<A<<","<<B<<"), BSD gives index bound "
		  <<(1+old_res1)<< " but index = "<<(1+res1)<<endl;
	    }
	  */
	  /*
	  res2 = try1(poly);
	  if(res1!=res2)
	    {
	      cout << A << " "<< B << "\t"<<flush;
	      cout<<"New try1 returns "<<res1<<", not "<<res2<<endl;
	    }
	  */
	}
      else
	{
	  //	  res1 = case2(A,B);
	  a=A; b=B; res1 = case2(a,b);
	  if(res1!=(res2=case2(A,B))) 
	    cout<<"bigint case2 returns "<<res1<<", not "<<res2<<endl;
	  /*  
	  res2 = try1(poly);
	  if(res1!=res2)
	    {
	      cout << A << " "<< B << "\t"<<flush;
	      cout<<"New try1 returns "<<res1<<", not "<<res2<<endl;
	    }
	  */
	}
#endif
#ifndef LOOP
      cout<<"Result from try1() = "<<res1<<"\t"<<endl;
#endif
      counter[res1]++;
      long c4 = (A+1)%3;
      long d32 = (A+B+1)%31;
      countab[d32][c4][res1]++;
#ifdef COMPARE_LEMMA7
      res2 = zpsol(zero,mpoly[3],mpoly[2],mpoly[1],mpoly[0],two,zero,0);
#ifndef LOOP      
      cout<<"Result from zpsol() = "<<res2<<endl;
#endif
      if(res1!=res2)
	cout<<"Results disagree for (A,B)=("<<A<<","<<B<<")"<<endl;
#endif
    }
#ifdef LOOP
  double total = counter[0]+counter[1];
  double r0=counter[0]/total;
  double r1=counter[1]/total;
  cout<<"Result = 0 for "<<counter[0]<<" cases, "<<r0<<" of total\n";
  cout<<"Result = 1 for "<<counter[1]<<" cases, "<<r1<<" of total\n";
  /*
  cout<<"\nd\%32\tc\%4\t#0\t#1\n";
  for(i=0; i<32; i++) for(j=0; j<4; j++) 
    {
      if((countab[i][j][0]>0)&&(countab[i][j][1]>0)) 
	cout<<i<<"\t"<<j<<"\t"<<countab[i][j][0]<<"\t"<<countab[i][j][1]<<endl;
    }
  */
#endif
}
