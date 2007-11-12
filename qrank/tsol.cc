#include "mquartic.h"
#include "mequiv.h"
#include "msoluble.h"
#include "samir.h"

int getquartic(quartic& g, int verbose)
{
  bigint a, b, c, d, e;
  Complex roots[4]; int nrr, type;
  
  if(verbose)  cout << "Enter quartic coefficients a,b,c,d,e ?" << endl;
  char ch; cin>>ch;
  if(ch=='(') cin>>a>>ch>>b>>ch>>c>>ch>>d>>ch>>e>>ch;
     else 
     {
       cin.putback(ch);
       cin >> a >> b >> c >> d >> e;
     }
     
  if (sign(a)==0&&sign(b)==0&&sign(c)==0&&sign(d)==0&&sign(e)==0)
     return 0;
    
  bigint I = 12*a*e - 3*b*d + c*c;
  bigint J = 72*a*c*e + 9*b*c*d - 27*a*d*d - 27*b*b*e - 2*pow(c,3);
  bigint disc = 4*pow(I,3)-J*J;
  type=0;
  g=quartic(a,b,c,d,e,(Complex*)roots,type,I,J,disc);
  return 1;
}

int main()
{
  initprimes("PRIMES",1);
  cout.precision(15);
  cin.flags( cin.flags() | ios::dec );  //force decimal input (bug fix)
  
  int verb; cout << "Verbose? "; cin >> verb;
  int method; cout << "Method? (0 for old, 1 for new Samir): "; cin>>method;
  int i, els=1, els1, els2;
  quartic g;

  while(getquartic(g, verb))
    {
      bigint I = g.getI(), J=g.getJ();
      cout << "I = " << I << ", J = " << J << "\n";
      bigint badp;
      bigintArray plist = pdivs(6*g.getdisc());

      cout << "Checking local solublity at primes " << plist << ":\n";
      for (i=0; i<plist.length; i++)
	{
	  if(method)
	    els1=new_qpsoluble(g,plist[i],verb);
	  else
	    els1=qpsoluble(g,plist[i]);

	  if(els1) cout << "Locally soluble at p = "<<plist[i]<<"\n";
	  else cout << "Not locally soluble at p = "<<plist[i]<<"\n";
	  //	  if(els1==els2)
	  //	    cout<<"New test agrees!\n";
	  //	  else
	  //	    cout<<"New test disagrees!\n";
	  els = els&els1;
	}
      if(!els) continue;
      
      cout << "Everywhere locally soluble.\n";

    }
}
