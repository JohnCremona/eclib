// FILE THILBERT.CC: test of Hilbert symbol functions

#include "marith.h"
#include "quadratic.h"
#include "conic.h"
#include "hilbert.h"

//#define AUTO

int main()
{
  bigint a,b,p,x0,y0,z0; long i;  int resp, res, checkres;

#ifdef AUTO
  long la,lb,abmax;
  cout<<"Enter max for a,b: "; cin>>abmax;
  for(la=-abmax; la<=abmax; la++)
  for(lb=-abmax; lb<=abmax; lb++)
    {
      a=la; b=lb;
      if(a*b==0) continue;
#else
  while(1) {
    cout<<"Enter nonzero a and b: "; cin>>a>>b;
    if(a*b==0) break;
#endif
    res=checkres=0;
    cout<<"("<<a<<","<<b<<")_p\n";
    resp = local_hilbert(a,b,0);
    cout<<0<<"\t"<<resp<<endl;
    res = res|resp;
    checkres = checkres^resp;

    resp = local_hilbert(a,b,2);
    cout<<2<<"\t"<<resp<<endl;
    res = res|resp;
    checkres = checkres^resp;

    bigintArray plist = pdivs(a);
    plist=merge(plist,pdivs(b));
    for (i=0; i<plist.length; i++)
      {
	p=plist[i]; if(p==2) continue;
	resp=local_hilbert(a,b,p);
	cout << p << "\t" << resp << endl;
	res = res|resp;
	checkres = checkres^resp;
      }
    cout<<"\nGlobal symbol = " << res << endl;
    cout<<"Check (should be 0) = " << checkres << endl;

    int gres = global_hilbert(a,b,plist);
    if(res==gres)
      cout<<"--agrees with single call to global_hilbert()\n";
    else
      cout<<"--DISAGREES with single call to global_hilbert()\n";

    quadratic q(a,0,b);
    int oldres = !solve_conic(q,1,x0,y0,z0,4);
    if(oldres==res)
      {
	cout<<"--agrees with solve_conic()\n";
      }
    else cout<<"--DISAGREES with solve_conic()\n";
  }
}
