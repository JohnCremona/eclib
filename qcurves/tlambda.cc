// tlambda.cc -- test program for of Silverman's "Lambda_bad" procedures

#include "points.h"
#include "lambda.h"   // for computing Silverman's set Lambda_bad

int main(){
  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);

  Curve E;

  while (cout<<"Input a curve: ", cin>>E, !E.isnull())
  {
    Curvedata C(E, 0);
//    cout<<"Curve "<< C <<endl;
    vector<bigint> badp = getbad_primes(C);
    CurveRed CR(C);
    CR.display(cout);

    long il,nl,iq;
    int verbose=0;
    vector<bigfloat> Lambda_bad = lambda_bad(CR,nl,verbose);

    if(verbose)
      {
	cout << "Lambda_bad has " << nl << " elements: ";
	for(il=0; il<nl; il++)
	  {if(il) cout<<", "; cout<<Lambda_bad[il];}
	cout<<endl;
      }


    Point P(C);
    Point Q(C);

    while 
      (
       cout<<"Input a point on curve in projective coords [X:Y:Z] ([0:0:0] to finish) ",
       cin >> P,
       !(!cin)&& P.isvalid()
       )
        {       
          cout << "Point " << P;
          int ord = order(P); 
          if(ord>0)cout<< " has order " << ord; 
          else     cout<< " has infinite order";
          cout << endl;
          cout << "Local heights:\n";
          bigfloat gh = 0;
	  bigint d = gcd(getZ(P),getX(P));
          vector<bigint> pdivsz=pdivs(d);
          for(iq=0; iq<pdivsz.size(); iq++)
            {
              bigint q = pdivsz[iq];
	      cout << q << ":\t\t"  << flush;
	      bigfloat ph = pheight(P,q);
	      gh+=ph;
	      cout << ph << endl;
            }
	  cout << "Sum so far =\t" << gh << endl;
	  bigfloat lxd = 2*log(I2bigfloat(d));
	  cout << "log(den(x(P))) = " << lxd << endl;
          for(iq=0; iq<badp.size(); iq++)
            {
              bigint pr = badp[iq];
	      if(div(pr,d)) continue;
              cout << pr << ":\t\t"  << flush;
              bigfloat  ph = pheight(P,pr);
              gh+=ph;
              cout << ph << endl;
            }
	  cout << "Sum so far =\t" << gh << endl;
          bigfloat rh = realheight(P);
          cout << "R:\t\t" << rh << endl;
	  gh += rh;
          cout << "Sum of local heights: " << gh << endl;
	  gh = height(P);
          cout << "global height of P:   " << gh << endl;
// Now attempt to reconstruct P from its height and real x-coordinate:

	  bigint XP = getX(P);
	  bigint ZP = getZ(P);
	  bigint dP = gcd(XP,ZP);
	  bigfloat xx = I2bigfloat(XP);
	  bigfloat xz = I2bigfloat(ZP);
	  bigfloat xp = xx/xz;
	  cout << "Attempting to reconstruct P from its global height and \n";
	  cout << "and real x-coordinate " << xp << endl;

	  if(make_point_from_x_and_ht(&C, Lambda_bad, xp, gh, &Q))
	    {
	      cout << "Success! Point = " << Q << endl;
	      cout << "Original point = " << P << endl;
	    }
	}
    bigfloat ht,xp;
    cout << "Enter height of a point: \n";
    cin>>ht;
    while(cout<<"Enter possible x-coord (0 to stop):",cin>>xp, !is_zero(xp))
      {
    if(make_point_from_x_and_ht(&C, Lambda_bad, xp, ht, &Q))
      {
	cout << "Point = " << Q << endl;
      }
    else
      {
	cout << "Failed to construct point\n";
      }
  }
  }
  cout<<endl;
}         // end main()

