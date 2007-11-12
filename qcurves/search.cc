//
// search.cc -- to find infinite points up to naive height
//
// same as findinf.cc apart from interface

#include "mwprocs.h"

int main()
{
  //  set_precision(30);
  set_precision("Enter number of decimal places");
  initprimes("PRIMES",0);

  bigfloat reg, ht_limit;
  bigint u,r,s,t;
  int verbose = 1, modopt=0, pp=0, change_flag; 
  long blength, rank;
  cerr<<"verbose (0/1)? ";             cin >>verbose;
  cerr<<"process points (0/1)? ";      cin >>pp;
  Curve E;

  while (1)
    {
      cerr<<"\nInput a curve: ";      
      cin >> E;
      if ( E.isnull() ) break;
      Curvedata C(E);
      if(verbose) cout << "Curve " << E << endl;
      Curvedata C_min = C.minimalize(u,r,s,t);
      change_flag = ((Curve)C_min) != E;
      if(change_flag&&verbose)
	{
	  cout<<"Searching on standard minimal curve "<<(Curve)C_min<<endl;
	  cout<<"(points found will be transferred back at end)"<<endl;
	}
      
      Point P(C);
      Point Q(C_min);

      mw mwbasis(&C_min, verbose, pp);

      cerr<<"enter search limit: ";      
      cin>>ht_limit;
      mwbasis.search(ht_limit, modopt, 0);
      
      if(pp)
	{
	  rank = mwbasis.getrank();
	  cout<<"Rank of points found is "<<rank<<endl;
	  PointArray b = mwbasis.getbasis();
	  for (long i=0; i<rank; i++)
	    { 
	      Q = b[i];
	      bigfloat hQ = height(Q);
	      cout << "\tGenerator "<<(i+1)<<" is "<<Q<<"; ";
	      cout << "height "<<hQ<<endl;
	      if(change_flag)
		{
		  P = shift(Q,&C,u,r,s,t,1);
		  cout<<"\t--maps back to "<<P<<" on curve "<<E<<endl;
		}
	    }
	  if(rank>1) cout<<"\t\tregulator is "<<mwbasis.regulator()<<endl;
	  cout<<endl;
// Pari output:
	  cout<<"[";
	  for(long i=0; i<rank; i++)
	    {
	      if(i) cout<<",";
	      output_pari(cout,b[i]);
	    }
	  cout<<"]\n";
	}
      else // just output the points found
	{
	  rank = mwbasis.getrank();
	  PointArray b = mwbasis.getbasis();
// Pari output:
	  cout<<"[";
	  for(long i=0; i<rank; i++)
	    {
	      if(i) cout<<",";
	      output_pari(cout,b[i]);
	    }
	  cout<<"]\n";
	}
    }
}

//end of file search.cc
