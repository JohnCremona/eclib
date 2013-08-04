// TMANIN_PROF.CC: Program for finding newforms 
#include <eclib/interface.h>
#include <eclib/timer.h>
#include <eclib/moddata.h>
#include <eclib/symb.h>
#include <eclib/cusp.h>
#include <eclib/homspace.h>
#include <eclib/oldforms.h>
#include <eclib/cperiods.h>
#include <eclib/newforms.h>

#define AUTOLOOP
#define LMFDB_ORDER       

int main( int argc, char** argv ) {

  long n, stopp, limit; 
  int output=1, verbose, sign=1, cuspidal=0, noldap=25;

  cout << "Program tmanin_prof.  Using METHOD = " << METHOD << " to find newforms" << endl;
#ifdef MODULAR
  cout << "MODULUS for linear algebra = " << MODULUS << endl;
#endif
  cout << "Verbose output? "; cin>>verbose;
  cout << "How many primes for Hecke eigenvalues? ";
  cin  >> stopp; cout << endl;
  cout << "Output Hecke eigenvalues to file? (0/1) ";  cin >> output;
  cout << "Sign? (-1/0/1) ";  cin >> sign;

  // Start timer
  timer profile("tmanin_runtimes.dat");
  profile.add("createfromscratch");

  cout<<"Enter first and last N: ";cin>>n>>limit;

  if( n <= 0 ) {
    cout << "Invalid level" << endl;
    exit( EXIT_FAILURE );
  }

  for( long level = n; level <= limit; level++ ) {  

      cout << "\n>>>Level " << level;
    
      stringstream l;
      l << "Level " << level << " ";
      profile.write(l.str());

      profile.start();

      if( verbose ) cout << endl; else cout << ":\t";

      newforms nf(level,verbose); 

      profile.start("createfromscratch");
  
      nf.createfromscratch(sign,noldap);
  
      profile.stop("createfromscratch");

#ifdef LMFDB_ORDER
      nf.sort();
      nf.make_projcoord();
#endif
      if( verbose > 1 ) nf.display();
      else cout << nf.n1ds << " newform(s) found."<<endl;

      if( verbose && nf.n1ds > 0 )
        cout << "Computing ap for primes up to " << prime_number(stopp) << endl;

      nf.addap( stopp );

      if( output ) nf.output_to_file();

      profile.stop();

      profile.show(1,"createfromscratch");
      profile.clearAll();

  }

  exit( EXIT_SUCCESS );
}
