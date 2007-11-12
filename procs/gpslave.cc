// parislave.cc: class for starting up a "slave" background gp process 
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
 
#ifdef USE_GP_FACTORING

#include <unistd.h>  // for unlink() (not needed on linux)
#include <sstream>
#include "gpslave.h"
#include "marith.h"

//#define DEBUG_PARI_SLAVE

parislave the_pari_slave;  // The one and only instance

// test if we have an executable gp to use...
int do_we_have_pari(string& gp_path)
{
  char* path_to_gp = getenv("PATH_TO_GP");
  if (path_to_gp==NULL) 
    gp_path = string("/usr/local/bin/gp");
  else
    gp_path = string(path_to_gp)+string("gp");
  string comm = string("[ -x ") + gp_path + string(" ]");
#ifdef DEBUG_PARI_SLAVE
  cout<<"In do_we_have_pari, command="<<comm<<endl;
#endif
  return !system(comm.c_str());
}

parislave::parislave()
{
  // test if we have an executable gp to use...
  string gpcommand;
  dummy = !do_we_have_pari(gpcommand);
  if(dummy) return;
  
  strcpy(gpinfilename,"/tmp/gp_in_XXXXXX");
  mkstemp(gpinfilename);
  string comm=string("rm ")+gpinfilename;  system(comm.c_str());
#ifdef DEBUG_PARI_SLAVE
  cout<<"Using file "<<gpinfilename<<" for gp input"<<endl;
#endif

  comm = string("mkfifo -m 600 ")+gpinfilename;
#ifdef DEBUG_PARI_SLAVE
  cout<<"About to run command "<<comm<< " via system() call"<<endl;
#endif
  system(comm.c_str());

  comm = gpcommand+" -q -f < "+gpinfilename;
#ifdef DEBUG_PARI_SLAVE
  cout<<"About to call popen( "<<comm<< " )"<<endl;
#endif
  gp = popen(comm.c_str(),"r");
  gpin.open(gpinfilename);
}

parislave::~parislave()
{
  if(dummy) return;
#ifdef DEBUG_PARI_SLAVE
  cout<< "Shutting down gp process..."<<endl;
#endif
  ofstream gpin(gpinfilename);
  gpin <<"\\q\n"<<flush;  // to shut down the gp process
#ifdef DEBUG_PARI_SLAVE
  cout<< "Deleting temp file "<<gpinfilename<<endl;
#endif
  /*
  string comm("rm ");  
  comm+=gpinfilename;
  system(comm.c_str());
  */
  unlink(gpinfilename);
}

int parislave::is_prime(const bigint& n)
{
  if(dummy)
    {
#ifdef DEBUG_PARI_SLAVE
      cout<<"dummy call to pari_slave::is_prime(): no gp running!"<<endl;
#endif
      return 1;
    }
#ifdef DEBUG_PARI_SLAVE
  cout << "Writing to gp's input stream: " 
       << "print(isprime(" << n << "))\n"<<flush;
#endif
  gpin << "print(isprime(" << n << "))\n"<<flush;
#ifdef DEBUG_PARI_SLAVE
  cout << "Reading from gp's output stream"<<endl;
#endif
  int ans;
  fscanf(gp,"%s",gpoutput);
#ifdef DEBUG_PARI_SLAVE
  cout<<"just read "<<gpoutput<<" from gp output"<<endl;
#endif
  istringstream gpout(gpoutput);
  gpout>>ans;
#ifdef DEBUG_PARI_SLAVE
      cout<<"Read ans="<<ans<<" from gp output"<<endl;
  cout<<"Finished reading from gp output"<<endl;
#endif
  return ans;
}

vector<bigint> parislave::factor(const bigint& n, int proof)
{
  vector<bigint> plist; bigint p;
  if(dummy)
    {
#ifdef DEBUG_PARI_SLAVE
      cout<<"Calling pdivs_trial()"<<endl;
#endif
      return pdivs_trial(n);
    }
#ifdef DEBUG_PARI_SLAVE
  cout << "Writing to gp's input stream: " 
       << "print(Vec(factor(" << n << ")[,1]))\n"<<flush;
#endif
  gpin << "print(Vec(factor(abs(" << n << "))[,1]))\n"<<flush;
#ifdef DEBUG_PARI_SLAVE
  cout << "Reading from gp's output stream"<<endl;
#endif
  int first=1; char c = ',';
  while(c!=']')
    {
      fscanf(gp,"%s",gpoutput);
#ifdef DEBUG_PARI_SLAVE
      cout<<"just read "<<gpoutput<<" from gp output"<<endl;
#endif
      istringstream gpout(gpoutput);
      if(first) gpout>>skipws>>c; // swallow leading "["
      first=0;
      gpout>>p;
#ifdef DEBUG_PARI_SLAVE
      cout<<"Reading p="<<p<<" from gp output"<<endl;
#endif
      plist.push_back(p);
      gpout>>skipws>>c>>skipws; // swallow ",", but it might turn out to be "]"
    }
#ifdef DEBUG_PARI_SLAVE
  cout<<"Finished reading from gp output"<<endl;
  if(proof) cout<<"Proving primality of factors found..."<<flush;
#endif
  if(proof)
    for(vector<bigint>::const_iterator pi=plist.begin(); pi!=plist.end(); pi++)
      {
	p =*pi;
	if(!is_prime(p))
	  {
	    cerr<<"WARNING:  pari's factor() returned p="<<p<<" for which isprim(p) FAILS!! Please report.";
	  }
      }
  return plist;
}
#endif // USE_GP_FACTORING
