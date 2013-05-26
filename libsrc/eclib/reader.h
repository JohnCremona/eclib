// reader.h: class for reading curves from file/tty
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2012 John Cremona
// 
// This file is part of the eclib package.
// 
// eclib is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// eclib is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with eclib; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////
 
#include <fstream>
#include <cstring> // for gcc >= 4.3
#define TTY "/dev/tty"

class CurveReader {
private:
  int tty_input;
  ifstream in;
public:
  CurveReader()
  {  
    tty_input=0;
    string filename;
    cerr << "Enter a filename for curve input (or tty): ";
    cin >> filename;
    cerr<<"filename entered is "<<filename<<endl;
    tty_input = (filename==string("tty"));
    if(tty_input)
      cerr<<"Taking input from terminal"<<endl;
    else
      {
	cerr<<"using filename "<<filename<<endl;
	in.open(filename.c_str());
	if(!in) 
	  {
	    cerr<<"Failed to open input file "<<filename<<endl;
	    abort();
	  }
      }
  }
  ~CurveReader() {if(!tty_input) in.close();}
  int operator>>(Curve& c)
  {
    if(tty_input)
      {
	cerr<<"Enter a curve (null to exit): "<<flush;
	cin>>ws;  if(cin.eof()) {cerr<<endl; return 0;}
	cin>>c;
      }
    else 
      {
	in>>ws;  if(in.eof()) return 0;
	in >> c;
      }
    return !c.isnull();
  }
};
