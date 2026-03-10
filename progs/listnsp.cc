// FILE LISTNSP.CC: Read and display precomputed Newspaces (d-dimensional newforms)
/////////////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2026 John Cremona
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

#define AUTOLOOP

#include "eclib/newspace.h"

int main()
{
  cout << "Program listnsp: read and display precomputed newforms of arbitrary dimension."
       << "\n\n (reading from directory " << newspaces_directory() << ")" << endl;
  long N;
#ifdef AUTOLOOP
  long min_n, max_n;
  cerr<<"Enter first and last level: ";
  cin >> min_n >> max_n;
  if (min_n<1) min_n=1;
  for (N=min_n; N<=max_n; N++) {
#else
  while (1) { cerr<<"Enter level: "; cin>>N; cout<<"\n"; if (N<1) break;
#endif
    cout << endl;
    cout << ">>>> Level " << N << endl;
    Newspace NS;
    int verbose = 0;
    NS.input_from_file(N, verbose);
    int nnf = NS.nforms();
    cout << nnf << (nnf==1? " newform" : " newforms") << " at level " << N;
    if (nnf)
      {
        cout << " of dimension";
        if (nnf>1)
          cout << "s " << NS.dimensions();
        else
          cout << " " << NS.dimensions()[0];
      }
    cout << endl;

    int show_aP = 1;     // 1: do display aP
    int show_AL = 1;     // 1: do display AL
    int show_traces = 1; // 0: do not display traces
    NS.display_newforms(show_aP, show_AL, show_traces);
    }     // end of level loop
  cout << endl;
  exit(0);
}   // end of main()
