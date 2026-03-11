// FILE MAKENSP.CC:  create newspace and d-dimensional newforms
//////////////////////////////////////////////////////////////////////////
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
#include "eclib/pari_init.h"

const scalar modulus(default_modulus<scalar>());

int main()
{
  cout << "Program makensp: constructing newforms of arbitrary dimension." << endl;
  eclib_pari_init();

  int verbose=1;
  cerr << "Verbose output? (0/1) ";
  cin >> verbose;

  int nap;
  cerr<<"Number of ap? ";
  cin>>nap;

  int output=1;
  cerr << "Output newspaces and newforms to files? (0/1) ";
  cin >> output;

  long N;
#ifdef AUTOLOOP
  long min_n, max_n;
  cerr<<"Enter first and last level: ";
  cin >> min_n >> max_n;
  if (min_n<1) min_n=1;
  for (N=min_n; N<=max_n; N++) {
#else
  while (1) {
    cerr<<"Enter level: "; cin>>N; cout<<"\n"; if (N<1) break;
#endif
    cout << endl;
    cout << ">>>> Level " << N << endl;
    homspace* H1 = get_homspace(N, modulus);
    int cdim = H1->h1cuspdim();
    if (verbose)
      {
        cout << "Constructed level " << N << " homspace of dimension " << H1->h1dim()
             << ", cuspidal dimension " << cdim
             << ", denominator " << H1->h1cdenom() << endl;
      }

    // The Newspace constructor looks for a suitable splitting
    // operator which is a linear combination of up to maxnp prime
    // Hecke operators with coefficients up to maxc:
    int maxnp = 7, maxc = 2;
    Newspace NS(H1, maxnp, maxc, verbose);
    if (!NS.ok())
      {
        cout << "Failed to find a splitting operator using linear combnations of " << maxnp
             << " operators with coefficients up to " << maxc
             << endl;
        continue; // to next level
      }
    int nnf = NS.nforms();
    if (verbose && nnf)
      cout << "Splitting using " << NS.splitopname() << endl;
    cout << nnf << " newform";
    if (nnf!=1) cout << "s";
    cout << " at level " << N;
    if (nnf)
      {
        cout << " of dimension";
        if (nnf>1)
          cout << "s " << NS.dimensions();
        else
          cout << " " << NS.dimensions()[0];
      }
    cout << endl;

    if (!nnf)
      {
        if (output)
          {
            if (verbose)
              cout << "Outputting newspace data for level " << N << endl;
            NS.output_to_file();
          }
        continue; // to next level
      }

    // We compute eigenvalues, coefficients and traces and then resort

    int inf=1;
    for (auto& F: NS.newforms)
      {
        if (verbose)
          cout << "Computing eigenvalues for newform #" << inf <<endl;
        F.compute_eigs_and_coefficients(nap, verbose);
        if (verbose)
          cout << endl;
        inf++;
      }

    // After computing eigenvalues, resort the newforms by their
    // sequence of traces.

    if (verbose)
      cout << "Resorting the newforms using traces" << endl;
    NS.sort_newforms();

    // Display complete newform data for newforms

    cout << endl;

    int show_aP = 1; // do display aP
    int show_AL = 1; // do display AL
    int show_traces = 1; // do display traces
    NS.display_newforms(show_aP, show_AL, show_traces);

    // Output newspace and newform data to files

    if (output)
      {
        if (verbose)
          cout << "Outputting Newspace data for level " << N
               << " to directory " << newspaces_directory() << endl;
        NS.output_to_file();
        if (verbose)
          cout << "Finished outputting Newspace data"<< endl;
      }
  }     // end of level loop
  cout << endl;
  exit(0);
}   // end of main()
