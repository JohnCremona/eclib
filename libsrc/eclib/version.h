// version.h: declaration of functions show_version(), eclib_version()
//////////////////////////////////////////////////////////////////////////
//
// Copyright 1990-2023 John Cremona
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
//

#if     !defined(_ECLIB_VERSION_H)
#define _ECLIB_VERSION_H      1       //flags that this file has been included

#include <string>
#include <eclib/templates.h>

// Return the current eclib version as a string, e.g. 'v20210317'

string eclib_version();

// Return the current eclib version as a triple [y,m,d], e.g. string, e.g. [2021,3,17]

vector<int> eclib_date();

// Display the current eclib version and compilation information

void show_version(ostream& os = cout);

// compare current eclib version date with a triple (y,m,d), returning
// +1 if the current version is later (more recent) than (y,m,d), 0 if
// equal and -1 if earlier (older).

int compare_eclib_version(int y, int m, int d);

#endif
