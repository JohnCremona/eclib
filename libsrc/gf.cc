// gf.cc:  interface for NTL's ZZ_p
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
 
#include "eclib/gf.h"

map<ZZ,ZZ_pContext> ZZ_pContextCache;

galois_field::galois_field(void)  //dummy
 :q(to_ZZ(2))
{
  ZZ_pContext_init_or_restore(q);
}

void ZZ_pContext_init_or_restore(const ZZ& q)
{
  auto t = ZZ_pContextCache.find(q);
  if(t==ZZ_pContextCache.end())
    {
      ZZ_p::init(q);
      ZZ_pContext c; c.save();
      ZZ_pContextCache[q]=c;
    }
  else
    {
      t->second.restore();
    }
}

galois_field::galois_field(const ZZ& qq)
 :q(qq)
{
  ZZ_pContext_init_or_restore(qq);
}
