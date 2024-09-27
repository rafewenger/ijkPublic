/*!
 *  @file ijkdual_constants.h
 *  @brief Constants used by ijkdual.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2018-2023 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#ifndef _IJKDUAL_CONSTANTS_H_
#define _IJKDUAL_CONSTANTS_H_

namespace IJKDUAL {

  const int DIM2(2);
  const int DIM3(3);
  const int DIM4(4);

  const int NUM_VERT_PER_TRIANGLE(3);
  const int NUM_VERT_PER_QUADRILATERAL(4);
  const int NUM_VERT_PER_PENTAGON(5);
  const int NUM_VERT_PER_HEXAGON(6);
  const int NUM_VERT_PER_SEPTAGON(7);
  const int NUM_VERT_PER_OCTAGON(8);

  const int TRI_ENCODING_BIT_SET_SIZE(16);
  const int INTERSECTS_FACET_BIT_SET_SIZE(16);
}

#endif
