/*!
 *  @file ijkMCtable_xitIO.h
 *  @brief I/O routines for .xit (XML Isosurface Table) file.
 *  - .xit is an xml format
 *  - Version 0.6
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2001-2024 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IJKXIO_
#define _IJKXIO_

#include <iostream>

#include "ijkMCtable.h"

namespace IJKXIO {

  // ***************************************************************
  // TYPES
  // ***************************************************************

  typedef enum { XIT_VERSION_1_0, XIT_VERSION_2_0,
                 XIT_VERSION_2_x,
                 UNKNOWN_XIT_VERSION } XIT_VERSION_TYPE;


  // ***************************************************************
  // READ .xit FILE
  // ***************************************************************

  /// @brief Read .xit file.
  void read_xit
    (std::istream & in, IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);


  // ***************************************************************
  // WRITE .xit FILE
  // ***************************************************************

  /// @brief Write .xit file.
  void write_xit
    (std::ostream & out, const XIT_VERSION_TYPE xit_version,
     const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  /// @brief Write version 2.0 .xit file.
  void write_xit_V2
    (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);


  // ***************************************************************
  // OLD VERSIONS
  // ***************************************************************

  /// @brief Read (old) version 1.0 .xit file.
  void read_xit_V1
    (std::istream & in, IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);

  /// @brief Write old version 1.0 .xit file.
  void write_xit_V1
    (std::ostream & out, const IJKMCUBE_TABLE::ISOSURFACE_TABLE & table);
}

#endif

