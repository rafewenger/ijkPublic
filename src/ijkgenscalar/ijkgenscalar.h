/*!
 *  @file ijkgenscalar.h
 *  @brief Generate a scalar field
 *  -Version v0.2.0
 */

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2018 Rephael Wenger

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

/*!
  \mainpage IJKGENSCALAR: Generate scalar grid.

  IJKGENSCALAR is a program for generating a regular scalar grid
  representing various scalar fields.  It can also be used to generate
  gradient vectors at the grid vertices.
*/

#ifndef _IJKGENSCALAR_
#define _IJKGENSCALAR_

#include <sstream>
#include <string>
#include <vector>

#include "ijk.tpp"
#include "ijkscalar_grid.tpp"
#include "ijkvector_grid.tpp"

#include "ijkgenscalar.tpp"


namespace IJKGENSCALAR {

  // ********************************************************
  // VERSION
  // ********************************************************

  inline const char * ProgramVersion()
    { return("0.2.0"); }

  // ********************************************************
  // TYPES
  // ********************************************************

  typedef size_t AXIS_SIZE_TYPE;
  typedef float SCALAR_TYPE;
  typedef float GRADIENT_COORD_TYPE;
  typedef float COORD_TYPE;
  typedef int GRID_COORD_TYPE;
  typedef float RADIUS_TYPE;
  typedef float DIFF_TYPE;
  typedef float ANGLE_TYPE;
  typedef int VERTEX_INDEX_TYPE;
  typedef VERTEX_INDEX_TYPE CUBE_INDEX_TYPE;
  typedef int NUM_TYPE;
  typedef unsigned long ISOTABLE_INDEX_TYPE;
  typedef unsigned int SEED_TYPE;
  typedef unsigned int MIN_MAX_TYPE;
  typedef unsigned int MAX_EXPONENT_TYPE;
  typedef unsigned long BOUNDARY_BITS;
  typedef IJK::GRID_NEIGHBORS<int, AXIS_SIZE_TYPE, int, int, int> BASE_GRID;
  typedef IJK::GRID_SPACING<COORD_TYPE, BASE_GRID> GRID;
  typedef IJK::SCALAR_GRID<GRID, SCALAR_TYPE> SCALAR_GRID;
  typedef IJK::VECTOR_GRID<GRID, NUM_TYPE, GRADIENT_COORD_TYPE> GRADIENT_GRID;
  typedef IJK::CUBE_FACE_INFO<NUM_TYPE,NUM_TYPE,NUM_TYPE> CUBE_TYPE;

  typedef FIELD_OBJECT_PROPERTIES_T
  <int, COORD_TYPE, GRADIENT_COORD_TYPE, RADIUS_TYPE, DIFF_TYPE, 
   ANGLE_TYPE, SCALAR_TYPE, NUM_TYPE>
  OBJECT_PROPERTIES;
  typedef IJKGENGEOM::GEOM_PARAM_T<OBJECT_PROPERTIES, AXIS_SIZE_TYPE, SEED_TYPE>
  GEOM_PARAM;
  typedef SET_SCALAR_DATA_T<CUBE_INDEX_TYPE, ISOTABLE_INDEX_TYPE,
                            VERTEX_INDEX_TYPE, NUM_TYPE> 
  SET_SCALAR_DATA_TYPE;
  typedef FIELD_PARAM_T
  <SET_SCALAR_DATA_TYPE, GEOM_PARAM, SCALAR_TYPE, MIN_MAX_TYPE,
   MAX_EXPONENT_TYPE, ISOTABLE_INDEX_TYPE> 
  FIELD_PARAM;

  typedef FIELD_INFO_T<OBJECT_PROPERTIES, SCALAR_GRID, GRADIENT_GRID> FIELD_INFO;
};

#endif
