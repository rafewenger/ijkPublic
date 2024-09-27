/// \file ijkmesh.h
/// compute mesh information

/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2018 Rephael Wenger

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

#ifndef _IJKMESH_H_
#define _IJKMESH_H_

#include "ijkcube.tpp"
#include "ijkmesh.tpp"
#include "ijkmesh_datastruct.tpp"

namespace IJKMESH {

  // **********************************************************************
  // Global constants
  // **********************************************************************

  const int DIM2(2);
  const int DIM3(3);
  const int NUM_VERT_PER_TETRAHEDRON(4);
  const int NUM_VERT_PER_HEXAHEDRON(8);
  const int TRI_ENCODING_BIT_SET_SIZE(16);


  // **********************************************************************
  // Types
  // **********************************************************************

  typedef float COORD_TYPE;
  typedef float COLOR_TYPE;


  // **********************************************************************
  // Datastructure types
  // **********************************************************************

  typedef typename IJK::CUBE_FACE_INFO<int,int,int> CUBE_TYPE;

  typedef typename IJK::POLYMESH<int,int> POLYMESH_TYPE;

  typedef typename IJK::VERTEX_POLY_INCIDENCE<int,int> 
  VERTEX_POLY_INCIDENCE_TYPE;
  typedef typename IJK::HEXMESH_ADJACENCY<int,int,int> HEXMESH_ADJACENCY_TYPE;
  typedef typename IJK::CUBE_FACE_INFO<int,int,int> CUBE_TYPE;


  // **********************************************************************
  // Arrays of measurements
  // **********************************************************************

  class MESH_MEASUREMENTS {

  public:
    std::vector<COORD_TYPE> scaled_jacobian_determinant;
    std::vector<COORD_TYPE> jacobian_shape;
  };


  // **********************************************************************
  // Data structure to replace coordinates.
  // **********************************************************************

  /// Class to replace coordinates.
  class REPLACE_COORD {
    
  public:

    /// Index [0..(dimension-1)] of coordinate to be replaced.
    int icoord;

    /// Replace coordinates that have value old_coord.
    COORD_TYPE old_coord;

    /// Replace value old_coord with value new_coord.
    COORD_TYPE new_coord;
  };

};

#endif
