/// \file ijkoriented_mesh.tpp
/// @brief ijk template classes for oriented polygonal mesh.
/// - Version 0.4.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 22018 Rephael Wenger

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

#ifndef _IJKORIENTED_MESH_
#define _IJKORIENTED_MESH_

#include "ijk.tpp"

#include <algorithm>
#include <vector>


namespace IJK {

  // *****************************************************************
  // Class ORIENTED_POLYGONAL_MESH_VERTEX
  // *****************************************************************

  template <typename ITYPE>
  class ORIENTED_POLYGONAL_MESH_VERTEX {

  public:

    /// Index of some half edge with vertex as origin.
    ITYPE incident_half_edge;
  };


  // *****************************************************************
  // Class ORIENTED_POLYGONAL_MESH_POLYGON
  // *****************************************************************

  template <typename ITYPE>
  class ORIENTED_POLYGONAL_MESH_POLYGON {

  public:

    /// Index of some half edge around polygon.
    ITYPE incident_half_edge;
  };


  // *****************************************************************
  // Class DIRECTED_HALF_EDGE
  // *****************************************************************

  template <typename VTYPE, typename ITYPE>
  class DIRECTED_HALF_EDGE {

  public:

    /// Vertex at origin of directed half edge.
    VTYPE origin;

    /// Index of next half edge around polygon.
    ITYPE next_half_edge;
  };

  // *****************************************************************
  // Class ORIENTED_POLYGONAL_MESH
  // *****************************************************************

  /// Oriented polygonal mesh
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename NTYPE>
  class ORIENTED_POLYGONAL_YMESH {

  protected:
    std::vector<VTYPE> vertex;
    std::vector<HALFE_TYPE> half_edge;
    std::vector<PTYPE> poly;

  public:

    /// Vertex type.
    typedef VTYPE VERTEX_TYPE;

    /// Half edge type
    typedef HALFE_TYPE HALF_EDGE_TYPE;

    /// Polygon type.
    typedef PTYPE VERTEX_TYPE;

    /// NUMBER type.
    typedef NTYPE NUMBER_TYPE;

  public:
    /// constructor
    ORIENTED_POLYGONAL_MESH(){};

    /// Number of Vertices.
    NTYPE NumVert() const               
    { return(vertex.size()); }

    /// Number of edges.
    NTYPE NumEdge() const               
    { return(half_edge.size()/2); }

    /// Number of polygons.
    NTYPE NumPoly() const               
    { return(poly.size()); }

    /// Opposite half edge.
    NTYPE OppositeHalfEdge(const NTYPE ihalfe) const
    { return(((ihalfe%2):(ihalfe-1),(ihalfe+1))); }

    /// Next half edge.
    NTYPE NextHalfEdge(const NTYPE ihalfe) const
    { return(half_edge[ihalfe].next_half_edge); }

    /// Origin.
    NTYPE Origin(const NTYPE ihalfe) const
    { return(half_edge[ihalfe].origin_vertex); }

    /// Destination.
    NTYPE Destination(const NTYPE ihalfe) const
    { return(Origin(OppositeHalfEdge(ihalfe))); }

    /// Set number of vertices.
    void SetNumVert();

    /// Add polygon.


  };

}

#endif
