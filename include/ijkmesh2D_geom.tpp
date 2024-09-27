/*!
 *  @file ijkmesh2D_geom.tpp
 *  @brief ijk template classes for 2D polygonal mesh geometric primitives.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2019-2020 Rephael Wenger

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

#ifndef _IJKMESH2D_GEOM_
#define _IJKMESH2D_GEOM_

#include "ijk.tpp"
#include "ijkmesh2D_datastruct.tpp"
#include "ijkcoord.tpp"

#include <algorithm>
#include <numeric>
#include <vector>


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // **************************************************
  /// @name COMPUTE EDGE MIDPOINTS
  // **************************************************

  ///@{

  /// Computer midpoint of line segment between two mesh vertices.
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  /// @pre Array midpoint_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE, typename ITYPE2>
  void compute_mesh2D_line_segment_midpoint
  (const DTYPE dimension, const COORD_TYPE0 vertex_coord[],
   const MESH_TYPE & mesh, const ITYPE iv0, const ITYPE2 iv1,
   COORD_TYPE1 midpoint_coord[])
  {
    const COORD_TYPE0 * v0coord = vertex_coord + iv0*dimension;
    const COORD_TYPE0 * v1coord = vertex_coord + iv1*dimension;

    IJK::compute_midpoint
      (dimension, v0coord, v1coord, midpoint_coord);
  }


  /// Computer midpoint of line segment between two mesh vertices.
  ///  - Version using C++ STL vector vertex_coord[].
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  /// @pre Array midpoint_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE, typename ITYPE2>
  void compute_mesh2D_line_segment_midpoint
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vertex_coord,
   const MESH_TYPE & mesh, const ITYPE iv0, const ITYPE2 iv1,
   COORD_TYPE1 midpoint_coord[])
  {
    compute_mesh2D_line_segment_midpoint
      (dimension, IJK::vector2pointer(vertex_coord), mesh, iv0, iv1,
       midpoint_coord);
  }


  /// Computer midpoint of line segment between two mesh vertices.
  ///  - Version using C++ STL vector vertex_coord[].
  ///  - Version using C++ STL vector midpoint_coord[].
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE, typename ITYPE2>
  void compute_mesh2D_line_segment_midpoint
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vertex_coord,
   const MESH_TYPE & mesh, const ITYPE iv0, const ITYPE2 iv1,
   std::vector<COORD_TYPE1> & midpoint_coord)
  {
    compute_mesh2D_line_segment_midpoint
      (dimension, vertex_coord, mesh, iv0, iv1, 
       IJK::vector2pointer(midpoint_coord));
  }


  /// Compute midpoint of mesh edge.
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  /// @pre Array midpoint_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_edge_midpoint
  (const DTYPE dimension, const COORD_TYPE0 vertex_coord[],
   const MESH_TYPE & mesh, const ITYPE ihalf_edge,
   COORD_TYPE1 midpoint_coord[])
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const VERTEX_INDEX_TYPE iv0 = mesh.FromVertexIndex(ihalf_edge);
    const VERTEX_INDEX_TYPE iv1 = mesh.ToVertexIndex(ihalf_edge);

    IJK::compute_mesh2D_line_segment_midpoint
      (dimension, vertex_coord, mesh, iv0, iv1, midpoint_coord);
  }


  /// Computer midpoint of mesh edge.
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  /// @pre Array midpoint_coord[] is preallocated to length at least dimension.
  ///  - Version using C++ STL vector vertex_coord[].
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_edge_midpoint
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vertex_coord,
   const MESH_TYPE & mesh, const ITYPE ihalf_edge,
   COORD_TYPE1 midpoint_coord[])
  {
    compute_mesh2D_edge_midpoint
      (dimension, IJK::vector2pointer(vertex_coord), mesh, 
       ihalf_edge, midpoint_coord);
  }


  /// Computer midpoint of mesh edge.
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  /// @pre Array midpoint_coord[] is preallocated to length at least dimension.
  /// - Version using C++ STL vector vertex_coord[].
  /// - Version using C++ STL vector midpoint_coord[].
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_edge_midpoint
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vertex_coord,
   const MESH_TYPE & mesh, const ITYPE ihalf_edge,
   std::vector<COORD_TYPE1> & midpoint_coord)
  {
    compute_mesh2D_edge_midpoint
      (dimension, vertex_coord, mesh, ihalf_edge, 
       IJK::vector2pointer(midpoint_coord));
  }


  /// Computer midpoint of midpoints of mesh edge.
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  /// @pre Array midpoint_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPEA, typename ITYPEB>
  void compute_mesh2D_midpoint_of_edge_midpoints
  (const DTYPE dimension, const COORD_TYPE0 vertex_coord[],
   const MESH_TYPE & mesh, 
   const ITYPEA ihalf_edgeA, const ITYPEB ihalf_edgeB,
   COORD_TYPE1 midpoint_coord[])
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const VERTEX_INDEX_TYPE ivA0 = mesh.FromVertexIndex(ihalf_edgeA);
    const VERTEX_INDEX_TYPE ivA1 = mesh.ToVertexIndex(ihalf_edgeA);
    const VERTEX_INDEX_TYPE ivB0 = mesh.FromVertexIndex(ihalf_edgeB);
    const VERTEX_INDEX_TYPE ivB1 = mesh.ToVertexIndex(ihalf_edgeB);
    const COORD_TYPE0 * vA0coord = vertex_coord + ivA0*dimension;
    const COORD_TYPE0 * vA1coord = vertex_coord + ivA1*dimension;
    const COORD_TYPE0 * vB0coord = vertex_coord + ivB0*dimension;
    const COORD_TYPE0 * vB1coord = vertex_coord + ivB1*dimension;

    IJK::copy_coord(dimension, vA0coord, midpoint_coord);
    IJK::add_coord(dimension, midpoint_coord, vA1coord, midpoint_coord);
    IJK::add_coord(dimension, midpoint_coord, vB0coord, midpoint_coord);
    IJK::add_coord(dimension, midpoint_coord, vB1coord, midpoint_coord);
    IJK::divide_coord(dimension, 4.0, midpoint_coord, midpoint_coord);
  }


  /// Computer midpoint of midpoints of mesh edge.
  ///  - Version using C++ STL vector vertex_coord[].
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  /// @pre Array midpoint_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPEA, typename ITYPEB>
  void compute_mesh2D_midpoint_of_edge_midpoints
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vertex_coord,
   const MESH_TYPE & mesh, 
   const ITYPEA ihalf_edgeA, const ITYPEB ihalf_edgeB,
   COORD_TYPE1 midpoint_coord[])
  {
    compute_mesh2D_midpoint_of_edge_midpoints
      (dimension, IJK::vector2pointer(vertex_coord), mesh, 
       ihalf_edgeA, ihalf_edgeB, midpoint_coord);
  }


  /// Computer midpoint of midpoints of mesh edge.
  /// @param[out] midpoint_coord[] Coordinates of midpoint.
  /// @pre Array midpoint_coord[] is preallocated to length at least dimension.
  /// - Version using C++ STL vector vertex_coord[].
  /// - Version using C++ STL vector midpoint_coord[].
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPEA, typename ITYPEB>
  void compute_mesh2D_midpoint_of_edge_midpoints
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vertex_coord,
   const MESH_TYPE & mesh, 
   const ITYPEA ihalf_edgeA, const ITYPEB ihalf_edgeB,
   std::vector<COORD_TYPE1> & midpoint_coord)
  {
    compute_mesh2D_midpoint_of_edge_midpoints
      (dimension, vertex_coord, mesh, ihalf_edgeA, ihalf_edgeB, 
       IJK::vector2pointer(midpoint_coord));
  }


  ///@}


  // **************************************************
  /// @name COMPUTE POLYGON CENTROIDS
  // **************************************************

  ///@{

  /// Computer centroid of polygon ipoly.
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_polygon_centroid
  (const DTYPE dimension, const COORD_TYPE0 vcoord[],
   const MESH_TYPE & mesh, const ITYPE ipoly, 
   COORD_TYPE1 centroid_coord[])
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const NUMBER_TYPE numv = mesh.NumPolygonVertices(ipoly);

    IJK::set_coord(dimension, 0, centroid_coord);

    for (NUMBER_TYPE j = 0; j < numv; j++) {

      const VERTEX_INDEX_TYPE jv = mesh.IndexOfVertexInPolygon(ipoly, j);

      IJK::add_coord
        (dimension, vcoord+jv*dimension, centroid_coord, centroid_coord);
    }

    if (numv > 0) {
      IJK::divide_coord(dimension, numv, centroid_coord, centroid_coord);
    }
  }


  /// Computer centroid of polygon ipoly.
  ///  - Version using C++ STL vector vertex_coord[].
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_polygon_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const MESH_TYPE & mesh, const ITYPE ipoly, 
   COORD_TYPE1 centroid_coord[])
  {
    compute_mesh2D_polygon_centroid
      (dimension, IJK::vector2pointer(vcoord), mesh, ipoly, centroid_coord);
  }

  /// Computer centroid of polygon ipoly.
  ///  - Version using C++ STL vector vertex_coord[] and centroid_coord[].
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_polygon_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const MESH_TYPE & mesh, const ITYPE ipoly, 
   std::vector<COORD_TYPE1> & centroid_coord)
  {
    centroid_coord.resize(dimension);

    compute_mesh2D_polygon_centroid
      (dimension, vcoord, mesh, ipoly, IJK::vector2pointer(centroid_coord));
  }


  /// Compute centroid of two polygons sharing an edge.
  /// @param ihalf_edgeA Index of one half edge of the shared edge.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_two_polygon_centroid
  (const DTYPE dimension, const COORD_TYPE0 vcoord[],
   const MESH_TYPE & mesh, const ITYPE ihalf_edgeA, 
   COORD_TYPE1 centroid_coord[])
  {
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const POLYGON_INDEX_TYPE ipolyA =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);
    const POLYGON_INDEX_TYPE ipolyB =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edgeB);

    const NUMBER_TYPE numvA = mesh.NumPolygonVertices(ipolyA);
    const NUMBER_TYPE numvB = mesh.NumPolygonVertices(ipolyB);

    IJK::set_coord(dimension, 0, centroid_coord);

    // Add coordinates of all vertices in ipolyA.
    for (NUMBER_TYPE j = 0; j < numvA; j++) {

      const VERTEX_INDEX_TYPE jv = mesh.IndexOfVertexInPolygon(ipolyA, j);

      IJK::add_coord
        (dimension, vcoord+jv*dimension, centroid_coord, centroid_coord);
    }

    // Add coordinates of vertices in ipolyA, 
    //   not including endpoints of ipolyB.
    HALF_EDGE_INDEX_TYPE ihalf_edge = 
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edgeB);
    ihalf_edge = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge);

    for (NUMBER_TYPE j = 0; j+2 < numvB; j++) {

      const VERTEX_INDEX_TYPE jv = mesh.FromVertexIndex(ihalf_edge);

      IJK::add_coord
        (dimension, vcoord+jv*dimension, centroid_coord, centroid_coord);

      ihalf_edge = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge);
    }

    if (numvA+numvB > 0) {
      IJK::divide_coord
        (dimension, numvA+numvB-2, centroid_coord, centroid_coord);
    }
  }


  /// Compute centroid of two polygons sharing an edge.
  /// - C++ STL vector type for array vcoord[].
  /// @param ihalf_edgeA Index of one half edge of the shared edge.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_two_polygon_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const MESH_TYPE & mesh, const ITYPE ihalf_edgeA, 
   COORD_TYPE1 centroid_coord[])
  {
    compute_mesh2D_two_polygon_centroid
      (dimension, IJK::vector2pointer(vcoord), mesh, 
       ihalf_edgeA, centroid_coord);
  }


  /// Compute centroid of two polygons sharing an edge.
  /// - C++ STL vector type for array vcoord[].
  /// - C++ STL vector type for array centroid_coord[].
  /// @param ihalf_edgeA Index of one half edge of the shared edge.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE>
  void compute_mesh2D_two_polygon_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const MESH_TYPE & mesh, const ITYPE ihalf_edgeA, 
   std::vector<COORD_TYPE1> & centroid_coord)
  {
    centroid_coord.resize(dimension);

    compute_mesh2D_two_polygon_centroid
      (dimension, vcoord, mesh, ihalf_edgeA, 
       IJK::vector2pointer(centroid_coord));
  }


  /// Computer centroid of vertices of polygon ipoly
  ///   and one additional point.
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename COORD_TYPE2,
            typename MESH_TYPE, typename ITYPE1>
  void compute_mesh2D_polygon_plus_point_centroid
  (const DTYPE dimension, const COORD_TYPE0 vcoord[],
   const COORD_TYPE1 point_coord[],
   const MESH_TYPE & mesh, const ITYPE1 ipoly,
   COORD_TYPE2 centroid_coord[])
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    const NUMBER_TYPE numv = mesh.NumPolygonVertices(ipoly);

    IJK::copy_coord(dimension, point_coord, centroid_coord);

    for (NUMBER_TYPE j = 0; j < numv; j++) {

      const VERTEX_INDEX_TYPE jv = mesh.IndexOfVertexInPolygon(ipoly, j);

      IJK::add_coord
        (dimension, vcoord+jv*dimension, centroid_coord, centroid_coord);
    }

    IJK::divide_coord(dimension, numv+1, centroid_coord, centroid_coord);
  }


  /// Computer centroid of vertices of polygon ipoly
  ///   and one additional point.
  ///  - Version using C++ STL vector vertex_coord[].
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename COORD_TYPE2,
            typename MESH_TYPE, typename ITYPE1>
  void compute_mesh2D_polygon_plus_point_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const COORD_TYPE1 point_coord[],
   const MESH_TYPE & mesh, const ITYPE1 ipoly,
   COORD_TYPE2 centroid_coord[])
  {
    compute_mesh2D_polygon_plus_point_centroid
      (dimension, IJK::vector2pointer(vcoord), point_coord,
       mesh, ipoly, centroid_coord);
  }


  /// Computer centroid of vertices of polygon ipoly
  ///   and one additional point.
  ///  - Version using C++ STL vector vertex_coord[].
  ///  - Version using C++ STL vector point_coord[].
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename COORD_TYPE2,
            typename MESH_TYPE, typename ITYPE1>
  void compute_mesh2D_polygon_plus_point_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const std::vector<COORD_TYPE1> & point_coord,
   const MESH_TYPE & mesh, const ITYPE1 ipoly,
   COORD_TYPE2 centroid_coord[])
  {
    compute_mesh2D_polygon_plus_point_centroid
      (dimension, vcoord, IJK::vector2pointer(point_coord),
       mesh, ipoly, centroid_coord);
  }


  /// Computer centroid of vertices of polygon ipoly
  ///   and one additional point.
  ///  - Version using C++ STL vector vertex_coord[].
  ///  - Version using C++ STL vector point_coord[].
  ///  - Version using C++ STL vector centroid_coord[].
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename COORD_TYPE2,
            typename MESH_TYPE, typename ITYPE1>
  void compute_mesh2D_polygon_plus_point_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const std::vector<COORD_TYPE1> & point_coord,
   const MESH_TYPE & mesh, const ITYPE1 ipoly,
   std::vector<COORD_TYPE2> & centroid_coord)
  {
    compute_mesh2D_polygon_plus_point_centroid
      (dimension, vcoord, point_coord, mesh, ipoly, 
       IJK::vector2pointer(centroid_coord));
  }


  /// Computer centroid of vertices of polygon ipoly
  ///   and one additional vertex.
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE1, typename ITYPE2>
  void compute_mesh2D_polygon_plus_vertex_centroid
  (const DTYPE dimension, const COORD_TYPE0 vcoord[],
   const MESH_TYPE & mesh, const ITYPE1 ipoly, ITYPE2 kv,
   COORD_TYPE1 centroid_coord[])
  {
    const COORD_TYPE0 * point_coord = vcoord + kv*dimension;

    compute_mesh2D_polygon_plus_point_centroid
      (dimension, vcoord, point_coord, mesh, ipoly, centroid_coord);
  }


  /// Computer centroid of vertices of polygon ipoly
  ///   and one additional vertex.
  ///  - Version using C++ STL vector vertex_coord[].
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE1, typename ITYPE2>
  void compute_mesh2D_polygon_plus_vertex_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const MESH_TYPE & mesh, const ITYPE1 ipoly, const ITYPE2 kv,
   COORD_TYPE1 centroid_coord[])
  {
    compute_mesh2D_polygon_plus_vertex_centroid
      (dimension, IJK::vector2pointer(vcoord), mesh, ipoly, kv,
       centroid_coord);
  }


  /// Computer centroid of vertices of polygon ipoly
  ///   and one additional vertex.
  ///  - Version using C++ STL vector vertex_coord[] and centroid_coord[].
  /// @param[out] centroid_coord[] Coordinates of centroid.
  /// @pre Array centroid_coord[] is preallocated to length at least dimension.
  template <typename DTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE, typename ITYPE1, typename ITYPE2>
  void compute_mesh2D_polygon_plus_vertex_centroid
  (const DTYPE dimension, const std::vector<COORD_TYPE0> & vcoord,
   const MESH_TYPE & mesh, const ITYPE1 ipoly, const ITYPE2 kv,
   std::vector<COORD_TYPE1> & centroid_coord)
  {
    centroid_coord.resize(dimension);

    compute_mesh2D_polygon_plus_vertex_centroid
      (dimension, vcoord, mesh, ipoly, kv,
       IJK::vector2pointer(centroid_coord));
  }


  ///@}


  // **************************************************
  /// @name COMPUTE EDGE LENGTHS
  // **************************************************

  ///@{

  /// Compute edge length squared.
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE>
  COORD_TYPE compute_edge_length_squared
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE ihalf_edge)
  {
    COORD_TYPE length_squared;

    const COORD_TYPE * v0coord = 
      vcoord + dimension*mesh.FromVertexIndex(ihalf_edge);
    const COORD_TYPE * v1coord = 
      vcoord + dimension*mesh.ToVertexIndex(ihalf_edge);
    IJK::compute_distance_squared
      (dimension, v0coord, v1coord, length_squared);

    return(length_squared);
  }


  /// Compute edge length squared.
  /// - C++ STL vector type for array vcoord[].
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE>
  COORD_TYPE compute_edge_length_squared
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord,
   const MESH_TYPE & mesh, const ITYPE ihalf_edge)
  {
    return(compute_edge_length_squared
           (dimension, IJK::vector2pointer(vcoord), mesh, ihalf_edge));
  }

  ///@}


  // **************************************************
  /// @name COMPUTE LONG/SHORT EDGES/QUADS
  // **************************************************

  ///@{

  /// Return index of longest half edge in polygon ipoly.
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE>
  typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
  compute_longest_polygon_half_edge
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE ipoly)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    COORD_TYPE length_squaredA(0.0), length_squaredB(0.0);

    const HALF_EDGE_INDEX_TYPE ihalf_edge0 = mesh.HalfEdgeIndex(ipoly, 0);
    length_squaredA = compute_edge_length_squared
      (dimension, vcoord, mesh, ihalf_edge0);

    HALF_EDGE_INDEX_TYPE ilongest_half_edge = ihalf_edge0;

    for (ITYPE i = 1; i < mesh.NumPolygonEdges(ipoly); i++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge = mesh.HalfEdgeIndex(ipoly, i);
      length_squaredB = compute_edge_length_squared
        (dimension, vcoord, mesh, ihalf_edge);

      if (length_squaredB > length_squaredA) {
        length_squaredA = length_squaredB;
        ilongest_half_edge = ihalf_edge;
      }
    }

    return(ilongest_half_edge);
  }


  /// Return index of longest half edge in polygon ipoly.
  /// - C++ STL vector type for array vcoord[].
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE>
  typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
  compute_longest_polygon_half_edge
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord, 
   const MESH_TYPE & mesh, const ITYPE ipoly)
  {
    return(compute_longest_polygon_half_edge
           (dimension, IJK::vector2pointer(vcoord), mesh, ipoly));
  }


  /// Compute indices of two longest half edges in polygon ipoly.
  template <typename DTYPE, typename COORD_TYPE, typename MESH_TYPE, 
            typename ITYPE0, typename ITYPE1, typename ITYPE2,
            typename LTYPE1, typename LTYPE2>
  void compute_two_longest_polygon_half_edges
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE0 ipoly,
   ITYPE1 & ilongest_half_edge, ITYPE2 & isecond_longest_half_edge,
   LTYPE1 & longest_length_squared, LTYPE2 & second_longest_length_squared)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    COORD_TYPE length_squared(0.0);

    const HALF_EDGE_INDEX_TYPE ihalf_edge0 = mesh.HalfEdgeIndex(ipoly, 0);
    longest_length_squared = compute_edge_length_squared
      (dimension, vcoord, mesh, ihalf_edge0);


    ilongest_half_edge = ihalf_edge0;
    isecond_longest_half_edge = ihalf_edge0;
    second_longest_length_squared = longest_length_squared;
    
    for (ITYPE0 i = 1; i < mesh.NumPolygonEdges(ipoly); i++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge = mesh.HalfEdgeIndex(ipoly, i);
      length_squared = compute_edge_length_squared
        (dimension, vcoord, mesh, ihalf_edge);

      if (length_squared > longest_length_squared) {
        second_longest_length_squared = longest_length_squared;
        longest_length_squared = length_squared;
        isecond_longest_half_edge = ilongest_half_edge;
        ilongest_half_edge = ihalf_edge;
      }
      else if (isecond_longest_half_edge == ilongest_half_edge ||
               length_squared > second_longest_length_squared) {
        second_longest_length_squared = length_squared;
        isecond_longest_half_edge = ihalf_edge;
      }
    }

  }


  /// Compute indices of two longest half edges in polygon ipoly.
  ///  - Version using C++ STL vector vertex_coord[].
  template <typename DTYPE, typename COORD_TYPE, typename MESH_TYPE, 
            typename ITYPE0, typename ITYPE1, typename ITYPE2,
            typename LTYPE1, typename LTYPE2>
  void compute_two_longest_polygon_half_edges
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord,
   const MESH_TYPE & mesh, const ITYPE0 ipoly,
   ITYPE1 & ilongest_half_edge, ITYPE2 & isecond_longest_half_edge,
   LTYPE1 & longest_length_squared, LTYPE2 & second_longest_length_squared)
  {
    compute_two_longest_polygon_half_edges
      (dimension, IJK::vector2pointer(vcoord), mesh, ipoly,
       ilongest_half_edge, isecond_longest_half_edge,
       longest_length_squared, second_longest_length_squared);
  }


  /// Compute indices of three longest half edges in polygon ipoly.
  template <typename DTYPE, typename COORD_TYPE, typename MESH_TYPE, 
            typename ITYPE0, typename ITYPE1, typename LTYPE>
  void compute_three_longest_polygon_half_edges
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE0 ipoly,
   ITYPE1 ihalf_edge_long[3], LTYPE length_squared[3])
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    const int THREE = 3;
    COORD_TYPE Lsquared(0.0);
    IJK::PROCEDURE_ERROR error
      ("IJK::compute_three_longest_polygon_half_edges");

    if (mesh.NumPolygonEdges(ipoly) == 0) {
      error.AddMessage("Polygon ", ipoly, " has no edges.");
      throw error;
    }

    // Set ihalf_edge_long[0].
    const ITYPE1 ihalf_edge0 = mesh.HalfEdgeIndex(ipoly, 0);
    ihalf_edge_long[0] = ihalf_edge0;
    length_squared[0] = 
      compute_edge_length_squared(dimension, vcoord, mesh, ihalf_edge0);
    
    // Initialize other ihalf_edge_long entries.
    for (ITYPE0 i = 1; i < THREE; i++) {
      ihalf_edge_long[i] = ihalf_edge0;
      length_squared[i] = length_squared[0];
    }

    for (ITYPE0 i = 1; i < mesh.NumPolygonEdges(ipoly); i++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge = mesh.HalfEdgeIndex(ipoly, i);
      Lsquared = compute_edge_length_squared
        (dimension, vcoord, mesh, ihalf_edge);

      int list_length = THREE;
      if (i < THREE) {
        ihalf_edge_long[i] = ihalf_edge; 
        length_squared[i] = Lsquared;
        list_length = i+1;
      }
      else {
        if (Lsquared > length_squared[THREE-1]) {
          // Replace ihalf_edge_long[2] with ihalf_edge.
          ihalf_edge_long[THREE-1] = ihalf_edge; 
          length_squared[THREE-1] = Lsquared;
        }
        else {
          // ihalf_edge is shorter than all three edges currently stored.
          continue;
        }
      }

      // Sort by swapping last half edge inserted in ihalf_edge_long[]
      int j = list_length-1;
      while (j > 0) {
        if (length_squared[j] > length_squared[j-1]) {
          std::swap(ihalf_edge_long[j], ihalf_edge_long[j-1]);
          std::swap(length_squared[j], length_squared[j-1]);
          j--;
        }
        else
          { break; }
      }
    }

  }


  /// Compute indices of three longest half edges in polygon ipoly.
  ///  - Version using C++ STL vector vertex_coord[].
  template <typename DTYPE, typename COORD_TYPE, typename MESH_TYPE, 
            typename ITYPE0, typename ITYPE1, typename LTYPE>
  void compute_three_longest_polygon_half_edges
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord,
   const MESH_TYPE & mesh, const ITYPE0 ipoly,
   ITYPE1 ihalf_edge_long[3], LTYPE length_squared[3])
  {
    compute_three_longest_polygon_half_edges
      (dimension, IJK::vector2pointer(vcoord), mesh, ipoly,
       ihalf_edge_long, length_squared);
  }


  /// Return index of shortest half edge in polygon ipoly.
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE, typename LTYPE>
  typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
  compute_shortest_polygon_half_edge
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE ipoly,
   LTYPE & shortest_length_squared)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    LTYPE length_squaredB;

    const HALF_EDGE_INDEX_TYPE ihalf_edge0 = mesh.HalfEdgeIndex(ipoly, 0);
    shortest_length_squared = compute_edge_length_squared
      (dimension, vcoord, mesh, ihalf_edge0);

    HALF_EDGE_INDEX_TYPE ishortest_half_edge = ihalf_edge0;
    
    for (ITYPE i = 0; i < mesh.NumPolygonEdges(ipoly); i++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge = mesh.HalfEdgeIndex(ipoly, i);
      length_squaredB = compute_edge_length_squared
      (dimension, vcoord, mesh, ihalf_edge);

      if (length_squaredB < shortest_length_squared) {
        shortest_length_squared = length_squaredB;
        ishortest_half_edge = ihalf_edge;
      }
    }

    return(ishortest_half_edge);
  }


  /// Return index of shortest half edge in polygon ipoly.
  /// - Version that does not return shortest_length_squared.
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE>
  typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
  compute_shortest_polygon_half_edge
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE ipoly)
  {
    COORD_TYPE shortest_length_squared;
    return(compute_shortest_polygon_half_edge
           (dimension, vcoord, mesh, ipoly, shortest_length_squared));
  }


  /// Return index of shortest half edge in polygon ipoly.
  /// - C++ STL vector type for array vcoord[].
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE, typename LTYPE>
  typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
  compute_shortest_polygon_half_edge
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord, 
   const MESH_TYPE & mesh, const ITYPE ipoly,
   LTYPE & shortest_length_squared)
  {
    return(compute_shortest_polygon_half_edge
           (dimension, IJK::vector2pointer(vcoord), mesh, ipoly,
            shortest_length_squared));
  }


  /// Return index of shortest half edge in polygon ipoly.
  /// - Version that does not return shortest_length_squared.
  /// - C++ STL vector type for array vcoord[].
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE>
  typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
  compute_shortest_polygon_half_edge
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord, 
   const MESH_TYPE & mesh, const ITYPE ipoly)
  {
    return(compute_shortest_polygon_half_edge
           (dimension, IJK::vector2pointer(vcoord), mesh, ipoly));
  }


  /// Get index of shortest non-split half edge in polygon ipoly.
  /// - Neither endpoint of the half edge is a split vertex.
  /// - Return false if no such edge is found.
  template <typename DTYPE, typename COORD_TYPE, typename MESH_TYPE, 
            typename ITYPE2, typename ITYPE3, typename LTYPE>
  bool compute_shortest_non_split_polygon_half_edge
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE2 ipoly,
   ITYPE3 & ishortest, LTYPE & shortest_length_squared)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::NUMBER_TYPE
      NUMBER_TYPE;

    LTYPE length_squaredB;
    bool flag_found(false);

    const HALF_EDGE_INDEX_TYPE ihalf_edge0 = mesh.HalfEdgeIndex(ipoly, 0);
    const COORD_TYPE * v0coord = 
      vcoord + dimension*mesh.FromVertexIndex(ihalf_edge0);
    const COORD_TYPE * v1coord = 
      vcoord + dimension*mesh.ToVertexIndex(ihalf_edge0);
    IJK::compute_distance_squared
      (dimension, v0coord, v1coord, shortest_length_squared);

    HALF_EDGE_INDEX_TYPE ishortest_half_edge = ihalf_edge0;
    if (!mesh.DoesFromVertexSplitEdge(ishortest_half_edge) &&
        !mesh.DoesToVertexSplitEdge(ishortest_half_edge))
      { flag_found = true; }
    
    for (NUMBER_TYPE i = 0; i < mesh.NumPolygonEdges(ipoly); i++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge = mesh.HalfEdgeIndex(ipoly, i);
      v0coord = vcoord + dimension*mesh.FromVertexIndex(ihalf_edge);
      v1coord = vcoord + dimension*mesh.ToVertexIndex(ihalf_edge);

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, length_squaredB);

      if (!mesh.DoesFromVertexSplitEdge(ihalf_edge) &&
          !mesh.DoesToVertexSplitEdge(ihalf_edge)) {

        if (!flag_found || (length_squaredB < shortest_length_squared)) {
          shortest_length_squared = length_squaredB;
          ishortest_half_edge = ihalf_edge;
          flag_found = true;
        }
      }

    }

    return(flag_found);
  }


  /// Get index of shortest non-split half edge in polygon ipoly.
  /// - Neither endpoint of the half edge is a split vertex.
  /// - Return false if no such edge is found.
  /// - C++ STL vector type for array vcoord[].
  template <typename DTYPE, typename COORD_TYPE, typename MESH_TYPE, 
            typename ITYPE2, typename ITYPE3, typename LTYPE>
  bool compute_shortest_non_split_polygon_half_edge
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord, 
   const MESH_TYPE & mesh, const ITYPE2 ipoly,
   ITYPE3 & ishortest, LTYPE & shortest_length_squared)
  {
    return(compute_shortest_non_split_polygon_half_edge
           (dimension, IJK::vector2pointer(vcoord), mesh, ipoly,
            ishortest, shortest_length_squared));
  }


  /// Return true if ihalf_edge0 is longest half edge in polygon ipoly.
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE>
  typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
  is_longest_polygon_half_edge
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE ihalf_edge0)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE 
      POLYGON_INDEX_TYPE;

    const POLYGON_INDEX_TYPE ipoly =
      mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    COORD_TYPE length_squared0(0.0), length_squared(0.0);

    const COORD_TYPE * v0coord = 
      vcoord + dimension*mesh.FromVertexIndex(ihalf_edge0);
    const COORD_TYPE * v1coord = 
      vcoord + dimension*mesh.ToVertexIndex(ihalf_edge0);
    IJK::compute_distance_squared
      (dimension, v0coord, v1coord, length_squared0);

    HALF_EDGE_INDEX_TYPE ihalf_edge = 
      mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);

    for (ITYPE i = i; i < mesh.NumPolygonEdges(ipoly); i++) {
      v0coord = vcoord + dimension*mesh.FromVertexIndex(ihalf_edge);
      v1coord = vcoord + dimension*mesh.ToVertexIndex(ihalf_edge);

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, length_squared);

      if (length_squared > length_squared0) { return(false); }
    }

    return(true);
  }


  /// Return true if ihalf_edge0 is longest half edge in polygon ipoly.
  /// - C++ STL vector type for array vcoord[].
  template <typename DTYPE, typename COORD_TYPE,
            typename MESH_TYPE, typename ITYPE>
  typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
  is_longest_polygon_half_edge
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord,
   const MESH_TYPE & mesh, const ITYPE ihalf_edge0)
  {
    return(is_longest_polygon_half_edge
           (dimension, IJK::vector2pointer(vcoord), mesh, ihalf_edge0)); 
  }


  /// Return true if distance between vertices is "long".
  /// - Return true if distance squared is greater than or equal to L.
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE0, typename ITYPE1, 
            typename LTYPE>
  bool is_long_distance(const DTYPE dimension,
                        const CTYPE vertex_coord[],
                        const ITYPE0 iv0, const ITYPE1 iv1,
                        const LTYPE L)
  {
    const CTYPE * v0coord = vertex_coord + iv0*dimension;
    const CTYPE * v1coord = vertex_coord + iv1*dimension;
    LTYPE distance_squared;

    IJK::compute_distance_squared
      (dimension, v0coord, v1coord, distance_squared);

    if (distance_squared >= L) { return(true); }

    return(false);
  }


  /// Return true if distance between vertices is "long".
  /// - Return true if distance squared is greater than or equal to L.
  /// - Version using C++ STL vector vertex_coord[].
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE0, typename ITYPE1, 
            typename LTYPE>
  bool is_long_distance(const DTYPE dimension,
                        const std::vector<CTYPE> & vertex_coord,
                        const ITYPE0 iv0, const ITYPE1 iv1,
                        const LTYPE L)
  {
    return(is_long_distance(dimension, IJK::vector2pointer(vertex_coord),
                            iv0, iv1, L));
  }


  /// Return true if quad is "long".
  /// @param[out] jlongest_edge_loc Location of longest quad edge.
  ///             Index is 0, 1, 2 or 3.
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE, typename ITYPE2, typename RTYPE>
  bool is_long_quad(const DTYPE dimension,
                    const CTYPE vertex_coord[],
                    const ITYPE quad_vert[], const RTYPE Rsquared,
                    ITYPE2 & jlongest_edge_loc)
  {
    const ITYPE2 NUM_VERT_PER_QUAD(4);
    CTYPE length_squared[NUM_VERT_PER_QUAD];

    jlongest_edge_loc = 0;
    for (ITYPE2 j0 = 0; j0 < NUM_VERT_PER_QUAD; j0++) {

      const ITYPE2 j1 = (j0+1)%NUM_VERT_PER_QUAD;
      const ITYPE iv0 = quad_vert[j0];
      const ITYPE iv1 = quad_vert[j1];
      const CTYPE * v0coord = vertex_coord + iv0*dimension;
      const CTYPE * v1coord = vertex_coord + iv1*dimension;

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, length_squared[j0]);

      if (length_squared[j0] > length_squared[jlongest_edge_loc])
        { jlongest_edge_loc = j0; }
    }

    const ITYPE2 j1 = (jlongest_edge_loc+1)%NUM_VERT_PER_QUAD;
    const ITYPE2 j2 = (jlongest_edge_loc+2)%NUM_VERT_PER_QUAD;
    const ITYPE2 j3 = (jlongest_edge_loc+3)%NUM_VERT_PER_QUAD;

    if (length_squared[j2] < length_squared[j1]*Rsquared)
      { return(false); }

    if (length_squared[j2] < length_squared[j3]*Rsquared)
      { return(false); }

    return(true);
  }


  /// Return true if quad is "long".
  /// @param[out] jlongest_edge_loc Location of longest quad edge.
  ///             Index is 0, 1, 2 or 3.
  /// - Version using C++ STL vector vertex_coord[].
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE, typename ITYPE2, typename RTYPE>
  bool is_long_quad(const DTYPE dimension,
                    const std::vector<CTYPE> & vertex_coord,
                    const ITYPE quad_vert[], const RTYPE Rsquared,
                    ITYPE2 & jlongest_edge_loc)
  {
    return(is_long_quad
           (dimension, IJK::vector2pointer(vertex_coord),
            quad_vert, Rsquared, jlongest_edge_loc));
  }


  /// Return true if shortest edge is between two long edges.
  /// @param[out] jlongest_edge_loc Location of longest quad edge.
  ///             Index is 0, 1, 2 or 3.
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE, typename ITYPE2, typename RTYPE>
  bool is_shortest_quad_edge_between_long_edges
  (const DTYPE dimension,
   const CTYPE vertex_coord[],
   const ITYPE quad_vert[], const RTYPE Rsquared,
   ITYPE2 & jshortest_edge_loc)
  {
    const ITYPE2 NUM_VERT_PER_QUAD(4);
    CTYPE length_squared[NUM_VERT_PER_QUAD];

    jshortest_edge_loc = 0;
    for (ITYPE2 j0 = 0; j0 < NUM_VERT_PER_QUAD; j0++) {

      const ITYPE2 j1 = (j0+1)%NUM_VERT_PER_QUAD;
      const ITYPE iv0 = quad_vert[j0];
      const ITYPE iv1 = quad_vert[j1];
      const CTYPE * v0coord = vertex_coord + iv0*dimension;
      const CTYPE * v1coord = vertex_coord + iv1*dimension;

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, length_squared[j0]);

      if (length_squared[j0] < length_squared[jshortest_edge_loc])
        { jshortest_edge_loc = j0; }
    }

    const ITYPE2 j1 = (jshortest_edge_loc+1)%NUM_VERT_PER_QUAD;
    const ITYPE2 j2 = (jshortest_edge_loc+2)%NUM_VERT_PER_QUAD;
    const ITYPE2 j3 = (jshortest_edge_loc+3)%NUM_VERT_PER_QUAD;

    const CTYPE shortest_length_squared = 
      length_squared[jshortest_edge_loc];

    if (shortest_length_squared*Rsquared > length_squared[j1])
      { return(false); }
    if (shortest_length_squared*Rsquared > length_squared[j3])
      { return(false); }

    return(true);
  }


  /// Return true if shortest edge is between two long edges.
  /// - Version using C++ STL vector vertex_coord[].
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE, typename ITYPE2,typename RTYPE>
  bool is_shortest_quad_edge_between_long_edges
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const ITYPE quad_vert[], const RTYPE R,
   ITYPE2 & jshortest_edge_loc)
  {
    return(is_shortest_quad_edge_between_long_edges
           (dimension, IJK::vector2pointer(vertex_coord),
            quad_vert, R, jshortest_edge_loc));
  }


  /// Return true if quad has a very short edge.
  /// @param[out] jshortest_edge_loc Location of longest quad edge.
  ///             Index is 0, 1, 2 or 3.
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE, typename ITYPE2, typename RTYPE>
  bool does_quad_have_very_short_edge
  (const DTYPE dimension,
   const CTYPE vertex_coord[],
   const ITYPE quad_vert[], const RTYPE Rsquared,
   ITYPE2 & jshortest_edge_loc)
  {
    const ITYPE2 NUM_VERT_PER_QUAD(4);
    CTYPE length_squared[NUM_VERT_PER_QUAD];

    jshortest_edge_loc = 0;
    for (ITYPE2 j0 = 0; j0 < NUM_VERT_PER_QUAD; j0++) {

      const ITYPE2 j1 = (j0+1)%NUM_VERT_PER_QUAD;
      const ITYPE iv0 = quad_vert[j0];
      const ITYPE iv1 = quad_vert[j1];
      const CTYPE * v0coord = vertex_coord + iv0*dimension;
      const CTYPE * v1coord = vertex_coord + iv1*dimension;

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, length_squared[j0]);

      if (length_squared[j0] < length_squared[jshortest_edge_loc])
        { jshortest_edge_loc = j0; }
    }


    const ITYPE2 j1 = (jshortest_edge_loc+1)%NUM_VERT_PER_QUAD;
    const ITYPE2 j2 = (jshortest_edge_loc+2)%NUM_VERT_PER_QUAD;
    const ITYPE2 j3 = (jshortest_edge_loc+3)%NUM_VERT_PER_QUAD;

    const CTYPE shortest_length_squared = 
      length_squared[jshortest_edge_loc];

    if (shortest_length_squared*Rsquared > length_squared[j1])
      { return(false); }
    if (shortest_length_squared*Rsquared > length_squared[j2])
      { return(false); }
    if (shortest_length_squared*Rsquared > length_squared[j3])
      { return(false); }

    return(true);
  }


  /// Return true if quad has a very short edge.
  /// - Version using C++ STL vector vertex_coord[].
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE, typename ITYPE2,typename RTYPE>
  bool does_quad_have_very_short_edge
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const ITYPE quad_vert[], const RTYPE R,
   ITYPE2 & jshortest_edge_loc)
  {
    return(does_quad_have_very_short_edge
           (dimension, IJK::vector2pointer(vertex_coord),
            quad_vert, R, jshortest_edge_loc));
  }


  /// Return true if two longest edges in polygon are adjacent
  ///   and both are greater than L.
  /// - If true, returns the long half edge that precedes the other.
  /// @param[out] ilong_half_edge0 Half edge ilong_half_edge0 and
  ///   IndexOfNextHalfInPolygon(ilong_half_edge0) are the two
  ///   adjacent long half edges.
  template <typename DTYPE, typename COORD_TYPE, typename MESH_TYPE, 
            typename ITYPE0, typename ITYPE2, typename RTYPE>
  bool has_two_adjacent_long_edges
  (const DTYPE dimension, const COORD_TYPE vcoord[],
   const MESH_TYPE & mesh, const ITYPE0 ipoly, const RTYPE Rsquared,
   ITYPE2 & ilong_half_edge0)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;

    HALF_EDGE_INDEX_TYPE ilongest_half_edge, isecond_longest_half_edge;
    COORD_TYPE longest_length_squared, second_longest_length_squared;
    COORD_TYPE shortest_length_squared;

    compute_two_longest_polygon_half_edges
      (dimension, vcoord, mesh, ipoly, 
       ilongest_half_edge, isecond_longest_half_edge,
       longest_length_squared, second_longest_length_squared);

    // Initialize.
    ilong_half_edge0 = ilongest_half_edge;

    if (mesh.IndexOfNextHalfEdgeInPolygon(ilongest_half_edge) ==
        isecond_longest_half_edge)
      { ilong_half_edge0 = ilongest_half_edge; }
    else if (mesh.IndexOfNextHalfEdgeInPolygon(isecond_longest_half_edge) ==
             ilongest_half_edge)
      { ilong_half_edge0 = isecond_longest_half_edge; }
    else {
      // Two longest edges are not adjacent.
      return(false);
    }

    compute_shortest_polygon_half_edge
      (dimension, vcoord, mesh, ipoly, shortest_length_squared);

    const COORD_TYPE L = shortest_length_squared*Rsquared;

    if (longest_length_squared < L) { return(false); }
    if (second_longest_length_squared < L) { return(false); }

    return(true);
  }


  /// Return true if two longest edges in polygon are adjacent
  ///   and both are greater than L.
  /// - If true, returns the long half edge that precedes the other.
  ///  - Version using C++ STL vector vertex_coord[].
  /// @param[out] ilong_half_edge0 Half edge ilong_half_edge0 and
  ///   IndexOfNextHalfInPolygon(ilong_half_edge0) are the two
  ///   adjacent long half edges.
  template <typename DTYPE, typename COORD_TYPE, typename MESH_TYPE, 
            typename ITYPE0, typename ITYPE2, typename RTYPE>
  bool has_two_adjacent_long_edges
  (const DTYPE dimension, const std::vector<COORD_TYPE> & vcoord,
   const MESH_TYPE & mesh, const ITYPE0 ipoly, const RTYPE Rsquared,
   ITYPE2 & ilong_half_edge0)
  {
    return(has_two_adjacent_long_edges
           (dimension, IJK::vector2pointer(vcoord), mesh,
            ipoly, Rsquared, ilong_half_edge0));
  }

  ///@}

}


#endif
