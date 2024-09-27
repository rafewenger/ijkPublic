/*!
 *  @file ijkdual_types.h
 *  @brief Type definitions for ijkdual.
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2023 Rephael Wenger

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

#ifndef _IJKDUAL_TYPES_H_
#define _IJKDUAL_TYPES_H_

#include <string>

#include "ijk.tpp"
#include "ijkscalar_grid.tpp"
#include "ijkmerge.tpp"

namespace IJKDUAL {


// **************************************************
// SCALAR TYPES
// **************************************************

  typedef int MESH_SIZE_TYPE;    ///< Mesh size type.
  typedef int GRID_SIZE_TYPE;    ///< Grid size type.
  typedef float SCALAR_TYPE;     ///< Scalar value type.
  typedef float COORD_TYPE;      ///< Isosurface vertex coordinate type.
  typedef int VERTEX_INDEX;      ///< Grid vertex index type.
  typedef float GRADIENT_TYPE;   ///< Gradient coordinate type.
  typedef int GRID_COORD_TYPE;   ///< Grid vertex coordinate type.
  typedef unsigned int AXIS_SIZE_TYPE;    ///< Axis size type.
  typedef int ISO_VERTEX_INDEX;  ///< Isosurface vertex index type.
  typedef int ISO_EDGE_INDEX;    ///< Isosurface edge index type.
  typedef unsigned int MERGE_INDEX;        ///< Merge index type.
  typedef unsigned char DIRECTION_TYPE;    ///< Direction type.

  /// @brief Boundary bits type.
  /// - Number of bits must be at least 2*dimension.
  typedef unsigned int BOUNDARY_BITS_TYPE;

  /// @brief Edge index type.
  /// - Vertex and edge indices must have the same type.
  typedef VERTEX_INDEX EDGE_INDEX;

  /// @brief Cube index type.
  /// - Vertex and cube indices must have the same type.
  typedef VERTEX_INDEX CUBE_INDEX;

  /// @brief Cube index type.
  /// - Vertex and cube indices must have the same type.
  typedef char CUBE_EDGE_INDEX;  

  /// @brief Isosurface vertex degree type.
  typedef unsigned int ISO_VERTEX_DEGREE;
  
  /// Facet index type.
  typedef unsigned char FACET_INDEX;

  /// Facet vertex index type.
  typedef unsigned char FACET_VERTEX_INDEX;

  /// @brief Type for integer representing a set of facets.
  /// - k'th bit is 1 if facet k is in the set.
  typedef unsigned int FACET_BITS_TYPE;

  /// Vertex degree type.
  typedef unsigned char DEGREE_TYPE;

  /// Isosurface lookup table index type.
  typedef unsigned long TABLE_INDEX;

  /// Random seed type. Used for testing.
  typedef unsigned long RANDOM_SEED_TYPE;


  // **************************************************
  // ARRAY TYPES
  // **************************************************

  typedef std::vector<COORD_TYPE>          /// Isosurface coordinate array.
    COORD_ARRAY;
  typedef std::vector<VERTEX_INDEX>        /// Grid vertex index array.
    VERTEX_INDEX_ARRAY; 
  typedef std::vector<ISO_VERTEX_INDEX>    /// Isosurface vertex index array.
    ISO_VERTEX_INDEX_ARRAY; 
  typedef std::vector<SCALAR_TYPE> SCALAR_ARRAY; ///< Scalar array.

  /// Array of facet vertex indices.
  typedef std::vector<FACET_VERTEX_INDEX> FACET_VERTEX_INDEX_ARRAY;

  /// Array of directions.
  typedef std::vector<DIRECTION_TYPE> DIRECTION_ARRAY;


  // **************************************************
  // ENUMERATED TYPES
  // **************************************************

  /*!
   *  @brief Mesh type.
   *  - SIMPLICIAL_COMPLEX: Mesh elements are triangles/tetrahedra/simplices.
   *  - CUBE_COMPLEX: Mesh elements are quads/cubes/hypercubes.
   *  - MIXED_MESH: Mesh elements are mixture of different polytopes,
   *    e.g. triangles and quadilaterals.
   *
   */
  typedef enum
    { SIMPLICIAL_COMPLEX, CUBE_COMPLEX, MIXED_MESH }
    MESH_TYPE;

  /*!
   *  @brief Isosurface quad/cube/hypercube representation.
   *  - Order of listing of quad/cube/hypercube vertices.
   *  - COORDINATE_VERTEX_ORDER: Coordinate vertex order,
   *    increasing x, then increasing y, then ....
   *    - 2D mesh; (0,0), (1,0), (0,1), (1,1).
   *    - 3D mesh; (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), 
   *      (0,1,1), (1.1.1).
   *  - CIRCULAR_VERTEX_ORDER: Order clockwise or counterclockwise
   *    around quad boundary (2D mesh only).
   *    - 2D mesh: (0,0), (1,0), (0,1), (1,1).
   */
  typedef enum { COORDINATE_VERTEX_ORDER, CIRCULAR_VERTEX_ORDER }
    CUBE_VERTEX_ORDER;

  /*!
   *  @brief Interpolation type.
   *  - LINEAR_INTERPOLATION: Determine the location of an isosurface vertex
   *    using linear interpolation on the endpoints of the cube, pyramid
   *    or simplex containing the isosurface vertex.
   *  - MULTILINEAR_INTERPOLATION: Determine the location of an 
   *    isosurface vertex using multilinear interpolation on the cube vertices.
   */
  typedef enum
    { LINEAR_INTERPOLATION, MULTILINEAR_INTERPOLATION }
    INTERPOLATION_TYPE;

  /*!
   *  @brief Isosurface vertex position method.
   *  - CUBE_CENTER: Position isosurface vertices at cube centers.
   *  - CENTROID_EDGE_ISO: Position isosurface vertices at the centroid
   *                       of the edge isosurface intersections.
   *  - DIAGONAL_INTERPOLATION: Position some isosurface vertices on diagonals
   *    using linear interpolation.  Position all other isosurface vertices
   *    at the centroid of the edge isosurface intersections.
   *  - RANDOM_ISOV_POS: Random position within cube. Used for testing.
   */
  typedef enum
    { CUBE_CENTER, CENTROID_EDGE_ISO, DIAGONAL_INTERPOLATION,
      IVOL_LIFTED02, QDUAL_INTERPOLATION, RANDOM_ISOV_POS } 
    VERTEX_POSITION_METHOD;

  /*!
   *  @brief Quadrilateral triangulation method.
   *  - UNDEFINED_TRI: Triangulation method not defined.
   *  - UNIFORM_TRI: Uniform triangulation. Split all quadrilaterals
   *    with diagonal (0,2).
   *  - SPLIT_MAX_ANGLE: Triangulate using diagonal incident 
   *    on maximum quadrilateral angle.
   *  - MAX_MIN_ANGLE: Triangulate using diagonal so that triangulation
   *    maximizes minimum angle.
   *  - TRI4_ALL_QUADS: Split each quadrilateral into four triangles,
   *    adding an additional vertex "in" each quadrilateral.
   *  - TRI4_MAX_MIN_ANGLE: Triangulate using either diagonals
   *    or splitting into four triangles, choosing the triangulation
   *    that maximizes the minimum angle.
   *  - TRI_DIAGONALS_OUTSIDE_ENVELOPE: Triangulate any quad that has
   *    a diagonal outside the envelope. (Other quads are not triangulated.)
   */
  typedef enum 
    { UNDEFINED_TRI, UNIFORM_TRI, SPLIT_MAX_ANGLE, MAX_MIN_ANGLE,
      TRI4_ALL_QUADS, TRI4_MAX_MIN_ANGLE, 
      TRI_DIAGONALS_OUTSIDE_ENVELOPE }
    QUAD_TRI_METHOD;

  /*!
   *  @brief Triangulation methods for diagonals outside the envelope.
   *  - Method for determining triangulation when exactly one diagonal
   *    is outside the envelope, and there is a choice of either
   *    using the other diagonal or splitting into four triangles.
   *  - ENVELOPE_UNDEFINED_TRI: Undefined envelope triangulation method.
   *  - ENVELOPE_ONLY_TRI2: If exactly one diagonal is outside the envelope,
   *    triangulate using the other diagonal. If both diagonals are
   *    outside the envelope, do nothing.
   *  - ENVELOPE_ONLY_TRI4: If either diagonal is outside the envelope,
   *    triangulate into four triangles.
   *  - ENVELOPE_PREFER_TRI2: If exactly one diagonal is outside the envelope,
   *    triangulate using the other diagonal. If both diagonals are outside
   *    the envelope, triangulate into four triangles.
   *  - ENVELOPE_MAX_MIN_ANGLE: If only one diagonal is outside the envelope,
   *    triangulate using the other diagonal or into four triangles,
   *    choosing the triangulation that maximizes the minimum angle.
   */
  typedef enum
    { ENVELOPE_UNDEFINED_TRI, ENVELOPE_ONLY_TRI2, ENVELOPE_ONLY_TRI4,
      ENVELOPE_PREFER_TRI2, ENVELOPE_MAX_MIN_ANGLE }
  ENVELOPE_QUAD_TRI_METHOD;

  /*!
   *  @brief Random distribution type.
   */
  typedef enum
    { UNIFORM_DISTRIBUTION, U_QUADRATIC_DISTRIBUTION,
      UNDEFINED_DISTRIBUTION } RANDOM_DISTRIBUTION;

  typedef enum { TRIANGLE, QUADRILATERAL, DEGENERATE_POLYGON }
    POLYGON_TYPE;


// **************************************************
// CLASSES
// **************************************************

  typedef IJK::BOX<VERTEX_INDEX> GRID_BOX;  ///< Grid box type.
}

#endif
