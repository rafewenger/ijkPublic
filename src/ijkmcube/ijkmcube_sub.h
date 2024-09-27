/*!
 *  @file ijkmcube_sub.h
 *  @brief Subroutines for ijkmcube.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2-23 Rephael Wenger

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

#ifndef _IJKMCUBE_SUB_
#define _IJKMCUBE_SUB_

#include <string>

#include "ijk.tpp"
#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"

namespace IJKMCUBE {

// **************************************************
// MARCHING CUBES SUBROUTINES
// **************************************************

  /*!
   *  @brief Merge identical edges in elist0[].
   *  @param elist0[] List of edges specified by endpoints.
   *    - (elist0[2*i],elist0[2*i+1]) = Edge i.
   *    - elist0[] may contain duplicate edges.
   *  @param[out] elist1[] New list of edges specified endpoints.
   *    - elist1[] contains edges in elist0[] with all duplicates removed.
   *    - (elist1[2*i],elist1[2*i+1]) = New edge i.
   *  @param[in,out] eindex[] Indices of edges.
   *    - Originally, eindex[j] is index of some edge e in elist0[].
   *    - eindex[j] is replaced by index of e in elist1[].
   *  @param param Merge parameters.
   */
  void merge_identical_edges
    (const std::vector<VERTEX_INDEX> & elist0,
     std::vector<VERTEX_INDEX> & elist1, 
     std::vector<ISO_VERTEX_INDEX> & eindex,
     const MERGE_EDGES_PARAMETERS & param);

  /*!
   *  @brief Compute position of isosurface vertices using linear interpolation.
   *  - Version where all vertices lie on regular grid edges.
   *  @param scalar_grid Scalar grid.
   *  @param isovalue Isosurface scalar value.
   *  @param isotable_type Type of isosurface table (BINARY, NEP, or IVOL).
   *  @param elist[] List of edges containing isosurface vertices.
   *    - Edges are listed by indices that encode index of lower endpoint
   *      and edge direction.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+k] = k'th coordinate of isosurface vertex i.
   *    @pre Array coord[] is preallocated to length at least
   *      dimension*elist.size().
   */
  void position_iso_vertices_linear
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const std::vector<EDGE_INDEX> & elist, COORD_TYPE * coord);

  /*!
   *  @brief Compute position of isosurface vertices using linear interpolation.
   *  - Version with intersected edges specified by pairs of endpoints
   *    in elist_endpoint[].
   *  @param scalar_grid Scalar grid.
   *  @param isovalue Isosurface scalar value.
   *  @param isotable_type Type of isosurface table (BINARY, NEP, or IVOL).
   *  @param elist_endpoint[] List of edges specified by pairs of edge endpoints.
   *    (elist_endpoint[2*i], elist_endpoint[2*i+1]) is edge i.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+k] = k'th coordinate of isosurface vertex i.
   *    @pre Array coord[] is preallocated to length at least
   *      dimension*elist.size().
   */
  void position_iso_vertices_linearB
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const std::vector<VERTEX_INDEX> & elist_endpoint, COORD_TYPE * coord);

  /*!
   *  @brief Compute position of isosurface vertices using linear interpolation.
   *  - Version with intersected edges specified by pairs of endpoints
   *    in elist_endpoint[].
   *  - Version where new, non-grid mesh vertices are stored
   *    in new_mesh_vertices.
   *  @param new_mesh_vertices Data structure storing new, non-grid mesh vertices.
   *  - If iv < scalar_grid.NumVertices(), then iv is a grid vertex.
   *  - If iv >= scalar_grid.NumVertices(), then iv is a
   *    new, non-grid mesh vertex, whose index in new_mesh_vertices
   *    is (iv-scalar_grid.NumVertices().
   *  @param elist_endpoint[] List of edges specified by pairs of edge endpoints.
   *    (elist_endpoint[2*i], elist_endpoint[2*i+1]) is edge i.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+k] = k'th coordinate of isosurface vertex i.
   *    @pre Array coord[] is preallocated to length at least
   *      dimension*elist.size().
   */
  void position_iso_vertices_linearB
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & elist_endpoint, COORD_TYPE * coord);

  /*!
   *  @brief Compute position of isosurface vertices using multilinear interpolation.
   *  - Apply multilinear interpolation on scalar values of cube vertices
   *    to determine scalar values of new mesh vertices.
   *  @param new_mesh_vertices Data structure storing new, non-grid mesh vertices.
   *  - If iv < scalar_grid.NumVertices(), then iv is a grid vertex.
   *  - If iv >= scalar_grid.NumVertices(), then iv is a
   *    new, non-grid mesh vertex, whose index in new_mesh_vertices
   *    is (iv-scalar_grid.NumVertices().
   *  @param elist_endpoint[] List of edges specified by pairs of edge endpoints.
   *    (elist_endpoint[2*i], elist_endpoint[2*i+1]) is edge i.
   */
  void position_iso_vertices_multilinear
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
     const ISOTABLE_TYPE isotable_type,
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & elist_endpoint, COORD_TYPE * coord);

  /*!
   *  @brief Merge identical isosurface vertices and compute their position.
   *  - Grid edges are represented by a single integer identifier.
   *  @param scalar_grid Scalar grid data.
   *  @param isovalue Isosurface scalar value.
   *  @param iso_simplices[] Array of isosurface simplex vertices.
   *         iso_simplices[dimension*is+k] =
   *           grid edge containing k'th vertex of simplex is.
   *  @param isotable_type Type of isosurface table (BINARY, NEP, or IVOL).
   *  @param isotable_type Type of isosurface table.
   *  @param simplex vert[] List of isosurface simplex vertices.
   *         simplex_vert[dimension*is+k] = k'th vertex of simplex is.
   *  @param vertex_coord[] List of isosurface vertex coordinates.
   *         vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv.
   *  @param merge_data Internal data structure for merging identical edges.
   *  @param[out] mcube_info Information about running time
   *    and grid cubes and edges.
   */
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     const ISOTABLE_TYPE isotable_type,
     std::vector<ISO_VERTEX_INDEX> & simplex_vert,
     std::vector<COORD_TYPE> & vertex_coord,
     MERGE_DATA & merge_data, MCUBE_INFO & mcube_info);

  /*!
   *  @brief Merge identical isosurface vertices and compute their position.
   *  - Version with intersected edges specified by pairs of endpoints
   *    in elist_endpoint[].
   *  @param scalar_grid Scalar grid data.
   *  @param isovalue Isosurface scalar value.
   *  @param elist_endpoint[] List of isosurface simplex vertices specified
   *    by edges intersected by the isosurface.
   *    - (elist_endpoint[2*i], elist_endpoint[2*i+1]) =
   *      endpoints of edge containing isosurface vertex i.
   *  @param isotable_type Type of isosurface table (BINARY, NEP, or IVOL).
   *  @param isotable_type Type of isosurface table.
   *  @param simplex vert[] List of isosurface simplex vertices.
   *         simplex_vert[dimension*is+k] = k'th vertex of simplex is.
   *  @param vertex_coord[] List of isosurface vertex coordinates.
   *         vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv.
   *  @param[out] mcube_info Information about running time
   *    and grid cubes and edges.
   */
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const std::vector<VERTEX_INDEX> & elist_endpoint,
     const ISOTABLE_TYPE isotable_type,
     std::vector<ISO_VERTEX_INDEX> & simplex_vert,
     std::vector<COORD_TYPE> & vertex_coord,
     const MERGE_EDGES_PARAMETERS &, MCUBE_INFO &);

  /*!
   *  @brief Merge identical isosurface vertices and compute their position.
   *  - Version with intersected edges specified by pairs of endpoints
   *    in elist_endpoint[].
   *  - Version where new, non-grid mesh vertices are stored
   *    in new_mesh_vertices.
   *  @param new_mesh_vertices Data structure storing new, non-grid mesh vertices.
   *  - If iv < scalar_grid.NumVertices(), then iv is a grid vertex.
   *  - If iv >= scalar_grid.NumVertices(), then iv is a
   *    new, non-grid mesh vertex, whose index in new_mesh_vertices
   *    is (iv-scalar_grid.NumVertices().
   */
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & elist_endpoint,
     const ISOTABLE_TYPE isotable_type,
     std::vector<ISO_VERTEX_INDEX> & simplex_vert,
     std::vector<COORD_TYPE> & vertex_coord,
     const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info);

  /*!
   *  @brief Merge identical isosurface vertices and compute their position.
   *  - Version that allows multilinear interpolation in computing
   *    scalar values of new mesh vertices.
   *  @param elist_endpoint[] List of isosurface simplex vertices specified
   *    by edges intersected by the isosurface.
   *    - (elist_endpoint[2*i], elist_endpoint[2*i+1]) =
   *      endpoints of edge containing isosurface vertex i.
   *  @param interpolation_type LINEAR_INTERPOLATION or MULTILINEAR_INTERPOLATION.
   *    - Type of interpolation used in computing scalar values
   *    of new mesh vertices.
   *
   */
  void merge_and_position_vertices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
     const MC_MESH_VERTEX_LIST & new_mesh_vertices,
     const std::vector<VERTEX_INDEX> & endpoint,
     const ISOTABLE_TYPE isotable_type,
     const INTERPOLATION_TYPE interpolation_type,
     std::vector<ISO_VERTEX_INDEX> & simplex_vert,
     std::vector<COORD_TYPE> & vertex_coord,
     const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info);

  /// @brief Convert isosurface vertex indices to pairs of endpoints
  ///   representing grid edges containing isosurface vertices.
  void convert_indices_to_endpoint_pairs
    (const MC_SCALAR_GRID_BASE & scalar_grid,
     const std::vector<ISO_VERTEX_INDEX> iso_vlist,
     const MERGE_DATA & merge_data,
     std::vector<VERTEX_INDEX> & elist_endpoint);


// **************************************************
// NEP MARCHING CUBES SUBROUTINES
// **************************************************

  /// Identify isosurface patches which lie in a facet.
  /// Set is_in_facet flag for each entry in isosurface lookup table.
  void set_in_facet_iso_patches(NEP_ISOSURFACE_TABLE & isotable);

// **************************************************
// INTERVAL VOLUME SUBROUTINES
// **************************************************

  /// Compute position of interval volume vertices using linear interpolation.
  void position_ivol_vertices_linear
    (const MC_SCALAR_GRID_BASE & scalar_grid,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord);

}

#endif
