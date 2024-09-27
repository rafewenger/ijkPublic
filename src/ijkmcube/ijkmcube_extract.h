/*!
 *  @file ijkmcube_extract.h
 *  @brief Subroutines for extracting isosurface mesh.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2013 Rephael Wenger

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

#ifndef _IJKMCUBE_EXTRACT_
#define _IJKMCUBE_EXTRACT_

#include <string>

#include "ijk.tpp"
#include "ijkmcube_types.h"
#include "ijkmcube_datastruct.h"

namespace IJKMCUBE {

// **************************************************
// EXTRACT ISOSURFACE MESH
// **************************************************

  /*!
   *  @brief Extract isosurface simplices from grid cubes.
   *  @param[out] iso_simplices[] List of isosurface simplex vertices
   *    - iso_simplices[dimension*is+k] =
   *      grid edge containing k'th vertex of simplex is.
   */
  void extract_iso_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

  /*!
   *  @brief Extract isosurface simplices and return grid cubes containing simplices.
   *  @param[out] iso_simplices[] List of isosurface simplex vertices
   *    - iso_simplices[dimension*is+k] =
   *      grid edge containing k'th vertex of simplex is.
   *  @param[out] cube_list[i] Cube containing i'th isosurface simplex.
   */
  void extract_iso_simplices_and_cube_info
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & cube_list,
     MCUBE_INFO & mcube_info);

  /*!
   *  @brief Extract isosurface simplices and return edges containing simplex vertices.
   *  @param[out] iso_simplices[] List of isosurface simplex vertices
   *    - iso_simplices[dimension*is+k] =
   *      grid edge containing k'th vertex of simplex is.
   *  @param[out] endpoint[] Endpoints of edges containing isosurface vertices.
   *    - (endpoint[2*i], endpoint[2*i+1]) =
   *      grid edge containing isosurface vertex i.
   */
  void extract_iso_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices with specified isosurface_topology.
  void extract_iso_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const POLY_ISOTABLE & poly_isotable,
     const SCALAR_TYPE isovalue,
     const ISOSURFACE_TOPOLOGY isosurface_topology,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint, 
     MC_MESH_VERTEX_LIST & new_mesh_vertices,
     MCUBE_INFO & mcube_info);

  /*!
   *  @brief Extract isosurface simplices using the asymptotic decider.
   *  - Returns list representing isosurface simplices
   *  @param scalar_grid Scalar grid.
   *  @param poly_isotable Isosurface lookup tables for polyhedra.
   *    - Isosurface lookup tables for cube, pyramid and simplex.
   *  @param isovalue Isosurface scalar value.
   *  @param[out] iso_simplices Isosurface simplex vertices.
   *    - iso_simplices[dimension*is+k] = k'th vertex of isosurface simplex is.
   *  @param[out] endpoint[] Endpoints of edges containing isosurface vertices.
   *    - (endpoint[2*i], endpoint[2*i+1]) = grid edge containing isosurface
   *      vertex i.
   */
  void extract_iso_simplices_adecider
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const POLY_ISOTABLE & poly_isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint, 
     MC_MESH_VERTEX_LIST & new_mesh_vertices,
     MCUBE_INFO & mcube_info);

  /// @brief Extract isosurface simplices.
  /// - Break cubes at saddle points to approximate linear topology.
  void extract_iso_simplices_saddle
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const POLY_ISOTABLE & poly_isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint, 
     MC_MESH_VERTEX_LIST & new_mesh_vertices,
     MCUBE_INFO & mcube_info);

  /// @brief Extract isosurface simplices resolving some cube ambiguities.
  /// - In ambiguous cubes with no ambiguous facets, resolve ambiguity
  ///   based on average scalar values of cube vertices.
  void extract_iso_simplices_cube_decider
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable, 
     const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info, 
     const SCALAR_TYPE isovalue, 
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices resolving some cube ambiguities.
  /// In ambiguous cubes with no ambiguous facets, resolve ambiguity
  ///   based on average scalar values of cube vertices.
  void extract_iso_simplices_cube_decider
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable, 
     const ISOSURFACE_TABLE_AMBIG_INFO & ambig_info, 
     const SCALAR_TYPE isovalue, 
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices from grid cubes in list.
  void extract_iso_simplices_from_list
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, 
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);

  /// Extract isosurface simplices using octree.
  void extract_iso_simplices_from_octree
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info,  const IJKOCTREE::OCTREE & octree);

  /// Extract isosurface simplices using partition into uniform regions.
  void extract_iso_simplices_from_minmax
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info, const MINMAX_REGIONS & minmax);

// **************************************************
// EXTRACT ISOSURFACE MESH USING NEP ISOTABLE
// **************************************************

  /*!
   *  @brief Extract isosurface simplices from grid cubes using nep lookup table.
   *  - Use negative, equals, positive (nep) isosurface lookup table.
   *  - @param nep_num_dup Parameter for processing duplicate isosurface triangles.
   *      - See \ref nepNumDupDoc "MC_DATA_FLAGS::nep_num_doc".
   */
  void extract_iso_simplices_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info);

  /// Extract nep isosurface vertices lying on the grid boundary.
  void extract_iso_simplices_nep_boundary
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int num_cube_facet_vertices,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, int & num_non_empty_cubes);

  /*!
   *  @brief Extract isosurface simplices using octree.
   *  - Use negative, equals, positive (nep) isosurface lookup table.
   *  - @param nep_num_dup Parameter for processing duplicate isosurface triangles.
   *      - See \ref nepNumDupDoc "MC_DATA_FLAGS::nep_num_doc".
   */
  void extract_iso_simplices_from_octree_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info,
     const IJKOCTREE::OCTREE & octree);

  /*!
   *  @brief Extract isosurface simplices using partition into uniform regions.
   *  - Use negative, equals, positive (nep) isosurface lookup table.
   *  - @param nep_num_dup Parameter for processing duplicate isosurface triangles.
   *      - See \ref nepNumDupDoc "MC_DATA_FLAGS::nep_num_doc".
   */
  void extract_iso_simplices_from_minmax_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices, MCUBE_INFO & mcube_info,
     const MINMAX_REGIONS & minmax);

  /*!
   *  @brief Extract isosurface simplices from list of cubes.
   *  - Use negative, equals, positive (nep) isosurface lookup table.
   *  - @param nep_num_dup Parameter for processing duplicate isosurface triangles.
   *      - See \ref nepNumDupDoc "MC_DATA_FLAGS::nep_num_doc".
   */
  void extract_iso_simplices_from_list_nep
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const NEP_ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue, const int nep_num_dup,
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     const bool extract_from_boundary,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     MCUBE_INFO & mcube_info);


// **************************************************
// EXTRACT INTERVAL VOLUME MESH
// **************************************************

  /// Extract interval volume simplices.
  void extract_ivol_simplices
    (const MC_SCALAR_GRID_BASE & scalar_grid, const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
     VERTEX_INDEX & num_non_empty_cubes);

  /// Extract interval volume simplices from list of cubes.
  void extract_ivol_simplices_from_list
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const ISOSURFACE_TABLE & isotable,
     const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
     const VERTEX_INDEX * vlist, const VERTEX_INDEX num_cubes, 
     std::vector<ISO_VERTEX_INDEX> & ivol_simplices,
     VERTEX_INDEX & num_non_empty_cubes);

// **************************************************
// EXTRACT ISOSURFACE MESH FROM MULTIRESOLUTION GRID
// **************************************************

  /// Extract isosurface simplices from multiresolution cubes, pyramids and tetrahedra.
  void multires_extract
    (const MC_SCALAR_GRID_BASE & scalar_grid, 
     const MULTIRES_GRID & multires_grid,
     const POLY_ISOTABLE & poly_isotable, const SCALAR_TYPE isovalue,
     std::vector<ISO_VERTEX_INDEX> & iso_simplices,
     std::vector<VERTEX_INDEX> & endpoint,
     MCUBE_INFO & mcube_info);
}

#endif
