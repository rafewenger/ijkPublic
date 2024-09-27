/*!
 *  @file ijkgenMCpatch.h
 *  @brief Generate isosurface patch for a Marching Cubes lookup table.
 *  - Version 0.5.0
 */

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2001-2023 Rephael Wenger

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

#ifndef _GENISOPATCH_
#define _GENISOPATCH_

#include <vector>

#include "ijkMCtable.h"

namespace IJKMCUBE_TABLE {

  /*!
   *  @brief Generate isosurface patch for Marching Cubes lookup table.
   *  @param polytope Mesh polytope.
   *  @param vertex_sign[iv] Vertex sign of vertex iv.  
   *    - -1 = negative; 0 = zero; 1 = positive
   *  @param[out] edge_list Return list of polytope edges
   *    - i'th simplex has vertices on edges
   *      (edge_list[d*i], edge_list[d*i+1], ..., edge_list[d*i+d-1]).
   *  @param[out] num_simplices Number of simplices in isosurface patch.
   */
  void gen_isopatch
    (const ISOTABLE_POLY_BASE & polytope,
     const int * const vertex_sign, 
     std::vector<ISOSURFACE_VERTEX_INDEX> & edge_list, int & num_simplices);

  /*!
   *  @overload
   *  @brief Generate isosurface patch for Marching Cubes lookup table.
   *  - Version that passes temporary arrays for faster processing.
   *  @param[out] temp_list[] Temporary array for point_list[].
   *  @pre Preallocated to size at least
   *    (polytope.NumEdge()+polytope.NumVertices().
   *  @param[out] temp_pcoord[] Temporary array for point coordinates.
   *    @pre Preallocated to size at least
   *    (polytope.NumEdges() + polytope.NumVertices())*polytope.Dimension().
   */
  void gen_isopatch
    (const ISOTABLE_POLY_BASE & polytope, 
     const int * const vertex_sign,
     std::vector<ISOSURFACE_VERTEX_INDEX> & edge_list, int & num_simplices,
     int * temp_plist, double * temp_pcoord);

  /*!
   *  @brief Generate isosurface patch for Marching Cubes edge groups lookup table.
   *  @pre Cube dimension is 3.
   */
  void gen_isopatch_edge_groups
    (const ISOTABLE_CUBE_3D & cube, 
     const int * const vertex_sign,
     std::vector<ISOSURFACE_VERTEX_INDEX> & edge_list, int & num_simplices);

  /*!
   *  @overload
   *  @brief Generate isosurface patch for Marching Cubes edge groups lookup table.
   *  - Version that passes temporary arrays for faster processing.
   *  @pre Cube dimension is 3.
   *  @param[out] temp_list[] Temporary array for point_list[].
   *  @pre Preallocated to size at least
   *    (cube.NumEdge()+cube.NumVertices().
   *  @param[out] temp_pcoord[] Temporary array for point coordinates.
   *  @pre Preallocated to size at least
   *    (cube.NumEdges() + cube.NumVertices())*cube.Dimension().
   */
  void gen_isopatch_edge_groups
    (const ISOTABLE_CUBE_3D & cube, 
     const int * const vertex_sign,
     std::vector<ISOSURFACE_VERTEX_INDEX> & edge_list,
     int & num_simplices,
     int * temp_plist, double * temp_pcoord);

  /*!
   *  @brief Generate isosurface patch for Marching Cubes 
   *    NEP (negative, equals, positive) lookup table.
   *  @param polytope Mesh polytope.
   *  @param vertex_sign[iv] Vertex sign of vertex iv.  
   *    - -1 = negative; 0 = zero; 1 = positive
   *  @param[out] isov_list List of isosurface vertices.
   *    - i'th simplex has vertices
   *      (isov_list[d*i], isov_list[d*i+1], ..., isov_list[d*i+d-1]).
   *  @param[out] isov_type[k] Isosurface vertex type corresponding 
   *    to isov_list[k].
   *  @param[out] num_simplices Number of simplices in isosurface patch.
   */
  void gen_isopatch_nep
    (const ISOTABLE_POLY_BASE & polytope, 
     const int * const vertex_sign,
     std::vector<ISOSURFACE_VERTEX_INDEX> & isov_list, 
     std::vector<ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE> & isov_type,
     int & num_simplices);

  /*!
   *  @overload
   *  @brief Generate isosurface patch for Marching Cubes
   *    NEP (negative, equals, positive) lookup table.
   *  - Version that passes temporary arrays for faster processing.
   *  @param[out] temp_plist[] Temporary array for point_list[].
   *  @pre Preallocated to size at least
   *    (polytope.NumEdge()+polytope.NumVertices().
   *  @param[out] temp_type[] Temporary array for isosurface vertex type.
   *  @pre Preallocated to size at least
   *    (polytope.NumEdge()+polytope.NumVertices().
   *  @param[out] temp_pcoord[] Temporary array for point coordinates.
   *  @pre Preallocated to size at least
   *    (polytope.NumEdges() + polytope.NumVertices())*polytope.Dimension().
   */
  void gen_isopatch_nep
    (const ISOTABLE_POLY_BASE & polytope, 
     const int * const vertex_sign,
     std::vector<ISOSURFACE_VERTEX_INDEX> & isov_list, 
     std::vector<ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE> & isov_type,
     int & num_simplices,
     int * temp_plist,
     ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE * temp_type, 
     double * temp_pcoord);

  /*!
   *  @brief Get isosurface patch boundary.
   *  - Note: Isosurface patch may have more than 
   *    one connected component.
   *  - Note: Boundary simplices have one less vertex
   *    than isosurface patch simplices.
   */
  /* OBSOLETE
   void get_isopatch_boundary
   (const ISOTABLE_POLY_BASE & polytope, 
    const int * const vertex_sign,
    std::vector<ISOSURFACE_VERTEX_INDEX> & boundary_simplex_vertex_list,
    int & num_boundary_simplices);
  */
  
}

#endif
