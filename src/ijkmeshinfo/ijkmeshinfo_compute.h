/*!
 *  @file ijkmeshinfo_compute.h
 *  @brief Compute angles, edge lengths, Jacobian determinants, etc.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2017-2024 Rephael Wenger

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

#include "ijkmeshinfo.h"

#ifndef _IJKMESHINFO_COMPUTE_
#define _IJKMESHINFO_COMPUTE_

namespace IJKMESHINFO {

  // **************************************************
  //! @name Compute angle routines
  // **************************************************

  //@{

  /// @brief Compute angle (in degrees) at k'th vertex of polygon ipoly.
  void compute_angle_in_degrees_at_polygon_vertex
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const int ipoly, const int k,
   const COORD_TYPE max_small_magnitude,
   ANGLE_TYPE & angle, bool & flag_identical);

  /*!
   *  @brief Compute min/max polygon angles.
   *  @pre mesh_dimension == 2.
   *  @param flag_internal If true, compute angles for interior polygons.
   *  @param num_poly_edges If num_poly_edges > 0, compute angles only
   *           for polygons with num_poly_edges.
   *  @param num_poly_edges If num_poly_edges = 0, compute angles 
   *           for all polygons.
   *  @param[out] min_angle Minimum polygon angle.
   *  @param[out] max_angle Maximum polygon angle.
   *  @param[out] poly_with_min_angle Index of polygon with min angle.
   *  @param[out] poly_with_max_angle Index of polygon with max angle.
   */
  void compute_min_max_polygon_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal, const int num_poly_edges,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
   int & poly_with_min_angle, int & poly_with_max_angle);

  /*!
   *  @overload
   *  @brief Compute min/max polygon angles. 
   *    (No restriction on number of poly edges.)
   *  - Version with num_poly_edges set to 0.
   */
  void compute_min_max_polygon_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
   int & poly_with_min_angle, int & poly_with_max_angle);

  /*!
   *  @overload
   *  @brief Compute min/max polygon angles. (No return poly.)
   *  - Version without return arguments poly_with_min_angle and
   *    poly_with_max_angle.
   */
  void compute_min_max_polygon_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal, const int num_poly_edges,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle);

  /*!
   *  @overload
   *  @brief Compute min/max polygon angles. 
   *    (No return poly and no restriction on number of poly edges.)
   *  - Version with num_poly_edges set to 0.
   *  - Version without return arguments poly_with_min_angle and
   *    poly_with_max_angle.
   */
  void compute_min_max_polygon_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal, ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle);

  /// @brief Compute min/max polygon angles.
  void compute_min_max_polygon_angles
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX poly_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle, int & num_angle);

  /*!
   *  @brief Compute min/max cos angles of a single polygon.
   *  - Note: acos(min_cos) is the max angle.<br>
   *          acos(max_cos) is the min angle.
   */
  void compute_min_max_cos_polygon_angles
  (const MESH_DATA & mesh_data, 
   const VERTEX_INDEX poly_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & cos_min, COORD_TYPE & cos_max, int & num_angle);

  /*!
   *  @brief Compute number of polygons with angles less than or equal 
   *    to angle_le[j] and greater than or equal to angle_le[j] for each j.
   *  @pre mesh_dimension == 2.
   *  @param[out] num_le[j] Number of polygons with angle less than
   *    or equal to angle_le[j].
   *  @param[out] num_ge[j] Number of polygons with angle greater than
   *    or equal to angle_ge[j].
   */
  void compute_num_polygon_angles
  (const int dimension,
   const POLYMESH_TYPE & polymesh,
   const std::vector<POLY_DATA> & poly_data,
   const COORD_TYPE * vertex_coord,
   const bool flag_internal, 
   const std::vector<ANGLE_TYPE> & angle_le, 
   const std::vector<ANGLE_TYPE> & angle_ge,
   const COORD_TYPE max_small_magnitude,
   std::vector<int> & num_le, std::vector<int> & num_ge);

  /*!
   *  @brief Compute cosine of min/max  angle over all tetrahedron facets.
   *  @pre Polytope ipoly in polymesh is a tetrahedron.
   */
  void compute_cos_min_max_tetrahedron_facet_angles
  (const int dimension,
   const VERTEX_INDEX tetrahedron_vert[], const COORD_TYPE * vertex_coord,
   COORD_TYPE & cos_min, COORD_TYPE & cos_max,
   int & num_angle);

  /*!
   *  @brief Compute min/max angle over all tetrahedron facets.
   *  @param num_vert Number of vertices.  
   *     This value is ignored since it should always be 4.
   *  @pre num_vert equals 4.
   */
  void compute_min_max_tetrahedron_facet_angles
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX tetrahedron_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
   int & num_angle);

  /*!
   *  @brief Compute min/max angle over all tetrahedra facets.
   *  @pre Tetrahedra are in 3D.
   */
  void compute_min_max_tetrahedra_facet_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
   int & poly_with_min_angle, int & poly_with_max_angle);

  /*!
   *  @overload
   *  @brief Compute min/max angle over all tetrahedra facets.
   *    (No return poly.)
   *  - Version which does not return poly_with_min_angle or
   *    poly_with_max_angle.
   */
  void compute_min_max_tetrahedra_facet_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle);

  /*!
   *  Compute number of tetrahedra facet angles less than or equal 
   *    to min_angle and  greater than or equal to max_angle.
   *  @pre mesh_dimension == 3.
   */
  void compute_num_tetrahedra_facet_angles
  (const int dimension, const int mesh_dimension, 
   const POLYMESH_TYPE & polymesh, const std::vector<POLY_DATA> & poly_data,
   const COORD_TYPE * vertex_coord,
   const bool flag_internal, 
   const ANGLE_TYPE min_angle, const ANGLE_TYPE max_angle,
   int & num_le, int & num_ge);

  //@}


  // **************************************************
  //! @name Compute dihedral angle routines
  // **************************************************

  //@{

  /// @brief Compute min/max dihedral angles of tetrahedra.
  /// @pre polymesh is a mesh of tetrahedra.
  void compute_min_max_tetmesh_dihedral_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle, 
   int & poly_with_min_angle, int & poly_with_max_angle);

  /*!
   *  @overload
   *  @brief Compute min/max dihedral angles of tetrahedra.
   *  - Version without return arguments poly_with_min_angle and
   *    poly_with_max_angle.
   */
  void compute_min_max_tetmesh_dihedral_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle);

  /*!
   *  @brief Compute number of dihedral angles less than or equal to min_angle and
   *    greater than or equal to max_angle.
   *  @pre mesh_dimension == 3.
   */
  void compute_num_tetmesh_dihedral_angles
  (const int dimension,
   const int mesh_dimension, 
   const POLYMESH_TYPE & polymesh,
   const std::vector<POLY_DATA> & poly_data,
   const COORD_TYPE * vertex_coord,
   const bool flag_internal, 
   const std::vector<ANGLE_TYPE> & angle_le, 
   const std::vector<ANGLE_TYPE> & angle_ge,
   std::vector<int> & num_le, std::vector<int> & num_ge);

  /// @brief Compute min/max dihedral angles of a tetrahedron.
  /// @param num_vert Number of vertices.
  ///    This value is ignored since it should always be 4.
  /// @pre num_vert equals 4.
  void compute_min_max_tetrahedron_dihedral_angles
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX poly_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle, int & num_val);

  //@}


  // **************************************************
  //! @name Compute edge length routines
  // **************************************************

  //@{

  /// @brief Compute edge length of edge ie.
  void compute_edge_length
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const int ie, COORD_TYPE & edge_length);

  /// @brief Compute min/max edge lengths of edges in mesh_data.edge_data.
  void compute_min_max_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal_edge,
   COORD_TYPE & min_edge_length, COORD_TYPE & max_edge_length,
   int & edge_with_min_length, int & edge_with_max_length);

  /// @brief Compute min/max edge lengths of edges in mesh_data.edge_data.
  /// - Version which does not return edges with min/max length.
  void compute_min_max_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal_edge,
   COORD_TYPE & min_edge_length, COORD_TYPE & max_edge_length);

  /// @brief Compute min/max edge lengths of polgyons.
  void compute_min_max_polygon_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_length, COORD_TYPE & max_length,
   int & poly_with_min_edge_length, int & poly_with_max_edge_length);

  /// @brief Compute min/max edge lengths of polygons.
  void compute_min_max_polygon_edge_lengths_select_poly_by_numv
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal, const int num_poly_vert,
   COORD_TYPE & min_length, COORD_TYPE & max_length,
   int & poly_with_min_edge_length, int & poly_with_max_edge_length);

  /// @brief Compute min/max edge lengths of polygons.
  /// - Version which does not return polygons with min/max edge lengths.
  void compute_min_max_polygon_edge_lengths_select_poly_by_numv
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal, const int num_poly_vert,
   COORD_TYPE & min_length, COORD_TYPE & max_length);

  /// @brief Compute min/max edge lengths of selected polygons
  ///   with number of edges equal to num_poly_edges.
  void compute_min_max_selected_polygon_edge_lengths
  (const MESH_DATA & mesh_data,
   const POLYMESH_TYPE & polymesh,
   const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   const int num_poly_edges,
   COORD_TYPE & min_length, COORD_TYPE & max_length,
   int & poly_with_min_edge_length, int & poly_with_max_edge_length);

  /// @brief Compute min/max edge lengths of one polygon.
  void compute_min_max_polygon_edge_lengths
  (const int dimension,
   const VERTEX_INDEX poly_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_length, COORD_TYPE & max_length, int & num_lengths);

  /// @brief Compute min/max edge lengths of one polygon.
  /// - Version using MESH_DATA.
  void compute_min_max_polygon_edge_lengths
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX poly_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_length, COORD_TYPE & max_length, int & num_lengths);

  /// @brief Compute min/max edge lengths of tetrahedra.
  void compute_min_max_tetrahedra_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_length, COORD_TYPE & max_length,
   int & poly_with_min_edge_length, int & poly_with_max_edge_length);

  /// @brief Compute min/max edge lengths of one tetrahedron.
  void compute_min_max_tetrahedron_edge_lengths
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX tet_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_length, COORD_TYPE & max_length, int & num_lengths);

  /// @brief Compute min/max edge lengths of a polymesh of simplices.
  void compute_min_max_simplices_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_length, COORD_TYPE & max_length,
   int & poly_with_min_edge_length, int & poly_with_max_edge_length);

  /// @brief Compute min/max edge lengths of one simplex.
  void compute_min_max_simplex_edge_lengths
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX simplex_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_length, COORD_TYPE & max_length,
   int & num_lengths);

  /// @brief Compute min/max edge lengths of hexahedra.
  void compute_min_max_hexahedra_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_length, COORD_TYPE & max_length,
   int & poly_with_min_edge_length, int & poly_with_max_edge_length);

  /// @brief Compute min/max edge lengths of one hexahedron.
  void compute_min_max_hexahedron_edge_lengths
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX hex_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_length, COORD_TYPE & max_length,
   int & num_lengths);

  //@}


  // **************************************************
  //! @name Compute edge length ratio routines
  // **************************************************

  //@{

  /// @brief Compute ratio of shortest to largest edge in a single polygon.
  void compute_polygon_edge_length_ratio
  (const int dimension,
   const VERTEX_INDEX poly_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & edge_length_ratio,
   int & num_lengths);

  /// @brief Compute ratio of shortest to largest edge in a single polygon.
  /// - Version using MESH_DATA.
  void compute_polygon_edge_length_ratio
  (const MESH_DATA & mesh_data, 
   const VERTEX_INDEX poly_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & edge_length_ratio,
   int & num_lengths);

  /// @brief Compute min edge length ratios of polygons.
  void compute_min_polygon_edge_length_ratios
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_ratio, int & poly_with_min_edge_length_ratio);

  /// @brief Compute min edge length ratios of polygons.
  void compute_min_polygon_edge_length_ratios_select_poly_by_numv
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal, const int num_poly_vert,
   COORD_TYPE & min_length, int & poly_with_min_edge_length);

  /// @brief Compute min edge length ratios of polygons.
  /// - Version which does not return polygon with min edge length ratio.
  void compute_min_polygon_edge_length_ratios_select_poly_by_numv
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal, const int num_poly_vert,
   COORD_TYPE & min_ratio);

  /// @brief Compute number of polygons with edge length ratio at most
  ///   max_edge_length_ratio.
  /// @pre mesh_dimension == 2.
  void compute_num_polygons_with_small_edge_length_ratios
  (const int dimension,
   const POLYMESH_TYPE & polymesh,
   const std::vector<POLY_DATA> & poly_data,
   const COORD_TYPE * vertex_coord,
   const bool flag_internal, 
   const COORD_TYPE ratio_bound, 
   int & num_le);

  //@}

  
  // **************************************************
  //! @name Compute Jacobian determinants
  // **************************************************

  //@{

  /// @brief Compute min/max Jacobian matrix determinants of hexahedra.
  /// @pre dimension = 3. 
  void compute_min_max_hexahedra_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant,
   int & poly_with_min_Jacobian_determinant, 
   int & poly_with_max_Jacobian_determinant);

  /// @brief Compute min/max Jacobian matrix determinants of hexahedra.
  /// - Version which does not return hexahedra containing min/max 
  ///   Jacobian determinants.
  void compute_min_max_hexahedra_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant);

  /// @brief Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
  /// - Eight Jacobian matrices at the eight hexahedron vertices and 
  ///     one Jacobian matrix at the hexahedron center.
  /// @pre dimension = 3. 
  void compute_min_max_hexahedron_Jacobian_determinants
  (const MESH_DATA & mesh_data, 
   const VERTEX_INDEX hex_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_Jacobian_determinant,
   COORD_TYPE & max_Jacobian_determinant,
   int & num_Jacobian_determinants);

  /// @brief Compute min/max Jacobian matrix determinants of hexahedra vertices.
  /// - Version which returns vertices with min/max Jacobian determinants.
  void compute_min_max_hex_vert_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const bool flag_internal_poly,
   const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
   int & poly_with_min_Jacobian_determinant, 
   int & poly_with_max_Jacobian_determinant,
   int & vert_with_min_Jacobian_determinant, 
   int & vert_with_max_Jacobian_determinant);

  /// @brief Compute min/max of the Jacobian matrix determinants of hexahedra vertices.
  /// - Version which does not return vertices or polytopes with min/max values.
  void compute_min_max_hex_vert_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const bool flag_internal_poly, const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant);

  /// @brief Compute min/max Jacobian matrix determinants of a vertex 
  ///   in a hexahedral mesh.
  /// - Version with input vertex-poly incidence data structure.
  void compute_min_max_hex_vert_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const VERTEX_INDEX iv0,
   const bool flag_internal_poly,
   const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant,
   int & poly_with_min_Jacobian_determinant, 
   int & poly_with_max_Jacobian_determinant,
   int & num_Jacobian_determinants);

  //@}


  // **************************************************
  //! @name Compute normalized Jacobian determinants
  // **************************************************

  //@{

  /// @brief Compute min/max normalized Jacobian matrix determinants of hexahedra.
  void compute_min_max_hexahedra_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
   int & poly_with_min_Jacobian_determinant, 
   int & poly_with_max_Jacobian_determinant);

  /// @brief Compute min/max normalized Jacobian matrix determinants of hexahedra.
  /// - Version which does not return hexahedra containing min/max 
  ///   normalized Jacobian determinants.
  void compute_min_max_hexahedra_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant);

  /// @brief Compute min/max of the nine normalized Jacobian matrix determinants 
  ///   of a hexahedron.
  /// @pre dimension = 3. 
  void compute_min_max_hexahedron_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX hex_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_Jacobian_determinant,
   COORD_TYPE & max_Jacobian_determinant,
   int & num_Jacobian_determinants);

  /// @brief Compute min/max normalized Jacobian matrix determinants 
  ///   of hexahedra vertices.
  /// - Version which returns vertices with min/max Jacobian determinants.
  void compute_min_max_hex_vert_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const bool flag_internal_poly,
   const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant,
   int & poly_with_min_Jacobian_determinant, 
   int & poly_with_max_Jacobian_determinant,
   int & vert_with_min_Jacobian_determinant, 
   int & vert_with_max_Jacobian_determinant);

  /// @brief Compute min/max of the normalized Jacobian matrix determinants 
  ///   of hexahedra vertices.
  /// - Version which does not return vertices or polytopes with min/max values.
  void compute_min_max_hex_vert_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const bool flag_internal_poly, const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant);

  /// @brief Compute min/max normalized Jacobian matrix determinants of a vertex 
  ///   in a hexahedral mesh.
  /// - Version with input vertex-poly incidence data structure.
  void compute_min_max_hex_vert_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const VERTEX_INDEX iv0,
   const bool flag_internal_poly,
   const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant,
   int & poly_with_min_Jacobian_determinant, 
   int & poly_with_max_Jacobian_determinant,
   int & num_Jacobian_determinants);

  /// @brief Compute min/max of the eight normalized Jacobian matrix determinants 
  ///   at the eight vertices of a hexahedron.
  /// @pre dimension = 3. 
  void compute_min_max_hex_vert_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX hex_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_Jacobian_determinant,
   COORD_TYPE & max_Jacobian_determinant,
   VERTEX_INDEX & vert_with_min_Jacobian_determinant,
   VERTEX_INDEX & vert_with_max_Jacobian_determinant,
   int & num_Jacobian_determinants);

  /// @brief Compute min/max normalized Jacobian matrix determinants of a vertex 
  ///   in a hexahedral mesh.
  /// - Version with input vertex-poly incidence data structure.
  void compute_min_max_hex_vert_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const VERTEX_INDEX iv0,
   const bool flag_internal_poly,
   const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant,
   int & poly_with_min_Jacobian_determinant, 
   int & poly_with_max_Jacobian_determinant,
   int & num_Jacobian_determinants);

  //@}
  

  // ***********************************************************************
  //! @name Compute shape metrics based on Jacobian matrices.
  // ***********************************************************************

  //@{
  
  /// @brief Compute min/max hexahedra shape metric based on Jacobian matrices.
  void compute_min_max_hexahedra_Jacobian_shape
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape,
   int & poly_with_min_Jacobian_shape, 
   int & poly_with_max_Jacobian_shape);

  /// @brief Compute min/max hexahedra shape metric based on Jacobian matrices.
  /// - Version which does not return hexahedra containing min/max shape.
  void compute_min_max_hexahedra_Jacobian_shape
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const bool flag_internal,
   COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape);

  /// @brief Compute min/max hexahedron shape metric based on Jacobian matrices.
  /// @pre dimension = 3. 
  void compute_min_max_hexahedron_Jacobian_shape
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX hex_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_Jacobian_determinant,
   COORD_TYPE & max_Jacobian_determinant,
   int & num_Jacobian_determinants);

  /// @brief Compute min/max of the eight Jacobian shape metrics
  ///   at the eight vertices of a hexahedron.
  /// @pre dimension = 3. 
  void compute_min_max_hex_vert_Jacobian_shape
  (const MESH_DATA & mesh_data,
   const VERTEX_INDEX hex_vert[], const int num_vert,
   const COORD_TYPE * vertex_coord,
   COORD_TYPE & min_Jacobian_shape,
   COORD_TYPE & max_Jacobian_shape,
   VERTEX_INDEX & vert_with_min_Jacobian_shape,
   VERTEX_INDEX & vert_with_max_Jacobian_shape,
   int & num_Jacobian_shapes);

  /// @brief Compute min/max shape metrics based on Jacobian matrices at hex vertices.
  /// - Version which returns vertices with min/max shape metric.
  void compute_min_max_hex_vert_Jacobian_shape
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const bool flag_internal_poly,
   const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_shape, 
   COORD_TYPE & max_Jacobian_shape,
   int & poly_with_min_Jacobian_shape, 
   int & poly_with_max_Jacobian_shape,
   int & vert_with_min_Jacobian_shape, 
   int & vert_with_max_Jacobian_shape);

  /// @brief Compute min/max shape metrics based on Jacobian matrices at hex vertices.
  /// - Version which does not return vertices or polytopes with min/max values.
  void compute_min_max_hex_vert_Jacobian_shape
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const bool flag_internal_poly, const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape);

  /// @brief Compute min/max shape metrics based on Jacobian matrices at hex vertices.
  /// - Version with input vertex-poly incidence data structure.
  void compute_min_max_hex_vert_Jacobian_shape
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
   const COORD_TYPE * vertex_coord,
   const VERTEX_INDEX iv0,
   const bool flag_internal_poly,
   const bool flag_internal_vert,
   COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape,
   int & poly_with_min_Jacobian_shape, 
   int & poly_with_max_Jacobian_shape,
   int & num_Jacobian_shapes);

  //@}


  // ***********************************************************************
  //! @name Compute polytope intersections.
  // ***********************************************************************

  //@{

  /// @brief Return true if x0 equals x1, or x0 or x1 equal zero.
  template <typename TYPEA>
  inline bool are_equal_or_zero(const TYPEA x0, const TYPEA x1)
  {
    if ((x0 == x1) || (x0 == 0) || (x1 == 0))
      { return(true); }
    else
      { return(false); }
  }
  
  /// @brief Compute min/max coordinates of triangle vertices.
  void compute_min_max_coord
  (const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
   const COORD_TYPE p2[DIM3],
   COORD_TYPE min_coord[DIM3], COORD_TYPE max_coord[DIM3]);

  /// @brief Compute orientation of (p1-p0), (p2-p0) and w.
  int compute_orientation
  (const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
   const COORD_TYPE p2[DIM3], const COORD_TYPE w[DIM3],
   const COORD_TYPE epsilon);
  
  /*!
   *  @brief Compute affine subspace spanning points p0, p1, p2.
   *  @param span_dimension Dimension of spanning subspace. 0, 1 or 2.
   *  @param w If span_dimension = 2, w is unit normal to the plane.
   *  @param w If span_dimension = 1, w is unit line direction.
   *  @param w If span_dimension = 0, w is p0.
   */
  void compute_span
  (const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
   const COORD_TYPE p2[DIM3], 
   int & span_dimension, COORD_TYPE w[DIM3]);

  /*!
   *  @brief Return true if triangle contains point.
   *  - Assumes plane containes point and triangle.
   *  @param normal[] Normal to plane containing point and triangle.
   */
  bool does_triangle_contain_point
  (const COORD_TYPE w0[DIM3], const COORD_TYPE w1[DIM3],
   const COORD_TYPE w2[DIM3], const COORD_TYPE p[DIM3],
   const COORD_TYPE normal[DIM3], const COORD_TYPE epsilon);

  /*!
   *  @brief Compute intersection of edge (v0,v1) and triangle (w0,w1,w2).
   *  - Use vectors (w1-w0) and (w2-w0)
   *  - Return true if edge intersects triangle.
   */
  bool intersect_edge_triangle
  (const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
   const COORD_TYPE w0[DIM3], const COORD_TYPE w1[DIM3],
   const COORD_TYPE w2[DIM3], const COORD_TYPE epsilon,
   COORD_TYPE intersection_point[DIM3]);
  
  /// @brief Return true if edges of triangle js1 intersect js2.
  bool intersect_triangle_edges_and_triangle
  (const VERTEX_INDEX triangle_vert[],
   const COORD_TYPE vertex_coord[],
   const VERTEX_INDEX vertex_merge[],
   const int js1, const int js2, const double epsilon,
   COORD_TYPE intersection_point[DIM3]);

  /*!
   *  @brief Return true if triangles intersect.
   *  - Computes intersection between every edge of one triangle
   *    with the other triangle.
   *  - Ignores edges of one triangle that have an endpoint
   *    at the same coordinates as some vertex of the other triangle.
   *  - Some algorithms do less comparisons, but comparing every edge
   *    of one triangle with the other triangle is possibly more robust.
   *  @param vertex_merge[iw] Representative index of set 
   *    of vertices with same coordinates as vertex iw.
   */
  bool intersect_triangle_triangle
  (const VERTEX_INDEX triangle_vert[],
   const COORD_TYPE vertex_coord[],
   const VERTEX_INDEX vertex_merge[],
   const int jt0, const int jt1, const double epsilon,
   COORD_TYPE intersection_point[DIM3]);

  /*!
   *  @overload
   *  @brief Return true if triangles intersect. (C++ STL vector.)
   *  - Version using C++ STL vector for vertex_merge[].
   */
  bool intersect_triangle_triangle
  (const VERTEX_INDEX simplex_vert[],
   const COORD_TYPE vertex_coord[],
   const std::vector<VERTEX_INDEX> & vertex_merge,
   const int js1, const int js2, const double epsilon,
   COORD_TYPE intersection_point[DIM3]);

  //@}

  
  // ***********************************************************************
  //! @name Compute mesh self intersections.
  // ***********************************************************************

  //@{

  /*!
   *  @brief Compute self intersections of triangle mesh. Simple (O(n^2)) algorithm.
   *  - Triangle mesh is embedded in R3.
   *  @param triangle_vert[i*3+j] j'th vertex of i'th triangle.
   *  @param num_triangles Number of triangles.
   *  @param vertex_coord[iw*DIM3+j] j'th coordinate of vertex iw.
   *  @param vertex_merge[iw] Representative index of set 
   *    of vertices with same coordinates as vertex iw.
   *  @param[out] intersecting_poly[] List of intersecting pairs of polytopes.
   *  @param[out] intersection_coord[] Intersection coordinates.
   *    - intersection_coord[i*dimension+j] =
   *      j'th coordinate of i'th pairs of intersecting polytopes.
   */
  void compute_trimesh_self_intersections_simple
  (const VERTEX_INDEX triangle_vert[],  const int num_triangles,
   const COORD_TYPE vertex_coord[],
   const COORD_TYPE selfI_epsilon,
   const std::vector<VERTEX_INDEX> & vertex_merge,
   INTERSECTING_POLY_ARRAY & intersecting_poly,
   std::vector<COORD_TYPE> & intersection_coord);

  //@}

}

#endif
