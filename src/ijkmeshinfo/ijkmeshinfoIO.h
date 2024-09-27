/*!
 *  @file ijkmeshinfoIO.h
 *  @brief IO routines for ijkmeshinfo.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Code
  Copyright (C) 2017-2022 Rephael Wenger

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

#ifndef _IJKMESHINFO_IO_
#define _IJKMESHINFO_IO_

#include <iostream>
#include <string>
#include <vector>

#include "ijkmeshinfo.h"


namespace IJKMESHINFO {


  // **************************************************
  // Class IO_INFO
  // **************************************************

  /// Information controlling input and output.
  class IO_INFO {

  protected:

    void Init();

  public:

    /// Output values less than or equal to angle_le[i].
    std::vector<ANGLE_TYPE> angle_le;

    /// Output values greater than or equal to angle_ge[i].
    std::vector<ANGLE_TYPE> angle_ge;

    /// Output values less than or equal to facet_angle_le.Value().
    IJK::SET_VALUE<ANGLE_TYPE> facet_angle_le;

    /// Output values less than or equal to facet_angle_ge.Value().
    IJK::SET_VALUE<ANGLE_TYPE> facet_angle_ge;

    /// Output values less than or equal to edge_length_ratio_le.Value().
    IJK::SET_VALUE<ANGLE_TYPE> edge_length_ratio_le;

    /// Limit on number of polytopes in output list.
    int max_num_poly_out;

    /// If true, output min angle.
    bool flag_output_min_angle;

    /// If true, output max angle.
    bool flag_output_max_angle;

    /// If true, output min edge length.
    bool flag_output_min_edge_length;

    /// If true, output max edge length.
    bool flag_output_max_edge_length;

    /// If true, output min ratio of min to max edge length in cell.
    bool flag_output_min_edge_length_ratio;

    /// If true, output min of the Jacobian determinants.
    bool flag_output_min_Jacobian_determinant;

    /// If true, output max of the Jacobian determinants.
    bool flag_output_max_Jacobian_determinant;

    /// If true, output min of the normalized Jacobian determinants.
    bool flag_output_min_normalized_Jacobian_determinant;

    /// If true, output max of the normalized Jacobian determinants.
    bool flag_output_max_normalized_Jacobian_determinant;

    /// If true, output min of the shape metrics based 
    ///   on the Jacobian matrices.
    bool flag_output_min_Jacobian_shape;

    /// If true, output max of the shape metrics based 
    ///   on the Jacobian matrices.
    bool flag_output_max_Jacobian_shape;

    /// If true, output all polytopes with minimum and maximum values.
    bool flag_output_all_min_max;

    /// If true, output general information about the mesh.
    bool flag_general_info;

    /// Compute Jacobian at vertices.
    IJK::BOOLEAN_SET_VALUE flag_vJacobian;

    /// Compute Jacobian at polytopes including poly center.
    IJK::BOOLEAN_SET_VALUE flag_pJacobian;

  public:
    /// Constructor
    IO_INFO():edge_length_ratio_le(0),
              flag_vJacobian(false),flag_pJacobian(true)
    { Init(); };
  };


  // **************************************************
  // POLYGON/POLYTOPE NAME ROUTINES
  // **************************************************

  void get_polygon_name
  (const int num_polygon_edges, std::string & polygon_name);


  // **************************************************
  // OUTPUT ROUTINES
  // **************************************************

  /// Output information about polytope ipoly.
  void output_poly_info
  (const int dimension, const POLYMESH_TYPE & polymesh, 
   const COORD_TYPE * vertex_coord, const int poly_index,
   const COORD_TYPE max_small_magnitude);

  /// Output degenerate polytopes.
  void output_degenerate_poly
  (const POLYMESH_TYPE & polymesh, const std::vector<POLY_DATA> & poly_data,
   const MESH_INFO & mesh_info);

  /// Output duplicate polytopes.
  void output_duplicate_poly
  (const POLYMESH_TYPE & polymesh, const std::vector<POLY_DATA> & poly_data,
   const MESH_INFO & mesh_info);

  /// Output polytopes with duplicate vertices.
  void output_duplicate_vertices
  (const int dimension, const POLYMESH_TYPE & polymesh,
   const COORD_TYPE vertex_coord[]);


  // **************************************************
  // OUTPUT POLYGON ANGLE ROUTINES
  // **************************************************

  /// Output the minimum and maximum angles of a set of polygons.
  void output_min_max_polygon_angle
  (const MESH_DATA & mesh_data,  const COORD_TYPE * vertex_coord, 
   const IO_INFO & io_info, 
   const bool flag_internal, const bool terse_flag,
   const int num_poly_edges);

  /// Output polygons with minimum angle.
  void output_polygons_with_min_angle
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const IO_INFO & io_info, const bool flag_internal,
   const bool flag_output_vertex_coord);

  /// Output polygons with maximum angle
  void output_polygons_with_max_angle
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const IO_INFO & io_info, const bool flag_internal,
   const bool flag_output_vertex_coord);

  // Output polygons with minimum and maximum angles
  void output_polygons_with_min_max_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const IO_INFO & io_info, const bool flag_internal,
   const bool flag_output_vertex_coord);

  /// Output polygons with angles at most angle_bound.
  void output_polygons_with_small_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const ANGLE_TYPE angle_bound);

  /// Output polygons with angles at most angle_bound[i] for each i.
  /// - Version with array angle_bound[].
  void output_polygons_with_small_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const std::vector<ANGLE_TYPE> & angle_bound);

  /// Output polygons with angles at least angle_bound.
  void output_polygons_with_large_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, ANGLE_TYPE angle_bound);

  /// Output polygons with angles at least angle_bound[i] for each i.
  /// - Version with array angle_bound[].
  void output_polygons_with_large_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const std::vector<ANGLE_TYPE> & angle_bound);


  // **************************************************
  // OUTPUT TETRAHEDRA FACET ANGLE ROUTINES
  // **************************************************

  /// Output the minimum and maximum angles of the set of triangles
  ///   which are the tetrahedra facets.
  void output_min_max_tetrahedra_facet_angle
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const IO_INFO & io_info, const bool flag_internal);

  /// Output number of tetrahedra facets with angles less than 
  ///   or greater than given values.
  void output_tetrahedra_facet_angle_count
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with minimum facet angle.
  void output_tetrahedra_with_min_facet_angle
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with maximum facet angle.
  void output_tetrahedra_with_max_facet_angle
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with minimum and maximum facet angles.
  void output_tetrahedra_with_min_max_facet_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with small facet angles.
  void output_tetrahedra_with_small_facet_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const ANGLE_TYPE angle_bound);

  /// Output tetrahedra with large facet angles.
  void output_tetrahedra_with_large_facet_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const ANGLE_TYPE angle_bound);


  // **************************************************
  // OUTPUT DIHEDRAL ANGLE ROUTINES
  // **************************************************

  void output_min_max_dihedral_angle
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output number of tetrahedra with angles less than or greater than
  ///   given values.
  void output_dihedral_angle_count
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with minimum dihedral angle.
  void output_tetrahedra_with_min_dihedral_angle
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with maximum dihedral angle.
  void output_tetrahedra_with_max_dihedral_angle
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with minimum and maximum dihedral angles.
  void output_tetrahedra_with_min_max_dihedral_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with small dihedral angles.
  void output_tetrahedra_with_small_dihedral_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const ANGLE_TYPE min_angle);

  /// Output tetrahedra with large dihedral angles.
  /// - Version with array angle_bound[].
  void output_tetrahedra_with_small_dihedral_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const std::vector<ANGLE_TYPE> & angle_bound);

  /// Output tetrahedra with large dihedral angles.
  void output_tetrahedra_with_large_dihedral_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const ANGLE_TYPE min_angle);

  /// Output tetrahedra with large dihedral angles.
  /// - Version with array angle_bound[].
  void output_tetrahedra_with_large_dihedral_angles
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const std::vector<ANGLE_TYPE> & angle_bound);


  // **************************************************
  // OUTPUT EDGE LENGTH ROUTINES
  // **************************************************

  /// Output min/max lengths of edges in array mesh_data.edge_data.
  void output_min_max_edge_lengths_using_edge_list
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal_edge);

  /// Output minimum/maximum polygon edge lengths.
  void output_min_max_polygon_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, 
   const bool flag_internal, const bool terse_flag,
   const int num_poly_edges);

  /// Output minimum/maximum tetrahedra edge lengths.
  void output_min_max_tetrahedra_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output minimum/maximum simplices edge lengths.
  void output_min_max_simplices_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output minimum/maximum hexahedra edge lengths.
  void output_min_max_hexahedra_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output polygons with minimum edge lengths.
  void output_polygons_with_min_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output polygons with max edge lengths.
  void output_polygons_with_max_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output polygons with minimum and maximum edge lengths.
  void output_polygons_with_min_max_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with minimum edge lengths.
  void output_tetrahedra_with_min_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output tetrahedra with maximum edge lengths.
  void output_tetrahedra_with_max_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output hexahedra with mininum edge lengths.
  void output_hexahedra_with_min_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output hexahedra with maximum edge lengths.
  void output_hexahedra_with_max_edge_lengths
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);


  // **************************************************
  // OUTPUT EDGE LENGTH RATIO ROUTINES
  // **************************************************

  /// Output minimum polygon edge length ratios.
  void output_min_polygon_edge_length_ratios
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, 
   const bool flag_internal, const bool terse_flag,
   const int num_poly_edges);

  // Output polygons with small ratios of shortest to longest edge lengths.
  void output_polygons_with_small_edge_length_ratios
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord, 
   const bool flag_internal, const COORD_TYPE ratio_bound);


  // **************************************************
  // OUTPUT JACOBIAN ROUTINES
  // **************************************************

  void output_min_max_hexahedra_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant);

  /// Output minimum and maximum hexahedra normalized 
  ///   Jacobian matrix determinants.
  void output_min_max_hexahedra_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant);

  /// Output minimum and maximum of the Jacobian matrix determinants
  ///   at the hex vertices (not at hex centers.)
  void output_min_max_hex_vert_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
   const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, 
   const bool flag_internal_poly,
   const bool flag_internal_vertex,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant);

  /// Output minimum and maximum of the normalized Jacobian matrix determinants
  ///   at the hex vertices (not at hex centers.)
  void output_min_max_hex_vert_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
   const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, 
   const bool flag_internal_poly,
   const bool flag_internal_vertex,
   COORD_TYPE & min_Jacobian_determinant, 
   COORD_TYPE & max_Jacobian_determinant);

  /// Output hexahedra with min Jacobian determinants.
  void output_hexahedra_with_min_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output hexahedra with max Jacobian determinants.
  void output_hexahedra_with_max_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output hexahedra with min normalized Jacobian determinants.
  void output_hexahedra_with_min_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output hexahedra with max normalized Jacobian determinants.
  void output_hexahedra_with_max_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output hex mesh vertices with min Jacobian determinants.
  void output_hex_mesh_vertices_with_min_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
   const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, 
   const bool flag_internal_poly,
   const bool flag_internal_vertex);

  /// Output hex mesh vertices with min normalized Jacobian determinants.
  void output_hex_mesh_vertices_with_min_normalized_Jacobian_determinants
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
   const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, 
   const bool flag_internal_poly,
   const bool flag_internal_vertex);


  // **************************************************
  // OUTPUT SHAPE METRIC BASED ON JACOBIAN MATRICES
  // **************************************************

  void output_min_max_hexahedra_Jacobian_shape
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal,
   COORD_TYPE & min_Jacobian_shape,
   COORD_TYPE & max_Jacobian_shape);

  void output_min_max_hex_vert_Jacobian_shape
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
   const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, 
   const bool flag_internal_poly,
   const bool flag_internal_vertex,
   COORD_TYPE & min_Jacobian_shape, 
   COORD_TYPE & max_Jacobian_shape);

  /// Output hexahedra with min Jacobian shape.
  void output_hexahedra_with_min_Jacobian_shape
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output hexahedra with max Jacobian shape.
  void output_hexahedra_with_max_Jacobian_shape
  (const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, const bool flag_internal);

  /// Output hex mesh vertices with min Jacobian shape metrics.
  void output_hex_mesh_vertices_with_min_Jacobian_shape
  (const MESH_DATA & mesh_data,
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
   const COORD_TYPE * vertex_coord,
   const IO_INFO & io_info, 
   const bool flag_internal_poly,
   const bool flag_internal_vertex);


  // **************************************************
  // PRINT ROUTINES
  // **************************************************

  void print_poly_index_and_vert
  (std::ostream & out, 
   const char * prefix, const char * separator, const char * suffix,
   const POLYMESH_TYPE & polymesh, const int ipoly);


  // **************************************************
  // USAGE/HELP MESSAGES
  // **************************************************

  void usage_msg();
  void usage_error();
  void help_msg();

};

#endif
