/*!
 *  @file ijkmeshinfo_compute.cpp
 *  @brief Compute angles, edge lengths, etc.
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

#define _USE_MATH_DEFINES

#include <cmath>

#include "ijk.tpp"
#include "ijkcoord.tpp"

#include "ijkmeshinfo_compute.h"
#include "ijkmeshinfo_compute.tpp"


// **************************************************
// COMPUTE ANGLE ROUTINES
// **************************************************

// Compute min/max polygon angles.
// @param flag_internal If true, compute angles for interior polygons.
// @param num_poly_vert If num_poly_vert > 0, compute angles only
//          for polygons with num_poly_vert.
// @param num_poly_edges If num_poly_vert = 0, compute angles 
//          for all polygons.
// @pre mesh_dimension == 2.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_vert,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  const int dimension = mesh_data.dimension;
  IJK::PROCEDURE_ERROR error("compute_min_max_polygon_angles");

  poly_with_min_angle = 0;
  poly_with_max_angle = 0;

  if (!check_mesh_dimension<DIM2>(mesh_data, error)) { throw error; }
  if (flag_internal) 
    { if (!check_boundary_facets(mesh_data, error)) { throw error; } }

  COORD_TYPE min_cos, max_cos;
  if (num_poly_vert == 0) {
    // Note: Poly with max angle has min_cos.
    //       Poly with min angle has max_cos.
    compute_min_max_plist_values
      (mesh_data, vertex_coord, flag_internal, 
       1, -1, min_cos, max_cos, poly_with_max_angle, poly_with_min_angle,
       compute_min_max_cos_polygon_angles);
  }
  else {
    // Note: Poly with max angle has min_cos.
    //       Poly with min angle has max_cos.
    compute_min_max_plist_values_select_poly_by_numv
      (mesh_data, vertex_coord, flag_internal, 1, -1, num_poly_vert,
       min_cos, max_cos, poly_with_max_angle, poly_with_min_angle,
       compute_min_max_cos_polygon_angles);
  }

  min_angle = std::acos(max_cos) * 180.0/M_PI;
  max_angle = std::acos(min_cos) * 180.0/M_PI;
}


// Compute min/max polygon angles.
// - Version with num_poly_edges set to 0.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  compute_min_max_polygon_angles
    (mesh_data, vertex_coord, flag_internal, 0, 
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute min/max polygon angles.
// - Version without return arguments poly_with_min_angle and
//   poly_with_max_angle.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_edges,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_polygon_angles
    (mesh_data, vertex_coord, flag_internal, num_poly_edges, 
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute min/max polygon angles.
// - Version without return arguments poly_with_min_angle and
//   poly_with_max_angle.
// - Version with num_poly_edges set to 0.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal, ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  compute_min_max_polygon_angles
    (mesh_data, vertex_coord, flag_internal, 0, min_angle, max_angle);
}


// Compute min/max polygon angles.
void IJKMESHINFO::compute_min_max_polygon_angles
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle, int & num_angle)
{
  const int dimension = mesh_data.dimension;

  // Initialize
  min_angle = 0;
  max_angle = 180;

  COORD_TYPE cos_min, cos_max;
  IJK::compute_cos_min_max_polygon_angles
    (dimension, vertex_coord, poly_vert, num_vert, 
     mesh_data.MaxSmallMagnitude(), cos_min, cos_max, num_angle);

  if (num_angle > 0) {
    min_angle = std::acos(cos_min) * 180.0/M_PI;
    max_angle = std::acos(cos_max) * 180.0/M_PI;
  }
}


// Compute min/max cosine of the angles of a single polygon.
// - Note: acos(min_cos) is the max angle.
//         acos(max_cos) is the min angle.
void IJKMESHINFO::compute_min_max_cos_polygon_angles
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_cos, COORD_TYPE & max_cos, int & num_angle)
{
  const int dimension = mesh_data.dimension;

  IJK::compute_cos_min_max_polygon_angles
    (dimension, vertex_coord, poly_vert, num_vert,
     mesh_data.MaxSmallMagnitude(), max_cos, min_cos, num_angle);
}


// Compute number of angles less than or equal to min_angle and
//   greater than or equal to max_angle.
// @pre mesh_dimension == 2.
void IJKMESHINFO::compute_num_polygon_angles
(const int dimension,
 const POLYMESH_TYPE & polymesh,
 const std::vector<POLY_DATA> & poly_data,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, 
 const std::vector<ANGLE_TYPE> & angle_le, 
 const std::vector<ANGLE_TYPE> & angle_ge,
 const COORD_TYPE max_small_magnitude,
 std::vector<int> & num_le, std::vector<int> & num_ge)
{
  std::vector<ANGLE_TYPE> cos_angle_le;
  std::vector<ANGLE_TYPE> cos_angle_ge;

  // Initialize to zero.
  IJK::init_vector(num_le, angle_le.size(), 0);
  IJK::init_vector(num_ge, angle_ge.size(), 0);

  // Set cos_angle_le and cos_angle_ge.
  for (int j = 0; j < angle_le.size(); j++) 
    { cos_angle_le.push_back(cos(angle_le[j]*M_PI/180.0)); }
  for (int j = 0; j < angle_ge.size(); j++) 
    { cos_angle_ge.push_back(cos(angle_ge[j]*M_PI/180.0)); }

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (flag_internal) {
      if (poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    int num_angle;
    IJK::compute_cos_min_max_polygon_angles
      (dimension, vertex_coord, polymesh.VertexList(ipoly), 
       polymesh.NumPolyVert(ipoly), max_small_magnitude,
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      for (int j = 0; j < angle_le.size(); j++) {
        if (cos_min_i >= cos_angle_le[j]) { num_le[j]++; }
      }

      for (int j = 0; j < angle_ge.size(); j++) {
        if (cos_max_i <= cos_angle_ge[j]) { num_ge[j]++; }
      }
    }
  }

}


// Compute number of dihedral angles less than or equal to min_angle and
//   greater than or equal to max_angle.
// @pre mesh_dimension == 3.
void IJKMESHINFO::compute_num_tetmesh_dihedral_angles
(const int dimension,
 const int mesh_dimension, 
 const POLYMESH_TYPE & polymesh,
 const std::vector<POLY_DATA> & poly_data,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, 
 const std::vector<ANGLE_TYPE> & angle_le, 
 const std::vector<ANGLE_TYPE> & angle_ge,
 std::vector<int> & num_le, std::vector<int> & num_ge)
{
  std::vector<ANGLE_TYPE> cos_angle_le;
  std::vector<ANGLE_TYPE> cos_angle_ge;

  // Initialize to zero.
  IJK::init_vector(num_le, angle_le.size(), 0);
  IJK::init_vector(num_ge, angle_ge.size(), 0);

  // Set cos_angle_le and cos_angle_ge.
  for (int j = 0; j < angle_le.size(); j++) 
    { cos_angle_le.push_back(cos(angle_le[j]*M_PI/180.0)); }
  for (int j = 0; j < angle_ge.size(); j++) 
    { cos_angle_ge.push_back(cos(angle_ge[j]*M_PI/180.0)); }

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (flag_internal) {
      if (poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    int num_angle;
    IJK::compute_cos_min_max_tetrahedron_dihedral_angles
      (dimension, polymesh.VertexList(ipoly), vertex_coord,
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      for (int j = 0; j < angle_le.size(); j++) {
        if (cos_min_i >= cos_angle_le[j]) { num_le[j]++; }
      }

      for (int j = 0; j < angle_ge.size(); j++) {
        if (cos_max_i <= cos_angle_ge[j]) { num_ge[j]++; }
      }
    }
  }

}


// Compute cosine of min/max angle over all tetrahedron facets.
void IJKMESHINFO::compute_cos_min_max_tetrahedron_facet_angles
(const int dimension,
 const VERTEX_INDEX tetrahedron_vert[], const COORD_TYPE * vertex_coord,
 COORD_TYPE & cos_min, COORD_TYPE & cos_max,
 int & num_angle)
{
  const int NUM_VERT_PER_TETRAHEDRON(4);

  cos_min = -1;
  cos_max = 1;
  num_angle = 0;

  for (int i0 = 0; i0 < NUM_VERT_PER_TETRAHEDRON; i0++) {
    const int iv0 = tetrahedron_vert[i0];

    const int i1 = (i0+1)%NUM_VERT_PER_TETRAHEDRON;
    const int i2 = (i0+2)%NUM_VERT_PER_TETRAHEDRON;
    const int i3 = (i0+3)%NUM_VERT_PER_TETRAHEDRON;

    const int iv1 = tetrahedron_vert[i1];
    const int iv2 = tetrahedron_vert[i2];
    const int iv3 = tetrahedron_vert[i3];

    double cos_angle;
    bool flag_duplicate_point;

    // Compute min/max angles of triangle (iv0, iv1, iv3).
    IJK::compute_cos_triangle_angle_coord_list
      (dimension, vertex_coord, iv0, iv1, iv2, cos_angle, 
       flag_duplicate_point);

    if (!flag_duplicate_point) {
      num_angle++;
      if (cos_angle > cos_min) { cos_min = cos_angle; }
      if (cos_angle < cos_max) { cos_max = cos_angle; }
    }
  }
}


// Compute min/max angle over all tetrahedron facets.
void IJKMESHINFO::compute_min_max_tetrahedron_facet_angles
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX tetrahedron_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & num_angle)
{
  const int dimension = mesh_data.dimension;

  // Initialize
  min_angle = 0;
  max_angle = 180;

  COORD_TYPE cos_min, cos_max;
  compute_cos_min_max_tetrahedron_facet_angles
    (dimension, tetrahedron_vert, vertex_coord, cos_min, cos_max, num_angle);

  if (num_angle > 0) {
    min_angle = std::acos(cos_min) * 180.0/M_PI;
    max_angle = std::acos(cos_max) * 180.0/M_PI;
  }
}


// Compute min/max angle over all tetrahedra facets.
// @pre Tetrahedra are in 3D.
void IJKMESHINFO::compute_min_max_tetrahedra_facet_angles
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  const int dimension = mesh_data.dimension;

  COORD_TYPE cos_min = -1;
  COORD_TYPE cos_max = 1;

  poly_with_min_angle = 0;
  poly_with_max_angle = 0;

  int num_angle;
  for (int ipoly = 0; ipoly < mesh_data.NumPoly(); ipoly++) {

    if (mesh_data.poly_data[ipoly].is_degenerate) { continue; }

    if (flag_internal) {
      if (mesh_data.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    compute_cos_min_max_tetrahedron_facet_angles
      (dimension, mesh_data.VertexList(ipoly), vertex_coord, 
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      if (cos_min_i > cos_min) { 
        cos_min = cos_min_i;
        poly_with_min_angle = ipoly;
      }

      if (cos_max_i < cos_max) {
        cos_max = cos_max_i;
        poly_with_max_angle = ipoly;
      }
    }
  }

  min_angle = std::acos(cos_min) * 180.0/M_PI;
  max_angle = std::acos(cos_max) * 180.0/M_PI;
}


// Compute min/max angle over all tetrahedra facets.
void IJKMESHINFO::compute_min_max_tetrahedra_facet_angles
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_tetrahedra_facet_angles
    (mesh_data, vertex_coord, flag_internal,
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute number of tetrahedra facet angles less than or equal 
//   to min_angle and  greater than or equal to max_angle.
// @pre mesh_dimension == 3.
void IJKMESHINFO::compute_num_tetrahedra_facet_angles
(const int dimension, const int mesh_dimension, 
 const POLYMESH_TYPE & polymesh, 
 const std::vector<POLY_DATA> & poly_data,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, 
 const ANGLE_TYPE min_angle, const ANGLE_TYPE max_angle,
 int & num_le, int & num_ge)
{
  // Initialize to zero.
  num_le = 0;
  num_ge = 0;

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  COORD_TYPE cos_min = cos(min_angle*M_PI/180.0);
  COORD_TYPE cos_max = cos(max_angle*M_PI/180.0);
  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (flag_internal) {
      if (poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    int num_angle;
    compute_cos_min_max_tetrahedron_facet_angles
      (dimension, polymesh.VertexList(ipoly), vertex_coord,
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      if (cos_min_i >= cos_min) { num_le++; }
      if (cos_max_i <= cos_max) { num_ge++; }
    }

  }

}


// **************************************************
// COMPUTE DIHEDRAL ANGLE ROUTINES
// **************************************************

// Compute min/max dihedral angles of tetrahedra.
void IJKMESHINFO::compute_min_max_tetmesh_dihedral_angles
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle,
 int & poly_with_min_angle, int & poly_with_max_angle)
{
  const int dimension = mesh_data.dimension;
  COORD_TYPE cos_min = -1;
  COORD_TYPE cos_max = 1;
  poly_with_min_angle = 0;
  poly_with_max_angle = 0;

  /* NOT YET IMPLEMENTED
  if (flag_internal && !are_boundary_facets_identified) {
    error.AddMessage("Programming error.  Need to compute boundary facets.");
    throw error;
  }
  */

  int num_angle;
  for (int ipoly = 0; ipoly < mesh_data.NumPoly(); ipoly++) {

    if (mesh_data.poly_data[ipoly].IsDegenerate()) { continue; }

    if (flag_internal) {
      if (mesh_data.poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE cos_min_i, cos_max_i;
    IJK::compute_cos_min_max_tetrahedron_dihedral_angles
      (dimension, mesh_data.VertexList(ipoly), vertex_coord,
       cos_min_i, cos_max_i, num_angle);

    if (num_angle > 0) {
      if (cos_min_i > cos_min) { 
        cos_min = cos_min_i; 
        poly_with_min_angle = ipoly;
      }

      if (cos_max_i < cos_max) {
        cos_max = cos_max_i;
        poly_with_max_angle = ipoly;
      }
    }
  }

  min_angle = std::acos(cos_min) * 180.0/M_PI;
  max_angle = std::acos(cos_max) * 180.0/M_PI;
}


void IJKMESHINFO::compute_min_max_tetmesh_dihedral_angles
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle)
{
  int poly_with_min_angle, poly_with_max_angle;

  compute_min_max_tetmesh_dihedral_angles
    (mesh_data, vertex_coord, flag_internal,
     min_angle, max_angle, poly_with_min_angle, poly_with_max_angle);
}


// Compute min/max dihedral angles of a tetrahedron.
void IJKMESHINFO::compute_min_max_tetrahedron_dihedral_angles
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 ANGLE_TYPE & min_angle, ANGLE_TYPE & max_angle, int & num_angle)
{
  const int dimension = mesh_data.dimension;

  // Initialize
  min_angle = 0;
  max_angle = 180;

  COORD_TYPE cos_min, cos_max;
  IJK::compute_cos_min_max_tetrahedron_dihedral_angles
    (dimension, poly_vert, vertex_coord, cos_min, cos_max, num_angle);

  if (num_angle > 0) {
    min_angle = std::acos(cos_min) * 180.0/M_PI;
    max_angle = std::acos(cos_max) * 180.0/M_PI;
  }
}


// **************************************************
// COMPUTE EDGE LENGTH ROUTINES
// **************************************************


// Compute edge length of edge ie.
void IJKMESHINFO::compute_edge_length
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const int ie, COORD_TYPE & edge_length)
{
  const int dimension = mesh_data.dimension;
  const VERTEX_INDEX iv0 = mesh_data.edge_data[ie].endpoint[0];
  const VERTEX_INDEX iv1 = mesh_data.edge_data[ie].endpoint[1];
  const COORD_TYPE * v0_coord = &(vertex_coord[iv0*dimension]);
  const COORD_TYPE * v1_coord = &(vertex_coord[iv1*dimension]);

  IJK::compute_distance(dimension, v0_coord, v1_coord, edge_length);
}


// Compute min/max edge lengths of polygons.
void IJKMESHINFO::compute_min_max_edge_lengths
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal_edge,
 COORD_TYPE & min_edge_length, COORD_TYPE & max_edge_length,
 int & edge_with_min_length, int & edge_with_max_length)
{
  COORD_TYPE edge_length;
  IJK::PROCEDURE_ERROR error("compute_min_max_edge_lengths");

  // Initialize.
  min_edge_length = 0;
  max_edge_length = 0;
  edge_with_min_length = 0;
  edge_with_max_length = 0;

  if (mesh_data.edge_data.size() == 0) { 
    if (mesh_data.NumPoly() > 0) {
      error.AddMessage
        ("Programming error. Edge array mesh_data.edge_data not set.");
      throw error;
    }

    // Else no polytopes or edges.
    return; 
  }


  bool is_set = false;
  for (int ie = 0; ie < mesh_data.edge_data.size(); ie++) {

    if (flag_internal_edge) {
      if (mesh_data.edge_data[ie].OnBoundary())
        { continue; } 
    }

    compute_edge_length(mesh_data, vertex_coord, ie, edge_length);

    if (!is_set || edge_length < min_edge_length) {
      min_edge_length = edge_length;
      edge_with_min_length = ie;
    }

    if (!is_set || edge_length > max_edge_length) {
      max_edge_length = edge_length;
      edge_with_max_length = ie;
    }

    is_set = true;
  }

}


// Compute min/max edge lengths of polygons.
void IJKMESHINFO::compute_min_max_edge_lengths
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal_edge,
 COORD_TYPE & min_edge_length, COORD_TYPE & max_edge_length)
{
  int edge_with_min_length, edge_with_max_length;

  compute_min_max_edge_lengths
    (mesh_data, vertex_coord, flag_internal_edge, 
     min_edge_length, max_edge_length, 
     edge_with_min_length, edge_with_max_length);
}

// Compute min/max edge lengths of polygons.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  compute_min_max_plist_values
    (mesh_data, vertex_coord, flag_internal,
     0, 0, min_length, max_length,
     poly_with_min_edge_length, poly_with_max_edge_length,
     compute_min_max_polygon_edge_lengths);
}


// Compute min/max edge lengths of polygons.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths_select_poly_by_numv
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_vert,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  if (num_poly_vert == 0) {
    compute_min_max_polygon_edge_lengths
      (mesh_data, vertex_coord, flag_internal, min_length, max_length,
       poly_with_min_edge_length, poly_with_max_edge_length);
  }
  else {
    compute_min_max_plist_values_select_poly_by_numv
      (mesh_data, vertex_coord, flag_internal,
       0, 0, num_poly_vert, min_length, max_length,
       poly_with_min_edge_length, poly_with_max_edge_length,
       compute_min_max_polygon_edge_lengths);
  }
}


// Compute min/max edge lengths of polygons.
// - Version which does not return polygons with min/max edge lengths.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths_select_poly_by_numv
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_vert,
 COORD_TYPE & min_length, COORD_TYPE & max_length)
{
  int poly_with_min_edge_length, poly_with_max_edge_length;

  compute_min_max_polygon_edge_lengths_select_poly_by_numv
    (mesh_data, vertex_coord, flag_internal, num_poly_vert,
     min_length, max_length, 
     poly_with_min_edge_length, poly_with_max_edge_length);
}

// Compute min/max edge lengths of one polygon.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths
(const int dimension,
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  COORD_TYPE edge_length_squared;

  // Initialize.
  min_length = 0;
  max_length = 0;

  COORD_TYPE min_edge_length_squared = 0;
  COORD_TYPE max_edge_length_squared = 0;
  num_lengths = 0;

  bool flag_set = false;
  for (int i0 = 0; i0 < num_vert; i0++) {
    const int i1 = (i0+1)%num_vert;

    const VERTEX_INDEX iv0 = poly_vert[i0];
    const VERTEX_INDEX iv1 = poly_vert[i1];

    const COORD_TYPE * v0coord = vertex_coord+iv0*dimension;
    const COORD_TYPE * v1coord = vertex_coord+iv1*dimension;

    IJK::compute_distance_squared
      (dimension, v0coord, v1coord, edge_length_squared);

    if (!flag_set || edge_length_squared < min_edge_length_squared) 
      { min_edge_length_squared = edge_length_squared; }
    if (!flag_set || edge_length_squared > max_edge_length_squared) 
      { max_edge_length_squared = edge_length_squared; }

    flag_set = true;
    num_lengths++;
  }

  if (num_lengths > 0) {
    min_length = std::sqrt(min_edge_length_squared);
    max_length = std::sqrt(max_edge_length_squared);
  }
}


// Compute min/max edge lengths of one polygon.
// - Version using MESH_DATA.
void IJKMESHINFO::compute_min_max_polygon_edge_lengths
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;

  compute_min_max_polygon_edge_lengths
    (dimension, poly_vert, num_vert, vertex_coord,
     min_length, max_length, num_lengths);
}


// Compute min/max edge lengths of tetrahedra.
void IJKMESHINFO::compute_min_max_tetrahedra_edge_lengths
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  const int dimension = mesh_data.dimension;

  compute_min_max_plist_values
    (mesh_data, vertex_coord, flag_internal,
     0, 0, min_length, max_length,
     poly_with_min_edge_length, poly_with_max_edge_length,
     compute_min_max_tetrahedron_edge_lengths);
}


// Compute min/max edge lengths of one tetrahedron.
void IJKMESHINFO::compute_min_max_tetrahedron_edge_lengths
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX tet_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord, 
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;
  const int NUM_VERT_PER_TETRAHEDRON(4);
  COORD_TYPE edge_length_squared;

  // Initialize.
  min_length = 0;
  max_length = 0;

  COORD_TYPE min_edge_length_squared = 0;
  COORD_TYPE max_edge_length_squared = 0;
  num_lengths = 0;

  bool flag_set = false;
  for (int i0 = 0; i0 < NUM_VERT_PER_TETRAHEDRON; i0++) {
    for (int i1 = i0+1; i1 < NUM_VERT_PER_TETRAHEDRON; i1++) {

      const VERTEX_INDEX iv0 = tet_vert[i0];
      const VERTEX_INDEX iv1 = tet_vert[i1];

      const COORD_TYPE * v0coord = vertex_coord+iv0*dimension;
      const COORD_TYPE * v1coord = vertex_coord+iv1*dimension;

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, edge_length_squared);

      if (!flag_set || edge_length_squared < min_edge_length_squared) 
        { min_edge_length_squared = edge_length_squared; }
      if (!flag_set || edge_length_squared > max_edge_length_squared) 
        { max_edge_length_squared = edge_length_squared; }

      flag_set = true;
      num_lengths++;
    }
  }

  if (num_lengths > 0) {
    min_length = std::sqrt(min_edge_length_squared);
    max_length = std::sqrt(max_edge_length_squared);
  }
}


// Compute min/max edge lengths of a polymesh of simplices.
void IJKMESHINFO::compute_min_max_simplices_edge_lengths
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  compute_min_max_plist_values
    (mesh_data, vertex_coord, flag_internal,
     0, 0, min_length, max_length,
     poly_with_min_edge_length, poly_with_max_edge_length,
     compute_min_max_simplex_edge_lengths);
}


// Compute min/max edge lengths of one simplex.
void IJKMESHINFO::compute_min_max_simplex_edge_lengths
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX simplex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;
  COORD_TYPE edge_length_squared;

  // Initialize.
  min_length = 0;
  max_length = 0;

  COORD_TYPE min_edge_length_squared = 0;
  COORD_TYPE max_edge_length_squared = 0;
  num_lengths = 0;

  bool flag_set = false;
  for (int i0 = 0; i0 < num_vert; i0++) {
    for (int i1 = i0+1; i1 < num_vert; i1++) {

      const VERTEX_INDEX iv0 = simplex_vert[i0];
      const VERTEX_INDEX iv1 = simplex_vert[i1];

      const COORD_TYPE * v0coord = vertex_coord+iv0*dimension;
      const COORD_TYPE * v1coord = vertex_coord+iv1*dimension;

      IJK::compute_distance_squared
        (dimension, v0coord, v1coord, edge_length_squared);

      if (!flag_set || edge_length_squared < min_edge_length_squared) 
        { min_edge_length_squared = edge_length_squared; }
      if (!flag_set || edge_length_squared > max_edge_length_squared) 
        { max_edge_length_squared = edge_length_squared; }

      flag_set = true;
      num_lengths++;
    }
  }

  if (num_lengths > 0) {
    min_length = std::sqrt(min_edge_length_squared);
    max_length = std::sqrt(max_edge_length_squared);
  }
}


// Compute min/max edge lengths of hexahedra.
void IJKMESHINFO::compute_min_max_hexahedra_edge_lengths
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & poly_with_min_edge_length, int & poly_with_max_edge_length)
{
  compute_min_max_plist_values
    (mesh_data, vertex_coord, flag_internal,
     0, 0, min_length, max_length,
     poly_with_min_edge_length, poly_with_max_edge_length,
     compute_min_max_hexahedron_edge_lengths);
}


// Compute min/max edge lengths of one hexahedron.
void IJKMESHINFO::compute_min_max_hexahedron_edge_lengths
(const MESH_DATA & mesh_data, 
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_length, COORD_TYPE & max_length,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;
  const static CUBE_TYPE cube(DIM3);
  COORD_TYPE edge_length_squared;

  // Initialize.
  min_length = 0;
  max_length = 0;

  COORD_TYPE min_edge_length_squared = 0;
  COORD_TYPE max_edge_length_squared = 0;
  num_lengths = 0;

  bool flag_set = false;
  for (int ie = 0; ie < cube.NumEdges(); ie++) {
    const int iend0 = cube.EdgeEndpoint(ie, 0);
    const int iend1 = cube.EdgeEndpoint(ie, 1);

    const VERTEX_INDEX iv0 = hex_vert[iend0];
    const VERTEX_INDEX iv1 = hex_vert[iend1];

    const COORD_TYPE * v0coord = vertex_coord+iv0*dimension;
    const COORD_TYPE * v1coord = vertex_coord+iv1*dimension;

    IJK::compute_distance_squared
      (dimension, v0coord, v1coord, edge_length_squared);

    if (!flag_set || edge_length_squared < min_edge_length_squared) 
      { min_edge_length_squared = edge_length_squared; }
    if (!flag_set || edge_length_squared > max_edge_length_squared) 
      { max_edge_length_squared = edge_length_squared; }

    flag_set = true;
    num_lengths++;
  }

  if (num_lengths > 0) {
    min_length = std::sqrt(min_edge_length_squared);
    max_length = std::sqrt(max_edge_length_squared);
  }

}


// **************************************************
// COMPUTE EDGE LENGTH RATIO ROUTINES
// **************************************************

// Compute ratio of shortest to largest edge in a single polygon.
void IJKMESHINFO::compute_polygon_edge_length_ratio
(const int dimension,
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & edge_length_ratio,
 int & num_lengths)
{
  COORD_TYPE min_length, max_length;

  compute_min_max_polygon_edge_lengths
    (dimension, poly_vert, num_vert, vertex_coord, 
     min_length, max_length, num_lengths);

  if (num_lengths > 0) {
    if (max_length == 0) {
      edge_length_ratio = 0.0;
    }
    else {
      edge_length_ratio = min_length/max_length;
    }
  }

}


// Compute ratio of shortest to largest edge in a single polygon.
void IJKMESHINFO::compute_polygon_edge_length_ratio
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX poly_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & edge_length_ratio,
 int & num_lengths)
{
  const int dimension = mesh_data.dimension;

  compute_polygon_edge_length_ratio
    (dimension, poly_vert, num_vert, vertex_coord,
     edge_length_ratio, num_lengths);
}

// Compute min edge length ratios of polygons.
void IJKMESHINFO::compute_min_polygon_edge_length_ratios
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_ratio, int & poly_with_min_edge_length_ratio)
{
  compute_min_plist_values
    (mesh_data, vertex_coord, flag_internal, 0, min_ratio,
     poly_with_min_edge_length_ratio,
     compute_polygon_edge_length_ratio);
}


// Compute min edge length ratios of polygons.
void IJKMESHINFO::compute_min_polygon_edge_length_ratios_select_poly_by_numv
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_vert,
 COORD_TYPE & min_ratio, int & poly_with_min_edge_length_ratio)
{
  if (num_poly_vert == 0) {
    compute_min_polygon_edge_length_ratios
      (mesh_data, vertex_coord, flag_internal, min_ratio,
       poly_with_min_edge_length_ratio);
  }
  else {
    compute_min_plist_values_select_poly_by_numv
      (mesh_data, vertex_coord, flag_internal,
       0, num_poly_vert, min_ratio,
       poly_with_min_edge_length_ratio,
       compute_polygon_edge_length_ratio);
  }
}


// Compute min edge length ratios of polygons.
// - Version which does not return polygons with min edge length ratio.
void IJKMESHINFO::compute_min_polygon_edge_length_ratios_select_poly_by_numv
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal, const int num_poly_vert,
 COORD_TYPE & min_ratio)
{
  int poly_with_min_edge_length_ratio;

  compute_min_polygon_edge_length_ratios_select_poly_by_numv
    (mesh_data, vertex_coord, flag_internal, num_poly_vert,
     min_ratio, poly_with_min_edge_length_ratio);
}


// Compute number of polygons with edge length ratio at most
//   max_edge_length_ratio.
// @pre mesh_dimension == 2.
void IJKMESHINFO::compute_num_polygons_with_small_edge_length_ratios
(const int dimension,
 const POLYMESH_TYPE & polymesh,
 const std::vector<POLY_DATA> & poly_data,
 const COORD_TYPE * vertex_coord,
 const bool flag_internal, 
 const COORD_TYPE ratio_bound, 
 int & num_le)
{
  // Initialize to zero.
  num_le = 0;

  for (int ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

    if (flag_internal) {
      if (poly_data[ipoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    COORD_TYPE edge_length_ratio;
    int num_lengths;
    const int * pvert = polymesh.VertexList(ipoly);
    const int num_poly_vert = polymesh.NumPolyVert(ipoly);
    compute_polygon_edge_length_ratio
      (dimension, pvert, num_poly_vert, vertex_coord,
       edge_length_ratio, num_lengths);

    if (num_lengths > 0 && edge_length_ratio <= ratio_bound) {
      num_le++;
    }
  }

}


// **************************************************
// COMPUTE JACOBIAN DETERMINANTS
// **************************************************

// Compute min/max Jacobian matrix determinants of hexahedra.
void IJKMESHINFO::compute_min_max_hexahedra_Jacobian_determinants
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant)
{
  compute_min_max_plist_values
    (mesh_data, vertex_coord, flag_internal,
     0, 0, min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant,
     compute_min_max_hexahedron_Jacobian_determinants);
}


// Compute min/max Jacobian matrix determinants of hexahedra.
// - Version which does not return hexahedra containing min/max 
//   Jacobian determinants.
void IJKMESHINFO::compute_min_max_hexahedra_Jacobian_determinants
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  int poly_with_min_Jacobian_determinant;
  int poly_with_max_Jacobian_determinant;

  compute_min_max_hexahedra_Jacobian_determinants
    (mesh_data, vertex_coord, flag_internal,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant);
}

// Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hexahedron_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_determinant,
 COORD_TYPE & max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);

  IJK::compute_min_max_hexahedron_Jacobian_determinant_3D
    (hex_vert, mesh_data.orientation, vertex_coord, cube,
     min_Jacobian_determinant, max_Jacobian_determinant);

  num_Jacobian_determinants = 9;
}


// Compute min/max Jacobian matrix determinants of hexahedra vertices.
// - Version which returns vertices with min/max Jacobian determinants.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant,
 int & vert_with_min_Jacobian_determinant, 
 int & vert_with_max_Jacobian_determinant)
{
  compute_min_max_vlist_values
    (mesh_data, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     0, 0, min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant, vert_with_max_Jacobian_determinant,
     compute_min_max_hex_vert_Jacobian_determinants);
}


// Compute min/max of the Jacobian matrix determinants of hexahedra vertices.
// - Version which does not return vertices or polytopes with min/max values.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly, const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  int poly_with_min_Jacobian_determinant;
  int poly_with_max_Jacobian_determinant;
  int vert_with_min_Jacobian_determinant;
  int vert_with_max_Jacobian_determinant;

  compute_min_max_hex_vert_Jacobian_determinants
    (mesh_data, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant,
     poly_with_max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant,
     vert_with_max_Jacobian_determinant);
}


// Compute min/max Jacobian matrix determinants of a vertex 
//   in a hexahedral mesh.
// - Version with input vertex-poly incidence data structure.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const VERTEX_INDEX iv0,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);
  const int NUM_CUBE_VERT = cube.NumVertices();
  const int num_incident_poly = 
    vertex_poly_incidence.NumIncidentPoly(iv0);
  COORD_TYPE Jacobian_determinant;

  // Initialize
  num_Jacobian_determinants = 0;
  min_Jacobian_determinant = 0.0;
  max_Jacobian_determinant = 0.0;
  poly_with_min_Jacobian_determinant = 0;
  poly_with_max_Jacobian_determinant = 0;

  if (flag_internal_vert) {
    if (mesh_data.vertex_data[iv0].OnBoundary()) { 
      // Vertex iv0 is not internal.
      return;
    }
  }

  for (int k = 0; k < num_incident_poly; k++) {

    const int kpoly = vertex_poly_incidence.IncidentPoly(iv0, k);
    const int kloc = vertex_poly_incidence.VertexLocInPolyVertexList(iv0, k);
    IJK::compute_Jacobian_determinant_at_hex_vertex_3D
      (mesh_data.VertexList(kpoly), mesh_data.orientation,
       vertex_coord, cube, kloc, Jacobian_determinant);

    if (flag_internal_poly) {
      if (mesh_data.poly_data[kpoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    if (num_Jacobian_determinants == 0) {
      min_Jacobian_determinant = Jacobian_determinant;
      max_Jacobian_determinant = Jacobian_determinant;
      poly_with_min_Jacobian_determinant = kpoly;
      poly_with_max_Jacobian_determinant = kpoly;
    }
    else {

      if (Jacobian_determinant < min_Jacobian_determinant) {
        min_Jacobian_determinant = Jacobian_determinant;
        poly_with_min_Jacobian_determinant = kpoly;
      }

      if (Jacobian_determinant > max_Jacobian_determinant) {
        max_Jacobian_determinant = Jacobian_determinant;
        poly_with_max_Jacobian_determinant = kpoly;
      }
    }

    num_Jacobian_determinants++;
  }

}


// **************************************************
// COMPUTE NORMALIZED JACOBIAN DETERMINANTS
// **************************************************

// Compute min/max normalized Jacobian matrix determinants of hexahedra.
void IJKMESHINFO::compute_min_max_hexahedra_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant)
{
  compute_min_max_plist_values
    (mesh_data, vertex_coord, flag_internal,
     0, 0, min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant,
     compute_min_max_hexahedron_normalized_Jacobian_determinants);
}


// Compute min/max normalized Jacobian matrix determinants of hexahedra.
// - Version which does not return hexahedra containing min/max 
//   normalized Jacobian determinants.
void IJKMESHINFO::compute_min_max_hexahedra_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  int poly_with_min_value;
  int poly_with_max_value;

  compute_min_max_hexahedra_normalized_Jacobian_determinants
    (mesh_data, vertex_coord, flag_internal,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_value, poly_with_max_value);
}


// Compute min/max of the nine normalized Jacobian matrix determinants 
//   of a hexahedron.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hexahedron_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_determinant,
 COORD_TYPE & max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);
  const COORD_TYPE max_small_magnitude(0.0);

  IJK::compute_min_max_hexahedron_normalized_Jacobian_determinant_3D
    (hex_vert, mesh_data.orientation, vertex_coord, cube, max_small_magnitude,
     min_Jacobian_determinant, max_Jacobian_determinant, 
     num_Jacobian_determinants);
}


// Compute min/max of the eight normalized Jacobian matrix determinants 
//   at the eight vertices of a hexahedron.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_determinant,
 COORD_TYPE & max_Jacobian_determinant,
 VERTEX_INDEX & vert_with_min_Jacobian_determinant,
 VERTEX_INDEX & vert_with_max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);

  const COORD_TYPE max_small_magnitude(0.0);

  IJK::compute_min_max_hex_vert_normalized_Jacobian_determinant_3D
    (hex_vert, mesh_data.orientation, vertex_coord, cube, max_small_magnitude,
     min_Jacobian_determinant, max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant, 
     vert_with_max_Jacobian_determinant,
     num_Jacobian_determinants);
}


// Compute min/max normalized Jacobian matrix determinants 
//   of hexahedra vertices.
// - Version which returns vertices with min/max Jacobian determinants.
void IJKMESHINFO::compute_min_max_hex_vert_normalized_Jacobian_determinants
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
 int & vert_with_max_Jacobian_determinant)
{
  compute_min_max_vlist_values
    (mesh_data, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     0, 0, min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant, poly_with_max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant, vert_with_max_Jacobian_determinant,
     compute_min_max_hex_vert_normalized_Jacobian_determinants);
}


// Compute min/max of the normalized Jacobian matrix determinants 
//   of hexahedra vertices.
// - Version which does not return vertices or polytopes with min/max values.
void IJKMESHINFO::compute_min_max_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly, const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant)
{
  int poly_with_min_Jacobian_determinant;
  int poly_with_max_Jacobian_determinant;
  int vert_with_min_Jacobian_determinant;
  int vert_with_max_Jacobian_determinant;

  compute_min_max_hex_vert_normalized_Jacobian_determinants
    (mesh_data, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_determinant, max_Jacobian_determinant,
     poly_with_min_Jacobian_determinant,
     poly_with_max_Jacobian_determinant,
     vert_with_min_Jacobian_determinant,
     vert_with_max_Jacobian_determinant);
}


// Compute min/max normalized Jacobian matrix determinants of a vertex 
//   in a hexahedral mesh.
// - Version with input vertex-poly incidence data structure.
void IJKMESHINFO::compute_min_max_hex_vert_normalized_Jacobian_determinants
(const MESH_DATA & mesh_data,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const VERTEX_INDEX iv0,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_determinant, COORD_TYPE & max_Jacobian_determinant,
 int & poly_with_min_Jacobian_determinant, 
 int & poly_with_max_Jacobian_determinant,
 int & num_Jacobian_determinants)
{
  const static CUBE_TYPE cube(DIM3);
  const int NUM_CUBE_VERT = cube.NumVertices();
  const int num_incident_poly = 
    vertex_poly_incidence.NumIncidentPoly(iv0);
  const COORD_TYPE max_small_magnitude(0.0);
  bool flag_zero;
  COORD_TYPE Jacobian_determinant;

  // Initialize
  num_Jacobian_determinants = 0;
  min_Jacobian_determinant = 0.0;
  max_Jacobian_determinant = 0.0;
  poly_with_min_Jacobian_determinant = 0;
  poly_with_max_Jacobian_determinant = 0;

  if (flag_internal_vert) {
    if (mesh_data.vertex_data[iv0].OnBoundary()) { 
      // Vertex iv0 is not internal.
      return;
    }
  }

  for (int k = 0; k < num_incident_poly; k++) {

    const int kpoly = vertex_poly_incidence.IncidentPoly(iv0, k);
    const int kloc = vertex_poly_incidence.VertexLocInPolyVertexList(iv0, k);
    IJK::compute_normalized_Jacobian_determinant_at_hex_vertex_3D
      (mesh_data.VertexList(kpoly), mesh_data.orientation,
       vertex_coord, cube, kloc, max_small_magnitude,
       Jacobian_determinant, flag_zero);

    if (flag_internal_poly) {
      if (mesh_data.poly_data[kpoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    if (flag_zero) { continue; }

    if (num_Jacobian_determinants == 0) {
      min_Jacobian_determinant = Jacobian_determinant;
      max_Jacobian_determinant = Jacobian_determinant;
      poly_with_min_Jacobian_determinant = kpoly;
      poly_with_max_Jacobian_determinant = kpoly;
    }
    else {

      if (Jacobian_determinant < min_Jacobian_determinant) {
        min_Jacobian_determinant = Jacobian_determinant;
        poly_with_min_Jacobian_determinant = kpoly;
      }

      if (Jacobian_determinant > max_Jacobian_determinant) {
        max_Jacobian_determinant = Jacobian_determinant;
        poly_with_max_Jacobian_determinant = kpoly;
      }
    }

    num_Jacobian_determinants++;
  }

}


// **************************************************
// COMPUTE SHAPE METRIC BASED ON JACOBIAN MATRICES
// **************************************************

// Compute min/max hexahedra shape metric based on Jacobian matrices.
void IJKMESHINFO::compute_min_max_hexahedra_Jacobian_shape
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape,
 int & poly_with_min_Jacobian_shape, 
 int & poly_with_max_Jacobian_shape)
{
 compute_min_max_plist_values
    (mesh_data, vertex_coord, flag_internal,
     0, 0, min_Jacobian_shape, max_Jacobian_shape,
     poly_with_min_Jacobian_shape, poly_with_max_Jacobian_shape,
     compute_min_max_hexahedron_Jacobian_shape);
}


// Compute min/max hexahedra shape metric based on Jacobian matrices.
// - Version which does not return hexahedra containing min/max shape.
void IJKMESHINFO::compute_min_max_hexahedra_Jacobian_shape
(const MESH_DATA & mesh_data, const COORD_TYPE * vertex_coord,
 const bool flag_internal,
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape)
{
  int poly_with_min_value;
  int poly_with_max_value;

  compute_min_max_hexahedra_Jacobian_shape
    (mesh_data, vertex_coord, flag_internal,
     min_Jacobian_shape, max_Jacobian_shape,
     poly_with_min_value, poly_with_max_value);
}


// Compute min/max hexahedron shape metric based on Jacobian matrices.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hexahedron_Jacobian_shape
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_shape,
 COORD_TYPE & max_Jacobian_shape,
 int & num_Jacobian_shapes)
{
  const static CUBE_TYPE cube(DIM3);
  const COORD_TYPE max_small_magnitude(0.0);

  IJK::compute_min_max_hexahedron_Jacobian_shape_3D
    (hex_vert, vertex_coord, cube, max_small_magnitude,
     min_Jacobian_shape, max_Jacobian_shape,
     num_Jacobian_shapes);
}


// Compute min/max of the eight Jacobian shape metrics
//   at the eight vertices of a hexahedron.
// @pre dimension = 3. 
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data,
 const VERTEX_INDEX hex_vert[], const int num_vert,
 const COORD_TYPE * vertex_coord,
 COORD_TYPE & min_Jacobian_shape,
 COORD_TYPE & max_Jacobian_shape,
 VERTEX_INDEX & vert_with_min_Jacobian_shape,
 VERTEX_INDEX & vert_with_max_Jacobian_shape,
 int & num_Jacobian_shapes)
{
  const static CUBE_TYPE cube(DIM3);

  const COORD_TYPE max_small_magnitude(0.0);

  IJK::compute_min_max_hex_vert_Jacobian_shape_3D
    (hex_vert, vertex_coord, cube, max_small_magnitude,
     min_Jacobian_shape, max_Jacobian_shape,
     vert_with_min_Jacobian_shape, 
     vert_with_max_Jacobian_shape,
     num_Jacobian_shapes);
}


// Compute min/max shape metrics based on Jacobian matrices at hex vertices.
// - Version which returns vertices with min/max shape metric.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_shape
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
 int & vert_with_max_Jacobian_shape)
{
  compute_min_max_vlist_values
    (mesh_data, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     0, 0, min_Jacobian_shape, max_Jacobian_shape,
     poly_with_min_Jacobian_shape, poly_with_max_Jacobian_shape,
     vert_with_min_Jacobian_shape, vert_with_max_Jacobian_shape,
     compute_min_max_hex_vert_Jacobian_shape);
}


// Compute min/max shape metrics based on Jacobian matrices at hex vertices.
// - Version which does not return vertices or polytopes with min/max values.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const bool flag_internal_poly, const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape)
{
  int poly_with_min_Jacobian_shape;
  int poly_with_max_Jacobian_shape;
  int vert_with_min_Jacobian_shape;
  int vert_with_max_Jacobian_shape;

  compute_min_max_hex_vert_Jacobian_shape
    (mesh_data, vertex_poly_incidence, vertex_coord, 
     flag_internal_poly, flag_internal_vert,
     min_Jacobian_shape, max_Jacobian_shape,
     poly_with_min_Jacobian_shape,
     poly_with_max_Jacobian_shape,
     vert_with_min_Jacobian_shape,
     vert_with_max_Jacobian_shape);
}


// Compute min/max shape metrics based on Jacobian matrices at hex vertices.
// - Version with input vertex-poly incidence data structure.
void IJKMESHINFO::compute_min_max_hex_vert_Jacobian_shape
(const MESH_DATA & mesh_data,
 const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence, 
 const COORD_TYPE * vertex_coord,
 const VERTEX_INDEX iv0,
 const bool flag_internal_poly,
 const bool flag_internal_vert,
 COORD_TYPE & min_Jacobian_shape, COORD_TYPE & max_Jacobian_shape,
 int & poly_with_min_Jacobian_shape, 
 int & poly_with_max_Jacobian_shape,
 int & num_Jacobian_shapes)
{
  const static CUBE_TYPE cube(DIM3);
  const int NUM_CUBE_VERT = cube.NumVertices();
  const int num_incident_poly = 
    vertex_poly_incidence.NumIncidentPoly(iv0);
  const COORD_TYPE max_small_magnitude(0.0);
  bool flag_zero;
  COORD_TYPE Jacobian_shape;

  // Initialize
  num_Jacobian_shapes = 0;
  min_Jacobian_shape = 0.0;
  max_Jacobian_shape = 0.0;
  poly_with_min_Jacobian_shape = 0;
  poly_with_max_Jacobian_shape = 0;

  if (flag_internal_vert) {
    if (mesh_data.vertex_data[iv0].OnBoundary()) { 
      // Vertex iv0 is not internal.
      return;
    }
  }

  for (int k = 0; k < num_incident_poly; k++) {

    const int kpoly = vertex_poly_incidence.IncidentPoly(iv0, k);
    const int kloc = vertex_poly_incidence.VertexLocInPolyVertexList(iv0, k);
    IJK::compute_Jacobian_shape_at_hex_vertex_3D
      (mesh_data.VertexList(kpoly), vertex_coord, cube, kloc, 
       max_small_magnitude, Jacobian_shape, flag_zero);


    if (flag_internal_poly) {
      if (mesh_data.poly_data[kpoly].ContainsBoundaryFacet()) 
        { continue; } 
    }

    if (flag_zero) { continue; }

    if (num_Jacobian_shapes == 0) {
      min_Jacobian_shape = Jacobian_shape;
      max_Jacobian_shape = Jacobian_shape;
      poly_with_min_Jacobian_shape = kpoly;
      poly_with_max_Jacobian_shape = kpoly;
    }
    else {

      if (Jacobian_shape < min_Jacobian_shape) {
        min_Jacobian_shape = Jacobian_shape;
        poly_with_min_Jacobian_shape = kpoly;
      }

      if (Jacobian_shape > max_Jacobian_shape) {
        max_Jacobian_shape = Jacobian_shape;
        poly_with_max_Jacobian_shape = kpoly;
      }
    }

    num_Jacobian_shapes++;
  }

}


// **************************************************
// INTERSECT POLYTOPES
// **************************************************

// Compute min/max coordinates of triangle vertices.
void IJKMESHINFO::compute_min_max_coord
(const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
 const COORD_TYPE p2[DIM3],
 COORD_TYPE min_coord[DIM3], COORD_TYPE max_coord[DIM3])
{
  COORD_TYPE minC, maxC;

  for (int d = 0; d < DIM3; d++) {
    minC = p0[d];
    maxC = p0[d];
    if (p1[d] < minC) { minC = p1[d]; }
    if (p1[d] > maxC) { maxC = p1[d]; }
    if (p2[d] < minC) { minC = p2[d]; }
    if (p2[d] > maxC) { maxC = p2[d]; }
    min_coord[d] = minC;
    max_coord[d] = maxC;
  }
}


// Compute orientation of (p1-p0), (p2-p0) and w.
int IJKMESHINFO::compute_orientation
(const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
 const COORD_TYPE p2[DIM3], const COORD_TYPE w[DIM3],
 const COORD_TYPE epsilon)
{
  COORD_TYPE v01[DIM3], v02[DIM3];
  COORD_TYPE det;

  using namespace IJK;
  
  subtract_coord_3D(p1, p0, v01);
  subtract_coord_3D(p2, p0, v02);

  determinant_3x3(v01, v02, w, det);
  if (std::fabs(det) < epsilon) { det = 0; }

  if (det < 0) { return(-1); }
  else if (det > 0) { return(1); }
  else { return(0); }
}


// Compute affine subspace spanning points p0, p1, p2.
// @param span_dimension Dimension of spanning subspace. 0, 1 or 2.
// @param w If span_dimension = 2, w is unit normal to the plane.
// @param w If span_dimension = 1, w is unit line direction.
// @param w If span_dimension = 0, w is p0.
void IJKMESHINFO::compute_span
(const COORD_TYPE p0[DIM3], const COORD_TYPE p1[DIM3],
 const COORD_TYPE p2[DIM3], 
 int & span_dimension, COORD_TYPE w[DIM3])
{
  COORD_TYPE v01[DIM3], v02[DIM3], v12[DIM3];
  COORD_TYPE normalized_v01[DIM3], normalized_v02[DIM3], 
    normalized_v12[DIM3];
  COORD_TYPE v01_mag, v02_mag, v12_mag;
  COORD_TYPE_PTR u0, u1, u2;
  COORD_TYPE u1_orth[DIM3], u2_orth[DIM3];
  COORD_TYPE normalized_u1_orth[DIM3], normalized_u2_orth[DIM3];
  COORD_TYPE u1_orth_mag, u2_orth_mag, w_mag;

  using namespace IJK;
  
  subtract_coord_3D(p1, p0, v01);
  subtract_coord_3D(p2, p0, v02);
  subtract_coord_3D(p2, p1, v12);

  normalize_vector_robust(DIM3, v01, 0, normalized_v01, v01_mag);
  normalize_vector_robust(DIM3, v02, 0, normalized_v02, v02_mag);
  normalize_vector_robust(DIM3, v12, 0, normalized_v12, v12_mag);

  if (v01_mag <= 0 && v02_mag <= 0 && v12_mag <= 0) {
    span_dimension = 0;
    copy_coord_3D(p0, w);
    return;
  }

  if (v01_mag >= v02_mag && v01_mag >= v12_mag) {
    u0 = normalized_v01;
    u1 = v02;
    u2 = v12;
  }
  else if (v02_mag >= v12_mag) {
    u0 = normalized_v02;
    u1 = v01;
    u2 = v12;
  }
  else {
    u0 = normalized_v12;
    u1 = v01;
    u2 = v02;
  }

  compute_orthogonal_vector(DIM3, u1, u0, u1_orth);
  compute_orthogonal_vector(DIM3, u2, u0, u2_orth);

  normalize_vector_robust(DIM3, u1_orth, 0, normalized_u1_orth, u1_orth_mag);
  normalize_vector_robust(DIM3, u2_orth, 0, normalized_u2_orth, u2_orth_mag);

  if (u1_orth_mag == 0 && u2_orth_mag == 0) {
    span_dimension = 1;
    copy_coord_3D(u0, w);
    return;
  }

  span_dimension = 2;
  if (u1_orth_mag >= u2_orth_mag) {
    compute_cross_product_3D(u0, normalized_u1_orth, w);
  }
  else {
    compute_cross_product_3D(u0, normalized_u2_orth, w);
  }

  // Renormalized, just in case.
  normalize_vector_robust(DIM3, w, 0, w, w_mag);
}


// Return true if triangle contains point.
// Assumes plane containes point and triangle.
// @param normal[] Normal to plane containing point and triangle.
bool IJKMESHINFO::does_triangle_contain_point
(const COORD_TYPE w0[DIM3], const COORD_TYPE w1[DIM3],
 const COORD_TYPE w2[DIM3], const COORD_TYPE p[DIM3],
 const COORD_TYPE normal[DIM3], const COORD_TYPE epsilon)
{
  COORD_TYPE min_coord[DIM3], max_coord[DIM3];
  int sign0, sign1, sign2;

  // Handles cases where w0, w1, w2 and p are collinear.
  compute_min_max_coord(w0, w1, w2, min_coord, max_coord);
  for (int d = 0; d < DIM3; d++) {
    if (p[d] + epsilon < min_coord[d]) { return(false); }
    if (p[d] - epsilon > max_coord[d]) { return(false); }
  }

  // Check that intersection_point is in triangle (w0,w1,w2).
  sign0 = compute_orientation(w0, w1, p, normal, epsilon);
  sign1 = compute_orientation(w1, w2, p, normal, epsilon);
  sign2 = compute_orientation(w2, w0, p, normal, epsilon);

  if (!are_equal_or_zero(sign0, sign1)) { return(false); }
  if (!are_equal_or_zero(sign1, sign2)) { return(false); }
  if (!are_equal_or_zero(sign0, sign2)) { return(false); }

  return(true);
}


// Compute intersection of edge (v0,v1) and triangle (w0,w1,w2).
// - Use vectors (w1-w0) and (w2-w0)
// - Return true if edge intersects triangle.
bool IJKMESHINFO::intersect_edge_triangle
(const COORD_TYPE v0[DIM3], const COORD_TYPE v1[DIM3],
 const COORD_TYPE w0[DIM3], const COORD_TYPE w1[DIM3],
 const COORD_TYPE w2[DIM3], const COORD_TYPE epsilon,
 COORD_TYPE intersection_point[DIM3])
{
  COORD_TYPE w[DIM3];
  COORD_TYPE v01[DIM3], normalized_v01[DIM3];
  COORD_TYPE v0w0[DIM3];
  COORD_TYPE u0[DIM3], u1[DIM3];
  COORD_TYPE w_mag, v01_mag;
  COORD_TYPE s0, s1, t;
  COORD_TYPE product0, product1;
  int span_dimension;

  using namespace IJK;
  
  // No intersection if line segment and triangle share a vertex.
  if (is_coord_equal_3D(v0, w0)) { return(false); };
  if (is_coord_equal_3D(v0, w1)) { return(false); };
  if (is_coord_equal_3D(v0, w2)) { return(false); };
  if (is_coord_equal_3D(v1, w0)) { return(false); };
  if (is_coord_equal_3D(v1, w1)) { return(false); };
  if (is_coord_equal_3D(v1, w2)) { return(false); };

  compute_span(w0, w1, w2, span_dimension, w);

  if (span_dimension < 2) { return(false); };

  // Check w is non-zero.
  compute_magnitude_3D(w, w_mag);
  if (w_mag <= epsilon) { return(false); };

  subtract_coord_3D(v1, v0, v01);
  subtract_coord_3D(w0, v0, v0w0);

  compute_inner_product_3D(v01, w, s0);
  compute_inner_product_3D(v0w0, w, s1);

  if (s0 == 0) { return(false); }

  if (fabs(s1) <= (1+epsilon)*fabs(s0)) 
    { t = s1/s0; }
  else
    { return(false); }

  if (t < -epsilon || t > 1+epsilon) 
    { return(false); }

  if (t*w_mag < -epsilon || t*w_mag > 1+epsilon)
    { return(false); }

  add_scaled_coord_3D(t, v01, v0, intersection_point);

  // Check that intersection_point is on line segment (v0,v1).
  subtract_coord_3D(intersection_point, v0, u0);
  subtract_coord_3D(intersection_point, v1, u1);
  compute_inner_product_3D(u0, v01, product0);
  compute_inner_product_3D(u1, v01, product1);

  if (product0 < epsilon || product1 > -epsilon) { return(false); }

  bool flag = does_triangle_contain_point
    (w0, w1, w2, intersection_point, w, epsilon);

  return(flag);
}


/*!
 *  @brief Return true if edges of triangle jt1 intersect jt2.
 */
bool IJKMESHINFO::intersect_triangle_edges_and_triangle
(const VERTEX_INDEX triangle_vert[],
 const COORD_TYPE vertex_coord[],
 const VERTEX_INDEX vertex_merge[],
 const int jt1, const int jt2, const double epsilon,
 COORD_TYPE intersection_point[DIM3])
{
  const int NUM_VERT_PER_TRIANGLE(3);

  const COORD_TYPE * w0 = 
    vertex_coord + DIM3*triangle_vert[jt2*DIM3];
  const COORD_TYPE * w1 = 
    vertex_coord + DIM3*triangle_vert[jt2*DIM3+1];
  const COORD_TYPE * w2 = 
    vertex_coord + DIM3*triangle_vert[jt2*DIM3+2];

  for (int k0 = 0; k0 < NUM_VERT_PER_TRIANGLE; k0++) {
    const COORD_TYPE * v0 = 
      vertex_coord + DIM3*triangle_vert[jt1*DIM3+k0];
    const int k1 = (k0+1)%DIM3;
    const COORD_TYPE * v1 = 
      vertex_coord + DIM3*triangle_vert[jt1*DIM3+k1];
        
    if (intersect_edge_triangle
        (v0, v1, w0, w1, w2, epsilon, intersection_point))
      { return(true); }
  }

  return false;
}


/*!
 *  @brief Return true if triangles intersect.
 *  - Computes intersection between every edge of one triangle
 *    with the other triangle.
 *  - Ignores edges of one triangle that have an endpoint
 *    at the same coordinates as some vertex of the other triangle.
 *  - Some algorithms do less comparisons, but comparing every edge
 *    of one triangle with the other triangle is possibly more robust.
 */
bool IJKMESHINFO::intersect_triangle_triangle
(const VERTEX_INDEX triangle_vert[],
 const COORD_TYPE vertex_coord[],
 const VERTEX_INDEX vertex_merge[],
 const int jt0, const int jt1, const double epsilon,
 COORD_TYPE intersection_point[DIM3])
{
  const int NUM_VERT_PER_TRIANGLE(3);

  if (intersect_triangle_edges_and_triangle
      (triangle_vert, vertex_coord, vertex_merge,
       jt1, jt0, epsilon, intersection_point))
    { return true; }

  if (intersect_triangle_edges_and_triangle
      (triangle_vert, vertex_coord, vertex_merge,
       jt0, jt1, epsilon, intersection_point))
    { return true; }

  for (int k0 = 0; k0 < NUM_VERT_PER_TRIANGLE; k0++) {
    for (int k1 = 0; k1 < NUM_VERT_PER_TRIANGLE; k1++) {
      const VERTEX_INDEX iv0 = 
        triangle_vert[jt0*NUM_VERT_PER_TRIANGLE+k0];
      const VERTEX_INDEX iv1 = 
        triangle_vert[jt1*NUM_VERT_PER_TRIANGLE+k1];

      if (iv0 == iv1) { continue; }

      const COORD_TYPE * v0 = vertex_coord + DIM3*iv0;
      const COORD_TYPE * v1 = vertex_coord + DIM3*iv1;

      if (IJK::is_coord_equal(DIM3, v0, v1)) {
        if (vertex_merge[iv0] != vertex_merge[iv1]) {
          // Distinct vertices iv0 and iv1 have same coordinates
          //   but are not merged.
          // Mesh intersects at v0.
          IJK::copy_coord(DIM3, v0, intersection_point);
          return true;
        }
      }
    }
  }

  return false;
}


bool IJKMESHINFO::intersect_triangle_triangle
(const VERTEX_INDEX triangle_vert[],
 const COORD_TYPE vertex_coord[],
 const std::vector<VERTEX_INDEX> & vertex_merge,
 const int jt1, const int jt2, const double epsilon,
 COORD_TYPE intersection_point[DIM3])
{
  return intersect_triangle_triangle
    (triangle_vert, vertex_coord, IJK::vector2pointer(vertex_merge),
     jt1, jt2, epsilon, intersection_point);
}


// **************************************************
// COMPUTE MESH SELF INTERSECTIONS
// **************************************************


// Compute self intersections of triangle mesh. Simple (O(n^2)) algorithm.
// - Triangle mesh is embedded in R3.
void IJKMESHINFO::compute_trimesh_self_intersections_simple
(const VERTEX_INDEX triangle_vert[],  const int num_triangles,
 const COORD_TYPE vertex_coord[],
 const COORD_TYPE selfI_epsilon,
 const std::vector<VERTEX_INDEX> & vertex_merge,
 INTERSECTING_POLY_ARRAY & intersecting_poly,
 std::vector<COORD_TYPE> & intersection_coord)
{
  COORD_TYPE intersection_point[DIM3];

  // Initialize
  intersecting_poly.clear();
  intersection_coord.clear();
  
  for (int jt1 = 0; jt1 < num_triangles; jt1++) {
    for (int jt2 = jt1+1; jt2 < num_triangles; jt2++) {

      if (intersect_triangle_triangle
          (triangle_vert, vertex_coord, vertex_merge,
           jt1, jt2, selfI_epsilon, intersection_point)) {
        intersecting_poly.push_back(std::make_pair(jt1,jt2));
        IJK::push_backIII(intersection_point, intersection_coord);
      }
    }
  }

}



