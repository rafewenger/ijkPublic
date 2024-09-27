/// \file ijkmcube_sub.cxx
/// Subroutines for ijkmcube

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2009 Rephael Wenger

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

#include <assert.h>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "ijkmcube_sub.h"
#include "ijkmcube_util.h"
#include "ijkoctree.h"
#include "ijkMCtable.h"

#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijkpoly.tpp"
#include "ijkmerge.tpp"
#include "ijkisopoly.tpp"
#include "ijkisocoord.tpp"

using namespace IJK;
using namespace IJKMCUBE;
using namespace IJKMCUBE_TABLE;

// **************************************************
// MARCHING CUBES SUBROUTINES
// **************************************************

// Call merge_identical_vertices and then position_vertices_linear.
// - Isosurface vertices lie on regular grid edges.
void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const std::vector<ISO_VERTEX_INDEX> & iso_simplices,
 const ISOTABLE_TYPE isotable_type,
 std::vector<ISO_VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 MERGE_DATA & merge_data, MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  std::vector<ISO_VERTEX_INDEX> iso_vlist;
  merge_identical(iso_simplices, iso_vlist, simplex_vert, merge_data);
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = iso_vlist.size();

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  if (numc > 0) {
    position_iso_vertices_linear
      (scalar_grid, isovalue, isotable_type, iso_vlist, &(vertex_coord[0]));
  }

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}


// Call merge_identical_edges and then position_iso_vertices_linearB.
//  - Version with intersected edges specified by pairs of endpoints
//    in elist_endpoint[].
void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const std::vector<VERTEX_INDEX> & elist_endpoint,
 const ISOTABLE_TYPE isotable_type,
 std::vector<ISO_VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_edges_parameter,
 MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  std::vector<VERTEX_INDEX> merged_endpoint;

  merge_identical_edges
    (elist_endpoint, merged_endpoint, simplex_vert, merge_edges_parameter);

  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = merged_endpoint.size()/2;

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  position_iso_vertices_linearB
    (scalar_grid, isovalue, isotable_type, merged_endpoint, 
     &(vertex_coord[0]));

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}


// Call merge_identical_edges and then position_iso_vertices_linearB
// - Version with intersected edges specified by pairs of endpoints
//   in elist_endpoint[].
// - Version where new, non-grid mesh vertices are stored
//   in new_mesh_vertices.
void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue, 
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & elist_endpoint,
 const ISOTABLE_TYPE isotable_type,
 std::vector<ISO_VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info)
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  std::vector<VERTEX_INDEX> merged_endpoint;
  merge_identical_edges
    (elist_endpoint, merged_endpoint, simplex_vert, merge_parameters);
			
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = merged_endpoint.size()/2;

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  position_iso_vertices_linearB
    (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
     merged_endpoint, &(vertex_coord[0]));

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}


// Call merge_identical_edges and then position_iso_vertices_linearB
// Input includes an array of edge endpoints.  Isosurface vertices lie on these edges.
// scalar_grid = scalar grid data
// isovalue = isosurface scalar value
// simplex_vert[] = vector of isosurface simplex vertices
//   simplex_vert[dimension*is+k] =
//     grid edge containing k'th vertex of simplex is.
// isotable_type = type of isosurface table
// endpoint[] = list of isosurface simplex vertices
//   (endpoint[2*i], endpoint[2*i+1]) = endoints
//    of edge containing isosurface vertex i.
// vertex_coord[] = list of isosurface vertex coordinates
//   vertex_coord[dimension*iv+k] = k'th coordinate of vertex iv
// merge_data = internal data structure for merging identical edges
// mcube_info = information about running time and grid cubes and edges
void IJKMCUBE::merge_and_position_vertices
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & endpoint,
 const ISOTABLE_TYPE isotable_type,
 const INTERPOLATION_TYPE interpolation_type,
 std::vector<ISO_VERTEX_INDEX> & simplex_vert,
 std::vector<COORD_TYPE> & vertex_coord,
 const MERGE_EDGES_PARAMETERS & merge_parameters, MCUBE_INFO & mcube_info)

{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  PROCEDURE_ERROR error("merge_and_position_vertices");

  clock_t t1 = clock();

  std::vector<VERTEX_INDEX> merged_endpoint;
  merge_identical_edges
    (endpoint, merged_endpoint, simplex_vert, merge_parameters);
			
  clock_t t2 = clock();

  // num_isov = number of isosurface vertices
  int num_isov = merged_endpoint.size()/2;

  int numc = num_isov*dimension;
  vertex_coord.resize(numc);

  switch(interpolation_type) {

  case LINEAR_INTERPOLATION:
    position_iso_vertices_linearB
      (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
       merged_endpoint, &(vertex_coord[0]));
    break;

  case MULTILINEAR_INTERPOLATION:
    position_iso_vertices_multilinear
      (scalar_grid, isovalue, isotable_type, new_mesh_vertices,
       merged_endpoint, &(vertex_coord[0]));
    break;

  default:
    error.AddMessage("Illegal interpolation type.");
    throw error;
  }

  clock_t t3 = clock();

  // store times
  mcube_info.time.merge = float(t2-t1)/CLOCKS_PER_SEC;
  mcube_info.time.position = float(t3-t2)/CLOCKS_PER_SEC;
}


// Merge identical vertices in list of isosurface vertices.
// - Return elist1, modifies eindex.
void IJKMCUBE::merge_identical_edges
(const std::vector<VERTEX_INDEX> & elist0,
 std::vector<VERTEX_INDEX> & elist1, 
 std::vector<ISO_VERTEX_INDEX> & eindex,
 const MERGE_EDGES_PARAMETERS & merge_edges_parameter)
{
  const int nume = elist0.size()/2;
  PROCEDURE_ERROR error("merge_identical_edges");

  if (elist0.size() == 0 && eindex.size() > 0) {
    error.AddMessage("Programming error.  elist0 has zero edges but eindex is not empty.");
    throw error;
  }

  if (!check_order(elist0, error)) { throw error; };

  // initialize
  elist1.clear();

  ARRAY<ISO_VERTEX_INDEX> elist0_map(nume);

  merge_pairs(elist0, elist1, elist0_map.Ptr(), 
              merge_edges_parameter);

  if (!check_pair_list
      (elist0, elist1, elist0_map.PtrConst(), error)) 
    { throw error; }

  for (unsigned int i = 0; i < eindex.size(); i++) 
    { eindex[i] = elist0_map[eindex[i]]; }
}


void IJKMCUBE::set_in_facet_iso_patches(NEP_ISOSURFACE_TABLE & isotable)
{
  const int numf = isotable.Polytope().NumFacets();
  const int numv_per_simplex = isotable.NumVerticesPerSimplex();

  isotable.SetIsInFacet(false);

  for (TABLE_INDEX it = 0; it < isotable.NumTableEntries(); it++) {
    int nums = isotable.NumSimplices(it);
    if (nums == 0)   // no isosurface patch for table entry it
      { continue; }

    for (int jf = 0; jf < numf; jf++) {
      bool in_facet = true;
      for (int is = 0; is < nums && in_facet; is++) {
	for (int k = 0; k < numv_per_simplex && in_facet; k++) {
	  ISOSURFACE_VERTEX_INDEX isov = isotable.SimplexVertex(it, is, k);
	  if (isotable.IsosurfaceVertex(isov).Type() ==
	      ISOSURFACE_VERTEX::VERTEX) {
	    int iv = isotable.IsosurfaceVertex(isov).Face();
	    if (!isotable.Polytope().IsVertexInFacet(jf, iv)) 
	      { in_facet = false; };
	  }
	  else {
	    in_facet = false;
	  }
	}
      }

      if (in_facet) {
	isotable.SetContainingFacet(it, jf);
	break;
      }
    }
  }

}


// **************************************************
// POSITION ISO VERTICES
// **************************************************

namespace {

  /// Get scalar values of polyhedron vertices.
  /// *** MAKE THIS A TEMPLATE ***
  inline void get_poly_scalar
  (const SCALAR_TYPE * scalar, const VERTEX_INDEX iv0,
   const VERTEX_INDEX * increment, const int num_vertices, 
   SCALAR_TYPE * poly_scalar)
  {
    for (int i = 0; i < num_vertices; i++) {
      VERTEX_INDEX iv1 = iv0 + increment[i];
      poly_scalar[i] = scalar[iv1];
    }
  }


  // Interpolate intersection of isosurface and mesh edge
  //   using multilinear interpolation on (hyper) cube vertices.
  void interpolate_using_multilinear
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const VERTEX_INDEX * increment, const int num_cube_vertices,
   const COORD_TYPE * coord0, const COORD_TYPE * coord1,
   COORD_TYPE * coord2)
  {
    const int dimension = scalar_grid.Dimension();
    VERTEX_INDEX iv_base;
    IJK::ARRAY<GRID_COORD_TYPE> base_coord(dimension);
    IJK::ARRAY<COORD_TYPE> coordA(dimension);
    IJK::ARRAY<COORD_TYPE> coordB(dimension);
    IJK::ARRAY<COORD_TYPE> coordC(dimension);
    IJK::ARRAY<SCALAR_TYPE> cube_scalar(num_cube_vertices);

    for (int d = 0; d < dimension; d++) {
      COORD_TYPE x = std::min(coord0[d], coord1[d]);
      base_coord[d] = GRID_COORD_TYPE(floor(x));
    }

    iv_base = scalar_grid.ComputeVertexIndex(base_coord.PtrConst());

    get_poly_scalar(scalar_grid.ScalarPtrConst(), iv_base, increment,
                    num_cube_vertices, cube_scalar.Ptr());


    subtract_coord(dimension, coord0, base_coord.PtrConst(), coordA.Ptr());
    subtract_coord(dimension, coord1, base_coord.PtrConst(), coordB.Ptr());

    const int NUM_ITER = 5;
    multilinear_interpolate_coord
    (dimension, coordA.PtrConst(), coordB.PtrConst(), isovalue, 
     num_cube_vertices, cube_scalar.PtrConst(), NUM_ITER, coordC.Ptr());

    add_coord(dimension, coordC.PtrConst(), base_coord.PtrConst(), coord2);
  }


  // Compute position of isosurface vertices using linear interpolation
  // scalar_grid = scalar grid data
  // vlist[] = list of isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  void position_iso_vertices_linear_binary
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & elist, COORD_TYPE * coord)
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    IJK::ARRAY<AXIS_SIZE_TYPE> axis_increment(dimension);
    const int nume = elist.size();
    IJK::ARRAY<GRID_COORD_TYPE> coord0(dimension);
    IJK::ARRAY<GRID_COORD_TYPE> coord1(dimension);
    IJK::ARRAY<COORD_TYPE> coord2(dimension);

    // compute_axis_increment
    compute_increment(scalar_grid, axis_increment.Ptr());

    const int num_isov_per_gridv = 
      get_num_iso_vertices_per_grid_vertex(dimension);


    for (int i = 0; i < nume; i++) {
      int d = elist[i]%num_isov_per_gridv;
      VERTEX_INDEX v0 = elist[i]/num_isov_per_gridv;
      VERTEX_INDEX v1 = v0+axis_increment[d];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      compute_coord(v0, dimension, axis_size, coord0.Ptr());
      compute_coord(v1, dimension, axis_size, coord1.Ptr());

      linear_interpolate_coord
        (dimension, s0, coord0.PtrConst(), 
         s1, coord1.PtrConst(), isovalue, coord2.Ptr());
      for (int d = 0; d < dimension; d++)
        coord[i*dimension+d] = coord2[d];
    }
  }

  void position_iso_vertices_linear_nep
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
  // compute position of isosurface vertices using linear interpolation
  //   vertices can have positive, negative or equals values
  // dimension = volume dimension
  // scalar_grid[iv] = scalar value at grid vertex iv
  // grid_length[i] = # grid vertices along axis i
  // vlist[] = list of isosurface vertices
  // coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   Precondition: array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    IJK::ARRAY<AXIS_SIZE_TYPE> axis_increment(dimension);
    const int numv = vlist.size();
    IJK::ARRAY<GRID_COORD_TYPE> coord0(dimension);
    IJK::ARRAY<GRID_COORD_TYPE> coord1(dimension);
    IJK::ARRAY<COORD_TYPE> coord2(dimension);
    const int VERTEX_OFFSET = dimension;

    // compute_axis_increment
    compute_increment(scalar_grid, axis_increment.Ptr());

    const int num_isov_per_gridv = 
      get_num_nep_iso_vertices_per_grid_vertex(dimension);

    for (int i = 0; i < numv; i++) {
      int k = vlist[i]%num_isov_per_gridv;
      VERTEX_INDEX v0 = vlist[i]/num_isov_per_gridv;

      if (k == VERTEX_OFFSET) {
        // isosurface vertex lies on grid vertex v0
        compute_coord(v0, dimension, axis_size, coord0.Ptr());
        for (int d = 0; d < dimension; d++)
          coord[i*dimension+d] = coord0[d];
      }
      else {
        // isosurface vertex lies on grid edge
        VERTEX_INDEX v1 = v0+axis_increment[k];

        SCALAR_TYPE s0 = scalar[v0];
        SCALAR_TYPE s1 = scalar[v1];

        compute_coord(v0, dimension, axis_size, coord0.Ptr());
        compute_coord(v1, dimension, axis_size, coord1.Ptr());

        linear_interpolate_coord
          (dimension, s0, coord0.PtrConst(), s1, coord1.PtrConst(), 
           isovalue, coord2.Ptr());
        for (int d = 0; d < dimension; d++)
          coord[i*dimension+d] = coord2[d];
      }
    }

  }


  // Compute position of isosurface vertices using linear interpolation.
  // - Version with intersected edges specified by pairs of endpoints
  //   in elist_endpoint[].
  // - Version where new, non-grid mesh vertices are stored
  //   in new_mesh_vertices.
  // @param new_mesh_vertices Data structure storing new, non-grid mesh vertices.
  //   - If iv < scalar_grid.NumVertices(), then iv is a grid vertex.
  //   - If iv >= scalar_grid.NumVertices(), then iv is a
  //     new, non-grid mesh vertex, whose index in new_mesh_vertices
  //     is (iv-scalar_grid.NumVertices().
  void position_iso_vertices_linear_binaryB
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const MC_MESH_VERTEX_LIST & new_mesh_vertices,
   const std::vector<VERTEX_INDEX> & elist_endpoint, COORD_TYPE * coord)
  {
    const int dimension = scalar_grid.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int nume = elist_endpoint.size()/2;
    IJK::ARRAY<COORD_TYPE> coord0(dimension);
    IJK::ARRAY<COORD_TYPE> coord1(dimension);
    IJK::ARRAY<COORD_TYPE> coord2(dimension);
    IJK::ARRAY<VERTEX_INDEX> axis_increment(dimension);

    const int num_cube_vertices = compute_num_cube_vertices(dimension);

    IJK::ARRAY<VERTEX_INDEX> increment(num_cube_vertices);
    compute_increment(scalar_grid, axis_increment.Ptr());
    compute_cube_vertex_increment
      (dimension, axis_increment.PtrConst(), increment.Ptr());

    for (int i = 0; i < nume; i++) {
      VERTEX_INDEX v0 = elist_endpoint[2*i];
      VERTEX_INDEX v1 = elist_endpoint[2*i+1];

      SCALAR_TYPE s0, s1;

      if (v0 < scalar_grid.NumVertices()) {
        s0 = scalar[v0];
        scalar_grid.ComputeCoord(v0, coord0.Ptr());
      }
      else {
        VERTEX_INDEX k = v0 - scalar_grid.NumVertices();
        new_mesh_vertices.CopyCoord(k, coord0.Ptr());
        s0 = new_mesh_vertices.Scalar(k);
      }

      if (v1 < scalar_grid.NumVertices()) {
        s1 = scalar[v1];
        scalar_grid.ComputeCoord(v1, coord1.Ptr());
      }
      else {
        VERTEX_INDEX k = v1 - scalar_grid.NumVertices();
        new_mesh_vertices.CopyCoord(k, coord1.Ptr());
        s1 = new_mesh_vertices.Scalar(k);
      }

      linear_interpolate_coord
        (dimension, s0, coord0.PtrConst(), s1, coord1.PtrConst(), 
         isovalue, coord2.Ptr());
      for (int d = 0; d < dimension; d++)
        { coord[i*dimension+d] = coord2[d]; }
    }
  }


  // Compute position of isosurface vertices using multilinear interpolation.
  // - Use multilinear interpolation to compute intersection of isosurface
  //   and mesh edge where one endpoint is not a regular grid edge
  //   and one endpoint is a regular grid edge.
  // @param elist_endpoint[] List of edges specified by pairs of edge endpoints.
  //   (elist_endpoint[2*i], elist_endpoint[2*i+1]) is edge i.
  // @param coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
  //   @pre Array coord[] is preallocated to length at least
  //                 dimension*vlist.size()
  void position_iso_vertices_multilinear_binary
  (const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
   const MC_MESH_VERTEX_LIST & new_mesh_vertices,
   const std::vector<VERTEX_INDEX> & elist_endpoint, COORD_TYPE * coord)

  {
    const int dimension = scalar_grid.Dimension();
    const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
    const int nume = elist_endpoint.size()/2;
    IJK::ARRAY<COORD_TYPE> coord0(dimension);
    IJK::ARRAY<COORD_TYPE> coord1(dimension);
    IJK::ARRAY<COORD_TYPE> coord2(dimension);
    IJK::ARRAY<VERTEX_INDEX> axis_increment(dimension);

    const int num_cube_vertices = compute_num_cube_vertices(dimension);

    IJK::ARRAY<VERTEX_INDEX> increment(num_cube_vertices);
    compute_increment(scalar_grid, axis_increment.Ptr());
    compute_cube_vertex_increment
      (dimension, axis_increment.PtrConst(), increment.Ptr());

    for (int i = 0; i < nume; i++) {
      VERTEX_INDEX v0 = elist_endpoint[2*i];
      VERTEX_INDEX v1 = elist_endpoint[2*i+1];

      SCALAR_TYPE s0, s1;
      bool is_new0 = false;
      bool is_new1 = false;

      if (v0 < scalar_grid.NumVertices()) {
        s0 = scalar[v0];
        scalar_grid.ComputeCoord(v0, coord0.Ptr());
      }
      else {
        VERTEX_INDEX k = v0 - scalar_grid.NumVertices();
        new_mesh_vertices.CopyCoord(k, coord0.Ptr());
        s0 = new_mesh_vertices.Scalar(k);
        is_new0 = true;
      }

      if (v1 < scalar_grid.NumVertices()) {
        s1 = scalar[v1];
        scalar_grid.ComputeCoord(v1, coord1.Ptr());
      }
      else {
        VERTEX_INDEX k = v1 - scalar_grid.NumVertices();
        new_mesh_vertices.CopyCoord(k, coord1.Ptr());
        s1 = new_mesh_vertices.Scalar(k);
        is_new1 = true;
      }

      if (is_new0 == is_new1) {
        linear_interpolate_coord
          (dimension, s0, coord0.PtrConst(), s1, coord1.PtrConst(), 
           isovalue, coord2.Ptr());
      }
      else {
        // *** WHY ISN'T MULTILINEAR INTERPOLATION APPLIED IN CREATING
        //     SCALAR VALUE OF NEW MESH VERTEX in new_mesh_vertices?
        interpolate_using_multilinear
          (scalar_grid, isovalue, increment.Ptr(), num_cube_vertices,
           coord0.PtrConst(), coord1.PtrConst(), coord2.Ptr());
      }

      for (int d = 0; d < dimension; d++)
        { coord[i*dimension+d] = coord2[d]; }
    }

  }

}


// Compute position of isosurface vertices using linear interpolation.
//  - Version where all vertices lie on regular grid edges.
void IJKMCUBE::position_iso_vertices_linear
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const std::vector<EDGE_INDEX> & elist, COORD_TYPE * coord)
{
  PROCEDURE_ERROR error("position_iso_vertices_linear");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_linear_binary(scalar_grid, isovalue, elist, coord);
    break;

  case NEP:
    position_iso_vertices_linear_nep(scalar_grid, isovalue, elist, coord);
    break;

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }
}


// Compute position of isosurface vertices using linear interpolation.
// - Version with intersected edges specified by pairs of endpoints
//   in elist_endpoint[].
void IJKMCUBE::position_iso_vertices_linearB
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const std::vector<VERTEX_INDEX> & elist_endpoint, COORD_TYPE * coord)
{
  PROCEDURE_ERROR error("position_iso_vertices_linear");

  switch(isotable_type) {

  case BINARY:
    compute_isov_coord_on_line_segment_linear
      (scalar_grid, isovalue, elist_endpoint, coord);
    break;

    /* NOT YET IMPLEMENTED
  case NEP:
    position_iso_vertices_linear_nep2(scalar_grid, isovalue, vlist, coord);
    break;
    */

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }

}


// Compute position of isosurface vertices using linear interpolation.
// - Version with intersected edges specified by pairs of endpoints
//   in elist_endpoint[].
void IJKMCUBE::position_iso_vertices_linearB
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord)
{
  PROCEDURE_ERROR error("position_iso_vertices_linearB");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_linear_binaryB
      (scalar_grid, isovalue, new_mesh_vertices, elist, coord);
    break;

    /* NOT YET IMPLEMENTED
  case NEP:
    position_iso_vertices_linear_nep2(scalar_grid, isovalue, vlist, coord);
    break;
    */

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }

}

void IJKMCUBE::position_iso_vertices_multilinear
(const MC_SCALAR_GRID_BASE & scalar_grid, const SCALAR_TYPE isovalue,
 const ISOTABLE_TYPE isotable_type,
 const MC_MESH_VERTEX_LIST & new_mesh_vertices,
 const std::vector<VERTEX_INDEX> & elist, COORD_TYPE * coord)
// Compute position of isosurface vertices using multilinear interpolation
// Input is an array of vertices representing pairs of edge endpoints.
// dimension = volume dimension
// scalar_grid[iv] = scalar value at grid vertex iv
// grid_length[i] = # grid vertices along axis i
// isotable_type = table type. Only BINARY or NEP are acceptable.
// elist[] = edge list. (elist[2*i], elist[2*i+1]) are endpoints of edge i.
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  PROCEDURE_ERROR error("position_iso_vertices_multilinear");

  switch(isotable_type) {

  case BINARY:
    position_iso_vertices_multilinear_binary
      (scalar_grid, isovalue, new_mesh_vertices, elist, coord);
    break;

    /* NOT YET IMPLEMENTED
  case NEP:
    position_iso_vertices_linear_nep2(scalar_grid, isovalue, vlist, coord);
    break;
    */

  default:
    error.AddMessage("Programming error.  Illegal isosurface isotable type.");
    break;
  }

}

// Convert isosurface vertex indices to pairs of endpoints
//   representing grid edges containing isosurface vertices.
void IJKMCUBE::convert_indices_to_endpoint_pairs
(const MC_SCALAR_GRID_BASE & scalar_grid,
 const std::vector<ISO_VERTEX_INDEX> iso_vlist,
 const MERGE_DATA & merge_data,
 std::vector<VERTEX_INDEX> & elist_endpoint)
{
  const int dimension = scalar_grid.Dimension();
  IJK::ARRAY<AXIS_SIZE_TYPE> axis_increment(dimension);

  elist_endpoint.clear();

  compute_increment(scalar_grid, axis_increment.Ptr());

  for (unsigned int i = 0; i < iso_vlist.size(); i++) {
    MERGE_INDEX isov = iso_vlist[i];
    VERTEX_INDEX v0 = merge_data.GetFirstEndpoint(isov);
    MERGE_INDEX dir = merge_data.GetEdgeDir(isov);
    VERTEX_INDEX v1 = v0 + axis_increment[dir];

    elist_endpoint.push_back(v0);
    elist_endpoint.push_back(v1);
  }

}


// **************************************************
// INTERVAL VOLUME SUBROUTINES
// **************************************************

void IJKMCUBE::position_ivol_vertices_linear
(const MC_SCALAR_GRID_BASE & scalar_grid,
 const SCALAR_TYPE isovalue0, const SCALAR_TYPE isovalue1,
 const std::vector<ISO_VERTEX_INDEX> & vlist, COORD_TYPE * coord)
// compute position of interval volume vertices using linear interpolation
// scalar_grid = scalar grid data
// vlist[] = list of isosurface vertices
// coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
//   Precondition: array coord[] is preallocated to length at least
//                 dimension*vlist.size()
{
  const int dimension = scalar_grid.Dimension();
  const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
  const SCALAR_TYPE * scalar = scalar_grid.ScalarPtrConst();
  IJK::ARRAY<AXIS_SIZE_TYPE> axis_increment(dimension);
  const int numv = vlist.size();
  IJK::ARRAY<GRID_COORD_TYPE> coord0(dimension);
  IJK::ARRAY<GRID_COORD_TYPE> coord1(dimension);
  IJK::ARRAY<COORD_TYPE> coord2(dimension);
  const int VERTEX_OFFSET = 2*dimension;

  // compute_axis_increment
  compute_increment(scalar_grid, axis_increment.Ptr());

  const int num_ivolv_per_gridv = 
    get_num_ivol_vertices_per_grid_vertex(dimension);


  for (int i = 0; i < numv; i++) {
    int k = vlist[i]%num_ivolv_per_gridv;
    VERTEX_INDEX v0 = vlist[i]/num_ivolv_per_gridv;

    if (k == VERTEX_OFFSET) {
      // isosurface vertex lies on grid vertex v0
      compute_coord(v0, dimension, axis_size, coord0.Ptr());
      for (int d = 0; d < dimension; d++)
        coord[i*dimension+d] = coord0[d];
    }
    else {
      // isosurface vertex lies on grid edge
      int d = k/2;
      VERTEX_INDEX v1 = v0+axis_increment[d];

      SCALAR_TYPE s0 = scalar[v0];
      SCALAR_TYPE s1 = scalar[v1];

      compute_coord(v0, dimension, axis_size, coord0.Ptr());
      compute_coord(v1, dimension, axis_size, coord1.Ptr());

      if (k % 2 == 0) {
        // linear interpolate using isovalue0
        linear_interpolate_coord
          (dimension, s0, coord0.PtrConst(), s1, coord1.PtrConst(), 
           isovalue0, coord2.Ptr());
      }
      else {
        // linear interpolate using isovalue1
        linear_interpolate_coord
          (dimension, s0, coord0.PtrConst(), s1, coord1.PtrConst(), 
           isovalue1, coord2.Ptr());
      }
      for (int d = 0; d < dimension; d++)
        coord[i*dimension+d] = coord2[d];
    }
  }

}


