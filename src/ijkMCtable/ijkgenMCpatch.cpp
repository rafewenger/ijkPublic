/*!
 *  @file ijkgenMCpatch.cxx
 *  @brief Generate isosurface patch for a Marching Cubes lookup table.
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

#include "ijkMCtable.h"
#include "ijkgenMCpatch.h"

#include "ijksimplex.tpp"

using namespace IJKMCUBE_TABLE;

// **************************************************
// Type Declarations
// **************************************************

typedef ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE ISOSURFACE_VERTEX_TYPE;

// **************************************************
// Function Declarations
// **************************************************

extern "C" void clarkson_convex_hull
(int dimension, double * point_coord, int num_points,
 int * convex_hull_dimension, int ** simplex_vert, int * num_simplices);
extern "C" void free_simplex_vertices(int * simplex_vert);

void add_simplex_to_isosurface
(const int dimension, const int * simplex_vert, const int is,
 const int * point_list, std::vector<ISOSURFACE_VERTEX_INDEX> & isov_list);

void add_simplex_to_isosurface
(const int dimension, const int * simplex_vert, const int is,
 const int * point_list, const ISOSURFACE_VERTEX_TYPE * point_type,
 std::vector<ISOSURFACE_VERTEX_INDEX> & isov_list,
 std::vector<ISOSURFACE_VERTEX_TYPE> & isov_type);

void add_triangle_vertices_to_triangle_vertex_list
(const int iv0, const int iv1, const int iv2,
 std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list);


// **************************************************
// inline functions
// **************************************************

/*!
 *  @brief Set bits in vertex_bits[].
 *  @param polytope Isosurface table polytope.
 *  @param point_list[] List of points.
 *    @pre Polytope vertices are listed first in point list
 *         followed by polytope edge midpoints.
 *  @param ip Set bits for point point_list[ip].
 *  @param num_vert_in_point_list Number of polytope vertices
 *    in point_list[].
 *  @param vertex_bits[out] Bit vector of polytope vertices.
 */
inline void set_vertex_bits
(const ISOTABLE_POLY_BASE & polytope,
 const int * point_list, const int ip,
 const int num_vert_in_point_list, long & vertex_bits)
{
  if (ip < num_vert_in_point_list) {
    // point_list[ip] is a polytope vertex
    vertex_bits  = vertex_bits | (1L << point_list[ip]);
  }
  else {
    // point_list[ip] is a polytope vertex
    int ie = point_list[ip];

    int iv0 = polytope.EdgeEndpoint(ie, 0);
    int iv1 = polytope.EdgeEndpoint(ie, 1);

    vertex_bits = vertex_bits | (1L << iv0);
    vertex_bits = vertex_bits | (1L << iv1);
  }
}


// **************************************************
// ijkgenpatch
// **************************************************


// Generate isosurface patch for standard Marching Cubes lookup table.
// - Version that passes temporary arrays for faster processing.
void IJKMCUBE_TABLE::gen_isopatch
(const ISOTABLE_POLY_BASE & polytope,
 const int * const vertex_sign,
 std::vector<ISOSURFACE_VERTEX_INDEX> & edge_list, int & num_simplices,
 int * temp_plist, double * temp_pcoord)
{
  const int dimension = polytope.Dimension();
  int * point_list = temp_plist;
  double * point_coord = temp_pcoord;

  num_simplices = 0;

  edge_list.clear();

  // add '+' vertices to point_list[]
  int nump = 0;
  for (int jv = 0; jv < polytope.NumVertices(); jv++) {
    if (vertex_sign[jv] >= 0) {
      // vertex jv is zero or positive
      for (int ic = 0; ic < dimension; ic++) {
        int j = nump*dimension+ic;
        point_coord[j] = polytope.VertexCoord(jv, ic);
      };

      point_list[nump] = jv;
      nump++;
    };
  };
  int num_pos_vertices = nump;

  // add midpoints of bipolar edges to point_list[]
  for (int ie = 0; ie < polytope.NumEdges(); ie++) {

    int iv0 = polytope.EdgeEndpoint(ie, 0);
    int iv1 = polytope.EdgeEndpoint(ie, 1);

    if ((vertex_sign[iv0] >= 0 && vertex_sign[iv1] < 0) ||
        (vertex_sign[iv0] < 0 && vertex_sign[iv1] >= 0)) {
      // edge connects positive and negative vertex
      // add edge midpoint to point_list[]
      for (int ic = 0; ic < dimension; ic++) {
        int j = nump*dimension+ic;
        point_coord[j] = polytope.MidpointCoord(ie, ic);
      };

      point_list[nump] = ie;
      nump++;
    };
  };

  assert(nump <= polytope.NumEdges() + polytope.NumVertices());

  int * simplex_vert;
  int num_hull_simplices = 0;
  int convex_hull_dimension = 0;
  clarkson_convex_hull(dimension, point_coord, nump,
                       &convex_hull_dimension, &simplex_vert,
                       &num_hull_simplices);

  assert(nump == 0 || convex_hull_dimension == dimension);

  for (int is = 0; is < num_hull_simplices; is++) {

    bool flag_isosurface_simplex = true;
    for (int ic = 0; ic < dimension; ic++)
      if (simplex_vert[is*dimension+ic] < num_pos_vertices)
        // ignore simplices sharing a polytope vertex
        flag_isosurface_simplex = false;

    // ignore simplices lying on polytope facets
    if (flag_isosurface_simplex) {
      long vertex_bits = 0L;
      for (int ic = 0; ic < dimension; ic++) {
        int ip = simplex_vert[is*dimension+ic];
        assert (ip >= num_pos_vertices);

        set_vertex_bits(polytope, point_list, ip, num_pos_vertices,
                        vertex_bits);
      };

      for (int jf = 0; jf < polytope.NumFacets(); jf++) {
        FACET f = polytope.Facet(jf);
        if ((vertex_bits & (~f)) == 0) {
          flag_isosurface_simplex = false;
          break;
        };
      };
    };

    if (flag_isosurface_simplex) {
      add_simplex_to_isosurface
        (dimension, simplex_vert, is, point_list, edge_list);
    };
  };

  assert(edge_list.size() % dimension == 0);

  num_simplices = edge_list.size()/dimension;

  free_simplex_vertices(simplex_vert);
}



// Generate isosurface patch for standard Marching Cubes lookup table.
void IJKMCUBE_TABLE::gen_isopatch
(const ISOTABLE_POLY_BASE & polytope,
 const int * const vertex_sign,
 std::vector<ISOSURFACE_VERTEX_INDEX> & edge_list, int & num_simplices)
{
  const int nume = polytope.NumEdges();
  const int numv = polytope.NumVertices();
  const int dimension = polytope.Dimension();
  int * temp_plist = new int[nume+numv];
  double * temp_pcoord = new double[(nume+numv)*dimension];

  gen_isopatch(polytope, vertex_sign, edge_list, num_simplices,
              temp_plist, temp_pcoord);

  delete [] temp_plist;
  delete [] temp_pcoord;
}


// **************************************************
// gen_isopatch_edge_groups
// **************************************************


// Local function.
namespace {

  /// Get number of adjacent vertices with same label.
  void get_matching_label_degree
  (const ISOTABLE_POLY_BASE & polytope,
   const int * const vertex_sign,
   std::vector<int> & degree)
  {
    std::fill(degree.begin(), degree.end(), 0);

    for (int ie = 0; ie < polytope.NumEdges(); ie++) {
      const int iend0 = polytope.EdgeEndpoint(ie, 0);
      const int iend1 = polytope.EdgeEndpoint(ie, 1);

      if (vertex_sign[iend0] == vertex_sign[iend1]) {
        degree[iend0]++;
        degree[iend1]++;
      }
    }
  }


  /// Get first vertex with given sign and degree.
  int get_vertex_with_sign_and_degree
  (const int num_poly_vertices,
   const int * vertex_sign,
   const std::vector<int> & degree,
   const int signX,
   const int degreeX)
  {
    for (int iv = 0; iv < num_poly_vertices; iv++) {
      if ((vertex_sign[iv] == signX) &&
          (degree[iv] == degreeX))
        { return iv; }
    }

    return (num_poly_vertices);
  }


  /// Return true if edge is bipolar.
  bool is_bipolar
  (const ISOTABLE_POLY_BASE & polytope,
   const int * const vertex_sign,
   const EDGE_INDEX ie)
  {
    const int iend0 = polytope.EdgeEndpoint(ie, 0);
    const int iend1 = polytope.EdgeEndpoint(ie, 1);

    if (vertex_sign[iend0] != vertex_sign[iend1])
      { return true; }
    else
      { return false; }
  }


  /// Get bipolar edge that is incident on two degree 1 vertices.
  int get_bipolar_edge_incident_on_two_degree1_vertices
  (const ISOTABLE_POLY_BASE & poly,
   const int * vertex_sign,
   const std::vector<int> & degree,
   IJK::ERROR & error)
  {
    const int DEGREE1(1);

    for (int ie = 0; ie < poly.NumEdges(); ie++) {
      if (is_bipolar(poly, vertex_sign, ie)) {
        const int iend0 = poly.EdgeEndpoint(ie, 0);
        const int iend1 = poly.EdgeEndpoint(ie, 1);
        if ((degree[iend0] == DEGREE1) &&
            (degree[iend1] == DEGREE1))
          { return ie; }
      }
    }

    error.AddMessage
      ("Programming error. No bipolar edges is incident on two degree 1 vertices.");
    throw error;
  }


  // Get bipolar edge incident on degree 0 and degree 2 vertices.
  int get_bipolar_edge_incident_on_degree0_and_degree2_vertices
  (const ISOTABLE_POLY_BASE & poly,
   const int * vertex_sign,
   const std::vector<int> & degree,
   IJK::ERROR & error)
  {
    const int DEGREE0(0);
    const int DEGREE2(2);

    for (int ie = 0; ie < poly.NumEdges(); ie++) {
      if (is_bipolar(poly, vertex_sign, ie)) {
        const int iend0 = poly.EdgeEndpoint(ie, 0);
        const int iend1 = poly.EdgeEndpoint(ie, 1);
        if ((degree[iend0] == DEGREE0) &&
            (degree[iend1] == DEGREE2))
          { return ie; }
        else if ((degree[iend0] == DEGREE2) &&
                 (degree[iend1] == DEGREE0))
          { return ie; }
      }
    }

    error.AddMessage
      ("Programming error. No bipolar edge is incident on degree 0 and degree 1 vertices.");
    throw error;
  }


  // Get bipolar edge incident on vertex iv.
  // @pre Some bipolar edge is incident on vertex iv.
  int get_bipolar_edge_incident_on_vertex
  (const ISOTABLE_POLY_BASE & polytope,
   const int * const vertex_sign,
   const int iv)
  {
    for (int ie = 0; ie < polytope.NumEdges(); ie++) {
      if ((iv == polytope.EdgeEndpoint(ie, 0)) ||
          (iv == polytope.EdgeEndpoint(ie, 1))) {
        if (is_bipolar(polytope, vertex_sign, ie))
          { return ie; }
      }
    }

    // Just in case.
    return 0;
  }


  // Return negative endpoint index (0 or 1) of bipolar edge ie.
  int get_negative_endpoint_index
  (const ISOTABLE_POLY_BASE & polytope,
   const int * const vertex_sign,
   const EDGE_INDEX ie)
  {
    const int iend0 = polytope.EdgeEndpoint(ie, 0);
    if (vertex_sign[iend0] == 0)
      { return 0; }
    else
      { return 1; }
  }


  void convert_edge_vertex_list_to_cycle
  (const int num_poly_edges,
   const std::vector<ISOSURFACE_VERTEX_INDEX> & edge_vertex_list,
   const ISOSURFACE_VERTEX_INDEX istart,
   std::vector<ISOSURFACE_VERTEX_INDEX> & cycle)
  {
    const int NUM_VERT_PER_EDGE(2);
    const int num_edges =
      edge_vertex_list.size()/NUM_VERT_PER_EDGE;
    std::vector<ISOSURFACE_VERTEX_INDEX> next_vertex(num_poly_edges,0);
    IJK::PROCEDURE_ERROR error("convert_edge_vertex_list_to_cycle");

    cycle.clear();
    for (int ie = 0; ie < num_edges; ie++) {
      const int iw0 = edge_vertex_list[2*ie];
      const int iw1 = edge_vertex_list[2*ie+1];
      next_vertex[iw0] = iw1;
    }

    int iw = istart;
    for (int i = 0; i < num_edges; i++) {
      cycle.push_back(iw);
      iw = next_vertex[iw];
    }

    if (iw != istart) {
      error.AddMessage
        ("Programming error. Edges in edge_vertex_list[] do not form a cycle.");
      error.AddMessage
        ("  Number of edges in edge_vertex_list[]: ", num_edges);
      error.AddMessage("  Starting vertex: ", istart);
      throw error;
    }
  }


  // Version that uses edge_vertex_list[0] as starting edge.
  void convert_edge_vertex_list_to_cycle
  (const int num_poly_edges,
   const std::vector<ISOSURFACE_VERTEX_INDEX> & edge_vertex_list,
   std::vector<ISOSURFACE_VERTEX_INDEX> & cycle)
  {
    // Initialize.
    cycle.clear();

    if (edge_vertex_list.size() == 0) {
      // Nothing to convert.
      return;
    }

    convert_edge_vertex_list_to_cycle
      (num_poly_edges, edge_vertex_list, edge_vertex_list[0], cycle);

  }


  /*!
   *  @brief Return true if isopatch is a pentagon and a triangle.
   */
  bool is_isopatch_pentagon_and_triangle
  (const std::vector<int> & num_edges_in_connected_component)
  {
    const int TWO_COMPONENTS(2);
    const int THREE_EDGES(3);
    const int FIVE_EDGES(5);

    if (num_edges_in_connected_component.size() !=
        TWO_COMPONENTS)
      { return false; }

    if ((num_edges_in_connected_component[0] == FIVE_EDGES) &&
        (num_edges_in_connected_component[1] == THREE_EDGES))
      { return true; }
    else if ((num_edges_in_connected_component[0] == THREE_EDGES) &&
             (num_edges_in_connected_component[1] == FIVE_EDGES))
      { return true; }
    else {
      return false;
    }
  }


  bool check_num_boundary_edges
  (const std::vector<ISOSURFACE_VERTEX_INDEX> & boundary_edge_vertex_list,
   const int expected_num_edges,
   IJK::ERROR & error)
  {
    if (boundary_edge_vertex_list.size() != 2*expected_num_edges) {
      error.AddMessage
        ("Programming error. Incorrect number of boundary edges.");
      error.AddMessage
        ("  Received ", boundary_edge_vertex_list.size()/2,
         " boundary edges.");
      error.AddMessage
        ("  Expected ", expected_num_edges, " boundary edges.");
      return false;
    }

    return true;
  }


  /*!
   *  @brief Add triangle isosurface patch to triangle_vertex_list.
   *  @pre boundary_edge_vertex_list[] has exactly 3 edges.
   *  @param[out] triangle_vertex_list[] List of triangle vertices.
   *    - May already contain some triangles.
   */
  void add_isopatch_triangle
  (const ISOTABLE_CUBE_3D & cube,
   std::vector<ISOSURFACE_VERTEX_INDEX> & boundary_edge_vertex_list,
   const std::vector<int> & degree,
   std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list)
  {
    const int NUM_TRIANGLE_EDGES(3);
    std::vector<ISOSURFACE_VERTEX_INDEX> cycle;
    IJK::PROCEDURE_ERROR error("add_isopatch_triangle");

    if (!check_num_boundary_edges
        (boundary_edge_vertex_list, NUM_TRIANGLE_EDGES, error))
      { throw error; }


    convert_edge_vertex_list_to_cycle
      (cube.NumEdges(), boundary_edge_vertex_list, cycle);

    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[1], cycle[2], triangle_vertex_list);
  }


  /*!
   *  @brief Generate pentagon isosurface patch for Marching Cubes
   *  edge groups lookup table,
   *  @pre boundary_edge_vertex_list[] has exactly 5 edges.
   */
  void gen_isopatch_pentagon
  (const ISOTABLE_CUBE_3D & cube,
   const int * const vertex_sign,
   std::vector<ISOSURFACE_VERTEX_INDEX> & boundary_edge_vertex_list,
   const std::vector<int> & degree,
   const bool flag_separate_neg,
   std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list)
  {
    const int DEGREE2(2);
    const int THREE(3);
    const int FOUR(4);
    const int NUM_PENTAGON_EDGES(5);
    const int num_cube_vertices = cube.NumVertices();
    std::vector<ISOSURFACE_VERTEX_INDEX> cycle;
    IJK::PROCEDURE_ERROR error("gen_isopatch_pentagon");

    // Initialize.
    triangle_vertex_list.clear();

    if (!check_num_boundary_edges
        (boundary_edge_vertex_list, NUM_PENTAGON_EDGES, error))
      { throw error; }

    const int num_positive_vertices =
      std::count(vertex_sign, vertex_sign+num_cube_vertices, 1);
    const int num_negative_vertices =
      num_cube_vertices - num_positive_vertices;

    int iv0;
    if (num_negative_vertices == THREE) {
      iv0 = get_vertex_with_sign_and_degree
        (num_cube_vertices, vertex_sign, degree, -1, DEGREE2);
    }
    else if (num_positive_vertices == THREE) {
      iv0 = get_vertex_with_sign_and_degree
        (num_cube_vertices, vertex_sign, degree, 1, DEGREE2);
    }
    else if (num_negative_vertices == FOUR) {
      if (flag_separate_neg) {
        iv0 = get_vertex_with_sign_and_degree
          (num_cube_vertices, vertex_sign, degree, -1, DEGREE2);
      }
      else {
        iv0 = get_vertex_with_sign_and_degree
          (num_cube_vertices, vertex_sign, degree, 1, DEGREE2);
      }
    }
    else {
      error.AddMessage
        ("Programming error. Unable to find degree 2 cube vertex.");
      throw error;
    }

    const int ie0 =
      get_bipolar_edge_incident_on_vertex(cube, vertex_sign, iv0);

    convert_edge_vertex_list_to_cycle
      (cube.NumEdges(), boundary_edge_vertex_list, ie0, cycle);

    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[1], cycle[2], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[2], cycle[3], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[3], cycle[4], triangle_vertex_list);
  }


  /*!
   *  @brief Generate pentagon isosurface patch for Marching Cubes
   *  edge groups lookup table,
   *  - Version that returns num_triangles.
   *  @pre boundary_edge_vertex_list[] has exactly 5 edges.
   */
  void gen_isopatch_pentagon
  (const ISOTABLE_CUBE_3D & cube,
   const int * const vertex_sign,
   std::vector<ISOSURFACE_VERTEX_INDEX> & boundary_edge_vertex_list,
   const std::vector<int> & degree,
   const bool flag_separate_neg,
   std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list,
   int & num_triangles)
  {
    gen_isopatch_pentagon
      (cube, vertex_sign, boundary_edge_vertex_list,
       degree, flag_separate_neg, triangle_vertex_list);
    num_triangles = 3;
  }


  /*!
   *  @brief Generate hexagon isosurface patch for Marching Cubes
   *  edge groups lookup table,
   *  @pre boundary_edge_vertex_list[] has exactly 6 edges.
   *  @pre No cube vertex has degree 3.
   */
  void gen_isopatch_hexagon
  (const ISOTABLE_CUBE_3D & cube,
   const int * const vertex_sign,
   std::vector<ISOSURFACE_VERTEX_INDEX> &
   boundary_edge_vertex_list,
   const std::vector<int> & degree,
   std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list,
   int & num_triangles)
  {
    const int DEGREE2(2);
    const int FOUR(4);
    const int NUM_HEXAGON_EDGES(6);
    const int num_cube_vertices = cube.NumVertices();
    std::vector<ISOSURFACE_VERTEX_INDEX> cycle;
    IJK::PROCEDURE_ERROR error("gen_isopatch_hexagon");

    // Initialize.
    num_triangles = 0;
    triangle_vertex_list.clear();

    const int num_positive_vertices =
      std::count(vertex_sign, vertex_sign+num_cube_vertices, 1);

    if (num_positive_vertices != FOUR) {
      error.AddMessage
        ("Programming error. Incorrect number of positive vertices.");
      error.AddMessage
        ("  num_positive_vertices: ", num_positive_vertices);
      error.AddMessage
        ("  Number of positive vertices should be ", FOUR, ".");
      throw error;
    }

    const int ie0 =
      get_bipolar_edge_incident_on_two_degree1_vertices
      (cube, vertex_sign, degree, error);

    if (!check_num_boundary_edges
        (boundary_edge_vertex_list, NUM_HEXAGON_EDGES, error))
      { throw error; }

    convert_edge_vertex_list_to_cycle
      (cube.NumEdges(), boundary_edge_vertex_list, ie0, cycle);

    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[1], cycle[2], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[2], cycle[3], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[3], cycle[4], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[4], cycle[5], triangle_vertex_list);
    num_triangles = 4;
  }


  /*!
   *  @brief Generate septagon isosurface patch for Marching Cubes
   *  edge groups lookup table,
   *  @pre boundary_edge_vertex_list[] has exactly 7 edges.
   *  @pre No cube vertex has degree 3.
   */
  void gen_isopatch_septagon
  (const ISOTABLE_CUBE_3D & cube,
   const int * const vertex_sign,
   std::vector<ISOSURFACE_VERTEX_INDEX> & boundary_edge_vertex_list,
   const std::vector<int> & degree,
   std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list,
   int & num_triangles)
  {
    const int DEGREE1(1);
    const int THREE(3);
    const int NUM_SEPTAGON_EDGES(7);
    const int num_cube_vertices = cube.NumVertices();
    std::vector<ISOSURFACE_VERTEX_INDEX> cycle;
    IJK::PROCEDURE_ERROR error("gen_isopatch_septagon");

    // Initialize.
    num_triangles = 0;
    triangle_vertex_list.clear();

    const int num_positive_vertices =
      std::count(vertex_sign, vertex_sign+num_cube_vertices, 1);
    const int num_negative_vertices =
      num_cube_vertices - num_positive_vertices;

    if (!check_num_boundary_edges
        (boundary_edge_vertex_list, NUM_SEPTAGON_EDGES, error))
      { throw error; }

    if ((num_positive_vertices != THREE) &&
        (num_negative_vertices != THREE)) {
      error.AddMessage
        ("Programming error. Incorrect number of positive/negative vertices.");
      error.AddMessage
        ("  num_positive_vertices: ", num_positive_vertices);
      error.AddMessage
        ("  num_negative_vertices: ", num_negative_vertices);
      error.AddMessage
        ("  Either number of positive or number of negative vertices");
      error.AddMessage("    should be ", THREE, ".");
      throw error;
    }

    const int ie0 =
      get_bipolar_edge_incident_on_degree0_and_degree2_vertices
      (cube, vertex_sign, degree, error);

    convert_edge_vertex_list_to_cycle
      (cube.NumEdges(), boundary_edge_vertex_list, ie0, cycle);

    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[1], cycle[3], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[3], cycle[4], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[0], cycle[4], cycle[6], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[1], cycle[2], cycle[3], triangle_vertex_list);
    add_triangle_vertices_to_triangle_vertex_list
      (cycle[4], cycle[5], cycle[6], triangle_vertex_list);
    num_triangles = 5;
  }


  /*!
   *  @brief Generate pentagon isosurface patch and triangle
   *  isosurface patch for Marching Cubes edge groups lookup table,
   *  @pre boundary_edge_vertex_list[] has exactly 5 edges.
   */
  void gen_isopatch_pentagon_and_triangle
  (const ISOTABLE_CUBE_3D & cube,
   const int * const vertex_sign,
   std::vector<ISOSURFACE_VERTEX_INDEX> & boundary_edge_vertex_list,
   const std::vector<int> & degree,
   const bool flag_separate_neg,
   std::vector<int> & boundary_edge_component,
   std::vector<int> & num_edges_in_connected_component,
   std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list,
   int & num_triangles)
  {
    const int TWO_COMPONENTS(2);
    const int NUM_VERT_PER_EDGE(2);
    const int NUM_TRIANGLE_EDGES(3);
    const int NUM_PENTAGON_EDGES(5);
    const int num_cube_vertices = cube.NumVertices();
    std::vector<ISOSURFACE_VERTEX_INDEX> cycle;
    IJK::PROCEDURE_ERROR error("gen_isopatch_pentagon_and_triangle");

    // Initialize.
    num_triangles = 0;
    triangle_vertex_list.clear();

    if (!check_num_boundary_edges
        (boundary_edge_vertex_list,
         NUM_PENTAGON_EDGES+NUM_TRIANGLE_EDGES, error))
      { throw error; }

    if (num_edges_in_connected_component.size() !=
        TWO_COMPONENTS) {
      error.AddMessage
        ("Programming error. Incorrect number of connected components.");
      error.AddMessage
        ("  Number of connected components: ",
         num_edges_in_connected_component.size());
      error.AddMessage
        ("  Expected number of connected components: ",
         TWO_COMPONENTS);
      throw error;
    }

    int icomp3, icomp5;
    if ((num_edges_in_connected_component[0] == NUM_PENTAGON_EDGES) &&
        (num_edges_in_connected_component[1] == NUM_TRIANGLE_EDGES)) {
      icomp5 = 0;
      icomp3 = 1;
    }
    else if ((num_edges_in_connected_component[0] == NUM_TRIANGLE_EDGES) &&
             (num_edges_in_connected_component[1] == NUM_PENTAGON_EDGES)) {
      icomp3 = 0;
      icomp5 = 1;
    }
    else {
      error.AddMessage
        ("Programming error. Incorrect number of edges in connected components.");
      error.AddMessage
        ("  Number of edges in component 0: ",
         num_edges_in_connected_component[0]);
      error.AddMessage
        ("  Number of edges in component 1: ",
         num_edges_in_connected_component[1]);
      error.AddMessage
        ("  Expected number of edges: ", NUM_TRIANGLE_EDGES,
         " in one component and ", NUM_PENTAGON_EDGES,
         " in the other.");
      throw error;
    }

    std::vector<ISOSURFACE_VERTEX_INDEX> triangle_edge_vertex_list;
    std::vector<ISOSURFACE_VERTEX_INDEX> pentagon_edge_vertex_list;

    IJK::get_simplices_in_connected_component
      (boundary_edge_vertex_list, NUM_VERT_PER_EDGE,
       boundary_edge_component, icomp3,
       triangle_edge_vertex_list);

    IJK::get_simplices_in_connected_component
      (boundary_edge_vertex_list, NUM_VERT_PER_EDGE,
       boundary_edge_component, icomp5,
       pentagon_edge_vertex_list);

    gen_isopatch_pentagon
      (cube, vertex_sign, pentagon_edge_vertex_list, degree,
       flag_separate_neg, triangle_vertex_list);

    add_isopatch_triangle
      (cube, triangle_edge_vertex_list,
       degree, triangle_vertex_list);

    num_triangles = 4;
  }


}




// Generate isosurface patch for Marching Cubes edge groups lookup table.
void IJKMCUBE_TABLE::gen_isopatch_edge_groups
(const ISOTABLE_CUBE_3D & cube,
 const int * const vertex_sign,
 std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list,
 int & num_triangles,
 int * temp_plist, double * temp_pcoord)
{
  const int DIM3(3);
  const int NUM_VERT_PER_TRIANGLE(3);
  const int NUM_VERT_PER_EDGE(2);
  const int ONE_COMPONENT(1);
  const int TWO_COMPONENTS(2);
  const int THREE_EDGES(3);
  const int FIVE_EDGES(5);
  const int SIX_EDGES(6);
  const int SEVEN_EDGES(7);
  const int DEGREE_THREE(3);
  const int num_poly_vert = cube.NumVertices();
  const bool flag_separate_neg = true;
  std::vector<ISOSURFACE_VERTEX_INDEX> boundary_edge_vertex_list;
  std::vector<int> boundary_edge_component;
  std::vector<int> triangle_containing_boundary_edge;
  std::vector<int> boundary_edge_swap_parity;
  std::vector<int> num_edges_in_connected_component;
  int num_components;
  std::vector<int> degree(num_poly_vert, 0);
  IJK::PROCEDURE_ERROR error("ijkgenpatch_edge_groups");

  // *** DEBUG ***
  /*
  using namespace std;
  cerr << "In " << __func__ << endl;
  */

  if (cube.Dimension() != DIM3) {
    error.AddMessage("Programming error. Dimension not three.");
    error.AddMessage
      ("  Routine ijkgenpatch_edge_groups() expects cube dimension 3.");
    error.AddMessage
      ("  Cube dimension: ", cube.Dimension());
    throw error;
  }

  // Get standard isosurface patch.
  gen_isopatch
    (cube, vertex_sign, triangle_vertex_list, num_triangles,
     temp_plist, temp_pcoord);

  // *** DEBUG ***
  /*
  using namespace std;
  cerr << "  Calling orient_all_simplices()." << endl;
  */

  IJK::orient_all_simplices
    (triangle_vertex_list, NUM_VERT_PER_TRIANGLE);

  // *** DEBUG ***
  /*
  using namespace std;
  cerr << "  Returned from gen_isopatch." << endl;
  cerr << "  Triangle vertex list:" << endl;
  IJK::print_list_as_int(cerr, "  ", triangle_vertex_list, "\n");
  cerr << "    num_triangles: " << num_triangles << endl;
  */

  IJK::get_simplex_boundary_facets
    (triangle_vertex_list, NUM_VERT_PER_TRIANGLE,
     boundary_edge_vertex_list,
     triangle_containing_boundary_edge,
     boundary_edge_swap_parity);

  IJK::reorient_simplices
    (boundary_edge_vertex_list, NUM_VERT_PER_EDGE,
     boundary_edge_swap_parity);

  const int num_boundary_edges =
    boundary_edge_vertex_list.size()/NUM_VERT_PER_EDGE;

  // *** DEBUG ***/
  /*
  using namespace std;
  cerr << "  Num boundary edges: "
       << triangle_containing_boundary_edge.size()
       << endl;
  */

  // *** DEBUG ***
  /*
  using namespace std;
  cerr << "  Boundary edge vertex list:" << endl;
  IJK::print_list_as_int(cerr, "  ", boundary_edge_vertex_list, "\n");
  */

  IJK::get_connected_components_in_simplicial_complex
    (boundary_edge_vertex_list, NUM_VERT_PER_EDGE,
     boundary_edge_component, num_components);

  IJK::get_num_simplices_in_each_connected_component
    (boundary_edge_component, num_boundary_edges,
     num_edges_in_connected_component, num_components);

  // *** DEBUG ***
  /*
  using namespace std;
  cerr << "  num_boundary_edges: " << num_boundary_edges << endl;
  cerr << "  num_components: " << num_components << endl;
  IJK::print_list(cerr, "  boundary_edge_component: ",
                  boundary_edge_component, "\n");
  IJK::print_list(cerr, "  num_edges_in_connected_component: ",
                  num_edges_in_connected_component, "\n");
  */


  // *** DEBUG ***
  /*
  using namespace std;
  cerr << "  Returned from get_num_simplices_in_each_connected_component()."
       << endl;
  cerr << "  num_edges_in_connected_component.size(): "
       << num_edges_in_connected_component.size() << endl;
  */

  get_matching_label_degree(cube, vertex_sign, degree);
  const auto degree3_iter = std::find
    (degree.begin(), degree.end(), DEGREE_THREE);

  if (num_components == ONE_COMPONENT) {

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  num_edges_in_connected_component[0]: "
         << num_edges_in_connected_component[0] << endl;
    */

    if (num_edges_in_connected_component[0] == FIVE_EDGES) {
      gen_isopatch_pentagon
        (cube, vertex_sign, boundary_edge_vertex_list, degree,
         flag_separate_neg, triangle_vertex_list, num_triangles);
    }
    else if (num_edges_in_connected_component[0] == SEVEN_EDGES) {
      gen_isopatch_septagon
        (cube, vertex_sign, boundary_edge_vertex_list, degree,
         triangle_vertex_list, num_triangles);
    }
    else if (num_edges_in_connected_component[0] == SIX_EDGES &&
             degree3_iter == degree.end()) {
      gen_isopatch_hexagon
        (cube, vertex_sign, boundary_edge_vertex_list, degree,
         triangle_vertex_list, num_triangles);
    }
    else {
      return;
    }
  }
  else if (is_isopatch_pentagon_and_triangle
           (num_edges_in_connected_component)) {
    gen_isopatch_pentagon_and_triangle
      (cube, vertex_sign, boundary_edge_vertex_list, degree,
       flag_separate_neg,
       boundary_edge_component, num_edges_in_connected_component,
       triangle_vertex_list, num_triangles);
  }

  return;
}


// Generate isosurface patch for Marching Cubes edge groups lookup table.
// @pre Cube dimension is 3.
// - Version without temporary arrays.
void IJKMCUBE_TABLE::gen_isopatch_edge_groups
(const ISOTABLE_CUBE_3D & cube,
 const int * const vertex_sign,
 std::vector<ISOSURFACE_VERTEX_INDEX> & edge_list, int & num_simplices)
{
  const int nume = cube.NumEdges();
  const int numv = cube.NumVertices();
  const int dimension = cube.Dimension();
  int * temp_plist = new int[nume+numv];
  double * temp_pcoord = new double[(nume+numv)*dimension];

  gen_isopatch_edge_groups
    (cube, vertex_sign, edge_list, num_simplices,
     temp_plist, temp_pcoord);

  delete [] temp_plist;
  delete [] temp_pcoord;
}


// **************************************************
// ijkgenpatch_nep
// **************************************************


// Generate isosurface patch for Marching Cubes
//   NEP (negative, equals, positive) lookup table.
// - Version that passes temporary arrays for faster processing.
void IJKMCUBE_TABLE::gen_isopatch_nep
(const ISOTABLE_POLY_BASE & polytope,
 const int * const vertex_sign,
 std::vector<ISOSURFACE_VERTEX_INDEX> & isov_list,
 std::vector<ISOSURFACE_VERTEX_TYPE> & isov_type,
 int & num_simplices,
 int * temp_plist, ISOSURFACE_VERTEX_TYPE * temp_type,
 double * temp_pcoord)
{
  const int dimension = polytope.Dimension();
  const int poly_numv = polytope.NumVertices();
  const int poly_nume = polytope.NumEdges();
  int * point_list = temp_plist;
  ISOSURFACE_VERTEX_TYPE * point_type = temp_type;
  double * point_coord = temp_pcoord;

  num_simplices = 0;

  isov_list.clear();
  isov_type.clear();

  if (poly_nume + poly_numv == 0) { return; };

  // add '+' or '0' vertices to point_list[]
  int nump = 0;
  for (int jv = 0; jv < poly_numv; jv++) {
    if (vertex_sign[jv] >= 0) {
      // vertex jv is zero or positive
      for (int ic = 0; ic < dimension; ic++) {
        int j = nump*dimension+ic;
        point_coord[j] = polytope.VertexCoord(jv, ic);
      };

      point_list[nump] = jv;
      point_type[nump] = ISOSURFACE_VERTEX::VERTEX;
      nump++;
    };
  };
  int num_vert_in_point_list = nump;

  // add midpoints of bipolar edges to point_list[]
  for (int ie = 0; ie < polytope.NumEdges(); ie++) {

    int iv0 = polytope.EdgeEndpoint(ie, 0);
    int iv1 = polytope.EdgeEndpoint(ie, 1);

    if ((vertex_sign[iv0] > 0 && vertex_sign[iv1] < 0) ||
        (vertex_sign[iv0] < 0 && vertex_sign[iv1] > 0)) {
      // edge connects positive and negative vertex
      // add edge midpoint to point_list[]
      for (int ic = 0; ic < dimension; ic++) {
        int j = nump*dimension+ic;
        point_coord[j] = polytope.MidpointCoord(ie, ic);
      };

      point_list[nump] = ie;
      point_type[nump] = ISOSURFACE_VERTEX::EDGE;
      nump++;
    };
  };

  assert(nump <= poly_nume + poly_numv);

  int * simplex_vert = NULL;
  int num_hull_simplices = 0;
  int convex_hull_dimension = 0;
  clarkson_convex_hull(dimension, point_coord, nump,
                       &convex_hull_dimension, &simplex_vert,
                       &num_hull_simplices);

  if (convex_hull_dimension == dimension) {

    for (int is = 0; is < num_hull_simplices; is++) {

      // ignore simplices lying on polytope facets
      long vertex_bits = 0L;
      for (int ic = 0; ic < dimension; ic++) {
        int ip = simplex_vert[is*dimension+ic];

        set_vertex_bits(polytope, point_list, ip, num_vert_in_point_list,
                        vertex_bits);
      };

      bool flag_isosurface_simplex = true;
      for (int jf = 0; jf < polytope.NumFacets(); jf++) {
        FACET f = polytope.Facet(jf);
        if ((vertex_bits & (~f)) == 0) {
          flag_isosurface_simplex = false;
          break;
        };
      };

      if (flag_isosurface_simplex) {

        // *** DEBUG ***
        /*
        using namespace std;
        int ip0 = simplex_vert[is*dimension];
        int ip1 = simplex_vert[is*dimension+1];
        int isov0 = point_list[ip0];
        int isov1 = point_list[ip1];
        cerr << "  Adding simplex: "
          << isov0 << " (type " << int(point_type[ip0])
          << "), " << isov1 << " (type " << int(point_type[ip1])
          << ")." << endl;
          */

        add_simplex_to_isosurface
          (dimension, simplex_vert, is, point_list, point_type,
           isov_list, isov_type);
      };
    };

  }
  else if (convex_hull_dimension+1 == dimension) {

    // check if points lie on polytope facet
    long vertex_bits = 0L;
    for (int ip = 0; ip < nump; ip++) {
      set_vertex_bits(polytope, point_list, ip, num_vert_in_point_list,
                      vertex_bits);
    };

    bool flag_in_facet =  false;
    for (int jf = 0; jf < polytope.NumFacets(); jf++) {
      FACET f = polytope.Facet(jf);
      if ((vertex_bits & (~f)) == 0) {
        flag_in_facet = true;
        break;
      };
    };

    if (flag_in_facet && (nump < poly_numv + poly_nume)) {
      // points lie on a polytope facet
      // construct triangulation of the convex hull of the points
      //   (not just the boundary of the convex hull)

      // add another vertex to increase the hull dimension
      int lastp = nump;
      nump++;

      for (int jv = 0; jv < poly_numv; jv++) {

        if (vertex_sign[jv] < 0) {
          // vertex jv is negative
          for (int ic = 0; ic < dimension; ic++) {
            int j = lastp*dimension+ic;
            point_coord[j] = polytope.VertexCoord(jv, ic);
          };

          // add vertex jv to the point list
          point_list[lastp] = jv;
          point_type[lastp] = ISOSURFACE_VERTEX::VERTEX;

          if (simplex_vert != NULL) {
            free_simplex_vertices(simplex_vert);
            simplex_vert = NULL;
          };

          clarkson_convex_hull(dimension, point_coord, nump,
                               &convex_hull_dimension, &simplex_vert,
                               &num_hull_simplices);

          if (convex_hull_dimension == dimension) {

            for (int is = 0; is < num_hull_simplices; is++) {

              // ignore simplices containing vertex point_list[lastp]
              bool flag_isosurface_simplex = true;
              for (int ic = 0; ic < dimension; ic++) {
                // ignore simplices incident on vertex point_list[lastp]
                int ip = simplex_vert[is*dimension+ic];
                if (ip == lastp) {
                  flag_isosurface_simplex = false;
                }
              };

              if (flag_isosurface_simplex) {
                add_simplex_to_isosurface
                  (dimension, simplex_vert, is, point_list, point_type,
                   isov_list, isov_type);
              };

            };

            break;
          }
        }
      }
    }
  }

  assert(isov_list.size() % dimension == 0);

  num_simplices = isov_list.size()/dimension;

  if (simplex_vert != NULL) {
    free_simplex_vertices(simplex_vert);
    simplex_vert = NULL;
  }
}


// Generate isosurface patch for Marching Cubes
//   NEP (negative, equals, positive) lookup table.
// - Version without temporary arrays.
void IJKMCUBE_TABLE::gen_isopatch_nep
(const ISOTABLE_POLY_BASE & polytope,
 const int * const vertex_sign,
 std::vector<ISOSURFACE_VERTEX_INDEX> & isov_list,
 std::vector<ISOSURFACE_VERTEX_TYPE> & isov_type,
 int & num_simplices)
{
  int nume = polytope.NumEdges();
  int numv = polytope.NumVertices();
  int dimension = polytope.Dimension();
  int * temp_plist = new int[nume+numv];
  ISOSURFACE_VERTEX_TYPE * temp_type = new ISOSURFACE_VERTEX_TYPE[nume+numv];
  double * temp_pcoord = new double[(nume+numv)*dimension];

  gen_isopatch_nep(polytope, vertex_sign, isov_list, isov_type, num_simplices,
                   temp_plist, temp_type, temp_pcoord);

  delete [] temp_pcoord;
  delete [] temp_type;
  delete [] temp_plist;
}


// **************************************************
// Subroutines
// **************************************************

// Add simplex "is" to isosurface.
// dimension = volume dimension
// simplex_vert[] = array of simplex vertices
// simplex_vert[is*dimension+j] = j'th vertex of simplex is
// point_list[] = list of points
// isov_list[] = array of isosurface simplex vertices
void add_simplex_to_isosurface
(const int dimension, const int * simplex_vert, const int is,
 const int * point_list, std::vector<ISOSURFACE_VERTEX_INDEX> & isov_list)
{
  for (int k = 0; k < dimension; k++) {
    int ip = simplex_vert[is*dimension+k];
    int isov = point_list[ip];
    isov_list.push_back(isov);
  }
}


// Add simplex "is" to isosurface.
// dimension = volume dimension
// simplex_vert[] = array of simplex vertices
// simplex_vert[is*dimension+j] = j'th vertex of simplex is
// point_list[] = list of points
// point_type[k] = type of point_list[k]
// isov_list[] = array of isosurface simplex vertices
// isov_type[k] = type of isov_list[k]
void add_simplex_to_isosurface
(const int dimension, const int * simplex_vert, const int is,
 const int * point_list, const ISOSURFACE_VERTEX_TYPE * point_type,
 std::vector<ISOSURFACE_VERTEX_INDEX> & isov_list,
 std::vector<ISOSURFACE_VERTEX_TYPE> & isov_type)
{
  for (int k = 0; k < dimension; k++) {
    int ip = simplex_vert[is*dimension+k];
    int isov = point_list[ip];
    ISOSURFACE_VERTEX_TYPE vtype = point_type[ip];

    isov_list.push_back(isov);
    isov_type.push_back(vtype);
  }
}


// Add triangle vertices to triangle_vertex_list.
void add_triangle_vertices_to_triangle_vertex_list
(const int iv0, const int iv1, const int iv2,
 std::vector<ISOSURFACE_VERTEX_INDEX> & triangle_vertex_list)
{
  triangle_vertex_list.push_back(iv0);
  triangle_vertex_list.push_back(iv1);
  triangle_vertex_list.push_back(iv2);
}
