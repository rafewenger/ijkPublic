/*!
 *  @file ijktriangulate_poly2D.tpp
 *  @brief ijk templates for triangulating 2D polygons.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2010-2024 Rephael Wenger

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

#ifndef _IJKTRIANGULATE_POLY2D_
#define _IJKTRIANGULATE_POLY2D_


#include "ijkmesh.tpp"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>


namespace IJK {

  // **************************************************
  //! @name TRIANGULATE POLYGONS
  // **************************************************

  ///@{

  /*!
   *  @brief Triangulate a polygon by adding diagonals from vertex poly_vert[0].
   *  - Creates fan triangulation from vertex poly_vert[0].
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - poly_vert[0] is a vertex of all new triangles.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename NTYPE, typename VTYPE0, typename VTYPE1>
  void triangulate_polygon_from_v0
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, 
   std::vector<VTYPE1> & tri_vert)
  {
    VTYPE0 v0 = poly_vert[0];
    for (NTYPE i = 1; i+1 < num_poly_vert; i++) {
      add_triangle_vertices(v0, poly_vert[i], poly_vert[i+1], tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a polygon by adding diagonals from vertex poly_vert[0].
   *  - Reverse triangle orientations.
   *  - Creates fan triangulation from vertex poly_vert[0].
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - poly_vert[0] is a vertex of all new triangles.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename NTYPE, typename VTYPE0, typename VTYPE1>
  void triangulate_polygon_from_v0_reverse_orient
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, 
   std::vector<VTYPE1> & tri_vert)
  {
    VTYPE0 v0 = poly_vert[0];
    for (NTYPE i = 1; i+1 < num_poly_vert; i++) {
      add_triangle_vertices(v0, poly_vert[i+1], poly_vert[i], tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a polygon by adding diagonals from vertex poly_vert[0].
   *  - Creates fan triangulation from vertex poly_vert[0].
   *  - Version with boolean flag_reverse_orient.
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_reverse_orient If true, reverse triangle orientations.
   */
  template <typename NTYPE, typename VTYPE0, typename VTYPE1>
  void triangulate_polygon_from_v0
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
   const bool flag_reverse_orient, std::vector<VTYPE1> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_polygon_from_v0_reverse_orient
	(num_poly_vert, poly_vert, tri_vert);
    }
    else {
      triangulate_polygon_from_v0
	(num_poly_vert, poly_vert, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a polygon using diagonals from vertex poly_vert[vloc0].
   *  - Creates fan triangulation from vertex poly_vert[vloc0].
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - poly_vert[vloc0] is a vertex of all new triangles.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename NTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_polygon_from_vertex
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
   const VTYPE1 vloc0,
   std::vector<VTYPE2> & tri_vert)
  {
    const VTYPE0 v0 = poly_vert[vloc0];
    NTYPE i1 = (vloc0+1)%num_poly_vert;
    NTYPE i2 = (i1+1)%num_poly_vert;
    while (i2 != NTYPE(vloc0)) {
      add_triangle_vertices(v0, poly_vert[i1], poly_vert[i2], tri_vert);
      i1 = i2;
      i2 = (i1+1)%num_poly_vert;
    }
  }


  /*!
   *  @brief Triangulate a polygon using diagonals from vertex poly_vert[vloc0].
   *  - Reverse triangle orientations.
   *  - Creates fan triangulation from vertex poly_vert[vloc0].
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - poly_vert[vloc0] is a vertex of all new triangles.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename NTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_polygon_from_vertex_reverse_orient
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
   const VTYPE1 vloc0,
   std::vector<VTYPE2> & tri_vert)
  {
    VTYPE0 v0 = poly_vert[vloc0];
    NTYPE i1 = (vloc0+1)%num_poly_vert;
    NTYPE i2 = (i1+1)%num_poly_vert;
    while (i2 != vloc0) {
      add_triangle_vertices(v0, poly_vert[i2], poly_vert[i1], tri_vert);
      i1 = i2;
      i2 = (i1+1)%num_poly_vert;
    }
  }


  /*!
   *  @brief Triangulate a polygon from an interior vertex w0.
   *  - Add a triangle between each polygon edge and vertex w0.
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2>
  void triangulate_polygon_from_interior_vertex
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, const VTYPE1 w0,
   std::vector<VTYPE2> & tri_vert)
  {
    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
      NTYPE i1 = (i0+1)%num_poly_vert;
      add_triangle_vertices(w0, poly_vert[i0], poly_vert[i1], tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a polygon from an interior vertex w0.
   *  - Add a triangle between each polygon edge and vertex w0.
   *  - Reverse triangle orientation.
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2>
  void triangulate_polygon_from_interior_vertex_reverse_orient
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, const VTYPE1 w0,
   std::vector<VTYPE2> & tri_vert)
  {
    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
      NTYPE i1 = (i0+1)%num_poly_vert;
      add_triangle_vertices(w0, poly_vert[i1], poly_vert[i0], tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a polygon from an interior vertex w0.
   *  - Add a triangle between each polygon edge and vertex w0.
   *  - Version using boolean flag_reverse_orient.
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_reverse_orient If true, reverse triangle orientations.
   */
  template <typename NTYPE, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2>
  void triangulate_polygon_from_interior_vertex
  (const NTYPE num_poly_vert, const VTYPE0 * poly_vert, const VTYPE1 w0,
   const bool flag_reverse_orient, std::vector<VTYPE2> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_polygon_from_interior_vertex_reverse_orient
        (num_poly_vert, poly_vert, w0, tri_vert);
    }
    else {
      triangulate_polygon_from_interior_vertex
        (num_poly_vert, poly_vert, w0, tri_vert);
    }
  }

  ///@}


  // **************************************************
  //! @name TRIANGULATE QUADRILATERALS
  // **************************************************

  ///@{

  /*!
   *  @brief Triangulate a quadrilateral using diagonal from vertex quad_vert[iloc0].
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename ITYPE>
  void triangulate_quad_from_vertex
  (const VTYPE0 * quad_vert, const ITYPE iloc0,
   std::vector<VTYPE1> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_QUAD(4);

    triangulate_polygon_from_vertex
      (NUM_VERT_PER_QUAD, quad_vert, iloc0, tri_vert);
  }


  /*!
   *  @brief Triangulate a quadrilateral using diagonal from vertex quad_vert[iloc0].
   *  - Version with boolean flag_reverse_orient.
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_reverse_orient If true, reverse triangle orientation.
   */
  template <typename VTYPE0, typename VTYPE1, typename ITYPE>
  void triangulate_quad_from_vertex
  (const VTYPE0 * quad_vert, const ITYPE iloc0,
   const bool flag_reverse_orient, std::vector<VTYPE1> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_QUAD(4);

    if (flag_reverse_orient) {
      triangulate_polygon_from_vertex_reverse_orient
        (NUM_VERT_PER_QUAD, quad_vert, iloc0, tri_vert);
    }
    else {
      triangulate_polygon_from_vertex
        (NUM_VERT_PER_QUAD, quad_vert, iloc0, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a quadrilateral using diagonal02.
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPEQ, typename VTYPET>
  void triangulate_quad_using_diagonal02
  (const VTYPEQ * quad_vert, std::vector<VTYPET> & tri_vert)
  {
    typedef typename std::vector<VTYPET>::size_type SIZE_TYPE;
    
    const VTYPEQ v0 = quad_vert[0];
    const VTYPEQ v1 = quad_vert[1];
    const VTYPEQ v2 = quad_vert[2];
    const VTYPEQ v3 = quad_vert[3];

    const SIZE_TYPE n = tri_vert.size();
    tri_vert.resize(n+6);

    // Add first triangle.
    tri_vert[n] = v0;
    tri_vert[n+1] = v1;
    tri_vert[n+2] = v2;

    // Add second triangle
    tri_vert[n+3] = v0;
    tri_vert[n+4] = v2;
    tri_vert[n+5] = v3;
  }

  
  /*!
   *  @brief Triangulate a list of quadrilaterals using diagonal (v0,v2).
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename NTYPEQ, typename VTYPE0, typename VTYPE1,
            typename NTYPE_TRI>
  void triangulate_quad_list_using_diagonal02
  (const VTYPE0 * quad_vert, const NTYPEQ num_quad,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad)
  {
    const int NUM_VERT_PER_QUAD(4);
    const NTYPEQ list_len = num_quad*NUM_VERT_PER_QUAD;
    num_triangulated_quad = 0;

    for (NTYPEQ k = 0; k < list_len; k += NUM_VERT_PER_QUAD) {
      triangulate_quad_using_diagonal02(quad_vert+k, tri_vert);
    }
    num_triangulated_quad = num_quad;
  }
    
    
  /*!
   *  @brief Triangulate a list of quadrilaterals using diagonal (v0,v2).
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  - Skip any quads that have been collapsed to a single vertex.
   *    - (Still triangulates quads that are collapsed to an edge
   *       or two edges.)
   */
  template <typename NTYPEQ, typename VTYPE0, typename VTYPE1,
            typename NTYPE_TRI>
  void triangulate_quad_list_using_diagonal02_skip_collapsed
  (const VTYPE0 * quad_vert, const NTYPEQ num_quad,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad)
  {
    const int NUM_VERT_PER_QUAD(4);
    const NTYPEQ list_len = num_quad*NUM_VERT_PER_QUAD;
    num_triangulated_quad = 0;

    for (NTYPEQ k = 0; k < list_len; k += NUM_VERT_PER_QUAD) {
      if (!is_quad_collapsed_to_one_vertex(quad_vert+k)) {
        triangulate_quad_using_diagonal02(quad_vert+k, tri_vert);
        num_triangulated_quad++;
      }
    }
  }
    

  /*!
   *  @brief Triangulate a list of quadrilaterals using diagonal (v0,v2).
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_skip_collapsed_quads If true, skip any quads
   *    that have been collapsed to a single vertex.
   *    - (Still triangulates quads that are collapsed to an edge
   *       or two edges.)
   *    - Default: true.
   */
  template <typename NTYPEQ, typename VTYPE0, typename VTYPE1,
            typename NTYPE_TRI>
  void triangulate_quad_list_using_diagonal02
  (const VTYPE0 * quad_vert, const NTYPEQ num_quad,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad,
   const bool flag_skip_collapsed_quads)
  {
    num_triangulated_quad = 0;
    
    if (flag_skip_collapsed_quads) {
      triangulate_quad_list_using_diagonal02_skip_collapsed
        (quad_vert, num_quad, tri_vert, num_triangulated_quad);
    }
    else {
      triangulate_quad_list_using_diagonal02_skip_collapsed
        (quad_vert, num_quad, tri_vert, num_triangulated_quad);
    }
  }


  /*!
   *  @overload
   *  @brief Triangulate a list of quadrilaterals 
   *    using diagonal (v0,v2). (C++ vector.)
   *  - C++ STL vector format for quad_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename NTYPE_TRI>
  void triangulate_quad_list_using_diagonal02
  (const std::vector<VTYPE0> quad_vert,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad,
   const bool flag_skip_collapsed_quads=true)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_QUAD(4);

    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    triangulate_quad_list_using_diagonal02
      (IJK::vector2pointer(quad_vert), num_quad,
       tri_vert, num_triangulated_quad, flag_skip_collapsed_quads);
  }


  /*!
   *  @brief Triangulate a quadrilateral from an interior vertex w0.
   *  - Add a triangle between each quad edge and vertex w0.
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_quad_from_interior_vertex
  (const VTYPE0 * quad_vert, const VTYPE1 w0, std::vector<VTYPE2> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD = 4;

    triangulate_polygon_from_interior_vertex
      (NUM_VERT_PER_QUAD, quad_vert, w0, tri_vert);
  }


  /*!
   *  @brief Triangulate a list of quadrilaterals from interior vertices.
   *  - For each quad, add a triangle between each quad edge and an interior vertex.
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  - Skip any quads that have been collapsed to a single vertex.
   *    - (Still triangulates quads that are collapsed to an edge
   *       or two edges.)
   *  @param first_vertex Index of first interior vertex.
   *    - Quad i is triangulated using interior vertex (first_vertex+i).
   */
  template <typename NTYPEQ, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2, typename NTYPETRIQ>
  void triangulate_quad_list_from_interior_vertices_skip_collapsed
  (const VTYPE0 * quad_vert, const NTYPEQ num_quad, 
   const VTYPE1 first_vertex, std::vector<VTYPE2> & tri_vert,
   NTYPETRIQ & num_triangulated_quad)
  {
    const int NUM_VERT_PER_QUAD(4);

    num_triangulated_quad = 0;
    
    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      const NTYPEQ k = iquad*NUM_VERT_PER_QUAD;

      if (!is_quad_collapsed_to_one_vertex(quad_vert+k)) {
        triangulate_polygon_from_interior_vertex
          (NUM_VERT_PER_QUAD, quad_vert+k, first_vertex+iquad, tri_vert);
        num_triangulated_quad++;
      }
    }
  }


  /*!
   *  @brief Triangulate a list of quadrilaterals from interior vertices.
   *  - For each quad, add a triangle between each quad edge and an interior vertex.
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  @param first_vertex Index of first interior vertex.
   *    - Quad i is triangulated using interior vertex (first_vertex+i).
   *  @param flag_skip_collapsed_quads If true, skip any quads
   *    that have been collapsed to a single vertex.
   *    - (Still triangulates quads that are collapsed to an edge
   *       or two edges.)
   *    - Default: true.
   */
  template <typename NTYPEQ, typename VTYPE0, typename VTYPE1, 
            typename VTYPE2, typename NTYPETRIQ>
  void triangulate_quad_list_from_interior_vertices
  (const VTYPE0 * quad_vert, const NTYPEQ num_quad, 
   const VTYPE1 first_vertex, std::vector<VTYPE2> & tri_vert,
   NTYPETRIQ & num_triangulated_quad,
   const bool flag_skip_collapsed_quads=true)
  {
    const int NUM_VERT_PER_QUAD(4);

    num_triangulated_quad = 0;
    
    if (flag_skip_collapsed_quads) {
      triangulate_quad_list_from_interior_vertices_skip_collapsed
        (quad_vert, num_quad, first_vertex,
         tri_vert, num_triangulated_quad);
    }
    else {
      for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {

        const NTYPEQ k = iquad*NUM_VERT_PER_QUAD;
        triangulate_polygon_from_interior_vertex
          (NUM_VERT_PER_QUAD, quad_vert+k, first_vertex+iquad,
           tri_vert);

        num_triangulated_quad++;
      }
    }
  }


  /*!
   *  @brief Triangulate a list of quadrilaterals 
   *    from interior vertices. (C++ vector.)
   *  - For each quad, add a triangle between each quad edge 
   *    and an interior vertex.
   *  - C++ STL vector format for quad_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename NTYPETRIQ>
  void triangulate_quad_list_from_interior_vertices
  (const std::vector<VTYPE0> & quad_vert,
   const VTYPE1 first_vertex, std::vector<VTYPE2> & tri_vert,
   NTYPETRIQ & num_triangulated_quad,
   const bool flag_skip_collapsed_quads=true)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_QUAD(4);

    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    triangulate_quad_list_from_interior_vertices
      (IJK::vector2pointer(quad_vert), num_quad, first_vertex,
       tri_vert, num_triangulated_quad);
  }


  /*!
   *  @brief Triangulate a quadrilateral with diagonal (0,2) or (1,3).
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_quad_using_diagonal
  (const VTYPE0 * quad_vert, const bool flag_diag02,
   std::vector<VTYPE1> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);

    if (flag_diag02) {
      triangulate_polygon_from_vertex(NUM_VERT_PER_QUAD, quad_vert, 0, tri_vert);
    }
    else {
      triangulate_polygon_from_vertex(NUM_VERT_PER_QUAD, quad_vert, 1, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a quadrilateral with diagonal (0,2) or (1,3).
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  - Reverse triangle_orientation.
   */
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_quad_using_diagonal_reverse_orient
  (const VTYPE0 * quad_vert, const bool flag_diag02,
   std::vector<VTYPE1> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);

    if (flag_diag02) {
      triangulate_polygon_from_vertex_reverse_orient
        (NUM_VERT_PER_QUAD, quad_vert, 0, tri_vert);
    }
    else {
      triangulate_polygon_from_vertex_reverse_orient
        (NUM_VERT_PER_QUAD, quad_vert, 1, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a quadrilateral with diagonal (0,2) or (1,3).
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_reverse_orient If true, reverse triangle orientation.
   */
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_quad_using_diagonal
  (const VTYPE0 * quad_vert, const bool flag_diag02,
   const bool flag_reverse_orient, std::vector<VTYPE1> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_quad_using_diagonal_reverse_orient
        (quad_vert, flag_diag02, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal(quad_vert, flag_diag02, tri_vert);
    }
  }


  ///@}


  // **************************************************
  //! @name TRIANGULATE PENTAGONS
  // **************************************************

  ///@{

  /*!
   *  @brief Triangulate a pentagon by adding a triangle between each polygon edge and vertex w0.
   *  - Pentagon vertices are listed in clockwise or counter-clockwise order
   *    around the pentagon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_reverse_orient If true, reverse triangle orientations.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_pentagon_with_vertex
  (const VTYPE0 * pentagon_vert, const VTYPE1 w0,
   const bool flag_reverse_orient, std::vector<VTYPE2> & tri_vert)
  {
    const int NUM_VERT_PER_PENTAGON(5);


      (NUM_VERT_PER_PENTAGON, pentagon_vert, w0, 
       flag_reverse_orient, tri_vert);
  }


  /*!
   *  @brief Triangulate a pentagon by adding triangles (v0,v1,v2), (v0, v2, v3) and (v0, v3, v4).
   *  - Pentagon vertices are listed in clockwise or counter-clockwise
   *    order around the pentagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPEB>
  void triangulate_pentagon
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(v0, v1, v2, tri_vert);
    add_triangle_vertices(v0, v2, v3, tri_vert);
    add_triangle_vertices(v0, v3, v4, tri_vert);
  }


  /*!
   *  @brief Triangulate a pentagon by adding triangles (v0,v1,v2), (v0, v2, v3) and (v0, v3, v4).
   *  - Version allowing reverse orientation.
   *  @param flag_reverse_orient If true, reverse orientation.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPEB>
  void triangulate_pentagon
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const bool flag_reverse_orient,
   std::vector<VTYPEB> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_pentagon(v0, v4, v3, v2, v1, tri_vert);
    }
    else {
      triangulate_pentagon(v0, v1, v2, v3, v4, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a pentagon by adding triangles incident on pentagon_vert[ivX].
   *  - Pentagon vertices are listed in clockwise or counter-clockwise
   *    order around the pentagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_pentagon
  (const VTYPE pentagon_vert[], const ITYPE ivX,
   std::vector<VTYPEB> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_PENTAGON(5);
    const ITYPE i1 = (ivX+1)%NUM_VERT_PER_PENTAGON;
    const ITYPE i2 = (ivX+2)%NUM_VERT_PER_PENTAGON;
    const ITYPE i3 = (ivX+3)%NUM_VERT_PER_PENTAGON;
    const ITYPE i4 = (ivX+4)%NUM_VERT_PER_PENTAGON;

    triangulate_pentagon
      (pentagon_vert[ivX], pentagon_vert[i1], pentagon_vert[i2],
       pentagon_vert[i3], pentagon_vert[i4], tri_vert);
  }


  /*!
   *  @brief Triangulate a pentagon by adding triangles incident on pentagon_vert[ivX].
   *  - Pentagon vertices are listed in clockwise or counter-clockwise
   *    order around the pentagon.
   *  - Add new triangles to vector tri_vert.
   *  - Reverse orientation.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_pentagon_reverse_orient
  (const VTYPE pentagon_vert[], const ITYPE ivX,
   std::vector<VTYPEB> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_PENTAGON(5);
    const ITYPE i1 = (ivX+1)%NUM_VERT_PER_PENTAGON;
    const ITYPE i2 = (ivX+2)%NUM_VERT_PER_PENTAGON;
    const ITYPE i3 = (ivX+3)%NUM_VERT_PER_PENTAGON;
    const ITYPE i4 = (ivX+4)%NUM_VERT_PER_PENTAGON;

    triangulate_pentagon
      (pentagon_vert[ivX], pentagon_vert[i4], pentagon_vert[i3],
       pentagon_vert[i2], pentagon_vert[i1], tri_vert);
  }


  /*!
   *  @brief Triangulate a pentagon by adding triangles incident on pentagon_vert[ivX].
   *  - Pentagon vertices are listed in clockwise or counter-clockwise
   *    order around the pentagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_pentagon
  (const VTYPE pentagon_vert[], const ITYPE ivX, 
   const bool flag_reverse_orient, std::vector<VTYPEB> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_pentagon_reverse_orient(pentagon_vert, ivX, tri_vert);
    }
    else {
      triangulate_pentagon(pentagon_vert, ivX, tri_vert);
    }
  }


  /// Triangulate a pentagon.
  /// - Triangles are specified by ears.
  /// @param ear0 Index of first ear in triangulation. In range [0..5].
  /// @param ear1 Index of second ear in triangulation. In range [0..5].
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, 
            typename EAR_TYPE0, typename EAR_TYPE1,
            typename VTYPEB>
  void triangulate_pentagon_by_ears
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   const EAR_TYPE0 ear0, const EAR_TYPE1 ear1,
   std::vector<VTYPEB> & tri_vert)
  {
    switch(ear0) {

    case 0:
      if (ear1 == 1 || ear1 == 3) 
        { triangulate_pentagon(v4, v0, v1, v2, v3, tri_vert); }
      else
        { triangulate_pentagon(v1, v2, v3, v4, v0, tri_vert); }
      break;

    case 1:
      if (ear1 == 2 || ear1 == 4) 
        { triangulate_pentagon(v0, v1, v2, v3, v4, tri_vert); }
      else
        { triangulate_pentagon(v2, v3, v4, v0, v1, tri_vert); }
      break;

    case 2:
      if (ear1 == 0 || ear1 == 3) 
        { triangulate_pentagon(v1, v2, v3, v4, v0, tri_vert); }
      else
        { triangulate_pentagon(v3, v4, v0, v1, v2, tri_vert); }
      break;

    case 3:
      if (ear1 == 0 || ear1 == 2) 
        { triangulate_pentagon(v4, v0, v1, v2, v3, tri_vert); }
      else
        { triangulate_pentagon(v2, v3, v4, v0, v1, tri_vert); }
      break;

    case 4:
    default:
      if (ear1 == 0 || ear1 == 2) 
        { triangulate_pentagon(v3, v4, v0, v1, v2, tri_vert); }
      else
        { triangulate_pentagon(v0, v1, v2, v3, v4, tri_vert); }
      break;
    }
  }


  /*!
   *  @brief Triangulate a pentagon by adding triangles 
   *    (w0,v0,v1), (w0,v1,v2), (w0,v2,v3), (w0,v3,v4) and (w0, v4, v0).
   *  - Pentagon vertices are listed in clockwise or counter-clockwise
   *    order around the pentagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename WTYPE, typename VTYPEB>
  void triangulate_pentagon_with_vertex
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const WTYPE w0,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(w0, v0, v1, tri_vert);
    add_triangle_vertices(w0, v1, v2, tri_vert);
    add_triangle_vertices(w0, v2, v3, tri_vert);
    add_triangle_vertices(w0, v3, v4, tri_vert);
    add_triangle_vertices(w0, v4, v0, tri_vert);
  }


  /*!
   *  @brief Triangulate a pentagon by adding triangles (w0,v1,v0), (w0,v2,v1), (w0,v3,v2), (w0,v4,v3) and (w0, v0, v4).
   *  - Pentagon vertices are listed in clockwise or counter-clockwise
   *    order around the pentagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename WTYPE, typename VTYPEB>
  void triangulate_pentagon_with_vertex_reverse_orient
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const WTYPE w0,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(w0, v1, v0, tri_vert);
    add_triangle_vertices(w0, v2, v1, tri_vert);
    add_triangle_vertices(w0, v3, v2, tri_vert);
    add_triangle_vertices(w0, v4, v3, tri_vert);
    add_triangle_vertices(w0, v0, v4, tri_vert);
  }


  /*!
   *  @brief Triangulate a pentagon by adding a triangle between each pentagon edge and vertex w0.
   *  - Pentagon vertices are listed in clockwise or counter-clockwise
   *    order around the pentagon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_reverse_orient If true, reverse triangle orientations.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename WTYPE, typename VTYPEB>
  void triangulate_pentagon_with_vertex
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const WTYPE w0,
   const bool flag_reverse_orient, std::vector<VTYPEB> & tri_vert)
  {

    if (flag_reverse_orient) {
      triangulate_pentagon_with_vertex_reverse_orient
        (v0, v1, v2, v3, v4, w0, tri_vert);
    }
    else {
      triangulate_pentagon_with_vertex
        (v0, v1, v2, v3, v4, w0, tri_vert);
    }
  }

  ///@}



  // **************************************************
  //! @name TRIANGULATE HEXAGONS
  // **************************************************

  ///@{

  /*!
   *  @brief Triangulate a hexagon by adding triangles (v0,v1,v2), (v0, v2, v3), (v0, v3, v4) and (v0, v4, v5).
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPEB>
  void triangulate_hexagon
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPE4 v4, const VTYPE5 v5, std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(v0, v1, v2, tri_vert);
    add_triangle_vertices(v0, v2, v3, tri_vert);
    add_triangle_vertices(v0, v3, v4, tri_vert);
    add_triangle_vertices(v0, v4, v5, tri_vert);
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangles (v0,v1,v2), (v0, v2, v3), (v0, v3, v4) and (v0, v4, v5).
   *  - Version allowing reverse orientation.
   *  @param flag_reverse_orient If true, reverse orientation.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPEB>
  void triangulate_hexagon
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const bool flag_reverse_orient, std::vector<VTYPEB> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_hexagon(v0, v5, v4, v3, v2, v1, tri_vert);
    }
    else {
      triangulate_hexagon(v0, v1, v2, v3, v4, v5, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangles incident on hexagon_vert[ivX].
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_hexagon
  (const VTYPE hexagon_vert[], const ITYPE ivX,
   std::vector<VTYPEB> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_HEXAGON(6);
    const ITYPE i1 = (ivX+1)%NUM_VERT_PER_HEXAGON;
    const ITYPE i2 = (ivX+2)%NUM_VERT_PER_HEXAGON;
    const ITYPE i3 = (ivX+3)%NUM_VERT_PER_HEXAGON;
    const ITYPE i4 = (ivX+4)%NUM_VERT_PER_HEXAGON;
    const ITYPE i5 = (ivX+5)%NUM_VERT_PER_HEXAGON;

    triangulate_hexagon
      (hexagon_vert[ivX], hexagon_vert[i1], hexagon_vert[i2],
       hexagon_vert[i3], hexagon_vert[i4], hexagon_vert[i5], tri_vert);
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangles incident on hexagon_vert[ivX].
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   *  - Reverse orientation.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_hexagon_reverse_orient
  (const VTYPE hexagon_vert[], const ITYPE ivX,
   std::vector<VTYPEB> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_HEXAGON(6);
    const ITYPE i1 = (ivX+1)%NUM_VERT_PER_HEXAGON;
    const ITYPE i2 = (ivX+2)%NUM_VERT_PER_HEXAGON;
    const ITYPE i3 = (ivX+3)%NUM_VERT_PER_HEXAGON;
    const ITYPE i4 = (ivX+4)%NUM_VERT_PER_HEXAGON;
    const ITYPE i5 = (ivX+5)%NUM_VERT_PER_HEXAGON;

    triangulate_hexagon
      (hexagon_vert[ivX], hexagon_vert[i5], hexagon_vert[i4], 
       hexagon_vert[i3], hexagon_vert[i2], hexagon_vert[i1], tri_vert);
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangles incident on hexagon_vert[ivX].
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_hexagon
  (const VTYPE hexagon_vert[], const ITYPE ivX, 
   const bool flag_reverse_orient, std::vector<VTYPEB> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_hexagon_reverse_orient(hexagon_vert, ivX, tri_vert);
    }
    else {
      triangulate_hexagon(hexagon_vert, ivX, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangle (iv0,iv2,iv4).
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPEB>
  void triangulate_hexagon_using_triangle_024
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(v0, v1, v2, tri_vert);
    add_triangle_vertices(v2, v3, v4, tri_vert);
    add_triangle_vertices(v4, v5, v0, tri_vert);
    add_triangle_vertices(v0, v2, v4, tri_vert);
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangle (iv4,iv2,iv0).
   *  - Triangles have reverse orientation of hexagon.
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPEB>
  void triangulate_hexagon_using_triangle_024_reverse_orient
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(v2, v1, v0, tri_vert);
    add_triangle_vertices(v4, v3, v2, tri_vert);
    add_triangle_vertices(v0, v5, v4, tri_vert);
    add_triangle_vertices(v4, v2, v0, tri_vert);
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangle (iv0,iv2,iv4).
   *  - Version with flag_reverse_orient.
   *  @param flag_reverse_orient
   *    If true, reverse orientation in creating triangles from the hexagon.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPEB>
  void triangulate_hexagon_using_triangle_024
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const bool flag_reverse_orient,
   std::vector<VTYPEB> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_hexagon_using_triangle_024_reverse_orient
        (v0, v1, v2, v3, v4, v5, tri_vert);
    }
    else {
      triangulate_hexagon_using_triangle_024
        (v0, v1, v2, v3, v4, v5, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a hexagon.
   *  - Triangles are specified by ears.
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   *  @param ear0 Index of first ear in triangulation. In range [0..5].
   *  @param ear1 Index of second ear in triangulation. In range [0..5].
   *  @param ear2 Index of third ear in triangulation. In range [0..5].
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename EAR_TYPE0, typename EAR_TYPE1, typename EAR_TYPE2,
            typename VTYPEB>
  void triangulate_hexagon_by_ears
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const EAR_TYPE0 ear0, const EAR_TYPE1 ear1, const EAR_TYPE2 ear2,
   std::vector<VTYPEB> & tri_vert)
  {
    EAR_TYPE1 ear1B;
    EAR_TYPE2 ear2B;

    if (ear1 > ear0) { ear1B = ear1-1; }
    else { ear1B = ear1; }
    if (ear2 > ear0) { ear2B = ear2-1; }
    else { ear2B = ear2; }

    switch(ear0) {

    case 0:
      add_triangle_vertices(v5, v0, v1, tri_vert);
      triangulate_pentagon_by_ears
        (v1, v2, v3, v4, v5, ear1B, ear2B, tri_vert);
      break;

    case 1:
      add_triangle_vertices(v0, v1, v2, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v2, v3, v4, v5, ear1B, ear2B, tri_vert);
      break;

    case 2:
      add_triangle_vertices(v1, v2, v3, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v1, v3, v4, v5, ear1B, ear2B, tri_vert);
      break;

    case 3:
      add_triangle_vertices(v2, v3, v4, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v1, v2, v4, v5, ear1B, ear2B, tri_vert);
      break;

    case 4:
      add_triangle_vertices(v3, v4, v5, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v1, v2, v3, v5, ear1B, ear2B, tri_vert);
      break;

    case 5:
    default:
      add_triangle_vertices(v4, v5, v0, tri_vert);
      triangulate_pentagon_by_ears
        (v0, v1, v2, v3, v4, ear1B, ear2B, tri_vert);
      break;
    }
  }


  /*!
   *  @brief Triangulate a hexagon.
   *  - Version with array of hexagon vertices.
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   */
  template <typename VTYPE,
            typename EAR_TYPE0, typename EAR_TYPE1, typename EAR_TYPE2,
            typename VTYPEB>
  void triangulate_hexagon_by_ears
  (const VTYPE hex_vert[],
   const EAR_TYPE0 ear0, const EAR_TYPE1 ear1, const EAR_TYPE2 ear2,
   std::vector<VTYPEB> & tri_vert)
  {
    triangulate_hexagon_by_ears
      (hex_vert[0], hex_vert[1], hex_vert[2], hex_vert[3], hex_vert[4],
       hex_vert[5], ear0, ear1, ear2, tri_vert);
  }


  /*!
   *  Triangulate a hexagon by cutting one ear and triangulating the remaining pentagon.
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ITYPE>
  void triangulate_hexagon_by_pentagon_tri5_and_triangle
  (const VTYPE0 hex_vert[], const ITYPE iv_ear, const VTYPE1 ivertX,
   const bool flag_reverse_orient, std::vector<VTYPE2> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_PENTAGON(5);
    const ITYPE NUM_VERT_PER_HEX(6);
    const ITYPE i0 = (iv_ear+NUM_VERT_PER_HEX-1)%NUM_VERT_PER_HEX;
    const ITYPE i1 = iv_ear;
    const ITYPE i2 = (iv_ear+1)%NUM_VERT_PER_HEX;
    const VTYPE0 pentagon_vert[NUM_VERT_PER_PENTAGON];

    add_triangle_vertices
      (hex_vert[i0], hex_vert[i1], hex_vert[i2], 
       flag_reverse_orient, tri_vert);

    ITYPE j = 0;
    for (ITYPE i = 0; i < NUM_VERT_PER_HEX; i++) {
      if (i != iv_ear) {
        pentagon_vert[j] = hex_vert[i];
        j++;
      }
    }

    if (j != NUM_VERT_PER_PENTAGON) {
      IJK::PROCEDURE_ERROR error
        ("triangulate_hexagon_by_pentagon_tri5_and_triangle");
      error.AddMessage
        ("Programming error.  Incorrect number of pentagon vertices.");
      error.AddMessage("  Added ", j, " vertices to pentagon.");
      throw error;
    }

    triangulate_pentagon_with_vertex
      (pentagon_vert, ivertX, flag_reverse_orient, tri_vert);
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangles containing w0.
   *  - Add triangulates (w0,v0,v1), (w0,v1,v2), (w0,v2,v3), (w0,v3,v4),
   *    (w0, v4, v5) and (w0, v5, v0).
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename WTYPE, typename VTYPEB>
  void triangulate_hexagon_with_vertex
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const WTYPE w0, std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(w0, v0, v1, tri_vert);
    add_triangle_vertices(w0, v1, v2, tri_vert);
    add_triangle_vertices(w0, v2, v3, tri_vert);
    add_triangle_vertices(w0, v3, v4, tri_vert);
    add_triangle_vertices(w0, v4, v5, tri_vert);
    add_triangle_vertices(w0, v5, v0, tri_vert);
  }


  /*!
   *  @brief Triangulate a hexagon by adding triangles containing w0.
   *  - Add triangles (w0,v1,v0), (w0,v2,v1), (w0,v3,v2), (w0,v4,v3),
   *    (w0, v3, v2) and (w0, v0, v5).
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename WTYPE, typename VTYPEB>
  void triangulate_hexagon_with_vertex_reverse_orient
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const WTYPE w0, std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(w0, v1, v0, tri_vert);
    add_triangle_vertices(w0, v2, v1, tri_vert);
    add_triangle_vertices(w0, v3, v2, tri_vert);
    add_triangle_vertices(w0, v4, v3, tri_vert);
    add_triangle_vertices(w0, v5, v4, tri_vert);
    add_triangle_vertices(w0, v0, v5, tri_vert);
  }


  /*!
   *  @brief Triangulate a hexagon by adding a triangle between each hexagon edge and vertex w0.
   *  - Hexagon vertices are listed in clockwise or counter-clockwise
   *    order around the hexagon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_reverse_orient If true, reverse triangle orientations.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename WTYPE, typename VTYPEB>
  void triangulate_hexagon_with_vertex
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const WTYPE w0, const bool flag_reverse_orient, 
   std::vector<VTYPEB> & tri_vert)
  {

    if (flag_reverse_orient) {
      triangulate_hexagon_with_vertex_reverse_orient
        (v0, v1, v2, v3, v4, v5, w0, tri_vert);
    }
    else {
      triangulate_hexagon_with_vertex
        (v0, v1, v2, v3, v4, v5, w0, tri_vert);
    }
  }

  ///@}


  // **************************************************
  //! @name TRIANGULATE SEPTAGONS
  // **************************************************

  ///@{

  /*!
   *  @brief Triangulate a septagon by adding triangles containing v0.
   *  - Add triangles (v0,v1,v2), (v0, v2, v3), (v0, v3, v4), (v0, v4, v5)
   *     and (v0, v5, v6).
   *  - Septagon vertices are listed in clockwise or counter-clockwise
   *    order around the septagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPE6, typename VTYPEB>
  void triangulate_septagon
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPE4 v4, const VTYPE5 v5, const VTYPE6 v6,
   std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(v0, v1, v2, tri_vert);
    add_triangle_vertices(v0, v2, v3, tri_vert);
    add_triangle_vertices(v0, v3, v4, tri_vert);
    add_triangle_vertices(v0, v4, v5, tri_vert);
    add_triangle_vertices(v0, v5, v6, tri_vert);
  }


  /*!
   *  @brief Triangulate a septagon by adding triangles incident on septagon_vert[ivX].
   *  - Septagon vertices are listed in clockwise or counter-clockwise
   *    order around the pentagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_septagon
  (const VTYPE septagon_vert[], const ITYPE ivX,
   std::vector<VTYPEB> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_SEPTAGON(7);
    const ITYPE i1 = (ivX+1)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i2 = (ivX+2)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i3 = (ivX+3)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i4 = (ivX+4)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i5 = (ivX+5)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i6 = (ivX+6)%NUM_VERT_PER_SEPTAGON;

    triangulate_septagon
      (septagon_vert[ivX], septagon_vert[i1], septagon_vert[i2],
       septagon_vert[i3], septagon_vert[i4], septagon_vert[i5],
       septagon_vert[i6], tri_vert);
  }


  /*!
   *  @brief Triangulate a septagon by adding triangles incident on septagon_vert[ivX].
   *  - Septagon vertices are listed in clockwise or counter-clockwise
   *    order around the septagon.
   *  - Add new triangles to vector tri_vert.
   *  - Reverse orientation.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_septagon_reverse_orient
  (const VTYPE septagon_vert[], const ITYPE ivX,
   std::vector<VTYPEB> & tri_vert)
  {
    const ITYPE NUM_VERT_PER_SEPTAGON(7);
    const ITYPE i1 = (ivX+1)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i2 = (ivX+2)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i3 = (ivX+3)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i4 = (ivX+4)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i5 = (ivX+5)%NUM_VERT_PER_SEPTAGON;
    const ITYPE i6 = (ivX+6)%NUM_VERT_PER_SEPTAGON;

    triangulate_septagon
      (septagon_vert[ivX], septagon_vert[i6], septagon_vert[i5],
       septagon_vert[i4], septagon_vert[i3], septagon_vert[i2], 
       septagon_vert[i1], tri_vert);
  }


  /*!
   *  @brief Triangulate a septagon by adding triangles incident on septagon_vert[ivX].
   *  - Septagon vertices are listed in clockwise or counter-clockwise
   *    order around the septagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE, typename ITYPE, typename VTYPEB>
  void triangulate_septagon
  (const VTYPE septagon_vert[], const ITYPE ivX, 
   const bool flag_reverse_orient, std::vector<VTYPEB> & tri_vert)
  {
    if (flag_reverse_orient) {
      triangulate_septagon_reverse_orient(septagon_vert, ivX, tri_vert);
    }
    else {
      triangulate_septagon(septagon_vert, ivX, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a septagon by adding triangles containing w0.
   *  - Add triangles (w0,v0,v1), (w0,v1,v2), (w0,v2,v3), (w0,v3,v4),
   *    (w0, v4, v5), (w0, v5, v6) and (w0, v6, v0).
   *  - Septagon vertices are listed in clockwise or counter-clockwise
   *    order around the septagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPE6, typename WTYPE, typename VTYPEB>
  void triangulate_septagon_with_vertex
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const VTYPE6 v6, const WTYPE w0, std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(w0, v0, v1, tri_vert);
    add_triangle_vertices(w0, v1, v2, tri_vert);
    add_triangle_vertices(w0, v2, v3, tri_vert);
    add_triangle_vertices(w0, v3, v4, tri_vert);
    add_triangle_vertices(w0, v4, v5, tri_vert);
    add_triangle_vertices(w0, v5, v6, tri_vert);
    add_triangle_vertices(w0, v6, v0, tri_vert);
  }


  /*!
   *  @brief Triangulate a septagon by adding triangles containing w0.
   *  - Add triangles (w0,v0,v1), (w0,v1,v2), (w0,v2,v3), (w0,v3,v4),
   *    (w0, v4, v5), (w0, v5, v6) and (w0, v6, v0).
   *  - Septagon vertices are listed in clockwise or counter-clockwise
   *    order around the septagon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPE6, typename WTYPE, typename VTYPEB>
  void triangulate_septagon_with_vertex_reverse_orient
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const VTYPE6 v6, const WTYPE w0, std::vector<VTYPEB> & tri_vert)
  {
    add_triangle_vertices(w0, v1, v0, tri_vert);
    add_triangle_vertices(w0, v2, v1, tri_vert);
    add_triangle_vertices(w0, v3, v2, tri_vert);
    add_triangle_vertices(w0, v4, v3, tri_vert);
    add_triangle_vertices(w0, v5, v4, tri_vert);
    add_triangle_vertices(w0, v6, v5, tri_vert);
    add_triangle_vertices(w0, v0, v6, tri_vert);
  }


  /*!
   *  @brief Triangulate a septagon by adding triangles containing w0.
   *  - Add a triangle between each septagon edge and vertex w0.
   *  - Septagon vertices are listed in clockwise or counter-clockwise
   *    order around the septagon.
   *  - Add new triangles to vector tri_vert.
   *  @param flag_reverse_orient If true, reverse triangle orientations.
   */
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPE5,
            typename VTYPE6, typename WTYPE, typename VTYPEB>
  void triangulate_septagon_with_vertex
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const VTYPE6 v6, const WTYPE w0, const bool flag_reverse_orient, 
   std::vector<VTYPEB> & tri_vert)
  {

    if (flag_reverse_orient) {
      triangulate_septagon_with_vertex_reverse_orient
        (v0, v1, v2, v3, v4, v5, v6, w0, tri_vert);
    }
    else {
      triangulate_septagon_with_vertex
        (v0, v1, v2, v3, v4, v5, v6, w0, tri_vert);
    }
  }

  ///@}


  // ***************************************************************
  //! @name TRI4 QUADS WITH DIAGONALS OUTSIDE THE ENVELOPE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Triangulate isosurface quadrilateral if it has 
   *    some diagonal outside the envelope.
   *  - Return true if isosurface quadrilateral is deleted.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If some diagonal is outside the envelope fromed by the quad edges
   *    and the dual grid edge, split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @param quad_vert[] Vertices of the quadrilateral.
   *  @param quad_info Quad information.
   *  @param ivert_in_quad_interior Index of vertex in quad interior.
   *    - Integer type.
   *    - Split quad into four triangles, with each triangle
   *      incident on vertex ivert_in_quad_interior.
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPEQ, typename VTYPEI, typename VTYPET,
            typename ISOQUAD_INFO_TYPE>
  bool tri4_dual_quad_DoutsideE
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   const VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE & quad_info,
   const VTYPEI ivert_in_quad_interior,
   std::vector<VTYPET> & tri_vert)
  {
    if (!quad_info.IsDiagonalInEnvelope(0) ||
        !quad_info.IsDiagonalInEnvelope(1)) {

      // Triangulate into four triangles using vertex ivert_in_quad_interior.
      triangulate_quad_from_interior_vertex
        (quad_vert, ivert_in_quad_interior, tri_vert);

      return true;
    }
    else {
      return false;
    }
  }


  /*!
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If some diagonal is outside the envelope fromed by the quad edges
   *    and the dual grid edge, split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   *  @param index_of_first_interior_vertex Index of vertex in interior of quad 0.
   *    - Integer type.
   *    - Index of vertex in interior of quad i is (index_of_first_interior_vertex+i).
   *    @pre Interior vertices have already been created and
   *    assigned coordinates in vertex_coord[].
   *  @param min_dist_allow_tri4 Allow triangulation into 4 triangles when each
   *    quad vertex is min_distance_allow_tri4 distance from one
   *    of the facets incident on the dual quad edge.
   *    - If some quad vertex is closer than min_distance_allow_tri4 to
   *    both neighboring facets incident on the dual quad edge,
   *    use a diagonal triangulation into two triangles.
   *  @param[out] quad_vert[] Isosurface quadrilateral vertices.
   *    - Quadrilateral i has vertices:
   *      (quad_vert[4*i], quad_vert[4*i+1], quad_vert[4*i+2], quad_vert[4*i+3]).
   *    - If quadrilateral i is deleted, then it's vertices are removed
   *      from quad_vert[].
   *  @param[out] quad_info[i] Information on i'th isosurface quadrilateral.
   *    - If quadrilateral i is deleted, then quad_info[i] is removed
   *      from quad_info[].
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPEQ, typename VTYPEI, typename VTYPET,
            typename ISOQUAD_INFO_TYPE, typename NTYPES>
  void tri4_dual_quad_list_with_diagonals_outside_envelope
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE quad_info[],
   const NTYPEQ num_quad,
   const VTYPEI index_of_first_interior_vertex,
   std::vector<VTYPET> & tri_vert,
   NTYPES & num_split)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const int NUM_VERT_PER_QUAD(4);
    const DTYPE dimension = grid.Dimension();
    const DTYPE DIM3(3);
    IJK::PROCEDURE_ERROR error
      ("tri4_dual_quad_list_with_diagonals_outside_envelope");

    // Initialize
    num_split = 0;

    if (num_quad > 0 && quad_info == NULL) {
      error.AddMessage
        ("Programming error.  Array quad_info[] is empty.");
      throw error;
    }

    if (dimension != DIM3) {
      error.AddMessage("Programming error. Vertices are not in dimension three.");
      error.AddMessage("Vertices in dimension: ", dimension);
      throw error;
    }

    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      const VTYPEI ivert_in_quad_interior =
        index_of_first_interior_vertex+iquad;
      VTYPEQ * quad_vert_iquad =
        quad_vert + iquad*NUM_VERT_PER_QUAD;

      const bool is_deleted =
        tri4_dual_quad_DoutsideE
        (grid, vertex_coord, quad_vert_iquad, quad_info[iquad],
         ivert_in_quad_interior, tri_vert);

      if (is_deleted) {
        collapse_quad_to_vertex(quad_vert_iquad);
        num_split++;
      }
    }
  }


  /*!
   *  @overload
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope.
   *  - Always triangulate into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - Version using C++ STL vectors vertex_coord[] and quad_vert[].
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ISOQUAD_INFO_TYPE, typename NTYPES>
  void tri4_dual_quad_list_with_diagonals_outside_envelope
  (const GRID_TYPE & grid,
   const std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE0> & quad_vert,
   const std::vector<ISOQUAD_INFO_TYPE> & quad_info,
   const VTYPE1 index_of_first_interior_vertex,
   std::vector<VTYPE2> & tri_vert,
   NTYPES & num_split)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    
    if (num_quad != quad_info.size()) {
      IJK::PROCEDURE_ERROR error
        ("tri4_dual_quad_list_with_diagonals_outside_envelope");
      error.AddMessage("Programming error. Incorrect size of quad_info[].");
      error.AddMessage
        ("Size of quad_info[] does not equal number of isosurface quadrilaterals.");
      error.AddMessage("  Size of quad info[]: ", quad_info.size(), ".");
      error.AddMessage("  Number of isosurface quadrilaterals: ",
                       num_quad, ".");
      throw error;
    }

    tri4_dual_quad_list_with_diagonals_outside_envelope
      (grid, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert),
       IJK::vector2pointer(quad_info), num_quad,
       index_of_first_interior_vertex,
       tri_vert, num_split);
  }


  /*!
   *  @brief Triangulate isosurface quadrilateral if it has 
   *    some diagonal outside the envelope, adding vertex.
   *  - Return true if isosurface quadrilateral is deleted.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If some diagonal is outside the envelope fromed by the quad edges
   *    and the dual grid edge, split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @param quad_vert[] Vertices of the quadrilateral.
   *  @param quad_info Quad information.
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPEQ, typename VTYPET,
            typename ISOQUAD_INFO_TYPE,
            typename STYPE, typename ADD_ISOV_FUNCTION>
  bool tri4_dual_quad_DoutsideE_add_vertex
  (const GRID_TYPE & grid,
   const VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE & quad_info,
   const STYPE isovalue,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPET> & tri_vert,
   const ADD_ISOV_FUNCTION & add_isov_function)
  {
    if (!quad_info.IsDiagonalInEnvelope(0) ||
        !quad_info.IsDiagonalInEnvelope(1)) {

      const VTYPEQ ivert_in_quad_interior =
        add_isov_function.AddDualIsov
        (grid.Dimension(), quad_vert, quad_info, isovalue,
         vertex_coord);
         
      // Triangulate into four triangles using vertex ivert_in_quad_interior.
      triangulate_quad_from_interior_vertex
        (quad_vert, ivert_in_quad_interior, tri_vert);

      return true;
    }
    else {
      return false;
    }
  }

  
  /*!
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope, adding vertices, as necessary.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If some diagonal is outside the envelope fromed by the quad edges
   *    and the dual grid edge, split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   *  @param min_dist_allow_tri4 Allow triangulation into 4 triangles when each
   *    quad vertex is min_distance_allow_tri4 distance from one
   *    of the facets incident on the dual quad edge.
   *    - If some quad vertex is closer than min_distance_allow_tri4 to
   *    both neighboring facets incident on the dual quad edge,
   *    use a diagonal triangulation into two triangles.
   *  @param[out] quad_vert[] Isosurface quadrilateral vertices.
   *    - Quadrilateral i has vertices:
   *      (quad_vert[4*i], quad_vert[4*i+1], quad_vert[4*i+2], quad_vert[4*i+3]).
   *    - If quadrilateral i is deleted, then it's vertices are removed
   *      from quad_vert[].
   *  @param[out] quad_info[i] Information on i'th isosurface quadrilateral.
   *    - If quadrilateral i is deleted, then quad_info[i] is removed
   *      from quad_info[].
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPEQ, typename VTYPET,
            typename ISOQUAD_INFO_TYPE,
            typename STYPE, typename NTYPES,
            typename ADD_ISOV_FUNCTION>
  void tri4_dual_quad_list_with_diagonals_outside_envelope_add_vertices
  (const GRID_TYPE & grid,
   VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE quad_info[],
   const NTYPEQ num_quad,
   const STYPE isovalue,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPET> & tri_vert,
   NTYPES & num_split,
   const ADD_ISOV_FUNCTION & add_isov_function)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const int NUM_VERT_PER_QUAD(4);
    const DTYPE dimension = grid.Dimension();
    const DTYPE DIM3(3);
    IJK::PROCEDURE_ERROR error
      ("tri4_dual_quad_list_with_diagonals_outside_envelope_add_vertices");

    // Initialize
    num_split = 0;

    if (num_quad > 0 && quad_info == NULL) {
      error.AddMessage
        ("Programming error.  Array quad_info[] is empty.");
      throw error;
    }

    if (dimension != DIM3) {
      error.AddMessage("Programming error. Vertices are not in dimension three.");
      error.AddMessage("Vertices in dimension: ", dimension);
      throw error;
    }

    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      VTYPEQ * quad_vert_iquad =
        quad_vert + iquad*NUM_VERT_PER_QUAD;

      const bool is_deleted =
        tri4_dual_quad_DoutsideE_add_vertex
        (grid, quad_vert_iquad, quad_info[iquad],
         isovalue, vertex_coord, tri_vert, add_isov_function);

      if (is_deleted) {
        collapse_quad_to_vertex(quad_vert_iquad);
        num_split++;
      }
    }
  }

    
  /*!
   *  @overload
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope, adding vertices, as necessary.
   *  - Always triangulate into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  - Version using C++ STL vectors vertex_coord[] and quad_vert[].
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPEQ, typename VTYPET,
            typename ISOQUAD_INFO_TYPE,
            typename STYPE, typename NTYPES,
            typename ADD_ISOV_FUNCTION>
  void tri4_dual_quad_list_with_diagonals_outside_envelope_add_vertices
  (const GRID_TYPE & grid,
   std::vector<VTYPEQ> & quad_vert,
   const std::vector<ISOQUAD_INFO_TYPE> & quad_info,
   const STYPE isovalue,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPET> & tri_vert,
   NTYPES & num_split,
   const ADD_ISOV_FUNCTION & add_isov_function)
  {
    typedef typename std::vector<VTYPEQ>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (num_quad != quad_info.size()) {
      IJK::PROCEDURE_ERROR error
        ("tri2_dual_quad_list_with_diagonals_outside_envelope_add_vertices");
      error.AddMessage("Programming error. Incorrect size of quad_info[].");
      error.AddMessage
        ("Size of quad_info[] does not equal number of isosurface quadrilaterals.");
      error.AddMessage("  Size of quad info: ", quad_info.size(), ".");
      error.AddMessage
        ("  Number of isosurface quadrilaterals: ", num_quad, ".");
      throw error;
    }

    tri4_dual_quad_list_with_diagonals_outside_envelope_add_vertices
      (grid, IJK::vector2pointer(quad_vert), 
       IJK::vector2pointer(quad_info), num_quad,
       isovalue, vertex_coord, tri_vert, num_split,
       add_isov_function);
  }

  //@}


  // ***************************************************************
  //! @name TRI2 OR TRI4 QUADS WITH DIAGONALS OUTSIDE THE ENVELOPE
  // ***************************************************************

  /*!
   *  @brief Triangulate isosurface quadrilateral if it has 
   *    some diagonal outside the envelope, prefer tri2.
   *  - Return number of triangles (0, 2 or 4) in triangulation 
   *    of isosurface quadrilateral. (0: no triangulation.)
   *  - Each quadrilateral is dual to a grid edge.
   *  - If exactly one diagonal is outside the envelope 
   *    formed by the quad edges and the dual grid edge, 
   *    split the quad into two triangles using the other diagonal.
   *  - If both diagonals are outside the envelope,
   *    split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @param quad_vert[] Vertices of the quadrilateral.
   *  @param quad_info Quad information.
   *  @param ivert_in_quad_interior Index of vertex in quad interior.
   *    - Integer type.
   *    - Split quad into four triangles, with each triangle
   *      incident on vertex ivert_in_quad_interior.
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPEQ, typename VTYPEI, typename VTYPET,
            typename ISOQUAD_INFO_TYPE>
  int tri2_or_tri4_dual_quad_DoutsideE_prefer_tri2
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   const VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE & quad_info,
   const VTYPEI ivert_in_quad_interior,
   std::vector<VTYPET> & tri_vert)
  {
    const int ZERO_TRIANGLES(0);
    const int TWO_TRIANGLES(2);
    const int FOUR_TRIANGLES(4);

    if (quad_info.IsDiagonalInEnvelope(0)) {
      if (quad_info.IsDiagonalInEnvelope(1)) {
        return ZERO_TRIANGLES;
      }
      else {
        // Only diagonal (v0,v2) is in the envelope.
        triangulate_quad_from_vertex(quad_vert, 0, tri_vert);
        return TWO_TRIANGLES;
      }
    }
    else if (quad_info.IsDiagonalInEnvelope(1)) {
      // Only diagonal (v1,v3) is in the envelope.
      triangulate_quad_from_vertex(quad_vert, 1, tri_vert);
      return TWO_TRIANGLES;
    }
    else {
      // Neither diagonal is in the envelope.
      // Triangulate into four triangles using vertex ivert_in_quad_interior.
      triangulate_quad_from_interior_vertex
        (quad_vert, ivert_in_quad_interior, tri_vert);

      return FOUR_TRIANGLES;
    }
  }


  /*!
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope, prefer tri2.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If exactly one diagonal is outside the envelope 
   *    formed by the quad edges and the dual grid edge, 
   *    split the quad into two triangles using the other diagonal.
   *  - If both diagonals are outside the envelope,
   *    split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   *  @param index_of_first_interior_vertex Index of vertex in interior of quad 0.
   *    - Integer type.
   *    - Index of vertex in interior of quad i is (index_of_first_interior_vertex+i).
   *    @pre Interior vertices have already been created and
   *    assigned coordinates in vertex_coord[].
   *  @param min_dist_allow_tri4 Allow triangulation into 4 triangles when each
   *    quad vertex is min_distance_allow_tri4 distance from one
   *    of the facets incident on the dual quad edge.
   *    - If some quad vertex is closer than min_distance_allow_tri4 to
   *    both neighboring facets incident on the dual quad edge,
   *    use a diagonal triangulation into two triangles.
   *  @param[out] quad_vert[] Isosurface quadrilateral vertices.
   *    - Quadrilateral i has vertices:
   *      (quad_vert[4*i], quad_vert[4*i+1], quad_vert[4*i+2], quad_vert[4*i+3]).
   *    - If quadrilateral i is deleted, then it's vertices are removed
   *      from quad_vert[].
   *  @param[out] quad_info[i] Information on i'th isosurface quadrilateral.
   *    - If quadrilateral i is deleted, then quad_info[i] is removed
   *      from quad_info[].
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPEQ, typename VTYPEI, typename VTYPET,
            typename ISOQUAD_INFO_TYPE,
            typename NTYPES2, typename NTYPES4>
  void tri2_or_tri4_dual_quad_list_DoutsideE_prefer_tri2
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE quad_info[],
   const NTYPEQ num_quad,
   const VTYPEI index_of_first_interior_vertex,
   std::vector<VTYPET> & tri_vert,
   NTYPES2 & num_tri2_split,
   NTYPES4 & num_tri4_split)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUAD(4);
    const int TWO_TRIANGLES(2);
    const int dimension = grid.Dimension();
    IJK::PROCEDURE_ERROR error
      ("tri4_dual_quad_list_DoutsideE_prefer_tri2");

    // Initialize
    num_tri2_split = 0;
    num_tri4_split = 0;

    if (num_quad > 0 && quad_info == NULL) {
      error.AddMessage
        ("Programming error.  Array quad_info[] is empty.");
      throw error;
    }

    if (dimension != DIM3) {
      error.AddMessage("Programming error. Vertices are not in dimension three.");
      error.AddMessage("Vertices in dimension: ", dimension);
      throw error;
    }

    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      const VTYPEI ivert_in_quad_interior =
        index_of_first_interior_vertex+iquad;
      VTYPEQ * quad_vert_iquad =
        quad_vert + iquad*NUM_VERT_PER_QUAD;

      const int num_tri =
        tri2_or_tri4_dual_quad_DoutsideE_prefer_tri2
        (grid, vertex_coord, quad_vert_iquad, quad_info[iquad],
         ivert_in_quad_interior, tri_vert);

      if (num_tri > 0) {
        collapse_quad_to_vertex(quad_vert_iquad);

        if (num_tri == TWO_TRIANGLES)
          { num_tri2_split++; }
        else
          { num_tri4_split++; }
      }
    }
  }


  /*!
   *  @overload
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope.
   *  - Triangulate into two triangles, is possible.
   *  - Triangulate into four triangles if both diagonals 
   *    are outside the envelope.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - Version using C++ STL vectors vertex_coord[] and quad_vert[].
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ISOQUAD_INFO_TYPE,
            typename NTYPES2, typename NTYPES4>
  void tri2_or_tri4_dual_quad_list_DoutsideE_prefer_tri2
  (const GRID_TYPE & grid,
   const std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE0> & quad_vert,
   const std::vector<ISOQUAD_INFO_TYPE> & quad_info,
   const VTYPE1 index_of_first_interior_vertex,
   std::vector<VTYPE2> & tri_vert,
   NTYPES2 & num_tri2_split,
   NTYPES4 & num_tri4_split)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (num_quad != quad_info.size()) {
      IJK::PROCEDURE_ERROR error
        ("tri4_dual_quad_list_DoutsideE_prefer_tri2");
      error.AddMessage("Programming error. Incorrect size of quad_info[].");
      error.AddMessage
        ("Size of quad_info[] does not equal number of isosurface quadrilaterals.");
      error.AddMessage("  Size of quad info[]: ", quad_info.size(), ".");
      error.AddMessage("  Number of isosurface quadrilaterals: ",
                       num_quad, ".");
      throw error;
    }

    tri2_or_tri4_dual_quad_list_DoutsideE_prefer_tri2
      (grid, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointerNC(quad_vert),
       IJK::vector2pointer(quad_info), num_quad,
       index_of_first_interior_vertex,
       tri_vert, num_tri2_split, num_tri4_split);
  }


  /*!
   *  @brief Triangulate isosurface quadrilateral if it has 
   *    some diagonal outside the envelope, prefer tri2.
   *    Add vertices, as necessary.
   *  - Return number of triangles (0, 2 or 4) in triangulation 
   *    of isosurface quadrilateral. (0: no triangulation.)
   *  - Each quadrilateral is dual to a grid edge.
   *  - If exactly one diagonal is outside the envelope 
   *    formed by the quad edges and the dual grid edge, 
   *    split the quad into two triangles using the other diagonal.
   *  - If both diagonals are outside the envelope,
   *    split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @param quad_vert[] Vertices of the quadrilateral.
   *  @param quad_info Quad information.
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPEQ, typename VTYPET,
            typename ISOQUAD_INFO_TYPE, typename STYPE,
            typename ADD_ISOV_FUNCTION>
  int tri2_or_tri4_dual_quad_DoutsideE_prefer_tri2_addV
  (const GRID_TYPE & grid,
   const VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE & quad_info,
   const STYPE isovalue,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPET> & tri_vert,
   const ADD_ISOV_FUNCTION & add_isov_function)
  {
    const int ZERO_TRIANGLES(0);
    const int TWO_TRIANGLES(2);
    const int FOUR_TRIANGLES(4);

    if (quad_info.IsDiagonalInEnvelope(0)) {
      if (quad_info.IsDiagonalInEnvelope(1)) {
        return ZERO_TRIANGLES;
      }
      else {
        // Only diagonal (v0,v2) is in the envelope.
        triangulate_quad_from_vertex(quad_vert, 0, tri_vert);
        return TWO_TRIANGLES;
      }
    }
    else if (quad_info.IsDiagonalInEnvelope(1)) {
      // Only diagonal (v1,v3) is in the envelope.
      triangulate_quad_from_vertex(quad_vert, 1, tri_vert);
      return TWO_TRIANGLES;
    }
    else {
      const VTYPEQ ivert_in_quad_interior =
        add_isov_function.AddDualIsov
        (grid.Dimension(), quad_vert, quad_info, isovalue,
         vertex_coord);
      
      // Neither diagonal is in the envelope.
      // Triangulate into four triangles using vertex ivert_in_quad_interior.
      triangulate_quad_from_interior_vertex
        (quad_vert, ivert_in_quad_interior, tri_vert);

      return FOUR_TRIANGLES;
    }
  }

  
  /*!
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope, prefer tri2. Add vertices, as necessary.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If exactly one diagonal is outside the envelope 
   *    formed by the quad edges and the dual grid edge, 
   *    split the quad into two triangles using the other diagonal.
   *  - If both diagonals are outside the envelope,
   *    split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   *  @param index_of_first_interior_vertex Index of vertex in interior of quad 0.
   *    - Integer type.
   *    - Index of vertex in interior of quad i is (index_of_first_interior_vertex+i).
   *    @pre Interior vertices have already been created and
   *    assigned coordinates in vertex_coord[].
   *  @param min_dist_allow_tri4 Allow triangulation into 4 triangles when each
   *    quad vertex is min_distance_allow_tri4 distance from one
   *    of the facets incident on the dual quad edge.
   *    - If some quad vertex is closer than min_distance_allow_tri4 to
   *    both neighboring facets incident on the dual quad edge,
   *    use a diagonal triangulation into two triangles.
   *  @param[out] quad_vert[] Isosurface quadrilateral vertices.
   *    - Quadrilateral i has vertices:
   *      (quad_vert[4*i], quad_vert[4*i+1], quad_vert[4*i+2], quad_vert[4*i+3]).
   *    - If quadrilateral i is deleted, then it's vertices are removed
   *      from quad_vert[].
   *  @param[out] quad_info[i] Information on i'th isosurface quadrilateral.
   *    - If quadrilateral i is deleted, then quad_info[i] is removed
   *      from quad_info[].
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPEQ, typename VTYPET,
            typename ISOQUAD_INFO_TYPE, typename STYPE,
            typename NTYPES2, typename NTYPES4,
            typename ADD_ISOV_FUNCTION>
  void tri2_or_tri4_dual_quad_list_DoutsideE_prefer_tri2_addV
  (const GRID_TYPE & grid,
   VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE quad_info[],
   const NTYPEQ num_quad,
   const STYPE isovalue,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPET> & tri_vert,
   NTYPES2 & num_tri2_split,
   NTYPES4 & num_tri4_split,
   const ADD_ISOV_FUNCTION & add_isov_function)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUAD(4);
    const int TWO_TRIANGLES(2);
    const int dimension = grid.Dimension();
    IJK::PROCEDURE_ERROR error
      ("tri4_dual_quad_list_DoutsideE_prefer_tri2_addV");

    // Initialize
    num_tri2_split = 0;
    num_tri4_split = 0;

    if (num_quad > 0 && quad_info == NULL) {
      error.AddMessage
        ("Programming error.  Array quad_info[] is empty.");
      throw error;
    }

    if (dimension != DIM3) {
      error.AddMessage("Programming error. Vertices are not in dimension three.");
      error.AddMessage("Vertices in dimension: ", dimension);
      throw error;
    }

    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      VTYPEQ * quad_vert_iquad =
        quad_vert + iquad*NUM_VERT_PER_QUAD;

      const int num_tri =
        tri2_or_tri4_dual_quad_DoutsideE_prefer_tri2_addV
        (grid, quad_vert_iquad, quad_info[iquad],
         isovalue, vertex_coord, tri_vert, add_isov_function);

      if (num_tri > 0) {
        collapse_quad_to_vertex(quad_vert_iquad);

        if (num_tri == TWO_TRIANGLES)
          { num_tri2_split++; }
        else
          { num_tri4_split++; }
      }
    }
  }

  
  /*!
   *  @overload
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope, adding vertices, when necessary
   *  - Triangulate into two triangles, is possible.
   *  - Triangulate into four triangles if both diagonals 
   *    are outside the envelope.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - Version using C++ STL vectors vertex_coord[] and quad_vert[].
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPEQ, typename VTYPET,
            typename ISOQUAD_INFO_TYPE, typename STYPE,
            typename NTYPES2, typename NTYPES4,
            typename ADD_ISOV_FUNCTION>
  void tri2_or_tri4_dual_quad_list_DoutsideE_prefer_tri2_addV
  (const GRID_TYPE & grid,
   std::vector<VTYPEQ> & quad_vert,
   const std::vector<ISOQUAD_INFO_TYPE> & quad_info,
   const STYPE isovalue,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPET> & tri_vert,
   NTYPES2 & num_tri2_split,
   NTYPES4 & num_tri4_split,
   const ADD_ISOV_FUNCTION & add_isov_function)
  {
    typedef typename std::vector<VTYPEQ>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (num_quad != quad_info.size()) {
      IJK::PROCEDURE_ERROR error
        ("tri4_dual_quad_list_DoutsideE_prefer_tri2_addV");
      error.AddMessage("Programming error. Incorrect size of quad_info[].");
      error.AddMessage
        ("Size of quad_info[] does not equal number of isosurface quadrilaterals.");
      error.AddMessage("  Size of quad info[]: ", quad_info.size(), ".");
      error.AddMessage("  Number of isosurface quadrilaterals: ",
                       num_quad, ".");
      throw error;
    }

    tri2_or_tri4_dual_quad_list_DoutsideE_prefer_tri2_addV
      (grid, IJK::vector2pointerNC(quad_vert),
       IJK::vector2pointer(quad_info), num_quad,
       isovalue, vertex_coord, tri_vert,
       num_tri2_split, num_tri4_split, add_isov_function);
  }


  // ***************************************************************
  //! @name TRI2 QUADS WITH ONE DIAGONAL OUTSIDE THE ENVELOPE
  // ***************************************************************

  /*!
   *  @brief Split isosurface quadrilaterals with exactly one diagonal
   *    outside the envelope.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If one diagonal is inside the envelope formed by the quad edges
   *    and the dual grid edge and one is outside the envelope,
   *    split using the diagonal inside the envelope.
   *  - Ignore quads that have both diagonals inside the envelope or
   *    both diagonals outside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   *  @param quad_info[i] Information on i'th isosurface quadrilateral.
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPEQ, typename VTYPET, typename ISOQUAD_INFO_TYPE,
            typename NTYPES>
  void tri2_dual_quad_list_with_exactly_one_diagonal_outside_envelope
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE quad_info[],
   const NTYPEQ num_quad,
   std::vector<VTYPET> & tri_vert,
   NTYPES & num_split)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const int NUM_VERT_PER_QUAD(4);
    const DTYPE dimension = grid.Dimension();
    const DTYPE DIM3(3);
    IJK::PROCEDURE_ERROR error
      ("tri2_dual_quad_list_with_one_diagonal_inside_envelope");

    // Initialize
    num_split = 0;

    if (num_quad > 0 && quad_info == NULL) {
      error.AddMessage
        ("Programming error.  Array quad_info[] is empty.");
      throw error;
    }

    if (dimension != DIM3) {
      error.AddMessage("Programming error. Vertices are not in dimension three.");
      error.AddMessage("Vertices in dimension: ", dimension);
      throw error;
    }

    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      VTYPEQ * quad_vert_iquad =
        quad_vert + iquad*NUM_VERT_PER_QUAD;

      if (quad_info[iquad].IsDiagonalInEnvelope(0)) {
        if (!quad_info[iquad].IsDiagonalInEnvelope(1)) {
          // Only diagonal (v0,v2) is in the envelope.
          triangulate_quad_from_vertex(quad_vert_iquad, 0, tri_vert);

          collapse_quad_to_vertex(quad_vert_iquad);
          num_split++;
        }
      }
      else if (quad_info[iquad].IsDiagonalInEnvelope(1)) {
        // Only diagonal (v1,v3) is in the envelope.
        triangulate_quad_from_vertex(quad_vert_iquad, 1, tri_vert);

        collapse_quad_to_vertex(quad_vert_iquad);
        num_split++;
      }

    }
  }


  /*!
   *  @brief Split isosurface quadrilaterals with exactly one diagonal 
   *    outside the envelope.
   *  - Split using the digaonal that is inside the envelope.
   *  - Version using C++ STL vectors vertex_coord[], quad_vert[] and quad_info[].
   *  - Ignore quads that have both diagonals inside the envelope or
   *    both diagonals outside the envelope.
   *  @param[out] num_split Number of quadrilaterals split.
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPEQ, typename VTYPET,
            typename ISOQUAD_INFO_TYPE,
            typename NTYPE>
  void tri2_dual_quad_list_with_exactly_one_diagonal_outside_envelope
  (const GRID_TYPE & grid,
   const std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPEQ> & quad_vert,
   const std::vector<ISOQUAD_INFO_TYPE> & quad_info,
   std::vector<VTYPET> & tri_vert,
   NTYPE & num_split)
  {
    typedef typename std::vector<VTYPEQ>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;


    if (num_quad != quad_info.size()) {
      IJK::PROCEDURE_ERROR error
        ("tri2_dual_quad_list_with_one_diagonal_inside_envelope");
      error.AddMessage
        ("Programming error. Incorrect size of quad_info[].");
      error.AddMessage
        ("  quad_info.size() = ", quad_info.size(), "");
      error.AddMessage
        ("  Expected ", num_quad, ", number of quads.");
      throw error;
    }

    tri2_dual_quad_list_with_exactly_one_diagonal_outside_envelope
      (grid, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointerNC(quad_vert),
       IJK::vector2pointer(quad_info), num_quad,
       tri_vert, num_split);
  }


  ///@}

}

#endif
