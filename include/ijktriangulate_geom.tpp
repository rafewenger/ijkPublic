/*!
 *  \file ijktriangulate_geom.tpp
 *  @brief DEPRECATED: ijk templates for triangulating polytopes (DEPRECATED).
 *  - Includes code based on geometry.
 *  - Version 0.4.0
 *  - Replaced by ijktriangulate_poly2D.tpp and ijktriagulate_poly3D.tpp,
 *    and ijktri2D_angle.tpp.
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2022 Rephael Wenger

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

#ifndef _IJKTRIANGULATE_GEOM_
#define _IJKTRIANGULATE_GEOM_

#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijktriangulate.tpp"
#include "ijktriangulate_poly2D.tpp"


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK_DEPRECATED {

  // ***************************************************************
  //! @name COMPUTE COS MIN TRIANGLE ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute the cosine of smallest angle of triangle (coord0, coord1, coord2).
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   *  - Version that returns POLY_TRIANGULATION_RESULT.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] poly_tri_result Data structure storing
   *    cos_min_triangulation_angle and flag_zero.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   poly_tri_result)
  {
    IJK::compute_cos_min_triangle_angle
      (dimension, coord0, coord1, coord2, max_small_magnitude,
       poly_tri_result.cos_min_triangulation_angle,
       poly_tri_result.flag_zero);

    poly_tri_result.num_triangles = 1;
    poly_tri_result.tri_vertex_index = 0;
    poly_tri_result.num_interior_tri_vertices = 0;
  }

  ///@}


  // ***************************************************************
  //! @name COMPUTE COS MIN TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute the cosine of the smallest triangle angle in polygon triangulation from vertex iv0.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param num_poly_vert Number of polygon vertices.
   *  @param vcoord[] Polygon vertex coordinates.
   *         vcoord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param iv0 Triangulation is from iv0.
   *  @pre iv0 is in range [0,num_poly_vert-1].
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename NTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename IVTYPE,typename COS_TYPE>
  void compute_cos_min_triangulation_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   const IVTYPE index_v0, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    NTYPE i1 = (index_v0+1)%num_poly_vert;
    NTYPE i2 = (i1+1)%num_poly_vert;
    while (i2 != index_v0) {
      
      const CTYPE * v0_coord = vcoord+poly_vert[index_v0]*dimension;
      const CTYPE * v1_coord = vcoord+poly_vert[i1]*dimension;
      const CTYPE * v2_coord = vcoord+poly_vert[i2]*dimension;

      IJK::compute_cos_min_triangle_angle
        (dimension, v0_coord, v1_coord, v2_coord, max_small_magnitude, 
         cosA, flagA);

      if (flagA) { 
        cos_min_angle = 0;
        flag_zero = true;
        return;
      }
      else {
        if (cosA > cos_min_angle) 
          { cos_min_angle = cosA; }
      }
      i1 = i2;
      i2 = (i1+1)%num_poly_vert;
    }

    return;
  }


  /*!  
   *  @brief Compute the cosine of the smallest triangle angle in polygon triangulation from vertex iv0.
   *  - C++ STL vector format for vcoord.
   */
  template <typename DTYPE, typename NTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename IVTYPE,typename COS_TYPE>
  void compute_cos_min_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vcoord, 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   const IVTYPE index_v0, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), num_poly_vert,
       poly_vert, max_small_magnitude, index_v0, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle in polygon triangulation from internal vertex vcoordX[].
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param poly_vcoord[] Pointers to polygon vertex coordinates.
   *         poly_vcoord[i] is a pointer to the i'th vertex coordinates.
   *  @param num_poly_vert Number of polygon vertices.
   *  @param vcoordX[] Internal vertex coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename CTYPEX,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangulation_angle_using_vcoordX
  (const DTYPE dimension, const CTYPE * poly_vcoord[], 
   const NTYPE num_poly_vert, const CTYPEX vcoordX[], 
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
      const NTYPE i1 = (i0+1)%num_poly_vert;
      
      IJK::compute_cos_min_triangle_angle
        (dimension, poly_vcoord[i0], poly_vcoord[i1], vcoordX, 
         max_small_magnitude, cosA, flagA);

      if (flagA) { 
        cos_min_angle = 0;
        flag_zero = true;
        return;
      }
      else {
        if (cosA > cos_min_angle) 
          { cos_min_angle = cosA; }
      }
    }
  }


  /*!
   *  @brief Compute triangulation vertex that maximizes min triangulation angle.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] List of polygon vertices.
   *  @param vcoord[] Polygon vertex coordinates.
   *         vcoord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] iv_max_min Index in poly_vert[] of vertex 
   *    which maximizes min triangulation angle.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if all triangulations have some triangle
   *    with two or three edges with length less than or equal 
   *    to max_small_magnitude.
   */
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename VTYPE, typename MTYPE, 
            typename IVTYPE, typename COS_TYPE>
  void compute_vertex_to_max_min_triangulation_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    // Initialize
    cos_min_angle = 1;
    flag_zero = true;
    iv_max_min = 0;

    for (VTYPE iv = 0; iv < num_poly_vert; iv++) {
      compute_cos_min_triangulation_angle
        (dimension, vcoord, num_poly_vert, poly_vert,
         max_small_magnitude, iv, cosA, flagA);

      if (!flagA && cosA < cos_min_angle) {
        flag_zero = false;
        iv_max_min = iv;
        cos_min_angle = cosA;
      }
    }

    return;
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle in polygon triangulation from vertex iv0.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param num_poly_vert Number of polygon vertices.
   *  @param vcoord[] Polygon vertex coordinates.
   *         vcoord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param iv0 Triangulation is from iv0.
   *  @pre iv0 is in range [0,num_poly_vert-1].
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename NTYPE, typename CTYPE,
            typename MTYPE, typename IVTYPE,typename COS_TYPE>
  void compute_cos_min_poly_v0triangulation_angle
  (const DTYPE dimension, const CTYPE * vcoord[], 
   const NTYPE num_poly_vert, const MTYPE max_small_magnitude, 
   const IVTYPE index_v0, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    NTYPE i1 = (index_v0+1)%num_poly_vert;
    NTYPE i2 = (i1+1)%num_poly_vert;
    while (i2 != index_v0) {
      
      const CTYPE * v0_coord = vcoord[index_v0];
      const CTYPE * v1_coord = vcoord[i1];
      const CTYPE * v2_coord = vcoord[i2];

      IJK::compute_cos_min_triangle_angle
        (dimension, v0_coord, v1_coord, v2_coord, max_small_magnitude, 
         cosA, flagA);

      if (flagA) { 
        cos_min_angle = 0;
        flag_zero = true;
        return;
      }
      else {
        if (cosA > cos_min_angle) 
          { cos_min_angle = cosA; }
      }
      i1 = i2;
      i2 = (i1+1)%num_poly_vert;
    }

    return;
  }


  /*!
   *  @brief Compute triangulation vertex that maximizes min triangulation angle.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] List of polygon vertices.
   *  @param vcoord[] Polygon vertex coordinates.
   *         vcoord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] iv_max_min Index in poly_vert[] of vertex 
   *    which maximizes min triangulation angle.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if all triangulations have some triangle
   *    with two or three edges with length less than or equal 
   *    to max_small_magnitude.
   */
  template <typename DTYPE, typename NTYPE, typename CTYPE, typename MTYPE, 
            typename IVTYPE, typename COS_TYPE>
  void compute_vertex_to_max_min_poly_triangulation_angle
  (const DTYPE dimension, const CTYPE * vcoord[], const NTYPE num_poly_vert, 
   const MTYPE max_small_magnitude, 
   IVTYPE & iv_max_min, COS_TYPE & cos_max_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    // Initialize
    cos_max_min_angle = 1;
    flag_zero = true;
    iv_max_min = 0;

    for (NTYPE iv = 0; iv < num_poly_vert; iv++) {
      compute_cos_min_poly_v0triangulation_angle
        (dimension, vcoord, num_poly_vert, max_small_magnitude, 
         iv, cosA, flagA);

      if (!flagA && cosA < cos_max_min_angle) {
        flag_zero = false;
        iv_max_min = iv;
        cos_max_min_angle = cosA;
      }
    }

    return;
  }

  ///@}


  // ***************************************************************
  //! @name COMPUTE COS MIN QUAD/TRIANGLE TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute the cosine of the smallest triangle angle in quadrilateral triangulation.
   *  - Quadrilateral is split by diagonal (coord0[],coord2[]) 
   *    into two triangles.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename COS_TYPE, typename MTYPE>
  void compute_cos_min_quad_tri02_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA, cosB;
    bool flagA, flagB;

    IJK::compute_cos_min_triangle_angle
      (dimension, coord0, coord1, coord2, max_small_magnitude, cosA, flagA);
    IJK::compute_cos_min_triangle_angle
      (dimension, coord0, coord3, coord2, max_small_magnitude, cosB, flagB);

    if (flagA || flagB) {
      cos_min_angle = 0;
      flag_zero = true;
    }
    else {
      if (cosA > cosB) { cos_min_angle = cosA; }
      else { cos_min_angle = cosB; }
      flag_zero = false;
    }
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle in quadrilateral triangulation.
   *  - Quadrilateral is split by diagonal (coord1[],coord3[]) 
   *    into two triangles.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle have two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, typename MTYPE>
  void compute_cos_min_quad_tri13_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3,
   const MTYPE max_small_magnitude, CTYPE4 & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_quad_tri02_angle
      (dimension, coord1, coord2, coord3, coord0, max_small_magnitude,
       cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle in quadrilateral triangulation.
   *  - Quadrilateral is split into four triangles formed from coordX[]
   *    and each of the four quadrilateral edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the four triangles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param coordX[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPEX, typename COS_TYPE,
            typename MTYPE>
  void compute_cos_min_quad_tri4_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3, const CTYPEX * coordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_QUAD(4);
    COS_TYPE cos[NUM_VERT_PER_QUAD];
    bool flag_tri_zero[NUM_VERT_PER_QUAD];

    IJK::compute_cos_min_triangle_angle
      (dimension, coord0, coord1, coordX, max_small_magnitude, 
       cos[0], flag_tri_zero[0]);
    IJK::compute_cos_min_triangle_angle
      (dimension, coord1, coord2, coordX, max_small_magnitude, 
       cos[1], flag_tri_zero[1]);
    IJK::compute_cos_min_triangle_angle
      (dimension, coord2, coord3, coordX, max_small_magnitude, 
       cos[2], flag_tri_zero[2]);
    IJK::compute_cos_min_triangle_angle
      (dimension, coord3, coord0, coordX, max_small_magnitude, 
       cos[3], flag_tri_zero[3]);

    bool is_cos_min_angle_set = false;
    for (int i = 0; i < NUM_VERT_PER_QUAD; i++) {

      if (flag_tri_zero[i]) {
        flag_zero = true;
        cos_min_angle = 0;
        return;
      }

      if (!is_cos_min_angle_set || cos[i] > cos_min_angle) {
        cos_min_angle = cos[i];
        is_cos_min_angle_set = true;
      }
    }

    flag_zero = false;
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle in splitting of atriangle.
   *  - Triangle is split into two smaller triangles.
   *  - Triangle (v0, v1, v2) is split by edge (vX, v1) 
   *    into (vX, v0, v1) and (vX, v1, v2).
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vcoord0[] Input coordinates.
   *  @param vcoord1[] Input coordinates.
   *  @param vcoord2[] Input coordinates.
   *  @param vcoordX[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_split_angle Cosine of min split triange angle.
   *  @param[out] flag_zero True, if some split triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPEX, typename COS_TYPE, typename MTYPE>
  void compute_cos_min_split_triangle_angle
  (const DTYPE dimension, const CTYPE0 * vcoord0, const CTYPE1 * vcoord1, 
   const CTYPE2 * vcoord2, const CTYPEX * vcoordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_quad_tri02_angle
      (dimension, vcoordX, vcoord0, vcoord1, vcoord2, max_small_magnitude,
       cos_min_angle, flag_zero);
  }

  ///@}


  // ***************************************************************
  //! @name COMPUTE COS MIN PENTAGON TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute the cosine of the smallest triangle angle in pentagon triangulation.
   *  - Pentagon is split into five triangles formed from vcoordX[]
   *    and each of the five pentagon edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the five triangles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vcoord0[] Coordinates of pentagon vertex 0.
   *  @param vcoord1[] Coordinates of pentagon vertex 1.
   *  @param vcoord2[] Coordinates of pentagon vertex 2.
   *  @param vcoord3[] Coordinates of pentagon vertex 3.
   *  @param vcoord4[] Coordinates of pentagon vertex 4.
   *  @param vcoordX[] Coordinates of vertex at center of star.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, typename CTYPEX,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_min_pentagon_tri5_angle
  (const DTYPE dimension, const CTYPE0 * vcoord0, const CTYPE1 * vcoord1, 
   const CTYPE2 * vcoord2, const CTYPE3 * vcoord3, const CTYPE4 * vcoord4,
   const CTYPEX * vcoordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    COS_TYPE cos[NUM_VERT_PER_PENTAGON];
    bool flag_tri_zero[NUM_VERT_PER_PENTAGON];

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoordX, max_small_magnitude,
       cos[0], flag_tri_zero[0]);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord2, vcoordX, max_small_magnitude,
       cos[1], flag_tri_zero[1]);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord2, vcoord3, vcoordX, max_small_magnitude,
       cos[2], flag_tri_zero[2]);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord3, vcoord4, vcoordX, max_small_magnitude,
       cos[3], flag_tri_zero[3]);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord4, vcoord0, vcoordX, max_small_magnitude,
       cos[4], flag_tri_zero[4]);

    bool is_cos_min_angle_set = false;
    for (int i = 0; i < NUM_VERT_PER_PENTAGON; i++) {

      if (flag_tri_zero[i]) {
        flag_zero = true;
        cos_min_angle = 0;
        return;
      }

      if (!is_cos_min_angle_set || cos[i] > cos_min_angle) {
        cos_min_angle = cos[i];
        is_cos_min_angle_set = true;
      }
    }

    flag_zero = false;
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle in pentagon triangulation.
   *  - Pentagon is split into five triangles formed from vcoordX[]
   *    and each of the five pentagon edges.
   *  - Version with pentagon vertices listed in array pentagon_vert[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vertex_coord[] Vertex coordinates.
   *    vertex_coord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param pentagon_vert[] List of pentagon vertices in clockwise
   *    or counter-clockwise order around the pentagon.
   *  @param ivX Triangulate into five triangles by connecting each
   *     pentagon edge to ivX.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE0, typename VTYPE1,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_min_pentagon_tri5_angle
  (const DTYPE dimension, const CTYPE vertex_coord[], 
   const VTYPE0 pentagon_vert[], const VTYPE1 ivX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const CTYPE * vcoord0 = vertex_coord + pentagon_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord + pentagon_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord + pentagon_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord + pentagon_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord + pentagon_vert[4]*dimension;
    const CTYPE * vcoordX = vertex_coord + ivX*dimension;

    compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoordX,
       max_small_magnitude, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle in pentagon triangulation.
   *  - Pentagon is split into five triangles formed from vcoordX[]
   *    and each of the five pentagon edges.
   *  - Version with pentagon vertices listed in array pentagon_vert[].
   *  - Version using C++ STL vector vertex_coord[].
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE0, typename VTYPE1,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_min_pentagon_tri5_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const VTYPE0 pentagon_vert[], const VTYPE1 ivX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_pentagon_tri5_angle
      (dimension, IJK::vector2pointer(vertex_coord), pentagon_vert, ivX,
       max_small_magnitude, cos_min_angle, flag_zero);
  }


  /// Compute min angle in triangulation of a pentagon with all triangles incident on vertex_coord0.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_pentagon_triangulation_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cos012, cos023, cos034;
    bool flag_zero_012, flag_zero_023, flag_zero_034;

    // Initialize.
    cos_min_angle = 1;
    flag_zero = true;

    IJK::compute_cos_min_triangle_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, 
       max_small_magnitude, cos012, flag_zero_012);
    if (flag_zero_012) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vertex_coord0, vertex_coord2, vertex_coord3, 
       max_small_magnitude, cos023, flag_zero_023);
    if (flag_zero_023) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vertex_coord0, vertex_coord3, vertex_coord4, 
       max_small_magnitude, cos034, flag_zero_034);
    if (flag_zero_034) { return; }

    cos_min_angle = cos012;
    if (cos_min_angle < cos023) { cos_min_angle = cos023; }
    if (cos_min_angle < cos034) { cos_min_angle = cos034; }

    flag_zero = false;
  }


  /// @brief Compute min angle in triangulation of a pentagon with all triangles incident on vertex_coord0.
  /// - Version returning data structure POLY_TRIANGULATION_RESULT.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_min_pentagon_triangulation_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const MTYPE max_small_magnitude,
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & pentagon_tri)
  {
    compute_cos_min_pentagon_triangulation_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, max_small_magnitude,
       pentagon_tri.cos_min_triangulation_angle, pentagon_tri.flag_zero);

    pentagon_tri.num_triangles = 3;
  }

  ///@}


  // ***************************************************************
  //! @name COMPUTE COS MAX MIN QUAD TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals.
   *  - Quadrilateral is either split by diagonal (coord0[], coord2[])
   *    or by diagonal (coord1[],coord3[]) into two triangles.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_diag02 If true, triangulation with min angle
   *       has diagonal (coord0[], coord2[]).  Otherwise, triangulation 
   *       with min angle has diagonal (coord1[], coord2[]).
   *  @param[out] flag_zero True, if some triangles in triangulation
   *    with min angle have two or three edges less than or equal
   *    to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename COS_TYPE, typename MTYPE>
  void compute_cos_max_min_quad_tri02_or_tri13_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, 
   bool & flag_diag02, bool & flag_zero)
  {
    COS_TYPE cos_min_diag02, cos_min_diag13;
    bool flag_zero_diag02, flag_zero_diag13;

    compute_cos_min_quad_tri02_angle
      (dimension, coord0, coord1, coord2, coord3, max_small_magnitude, 
       cos_min_diag02, flag_zero_diag02);
    compute_cos_min_quad_tri13_angle
      (dimension, coord0, coord1, coord2, coord3, max_small_magnitude, 
       cos_min_diag13, flag_zero_diag13);

    int index_selected;
    IJK::select_min
      (cos_min_diag02, flag_zero_diag02, 
       cos_min_diag13, flag_zero_diag13, 
       cos_min_angle, index_selected, flag_zero);

    if (index_selected == 0) { flag_diag02 = true; }
    else { flag_diag02 = false; }
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals.
   *  - Version with list of pointers to input coordinates 
   *    in array coord_list[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord_list[] Array of pointers to 4 input coordinates.
   */
  template <typename DTYPE, typename CTYPE,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_max_min_quad_tri02_or_tri13_angle
  (const DTYPE dimension, const CTYPE * coord_list[],
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, 
   bool & flag_diag02, bool & flag_zero)
  {
    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, coord_list[0], coord_list[1], coord_list[2], coord_list[3],
       max_small_magnitude, cos_min_angle, flag_diag02, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals.
   *  - Version with quad vertices and array of vertex coordinates.
   *  - Quad vertices are:
   *       (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3]),
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_max_min_quad_tri02_or_tri13_angle
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const VTYPE quad_vert, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_diag02, bool & flag_zero)
  {
    const CTYPE * vcoord0 = vertex_coord + dimension*quad_vert[0];
    const CTYPE * vcoord1 = vertex_coord + dimension*quad_vert[1];
    const CTYPE * vcoord2 = vertex_coord + dimension*quad_vert[2];
    const CTYPE * vcoord3 = vertex_coord + dimension*quad_vert[3];

    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
       max_small_magnitude, cos_min_angle, flag_diag02, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals.
   *  - Version with quad vertices and array of vertex coordinates.
   *  - Version using C++ STL vector vertex_coord[].
   *  - Quad vertices are:
   *       (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3]),
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_max_min_quad_tri02_or_tri13_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const VTYPE quad_vert, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_diag02, bool & flag_zero)
  {
    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, IJK::vector2pointer(vertex_coord), quad_vert,
       max_small_magnitude, cos_min_angle, flag_diag02, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals.
   *  - *** DEPRECATED ***
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  @param[out] quad_tri Quad triangulation result.
   *  - Returns quad_tri.tri_vertex_index set to 0 or 1.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename MTYPE,
            typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3,
   const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {
    COS_TYPE cos_min_diag02, cos_min_diag13;
    bool flag_zero_diag02, flag_zero_diag13;

    compute_cos_min_quad_tri02_angle
      (dimension, coord0, coord1, coord2, coord3, max_small_magnitude, 
       cos_min_diag02, flag_zero_diag02);
    compute_cos_min_quad_tri13_angle
      (dimension, coord0, coord1, coord2, coord3, max_small_magnitude, 
       cos_min_diag13, flag_zero_diag13);

    int index_selected;
    IJK::select_min
      (cos_min_diag02, flag_zero_diag02, 
       cos_min_diag13, flag_zero_diag13, 
       quad_tri.cos_min_triangulation_angle, index_selected, 
       quad_tri.flag_zero);

    quad_tri.num_triangles = 2;
    quad_tri.num_interior_tri_vertices = 0;
    quad_tri.tri_vertex_index = index_selected;
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals.
   *  - *** DEPRECATED ***
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param[out] quad_tri Quad triangulation result.
   */
  template <typename DTYPE, typename CTYPE, typename MTYPE,
            typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_angle
  (const DTYPE dimension, const CTYPE * vcoord[4],
   const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {
    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord[0], vcoord[1], vcoord[2], vcoord[3], 
       max_small_magnitude, quad_tri);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals.
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param[out] quad_tri Quad triangulation result.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_angle
  (const DTYPE dimension, const CTYPE * vertex_coord, 
   const VTYPE quad_vert[], const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {

    const CTYPE * vcoord0 = vertex_coord + dimension*quad_vert[0];
    const CTYPE * vcoord1 = vertex_coord + dimension*quad_vert[1];
    const CTYPE * vcoord2 = vertex_coord + dimension*quad_vert[2];
    const CTYPE * vcoord3 = vertex_coord + dimension*quad_vert[3];

    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
       max_small_magnitude, quad_tri);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals.
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param[out] quad_tri Quad triangulation result.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const VTYPE quad_vert[], const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);

    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vertex_coord_ptr, quad_vert, max_small_magnitude, 
       quad_tri);
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle in quadrilateral triangulation.
   *  - Quadrilateral is split into four triangles formed from coordX[]
   *    and each of the four quadrilateral edges.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Quadrilateral is split into four triangles formed from coordX[]
   *    and each of the four quadrilateral edges.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_quad_tri4_angle
  (const DTYPE dimension, const CTYPE * vcoord[4], const CTYPEX coordX[],
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_quad_tri4_angle
      (dimension, vcoord[0], vcoord[1], vcoord[2], vcoord[3], coordX,
       max_small_magnitude, cos_min_angle, flag_zero);
  }

  /*!
   *  @brief Compute the cosine of the smallest triangle angle in quadrilateral triangulation.
   *  - Quadrilateral is split into four triangles formed from coordX[]
   *    and each of the four quadrilateral edges.
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Quadrilateral is split into four triangles formed from coordX[]
   *    and each of the four quadrilateral edges.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPE, typename MTYPE, typename COS_TYPE>
  void compute_cos_min_quad_tri4_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const CTYPEX coordX[], const VTYPE quad_vert[],
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord[NUM_VERT_PER_QUAD] =
      { vertex_coord_ptr+quad_vert[0]*dimension,
        vertex_coord_ptr+quad_vert[1]*dimension,
        vertex_coord_ptr+quad_vert[2]*dimension,
        vertex_coord_ptr+quad_vert[3]*dimension };

    compute_cos_min_quad_tri4_angle
      (dimension, vcoord, coordX, max_small_magnitude, 
       cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals
   *    or triangulation into four triangles.
   *  - Quadrilateral is either split by diagonal (coord0[], coord2[])
   *    or by diagonal (coord1[],coord3[]) into two triangles or
   *    split into 4 triangles all incident on coordX[].
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param coordX[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_tri4 If true, split into 4 triangles.
   *  @param[out] flag_diag02 If true, triangulation with min angle
   *       has diagonal (coord0[], coord2[]).  Otherwise, triangulation 
   *       with min angle has diagonal (coord1[], coord2[]).
   *  @param[out] flag_zero True, if some triangles in triangulation
   *    with min angle have two or three edges less than or equal
   *    to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPEX, typename COS_TYPE,
            typename MTYPE>
  void compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, const CTYPE0 * vcoord0, const CTYPE1 * vcoord1, 
   const CTYPE2 * vcoord2, const CTYPE3 * vcoord3, const CTYPEX * vcoordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, 
   bool & flag_tri4, bool & flag_diag02, bool & flag_zero)
  {
    COS_TYPE cos_min_diag, cos_min_tri4;
    bool flag_zero_diag, flag_zero_tri4;

    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
       cos_min_diag, flag_diag02, flag_zero_diag);
    compute_cos_min_quad_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoordX,
       max_small_magnitude, cos_min_tri4, flag_zero_tri4);

    int index_selected;
    IJK::select_min
      (cos_min_tri4, flag_zero_tri4, 
       cos_min_diag, flag_zero_diag, 
       cos_min_angle, index_selected, flag_zero);

    if (index_selected == 0) { flag_tri4 = true; }
    else { flag_tri4 = false; }
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals
   *    or triangulation into four triangles.
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  @param[out] quad_tri_result Quad triangulation result.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPEX, typename MTYPE,
            typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3, const CTYPEX * vcoordX,
   const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    bool flag_tri4, flag_diag02;

    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, coord0, coord1, coord2, coord3, vcoordX,
       max_small_magnitude, quad_tri_result.cos_min_triangulation_angle,
       flag_tri4, flag_diag02, quad_tri_result.flag_zero);

    if (flag_diag02) 
      { quad_tri_result.tri_vertex_index = 0; }
    else
      { quad_tri_result.tri_vertex_index = 1; }

    if (flag_tri4) {
      quad_tri_result.num_interior_tri_vertices = 1;
      quad_tri_result.tri_vertex_index = 4;
      quad_tri_result.num_triangles = 4;
    }
    else {
      quad_tri_result.num_interior_tri_vertices = 0;
      quad_tri_result.num_triangles = 2;
    }
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals
   *    or triangulation into four triangles.
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  - Version with array of coordinates.
   *  @param[out] quad_tri_result Quad triangulation result.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, const CTYPE * quad_coord[4],
   const CTYPEX * vcoordX, const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> &
   quad_tri_result)
  {
    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, quad_coord[0], quad_coord[1], 
       quad_coord[2], quad_coord[3], vcoordX, max_small_magnitude,
       quad_tri_result);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals
   *    or triangulation into four triangles.
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  - Version with quad vertices and array of vertex coordinates.
   *  - Quad vertices are:
   *       (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3]),
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename VTYPE, typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, const CTYPE * vertex_coord, const CTYPEX * vcoordX,
   const VTYPE quad_vert, const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    const CTYPE * vcoord0 = vertex_coord + dimension*quad_vert[0];
    const CTYPE * vcoord1 = vertex_coord + dimension*quad_vert[1];
    const CTYPE * vcoord2 = vertex_coord + dimension*quad_vert[2];
    const CTYPE * vcoord3 = vertex_coord + dimension*quad_vert[3];

    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoordX, 
       max_small_magnitude, quad_tri_result);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals
   *    or triangulation into four triangles.
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  - Version with quad vertices and array of vertex coordinates.
   *  - Version using C++ STL vector vertex_coord[].
   *  - Quad vertices are:
   *       (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3]),
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, typename VTYPE, 
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX, const VTYPE quad_vert, 
   const MTYPE max_small_magnitude, 
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, IJK::vector2pointer(vertex_coord), vcoordX, quad_vert, 
       max_small_magnitude, quad_tri_result);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle.
   *  - Max min angle over the two triangulations using the quad diagonals
   *    or triangulation into four triangles.
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  - Version with quad vertices and array of vertex coordinates.
   *  - Version using C++ STL vector vertex_coord[].
   *  - Version including block_diagonal.
   *  - Quad vertices are:
   *       (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3]),
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, typename VTYPE, 
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX, const VTYPE quad_vert, 
   const MTYPE max_small_magnitude, 
   const bool block_diagonal[2],
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    if (block_diagonal[0] || block_diagonal[1]) {
      compute_cos_min_quad_tri4_angle
        (dimension, vertex_coord, vcoordX, quad_vert, max_small_magnitude, 
         quad_tri_result.cos_min_triangulation_angle, quad_tri_result.flag_zero);
      quad_tri_result.num_triangles = 4;
      quad_tri_result.num_interior_tri_vertices = 1;
      quad_tri_result.tri_vertex_index = 4;
    }
    else {
      compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
        (dimension, vertex_coord, vcoordX, quad_vert, max_small_magnitude, 
         quad_tri_result);
    }
  }

  ///@}


  // ***************************************************************
  //! @name COMPUTE COS MAX MIN PENTAGON TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*! 
   *  @brief Compute vertex to maximize min triangulation angle of a pentagon.
   *  - Faster than calling compute_vertex_to_max_min_triangulation_angle.
   *  @param vertex_coord0 Coordinates of vertex 0.
   *  @param vertex_coord1 Coordinates of vertex 1.
   *  @param vertex_coord2 Coordinates of vertex 2.
   *  @param vertex_coord3 Coordinates of vertex 3.
   *  @param vertex_coord4 Coordinates of vertex 4.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *    equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] iv_max_min Index in pentagon_vert[] of vertex 
   *    which maximizes min triangulation angle.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if all triangulations have some triangle
   *    with two or three edges with length less than or equal 
   *    to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE, typename MTYPE, 
            typename IVTYPE, typename COS_TYPE>
  void compute_vertex_to_max_min_pentagon_triangulation_angle
  (const DTYPE dimension, const CTYPE vertex_coord0[], 
   const CTYPE vertex_coord1[], const CTYPE vertex_coord2[], 
   const CTYPE vertex_coord3[], const CTYPE vertex_coord4[],
   const MTYPE max_small_magnitude,
   IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    const CTYPE * vcoord[NUM_VERT_PER_PENTAGON] =
      { vertex_coord0, vertex_coord1, vertex_coord2, 
        vertex_coord3, vertex_coord4 };
    CTYPE cos_ear_triangle[NUM_VERT_PER_PENTAGON];
    bool flag_zero_ear[NUM_VERT_PER_PENTAGON];
    CTYPE cos_center_triangle[NUM_VERT_PER_PENTAGON];
    bool flag_zero_center_triangle[NUM_VERT_PER_PENTAGON];

    // Initialize
    iv_max_min = 0;
    cos_min_angle = 1;
    flag_zero = true;

    for (int i0 = 0; i0 < NUM_VERT_PER_PENTAGON; i0++) {
      const int i1 = (i0+1)%NUM_VERT_PER_PENTAGON;
      const int i2 = (i0+2)%NUM_VERT_PER_PENTAGON;
      const int i4 = (i0+4)%NUM_VERT_PER_PENTAGON;

      IJK::compute_cos_min_triangle_angle
        (dimension, vcoord[i0], vcoord[i1], vcoord[i2], max_small_magnitude, 
         cos_ear_triangle[i1], flag_zero_ear[i1]);

      IJK::compute_cos_min_triangle_angle
        (dimension, vcoord[i0], vcoord[i2], vcoord[i4], max_small_magnitude, 
         cos_center_triangle[i2], flag_zero_center_triangle[i2]);
    }

    bool flag_found = false;
    for (int i0 = 0; i0 < NUM_VERT_PER_PENTAGON; i0++) {
      const int i1 = (i0+1)%NUM_VERT_PER_PENTAGON;
      const int i2 = (i0+2)%NUM_VERT_PER_PENTAGON;

      if (flag_zero_ear[i0] || flag_zero_ear[i2] ||
          flag_zero_center_triangle[i1]) { continue; }

      COS_TYPE cos_min_i1 = cos_center_triangle[i1];
      if (cos_ear_triangle[i0] > cos_min_i1)
        { cos_min_i1 = cos_ear_triangle[i0]; }
      if (cos_ear_triangle[i2] > cos_min_i1)
        { cos_min_i1 = cos_ear_triangle[i2]; }

      if (flag_found) {
        if (cos_min_angle > cos_min_i1) {
          iv_max_min = i1;
          cos_min_angle = cos_min_i1;
        }
      }
      else {
        iv_max_min = i1;
        cos_min_angle = cos_min_i1;
        flag_found = true;
      }
    }

    flag_zero = !flag_found;
  }


  /*!
   *  @brief Compute vertex to maximize min triangulation angle of a pentagon.
   *  - Faster than calling compute_vertex_to_max_min_triangulation_angle.
   *  @param vertex_coord[] List of vertex coordinates.
   *    vertex_coord[i*3+j] is j'th coordinate of i'th vertex.
   *  @param pentagon_vert[] List of pentagon vertices in clockwise
   *    or counter-clockwise order around the pentagon.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *     equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] iv_max_min Index in pentagon_vert[] of vertex 
   *     which maximizes min triangulation angle.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if all triangulations have some triangle
   *     with two or three edges with length less than or equal 
   *     to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE, typename MTYPE,
            typename IVTYPE, typename COS_TYPE>
  void compute_vertex_to_max_min_pentagon_triangulation_angle
  (const DTYPE dimension, 
   const CTYPE vertex_coord[],
   const VTYPE pentagon_vert[],
   const MTYPE max_small_magnitude,
   IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    const CTYPE * vcoord[NUM_VERT_PER_PENTAGON] =
      { vertex_coord+pentagon_vert[0]*dimension,
        vertex_coord+pentagon_vert[1]*dimension,
        vertex_coord+pentagon_vert[2]*dimension,
        vertex_coord+pentagon_vert[3]*dimension,
        vertex_coord+pentagon_vert[4]*dimension };

    compute_vertex_to_max_min_pentagon_triangulation_angle
      (dimension, vcoord[0], vcoord[1], vcoord[2], vcoord[3], vcoord[4],
       max_small_magnitude, iv_max_min, cos_min_angle, flag_zero);
  }


  /*! 
   *  @brief Compute vertex to maximize min triangulation angle of a pentagon.
   *  - Version using C++ STL vector for array vertex_coord[].
   *  - Faster than calling compute_vertex_to_max_min_triangulation_angle.
   *  @param vertex_coord[] List of vertex coordinates.
   *    vertex_coord[i*3+j] is j'th coordinate of i'th vertex.
   *  @param pentagon_vert[] List of pentagon vertices in clockwise
   *    or counter-clockwise order around the pentagon.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *     equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] iv_max_min Index in pentagon_vert[] of vertex 
   *     which maximizes min triangulation angle.
   *   @param[out] cos_min_angle Cosine of min triange angle.
   *   @param[out] flag_zero True, if all triangulations have some triangle
   *    with two or three edges with length less than or equal 
   *    to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE, typename MTYPE,
            typename IVTYPE, typename COS_TYPE>
  void compute_vertex_to_max_min_pentagon_triangulation_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE pentagon_vert[],
   const MTYPE max_small_magnitude,
   IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_vertex_to_max_min_pentagon_triangulation_angle
      (dimension, IJK::vector2pointer(vertex_coord), pentagon_vert,
       max_small_magnitude, iv_max_min, cos_min_angle, flag_zero);
  }


  /*! 
   *  @brief Compute max min triangulation angle of a pentagon.
   *  - Max min angle over triangulations into 3 triangles
   *    incident on a pentagon vertex or into 5 triangles 
   *    incident on vcoordX[].
   *  @param vcoordX[] Coordinates of vertex at center of star.
   *  @param[out] iv_max_min Index 0,1,2,3 or 4, of pentagon vertex
   *    which maximizes the min angle in triangulation into 3 triangles.
   */
  template <typename DTYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename CTYPE3, typename CTYPE4, typename CTYPEX,
            typename MTYPE, typename COS_TYPE, typename NTYPE>
  void compute_cos_max_min_pentagon_tri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPEX vcoordX[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_tri5, 
   NTYPE & iv_max_min, bool & flag_zero)
  {
    COS_TYPE cos_min_tri3_angle, cos_min_tri5_angle;
    bool flag_tri3_zero, flag_tri5_zero;

    compute_vertex_to_max_min_pentagon_triangulation_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vertex_coord4, max_small_magnitude, iv_max_min,
       cos_min_tri3_angle, flag_tri3_zero);
    compute_cos_min_pentagon_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vertex_coord4, vcoordX,
       max_small_magnitude, cos_min_tri5_angle, flag_tri5_zero);

    if (flag_tri3_zero) { flag_tri5 = true; }
    else if (flag_tri5_zero) { flag_tri5 = false; }
    else if (cos_min_tri5_angle <= cos_min_tri3_angle)
      { flag_tri5 = true; }
    else 
      { flag_tri5 = false; }

    if (flag_tri5) {
      flag_tri5 = true;
      flag_zero = flag_tri5_zero;
      cos_min_angle = cos_min_tri5_angle;
    }
    else {
      flag_tri5 = false;
      flag_zero = flag_tri3_zero;
      cos_min_angle = cos_min_tri3_angle;
    }

  }


  /*!
   *  @brief Compute max min triangulation angle of a pentagon.
   *  - Max min angle over triangulations into 3 triangles
   *    incident on a pentagon vertex or into 5 triangles 
   *    incident on vcoordX[].
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  @param vcoordX[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename CTYPE3, typename CTYPE4, typename CTYPEX,
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_pentagon_tri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPEX vcoordX[],
   const MTYPE max_small_magnitude,
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> &
   pentagon_tri_result)
  {
    bool flag_tri5;

    compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, vcoordX, max_small_magnitude,
       pentagon_tri_result.cos_min_triangulation_angle, flag_tri5,
       pentagon_tri_result.tri_vertex_index, pentagon_tri_result.flag_zero);

    if (flag_tri5) {
      pentagon_tri_result.num_triangles = 5;
      pentagon_tri_result.num_interior_tri_vertices = 1;
    }
    else {
      pentagon_tri_result.num_triangles = 3;
      pentagon_tri_result.num_interior_tri_vertices = 0;
    }

  }


  /*! 
   *  @brief Compute max min triangulation angle of a pentagon.
   *  - Max min angle over triangulations into 3 triangles
   *    incident on a pentagon vertex or into 5 triangles 
   *    incident on vcoordX[].
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  - Version with pentagon vertices and array of vertex coordinates.
   *  - Pentagon vertices are:
   *       (pentagon_vert[0], pentagon_vert[1], ..., pentagon_vert[4])
   *    listed in clockwise/counter-clockwise order around the pentagon.
   * @param vcoordX[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, typename VTYPE,
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_pentagon_tri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const CTYPEX vcoordX[],
   const VTYPE pentagon_vert[],
   const MTYPE max_small_magnitude,
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & pentagon_tri_result)
  {
    const CTYPE * vcoord0 = vertex_coord + dimension*pentagon_vert[0];
    const CTYPE * vcoord1 = vertex_coord + dimension*pentagon_vert[1];
    const CTYPE * vcoord2 = vertex_coord + dimension*pentagon_vert[2];
    const CTYPE * vcoord3 = vertex_coord + dimension*pentagon_vert[3];
    const CTYPE * vcoord4 = vertex_coord + dimension*pentagon_vert[4];

    compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoordX, 
       max_small_magnitude, pentagon_tri_result);
  }


  /*! 
   *  @brief Compute max min triangulation angle of a pentagon.
   *  - Max min angle over triangulations into 3 triangles
   *    incident on a pentagon vertex or into 5 triangles 
   *    incident on vcoordX[].
   *  - Version using C++ STL vector vertex_coord[].
   *  - Pentagon vertices are:
   *       (pentagon_vert[0], pentagon_vert[1], ..., pentagon_vert[4])
   *    listed in clockwise/counter-clockwise order around the pentagon.
   *  @param vcoordX[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, typename VTYPE,
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_pentagon_tri3_or_tri5_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX vcoordX[],
   const VTYPE pentagon_vert[],
   const MTYPE max_small_magnitude,
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> &
   pentagon_tri_result)
  {
    compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, IJK::vector2pointer(vertex_coord), vcoordX, pentagon_vert,
       max_small_magnitude, pentagon_tri_result);
  }

  ///@}


  // ****************************************************************
  //! @name COMPUTE COS MIN HEXAGON ANGLE
  // ****************************************************************

  /// @brief Compute min angle in triangulation of hexagon
  /// - All triangles are incident on vertex_coord0.
  template <typename DTYPE, typename CTYPE, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_hexagon_triangulation_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord0,
   const CTYPE * vertex_coord1,
   const CTYPE * vertex_coord2,
   const CTYPE * vertex_coord3,
   const CTYPE * vertex_coord4,
   const CTYPE * vertex_coord5,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cos012, cos023, cos034, cos045;
    bool flag_zero_012, flag_zero_023, flag_zero_034, flag_zero_045;

    // Initialize.
    cos_min_angle = 1;
    flag_zero = true;

    IJK::compute_cos_min_triangle_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, 
       max_small_magnitude, cos012, flag_zero_012);
    if (flag_zero_012) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vertex_coord0, vertex_coord2, vertex_coord3, 
       max_small_magnitude, cos023, flag_zero_023);
    if (flag_zero_023) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vertex_coord0, vertex_coord3, vertex_coord4, 
       max_small_magnitude, cos034, flag_zero_034);
    if (flag_zero_034) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vertex_coord0, vertex_coord4, vertex_coord5, 
       max_small_magnitude, cos045, flag_zero_045);
    if (flag_zero_045) { return; }

    cos_min_angle = cos012;
    if (cos_min_angle < cos023) { cos_min_angle = cos023; }
    if (cos_min_angle < cos034) { cos_min_angle = cos034; }
    if (cos_min_angle < cos045) { cos_min_angle = cos045; }

    flag_zero = false;
  }


  // ****************************************************************
  //! @name COMPUTE COS MAX MIN HEXAGON TRIANGULATION ANGLE
  // ****************************************************************

  ///@{

  /*!
   *  @brief Compute vertex to maximize min triangulation angle of a hexagon.
   *  - Faster than calling compute_vertex_to_max_min_triangulation_angle.
   *  @param vertex_coord0 Coordinates of vertex 0.
   *  @param vertex_coord1 Coordinates of vertex 1.
   *  @param vertex_coord2 Coordinates of vertex 2.
   *  @param vertex_coord3 Coordinates of vertex 3.
   *  @param vertex_coord4 Coordinates of vertex 4.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *     equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] iv_max_min Index in pentagon_vert[] of vertex 
   *     which maximizes min triangulation angle.
   *   @param[out] cos_min_angle Cosine of min triange angle.
   *   @param[out] flag_zero True, if all triangulations have some triangle
   *     with two or three edges with length less than or equal 
   *     to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE, typename MTYPE, 
            typename IVTYPE, typename COS_TYPE>
  void compute_vertex_to_max_min_hexagon_triangulation_angle
  (const DTYPE dimension, 
   const CTYPE vertex_coord0[], const CTYPE vertex_coord1[], 
   const CTYPE vertex_coord2[], const CTYPE vertex_coord3[], 
   const CTYPE vertex_coord4[], const CTYPE vertex_coord5[],
   const MTYPE max_small_magnitude,
   IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_HEXAGON(6);
    const CTYPE * vcoord[NUM_VERT_PER_HEXAGON] =
      { vertex_coord0, vertex_coord1, vertex_coord2, 
        vertex_coord3, vertex_coord4, vertex_coord5 };

    compute_vertex_to_max_min_poly_triangulation_angle
      (dimension, vcoord, NUM_VERT_PER_HEXAGON, max_small_magnitude, 
       iv_max_min, cos_min_angle, flag_zero);
  }
   

  /*!
   *  @brief Compute vertex to maximize min triangulation angle of a hexagon.
   *  - Max min angle over 4 triangles incident on a hexagon vertex
   *    or 6 triangles incident on vcoordX[].
   *  @param vcoordX[] Coordinates of vertex at center of star.
   *  @param[out] iv_max_min Index 0,1,2,3, 4 or 5, of hexagon vertex
   *    which maximizes the min angle in triangulation into 4 triangles.
   *    - If (iv_max_min == 6), then triangulate from vcoordX[].
   */
  template <typename DTYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename CTYPE3, typename CTYPE4, typename CTYPE5,
            typename CTYPEX, typename MTYPE, typename COS_TYPE, typename NTYPE>
  void compute_cos_max_min_hexagon_tri4_or_tri6_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPE5 vertex_coord5[],
   const CTYPEX vcoordX[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_max_min_angle, 
   NTYPE & num_triangles,
   NTYPE & iv_max_min, bool & flag_zero)
  {
    COS_TYPE cos_max_min_tri4_angle, cos_min_tri6_angle;
    bool flag_zero_tri4, flag_zero_tri6;

    compute_vertex_to_max_min_hexagon_triangulation_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vertex_coord4, vertex_coord5, max_small_magnitude, iv_max_min,
       cos_max_min_tri4_angle, flag_zero_tri4);
    IJK::compute_cos_min_hexagon_tri6_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vertex_coord4, vertex_coord5, vcoordX,
       max_small_magnitude, cos_min_tri6_angle, flag_zero_tri6);


    // Initialize
    num_triangles = 4;

    NTYPE index;
    IJK::select_min(cos_max_min_tri4_angle, flag_zero_tri4,
               cos_min_tri6_angle, flag_zero_tri6,
               cos_max_min_angle, index, flag_zero);

    if (index == 1) {
      num_triangles = 6;

      // If triangulation from vcoordX[], set iv_max_min to 6.
      iv_max_min = 6;
    }
  }


  /*!
   *  @brief Compute vertex to maximize min triangulation angle of a hexagon.
   *  - Max min angle over 4 triangles incident on a hexagon vertex
   *    or 6 triangles incident on vcoordX[].
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  @param vcoordX[] Coordinates of vertex at center of star.
   *  @param[out] iv_max_min Index 0,1,2,3, 4 or 5, of hexagon vertex
   *    which maximizes the min angle in triangulation into 4 or 6 triangles.
   */
  template <typename DTYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename CTYPE3, typename CTYPE4, typename CTYPE5,
            typename CTYPEX, typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_hexagon_tri4_or_tri6_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPE5 vertex_coord5[],
   const CTYPEX vcoordX[],
   const MTYPE max_small_magnitude,
   IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> &
   hexagon_tri_result)
  {
    compute_cos_max_min_hexagon_tri4_or_tri6_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, vertex_coord5, vcoordX,
       max_small_magnitude, hexagon_tri_result.cos_min_triangulation_angle,
       hexagon_tri_result.num_triangles, hexagon_tri_result.tri_vertex_index,
       hexagon_tri_result.flag_zero);

    if (hexagon_tri_result.num_triangles == 6) {
      hexagon_tri_result.num_interior_tri_vertices = 1;
    }
    else {
      hexagon_tri_result.num_interior_tri_vertices = 0;
    }

  }

  ///@}


  // **************************************************
  //! @name COMPUTE DISTANCE
  // **************************************************

  ///@{

  /*!
   *  @brief Compute the distance from pcoord[] to the closest grid facet
   *    incident on edge.
   *  @tparam GRID_TYPE Must have function SpacingPtrConst().
   *  @param grid Regular grid.
   */
  template <typename GRID_TYPE, typename CTYPE, 
            typename VTYPE, typename DIR_TYPE, typename DIST_TYPE>
  void compute_unscaled_distance_to_closest_facet
  (const GRID_TYPE & grid, const CTYPE pcoord[],
   const VTYPE endpoint0, const DIR_TYPE edge_dir, DIST_TYPE & distance)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> end0_coord(dimension);

    grid.ComputeCoord(endpoint0, end0_coord.Ptr());

    // Initialize
    distance = 0;

    bool is_distance_set = false;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != edge_dir) {
        DIST_TYPE distance2 = end0_coord[d]-pcoord[d]/grid.Spacing(d);
        if (distance2 < 0) { distance2 = -distance2; }
        if (!is_distance_set || distance2 < distance) {
          distance = distance2; 
          is_distance_set = true;
        }
      }
    }
  }

  /*!
   *  @brief Compute the distance from pcoord[] to the closest and farthest grid facets incident on edge.
   *  @tparam GRID_TYPE Must have function SpacingPtrConst().
   *  @param grid Regular grid.
   */
  template <typename GRID_TYPE, typename CTYPE, 
            typename VTYPE, typename DIR_TYPE, typename DIST_TYPE>
  void compute_unscaled_distances_to_facets
  (const GRID_TYPE & grid, const CTYPE pcoord[],
   const VTYPE endpoint0, const DIR_TYPE edge_dir, 
   DIST_TYPE & closest_distance, DIST_TYPE & farthest_distance)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> end0_coord(dimension);

    grid.ComputeCoord(endpoint0, end0_coord.Ptr());

    // Initialize
    closest_distance = 0;
    farthest_distance = 1;

    bool is_distance_set = false;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != edge_dir) {
        DIST_TYPE distance2 = end0_coord[d]-pcoord[d]/grid.Spacing(d);
        if (distance2 < 0) { distance2 = -distance2; }
        if (!is_distance_set) {
          closest_distance = distance2;
          farthest_distance = distance2;
        }
        else {
          if (distance2 < closest_distance) 
            { closest_distance = distance2; }
          if (distance2 > farthest_distance) 
            { farthest_distance = distance2; }
        }
        is_distance_set = true;
      }
    }
  }


  /*!
   *  @brief Compute the minimum distance from the four quad vertices to the closest grid facet incident on edge (endpoint0, endpoint1).
   *  @tparam GRID_TYPE Must have function SpacingPtrConst().
   *  @param grid Regular grid.
   */
  template <typename GRID_TYPE, typename CTYPE, 
            typename VTYPE0, typename VTYPE1, typename DIR_TYPE,
            typename DIST_TYPE>
  void compute_unscaled_min_distance_to_closest_facet
  (const GRID_TYPE & grid, const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const VTYPE1 endpoint0, const DIR_TYPE edge_dir, DIST_TYPE & min_distance)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE dimension = grid.Dimension();
    const NTYPE NUM_VERT_PER_QUAD(4);

    const CTYPE * pcoord0 = vertex_coord+dimension*quad_vert[0];
    compute_unscaled_distance_to_closest_facet
      (grid, pcoord0, endpoint0, edge_dir, min_distance);

    for (NTYPE i = 1; i < NUM_VERT_PER_QUAD; i++) {
      const CTYPE * pcoord = vertex_coord+dimension*quad_vert[i];

      DIST_TYPE distance;
      compute_unscaled_distance_to_closest_facet
        (grid, pcoord, endpoint0, edge_dir, distance);

      if (distance < min_distance) 
        { min_distance = distance; }
    }
  }


  /*!
   *  @brief Compute the minimum distance from the four quad vertices 
   *    to the closest and farthest grid facets incident on an edge.
   *  @tparam GRID_TYPE Must have function SpacingPtrConst().
   *  @param grid Regular grid.
   *  @param[out] min_distance Minimum distance from any quad vert
   *     to the facets incident on the edge.
   *  @param[out] min_max_distance Minimum of the maximum distance 
   *     from a quad vert to the facets incident on the edge.
   */
  template <typename GRID_TYPE, typename CTYPE, 
            typename VTYPE0, typename VTYPE1, typename DIR_TYPE,
            typename DIST_TYPE>
  void compute_unscaled_distances_from_quad_vert_to_facets
  (const GRID_TYPE & grid, const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const VTYPE1 endpoint0, const DIR_TYPE edge_dir, 
   DIST_TYPE & min_distance,
   DIST_TYPE & min_max_distance)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE dimension = grid.Dimension();
    const NTYPE NUM_VERT_PER_QUAD(4);

    const CTYPE * pcoord0 = vertex_coord+dimension*quad_vert[0];
    compute_unscaled_distances_to_facets
      (grid, pcoord0, endpoint0, edge_dir, min_distance, min_max_distance);

    for (NTYPE i = 1; i < NUM_VERT_PER_QUAD; i++) {
      const CTYPE * pcoord = vertex_coord+dimension*quad_vert[i];

      DIST_TYPE closest_distance, farthest_distance;
      compute_unscaled_distances_to_facets
        (grid, pcoord, endpoint0, edge_dir, 
         closest_distance, farthest_distance);

      if (closest_distance < min_distance)
        { min_distance = closest_distance; }
      if (farthest_distance < min_max_distance)
        { min_max_distance = farthest_distance; }
    }
  }


  /*!
   *  @brief Return true if some quad vertex is close to one of the facets incident on edge (endpoint0, endpoint1).
   *  @tparam GRID_TYPE Must have function SpacingPtrConst()
   *    and EdgeDirecton().
   *  @param grid Regular grid.
   */
  template <typename GRID_TYPE, typename CTYPE, 
            typename VTYPE0, typename VTYPE1, 
            typename DIR_TYPE,
            typename DIST_TYPE>
  bool is_some_quad_vert_close_to_facet
  (const GRID_TYPE & grid, const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[], const VTYPE1 endpoint0, const DIR_TYPE edge_dir,
   const DIST_TYPE min_distance)
  {
    CTYPE distance;
    compute_unscaled_min_distance_to_closest_facet
      (grid, vertex_coord, quad_vert, endpoint0, edge_dir, distance);

    if (distance < min_distance) { return(true); }
    else { return(false); }
  }

  ///@}


  // **************************************************
  /// @name TRIANGULATE POLYGON
  // **************************************************

  ///@{

  /*!
   *  @brief Triangulate a single polygon splitting maximum polygon angle.
   *  - Add diagonals from the polygon vertex with the maximum angle.
   *  - Polygon vertices are listed in clockwise or counter-clockwise order
   *    around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  - If some polygon edge has length (near) zero, the triangulation
   *    will be from an arbitrary vertex.
   */
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_polygon_split_max_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const NTYPE num_poly_vert,
   const VTYPE0 * poly_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    CTYPE min_cos;
    bool flag_zero;

    if (num_poly_vert < 3) { return; };

    NTYPE k_min_cos = 0;        // Index of the vertex with min_cos
    NTYPE iv0 = poly_vert[num_poly_vert-1];
    NTYPE iv1 = poly_vert[0];
    NTYPE iv2 = poly_vert[1];

    IJK::compute_cos_triangle_angle_coord_list
      (dimension, vertex_coord, iv0, iv1, iv2, max_small_magnitude,
       min_cos, flag_zero);

    for (NTYPE i = 1; i < num_poly_vert; i++) {
      iv0 = poly_vert[i-1];
      iv1 = poly_vert[i];
      iv2 = poly_vert[((i+1)%num_poly_vert)];
      CTYPE cos_v0;

      IJK::compute_cos_triangle_angle_coord_list
        (dimension, vertex_coord, iv0, iv1, iv2, max_small_magnitude,
         cos_v0, flag_zero);

      if (cos_v0 < min_cos) {
        k_min_cos = i;
        min_cos = cos_v0;
      }
    }

    IJK::triangulate_polygon_from_vertex(num_poly_vert, poly_vert, k_min_cos, tri_vert);
  }


  /// @brief Triangulate list of polygons, splitting max angles.
  /// - Add diagonals from the polygon vertex with the maximum angle.
  template <typename DTYPE, typename CTYPE,
            typename NTYPE0, typename NTYPE1, 
            typename VTYPE0, typename VTYPE1,
            typename MTYPE, typename ITYPE>
  void triangulate_polygon_list_split_max_angle
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const NTYPE0 * num_poly_vert, const VTYPE0 * poly_vert,
   const ITYPE * first_poly_vert, const NTYPE1 num_poly,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    for (NTYPE1 ipoly = 0; ipoly < num_poly; ipoly++) {
      triangulate_polygon_split_max_angle
        (dimension, vertex_coord, num_poly_vert[ipoly], 
         poly_vert+first_poly_vert[ipoly], max_small_magnitude, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate list of polygons, splitting max angles.
   *  - C++ STL vector format for vertex_coord[], num_poly_vert[],
   *    poly_vert[], and first_poly_vert[].
   */
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE, typename ITYPE>
  void triangulate_polygon_list_split_max_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const std::vector<NTYPE> & num_poly_vert, 
   const std::vector<VTYPE0> & poly_vert,
   const std::vector<ITYPE> & first_poly_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    triangulate_polygon_list_split_max_angle
      (dimension, IJK::vector2pointer(vertex_coord), 
       IJK::vector2pointer(num_poly_vert), IJK::vector2pointer(poly_vert),
       IJK::vector2pointer(first_poly_vert), num_poly_vert.size(),
       max_small_magnitude, tri_vert);
  }

  ///@}


  // **************************************************
  /// @name TRIANGULATE QUADRILATERAL
  // **************************************************

  ///@{

  /// Triangulate a list of quadrilaterals, splitting max angles.
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_split_max_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const NTYPE num_quad,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_QUAD (4);

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      NTYPE k = iquad*NUM_VERT_PER_QUAD;
      triangulate_polygon_split_max_angle
        (dimension, vertex_coord, NUM_VERT_PER_QUAD, quad_vert+k, 
         max_small_magnitude, tri_vert);
    }
  }


  /// @brief Triangulate a list of quadrilaterals, splitting max angles.
  /// - C++ STL vector format for vertex_coord and quad_vert.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_split_max_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE0> & quad_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);

    SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    triangulate_quad_split_max_angle
      (dimension, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert), num_quad, max_small_magnitude,
       tri_vert);

  }


  /*!
   *  @brief Triangulate a single quad, maximizing minimum triangle angle.
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    const CTYPE * vcoord0 = vertex_coord+quad_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+quad_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+quad_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+quad_vert[3]*dimension;
    CTYPE cos_min_angle;
    bool flag_diag02;
    bool flag_zero;

    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
       cos_min_angle, flag_diag02, flag_zero);

    triangulate_quad_using_diagonal(quad_vert, flag_diag02, tri_vert);
  }


  /// Triangulate a list of quadrilaterals, maximizing minimum triangle angle.
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE>
  void triangulate_quad_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const NTYPE num_quad,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_QUAD(4);

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      NTYPE k = iquad*NUM_VERT_PER_QUAD;
      triangulate_quad_max_min_angle
        (dimension, vertex_coord, quad_vert+k, max_small_magnitude, tri_vert);
    }
  }


  /// @brief Triangulate a list of quadrilaterals, maximizing minimum triangle angle.
  /// - C++ STL vector format for vertex_coord and quad_vert.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_max_min_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE0> & quad_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    triangulate_quad_max_min_angle
      (dimension, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert), num_quad, max_small_magnitude,
       tri_vert);
  }


  /*!
   *  @brief Triangulate polygon, maximizing minimum angle.
   *  - Max min angle over all triantulations with triangles incident
   *    on one vertex.
   *  - Add new triangles to vector tri_vert.
   *  @param[out] index_tri_vert Index of triangulation vertex.
   *     All triangulation triangles are incident on poly_vert[index_tri_vert].
   */
  template <typename DTYPE, typename NTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE,
            typename COS_TYPE, typename IVTYPE>
  void triangulate_polygon_max_min_angle
  (const DTYPE dimension, const CTYPE * vertex_coord, 
   const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert,
   IVTYPE & index_tri_vert,
   COS_TYPE & cos_min_angle)
  {
    bool flagA;

    compute_vertex_to_max_min_triangulation_angle
      (dimension, vertex_coord, num_poly_vert, poly_vert,
       max_small_magnitude, index_tri_vert, cos_min_angle, flagA);

    if (flagA) {
      index_tri_vert = 0;
      cos_min_angle = 0;
    }

    triangulate_polygon_from_vertex(num_poly_vert, poly_vert, index_tri_vert, tri_vert);
  }


  /*! 
   *  @brief Triangulate polygon, maximizing minimum angle.
   *  - Max min angle over all triantulations with triangles incident
   *    on one vertex.
   *  - Version that does not return index_tri_vert or cos_min_angle.
   */
  template <typename DTYPE, typename NTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_polygon_max_min_angle
  (const DTYPE dimension, const CTYPE * vertex_coord, 
   const NTYPE num_poly_vert, const VTYPE0 * poly_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    NTYPE index_tri_vert;
    CTYPE cos_min_angle;

    triangulate_polygon_max_min_angle
      (dimension, vertex_coord, num_poly_vert, poly_vert,
       max_small_magnitude, tri_vert, index_tri_vert, cos_min_angle);
  }


  /*!
   *  @brief Triangulate polygon, maximizing minimum angle.
   *  - Max min angle over all triantulations with triangles incident
   *    on one vertex.
   *  - Version that does not return index_tri_vert or cos_min_angle.
   *  - C++ STL vector format for vertex_coord and poly_vert.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_polygon_max_min_angle
  (const DTYPE dimension, 
   const std::vector<CTYPE> & vertex_coord, 
   const std::vector<VTYPE0> & poly_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    triangulate_polygon_max_min_angle
      (dimension, IJK::vector2pointer(vertex_coord),
       poly_vert.size(), IJK::vector2pointer(poly_vert),
       max_small_magnitude, tri_vert);
  }


  /// Return true if vertex iw separates iv0 and iv1 from iv2 and iv3 in direction d.
  template <typename CTYPE, typename VTYPE0, typename VTYPE1, 
            typename DIR_TYPE>
  bool does_vertex_separate
  (const CTYPE vertex_coord[],
   const VTYPE0 iw, const VTYPE1 iv0, const VTYPE1 iv1,
   const VTYPE1 iv2, const VTYPE1 iv3, const DIR_TYPE d)
  {
    const DIR_TYPE DIM3(3);
    const CTYPE wc = vertex_coord[DIM3*iw+d];
    const CTYPE vc0 = vertex_coord[DIM3*iv0+d];
    const CTYPE vc1 = vertex_coord[DIM3*iv1+d];
    const CTYPE vc2 = vertex_coord[DIM3*iv2+d];
    const CTYPE vc3 = vertex_coord[DIM3*iv3+d];

    if (vc0 <= wc && vc1 <= wc && vc2 >= wc && vc3 >= wc)
      { return(true); }

    if (vc0 >= wc && vc1 >= wc && vc2 <= wc && vc3 <= wc)
      { return(true); }

    return(false); 
  }


  /// Return true if quad vertices surround the dual edge.
  template <typename CTYPE, typename VTYPE0, typename VTYPE1,
            typename ISOPOLY_INFO_TYPE>
  bool does_quad_surround_dual_edge
  (const CTYPE vertex_coord[], 
   const VTYPE0 quad_vert[],
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const VTYPE1 vertex_on_dual_edge)
  {
    typedef typename ISOPOLY_INFO_TYPE::DIRECTION_TYPE DIR_TYPE;

    const DIR_TYPE DIM3(3);
    const DIR_TYPE edge_dir = isopoly_info.GridEdgeDirection();

    const DIR_TYPE d1 = (edge_dir+1)%DIM3;
    const DIR_TYPE d2 = (edge_dir+2)%DIM3;

    if (does_vertex_separate
        (vertex_coord, vertex_on_dual_edge, 
         quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3], d1)) {

      if (does_vertex_separate
          (vertex_coord, vertex_on_dual_edge, 
           quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[0], d2)) {

        return(true);
      }
    }

    if (does_vertex_separate
        (vertex_coord, vertex_on_dual_edge, 
         quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3], d2)) {

      if (does_vertex_separate
          (vertex_coord, vertex_on_dual_edge, 
           quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[0], d1)) {

        return(true);
      }
    }

    return(false);
  }


  /*!
   *  @brief Triangulate a set of quadrilaterals.
   *  - Add 4 triangles when all quad vertices are min distance from facets
   *    incident on dual quad edge.
   *    <br> Otherwise, use diagonal which minimizes quad angle.
   *  - Split quad i into four triangles by adding a triangle between
   *    each quad edge and vertex (first_vertex+i).
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise 
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ISOPOLY_INFO_TYPE, typename MTYPE, typename DIST_TYPE>
  void triangulate_quad_tri4_by_distance
  (const GRID_TYPE & grid, const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[], const NTYPE num_quad,
   const ISOPOLY_INFO_TYPE isopoly_info[], const VTYPE1 first_vertex, 
   const DIST_TYPE min_distance,  const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    const NTYPE NUM_VERT_PER_QUAD(4);

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {

      const NTYPE k = iquad*NUM_VERT_PER_QUAD;
      const VTYPE0 iv0 = isopoly_info[iquad].GridEdgeEndpoint0();
      const DTYPE edge_direction = 
	isopoly_info[iquad].GridEdgeDirection();

      if (is_some_quad_vert_close_to_facet
          (grid, vertex_coord, quad_vert+k, iv0, 
           edge_direction, min_distance)) {

        // Split by some diagonal into two triangles.
        triangulate_quad_max_min_angle
          (dimension, vertex_coord, quad_vert+k, max_small_magnitude,
           tri_vert);
      }
      else {
        // Split into four triangles.
        triangulate_polygon_from_interior_vertex
          (NUM_VERT_PER_QUAD, quad_vert+k, first_vertex+iquad, tri_vert);
      }

    }
  }


  /*!
   *  @brief Triangulate a set of quadrilaterals.
   *  - Add 4 triangles when all quad vertices are min distance from facets
   *    incident on dual quad edge.
   *    <br> Otherwise, use diagonal which minimizes quad angle.
   *  - Split quad i into four triangles by adding a triangle between
   *    each quad edge and vertex (first_vertex+i).
   *  - Quadrilateral vertices are listed in clockwise or counter-clockwise
   *    order around the polygon.
   *  - Add new triangles to vector tri_vert.
   *  - C++ STL vector format for vertex_coord and quad_vert and dual_edge[].
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename GRID_EDGE_TYPE, typename MTYPE, typename DIST_TYPE>
  void triangulate_quad_tri4_by_distance
  (const GRID_TYPE & grid, const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE0> & quad_vert, 
   const std::vector<GRID_EDGE_TYPE> & dual_edge, 
   const VTYPE1 first_vertex, 
   const DIST_TYPE min_distance,  const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    triangulate_quad_tri4_by_distance
      (grid, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert), num_quad,
       IJK::vector2pointer(dual_edge), first_vertex,
       min_distance, max_small_magnitude, tri_vert);
  }


  /*!
   *  @brief Triangulate a single quad into two or four triangles.
   *  - Minimize maximum triangle angle.
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_quad_tri4_max_min_angle
  (const DTYPE dimension,
   const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const VTYPE1 ivert_on_dual_edge,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const CTYPE * vcoord0 = vertex_coord+quad_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+quad_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+quad_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+quad_vert[3]*dimension;
    const CTYPE * vcoordX = vertex_coord+ivert_on_dual_edge*dimension;
    CTYPE cos_min_angle;
    bool flag_tri4, flag_diag02, flag_zero;

    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoordX,
       max_small_magnitude, cos_min_angle, flag_tri4, flag_diag02, flag_zero);

    if (flag_tri4) {
      triangulate_polygon_from_interior_vertex
        (NUM_VERT_PER_QUAD, quad_vert, ivert_on_dual_edge, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal(quad_vert, flag_diag02, tri_vert);
    }
  }


  /*! 
   *  @brief Triangulate a list of isosurface quadrilaterals.
   *  - Each quadrilateral is dual to a grid edge.
   *  @param min_dist_allow_tri4 Allow triangulation into 4 triangles when all
   *    quad vertices are min_distance_allow_tri4 distance from facets
   *    incident on dual quad edges.
   *    - If some quad vertex is closer than min_distance_allow_tri4 to
   *    a facet incident on the dual quad edge, use a diagonal triangulation
   *    into two triangles.
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ISOPOLY_INFO_TYPE,
            typename DIST_TYPE, typename MTYPE>
  void triangulate_dual_quad_tri4_max_min_angle
  (const GRID_TYPE & grid, 
   const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const NTYPE num_quad,
   const ISOPOLY_INFO_TYPE isopoly_info[],
   const VTYPE1 first_vertex_on_dual_edge, 
   const DIST_TYPE min_distance_allow_tri4,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    const DTYPE DIM3(3);
    const NTYPE NUM_VERT_PER_QUAD(4);
    CTYPE min_distance, min_max_distance;
    IJK::PROCEDURE_ERROR error("triangulate_dual_quad_tri4_max_min_angle");

    if (num_quad > 0 && isopoly_info == NULL) {
      error.AddMessage
        ("Programming error.  Array isopoly_info[] is empty.");
      throw error;
    }

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      const VTYPE2 iw = first_vertex_on_dual_edge+iquad;
      const NTYPE k = iquad*NUM_VERT_PER_QUAD;

      if (dimension == DIM3 &&
          !does_quad_surround_dual_edge
          (vertex_coord, quad_vert+k, isopoly_info[iquad], iw)) {

        // Split by some diagonal into two triangles.
        triangulate_quad_max_min_angle
          (dimension, vertex_coord, quad_vert+k, max_small_magnitude,
           tri_vert);
      }
      else {
	const VTYPE0 iv0 = isopoly_info[iquad].GridEdgeEndpoint0();
	const DTYPE edge_direction = 
	  isopoly_info[iquad].GridEdgeDirection();

        compute_unscaled_distances_from_quad_vert_to_facets
          (grid, vertex_coord, quad_vert+k, iv0, edge_direction, 
           min_distance, min_max_distance);

        if (min_max_distance > min_distance_allow_tri4) {
          // Split into two or four triangles, whichever maximizes min angle.
          triangulate_quad_tri4_max_min_angle
            (dimension, vertex_coord, quad_vert+k, iw, max_small_magnitude, 
             tri_vert);
        }
        else {
          // Split by some diagonal into two triangles.
          triangulate_quad_max_min_angle
            (dimension, vertex_coord, quad_vert+k, max_small_magnitude,
             tri_vert);
        }
      }
    }
  }


  /*!
   *  @brief Triangulate a list of isosurface quadrilaterals.
   *  - Each quadrilateral is dual to a grid edge.
   *  - C++ STL vector format for vertex_coord and quad_vert.
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2, 
            typename GRID_EDGE_TYPE,
            typename DIST_TYPE, typename MTYPE>
  void triangulate_dual_quad_tri4_max_min_angle
  (const GRID_TYPE & grid, 
   const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE0> & quad_vert,
   const std::vector<GRID_EDGE_TYPE> & dual_edge,
   const VTYPE1 first_vertex_on_dual_edge, 
   const DIST_TYPE min_distance_allow_tri4,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (num_quad != dual_edge.size()) {
      IJK::PROCEDURE_ERROR error("triangulate_dual_quad_tri4_max_min_angle");
      error.AddMessage
        ("Programming error. Number of dual edges does not equal number of isosurface quadrilaterals.");
      error.AddMessage("  Number of dual edges: ", dual_edge.size(), ".");
      error.AddMessage("  Number of isosurface quadrilaterals: ", 
                       num_quad, ".");
      throw error;
    }

    triangulate_dual_quad_tri4_max_min_angle
      (grid, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert), num_quad, 
       IJK::vector2pointer(dual_edge),
       first_vertex_on_dual_edge,
       min_distance_allow_tri4, max_small_magnitude,
       tri_vert);
  }


  /// Triangulate a list of isosurface quadrilaterals.
  template <typename GRID_TYPE, typename CTYPE, typename NTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_quad_tri4_max_min_angle
  (const GRID_TYPE & grid, 
   const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const NTYPE num_quad,
   const VTYPE1 first_vertex_on_dual_edge, 
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    const NTYPE NUM_VERT_PER_QUAD(4);
    IJK::PROCEDURE_ERROR error("triangulate_quad_tri4_max_min_angle");

    if (num_quad > 0 && quad_vert == NULL) {
      error.AddMessage
        ("Programming error.  Array quad_vert[] is empty.");
      throw error;
    }

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      const VTYPE2 iw = first_vertex_on_dual_edge+iquad;
      const NTYPE k = iquad*NUM_VERT_PER_QUAD;

      // Split into two or four triangles, whichever maximizes min angle.
      triangulate_quad_tri4_max_min_angle
        (dimension, vertex_coord, quad_vert+k, iw, max_small_magnitude, 
         tri_vert);
    }
  }


  /// @brief Triangulate a list of isosurface quadrilaterals.
  /// - C++ STL vector format for vertex_coord and quad_vert.
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2, 
            typename MTYPE>
  void triangulate_quad_tri4_max_min_angle
  (const GRID_TYPE & grid, 
   const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE0> & quad_vert,
   const VTYPE1 first_additional_vertex,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    triangulate_quad_tri4_max_min_angle
      (grid, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert), num_quad, 
       first_additional_vertex,
       max_small_magnitude,
       tri_vert);
  }


  /*!
   *  @brief Triangulate a single quad into four triangles or two triangles using diag 02.
   *  - Minimize maximum triangle angle.
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_quad_tri4_or_tri02_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivert_on_dual_edge,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const CTYPE * vcoord0 = vertex_coord+quad_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+quad_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+quad_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+quad_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord+ivert_on_dual_edge*dimension;
    CTYPE cos_min_diag02, cos_min_tri4;
    bool flag_zero_diag02, flag_zero_tri4;

    compute_cos_min_quad_tri02_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
       cos_min_diag02, flag_zero_diag02);
    compute_cos_min_quad_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4,
       max_small_magnitude, cos_min_tri4, flag_zero_tri4);

    bool flag_tri4;

    if (flag_zero_tri4) { 
      flag_tri4 = false;
    }
    else if (flag_zero_diag02) {
      flag_tri4 = true;
    }
    else {
      if (cos_min_diag02 < cos_min_tri4) 
        { flag_tri4 = false; }
      else
        { flag_tri4 = true; }
    }

    if (flag_tri4) {
      triangulate_polygon_from_interior_vertex
        (NUM_VERT_PER_QUAD, quad_vert, ivert_on_dual_edge, tri_vert);
    }
    else {
      triangulate_polygon_from_vertex(NUM_VERT_PER_QUAD, quad_vert, 0, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a single quad into four triangles or two triangles
   *    using diag 02.
   *  - Minimize maximum triangle angle.
   *  - C++ STL vector format for vertex_coord.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_quad_tri4_or_tri02_max_min_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivert_on_dual_edge,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    triangulate_quad_tri4_or_tri02_max_min_angle
      (dimension, IJK::vector2pointer(vertex_coord), quad_vert,
       ivert_on_dual_edge, max_small_magnitude, tri_vert);
  }


  /*!
   *  @brief Triangulate a single quad into four triangles or two triangles using diag 13.
   *  - Minimize maximum triangle angle.
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_quad_tri4_or_tri13_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivert_on_dual_edge,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const CTYPE * vcoord0 = vertex_coord+quad_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+quad_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+quad_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+quad_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord+ivert_on_dual_edge*dimension;
    CTYPE cos_min_diag13, cos_min_tri4;
    bool flag_zero_diag13, flag_zero_tri4;

    compute_cos_min_quad_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, max_small_magnitude, 
       cos_min_diag13, flag_zero_diag13);
    compute_cos_min_quad_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4,
       max_small_magnitude, cos_min_tri4, flag_zero_tri4);

    bool flag_tri4;


    if (flag_zero_tri4) { 
      flag_tri4 = false;
    }
    else if (flag_zero_diag13) {
      flag_tri4 = true;
    }
    else {
      if (cos_min_diag13 < cos_min_tri4) 
        { flag_tri4 = false; }
      else
        { flag_tri4 = true; }
    }

    if (flag_tri4) {
      triangulate_polygon_from_interior_vertex
        (NUM_VERT_PER_QUAD, quad_vert, ivert_on_dual_edge, tri_vert);
    }
    else {
      triangulate_polygon_from_vertex(NUM_VERT_PER_QUAD, quad_vert, 1, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a single quad into four triangles or two triangles
   *    using diag 13.
   *  - Minimize maximum triangle angle.
   *  - C++ STL vector format for vertex_coord.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_quad_tri4_or_tri13_max_min_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivert_on_dual_edge,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    triangulate_quad_tri4_or_tri13_max_min_angle
      (dimension, IJK::vector2pointer(vertex_coord), quad_vert,
       ivert_on_dual_edge, max_small_magnitude, tri_vert);
  }

  ///@}


  // **************************************************
  /// @name TRIANGULATE PENTAGON
  // **************************************************

  ///@{

  /*!
   *  @brief Triangulate a single pentagon into five triangles.
   *  - Add new coord vcoordX to vertex_coord[].
   *  - Pentagon vertices are listed in clockwise or counter-clockwise order
   *    around the pentagon.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename CTYPE0, typename CTYPE1>
 void triangulate_pentagon_tri5
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   const CTYPE0 vcoordX[],
   const bool flag_reverse_orient,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<VTYPE5> & tri_vert)
  {
    const VTYPE5 iv_new = IJK::insert_coord(dimension, vcoordX, vertex_coord);
    IJK::triangulate_pentagon_with_vertex
      (v0, v1, v2, v3, v4, iv_new, flag_reverse_orient, tri_vert);
  }


  /*!
   *  @brief Triangulate a single pentagon into five triangles.
   *  - Version where pentagon vertices are in an array.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE, typename CTYPE0, typename CTYPE1,
            typename TRI_VTYPE>
 void triangulate_pentagon_tri5
  (const DTYPE dimension,
   const VTYPE pentagon_vert[],
   const CTYPE0 vcoordX[],
   const bool flag_reverse_orient,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    triangulate_pentagon_tri5
      (dimension, pentagon_vert[0], pentagon_vert[1], pentagon_vert[2],
       pentagon_vert[3], pentagon_vert[4], vcoordX, flag_reverse_orient,
       vertex_coord, tri_vert);
  }


  /*!
   *  @brief Triangulate a single pentagon into three or five triangles based on flag_tri5.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Pentagon vertices are listed in clockwise or counter-clockwise order
   *    around the pentagon.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   *  @param iv_index Index of vertex to be used in triangulation
   *    into three triangles.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename NTYPE,
            typename CTYPE0, typename CTYPE1>
 void triangulate_pentagon_tri3_or_tri5
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   const CTYPE0 vcoordX[],
   const bool flag_tri5,
   const NTYPE iv_index,
   const bool flag_reverse_orient,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<VTYPE5> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_PENTAGON(5);
    const VTYPE0 pentagon_vert[NUM_VERT_PER_PENTAGON] = 
      { v0, v1, v2, v3, v4 };

    if (flag_tri5) {
      triangulate_pentagon_tri5
        (dimension, v0, v1, v2, v3, v4, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      triangulate_pentagon
        (pentagon_vert, iv_index, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a single pentagon into three or five triangles.
   *  - Version using data structure POLY_TRIANGULATION_RESULT.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename CTYPE0, typename CTYPE1,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
 void triangulate_pentagon_tri3_or_tri5
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   const CTYPE0 vcoordX[],
   const bool flag_reverse_orient,
   const IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & poly_tri,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<VTYPE5> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_PENTAGON(5);
    const VTYPE0 pentagon_vert[NUM_VERT_PER_PENTAGON] = 
      { v0, v1, v2, v3, v4 };

    if (poly_tri.num_triangles == 5) {
      triangulate_pentagon_tri5
        (dimension, v0, v1, v2, v3, v4, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      IJK::triangulate_pentagon
        (pentagon_vert, poly_tri.tri_vertex_index, flag_reverse_orient, 
         tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a single pentagon into three or five triangles based on flag_tri5.
   *  - Triangulate from v0 in triangulation into 3 triangles.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Pentagon vertices are listed in clockwise or counter-clockwise order
   *    around the pentagon.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename CTYPE0, typename CTYPE1>
 void triangulate_pentagon_v0tri3_or_tri5
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   const CTYPE0 vcoordX[],
   const bool flag_tri5,
   const bool flag_reverse_orient,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<VTYPE5> & tri_vert)
  {
    if (flag_tri5) {
      triangulate_pentagon_tri5
        (dimension, v0, v1, v2, v3, v4, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      IJK::triangulate_pentagon
        (v0, v1, v2, v3, v4,flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a single pentagon into three or five triangles.
   *  - Version using data structure POLY_TRIANGULATION_RESULT.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename CTYPE0, typename CTYPE1,
            typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
 void triangulate_pentagon_v0tri3_or_tri5
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   const CTYPE0 vcoordX[],
   const IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> &
   pentagon_tri,
   const bool flag_reverse_orient,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<VTYPE5> & tri_vert)
  {
    if (pentagon_tri.num_triangles == 5) {
      triangulate_pentagon_tri5
        (dimension, v0, v1, v2, v3, v4, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      IJK::triangulate_pentagon
        (v0, v1, v2, v3, v4, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a pentagon.
   *  - Maximize minimum triangle angle.
   *  - Pentagon vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename CTYPE, typename VTYPE, typename MTYPE,
            typename TRI_VTYPE>
  void triangulate_pentagon_max_min_angle
  (const CTYPE vertex_coord[],
   const VTYPE pentagon_vert[],
   const MTYPE max_small_magnitude,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const int DIM3(3);
    int iv;
    CTYPE cos_min_angle;
    bool flag_zero;

    compute_vertex_to_max_min_pentagon_triangulation_angle
      (DIM3, vertex_coord, pentagon_vert, max_small_magnitude, iv,
       cos_min_angle, flag_zero);

    if (flag_zero) {
      IJK::triangulate_pentagon(pentagon_vert, 0, tri_vert);
    }
    else {
      IJK::triangulate_pentagon(pentagon_vert, iv, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a pentagon into three or five triangles.
   *  - Maximize minimum triangle angle.
   *  - Pentagon vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  @param ivX Triangulate into five triangles by connecting each
   *      pentagon edge to ivX.
   */
  template <typename DTYPE, typename VTYPE0, typename CTYPE0, typename CTYPE1,
            typename MTYPE, typename TRI_VTYPE>
  void triangulate_pentagon_tri5_max_min_angle
  (const DTYPE dimension,
   const VTYPE0 pentagon_vert[],
   const CTYPE0 vcoordX[],
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const CTYPE1 * vcoord_ptr = &(vertex_coord[0]);
    const CTYPE1 * vcoord0 = vcoord_ptr + pentagon_vert[0]*dimension;
    const CTYPE1 * vcoord1 = vcoord_ptr + pentagon_vert[1]*dimension;
    const CTYPE1 * vcoord2 = vcoord_ptr + pentagon_vert[2]*dimension;
    const CTYPE1 * vcoord3 = vcoord_ptr + pentagon_vert[3]*dimension;
    const CTYPE1 * vcoord4 = vcoord_ptr + pentagon_vert[4]*dimension;
    int iv_max_min;
    CTYPE1 cos_min_angle;
    bool flag_zero;
    bool flag_tri5;

    compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoordX, 
       max_small_magnitude, cos_min_angle, flag_tri5, iv_max_min, flag_zero);

    if (flag_zero) {
      IJK::triangulate_pentagon(pentagon_vert, 0, flag_reverse_orient, tri_vert);
    }
    else if (flag_tri5) {
      triangulate_pentagon_tri5
        (dimension, pentagon_vert, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      IJK::triangulate_pentagon
        (pentagon_vert, iv_max_min, flag_reverse_orient, tri_vert);
    }

  }


  template <typename CTYPE, typename VTYPE, typename MTYPE,
            typename TRI_VTYPE>
  void triangulate_pentagon_max_min_angle
  (const std::vector<CTYPE> & vertex_coord,
   const VTYPE pentagon_vert[],
   const MTYPE max_small_magnitude,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    triangulate_pentagon_max_min_angle
      (IJK::vector2pointer(vertex_coord), pentagon_vert, 
       max_small_magnitude, tri_vert);
  }


  template <typename CTYPE, typename VTYPE, typename NTYPE,
            typename MTYPE, typename TRI_VTYPE>
  void triangulate_pentagon_list_max_min_angle
  (const std::vector<CTYPE> & vertex_coord,
   const VTYPE pentagon_vert[],
   const NTYPE num_pentagon,
   const MTYPE max_small_magnitude,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const int NUM_VERT_PER_PENTAGON(5);

    for (int i = 0; i < num_pentagon; i++) {
      triangulate_pentagon_max_min_angle
        (vertex_coord, pentagon_vert+i*NUM_VERT_PER_PENTAGON, 
         max_small_magnitude, tri_vert);
    }
  }


  template <typename CTYPE, typename VTYPE,
            typename MTYPE, typename TRI_VTYPE>
  void triangulate_pentagon_list_max_min_angle
  (const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE> & pentagon_vert,
   const MTYPE max_small_magnitude,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    const int num_pentagon = pentagon_vert.size()/NUM_VERT_PER_PENTAGON;

    triangulate_pentagon_list_max_min_angle
      (vertex_coord, IJK::vector2pointer(pentagon_vert), num_pentagon,
       max_small_magnitude, tri_vert);
  }


  /*!
   *  @brief Triangulate list of pentagons into three or five triangles each.
   *  - Add vertex at centroid of pentagon for triangulations
   *    into 5 triangles.
   */
  template <typename DTYPE, typename VTYPE, typename NTYPE, typename MTYPE, 
            typename CTYPE, typename TRI_VTYPE>
  void triangulate_pentagon_list_tri5_centroid_max_min_angle
  (const DTYPE dimension,
   const VTYPE pentagon_vert[],
   const NTYPE num_pentagon,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    IJK::ARRAY<CTYPE> pentagon_centroid(dimension);

    for (int i = 0; i < num_pentagon; i++) {

      const VTYPE * first_pentagon_vert =  
        pentagon_vert + i*NUM_VERT_PER_PENTAGON;

      IJK::compute_pentagon_centroid
        (dimension, first_pentagon_vert, vertex_coord,
         pentagon_centroid.Ptr());

      triangulate_pentagon_tri5_max_min_angle
        (dimension, pentagon_vert+i*NUM_VERT_PER_PENTAGON,
         pentagon_centroid.PtrConst(), false,
         max_small_magnitude, vertex_coord, tri_vert);
    }
  }


  /// @brief Triangulate list of pentagons into three or five triangles each.
  /// - Version using C++ STL vector to store array pentagon_vert[].
  template <typename DTYPE, typename VTYPE, typename CTYPE, 
            typename MTYPE, typename TRI_VTYPE>
  void triangulate_pentagon_list_tri5_centroid_max_min_angle
  (const DTYPE dimension,
   const std::vector<VTYPE> & pentagon_vert,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    const int num_pentagon = pentagon_vert.size()/NUM_VERT_PER_PENTAGON;

    triangulate_pentagon_list_tri5_centroid_max_min_angle
      (dimension, IJK::vector2pointer(pentagon_vert), num_pentagon,
       max_small_magnitude, vertex_coord, tri_vert);
  }

  ///@}


  // ***************************************************************
  /// @name TRIANGULATE HEXAGON
  // ***************************************************************

  ///@{

  /*!
   *  @brief Triangulate hexagon.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Split hexagon with edges (v0,v3) or (v1,v4) into two quadrilaterals
   *    and triangulate.
   *  @param[in] hex_vert[] Six hexagon vertices.
   *  - Six vertices are listed in clockwise or counter-clockwise order
   *      around the hexagon.
   *  @param[in] ivertX Index of possible vertex to use to triangulate
   *    hexagon.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_hexagon_max_min_angle_split03_or14
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * hex_vert,
   const VTYPE1 ivertX,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_HEX(6);
    const CTYPE * vcoord0 = vertex_coord+hex_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+hex_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+hex_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+hex_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord+hex_vert[4]*dimension;
    const CTYPE * vcoord5 = vertex_coord+hex_vert[5]*dimension;
    const CTYPE * vcoordX = vertex_coord+ivertX*dimension;
    const VTYPE0 quadA_vert[NUM_VERT_PER_QUAD] = 
      { hex_vert[0], hex_vert[1], hex_vert[4], hex_vert[5] };
    const VTYPE0 * quadB_vert = &(hex_vert[1]);
    const VTYPE0 quadC_vert[NUM_VERT_PER_QUAD] = 
      { hex_vert[5], hex_vert[0], hex_vert[3], hex_vert[4] };
    const VTYPE0 * quadD_vert = &(hex_vert[0]);
    IJK::POLY_TRIANGULATION_RESULT<16,CTYPE,int> quadA_tri_result;
    IJK::POLY_TRIANGULATION_RESULT<16,CTYPE,int> quadB_tri_result;
    IJK::POLY_TRIANGULATION_RESULT<16,CTYPE,int> quadC_tri_result;
    IJK::POLY_TRIANGULATION_RESULT<16,CTYPE,int> quadD_tri_result;
    CTYPE cos_min_hex;
    bool flag_zero_hex;

    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vertex_coord, quadA_vert, max_small_magnitude,
       quadA_tri_result);
    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vertex_coord, quadB_vert, max_small_magnitude,
       quadB_tri_result);
    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vertex_coord, quadC_vert, max_small_magnitude,
       quadC_tri_result);
    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vertex_coord, quadD_vert, max_small_magnitude,
       quadD_tri_result);

    const CTYPE cos_max_min_quadsAB = 
      std::max(quadA_tri_result.cos_min_triangulation_angle,
               quadB_tri_result.cos_min_triangulation_angle);
    const CTYPE cos_max_min_quadsCD = 
      std::max(quadC_tri_result.cos_min_triangulation_angle,
               quadD_tri_result.cos_min_triangulation_angle);

    const bool flag_zeroAB = (quadA_tri_result.flag_zero || 
                              quadB_tri_result.flag_zero);
    const bool flag_zeroCD = (quadC_tri_result.flag_zero || 
                              quadD_tri_result.flag_zero);
       

    IJK::compute_cos_min_hexagon_tri6_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoord5, 
       vcoordX, max_small_magnitude, cos_min_hex, flag_zero_hex);

    int index;
    CTYPE cos_min_angle;
    bool flag_zero;
    IJK::select_minIII
      (cos_max_min_quadsAB, flag_zeroAB,
       cos_max_min_quadsCD, flag_zeroCD,
       cos_min_hex, flag_zero_hex,
       cos_min_angle, index, flag_zero);
    
    if (index == 0) {
      // Triangulate quadA and quadB separately using quad diagonals.
      IJK::triangulate_quad_from_vertex
        (quadA_vert, quadA_tri_result.tri_vertex_index, 
         flag_reverse_orient, tri_vert);
      IJK::triangulate_quad_from_vertex
        (quadB_vert, quadB_tri_result.tri_vertex_index, 
         flag_reverse_orient, tri_vert);
    }
    else if (index == 1) {
      // Triangulate quadC and quadD separately using quad diagonals.
      IJK::triangulate_quad_from_vertex
        (quadC_vert, quadC_tri_result.tri_vertex_index, 
         flag_reverse_orient, tri_vert);
      IJK::triangulate_quad_from_vertex
        (quadD_vert, quadD_tri_result.tri_vertex_index, 
         flag_reverse_orient, tri_vert);
    }
    else if (index == 2) {
      // Use triangulation into 6 triangles incident on vcoordX.
      IJK::triangulate_polygon_from_interior_vertex
        (NUM_VERT_PER_HEX, hex_vert, ivertX, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate hexagon.
   *  - Maximize minimum triangle angle.
   *  - C++ STL vector format for vertex_coord.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_hexagon_max_min_angle_split03_or14
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivertX,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    triangulate_hexagon_max_min_angle_split03_or14
      (dimension, IJK::vector2pointer(vertex_coord), quad_vert,
       ivertX, flag_reverse_orient, max_small_magnitude, tri_vert);
  }


  /*!
   *  @brief Triangulate a single hexagon into six triangles.
   *  - Add new coord vcoordX to vertex_coord[].
   *  - Hexagon vertices are listed in clockwise or counter-clockwise order
   *    around the hexagon.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename VTYPE6,
            typename CTYPE0, typename CTYPE1>
 void triangulate_hexagon_tri6
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const CTYPE0 vcoordX[],
   const bool flag_reverse_orient,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<VTYPE6> & tri_vert)
  {
    const VTYPE6 iv_new = IJK::insert_coord(dimension, vcoordX, vertex_coord);
    IJK::triangulate_hexagon_with_vertex
      (v0, v1, v2, v3, v4, v5, iv_new, flag_reverse_orient, tri_vert);
  }


  /*! 
   *  @brief Triangulate a single hexagon into four or six triangles.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Hexagon vertices are listed in clockwise or counter-clockwise order
   *    around the hexagon.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   *  @param iv_index Index of vertex to be used in triangulation
   *    into three triangles.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename TRI_VTYPE, typename NTYPE,
            typename CTYPE0, typename CTYPE1>
 void triangulate_hexagon_tri4_or_tri6
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const CTYPE0 vcoordX[],
   const bool flag_reverse_orient,
   const NTYPE num_triangles,
   const NTYPE iv_index,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_HEXAGON(6);
    const TRI_VTYPE hexagon_vert[NUM_VERT_PER_HEXAGON] = 
      { v0, v1, v2, v3, v4, v5 };

    if (num_triangles == 6) {
      triangulate_hexagon_tri6
        (dimension, v0, v1, v2, v3, v4, v5, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      IJK::triangulate_hexagon
        (hexagon_vert, iv_index, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *   @brief Triangulate a single hexagon into four or six triangles.
   *  - Version using data structure POLY_TRIANGULATION_RESULT.
   *  @param vcoordX[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename TRI_VTYPE,
            typename CTYPE0, typename CTYPE1,
            typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void triangulate_hexagon_tri4_or_tri6
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4, const VTYPE5 v5,
   const CTYPE0 vcoordX[],
   const bool flag_reverse_orient,
   const IJK::POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> &
   poly_tri,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    triangulate_hexagon_tri4_or_tri6
      (dimension, v0, v1, v2, v3, v4, v5, vcoordX, flag_reverse_orient,
       poly_tri.num_triangles, poly_tri.tri_vertex_index,
       vertex_coord, tri_vert);
  }

  ///@}


  // ***************************************************************
  /// @name TRIANGULATE SEPTAGON
  // ***************************************************************

  ///@{

  /// Compute min angle in triangulation of septagon
  ///   with all triangles incident on vertex_coord0.
  /// *** PROBABLY SHOULD BE "const CTYPE vcoord0[], ..."
  template <typename DTYPE, typename CTYPE, typename MTYPE, 
            typename COS_TYPE>
  void compute_cos_min_sept_triangulation_angle
  (const DTYPE dimension,
   const CTYPE vcoord0,
   const CTYPE vcoord1,
   const CTYPE vcoord2,
   const CTYPE vcoord3,
   const CTYPE vcoord4,
   const CTYPE vcoord5,
   const CTYPE vcoord6,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cos012, cos023, cos034, cos045, cos056;
    bool flag_zero_012, flag_zero_023, flag_zero_034;
    bool flag_zero_045, flag_zero_056;

    // Initialize.
    cos_min_angle = 1;
    flag_zero = true;

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoord2, max_small_magnitude,
       cos012, flag_zero_012);
    if (flag_zero_012) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord2, vcoord3, max_small_magnitude,
       cos023, flag_zero_023);
    if (flag_zero_023) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord3, vcoord4, max_small_magnitude,
       cos034, flag_zero_034);
    if (flag_zero_034) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord4, vcoord5, max_small_magnitude,
       cos045, flag_zero_045);
    if (flag_zero_045) { return; }

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord5, vcoord6, max_small_magnitude,
       cos056, flag_zero_056);
    if (flag_zero_056) { return; }

    cos_min_angle = cos012;
    if (cos_min_angle < cos023) { cos_min_angle = cos023; }
    if (cos_min_angle < cos034) { cos_min_angle = cos034; }
    if (cos_min_angle < cos045) { cos_min_angle = cos045; }
    if (cos_min_angle < cos056) { cos_min_angle = cos056; }

    flag_zero = false;
  }

  ///@}



  // **************************************************
  //! @name ENVELOPE
  // **************************************************

  ///@{

  /*!
   *  @brief Return true if diagonal (w0,w2) is in envelope.
   *  @param v0[] Grid vertex v0.
   *  @param v1[] Grid vertex v1. (v0,v1) is a grid edge.
   *  @param w0[] Isosurface vertex w0.
   *  @param w1[] Isosurface vertex w1.
   *  @param w2[] Isosurface vertex w2.
   *  @param w3[] Isosurface vertex w3.
   *         (w0,w1,w2,w3) is an isosurface quadrilateral q dual to (v0,v1).
   *         Vertices (w0,w1,w2,w3) are listed in order around q.
   *  @pre  Orientation of (w0,w1,v0,v1) is positive.
   *  @param epsilon Tolerance.  
   *         Determinants less than epsilon are considered equivalent to 0.
   */
  template <typename VCOORD_TYPE, typename WCOORD_TYPE,
            typename EPSILON_TYPE>
  bool is_in_envelope_3D
  (const VCOORD_TYPE v0[3], const VCOORD_TYPE v1[3],
   const WCOORD_TYPE w0[3], const WCOORD_TYPE w1[3],
   const WCOORD_TYPE w2[3], const WCOORD_TYPE w3[3],
   const EPSILON_TYPE epsilon)
  {
    double D;

    IJK::determinant_point_3D(w0, w2, v1, w1, D);
    if (D < epsilon) { return(false); }
    IJK::determinant_point_3D(w0, w2, v0, w1, D);
    if (-D < epsilon) { return(false); }
    IJK::determinant_point_3D(w0, w2, v1, w3, D);
    if (-D < epsilon) { return(false); }
    IJK::determinant_point_3D(w0, w2, v0, w3, D);
    if (D < epsilon) { return(false); }

    return(true);
  }

  ///@}

}

#endif
