/*!
 *  @file ijktri2D_angle.tpp
 *  @brief ijk templates for computing angles in 2D meshes.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2023 Rephael Wenger

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

#ifndef _IJKTRI2D_ANGLE_
#define _IJKTRI2D_ANGLE_

#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkisocoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijktri2D_info.tpp"
#include "ijktriangulate_poly2D.tpp"


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // ***************************************************************
  //! @name COMPUTE COS MIN TRIANGLE ANGLE
  // ***************************************************************

  //@{

  /*!
   *  @brief Compute the cosine of smallest angle of triangle 
   *    (coord0, coord1, coord2).
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
   IJK::POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   poly_tri_result)
  {
    IJK::compute_cos_min_triangle_angle
      (dimension, coord0, coord1, coord2, max_small_magnitude,
       poly_tri_result.cos_min_triangulation_angle,
       poly_tri_result.flag_zero);

    poly_tri_result.num_triangles = 1;
    poly_tri_result.tri_vertex_index = 0;
    poly_tri_result.num_interior_tri_vertices = 0;
    poly_tri_result.triangulation_encoding.Clear();
  }

  //@}


  // ***************************************************************
  //! @name COMPUTE COS MIN TRIANGULATION ANGLE
  // ***************************************************************

  //@{

  /*!
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in a fan triangulation.
   *  - Fan triangulation is from vertex iv0.
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
  void compute_cos_min_fan_triangulation_angle
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

      compute_cos_min_triangle_angle
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
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in fan triangulation.
   *  - Fan triangulation is from vertex iv0.
   *  - Version using C++ STL vector vertex_coord[].
   */
  template <typename DTYPE, typename NTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename IVTYPE,typename COS_TYPE>
  void compute_cos_min_fan_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vcoord, 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   const IVTYPE index_v0, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_fan_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), num_poly_vert,
       poly_vert, max_small_magnitude, index_v0, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in fan triangulation.
   *  - Fan triangulation is from vertex iv0.
   *  - Version using list of pointers to vertex coordinates.
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
  void compute_cos_min_fan_triangulation_angleP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE * const vcoord_ptr[], const MTYPE max_small_magnitude, 
   const IVTYPE index_v0, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    NTYPE i1 = (index_v0+1)%num_poly_vert;
    NTYPE i2 = (i1+1)%num_poly_vert;
    while (i2 != index_v0) {
      
      compute_cos_min_triangle_angle
        (dimension, vcoord_ptr[index_v0], vcoord_ptr[i1], 
         vcoord_ptr[i2], max_small_magnitude, cosA, flagA);

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
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in triangulation from interior vertex.
   *  - Triangulation is from interior vertex at vcoordX[].
   *  - Version using list of pointers to vertex coordinates.
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
  void compute_cos_min_triangulation_angle_using_vcoordXP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE * const vcoord_ptr[], const CTYPEX vcoordX[], 
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
      const NTYPE i1 = (i0+1)%num_poly_vert;
      
      compute_cos_min_triangle_angle
        (dimension, vcoord_ptr[i0], vcoord_ptr[i1], vcoordX, 
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
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in triangulation from interior vertex.
   *  - Triangulation is from interior vertex at vcoordX[].
   *  - Version using vcoord[] and poly_vert[].
   *  @param vcoord[] Array of coordinates.  
   *  - Not all coordinates are necessarily polygon vertex coordinates.
   *  @param vcoordX[] Internal vertex coordinates.
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] Array of polygon vertices indexing vcoord[].
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename NTYPE, typename VTYPE,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangulation_angle_using_vcoordX
  (const DTYPE dimension, const CTYPE vcoord[], const CTYPEX vcoordX[],
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
      const NTYPE i1 = (i0+1)%num_poly_vert;

      const CTYPE * v0_coord = vcoord+poly_vert[i0]*dimension;
      const CTYPE * v1_coord = vcoord+poly_vert[i1]*dimension;
      
      compute_cos_min_triangle_angle
        (dimension, v0_coord, v1_coord, vcoordX, max_small_magnitude, 
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
    }
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in triangulation of a polygon minus two edges.
   *  - Triangulation is from interior vertex at vcoordX[] of a polygon
   *    minus two edges incident on iloc1.
   *  - Version using vcoord[] and poly_vert[].
   *  @param vcoord[] Array of coordinates.  
   *  @param vcoordX[] Internal vertex coordinates.
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] Array of polygon vertices indexing vcoord[].
   *  @param iloc1 Edges incident on vertex at iloc1 are not part
   *         of the triangulation.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename NTYPE, typename VTYPE, typename ITYPE,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangulationX_angle_on_polygon_minus_two_edges
  (const DTYPE dimension, const CTYPE vcoord[], const CTYPEX vcoordX[],
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const ITYPE iloc1, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const ITYPE iloc0 = (iloc1 + num_poly_vert - 1)%num_poly_vert;
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {

      if (i0 == iloc0 || i0 == iloc1) {
        // Skip triangles with edges incident on poly_vert[iloc1].
        continue;
      }

      const NTYPE i1 = (i0+1)%num_poly_vert;

      const CTYPE * v0_coord = vcoord+poly_vert[i0]*dimension;
      const CTYPE * v1_coord = vcoord+poly_vert[i1]*dimension;
      
      compute_cos_min_triangle_angle
        (dimension, v0_coord, v1_coord, vcoordX, max_small_magnitude, 
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
    }
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in triangulation of a polygon minus two edges.
   *  - Triangulation is from interior vertex at vcoordX[] of a polygon
   *    minus two edges incident on iloc1.
   *  - Version using list of pointers to vertex coordinates.
   *  @param vcoord[] Array of coordinates.  
   *  @param vcoordX[] Internal vertex coordinates.
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] Array of polygon vertices indexing vcoord[].
   *  @param iloc1 Edges incident on vertex at iloc1 are not part
   *         of the triangulation.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename NTYPE, typename ITYPE,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangulationX_angle_on_polygon_minus_two_edgesP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE * const vcoord_ptr[], const CTYPEX vcoordX[],
   const ITYPE iloc1, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const ITYPE iloc0 = (iloc1 + num_poly_vert - 1)%num_poly_vert;
    COS_TYPE cosA;
    bool flagA;

    cos_min_angle = -1;
    flag_zero = false;

    for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {

      if (i0 == iloc0 || i0 == iloc1) {
        // Skip triangles with edges incident on poly_vert[iloc1].
        continue;
      }

      const NTYPE i1 = (i0+1)%num_poly_vert;
      
      compute_cos_min_triangle_angle
        (dimension, vcoord_ptr[i0], vcoord_ptr[i1], vcoordX, 
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
   *  @brief Compute the cosine of the smallest angle 
   *    in the triangulation of polygon minus ear.
   *  - Polygon vertex iloc_ear is cut off in an ear
   *    and all other vertices joined to vcoordX[].
   *  - Version using vcoord[] and poly_vert[].
   *  @param vcoord[] Array of coordinates.  
   *  - Not all coordinates are necessarily polygon vertex coordinates.
   *  @param vcoordX[] Internal vertex coordinates.
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] Array of polygon vertices indexing vcoord[].
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename NTYPE, typename VTYPE,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangulation_angle_using_cut_ear_vcoordX
  (const DTYPE dimension, const CTYPE vcoord[], const CTYPEX vcoordX[],
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const NTYPE iloc_ear, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const NTYPE iloc0 = (iloc_ear+num_poly_vert-1)%num_poly_vert;
    const NTYPE iloc2 = (iloc_ear+1)%num_poly_vert;
    const VTYPE iv0 = poly_vert[iloc0];
    const VTYPE iv1 = poly_vert[iloc_ear];
    const VTYPE iv2 = poly_vert[iloc2];
    const CTYPE * vcoord0 = vcoord + iv0*dimension;
    const CTYPE * vcoord1 = vcoord + iv1*dimension;
    const CTYPE * vcoord2 = vcoord + iv2*dimension;
    COS_TYPE cos_triangleA, cos_triangleB;
    bool flagA_zero, flagB_zero;

    cos_min_angle = -1;
    flag_zero = false;

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoord2, max_small_magnitude,
       cos_triangleA, flagA_zero);

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoordX, vcoord2, max_small_magnitude,
       cos_triangleB, flagB_zero);

    compute_cos_min_triangulationX_angle_on_polygon_minus_two_edges
      (dimension, vcoord, vcoordX, num_poly_vert, poly_vert,
       iloc_ear, max_small_magnitude, cos_min_angle, flag_zero);

    flag_zero = (flag_zero || flagA_zero || flagB_zero);
    cos_min_angle = std::max(cos_min_angle, cos_triangleA);
    cos_min_angle = std::max(cos_min_angle, cos_triangleB);
  }


  /*!
   *  @brief Compute the cosine of the smallest angle 
   *    in the triangulation of polygon minus ear.
   *  - Polygon vertex iloc_ear is cut off in an ear
   *    and all other vertices joined to vcoordX[].
   *  - Version using list of pointers to vertex coordinates.
   *  @param vcoord[] Array of coordinates.  
   *  - Not all coordinates are necessarily polygon vertex coordinates.
   *  @param vcoordX[] Internal vertex coordinates.
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] Array of polygon vertices indexing vcoord[].
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename NTYPE, typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangulation_angle_using_cut_ear_vcoordXP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE * const vcoord_ptr[], const CTYPEX vcoordX[],   
   const NTYPE iloc_ear, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const NTYPE iloc0 = (iloc_ear+num_poly_vert-1)%num_poly_vert;
    const NTYPE iloc2 = (iloc_ear+1)%num_poly_vert;
    COS_TYPE cos_triangleA, cos_triangleB;
    bool flagA_zero, flagB_zero;

    cos_min_angle = -1;
    flag_zero = false;

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord_ptr[iloc0], vcoord_ptr[iloc_ear], 
       vcoord_ptr[iloc2], max_small_magnitude, 
       cos_triangleA, flagA_zero);

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord_ptr[iloc0], vcoordX, vcoord_ptr[iloc2], 
       max_small_magnitude, cos_triangleB, flagB_zero);

    compute_cos_min_triangulationX_angle_on_polygon_minus_two_edgesP
      (dimension, num_poly_vert, vcoord_ptr, vcoordX,
       iloc_ear, max_small_magnitude, cos_min_angle, flag_zero);

    flag_zero = (flag_zero || flagA_zero || flagB_zero);
    cos_min_angle = std::max(cos_min_angle, cos_triangleA);
    cos_min_angle = std::max(cos_min_angle, cos_triangleB);

  }


  /*!
   *  @brief Compute the cosine of the smallest angle 
   *    in the triangulation  of polygon minus ear.
   *  - Polygon vertex iloc_ear is cut off in an ear
   *    and all other vertices joined to vcoordX[].
   *  - Version using list of pointers to vertex coordinates.
   *  - Version trying two interior coordinates.
   *  @param vcoord[] Array of coordinates.  
   *  - Not all coordinates are necessarily polygon vertex coordinates.
   *  @param vcoordX[] Internal vertex coordinates.
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] Array of polygon vertices indexing vcoord[].
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename NTYPE, typename MTYPE, typename COS_TYPE,
            typename ITYPE>
  void compute_cos_min_triangulation_angle_using_cut_ear_vcoordXPx2
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE * const vcoord_ptr[], 
   const CTYPEX vcoordX0[], const CTYPEX vcoordX1[],
   const NTYPE iloc_ear, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, 
   ITYPE & index_selected_IV, bool & flag_zero)
  {
    COS_TYPE cosB;
    bool flagB_zero;

    index_selected_IV = 0;
    compute_cos_min_triangulation_angle_using_cut_ear_vcoordXP
      (dimension, num_poly_vert, vcoord_ptr,
       vcoordX0, iloc_ear, max_small_magnitude,
       cos_min_angle, flag_zero);

    compute_cos_min_triangulation_angle_using_cut_ear_vcoordXP
      (dimension, num_poly_vert, vcoord_ptr,
       vcoordX1, iloc_ear, max_small_magnitude,
       cosB, flagB_zero);

    if (flag_zero || (!flagB_zero && cosB < cos_min_angle)) {
      cos_min_angle = cosB;
      flag_zero = flagB_zero;
      index_selected_IV = 1;
    }

  }


  /*!
   *  @brief Compute cosine of the min angle of a triangle 
   *    in a triangulation.
   *  - Polygon triangulation is given by tri_encoding.
   */
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPE, typename MTYPE, typename COS_TYPE,
            int BIT_SET_SIZE>
  void compute_cos_min_encoded_triangulation_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> & tri_encoding,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_TRIANGLE(3);
    VTYPE tri_vert1[NUM_VERT_PER_TRIANGLE];
    const CTYPE * tri_vcoord[NUM_VERT_PER_TRIANGLE];
    std::vector<VTYPE> tri_list;
    COS_TYPE cosX;
    bool flagX_zero;

    // Initialize
    cos_min_angle = -1.0;
    flag_zero = false;

    gen_polygon_triangulation_list_based_on_encoding
      (num_poly_vert, tri_encoding, tri_list);

    const NTYPE num_triangles = 
      tri_list.size()/NUM_VERT_PER_TRIANGLE;


    for (NTYPE i = 0; i < num_triangles; i++) {

      for (NTYPE k = 0; k < NUM_VERT_PER_TRIANGLE; k++) {
        const NTYPE iloc = tri_list[k + NUM_VERT_PER_TRIANGLE*i]; 
        tri_vert1[k] = poly_vert[iloc];
      }

      IJK::compute_cos_min_triangle_angle
        (dimension, vcoord, tri_vert1, max_small_magnitude, cosX, flagX_zero);

      flag_zero = (flagX_zero || flag_zero);
      cos_min_angle = std::max(cos_min_angle, cosX);
    }

  }


  /*!
   *  @brief Compute cosine of the min angle of a triangle 
   *    in a triangulation.
   *  - Polygon triangulation is given by tri_encoding.
   *  - Version using list of pointers to vertex coordinates.
   */
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename MTYPE, typename COS_TYPE,
            int BIT_SET_SIZE>
  void compute_cos_min_encoded_triangulation_angleP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE * const vcoord_ptr[], 
   const POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> & tri_encoding,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_TRIANGLE(3);
    const CTYPE * tri_vcoord[NUM_VERT_PER_TRIANGLE];
    std::vector<NTYPE> tri_list;
    COS_TYPE cosX;
    bool flagX_zero;

    // Initialize
    cos_min_angle = -1.0;
    flag_zero = false;

    gen_polygon_triangulation_list_based_on_encoding
      (num_poly_vert, tri_encoding, tri_list);

    const NTYPE num_triangles = 
      tri_list.size()/NUM_VERT_PER_TRIANGLE;

    for (NTYPE i = 0; i < num_triangles; i++) {

      for (NTYPE k = 0; k < NUM_VERT_PER_TRIANGLE; k++) {
        const NTYPE iloc = tri_list[k + NUM_VERT_PER_TRIANGLE*i]; 
        tri_vcoord[k] = vcoord_ptr[iloc];
      }

      IJK::compute_cos_min_triangle_angle
        (dimension, tri_vcoord, max_small_magnitude, cosX, flagX_zero);

      flag_zero = (flagX_zero || flag_zero);
      cos_min_angle = std::max(cos_min_angle, cosX);
    }

  }

  //@}


  // ***************************************************************
  //! @name COMPUTE COS MIN PENTAGON TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute the cosine of the smallest triangle angle
   *    in the triangulation of a pentagon into five triangles.
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

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoordX, max_small_magnitude,
       cos[0], flag_tri_zero[0]);
    compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord2, vcoordX, max_small_magnitude,
       cos[1], flag_tri_zero[1]);
    compute_cos_min_triangle_angle
      (dimension, vcoord2, vcoord3, vcoordX, max_small_magnitude,
       cos[2], flag_tri_zero[2]);
    compute_cos_min_triangle_angle
      (dimension, vcoord3, vcoord4, vcoordX, max_small_magnitude,
       cos[3], flag_tri_zero[3]);
    compute_cos_min_triangle_angle
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
   *  @brief Compute the cosine of the smallest triangle angle
   *    in the triangulation of a pentagon into five triangles.
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
   *  @brief Compute the cosine of the smallest triangle angle
   *    in the triangulation of a pentagon into five triangles.
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


  /// @brief Compute min angle in triangulation of a pentagon
  ///   with all triangles incident on vertex_coord0.
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


  /// @brief Compute min angle in triangulation of a pentagon
  ///   with all triangles incident on vertex_coord0.
  /// - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & pentagon_tri)
  {
    compute_cos_min_pentagon_triangulation_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, max_small_magnitude,
       pentagon_tri.cos_min_triangulation_angle, pentagon_tri.flag_zero);

    pentagon_tri.num_triangles = 3;
  }

  ///@}

  
  // ***************************************************************
  //! @name SPLIT MAX ANGLE
  // ***************************************************************

  //@{


  /*!
   *  @brief Triangulate a single polygon by splitting maximum polygon angle.
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

    compute_cos_triangle_angle_coord_list
      (dimension, vertex_coord, iv0, iv1, iv2, max_small_magnitude,
       min_cos, flag_zero);

    for (NTYPE i = 1; i < num_poly_vert; i++) {
      iv0 = poly_vert[i-1];
      iv1 = poly_vert[i];
      iv2 = poly_vert[((i+1)%num_poly_vert)];
      CTYPE cos_v0;

      compute_cos_triangle_angle_coord_list
        (dimension, vertex_coord, iv0, iv1, iv2, max_small_magnitude,
         cos_v0, flag_zero);

      if (cos_v0 < min_cos) {
        k_min_cos = i;
        min_cos = cos_v0;
      }
    }

    triangulate_polygon_from_vertex(num_poly_vert, poly_vert, k_min_cos, tri_vert);
  }


  /// @brief Triangulate each polygon in list by splitting maximum polygon angle.
  template <typename DTYPE, typename CTYPE,
            typename NTYPEV, typename NTYPEP, 
            typename VTYPE0, typename VTYPE1,
            typename MTYPE, typename ITYPE,
            typename NTYPE_TRI>
  void triangulate_polygon_list_split_max_angle
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const NTYPEV * num_poly_vert, const VTYPE0 * poly_vert,
   const ITYPE * first_poly_vert, const NTYPEP num_poly,
   const MTYPE max_small_magnitude,   
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_poly)
  {
    num_triangulated_poly = 0;
    
    for (NTYPEP ipoly = 0; ipoly < num_poly; ipoly++) {
      triangulate_polygon_split_max_angle
        (dimension, vertex_coord, num_poly_vert[ipoly], 
         poly_vert+first_poly_vert[ipoly], max_small_magnitude, tri_vert);

      num_triangulated_poly++;
    }
  }


  /*!
   *  @brief Triangulate each polygon in list by splitting maximum polygon angle.
   *  - Version using C++ STL vectors vertex_coord[], num_poly_vert[],
   *    poly_vert[] and first_poly_vert[].
   */
  template <typename DTYPE, typename CTYPE, typename NTYPEV,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE, typename ITYPE,
            typename NTYPE_TRI>
  void triangulate_polygon_list_split_max_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const std::vector<NTYPEV> & num_poly_vert, 
   const std::vector<VTYPE0> & poly_vert,
   const std::vector<ITYPE> & first_poly_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_poly)
  {
    triangulate_polygon_list_split_max_angle
      (dimension, IJK::vector2pointer(vertex_coord), 
       IJK::vector2pointer(num_poly_vert), IJK::vector2pointer(poly_vert),
       IJK::vector2pointer(first_poly_vert), num_poly_vert.size(),
       max_small_magnitude, tri_vert, num_triangulated_poly);
  }

  
  /*!
   *  @brief Triangulate each polygon in list by splitting maximum polygon angle.
   *  - Version using C++ STL vectors vertex_coord[], num_poly_vert[],
   *    poly_vert[] and first_poly_vert[].
   *  - Version that does NOT return num_triangulate_poly.
   */
  template <typename DTYPE, typename CTYPE, typename NTYPEV,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE, typename ITYPE>
  void triangulate_polygon_list_split_max_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const std::vector<NTYPEV> & num_poly_vert, 
   const std::vector<VTYPE0> & poly_vert,
   const std::vector<ITYPE> & first_poly_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert)
  {
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;
    SIZE_TYPE num_triangulated_poly;

    triangulate_polygon_list_split_max_angle
      (dimension, vertex_coord, num_poly_vert, poly_vert,
       first_poly_vert, max_small_magnitude,
       tri_vert, num_triangulated_poly);
  }
  
  //@}
  

  // ***************************************************************
  //! @name TRIANGULATE QUADRILATERALS - SPLIT MAX ANGLE
  // ***************************************************************

  //@{

  /// @brief Triangulate each quadrilateral in list by splitting maximum quad angle.
  template <typename DTYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPE0, typename VTYPE1, typename MTYPE,
            typename NTYPE_TRI>
  void triangulate_quad_split_max_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const NTYPEQ num_quad,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad)
  {
    const int NUM_VERT_PER_QUAD (4);

    num_triangulated_quad = 0;
    
    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      const NTYPEQ k = iquad*NUM_VERT_PER_QUAD;
      triangulate_polygon_split_max_angle
        (dimension, vertex_coord, NUM_VERT_PER_QUAD, quad_vert+k, 
         max_small_magnitude, tri_vert);
      num_triangulated_quad++;
    }
  }


  /*!
   *  @brief Triangulate each quadrilateral in list by splitting maximum quad angle.
   *  - Version using C++ STL vectors vertex_coord[] and quad_vert[].
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE,
            typename NTYPE_TRI>
  void triangulate_quad_split_max_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE0> & quad_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_QUAD(4);

    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    triangulate_quad_split_max_angle
      (dimension, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert), num_quad, max_small_magnitude,
       tri_vert, num_triangulated_quad);
  }

  //@}

}

#endif
