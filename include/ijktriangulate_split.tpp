/// \file ijktriangulate_split.tpp
/// @brief ijk templates for triangulating split polygons
///   or polygons with split edges.
/// - Includes code based on geometry.
/// - Version 0.4.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2021 Rephael Wenger

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

#ifndef _IJKTRIANGULATE_SPLITE_
#define _IJKTRIANGULATE_SPLITE_

#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijktri2D_info.tpp"


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // ***************************************************************
  //! @name COMPUTE COS MAX MIN SUBQUAD TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    when the quad is split into two subquads
   *      (v0, v1, vX12, vX03) and (vX03, vX12, v2, v3).
   *  - Quad vertices have coordinates (vcoord0, vcoord1, vcoord2, vcoord3)
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vcoord0[] Coordinates of quad vertex 0.
   *  @param vcoord1[] Coordinates of quad vertex 1.
   *  @param vcoord2[] Coordinates of quad vertex 2.
   *  @param vcoord3[] Coordinates of quad vertex 3.
   *  - Note: upper_coord1[] precedes upper_coord0[] in the parameter list.
   *  @param vcoordX03[] Coordinate of point on edge (v0, v3).
   *  @param vcoordX12[] Coordinate of point on edge (v1, v2).
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_diag02[0] If true, quad (v0, v1, vX12, vX03)
   *       with max min triangulation angle has diagonal (v0, vX12).
   *       Otherwise, quad (v0, v1, vX12, vX01) has diagonal (v1, vX03).
   *  @param[out] flag_diag02[1] If true, quad (vX03, vX12, v2, v3)
   *       with max min triangulation angle has diagonal (vX03, v2).
   *       Otherwise, quad (vX03, vX12, v2, v3) has diagonal (vX12, v3).
   *  @param[out] flag_zero True, if some triangles in triangulation
   *    which maximizes min angle have two or three edges less than or equal
   *    to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPEX,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12
  (const DTYPE dimension, 
   const CTYPE0 vcoord0[], const CTYPE1 vcoord1[], 
   const CTYPE2 vcoord2[], const CTYPE3 vcoord3[],
   const CTYPEX vcoordX03[], const CTYPEX vcoordX12[],
   const MTYPE max_small_magnitude, COS_TYPE & cos_max_min_angle, 
   bool flag_diag02[2], bool & flag_zero)
  {
    COS_TYPE cos_angle[2];
    bool flag[2];

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord1, vcoordX12, vcoordX03, max_small_magnitude,
       cos_angle[0], flag_diag02[0], flag[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoordX03, vcoordX12, vcoord2, vcoord3, max_small_magnitude,
       cos_angle[1], flag_diag02[1], flag[1]);

    cos_max_min_angle = std::max(cos_angle[0], cos_angle[1]);
    flag_zero = (flag[0] || flag[1]);
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    when the quad is split into two subquads
   *      (v0, v1, vX12, vX03) and (vX03, vX12, v2, v3).
   *  - Version returning data structure POLY_TRIANGULATION_RESULTX2.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPEX,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12
  (const DTYPE dimension, 
   const CTYPE0 vcoord0[], const CTYPE1 vcoord1[], 
   const CTYPE2 vcoord2[], const CTYPE3 vcoord3[],
   const CTYPEX vcoordX03[], const CTYPEX vcoordX12[],
   const MTYPE max_small_magnitude, 
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    bool flag_diag02[2], flag_zero;

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
       vcoordX03, vcoordX12, max_small_magnitude, 
       quad_tri_result.cos_min_triangulation_angle, flag_diag02, flag_zero);

    quad_tri_result.num_triangles = 4;
    quad_tri_result.num_split_edges = 2;
    quad_tri_result.num_tri_ears = 2;

    if (flag_diag02[0]) 
      { quad_tri_result.ear_list[0] = 1; }
    else 
      { quad_tri_result.ear_list[0] = 0; }

    if (flag_diag02[1]) 
      { quad_tri_result.ear_list[1] = 3; }
    else
      { quad_tri_result.ear_list[1] = 2; }
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    when the quad is split into two quads at vcoordX03 and vcoordX12.
   *  - Version using lower_vcoord[] and upper_vcoord[].
   *  - Version returning data structure POLY_TRIANGULATION_RESULTX2_OLD.
   */
  template <typename DTYPE, typename CTYPEL, typename CTYPEU, typename CTYPEX,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12_LU
  (const DTYPE dimension, 
   const CTYPEL * lower_vcoord[2], const CTYPEU upper_vcoord[2],
   const CTYPEX vcoordX03[], const CTYPEX vcoordX12[],
   const MTYPE max_small_magnitude, 
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    // Note: upper_vcoord[1] precedes upper_vcoord[0] in the parameter list.
    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12
      (dimension, lower_vcoord[0], lower_vcoord[1], upper_vcoord[1], 
       upper_vcoord[0], vcoordX03, vcoordX12, 
       max_small_magnitude, quad_tri_result);
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    when the quad is split into two subquads
   *      (v0, vX01, vX23, v3) and (vX01, v1, v2, vX23).
   *  - Quad vertices have coordinates (vcoord0, vcoord1, vcoord2, vcoord3)
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vcoord0[] Coordinates of quad vertex 0.
   *  @param vcoord1[] Coordinates of quad vertex 1.
   *  @param vcoord2[] Coordinates of quad vertex 2.
   *  @param vcoord3[] Coordinates of quad vertex 3.
   *  - Note: upper_coord1[] precedes upper_coord0[] in the parameter list.
   *  @param vcoordX01[] Coordinate of point on edge (v0, v1).
   *  @param vcoordX23[] Coordinate of point on edge (v2, v3).
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_diag02[0] If true, quad (v0, vX01, vX23, v3)
   *       with max min triangulation angle has diagonal (v0, vX23).
   *       Otherwise, quad (v0, vX01, vX23, v3) has diagonal (vX01, v3).
   *  @param[out] flag_diag02[1] If true, quad (vX01, v1, v2, vX23)
   *       with max min triangulation angle has diagonal (vX01, v2).
   *       Otherwise, quad (vX01, v1, v2, vX23) has diagonal (v1, vX23).
   *  @param[out] flag_zero True, if some triangles in triangulation
   *    which maximizes min angle have two or three edges less than or equal
   *    to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPEX,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_max_min_quad_split_into_two_subquads_vX01_vX23
  (const DTYPE dimension, 
   const CTYPE0 vcoord0[], const CTYPE1 vcoord1[], 
   const CTYPE2 vcoord2[], const CTYPE3 vcoord3[],
   const CTYPEX vcoordX01[], const CTYPEX vcoordX23[],
   const MTYPE max_small_magnitude, COS_TYPE & cos_max_min_angle, 
   bool flag_diag02[2], bool & flag_zero)
  {
    COS_TYPE cos_angle[2];
    bool flag[2];

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoordX01, vcoordX23, vcoord3, max_small_magnitude,
       cos_angle[0], flag_diag02[0], flag[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoordX01, vcoord1, vcoord2, vcoordX23, max_small_magnitude,
       cos_angle[1], flag_diag02[1], flag[1]);

    cos_max_min_angle = std::max(cos_angle[0], cos_angle[1]);
    flag_zero = (flag[0] || flag[1]);
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    when the quad is split into two subquads
   *      (v0, vX01, vX23, v3) and (vX01, v1, v2, vX23).
   *  - Version returning data structure POLY_TRIANGULATION_RESULTX2_OLD.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPEX,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split_into_two_subquads_vX01_vX23
  (const DTYPE dimension, 
   const CTYPE0 vcoord0[], const CTYPE1 vcoord1[], 
   const CTYPE2 vcoord2[], const CTYPE3 vcoord3[],
   const CTYPEX vcoordX01[], const CTYPEX vcoordX23[],
   const MTYPE max_small_magnitude, 
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    bool flag_diag02[2], flag_zero;

    compute_cos_max_min_quad_split_into_two_subquads_vX01_vX23
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
       vcoordX01, vcoordX23, max_small_magnitude, 
       quad_tri_result.cos_min_triangulation_angle, flag_diag02, flag_zero);

    quad_tri_result.num_triangles = 4;
    quad_tri_result.num_split_edges = 2;
    quad_tri_result.num_tri_ears = 2;

    if (flag_diag02[0]) 
      { quad_tri_result.ear_list[0] = 3; }
    else 
      { quad_tri_result.ear_list[0] = 0; }

    if (flag_diag02[1]) 
      { quad_tri_result.ear_list[1] = 1; }
    else
      { quad_tri_result.ear_list[1] = 2; }
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    when the quad is split into two subquads by either a vertical
   *    or a horizontal edge.
   *  - Quad vertices have coordinates (vcoord0, vcoord1, vcoord2, vcoord3)
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vcoord0[] Coordinates of quad vertex 0.
   *  @param vcoord1[] Coordinates of quad vertex 1.
   *  @param vcoord2[] Coordinates of quad vertex 2.
   *  @param vcoord3[] Coordinates of quad vertex 3.
   *  @param vcoordX01_or_X03[] Coordinate of point either on edge (v0, v1)
   *       or on edge (v0, v3).
   *  @param vcoordX23_or_X12[] Coordinate of point either on edge (v2, v3)
   *       or on edge (v1, v2).           
   *  @param max_small_magnitude Vectors with magnitude less than or
   *       equal to max_small_magnitude are set to 0.
   *  @param split_edge_index Index of one of the split edges.
   *       One split edge is (split_edge_index, (split_edge_index+1) mod 4).
   *       The other split edge is the parallel edge:
   *       ((split_edge_index+2) mod 4, (split_edge_index+3) mod 4).  
   *  @pre max_small_magnitude >= 0.
   *  @param[out] quad_tri_result Quad triangulation result.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPEX,
            typename ITYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split_into_two_subquads
  (const DTYPE dimension, 
   const CTYPE0 vcoord0[], const CTYPE1 vcoord1[], 
   const CTYPE2 vcoord2[], const CTYPE3 vcoord3[],
   const CTYPEX vcoordX01_or_X03[], const CTYPEX vcoordX23_or_X12[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    if ((split_edge_index == 0) || (split_edge_index == 2)) {
      compute_cos_max_min_quad_split_into_two_subquads_vX01_vX23
        (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
         vcoordX01_or_X03, vcoordX23_or_X12, max_small_magnitude,
         quad_tri_result);
    }
    else {
      compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12
        (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
         vcoordX01_or_X03, vcoordX23_or_X12, max_small_magnitude,
         quad_tri_result);
    }
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    when the quad is split into two subquads by either a vertical
   *    or a horizontal edge.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename ITYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split_into_two_subquads
  (const DTYPE dimension, 
   const CTYPE * vcoord[4],
   const CTYPEX vcoordX01_or_X03[], const CTYPEX vcoordX23_or_X12[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    compute_cos_max_min_quad_split_into_two_subquads
      (dimension, vcoord[0], vcoord[1], vcoord[2], vcoord[3],
       vcoordX01_or_X03, vcoordX23_or_X12, split_edge_index,
       max_small_magnitude, quad_tri_result);
  }


  /*!
   *  Compute the cosine of the max min quad triangulation angle
   *    when the quad is split into two subquads by either a vertical
   *    or a horizontal edge.
   * - Version with list of vertices quad_vert[].
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename VTYPE, typename ITYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split_into_two_subquads
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX vcoordX01_or_X03[], const CTYPEX vcoordX23_or_X12[],
   const VTYPE quad_vert[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord[NUM_VERT_PER_QUAD] =
      { vertex_coord_ptr+quad_vert[0]*dimension,
        vertex_coord_ptr+quad_vert[1]*dimension,
        vertex_coord_ptr+quad_vert[2]*dimension,
        vertex_coord_ptr+quad_vert[3]*dimension };

    compute_cos_max_min_quad_split_into_two_subquads
      (dimension, vcoord, vcoordX01_or_X03, vcoordX23_or_X12, 
       split_edge_index, max_small_magnitude, quad_tri_result);
  }

  ///@}


  // ***************************************************************
  //! @name COMPUTE COS MAX MIN PENTAGON TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  Compute max of min quad triangulation angle where quad is split
   *    into 3 triangles.
   *  - Vertex is added between v0 and v3.
   *  - Triangulate with 3 triangles incident on v1, v2 or vcoordX03[].
   *  @param vcoordX03[] Coordinates of vertex on edge (v0,v3).
   *  @param[out] tri_vertex_index 0,1,2,3,4 of vertex used in
   *     triangulation into three triangles.
   *    - If (tri_vertex_index == 4), then triangulate from vcoordX03[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, typename CTYPEX,
            typename MTYPE, typename COS_TYPE, typename ITYPE>
  void compute_cos_max_min_quad_tri3_angle_vX03
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX03[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, ITYPE & tri_vertex_index, bool & flag_zero)
  {
    COS_TYPE cos_min_v1_angle, cos_min_v2_angle, cos_min_vX03_angle;
    bool flag_zero_v1, flag_zero_v2, flag_zero_vX03;

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vcoordX03, vertex_coord0, vertex_coord1, vertex_coord2, 
       vertex_coord3, max_small_magnitude, cos_min_vX03_angle, 
       flag_zero_vX03);
    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vertex_coord1, vertex_coord2, vertex_coord3, vcoordX03,
       vertex_coord0, max_small_magnitude, cos_min_v1_angle, 
       flag_zero_v1);
    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vertex_coord2, vertex_coord3, vcoordX03, vertex_coord0, 
       vertex_coord1, max_small_magnitude, cos_min_v2_angle, 
       flag_zero_v2);

    int index_selected;
    select_minIII
      (cos_min_v1_angle, flag_zero_v1, 
       cos_min_v2_angle, flag_zero_v2, 
       cos_min_vX03_angle, flag_zero_vX03, 
       cos_min_angle, index_selected, flag_zero);

    if (index_selected == 0)
      { tri_vertex_index = 1; }
    else if (index_selected == 1)
      { tri_vertex_index = 2; }
    else
      { tri_vertex_index = 4; }
  }


  /*!
   *  Compute max of min quad triangulation angle where quad is split
   *    into 3 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate with 3 triangles incident on vcoordX03[]
   *      or (split_edge_index+2) mod 4 or (split_edge_index+3) mod 4.
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param[out] tri_vertex_index 0,1,2,3,4 of vertex used in
   *     triangulation into three triangles.
   *    - If (tri_vertex_index == 4), then triangulate from vcoordX[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, typename CTYPEX,
            typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_max_min_quad_tri3_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX[],
   const ITYPE0 split_edge_index,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_max_min_angle, ITYPE1 & tri_vertex_index, bool & flag_zero)
  {
    const int NUM_VERT_PER_QUAD(4);

    switch (split_edge_index) {

    case 0:
      // Edge (v0,v1) is split.
      compute_cos_max_min_quad_tri3_angle_vX03
        (dimension, vertex_coord1, vertex_coord2,
         vertex_coord3, vertex_coord0, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      if (tri_vertex_index < NUM_VERT_PER_QUAD) {
        tri_vertex_index = (tri_vertex_index+1)%NUM_VERT_PER_QUAD;
      }
      break;

    case 1:
      // Edge (v1,v2) is split.
      compute_cos_max_min_quad_tri3_angle_vX03
        (dimension, vertex_coord2, vertex_coord3, vertex_coord0, 
         vertex_coord1, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      if (tri_vertex_index < NUM_VERT_PER_QUAD) {
        tri_vertex_index = (tri_vertex_index+2)%NUM_VERT_PER_QUAD;
      }
      break;

    case 2:
      // Edge (v2,v3) is split.
      compute_cos_max_min_quad_tri3_angle_vX03
        (dimension, vertex_coord3, vertex_coord0, vertex_coord1, 
         vertex_coord2, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      if (tri_vertex_index < NUM_VERT_PER_QUAD) {
        tri_vertex_index = (tri_vertex_index+3)%NUM_VERT_PER_QUAD;
      }
      break;

    case 3:
    default:
      // Edge (v3,v0) is split.
      compute_cos_max_min_quad_tri3_angle_vX03
        (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
         vertex_coord3, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      return;
    }

  }


  /*!
   *  Compute max of min quad triangulation angle where quad is split
   *    into 5 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate with 5 triangles incident on vcoordY[].
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param[out] tri_vertex_index 0,1,2,3,4 of vertex used in
   *     triangulation into three triangles.
   *    - If (tri_vertex_index == 4), then triangulate from vcoordX[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename ITYPE, typename MTYPE, typename COS_TYPE>
  void compute_cos_min_quad_tri5_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_max_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_QUAD(4);

    switch (split_edge_index) {

    case 0:
      // Edge (v0,v1) is split.
      IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
        (dimension, vertex_coord0, vcoordX, vertex_coord1,
         vertex_coord2, vertex_coord3, vcoordY, max_small_magnitude,
         cos_max_min_angle, flag_zero);
      break;

    case 1:
      // Edge (v1,v2) is split.
      IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
        (dimension, vertex_coord0, vertex_coord1, vcoordX,
         vertex_coord2, vertex_coord3, vcoordY, max_small_magnitude,
         cos_max_min_angle, flag_zero);
      break;

    case 2:
      // Edge (v2,v3) is split.
      IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
        (dimension, vertex_coord0, vertex_coord1, vertex_coord2, 
         vcoordX, vertex_coord3, vcoordY, max_small_magnitude,
         cos_max_min_angle, flag_zero);
      break;

    case 3:
    default:
      // Edge (v3,v0) is split.
      IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
        (dimension, vertex_coord0, vertex_coord1, vertex_coord2, 
         vertex_coord3, vcoordX, vcoordY, max_small_magnitude,
         cos_max_min_angle, flag_zero);
      break;
    }

  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added between v0 and v3.
   *  - Triangulate either with 3 triangles incident on v1, v2 or vcoordX03[],
   *      or with 5 triangles incident on vcoordY[].
   *  @param vcoordX03[] Coordinates of vertex on edge (v0,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param[out] tri_vertex_index 0,1,2,3,4 of vertex used in
   *    triangulation into three triangles.
   *    - If (tri_vertex_index == 4), then triangulate from vcoordX03[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, typename COS_TYPE, typename ITYPE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX03[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, 
   bool & flag_tri5,
   ITYPE & tri_vertex_index,
   bool & flag_zero)
  {
    COS_TYPE cos_min_tri3_angle, cos_min_tri5_angle;
    bool flag_tri3_zero, flag_tri5_zero;

    compute_cos_max_min_quad_tri3_angle_vX03
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX03, max_small_magnitude, cos_min_tri3_angle, tri_vertex_index,
       flag_tri3_zero);
    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX03, vcoordY, max_small_magnitude,
       cos_min_tri5_angle, flag_tri5_zero);

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

      // (tri_vertex_index == 5) represents that all the triangles
      //   are incident on vcoordY[].
      tri_vertex_index = 5;
    }
    else {
      flag_tri5 = false;
      flag_zero = flag_tri3_zero;
      cos_min_angle = cos_min_tri3_angle;
    }

  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added between v0 and v3.
   *  - Triangulate either with 3 triangles incident on v1, v2 or vcoordX03[],
   *     or with 5 triangles incident on vcoordY[].
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  @param vcoordX03[] Coordinates of vertex on edge (v0,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX03[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {
    bool flag_tri5;
    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX03, vcoordY, max_small_magnitude,
       quad_tri.cos_min_triangulation_angle, flag_tri5, 
       quad_tri.tri_vertex_index, quad_tri.flag_zero);

    if (flag_tri5) {
      quad_tri.num_interior_tri_vertices = 1;
      quad_tri.num_triangles = 5;

      // (tri_vertex_index == 5) represents that all the triangles
      //   are incident on vcoordY[].
      quad_tri.tri_vertex_index = 5;
    }
    else {
      quad_tri.num_interior_tri_vertices = 0;
      quad_tri.num_triangles = 3;
    }
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex vX is added between v0 and v3.
   *  - Triangulate either with 3 triangles incident on vcoordX03[],
   *      or with 5 triangles incident on vcoordY[].
   *  - Note: Difference between this routine and
   *      compute_cos_max_min_quad_tri3_or_tri5_angle_vX03 (not _vXtri3_)
   *      is that the other routine allows triangulation 
   *      into 3 triangles from vertices other than vcoordX03[].
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename CTYPE3, typename CTYPEX, typename CTYPEY,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX03[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_tri5, bool & flag_zero)
  {
    COS_TYPE cos_min_tri3_angle, cos_min_tri5_angle;
    bool flag_tri3_zero, flag_tri5_zero;

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vcoordX03, vertex_coord0, vertex_coord1, 
       vertex_coord2, vertex_coord3, max_small_magnitude, 
       cos_min_tri3_angle, flag_tri3_zero);
    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX03, vcoordY, max_small_magnitude, 
       cos_min_tri5_angle, flag_tri5_zero);

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
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added between v0 and v3.
   *  - Triangulate either with 3 triangles incident on vcoordX03[],
   *      or with 5 triangles incident on vcoordY[].
   *  @param vcoordX03[] Coordinates of vertex on edge (v0,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param[out] quad_tri Triangulation information.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX03[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {
    bool flag_tri5;
    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX03, vcoordY, max_small_magnitude,
       quad_tri.cos_min_triangulation_angle, flag_tri5, quad_tri.flag_zero);

    if (flag_tri5) {
      quad_tri.num_interior_tri_vertices = 1;
      quad_tri.num_triangles = 5;

      // (tri_vertex_index == 5) represents that all the triangles
      //   are incident on vcoordY[].
      quad_tri.tri_vertex_index = 5;
    }
    else {
      quad_tri.num_interior_tri_vertices = 0;
      quad_tri.num_triangles = 3;
      quad_tri.tri_vertex_index = 4;
    }

  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate either with 3 triangles incident on vcoordX[],
   *      or (split_edge_index+2) mod 4 or (split_edge_index+3) mod 4
   *      or with 5 triangles incident on vcoordY[].
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param[out] tri_vertex_index 0,1,2,3,4 of vertex used in
   *     triangulation into three triangles.
   *    If (tri_vertex_index == 4), then triangulate from vcoordX[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE0 split_edge_index,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, 
   bool & flag_tri5,
   ITYPE1 & tri_vertex_index,
   bool & flag_zero)
  {
    COS_TYPE cos_min_tri3_angle, cos_min_tri5_angle;
    bool flag_tri3_zero, flag_tri5_zero;

    compute_cos_max_min_quad_tri3_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX, split_edge_index, max_small_magnitude, 
       cos_min_tri3_angle, tri_vertex_index, flag_tri3_zero);

    compute_cos_min_quad_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX, vcoordY, split_edge_index, max_small_magnitude,
       cos_min_tri5_angle, flag_tri5_zero);

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

      // (tri_vertex_index == 5) represents that all the triangles
      //   are incident on vcoordY[].
      tri_vertex_index = 5;
    }
    else {
      flag_tri5 = false;
      flag_zero = flag_tri3_zero;
      cos_min_angle = cos_min_tri3_angle;
    }

  }


  /*! 
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate either with 3 triangles incident on vcoordX[],
   *      or (split_edge_index+2) mod 4 or (split_edge_index+3) mod 4
   *      or with 5 triangles incident on vcoordY[].
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {
    bool flag_tri5;
    compute_cos_max_min_quad_tri3_or_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX, vcoordY, split_edge_index, max_small_magnitude,
       quad_tri.cos_min_triangulation_angle, flag_tri5, quad_tri.tri_vertex_index,
       quad_tri.flag_zero);

    if (flag_tri5) {
      quad_tri.num_interior_tri_vertices = 1;
      quad_tri.num_triangles = 5;

      // (tri_vertex_index == 5) represents that all the triangles
      //   are incident on vcoordY[].
      quad_tri.tri_vertex_index = 5;
    }
    else {
      quad_tri.num_interior_tri_vertices = 0;
      quad_tri.num_triangles = 3;
    }
  }


  /*! 
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   * @param vcoordX[] Coordinates of vertex on split edge.
   * @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, typename CTYPEY,
            typename ITYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE * vcoord[4],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    compute_cos_max_min_quad_tri3_or_tri5_angle
      (dimension, vcoord[0], vcoord[1], vcoord[2], vcoord[3],
       vcoordX, vcoordY, split_edge_index, max_small_magnitude, 
       quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Version with list of vertices quad_vert[].
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, typename CTYPEY,
            typename VTYPE, typename ITYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const VTYPE quad_vert[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord[NUM_VERT_PER_QUAD] =
      { vertex_coord_ptr+quad_vert[0]*dimension,
        vertex_coord_ptr+quad_vert[1]*dimension,
        vertex_coord_ptr+quad_vert[2]*dimension,
        vertex_coord_ptr+quad_vert[3]*dimension };

    compute_cos_max_min_quad_tri3_or_tri5_angle
      (dimension, vcoord, vcoordX, vcoordY, split_edge_index, 
       max_small_magnitude, quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate either with 3 triangles incident on vcoordX[],
   *      or with 5 triangles incident on vcoordY[].
   * @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename CTYPE3, typename CTYPEX, typename CTYPEY,
            typename ITYPE, typename MTYPE, typename COS_TYPE>
  void compute_cos_max_min_quad_vXtri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_tri5, bool & flag_zero)
  {
    COS_TYPE cos_min_tri3_angle, cos_min_tri5_angle;
    bool flag_tri3_zero, flag_tri5_zero;

    switch (split_edge_index) {

    case 0:
      // Edge (v0,v1) is split.
      IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
        (dimension, vcoordX, vertex_coord1, vertex_coord2, 
         vertex_coord3, vertex_coord0, max_small_magnitude, 
         cos_min_tri3_angle, flag_tri3_zero);
      break;

    case 1:
      // Edge (v1,v2) is split.
      IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
        (dimension, vcoordX, vertex_coord2, vertex_coord3, 
         vertex_coord0, vertex_coord1, max_small_magnitude, 
         cos_min_tri3_angle, flag_tri3_zero);
      break;

    case 2:
      // Edge (v2,v3) is split.
      IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
        (dimension, vcoordX, vertex_coord3, vertex_coord0, 
         vertex_coord1, vertex_coord2, max_small_magnitude, 
         cos_min_tri3_angle, flag_tri3_zero);
      break;

    default:
    case 3:
      // Edge (v3,v0) is split.
      IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
        (dimension, vcoordX, vertex_coord0, vertex_coord1, 
         vertex_coord2, vertex_coord3, max_small_magnitude, 
         cos_min_tri3_angle, flag_tri3_zero);
      break;

    }

    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX, vcoordY, max_small_magnitude, 
       cos_min_tri5_angle, flag_tri5_zero);

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
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate either with 3 triangles incident on vcoordX[],
   *      or with 5 triangles incident on vcoordY[].
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param[out] quad_tri Triangulation information.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_vXtri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri)
  {
    bool flag_tri5;
    compute_cos_max_min_quad_vXtri3_or_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vcoordX, vcoordY, split_edge_index, max_small_magnitude,
       quad_tri.cos_min_triangulation_angle, flag_tri5, quad_tri.flag_zero);

    if (flag_tri5) {
      quad_tri.num_interior_tri_vertices = 1;
      quad_tri.num_triangles = 5;

      // (tri_vertex_index == 5) represents that all the triangles
      //   are incident on vcoordY[].
      quad_tri.tri_vertex_index = 5;
    }
    else {
      quad_tri.num_interior_tri_vertices = 0;
      quad_tri.num_triangles = 3;
      quad_tri.tri_vertex_index = 4;
    }

  }


  /*! 
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added between v0 and v3.
   *  - Triangulate either with 3 triangles incident on v1, v2 
   *      or vcoordX03[], or with 5 triangles incident on vcoordY[].
   *  @param vcoordX03[] Coordinates of vertex on edge (v0,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param block_diagonal[0] If true, do not allow triangulation 
   *    with diagonal (v0,v2).
   *  @param block_diagonal[1] If true, do not allow triangulation 
   *    with diagonal (v1,v3).
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX03[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {

    if (block_diagonal[0] || block_diagonal[1]) {
      compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
        (dimension, vertex_coord0, vertex_coord1, 
         vertex_coord2, vertex_coord3, vcoordX03, vcoordY, 
         max_small_magnitude, quad_tri);
    }
    else {
      compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
        (dimension, vertex_coord0, vertex_coord1, 
         vertex_coord2, vertex_coord3, vcoordX03, vcoordY, 
         max_small_magnitude, quad_tri);
    }
  }


  /*! 
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate either with 3 triangles incident on vcoordX[],
   *      or (split_edge_index+2) mod 4 or (split_edge_index+3) mod 4
   *      or with 5 triangles incident on vcoordY[].
   *  @param vcoordX[] Coordinates of vertex on edge (v0,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param block_diagonal[0] If true, do not allow triangulation 
   *    with diagonal (v0,v2).
   *  @param block_diagonal[1] If true, do not allow triangulation 
   *    with diagonal (v1,v3).
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {

    if (block_diagonal[0] || block_diagonal[1]) {
      compute_cos_max_min_quad_vXtri3_or_tri5_angle
        (dimension, vertex_coord0, vertex_coord1, 
         vertex_coord2, vertex_coord3, vcoordX, vcoordY, 
         split_edge_index, max_small_magnitude, quad_tri_result);
    }
    else {
      compute_cos_max_min_quad_tri3_or_tri5_angle
        (dimension, vertex_coord0, vertex_coord1, 
         vertex_coord2, vertex_coord3, vcoordX, vcoordY, 
         split_edge_index, max_small_magnitude, quad_tri_result);
    }
  }


  /*! 
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate either with 3 triangles incident on vcoordX[],
   *      or (split_edge_index+2) mod 4 or (split_edge_index+3) mod 4
   *      or with 5 triangles incident on vcoordY[].
   *  @param vcoordX[] Coordinates of vertex on edge (v0,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param block_diagonal[0] If true, do not allow triangulation 
   *    with diagonal (v0,v2).
   *  @param block_diagonal[1] If true, do not allow triangulation 
   *    with diagonal (v1,v3).
   */
  template <typename DTYPE, typename CTYPE,
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle
  (const DTYPE dimension,
   const CTYPE vcoord[4], const CTYPEX vcoordX[], const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    compute_cos_max_min_quad_tri3_or_tri5_angle
      (dimension, vcoord[0], vcoord[1], vcoord[2], vcoord[3],
       vcoordX, vcoordY, split_edge_index, max_small_magnitude,
       block_diagonal, quad_tri_result);
  }

  /*! 
   *  Compute max of min quad triangulation angle where
   *    quad is split into 3 or 5 triangles.
   *  - Version with quad vertices and array of vertex coordinates.
   *  - Quad vertices are:
   *       (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3]),
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Triangulate either with 3 triangles incident on vcoordX[],
   *      or (split_edge_index+2) mod 4 or (split_edge_index+3) mod 4
   *      or with 5 triangles incident on vcoordY[].
   */
  template <typename DTYPE, typename CTYPE, 
            typename CTYPEX, typename CTYPEY, 
            typename VTYPE, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri3_or_tri5_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const CTYPEX vcoordX[], const CTYPEY vcoordY[],
   const VTYPE quad_vert[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord[NUM_VERT_PER_QUAD] =
      { vertex_coord_ptr+quad_vert[0]*dimension,
        vertex_coord_ptr+quad_vert[1]*dimension,
        vertex_coord_ptr+quad_vert[2]*dimension,
        vertex_coord_ptr+quad_vert[3]*dimension };

    compute_cos_max_min_quad_tri3_or_tri5_angle
      (dimension, vcoord, vcoordX, vcoordY, split_edge_index, 
       max_small_magnitude, block_diagonal, quad_tri_result);
  }
 

  /*!
   *  Compute the cosine of the smallest triangle angle
   *    in the triangulation of a pentagon into five triangles and
   *    triangle (vcoord0[], vcoord4[], vcoordY[]).
   *  - Quadrilateral is split into five triangles formed from vcoordY[]
   *    and each of the five pentagon edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the five triangles and triangle (vcoord0[], vcoord4[], vcoordY[]).
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param vcoord0[] Coordinates of pentagon vertex 0.
   *  @param vcoord1[] Coordinates of pentagon vertex 1.
   *  @param vcoord2[] Coordinates of pentagon vertex 2.
   *  @param vcoord3[] Coordinates of pentagon vertex 3.
   *  @param vcoord4[] Coordinates of pentagon vertex 4.
   *  @param vcoordX[] Triangle coordinate.
   *  @param vcoord_star_center[] Coordinates of vertex at center of star.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, 
            typename CTYPE5, typename CTYPE6,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_min_pentagon_tri5_and_triangle_angle
  (const DTYPE dimension, const CTYPE0 * vcoord0, const CTYPE1 * vcoord1, 
   const CTYPE2 * vcoord2, const CTYPE3 * vcoord3, const CTYPE4 * vcoord4,
   const CTYPE5 * vcoordX,
   const CTYPE6 * vcoord_star_center,
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cos_min_pentagon_angle;
    COS_TYPE cos_min_triangle_angle;
    bool flag_pentagon_zero;
    bool flag_triangle_zero;
    
    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4,
       vcoord_star_center, max_small_magnitude, 
       cos_min_pentagon_angle, flag_pentagon_zero);
    comput_cos_min_triangle_angle
      (dimension, vcoord0, vcoord4, vcoordX, max_small_magnitude,
       cos_min_triangle_angle, flag_triangle_zero);

    flag_zero = flag_pentagon_zero && flag_triangle_zero;

    if (cos_min_pentagon_angle > cos_min_triangle_angle) 
      { cos_min_angle = cos_min_pentagon_angle; }
    else
      { cos_min_angle = cos_min_triangle_angle; }
  }


  /// Compute max of min triangulation of pentagon, triangle, pentagon.
  /// - Pentagon A is (v0, v1, v2, v6, v7).
  /// - Triangle B is (v1, v3, v2).
  /// - Pentagon C is (v2, v3, v4, v5, v6).
  template <typename DTYPE, typename CTYPE, typename CTYPEA,
            typename CTYPEC, typename MTYPE, 
            typename COS_TYPE, typename NTYPE0, typename NTYPE1>
  void compute_cos_max_min_pentagon_triangle_pentagon_angle
  (const DTYPE dimension,
   const CTYPE vcoord0[],
   const CTYPE vcoord1[],
   const CTYPE vcoord2[],
   const CTYPE vcoord3[],
   const CTYPE vcoord4[],
   const CTYPE vcoord5[],
   const CTYPE vcoord6[],
   const CTYPE vcoord7[],
   const CTYPEA vcoordA[],
   const CTYPEC vcoordC[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, 
   bool & flag_tri5_pentagonA, 
   bool & flag_tri5_pentagonC, 
   NTYPE0 & ivA_max_min, 
   NTYPE1 & ivC_max_min,
   bool & flag_zero)
  {
    COS_TYPE cos_min_triangleB, cos_min_pentagonA, cos_min_pentagonC;
    bool flagA_zero, flagB_zero, flagC_zero;

    // Compute cos min triangle angle.
    compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord3, vcoord2, max_small_magnitude, 
       cos_min_triangleB, flagB_zero);

    // Compute cos min pentagon angles.
    compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord6, vcoord7, vcoordA,
       max_small_magnitude, cos_min_pentagonA, flag_tri5_pentagonA, 
       ivA_max_min, flagA_zero);

    compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vcoord2, vcoord3, vcoord4, vcoord5, vcoord6, vcoordC,
       max_small_magnitude, cos_min_pentagonC, flag_tri5_pentagonC, 
       ivC_max_min, flagC_zero);

    flag_zero = (flagA_zero || flagB_zero || flagC_zero);
    cos_min_angle = std::max(cos_min_pentagonA, cos_min_pentagonC);
    cos_min_angle = std::max(cos_min_angle, cos_min_triangleB);
  }


  ///@}


  // **********************************************************************
  //! @name COMPUTE COS MIN TRI ANGLE - SPLIT EDGE QUADS
  // **********************************************************************

  /*!
   *  Compute min quad triangulation angle where adjacent
   *    quad edges (v0,v1) and (v1,v2) are_split.
   *  - Split into 6 triangles around internal vertex at vcoordY[].
   *  - Vertices are added on split edges (v0,v1) and (v1,v2).
   *  @param vcoordX01[] Coordinates of vertex on edge (v0,v1).
   *  @param vcoordX12[] Coordinates of vertex on edge (v1,v2).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_min_quad_splitL_tri6_angle_vX01_vX12
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    compute_cos_min_hexagon_tri6_angle
      (dimension, vertex_coord0, vcoordX01, vertex_coord1, vcoordX12,
       vertex_coord2, vertex_coord3, vcoordY, max_small_magnitude, 
       quad_tri_result.cos_min_triangulation_angle, 
       quad_tri_result.flag_zero);

    quad_tri_result.num_interior_tri_vertices = 1;
    quad_tri_result.num_triangles = 6;
    quad_tri_result.tri_vertex_index = 6;

  }


  /*!
   *  Compute min quad triangulation angle where
   *    quad edges (v0,v1), (v1,v2) and (v2,v3) are_split.
   *  - Split into 7 triangles around internal vertex at vcoordY[].
   *  - Vertices are added on split edges (v0,v1), (v1,v2) and (v2,v3).
   *  @param vcoordX01[] Coordinates of vertex on edge (v0,v1).
   *  @param vcoordX12[] Coordinates of vertex on edge (v1,v2).
   *  @param vcoordX23[] Coordinates of vertex on edge (v2,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_min_quad_split3e_tri7_angle_vX01_vX12_vX23
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEX vcoordX23[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    COS_TYPE cos_min_pentagon_tri5, cos_min_tri7;
    COS_TYPE cos_min_triangle1_angle, cos_min_triangle2_angle;
    COS_TYPE cos_min_2ear_angle;
    COS_TYPE cos_min_angle;
    bool flag_zero_pentagon_tri5, flag_zero_tri7;
    bool flag_zero_triangle1, flag_zero_triangle2;
    bool flag_zero_2ear;
    bool flag_zero;

    compute_cos_min_septagon_tri7_angle
      (dimension, vertex_coord0, vcoordX01, vertex_coord1, vcoordX12,
       vertex_coord2, vcoordX23, vertex_coord3, vcoordY, 
       max_small_magnitude, cos_min_tri7, flag_zero_tri7);

    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vertex_coord0, vcoordX01, vcoordX12, vcoordX23, 
       vertex_coord3, vcoordY, max_small_magnitude, 
       cos_min_pentagon_tri5, flag_zero_pentagon_tri5);

    compute_cos_min_triangle_angle
      (dimension, vcoordX01, vertex_coord1, vcoordX12, max_small_magnitude,
       cos_min_triangle1_angle, flag_zero_triangle1);

    compute_cos_min_triangle_angle
      (dimension, vcoordX12, vertex_coord2, vcoordX23, max_small_magnitude,
       cos_min_triangle2_angle, flag_zero_triangle2);

    cos_min_2ear_angle = 
      std::max(cos_min_pentagon_tri5, cos_min_triangle1_angle);
    cos_min_2ear_angle = 
      std::max(cos_min_2ear_angle, cos_min_triangle2_angle);
    flag_zero_2ear = (flag_zero_pentagon_tri5 || flag_zero_triangle1 ||
                      flag_zero_triangle2);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "cos_min_tri7: " << cos_min_tri7 << endl;
    cerr << "cos_min_pentagon_tri5: " << cos_min_pentagon_tri5 << endl;
    cerr << "cos_min_triangle1_angle: " << cos_min_triangle1_angle << endl;
    cerr << "cos_min_triangle2_angle: " << cos_min_triangle2_angle << endl;
    cerr << "cos_min_2ear_angle: " << cos_min_2ear_angle << endl;
    */

    int index;
    select_min
      (cos_min_2ear_angle, flag_zero_2ear,
       cos_min_tri7, flag_zero_tri7,
       cos_min_angle, index, flag_zero);

    quad_tri_result.cos_min_triangulation_angle = cos_min_angle;
    quad_tri_result.flag_zero = flag_zero;
    quad_tri_result.num_split_edges = 3;

    switch(index) {

    case 0:
      quad_tri_result.num_interior_tri_vertices = 1;
      quad_tri_result.num_triangles = 7;
      quad_tri_result.tri_vertex_index = 7;
      quad_tri_result.num_tri_ears = 2;
      quad_tri_result.ear_list[0] = 1;
      quad_tri_result.ear_list[1] = 2;
      break;

    default:
      quad_tri_result.num_interior_tri_vertices = 1;
      quad_tri_result.num_triangles = 7;
      quad_tri_result.tri_vertex_index = 7;
      quad_tri_result.num_tri_ears = 0;
      break;
    }
  }



  // **********************************************************************
  //! @name COMPUTE COS MAX MIN TRI ANGLE - SPLIT EDGE QUADS OR TRIANGLES
  // **********************************************************************

  ///@{

  /*!
   *  Compute max of min quad triangulation angle where quad edges
   *    can be split once or twice.
   *  - Vertex is added between v0 and v3 or between v2 and v3.
   *  @param vcoordX03[] Coordinates of vertex on edge (v0,v3).
   *  @param vcoordX23[] Coordinates of vertex on edge (v2,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param[out] tri_vertex_index 0,1,2,3,4,5 of vertex used 
   *    in triangulation.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, typename COS_TYPE, typename NTYPE>
  void compute_cos_max_min_quad_split2_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX03[],
   const CTYPEX vcoordX23[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   COS_TYPE & cos_min_angle,
   NTYPE & num_triangles,
   NTYPE & tri_vertex_index,
   NTYPE & num_tri_ears,
   bool & flag_zero)
  {
    COS_TYPE cos_min_tri6_angle, cos_min_triX_angle;
    COS_TYPE cos_min_v0tri3_angle, cos_min_v2tri3_angle;
    COS_TYPE cos_min_tri5_angle, cos_min_triX_tri5_angle;
    bool flag_zero_tri6;
    bool flag_zero_triX, flag_zero_tri5, flag_zero_triX_tri5;
    bool flag_zero_v0tri3, flag_zero_v2tri3;
    

    compute_cos_min_hexagon_tri6_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vcoordX23, 
       vertex_coord3, vcoordX03, vcoordY, max_small_magnitude, 
       cos_min_tri6_angle, flag_zero_tri6);

    compute_cos_min_triangle_angle
      (dimension, vertex_coord3, vcoordX03, vcoordX23, max_small_magnitude, 
       cos_min_triX_angle, flag_zero_triX);

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vcoordX23, vcoordX03, max_small_magnitude, 
       cos_min_v0tri3_angle, flag_zero_v0tri3);

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vertex_coord2, vcoordX23, vcoordX03,
       vertex_coord0, vertex_coord1, max_small_magnitude, 
       cos_min_v2tri3_angle, flag_zero_v2tri3);

    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vcoordX23, vcoordX03, vcoordY, max_small_magnitude, 
       cos_min_tri5_angle, flag_zero_tri5);

    cos_min_v0tri3_angle = std::max(cos_min_v0tri3_angle, cos_min_triX_angle);
    cos_min_v2tri3_angle = std::max(cos_min_v2tri3_angle, cos_min_triX_angle);
    cos_min_triX_tri5_angle = std::max(cos_min_tri5_angle, cos_min_triX_angle);
    flag_zero_v0tri3 = (flag_zero_v0tri3 || flag_zero_triX);
    flag_zero_v2tri3 = (flag_zero_v2tri3 || flag_zero_triX);
    flag_zero_triX_tri5 = (flag_zero_tri5 || flag_zero_triX);

    NTYPE index;
    if (block_diagonal[0]) {
      select_min
        (cos_min_tri6_angle, flag_zero_tri6,
         cos_min_triX_tri5_angle, flag_zero_triX_tri5,
         cos_min_angle, index, flag_zero);
    }
    else {
      select_minIV
        (cos_min_tri6_angle, flag_zero_tri6,
         cos_min_triX_tri5_angle, flag_zero_triX_tri5,
         cos_min_v0tri3_angle, flag_zero_v0tri3,
         cos_min_v2tri3_angle, flag_zero_v2tri3,
         cos_min_angle, index, flag_zero);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "  vcoord2: ", vertex_coord2, "\n");
    IJK::print_coord3D(cerr, "  vcoordX23: ", vcoordX23, "\n");
    IJK::print_coord3D(cerr, "  vcoordX03: ", vcoordX03, "\n");
    IJK::print_coord3D(cerr, "  vcoord0: ", vertex_coord0, "\n");
    IJK::print_coord3D(cerr, "  vcoord1: ", vertex_coord1, "\n");
    cerr << "  cos_min_triX_angle: " << cos_min_triX_angle << endl;
    cerr << "  cos_min_tri6_angle: " << cos_min_tri6_angle << endl;
    cerr << "  cos_min_v0tri3_angle: " << cos_min_v0tri3_angle << endl;
    cerr << "  cos_min_v2tri3_angle: " << cos_min_v2tri3_angle << endl;
    cerr << "  cos_min_triX_tri5_angle: "
         << cos_min_triX_tri5_angle << endl;
    cerr << "  index: " << index << endl;
    */


    switch(index) {

    case 1:
      num_triangles = 6;
      tri_vertex_index = 0;
      num_tri_ears = 1;
      break;

    case 2:
      num_triangles = 4;
      tri_vertex_index = 0;
      num_tri_ears = 2;
      break;

    case 3:
      num_triangles = 4;
      tri_vertex_index = 2;
      num_tri_ears = 2;
      break;

    default:
      // Triangulate with six triangles incident on vcoordY[].
      num_triangles = 6;
      tri_vertex_index = 0;
      num_tri_ears = 0;
    }

  }


  /*!
   *  Compute max of min quad triangulation angle where quad edges
   *    can be split once or twice.
   *  - Vertex is added between v0 and v3 or between v2 and v3.
   *  @param vcoordX03[] Coordinates of vertex on edge (v0,v3).
   *  @param vcoordX23[] Coordinates of vertex on edge (v2,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split2_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX03[],
   const CTYPEX vcoordX23[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    compute_cos_max_min_quad_split2_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, 
       vertex_coord3, vcoordX03, vcoordX23, vcoordY, max_small_magnitude,
       block_diagonal, 
       quad_tri_result.cos_min_triangulation_angle, 
       quad_tri_result.num_triangles,
       quad_tri_result.tri_vertex_index, 
       quad_tri_result.num_tri_ears,
       quad_tri_result.flag_zero);

    if (quad_tri_result.num_triangles == 5 || 
        quad_tri_result.num_triangles == 6) 
      { quad_tri_result.num_interior_tri_vertices = 1; }
    else 
      { quad_tri_result.num_interior_tri_vertices = 0; }

    quad_tri_result.num_split_edges = 2;
  }


  /*!
   *  Compute max of min quad triangulation angle where adjacent
   *    quad edges (v0,v1) and (v1,v2) are_split.
   *  - Split into 6 triangles around internal vertex at vcoordY[].
   *  - Vertices are added on split edges (v0,v1) and (v1,v2).
   *  @param vcoordX01[] Coordinates of vertex on edge (v0,v1).
   *  @param vcoordX12[] Coordinates of vertex on edge (v1,v2).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_splitL_tri6_angle_vX01_vX12
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    COS_TYPE cos_min_tri6_angle, cos_min_triX_angle;
    COS_TYPE cos_min_tri5_angle, cos_min_triX_tri5_angle;
    bool flag_zero_tri6;
    bool flag_zero_triX, flag_zero_tri5, flag_zero_triX_tri5;

    compute_cos_min_hexagon_tri6_angle
      (dimension, vertex_coord0, vcoordX01, vertex_coord1, vcoordX12,
       vertex_coord2, vertex_coord3, vcoordY, max_small_magnitude, 
       cos_min_tri6_angle, flag_zero_tri6);

    compute_cos_min_triangle_angle
      (dimension, vertex_coord1, vcoordX12, vcoordX01, max_small_magnitude, 
       cos_min_triX_angle, flag_zero_triX);

    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vertex_coord0, vcoordX01, vcoordX12, vertex_coord2, 
       vertex_coord3, vcoordY, max_small_magnitude, 
       cos_min_tri5_angle, flag_zero_tri5);

    cos_min_triX_tri5_angle = std::max(cos_min_tri5_angle, cos_min_triX_angle);
    flag_zero_triX_tri5 = (flag_zero_tri5 || flag_zero_triX);

    NTYPE index;
    select_min
      (cos_min_tri6_angle, flag_zero_tri6,
       cos_min_triX_tri5_angle, flag_zero_triX_tri5,
       quad_tri_result.cos_min_triangulation_angle, index, 
       quad_tri_result.flag_zero);

    quad_tri_result.num_interior_tri_vertices = 1;
    quad_tri_result.num_triangles = 6;
    quad_tri_result.tri_vertex_index = 6;

    switch(index) {

    case 1:
      // Triangulate with five triangles incident on vcoordY[]
      //   and triangle (vcoordX01, v1, vcoordX12).
      quad_tri_result.num_tri_ears = 1;
      quad_tri_result.ear_list[0] = 1;
      break;

    case 0:
    default:
      // Triangulate with six triangles incident on vcoordY[].
      quad_tri_result.num_tri_ears = 0;
    }
  }


  /*!
   *  Compute max of min quad triangulation angle where adjacent
   *    quad edges (v0,v1) and (v1,v2) are split.
   *  - Vertices are added on split edges (v0,v1) and (v1,v2).
   *  - Split into 5 triangles.
   *  @param vcoordX01[] Coordinates of vertex on edge (v0,v1).
   *  @param vcoordX12[] Coordinates of vertex on edge (v1,v2).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_splitL_tri5_angle_vX01_vX12
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    COS_TYPE cos_min_triX_angle;
    COS_TYPE cos_min_v0tri3_angle, cos_min_v2tri3_angle;
    bool flag_zero_triX;
    bool flag_zero_v0tri3, flag_zero_v2tri3;

    compute_cos_min_triangle_angle
      (dimension, vertex_coord1, vcoordX12, vcoordX01, max_small_magnitude, 
       cos_min_triX_angle, flag_zero_triX);

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vertex_coord0, vcoordX01, vcoordX12, 
       vertex_coord2, vertex_coord3, max_small_magnitude, 
       cos_min_v0tri3_angle, flag_zero_v0tri3);

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vertex_coord2, vertex_coord3, vertex_coord0, 
       vcoordX01, vcoordX12, max_small_magnitude, 
       cos_min_v2tri3_angle, flag_zero_v2tri3);

    cos_min_v0tri3_angle = std::max(cos_min_v0tri3_angle, cos_min_triX_angle);
    cos_min_v2tri3_angle = std::max(cos_min_v2tri3_angle, cos_min_triX_angle);
    flag_zero_v0tri3 = (flag_zero_v0tri3 || flag_zero_triX);
    flag_zero_v2tri3 = (flag_zero_v2tri3 || flag_zero_triX);

    NTYPE index;
    select_min
      (cos_min_v0tri3_angle, flag_zero_v0tri3,
       cos_min_v2tri3_angle, flag_zero_v2tri3,
       quad_tri_result.cos_min_triangulation_angle, index, 
       quad_tri_result.flag_zero);


    quad_tri_result.num_interior_tri_vertices = 0;
    quad_tri_result.num_triangles = 4;
    quad_tri_result.num_tri_ears = 2;
    quad_tri_result.ear_list[0] = 1;
    quad_tri_result.ear_list[1] = 3;
    
    switch(index) {

    case 0:
    default:
      quad_tri_result.tri_vertex_index = 0;
      break;

    case 1:
      quad_tri_result.tri_vertex_index = 2;
      break;
    }
  }


  /*!
   *  Compute max of min quad triangulation angle where adjacent
   *    quad edges (v0,v1) and (v1,v2) are split.
   *  - Vertices are added on split edges (v0,v1) and (v1,v2).
   *  - Consider triangulation into five and into six triangles.
   *  @param vcoordX01[] Coordinates of vertex on edge (v0,v1).
   *  @param vcoordX12[] Coordinates of vertex on edge (v1,v2).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_splitL_angle_vX01_vX12
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quad_tri5_result;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quad_tri6_result;
    COS_TYPE cos_max_min_angle;
    bool flag_zero;

    compute_cos_max_min_quad_splitL_tri5_angle_vX01_vX12
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vcoordX01, vcoordX12, vcoordY, max_small_magnitude,
       quad_tri5_result);

    compute_cos_max_min_quad_splitL_tri6_angle_vX01_vX12
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vcoordX01, vcoordX12, vcoordY, max_small_magnitude,
       quad_tri6_result);

    NTYPE index;
    select_min
      (quad_tri5_result.cos_min_triangulation_angle, quad_tri5_result.flag_zero,
       quad_tri6_result.cos_min_triangulation_angle, quad_tri6_result.flag_zero,
       cos_max_min_angle, index, flag_zero);

    if (index == 0) 
      { quad_tri_result = quad_tri5_result; }
    else
      { quad_tri_result = quad_tri6_result; }
  }


  /*!
   *  Compute max of min quad triangulation angle where two adjacent
   *    quad edges are split.
   *  - Vertices are added on split edges.
   *  @param vcoordXA[] Coordinates of vertex on first split edge.
   *  @param vcoordXB[] Coordinates of vertex on second split edge.
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param split_edgeA_index Index of first split edge.
   *    - Second split edge has index (split_edgeA_index+1) mod 4.
   *    - Split edge A is (split_edgeA_index, (split_edgeA_index+1) mod 4).
   *    - Split edge B is 
   *        ((split_edgeA_index+2) mod 4, (split_edgeA_index+3) mod 4).
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename ITYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_splitL_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordXA[],
   const CTYPEX vcoordXB[],
   const CTYPEY vcoordY[],
   const ITYPE split_edgeA_index,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    const ITYPE NUM_VERT_PER_QUAD(4);
    const bool block_diagonal_0_or_1 = 
      (block_diagonal[0] || block_diagonal[1]);

    switch (split_edgeA_index) {

    case 0:
    default:
      // Edges (v0,v1) and (v1,v2) are split.
      if (block_diagonal_0_or_1) {
        compute_cos_max_min_quad_splitL_tri6_angle_vX01_vX12
          (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
           vertex_coord3, vcoordXA, vcoordXB, vcoordY,
           max_small_magnitude, quad_tri_result);
      }
      else {
        compute_cos_max_min_quad_splitL_angle_vX01_vX12
          (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
           vertex_coord3, vcoordXA, vcoordXB, vcoordY,
           max_small_magnitude, quad_tri_result);
      }
      break;


    case 1:
      // Edges (v1,v2) and (v2,v3) are split.
      if (block_diagonal_0_or_1) {
        compute_cos_max_min_quad_splitL_tri6_angle_vX01_vX12
          (dimension, vertex_coord1, vertex_coord2,
           vertex_coord3, vertex_coord0, vcoordXA, vcoordXB, vcoordY,
           max_small_magnitude, quad_tri_result);
      }
      else {
        compute_cos_max_min_quad_splitL_angle_vX01_vX12
          (dimension, vertex_coord1, vertex_coord2,
           vertex_coord3, vertex_coord1, vcoordXA, vcoordXB, vcoordY,
           max_small_magnitude, quad_tri_result);
      }
      break;

    case 2:
      // Edges (v2,v3) and (v3,v0) are split.
      if (block_diagonal_0_or_1) {
        compute_cos_max_min_quad_splitL_tri6_angle_vX01_vX12
          (dimension, vertex_coord2, vertex_coord3, vertex_coord0, 
           vertex_coord1, vcoordXA, vcoordXB, vcoordY,
           max_small_magnitude, quad_tri_result);
      }
      else {
        compute_cos_max_min_quad_splitL_angle_vX01_vX12
          (dimension, vertex_coord2, vertex_coord3, vertex_coord0,
           vertex_coord1, vcoordXA, vcoordXB, vcoordY,
           max_small_magnitude, quad_tri_result);
      }
      break;

    case 3:
      // Edges (v3,v0) and (v0,v2) are split.
      if (block_diagonal_0_or_1) {
        compute_cos_max_min_quad_splitL_tri6_angle_vX01_vX12
          (dimension, vertex_coord3, vertex_coord0, vertex_coord1,
           vertex_coord2, vcoordXA, vcoordXB, vcoordY,
           max_small_magnitude, quad_tri_result);
      }
      else {
        compute_cos_max_min_quad_splitL_angle_vX01_vX12
          (dimension, vertex_coord3, vertex_coord0, vertex_coord1,
           vertex_coord2, vcoordXA, vcoordXB, vcoordY,
           max_small_magnitude, quad_tri_result);
      }
      break;
    }

    // Adjust quad_tri_result.tri_vertex_index based on split_edgeA_index.
    if (quad_tri_result.tri_vertex_index < NUM_VERT_PER_QUAD) {
      quad_tri_result.tri_vertex_index = 
        (quad_tri_result.tri_vertex_index+split_edgeA_index)%NUM_VERT_PER_QUAD;
    }

  }


  /*!
   *  Compute max of min quad triangulation angle where two adjacent
   *    quad edges are split.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename DTYPE, typename CTYPE, 
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_splitL_angle
  (const DTYPE dimension,
   const CTYPE vertex_coord[4],
   const CTYPEX vcoordXA[],
   const CTYPEX vcoordXB[],
   const CTYPEY vcoordY[],
   const ITYPE split_edgeA_index,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    compute_cos_max_min_quad_splitL_angle
      (dimension, vertex_coord[0], vertex_coord[1], vertex_coord[2],
       vertex_coord[3], vcoordXA, vcoordXB, vcoordY, split_edgeA_index,
       max_small_magnitude, block_diagonal, quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where two adjacent
   *    quad edges are split.
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename DTYPE, typename CTYPE, 
            typename CTYPEX, typename CTYPEY, 
            typename VTYPE, typename ITYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_splitL_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const CTYPEX vcoordXA[], const CTYPEX vcoordXB[], const CTYPEY vcoordY[],
   const VTYPE quad_vert[],
   const ITYPE split_edgeA_index,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[2],
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord[NUM_VERT_PER_QUAD] =
      { vertex_coord_ptr+quad_vert[0]*dimension,
        vertex_coord_ptr+quad_vert[1]*dimension,
        vertex_coord_ptr+quad_vert[2]*dimension,
        vertex_coord_ptr+quad_vert[3]*dimension };

    compute_cos_max_min_quad_splitL_angle
      (dimension, vcoord, vcoordXA, vcoordXB, vcoordY, split_edgeA_index,
       max_small_magnitude, block_diagonal, quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    quad edges (v0,v1), (v1,v2) and (v2,v3) are_split.
   *  - Split into 7 triangles around internal vertex at vcoordY[].
   *  - Vertices are added on split edges (v0,v1), (v1,v2) and (v2,v3).
   *  @param vcoordX01[] Coordinates of vertex on edge (v0,v1).
   *  @param vcoordX12[] Coordinates of vertex on edge (v1,v2).
   *  @param vcoordX23[] Coordinates of vertex on edge (v2,v3).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split3e_tri7_angle_vX01_vX12_vX23
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEX vcoordX23[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE, NTYPE> & 
   quad_tri_result)
  {
    COS_TYPE cos_min_tri7_angle;
    bool flag_zero_tri7;

    compute_cos_min_septagon_tri7_angle
      (dimension, vertex_coord0, vcoordX01, vertex_coord1, vcoordX12,
       vertex_coord2, vcoordX23, vertex_coord3, vcoordY, 
       max_small_magnitude, cos_min_tri7_angle, flag_zero_tri7);

    quad_tri_result.cos_min_triangulation_angle = cos_min_tri7_angle;
    quad_tri_result.flag_zero = flag_zero_tri7;
    quad_tri_result.num_interior_tri_vertices = 1;
    quad_tri_result.num_triangles = 7;
    quad_tri_result.tri_vertex_index = 7;
    quad_tri_result.num_tri_ears = 0;
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    three quad edges are_split.
   *  - Split into 7 triangles around internal vertex at vcoordY[].
   *  - Vertices are added on split edges.
   *  @param vcoordXA[] Coordinates of vertex on first split edge.
   *  @param vcoordXB[] Coordinates of vertex on second split edge.
   *  @param vcoordXC[] Coordinates of vertex on third split edge.
   *  @param vcoordY[] Coordinates of vertex at center of star.
   *  @param splitA_edge_index Index of first split edge.
   *  - Second split edge is (splitA_edge_index+1) mod 4.
   *  - Third split edge is (splitA_edge_index+2) mod 4.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split3e_tri7_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordXA[],
   const CTYPEX vcoordXB[],
   const CTYPEX vcoordXC[],
   const CTYPEY vcoordY[],
   const ITYPE split_edgeA_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    const ITYPE NUM_VERT_PER_QUAD(4);

    switch (split_edgeA_index) {

    case 0:
    default:
      compute_cos_max_min_quad_split3e_tri7_angle_vX01_vX12_vX23
        (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
         vertex_coord3, vcoordXA, vcoordXB, vcoordXC, vcoordY,
         max_small_magnitude, quad_tri_result);
      break;

    case 1:
      compute_cos_max_min_quad_split3e_tri7_angle_vX01_vX12_vX23
        (dimension, vertex_coord1, vertex_coord2, vertex_coord3,
         vertex_coord0, vcoordXA, vcoordXB, vcoordXC, vcoordY,
         max_small_magnitude, quad_tri_result);
      break;

    case 2:
      compute_cos_max_min_quad_split3e_tri7_angle_vX01_vX12_vX23
        (dimension, vertex_coord2, vertex_coord3, vertex_coord0,
         vertex_coord1, vcoordXA, vcoordXB, vcoordXC, vcoordY,
         max_small_magnitude, quad_tri_result);
      break;

    case 3:
      compute_cos_max_min_quad_split3e_tri7_angle_vX01_vX12_vX23
        (dimension, vertex_coord3, vertex_coord0, vertex_coord1,
         vertex_coord2, vcoordXA, vcoordXB, vcoordXC, vcoordY,
         max_small_magnitude, quad_tri_result);
      break;
    }

    // Adjust quad_tri_result.tri_vertex_index based on split_edgeA_index.
    if (quad_tri_result.tri_vertex_index < NUM_VERT_PER_QUAD) {
      quad_tri_result.tri_vertex_index = 
        (quad_tri_result.tri_vertex_index+split_edgeA_index)%NUM_VERT_PER_QUAD;
    }
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    three quad edges are_split.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename DTYPE, typename CTYPE,
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split3e_tri7_angle
  (const DTYPE dimension,
   const CTYPE vertex_coord[4],
   const CTYPEX vcoordXA[],
   const CTYPEX vcoordXB[],
   const CTYPEX vcoordXC[],
   const CTYPEY vcoordY[],
   const ITYPE split_edgeA_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    compute_cos_max_min_quad_split3e_tri7_angle
      (dimension, vertex_coord[0], vertex_coord[1], vertex_coord[2],
       vertex_coord[3], vcoordXA, vcoordXB, vcoordXC, vcoordY,
       split_edgeA_index, max_small_magnitude, quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    three quad edges are_split.
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename DTYPE, typename CTYPE,
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split3e_tri7_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const CTYPEX vcoordXA[], const CTYPEX vcoordXB[], 
   const CTYPEX vcoordXC[], const CTYPEY vcoordY[],
   const VTYPE quad_vert[], const ITYPE split_edgeA_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord[NUM_VERT_PER_QUAD] =
      { vertex_coord_ptr+quad_vert[0]*dimension,
        vertex_coord_ptr+quad_vert[1]*dimension,
        vertex_coord_ptr+quad_vert[2]*dimension,
        vertex_coord_ptr+quad_vert[3]*dimension };

    compute_cos_max_min_quad_split3e_tri7_angle
      (dimension, vcoord, vcoordXA, vcoordXB, vcoordXC, vcoordY, 
       split_edgeA_index, max_small_magnitude, quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    three quad edges are_split.
   *  - Version using array quad_vert[] of quad vertices.
   *  - Version storing vcoordXA[], vcoordXB[] and vcoordXC[] 
   *      in a single array.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param vcoordX[] Array of vertex coordinates.
   *  - vcoordXA[] starts at is vcoordX.
   *  - vcoordXB[] starts at vcoordX+dimension.
   *  - vcoordXC[] starts at vcoordX+2*dimension.
   */
  template <typename DTYPE, typename CTYPE,
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split3e_tri7_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const CTYPEX vcoordX[], const CTYPEY vcoordY[],
   const VTYPE quad_vert[], const ITYPE split_edgeA_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    const CTYPEX * vcoordXA = vcoordX;
    const CTYPEX * vcoordXB = vcoordX+dimension;
    const CTYPEX * vcoordXC = vcoordX+2*dimension;

    compute_cos_max_min_quad_split3e_tri7_angle
      (dimension, vertex_coord, vcoordXA, vcoordXB, vcoordXC,
       vcoordY, quad_vert, split_edgeA_index, max_small_magnitude,
       quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    all quad edges are split.
   *  - Split into 8 triangles around internal vertex at vcoordY[].
   *  - Vertices are added on each quad edge.
   *  @param vcoordX01[] Coordinates of vertex on edge (v0,v1).
   *  @param vcoordX12[] Coordinates of vertex on edge (v1,v2).
   *  @param vcoordX23[] Coordinates of vertex on edge (v2,v3).
   *  @param vcoordX30[] Coordinates of vertex on edge (v3,v0).
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split4e_tri8_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEX vcoordX23[],
   const CTYPEX vcoordX30[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    COS_TYPE cos_min_tri8_angle;
    bool flag_zero_tri8;

    compute_cos_min_octagon_tri8_angle
      (dimension, vertex_coord0, vcoordX01, vertex_coord1, vcoordX12,
       vertex_coord2, vcoordX23, vertex_coord3, vcoordX30, vcoordY, 
       max_small_magnitude, cos_min_tri8_angle, flag_zero_tri8);

    quad_tri_result.cos_min_triangulation_angle = cos_min_tri8_angle;
    quad_tri_result.flag_zero = flag_zero_tri8;
    quad_tri_result.num_interior_tri_vertices = 1;
    quad_tri_result.num_triangles = 8;
    quad_tri_result.tri_vertex_index = 8;
    quad_tri_result.num_tri_ears = 0;
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    all quad edges are split.
   *  - Split into 8 triangles around internal vertex at vcoordY[].
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename DTYPE, 
            typename CTYPE, typename CTYPEX, typename CTYPEY,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split4e_tri8_angle
  (const DTYPE dimension,
   const CTYPE vertex_coord[4],
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEX vcoordX23[],
   const CTYPEX vcoordX30[],
   const CTYPEY vcoordY[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    compute_cos_max_min_quad_split4e_tri8_angle
      (dimension, vertex_coord[0], vertex_coord[1], vertex_coord[2],
       vertex_coord[3], vcoordX01, vcoordX12, vcoordX23, vcoordX30,
       vcoordY, max_small_magnitude, quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    all quad edges are split.
   *  - Split into 8 triangles around internal vertex at vcoordY[].
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename DTYPE, 
            typename CTYPE, typename CTYPEX, typename CTYPEY,
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split4e_tri8_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord, 
   const CTYPEX vcoordX01[],
   const CTYPEX vcoordX12[],
   const CTYPEX vcoordX23[],
   const CTYPEX vcoordX30[],
   const CTYPEY vcoordY[],
   const VTYPE quad_vert[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord[NUM_VERT_PER_QUAD] =
      { vertex_coord_ptr+quad_vert[0]*dimension,
        vertex_coord_ptr+quad_vert[1]*dimension,
        vertex_coord_ptr+quad_vert[2]*dimension,
        vertex_coord_ptr+quad_vert[3]*dimension };

    compute_cos_max_min_quad_split4e_tri8_angle
      (dimension, vcoord, vcoordX01, vcoordX12, vcoordX23, vcoordX30,
       vcoordY, max_small_magnitude, quad_tri_result);
  }


  /*!
   *  Compute max of min quad triangulation angle where
   *    all quad edges are split.
   *  - Split into 8 triangles around internal vertex at vcoordY[].
   *  - Version storing vcoordXA[], vcoordXB[] and vcoordXC[] 
   *      in a single array.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param vcoordX[] Array of vertex coordinates.
   *  - vcoordXA[] starts at is vcoordX.
   *  - vcoordXB[] starts at vcoordX+dimension.
   *  - vcoordXC[] starts at vcoordX+2*dimension.
   *  - vcoordXD[] starts at vcoordX+3*dimension.
   */
  template <typename DTYPE, 
            typename CTYPE, typename CTYPEX, typename CTYPEY,
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_split4e_tri8_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord, 
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const VTYPE quad_vert[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    const CTYPEX * vcoordXA = vcoordX;
    const CTYPEX * vcoordXB = vcoordX+dimension;
    const CTYPEX * vcoordXC = vcoordX+2*dimension;
    const CTYPEX * vcoordXD = vcoordX+3*dimension;

    compute_cos_max_min_quad_split4e_tri8_angle
      (dimension, vertex_coord, vcoordXA, vcoordXB, vcoordXC, vcoordXD,
       vcoordY, quad_vert, max_small_magnitude, quad_tri_result);
  }

  ///@}


  // **********************************************************************
  //! @name COMPUTE COS MAX MIN TRI ANGLE - SPLIT EDGE PENTAGON
  // **********************************************************************

  ///@{

  /*! Compute max of min pentagon triangulation angle where pentagon is split
   *    into 4 triangles.
   *  - Vertex is added between v0 and v4.
   *  - Triangulate with 3 triangles incident on v1, v2, v3 or vcoordX04[].
   *  @param vcoordX04[] Coordinates of vertex on edge (v0,v4).
   *  @param[out] tri_vertex_index 0,1,2,3,4,5 of vertex used in
   *     triangulation into three triangles.
   *  If (tri_vertex_index == 5), then triangulate from vcoordX04[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename CTYPEX, typename MTYPE, 
            typename COS_TYPE, typename ITYPE>
  void compute_cos_max_min_pentagon_tri4_angle_vX04
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPEX vcoordX04[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, ITYPE & tri_vertex_index, bool & flag_zero)
  {
    COS_TYPE cos_min_v1_angle, cos_min_v2_angle, cos_min_v3_angle;
    COS_TYPE cos_min_vX04_angle;
    bool flag_zero_v1, flag_zero_v2, flag_zero_v3, flag_zero_vX04;

    compute_cos_min_hexagon_triangulation_angle
      (dimension, vcoordX04, vertex_coord0, vertex_coord1, vertex_coord2, 
       vertex_coord3, vertex_coord4, max_small_magnitude, cos_min_vX04_angle, 
       flag_zero_vX04);
    compute_cos_min_hexagon_triangulation_angle
      (dimension, vertex_coord1, vertex_coord2, vertex_coord3, 
       vertex_coord4, vcoordX04, vertex_coord0, 
       max_small_magnitude, cos_min_v1_angle, flag_zero_v1);
    compute_cos_min_hexagon_triangulation_angle
      (dimension, vertex_coord2, vertex_coord3, vertex_coord4,
       vcoordX04, vertex_coord0, vertex_coord1, 
       max_small_magnitude, cos_min_v2_angle, flag_zero_v2);
    compute_cos_min_hexagon_triangulation_angle
      (dimension, vertex_coord3, vertex_coord4, vcoordX04, 
       vertex_coord0, vertex_coord1, vertex_coord2,
       max_small_magnitude, cos_min_v3_angle, flag_zero_v3);

    int index_selected;
    select_minIV
      (cos_min_v1_angle, flag_zero_v1, 
       cos_min_v2_angle, flag_zero_v2, 
       cos_min_v3_angle, flag_zero_v3, 
       cos_min_vX04_angle, flag_zero_vX04,
       cos_min_angle, index_selected, flag_zero);

    if (index_selected == 0)
      { tri_vertex_index = 1; }
    else if (index_selected == 1)
      { tri_vertex_index = 2; }
    else if (index_selected == 2)
      { tri_vertex_index = 3; }
    else
      { tri_vertex_index = 5; }
  }


  /// Compute max of min pentagon triangulation angle where pentagon is split
  ///   into 4 triangles.
  /// - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 5).
  /// - Triangulate with 4 triangles incident on vcoordX05[]
  ///     or (split_edge_index+2) mod 5 or (split_edge_index+3) mod 5
  ///     or (split_edge_index+4) mod 5.
  /// @param vcoordX[] Coordinates of vertex on split edge.
  /// @param[out] tri_vertex_index 0,1,2,3,4,5 of vertex used in
  ///    triangulation into four triangles.
  ///   If (tri_vertex_index == 5), then triangulate from vcoordX[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename CTYPEX, typename ITYPE0, typename ITYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_max_min_pentagon_tri4_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPEX vcoordX[],
   const ITYPE0 split_edge_index,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_max_min_angle, ITYPE1 & tri_vertex_index, bool & flag_zero)
  {
    const int NUM_VERT_PER_PENTAGON(5);

    switch (split_edge_index) {

    case 0:
      // Edge (v0,v1) is split.
      compute_cos_max_min_pentagon_tri4_angle_vX04
        (dimension, vertex_coord1, vertex_coord2, vertex_coord3, 
         vertex_coord4, vertex_coord0, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      if (tri_vertex_index < NUM_VERT_PER_PENTAGON) {
        tri_vertex_index = (tri_vertex_index+1)%NUM_VERT_PER_PENTAGON;
      }
      break;

    case 1:
      // Edge (v1,v2) is split.
      compute_cos_max_min_pentagon_tri4_angle_vX04
        (dimension, vertex_coord2, vertex_coord3, vertex_coord4,
         vertex_coord0, vertex_coord1, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      if (tri_vertex_index < NUM_VERT_PER_PENTAGON) {
        tri_vertex_index = (tri_vertex_index+2)%NUM_VERT_PER_PENTAGON;
      }
      break;

    case 2:
      // Edge (v2,v3) is split.
      compute_cos_max_min_pentagon_tri4_angle_vX04
        (dimension, vertex_coord3, vertex_coord4, vertex_coord0, 
         vertex_coord1, vertex_coord2, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      if (tri_vertex_index < NUM_VERT_PER_PENTAGON) {
        tri_vertex_index = (tri_vertex_index+3)%NUM_VERT_PER_PENTAGON;
      }
      break;

    case 3:
      // Edge (v3,v4) is split.
      compute_cos_max_min_pentagon_tri4_angle_vX04
        (dimension, vertex_coord4, vertex_coord0, vertex_coord1, 
         vertex_coord2, vertex_coord3, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      if (tri_vertex_index < NUM_VERT_PER_PENTAGON) {
        tri_vertex_index = (tri_vertex_index+4)%NUM_VERT_PER_PENTAGON;
      }
      break;

    case 4:
    default:
      // Edge (v4,v0) is split.
      compute_cos_max_min_pentagon_tri4_angle_vX04
        (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
         vertex_coord3, vertex_coord4, vcoordX, max_small_magnitude,
         cos_max_min_angle, tri_vertex_index, flag_zero);
      return;
    }

  }


  /*!
   *  Compute max of min pentagon triangulation angle where pentagon is split
   *    into 4 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 5).
   *  - Triangulate with 4 triangles incident on vcoordX05[]
   *      or (split_edge_index+2) mod 5 or (split_edge_index+3) mod 5
   *      or (split_edge_index+4) mod 5.
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param[out] tri_vertex_index 0,1,2,3,4,5 of vertex used in
   *     triangulation into four triangles.
   *  If (tri_vertex_index == 5), then triangulate from vcoordX[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename CTYPEX, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_pentagon_tri4_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPEX vcoordX[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   pentagon_tri_result)
  {
    COS_TYPE cos_max_min_angle;
    NTYPE tri_vertex_index;
    bool flag_zero;

    compute_cos_max_min_pentagon_tri4_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, vcoordX, split_edge_index,
       max_small_magnitude, cos_max_min_angle, tri_vertex_index,
       flag_zero);

    pentagon_tri_result.cos_min_triangulation_angle = cos_max_min_angle;
    pentagon_tri_result.tri_vertex_index = tri_vertex_index;
    pentagon_tri_result.flag_zero = flag_zero;
  }


  /*!
   *  Compute max of min pentagon triangulation angle where pentagon is split
   *    into 6 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 5).
   *  - Triangulate with 6 triangles incident on vcoordY[].
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param[out] tri_vertex_index 0,1,2,3,4,5 of vertex used in
   *     triangulation into three triangles.
   *    - If (tri_vertex_index == 5), then triangulate from vcoordX[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename CTYPEX, typename CTYPEY,
            typename ITYPE, typename MTYPE, typename COS_TYPE>
  void compute_cos_min_pentagon_tri6_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_max_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_PENTAGON(5);

    switch (split_edge_index) {

    case 0:
      // Edge (v0,v1) is split.
      compute_cos_min_hexagon_tri6_angle
        (dimension, vertex_coord0, vcoordX, vertex_coord1,
         vertex_coord2, vertex_coord3, vertex_coord4, vcoordY, 
         max_small_magnitude, cos_max_min_angle, flag_zero);
      break;

    case 1:
      // Edge (v1,v2) is split.
      compute_cos_min_hexagon_tri6_angle
        (dimension, vertex_coord0, vertex_coord1, vcoordX,
         vertex_coord2, vertex_coord3, vertex_coord4, vcoordY, 
         max_small_magnitude, cos_max_min_angle, flag_zero);
      break;

    case 2:
      // Edge (v2,v3) is split.
      compute_cos_min_hexagon_tri6_angle
        (dimension, vertex_coord0, vertex_coord1, vertex_coord2, 
         vcoordX, vertex_coord3, vertex_coord4, vcoordY, 
         max_small_magnitude, cos_max_min_angle, flag_zero);
      break;

    case 3:
      // Edge (v3,v4) is split.
      compute_cos_min_hexagon_tri6_angle
        (dimension, vertex_coord0, vertex_coord1, vertex_coord2, 
         vertex_coord3, vcoordX, vertex_coord4, vcoordY, 
         max_small_magnitude, cos_max_min_angle, flag_zero);
      break;

    case 4:
    default:
      // Edge (v4,v0) is split.
      compute_cos_min_hexagon_tri6_angle
        (dimension, vertex_coord0, vertex_coord1, vertex_coord2, 
         vertex_coord3, vertex_coord4, vcoordX, vcoordY, 
         max_small_magnitude, cos_max_min_angle, flag_zero);
      break;
    }

  }


  /*!
   *  Compute max of min pentagon triangulation angle where pentagon is split
   *    into 6 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 5).
   *  - Triangulate with 6 triangles incident on vcoordY[].
   *  - Version returning data structure POLY_TRIANGULATION_RESULT.
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param[out] tri_vertex_index 0,1,2,3,4,5 of vertex used in
   *     triangulation into three triangles.
   *    - If (tri_vertex_index == 5), then triangulate from vcoordX[].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename CTYPEX, typename CTYPEY,
            typename ITYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_min_pentagon_tri6_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   pentagon_tri_result)
  {
    compute_cos_min_pentagon_tri6_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, vcoordX, vcoordY,
       split_edge_index, max_small_magnitude, 
       pentagon_tri_result.cos_min_triangulation_angle, pentagon_tri_result.flag_zero);

    pentagon_tri_result.num_interior_tri_vertices = 1;
    pentagon_tri_result.tri_vertex_index = 6;
  }


  /*!
   *  Compute max of min pentagon triangulation angle where pentagon
   *    is split into 4 or 6 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 5).
   *  - Triangulate either with 4 triangles incident on vcoordX[],
   *      or (split_edge_index+2) mod 5 or (split_edge_index+3) mod 5
   *      or (split_edge_index+3) mod 5
   *      or with 6 triangles incident on vcoordY[].
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_pentagon_tri4_or_tri6_angle
  (const DTYPE dimension,
   const CTYPE0 vertex_coord0[],
   const CTYPE1 vertex_coord1[],
   const CTYPE2 vertex_coord2[],
   const CTYPE3 vertex_coord3[],
   const CTYPE4 vertex_coord4[],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   pentagon_tri_result)
  {
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      pentagon_tri4_result;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      pentagon_tri6_result;

    compute_cos_max_min_pentagon_tri4_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, vcoordX, split_edge_index,
       max_small_magnitude, pentagon_tri4_result);

    compute_cos_min_pentagon_tri6_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, vcoordX, vcoordY, split_edge_index,
       max_small_magnitude, pentagon_tri6_result);

    if (pentagon_tri6_result.cos_min_triangulation_angle <=
        pentagon_tri4_result.cos_min_triangulation_angle) 
      { pentagon_tri_result = pentagon_tri6_result; }
    else 
      { pentagon_tri_result = pentagon_tri4_result; }
  }


  /*!
   *  Compute max of min pentagon triangulation angle where pentagon
   *    is split into 4 or 6 triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 5).
   *  - Version using array vcoord[] of pentagon vertex coordinates.
   *  - Pengaton vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3], vcoord[4])
   *    listed in clockwise/counter-clockwise order around the pentagon.
   *  @param vcoordX[] Coordinates of vertex on split edge.
   *  @param vcoordY[] Coordinates of vertex at center of star.
   */
  template <typename DTYPE, typename CTYPE,
            typename CTYPEX, typename CTYPEY, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_pentagon_tri4_or_tri6_angle
  (const DTYPE dimension, const CTYPE * vcoord[5],
   const CTYPEX vcoordX[],
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   pentagon_tri_result)
  {
    compute_cos_max_min_pentagon_tri4_or_tri6_angle
      (dimension, vcoord[0], vcoord[1], vcoord[2], vcoord[3], vcoord[4],
       vcoordX, vcoordY, split_edge_index, max_small_magnitude,
       pentagon_tri_result);
  }

  /*! 
   *  Compute max of min pentagon triangulation angle where
   *    pentagon is split into 4 or 6 triangles.
   *  - Version with pentagon vertices and array of vertex coordinates.
   *  - Pentagon vertices are:
   *       (pentagon_vert[0], pentagon_vert[1], ..., pentagon_vert[4]),
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 5).
   *  - Triangulate either with 4 triangles incident on vcoordX[],
   *      or with 6 triangles incident on vcoordY[].
   */
  template <typename DTYPE, typename CTYPE, 
            typename CTYPEX, typename CTYPEY, 
            typename VTYPE, typename ITYPE,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_pentagon_tri4_or_tri6_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord, 
   const CTYPEX vcoordX[], const CTYPEY vcoordY[],
   const VTYPE pentagon_vert[5],
   const ITYPE split_edge_index,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   pentagon_tri_result)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord[NUM_VERT_PER_PENTAGON] =
      { vertex_coord_ptr+pentagon_vert[0]*dimension,
        vertex_coord_ptr+pentagon_vert[1]*dimension,
        vertex_coord_ptr+pentagon_vert[2]*dimension,
        vertex_coord_ptr+pentagon_vert[3]*dimension,
        vertex_coord_ptr+pentagon_vert[4]*dimension };

    compute_cos_max_min_pentagon_tri4_or_tri6_angle
      (dimension, vcoord, vcoordX, vcoordY, split_edge_index, 
       max_small_magnitude, pentagon_tri_result);
  }


  // ****************************************************************
  //! @name TRIANGULATE QUADS WITH SPLIT EDGES
  // ****************************************************************

  ///@{

  /*!
   *  Triangulate quad into three or five triangles.
   *  - Edge (v0,v3) is split by vertex ivX03.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *      around the quadrilateral.
   * - Add new triangles to vector tri_vert.
   * @param ivX03 Index of vertex splitting edge (v0,v3).
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPEX, 
            typename NTYPE, typename ITYPE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
 void triangulate_quad_tri3_or_tri5_vX03
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX ivX03, const CTYPEY vcoordY[],
   const NTYPE num_triangles,
   const ITYPE tri_vertex_index,
   const bool flag_reverse_orient,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_PENTAGON(5);
    const VTYPE0 pentagon_vert[NUM_VERT_PER_PENTAGON] = 
      { v0, v1, v2, v3, ivX03 };

    if (num_triangles == 5) {
      IJK_DEPRECATED::triangulate_pentagon_tri5
        (dimension, v0, v1, v2, v3, ivX03, vcoordY, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      IJK::triangulate_pentagon
        (pentagon_vert, tri_vertex_index, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  Triangulate quad into three or five triangles.
   *  - Edge (v0,v3) is split by vertex ivX03.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *      around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  - Version using data structure POLY_TRIANGULATION_RESULT.
   *  @param ivX03 Index of vertex splitting edge (v0,v3).
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPEX, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
 void triangulate_quad_tri3_or_tri5_vX03
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX ivX03, const CTYPEY vcoordY[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, v0, v1, v2, v3, ivX03, vcoordY,
       quad_tri.num_triangles, quad_tri.tri_vertex_index,
       flag_reverse_orient, vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate quad into three or five triangles.
   *  - Edge (v0,v3) is split by vertex ivX03.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *      around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  - Version using data structure POLY_TRIANGULATION_RESULT.
   *  - Version using array quad_vert[].
   *  @param ivX03 Index of vertex splitting edge (v0,v3).
   */
  template <typename DTYPE, typename VTYPE, typename VTYPEX, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
 void triangulate_quad_tri3_or_tri5_vX03
  (const DTYPE dimension, const VTYPE quad_vert[4],
   const VTYPEX ivX03, const CTYPEY vcoordY[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3],
       ivX03, vcoordY, flag_reverse_orient, quad_tri_result, 
       vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate quad into three or five triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *      around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  @param iv_split Index of vertex splitting edge.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPEX, 
            typename ITYPE0, typename ITYPE1, typename NTYPE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
 void triangulate_quad_tri3_or_tri5
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX iv_split, 
   const CTYPEY vcoordY[],
   const ITYPE0 split_edge_index,
   const NTYPE num_triangles,
   const ITYPE1 tri_vertex_index,
   const bool flag_reverse_orient,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);
    ITYPE1 tri_vertex_indexB = tri_vertex_index;

    switch (split_edge_index) {

    case 0:
      if (tri_vertex_index < NUM_VERT_PER_QUAD) 
        { tri_vertex_indexB = (tri_vertex_index+3)%NUM_VERT_PER_QUAD; }
      triangulate_quad_tri3_or_tri5_vX03
        (dimension, v1, v2, v3, v0, iv_split, vcoordY, 
         num_triangles, tri_vertex_indexB,
         flag_reverse_orient, vertex_coord, tri_vert);
      break;

    case 1:
      if (tri_vertex_index < NUM_VERT_PER_QUAD) 
        { tri_vertex_indexB = (tri_vertex_index+2)%NUM_VERT_PER_QUAD; }
      triangulate_quad_tri3_or_tri5_vX03
        (dimension, v2, v3, v0, v1, iv_split, vcoordY, 
         num_triangles, tri_vertex_indexB,
         flag_reverse_orient, vertex_coord, tri_vert);
      break;

    case 2:
      if (tri_vertex_index < NUM_VERT_PER_QUAD) 
        { tri_vertex_indexB = (tri_vertex_index+1)%NUM_VERT_PER_QUAD; }
      triangulate_quad_tri3_or_tri5_vX03
        (dimension, v3, v0, v1, v2, iv_split, vcoordY, 
         num_triangles, tri_vertex_indexB,
         flag_reverse_orient, vertex_coord, tri_vert);
      break;

    default:
    case 3:
      triangulate_quad_tri3_or_tri5_vX03
        (dimension, v0, v1, v2, v3, iv_split, vcoordY, 
         num_triangles, tri_vertex_index,
         flag_reverse_orient, vertex_coord, tri_vert);
      break;
    }

  }


  /*!
   *  Triangulate quad into three or five triangles.
   *  - Vertex is added on edge (split_edge_index,(split_edge_index+1) mod 4).
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *      around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  - Version using data structure POLY_TRIANGULATION_RESULT.
   *  @param iv_split Index of vertex splitting edge.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPEX, 
            typename ITYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
 void triangulate_quad_tri3_or_tri5
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX iv_split, 
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri,
   const bool flag_reverse_orient,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    triangulate_quad_tri3_or_tri5
      (dimension, v0, v1, v2, v3, iv_split, vcoordY, split_edge_index,
       quad_tri.num_triangles, quad_tri.tri_vertex_index,
       flag_reverse_orient,vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate quad into three or five triangles.
   *  - Version using data structure POLY_TRIANGULATION_RESULT.
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  @param iv_split Index of vertex splitting edge.
   */
  template <typename DTYPE, typename VTYPE, typename VTYPEX, 
            typename ITYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
 void triangulate_quad_tri3_or_tri5
  (const DTYPE dimension, const VTYPE quad_vert[4],
   const VTYPEX iv_split, 
   const CTYPEY vcoordY[],
   const ITYPE split_edge_index,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri,
   const bool flag_reverse_orient,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    triangulate_quad_tri3_or_tri5
      (dimension, quad_vert[0], quad_vert[1], quad_vert[2],
       quad_vert[3], iv_split, vcoordY, split_edge_index, quad_tri,
       flag_reverse_orient, vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate quad with split edges (v0,v3) and (v2,v3).
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *      around the quadrilateral.
   * - Add new triangles to vector tri_vert.
   *  @param vcoordY[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPEX, 
            typename NTYPE, typename ITYPE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
  void triangulate_quad_with_split_v03_v23
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX ivX03, const VTYPEX ivX23, const CTYPEY vcoordY[],
   const NTYPE num_triangles,
   const NTYPE num_tri_ears,
   const ITYPE tri_vertex_index,
   const bool flag_reverse_orient,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    const NTYPE NUM_VERT_PER_PENTAGON(5);
    const VTYPE0 pentagon_vert[NUM_VERT_PER_PENTAGON] = 
      { v0, v1, v2, ivX23, ivX03 };

    if (num_triangles == 6) {

      if (num_tri_ears == 0) {
        IJK_DEPRECATED::triangulate_hexagon_tri6
          (dimension, v0, v1, v2, ivX23, v3, ivX03, vcoordY, 
           flag_reverse_orient, vertex_coord, tri_vert);
      }
      else {
        add_triangle_vertices(ivX23, v3, ivX03, flag_reverse_orient, tri_vert);
        IJK_DEPRECATED::triangulate_pentagon_tri5
          (dimension, v0, v1, v2, ivX23, ivX03, vcoordY, 
           flag_reverse_orient, vertex_coord, tri_vert);
      }
    }
    else {
      add_triangle_vertices(ivX23, v3, ivX03, flag_reverse_orient, tri_vert);
      IJK::triangulate_pentagon
        (pentagon_vert, tri_vertex_index, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  Triangulate quad with split edges (v0,v3) and (v2,v3).
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *      around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  - Version using data structure POLY_TRIANGULATION_RESULT.
   *  @param vcoordY[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPEX, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
  void triangulate_quad_with_split_v03_v23
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX ivX03, const VTYPEX ivX23, const CTYPEY vcoordY[],
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri,
   const bool flag_reverse_orient,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    triangulate_quad_with_split_v03_v23
      (dimension, v0, v1, v2, v3, ivX03, ivX23, vcoordY,
       quad_tri.num_triangles, quad_tri.num_tri_ears, 
       quad_tri.tri_vertex_index, flag_reverse_orient, 
       vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate quad with two adjacent split edges.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *      around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordY[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPEX, 
            typename ITYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
  void triangulate_quad_with_two_adjacent_split_edges
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX ivXA, const VTYPEX ivXB, const CTYPEY vcoordY[], 
   const ITYPE split_edge_index,  
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri,
   const bool flag_reverse_orient,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);
    const NTYPE num_triangles = quad_tri.num_triangles;
    const NTYPE num_tri_ears = quad_tri.num_tri_ears;
    ITYPE tri_vertex_index = quad_tri.tri_vertex_index;

    switch(split_edge_index) {
      
    case 0:
      if (tri_vertex_index < NUM_VERT_PER_QUAD)
        { tri_vertex_index = (tri_vertex_index+3)%NUM_VERT_PER_QUAD; }
      triangulate_quad_with_split_v03_v23
        (dimension, v1, v2, v3, v0, ivXA, ivXB, vcoordY,
         num_triangles, num_tri_ears, tri_vertex_index, 
         flag_reverse_orient, vertex_coord, tri_vert);
      break;

    case 1:
      if (tri_vertex_index < NUM_VERT_PER_QUAD)
        { tri_vertex_index = (tri_vertex_index+2)%NUM_VERT_PER_QUAD; }
      triangulate_quad_with_split_v03_v23
        (dimension, v2, v3, v0, v1, ivXA, ivXB, vcoordY,
         num_triangles, num_tri_ears, tri_vertex_index, 
         flag_reverse_orient, vertex_coord, tri_vert);
      break;

    case 2:
      if (tri_vertex_index < NUM_VERT_PER_QUAD)
        { tri_vertex_index = (tri_vertex_index+1)%NUM_VERT_PER_QUAD; }
      triangulate_quad_with_split_v03_v23
        (dimension, v3, v0, v1, v2, ivXA, ivXB, vcoordY,
         num_triangles, num_tri_ears, tri_vertex_index, 
         flag_reverse_orient, vertex_coord, tri_vert);
      break;

    default:
    case 3:
    triangulate_quad_with_split_v03_v23
      (dimension, v0, v1, v2, v3, ivXA, ivXB, vcoordY,
       num_triangles, num_tri_ears, tri_vertex_index, 
       flag_reverse_orient, vertex_coord, tri_vert);
      break;
    }

  }


  /*!
   *  Triangulate quad with two adjacent split edges.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordY[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE, typename VTYPEX, 
            typename ITYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename CTYPEY, typename CTYPE2, typename VTYPET>
  void triangulate_quad_with_two_adjacent_split_edges
  (const DTYPE dimension, const VTYPE quad_vert[4],
   const VTYPEX ivXA, const VTYPEX ivXB, const CTYPEY vcoordY[],
   const ITYPE split_edge_index,  
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri,
   const bool flag_reverse_orient,
   std::vector<CTYPE2> & vertex_coord,
   std::vector<VTYPET> & tri_vert)
  {
    triangulate_quad_with_two_adjacent_split_edges
      (dimension, quad_vert[0], quad_vert[1], quad_vert[2],
       quad_vert[3], ivXA, ivXB, vcoordY, split_edge_index, 
       quad_tri, flag_reverse_orient, vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate quad with three split edges into 7 triangles.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are (v0, v1, v3, v4)
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordY[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3,
            typename VTYPEX, typename VTYPEY, typename VTYPET>
  void triangulate_quad_split3e_tri7_vX01_vX12_vX23
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX ivX01, const VTYPEX ivX12, const VTYPEX ivX23,
   const VTYPEY ivY, const bool flag_reverse_orient,
   std::vector<VTYPET> & tri_vert)
  {
    const int NUM_VERT_PER_SEPTAGON(7);
    const VTYPE0 septagon_vert[NUM_VERT_PER_SEPTAGON] = 
      { v0, ivX01, v1, ivX12, v2, ivX23, v3 };

    triangulate_polygon_with_vertex
      (NUM_VERT_PER_SEPTAGON, septagon_vert, ivY, flag_reverse_orient, 
       tri_vert);
  }


  /*!
   *  Triangulate quad with three split edges into 7 triangles.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Quad vertices are (v0, v1, v3, v4)
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordY[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3,
            typename VTYPEX, typename VTYPEY,
            typename ITYPE, typename VTYPET>
  void triangulate_quad_split3e_tri7
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3, 
   const VTYPEX ivXA, const VTYPEX ivXB, const VTYPEX ivXC,
   const VTYPEY ivY, 
   const ITYPE split_edge_index,
   const bool flag_reverse_orient,
   std::vector<VTYPET> & tri_vert)
  {
    switch (split_edge_index) {

    case 1:
      triangulate_quad_split3e_tri7_vX01_vX12_vX23
        (dimension, v1, v2, v3, v0, ivXA, ivXB, ivXC, ivY, 
         flag_reverse_orient, tri_vert);
      break;

    case 2:
      triangulate_quad_split3e_tri7_vX01_vX12_vX23
        (dimension, v2, v3, v0, v1, ivXA, ivXB, ivXC, ivY, 
         flag_reverse_orient, tri_vert);
      break;

    case 3:
      triangulate_quad_split3e_tri7_vX01_vX12_vX23
        (dimension, v3, v0, v1, v3, ivXA, ivXB, ivXC, ivY, 
         flag_reverse_orient, tri_vert);
      break;

    default:
    case 0:
      triangulate_quad_split3e_tri7_vX01_vX12_vX23
        (dimension, v0, v1, v2, v3, ivXA, ivXB, ivXC, ivY, 
         flag_reverse_orient, tri_vert);
      break;
    }
  }


  /*!
   *  Triangulate quad with three split edges into 7 triangles.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Version using array quad_vert[] of quad vertices.
   *  - Quad vertices are
   *    (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordY[] Coordinates of vertex used in additional triangulation.
   */
  template <typename DTYPE, typename VTYPE,
            typename VTYPEX, typename VTYPEY,
            typename ITYPE, typename VTYPET>
  void triangulate_quad_split3e_tri7
  (const DTYPE dimension, const VTYPE quad_vert[4],
   const VTYPEX ivXA, const VTYPEX ivXB, const VTYPEX ivXC,
   const VTYPEY ivY, 
   const ITYPE split_edge_index,
   const bool flag_reverse_orient,
   std::vector<VTYPET> & tri_vert)
  {
    triangulate_quad_split3e_tri7
      (dimension, quad_vert[0], quad_vert[1], quad_vert[2],
       quad_vert[3], ivXA, ivXB, ivXC, ivY, split_edge_index,
       flag_reverse_orient, tri_vert);
  }

  ///@}

}

#endif
