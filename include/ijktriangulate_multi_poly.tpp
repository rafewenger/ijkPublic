/// \file ijktriangulate_multi_poly.tpp
/// @brief ijk templates for triangulating multiple polytopes.
/// *** ??? DEPRECATED ??? ***
/// - Version 0.4.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2019-2021 Rephael Wenger

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


#ifndef _IJKTRIANGULATE_MULTI_POLY_
#define _IJKTRIANGULATE_MULTI_POLY_

#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijktriangulate_geom.tpp"
#include "ijktriangulate_split.tpp"


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // *********************************************************************
  //! @name COMPUTE COS MAX MIN TRIANGULATION ANGLE OF TWO QUADS
  // *********************************************************************

  ///@{

  /*!
   *  Compute the cosine of the max min triangulation angle of two quads
   *    over the triangulations using the two quad diagonals.
   *  - Quadrilateral A vertices are at coordinates: 
   *      (lower_vcoord0[], lower_vcoord1[], 
   *       upper_vcoord1[], upper_vcoord0[])
   *  - Quadrilateral B vertices are at coordinates: 
   *      (lower_vcoord1[], lower_vcoord2[], 
   *       upper_vcoord2[], upper_vcoord1[])
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_quadA_diag02 If true, triangulation with min angle
   *       has diagonal (lower_vcoord0[], upper_vcoord1[]).  
   *       Otherwise, triangulation with min angle has 
   *       diagonal (lower_vcoord1[], upper_vcoord0[]).
   *  @param[out] flag_quadB_diag02 If true, triangulation with min angle
   *       has diagonal (lower_vcoord1[], upper_vcoord2[]).  
   *       Otherwise, triangulation with min angle has 
   *       diagonal (lower_vcoord2[], upper_vcoord1[]).
   *  @param[out] flag_zero True, if some triangles in triangulation
   *    with min angle have two or three edges less than or equal
   *    to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE_L0, typename CTYPE_L1, typename CTYPE_L2, 
            typename CTYPE_U0, typename CTYPE_U1, typename CTYPE_U2, 
            typename COS_TYPE, typename MTYPE>
  void compute_cos_max_min_two_quads_tri02_or_tri13_angle
  (const DTYPE dimension, 
   const CTYPE_L0 * lower_vcoord0, const CTYPE_L1 * lower_vcoord1, 
   const CTYPE_L2 * lower_vcoord2, const CTYPE_U0 * upper_vcoord0, 
   const CTYPE_U1 * upper_vcoord1, const CTYPE_U2 * upper_vcoord2,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, 
   bool & flag_quadA_diag02, bool & flag_quadB_diag02, bool & flag_zero)
  {
    COS_TYPE cos_min_quadA, cos_min_quadB;
    bool flag_quadA_zero, flag_quadB_zero;

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord0, lower_vcoord1, upper_vcoord1, upper_vcoord0, 
       max_small_magnitude, cos_min_quadA, flag_quadA_diag02, flag_quadA_zero);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord1, lower_vcoord2, upper_vcoord2, upper_vcoord1,
       max_small_magnitude, cos_min_quadB, flag_quadB_diag02, flag_quadB_zero);

    flag_zero = (flag_quadA_zero && flag_quadB_zero);
    cos_min_angle = std::max(cos_min_quadA, cos_min_quadB);
  }


  /*!
   *  Compute the cosine of the max min triangulation angle of two quads
   *    over the triangulations using the two quad diagonals.
   *  - Version with vertex coordinates stored in arrays
   *      lower_vcoord[] and upper_vcoord[].
   *  - Quadrilateral A vertices are at coordinates: 
   *      (lower_vcoord[0][], lower_vcoord[1][], 
   *       upper_vcoord[1][], upper_vcoord[0][])
   *  - Quadrilateral B vertices are at coordinates: 
   *      (lower_vcoord[1][], lower_vcoord[2][], 
   *       upper_vcoord[2][], upper_vcoord[1][])
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_quadA_diag02 If true, triangulation with min angle
   *       has diagonal (lower_vcoord[0][], upper_vcoord[1][]).  
   *       Otherwise, triangulation with min angle has 
   *       diagonal (lower_vcoord[1][], upper_vcoord[0][]).
   *  @param[out] flag_quadB_diag02 If true, triangulation with min angle
   *       has diagonal (lower_vcoord[1][], upper_vcoord[2][]).  
   *       Otherwise, triangulation with min angle has 
   *       diagonal (lower_vcoord[2][], upper_vcoord[1][]).
   *  @param[out] flag_zero True, if some triangles in triangulation
   *    with min angle have two or three edges less than or equal
   *    to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPEL, typename CTYPEU,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_max_min_two_quads_tri02_or_tri13_angle
  (const DTYPE dimension, 
   const CTYPEL * lower_vcoord[3], const CTYPEU * upper_vcoord[3],
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, 
   bool flag_diag02[2], bool & flag_zero)
  {
    compute_cos_max_min_two_quads_tri02_or_tri13_angle
      (dimension, lower_vcoord[0], lower_vcoord[1], lower_vcoord[2],
       upper_vcoord[0], upper_vcoord[1], upper_vcoord[2], max_small_magnitude,
       cos_min_angle, flag_diag02[0], flag_diag02[1], flag_zero);
  }


  /*!
   *  Compute the cosine of the max min triangulation angle of two quads.
   *  - Quadrilateral A vertices are at coordinates: 
   *      (lower_vcoord0[], lower_vcoord1[], 
   *       upper_vcoord1[], upper_vcoord0[])
   *  - Quadrilateral B vertices are at coordinates: 
   *      (lower_vcoord1[], lower_vcoord2[], 
   *       upper_vcoord2[], upper_vcoord1[])
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of the minimum angle in the two quads.
   *  @param[out] flag_zero True if for some quad all triangulations
   *           have zero degree edges.
   */
  template <typename DTYPE, 
            typename CTYPE_L0, typename CTYPE_L1, typename CTYPE_L2, 
            typename CTYPE_U0, typename CTYPE_U1, typename CTYPE_U2, 
            typename CTYPEX, typename MTYPE, 
            typename COS_TYPE0, typename COS_TYPE1, 
            typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_two_quads_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, 
   const CTYPE_L0 * lower_vcoord0, const CTYPE_L1 * lower_vcoord1, 
   const CTYPE_L2 * lower_vcoord2, const CTYPE_U0 * upper_vcoord0, 
   const CTYPE_U1 * upper_vcoord1, const CTYPE_U2 * upper_vcoord2,
   const CTYPEX * quadA_vcoordX, const CTYPEX * quadB_vcoordX,
   const MTYPE max_small_magnitude, 
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE0,NTYPE> 
   quad_tri_result[2],
   COS_TYPE1 & cos_min_angle,
   bool & flag_zero)
  {
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, lower_vcoord0, lower_vcoord1, upper_vcoord1, upper_vcoord0,
       quadA_vcoordX, max_small_magnitude, quad_tri_result[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, lower_vcoord1, lower_vcoord2, upper_vcoord2, upper_vcoord1,
       quadB_vcoordX, max_small_magnitude, quad_tri_result[1]);

    flag_zero = (quad_tri_result[0].flag_zero &&
                 quad_tri_result[1].flag_zero);
    cos_min_angle = std::max(quad_tri_result[0].cos_min_triangulation_angle,
                             quad_tri_result[1].cos_min_triangulation_angle);
  }



  /*!
   *  Compute the cosine of the max min triangulation angle of two quads.
   *  - Version with vertex coordinates stored in arrays
   *      lower_vcoord[] and upper_vcoord[].
   *  - Quadrilateral A vertices are at coordinates: 
   *      (lower_vcoord[0][], lower_vcoord[1][], 
   *       upper_vcoord[1][], upper_vcoord[0][])
   *  - Quadrilateral B vertices are at coordinates: 
   *      (lower_vcoord[1][], lower_vcoord[2][], 
   *       upper_vcoord[2][], upper_vcoord[1][])
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of the minimum angle in the two quads.
   *  @param[out] flag_zero True if for some quad all triangulations
   *           have zero degree edges.
   */
  template <typename DTYPE, 
            typename CTYPEL, typename CTYPEU, typename CTYPEX, typename MTYPE, 
            typename COS_TYPE0, 
            typename COS_TYPE1, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_two_quads_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, 
   const CTYPEL * lower_vcoord[3], const CTYPEU * upper_vcoord[3],
   const CTYPEX * quadA_vcoordX, const CTYPEX * quadB_vcoordX,
   const MTYPE max_small_magnitude, 
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE0,NTYPE> 
   quad_tri_result[2],
   COS_TYPE1 & cos_min_angle,
   bool & flag_zero)
  {
    compute_cos_max_min_two_quads_tri02_or_tri13_or_tri4_angle
      (dimension, lower_vcoord[0], lower_vcoord[1], lower_vcoord[2],
       upper_vcoord[0], upper_vcoord[1], upper_vcoord[2],
       quadA_vcoordX, quadB_vcoordX, max_small_magnitude, 
       quad_tri_result, cos_min_angle, flag_zero);
  }


  /*!
   *  Compute the cosine of the max min triangulation angle of two quads
   *    over the triangulations using a point between lower_v1 and upper_v1.
   *  - Quadrilateral A vertices are at coordinates: 
   *      (lower_vcoord0[], lower_vcoord1[], 
   *       upper_vcoord1[], upper_vcoord0[])
   *  - Quadrilateral B vertices are at coordinates: 
   *      (lower_vcoord1[], lower_vcoord2[], 
   *       upper_vcoord2[], upper_vcoord1[])
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, 
            typename CTYPE_L0, typename CTYPE_L1, typename CTYPE_L2, 
            typename CTYPE_U0, typename CTYPE_U1, typename CTYPE_U2, 
            typename CTYPEX, typename CTYPEA, typename CTYPEB,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_two_quads_X11_vXtri3_or_tri5_angle
  (const DTYPE dimension, 
   const CTYPE_L0 * lower_vcoord0, const CTYPE_L1 * lower_vcoord1, 
   const CTYPE_L2 * lower_vcoord2, const CTYPE_U0 * upper_vcoord0, 
   const CTYPE_U1 * upper_vcoord1, const CTYPE_U2 * upper_vcoord2,
   const CTYPEX * vcoordX11, const CTYPEA * vcoordA, const CTYPEB * vcoordB,
   const MTYPE max_small_magnitude, 
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[2])
  {
    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, upper_vcoord1, upper_vcoord0, 
       lower_vcoord0, lower_vcoord1, vcoordX11, vcoordA, 
       max_small_magnitude, quad_tri[0]);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, lower_vcoord1, lower_vcoord2, 
       upper_vcoord2, upper_vcoord1, vcoordX11, vcoordB, 
       max_small_magnitude, quad_tri[1]);
  }


  /*!
   *  Compute the cosine of the max min triangulation angle of two quads
   *    over the triangulations using a point between v1 and v4.
   *  - Quadrilateral A vertices are at coordinates: 
   *      (lower_vcoord[0][], lower_vcoord[1][], 
   *       upper_vcoord[1][], upper_vcoord[0][])
   *  - Quadrilateral B vertices are at coordinates: 
   *      (lower_vcoord[1][], lower_vcoord[2][], 
   *       upper_vcoord[2][], upper_vcoord[1][])
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of all the angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPEL, typename CTYPEU,
            typename CTYPEX, typename CTYPEA, typename CTYPEB,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_two_quads_X11_vXtri3_or_tri5_angle
  (const DTYPE dimension,  
   const CTYPEL * lower_vcoord[3], const CTYPEU * upper_vcoord[3],
   const CTYPEX * vcoordX11, const CTYPEA * vcoordA, const CTYPEB * vcoordB,
   const MTYPE max_small_magnitude, 
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[2])
  {
    compute_cos_max_min_two_quads_X11_vXtri3_or_tri5_angle
      (dimension, lower_vcoord[0], lower_vcoord[1], lower_vcoord[2],
       upper_vcoord[0], upper_vcoord[1], upper_vcoord[2],
       vcoordX11, vcoordA, vcoordB, max_small_magnitude, quad_tri);
  }

  ///@}


  // ****************************************************************
  //! @name COMPUTE COS MAX MIN TRIANGULATION ANGLE OF THREE QUADS
  // ****************************************************************

  ///@{

  /*! 
   *  Compute the cos of the max min angle of three quadrilaterals
   *    when each is triangulated using a single diagonal.
   *  - Quad A has vertex coordinates:
   *    (lower_vcoord0, lower_vcoord1, upper_vcoord1, upper_vcoord0)
   *  - Quad B has vertex coordinates:
   *    (lower_vcoord1, lower_vcoord2, upper_vcoord2, upper_vcoord1)
   *  - Quad C has vertex coordinates:
   *    (lower_vcoord2, lower_vcoord3, upper_vcoord3, upper_vcoord2)
   */
  template <typename DTYPE, typename CTYPEL, typename CTYPEU,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_three_quads_tri02_or_tri13_angle
  (const DTYPE dimension,
   const CTYPEL lower_vcoord0[], const CTYPEL lower_vcoord1[],
   const CTYPEL lower_vcoord2[], const CTYPEL lower_vcoord3[],
   const CTYPEU upper_vcoord0[], const CTYPEU upper_vcoord1[],
   const CTYPEU upper_vcoord2[], const CTYPEU upper_vcoord3[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2[3])
  {
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord0, lower_vcoord1, upper_vcoord1, upper_vcoord0, 
       max_small_magnitude, quad_tri2[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord1, lower_vcoord2, upper_vcoord2, upper_vcoord1, 
       max_small_magnitude, quad_tri2[1]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord2, lower_vcoord3, upper_vcoord3, upper_vcoord2, 
       max_small_magnitude, quad_tri2[2]);
  }


  /*! 
   *  Compute the cos of the max min angle of three quadrilaterals
   *    when each is triangulated using a single diagonal.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad A has vertex coordinates:
   *    (lower_vcoord[0], lower_vcoord[1], upper_vcoord[1], upper_vcoord[0])
   *  - Quad B has vertex coordinates:
   *    (lower_vcoord[1], lower_vcoord[2], upper_vcoord[2], upper_vcoord[1])
   *  - Quad C has vertex coordinates:
   *    (lower_vcoord[2], lower_vcoord[3], upper_vcoord[3], upper_vcoord[2])
   */
  template <typename DTYPE, typename CTYPEL, typename CTYPEU,
            typename MTYPE,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_three_quads_tri02_or_tri13_angle
  (const DTYPE dimension,
   const CTYPEL * lower_vcoord[4], const CTYPEU * upper_vcoord[4],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2[3])
  {
    compute_cos_max_min_three_quads_tri02_or_tri13_angle
      (dimension, 
       lower_vcoord[0], lower_vcoord[1], lower_vcoord[2], lower_vcoord[3],
       upper_vcoord[0], upper_vcoord[1], upper_vcoord[2], upper_vcoord[3],
       max_small_magnitude, quad_tri2);
  }


  /*! 
   *  Compute the cos of the max min angle of three quadrilaterals
   *    when each is triangulated using a single diagonal.
   *  - Version where input is vertex_coord[] and lists of lower and upper
   *   quadrilateral vertices.
   *  - Quad A has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad C has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPEL, typename VTYPEU, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_three_quads_tri02_or_tri13_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2[3])
  {
    const int NUM_QUADS(3);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension };

    compute_cos_max_min_three_quads_tri02_or_tri13_angle
      (dimension, lower_vcoord, upper_vcoord, max_small_magnitude, quad_tri2);
  }


  /*! 
   *  Compute the cosine of the max min triangulation angle of three quads.
   *  - Quad A has vertex coordinates:
   *    (lower_vcoord0, lower_vcoord1, upper_vcoord1, upper_vcoord0)
   *  - Quad B has vertex coordinates:
   *    (lower_vcoord1, lower_vcoord2, upper_vcoord2, upper_vcoord1)
   *  - Quad C has vertex coordinates:
   *    (lower_vcoord2, lower_vcoord3, upper_vcoord3, upper_vcoord2)
   */
  template <typename DTYPE, typename CTYPEL, typename CTYPEU,
            typename CTYPEX,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_three_quads_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension,
   const CTYPEL lower_vcoord0[], const CTYPEL lower_vcoord1[],
   const CTYPEL lower_vcoord2[], const CTYPEL lower_vcoord3[],
   const CTYPEU upper_vcoord0[], const CTYPEU upper_vcoord1[],
   const CTYPEU upper_vcoord2[], const CTYPEU upper_vcoord3[],
   const CTYPEX * quadA_vcoordX, const CTYPEX * quadB_vcoordX,
   const CTYPEX * quadC_vcoordX,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri_result[3])
  {
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, lower_vcoord0, lower_vcoord1, upper_vcoord1, upper_vcoord0, 
       quadA_vcoordX, max_small_magnitude, quad_tri_result[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, lower_vcoord1, lower_vcoord2, upper_vcoord2, upper_vcoord1, 
       quadB_vcoordX, max_small_magnitude, quad_tri_result[1]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, lower_vcoord2, lower_vcoord3, upper_vcoord3, upper_vcoord2, 
       quadC_vcoordX, max_small_magnitude, quad_tri_result[2]);
  }


  /*! Compute the cosine of the max min triangulation angle of three quads.
   *  - Version using array vcoord[] of quad vertex coordinates.
   *  - Quad vertices have coordinates 
   *    (vcoord[0], vcoord[1], vcoord[2], vcoord[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename DTYPE, typename CTYPEL, typename CTYPEU,
            typename CTYPEX,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_three_quads_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension,
   const CTYPEL * lower_vcoord[4], const CTYPEU * upper_vcoord[4],
   const CTYPEX * quad_vcoordX[3],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   quad_tri_result[3])
  {
    compute_cos_max_min_three_quads_tri02_or_tri13_or_tri4_angle
      (dimension, 
       lower_vcoord[0], lower_vcoord[1], lower_vcoord[2], lower_vcoord[3],
       upper_vcoord[0], upper_vcoord[1], upper_vcoord[2], upper_vcoord[3],
       quad_vcoordX[0], quad_vcoordX[1], quad_vcoordX[2],
       max_small_magnitude, quad_tri_result);
  }


  /*! Compute the cosine of the max min triangulation angle of three quads.
   *  - Version where input is vertex_coord[] and lists of lower and upper
   *   quadrilateral vertices.
   *  - Quad A has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices 
   *     (lower_quad_vert[1], lower_quad_vert[2], 
   *      upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad C has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename VTYPEL, typename VTYPEU,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_three_quads_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const CTYPEX * quadA_vcoordX,
   const CTYPEX * quadB_vcoordX,
   const CTYPEX * quadC_vcoordX,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   quad_tri_result[3])
  {
    const int NUM_QUADS(3);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension };

    compute_cos_max_min_three_quads_tri02_or_tri13_or_tri4_angle
      (dimension, 
       lower_vcoord[0], lower_vcoord[1], lower_vcoord[2], lower_vcoord[3],
       upper_vcoord[0], upper_vcoord[1], upper_vcoord[2], upper_vcoord[3],
       quadA_vcoordX, quadB_vcoordX, quadC_vcoordX,
       max_small_magnitude, quad_tri_result);
  }


  /*! 
   *  Compute the cosine of the max min triangulation angle of three quads.
   *  - Consider triangulation of quad A and quad C into 5 triangles.
   *  - Consider triangulation of quad A into three triangles 
   *    by starring from vcoordX16 or iv0 or iv7.
   *  - Consider triangulation of quad C into three triangles 
   *    by starring from vcoordX25 or iv3 or iv4.
   *  - Version returning data structure POLY_TRIANGULATION_RESULTX2.
   *  - Quad A has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad C has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   *  @param quad_vcoordX[i] Coordinates of vertex inside quad i.
   *  @param block_diagonal[i][0] 
   *    Do not allow triangulation with diagonal (v0,v2) of quad i.
   *  @param block_diagonal[i][1] 
   *    Do not allow triangulation with diagonal (v1,v3) of quad i.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename CTYPEX2,
            typename VTYPEL, typename VTYPEU, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_three_quads_angle_allow_tri5
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX11,
   const CTYPEX * vcoordX22,
   const CTYPEX2 * quad_vcoordX[],
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const MTYPE max_small_magnitude,
   const bool block_diagonal[3][2],
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   quad_tri_result[3],
   COS_TYPE & cos_min_angle,
   bool & flag_zero)
  {
    const int NUM_QUADS(3);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension };
    const CTYPE ** lower_quadB_vcoord = lower_vcoord+1;
    const CTYPE ** upper_quadB_vcoord = upper_vcoord+1;
    const CTYPEX * quadA_vcoord = quad_vcoordX[0];
    const CTYPEX * quadC_vcoord = quad_vcoordX[2];
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quad_tri2_or_tri4[NUM_QUADS];
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadA_triX11, quadC_triX22;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadB_triX11, quadB_triX22; 
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      subquadsB_tri_result; 
    COS_TYPE cos_min_three_quads, cos_min_quadX; 
    COS_TYPE cos_min_triX11, cos_min_triX22;
    COS_TYPE cos_min_triX11_triX22_quadX;
    bool flag_zero_three_quads,flag_zero_quadX;
    bool flag_zero_triX11, flag_zero_triX22;

    compute_cos_max_min_three_quads_tri02_or_tri13_or_tri4_angle
      (dimension, lower_vcoord, upper_vcoord, quad_vcoordX,
      max_small_magnitude, quad_tri2_or_tri4);

    /* OBSOLETE/ERROR? IS THIS NEEDED FOR ANYTHING?
    const COS_TYPE cos_min_quadAC =
      std::max(quad_tri2_or_tri4[0].cos_min_triangulation_angle,
               quad_tri2_or_tri4[2].cos_min_triangulation_angle);
    */

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vcoordX11, lower_vcoord[1], lower_vcoord[2], 
       upper_vcoord[2], upper_vcoord[1], max_small_magnitude, quadB_triX11);
    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vcoordX22, upper_vcoord[2], upper_vcoord[1], 
       lower_vcoord[1], lower_vcoord[2], max_small_magnitude, quadB_triX22);

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12_LU
      (dimension, lower_quadB_vcoord, upper_quadB_vcoord, vcoordX11, vcoordX22,
       max_small_magnitude, subquadsB_tri_result);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, upper_vcoord[1], upper_vcoord[0], 
       lower_vcoord[0], lower_vcoord[1], vcoordX11, quadA_vcoord, 
       max_small_magnitude, block_diagonal[0], quadA_triX11);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[2], lower_vcoord[3], 
       upper_vcoord[3], upper_vcoord[2], vcoordX22, quadC_vcoord, 
       max_small_magnitude, block_diagonal[2], quadC_triX22);

    cos_min_three_quads = 
      std::max(quad_tri2_or_tri4[0].cos_min_triangulation_angle,
               quad_tri2_or_tri4[1].cos_min_triangulation_angle);
    cos_min_three_quads = 
      std::max(cos_min_three_quads,
               quad_tri2_or_tri4[2].cos_min_triangulation_angle);
    flag_zero_three_quads =
      (quad_tri2_or_tri4[0].flag_zero || quad_tri2_or_tri4[1].flag_zero ||
       quad_tri2_or_tri4[2].flag_zero);

    cos_min_triX11 = std::max(quadA_triX11.cos_min_triangulation_angle,
                              quadB_triX11.cos_min_triangulation_angle);
    cos_min_triX11 = std::max(cos_min_triX11, 
                              quad_tri2_or_tri4[2].cos_min_triangulation_angle);
    cos_min_triX22 = std::max(quadC_triX22.cos_min_triangulation_angle,
                              quadB_triX22.cos_min_triangulation_angle);
    cos_min_triX22 = std::max(cos_min_triX22,
                              quad_tri2_or_tri4[0].cos_min_triangulation_angle);

    flag_zero_triX11 = (quadA_triX11.flag_zero || quadB_triX11.flag_zero || 
                        quad_tri2_or_tri4[2].flag_zero);
    flag_zero_triX22 = (quadC_triX22.flag_zero || quadB_triX22.flag_zero ||
                        quad_tri2_or_tri4[0].flag_zero);

    cos_min_quadX = std::max(subquadsB_tri_result.cos_min_triangulation_angle,
                             quadA_triX11.cos_min_triangulation_angle);
    cos_min_quadX = std::max(cos_min_quadX, 
                             quadC_triX22.cos_min_triangulation_angle);
    flag_zero_quadX = (subquadsB_tri_result.flag_zero ||
                       quadA_triX11.flag_zero || quadC_triX22.flag_zero);

    cos_min_triX11_triX22_quadX = std::min(cos_min_triX11, cos_min_triX22);
    cos_min_triX11_triX22_quadX = 
      std::min(cos_min_triX11_triX22_quadX, cos_min_quadX);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  quad_tri2_or_tri4[0].cos_min_triangulation_angle: "
         << quad_tri2_or_tri4[0].cos_min_triangulation_angle << endl;
    cerr << "  quad_tri2_or_tri4[1].cos_min_triangulation_angle: "
         << quad_tri2_or_tri4[1].cos_min_triangulation_angle << endl;
    cerr << "  quad_tri2_or_tri4[2].cos_min_triangulation_angle: "
         << quad_tri2_or_tri4[2].cos_min_triangulation_angle << endl;
    cerr << "  quadA_triX11.cos_min_triangulation_angle: "
         << quadA_triX11.cos_min_triangulation_angle << endl;
    cerr << "  cos_min_three_quads: "
         << cos_min_three_quads << endl;
    cerr << "  quadA_triX11.num_triangles: " 
         << quadA_triX11.num_triangles << endl;
    cerr << "  quadA_triX11.tri_vertex_index: "
         << quadA_triX11.tri_vertex_index << endl;
    cerr << "  quadC_triX22.cos_min_triangulation_angle: "
         << quadC_triX22.cos_min_triangulation_angle << endl;
    cerr << "  quadC_triX22.num_triangles: " 
         << quadC_triX22.num_triangles << endl;
    cerr << "  quadC_triX22.tri_vertex_index: "
         << quadC_triX22.tri_vertex_index << endl;
    cerr << "  subquadsB_tri_result.cos_min_triangulation_angle: "
         << subquadsB_tri_result.cos_min_triangulation_angle << endl;
    cerr << "  cos_min_quadX: " << cos_min_quadX << endl;
    cerr << "  cos_min_triX11_triX22_quadX: "
         << cos_min_triX11_triX22_quadX << endl;
    cerr << "  block_diagonal[0][0]: " << int(block_diagonal[0][0]) << endl;
    cerr << "  block_diagonal[0][1]: " << int(block_diagonal[0][1]) << endl;
    cerr << "  block_diagonal[2][0]: " << int(block_diagonal[2][0]) << endl;
    cerr << "  block_diagonal[2][1]: " << int(block_diagonal[2][1]) << endl;
    */

    if (flag_zero_triX11 || flag_zero_triX22 || flag_zero_quadX) {
      quad_tri_result[0].Copy0Split(quad_tri2_or_tri4[0]);
      quad_tri_result[1].Copy0Split(quad_tri2_or_tri4[1]);
      quad_tri_result[2].Copy0Split(quad_tri2_or_tri4[2]);
      flag_zero = flag_zero_three_quads;
      cos_min_angle = cos_min_three_quads;
    }
    else if (!flag_zero_three_quads &&
             (cos_min_three_quads <= cos_min_triX11_triX22_quadX)) {
      quad_tri_result[0].Copy0Split(quad_tri2_or_tri4[0]);
      quad_tri_result[1].Copy0Split(quad_tri2_or_tri4[1]);
      quad_tri_result[2].Copy0Split(quad_tri2_or_tri4[2]);
      flag_zero = flag_zero_three_quads;
      cos_min_angle = cos_min_three_quads;
    }
    else if (cos_min_quadX <= cos_min_triX11 &&
             cos_min_quadX <= cos_min_triX22) {
      quad_tri_result[0].Copy1Split(quadA_triX11);
      quad_tri_result[1].Copy2Split(subquadsB_tri_result);
      quad_tri_result[2].Copy1Split(quadC_triX22);
      flag_zero = flag_zero_quadX;
      cos_min_angle = cos_min_quadX;
    }
    else if (cos_min_triX11 <= cos_min_triX22) {
      quad_tri_result[0].Copy1Split(quadA_triX11);
      quad_tri_result[1].Copy1Split(quadB_triX11);
      quad_tri_result[2].Copy0Split(quad_tri2_or_tri4[2]);
      flag_zero = flag_zero_triX11;
      cos_min_angle = cos_min_triX11;
    }
    else {
      quad_tri_result[0].Copy0Split(quad_tri2_or_tri4[0]);
      quad_tri_result[1].Copy1Split(quadB_triX22);
      quad_tri_result[2].Copy1Split(quadC_triX22);
      flag_zero = flag_zero_triX22;
      cos_min_angle = cos_min_triX22;
    }

  }


  /*!
   *  Compute the cosine of the max min triangulation angle of three quads.
   *  - Consider triangulation of quad A and quad C into 5 triangles.
   *  - Consider triangulation of quad A into three triangles 
   *    by starring from vcoordX11 or lower_quad_vert[0] or upper_quad_vert[0].
   *  - Consider triangulation of quad C into three triangles 
   *    by starring from vcoordX31 or lower_quad_vert[5] or upper_quad_vert[2].
   *  @param block_diagonal[i][0] 
   *    Do not allow triangulation with diagonal (v0,v2) of quad i.
   *  @param block_diagonal[i][1] 
   *    Do not allow triangulation with diagonal (v1,v3) of quad i.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_three_quadsL_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX11,
   const CTYPEX * vcoordX31,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadB_vcoord,
   const CTYPEX * quadC_vcoord,
   const VTYPE * lower_quad_vert,
   const VTYPE * upper_quad_vert,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[3][2],
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[3],
   COS_TYPE & cos_min_angle,
   bool & flag_zero)
  {
    const int NUM_QUADS(3);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoordL0 = vertex_coord_ptr+lower_quad_vert[0]*dimension;
    const CTYPE * vcoordL1 = vertex_coord_ptr+lower_quad_vert[1]*dimension;
    const CTYPE * vcoordL2 = vertex_coord_ptr+lower_quad_vert[2]*dimension;
    const CTYPE * vcoordL3 = vertex_coord_ptr+lower_quad_vert[3]*dimension;
    const CTYPE * vcoordL4 = vertex_coord_ptr+lower_quad_vert[4]*dimension;
    const CTYPE * vcoordU0 = vertex_coord_ptr+upper_quad_vert[0]*dimension;
    const CTYPE * vcoordU1 = vertex_coord_ptr+upper_quad_vert[1]*dimension;
    const CTYPE * vcoordU2 = vertex_coord_ptr+upper_quad_vert[2]*dimension;
    COS_TYPE cos_min_tri2, cos_min_triAB, cos_min_triBC, cos_min_triABC;
    COS_TYPE cos_min_triAB_or_triBC, cos_min_triAB_or_triBC_or_triABC;
    bool flag_zero_AB, flag_zero_BC, flag_zero_ABC, flag_zero_tri2;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2_or_tri4[NUM_QUADS];
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadA_triX11, quadC_triX31;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadB_triX11, quadB_triX31; 
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quadB_split2_tri;

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoordL0, vcoordL1, vcoordU1, vcoordU0, 
       quadA_vcoord, max_small_magnitude, quad_tri2_or_tri4[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoordL1, vcoordL2, vcoordL3, vcoordU1, 
       quadB_vcoord, max_small_magnitude, quad_tri2_or_tri4[1]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoordL3, vcoordL4, vcoordU2, vcoordU1, 
       quadC_vcoord, max_small_magnitude, quad_tri2_or_tri4[2]);
    
    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, vcoordU1, vcoordU0, vcoordL0, vcoordL1, vcoordX11,
       quadA_vcoord, max_small_magnitude, block_diagonal[0], quadA_triX11);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, vcoordL3, vcoordL4, vcoordU2, vcoordU1, vcoordX31,
       quadC_vcoord, max_small_magnitude, block_diagonal[2], quadC_triX31);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, vcoordL1, vcoordL2, vcoordL3, vcoordU1, vcoordX11,
       quadB_vcoord, max_small_magnitude, block_diagonal[1], quadB_triX11);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, vcoordU1, vcoordL1, vcoordL2, vcoordL3, vcoordX31,
       quadB_vcoord, max_small_magnitude, block_diagonal[1], quadB_triX31);

    // Compute the max min triangulation which splits two edges of quadB.
    compute_cos_max_min_quad_split2_angle
    (dimension, vcoordL1, vcoordL2, vcoordL3, vcoordU1,  
     vcoordX11, vcoordX31, quadB_vcoord, 
     max_small_magnitude, block_diagonal[1], quadB_split2_tri);

    // cos_min_tri2 is the min angle when quads A, B and C are
    //   triangulated separately, each into two triangles.
    cos_min_tri2 = std::max(quad_tri2_or_tri4[0].cos_min_triangulation_angle,
                            quad_tri2_or_tri4[1].cos_min_triangulation_angle);
    cos_min_tri2 = std::max(cos_min_tri2,
                            quad_tri2_or_tri4[2].cos_min_triangulation_angle);
    flag_zero_tri2 = (quad_tri2_or_tri4[0].flag_zero || 
                      quad_tri2_or_tri4[1].flag_zero ||
                      quad_tri2_or_tri4[2].flag_zero);

    cos_min_triAB = std::max(quadA_triX11.cos_min_triangulation_angle,
                             quadB_triX11.cos_min_triangulation_angle);
    cos_min_triAB = std::max(cos_min_triAB,
                             quad_tri2_or_tri4[2].cos_min_triangulation_angle);
    flag_zero_AB = quadA_triX11.flag_zero || quadB_triX11.flag_zero || 
      quad_tri2_or_tri4[2].flag_zero;

    cos_min_triBC = std::max(quadC_triX31.cos_min_triangulation_angle,
                             quadB_triX31.cos_min_triangulation_angle);
    cos_min_triBC = std::max(cos_min_triBC,
                             quad_tri2_or_tri4[0].cos_min_triangulation_angle);
    flag_zero_BC = quadC_triX31.flag_zero || quadB_triX31.flag_zero || 
      quad_tri2_or_tri4[2].flag_zero;

    cos_min_triABC = std::max(quadA_triX11.cos_min_triangulation_angle,
                              quadC_triX31.cos_min_triangulation_angle);
    cos_min_triABC = std::max(cos_min_triABC,
                              quadB_split2_tri.cos_min_triangulation_angle);
    flag_zero_ABC = quadA_triX11.flag_zero || quadC_triX31.flag_zero || 
      quadB_split2_tri.flag_zero;

    cos_min_triAB_or_triBC = std::min(cos_min_triAB, cos_min_triBC);
    cos_min_triAB_or_triBC_or_triABC =
      std::min(cos_min_triAB_or_triBC, cos_min_triABC);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "cos_min_tri2: " << cos_min_tri2 << endl;
    cerr << "cos_min_triABC: " << cos_min_triABC << endl;
    */


    if (flag_zero_AB || flag_zero_BC || flag_zero_ABC ||
        cos_min_tri2 <= cos_min_triAB_or_triBC_or_triABC) {
      quad_tri[0].Copy0Split(quad_tri2_or_tri4[0]);
      quad_tri[1].Copy0Split(quad_tri2_or_tri4[1]);
      quad_tri[2].Copy0Split(quad_tri2_or_tri4[2]);
      flag_zero = flag_zero_tri2;
      cos_min_angle = cos_min_tri2;
    }
    else if (cos_min_triABC < cos_min_triAB_or_triBC) {
      quad_tri[0].Copy1Split(quadA_triX11);
      quad_tri[1].Copy2Split(quadB_split2_tri);
      quad_tri[2].Copy1Split(quadC_triX31);
      flag_zero = flag_zero_ABC;
      cos_min_angle = cos_min_triABC;
    }
    else if (cos_min_triAB <= cos_min_triBC) {
      quad_tri[0].Copy1Split(quadA_triX11);
      quad_tri[1].Copy1Split(quadB_triX11);
      quad_tri[2].Copy0Split(quad_tri2_or_tri4[2]);
      flag_zero = flag_zero_AB;
      cos_min_angle = cos_min_triAB;
    }
    else {
      // cos_min_triBC < cos_min_triAB)
      quad_tri[0].Copy0Split(quad_tri2_or_tri4[0]);
      quad_tri[1].Copy1Split(quadB_triX31);
      quad_tri[2].Copy1Split(quadC_triX31);
      flag_zero = flag_zero_BC;
      cos_min_angle = cos_min_triBC;
    }

  }

  ///@}


  // ****************************************************************
  //! @name COMPUTE COS ANGLE QUAD-TRIANGLE-QUAD
  // ****************************************************************

  ///@{

  /// Compute the cosine of the max min triangulation angle
  ///   of quad-triangle-quad.
  /// @param block_diagonal[i][0] 
  ///   Do not allow triangulation with diagonal (v0,v2) of quad i.
  /// @param block_diagonal[i][1] 
  ///   Do not allow triangulation with diagonal (v1,v3) of quad i.
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename ITYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_triangle_quad_angle_2pent
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * pentagonA_vcoord,
   const CTYPEX * quadB_vcoord,
   const CTYPEX * pentagonC_vcoord,
   const VTYPE * poly_vert,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_max_min_angle,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,ITYPE> poly_tri[3],
   bool & flag_zero)
  {
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord0 = vertex_coord_ptr+poly_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+poly_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord_ptr+poly_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    const CTYPE * vcoord5 = vertex_coord_ptr+poly_vert[5]*dimension;
    const CTYPE * vcoord6 = vertex_coord_ptr+poly_vert[6]*dimension;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,ITYPE> 
      quadA_tri2, quadC_tri2;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,ITYPE> pentagonA_tri;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,ITYPE> pentagonC_tri;
    COS_TYPE cos_min_triB_angle, cos_min_triX_angle;
    COS_TYPE cos_min_tri2, cos_min_pent2;
    bool flag_zero_triB, flag_zero_triX;
    bool flag_zero_tri2, flag_zero_pent2;

    // Compute cos max min angle of quad A triangulation into two triangles.
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord5, vcoord6, max_small_magnitude, 
       quadA_tri2);

    // Compute cos max min angle of quad C triangulation into two triangles.
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord2, vcoord3, vcoord4, vcoord5, max_small_magnitude, 
       quadC_tri2);

    // Compute cos min angle of triangle B.
    compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord2, vcoord5, max_small_magnitude,
       cos_min_triB_angle, flag_zero_triB);

    // Compute cos min angle of triangle (v1, v2, quadB_vcoord).
    compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord2, quadB_vcoord, max_small_magnitude,
       cos_min_triX_angle, flag_zero_triX);

    // Compute cos max min angle of pentagon A triangulation.
    IJK_DEPRECATED::compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vcoord0, vcoord1, quadB_vcoord, vcoord5, vcoord6,
       pentagonA_vcoord, max_small_magnitude, pentagonA_tri);

    // Compute cos max min angle of pentagon C triangulation.
    IJK_DEPRECATED::compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vcoord2, vcoord3, vcoord4, vcoord5, quadB_vcoord,
       pentagonC_vcoord, max_small_magnitude, pentagonC_tri);

    // cos_min_tri2 is the min angle when quads A and C are
    //   triangulated separately, each into two triangles.
    cos_min_tri2 = std::max(quadA_tri2.cos_min_triangulation_angle, 
                            quadC_tri2.cos_min_triangulation_angle);
    cos_min_tri2 = std::max(cos_min_tri2, cos_min_triB_angle);
    flag_zero_tri2 = (quadA_tri2.flag_zero || quadC_tri2.flag_zero);

    // cos_min_pent2 is the min angle when triangle B is replaced
    //   by triangle (v1, v2, quadB_vcoord) and pentagons A and C
    //   are triangulated.
    cos_min_pent2 = std::max(pentagonA_tri.cos_min_triangulation_angle, 
                             pentagonC_tri.cos_min_triangulation_angle);
    cos_min_pent2 = std::max(cos_min_pent2, cos_min_triX_angle);
    flag_zero_pent2 = (pentagonA_tri.flag_zero || pentagonC_tri.flag_zero ||
                       flag_zero_triX);

    // Triangle B is either retained or replaced by another triangle.
    poly_tri[1].num_triangles = 1;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "pentagonA_tri.cos_min_triangulation_angle: " 
         << pentagonA_tri.cos_min_triangulation_angle << endl;
    cerr << "pentagonC_tri.cos_min_triangulation_angle: " 
         << pentagonC_tri.cos_min_triangulation_angle << endl;
    cerr << "pentagonA_tri.tri_vertex_index: "
         << pentagonA_tri.tri_vertex_index << endl;
    cerr << "pentagonC_tri.tri_vertex_index: "
         << pentagonC_tri.tri_vertex_index << endl;
    cerr << "pentagonA_tri.num_triangles: "
         << pentagonA_tri.num_triangles << endl;
    cerr << "pentagonC_tri.num_triangles: "
         << pentagonC_tri.num_triangles << endl;
    cerr << "cos_min_pent2: " << cos_min_pent2 << endl;
    cerr << "cos_min_tri2: " << cos_min_tri2 << endl;
    */


    if (flag_zero_pent2 || (cos_min_tri2 <= cos_min_pent2)) {
      poly_tri[0] = quadA_tri2;
      poly_tri[2] = quadC_tri2;
      cos_max_min_angle = cos_min_tri2;
      flag_zero = flag_zero_tri2;
    }
    else {
      poly_tri[0] = pentagonA_tri;
      poly_tri[2] = pentagonC_tri;
      cos_max_min_angle = cos_min_pent2;
      flag_zero = flag_zero_pent2;
    }
  }

  ///@}


  // ****************************************************************
  //! @name COMPUTE COS ANGLE QUAD-TRIANGLE-PENTAGON
  // ****************************************************************

  ///@{

  /// Compute the cosine of the max min triangulation angle 
  ///   of a quad-triangle-pentagon.
  /// - Consider triangulation of quad into 5 triangles.
  /// - Consider triangulation of pentagon into 6 triangles.
  /// @param vcoordX[] Coordinates of three new points on edges.
  ///      (Usually midpoints of those edges.)
  /// - First point is on edge (lower_quad_vert[1], upper_quad_vert[1]).
  /// - Second point is on edge (lower_quad_vert[1], upper_quad_vert[2]).
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPEL, typename VTYPER,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_triangle_pentagon_angle_allow_tri5_and_tri6
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * poly_vcoord,
   const VTYPEL lower_poly_vert[],
   const VTYPER upper_poly_vert[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_max_min_angle,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   poly_tri_result[3],
   bool & flag_zero)
  {
    const int NUM_LOWER_POLY_VERT(4);
    const int NUM_UPPER_POLY_VERT(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_POLY_VERT] =
      { vertex_coord_ptr+lower_poly_vert[0]*dimension,
        vertex_coord_ptr+lower_poly_vert[1]*dimension,
        vertex_coord_ptr+lower_poly_vert[2]*dimension,
        vertex_coord_ptr+lower_poly_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_POLY_VERT] =
      { vertex_coord_ptr+upper_poly_vert[0]*dimension,
        vertex_coord_ptr+upper_poly_vert[1]*dimension,
        vertex_coord_ptr+upper_poly_vert[2]*dimension,
        vertex_coord_ptr+upper_poly_vert[3]*dimension };
    const CTYPE * vcoordA = poly_vcoord;
    const CTYPE * vcoordB = poly_vcoord+dimension;
    const CTYPE * vcoordC = poly_vcoord+2*dimension;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadA_tri2;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> pentagonC_tri;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      polyA_triX, polyC_triX;
    COS_TYPE cos_min_triangleB_angle, cos_min_triangleX_angle;
    COS_TYPE cos_max_min_separate_tri, cos_max_min_triX;
    bool flag_zero_separate_tri, flag_zero_triX;
    bool flag_zero_triangleB, flag_zero_triangleX;

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord[0], lower_vcoord[1], 
       upper_vcoord[1], upper_vcoord[0],
       max_small_magnitude, quadA_tri2);

    compute_cos_min_triangle_angle
      (dimension, lower_vcoord[1], upper_vcoord[2], upper_vcoord[1],
       max_small_magnitude, cos_min_triangleB_angle, flag_zero_triangleB);

    compute_cos_min_triangle_angle
      (dimension, vcoordB, upper_vcoord[2], upper_vcoord[1],
       max_small_magnitude, cos_min_triangleX_angle, flag_zero_triangleX);

    // Compute cos max min angle of pentagon C triangulation.
    IJK_DEPRECATED::compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, lower_vcoord[1], lower_vcoord[2], lower_vcoord[3], 
       upper_vcoord[3], upper_vcoord[2], vcoordC, max_small_magnitude, 
       pentagonC_tri);

    cos_max_min_separate_tri = 
      std::max(quadA_tri2.cos_min_triangulation_angle,
               cos_min_triangleB_angle);
    cos_max_min_separate_tri = 
      std::max(cos_max_min_separate_tri,
               pentagonC_tri.cos_min_triangulation_angle);
    flag_zero_separate_tri =(quadA_tri2.flag_zero || flag_zero_triangleB ||
                             pentagonC_tri.flag_zero);

    // Initialize
    poly_tri_result[0] = quadA_tri2;
    poly_tri_result[1].num_triangles = 1;
    poly_tri_result[2].Copy(pentagonC_tri);
    flag_zero = flag_zero_separate_tri;
    cos_max_min_angle = cos_max_min_separate_tri;

    // Compute cos max min angle of triangulation of quad A with 
    //   additional vertex at vcoordB.
    IJK_DEPRECATED::compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, upper_vcoord[1], upper_vcoord[0], 
       lower_vcoord[0], lower_vcoord[1], 
       vcoordB, vcoordA, max_small_magnitude, polyA_triX);

    // Compute cos max min angle of triangulation of pentagon C with
    //   additional vertex at vcoordB.
    IJK_DEPRECATED::compute_cos_max_min_hexagon_tri4_or_tri6_angle
      (dimension, lower_vcoord[1], lower_vcoord[2], lower_vcoord[3],
       upper_vcoord[3], upper_vcoord[2], vcoordB, vcoordC,
       max_small_magnitude, polyC_triX);

    cos_max_min_triX = std::max(cos_min_triangleX_angle,
                                polyA_triX.cos_min_triangulation_angle);
    cos_max_min_triX = std::max(cos_max_min_triX,
                                polyC_triX.cos_min_triangulation_angle);

    flag_zero_triX = (flag_zero_triangleX || polyA_triX.flag_zero ||
                      polyC_triX.flag_zero);
    

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "  vcoordB: ", vcoordB, "\n");
    cerr << "  cos_max_min_separate_tri: "
         << cos_max_min_separate_tri << endl;
    cerr << "  cos_max_min_triX: "
         << cos_max_min_triX << endl;
    cerr << "  cos_min_triangleX_angle: "
         << cos_min_triangleX_angle << endl;
    cerr << "  polyA_triX.cos_min_triangulation_angle: "
         << polyA_triX.cos_min_triangulation_angle << endl;
    cerr << "  polyC_triX.cos_min_triangulation_angle: "
         << polyC_triX.cos_min_triangulation_angle << endl;
    cerr << "  polyC_triX.num_triangles: "
         << polyC_triX.num_triangles << endl;
    cerr << "  polyC_triX.triangle_vertex_index: "
         << polyC_triX.tri_vertex_index << endl;
    */

    int index;
    select_min(cos_max_min_separate_tri, flag_zero_separate_tri,
               cos_max_min_triX, flag_zero_triX,
               cos_max_min_angle, index, flag_zero);

    if (index == 1) {
      poly_tri_result[0] = polyA_triX;
      poly_tri_result[2]= polyC_triX;
    }
    // Otherwise, use default values.

  }

  ///@}


  // ****************************************************************
  //! @name COMPUTE COS ANGLE PENTAGON-QUAD-QUAD
  // ****************************************************************

  ///@{

  /// Compute the cosine of the max min triangulation angle 
  ///   of a pentagoin-quad-quad.
  /// - Consider triangulation of quad into 5 triangles.
  /// - Consider triangulation of pentagon into 6 triangles.
  /// @param vcoordX[] Coordinates of three new points on edges.
  ///      (Usually midpoints of those edges.)
  /// - First point is on edge (lower_quad_vert[1], upper_quad_vert[2]).
  /// - Second point is on edge (lower_quad_vert[2], upper_quad_vert[3]).
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPEL, typename VTYPER,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_pentagon_quad_quad_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX12,
   const CTYPEX * vcoordX23,
   const CTYPEX * poly_vcoord,
   const VTYPEL lower_poly_vert[],
   const VTYPER upper_poly_vert[],
   const MTYPE max_small_magnitude,
   const bool block_diagonal[3][2],
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   poly_tri_result[3],
   COS_TYPE & cos_max_min_angle,
   bool & flag_zero)
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_PENTAGON(5);
    const int NUM_LOWER_POLY_VERT(4);
    const int NUM_UPPER_POLY_VERT(5);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_POLY_VERT] =
      { vertex_coord_ptr+lower_poly_vert[0]*dimension,
        vertex_coord_ptr+lower_poly_vert[1]*dimension,
        vertex_coord_ptr+lower_poly_vert[2]*dimension,
        vertex_coord_ptr+lower_poly_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_POLY_VERT] =
      { vertex_coord_ptr+upper_poly_vert[0]*dimension,
        vertex_coord_ptr+upper_poly_vert[1]*dimension,
        vertex_coord_ptr+upper_poly_vert[2]*dimension,
        vertex_coord_ptr+upper_poly_vert[3]*dimension,
        vertex_coord_ptr+upper_poly_vert[4]*dimension };
    const CTYPE ** lower_quadB_vcoord = lower_vcoord+1;
    const CTYPE ** upper_quadB_vcoord = upper_vcoord+2;
    const VTYPEL pentagonA_vert[NUM_VERT_PER_PENTAGON] = 
      { lower_poly_vert[0], lower_poly_vert[1], 
        upper_poly_vert[2], upper_poly_vert[1], upper_poly_vert[0] };
    const VTYPEL quadB_vert[NUM_VERT_PER_QUAD] = 
      { lower_poly_vert[1], lower_poly_vert[2], 
        upper_poly_vert[3], upper_poly_vert[2] };
    const VTYPEL quadC_vert[NUM_VERT_PER_QUAD] = 
      { lower_poly_vert[2], lower_poly_vert[3],
        upper_poly_vert[4], upper_poly_vert[3] };
    const CTYPE * vcoordA = poly_vcoord;
    const CTYPE * vcoordB = poly_vcoord+dimension;
    const CTYPE * vcoordC = poly_vcoord+2*dimension;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      pentagonA_tri3_or_tri5;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadB_tri2_or_tri4;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadC_tri2_or_tri4;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadB_triX12, quadB_triX23;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      subquadsB_tri_result;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      pentagonA_triX12, quadC_triX23;

    // Compute separate triangulations of poly A, B, and C.
    IJK_DEPRECATED::compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vertex_coord, vcoordA, pentagonA_vert,
       max_small_magnitude, pentagonA_tri3_or_tri5);

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vertex_coord, vcoordB, quadB_vert,
       max_small_magnitude, quadB_tri2_or_tri4);

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vertex_coord, vcoordC, quadC_vert,
       max_small_magnitude, quadC_tri2_or_tri4);

    // Compute cos max min angle of triangulation of pentagon A
    //   with additional vertex at vcoordX12.
    IJK_DEPRECATED::compute_cos_max_min_hexagon_tri4_or_tri6_angle
      (dimension, upper_vcoord[2], upper_vcoord[1], upper_vcoord[0],
       lower_vcoord[0], lower_vcoord[1], vcoordX12, vcoordA,
       max_small_magnitude, pentagonA_triX12);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[1], lower_vcoord[2], 
       upper_vcoord[3], upper_vcoord[2], vcoordX12, vcoordB, 
       max_small_magnitude, block_diagonal[1], quadB_triX12);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, upper_vcoord[3], upper_vcoord[2], 
       lower_vcoord[1], lower_vcoord[2], vcoordX23, vcoordB, 
       max_small_magnitude, block_diagonal[1], quadB_triX23);

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12_LU
      (dimension, lower_quadB_vcoord, upper_quadB_vcoord, 
       vcoordX12, vcoordX23,
       max_small_magnitude, subquadsB_tri_result);

    compute_cos_max_min_quad_tri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[2], lower_vcoord[3], 
       upper_vcoord[4], upper_vcoord[3], vcoordX23, vcoordC, 
       max_small_magnitude, block_diagonal[2], quadC_triX23);

    COS_TYPE cos_max_min_separate_tri = 
      std::max(pentagonA_tri3_or_tri5.cos_min_triangulation_angle,
               quadB_tri2_or_tri4.cos_min_triangulation_angle);
    cos_max_min_separate_tri = 
      std::max(cos_max_min_separate_tri,
               quadC_tri2_or_tri4.cos_min_triangulation_angle);
    bool flag_zero_separate_tri =
      (pentagonA_tri3_or_tri5.flag_zero || quadB_tri2_or_tri4.flag_zero ||
       quadC_tri2_or_tri4.flag_zero);

    COS_TYPE cos_max_min_three_poly =
      std::max(pentagonA_triX12.cos_min_triangulation_angle,
               subquadsB_tri_result.cos_min_triangulation_angle);
    cos_max_min_three_poly =
      std::max(cos_max_min_three_poly, quadC_triX23.cos_min_triangulation_angle);
    bool flag_zero_three_poly =
      (pentagonA_triX12.flag_zero || subquadsB_tri_result.flag_zero ||
       quadC_triX23.flag_zero);

    COS_TYPE cos_max_min_triX12 = 
      std::max(pentagonA_triX12.cos_min_triangulation_angle,
               quadB_triX12.cos_min_triangulation_angle);
    cos_max_min_triX12 =
      std::max(cos_max_min_triX12, quadC_tri2_or_tri4.cos_min_triangulation_angle);
    bool flag_zero_triX12 = 
      (pentagonA_triX12.flag_zero || quadB_triX12.flag_zero ||
       quadC_tri2_or_tri4.flag_zero);

    COS_TYPE cos_max_min_triX23 = 
      std::max(pentagonA_tri3_or_tri5.cos_min_triangulation_angle,
               quadB_triX23.cos_min_triangulation_angle);
    cos_max_min_triX23 =
      std::max(cos_max_min_triX23, quadC_triX23.cos_min_triangulation_angle);
    bool flag_zero_triX23 = 
      (pentagonA_tri3_or_tri5.flag_zero || quadB_triX23.flag_zero ||
       quadC_triX23.flag_zero);

    int index;
    select_minIV(cos_max_min_separate_tri, flag_zero_separate_tri,
                 cos_max_min_triX12, flag_zero_triX12,
                 cos_max_min_triX23, flag_zero_triX23,
                 cos_max_min_three_poly, flag_zero_three_poly,
                 cos_max_min_angle, index, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  cos_max_min_separate_tri: "
         << cos_max_min_separate_tri << endl;
    cerr << "  cos_max_min_three_poly: "
         << cos_max_min_three_poly << endl;
    cerr << "  cos_max_min_triX12: "
         << cos_max_min_triX12 << endl;
    cerr << "  cos_max_min_triX23: "
         << cos_max_min_triX23 << endl;
    cerr << "  index: " << index << endl;
    cerr << endl;
    cerr << "  pentagonA_tri3_or_tri5.cos_min_triangulation_angle: "
         << pentagonA_tri3_or_tri5.cos_min_triangulation_angle << endl;
    cerr << "  quadB_tri2_or_tri4.cos_min_triangulation_angle: "
         << quadB_tri2_or_tri4.cos_min_triangulation_angle << endl;
    cerr << "  quadC_tri2_or_tri4.cos_min_triangulation_angle: "
         << quadC_tri2_or_tri4.cos_min_triangulation_angle << endl;
    cerr << "  pentagonA_triX12.cos_min_triangulation_angle: "
         << pentagonA_triX12.cos_min_triangulation_angle << endl;
    cerr << "  subquadsB_tri_result.cos_min_triangulation_angle: "
         << subquadsB_tri_result.cos_min_triangulation_angle << endl;
    cerr << "  quadC_triX23.cos_min_triangulation_angle: "
         << quadC_triX23.cos_min_triangulation_angle << endl;
    cerr << "  block_diagonal[2]: "
         << int(block_diagonal[2][0]) << " "
         << int(block_diagonal[2][1]) << endl;
    */


    switch(index) {

    case 1: 
      poly_tri_result[0].Copy1Split(pentagonA_triX12);
      poly_tri_result[1].Copy1Split(quadB_triX12);
      poly_tri_result[2].Copy0Split(quadC_tri2_or_tri4);
      break;
      break;

    case 2: 
      poly_tri_result[0].Copy0Split(pentagonA_tri3_or_tri5);
      poly_tri_result[1].Copy1Split(quadB_triX23);
      poly_tri_result[2].Copy1Split(quadC_triX23);
      break;

    case 3:
      poly_tri_result[0].Copy1Split(pentagonA_triX12);
      poly_tri_result[1].Copy2Split(subquadsB_tri_result);
      poly_tri_result[2].Copy1Split(quadC_triX23);
      break;

    default:
      // Triangulate all polygons separately.
      poly_tri_result[0].Copy0Split(pentagonA_tri3_or_tri5);
      poly_tri_result[1].Copy0Split(quadB_tri2_or_tri4);
      poly_tri_result[2].Copy0Split(quadC_tri2_or_tri4);
    }

  }

  ///@}


  // ****************************************************************
  /// @name COMPUTE COS ANGLE QUAD-QUAD-TRIANGLE
  // ****************************************************************

  ///@{

  /// Compute the cosine of the max min triangulation angle 
  ///   of quad-quad-triangle.
  /// - Consider triangulation of quad A into 5 triangles.
  /// @param poly_vert[] Seven vertices shared by the two quads and triangle.
  /// - Quad A has vertices 
  ///     (poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6])
  /// - Quad B has vertices 
  ///     (poly_vert[1], poly_vert[2], poly_vert[4], poly_vert[5])
  /// - Triangle C has vertices 
  ///     (poly_vert[2], poly_vert[3], poly_vert[4])
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_quad_quad_triangle_angle_allow_tri5
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX15,
   const CTYPEX * vcoordX24,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadB_vcoord,
   const VTYPE * poly_vert,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> poly_tri[3],
   COS_TYPE & cos_min_angle,
   bool & flag_zero)
  {
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord0 = vertex_coord_ptr+poly_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+poly_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord_ptr+poly_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    const CTYPE * vcoord5 = vertex_coord_ptr+poly_vert[5]*dimension;
    const CTYPE * vcoord6 = vertex_coord_ptr+poly_vert[6]*dimension;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadA_tri_result;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadB_tri_result;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadA_triX15;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadB_triX15, quadB_triX24;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      subquadsB_tri_result;
    COS_TYPE cos_min_triangleC_angle;
    COS_TYPE cos_min_splitC_angle;
    COS_TYPE cos_min_quadX;
    COS_TYPE cos_min_triX15, cos_min_triX24;
    bool flag_zero_triangleC, flag_zero_splitC;

    // Initialize
    cos_min_angle = 1.0;
    flag_zero = false;

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord5, vcoord6, quadA_vcoord,
       max_small_magnitude, quadA_tri_result);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord1, vcoord2, vcoord4, vcoord5, quadB_vcoord,
       max_small_magnitude, quadB_tri_result);

    compute_cos_min_triangle_angle
      (dimension, vcoord2, vcoord3, vcoord4, max_small_magnitude, 
       cos_min_triangleC_angle, flag_zero_triangleC);

    // Min angle is max cos.
    COS_TYPE cos_min_quadA_angle = quadA_tri_result.cos_min_triangulation_angle;
    COS_TYPE cos_min_quadB_angle = quadB_tri_result.cos_min_triangulation_angle;
    COS_TYPE cos_min_AB_angle = 
      std::max(cos_min_quadA_angle, cos_min_quadB_angle);
    COS_TYPE cos_min_BC_angle = 
      std::max(cos_min_quadB_angle, cos_min_triangleC_angle);
    COS_TYPE cos_min_ABC_angle = 
      std::max(cos_min_AB_angle, cos_min_triangleC_angle);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "quad-quad-triangle vert: ";
    IJK::print_list(cerr, poly_vert, 7);
    cerr << endl;
    */

    // Compute cos max min angle triangulation of quadA.
    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vcoord5, vcoord6, vcoord0, vcoord1, 
       vcoordX15, quadA_vcoord, max_small_magnitude, quadA_triX15);

    // Compute cos min angle in triangulations of quadB.
    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vcoordX15, vcoord1, vcoord2, vcoord4, vcoord5,
       max_small_magnitude, quadB_triX15);

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vcoordX24, vcoord4, vcoord5, vcoord1, vcoord2,
       max_small_magnitude, quadB_triX24);

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12
      (dimension, vcoord1, vcoord2, vcoord4, vcoord5, vcoordX15, vcoordX24,
       max_small_magnitude, subquadsB_tri_result);

    // Compute cos min angle in triangulations of triangle C.
    IJK_DEPRECATED::compute_cos_min_split_triangle_angle
      (dimension, vcoord2, vcoord3, vcoord4, vcoordX24, max_small_magnitude,
       cos_min_splitC_angle, flag_zero_splitC);

    cos_min_triX15 = std::max(quadA_triX15.cos_min_triangulation_angle,
                              quadB_triX15.cos_min_triangulation_angle);
    cos_min_triX24 = std::max(cos_min_splitC_angle, 
                              quadB_triX24.cos_min_triangulation_angle);
    cos_min_triX24 = 
      std::max(cos_min_triX24, quadA_tri_result.cos_min_triangulation_angle);

    cos_min_quadX = std::max(subquadsB_tri_result.cos_min_triangulation_angle,
                             quadA_triX15.cos_min_triangulation_angle);
    cos_min_quadX = std::max(cos_min_quadX, cos_min_splitC_angle);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  cos_min_quadA_angle: "
         << cos_min_quadA_angle << endl;
    cerr << "  cos_min_quadB_angle: "
         << cos_min_quadB_angle << endl;
    cerr << "  cos_min_triangleC_angle: "
         << cos_min_triangleC_angle << endl;
    cerr << "  cos_min_splitC_angle: " << cos_min_splitC_angle << endl;
    cerr << "  quadB_triX24.cos_min_triangulation_angle: "
         << quadB_triX24.cos_min_triangulation_angle << endl;
    cerr << "  cos_min_ABC_angle: " << cos_min_ABC_angle << endl;
    cerr << "  cos_min_triX15: " << cos_min_triX15 << endl;
    cerr << "  cos_min_triX24: " << cos_min_triX24 << endl;
    cerr << "  cos_min_quadX: " << cos_min_quadX << endl;
    */


    if (cos_min_ABC_angle <= cos_min_triX15 &&
        cos_min_ABC_angle <= cos_min_triX24 &&
        cos_min_ABC_angle <= cos_min_quadX) {
      // Triangulate quads A and B separately.
      poly_tri[0].Copy0Split(quadA_tri_result);
      poly_tri[1].Copy0Split(quadB_tri_result);
      poly_tri[2].num_triangles = 1;
    }
    else if (cos_min_quadX < cos_min_triX15 &&
             cos_min_quadX < cos_min_triX24) {
      // Triangulate quads and split triangle C.
      poly_tri[0].Copy1Split(quadA_triX15);
      poly_tri[1].Copy2Split(subquadsB_tri_result);
      poly_tri[2].num_triangles = 2;
      poly_tri[2].num_split_edges = 1;
    }
    else if (cos_min_triX15 < cos_min_triX24) {
      // Triangulate quads.  Don't split triangle C.
      poly_tri[0].Copy1Split(quadA_triX15);
      poly_tri[1].Copy1Split(quadB_triX15);
      poly_tri[2].num_triangles = 1;
    }
    else if (cos_min_triX24 < cos_min_BC_angle) {
      // Triangulate quadB and split triangle C.
      poly_tri[0].Copy0Split(quadA_tri_result);
      poly_tri[1].Copy1Split(quadB_triX24);
      poly_tri[2].num_triangles = 2;
      poly_tri[2].num_split_edges = 1;
    }
    else {
      // Triangulate quads A and B separately.
      poly_tri[0].Copy0Split(quadA_tri_result);
      poly_tri[1].Copy0Split(quadB_tri_result);
      poly_tri[2].num_triangles = 1;
    }

  }

  ///@}


  // ****************************************************************
  //! @name COMPUTE COS MAX MIN TRIANGULATION ANGLE OF FOUR QUADS
  // ****************************************************************

  ///@{

  /*!
   *  Compute the cos of the max min angle of four quadrilaterals
   *    when each is triangulated using a single diagonal.
   *  - Quad A has vertex coordinates:
   *      (lower_vcoord0, lower_vcoord1, upper_vcoord1, upper_vcoord0)
   *  - Quad B has vertex coordinates:
   *      (lower_vcoord1, lower_vcoord2, upper_vcoord2, upper_vcoord1)
   *  - Quad C has vertex coordinates:
   *      (lower_vcoord2, lower_vcoord3, upper_vcoord3, upper_vcoord2)
   *  - Quad D has vertex coordinates:
   *      (lower_vcoord3, lower_vcoord4, upper_vcoord4, upper_vcoord3)
   */
  template <typename DTYPE, typename CTYPEL, typename CTYPEU,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quads_tri02_or_tri13_angle
  (const DTYPE dimension,
   const CTYPEL lower_vcoord0[], const CTYPEL lower_vcoord1[],
   const CTYPEL lower_vcoord2[], const CTYPEL lower_vcoord3[],
   const CTYPEL lower_vcoord4[],
   const CTYPEU upper_vcoord0[], const CTYPEU upper_vcoord1[],
   const CTYPEU upper_vcoord2[], const CTYPEU upper_vcoord3[],
   const CTYPEU upper_vcoord4[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2[4])
  {
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord0, lower_vcoord1, upper_vcoord1, upper_vcoord0, 
       max_small_magnitude, quad_tri2[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord1, lower_vcoord2, upper_vcoord2, upper_vcoord1, 
       max_small_magnitude, quad_tri2[1]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord2, lower_vcoord3, upper_vcoord3, upper_vcoord2, 
       max_small_magnitude, quad_tri2[2]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, lower_vcoord3, lower_vcoord4, upper_vcoord4, upper_vcoord3, 
       max_small_magnitude, quad_tri2[3]);
  }


  /// Compute the cos of the max min angle of four quadrilaterals
  ///   when each is triangulated using a single diagonal.
  template <typename DTYPE, typename CTYPEL, typename CTYPEU,
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quads_tri02_or_tri13_angle
  (const DTYPE dimension,
   const CTYPEL * lower_vcoord[5], const CTYPEU * upper_vcoord[5],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2[4])
  {
    compute_cos_max_min_four_quads_tri02_or_tri13_angle
      (dimension, lower_vcoord[0], lower_vcoord[1], lower_vcoord[2],
       lower_vcoord[3], lower_vcoord[4], upper_vcoord[0], upper_vcoord[1],
       upper_vcoord[2], upper_vcoord[3], upper_vcoord[4],
       max_small_magnitude, quad_tri2);
  }
 

  /// Compute the cos of the max min angle of four quadrilaterals
  ///   when each is triangulated using a single diagonal.
  /// - Quad A has vertices 
  ///     (lower_quad_vert[0], lower_quad_vert[1], 
  ///      upper_quad_vert[1], upper_quad_vert[0])
  /// - Quad B has vertices 
  ///     (lower_quad_vert[1], lower_quad_vert[2], 
  ///      upper_quad_vert[2], upper_quad_vert[1])
  /// - Quad C has vertices 
  ///     (lower_quad_vert[2], lower_quad_vert[3], 
  ///      upper_quad_vert[3], upper_quad_vert[2])
  /// - Quad D has vertices 
  ///     (lower_quad_vert[3], lower_quad_vert[4], 
  ///      upper_quad_vert[4], upper_quad_vert[3])
  template <typename DTYPE, typename CTYPE,
            typename VTYPEL, typename VTYPEU, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quads_tri02_or_tri13_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2[4])
  {
    const int NUM_QUADS(4);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension,
        vertex_coord_ptr+lower_quad_vert[4]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension,
        vertex_coord_ptr+upper_quad_vert[4]*dimension };

    compute_cos_max_min_four_quads_tri02_or_tri13_angle
      (dimension, lower_vcoord, upper_vcoord, max_small_magnitude, quad_tri2);
  }


  /*!
   *  Compute the cos of the max min angle of a quad and three adjacent quads
   *    when each is triangulated using a single diagonal.
   *  - Quad A has vertex coordinates: 
   *      (vcoord0[], vcoord3[], vcoord6[], vcoord9[]).
   *  - Quad B has vertex coordinates: 
   *      (vcoord0[], vcoord1[], vcoord2[], vcoord3[]).
   *  - Quad C has vertex coordinates:
   *      (vcoord3[], vcoord4[], vcoord5[], vcoord6[]).
   *  - Quad D has vertex coordinates:
   *      (vcoord6[], vcoord7[], vcoord8[], vcoord9[]).
   */
  template <typename DTYPE, typename CTYPE, typename MTYPE, 
            typename COS_TYPE0, typename COS_TYPE1, 
            typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quadsT_tri02_or_tri13_angle
  (const DTYPE dimension,
   const CTYPE vcoord0[], const CTYPE vcoord1[], const CTYPE vcoord2[],
   const CTYPE vcoord3[], const CTYPE vcoord4[], const CTYPE vcoord5[],
   const CTYPE vcoord6[], const CTYPE vcoord7[], const CTYPE vcoord8[],
   const CTYPE vcoord9[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE0,NTYPE> quad_tri2[4],
   COS_TYPE1 & cos_min_angle,
   bool & flag_zero)
  {
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord3, vcoord6, vcoord9,
       max_small_magnitude, quad_tri2[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
       max_small_magnitude, quad_tri2[1]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord3, vcoord4, vcoord5, vcoord6,
       max_small_magnitude, quad_tri2[2]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord6, vcoord7, vcoord8, vcoord9,
       max_small_magnitude, quad_tri2[3]);


    cos_min_angle = std::max(quad_tri2[0].cos_min_triangulation_angle,
                             quad_tri2[1].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             quad_tri2[2].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             quad_tri2[3].cos_min_triangulation_angle);

    flag_zero = (quad_tri2[0].flag_zero || quad_tri2[1].flag_zero ||
                 quad_tri2[2].flag_zero || quad_tri2[3].flag_zero);
  }


  ///  Compute the cos of the max min angle of a quad and three adjacent quads
  ///    when each is triangulated using a single diagonal.
  /// - Version using array vlist_coord[].
  template <typename DTYPE, typename CTYPE, typename MTYPE, 
            typename COS_TYPE0, typename COS_TYPE1, 
            typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quadsT_tri02_or_tri13_angle
  (const DTYPE dimension, const CTYPE * vlist_coord[10],   
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE0,NTYPE> 
   quad_tri_result[4],
   COS_TYPE1 & cos_min_angle,
   bool & flag_zero)
  {
    compute_cos_max_min_four_quadsT_tri02_or_tri13_angle
      (dimension, vlist_coord[0], vlist_coord[1], vlist_coord[2], 
       vlist_coord[3], vlist_coord[4], vlist_coord[5], vlist_coord[6], 
       vlist_coord[7], vlist_coord[8], vlist_coord[9],
       max_small_magnitude, quad_tri_result, cos_min_angle, flag_zero);
  }


  /*!
   *  Compute the cos of the max min angle of a quad and three adjacent quads
   *    when each is triangulated using a single diagonal.
   *  - Version using arrays vertex_coord[] and quad_vlist[].
   *  - Quad A has vertices 
   *      (quad_vlist[0], quad_vlist[3], quad_vlist[6], quad_vlist[9]).
   *  - Quad B has vertices 
   *      (quad_vlist[0], quad_vlist[1], quad_vlist[2], quad_vlist[3]).
   *  - Quad C has vertices 
   *      (quad_vlist[3], quad_vlist[4], quad_vlist[5], quad_vlist[6]).
   *  - Quad D has vertices 
   *      (quad_vlist[6], quad_vlist[7], quad_vlist[8], quad_vlist[9]).
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE, typename MTYPE, 
            typename COS_TYPE0, typename COS_TYPE1, 
            typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quadsT_tri02_or_tri13_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE quad_vlist[10],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE0,NTYPE> quad_tri2[4],
   COS_TYPE1 & cos_min_angle,
   bool & flag_zero)
  {
    const int NUM_QUADS(4);
    const int VLIST_LENGTH(10);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vlist_coord[VLIST_LENGTH] =
      { vertex_coord_ptr+quad_vlist[0]*dimension,
        vertex_coord_ptr+quad_vlist[1]*dimension,
        vertex_coord_ptr+quad_vlist[2]*dimension,
        vertex_coord_ptr+quad_vlist[3]*dimension,
        vertex_coord_ptr+quad_vlist[4]*dimension,
        vertex_coord_ptr+quad_vlist[5]*dimension,
        vertex_coord_ptr+quad_vlist[6]*dimension,
        vertex_coord_ptr+quad_vlist[7]*dimension,
        vertex_coord_ptr+quad_vlist[8]*dimension,
        vertex_coord_ptr+quad_vlist[9]*dimension };

    compute_cos_max_min_four_quadsT_tri02_or_tri13_angle
      (dimension, vlist_coord, max_small_magnitude, quad_tri2,
       cos_min_angle, flag_zero);
  }


  /*!
   *  Compute the cosine of the max min triangulation angle of a quad
   *    and three adjacent quads.
   *  - Quad A has vertex coordinates: 
   *      (vcoord0[], vcoord3[], vcoord6[], vcoord9[]).
   *  - Quad B has vertex coordinates: 
   *      (vcoord0[], vcoord1[], vcoord2[], vcoord3[]).
   *  - Quad C has vertex coordinates:
   *      (vcoord3[], vcoord4[], vcoord5[], vcoord6[]).
   *  - Quad D has vertex coordinates:
   *      (vcoord6[], vcoord7[], vcoord8[], vcoord9[]).
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename MTYPE, typename COS_TYPE0, typename COS_TYPE1, 
            typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quadsT_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension,
   const CTYPE vcoord0[], const CTYPE vcoord1[], const CTYPE vcoord2[],
   const CTYPE vcoord3[], const CTYPE vcoord4[], const CTYPE vcoord5[],
   const CTYPE vcoord6[], const CTYPE vcoord7[], const CTYPE vcoord8[],
   const CTYPE vcoord9[],
   const CTYPEX * quadA_vcoordX, const CTYPEX * quadB_vcoordX,
   const CTYPEX * quadC_vcoordX, const CTYPEX * quadD_vcoordX,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE0,NTYPE> quad_tri_result[4],
   COS_TYPE1 & cos_min_angle,
   bool & flag_zero)
  {
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord0, vcoord3, vcoord6, vcoord9,
       quadA_vcoordX, max_small_magnitude, quad_tri_result[0]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
       quadB_vcoordX, max_small_magnitude, quad_tri_result[1]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord3, vcoord4, vcoord5, vcoord6,
       quadC_vcoordX, max_small_magnitude, quad_tri_result[2]);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord6, vcoord7, vcoord8, vcoord9,
       quadD_vcoordX, max_small_magnitude, quad_tri_result[3]);

    cos_min_angle = std::max(quad_tri_result[0].cos_min_triangulation_angle,
                             quad_tri_result[1].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             quad_tri_result[2].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             quad_tri_result[3].cos_min_triangulation_angle);

    flag_zero = (quad_tri_result[0].flag_zero || quad_tri_result[1].flag_zero ||
                 quad_tri_result[2].flag_zero || quad_tri_result[3].flag_zero);
  }


  ///  Compute the cos of the max min angle of a quad and three adjacent quads
  ///    when each is triangulated using a single diagonal.
  /// - Version using array vlist_coord[].
  template <typename DTYPE, typename CTYPE, typename CTYPEX, typename MTYPE, 
            typename COS_TYPE0, typename COS_TYPE1, 
            typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quadsT_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, const CTYPE * vlist_coord[10],   
   const CTYPEX * quad_vcoordX[4],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE0,NTYPE> 
   quad_tri_result[4],
   COS_TYPE1 & cos_min_angle,
   bool & flag_zero)
  {
    compute_cos_max_min_four_quadsT_tri02_or_tri13_or_tri4_angle
      (dimension, vlist_coord[0], vlist_coord[1], vlist_coord[2], 
       vlist_coord[3], vlist_coord[4], vlist_coord[5], vlist_coord[6], 
       vlist_coord[7], vlist_coord[8], vlist_coord[9],
       quad_vcoordX[0], quad_vcoordX[1], quad_vcoordX[2], quad_vcoordX[3],
       max_small_magnitude, quad_tri_result, cos_min_angle, flag_zero);
  }


  /*!
   *  Compute the cosine of the max min triangulation angle of a quad
   *    and three adjacent quads.
   *  - Version using arrays vertex_coord[] and quad_vlist[].
   *  - Quad A has vertices 
   *      (quad_vlist[0], quad_vlist[3], quad_vlist[6], quad_vlist[9]).
   *  - Quad B has vertices 
   *      (quad_vlist[0], quad_vlist[1], quad_vlist[2], quad_vlist[3]).
   *  - Quad C has vertices 
   *      (quad_vlist[3], quad_vlist[4], quad_vlist[5], quad_vlist[6]).
   *  - Quad D has vertices 
   *      (quad_vlist[6], quad_vlist[7], quad_vlist[8], quad_vlist[9]).
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX,
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE0, typename COS_TYPE1, 
            typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quadsT_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE quad_vlist[10],
   const CTYPEX * quadA_vcoordX,
   const CTYPEX * quadB_vcoordX,
   const CTYPEX * quadC_vcoordX,
   const CTYPEX * quadD_vcoordX,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE0,NTYPE> quad_tri_result[4],
   COS_TYPE1 & cos_min_angle,
   bool & flag_zero)
  {
    const int VLIST_LENGTH(10);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vlist_coord[VLIST_LENGTH] =
      { vertex_coord_ptr+quad_vlist[0]*dimension,
        vertex_coord_ptr+quad_vlist[1]*dimension,
        vertex_coord_ptr+quad_vlist[2]*dimension,
        vertex_coord_ptr+quad_vlist[3]*dimension,
        vertex_coord_ptr+quad_vlist[4]*dimension,
        vertex_coord_ptr+quad_vlist[5]*dimension,
        vertex_coord_ptr+quad_vlist[6]*dimension,
        vertex_coord_ptr+quad_vlist[7]*dimension,
        vertex_coord_ptr+quad_vlist[8]*dimension,
        vertex_coord_ptr+quad_vlist[9]*dimension };


    compute_cos_max_min_four_quadsT_tri02_or_tri13_or_tri4_angle
      (dimension, vlist_coord[0], vlist_coord[1], vlist_coord[2], 
       vlist_coord[3], vlist_coord[4], vlist_coord[5], vlist_coord[6], 
       vlist_coord[7], vlist_coord[8], vlist_coord[9],
       quadA_vcoordX, quadB_vcoordX, quadC_vcoordX, quadD_vcoordX,
       max_small_magnitude, quad_tri_result, cos_min_angle, flag_zero);
  }


  /*!
   *  Compute the cosine of the max min triangulation angle of four quads.
   *  - Consider triangulation of quads into 5 triangles.
   *  - Quad A has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad C has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   *  - Quad D has vertices 
   *      (lower_quad_vert[3], lower_quad_vert[4], 
   *       upper_quad_vert[4], upper_quad_vert[3])
   *  @param vcoordX[] Coordinates of three new points on edges.
   *       (Usually midpoints of those edges.)
   *  - First point is on edge (lower_quad_vert[1], upper_quad_vert[1]).
   *  - Second point is on edge (lower_quad_vert[2], upper_quad_vert[2]).
   *  - Third point is on edge (lower_quad_vert[3], upper_quad_vert[3]).
  */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPEL, typename VTYPER, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quads_angle_allow_tri5
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX,
   const CTYPEX * quad_vcoord,
   const VTYPEL lower_quad_vert[],
   const VTYPER upper_quad_vert[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[4],
   bool & flag_zero)
  {
    const int NUM_QUADS(4);    
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension,
        vertex_coord_ptr+lower_quad_vert[4]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension,
        vertex_coord_ptr+upper_quad_vert[4]*dimension };
    const CTYPE ** lower_quadB_vcoord = lower_vcoord+1;
    const CTYPE ** upper_quadB_vcoord = upper_vcoord+1;
    const CTYPE ** lower_quadC_vcoord = lower_vcoord+2;
    const CTYPE ** upper_quadC_vcoord = upper_vcoord+2;
    const CTYPE * vcoordX11 = vcoordX;
    const CTYPE * vcoordX22 = vcoordX+dimension;
    const CTYPE * vcoordX33 = vcoordX+2*dimension;
    const CTYPE * quadA_vcoord = quad_vcoord;
    const CTYPE * quadB_vcoord = quad_vcoord+dimension;
    const CTYPE * quadC_vcoord = quad_vcoord+2*dimension;
    const CTYPE * quadD_vcoord = quad_vcoord+3*dimension;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2[NUM_QUADS];
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadA_triX11, quadD_triX33;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadB_triX11, quadB_triX22;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadC_triX33, quadC_triX22;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> subquadsB_tri_result;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> subquadsC_tri_result;
    COS_TYPE cos_min_quadABCD;
    COS_TYPE cos_min_triBC, cos_min_triABC, cos_min_triBCD, cos_min_triABCD;
    COS_TYPE cos_min_triAB_triCD;
    bool flag_zero_quadABCD, flag_zero_triBC;
    bool flag_zero_triABC, flag_zero_triBCD, flag_zero_triABCD;
    bool flag_zero_triAB_triCD;

    compute_cos_max_min_four_quads_tri02_or_tri13_angle
      (dimension, vertex_coord, lower_quad_vert, upper_quad_vert,
       max_small_magnitude, quad_tri2);

    cos_min_quadABCD = std::max(quad_tri2[0].cos_min_triangulation_angle,
                                quad_tri2[1].cos_min_triangulation_angle);
    cos_min_quadABCD = std::max(cos_min_quadABCD,
                                quad_tri2[2].cos_min_triangulation_angle);
    cos_min_quadABCD = std::max(cos_min_quadABCD,
                                quad_tri2[3].cos_min_triangulation_angle);

    flag_zero_quadABCD =(quad_tri2[0].flag_zero || quad_tri2[1].flag_zero ||
                         quad_tri2[2].flag_zero || quad_tri2[3].flag_zero);

    // Initialize
    quad_tri[0].Copy0Split(quad_tri2[0]);
    quad_tri[1].Copy0Split(quad_tri2[1]);
    quad_tri[2].Copy0Split(quad_tri2[2]);
    quad_tri[3].Copy0Split(quad_tri2[3]);
    flag_zero = flag_zero_quadABCD;
    cos_min_angle = cos_min_quadABCD;

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, upper_vcoord[2], upper_vcoord[1], 
       lower_vcoord[1], lower_vcoord[2],
       vcoordX22, quadB_vcoord, max_small_magnitude, quadB_triX22);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[2], lower_vcoord[3], 
       upper_vcoord[3], upper_vcoord[2],
       vcoordX22, quadC_vcoord, max_small_magnitude, quadC_triX22);

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12_LU
      (dimension, lower_quadB_vcoord, upper_quadB_vcoord, vcoordX11, vcoordX22,
       max_small_magnitude, subquadsB_tri_result);

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12_LU
      (dimension, lower_quadC_vcoord, upper_quadC_vcoord, vcoordX22, vcoordX33,
       max_small_magnitude, subquadsC_tri_result);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, upper_vcoord[1], upper_vcoord[0], 
       lower_vcoord[0], lower_vcoord[1], 
       vcoordX11, quadA_vcoord, max_small_magnitude, quadA_triX11);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[3], lower_vcoord[4], 
       upper_vcoord[4], upper_vcoord[3],
       vcoordX33, quadD_vcoord, max_small_magnitude, quadD_triX33);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[1], lower_vcoord[2], 
       upper_vcoord[2], upper_vcoord[1],
       vcoordX11, quadB_vcoord, max_small_magnitude, quadB_triX11);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, upper_vcoord[3], upper_vcoord[2], 
       lower_vcoord[2], lower_vcoord[3], 
       vcoordX33, quadC_vcoord, max_small_magnitude, quadC_triX33);

    // Ignore quads A and D when only triangulating BC.
    cos_min_triBC = std::max(quadB_triX22.cos_min_triangulation_angle,
                             quadC_triX22.cos_min_triangulation_angle);
    flag_zero_triBC = (quadB_triX22.flag_zero || quadC_triX22.flag_zero);

    // Ignore quad D when only triangulating ABC.
    cos_min_triABC = std::max(quadA_triX11.cos_min_triangulation_angle,
                              quadC_triX22.cos_min_triangulation_angle);
    cos_min_triABC = std::max(cos_min_triABC, 
                              subquadsB_tri_result.cos_min_triangulation_angle);
    flag_zero_triABC = (quadA_triX11.flag_zero || quadC_triX22.flag_zero ||
                        subquadsB_tri_result.flag_zero);

    // Ignore quad A when only triangulating BCD.
    cos_min_triBCD = std::max(quadB_triX22.cos_min_triangulation_angle,
                              quadD_triX33.cos_min_triangulation_angle);
    cos_min_triBCD = std::max(cos_min_triBCD,
                              subquadsC_tri_result.cos_min_triangulation_angle);
    flag_zero_triBCD = (quadB_triX22.flag_zero || quadD_triX33.flag_zero ||
                        subquadsC_tri_result.flag_zero);

    // Triangulation of ABCD where edges (lower_v1,upper_v1), 
    //   (lower_v2,upper_v2) and (lower_v3,upper_v3) are split.
    cos_min_triABCD = std::max(quadA_triX11.cos_min_triangulation_angle,
                               quadD_triX33.cos_min_triangulation_angle);
    cos_min_triABCD = std::max(cos_min_triABCD, 
                               subquadsB_tri_result.cos_min_triangulation_angle);
    cos_min_triABCD = std::max(cos_min_triABCD, 
                               subquadsC_tri_result.cos_min_triangulation_angle);
    flag_zero_triABCD = (quadA_triX11.flag_zero || quadD_triX33.flag_zero ||
                         subquadsB_tri_result.flag_zero ||
                         subquadsC_tri_result.flag_zero);

    // Consider separate triangulations of AB and of CD.
    const COS_TYPE cos_min_triAB =
      std::max(quadA_triX11.cos_min_triangulation_angle,
               quadB_triX11.cos_min_triangulation_angle);
    const COS_TYPE cos_min_triCD =
      std::max(quadD_triX33.cos_min_triangulation_angle,
               quadC_triX33.cos_min_triangulation_angle);
    cos_min_triAB_triCD = std::max(cos_min_triAB, cos_min_triCD);
    flag_zero_triAB_triCD =
      (quadA_triX11.flag_zero || quadB_triX11.flag_zero ||
       quadC_triX33.flag_zero || quadD_triX33.flag_zero);

    int index;
    select_minVI
      (cos_min_quadABCD, flag_zero_quadABCD,
       cos_min_triABCD, flag_zero_triABCD,
       cos_min_triABC, flag_zero_triABC,
       cos_min_triBCD, flag_zero_triBCD,
       cos_min_triBC, flag_zero_triBC,
       cos_min_triAB_triCD, flag_zero_triAB_triCD,
       cos_min_angle, index, flag_zero);


    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  quad_tri2[0].cos_min_triangulation_angle: "
         << quad_tri2[0].cos_min_triangulation_angle << endl;
    cerr << "  quad_tri2[1].cos_min_triangulation_angle: "
         << quad_tri2[1].cos_min_triangulation_angle << endl;
    cerr << "  quad_tri2[2].cos_min_triangulation_angle: "
         << quad_tri2[2].cos_min_triangulation_angle << endl;
    cerr << "  quad_tri2[3].cos_min_triangulation_angle: "
         << quad_tri2[3].cos_min_triangulation_angle << endl;
    cerr << "  quadA_triX11.cos_min_triangulation_angle: "
         << quadA_triX11.cos_min_triangulation_angle << endl;
    cerr << "  quadD_triX33.cos_min_triangulation_angle: "
         << quadD_triX33.cos_min_triangulation_angle << endl;
    cerr << "  cos_min_quadABCD: " << cos_min_quadABCD << endl;
    cerr << "  subquadsB_tri_result.cos_min_triangulation_angle: "
         << subquadsB_tri_result.cos_min_triangulation_angle << endl;
    cerr << "  subquadsC_tri_result.cos_min_triangulation_angle: "
         << subquadsC_tri_result.cos_min_triangulation_angle << endl;
    cerr << "  cos_min_triABCD: " << cos_min_triABCD << endl;
    cerr << "  quadB_triX22.cos_min_triangulation_angle: "
         << quadB_triX22.cos_min_triangulation_angle << endl;
    cerr << "  quadC_triX22.cos_min_triangulation_angle: "
         << quadC_triX22.cos_min_triangulation_angle << endl;
    */


    switch (index) {

    case 1:
      // Triangulate ABCD.
      quad_tri[0].Copy1Split(quadA_triX11);
      quad_tri[1].Copy2Split(subquadsB_tri_result);
      quad_tri[2].Copy2Split(subquadsC_tri_result);
      quad_tri[3].Copy1Split(quadD_triX33);
      flag_zero = flag_zero_triABCD;
      cos_min_angle = cos_min_triABCD;
      break;

    case 2:
      // Triangulate ABC.
      quad_tri[0].Copy1Split(quadA_triX11);
      quad_tri[1].Copy2Split(subquadsB_tri_result);
      quad_tri[2].Copy1Split(quadC_triX22);
      flag_zero = flag_zero_triABC;
      cos_min_angle = cos_min_triABC;
      break;

    case 3:
      // Triangulate BCD
      quad_tri[1].Copy1Split(quadB_triX22);
      quad_tri[2].Copy2Split(subquadsC_tri_result);
      quad_tri[3].Copy1Split(quadD_triX33);
      flag_zero = flag_zero_triBCD;
      cos_min_angle = cos_min_triBCD;
      break;

    case 4:
      // Triangulate BC.
      quad_tri[1].Copy1Split(quadB_triX22);
      quad_tri[2].Copy1Split(quadC_triX22);
      flag_zero = flag_zero_triBC;
      cos_min_angle = cos_min_triBC;
      break;

    case 5:
      // Triangulate AB and CD separately.
      quad_tri[0].Copy1Split(quadA_triX11);
      quad_tri[1].Copy1Split(quadB_triX11);
      quad_tri[2].Copy1Split(quadC_triX33);
      quad_tri[3].Copy1Split(quadD_triX33);
      flag_zero = flag_zero_triAB_triCD;
      cos_min_angle = cos_min_triAB_triCD;
      break;

    default:
      // Triangulate all quads separately into 2 triangles.
      // Leave initial values.
      break;
    }
  }


  /*!
   *  Compute the cosine of the max min triangulation angle of quad A
   *    and three adjacent quads (T shape).
   *  - Consider triangulation of quads into 5 triangles.
   *  - Quad A has vertices 
   *      (quad_vlist[0], quad_vlist[3], quad_vlist[6], quad_vlist[9]).
   *  - Quad B has vertices 
   *      (quad_vlist[0], quad_vlist[1], quad_vlist[2], quad_vlist[3]).
   *  - Quad C has vertices 
   *      (quad_vlist[3], quad_vlist[4], quad_vlist[5], quad_vlist[6]).
   *  - Quad D has vertices 
   *      (quad_vlist[6], quad_vlist[7], quad_vlist[8], quad_vlist[9]).
   *  @param vcoordX[] Coordinates of three new points on edges.
   *       (Usually midpoints of those edges.)
   *  - First point is on edge (quad_vlist[0], quad_vlist[3]).
   *  - Second point is on edge (quad_vlist[3], quad_vlist[6]).
   *  - Third point is on edge (quad_vlist[6], quad_vlist[9]).
  */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_quadsT_angle_allow_tri5
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX,
   const CTYPEX * quad_vcoord,
   const VTYPE quad_vlist[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[4],
   COS_TYPE & cos_min_angle,
   bool & flag_zero)
  {
    const int NUM_QUADS(4);    
    const int VLIST_LENGTH(10);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vlist_coord[VLIST_LENGTH] =
      { vertex_coord_ptr+quad_vlist[0]*dimension,
        vertex_coord_ptr+quad_vlist[1]*dimension,
        vertex_coord_ptr+quad_vlist[2]*dimension,
        vertex_coord_ptr+quad_vlist[3]*dimension,
        vertex_coord_ptr+quad_vlist[4]*dimension,
        vertex_coord_ptr+quad_vlist[5]*dimension,
        vertex_coord_ptr+quad_vlist[6]*dimension,
        vertex_coord_ptr+quad_vlist[7]*dimension,
        vertex_coord_ptr+quad_vlist[8]*dimension,
        vertex_coord_ptr+quad_vlist[9]*dimension };
    const CTYPE * vcoordX03 = vcoordX;
    const CTYPE * vcoordX36 = vcoordX+dimension;
    const CTYPE * vcoordX69 = vcoordX+2*dimension;
    const CTYPE * quadA_vcoord = quad_vcoord;
    const CTYPE * quadB_vcoord = quad_vcoord+dimension;
    const CTYPE * quadC_vcoord = quad_vcoord+2*dimension;
    const CTYPE * quadD_vcoord = quad_vcoord+3*dimension;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quad_tri2_or_tri4[NUM_QUADS];
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadB_triX03, quadC_triX36, quadD_triX69;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadA_triX03, quadA_triX36, quadA_triX69;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      quadA_triX_3edges, quadA_triX_03_69, quadA_triX_03_36, quadA_triX_36_69;
    COS_TYPE cos_min_quadABCD, cos_min_triABCD;
    bool flag_zero_quadABCD, flag_zero_triABCD;

    compute_cos_max_min_four_quadsT_tri02_or_tri13_or_tri4_angle
      (dimension, vertex_coord, quad_vlist, quadA_vcoord, quadB_vcoord,
       quadC_vcoord, quadD_vcoord, max_small_magnitude, quad_tri2_or_tri4,
       cos_min_quadABCD, flag_zero_quadABCD);

    // Initialize
    quad_tri[0].Copy0Split(quad_tri2_or_tri4[0]);
    quad_tri[1].Copy0Split(quad_tri2_or_tri4[1]);
    quad_tri[2].Copy0Split(quad_tri2_or_tri4[2]);
    quad_tri[3].Copy0Split(quad_tri2_or_tri4[3]);
    flag_zero = flag_zero_quadABCD;
    cos_min_angle = cos_min_quadABCD;


    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vlist_coord[0], vlist_coord[1], 
       vlist_coord[2], vlist_coord[3], vcoordX03, quadB_vcoord,
       max_small_magnitude, quadB_triX03);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vlist_coord[3], vlist_coord[4], 
       vlist_coord[5], vlist_coord[6], vcoordX36, quadC_vcoord, 
       max_small_magnitude, quadC_triX36);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vlist_coord[6], vlist_coord[7], 
       vlist_coord[8], vlist_coord[9], vcoordX69, quadD_vcoord, 
       max_small_magnitude, quadD_triX69);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vlist_coord[3], vlist_coord[6], 
       vlist_coord[9], vlist_coord[0], vcoordX03, quadA_vcoord, 
       max_small_magnitude, quadA_triX03);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vlist_coord[6], vlist_coord[9], 
       vlist_coord[0], vlist_coord[3], vcoordX36, quadA_vcoord, 
       max_small_magnitude, quadA_triX36);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vlist_coord[9], vlist_coord[0], 
       vlist_coord[3], vlist_coord[6], vcoordX69, quadA_vcoord, 
       max_small_magnitude, quadA_triX69);

    compute_cos_max_min_quad_splitL_tri6_angle_vX01_vX12
      (dimension, vlist_coord[0], vlist_coord[3], 
       vlist_coord[6], vlist_coord[9], vcoordX03, vcoordX36, quadA_vcoord, 
       max_small_magnitude, quadA_triX_03_36);

    compute_cos_max_min_quad_splitL_tri6_angle_vX01_vX12
      (dimension, vlist_coord[3], vlist_coord[6], 
       vlist_coord[9], vlist_coord[0], vcoordX36, vcoordX69, quadA_vcoord, 
       max_small_magnitude, quadA_triX_36_69);

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12
      (dimension, vlist_coord[9], vlist_coord[0], 
       vlist_coord[3], vlist_coord[6], vcoordX69, vcoordX03,
       max_small_magnitude, quadA_triX_03_69);

    compute_cos_min_quad_split3e_tri7_angle_vX01_vX12_vX23
      (dimension, vlist_coord[0], vlist_coord[3], 
       vlist_coord[6], vlist_coord[9],vcoordX03, vcoordX36, vcoordX69, 
       quadA_vcoord, max_small_magnitude, quadA_triX_3edges);


    // Triangulate all quads.
    cos_min_triABCD = quadA_triX_3edges.cos_min_triangulation_angle;
    cos_min_triABCD = std::max(cos_min_triABCD, quadB_triX03.cos_min_triangulation_angle);
    cos_min_triABCD = std::max(cos_min_triABCD, quadC_triX36.cos_min_triangulation_angle);
    cos_min_triABCD = std::max(cos_min_triABCD, quadD_triX69.cos_min_triangulation_angle);

    flag_zero_triABCD = (quadA_triX_3edges.flag_zero || 
                         quadB_triX03.flag_zero || quadC_triX36.flag_zero ||
                         quadD_triX69.flag_zero);

    // Triangulate only quads A and B.
    COS_TYPE cos_min_triAB = std::max(quadA_triX03.cos_min_triangulation_angle,
                                      quadB_triX03.cos_min_triangulation_angle);
    bool flag_zero_triAB = (quadA_triX03.flag_zero || quadB_triX03.flag_zero);

    // Triangulate only quads A and C.
    COS_TYPE cos_min_triAC = std::max(quadA_triX36.cos_min_triangulation_angle,
                                      quadC_triX36.cos_min_triangulation_angle);
    bool flag_zero_triAC = (quadA_triX36.flag_zero || quadC_triX36.flag_zero);

    // Triangulate only quads A and D.
    COS_TYPE cos_min_triAD = std::max(quadA_triX69.cos_min_triangulation_angle,
                                      quadD_triX69.cos_min_triangulation_angle);
    bool flag_zero_triAD = (quadA_triX69.flag_zero || quadD_triX69.flag_zero);

    // Triangulate A, B and C.
    COS_TYPE cos_min_triABC = quadA_triX_03_36.cos_min_triangulation_angle;
    cos_min_triABC = std::max(cos_min_triABC, quadB_triX03.cos_min_triangulation_angle);
    cos_min_triABC = std::max(cos_min_triABC, quadC_triX36.cos_min_triangulation_angle);
    bool flag_zero_triABC = (quadA_triX_03_36.flag_zero || 
                             quadB_triX03.flag_zero || quadC_triX36.flag_zero);

    // Triangulate A, C and D.
    COS_TYPE cos_min_triACD = quadA_triX_36_69.cos_min_triangulation_angle;
    cos_min_triACD = std::max(cos_min_triACD, quadC_triX36.cos_min_triangulation_angle);
    cos_min_triACD = std::max(cos_min_triACD, quadD_triX69.cos_min_triangulation_angle);
    bool flag_zero_triACD = (quadA_triX_36_69.flag_zero || 
                             quadC_triX36.flag_zero || quadD_triX69.flag_zero);

    // Triangulate A, B and D.
    COS_TYPE cos_min_triABD = quadA_triX_03_69.cos_min_triangulation_angle;
    cos_min_triABD = std::max(cos_min_triABD, quadB_triX03.cos_min_triangulation_angle);
    cos_min_triABD = std::max(cos_min_triABD, quadD_triX69.cos_min_triangulation_angle);
    bool flag_zero_triABD = (quadA_triX_03_69.flag_zero || 
                             quadB_triX03.flag_zero || quadD_triX69.flag_zero);


    int index;
    select_minVIII
      (cos_min_quadABCD, flag_zero_quadABCD,
       cos_min_triAB, flag_zero_triAB,
       cos_min_triAC, flag_zero_triAC,
       cos_min_triAD, flag_zero_triAD,
       cos_min_triABC, flag_zero_triABC,
       cos_min_triABD, flag_zero_triABD,
       cos_min_triACD, flag_zero_triACD,
       cos_min_triABCD, flag_zero_triABCD,
       cos_min_angle, index, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "cos_min_tri2: "
         << quad_tri2_or_tri4[0].cos_min_triangulation_angle << " "
         << quad_tri2_or_tri4[1].cos_min_triangulation_angle << " "
         << quad_tri2_or_tri4[2].cos_min_triangulation_angle << " "
         << quad_tri2_or_tri4[3].cos_min_triangulation_angle << " " << endl;
    cerr << "quadA_triX_03_69.cos_min_triangulation_angle: "
         << quadA_triX_03_69.cos_min_triangulation_angle << endl;
    cerr << "quadB_triX03.cos_min_triangulation_angle: "
         << quadB_triX03.cos_min_triangulation_angle << endl;
    cerr << "quadC_triX36.cos_min_triangulation_angle: "
         << quadC_triX36.cos_min_triangulation_angle << endl;
    cerr << "quadD_triX69.cos_min_triangulation_angle: "
         << quadD_triX69.cos_min_triangulation_angle << endl;
    cerr << "cos_min_quadABCD: " << cos_min_quadABCD << endl;
    cerr << "cos_min_triAB: " << cos_min_triAB << endl;
    cerr << "cos_min_triAC: " << cos_min_triAC << endl;
    cerr << "cos_min_triAD: " << cos_min_triAD << endl;
    cerr << "cos_min_triABC: " << cos_min_triABC << endl;
    cerr << "cos_min_triABD: " << cos_min_triABD << endl;
    cerr << "cos_min_triACD: " << cos_min_triACD << endl;
    cerr << "cos_min_triABCD: " << cos_min_triABCD << endl;
    cerr << "  switch index: " << index << endl;
    */

    switch (index) {

    case 1:
      // Triangulate AB.
      quad_tri[0].Copy1Split(quadA_triX03);
      quad_tri[1].Copy1Split(quadB_triX03);
      cos_min_angle = cos_min_triAB;
      flag_zero = flag_zero_triAB;
      break;

    case 2:
      // Triangulate AC.
      quad_tri[0].Copy1Split(quadA_triX36);
      quad_tri[2].Copy1Split(quadC_triX36);
      cos_min_angle = cos_min_triAC;
      flag_zero = flag_zero_triAC;
      break;

    case 3:
      // Triangulate AD.
      quad_tri[0].Copy1Split(quadA_triX69);
      quad_tri[3].Copy1Split(quadD_triX69);
      cos_min_angle = cos_min_triAD;
      flag_zero = flag_zero_triAD;
      break;

    case 4:
      // Triangulate ABC.
      quad_tri[0].Copy2Split(quadA_triX_03_36);
      quad_tri[1].Copy1Split(quadB_triX03);
      quad_tri[2].Copy1Split(quadC_triX36);
      flag_zero = flag_zero_triABC;
      cos_min_angle = cos_min_triABC;
      break;

    case 5:
      // Triangulate ABD.
      quad_tri[0].Copy2Split(quadA_triX_03_69);
      quad_tri[1].Copy1Split(quadB_triX03);
      quad_tri[3].Copy1Split(quadD_triX69);
      flag_zero = flag_zero_triABD;
      cos_min_angle = cos_min_triABD;
      break;

    case 6:
      // Triangulate ACD.
      quad_tri[0].Copy2Split(quadA_triX_36_69);
      quad_tri[2].Copy1Split(quadC_triX36);
      quad_tri[3].Copy1Split(quadD_triX69);
      flag_zero = flag_zero_triACD;
      cos_min_angle = cos_min_triACD;
      break;

    case 7:
      // Triangulate ABCD.
      quad_tri[0].Copy3Split(quadA_triX_3edges);
      quad_tri[1].Copy1Split(quadB_triX03);
      quad_tri[2].Copy1Split(quadC_triX36);
      quad_tri[3].Copy1Split(quadD_triX69);
      flag_zero = flag_zero_triACD;
      cos_min_angle = cos_min_triACD;
      break;

    default:
      // Triangulate all quads separately into 2 triangles.
      // Leave initial values.
      break;
    }

  }


  /*!
   *  Compute the cosine of the max min triangulation angle of
   *    four polygons in T shape, a quad and an adjacent
   *    triangle, quad and triangle.
   *    and three adjacent quads (T shape).
   *  - Consider triangulation of quads into 5 triangles.
   *  - Quad A has vertices 
   *      (quad_vlist[0], quad_vlist[3], quad_vlist[6], quad_vlist[9]).
   *  - Quad B has vertices 
   *      (quad_vlist[0], quad_vlist[1], quad_vlist[2], quad_vlist[3]).
   *  - Quad C has vertices 
   *      (quad_vlist[3], quad_vlist[4], quad_vlist[5], quad_vlist[6]).
   *  - Quad D has vertices 
   *      (quad_vlist[6], quad_vlist[7], quad_vlist[8], quad_vlist[9]).
   *  @param vcoordX[] Coordinates of three new points on edges.
   *       (Usually midpoints of those edges.)
   *  - First point is on edge (quad_vlist[0], quad_vlist[3]).
   *  - Second point is on edge (quad_vlist[3], quad_vlist[6]).
   *  - Third point is on edge (quad_vlist[6], quad_vlist[9]).
  */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_four_polyT_QTQT_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX,
   const CTYPEX quadA_vcoordX[],
   const CTYPEX quadC_vcoordX[],
   const VTYPE poly_vlist[],
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[4],
   COS_TYPE & cos_min_angle,
   bool & flag_zero)
  {
    const int NUM_VERT_PER_TRIANGLE(3);
    const int NUM_VERT_PER_QUAD(4);
    const int VLIST_LENGTH(8);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vlist_coord[VLIST_LENGTH] =
      { vertex_coord_ptr+poly_vlist[0]*dimension,
        vertex_coord_ptr+poly_vlist[1]*dimension,
        vertex_coord_ptr+poly_vlist[2]*dimension,
        vertex_coord_ptr+poly_vlist[3]*dimension,
        vertex_coord_ptr+poly_vlist[4]*dimension,
        vertex_coord_ptr+poly_vlist[5]*dimension,
        vertex_coord_ptr+poly_vlist[6]*dimension,
        vertex_coord_ptr+poly_vlist[7]*dimension };
    const VTYPE quadA_vert[NUM_VERT_PER_QUAD] = 
      { poly_vlist[0], poly_vlist[2], poly_vlist[5], poly_vlist[7] };
    const VTYPE triangleB_vert[NUM_VERT_PER_TRIANGLE] = 
      { poly_vlist[0], poly_vlist[1], poly_vlist[2] };
    const VTYPE quadC_vert[NUM_VERT_PER_QUAD] = 
      { poly_vlist[2], poly_vlist[3], poly_vlist[4], poly_vlist[5] };
    const VTYPE triangleD_vert[NUM_VERT_PER_TRIANGLE] = 
      { poly_vlist[5], poly_vlist[6], poly_vlist[7] };
    const CTYPE * vcoordX02 = vcoordX;
    const CTYPE * vcoordX25 = vcoordX+dimension;
    const CTYPE * vcoordX57 = vcoordX+2*dimension;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadA_tri2_or_tri4;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadC_tri2_or_tri4;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadC_triX25;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quadA_triX_3edges;
    COS_TYPE cos_min_triangleB_angle, cos_min_triangleD_angle;
    COS_TYPE cos_min_splitB_angle, cos_min_splitD_angle;
    COS_TYPE cos_min_no_split_edges;
    COS_TYPE cos_min_split3e;
    bool flagB_zero, flagD_zero, flag_zero_splitB, flag_zero_splitD;
    bool flag_zero_no_split_edges, flag_zero_split3e;

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vertex_coord, quadA_vcoordX, quadA_vert,
       max_small_magnitude, quadA_tri2_or_tri4);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vertex_coord, quadC_vcoordX, quadC_vert,
       max_small_magnitude, quadC_tri2_or_tri4);

    compute_cos_min_triangle_angle
      (dimension, vertex_coord, triangleB_vert, max_small_magnitude,
       cos_min_triangleB_angle, flagB_zero);
    compute_cos_min_triangle_angle
      (dimension, vertex_coord, triangleD_vert, max_small_magnitude,
       cos_min_triangleD_angle, flagD_zero);

    cos_min_no_split_edges = 
      std::max(quadA_tri2_or_tri4.cos_min_triangulation_angle,
               quadC_tri2_or_tri4.cos_min_triangulation_angle);
    cos_min_no_split_edges = 
      std::max(cos_min_no_split_edges, cos_min_triangleB_angle);
    cos_min_no_split_edges = 
      std::max(cos_min_no_split_edges, cos_min_triangleD_angle);
    flag_zero_no_split_edges =
      (quadA_tri2_or_tri4.flag_zero || quadC_tri2_or_tri4.flag_zero || 
       flagB_zero || flagD_zero);

    // Initialize
    quad_tri[0].Copy0Split(quadA_tri2_or_tri4);
    quad_tri[1].num_triangles = 1;
    quad_tri[2].Copy0Split(quadC_tri2_or_tri4);
    quad_tri[3].num_triangles = 1;
    cos_min_angle = cos_min_no_split_edges;
    flag_zero = flag_zero_no_split_edges;

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, vlist_coord[2], vlist_coord[3], 
       vlist_coord[4], vlist_coord[5], vcoordX25, quadC_vcoordX,
       max_small_magnitude, quadC_triX25);

    compute_cos_min_quad_split3e_tri7_angle_vX01_vX12_vX23
      (dimension, vlist_coord[0], vlist_coord[2], 
       vlist_coord[5], vlist_coord[7],vcoordX02, vcoordX25, vcoordX57,
       quadA_vcoordX, max_small_magnitude, quadA_triX_3edges);

    // Compute cos min angle in triangulations of triangles B and C.
    IJK_DEPRECATED::compute_cos_min_split_triangle_angle
      (dimension, vlist_coord[0], vlist_coord[1], vlist_coord[2], 
       vcoordX02, max_small_magnitude, 
       cos_min_splitB_angle, flag_zero_splitB);

    IJK_DEPRECATED::compute_cos_min_split_triangle_angle
      (dimension, vlist_coord[5], vlist_coord[6], vlist_coord[7], 
       vcoordX57, max_small_magnitude, 
       cos_min_splitD_angle, flag_zero_splitD);

    cos_min_split3e =
      std::max(quadA_triX_3edges.cos_min_triangulation_angle,
               quadC_triX25.cos_min_triangulation_angle);
    cos_min_split3e =
      std::max(cos_min_split3e, cos_min_splitB_angle);
    cos_min_split3e =
      std::max(cos_min_split3e, cos_min_splitD_angle);

    flag_zero_split3e =
      (quadA_triX_3edges.flag_zero && quadC_triX25.flag_zero &&
       flag_zero_splitB && flag_zero_splitD);

    int index;
    select_min
      (cos_min_no_split_edges, flag_zero_no_split_edges,
       cos_min_split3e, flag_zero_split3e,
       cos_min_angle, index, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "cos_min_no_split_edges: " << cos_min_no_split_edges << endl;
    cerr << "cos_min_triangleB_angle: " << cos_min_triangleB_angle << endl;
    cerr << "cos_min_triangleD_angle: " << cos_min_triangleD_angle << endl;
    cerr << "quad_AtriX_3edges.cos_min_triangulation_angle: "
         << quadA_triX_3edges.cos_min_triangulation_angle << endl;
    cerr << "quadC_triX25.cos_min_triangulation_angle: "
         << quadC_triX25.cos_min_triangulation_angle << endl;
    cerr << "cos_min_splitB_angle: " << cos_min_splitB_angle << endl;
    cerr << "cos_min_splitD_angle: " << cos_min_splitD_angle << endl;
    cerr << "cos_min_split3e: " << cos_min_split3e << endl;
    */

    switch (index) {

    case 1:
      // Triangulate by splitting all quads and triangles.
      quad_tri[0].Copy3Split(quadA_triX_3edges);
      quad_tri[2].Copy1Split(quadC_triX25);
      quad_tri[1].num_triangles = 2;
      quad_tri[3].num_triangles = 2;
      flag_zero = flag_zero_split3e;
      cos_min_angle = cos_min_split3e;
      break;

    default:
      // Triangulate all quads separately.
      // Leave initial values.
      break;
    }

  }

  ///@}


  // ****************************************************************
  //! @name COMPUTE COS ANGLE TRIANGLE-QUADx3 and QUADx3-TRIANGLE
  // ****************************************************************

  ///@{

  /*!
   *  Compute the cosine of the max min triangulation angle 
   *    of a triangle followed by three quads.
   *  - Consider triangulation of quads into 5 triangles.
   *  - Triangle A has vertices:
   *      (triangle_vert, lower_quad_vert[0], upper_quad_vert[1]).
   *  - Quad B has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad C has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad D has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   *  @param vcoordX[] Coordinates of three new points on edges.
   *       (Usually midpoints of those edges.)
   *  - First point is on edge (lower_quad_vert[0], upper_quad_vert[0]).
   *  - Second point is on edge (lower_quad_vert[1], upper_quad_vert[1]).
   *  - Third point is on edge (lower_quad_vert[2], upper_quad_vert[2]).
  */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename VTYPEL, typename VTYPER, typename VTYPET, 
            typename MTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE>
  void compute_cos_max_min_triangle_quadx3_angle_allow_tri5
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const CTYPEX * vcoordX,
   const CTYPEX * quad_vcoord,
   const VTYPEL lower_quad_vert[],
   const VTYPER upper_quad_vert[],
   const VTYPET triangle_vert,
   const MTYPE max_small_magnitude,
   POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   poly_tri_result[4],
   COS_TYPE & cos_min_angle,
   bool & flag_use_vcoordX11,
   bool & flag_zero)
  {
    const int NUM_QUADS(3);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension };
    const CTYPE * triangle_vcoord =
      vertex_coord_ptr+triangle_vert*dimension;
    const CTYPE ** lower_quadB_vcoord = lower_vcoord;
    const CTYPE ** upper_quadB_vcoord = upper_vcoord;
    const CTYPE ** lower_quadC_vcoord = lower_vcoord+1;
    const CTYPE ** upper_quadC_vcoord = upper_vcoord+1;
    const CTYPE * vcoordX00 = vcoordX;
    const CTYPE * vcoordX11 = vcoordX+dimension;
    const CTYPE * vcoordX22 = vcoordX+2*dimension;
    const CTYPE * quadB_vcoord = quad_vcoord;
    const CTYPE * quadC_vcoord = quad_vcoord+dimension;
    const CTYPE * quadD_vcoord = quad_vcoord+2*dimension;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri2[NUM_QUADS];
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadB_triX00, quadB_triX11;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadC_triX22, quadC_triX11;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> subquadsB_tri_result;
    POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> subquadsC_tri_result;
    POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quadD_triX22;
    COS_TYPE cos_max_min_tri2;
    COS_TYPE cos_min_triangleA_angle, cos_min_splitA_angle;
    COS_TYPE cos_min_triBC, cos_min_triABC, cos_min_triBCD, cos_min_triABCD;
    COS_TYPE cos_min_triAB_triCD;
    bool flag_zero_tri2;
    bool flag_zero_triangleA, flag_zero_splitA;
    bool flag_zero_triBC;
    bool flag_zero_triABC, flag_zero_triBCD, flag_zero_triABCD;
    bool flag_zero_triAB_triCD;

    compute_cos_max_min_three_quads_tri02_or_tri13_angle
      (dimension, vertex_coord, lower_quad_vert, upper_quad_vert,
       max_small_magnitude, quad_tri2);

    compute_cos_min_triangle_angle
      (dimension, lower_vcoord[0], upper_vcoord[0], triangle_vcoord,
       max_small_magnitude, cos_min_triangleA_angle, flag_zero_triangleA);

    cos_max_min_tri2 = std::max(quad_tri2[0].cos_min_triangulation_angle,
                                quad_tri2[1].cos_min_triangulation_angle);
    cos_max_min_tri2 = std::max(cos_max_min_tri2,
                                quad_tri2[2].cos_min_triangulation_angle);
    cos_max_min_tri2 = std::max(cos_max_min_tri2,
                                cos_min_triangleA_angle);
    flag_zero_tri2 =(quad_tri2[0].flag_zero || quad_tri2[1].flag_zero ||
                     quad_tri2[2].flag_zero || flag_zero_triangleA);

    // Initialize
    poly_tri_result[0].num_triangles = 1;
    poly_tri_result[1].Copy(quad_tri2[0]);
    poly_tri_result[2].Copy(quad_tri2[1]);
    poly_tri_result[3].Copy(quad_tri2[2]);
    flag_zero = flag_zero_tri2;
    cos_min_angle = cos_max_min_tri2;
    flag_use_vcoordX11 = false;


    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, upper_vcoord[1], upper_vcoord[0], 
       lower_vcoord[0], lower_vcoord[1],
       vcoordX11, quadB_vcoord, max_small_magnitude, quadB_triX11);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[1], lower_vcoord[2], 
       upper_vcoord[2], upper_vcoord[1],
       vcoordX11, quadC_vcoord, max_small_magnitude, quadC_triX11);

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12_LU
      (dimension, lower_quadB_vcoord, upper_quadB_vcoord, vcoordX00, vcoordX11,
       max_small_magnitude, subquadsB_tri_result);

    compute_cos_max_min_quad_split_into_two_subquads_vX03_vX12_LU
      (dimension, lower_quadC_vcoord, upper_quadC_vcoord, vcoordX11, vcoordX22,
       max_small_magnitude, subquadsC_tri_result);

    
    // Compute cos min angle in triangulations of triangle A.
    IJK_DEPRECATED::compute_cos_min_split_triangle_angle
      (dimension, lower_vcoord[0], triangle_vcoord, upper_vcoord[0], 
       vcoordX00, max_small_magnitude, 
       cos_min_splitA_angle, flag_zero_splitA);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[2], lower_vcoord[3], 
       upper_vcoord[3], upper_vcoord[2],
       vcoordX22, quadD_vcoord, max_small_magnitude, quadD_triX22);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, lower_vcoord[0], lower_vcoord[1], 
       upper_vcoord[1], upper_vcoord[0],
       vcoordX00, quadB_vcoord, max_small_magnitude, quadB_triX00);

    compute_cos_max_min_quad_vXtri3_or_tri5_angle_vX03
      (dimension, upper_vcoord[2], upper_vcoord[1], 
       lower_vcoord[1], lower_vcoord[2], 
       vcoordX22, quadC_vcoord, max_small_magnitude, quadC_triX22);

    // Ignore quads A and D when only triangulating BC.
    cos_min_triBC = std::max(quadB_triX11.cos_min_triangulation_angle,
                             quadC_triX11.cos_min_triangulation_angle);
    flag_zero_triBC = (quadB_triX11.flag_zero || quadC_triX11.flag_zero);

    // Ignore quad D when only triangulating ABC.
    cos_min_triABC = std::max(cos_min_splitA_angle,
                              quadC_triX22.cos_min_triangulation_angle);
    cos_min_triABC = std::max(cos_min_triABC, 
                              subquadsB_tri_result.cos_min_triangulation_angle);
    flag_zero_triABC = (flag_zero_splitA || quadC_triX22.flag_zero ||
                        subquadsB_tri_result.flag_zero);

    // Ignore quad A when only triangulating BCD.
    cos_min_triBCD = std::max(quadB_triX11.cos_min_triangulation_angle,
                              quadD_triX22.cos_min_triangulation_angle);
    cos_min_triBCD = std::max(cos_min_triBCD, 
                              subquadsC_tri_result.cos_min_triangulation_angle);
    flag_zero_triBCD = (quadB_triX11.flag_zero || quadD_triX22.flag_zero ||
                        subquadsC_tri_result.flag_zero);

    // Triangulations of ABCD where edges (lower_v0,upper_v0)
    //   (lower_v1,upper_v1) and (lower_v2,upper_v2) are split.
    cos_min_triABCD = std::max(cos_min_splitA_angle,
                               quadD_triX22.cos_min_triangulation_angle);
    cos_min_triABCD = std::max(cos_min_triABCD, 
                               subquadsB_tri_result.cos_min_triangulation_angle);
    cos_min_triABCD = std::max(cos_min_triABCD, 
                               subquadsC_tri_result.cos_min_triangulation_angle);
    flag_zero_triABCD = (flag_zero_splitA || quadD_triX22.flag_zero ||
                         subquadsB_tri_result.flag_zero ||
                         subquadsC_tri_result.flag_zero);

    // Consider separate triangulations of AB and of CD.
    const COS_TYPE cos_min_triAB =
      std::max(cos_min_splitA_angle,
               quadB_triX11.cos_min_triangulation_angle);
    const COS_TYPE cos_min_triCD =
      std::max(quadD_triX22.cos_min_triangulation_angle,
               quadC_triX22.cos_min_triangulation_angle);
    cos_min_triAB_triCD = std::max(cos_min_triAB, cos_min_triCD);
    flag_zero_triAB_triCD =
      (flag_zero_splitA || quadB_triX00.flag_zero ||
       quadC_triX22.flag_zero || quadD_triX22.flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  cos_min_triangleA_angle: "
         << cos_min_triangleA_angle << endl;
    cerr << "  quad_tri2[0].cos_min_triangulation_angle: "
         << quad_tri2[0].cos_min_triangulation_angle << endl;
    cerr << "  quad_tri2[1].cos_min_triangulation_angle: "
         << quad_tri2[1].cos_min_triangulation_angle << endl;
    cerr << "  quad_tri2[2].cos_min_triangulation_angle: "
         << quad_tri2[2].cos_min_triangulation_angle << endl;
    cerr << "  quadB_triX11.cos_min_triangulation_angle: "
         << quadB_triX11.cos_min_triangulation_angle << endl;
    cerr << "  quadB_triX11.tri_vertex_index: "
         << quadB_triX11.tri_vertex_index << endl;
    cerr << "  quadC_triX11.cos_min_triangulation_angle: "
         << quadC_triX11.cos_min_triangulation_angle << endl;
    cerr << "  quadC_triX11.tri_vertex_index: "
         << quadC_triX11.tri_vertex_index << endl;
    cerr << "  quadC_triX11.num_internal_tri_vertices: "
         << int(quadC_triX11.num_internal_tri_vertices) << endl;
    cerr << "  quadD_triX33.cos_min_triangulation_angle: "
         << quadD_triX22.cos_min_triangulation_angle << endl;
    cerr << "  cos_min_splitA_angle: "
         << cos_min_splitA_angle << endl;
    cerr << "  cos_max_min_tri2: " << cos_max_min_tri2 << endl;
    cerr << "  subquadsB_tri_result.cos_min_triangulation_angle: "
         << subquadsB_tri_result.cos_min_triangulation_angle << endl;
    cerr << "  subquadsC_tri_result.cos_min_triangulation_angle: "
         << subquadsC_tri_result.cos_min_triangulation_angle << endl;
    cerr << "  cos_min_triABCD: " << cos_min_triABCD << endl;
    */

    int index;
    select_minVI
      (cos_max_min_tri2, flag_zero_tri2,
       cos_min_triABCD, flag_zero_triABCD,
       cos_min_triABC, flag_zero_triABC,
       cos_min_triBCD, flag_zero_triBCD,
       cos_min_triBC, flag_zero_triBC,
       cos_min_triAB_triCD, flag_zero_triAB_triCD,
       cos_min_angle, index, flag_zero);

    switch (index) {

    case 1:
      // Triangulate ABCD.
      poly_tri_result[0].num_triangles = 2;
      poly_tri_result[1] = subquadsB_tri_result;
      poly_tri_result[2] = subquadsC_tri_result;
      poly_tri_result[3].Copy(quadD_triX22);
      flag_zero = flag_zero_triABCD;
      cos_min_angle = cos_min_triABCD;
      flag_use_vcoordX11 = true;
      break;

    case 2:
      // Triangulate ABC.
      poly_tri_result[0].num_triangles = 2;
      poly_tri_result[1] = subquadsB_tri_result;
      poly_tri_result[2].Copy(quadC_triX11);
      flag_zero = flag_zero_triABC;
      cos_min_angle = cos_min_triABC;
      flag_use_vcoordX11 = true;
      break;

    case 3:
      // Triangulate BCD
      poly_tri_result[1].Copy(quadB_triX11);
      poly_tri_result[2] = subquadsC_tri_result;
      poly_tri_result[3].Copy(quadD_triX22);
      flag_zero = flag_zero_triBCD;
      cos_min_angle = cos_min_triBCD;
      flag_use_vcoordX11 = true;
      break;

    case 4:
      // Triangulate BC.
      poly_tri_result[1].Copy(quadB_triX11);
      poly_tri_result[2].Copy(quadC_triX11);
      flag_zero = flag_zero_triBC;
      cos_min_angle = cos_min_triBC;
      flag_use_vcoordX11 = true;
      break;

    case 5:
      // Triangulate AB and CD separately.
      poly_tri_result[0].num_triangles = 2;
      poly_tri_result[1].Copy(quadB_triX00);
      poly_tri_result[2].Copy(quadC_triX22);
      poly_tri_result[3].Copy(quadD_triX22);
      flag_zero = flag_zero_triAB_triCD;
      cos_min_angle = cos_min_triAB_triCD;
      flag_use_vcoordX11 = false;
      break;

    default:
      // Triangulate all quads separately into 2 triangles.
      // Leave initial values.
      break;
    }

  }

  ///@}


  // ****************************************************************
  /// @name TRIANGULATE TWO QUADRILATERALS
  // ****************************************************************

  ///@{

  /*!
   *  Triangulate two quadrilaterals sharing an edge. (OLD VERSION.)
   *  - Add new triangles to vector tri_vert.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Allow triangulation of quads into 5 triangles.
   *  @param[in] quad_vert[] Six vertices shared by the two quadrilaterals.
   *  - Six vertices are listed in clockwise or counter-clockwise order
   *      around the two quadrilaterals.
   *  - Quadrilaterals are: 
   *      (quad_vert[0], quad_vert[1], quad_vert[4], quad_vert[5])
   *      (quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[4])
   */
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename TRI_VTYPE, typename CTYPE, typename CTYPEX>
  void triangulate_two_quads_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const bool flag_quad_tri5[2],
   const CTYPEX * vcoordX11,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadB_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const TRI_VTYPE ivertX11 = 
      IJK::insert_coord(dimension, vcoordX11, vertex_coord);

    // Triangulate quadA with additional vertex at ivertX11.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX11, upper_quad_vert[1], upper_quad_vert[0], 
       lower_quad_vert[0], lower_quad_vert[1], quadA_vcoord, flag_quad_tri5[0],
       flag_reverse_orient, vertex_coord, tri_vert);

    // Triangulate quadB with additional vertex at ivertX11.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX11, lower_quad_vert[1], lower_quad_vert[2], 
       upper_quad_vert[2], upper_quad_vert[1], quadB_vcoord, flag_quad_tri5[1],
       flag_reverse_orient, vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate two quadrilaterals sharing an edge. (NEW VERSION.)
   *  - Version using POLY_TRIANGULATION_RESULT.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Allow triangulation of quads into 5 triangles.
   *  @param[in] quad_vert[] Six vertices shared by the two quadrilaterals.
   *  - Six vertices are listed in clockwise or counter-clockwise order
   *      around the two quadrilaterals.
   *  - Quadrilaterals are: 
   *      (quad_vert[0], quad_vert[1], quad_vert[4], quad_vert[5])
   *      (quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[4])
  */
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename TRI_VTYPE, typename CTYPE, typename CTYPEX>
  void triangulate_two_quads_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quadA_tri_result,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quadB_tri_result,
   const CTYPEX * vcoordX11,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadB_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const TRI_VTYPE ivertX11 = 
      IJK::insert_coord(dimension, vcoordX11, vertex_coord);

    // Triangulate quadA with additional vertex at ivertX11.
    IJK_DEPRECATED::triangulate_pentagon_tri3_or_tri5
      (dimension, upper_quad_vert[1], upper_quad_vert[0], 
       lower_quad_vert[0], lower_quad_vert[1], ivertX11,
       quadA_vcoord, flag_reverse_orient, quadA_tri_result,
       vertex_coord, tri_vert);

    // Triangulate quadB with additional vertex at ivertX11.
    IJK_DEPRECATED::triangulate_pentagon_tri3_or_tri5
      (dimension, lower_quad_vert[1], lower_quad_vert[2], 
       upper_quad_vert[2], upper_quad_vert[1], ivertX11,
       quadB_vcoord, flag_reverse_orient, quadB_tri_result,
       vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate two quadrilaterals sharing an edge. (NEW VERSION.)
   *  - Version using array quad_tri_result[2].
   *  - Add new triangles to vector tri_vert.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Allow triangulation of quads into 5 triangles.
   *  @param[in] quad_vert[] Six vertices shared by the two quadrilaterals.
   *  - Six vertices are listed in clockwise or counter-clockwise order
   *      around the two quadrilaterals.
   *  - Quadrilaterals are: 
   *      (quad_vert[0], quad_vert[1], quad_vert[4], quad_vert[5])
   *      (quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[4])
  */
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename TRI_VTYPE, typename CTYPE, typename CTYPEX>
  void triangulate_two_quads_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri_result[2],
   const CTYPEX * vcoordX11,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadB_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    triangulate_two_quads_allow_tri5
      (dimension, lower_quad_vert, upper_quad_vert, flag_reverse_orient,
       quad_tri_result[0], quad_tri_result[1], vcoordX11,
       quadA_vcoord, quadB_vcoord, vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate two quadrilaterals sharing an edge.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coord, if necessary, to vertex_coord[].
   *  - Consider triangulation of quads into 5 triangles.
   *  - Only triangulate quads where merging with quads improves triangulation.
   *  - Quadrilaterals are: 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], quad_vert[1])
   *  @param[out] flag_split True if two quads are split into triangles
   *      which are added to tri_vert[].
  */
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename TRI_VTYPE, typename CTYPE, typename MTYPE>
  void triangulate_two_quads_max_min_angle_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert,
   bool & flag_split)
  {
    const int NUM_QUADS(2);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension };
    IJK::ARRAY<CTYPE> vcoordX(dimension);
    IJK::ARRAY<CTYPE> quadA_centroid(dimension), quadB_centroid(dimension);
    const TRI_VTYPE quadA_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[0], lower_quad_vert[1], 
        upper_quad_vert[1], upper_quad_vert[0] };
    const TRI_VTYPE quadB_vert[NUM_VERT_PER_QUAD] =
      { lower_quad_vert[1], lower_quad_vert[2],
        upper_quad_vert[2], upper_quad_vert[1] };
    CTYPE cos_min_two_quads, cos_min_vertX;
    POLY_TRIANGULATION_RESULT<16,CTYPE,int> quad_tri_result[2];
    POLY_TRIANGULATION_RESULT<16,CTYPE,int> quad_tri_X11_result[2];
    bool flag_zero_two_quads, flag_zero_vertX;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  Lower quad vert: ";
    IJK::print_list(cerr, lower_quad_vert, NUM_LOWER_QUAD_VERT);
    cerr << "\n";
    cerr << "  Upper quad vert: ";
    IJK::print_list(cerr, upper_quad_vert, NUM_UPPER_QUAD_VERT);
    cerr << "\n";
    */

    // Initialize to no split.
    flag_split = false;

    IJK::compute_midpoint
      (dimension, lower_vcoord[1], upper_vcoord[1], vcoordX.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadA_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadA_centroid.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadB_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadB_centroid.Ptr());

    compute_cos_max_min_two_quads_tri02_or_tri13_or_tri4_angle
      (dimension, lower_vcoord, upper_vcoord, 
       quadA_centroid.PtrConst(), quadB_centroid.PtrConst(),  
       max_small_magnitude,
       quad_tri_result, cos_min_two_quads, flag_zero_two_quads);

    compute_cos_max_min_two_quads_X11_vXtri3_or_tri5_angle
      (dimension, lower_vcoord, upper_vcoord, vcoordX.PtrConst(), 
       quadA_centroid.PtrConst(), quadB_centroid.PtrConst(),  
       max_small_magnitude, quad_tri_X11_result);

    cos_min_vertX = std::max(quad_tri_X11_result[0].cos_min_triangulation_angle,
                             quad_tri_X11_result[1].cos_min_triangulation_angle);
    flag_zero_vertX = 
      (quad_tri_X11_result[0].flag_zero || quad_tri_X11_result[1].flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "cos_min_two_quads: " << cos_min_two_quads << endl;
    cerr << "cos_min_vertX: " << cos_min_vertX << endl;
    */

    if (!flag_zero_vertX) {

      if (flag_zero_two_quads || (cos_min_vertX < cos_min_two_quads)) {
        triangulate_two_quads_allow_tri5
          (dimension, lower_quad_vert, upper_quad_vert, flag_reverse_orient, 
           quad_tri_X11_result, vcoordX.PtrConst(), 
           quadA_centroid.PtrConst(), quadB_centroid.PtrConst(),
           vertex_coord, tri_vert);

        flag_split = true;
      }
    }

  }


  /*! 
   *  Triangulate subquads of a quad where the subquads are:
   *    (lower_quad_vert[0], lower_quad_vert[1], ivertX11, ivertX00) and
   *    (ivertX00, ivertX11, upper_quad_vert[1], upper_quad_vert[0]).
   *  - Quad vertices are:
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Procedure suffix LU indicates that quad vertices are listed
   *      as array of lower vertices and array of upper quad vertices.
   */
  template <typename VTYPEL, typename VTYPEU, typename VTYPEX,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, typename TRI_VTYPE>
  void triangulate_subquads_using_diagonals_LU
  (const VTYPEL lower_quad_vert[2], const VTYPEU upper_quad_vert[2], 
   const VTYPEX ivertX00, const VTYPEX ivertX11,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result,
   const bool flag_reverse_orient,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);
    const TRI_VTYPE subquadA_vert[NUM_VERT_PER_QUAD] =
      { lower_quad_vert[0], lower_quad_vert[1], ivertX11, ivertX00 };
    const TRI_VTYPE subquadB_vert[NUM_VERT_PER_QUAD] =
      { ivertX00, ivertX11, upper_quad_vert[1], upper_quad_vert[0] };

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_quad_vertices(cerr, "  Subquad A: ", subquadA_vert, "\n");
    IJK::print_quad_vertices(cerr, "  Subquad B: ", subquadB_vert, "\n");
    cerr << "  quad_tri_result.ear_list: "
         << quad_tri_result.ear_list[0] << " " 
         << quad_tri_result.ear_list[1] << endl;
    */


    if (quad_tri_result.ear_list[0] == 0) {
      triangulate_quad_using_diagonal
        (subquadA_vert, false, flag_reverse_orient, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal
        (subquadA_vert, true, flag_reverse_orient, tri_vert);
    }

    if (quad_tri_result.ear_list[1] == 2) {
      triangulate_quad_using_diagonal
        (subquadB_vert, false, flag_reverse_orient, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal
        (subquadB_vert, true, flag_reverse_orient, tri_vert);
    }
  }


  /*! 
   *  Triangulate subquads of a quad where the subquads are:
   *    (quad_vert[0], quad_vert[1], ivertX12, ivertX03) and
   *    (ivertX03, ivertX12, quad_vert[2], quad_vert[3]).
   *  - Quad vertices are:
   *      (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename VTYPE, typename VTYPEX,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, typename TRI_VTYPE>
  void triangulate_subquads_using_diagonals_vX03_vX12
  (const VTYPE quad_vert[4], const VTYPEX ivertX03, const VTYPEX ivertX12,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result,
   const bool flag_reverse_orient,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);
    const TRI_VTYPE subquadA_vert[NUM_VERT_PER_QUAD] =
      { quad_vert[0], quad_vert[1], ivertX12, ivertX03 };
    const TRI_VTYPE subquadB_vert[NUM_VERT_PER_QUAD] =
      { ivertX03, ivertX12, quad_vert[2], quad_vert[3] };

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_quad_vertices(cerr, "  Subquad A: ", subquadA_vert, "\n");
    IJK::print_quad_vertices(cerr, "  Subquad B: ", subquadB_vert, "\n");
    cerr << "  quad_tri_result.ear_list: "
         << quad_tri_result.ear_list[0] << " " 
         << quad_tri_result.ear_list[1] << endl;
    */


    if (quad_tri_result.ear_list[0] == 0) {
      triangulate_quad_using_diagonal
        (subquadA_vert, false, flag_reverse_orient, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal
        (subquadA_vert, true, flag_reverse_orient, tri_vert);
    }

    if (quad_tri_result.ear_list[1] == 2) {
      triangulate_quad_using_diagonal
        (subquadB_vert, false, flag_reverse_orient, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal
        (subquadB_vert, true, flag_reverse_orient, tri_vert);
    }
  }


  /*! 
   *  Triangulate subquads of a quad where the subquads are:
   *    (quad_vert[0], ivertX01, ivertX23, quad_vert[3]) and
   *    (ivertX01, quad_vert[1], quad_vert[2], ivertX23).
   *  - Quad vertices are:
   *      (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   */
  template <typename VTYPE, typename VTYPEX,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, typename TRI_VTYPE>
  void triangulate_subquads_using_diagonals_vX01_vX23
  (const VTYPE quad_vert[4], const VTYPEX ivertX01, const VTYPEX ivertX23,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result,
   const bool flag_reverse_orient,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);
    const TRI_VTYPE subquadA_vert[NUM_VERT_PER_QUAD] =
      { quad_vert[0], ivertX01, ivertX23, quad_vert[3] };
    const TRI_VTYPE subquadB_vert[NUM_VERT_PER_QUAD] =
      { ivertX01, quad_vert[1], quad_vert[2], ivertX23 };

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_quad_vertices(cerr, "  Subquad A: ", subquadA_vert, "\n");
    IJK::print_quad_vertices(cerr, "  Subquad B: ", subquadB_vert, "\n");
    cerr << "  quad_tri_result.ear_list: "
         << quad_tri_result.ear_list[0] << " " 
         << quad_tri_result.ear_list[1] << endl;
    */

    if (quad_tri_result.ear_list[0] == 0) {
      triangulate_quad_using_diagonal
        (subquadA_vert, false, flag_reverse_orient, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal
        (subquadA_vert, true, flag_reverse_orient, tri_vert);
    }

    if (quad_tri_result.ear_list[1] == 2) {
      triangulate_quad_using_diagonal
        (subquadB_vert, false, flag_reverse_orient, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal
        (subquadB_vert, true, flag_reverse_orient, tri_vert);
    }
  }


  /*! 
   *  Triangulate subquads of a quad.
   *  - Quad vertices are:
   *      (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])
   *    listed in clockwise/counter-clockwise order around the quad.
   *  - Subquads are determined by split_edge_index.
   */
  template <typename VTYPE, typename VTYPEX, typename ITYPE,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, typename TRI_VTYPE>
  void triangulate_subquads_using_diagonals
  (const VTYPE quad_vert[4], 
   const VTYPEX ivertX01_or_X03, const VTYPEX ivertX23_or_X12,
   const ITYPE split_edge_index,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result,
   const bool flag_reverse_orient,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    if (split_edge_index == 0 || split_edge_index == 2) {
      triangulate_subquads_using_diagonals_vX01_vX23
        (quad_vert, ivertX01_or_X03, ivertX23_or_X12,
         quad_tri_result, flag_reverse_orient, tri_vert);
    }
    else {
      triangulate_subquads_using_diagonals_vX03_vX12
        (quad_vert, ivertX01_or_X03, ivertX23_or_X12,
         quad_tri_result, flag_reverse_orient, tri_vert);
    }

  }


  ///@}


  // **************************************************
  /// @name TRIANGULATE THREE QUADRILATERALS
  // **************************************************


  ///@{

  /*!
   *  Triangulate three quadrilaterals.
   *  - Break middle quad into two quads and triangulate each.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords to vertex_coord[].
   *  - Quad A has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad C has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   */
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename TRI_VTYPE, typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_three_quads_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_triA_result,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_triB_result,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_triC_result,
   const CTYPEX * vcoordX11,
   const CTYPEX * vcoordX22,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadC_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const TRI_VTYPE * lower_quadB_vert = lower_quad_vert+1;
    const TRI_VTYPE * upper_quadB_vert = upper_quad_vert+1;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  lower_quad_vert: ";
    IJK::print_list(cerr, lower_quad_vert, 4);
    cerr << endl;
    cerr << "  upper_quad_vert: ";
    IJK::print_list(cerr, upper_quad_vert, 4);
    cerr << endl;
    */

    const TRI_VTYPE ivertX11 = 
      IJK::insert_coord(dimension, vcoordX11, vertex_coord);
    const TRI_VTYPE ivertX22 = 
      IJK::insert_coord(dimension, vcoordX22, vertex_coord);

    // Triangulate subquads of quad B.
    triangulate_subquads_using_diagonals_LU
      (lower_quadB_vert, upper_quadB_vert, ivertX11, ivertX22,
       quad_triB_result, flag_reverse_orient, tri_vert);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "quad_triA_result.num_triangles: "
         << quad_triA_result.num_triangles << endl;
    cerr << "quad_tri_resultC.num_triangles: "
         << quad_triC_result.num_triangles << endl;
    */

    // Triangulate quadA with additional vertex at ivertX11.
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, upper_quad_vert[1], upper_quad_vert[0], 
       lower_quad_vert[0], lower_quad_vert[1], ivertX11, quadA_vcoord, 
       flag_reverse_orient, quad_triA_result, vertex_coord, tri_vert);

    // Triangulate quadC with additional vertex at ivertX22.
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, lower_quad_vert[2], lower_quad_vert[3], 
       upper_quad_vert[3], upper_quad_vert[2], ivertX22, quadC_vcoord, 
       flag_reverse_orient, quad_triC_result, vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate three quadrilaterals.
   *  - Version using array quad_tri_result.
   */
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename TRI_VTYPE, typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_three_quads_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri_result[3],
   const CTYPEX * vcoordX11,
   const CTYPEX * vcoordX22,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadC_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    triangulate_three_quads_allow_tri5
      (dimension, lower_quad_vert, upper_quad_vert, flag_reverse_orient,
       quad_tri_result[0], quad_tri_result[1], quad_tri_result[2],
       vcoordX11, vcoordX22, quadA_vcoord, quadC_vcoord,
       vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate three quadrilaterals where all three quads 
   *    have a common vertex (L shape.)
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords to vertex_coord[].
   *  @param lower_quad_vert[] Five vertices shared by the three quadrilaterals.
   *  @param upper_quad_vert[] Three vertices shared by the 
   *           three quadrilaterals.
   *  - Quad A has vertices (lower_quad_vert[0], lower_quad_vert[1], 
   *                         upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices (lower_quad_vert[1], lower_quad_vert[2], 
   *                         lower_quad_vert[3], upper_quad_vert[1])
   *  - Quad C has vertices (lower_quad_vert[3], lower_quad_vert[4], 
   *                         upper_quad_vert[2], upper_quad_vert[1])
  */
  template <typename DTYPE, typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX, 
            typename VTYPE0, typename VTYPE1>
  void triangulate_three_quadsL
  (const DTYPE dimension,
   const VTYPE0 * lower_quad_vert,
   const VTYPE0 * upper_quad_vert,
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[3],
   const CTYPEX * vcoordX11,
   const CTYPEX * vcoordX31,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadB_vcoord,
   const CTYPEX * quadC_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert)
  {
    const VTYPE0 ivertX11 = 
      IJK::insert_coord(dimension, vcoordX11, vertex_coord);
    const VTYPE0 ivertX31 = 
      IJK::insert_coord(dimension, vcoordX31, vertex_coord);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  lower_quad_vert: ";
    IJK::print_list(cerr, lower_quad_vert, 5);
    cerr << endl;
    cerr << "  upper_quad_vert: ";
    IJK::print_list(cerr, upper_quad_vert, 3);
    cerr << endl;
    cerr << "  quad_tri[1].num_triangles: " 
         << quad_tri[1].num_triangles << endl;
    cerr << "  quad_tri[1].num_tri_ears: " 
         << quad_tri[1].num_tri_ears << endl;
    */

    // Triangulate quadB.
    triangulate_quad_with_split_v03_v23
      (dimension, lower_quad_vert[1], lower_quad_vert[2],
       lower_quad_vert[3], upper_quad_vert[1],
       ivertX11, ivertX31, quadB_vcoord, quad_tri[1],
       flag_reverse_orient, vertex_coord, tri_vert);

    // Triangulate quadA with additional vertex at ivertX11.
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, upper_quad_vert[1], upper_quad_vert[0], 
       lower_quad_vert[0], lower_quad_vert[1], ivertX11, quadA_vcoord, 
       flag_reverse_orient, quad_tri[0], vertex_coord, tri_vert);

    // Triangulate quadC with additional vertex at ivertX31.
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, lower_quad_vert[3], lower_quad_vert[4], 
       upper_quad_vert[2], upper_quad_vert[1], ivertX31, quadC_vcoord, 
       flag_reverse_orient, quad_tri[2], vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate three quadrilaterals.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  - Consider triangulation of quad A and quad C into 5 triangles.
   *  - Consider triangulation of quad A into three triangles 
   *    by starring from vcoordX16 or iv0 or iv7.
   *  - Consider triangulation of quad C into three triangles 
   *    by starring from vcoordX25 or iv3 or iv4.
   *  - Only triangulate quads where merging with other quads 
   *    improves triangulation.
   *  - Quad A has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad C has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   *  @param block_diagonal[i][0] 
   *    Do not allow triangulation with diagonal (v0,v2) of quad i.
   *  @param block_diagonal[i][1] 
   *    Do not allow triangulation with diagonal (v1,v3) of quad i.
   *  @param[out] flag_split[] flag_split[i] is true if quad i is
   *      is split into triangles which have been added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE,
            typename VTYPEL, typename VTYPEU, typename TRI_VTYPE,
            typename MTYPE>
  void triangulate_three_quads_max_min_angle_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[3][2],
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert,
   bool flag_split[3])
  {
    const int NUM_QUADS(3);
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension };

    const TRI_VTYPE quadA_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[0], lower_quad_vert[1], 
        upper_quad_vert[1], upper_quad_vert[0] };
    const TRI_VTYPE quadB_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[1], lower_quad_vert[2], 
        upper_quad_vert[2], upper_quad_vert[1] };
    const TRI_VTYPE quadC_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[2], lower_quad_vert[3], 
        upper_quad_vert[3], upper_quad_vert[2] };

    IJK::ARRAY<CTYPE> vcoordX11(dimension);
    IJK::ARRAY<CTYPE> vcoordX22(dimension);
    IJK::ARRAY<CTYPE> quad_centroid_coord(NUM_QUADS*dimension);
    CTYPE * quadA_centroid = quad_centroid_coord.Ptr();
    CTYPE * quadB_centroid = (quad_centroid_coord.Ptr()+dimension);
    CTYPE * quadC_centroid = (quad_centroid_coord.Ptr()+2*dimension);
    const CTYPE * quad_centroid[NUM_QUADS];
    CTYPE cos_min_angle;
    POLY_TRIANGULATION_RESULTX2<16,CTYPE,int> quad_tri_result[NUM_QUADS];
    bool flag_zero;


    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  lower quad vert: ";
    IJK::print_list(cerr, lower_quad_vert, 4);
    cerr << endl;
    cerr << "  upper quad vert: ";
    IJK::print_list(cerr, upper_quad_vert, 4);
    cerr << endl;
    */

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = false;

    // Set quad_centroid[].
    quad_centroid[0] = quadA_centroid;
    quad_centroid[1] = quadB_centroid;
    quad_centroid[2] = quadC_centroid;
    
    IJK::compute_midpoint
      (dimension, lower_vcoord[1], upper_vcoord[1], vcoordX11.Ptr());
    IJK::compute_midpoint
      (dimension, lower_vcoord[2], upper_vcoord[2], vcoordX22.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadA_vert, NUM_VERT_PER_QUAD, vertex_coord, quadA_centroid);
    IJK::compute_coord_centroid
      (dimension, quadB_vert, NUM_VERT_PER_QUAD, vertex_coord, quadB_centroid);
    IJK::compute_coord_centroid
      (dimension, quadC_vert, NUM_VERT_PER_QUAD, vertex_coord, quadC_centroid);

    compute_cos_max_min_three_quads_angle_allow_tri5
      (dimension, vertex_coord, vcoordX11.PtrConst(), vcoordX22.PtrConst(),
       quad_centroid,
       lower_quad_vert, upper_quad_vert, max_small_magnitude, 
       block_diagonal, quad_tri_result, cos_min_angle, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "quad_tri_result[0].cos_min_triangulation_angle: "
         << quad_tri_result[0].cos_min_triangulation_angle << endl;
    cerr << "quad_tri_result[1].cos_min_triangulation_angle: "
         << quad_tri_result[1].cos_min_triangulation_angle << endl;
    cerr << "quad_tri_result[2].cos_min_triangulation_angle: "
         << quad_tri_result[2].cos_min_triangulation_angle << endl;
    */

    if (quad_tri_result[1].num_split_edges == 2) {

      triangulate_three_quads_allow_tri5
        (dimension, lower_quad_vert, upper_quad_vert, flag_reverse_orient, 
         quad_tri_result, vcoordX11.PtrConst(), vcoordX22.PtrConst(), 
         quadA_centroid, quadC_centroid,
         vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = true;
    }
    else if (quad_tri_result[1].num_split_edges == 1) {

      if (quad_tri_result[0].num_split_edges > 0 &&
          quad_tri_result[2].num_split_edges == 0) {

        const TRI_VTYPE ivertX11 = 
          IJK::insert_coord(dimension, vcoordX11.PtrConst(), vertex_coord);

        // Triangulate quadA with additional vertex at ivertX11.
        triangulate_quad_tri3_or_tri5_vX03
          (dimension, upper_quad_vert[1], upper_quad_vert[0], 
           lower_quad_vert[0], lower_quad_vert[1], ivertX11, 
           quadA_centroid, flag_reverse_orient, 
           quad_tri_result[0], vertex_coord, tri_vert);

        // Triangulate quadB with additional vertex at ivertX111.
        IJK::triangulate_pentagon
          (ivertX11, lower_quad_vert[1], lower_quad_vert[2], 
           upper_quad_vert[2], upper_quad_vert[1], 
           flag_reverse_orient, tri_vert);

        flag_split[0] = flag_split[1] = true;
      }
      else if (quad_tri_result[0].num_split_edges == 0 &&
               quad_tri_result[2].num_split_edges == 1) {

        // Use triangulation around X22
        const TRI_VTYPE ivertX22 = IJK::insert_coord
          (dimension, vcoordX22.PtrConst(), vertex_coord);

        // Triangulate quadB with additional vertex at ivertX22.
        IJK::triangulate_pentagon
          (ivertX22, upper_quad_vert[2], upper_quad_vert[1], 
           lower_quad_vert[1], lower_quad_vert[2], 
           flag_reverse_orient, tri_vert);

        // Triangulate quadC with additional vertex at ivertX22.
        triangulate_quad_tri3_or_tri5_vX03
          (dimension, lower_quad_vert[2], lower_quad_vert[3], 
           upper_quad_vert[3], upper_quad_vert[2], 
           ivertX22, quadC_centroid, flag_reverse_orient, 
           quad_tri_result[2], vertex_coord, tri_vert);

        flag_split[1] = flag_split[2] = true;
      }
    }

  }


  /*!
   *  Triangulate three quadrilaterals where all three quads 
   *    have a common vertex (L shape.)
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  - Consider triangulation of quad A and quad C into 5 triangles.
   *  - Consider triangulation of quad A into three triangles iv0 or iv7.
   *  - Consider triangulation of quad C into three triangles iv3 or iv4.
   *  - Only triangulate quads where merging with other quads 
   *    improves triangulation.
   *  @param lower_quad_vert[] Five vertices shared by the three quadrilaterals.
   *  @param upper_quad_vert[] Three vertices shared by the 
   *           three quadrilaterals.
   *  - Quad A has vertices (lower_quad_vert[0], lower_quad_vert[1], 
   *                         upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices (lower_quad_vert[1], lower_quad_vert[2], 
   *                         lower_quad_vert[3], upper_quad_vert[1])
   *  - Quad C has vertices (lower_quad_vert[3], lower_quad_vert[4], 
   *                         upper_quad_vert[2], upper_quad_vert[1])
   *  @param block_diagonal[i][0] 
   *    Do not allow triangulation with diagonal (v0,v2) of quad i.
   *  @param block_diagonal[i][1] 
   *    Do not allow triangulation with diagonal (v1,v3) of quad i.
   *  @param[out] flag_split[] flag_split[i] is true if quad i is
   *      is split into triangles which have been added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2, typename MTYPE>
  void triangulate_three_quadsL_max_min_angle
  (const DTYPE dimension,
   const VTYPE0 * lower_quad_vert,
   const VTYPE1 * upper_quad_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[3][2],
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE2> & tri_vert,
   bool flag_split[3])
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_QUADS(3);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoordL1 = vertex_coord_ptr+lower_quad_vert[1]*dimension;
    const CTYPE * vcoordL3 = vertex_coord_ptr+lower_quad_vert[3]*dimension;
    const CTYPE * vcoordU1 = vertex_coord_ptr+upper_quad_vert[1]*dimension;
    const VTYPE0 quadA_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[0], lower_quad_vert[1], 
        upper_quad_vert[1], upper_quad_vert[0] };
    const VTYPE0 quadB_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[1], lower_quad_vert[2], 
        lower_quad_vert[3], upper_quad_vert[1] };
    const VTYPE0 quadC_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[3], lower_quad_vert[4], 
        upper_quad_vert[2], upper_quad_vert[1] };
    IJK::ARRAY<CTYPE> vcoordX11(dimension);
    IJK::ARRAY<CTYPE> vcoordX31(dimension);
    IJK::ARRAY<CTYPE> quadA_centroid(dimension);
    IJK::ARRAY<CTYPE>quadB_centroid(dimension);
    IJK::ARRAY<CTYPE> quadC_centroid(dimension);
    CTYPE cos_min_angle;
    POLY_TRIANGULATION_RESULTX2<16,CTYPE,int> quad_tri[NUM_QUADS];
    bool flag_zero;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << endl;
    cerr << "In " << __func__ << endl;
    cerr << "  lower_quad_vert: ";
    IJK::print_list(cerr, lower_quad_vert, 5);
    cerr << endl;
    cerr << "  upper_quad_vert: ";
    IJK::print_list(cerr, upper_quad_vert, 3);
    cerr << endl;
    */

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = false;

    IJK::compute_midpoint(dimension, vcoordL1, vcoordU1, vcoordX11.Ptr());
    IJK::compute_midpoint(dimension, vcoordL3, vcoordU1, vcoordX31.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadA_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadA_centroid.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadB_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadB_centroid.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadC_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadC_centroid.Ptr());

    compute_cos_max_min_three_quadsL_angle
      (dimension, vertex_coord, vcoordX11.PtrConst(), vcoordX31.PtrConst(),
       quadA_centroid.PtrConst(), quadB_centroid.PtrConst(), 
       quadC_centroid.PtrConst(),
       lower_quad_vert, upper_quad_vert, max_small_magnitude, block_diagonal,
       quad_tri, cos_min_angle, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "cos_min_angle: " << cos_min_angle << endl;
    */

    if (quad_tri[1].num_split_edges == 2) {

      triangulate_three_quadsL
        (dimension, lower_quad_vert, upper_quad_vert,
         flag_reverse_orient, quad_tri,
         vcoordX11.PtrConst(), vcoordX31.PtrConst(), 
         quadA_centroid.PtrConst(), quadB_centroid.PtrConst(),
         quadC_centroid.PtrConst(),
         vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = true;
    }
    else if (quad_tri[1].num_split_edges == 1) {

      if (quad_tri[0].num_split_edges == 1 &&
          quad_tri[2].num_split_edges == 0) {

        const VTYPE0 ivertX11 = 
          IJK::insert_coord(dimension, vcoordX11.PtrConst(), vertex_coord);

        // Triangulate quadA with additional vertex at ivertX11.
        triangulate_quad_tri3_or_tri5_vX03
          (dimension, upper_quad_vert[1], upper_quad_vert[0], 
           lower_quad_vert[0], lower_quad_vert[1], 
           ivertX11, quadA_centroid.PtrConst(), 
           flag_reverse_orient, quad_tri[0], vertex_coord, tri_vert);


        // Triangulate quadB with additional vertex at ivertX11.
        triangulate_quad_tri3_or_tri5_vX03
          (dimension, lower_quad_vert[1], lower_quad_vert[2], 
           lower_quad_vert[3], upper_quad_vert[1],
           ivertX11, quadB_centroid.PtrConst(), 
           flag_reverse_orient, quad_tri[1], vertex_coord, tri_vert);

        flag_split[0] = flag_split[1] = true;
      }
      else if (quad_tri[0].num_split_edges == 0 &&
               quad_tri[2].num_split_edges == 1) {

        // Use triangulation around X31
        const VTYPE0 ivertX31 = IJK::insert_coord
          (dimension, vcoordX31.PtrConst(), vertex_coord);

        // Triangulate quadB with additional vertex at ivertX31.
        triangulate_quad_tri3_or_tri5_vX03
          (dimension, upper_quad_vert[1], lower_quad_vert[1], 
           lower_quad_vert[2], lower_quad_vert[3], 
           ivertX31, quadB_centroid.PtrConst(), 
           flag_reverse_orient, quad_tri[1], vertex_coord, tri_vert);


        // Triangulate quadC with additional vertex at ivertX31.
        triangulate_quad_tri3_or_tri5_vX03
          (dimension, lower_quad_vert[3], lower_quad_vert[4], 
           upper_quad_vert[2], upper_quad_vert[1],
           ivertX31, quadC_centroid.PtrConst(), 
           flag_reverse_orient, quad_tri[2], vertex_coord, tri_vert);


        flag_split[1] = flag_split[2] = true;
      }
    }

  }

  ///@}


  // **************************************************
  /// @name TRIANGULATE FOUR QUADRILATERALS
  // **************************************************

  ///@{

  /*!
   *  Triangulate four quadrilaterals.
   *  - Break each of the middle two quads into two smaller quads 
   *      and triangulate each.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords to vertex_coord[].
   *  - Quad A has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad C has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   *  - Quad D has vertices 
   *      (lower_quad_vert[3], lower_quad_vert[4], 
   *       upper_quad_vert[4], upper_quad_vert[3])
   *  @param vcoordX[] Coordinates of three new points on edges.
   *       (Usually midpoints of those edges.)
   *  - First point is on edge (lower_quad_vert[1], upper_quad_vert[1]).
   *  - Second point is on edge (lower_quad_vert[2], upper_quad_vert[2]).
   *  - Third point is on edge (lower_quad_vert[3], upper_quad_vert[3]).
  */
  template <typename DTYPE, 
            typename VTYPEL, typename VTYPEU, typename TRI_VTYPE,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_four_quads_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[4],
   const CTYPEX * vcoordX,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadD_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    TRI_VTYPE ivertX[3];
    const TRI_VTYPE * lower_quadB_vert = lower_quad_vert+1;
    const TRI_VTYPE * upper_quadB_vert = upper_quad_vert+1;
    const TRI_VTYPE * lower_quadC_vert = lower_quad_vert+2;
    const TRI_VTYPE * upper_quadC_vert = upper_quad_vert+2;

    ivertX[0] = 
      IJK::insert_coord(dimension, vcoordX, vertex_coord);
    ivertX[1] = 
      IJK::insert_coord(dimension, vcoordX+dimension, vertex_coord);
    ivertX[2] = 
      IJK::insert_coord(dimension, vcoordX+2*dimension, vertex_coord);

    // Triangulate subquads of quad B.
    triangulate_subquads_using_diagonals_LU
      (lower_quadB_vert, upper_quadB_vert, ivertX[0], ivertX[1],
       quad_tri[1], flag_reverse_orient, tri_vert);

    // Triangulate subquads of quad C.
    triangulate_subquads_using_diagonals_LU
      (lower_quadC_vert, upper_quadC_vert, ivertX[1], ivertX[2],
       quad_tri[2], flag_reverse_orient, tri_vert);

    // Triangulate quadA with additional vertex at ivertX11.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX[0], upper_quad_vert[1], upper_quad_vert[0], 
       lower_quad_vert[0], lower_quad_vert[1], quadA_vcoord, quad_tri[0],
       flag_reverse_orient, vertex_coord, tri_vert);

    // Triangulate quadD with additional vertex at ivertX22.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX[2], lower_quad_vert[3], lower_quad_vert[4], 
       upper_quad_vert[4], upper_quad_vert[3], quadD_vcoord, quad_tri[3],
       flag_reverse_orient, vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate four quadrilaterals.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  - Consider triangulation of quads into 5 triangles.
   *  - Only triangulate quads where merging with other quads 
   *    improves triangulation.
   *  - Quad A has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad B has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad C has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   *  - Quad D has vertices 
   *      (lower_quad_vert[3], lower_quad_vert[4], 
   *       upper_quad_vert[4], upper_quad_vert[3])
   *  @param[out] flag_split[] flag_split[i] is true if quad i is
   *      is split into triangles which have been added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE,
            typename VTYPEL, typename VTYPEU, typename TRI_VTYPE,
            typename MTYPE>
  void triangulate_four_quads_max_min_angle_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert,
   bool flag_split[4])
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_QUADS(4);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension,
        vertex_coord_ptr+lower_quad_vert[4]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension,
        vertex_coord_ptr+upper_quad_vert[4]*dimension };
    const TRI_VTYPE quadA_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[0], lower_quad_vert[1], 
        upper_quad_vert[1], upper_quad_vert[0] };
    const TRI_VTYPE quadB_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[1], lower_quad_vert[2], 
        upper_quad_vert[2], upper_quad_vert[1] };
    const TRI_VTYPE quadC_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[2], lower_quad_vert[3], 
        upper_quad_vert[3], upper_quad_vert[2] };
    const TRI_VTYPE quadD_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[3], lower_quad_vert[4], 
        upper_quad_vert[4], upper_quad_vert[3] };
    IJK::ARRAY<CTYPE> quad_centroid(NUM_QUADS*dimension);
    CTYPE cos_min_angle;
    POLY_TRIANGULATION_RESULTX2<16,CTYPE,int> quad_tri_result[NUM_QUADS];
    bool flag_zero;

    // Coordinates of edge midpoints.
    // - First point is on edge (lower_quad_vert[1], upper_quad_vert[1]).
    // - Second point is on edge (lower_quad_vert[2], upper_quad_vert[2]).
    // - Third point is on edge (lower_quad_vert[3], upper_quad_vert[3]).
    IJK::ARRAY<CTYPE> vcoordX(3*dimension);

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = false;

    CTYPE * vcoordX11 = vcoordX.Ptr();
    CTYPE * vcoordX22 = vcoordX.Ptr()+dimension;
    CTYPE * vcoordX33 = vcoordX.Ptr()+2*dimension;

    IJK::compute_midpoint
      (dimension, lower_vcoord[1], upper_vcoord[1], vcoordX11);
    IJK::compute_midpoint
      (dimension, lower_vcoord[2], upper_vcoord[2], vcoordX22);
    IJK::compute_midpoint
      (dimension, lower_vcoord[3], upper_vcoord[3], vcoordX33);

    CTYPE * quadA_centroid = quad_centroid.Ptr();
    CTYPE * quadB_centroid = quad_centroid.Ptr()+dimension;
    CTYPE * quadC_centroid = quad_centroid.Ptr()+2*dimension;
    CTYPE * quadD_centroid = quad_centroid.Ptr()+3*dimension;

    IJK::compute_coord_centroid
      (dimension, quadA_vert, NUM_VERT_PER_QUAD, vertex_coord, quadA_centroid);
    IJK::compute_coord_centroid
      (dimension, quadB_vert, NUM_VERT_PER_QUAD, vertex_coord, quadB_centroid);
    IJK::compute_coord_centroid
      (dimension, quadC_vert, NUM_VERT_PER_QUAD, vertex_coord, quadC_centroid);
    IJK::compute_coord_centroid
      (dimension, quadD_vert, NUM_VERT_PER_QUAD, vertex_coord, quadD_centroid);

    compute_cos_max_min_four_quads_angle_allow_tri5
    (dimension, vertex_coord, vcoordX.PtrConst(), quad_centroid.PtrConst(),
     lower_quad_vert, upper_quad_vert, max_small_magnitude, cos_min_angle, 
     quad_tri_result, flag_zero);


    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  quad_tri_result[1].num_triangles: " 
         << quad_tri_result[1].num_triangles << endl;
    cerr << "  quad_tri_result[2].num_triangles: " 
         << quad_tri_result[2].num_triangles << endl;
    */


    // *** DEBUG ***
    /*
    using namespace std;
    if (quad_tri_result[1].num_triangles == 4 ||
        quad_tri_result[2].num_triangles == 4) {
      cerr << "In " << __func__ << endl;
      cerr << "  lower quad vert: ";
      IJK::print_list(cerr, lower_quad_vert, 5);
      cerr << endl;
      cerr << "  upper quad vert: ";
      IJK::print_list(cerr, upper_quad_vert, 5);
      cerr << endl;
      cerr << "  cos_min_angle: " << cos_min_angle << endl;
      cerr << "  quad_tri_result[0].num_triangles: " 
           << quad_tri_result[0].num_triangles << endl;
      cerr << "  quad_tri_result[0].tri_vertex_index: "
           << quad_tri_result[0].tri_vertex_index << endl;
      cerr << "  quad_tri_result[1].num_triangles: " 
           << quad_tri_result[1].num_triangles << endl;  \
      cerr << "  quad_tri_result[2].num_triangles: " 
           << quad_tri_result[2].num_triangles << endl;
      cerr << "  quad_tri_result[3].num_triangles: " 
           << quad_tri_result[3].num_triangles << endl;
      cerr << "  quad_tri_result[3].tri_vertex_index: "
           << quad_tri_result[3].tri_vertex_index << endl;
      cerr << "  quad_tri_result[1].ear_list: "
           << quad_tri_result[1].ear_list[0] << "  "
           << quad_tri_result[1].ear_list[1] << endl;
      cerr << "  quad_tri_result[2].ear_list: "
           << quad_tri_result[2].ear_list[0] << "  "
           << quad_tri_result[2].ear_list[1] << endl;
    }
    */



    if (quad_tri_result[1].num_split_edges == 2 &&
        quad_tri_result[2].num_split_edges == 2) {

      triangulate_four_quads_allow_tri5
        (dimension, lower_quad_vert, upper_quad_vert, flag_reverse_orient, 
         quad_tri_result, vcoordX.PtrConst(), quadA_centroid, quadD_centroid, 
         vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = true;
    }
    else if (quad_tri_result[0].num_split_edges == 1 &&
             quad_tri_result[1].num_split_edges == 1 &&
             quad_tri_result[2].num_split_edges == 1 &&
             quad_tri_result[3].num_split_edges == 1) {

      const VTYPEL * lower_quadAB_vert = lower_quad_vert;
      const VTYPEU * upper_quadAB_vert = upper_quad_vert;
      const VTYPEL * lower_quadCD_vert = lower_quad_vert+2;
      const VTYPEU * upper_quadCD_vert = upper_quad_vert+2;

      bool flag_quadAB_tri5[2];
      bool flag_quadCD_tri5[2];

      flag_quadAB_tri5[0] = (quad_tri_result[0].num_triangles == 5);
      flag_quadAB_tri5[1] = (quad_tri_result[1].num_triangles == 5);
      flag_quadCD_tri5[0] = (quad_tri_result[2].num_triangles == 5);
      flag_quadCD_tri5[1] = (quad_tri_result[3].num_triangles == 5);

      triangulate_two_quads_allow_tri5
        (dimension, lower_quadAB_vert, upper_quadAB_vert, 
         flag_reverse_orient, flag_quadAB_tri5, vcoordX11, 
         quadA_centroid, quadB_centroid, vertex_coord, tri_vert);

      triangulate_two_quads_allow_tri5
        (dimension, lower_quadCD_vert, upper_quadCD_vert,
         flag_reverse_orient, flag_quadCD_tri5, vcoordX33, 
         quadC_centroid, quadD_centroid, vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = true;
    }
    else if (quad_tri_result[0].num_split_edges == 1 &&
             quad_tri_result[1].num_split_edges == 2 &&
             quad_tri_result[2].num_split_edges == 1 &&
             quad_tri_result[3].num_split_edges == 0) {

      // **** DEBUG ***
      /*
      using namespace std;
      cerr << "  Calling triangulate_three_quads_allow_tri5 on ABC." << endl;
      */

      triangulate_three_quads_allow_tri5
        (dimension, lower_quad_vert, upper_quad_vert, flag_reverse_orient, 
         quad_tri_result, vcoordX11, vcoordX22, quadA_centroid, quadC_centroid,
         vertex_coord, tri_vert);


      flag_split[0] = flag_split[1] = flag_split[2] = true;
      flag_split[3] = false;
    }
    else if (quad_tri_result[0].num_split_edges == 0 &&
             quad_tri_result[1].num_split_edges == 1 &&
             quad_tri_result[2].num_split_edges == 2 &&
             quad_tri_result[3].num_split_edges == 1) {

      // **** DEBUG ***
      /*
      using namespace std;
      cerr << "  Calling triangulate_three_quads_allow_tri5 on BCD." << endl;
      */

      const VTYPEL * lower_quadBCD_vert = lower_quad_vert+1;
      const VTYPEU * upper_quadBCD_vert = upper_quad_vert+1;

      triangulate_three_quads_allow_tri5
        (dimension, lower_quadBCD_vert, upper_quadBCD_vert, 
         flag_reverse_orient, quad_tri_result+1, vcoordX22, vcoordX33, 
         quadB_centroid, quadD_centroid, vertex_coord, tri_vert);

      flag_split[0] = false;
      flag_split[1] = flag_split[2] = flag_split[3] = true;
    }
    else if (quad_tri_result[1].num_split_edges == 1 &&
             quad_tri_result[2].num_split_edges == 1) {

      // **** DEBUG ***
      /*
      using namespace std;
      cerr << "  Calling triangulate_two_quads_allow_tri5 on BC." << endl;
      */

      const VTYPEL * lower_quadBC_vert = lower_quad_vert+1;
      const VTYPEU * upper_quadBC_vert = upper_quad_vert+1;

      bool flag_quad_tri5[2];

      flag_quad_tri5[0] = (quad_tri_result[1].num_triangles == 5);
      flag_quad_tri5[1] = (quad_tri_result[2].num_triangles == 5);

      triangulate_two_quads_allow_tri5
        (dimension, lower_quadBC_vert, upper_quadBC_vert, 
         flag_reverse_orient, flag_quad_tri5, vcoordX22, 
         quadB_centroid, quadC_centroid, vertex_coord, tri_vert);

      flag_split[0] = flag_split[3] = false;
      flag_split[1] = flag_split[2] = true;
    }

  }


  /*!
   *  Triangulate quad and three neighboring quads (T shape).
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords to vertex_coord[].
   *  @param vcoordX[] Coordinates of three new points on edges.
   *       (Usually midpoints of those edges.)
   *  - First point is on edge (lower_quad_vert[1], upper_quad_vert[1]).
   *  - Second point is on edge (lower_quad_vert[2], upper_quad_vert[2]).
   *  - Third point is on edge (lower_quad_vert[3], upper_quad_vert[3]).
  */
  template <typename DTYPE, typename VTYPE, typename TRI_VTYPE,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_four_quadsT_allow_tri5
  (const DTYPE dimension,
   const VTYPE quad_vlist[10],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri[4],
   const CTYPEX * vcoordX,
   const CTYPEX * quad_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const CTYPEX * quadA_vcoordX = quad_vcoord;
    const CTYPEX * quadB_vcoordX = quad_vcoord+dimension;
    const CTYPEX * quadC_vcoordX = quad_vcoord+2*dimension;
    const CTYPEX * quadD_vcoordX = quad_vcoord+3*dimension;
    TRI_VTYPE ivertX[3];
    TRI_VTYPE iv_quadA;

    ivertX[0] = 
      IJK::insert_coord(dimension, vcoordX, vertex_coord);
    ivertX[1] = 
      IJK::insert_coord(dimension, vcoordX+dimension, vertex_coord);
    ivertX[2] = 
      IJK::insert_coord(dimension, vcoordX+2*dimension, vertex_coord);
    iv_quadA =
      IJK::insert_coord(dimension, quadA_vcoordX, vertex_coord);

    // Triangulate quadB with additional vertex at quadB_vcoordX.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX[0], quad_vlist[0], quad_vlist[1], 
       quad_vlist[2], quad_vlist[3], quadB_vcoordX, quad_tri[1],
       flag_reverse_orient, vertex_coord, tri_vert);

    // Triangulate quadC with additional vertex at quadC_vcoordX.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX[1], quad_vlist[3], quad_vlist[4], 
       quad_vlist[5], quad_vlist[6], quadC_vcoordX, quad_tri[2],
       flag_reverse_orient, vertex_coord, tri_vert);

    // Triangulate quadD with additional vertex at quadD_vcoordX.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX[2], quad_vlist[6], quad_vlist[7], 
       quad_vlist[8], quad_vlist[9], quadD_vcoordX, quad_tri[3],
       flag_reverse_orient, vertex_coord, tri_vert);

    // Triangulate quadA with additional vertex at quadA_vcoord.
    triangulate_septagon_with_vertex
      (quad_vlist[0], ivertX[0], quad_vlist[3], ivertX[1],
       quad_vlist[6], ivertX[2], quad_vlist[9], iv_quadA, 
       flag_reverse_orient, tri_vert);
  }


  /*!
   *  Triangulate quad and three neighboring quads (T shape).
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  - Consider triangulation of central quad into 6 triangles.
   *  - Consider triangulation of neighboring quads into 5 triangles.
   *  - Only triangulate quads where merging with other quads 
   *    improves triangulation.
   *  - Quad A has vertices 
   *      (quad_vlist[0], quad_vlist[3], quad_vlist[6], quad_vlist[9]).
   *  - Quad B has vertices 
   *      (quad_vlist[0], quad_vlist[1], quad_vlist[2], quad_vlist[3]).
   *  - Quad C has vertices 
   *      (quad_vlist[3], quad_vlist[4], quad_vlist[5], quad_vlist[6]).
   *  - Quad D has vertices 
   *      (quad_vlist[6], quad_vlist[7], quad_vlist[8], quad_vlist[9]).
   *  @param[out] flag_split[] flag_split[i] is true if quad i is split
   *      into triangles which have been added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE, typename VTYPE, typename TRI_VTYPE,
            typename MTYPE>
  void triangulate_four_quadsT_max_min_angle_allow_tri5
  (const DTYPE dimension,
   const VTYPE quad_vlist[10],
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert,
   bool flag_split[4])
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_QUADS(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * quad_vlist_coord[10] =
      { vertex_coord_ptr+quad_vlist[0]*dimension,
        vertex_coord_ptr+quad_vlist[1]*dimension,
        vertex_coord_ptr+quad_vlist[2]*dimension,
        vertex_coord_ptr+quad_vlist[3]*dimension,
        vertex_coord_ptr+quad_vlist[4]*dimension,
        vertex_coord_ptr+quad_vlist[5]*dimension,
        vertex_coord_ptr+quad_vlist[6]*dimension,
        vertex_coord_ptr+quad_vlist[7]*dimension,
        vertex_coord_ptr+quad_vlist[8]*dimension,
        vertex_coord_ptr+quad_vlist[9]*dimension };
    const TRI_VTYPE quadA_vert[NUM_VERT_PER_QUAD] = 
      { quad_vlist[0], quad_vlist[3], quad_vlist[6], quad_vlist[9] };
    const TRI_VTYPE quadB_vert[NUM_VERT_PER_QUAD] = 
      { quad_vlist[0], quad_vlist[1], quad_vlist[2], quad_vlist[3] };
    const TRI_VTYPE quadC_vert[NUM_VERT_PER_QUAD] = 
      { quad_vlist[3], quad_vlist[4], quad_vlist[5], quad_vlist[6] };
    const TRI_VTYPE quadD_vert[NUM_VERT_PER_QUAD] = 
      { quad_vlist[6], quad_vlist[7], quad_vlist[8], quad_vlist[9] };
    IJK::ARRAY<CTYPE> quad_centroid(NUM_QUADS*dimension);
    CTYPE cos_min_angle;
    POLY_TRIANGULATION_RESULTX2<16,CTYPE,int> 
      quad_tri_result[NUM_QUADS];
    bool flag_zero;
    bool flag_quad_tri5[2];

    // Coordinates of edge midpoints.
    // - First point is on edge (quad_vlist[0], quad_vlist[3]).
    // - Second point is on edge (quad_vlist[3], quad_vlist[6]).
    // - Third point is on edge (quad_vlist[6], quad_vlist[9]).
    IJK::ARRAY<CTYPE> vcoordX(3*dimension);

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = false;

    CTYPE * vcoordX03 = vcoordX.Ptr();
    CTYPE * vcoordX36 = vcoordX.Ptr()+dimension;
    CTYPE * vcoordX69 = vcoordX.Ptr()+2*dimension;

    IJK::compute_midpoint
      (dimension, quad_vlist_coord[0], quad_vlist_coord[3], vcoordX03);
    IJK::compute_midpoint
      (dimension, quad_vlist_coord[3], quad_vlist_coord[6], vcoordX36);
    IJK::compute_midpoint
      (dimension, quad_vlist_coord[6], quad_vlist_coord[9], vcoordX69);

    CTYPE * quadA_centroid = quad_centroid.Ptr();
    CTYPE * quadB_centroid = quad_centroid.Ptr()+dimension;
    CTYPE * quadC_centroid = quad_centroid.Ptr()+2*dimension;
    CTYPE * quadD_centroid = quad_centroid.Ptr()+3*dimension;

    IJK::compute_coord_centroid
      (dimension, quadA_vert, NUM_VERT_PER_QUAD, vertex_coord, quadA_centroid);
    IJK::compute_coord_centroid
      (dimension, quadB_vert, NUM_VERT_PER_QUAD, vertex_coord, quadB_centroid);
    IJK::compute_coord_centroid
      (dimension, quadC_vert, NUM_VERT_PER_QUAD, vertex_coord, quadC_centroid);
    IJK::compute_coord_centroid
      (dimension, quadD_vert, NUM_VERT_PER_QUAD, vertex_coord, quadD_centroid);

    compute_cos_max_min_four_quadsT_angle_allow_tri5
    (dimension, vertex_coord, vcoordX.PtrConst(), quad_centroid.PtrConst(),
     quad_vlist, max_small_magnitude, quad_tri_result,
     cos_min_angle, flag_zero);

    if (quad_tri_result[0].num_split_edges == 3) {
      triangulate_four_quadsT_allow_tri5
        (dimension, quad_vlist, flag_reverse_orient, 
         quad_tri_result, vcoordX.PtrConst(), quad_centroid.PtrConst(),
         vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = true;
    }
    else if (quad_tri_result[0].num_split_edges > 0) {

      if (quad_tri_result[1].num_split_edges  > 0 &&
          quad_tri_result[2].num_split_edges  == 0 &&
          quad_tri_result[3].num_split_edges  == 0) {

        const VTYPE lower_quad_vert[3] = 
          { quad_vlist[9], quad_vlist[0], quad_vlist[1] };
        const VTYPE upper_quad_vert[3] = 
          { quad_vlist[6], quad_vlist[3], quad_vlist[2] };
        flag_quad_tri5[0] = (quad_tri_result[0].num_triangles == 5);
        flag_quad_tri5[1] = (quad_tri_result[1].num_triangles == 5);

        triangulate_two_quads_allow_tri5
          (dimension, lower_quad_vert, upper_quad_vert,
           flag_reverse_orient, flag_quad_tri5, vcoordX03, 
           quadA_centroid, quadB_centroid, vertex_coord, tri_vert);

        flag_split[0] = flag_split[1] = true;
        flag_split[2] = flag_split[3] = false;
      }
      else if (quad_tri_result[1].num_split_edges == 0 &&
               quad_tri_result[2].num_split_edges  > 0 &&
               quad_tri_result[3].num_split_edges  == 0) {

        // Currently no split.
        flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = false;
      }
      else if (quad_tri_result[1].num_split_edges == 0 &&
               quad_tri_result[2].num_split_edges == 0 &&
               quad_tri_result[3].num_split_edges  > 0) {

        const VTYPE lower_quad_vert[3] = 
          { quad_vlist[3], quad_vlist[6], quad_vlist[7] };
        const VTYPE upper_quad_vert[3] = 
          { quad_vlist[0], quad_vlist[9], quad_vlist[8] };
        flag_quad_tri5[0] = (quad_tri_result[0].num_triangles == 5);
        flag_quad_tri5[1] = (quad_tri_result[3].num_triangles == 5);

        triangulate_two_quads_allow_tri5
          (dimension, lower_quad_vert, upper_quad_vert,
           flag_reverse_orient, flag_quad_tri5, vcoordX69, 
           quadA_centroid, quadD_centroid, vertex_coord, tri_vert);

        flag_split[0] = flag_split[3] = true;
        flag_split[1] = flag_split[2] = false;
      }
      else if (quad_tri_result[1].num_split_edges > 0 &&
               quad_tri_result[2].num_split_edges == 0 &&
               quad_tri_result[3].num_split_edges > 0) {

        const VTYPE lower_quad_vert[4] = 
          { quad_vlist[8], quad_vlist[9], quad_vlist[0], quad_vlist[1] };
        const VTYPE upper_quad_vert[4] = 
          { quad_vlist[7], quad_vlist[6], quad_vlist[3], quad_vlist[2] };

        triangulate_three_quads_allow_tri5
          (dimension, lower_quad_vert, upper_quad_vert, flag_reverse_orient, 
           quad_tri_result[3], quad_tri_result[0], quad_tri_result[1], 
           vcoordX69, vcoordX03, quadD_centroid, quadB_centroid,
           vertex_coord, tri_vert);

        flag_split[0] = flag_split[1] = flag_split[3] = true;
        flag_split[2] = false;
      }
      else {
        flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = false;
      }
    }
    else {
      flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = false;
    }

  }


  /*!
   *  Triangulate four polygons in T shape, a quad and an adjacent
   *    triangle, quad and triangle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords to vertex_coord[].
   *  @param vcoordX[] Coordinates of three new points on edges.
   *       (Usually midpoints of those edges.)
  */
  template <typename DTYPE, typename VTYPE, typename TRI_VTYPE,
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_four_polyT_QTQT_max_min_angle
  (const DTYPE dimension,
   const VTYPE poly_vlist[8],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   poly_tri[4],
   const CTYPEX * vcoordX,
   const CTYPEX quadA_vcoordX[],
   const CTYPEX quadC_vcoordX[],
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    TRI_VTYPE ivertX[3];
    TRI_VTYPE iv_quadA;

    ivertX[0] = 
      IJK::insert_coord(dimension, vcoordX, vertex_coord);
    ivertX[1] = 
      IJK::insert_coord(dimension, vcoordX+dimension, vertex_coord);
    ivertX[2] = 
      IJK::insert_coord(dimension, vcoordX+2*dimension, vertex_coord);
    iv_quadA =
      IJK::insert_coord(dimension, quadA_vcoordX, vertex_coord);

    if (poly_tri[0].num_tri_ears == 0) {

      // Triangulate quadA with additional vertex at quadA_vcoord.
      triangulate_septagon_with_vertex
        (poly_vlist[0], ivertX[0], poly_vlist[2], ivertX[1],
         poly_vlist[5], ivertX[2], poly_vlist[7], iv_quadA, 
         flag_reverse_orient, tri_vert);
    }
    else {
      // Triangulate quadA with additional vertex at quadA_vcoord
      //   and ears at poly_vlist[1] and poly_vlist[2].
      IJK::triangulate_pentagon_with_vertex
        (poly_vlist[0], ivertX[0], ivertX[1], ivertX[2], poly_vlist[7], 
         iv_quadA, flag_reverse_orient, tri_vert);

      // Add ears.
      add_triangle_vertices
        (ivertX[0], poly_vlist[2], ivertX[1], flag_reverse_orient, tri_vert);
      add_triangle_vertices
        (ivertX[1], poly_vlist[5], ivertX[2], flag_reverse_orient, tri_vert);
    }

    // Triangulate quadC with additional vertex at quadC_vcoordX.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX[1], poly_vlist[2], poly_vlist[3], 
       poly_vlist[4], poly_vlist[5], quadC_vcoordX, poly_tri[2],
       flag_reverse_orient, vertex_coord, tri_vert);

    // Split triangleB into two triangles.
    add_triangle_vertices
      (ivertX[0], poly_vlist[0], poly_vlist[1], flag_reverse_orient, tri_vert);
    add_triangle_vertices
      (ivertX[0], poly_vlist[1], poly_vlist[2], flag_reverse_orient, tri_vert);

    // Split triangleD into two triangles.
    add_triangle_vertices
      (ivertX[2], poly_vlist[5], poly_vlist[6], flag_reverse_orient, tri_vert);
    add_triangle_vertices
      (ivertX[2], poly_vlist[6], poly_vlist[7], flag_reverse_orient, tri_vert);
  }


  /*!
   *  Triangulate four polygons in T shape, a quad and an adjacent
   *    triangle, quad and triangle.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  - Consider triangulation of central quad into 6 triangles.
   *  - Consider triangulation of neighboring quads into 5 triangles.
   *  - Only triangulate quads where merging with other quads 
   *    improves triangulation.
   *  - Quad A has vertices 
   *      (poly_vlist[0], poly_vlist[2], poly_vlist[5], poly_vlist[7]).
   *  - Triangle B has vertices 
   *      (poly_vlist[0], poly_vlist[1], poly_vlist[2]).
   *  - Quad C has vertices 
   *      (poly_vlist[2], poly_vlist[3], poly_vlist[4], poly_vlist[5]).
   *  - Triangle D has vertices 
   *      (poly_vlist[5], poly_vlist[6], poly_vlist[7]).
   *  @param[out] flag_split[] flag_split[i] is true if poly i is split
   *      into triangles which have been added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE, typename VTYPE, typename TRI_VTYPE,
            typename MTYPE>
  void triangulate_four_polyT_QTQT_max_min_angle
  (const DTYPE dimension,
   const VTYPE poly_vlist[8],
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert,
   bool flag_split[4])
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_POLY(4);
    const int VLIST_LENGTH(8);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * poly_vlist_coord[VLIST_LENGTH] =
      { vertex_coord_ptr+poly_vlist[0]*dimension,
        vertex_coord_ptr+poly_vlist[1]*dimension,
        vertex_coord_ptr+poly_vlist[2]*dimension,
        vertex_coord_ptr+poly_vlist[3]*dimension,
        vertex_coord_ptr+poly_vlist[4]*dimension,
        vertex_coord_ptr+poly_vlist[5]*dimension,
        vertex_coord_ptr+poly_vlist[6]*dimension,
        vertex_coord_ptr+poly_vlist[7]*dimension };
    const TRI_VTYPE quadA_vert[NUM_VERT_PER_QUAD] = 
      { poly_vlist[0], poly_vlist[2], poly_vlist[5], poly_vlist[7] };
    const TRI_VTYPE quadC_vert[NUM_VERT_PER_QUAD] = 
      { poly_vlist[2], poly_vlist[3], poly_vlist[4], poly_vlist[5] };
    IJK::ARRAY<CTYPE> quadA_centroid(dimension);
    IJK::ARRAY<CTYPE> quadC_centroid(dimension);
    CTYPE cos_min_angle;
    POLY_TRIANGULATION_RESULTX2<16,CTYPE,int> poly_tri_result[NUM_POLY];
    bool flag_zero;

    // Coordinates of edge midpoints.
    // - First point is on edge (poly_vlist[0], poly_vlist[2]).
    // - Second point is on edge (poly_vlist[2], poly_vlist[5]).
    // - Third point is on edge (poly_vlist[5], poly_vlist[7]).
    IJK::ARRAY<CTYPE> vcoordX(3*dimension);

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = false;

    CTYPE * vcoordX02 = vcoordX.Ptr();
    CTYPE * vcoordX25 = vcoordX.Ptr()+dimension;
    CTYPE * vcoordX57 = vcoordX.Ptr()+2*dimension;

    IJK::compute_midpoint
      (dimension, poly_vlist_coord[0], poly_vlist_coord[2], vcoordX02);
    IJK::compute_midpoint
      (dimension, poly_vlist_coord[2], poly_vlist_coord[5], vcoordX25);
    IJK::compute_midpoint
      (dimension, poly_vlist_coord[5], poly_vlist_coord[7], vcoordX57);

    IJK::compute_coord_centroid
      (dimension, quadA_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadA_centroid.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadC_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadC_centroid.Ptr());

    compute_cos_max_min_four_polyT_QTQT_angle
      (dimension, vertex_coord, vcoordX.PtrConst(), 
       quadA_centroid.PtrConst(), quadC_centroid.PtrConst(), poly_vlist, 
       max_small_magnitude, poly_tri_result, cos_min_angle, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  poly vert: ";
    IJK::print_list(cerr, poly_vlist, VLIST_LENGTH);
    cerr << endl;
    cerr << "  poly_tri_result[0].num_triangles: "
         << poly_tri_result[0].num_triangles << endl;
    cerr << "  poly_tri_result[1].num_triangles: "
         << poly_tri_result[1].num_triangles << endl;
    cerr << "  poly_tri_result[2].num_triangles: "
         << poly_tri_result[2].num_triangles << endl;
    cerr << "  poly_tri_result[3].num_triangles: "
         << poly_tri_result[3].num_triangles << endl;
    cerr << "cos_min_angle: " << cos_min_angle << endl;
    */

    if (poly_tri_result[0].num_split_edges == 3) {

      triangulate_four_polyT_QTQT_max_min_angle
        (dimension, poly_vlist, flag_reverse_orient, 
         poly_tri_result, vcoordX.PtrConst(), 
         quadA_centroid.PtrConst(), quadC_centroid.PtrConst(),
         vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = true;
    }
    else {
      flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = false;
    }
  }

  ///@}


  // **************************************************
  /// @name TRIANGULATE QUAD-TRIANGLE-QUAD
  // **************************************************


  ///@{

  /*!
   *  Triangulate quad-triangle-quad.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  @param poly_vert[] Seven vertices shared by the triangle and two quads.
   *  - Quad A has vertices 
   *      (poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6])
   *  - Triangle B has vertices (poly_vert[1], poly_vert[4], poly_vert[5])
   *  - Quad C has vertices 
   *      (poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4])
   *  @param ivertX15 Index of vertex splitting (poly_vert[1], poly_vert[5]).
   *  @param ivertX14 Index of vertex splitting (poly_vert[1], poly_vert[4]).
   *   - Note: ivertX15 precedes ivertX14 in the argument list because
   *     ivertX15 is on the boundary of quadA.
   *  @param[out] is_triangleB_split If true, triangleB is split into
   *     new triangles which are added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, 
            typename VTYPE2, typename VTYPE3,
            typename MTYPE>
  void triangulate_quad_triangle_quad_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * poly_vert,
   const VTYPE1 ivertX15,
   const VTYPE2 ivertX14,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE3> & tri_vert,
   bool & is_triangleB_split)
  {
    const int NUM_VERT_PER_TRIANGLE(3);
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vcoord0 = vertex_coord+poly_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+poly_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+poly_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord+poly_vert[4]*dimension;
    const CTYPE * vcoord5 = vertex_coord+poly_vert[5]*dimension;
    const CTYPE * vcoord6 = vertex_coord+poly_vert[6]*dimension;
    const CTYPE * vcoordX14 = vertex_coord+ivertX14*dimension;
    const CTYPE * vcoordX15 = vertex_coord+ivertX15*dimension;
    const VTYPE0 quadA_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6] };
    const VTYPE0 triB_vert[NUM_VERT_PER_TRIANGLE] = 
      { poly_vert[1], poly_vert[4], poly_vert[5] };
    const VTYPE0 quadC_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4] };
    CTYPE cos_min_quadA, cos_min_triB, cos_min_quadC;
    CTYPE cos_min_pentagonAB, cos_min_pentagonBC;
    CTYPE cos_min_no_edge_split;
    CTYPE cos_min_pentagonAB_quadC, cos_min_pentagonBC_quadA;
    CTYPE cos_min_quadX_upper, cos_min_quadX;
    CTYPE cos_pentagon_triX14, cos_pentagon_triX15;
    bool flag_diag02_quadA, flag_diag02_quadC;
    bool flag_diag02_quadX_upper;
    bool flag_zero_quadA, flag_zero_triB, flag_zero_quadC;
    bool flag_zero_no_edge_split;
    bool flag_zero_pentagonAB, flag_zero_pentagonBC;
    bool flag_zero_pentagonAB_quadC, flag_zero_pentagonBC_quadA;
    bool flag_pentagon_quadX_upper;
    bool flag_zero_pentagon_triX14, flag_zero_pentagon_triX15;
    bool flag_zero_quadX_upper, flag_zero_quadX;

    // Initialize
    is_triangleB_split = false;

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord0, vcoord1, vcoord5, vcoord6, max_small_magnitude, 
       cos_min_quadA, flag_diag02_quadA, flag_zero_quadA);
    compute_cos_min_triangle_angle
      (dimension, vcoord1, vcoord4, vcoord5, max_small_magnitude, 
       cos_min_triB, flag_zero_triB);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord1, vcoord2, vcoord3, vcoord4, max_small_magnitude, 
       cos_min_quadC, flag_diag02_quadC, flag_zero_quadC);
    cos_min_no_edge_split = std::max(cos_min_quadA, cos_min_quadC);
    cos_min_no_edge_split = std::max(cos_min_no_edge_split, cos_min_triB);
    flag_zero_no_edge_split = flag_zero_quadA || flag_zero_triB || 
      flag_zero_quadC;

    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord4, vcoord5, vcoord6, vcoordX15, 
       max_small_magnitude, cos_min_pentagonAB, flag_zero_pentagonAB);
    cos_min_pentagonAB_quadC = cos_min_pentagonAB;
    if (cos_min_pentagonAB_quadC < cos_min_quadC)
      { cos_min_pentagonAB_quadC = cos_min_quadC; }
    flag_zero_pentagonAB_quadC = flag_zero_pentagonAB || flag_zero_quadC;

    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord1, vcoord2, vcoord3, vcoord4, vcoord5, vcoordX14, 
       max_small_magnitude, cos_min_pentagonBC, flag_zero_pentagonBC);
    cos_min_pentagonBC_quadA = cos_min_pentagonBC;
    if (cos_min_pentagonBC_quadA < cos_min_quadA)
      { cos_min_pentagonBC_quadA = cos_min_quadA; }
    flag_zero_pentagonBC_quadA = flag_zero_pentagonBC || flag_zero_quadA;

    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoordX15, vcoordX14, vcoord4, vcoord5, max_small_magnitude, 
       cos_min_quadX_upper, flag_diag02_quadX_upper, 
       flag_zero_quadX_upper);

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vcoordX15, vcoord5, vcoord6, vcoord0, vcoord1, 
       max_small_magnitude, cos_pentagon_triX15, flag_zero_pentagon_triX15);

    IJK_DEPRECATED::compute_cos_min_pentagon_triangulation_angle
      (dimension, vcoordX14, vcoord1, vcoord2, vcoord3, vcoord4, 
       max_small_magnitude, cos_pentagon_triX14, flag_zero_pentagon_triX14);
 
    cos_min_quadX = cos_min_quadX_upper;
    if (cos_min_quadX < cos_pentagon_triX15)
      { cos_min_quadX = cos_pentagon_triX15; }
    if (cos_min_quadX < cos_pentagon_triX14)
      { cos_min_quadX = cos_pentagon_triX14; }
    flag_zero_quadX = flag_zero_quadX_upper|| flag_zero_pentagon_triX15 || 
      flag_zero_pentagon_triX14;

    // Default is do not split edges (1,4) or (1,5).
    bool flag_tri_no_edge_split(false);
    bool flag_tri_pentagonAB_quadC(false), flag_tri_pentagonBC_quadA(false);
    bool flag_tri_quadX(false);

    if (flag_zero_pentagonAB_quadC || flag_zero_pentagonBC_quadA ||
        flag_zero_quadX) { flag_tri_no_edge_split = true; }
    else {
      if (cos_min_quadX <= cos_min_pentagonAB_quadC &&
          cos_min_quadX <= cos_min_pentagonBC_quadA)
        { flag_tri_quadX = true; }
      else if (cos_min_pentagonAB_quadC < cos_min_pentagonBC_quadA)
        { flag_tri_pentagonAB_quadC = true;}
      else
        { flag_tri_pentagonBC_quadA = true;}

      if (!flag_zero_no_edge_split) {
        if (cos_min_no_edge_split <= cos_min_quadX &&
            cos_min_no_edge_split <= cos_min_pentagonAB_quadC &&
            cos_min_no_edge_split <= cos_min_pentagonBC_quadA) {
          flag_tri_quadX = false;
          flag_tri_pentagonAB_quadC = false;
          flag_tri_pentagonBC_quadA = false;
          flag_tri_no_edge_split = true;
        }
      }
    }

    if (flag_tri_no_edge_split) {
      // Triangulate quadA and quadC separately using quad diagonals.
      triangulate_quad_using_diagonal
        (quadA_vert, flag_diag02_quadA, flag_reverse_orient, tri_vert);
      triangulate_quad_using_diagonal
        (quadC_vert, flag_diag02_quadC, flag_reverse_orient, tri_vert);
    }
    else if (flag_tri_quadX) {
      const VTYPE0 quadX_upper_vert[NUM_VERT_PER_QUAD] = 
        { ivertX15, ivertX14, poly_vert[4], poly_vert[5] };

      // Triangulate quadX_upper
      triangulate_quad_using_diagonal
        (quadX_upper_vert, flag_diag02_quadX_upper, flag_reverse_orient, 
         tri_vert);

      // Add triangle below quadX_upper
      add_triangle_vertices
        (poly_vert[1], ivertX14, ivertX15, flag_reverse_orient, tri_vert);

      // Triangulate quadA with additional vertex at ivertX16.
      IJK::triangulate_pentagon
        (ivertX15, poly_vert[5], poly_vert[6], poly_vert[0], poly_vert[1],
         flag_reverse_orient, tri_vert);

      // Triangulate quadC with additional vertex at ivertX14.
      IJK::triangulate_pentagon
        (ivertX14, poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4],
         flag_reverse_orient, tri_vert);

      is_triangleB_split = true;
    }
    else if (flag_tri_pentagonAB_quadC) {

      // Triangulate pentagon AB with additional vertex at ivertX15.
      IJK::triangulate_pentagon_with_vertex
        (poly_vert[0], poly_vert[1], poly_vert[4], poly_vert[5], poly_vert[6],
         ivertX15, flag_reverse_orient, tri_vert);

      // Triangulate quadC using quad diagonals.
      triangulate_quad_using_diagonal
        (quadC_vert, flag_diag02_quadC, flag_reverse_orient, tri_vert);

      is_triangleB_split = true;
    }
    else {
      // Use triangulation pentagonBC quadA

      // Triangulate quadA using quad diagonals.
      triangulate_quad_using_diagonal
        (quadA_vert, flag_diag02_quadA, flag_reverse_orient, tri_vert);

      // Triangulate pentagonBC with additional vertex at ivertX14.
      IJK::triangulate_pentagon_with_vertex
        (poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4], poly_vert[5],
         flag_reverse_orient, tri_vert);

      is_triangleB_split = true;
    }
    
  }


  /// Triangulate quad, triangle, quad.
  /// - Maximize minimum triangle angle.
  /// - C++ STL vector format for vertex_coord.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, 
            typename VTYPE2, typename VTYPE3,
            typename MTYPE>
  void triangulate_quad_triangle_quad_max_min_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE0 * poly_vert,
   const VTYPE1 ivertX15,
   const VTYPE2 ivertX14,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE3> & tri_vert,
   bool & is_triangleB_split)
  {
    triangulate_quad_triangle_quad_max_min_angle
      (dimension, IJK::vector2pointer(vertex_coord), poly_vert,
       ivertX15, ivertX14, flag_reverse_orient, max_small_magnitude, 
       tri_vert, is_triangleB_split);
  }


  /*!
   *  Triangulate quad, triangle, quad.
   *  - Maximize minimum triangle angle.
   *  - Allow replacement of quad, triangle, quad with quad, pentagon, quad.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  @param poly_vert[] Seven vertices shared by the two quadrilaterals 
   *      and triangle.
   *  - Quad A has vertices 
   *      (poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6])
   *  - Triangle B has vertices 
   *      (poly_vert[1], poly_vert[2], poly_vert[5])
   *  - Quad C has vertices 
   *      (poly_vert[2], poly_vert[3], poly_vert[4], poly_vert[5])
   *  @param[out] flag_replace If true, then all polygons are replaced.
  */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_triangle_quad_max_min_angle_allow_2pent
  (const DTYPE dimension,
   const VTYPE0 * poly_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert,
   bool & flag_replace)
  {
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord0 = vertex_coord_ptr+poly_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+poly_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord_ptr+poly_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    const CTYPE * vcoord5 = vertex_coord_ptr+poly_vert[5]*dimension;
    const CTYPE * vcoord6 = vertex_coord_ptr+poly_vert[6]*dimension;
    IJK::ARRAY<CTYPE> coordB(dimension);
    IJK::ARRAY<CTYPE> pentagonA_centroid(dimension);
    IJK::ARRAY<CTYPE> pentagonC_centroid(dimension);
    IJK::POLY_TRIANGULATION_RESULT<16,CTYPE,int> poly_tri[3];
    CTYPE cos_max_min_angle;
    bool flag_zero;

    IJK::compute_midpoint_of_midpoints
      (dimension, vcoord1, vcoord5, vcoord2, vcoord5, coordB.Ptr());
    IJK::compute_pentagon_centroid
      (dimension, vcoord0, vcoord1, coordB.PtrConst(), vcoord5, vcoord6,
       pentagonA_centroid.Ptr());
    IJK::compute_pentagon_centroid
      (dimension, vcoord2, vcoord3, vcoord4, vcoord5, coordB.PtrConst(),
       pentagonC_centroid.Ptr());

    compute_cos_max_min_quad_triangle_quad_angle_2pent
      (dimension, vertex_coord, pentagonA_centroid.PtrConst(),
       coordB.PtrConst(), pentagonC_centroid.PtrConst(), poly_vert,
       max_small_magnitude, cos_max_min_angle, poly_tri, flag_zero);


    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  poly_tri[0].num_triangles: " << poly_tri[0].num_triangles << endl;
    cerr << "  poly_tri[1].num_triangles: " << poly_tri[1].num_triangles << endl;
    cerr << "  poly_tri[2].num_triangles: " << poly_tri[2].num_triangles << endl;
    IJK::print_coord3D(cerr, "coordB: ", coordB.PtrConst(), "\n");
    */

    // *** SHOULD REPLACE WITH num_split_edges.  ***
    if (poly_tri[0].num_triangles == 2 || poly_tri[2].num_triangles == 2) {
      flag_replace = false;
    }
    else {
      flag_replace = true;

      const VTYPE0 isovB = 
        IJK::insert_coord(dimension, coordB.PtrConst(), vertex_coord);

      // Add triangle.
      add_triangle_vertices
        (poly_vert[1], poly_vert[2], isovB, flag_reverse_orient, tri_vert);

      IJK_DEPRECATED::triangulate_pentagon_tri3_or_tri5
        (dimension, poly_vert[0], poly_vert[1], isovB, poly_vert[5],
         poly_vert[6], pentagonA_centroid.PtrConst(), flag_reverse_orient,
         poly_tri[0], vertex_coord, tri_vert);

      IJK_DEPRECATED::triangulate_pentagon_tri3_or_tri5
        (dimension, poly_vert[2], poly_vert[3], poly_vert[4], poly_vert[5], 
         isovB, pentagonC_centroid.PtrConst(), flag_reverse_orient,
         poly_tri[2], vertex_coord, tri_vert);
    }
  }


  /*!
   *  Triangulate quad, triangle, quad.
   *  - Maximize minimum triangle angle.
   *  - Allow replacement of quad, triangle, quad with quad, pentagon, quad.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  @param poly_vert[] Seven vertices shared by the two quadrilaterals 
   *      and triangle.
   *  - Quad A has vertices 
   *      (poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6])
   *  - Triangle B has vertices 
   *      (poly_vert[1], poly_vert[2], poly_vert[5])
   *  - Quad C has vertices 
   *      (poly_vert[2], poly_vert[3], poly_vert[4], poly_vert[5])
  */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_triangle_quad_max_min_angle_allow_2pent_and_sept
  (const DTYPE dimension,
   const VTYPE0 * poly_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_PENTAGON(5);
    const int NUM_VERT_PER_SEPTAGON(7);
    const int NUM_POLY(3);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord0 = vertex_coord_ptr+poly_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+poly_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord_ptr+poly_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    const CTYPE * vcoord5 = vertex_coord_ptr+poly_vert[5]*dimension;
    const CTYPE * vcoord6 = vertex_coord_ptr+poly_vert[6]*dimension;
    const VTYPE0 quadA_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6] };
    const VTYPE0 triangleB_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[1], poly_vert[2], poly_vert[5] };
    const VTYPE0 quadC_vert[NUM_VERT_PER_QUAD] = 
            { poly_vert[2], poly_vert[3], poly_vert[4], poly_vert[5] };
    IJK::ARRAY<CTYPE> coordB(dimension);
    IJK::ARRAY<CTYPE> pentagonA_centroid(dimension);
    IJK::ARRAY<CTYPE> pentagonC_centroid(dimension);
    CTYPE cos_min_ABC_angle, cos_min_septA_angle, cos_min_septC_angle;
    VTYPE0 ivA_max_min, ivC_max_min;
    bool flagA_tri5, flagC_tri5, flag_zero_ABC;
    bool flag_zero_septA, flag_zero_septC;
    bool flag_tri_septA, flag_tri_septC;

    IJK::compute_midpoint_of_midpoints
      (dimension, vcoord1, vcoord5, vcoord2, vcoord5, coordB.Ptr());
    IJK::compute_pentagon_centroid
      (dimension, vcoord0, vcoord1, coordB.PtrConst(), vcoord5, vcoord6,
       pentagonA_centroid.Ptr());
    IJK::compute_pentagon_centroid
      (dimension, vcoord2, vcoord3, vcoord4, vcoord5, coordB.PtrConst(),
       pentagonC_centroid.Ptr());

    compute_cos_max_min_pentagon_triangle_pentagon_angle
      (dimension, vcoord0, vcoord1, coordB.PtrConst(), vcoord2, vcoord3, 
       vcoord4, vcoord5, vcoord6, 
       pentagonA_centroid.PtrConst(), pentagonC_centroid.PtrConst(),
       max_small_magnitude, cos_min_ABC_angle,
       flagA_tri5, flagC_tri5, ivA_max_min, ivC_max_min, flag_zero_ABC);

    // Compute cos min angle for triangulation from vertex v0.
    compute_cos_min_sept_triangulation_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoord5,
       vcoord6, max_small_magnitude, cos_min_septA_angle, flag_zero_septA);

    // Compute cos min angle for triangulation from vertex v3.
    compute_cos_min_sept_triangulation_angle
      (dimension, vcoord3, vcoord4, vcoord5, vcoord6, vcoord0, vcoord1,
       vcoord2, max_small_magnitude, cos_min_septC_angle, flag_zero_septC);

    int index_selected;
    bool flag_zero;
    CTYPE cos_min_angle;
    select_minIII
      (cos_min_ABC_angle, flag_zero_ABC,
       cos_min_septA_angle, flag_zero_septA,
       cos_min_septC_angle, flag_zero_septC,
       cos_min_angle, index_selected, flag_zero);

    if ((poly_vert[0] == poly_vert[3]) || (poly_vert[4] == poly_vert[6])) {
      // Triangulating from poly_vert[0] or poly_vert[3] could create non-manifold.
      // Allow only quad-triangle-quad triangulation.
      index_selected = 0;
    }

    if (index_selected == 0 || flag_zero) {

      const VTYPE0 isovB = 
        IJK::insert_coord(dimension, coordB.PtrConst(), vertex_coord);

      // Add triangle.
      add_triangle_vertices
        (poly_vert[1], poly_vert[2], isovB, flag_reverse_orient, tri_vert);

      const VTYPE0 pentagonA_vert[NUM_VERT_PER_PENTAGON] =
        { poly_vert[0], poly_vert[1], isovB, poly_vert[5], poly_vert[6] };

      const VTYPE0 pentagonC_vert[NUM_VERT_PER_PENTAGON] =
        { isovB, poly_vert[2], poly_vert[3], poly_vert[4], poly_vert[5] };

      IJK_DEPRECATED::triangulate_pentagon_tri3_or_tri5
        (dimension, poly_vert[0], poly_vert[1], isovB, poly_vert[5],
         poly_vert[6], pentagonA_centroid.PtrConst(), flagA_tri5,
         ivA_max_min, flag_reverse_orient, vertex_coord, tri_vert);

      IJK_DEPRECATED::triangulate_pentagon_tri3_or_tri5
        (dimension, isovB, poly_vert[2], poly_vert[3], poly_vert[4],
         poly_vert[5], pentagonC_centroid.PtrConst(), flagC_tri5,
         ivC_max_min, flag_reverse_orient, vertex_coord, tri_vert);
    }
    else if (index_selected == 1) {
      triangulate_septagon(poly_vert, 0, flag_reverse_orient, tri_vert);
    }
    else {
      triangulate_septagon(poly_vert, 3, flag_reverse_orient, tri_vert);
    }

  }


  ///@}


  // **************************************************
  /// @name TRIANGULATE TRIANGLE-QUAD
  // **************************************************

  ///@{

  /// Triangulate triangle-quad.
  /// - Add new triangles to vector tri_vert.
  /// - Add new coords to vertex_coord[].
  /// - Triangle has vertices:
  ///     (triangle_vert, lower_quad_vert[0], upper_quad_vert[1]).
  /// - Quad has vertices 
  ///     (lower_quad_vert[0], lower_quad_vert[1], 
  ///      upper_quad_vert[1], upper_quad_vert[0])
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename VTYPET, typename TRI_VTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_triangle_quad_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const VTYPET triangle_vert,
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result,
   const CTYPEX * vcoordX00,
   const CTYPEX * quad_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const TRI_VTYPE ivertX00 = 
      IJK::insert_coord(dimension, vcoordX00, vertex_coord);

    // Split triangleA into two triangles.
    add_triangle_vertices
      (triangle_vert, lower_quad_vert[0], ivertX00,
       flag_reverse_orient, tri_vert);
    add_triangle_vertices
      (upper_quad_vert[0], triangle_vert, ivertX00,
       flag_reverse_orient, tri_vert);

    // Triangulate quad with additional vertex at ivertX00.
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, lower_quad_vert[0], lower_quad_vert[1], 
       upper_quad_vert[1], upper_quad_vert[0], ivertX00, quad_vcoord, 
       flag_reverse_orient, quad_tri_result, vertex_coord, tri_vert);
  }


  /// Triangulate triangle-quad.
  /// - Maximize minimum triangle angle.
  /// - Add new triangles to vector tri_vert.
  /// - Triangle A or triangulation of triangle A is added to tri_vert.
  /// @param vcoordX14 Coordinates of vertex on edge [v1,v4].
  /// @param poly_vert[] Five vertices shared by the triangle and quad.
  /// - Triangle A has vertices   
  ///     (poly_vert[0], poly_vert[1], poly_vert[4])
  /// - Quad B has vertices 
  ///     (poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4])
  /// @param[out] flag_split True if triangle and quad are split into
  ///     triangles which are added to tri_vert.
  template <typename DTYPE, typename CTYPEX, typename CTYPE,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE>
  void triangulate_triangle_quad_max_min_angle
  (const DTYPE dimension,
   const CTYPEX * vcoordX14,
   const VTYPE0 * poly_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert,
   bool & flag_split)
  {
    const int NUM_VERT_PER_TRIANGLE(3);
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord0 = vertex_coord_ptr+poly_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+poly_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord_ptr+poly_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    CTYPE cos_min_triA, cos_min_quadB, cos_min_pentagon;
    CTYPE cos_min_polyAB;
    bool flag_diag02_quadB;
    bool flag_zero_triA, flag_zero_quadB, flag_zero_pentagon;

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoord4, max_small_magnitude, 
       cos_min_triA, flag_zero_triA);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord1, vcoord2, vcoord3, vcoord4, max_small_magnitude, 
       cos_min_quadB, flag_diag02_quadB, flag_zero_quadB);
    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoordX14, 
       max_small_magnitude, cos_min_pentagon, flag_zero_pentagon);

    cos_min_polyAB = std::max(cos_min_triA, cos_min_quadB);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "cos_min_pentagon: " << cos_min_pentagon << endl;
    cerr << "cos_min_polyAB: " << cos_min_polyAB << endl;
    */

    if (!flag_zero_pentagon && (cos_min_pentagon < cos_min_polyAB)) {

      // Triangulate pentagon with additional vertex at vcoordX14.
      IJK_DEPRECATED::triangulate_pentagon_tri5
        (dimension, poly_vert, vcoordX14, flag_reverse_orient, 
         vertex_coord, tri_vert);
      flag_split = true;
    }
    else {
      flag_split = false;
    }
  }


  /// Triangulate triangle-quad.
  /// - Maximize minimum triangle angle.
  /// - Add new triangles to vector tri_vert.
  /// - Triangle A or triangulation of triangle A is added to tri_vert.
  /// @param poly_vert[] Five vertices shared by the triangle and quad.
  /// - Triangle A has vertices   
  ///     (poly_vert[0], poly_vert[1], poly_vert[4])
  /// - Quad B has vertices 
  ///     (poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4])
  /// @param[out] flag_split True if triangle and quad are split into
  ///     triangles which are added to tri_vert.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE>
  void triangulate_triangle_quad_max_min_angle
  (const DTYPE dimension,
   const VTYPE0 * poly_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert,
   bool & flag_split)
  {
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    IJK::ARRAY<CTYPE> vcoordX14(dimension);

    IJK::compute_midpoint(dimension, vcoord1, vcoord4, vcoordX14.Ptr());

    triangulate_triangle_quad_max_min_angle
      (dimension, vcoordX14.Ptr(), poly_vert, flag_reverse_orient, 
       max_small_magnitude, vertex_coord, tri_vert, flag_split);
  }


  /// Triangulate triangle-quad.
  /// - Maximize minimum triangle angle.
  /// - Add new triangles to vector tri_vert.
  /// - Triangle A or triangulation of triangle A is added to tri_vert.
  /// @param vcoordX14 Coordinates of vertex on edge [v1,v4].
  /// @param quad_centroid Centroid of quad B.
  /// @param poly_vert[] Five vertices shared by the triangle and quad.
  /// - Triangle A has vertices   
  ///     (poly_vert[0], poly_vert[1], poly_vert[4])
  /// - Quad B has vertices 
  ///     (poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4])
  /// @param[out] flag_split True if triangle and quad are split into
  ///     triangles which are added to tri_vert.
  template <typename DTYPE, typename CTYPEX, typename CTYPE,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE>
  void triangulate_triangle_quad_max_min_angle_allow_tri4
  (const DTYPE dimension,
   const CTYPEX * vcoordX14,
   const CTYPEX * quad_centroid,
   const VTYPE0 * poly_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert,
   bool & flag_split)
  {
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord0 = vertex_coord_ptr+poly_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+poly_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord_ptr+poly_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    CTYPE cos_min_triA, cos_min_quadB, cos_min_pentagon;
    CTYPE cos_min_polyAB;
    bool flag_tri4, flag_diag02_quadB;
    bool flag_zero_triA, flag_zero_quadB, flag_zero_pentagon;

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoord4, max_small_magnitude, 
       cos_min_triA, flag_zero_triA);
    IJK_DEPRECATED::compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord1, vcoord2, vcoord3, vcoord4, quad_centroid,
       max_small_magnitude, cos_min_quadB, flag_tri4, 
       flag_diag02_quadB, flag_zero_quadB);
    IJK_DEPRECATED::compute_cos_min_pentagon_tri5_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoordX14, 
       max_small_magnitude, cos_min_pentagon, flag_zero_pentagon);

    cos_min_polyAB = std::max(cos_min_triA, cos_min_quadB);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "cos_min_pentagon: " << cos_min_pentagon << endl;
    cerr << "cos_min_triA: " << cos_min_triA << endl;
    cerr << "cos_min_quadB: " << cos_min_quadB << endl;
    cerr << "cos_min_polyAB: " << cos_min_polyAB << endl;
    */

    if (!flag_zero_pentagon && (cos_min_pentagon < cos_min_polyAB)) {

      // Triangulate pentagon with additional vertex at vcoordX14.
      IJK_DEPRECATED::triangulate_pentagon_tri5
        (dimension, poly_vert, vcoordX14, flag_reverse_orient, 
         vertex_coord, tri_vert);
      flag_split = true;
    }
    else {
      flag_split = false;
    }
  }


  /// Triangulate triangle-quad.
  /// - Maximize minimum triangle angle.
  /// - Add new triangles to vector tri_vert.
  /// - Triangle A or triangulation of triangle A is added to tri_vert.
  /// - Allow quad to be split into 4 triangles.
  /// @param poly_vert[] Five vertices shared by the triangle and quad.
  /// - Triangle A has vertices   
  ///     (poly_vert[0], poly_vert[1], poly_vert[4])
  /// - Quad B has vertices 
  ///     (poly_vert[1], poly_vert[2], poly_vert[3], poly_vert[4])
  /// @param[out] flag_split True if triangle and quad are split into
  ///     triangles which are added to tri_vert.
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE>
  void triangulate_triangle_quad_max_min_angle_allow_tri4
  (const DTYPE dimension,
   const VTYPE0 * poly_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert,
   bool & flag_split)
  {
    const int NUM_VERT_PER_QUAD(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    const VTYPE0 * quad_vert = poly_vert+1;
    IJK::ARRAY<CTYPE> vcoordX14(dimension);
    IJK::ARRAY<CTYPE> quad_centroid(dimension);

    IJK::compute_midpoint(dimension, vcoord1, vcoord4, vcoordX14.Ptr());
    IJK::compute_coord_centroid
      (dimension, quad_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quad_centroid.Ptr());

    triangulate_triangle_quad_max_min_angle_allow_tri4
      (dimension, vcoordX14.Ptr(), quad_centroid.Ptr(), poly_vert,
       flag_reverse_orient,  max_small_magnitude, vertex_coord, tri_vert, 
       flag_split);
  }

  ///@}


  // ****************************************************************
  /// @name TRIANGULATE QUAD-QUAD-TRIANGLE or TRIANGLE-QUAD-QUAD
  // ****************************************************************

  ///@{

  /// Triangulate triangle-quad-quad.
  /// - Break middle quad into two quads and triangulate each.
  /// - Add new triangles to vector tri_vert.
  /// - Add new coords to vertex_coord[].
  /// - triangle A has vertices 
  ///     (lower_quad_vert[0], upper_quad_vert[1], triangle_vert)
  /// - Quad B has vertices 
  ///     (lower_quad_vert[1], lower_quad_vert[2], 
  ///      upper_quad_vert[2], upper_quad_vert[1])
  /// - Quad C has vertices 
  ///     (lower_quad_vert[2], lower_quad_vert[3], 
  ///      upper_quad_vert[3], upper_quad_vert[2])
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename VTYPET, typename TRI_VTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_triangle_quad_quad_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const VTYPET triangle_vert,
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> quad_tri_result[2],
   const CTYPEX * vcoordX00,
   const CTYPEX * vcoordX11,
   const CTYPEX * quadC_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const TRI_VTYPE * lower_quadB_vert = lower_quad_vert;
    const TRI_VTYPE * upper_quadB_vert = upper_quad_vert;

    const TRI_VTYPE ivertX00 = 
      IJK::insert_coord(dimension, vcoordX00, vertex_coord);
    const TRI_VTYPE ivertX11 = 
      IJK::insert_coord(dimension, vcoordX11, vertex_coord);

    // Triangulate subquads of quad B.
    triangulate_subquads_using_diagonals_LU
      (lower_quadB_vert, upper_quadB_vert, ivertX00, ivertX11,
       quad_tri_result[0], flag_reverse_orient, tri_vert);

    // Split triangle A into two triangles.
    add_triangle_vertices
      (triangle_vert, lower_quad_vert[0], ivertX00,
       flag_reverse_orient, tri_vert);
    add_triangle_vertices
      (upper_quad_vert[0], triangle_vert, ivertX00,
       flag_reverse_orient, tri_vert);

    // Triangulate quadC with additional vertex at ivertX11.
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, lower_quad_vert[1], lower_quad_vert[2], 
       upper_quad_vert[2], upper_quad_vert[1], ivertX11, quadC_vcoord, 
       flag_reverse_orient, quad_tri_result[1], vertex_coord, tri_vert);
  }


  /// Triangulate quad-quad-triangle, allow tri5.
  /// - Break middle quad into two quads and triangulate each.
  /// - Add new triangles to vector tri_vert.
  /// - Add new coords to vertex_coord[].
  /// @param poly_vert[] Seven vertices shared by the two quads and triangle.
  /// - Quad A has vertices 
  ///     (poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6])
  /// - Quad B has vertices 
  ///     (poly_vert[1], poly_vert[2], poly_vert[4], poly_vert[5])
  /// - Triangle C has vertices 
  ///     (poly_vert[2], poly_vert[3], poly_vert[4])
  template <typename DTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE,
            typename CTYPE, typename CTYPEX, 
            typename VTYPE0, typename VTYPE1>
  void triangulate_quad_quad_triangle_allow_tri5
  (const DTYPE dimension,
   const VTYPE0 * poly_vert,
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> poly_tri[3],
   const CTYPEX * vcoordX15,
   const CTYPEX * vcoordX24,
   const CTYPEX * quadA_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert)
  {
    const int NUM_VERT_PER_QUAD(4);
    const VTYPE0 quadB_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[1], poly_vert[2], poly_vert[4], poly_vert[5] };

    const VTYPE0 ivertX15 = 
      IJK::insert_coord(dimension, vcoordX15, vertex_coord);
    const VTYPE0 ivertX24 = 
      IJK::insert_coord(dimension, vcoordX24, vertex_coord);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    IJK::print_quad_vertices(cerr, "  Quad B vert: ", quadB_vert, "\n");
    */

    triangulate_subquads_using_diagonals_vX03_vX12
      (quadB_vert, ivertX15, ivertX24, poly_tri[1], flag_reverse_orient,
       tri_vert);

    // Triangulate quadA with additional vertex at ivertX15.
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX15, poly_vert[5], poly_vert[6], poly_vert[0], 
       poly_vert[1], quadA_vcoord, poly_tri[0],
       flag_reverse_orient, vertex_coord, tri_vert);

    // Split triangleC into two triangles.
    add_triangle_vertices
      (ivertX24, poly_vert[2], poly_vert[3], flag_reverse_orient, tri_vert);
    add_triangle_vertices
      (ivertX24, poly_vert[3], poly_vert[4], flag_reverse_orient, tri_vert);
  }


  /// Triangulate quad-quad-triangle.
  /// - Maximize minimum triangle angle.
  /// - Add new triangles to vector tri_vert.
  /// - Add new coords, if necessary, to vertex_coord[].
  /// - Consider triangulation of quad A into 5 triangles.
  /// - Only triangulate quads where merging with other quads 
  ///   improves triangulation.
  /// @param poly_vert[] Seven vertices shared by the two quads and triangle.
  /// - Quad A has vertices 
  ///     (poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6])
  /// - Quad B has vertices 
  ///     (poly_vert[1], poly_vert[2], poly_vert[4], poly_vert[5])
  /// - Triangle C has vertices 
  ///     (poly_vert[2], poly_vert[3], poly_vert[4])
  /// @param[out] flag_split[] flag_split[i] is true if quad i is
  ///     is split into triangles which have been added to tri_vert[].
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_quad_triangle_max_min_angle_allow_tri5
  (const DTYPE dimension,
   const VTYPE0 * poly_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<VTYPE1> & tri_vert,
   bool flag_split[3])
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_POLY(3);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * vcoord1 = vertex_coord_ptr+poly_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+poly_vert[2]*dimension;
    const CTYPE * vcoord4 = vertex_coord_ptr+poly_vert[4]*dimension;
    const CTYPE * vcoord5 = vertex_coord_ptr+poly_vert[5]*dimension;
    const VTYPE0 quadA_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[0], poly_vert[1], poly_vert[5], poly_vert[6] };
    const VTYPE0 quadB_vert[NUM_VERT_PER_QUAD] = 
      { poly_vert[1], poly_vert[2], poly_vert[4], poly_vert[5] };
    IJK::ARRAY<CTYPE> vcoordX15(dimension);
    IJK::ARRAY<CTYPE> vcoordX24(dimension);
    IJK::ARRAY<CTYPE> quadA_centroid(dimension);
    IJK::ARRAY<CTYPE> quadB_centroid(dimension);
    POLY_TRIANGULATION_RESULTX2<16,CTYPE,int> poly_tri[NUM_POLY];
    CTYPE cos_min_angle;
    bool flag_zero;

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = false;

    IJK::compute_midpoint(dimension, vcoord1, vcoord5, vcoordX15.Ptr());
    IJK::compute_midpoint(dimension, vcoord2, vcoord4, vcoordX24.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadA_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadA_centroid.Ptr());
    IJK::compute_coord_centroid
      (dimension, quadB_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadB_centroid.Ptr());

    compute_cos_max_min_quad_quad_triangle_angle_allow_tri5
      (dimension, vertex_coord, vcoordX15.PtrConst(), vcoordX24.PtrConst(),
       quadA_centroid.PtrConst(), quadB_centroid.PtrConst(),
       poly_vert, max_small_magnitude, poly_tri, 
       cos_min_angle, flag_zero);


    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "poly_tri[0].num_triangles: " << poly_tri[0].num_triangles << endl;
    cerr << "poly_tri[1].num_triangles: " << poly_tri[1].num_triangles << endl;
    cerr << "poly_tri[2].num_triangles: " << poly_tri[2].num_triangles << endl;
    cerr << "poly_tri[0].num_split_edges: " 
         << poly_tri[0].num_split_edges << endl;
    cerr << "poly_tri[1].num_split_edges: " 
         << poly_tri[1].num_split_edges << endl;
    cerr << "poly_tri[2].num_split_edges: " 
         << poly_tri[2].num_split_edges << endl;
    */

    if (poly_tri[1].num_split_edges == 2) {
      triangulate_quad_quad_triangle_allow_tri5
        (dimension, poly_vert, flag_reverse_orient, poly_tri,
         vcoordX15.PtrConst(), vcoordX24.PtrConst(), 
         quadA_centroid.PtrConst(), vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = true;
    }
    else if (poly_tri[1].num_split_edges == 1) {

      if (poly_tri[2].num_split_edges == 1) {

        // Use triangulation around X25
        const VTYPE0 ivertX24 = IJK::insert_coord
          (dimension, vcoordX24.PtrConst(), vertex_coord);

        // Triangulate quadB with additional vertex at ivertX24.
        IJK::triangulate_pentagon
          (ivertX24, poly_vert[4], poly_vert[5], poly_vert[1], poly_vert[2],
           flag_reverse_orient, tri_vert);

        // Split triangleC into two triangles.
        add_triangle_vertices
          (ivertX24, poly_vert[2], poly_vert[3], 
           flag_reverse_orient, tri_vert);
        add_triangle_vertices
          (ivertX24, poly_vert[3], poly_vert[4], 
           flag_reverse_orient, tri_vert);

        flag_split[1] = flag_split[2] = true;
      }
      else if (poly_tri[0].num_split_edges == 1) {
        const VTYPE0 ivertX15 = 
          IJK::insert_coord(dimension, vcoordX15.PtrConst(), vertex_coord);

        // Triangulate quadA with additional vertex at ivertX15.
        IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
          (dimension, ivertX15, poly_vert[5], poly_vert[6], poly_vert[0], 
           poly_vert[1], quadA_centroid.PtrConst(), poly_tri[0],
           flag_reverse_orient, vertex_coord, tri_vert);

        // Triangulate quadB with additional vertex at ivertX15.
        IJK::triangulate_pentagon
          (ivertX15, poly_vert[1], poly_vert[2], poly_vert[4], poly_vert[5],
           flag_reverse_orient, tri_vert);

        flag_split[0] = flag_split[1] = true;
      }
    }

  }

  ///@}


  // ****************************************************************
  /// @name TRIANGULATE TRIANGLE-QUADx3 and QUADx3-TRIANGLE
  // ****************************************************************

  ///@{

  /*!
   *  Triangulate triangle-quad-quad-quad.
   *  - Break each of the middle two quads into two smaller quads 
   *      and triangulate each.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords to vertex_coord[].
   *  - triangle A has vertices 
   *      (lower_quad_vert[0], upper_quad_vert[1], triangle_vert)
   *  - Quad B has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad C has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad D has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   *  @param vcoordX[] Coordinates of three new points on edges.
   *       (Usually midpoints of those edges.)
   *  - First point is on edge (lower_quad_vert[0], upper_quad_vert[0]).
   *  - Second point is on edge (lower_quad_vert[1], upper_quad_vert[1]).
   *  - Third point is on edge (lower_quad_vert[2], upper_quad_vert[2]).
  */
  template <typename DTYPE, 
            typename VTYPEL, typename VTYPEU, typename VTYPET,
            typename TRI_VTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_triangle_quadx3_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const VTYPET triangle_vert,
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   poly_tri[4],
   const CTYPEX * vcoordX,
   const CTYPEX * quadD_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    TRI_VTYPE ivertX[3];
    const TRI_VTYPE * lower_quadB_vert = lower_quad_vert;
    const TRI_VTYPE * upper_quadB_vert = upper_quad_vert;
    const TRI_VTYPE * lower_quadC_vert = lower_quad_vert+1;
    const TRI_VTYPE * upper_quadC_vert = upper_quad_vert+1;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  lower_quad_vert: ";
    IJK::print_list(cerr, lower_quad_vert, 4);
    cerr << endl;
    cerr << "  upper_quad_vert: ";
    IJK::print_list(cerr, upper_quad_vert, 4);
    cerr << endl;
    */

    ivertX[0] = 
      IJK::insert_coord(dimension, vcoordX, vertex_coord);
    ivertX[1] = 
      IJK::insert_coord(dimension, vcoordX+dimension, vertex_coord);
    ivertX[2] = 
      IJK::insert_coord(dimension, vcoordX+2*dimension, vertex_coord);

    // Triangulate subquads of quad B.
    triangulate_subquads_using_diagonals_LU
      (lower_quadB_vert, upper_quadB_vert, ivertX[0], ivertX[1],
       poly_tri[1], flag_reverse_orient, tri_vert);

    // Triangulate subquads of quad C.
    triangulate_subquads_using_diagonals_LU
      (lower_quadC_vert, upper_quadC_vert, ivertX[1], ivertX[2],
       poly_tri[2], flag_reverse_orient, tri_vert);

    // Split triangle A into two triangles.
    add_triangle_vertices
      (triangle_vert, lower_quad_vert[0], ivertX[0],
       flag_reverse_orient, tri_vert);
    add_triangle_vertices
      (upper_quad_vert[0], triangle_vert, ivertX[0],
       flag_reverse_orient, tri_vert);

    // Triangulate quadD with additional vertex at ivertX[2].
    IJK_DEPRECATED::triangulate_pentagon_v0tri3_or_tri5
      (dimension, ivertX[2], lower_quad_vert[2], lower_quad_vert[3], 
       upper_quad_vert[3], upper_quad_vert[2], quadD_vcoord, poly_tri[3],
       flag_reverse_orient, vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate triangle-quad-quad-quad.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  - Consider triangulation of quads into 5 triangles.
   *  - Only triangulate quads where merging with other quads 
   *    improves triangulation.
   *  - triangle A has vertices 
   *      (lower_quad_vert[0], upper_quad_vert[1], triangle_vert)
   *  - Quad B has vertices 
   *      (lower_quad_vert[0], lower_quad_vert[1], 
   *       upper_quad_vert[1], upper_quad_vert[0])
   *  - Quad C has vertices 
   *      (lower_quad_vert[1], lower_quad_vert[2], 
   *       upper_quad_vert[2], upper_quad_vert[1])
   *  - Quad D has vertices 
   *      (lower_quad_vert[2], lower_quad_vert[3], 
   *       upper_quad_vert[3], upper_quad_vert[2])
   *  @param[out] flag_split[] flag_split[i] is true if poly i is
   *      is split into triangles which have been added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE,
            typename VTYPEL, typename VTYPEU, typename VTYPET,
            typename TRI_VTYPE, typename MTYPE>
  void triangulate_triangle_quadx3_max_min_angle_allow_tri5
  (const DTYPE dimension,
   const VTYPEL lower_quad_vert[],
   const VTYPEU upper_quad_vert[],
   const VTYPET triangle_vert,
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert,
   bool flag_split[4])
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_QUADS(3);
    const int NUM_POLY(NUM_QUADS+1);
    const int NUM_LOWER_QUAD_VERT(NUM_QUADS+1);
    const int NUM_UPPER_QUAD_VERT(NUM_QUADS+1);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_QUAD_VERT] =
      { vertex_coord_ptr+lower_quad_vert[0]*dimension,
        vertex_coord_ptr+lower_quad_vert[1]*dimension,
        vertex_coord_ptr+lower_quad_vert[2]*dimension,
        vertex_coord_ptr+lower_quad_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_QUAD_VERT] =
      { vertex_coord_ptr+upper_quad_vert[0]*dimension,
        vertex_coord_ptr+upper_quad_vert[1]*dimension,
        vertex_coord_ptr+upper_quad_vert[2]*dimension,
        vertex_coord_ptr+upper_quad_vert[3]*dimension };
    const TRI_VTYPE quadB_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[0], lower_quad_vert[1], 
        upper_quad_vert[1], upper_quad_vert[0] };
    const TRI_VTYPE quadC_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[1], lower_quad_vert[2], 
        upper_quad_vert[2], upper_quad_vert[1] };
    const TRI_VTYPE quadD_vert[NUM_VERT_PER_QUAD] = 
      { lower_quad_vert[2], lower_quad_vert[3], 
        upper_quad_vert[3], upper_quad_vert[2] };
    IJK::ARRAY<CTYPE> quad_centroid(NUM_QUADS*dimension);
    CTYPE cos_min_angle;
    POLY_TRIANGULATION_RESULTX2<16,CTYPE,int> poly_tri_result[NUM_POLY];
    bool flag_use_vcoordX11;
    bool flag_zero;

    // Coordinates of edge midpoints.
    // - First point is on edge (lower_quad_vert[0], upper_quad_vert[0]).
    // - Second point is on edge (lower_quad_vert[1], upper_quad_vert[1]).
    // - Third point is on edge (lower_quad_vert[2], upper_quad_vert[2]).
    IJK::ARRAY<CTYPE> vcoordX(3*dimension);

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = false;

    CTYPE * vcoordX00 = vcoordX.Ptr();
    CTYPE * vcoordX11 = vcoordX.Ptr()+dimension;
    CTYPE * vcoordX22 = vcoordX.Ptr()+2*dimension;

    IJK::compute_midpoint
      (dimension, lower_vcoord[0], upper_vcoord[0], vcoordX00);
    IJK::compute_midpoint
      (dimension, lower_vcoord[1], upper_vcoord[1], vcoordX11);
    IJK::compute_midpoint
      (dimension, lower_vcoord[2], upper_vcoord[2], vcoordX22);

    CTYPE * quadB_centroid = quad_centroid.Ptr();
    CTYPE * quadC_centroid = quad_centroid.Ptr()+dimension;
    CTYPE * quadD_centroid = quad_centroid.Ptr()+2*dimension;

    IJK::compute_coord_centroid
      (dimension, quadB_vert, NUM_VERT_PER_QUAD, vertex_coord, quadB_centroid);
    IJK::compute_coord_centroid
      (dimension, quadC_vert, NUM_VERT_PER_QUAD, vertex_coord, quadC_centroid);
    IJK::compute_coord_centroid
      (dimension, quadD_vert, NUM_VERT_PER_QUAD, vertex_coord, quadD_centroid);

    compute_cos_max_min_triangle_quadx3_angle_allow_tri5
    (dimension, vertex_coord, vcoordX.PtrConst(), quad_centroid.PtrConst(),
     lower_quad_vert, upper_quad_vert, triangle_vert,
     max_small_magnitude, poly_tri_result, cos_min_angle, 
     flag_use_vcoordX11, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
      cerr << "In " << __func__ << endl;
      cerr << "  lower quad vert: ";
      IJK::print_list(cerr, lower_quad_vert, NUM_LOWER_QUAD_VERT);
      cerr << endl;
      cerr << "  upper quad vert: ";
      IJK::print_list(cerr, upper_quad_vert, NUM_UPPER_QUAD_VERT);
      cerr << endl;
      cerr << "  triangle_vert: " << triangle_vert << endl;
      cerr << "  cos_min_angle: " << cos_min_angle << endl;
      cerr << "  poly_tri_result[0].num_triangles: " 
           << poly_tri_result[0].num_triangles << endl;
      cerr << "  poly_tri_result[1].num_triangles: " 
           << poly_tri_result[1].num_triangles << endl;
      cerr << "  poly_tri_result[1].tri_vertex_index: "
           << poly_tri_result[1].tri_vertex_index << endl;
      cerr << "  poly_tri_result[2].num_triangles: " 
           << poly_tri_result[2].num_triangles << endl;
      cerr << "  poly_tri_result[2].tri_vertex_index: "
           << poly_tri_result[2].tri_vertex_index << endl;
      cerr << "  poly_tri_result[2].num_internal_tri_vertices: "
         << int(poly_tri_result[2].num_internal_tri_vertices) << endl;
      cerr << "  poly_tri_result[3].num_triangles: " 
           << poly_tri_result[3].num_triangles << endl;
      cerr << "  poly_tri_result[3].tri_vertex_index: "
           << poly_tri_result[3].tri_vertex_index << endl;
      cerr << "  poly_tri_result[1].ear_list: "
           << poly_tri_result[1].ear_list[0] << "  "
           << poly_tri_result[1].ear_list[1] << endl;
      cerr << "  poly_tri_result[2].ear_list: "
           << poly_tri_result[2].ear_list[0] << "  "
           << poly_tri_result[2].ear_list[1] << endl;
      cerr << "  flag_use_vcoordX11: "
           << int(flag_use_vcoordX11) << endl;
    */

    if (poly_tri_result[0].num_triangles == 2 && 
        poly_tri_result[1].num_triangles > 2 &&
        poly_tri_result[2].num_triangles > 2 && 
        poly_tri_result[3].num_triangles > 2) {

      if (flag_use_vcoordX11) {

        // *** DEBUG ***
        /*
        using namespace std;
        cerr << "Calling triangulate_triangle_quadx3_allow_tri5." << endl;
        */

        triangulate_triangle_quadx3_allow_tri5
          (dimension, lower_quad_vert, upper_quad_vert, triangle_vert,
           flag_reverse_orient, poly_tri_result, vcoordX.PtrConst(), 
           quadD_centroid, vertex_coord, tri_vert);
      }
      else {
        const VTYPEL * lower_quadAB_vert = lower_quad_vert;
        const VTYPEU * upper_quadAB_vert = upper_quad_vert;
        const VTYPEL * lower_quadCD_vert = lower_quad_vert+1;
        const VTYPEU * upper_quadCD_vert = upper_quad_vert+1;

        bool flag_quadCD_tri5[2];
        flag_quadCD_tri5[0] = (poly_tri_result[2].num_triangles == 5);
        flag_quadCD_tri5[1] = (poly_tri_result[3].num_triangles == 5);

        // *** DEBUG ***
        /*
        using namespace std;
        cerr << "Calling triangulate_triangle_quad_allow_tri5 and" 
             << endl;
        cerr << "  triangulate_two_quads_allow_tri5." << endl;
        */

        triangulate_triangle_quad_allow_tri5
          (dimension, lower_quadAB_vert, upper_quadAB_vert, triangle_vert,
           flag_reverse_orient, poly_tri_result[1], vcoordX00,
           quadB_centroid, vertex_coord, tri_vert);

        triangulate_two_quads_allow_tri5
          (dimension, lower_quadCD_vert, upper_quadCD_vert,
           flag_reverse_orient, flag_quadCD_tri5, vcoordX22, 
           quadC_centroid, quadD_centroid, vertex_coord, tri_vert);
      }

      flag_split[0] = flag_split[1] = flag_split[2] = flag_split[3] = true;
    }
    else if (poly_tri_result[0].num_triangles == 2 && 
             poly_tri_result[1].num_triangles > 2 &&
             poly_tri_result[2].num_triangles > 2 && 
             poly_tri_result[3].num_triangles == 2) {

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "Calling triangulate_triangle_quad_quad_allow_tri5." << endl;
      */

      triangulate_triangle_quad_quad_allow_tri5
        (dimension, lower_quad_vert, upper_quad_vert, triangle_vert,
         flag_reverse_orient, poly_tri_result+1, vcoordX00, vcoordX11, 
         quadC_centroid, vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = true;
      flag_split[3] = false;
    }
    else if (poly_tri_result[0].num_triangles == 1 && 
             poly_tri_result[1].num_triangles > 2 &&
             poly_tri_result[2].num_triangles > 2 &&
             poly_tri_result[3].num_triangles > 2) {

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "Calling triangulate_three_quads_allow_tri5." << endl;
      */

      triangulate_three_quads_allow_tri5
        (dimension, lower_quad_vert, upper_quad_vert, 
         flag_reverse_orient, poly_tri_result+1, vcoordX11, vcoordX22, 
         quadB_centroid, quadD_centroid, vertex_coord, tri_vert);

      flag_split[0] = false;
      flag_split[1] = flag_split[2] = flag_split[3] = true;
    }
    else if (poly_tri_result[0].num_triangles == 1 &&
             poly_tri_result[1].num_triangles > 2 && 
             poly_tri_result[2].num_triangles > 2) {

      bool flag_quad_tri5[2];

      flag_quad_tri5[0] = (poly_tri_result[1].num_triangles == 5);
      flag_quad_tri5[1] = (poly_tri_result[2].num_triangles == 5);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "Calling triangulate_two_quads_allow_tri5." << endl;
      */

      triangulate_two_quads_allow_tri5
        (dimension, lower_quad_vert, upper_quad_vert, 
         flag_reverse_orient, flag_quad_tri5, vcoordX11, 
         quadB_centroid, quadC_centroid, vertex_coord, tri_vert);

      flag_split[0] = flag_split[3] = false;
      flag_split[1] = flag_split[2] = true;
    }

  }

  ///@}


  // ****************************************************************
  /// @name TRIANGULATE QUAD-TRIANGLE-PENTAGON
  // ****************************************************************

  ///@{

  /// Triangulate quad-triangle-pentagon.
  /// - Add new triangles to vector tri_vert.
  /// - Add new coords to vertex_coord[].
  template <typename DTYPE, 
            typename VTYPEL, typename VTYPEU, typename TRI_VTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_quad_triangle_pentagon_allow_tri5_and_tri6
  (const DTYPE dimension,
   const VTYPEL lower_poly_vert[],
   const VTYPEU upper_poly_vert[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> poly_tri[3],
   const CTYPEX * poly_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const CTYPE * vcoordA = poly_vcoord;
    const CTYPE * vcoordB = poly_vcoord+dimension;
    const CTYPE * vcoordC = poly_vcoord+2*dimension;
    TRI_VTYPE ivB;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  lower_quad_vert: ";
    IJK::print_list(cerr, lower_poly_vert, 4);
    cerr << endl;
    cerr << "  upper_quad_vert: ";
    IJK::print_list(cerr, upper_poly_vert, 4);
    cerr << endl;
    */

    ivB = IJK::insert_coord(dimension, vcoordB, vertex_coord);

    // Add triangle B with vertex at vcoordB.
    add_triangle_vertices
      (upper_poly_vert[2], upper_poly_vert[1], ivB, flag_reverse_orient, 
       tri_vert);

    // Triangulate quadA with additional vertex at vcoordB.
    IJK_DEPRECATED::triangulate_pentagon_tri3_or_tri5
      (dimension, upper_poly_vert[1], upper_poly_vert[0], 
       lower_poly_vert[0], lower_poly_vert[1], ivB, vcoordA,
       flag_reverse_orient, poly_tri[0], vertex_coord, tri_vert);

    // Triangulate pentagonC with additional vertex at vcoordB.
    IJK_DEPRECATED::triangulate_hexagon_tri4_or_tri6
      (dimension, lower_poly_vert[1], lower_poly_vert[2], lower_poly_vert[3],
       upper_poly_vert[3], upper_poly_vert[2], ivB, vcoordC, 
       flag_reverse_orient, poly_tri[2], vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate quad-triangle-pentagon.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  - Consider triangulation of quad into 5 triangles.
   *  - Consider triangulation of pentagon into 6 triangles.
   *  - Only triangulate quads where merging with other quads 
   *    improves triangulation.
   *  - Quad A has vertices 
   *      (lower_poly_vert[0], lower_poly_vert[1], 
   *       upper_poly_vert[1], upper_poly_vert[0])
   *  - Triangle B has vertices 
   *      (lower_poly_vert[1], upper_poly_vert[2], upper_poly_vert[1])
   *  - Pentagon C has vertices 
   *      (lower_poly_vert[1], lower_poly_vert[2], lower_poly_vert[3],
   *       upper_poly_vert[3], upper_poly_vert[2])
   *  @param[out] flag_split[] flag_split[i] is true if poly i is
   *      is split into triangles which have been added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE, typename VTYPEL, typename VTYPEU, 
            typename TRI_VTYPE, typename MTYPE>
  void triangulate_quad_triangle_pentagon_max_min_angle_allow_tri5_and_tri6
  (const DTYPE dimension,
   const VTYPEL lower_poly_vert[],
   const VTYPEU upper_poly_vert[],
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert,
   bool flag_split[3])
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_PENTAGON(5);
    const int NUM_POLY(3);
    const int NUM_LOWER_POLY_VERT(4);
    const int NUM_UPPER_POLY_VERT(4);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_POLY_VERT] =
      { vertex_coord_ptr+lower_poly_vert[0]*dimension,
        vertex_coord_ptr+lower_poly_vert[1]*dimension,
        vertex_coord_ptr+lower_poly_vert[2]*dimension,
        vertex_coord_ptr+lower_poly_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_POLY_VERT] =
      { vertex_coord_ptr+upper_poly_vert[0]*dimension,
        vertex_coord_ptr+upper_poly_vert[1]*dimension,
        vertex_coord_ptr+upper_poly_vert[2]*dimension,
        vertex_coord_ptr+upper_poly_vert[3]*dimension };
    const TRI_VTYPE quadA_vert[NUM_VERT_PER_QUAD] = 
      { lower_poly_vert[0], lower_poly_vert[1], 
        upper_poly_vert[1], upper_poly_vert[0] };
    const TRI_VTYPE pentagonC_vert[NUM_VERT_PER_PENTAGON] = 
      { lower_poly_vert[1], lower_poly_vert[2], lower_poly_vert[2],
        upper_poly_vert[3], upper_poly_vert[2] };
    IJK::ARRAY<CTYPE> poly_central_coord(NUM_POLY*dimension);
    CTYPE cos_min_angle;
    POLY_TRIANGULATION_RESULT<16,CTYPE,int> poly_tri_result[NUM_POLY];
    bool flag_zero;

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = false;

    CTYPE * quadA_centroid = poly_central_coord.Ptr();
    CTYPE * coordB = poly_central_coord.Ptr()+dimension;
    CTYPE * pentagonC_centroid = poly_central_coord.Ptr()+2*dimension;

    IJK::compute_coord_centroid
      (dimension, quadA_vert, NUM_VERT_PER_QUAD, vertex_coord, quadA_centroid);
    IJK::compute_coord_centroid
      (dimension, pentagonC_vert, NUM_VERT_PER_PENTAGON, vertex_coord, 
       pentagonC_centroid);

    IJK::compute_midpoint_of_midpoints
      (dimension, lower_vcoord[1], upper_vcoord[1], 
       lower_vcoord[1], upper_vcoord[2], coordB);

    compute_cos_max_min_quad_triangle_pentagon_angle_allow_tri5_and_tri6
    (dimension, vertex_coord, poly_central_coord.PtrConst(),
     lower_poly_vert, upper_poly_vert, max_small_magnitude, 
     cos_min_angle, poly_tri_result, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
      cerr << "In " << __func__ << endl;
      cerr << "  lower poly vert: ";
      IJK::print_list(cerr, lower_poly_vert, NUM_LOWER_POLY_VERT);
      cerr << endl;
      cerr << "  upper poly vert: ";
      IJK::print_list(cerr, upper_poly_vert, NUM_UPPER_POLY_VERT);
      cerr << endl;
      cerr << "  cos_min_angle: " << cos_min_angle << endl;
      cerr << "  poly_tri_result[0].num_triangles: " 
           << poly_tri_result[0].num_triangles << endl;
      cerr << "  poly_tri_result[1].num_triangles: " 
           << poly_tri_result[1].num_triangles << endl;
      cerr << "  poly_tri_result[2].num_triangles: " 
           << poly_tri_result[2].num_triangles << endl;
      cerr << "  poly_tri_result[2].tri_vertex_index: "
           << poly_tri_result[2].tri_vertex_index << endl;
    */

    // *** SHOULD USE POLY_TRIANGULATION_RESULTX2 and .num_split_edges.

    if (poly_tri_result[0].num_triangles > 2) {

      triangulate_quad_triangle_pentagon_allow_tri5_and_tri6
          (dimension, lower_poly_vert, upper_poly_vert,
           flag_reverse_orient, poly_tri_result, 
           poly_central_coord.PtrConst(), vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = true;
    }
    else {
      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "Triangulating only the pentagon." << endl;
      */

      // Triangulate only the pentagon.
      IJK_DEPRECATED::triangulate_pentagon_tri3_or_tri5
        (dimension, lower_poly_vert[1], lower_poly_vert[2], 
         lower_poly_vert[3], upper_poly_vert[3], upper_poly_vert[2],
         pentagonC_centroid, flag_reverse_orient, poly_tri_result[2], 
         vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = false;
      flag_split[2] = true;
    }

  }

  ///@}


  // ****************************************************************
  /// @name TRIANGULATE PENTAGON-QUAD-QUAD
  // ****************************************************************

  ///@{

  /*!
   *  Triangulate pentagon-quad-quad.
   */
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename TRI_VTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_pentagon_quad_quad
  (const DTYPE dimension,
   const VTYPEL lower_poly_vert[],
   const VTYPEU upper_poly_vert[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> 
   poly_tri_result[3],
   const CTYPEX * vcoordX12,
   const CTYPEX * vcoordX23,
   const CTYPEX * pentagonA_vcoord,
   const CTYPEX * quadB_vcoord,
   const CTYPEX * quadC_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const TRI_VTYPE * lower_quadB_vert = lower_poly_vert+1;
    const TRI_VTYPE * upper_quadB_vert = upper_poly_vert+2;

    const TRI_VTYPE ivertX12 = 
      IJK::insert_coord(dimension, vcoordX12, vertex_coord);
    const TRI_VTYPE ivertX23 = 
      IJK::insert_coord(dimension, vcoordX23, vertex_coord);

    // Triangulate pentagonA with additional vertex at ivertX12.
    IJK_DEPRECATED::triangulate_hexagon_tri4_or_tri6
      (dimension, upper_poly_vert[2], upper_poly_vert[1], upper_poly_vert[0],
       lower_poly_vert[0], lower_poly_vert[1], ivertX12, pentagonA_vcoord, 
       flag_reverse_orient, poly_tri_result[0], vertex_coord, tri_vert);

    // Triangulate subquads of quad B.
    triangulate_subquads_using_diagonals_LU
      (lower_quadB_vert, upper_quadB_vert, ivertX12, ivertX23,
       poly_tri_result[1], flag_reverse_orient, tri_vert);

    // Triangulate quadC with additional vertex at ivertX23.
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, lower_poly_vert[2], lower_poly_vert[3], 
       upper_poly_vert[4], upper_poly_vert[3], ivertX23, quadC_vcoord, 
       flag_reverse_orient, poly_tri_result[2], vertex_coord, tri_vert);
  }


  /*!
   *  Triangulate pentagon-quad-quad.
   */
  template <typename DTYPE, typename VTYPEL, typename VTYPEU, 
            typename TRI_VTYPE, 
            typename COS_TYPE, typename NTYPE, int BIT_SET_SIZE, 
            typename CTYPE, typename CTYPEX>
  void triangulate_pentagon_quad
  (const DTYPE dimension,
   const VTYPEL lower_poly_vert[],
   const VTYPEU upper_poly_vert[],
   const bool flag_reverse_orient,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   pentagonA_tri_result,
   const POLY_TRIANGULATION_RESULTX2<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quadB_tri_result,
   const CTYPEX * vcoordX12,
   const CTYPEX * quadA_vcoord,
   const CTYPEX * quadB_vcoord,
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert)
  {
    const TRI_VTYPE ivertX12 = 
      IJK::insert_coord(dimension, vcoordX12, vertex_coord);

    // Triangulate pentagonA with additional vertex at ivertX12.
    IJK_DEPRECATED::triangulate_hexagon_tri4_or_tri6
      (dimension, upper_poly_vert[2], upper_poly_vert[1], upper_poly_vert[0],
       lower_poly_vert[0], lower_poly_vert[1], ivertX12, quadA_vcoord, 
       flag_reverse_orient, pentagonA_tri_result, vertex_coord, tri_vert);

    // Triangulate quadB with additional vertex at ivertX12.
    triangulate_quad_tri3_or_tri5_vX03
      (dimension, lower_poly_vert[1], lower_poly_vert[2], 
       upper_poly_vert[3], upper_poly_vert[2], ivertX12, quadB_vcoord, 
       flag_reverse_orient, quadB_tri_result, vertex_coord, tri_vert);
  }

  /*!
   *  Triangulate pentagon-quad-quad.
   *  - Maximize minimum triangle angle.
   *  - Add new triangles to vector tri_vert.
   *  - Add new coords, if necessary, to vertex_coord[].
   *  - Consider triangulation of quad into 5 triangles.
   *  - Consider triangulation of pentagon into 6 triangles.
   *  - Only triangulate quads where merging with other quads 
   *    improves triangulation.
   *  - Quad A has vertices 
   *      (lower_poly_vert[0], lower_poly_vert[1], 
   *       upper_poly_vert[1], upper_poly_vert[0])
   *  - Triangle B has vertices 
   *      (lower_poly_vert[1], upper_poly_vert[2], upper_poly_vert[1])
   *  - Pentagon C has vertices 
   *      (lower_poly_vert[1], lower_poly_vert[2], lower_poly_vert[3],
   *       upper_poly_vert[3], upper_poly_vert[2])
   *  @param[out] flag_split[] flag_split[i] is true if poly i is
   *      is split into triangles which have been added to tri_vert[].
  */
  template <typename DTYPE, typename CTYPE, typename VTYPEL, typename VTYPEU, 
            typename TRI_VTYPE, typename MTYPE>
  void triangulate_pentagon_quad_quad_max_min_angle
  (const DTYPE dimension,
   const VTYPEL lower_poly_vert[],
   const VTYPEU upper_poly_vert[],
   const bool flag_reverse_orient,
   const MTYPE max_small_magnitude,
   const bool block_diagonal[3][2],
   std::vector<CTYPE> & vertex_coord,
   std::vector<TRI_VTYPE> & tri_vert,
   bool flag_split[3])
  {
    const int NUM_VERT_PER_QUAD(4);
    const int NUM_VERT_PER_PENTAGON(5);
    const int NUM_POLY(3);
    const int NUM_LOWER_POLY_VERT(4);
    const int NUM_UPPER_POLY_VERT(5);
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);
    const CTYPE * lower_vcoord[NUM_LOWER_POLY_VERT] =
      { vertex_coord_ptr+lower_poly_vert[0]*dimension,
        vertex_coord_ptr+lower_poly_vert[1]*dimension,
        vertex_coord_ptr+lower_poly_vert[2]*dimension,
        vertex_coord_ptr+lower_poly_vert[3]*dimension };
    const CTYPE * upper_vcoord[NUM_UPPER_POLY_VERT] =
      { vertex_coord_ptr+upper_poly_vert[0]*dimension,
        vertex_coord_ptr+upper_poly_vert[1]*dimension,
        vertex_coord_ptr+upper_poly_vert[2]*dimension,
        vertex_coord_ptr+upper_poly_vert[3]*dimension,
        vertex_coord_ptr+upper_poly_vert[4]*dimension };
    const TRI_VTYPE pentagonA_vert[NUM_VERT_PER_PENTAGON] = 
      { lower_poly_vert[0], lower_poly_vert[1], 
        upper_poly_vert[2], upper_poly_vert[1], upper_poly_vert[0] };
    const TRI_VTYPE quadB_vert[NUM_VERT_PER_QUAD] = 
      { upper_poly_vert[3], upper_poly_vert[2], 
        lower_poly_vert[1], lower_poly_vert[2] };
    const TRI_VTYPE quadC_vert[NUM_VERT_PER_QUAD] = 
      { lower_poly_vert[2], lower_poly_vert[3], 
        upper_poly_vert[4], upper_poly_vert[3] };
        
    IJK::ARRAY<CTYPE> vcoordX(2*dimension);
    IJK::ARRAY<CTYPE> poly_central_coord(NUM_POLY*dimension);
    POLY_TRIANGULATION_RESULTX2<16,CTYPE,int> poly_tri_result[NUM_POLY];
    CTYPE cos_max_min_angle;
    bool flag_zero;

    // Initialize.
    flag_split[0] = flag_split[1] = flag_split[2] = false;

    CTYPE * const vcoordX12 = vcoordX.Ptr();
    CTYPE * const vcoordX23 = vcoordX.Ptr()+dimension;
    CTYPE * const pentagonA_centroid = poly_central_coord.Ptr();
    CTYPE * const quadB_centroid = poly_central_coord.Ptr()+dimension;
    CTYPE * const quadC_centroid = poly_central_coord.Ptr()+2*dimension;

    IJK::compute_midpoint
      (dimension, lower_vcoord[1], upper_vcoord[2], vcoordX12);
    IJK::compute_midpoint
      (dimension, lower_vcoord[2], upper_vcoord[3], vcoordX23);

    IJK::compute_coord_centroid
      (dimension, pentagonA_vert, NUM_VERT_PER_PENTAGON, vertex_coord, 
       pentagonA_centroid);
    IJK::compute_coord_centroid
      (dimension, quadB_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadB_centroid);
    IJK::compute_coord_centroid
      (dimension, quadC_vert, NUM_VERT_PER_QUAD, vertex_coord, 
       quadC_centroid);

    compute_cos_max_min_pentagon_quad_quad_angle
      (dimension, vertex_coord, vcoordX12, vcoordX23,
       poly_central_coord.PtrConst(), lower_poly_vert, upper_poly_vert, 
       max_small_magnitude, block_diagonal, poly_tri_result,
       cos_max_min_angle, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  lower_poly_vert: ";
    IJK::print_list(cerr, lower_poly_vert, NUM_LOWER_POLY_VERT);
    IJK::print_list(cerr, upper_poly_vert, NUM_UPPER_POLY_VERT);
    cerr << "\n";
    */

    // Initialize
    flag_split[0] = flag_split[1] = flag_split[2];

    if (poly_tri_result[1].num_split_edges == 2) {
      triangulate_pentagon_quad_quad
        (dimension, lower_poly_vert, upper_poly_vert, 
         flag_reverse_orient, poly_tri_result,
         vcoordX12, vcoordX23, 
         pentagonA_centroid, quadB_centroid, quadC_centroid,
         vertex_coord, tri_vert);

      flag_split[0] = flag_split[1] = flag_split[2] = true;
    }
    else if (poly_tri_result[1].num_split_edges == 1) {

      if (poly_tri_result[0].num_split_edges == 1 &&
          poly_tri_result[2].num_split_edges == 0) {
        triangulate_pentagon_quad
          (dimension, lower_poly_vert, upper_poly_vert, 
           flag_reverse_orient, poly_tri_result[0], poly_tri_result[1],
           vcoordX12, pentagonA_centroid, quadB_centroid,
           vertex_coord, tri_vert);

        flag_split[0] = flag_split[1] = true;
      }
      else if (poly_tri_result[0].num_split_edges == 0 &&
               poly_tri_result[2].num_split_edges == 1) {
        const VTYPEL * lower_quadBC_vert = lower_poly_vert+1;
        const VTYPEU * upper_quadBC_vert = upper_poly_vert+2;
        
        triangulate_two_quads_allow_tri5
          (dimension, lower_quadBC_vert, upper_quadBC_vert, 
           flag_reverse_orient, poly_tri_result[1], poly_tri_result[2],
           vcoordX23, quadB_centroid, quadC_centroid, 
           vertex_coord, tri_vert);

        flag_split[1] = flag_split[2] = true;
      }
    }

  }

  ///@}

}


#endif
