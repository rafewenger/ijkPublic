/*!
 *  @file ijktri2D_max_min_angle.tpp
 *  @brief ijk templates for triangulations that maximize min angles in 2D meshes.
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

#ifndef _IJKTRI2D_MAX_MIN_ANGLE_
#define _IJKTRI2D_MAX_MIN_ANGLE_

#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkisocoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijktri2D_angle.tpp"
#include "ijktri2D_info.tpp"
#include "ijktriangulate_poly2D.tpp"


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // ***************************************************************
  //! @name MAXIMIZE MIN TRIANGULATION ANGLE
  // ***************************************************************

  //@{

  /*!
   *  @brief Compute fan triangulation vertex which maximizes 
   *    min triangulation angle.
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
  void compute_fan_triangulation_to_max_min_angle
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
      compute_cos_min_fan_triangulation_angle
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
   *  @brief Compute fan triangulation vertex which maximizes 
   *    min triangulation angle.
   *  - Version that returns POLYGON_TRIANGULATION_RESULT.
   */
  template <typename DTYPE, typename CTYPE, 
            typename NTYPE1, typename NTYPE2,
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_fan_triangulation_to_max_min_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const NTYPE1 num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude, 
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE, NTYPE2> & 
   poly_tri_result)
  {
    NTYPE1 iv_max_min;
    COS_TYPE cos_min_angle;
    bool flag_zero;

    poly_tri_result.Clear();

    compute_fan_triangulation_to_max_min_angle
      (dimension, vcoord, num_poly_vert, poly_vert,
       max_small_magnitude, iv_max_min, cos_min_angle, flag_zero);

    poly_tri_result.flag_zero = flag_zero;
    if (!flag_zero) {
      poly_tri_result.SetNoInterior
        (cos_min_angle, num_poly_vert-2, iv_max_min);
      poly_tri_result.triangulation_encoding.SetFan
        (iv_max_min, num_poly_vert);
    }
  }


  /*!
   *  @brief Compute fan triangulation vertex which maximizes 
   *    min triangulation angle.
   *  - Version using list of pointers to vertex coordinates.
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
            typename CTYPE, typename MTYPE, 
            typename IVTYPE, typename COS_TYPE>
  void compute_fan_triangulation_to_max_min_angleP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE vcoord_ptr[], const MTYPE max_small_magnitude, 
   IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    // Initialize
    cos_min_angle = 1;
    flag_zero = true;
    iv_max_min = 0;

    for (IVTYPE iv = 0; iv < num_poly_vert; iv++) {
      compute_cos_min_fan_triangulation_angleP
        (dimension, num_poly_vert, vcoord_ptr, 
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
   *  @brief Compute fan triangulation vertex which maximizes 
   *    min triangulation angle.
   *  - Version that returns POLYGON_TRIANGULATION_RESULT.
   *  - Version using list of pointers to vertex coordinates.
   */
  template <typename DTYPE, typename CTYPE, 
            typename NTYPE1, typename NTYPE2, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_fan_triangulation_to_max_min_angleP
  (const DTYPE dimension, const NTYPE1 num_poly_vert,
   const CTYPE vcoord_ptr[], const MTYPE max_small_magnitude, 
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE, NTYPE2> & 
   poly_tri_result)
  {
    NTYPE1 iv_max_min;
    COS_TYPE cos_min_angle;
    bool flag_zero;

    poly_tri_result.Clear();

    compute_fan_triangulation_to_max_min_angleP
      (dimension, num_poly_vert, vcoord_ptr, 
       max_small_magnitude, iv_max_min, cos_min_angle, flag_zero);

    poly_tri_result.flag_zero = flag_zero;
    if (!flag_zero) {
      poly_tri_result.SetNoInterior
        (cos_min_angle, num_poly_vert-2, iv_max_min);
      poly_tri_result.triangulation_encoding.SetFan
        (iv_max_min, num_poly_vert);
    }
  }


  /*!
   *  @brief Compute fan triangulation vertex which maximizes 
   *    min triangulation angle.
   *  - Version using list of pointers to vertex coordinates.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param num_poly_vert Number of polygon vertices.
   *  @param poly_vert[] List of polygon vertices.
   *  @param vcoord[] Polygon vertex coordinates.
   *         vcoord[i*dimension+j] is j'th coordinate of i'th vertex.
   *  @param flag_not_ear[i] If true, do not consider any triangulation
   *    that has vertex i as an ear.
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
            typename CTYPE, typename BIT_SET_TYPE, typename MTYPE, 
            typename IVTYPE, typename COS_TYPE>
  void compute_fan_triangulation_avoid_ears_to_max_min_angleP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE vcoord_ptr[], const BIT_SET_TYPE & flag_not_ear,
   const MTYPE max_small_magnitude, 
   IVTYPE & iv_max_min, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    COS_TYPE cosA;
    bool flagA;

    // Initialize
    cos_min_angle = 1;
    flag_zero = true;
    iv_max_min = 0;

    for (IVTYPE iv1 = 0; iv1 < num_poly_vert; iv1++) {

      const IVTYPE iv0 = (iv1+(num_poly_vert-1))%num_poly_vert;
      const IVTYPE iv2 = (iv1+1)%num_poly_vert;

      if (flag_not_ear[iv0] || flag_not_ear[iv2]) {
        // Creating fan triangulation at iv1 will 
        //   cause iv0 and iv2 to be ears.
        // Skip fan triangulation at iv1.
        continue;
      }

      compute_cos_min_fan_triangulation_angleP
        (dimension, num_poly_vert, vcoord_ptr, 
         max_small_magnitude, iv1, cosA, flagA);

      if (!flagA && cosA < cos_min_angle) {
        flag_zero = false;
        iv_max_min = iv1;
        cos_min_angle = cosA;
      }
    }

    return;
  }


  /*!
   *  @brief Compute fan triangulation vertex which maximizes 
   *    min triangulation angle.
   *  - Consider only triangulations that do not create ears at flagged vertices.
   *  - Version that returns POLYGON_TRIANGULATION_RESULT.
   *  - Version using list of pointers to vertex coordinates.
   */
  template <typename DTYPE, typename CTYPE, 
            typename NTYPE1, typename NTYPE2, 
	    typename BIT_SET_TYPE, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_fan_triangulation_avoid_ears_to_max_min_angleP
  (const DTYPE dimension, const NTYPE1 num_poly_vert,
   const CTYPE vcoord_ptr[], const BIT_SET_TYPE & flag_no_ears,
   const MTYPE max_small_magnitude, 
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE, NTYPE2> & 
   poly_tri_result)
  {
    NTYPE1 iv_max_min;
    COS_TYPE cos_min_angle;
    bool flag_zero;

    poly_tri_result.Clear();

    compute_fan_triangulation_avoid_ears_to_max_min_angleP
      (dimension, num_poly_vert, vcoord_ptr, flag_no_ears,
       max_small_magnitude, iv_max_min, cos_min_angle, flag_zero);

    poly_tri_result.flag_zero = flag_zero;
    if (!flag_zero) {
      poly_tri_result.SetNoInterior
        (cos_min_angle, num_poly_vert-2, iv_max_min);
      poly_tri_result.triangulation_encoding.SetFan
        (iv_max_min, num_poly_vert);
    }
  }


  /*!
   *  @brief Compute triangulation that maximizes 
   *    the minimum triangle angle.
   *  - Consider all possible polygon triangulations.
   *  @param[out] cos_min_angle Cosine of the minimum triangle angle
   *    over all triangles in the triangulation.
   *  - Note: Cosine of the minimum triangle angle is the maximum
   *    cosize over all angles of all triangles in the triangulation.
   *  @param[out] flag_zero True if all triangulations have
   *    some edge with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE,
            typename NTYPE, typename VTYPE, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_over_all_to_max_min_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const NTYPE num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero,
   POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> & tri_encoding)
  {
    POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> tri_encodingX;
    COS_TYPE cosX;
    bool flagX_zero;

    // Initialize
    cos_min_angle = -1.0;
    flag_zero = false;

    tri_encodingX.FirstEncoding(num_poly_vert);
    compute_cos_min_encoded_triangulation_angle
      (dimension, vcoord, num_poly_vert, poly_vert,
       tri_encodingX, max_small_magnitude, cos_min_angle, flag_zero);
    tri_encoding = tri_encodingX;

    while(tri_encodingX.NextEncoding(num_poly_vert)) {

      compute_cos_min_encoded_triangulation_angle
        (dimension, vcoord, num_poly_vert, poly_vert,
         tri_encodingX, max_small_magnitude, cosX, flagX_zero);

      if (!flagX_zero) {

        if (flag_zero || (cos_min_angle > cosX)) {
          cos_min_angle = cosX;
          flag_zero = flagX_zero;
          tri_encoding = tri_encodingX;
        }
      }
    }

  }


  /*!
   *  @brief Compute triangulation that maximizes 
   *    the minimum triangle angle.
   *  - Consider all possible polygon triangulations.
   *  - Version that returns POLYGON_TRIANGULATION_RESULT.
   */
  template <typename DTYPE, typename CTYPE,
            typename NTYPE1, typename NTYPE2,
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_over_all_to_max_min_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const NTYPE1 num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude,
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE2> &
   poly_tri_result)
  {
    POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> tri_encoding;
    COS_TYPE cos_min_angle;
    bool flag_zero;

    compute_triangulation_over_all_to_max_min_angle
      (dimension, vcoord, num_poly_vert, poly_vert, max_small_magnitude,
       cos_min_angle, flag_zero, tri_encoding);

    poly_tri_result.Clear();
    poly_tri_result.flag_zero = flag_zero;
 
    if (!flag_zero) {
      poly_tri_result.cos_min_triangulation_angle = cos_min_angle;
      poly_tri_result.triangulation_encoding = tri_encoding;
      poly_tri_result.num_triangles = num_poly_vert - 2;
    }
  }


  /*!
   *  @brief Compute triangulation that maximizes 
   *    the minimum triangle angle.
   *  - Consider all possible polygon triangulations.
   *  - Version using list of pointers to vertex coordinates.
   *  @param[out] cos_min_angle Cosine of the minimum triangle angle
   *    over all triangles in the triangulation.
   *  @param[out] flag_zero True if all triangulations have
   *    some edge with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE,
            typename NTYPE, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_over_all_to_max_min_angleP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPE * const vcoord_ptr[], 
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero,
   POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> & tri_encoding)
  {
    POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> tri_encodingX;
    COS_TYPE cosX;
    bool flagX_zero;

    // Initialize
    cos_min_angle = -1.0;
    flag_zero = false;

    tri_encodingX.FirstEncoding(num_poly_vert);
    compute_cos_min_encoded_triangulation_angleP
      (dimension, num_poly_vert, vcoord_ptr, 
       tri_encodingX, max_small_magnitude, cos_min_angle, flag_zero);
    tri_encoding = tri_encodingX;

    // *** DEBUG ***
    /*
    using namespace std;
    tri_encoding.PrintTriangulationRepresentation
      (cerr, num_poly_vert, "first_tri_encoding (B): ", ".  ");
    cerr << " cos_min_angle: " << cos_min_angle << endl;
    */

    while(tri_encodingX.NextEncoding(num_poly_vert)) {

      compute_cos_min_encoded_triangulation_angleP
        (dimension, num_poly_vert, vcoord_ptr,
         tri_encodingX, max_small_magnitude, cosX, flagX_zero);

      if (!flagX_zero) {

        // *** DEBUG ***
        /*
        using namespace std;
        tri_encodingX.PrintTriangulationRepresentation
          (cerr, num_poly_vert, "  Encoding: ", "");
        cerr << "  cosX: " << cosX << endl;
        */

        if (flag_zero || (cos_min_angle > cosX)) {
          cos_min_angle = cosX;
          flag_zero = flagX_zero;
          tri_encoding = tri_encodingX;
        }
      }
    }

  }


  /*!
   *  @brief Compute triangulation that maximizes 
   *    the minimum triangle angle.
   *  - Consider all possible polygon triangulations.
   *  - Version using list of pointers to vertex coordinates.
   *  - Version that returns POLYGON_TRIANGULATION_RESULT.
   */
  template <typename DTYPE, typename CTYPE,
            typename NTYPE1, typename NTYPE2, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_over_all_to_max_min_angleP
  (const DTYPE dimension, const NTYPE1 num_poly_vert, 
   const CTYPE vcoord_ptr[], 
   const MTYPE max_small_magnitude,
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE2> &
   poly_tri_result)
  {
    POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> tri_encoding;
    COS_TYPE cos_min_angle;
    bool flag_zero;

    compute_triangulation_over_all_to_max_min_angleP
      (dimension, num_poly_vert, vcoord_ptr, max_small_magnitude,
       cos_min_angle, flag_zero, tri_encoding);

    poly_tri_result.Clear();
    poly_tri_result.flag_zero = flag_zero;
 
    if (!flag_zero) {
      poly_tri_result.cos_min_triangulation_angle = cos_min_angle;
      poly_tri_result.triangulation_encoding = tri_encoding;
      poly_tri_result.num_triangles = num_poly_vert - 2;
    }
  }


  /*!
   *  @brief Compute triangulation that maximizes 
   *    the minimum triangle angle.
   *  - Consider only triangulations that do not create ears at flagged vertices.
   *  - Version using list of pointers to vertex coordinates.
   *  @param flag_not_ear[i] If true, do not consider any triangulation
   *    that has vertex i as an ear.
   *  @param[out] cos_min_angle Cosine of the minimum triangle angle
   *    over all triangles in the triangulation.
   *  @param[out] flag_zero True if all triangulations have
   *    some edge with length less than or equal to max_small_magnitude.
   *  @param[out] flag_found True if some triangulation was found.
   */
  template <typename DTYPE, typename CTYPE,
            typename NTYPE, typename BIT_SET_TYPE, typename MTYPE,
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_avoid_ears_to_max_min_angleP
  (const DTYPE dimension, const NTYPE num_poly_vert,
   const CTYPE * const vcoord_ptr[],
   const BIT_SET_TYPE & flag_not_ear,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle, bool & flag_zero,
   POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> & tri_encoding,
   bool & flag_found)
  {
    POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> tri_encodingX;
    COS_TYPE cosX;
    bool flagX_zero;
    typename POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE>::BIT_SET_TYPE flag_ears;

    // Initialize
    cos_min_angle = -1.0;
    flag_zero = false;
    flag_found = false;
    tri_encoding.Clear();

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  In " << __func__ << endl;
    cerr << "  flag_not_ear: ";
    for (int i = 0; i < num_poly_vert; i++) 
      { cerr << " " << flag_not_ear[i]; }
    cerr << endl;
    */

    tri_encodingX.FirstEncoding(num_poly_vert);
    tri_encodingX.GetFlagEars(num_poly_vert, flag_ears);
    if (!does_some_true_bit_match(flag_not_ear, flag_ears, num_poly_vert)) {
      compute_cos_min_encoded_triangulation_angleP
        (dimension, num_poly_vert, vcoord_ptr,
         tri_encodingX, max_small_magnitude, cos_min_angle, flag_zero);
      tri_encoding = tri_encodingX;
      flag_found = true;

      // *** DEBUG ***
      /*
      using namespace std;
      tri_encodingX.PrintTriangulationRepresentation
      (cerr, num_poly_vert, "    Encoding: ", "\n");
      */
    }

    // *** DEBUG ***
    /*
    using namespace std;
    tri_encoding.PrintTriangulationRepresentation
      (cerr, num_poly_vert, "first_tri_encoding (B): ", ".  ");
    cerr << " cos_min_angle: " << cos_min_angle << endl;
    */

    while(tri_encodingX.NextEncoding(num_poly_vert)) {

      tri_encodingX.GetFlagEars(num_poly_vert, flag_ears);
      if (!does_some_true_bit_match(flag_not_ear, flag_ears, num_poly_vert)) {

        compute_cos_min_encoded_triangulation_angleP
          (dimension, num_poly_vert, vcoord_ptr,
           tri_encodingX, max_small_magnitude, cosX, flagX_zero);

        if (!flag_found || !flagX_zero) {

          // *** DEBUG ***
          /*
          using namespace std;
          tri_encodingX.PrintTriangulationRepresentation
            (cerr, num_poly_vert, "  Encoding: ", "");
          cerr << "  cosX: " << cosX << endl;
          */

          if (!flag_found || (flag_zero || (cos_min_angle > cosX))) {
            cos_min_angle = cosX;
            flag_zero = flagX_zero;
            tri_encoding = tri_encodingX;
          }
          flag_found = true;
        }
      }
    }

  }


  /*!
   *  @brief Compute triangulation that maximizes 
   *    the minimum triangle angle.
   *  - Consider only triangulations that do not create ears at flagged vertices.
   *  - Version using list of pointers to vertex coordinates.
   *  - Version that returns POLYGON_TRIANGULATION_RESULT.
   *  @param flag_not_ear[i] If true, do not consider any triangulation
   *    that has vertex i as an ear.
   *  @param[out] poly_tri_result Triangulation that maximizes minimum triangle angle.
   *     - Note: If every triangulation has some flagged
   *       vertex as an ear, then poly_tri_result will have 0 triangles.
   */
  template <typename DTYPE, typename CTYPE,
            typename NTYPE1, typename NTYPE2, typename BIT_SET_TYPE,
            typename MTYPE, typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_avoid_ears_to_max_min_angleP
  (const DTYPE dimension, const NTYPE1 num_poly_vert,
   const CTYPE vcoord_ptr[],
   const BIT_SET_TYPE & flag_not_ear,
   const MTYPE max_small_magnitude,
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE2> &
   poly_tri_result)
  {
    POLYGON_TRIANGULATION_ENCODING<BIT_SET_SIZE> tri_encoding;
    COS_TYPE cos_min_angle;
    bool flag_zero, flag_found;

    compute_triangulation_avoid_ears_to_max_min_angleP
      (dimension, num_poly_vert, vcoord_ptr, flag_not_ear, max_small_magnitude,
       cos_min_angle, flag_zero, tri_encoding, flag_found);

    poly_tri_result.Clear();
    if (flag_found) {

      poly_tri_result.flag_zero = flag_zero;

      if (!flag_zero) {
        poly_tri_result.cos_min_triangulation_angle = cos_min_angle;
        poly_tri_result.triangulation_encoding = tri_encoding;
        poly_tri_result.num_triangles = num_poly_vert - 2;
      }
    }
  }


  /*!
   *  @brief Compute triangulation that maximizes min angle 
   *    using interior point.
   *  - Triangulation is from point at vcoordX0[] or vcoordX1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param poly_vcoord[] Pointers to polygon vertex coordinates.
   *         poly_vcoord[i] is a pointer to the i'th vertex coordinates.
   *  @param num_poly_vert Number of polygon vertices.
   *  @param vcoordX0[] Internal vertex coordinates 0.
   *  @param vcoordX1[] Internal vertex coordinates 1.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *    - Note: Cosine of the smallest angle is the largest cosine of an angle.
   *  @param[out] index_selected_IV Index of selected internal
   *           vertex_coordinates. Either 0 or 1.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename NTYPE, 
            typename CTYPEV, typename CTYPEX0, typename CTYPEX1,
            typename MTYPE, typename COS_TYPE, typename ITYPES>
  void compute_triangulation_from_vcoordx2_to_max_min_angleP
  (const DTYPE dimension, const NTYPE num_poly_vert, 
   const CTYPEV * const vcoord_ptr[], 
   const CTYPEX0 vcoordX0[], 
   const CTYPEX1 vcoordX1[], 
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, 
   ITYPES & index_selected_IV,
   bool & flag_zero)
  {
    COS_TYPE cosB;
    bool flagB_zero;

    index_selected_IV = 0;
    compute_cos_min_triangulation_angle_using_vcoordXP
      (dimension, num_poly_vert, vcoord_ptr, vcoordX0,
       max_small_magnitude, cos_min_angle, flag_zero);

    compute_cos_min_triangulation_angle_using_vcoordXP
      (dimension, num_poly_vert, vcoord_ptr, vcoordX1,
       max_small_magnitude, cosB, flagB_zero);

    if (flag_zero || (!flagB_zero && cosB < cos_min_angle)) {
      cos_min_angle = cosB;
      flag_zero = flagB_zero;
      index_selected_IV = 1;
    }
  }


  /*!
   *  @brief Compute triangulation to max min angle.
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE1, typename NTYPE2,
            typename VTYPE, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_to_max_min_angle
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 interior_vcoord[], 
   const NTYPE1 num_poly_vert, const VTYPE poly_vert[],
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE2> & 
   poly_tri_result)
  {
    COS_TYPE cos_min_angle;
    bool flag_zero;

    poly_tri_result.Clear();

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  In " << __func__ << endl;
    IJK::print_list(cerr, "  Poly vert: ", poly_vert, num_poly_vert, "\n");
    */

    if (triangulation_settings.flag_all_triangulations_from_polygon_vertices) {
      compute_triangulation_over_all_to_max_min_angle
        (dimension, vcoord, num_poly_vert, poly_vert,
         max_small_magnitude, poly_tri_result);
    }
    else {
      compute_fan_triangulation_to_max_min_angle
        (dimension, vcoord, num_poly_vert, poly_vert, 
         max_small_magnitude, poly_tri_result);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "  After computer_max_min_fan_triangulation_angle" << endl;
      poly_tri_result.Print(cerr, "    ", num_poly_vert);
      */
    }

    if (triangulation_settings.flag_triangulate_from_interior_vertex) {

      compute_cos_min_triangulation_angle_using_vcoordX
        (dimension, vcoord, interior_vcoord, num_poly_vert, poly_vert,
         max_small_magnitude, cos_min_angle, flag_zero);

      if (!flag_zero) {
        if (poly_tri_result.cos_min_triangulation_angle > cos_min_angle) {
          poly_tri_result.SetInterior(cos_min_angle, num_poly_vert, 1);
        }
      }
    }

    if (triangulation_settings.flag_split_ear_triangulate_from_interior_vertex) {
      NTYPE1 iloc_min, iloc_max;
      COS_TYPE cos_min, cos_max;
      NTYPE1 num_angle;

      IJK::compute_cos_min_max_polygon_angles
        (dimension, vcoord, poly_vert, num_poly_vert, max_small_magnitude,
         cos_min, cos_max, iloc_min, iloc_max, num_angle);

      if (num_angle == num_poly_vert) {

        compute_cos_min_triangulation_angle_using_cut_ear_vcoordX
          (dimension, vcoord, interior_vcoord, num_poly_vert, poly_vert,
           iloc_min, max_small_magnitude, cos_min_angle, flag_zero);

        if (!flag_zero) {
          if (poly_tri_result.cos_min_triangulation_angle > cos_min_angle) {
            poly_tri_result.SetInterior(cos_min_angle, num_poly_vert, 1);
            poly_tri_result.SetEarFlag(iloc_min, true);
          }
        }
      }
    }

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_list(cerr, "  Poly vert: ", poly_vert, num_poly_vert, "\n");
    poly_tri_result.Print(cerr, "    ", num_poly_vert);
    */
  }


  /*!
   *  @brief Compute triangulation to max min angle.
   *  - Version using list of pointers to vertex coordinates.
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE1, typename NTYPE2,
            typename BIT_SET_TYPE, typename MTYPE,
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_to_max_min_angleP
  (const DTYPE dimension, const NTYPE1 num_poly_vert, 
   const CTYPE0 * const vcoord_ptr[], const CTYPE1 interior_vcoord[],
   const BIT_SET_TYPE & flag_not_ear,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE2> & 
   poly_tri_result)
  {
    const NTYPE1 NUM_VERT_IN_HEXAGON(6);
    COS_TYPE cos_min_angle;
    bool flag_zero;

    poly_tri_result.Clear();

    const bool flag_avoid_ears =
      is_some_bit_true(flag_not_ear, num_poly_vert);


    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    triangulation_settings.Print(cerr, "  ");
    */

    if (triangulation_settings.flag_all_triangulations_from_polygon_vertices ||
        (triangulation_settings.flag_all_triangulations_from_hexagon_vertices &&
         num_poly_vert == NUM_VERT_IN_HEXAGON)) {

      if (flag_avoid_ears) {
        compute_triangulation_avoid_ears_to_max_min_angleP
          (dimension, num_poly_vert, vcoord_ptr, flag_not_ear,
           max_small_magnitude, poly_tri_result);
      }
      else {
        compute_triangulation_over_all_to_max_min_angleP
          (dimension, num_poly_vert, vcoord_ptr,
           max_small_magnitude, poly_tri_result);
      }
    }
    else {
      if (flag_avoid_ears) {
        compute_fan_triangulation_avoid_ears_to_max_min_angleP
        (dimension, num_poly_vert, vcoord_ptr, flag_not_ear,
         max_small_magnitude, poly_tri_result);
      }
      else {
        compute_fan_triangulation_to_max_min_angleP
        (dimension, num_poly_vert, vcoord_ptr,
         max_small_magnitude, poly_tri_result);
      }
    }

    if (triangulation_settings.flag_triangulate_from_interior_vertex) {

      compute_cos_min_triangulation_angle_using_vcoordXP
        (dimension, num_poly_vert, vcoord_ptr, interior_vcoord,
         max_small_magnitude, cos_min_angle, flag_zero);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "  cos_min_angle (I): " << cos_min_angle << endl;
      */

      if (!flag_zero) {
        if (poly_tri_result.cos_min_triangulation_angle > cos_min_angle) {
          poly_tri_result.SetInterior(cos_min_angle, num_poly_vert, 1);
        }
      }
    }

    if (triangulation_settings.flag_split_ear_triangulate_from_interior_vertex) {
      NTYPE1 iloc_min, iloc_max;
      COS_TYPE cos_min, cos_max;
      NTYPE1 num_angle;

      IJK::compute_cos_min_max_polygon_angles
        (dimension, num_poly_vert, vcoord_ptr, max_small_magnitude,
         cos_min, cos_max, iloc_min, iloc_max, num_angle);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "  num_angle: " << num_angle << endl;
      cerr << "  iloc_min: " << iloc_min << endl;
      cerr << "  num_poly_vert: " << num_poly_vert << endl;
      IJK::print_coord3D(cerr, "  vcoord[0]: ",
                         vcoord_ptr[0], "\n");
      IJK::print_coord3D(cerr, "  vcoord[1]: ",
                         vcoord_ptr[1], "\n");
      IJK::print_coord3D(cerr, "  coord of min angle: ",
                         vcoord_ptr[iloc_min], "\n");
      IJK::print_coord3D(cerr, "  interior: ", interior_vcoord, "\n");
      */
      
      if (num_angle == num_poly_vert) {

        compute_cos_min_triangulation_angle_using_cut_ear_vcoordXP
          (dimension, num_poly_vert, vcoord_ptr, interior_vcoord,
           iloc_min, max_small_magnitude, cos_min_angle, flag_zero);

        // *** DEBUG ***
        /*
        using namespace std;
        cerr << "  cut_ear cos_min_angle: " << cos_min_angle << endl;
        */

        if (!flag_zero) {
          if (poly_tri_result.cos_min_triangulation_angle > cos_min_angle) {
            poly_tri_result.SetInterior(cos_min_angle, num_poly_vert, 1);
            poly_tri_result.SetEarFlag(iloc_min, true);
          }
        }
      }
    }

  }

  
  /*!
   *  @brief Compute triangulation to max min angle.
   *  - Version using list of pointers to vertex coordinates.
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE1, typename NTYPE2, typename MTYPE,
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_to_max_min_angleP
  (const DTYPE dimension, const NTYPE1 num_poly_vert,
   const CTYPE0 * const vcoord_ptr[], const CTYPE1 interior_vcoord[],
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE2> &
   poly_tri_result)
  {
    std::bitset<BIT_SET_SIZE> flag_not_ear;

    compute_triangulation_to_max_min_angleP
    (dimension, num_poly_vert, vcoord_ptr, interior_vcoord, flag_not_ear,
     max_small_magnitude, triangulation_settings, poly_tri_result);
  }

  
  /*!
   *  @brief Compute triangulation to max min angle.
   *  - Version using list of pointers to vertex coordinates.
   *  - Version using C++ STL vector poly_vcoord_ptr[].
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, typename CTYPEV, typename CTYPEI,
            typename NTYPE, typename BIT_SET_TYPE, typename MTYPE,
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_to_max_min_angleP
  (const DTYPE dimension,
   const std::vector<CTYPEV> & poly_vcoord_ptr,
   const CTYPEI interior_vcoord[],
   const BIT_SET_TYPE & flag_not_ear,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE> &
   poly_tri_result)
  {
    compute_triangulation_to_max_min_angleP
      (dimension, poly_vcoord_ptr.size(),
       IJK::vector2pointer(poly_vcoord_ptr),
       interior_vcoord, flag_not_ear, max_small_magnitude,
       triangulation_settings, poly_tri_result);
  }


  /*!
   *  @brief Compute triangulation to max min angle.
   *  - Version using list of pointers to vertex coordinates.
   *  - Version using C++ STL vector poly_vcoord_ptr[].
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, typename CTYPEV, typename CTYPEI,
            typename NTYPE, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE>
  void compute_triangulation_to_max_min_angleP
  (const DTYPE dimension, 
   const std::vector<CTYPEV> & poly_vcoord_ptr, 
   const CTYPEI interior_vcoord[], 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   poly_tri_result)
  {
    compute_triangulation_to_max_min_angleP
      (dimension, poly_vcoord_ptr.size(), 
       IJK::vector2pointer(poly_vcoord_ptr),
       interior_vcoord, max_small_magnitude, triangulation_settings,
       poly_tri_result);
  }


  /*!
   *  @brief Compute cosine of the min angle in the polygon 
   *    triangulation that maximizes the min angle.
   *  - Version using list of pointers to vertex coordinates.
   *  - Version trying two interior coordinates.
   *  @param triangulation_settings
   *    Indicates types of triangulations to consider.
   */
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEI0, typename CTYPEI1,
            typename NTYPE1, typename NTYPE2, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE,
            typename ITYPES>
  void compute_cos_max_min_polygon_triangulation_angleP_IVx2
  (const DTYPE dimension, const NTYPE1 num_poly_vert, 
   const CTYPEV * const vcoord_ptr[], 
   const CTYPEI0 interior_vcoord0[], 
   const CTYPEI1 interior_vcoord1[], 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE2> & 
   poly_tri_result,
   ITYPES & index_selected_IV)
  {
    const NTYPE1 NUM_VERT_IN_HEXAGON(6);
    COS_TYPE cos_min_angle;
    bool flag_zero;
    

    poly_tri_result.Clear();

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */

    if (triangulation_settings.flag_all_triangulations_from_polygon_vertices ||
        (triangulation_settings.flag_all_triangulations_from_hexagon_vertices &&
         num_poly_vert == NUM_VERT_IN_HEXAGON)) {
      compute_triangulation_over_all_to_max_min_angleP
        (dimension, num_poly_vert, vcoord_ptr, 
         max_small_magnitude, poly_tri_result);
    }
    else {
      compute_fan_triangulation_to_max_min_angleP
        (dimension, num_poly_vert, vcoord_ptr,
         max_small_magnitude, poly_tri_result);
    }

    if (triangulation_settings.flag_triangulate_from_interior_vertex) {

      compute_triangulation_from_vcoordx2_to_max_min_angleP
        (dimension, num_poly_vert, vcoord_ptr, 
         interior_vcoord0, interior_vcoord1, max_small_magnitude,
         cos_min_angle, index_selected_IV, flag_zero);
         
      if (!flag_zero) {
        if (poly_tri_result.cos_min_triangulation_angle > cos_min_angle) {
          poly_tri_result.SetInterior(cos_min_angle, num_poly_vert, 1);
        }
      }
    }

    if (triangulation_settings.flag_split_ear_triangulate_from_interior_vertex) {
      NTYPE1 iloc_min, iloc_max;
      COS_TYPE cos_min, cos_max;
      NTYPE1 num_angle;

      IJK::compute_cos_min_max_polygon_angles
        (dimension, num_poly_vert, vcoord_ptr, max_small_magnitude,
         cos_min, cos_max, iloc_min, iloc_max, num_angle);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "  num_angle: " << num_angle << endl;
      cerr << "  iloc_min: " << iloc_min << endl;
      cerr << "  num_poly_vert: " << num_poly_vert << endl;
      IJK::print_coord3D(cerr, "  vcoord[0]: ",
                         vcoord_ptr[0], "\n");
      IJK::print_coord3D(cerr, "  vcoord[1]: ",
                         vcoord_ptr[1], "\n");
      IJK::print_coord3D(cerr, "  coord of min angle: ",
                         vcoord_ptr[iloc_min], "\n");
      */
      
      if (num_angle == num_poly_vert) {

        ITYPES index_selectedB;
        compute_cos_min_triangulation_angle_using_cut_ear_vcoordXPx2
          (dimension, num_poly_vert, vcoord_ptr, 
           interior_vcoord0, interior_vcoord1,
           iloc_min, max_small_magnitude, 
           cos_min_angle, index_selectedB, flag_zero);

        if (!flag_zero) {
          if (poly_tri_result.cos_min_triangulation_angle > cos_min_angle) {
            poly_tri_result.SetInterior(cos_min_angle, num_poly_vert, 1);
            poly_tri_result.SetEarFlag(iloc_min, true);
            index_selected_IV = index_selectedB;
          }
        }
      }
    }

  }


  /// @brief Compute cosine of the min angle in the polygon triangulation 
  ///   that maximizes the min angle.
  /// - Version using list of pointers to vertex coordinates.
  /// - Version trying two interior coordinates[].
  /// @param triangulation_settings 
  ///   Indicates types of triangulations to consider.
  template <typename DTYPE, typename CTYPEV, 
            typename CTYPEI0, typename CTYPEI1,
            typename NTYPE1, typename NTYPE2, typename MTYPE, 
            typename COS_TYPE, int BIT_SET_SIZE,
            typename CTYPES>
  void compute_cos_max_min_polygon_triangulation_angleP_IVx2
  (const DTYPE dimension, const NTYPE1 num_poly_vert, 
   const CTYPEV * const vcoord_ptr[], 
   const CTYPEI0 interior_vcoord0[], 
   const CTYPEI1 interior_vcoord1[], 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE2> & 
   poly_tri_result,
   std::vector<CTYPES> & selected_interior_vcoord)
  {
    // Resize selected_interior_vcoord[] to dimension, just in case.
    selected_interior_vcoord.resize(dimension);

    NTYPE1 index_selected_IV;

    compute_cos_max_min_polygon_triangulation_angleP_IVx2
      (dimension, num_poly_vert, vcoord_ptr,
       interior_vcoord0, interior_vcoord1, max_small_magnitude,
       triangulation_settings, poly_tri_result, index_selected_IV);

    if (index_selected_IV == 0) {
      std::copy(interior_vcoord0, interior_vcoord0+dimension,
                selected_interior_vcoord.begin());
    }
    else {
      std::copy(interior_vcoord1, interior_vcoord1+dimension,
                selected_interior_vcoord.begin());
    }
  }


  // ***************************************************************
  //! @name COMPUTE COS MIN QUAD/TRIANGLE TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in the triangulation of a quadrilateral.
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

    compute_cos_min_triangle_angle
      (dimension, coord0, coord1, coord2, max_small_magnitude, cosA, flagA);
    compute_cos_min_triangle_angle
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
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in the triangulation of a quadrilateral.
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
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in the triangulation of a quadrilateral.
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

    compute_cos_min_triangle_angle
      (dimension, coord0, coord1, coordX, max_small_magnitude, 
       cos[0], flag_tri_zero[0]);
    compute_cos_min_triangle_angle
      (dimension, coord1, coord2, coordX, max_small_magnitude, 
       cos[1], flag_tri_zero[1]);
    compute_cos_min_triangle_angle
      (dimension, coord2, coord3, coordX, max_small_magnitude, 
       cos[2], flag_tri_zero[2]);
    compute_cos_min_triangle_angle
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
   *  @brief Compute the cosine of the smallest triangle angle 
   *    in the splitting of a triangle into two smaller triangles.
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
  //! @name COMPUTE COS MAX MIN QUAD TRIANGULATION ANGLE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
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
    select_min
      (cos_min_diag02, flag_zero_diag02, 
       cos_min_diag13, flag_zero_diag13, 
       cos_min_angle, index_selected, flag_zero);

    if (index_selected == 0) { flag_diag02 = true; }
    else { flag_diag02 = false; }

    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_coord3D(cerr, "    coord0: ", coord0, "\n");
    IJK::print_coord3D(cerr, "    coord1: ", coord1, "\n");
    IJK::print_coord3D(cerr, "    coord2: ", coord2, "\n");
    IJK::print_coord3D(cerr, "    coord3: ", coord3, "\n");
    cerr << "    cos_min_diag02: " << cos_min_diag02 << endl;
    cerr << "    cos_min_diag13: " << cos_min_diag13 << endl;
    */
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
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
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
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
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
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
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
   *  - *** DEPRECATED ***
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
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
    select_min
      (cos_min_diag02, flag_zero_diag02, 
       cos_min_diag13, flag_zero_diag13, 
       quad_tri.cos_min_triangulation_angle, index_selected, 
       quad_tri.flag_zero);

    quad_tri.num_triangles = 2;
    quad_tri.num_interior_tri_vertices = 0;
    quad_tri.tri_vertex_index = index_selected;
    quad_tri.triangulation_encoding.Clear();
    quad_tri.triangulation_encoding.Set(1-index_selected, true);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
   *  - *** DEPRECATED ***
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {
    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vcoord[0], vcoord[1], vcoord[2], vcoord[3], 
       max_small_magnitude, quad_tri);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
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
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals.
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri)
  {
    const CTYPE * vertex_coord_ptr = &(vertex_coord[0]);

    compute_cos_max_min_quad_tri02_or_tri13_angle
      (dimension, vertex_coord_ptr, quad_vert, max_small_magnitude, 
       quad_tri);
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle
   *    in the triangulation of a quadrilateral.
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
   *  @brief Compute the cosine of the smallest triangle angle
   *    in the triangulation of a quadrilateral.
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
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals or
   *    the triangulation into four triangles.
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

    // Prefer diag triangulation over tri4 triangulation.
    select_min
      (cos_min_diag, flag_zero_diag,
       cos_min_tri4, flag_zero_tri4, 
       cos_min_angle, index_selected, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr.precision(24);
    cerr << "index_selected: " << index_selected << endl;
    cerr << "    cos_min_tri4: " << cos_min_tri4 << endl;
    cerr << "    cos_min_diag: " << cos_min_diag << endl;
    cerr << "    flag_diag02: " << int(flag_diag02) << endl;
    */

    if (index_selected == 1) { flag_tri4 = true; }
    else { flag_tri4 = false; }
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals or
   *    the triangulation into four triangles.
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    bool flag_tri4, flag_diag02;

    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, coord0, coord1, coord2, coord3, vcoordX,
       max_small_magnitude, quad_tri_result.cos_min_triangulation_angle,
       flag_tri4, flag_diag02, quad_tri_result.flag_zero);

    quad_tri_result.triangulation_encoding.Clear();

    if (flag_diag02) {
      quad_tri_result.tri_vertex_index = 0; 
      quad_tri_result.triangulation_encoding.Set(1, true);
    }
    else {
      quad_tri_result.tri_vertex_index = 1; 
      quad_tri_result.triangulation_encoding.Set(0, true);
    }

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
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals or
   *    the triangulation into four triangles.
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
   *  - Version with array of coordinates.
   *  @param[out] quad_tri_result Quad triangulation result.
   */
  template <typename DTYPE, typename CTYPE, typename CTYPEX, 
            typename MTYPE, typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
  void compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
  (const DTYPE dimension, const CTYPE * quad_coord[4],
   const CTYPEX * vcoordX, const MTYPE max_small_magnitude, 
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   quad_tri_result)
  {
    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, quad_coord[0], quad_coord[1], 
       quad_coord[2], quad_coord[3], vcoordX, max_small_magnitude,
       quad_tri_result);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals or
   *    the triangulation into four triangles.
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
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
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals or
   *    the triangulation into four triangles.
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
  {
    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, IJK::vector2pointer(vertex_coord), vcoordX, quad_vert, 
       max_small_magnitude, quad_tri_result);
  }


  /*!
   *  @brief Compute the cosine of the max min quad triangulation angle
   *    over the two triangulations using the quad diagonals or
   *    the triangulation into four triangles.
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & quad_tri_result)
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

      compute_cos_min_triangle_angle
        (dimension, vcoord[i0], vcoord[i1], vcoord[i2], max_small_magnitude, 
         cos_ear_triangle[i1], flag_zero_ear[i1]);

      compute_cos_min_triangle_angle
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
   *  - Version using C++ STL vector vertex_coord[].
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
   *  @brief Compute min angle in triangulation of a pentagon
   *    either with 3 triangles incident on a pentagon vertex
   *    or with 5 triangles incident on vcoordX[].
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
   *  @brief Compute min angle in triangulation of a pentagon
   *    either with 3 triangles incident on a pentagon vertex
   *    or with 5 triangles incident on vcoordX[].
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   pentagon_tri_result)
  {
    const NTYPE NUM_PENTAGON_VERTICES(5);
    bool flag_tri5;

    compute_cos_max_min_pentagon_tri3_or_tri5_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, vcoordX, max_small_magnitude,
       pentagon_tri_result.cos_min_triangulation_angle, flag_tri5,
       pentagon_tri_result.tri_vertex_index, pentagon_tri_result.flag_zero);

    pentagon_tri_result.triangulation_encoding.Clear();

    if (flag_tri5) {
      pentagon_tri_result.num_triangles = 5;
      pentagon_tri_result.num_interior_tri_vertices = 1;
    }
    else {
      pentagon_tri_result.num_triangles = 3;
      pentagon_tri_result.num_interior_tri_vertices = 0;
      pentagon_tri_result.triangulation_encoding.SetFan
        (pentagon_tri_result.tri_vertex_index, NUM_PENTAGON_VERTICES);
    }

  }


  /*! @brief Compute min angle in triangulation of a pentagon
   *    either with 3 triangles incident on a pentagon vertex
   *    or with 5 triangles incident on vcoordX[].
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & pentagon_tri_result)
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


  /*! @brief Compute min angle in triangulation of a pentagon
   *    either with 3 triangles incident on a pentagon vertex
   *    or with 5 triangles incident on vcoordX[].
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
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
  ///   with all triangles incident on vertex_coord0.
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
   *  @brief Compute min angle in triangulation of a hexagon
   *    either with 4 triangles incident on a hexagon vertex
   *    or with 6 triangles incident on vcoordX[].
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
    compute_cos_min_hexagon_tri6_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2, vertex_coord3,
       vertex_coord4, vertex_coord5, vcoordX,
       max_small_magnitude, cos_min_tri6_angle, flag_zero_tri6);


    // Initialize
    num_triangles = 4;

    NTYPE index;
    select_min(cos_max_min_tri4_angle, flag_zero_tri4,
               cos_min_tri6_angle, flag_zero_tri6,
               cos_max_min_angle, index, flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    IJK::print_coord3D(cerr, "  vcoordX: ", vcoordX, "\n");
    cerr << "  cos_min_tri6_angle: " << cos_min_tri6_angle << endl;
    cerr << "  cos_max_min_tri4_angle: " << cos_max_min_tri4_angle << endl;
    cerr << "  iv_max_min: " << iv_max_min << endl;
    */


    if (index == 1) {
      num_triangles = 6;

      // If triangulation from vcoordX[], set iv_max_min to 6.
      iv_max_min = 6;
    }
  }


  /*!
   *  @brief Compute min angle in triangulation of a hexagon
   *    either with 4 triangles incident on a hexagon vertex
   *    or with 6 triangles incident on vcoordX[].
   *  - Version returning data structure POLYGON_TRIANGULATION_RESULT.
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
   POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   hexagon_tri_result)
  {
    const NTYPE NUM_HEX_VERTICES(6);

    compute_cos_max_min_hexagon_tri4_or_tri6_angle
      (dimension, vertex_coord0, vertex_coord1, vertex_coord2,
       vertex_coord3, vertex_coord4, vertex_coord5, vcoordX,
       max_small_magnitude, hexagon_tri_result.cos_min_triangulation_angle,
       hexagon_tri_result.num_triangles, hexagon_tri_result.tri_vertex_index,
       hexagon_tri_result.flag_zero);

    hexagon_tri_result.triangulation_encoding.Clear();
    if (hexagon_tri_result.num_triangles == 6) {
      hexagon_tri_result.num_interior_tri_vertices = 1;
    }
    else {
      hexagon_tri_result.num_interior_tri_vertices = 0;
      hexagon_tri_result.triangulation_encoding.SetFan
        (hexagon_tri_result.tri_vertex_index, NUM_HEX_VERTICES);
    }

  }

  ///@}


  // ***************************************************************
  //! @name TRIANGULATE POLYGON
  // ***************************************************************

  ///@{

  /*!
   *  @brief Triangulate a single polygon to maximize the minimum triangle angle.
   *  - All triangulation triangles are incident on one vertex,
   *    i.e. all triangulations are fan triangulations.
   *  - Add new triangles to vector tri_vert.
   *  @param[out] index_tri_vert Index of triangulation vertex.
   *    - All triangulation triangles are incident on poly_vert[index_tri_vert].
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

    compute_fan_triangulation_to_max_min_angle
      (dimension, vertex_coord, num_poly_vert, poly_vert,
       max_small_magnitude, index_tri_vert, cos_min_angle, flagA);

    if (flagA) {
      index_tri_vert = 0;
      cos_min_angle = 0;
    }

    triangulate_polygon_from_vertex
      (num_poly_vert, poly_vert, index_tri_vert, tri_vert);
  }


  /*!
   *  @brief Triangulate a single polygon to maximize the minimum triangle angle.
   *  - All triangulation triangles are incident on one vertex.
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
   *  @brief Triangulate a single polygon to maximize the minimum triangle angle.
   *  - All triangulation triangles are incident on one vertex.
   *  - Version that does not return index_tri_vert or cos_min_angle.
   *  - Version using C++ STL vectors vertex_coord[] and poly_vert[].
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

  ///@}


  // ***************************************************************
  //! @name TRIANGULATE QUADRILATERAL
  // ***************************************************************

  ///@{

  /*!
   *  @brief Triangulate a single quadrilateral to maximize the minimum triangle angle.
   *  - Triangulate using a quadrilateral diagonal, splitting the quadrilateral
   *    into two triangles.
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename MTYPE>
  void triangulate_quad_tri2_max_min_angle
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


    /*!
   *  @brief Triangulate each quadrilateral in list to maximize 
   *    the minimum triangle angle.
   *  - Triangulate using a quadrilateral diagonal, splitting 
   *    the quadrilateral into two triangles.
   *  - For each quad, choose the triangulation that maximizes 
   *    the minimum triangle angle.
   *  - Skipped quadrilaterals whose vertices are all identical.
   */
  template <typename DTYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE, typename NTYPE_TRI>
  void triangulate_quad_tri2_max_min_angle_skip_collapsed
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const NTYPEQ num_quad,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad)
  {
    const int NUM_VERT_PER_QUAD(4);

    num_triangulated_quad = 0;
    
    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      const NTYPEQ k = iquad*NUM_VERT_PER_QUAD;
      if (!is_quad_collapsed_to_one_vertex(quad_vert+k)) {
        triangulate_quad_tri2_max_min_angle
          (dimension, vertex_coord, quad_vert+k, max_small_magnitude,
           tri_vert);
        num_triangulated_quad++;
      }
    }
  }
  
    
  /*!
   *  @brief Triangulate each quadrilateral in list to maximize 
   *    the minimum triangle angle.
   *  - Triangulate using a quadrilateral diagonal, splitting 
   *    the quadrilateral into two triangles.
   *  - For each quad, choose the triangulation that maximizes 
   *    the minimum triangle angle.
   *  @param flag_skip_collapsed_quads If true, skip any quads
   *    that have been collapsed to a single vertex.
   *    - (Still triangulates quads that are collapsed to an edge
   *       or two edges.)
   *    - Default: true.
   */
  template <typename DTYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE, typename NTYPE_TRI>
  void triangulate_quad_tri2_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const NTYPEQ num_quad,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad,
   const bool flag_skip_collapsed_quads=true)
  {
    const int NUM_VERT_PER_QUAD(4);

    num_triangulated_quad = 0;
    
    if (flag_skip_collapsed_quads) {
      triangulate_quad_tri2_max_min_angle_skip_collapsed
        (dimension, vertex_coord, quad_vert, num_quad,
         max_small_magnitude, tri_vert, num_triangulated_quad);
    }
    else {
      for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
        const NTYPEQ k = iquad*NUM_VERT_PER_QUAD;
        triangulate_quad_tri2_max_min_angle
          (dimension, vertex_coord, quad_vert+k, max_small_magnitude,
           tri_vert);
        num_triangulated_quad++;
      }
    }
  }


  /*!
   *  @overload
   *  @brief Triangulate each quadrilateral in list to maximize 
   *    the minimum triangle angle.
   *  - For each quad, choose the triangulation that maximizes 
   *    the minimum triangle angle.
   *  - Version using C++ STL vectors vertex_coord[] and quad_vert[].
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1,
            typename MTYPE, typename NTYPE_TRI>
  void triangulate_quad_tri2_max_min_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE0> & quad_vert,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE1> & tri_vert,
   NTYPE_TRI & num_triangulated_quad,
   const bool flag_skip_collapsed_quads=true)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    triangulate_quad_tri2_max_min_angle
      (dimension, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert), num_quad, max_small_magnitude,
       tri_vert, num_triangulated_quad, flag_skip_collapsed_quads);
  }


  /*!
   *  @brief Return true if axis parallel plane separates iv0 and iv1 
   *    from iv2 and iv3.
   *  @param point_coord[] Coordinates of point on plane.
   *  @param normal_direction Axis direction (0, 1, or 2) normal to plane.
   *  @param vertex_coord[] Array of vertex coordinates.
   *    - Vertex i has coordinates:
   *      (vertex_coord[3*i], vertex_coord[3*i+1], vertex_coord[3*i+2])
   */
  template <typename CTYPEV, typename CTYPEP,
	    typename VTYPE, typename DIR_TYPE>
  bool does_axis_parallel_plane_separate_3D
  (const CTYPEP point_coord[],
   const DIR_TYPE normal_direction,
   const CTYPEV vertex_coord[],
   const VTYPE iv0, const VTYPE iv1,
   const VTYPE iv2, const VTYPE iv3)
  {
    const DIR_TYPE DIM3(3);
    const CTYPEV wc = point_coord[normal_direction];
    const CTYPEV vc0 = vertex_coord[DIM3*iv0+normal_direction];
    const CTYPEV vc1 = vertex_coord[DIM3*iv1+normal_direction];
    const CTYPEV vc2 = vertex_coord[DIM3*iv2+normal_direction];
    const CTYPEV vc3 = vertex_coord[DIM3*iv3+normal_direction];

    if (vc0 <= wc && vc1 <= wc && vc2 >= wc && vc3 >= wc)
      { return(true); }

    if (vc0 >= wc && vc1 >= wc && vc2 <= wc && vc3 <= wc)
      { return(true); }

    return(false);
  }


  /*!
   *  @brief Return true if quad vertices surround the line.
   *  @param point_coord[] Coordinates of point on line.
   *  @param line_direction Line direction: 0, 1, or 2.
   */
  template <typename CTYPEP, typename CTYPEV, typename DIR_TYPE,
	    typename VTYPE>
  bool does_quad_surround_line
  (const CTYPEP point_coord[],
   const DIR_TYPE line_direction,
   const CTYPEV vertex_coord[],
   const VTYPE quad_vert[])
  {
    const int DIM3(3);

    const DIR_TYPE d1 = (line_direction+1)%DIM3;
    const DIR_TYPE d2 = (line_direction+2)%DIM3;

    if (does_axis_parallel_plane_separate_3D
        (point_coord, d1, vertex_coord,
         quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])) {

      if (does_axis_parallel_plane_separate_3D
          (point_coord, d2, vertex_coord,
           quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[0])) {

        return(true);
      }
    }

    if (does_axis_parallel_plane_separate_3D
        (point_coord, d2, vertex_coord,
         quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3])) {

      if (does_axis_parallel_plane_separate_3D
          (point_coord, d1, vertex_coord,
           quad_vert[1], quad_vert[2], quad_vert[3], quad_vert[0])) {

        return(true);
      }
    }

    return(false);
  }


  /// @brief Return true if quad vertices surround the dual edge.
  template <typename CTYPE, typename VTYPE0, typename VTYPE1,
            typename GRID_EDGE_TYPE>
  bool does_quad_surround_dual_edge
  (const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const GRID_EDGE_TYPE & dual_edge,
   const VTYPE1 vertex_on_dual_edge)
  {
    typedef typename GRID_EDGE_TYPE::DIRECTION_TYPE DIR_TYPE;

    const DIR_TYPE DIM3(3);
    const DIR_TYPE edge_dir = dual_edge.GridEdgeDirection();
    const CTYPE * point_coord =
      vertex_coord+(vertex_on_dual_edge*DIM3);

    return does_quad_surround_line
      (point_coord, edge_dir, vertex_coord, quad_vert);
  }


  /*!
   *  @brief Triangulate a single quad into two or four triangles.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   *  @param ivert_in_quad_interior Index of isosurface vertex in quad "interior".
   *    - If quad is triangulated into four triangles, then each triangle
   *      is incident on vertex ivert_in_quad_interior.
   *    - If quad is triangulated into two triangles, then vertex
   *      ivert_in_quad_interior is ignored.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  void triangulate_quad_tri2_or_tri4_max_min_angle
  (const DTYPE dimension,
   const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const VTYPE1 ivert_in_quad_interior,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert,
   bool & flag_tri4)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);
    const CTYPE * vcoord0 = vertex_coord+quad_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+quad_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+quad_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+quad_vert[3]*dimension;
    const CTYPE * vcoordX = vertex_coord+ivert_in_quad_interior*dimension;
    CTYPE cos_min_angle;
    bool flag_diag02, flag_zero;


    // *** DEBUG ***
    /*
    using namespace std;
    IJK::print_quad_vertices(cerr, "Triangulating quad: ", quad_vert, "\n");
    IJK::print_coord3D(cerr, "  vcoord0: ", vcoord0, "\n");
    IJK::print_coord3D(cerr, "  vcoord1: ", vcoord1, "\n");
    IJK::print_coord3D(cerr, "  vcoord2: ", vcoord2, "\n");
    IJK::print_coord3D(cerr, "  vcoord3: ", vcoord3, "\n");
    IJK::print_coord3D(cerr, "  vcoordX: ", vcoordX, "\n");
    */

    compute_cos_max_min_quad_tri02_or_tri13_or_tri4_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoordX,
       max_small_magnitude, cos_min_angle,
       flag_tri4, flag_diag02, flag_zero);

    if (flag_tri4) {
      triangulate_polygon_from_interior_vertex
        (NUM_VERT_PER_QUAD, quad_vert, ivert_in_quad_interior, tri_vert);
    }
    else {
      triangulate_quad_using_diagonal(quad_vert, flag_diag02, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a single quad into four triangles 
   *    or two triangles using diagonal (v0,v2).
   *  - Return number of triangles (2 or 4) in the triangulation.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  int triangulate_quad_tri4_or_tri02_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivert_in_quad_interior,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    const int TWO_TRIANGLES(2);
    const int FOUR_TRIANGLES(4);
    const CTYPE * vcoord0 = vertex_coord+quad_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+quad_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+quad_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+quad_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord+ivert_in_quad_interior*dimension;
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
      if (cos_min_diag02 <= cos_min_tri4)
        { flag_tri4 = false; }
      else
        { flag_tri4 = true; }
    }

    if (flag_tri4) {
      triangulate_quad_from_interior_vertex
        (quad_vert, ivert_in_quad_interior, tri_vert);

      return FOUR_TRIANGLES;
    }
    else {
      triangulate_quad_from_vertex(quad_vert, 0, tri_vert);

      return TWO_TRIANGLES;
    }
  }


  /*!
   *  @brief Triangulate a single quad into four triangles 
   *    or two triangles using diagonal (v0,v2).
   *  - Return number of triangles (2 or 4) in the triangulation.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Version using C++ STL vector vertex_coord[].
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  int triangulate_quad_tri4_or_tri02_max_min_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivert_in_quad_interior,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    return triangulate_quad_tri4_or_tri02_max_min_angle
	    (dimension, IJK::vector2pointer(vertex_coord), quad_vert,
	     ivert_in_quad_interior, max_small_magnitude, tri_vert);
  }


  /*!
   *  @brief Triangulate a single quad into four triangles 
   *    or two triangles using diagonal (v1,v3).
   *  - Return number of triangles (2 or 4) in the triangulation.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Quad vertices are listed in clockwise or counter-clockwise order
   *    around the quadrilateral.
   *  - Add new triangles to vector tri_vert.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  int triangulate_quad_tri4_or_tri13_max_min_angle
  (const DTYPE dimension,
   const CTYPE * vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivert_in_quad_interior,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    const int TWO_TRIANGLES(2);
    const int FOUR_TRIANGLES(4);
    const CTYPE * vcoord0 = vertex_coord+quad_vert[0]*dimension;
    const CTYPE * vcoord1 = vertex_coord+quad_vert[1]*dimension;
    const CTYPE * vcoord2 = vertex_coord+quad_vert[2]*dimension;
    const CTYPE * vcoord3 = vertex_coord+quad_vert[3]*dimension;
    const CTYPE * vcoord4 = vertex_coord+ivert_in_quad_interior*dimension;
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
      triangulate_quad_from_interior_vertex
        (quad_vert, ivert_in_quad_interior, tri_vert);

      return FOUR_TRIANGLES;
    }
    else {
      triangulate_quad_from_vertex(quad_vert, 1, tri_vert);

      return TWO_TRIANGLES;
    }
  }


  /*!
   *  @brief Triangulate a single quad into four triangles 
   *    or two triangles using diagonal (v1,v3).
   *  - Return number of triangles (2 or 4) in the triangulation.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Version using C++ STL vector vertex_coord[].
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE>
  int triangulate_quad_tri4_or_tri13_max_min_angle
  (const DTYPE dimension,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE0 * quad_vert,
   const VTYPE1 ivert_in_quad_interior,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert)
  {
    return triangulate_quad_tri4_or_tri13_max_min_angle
	    (dimension, IJK::vector2pointer(vertex_coord), quad_vert,
	     ivert_in_quad_interior, max_small_magnitude, tri_vert);
  }

  ///@}


  // ***************************************************************
  //! @name TRIANGULATE A LIST OF QUADRILATERALS
  // ***************************************************************

  ///@{


  /*!
   *  @brief Triangulate a list of isosurface quadrilaterals, 
   *    splitting each into two or four triangles.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Skipped quadrilaterals whose vertices are all identical.
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE, typename NTYPE_TRI2, typename NTYPE_TRI4>
  void triangulate_quad_tri2_or_tri4_max_min_angle_skip_collapsed
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const NTYPEQ num_quad,
   const VTYPE1 index_of_first_interior_vertex,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert,
   NTYPE_TRI2 & num_tri2,
   NTYPE_TRI4 & num_tri4)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    const int NUM_VERT_PER_QUAD(4);
    IJK::PROCEDURE_ERROR error
      ("triangulate_quad_tri2_or_tri4_max_min_angle_skip_collaped");

    if (num_quad > 0 && quad_vert == NULL) {
      error.AddMessage
        ("Programming error.  Array quad_vert[] is empty.");
      throw error;
    }

    num_tri2 = 0;
    num_tri4 = 0;
    
    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      const VTYPE2 iw = index_of_first_interior_vertex+iquad;
      const NTYPEQ k = iquad*NUM_VERT_PER_QUAD;
      bool flag_tri4;

      if (!is_quad_collapsed_to_one_vertex(quad_vert+k)) {
        // Split into two or four triangles, whichever maximizes min angle.
        triangulate_quad_tri2_or_tri4_max_min_angle
          (dimension, vertex_coord, quad_vert+k, iw, max_small_magnitude,
           tri_vert, flag_tri4);

        if (flag_tri4) { num_tri4++; }
        else { num_tri2++; }
      }
    }
  }

  
  /*!
   *  @brief Triangulate a list of isosurface quadrilaterals, 
   *    splitting each into two or four triangles.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   */
  
  template <typename GRID_TYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE,
            typename NTYPE_TRI2, typename NTYPE_TRI4>
  void triangulate_quad_tri2_or_tri4_max_min_angle
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   const VTYPE0 quad_vert[],
   const NTYPEQ num_quad,
   const VTYPE1 index_of_first_interior_vertex,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert,
   NTYPE_TRI2 & num_tri2,
   NTYPE_TRI4 & num_tri4,
   const bool flag_skip_collapsed=true)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    const int NUM_VERT_PER_QUAD(4);
    IJK::PROCEDURE_ERROR error("triangulate_quad_tri2_or_tri4_max_min_angle");

    if (num_quad > 0 && quad_vert == NULL) {
      error.AddMessage
        ("Programming error.  Array quad_vert[] is empty.");
      throw error;
    }

    if (flag_skip_collapsed) {
      triangulate_quad_tri2_or_tri4_max_min_angle_skip_collapsed
        (grid, vertex_coord, quad_vert, num_quad,
         index_of_first_interior_vertex, max_small_magnitude,
         tri_vert, num_tri2, num_tri4);
    }
    else {
      for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
        const VTYPE2 iw = index_of_first_interior_vertex+iquad;
        const NTYPEQ k = iquad*NUM_VERT_PER_QUAD;
        bool flag_tri4;

        // Split into two or four triangles, whichever maximizes min angle.
        triangulate_quad_tri2_or_tri4_max_min_angle
          (dimension, vertex_coord, quad_vert+k, iw, max_small_magnitude,
           tri_vert, flag_tri4);
        
        if (flag_tri4) { num_tri4++; }
        else { num_tri2++; }
      }
    }
  }


  /*!
   *  @brief Triangulate a list of isosurface quadrilaterals.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Version using C++ STL vectors vertex_coord[] and quad_vert[].
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE,
            typename NTYPE_TRI2, typename NTYPE_TRI4>
  void triangulate_quad_tri2_or_tri4_max_min_angle
  (const GRID_TYPE & grid,
   const std::vector<CTYPE> & vertex_coord,
   const std::vector<VTYPE0> & quad_vert,
   const VTYPE1 first_additional_vertex,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE2> & tri_vert,
   NTYPE_TRI2 & num_tri2,
   NTYPE_TRI4 & num_tri4,
   const bool flag_skip_collapsed=true)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    triangulate_quad_tri2_or_tri4_max_min_angle
      (grid, IJK::vector2pointer(vertex_coord),
       IJK::vector2pointer(quad_vert), num_quad,
       first_additional_vertex, max_small_magnitude,
       tri_vert, num_tri2, num_tri4, flag_skip_collapsed);
  }

  ///@}


  // ***************************************************************
  //! @name TRIANGULATE PENTAGON
  // ***************************************************************

  ///@{

  // *** SHOULD RENAME AS triangulate_pentagon_with_interior_vertex
  //     and move to ijktriangulate_poly2D.tpp.
  // *** SHOULD CREATE VERSION WITHOUT flag_reverse_orient ***
  /*!
   *  @brief Triangulate a single pentagon into five triangles.
   *  - Add new coord vcoordX to vertex_coord[].
   *  - Pentagon vertices are listed in clockwise or counter-clockwise order
   *    around the pentagon.
   *  - Add new triangles to vector tri_vert.
   *  @param vcoordX[] Coordinates of interior vertex used in triangulation.
   *    - Every triangle in the triangulation is incident on vcoordX[].
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
    triangulate_pentagon_with_vertex
      (v0, v1, v2, v3, v4, iv_new, flag_reverse_orient, tri_vert);
  }


  /*!
   *  @brief Triangulate a single pentagon into five triangles.
   *  - Version where pentagon vertices are in an array.
   *  @param vcoordX[] Coordinates of interior vertex used in triangulation.
   *    - Every triangle in the triangulation is incident on vcoordX[].
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
   *  @param vcoordX[] Coordinates of vertex used in triangulation into five triangles.
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
    const int NUM_VERT_PER_PENTAGON(5);
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
   *  @brief Triangulate a single pentagon into three or five triangles based on poly_tri.
   *  - Version using data structure POLYGON_TRIANGULATION_RESULT.
   *  @param vcoordX[] Coordinates of vertex used in triangulation into five triangles.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3, typename VTYPE4,
            typename VTYPE5, typename CTYPE0, typename CTYPE1,
            typename COS_TYPE, typename NTYPE,
            int BIT_SET_SIZE>
 void triangulate_pentagon_tri3_or_tri5
  (const DTYPE dimension,
   const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const VTYPE3 v3, const VTYPE4 v4,
   const CTYPE0 vcoordX[],
   const bool flag_reverse_orient,
   const POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
   poly_tri,
   std::vector<CTYPE1> & vertex_coord,
   std::vector<VTYPE5> & tri_vert)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    const VTYPE0 pentagon_vert[NUM_VERT_PER_PENTAGON] = 
      { v0, v1, v2, v3, v4 };

    if (poly_tri.num_triangles == 5) {
      triangulate_pentagon_tri5
        (dimension, v0, v1, v2, v3, v4, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      triangulate_pentagon
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
   *  @param vcoordX[] Coordinates of vertex used in triangulation into five triangles.
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
      triangulate_pentagon
        (v0, v1, v2, v3, v4,flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a single pentagon into three or five triangles based on poly_tri.
   *  - Triangulate from v0 in triangulation into 3 triangles.
   *  - Version using data structure POLYGON_TRIANGULATION_RESULT.
   *  @param vcoordX[] Coordinates of vertex used in triangulation into five triangles.
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
   const POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
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
      triangulate_pentagon
        (v0, v1, v2, v3, v4, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a pentagon into three triangles.
   *  - Maximize minimum triangle angle.
   *  - Consider all five possible triangulations.
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
      triangulate_pentagon(pentagon_vert, 0, tri_vert);
    }
    else {
      triangulate_pentagon(pentagon_vert, iv, tri_vert);
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
      triangulate_pentagon(pentagon_vert, 0, flag_reverse_orient, tri_vert);
    }
    else if (flag_tri5) {
      triangulate_pentagon_tri5
        (dimension, pentagon_vert, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      triangulate_pentagon
        (pentagon_vert, iv_max_min, flag_reverse_orient, tri_vert);
    }

  }


  /*!
   *  @brief Triangulate a pentagon into three or five triangles.
   *  - Maximize minimum triangle angle.
   *  - Version using C++ STL vector vertex_coord[].
   */
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


  /*!
   *  @brief Triangulate each pentagonin list into three triangles.
   *  - Maximize minimum triangle angle.
   *  - For each pentagon, consider five different possible triangulations.
   */
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


  /*!
   *  @brief Triangulate each pentagon in list into three triangles.
   *  - Maximize minimum triangle angle.
   *  - Version using C++ STL vector pentagon_vert[].
   *  - For each pentagon, consider five different possible triangulations.
   */
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
   *  @brief Triangulate each pentagon in list into three or five triangles.
   *  - For triangulations into five triangles, add vertex at centroid
   *    of pentagon.
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


  /*!
   *  @brief Triangulate each pentagon in list into three or five triangles.
   *  - Version using C++ STL vector pentagon_vert[].
   *  - For triangulations into five triangles, add vertex at centroid
   *    of pentagon.
   */
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
  //! @name TRIANGULATE HEXAGON
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
            typename MTYPE, int BIT_SET_SIZE>
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
    POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,CTYPE,int> quadA_tri_result;
    POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,CTYPE,int> quadB_tri_result;
    POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,CTYPE,int> quadC_tri_result;
    POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,CTYPE,int> quadD_tri_result;
    CTYPE cos_min_hex;
    bool flag_zero_hex;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  hex_vert: ";
    IJK::print_list(cerr, hex_vert, NUM_VERT_PER_HEX);
    cerr << endl;
    */

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
       

    compute_cos_min_hexagon_tri6_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3, vcoord4, vcoord5, 
       vcoordX, max_small_magnitude, cos_min_hex, flag_zero_hex);

    int index;
    CTYPE cos_min_angle;
    bool flag_zero;
    select_minIII
      (cos_max_min_quadsAB, flag_zeroAB,
       cos_max_min_quadsCD, flag_zeroCD,
       cos_min_hex, flag_zero_hex,
       cos_min_angle, index, flag_zero);
    
    if (index == 0) {
      // Triangulate quadA and quadB separately using quad diagonals.
      triangulate_quad
        (quadA_vert, quadA_tri_result.tri_vertex_index, 
         flag_reverse_orient, tri_vert);
      triangulate_quad
        (quadB_vert, quadB_tri_result.tri_vertex_index, 
         flag_reverse_orient, tri_vert);
    }
    else if (index == 1) {
      // Triangulate quadC and quadD separately using quad diagonals.
      triangulate_quad
        (quadC_vert, quadC_tri_result.tri_vertex_index, 
         flag_reverse_orient, tri_vert);
      triangulate_quad
        (quadD_vert, quadD_tri_result.tri_vertex_index, 
         flag_reverse_orient, tri_vert);
    }
    else if (index == 2) {
      // Use triangulation into 6 triangles incident on vcoordX.
      triangulate_polygon_with_vertex
        (NUM_VERT_PER_HEX, hex_vert, ivertX, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate hexagon.
   *  - Maximize minimum triangle angle.
   *  - Version using C++ STL vector vertex_coord[].
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
    triangulate_hexagon_with_vertex
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
    const int NUM_VERT_PER_HEXAGON(6);
    const TRI_VTYPE hexagon_vert[NUM_VERT_PER_HEXAGON] = 
      { v0, v1, v2, v3, v4, v5 };

    if (num_triangles == 6) {
      triangulate_hexagon_tri6
        (dimension, v0, v1, v2, v3, v4, v5, vcoordX, flag_reverse_orient, 
         vertex_coord, tri_vert);
    }
    else {
      triangulate_hexagon
        (hexagon_vert, iv_index, flag_reverse_orient, tri_vert);
    }
  }


  /*!
   *  @brief Triangulate a single hexagon into four or six triangles.
   *  - Version using data structure POLYGON_TRIANGULATION_RESULT.
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
   const POLYGON_TRIANGULATION_RESULT<BIT_SET_SIZE,COS_TYPE,NTYPE> & 
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
  //! @name TRIANGULATE SEPTAGON
  // ***************************************************************

  ///@{

  /*!
   *  @brief Compute min angle in triangulation of septagon
   *    with all triangles incident on vertex_coord0.
   *  *** PROBABLY SHOULD BE "const CTYPE vcoord0[], ..."
   */
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


  // ***************************************************************
  //! @name TRIANGULATE QUADS WITH DIAGONALS OUTSIDE THE ENVELOPE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Triangulate isosurface quadrilateral if it has some diagonal outside the envelope.
   *  - Return number of triangles (0, 2 or 4) in triangulation 
   *    of isosurface quadrilateral. (0: no triangulation.)
   *  - If the isosurface quadrilateral is not triangulated
   *    (all diagonals are inside the envelope), then return 0.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Quadrilateral is dual to a grid edge.
   *  - Only use diagonals that are inside the envelope formed
   *    by the quad edges and the dual grid edge.
   *  - If no diagonals are inside the envelope, triangulate into four triangles.
   *  - Checks that the quad edges "surround" the dual grid edge.
   *  - If the quad edges do not "surround" the dual grid edge
   *    and at least one diagonal is within the envelope,
   *    then split into two triangles using a single quad diagonal.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   *  @param quad_vert[] Vertices of the quadrilateral.
   *  @param quad_info Quad information.
   *  @param ivert_in_quad_interior Index of vertex in quad interior.
   *    - Integer type.
   *    - If quad is triangulated into four triangles, then each triangle
   *      is incident on vertex ivert_in_quad_interior.
   *  @param min_dist_allow_tri4 Allow triangulation into 4 triangles when each
   *    quad vertex is min_distance_allow_tri4 distance from one
   *    of the facets incident on the dual quad edge.
   *    - If some quad vertex is closer than min_distance_allow_tri4 to
   *    both neighboring facets incident on the dual quad edge,
   *    use a diagonal triangulation into two triangles.
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPEQ, typename VTYPEI, typename VTYPET,
            typename ISOQUAD_INFO_TYPE, typename MTYPE>
  int tri2_or_tri4_dual_quad_DoutsideE_3D
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   const VTYPEQ quad_vert[],
   const ISOQUAD_INFO_TYPE & quad_info,
   const VTYPEI ivert_in_quad_interior,
   const MTYPE max_small_magnitude,
   std::vector<VTYPET> & tri_vert)
  {
    const int DIM3(3);
    const int ZERO_TRIANGLES(0);
    const int FOUR_TRIANGLES(4);

    if (quad_info.IsDiagonalInEnvelope(0)) {

      if (quad_info.IsDiagonalInEnvelope(1)) {
        // Both diagonals are in the envelope.
        // Don't triangulate.

        return ZERO_TRIANGLES;
      }
      else {
        // Only diagonal (v0,v2) is in the envelope.
        const int num_triangles =
          triangulate_quad_tri4_or_tri02_max_min_angle
          (DIM3, vertex_coord, quad_vert, ivert_in_quad_interior,
           max_small_magnitude, tri_vert);

        return num_triangles;
      }
    }
    else if (quad_info.IsDiagonalInEnvelope(1)) {
      // Only diagonal (v1,v3) is in the envelope.
      const int num_triangles =
        triangulate_quad_tri4_or_tri13_max_min_angle
        (DIM3, vertex_coord, quad_vert, ivert_in_quad_interior,
         max_small_magnitude, tri_vert);

      return num_triangles;
    }
    else {
      // Neither diagonal is in envelope.
      // Triangulate into four triangles using vertex ivert_in_quad_interior.
      triangulate_quad_from_interior_vertex
        (quad_vert, ivert_in_quad_interior, tri_vert);
      return FOUR_TRIANGLES;
    }

  }

    
  /*!
   *  @brief Triangulate isosurface quadrilaterals with some diagonal 
   *    outside the envelope.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Each quadrilateral is dual to a grid edge.
   *  - Only use diagonals that are inside the envelope formed
   *    by the quad edges and the dual grid edge.
   *  - If no diagonals are inside the envelope, triangulate into four triangles.
   *  - Checks that the quad edges "surround" the dual grid edge.
   *  - If the quad edges do not "surround" the dual grid edge
   *    and at least one diagonal is within the envelope,
   *    then split into two triangles using a single quad diagonal.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  - Collapse any triangulated quadrilaterals to a single vertex
   *    to mark their replacement by triangles.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   *  @param quad_info[i] Information on i'th isosurface quadrilateral.
   *  @param index_of_first_interior_vertex Index of vertex in interior of quad 0.
   *    - Integer type.
   *    - Index of vertex in interior of quad i is (index_of_first_interior_vertex+i).
   *    @pre Interior vertices have already been created and
   *    assigned coordinates in vertex_coord[].
   */
  template <typename GRID_TYPE, typename CTYPE, typename NTYPEQ,
            typename VTYPEQ, typename VTYPEI, typename VTYPET,
            typename ISOQUAD_INFO_TYPE, typename MTYPE,
            typename NTYPES2, typename NTYPES4>
  void tri2_or_tri4_dual_quad_list_with_diagonals_outside_envelope
  (const GRID_TYPE & grid,
   const CTYPE vertex_coord[],
   const VTYPEI index_of_first_interior_vertex,
   const MTYPE max_small_magnitude,
   VTYPEQ quad_vert[],
   const NTYPEQ num_quad,
   const ISOQUAD_INFO_TYPE quad_info[],
   std::vector<VTYPET> & tri_vert,
   NTYPES2 & num_tri2_split, NTYPES4 & num_tri4_split)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const int NUM_VERT_PER_QUAD(4);
    const int ZERO_TRIANGLES(0);
    const int TWO_TRIANGLES(2);
    const DTYPE dimension = grid.Dimension();
    const DTYPE DIM3(3);
    IJK::PROCEDURE_ERROR error
      ("tri2_or_tri4_dual_quad_list_with_diagonals_outside_envelope");

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

      const int num_triangles =
        tri2_or_tri4_dual_quad_DoutsideE_3D
        (grid, vertex_coord, quad_vert_iquad, quad_info[iquad],
         ivert_in_quad_interior, max_small_magnitude, tri_vert);

      if (num_triangles > ZERO_TRIANGLES) {

        collapse_quad_to_vertex(quad_vert_iquad);

        if (num_triangles == TWO_TRIANGLES)
          { num_tri2_split++; }
        else {
          // Isosurface quadrilateral was split into four triangles.
          num_tri4_split++;
        }
      }
    }
  }


  /*!
   *  @overload
   *  @brief Triangulate a list of isosurface quadrilaterals, (C++ vector.)
   *    splitting each into two or four triangles.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Only use diagonals that are inside the envelope formed
   *    by the quad edges and the dual grid edge.
   *  - Version using C++ STL vectors vertex_coord[] and quad_vert[].
   */
  template <typename GRID_TYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ISOQUAD_INFO_TYPE, typename MTYPE,
            typename NTYPES2, typename NTYPES4>
  void tri2_or_tri4_dual_quad_list_with_diagonals_outside_envelope
  (const GRID_TYPE & grid,
   const std::vector<CTYPE> & vertex_coord,
   const VTYPE1 index_of_first_interior_vertex,
   const MTYPE max_small_magnitude,
   std::vector<VTYPE0> & quad_vert,
   const std::vector<ISOQUAD_INFO_TYPE> & quad_info,
   std::vector<VTYPE2> & tri_vert,
   NTYPES2 & num_tri2_split, NTYPES4 & num_tri4_split)
  {
    typedef typename std::vector<VTYPE0>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_QUAD(4);
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (num_quad != quad_info.size()) {
      IJK::PROCEDURE_ERROR error
        ("tri2_or_tri4_dual_quad_list_with_diagonals_outside_envelope");
      error.AddMessage
        ("Programming error. Size of quad_info[] incorrect.");
      error.AddMessage
        ("  Size of quad_info[] does not equal number of isosurface quadrilaterals.");
      error.AddMessage("  Size of quad info: ", quad_info.size(), ".");
      error.AddMessage("  Number of isosurface quadrilaterals: ",
                       num_quad, ".");
      throw error;
    }

    tri2_or_tri4_dual_quad_list_with_diagonals_outside_envelope
      (grid, IJK::vector2pointer(vertex_coord),
       index_of_first_interior_vertex, max_small_magnitude,
       IJK::vector2pointer(quad_vert), num_quad,
       IJK::vector2pointer(quad_info),
       tri_vert, num_tri2_split, num_tri4_split);
  }

  ///@}

}

#endif
