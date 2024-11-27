/*!
 *  @file ijkcoord.tpp
 *  @brief ijk templates for coordinate arithmetic
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009-2024 Rephael Wenger

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

#ifndef _IJKCOORD_
#define _IJKCOORD_

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "ijk.tpp"


namespace IJK {

  // **************************************************
  //! @name Basic coordinate operations
  // **************************************************

  //@{

  /*!
   *  @brief Set all vertex coordinates to \a c.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param c Scalar constant.  Set vertex coordinates to \a c.
   *  @param[out] coord[] Output coordinates.
   */
  template <typename DTYPE, typename STYPE, typename CTYPE>
  void set_coord(const DTYPE dimension, const STYPE c, CTYPE coord[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord[d] = CTYPE(c); };
  }


  /*!
   *  @brief Copy \a coord0[] to \a coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param[out] coord1[] Output coordinates.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  void copy_coord(const DTYPE dimension, 
                  const CTYPE0 coord0[], CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = coord0[d]; };
  }


  /*!
   *  @brief Copy \a coord0[] to \a coord1[]. \n
   *  - Faster algorithm when \a coord0[] and \a coord1[] are both type float.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param[out] coord1[] Output coordinates.
   */
  template <typename DTYPE>
  void copy_coord(const DTYPE dimension, 
                  const float * coord0[], const float * coord1[])
  {
    std::copy(coord0, coord0+dimension, coord1);
  }


  /*!
   *  @brief Copy \a coord0[] to \a coord1[]. \n
   *  - Faster algorithm when \a coord0[] and \a coord1[] are both type double.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param[out] coord1[] Output coordinates.
   */
  template <typename DTYPE>
  void copy_coord(const DTYPE dimension, 
                  const double * coord0[], const double * coord1[])
  {
    std::copy(coord0, coord0+dimension, coord1);
  }


  /*!
   *  @brief Add \a coord0[] to \a coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param[out] coord2[] Output coordinate.
   *    - coord2[] = (\a coord0[] + \a coord1[]).
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void add_coord(const DTYPE dimension, const CTYPE0 coord0[],
                 const CTYPE1 coord1[], CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { coord2[d] = coord0[d] + coord1[d]; };
  }


  /*!
   *  @brief Subtract \a coord1[] from \a coord0[].
   *  @param dimension = Coordinate dimension (= number of coordinates.)
   *  @param coord0 = Input coordinates.
   *  @param coord1 = Input coordinates.
   *  @param[out] coord2 = Output coordinate equal to 
   *     (\a coord0[] - \a coord1[]).
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void subtract_coord
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord2[d] = coord0[d] - coord1[d]; };
  }


  /*!
   *  @brief Multiply \a coord0[] by the scalar \a s.
   *  @param dimension = Coordinate dimension (= number of coordinates.)
   *  @param coord0 = Input coordinates.
   *  @param s = Scalar.
   *  @param[out] coord1 = Output coordinate equal to (\a s * \a coord0[]).
   */
  template <typename DTYPE, typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void multiply_coord
  (const DTYPE dimension, const STYPE s, const CTYPE0 coord0[],
   CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = s*coord0[d]; };
  }


  /*!
   *  @brief Divide \a coord0[] by the scalar \a s.
   *  @param dimension = Coordinate dimension (= number of coordinates.)
   *  @param coord0 = Input coordinates.
   *  @param s = Non-zero scalar.  
   *  @pre Scalar s is not zero.
   *  @param[out] coord1 = Output coordinate equal to (\a s * \a coord0[]).
   */
  template <typename DTYPE, typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void divide_coord
  (const DTYPE dimension, const STYPE s, const CTYPE0 coord0[],
   CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = coord0[d]/s; };
  }

  /*!
   * Scale each \a coord0[d] by \a scale[d].
   * @param dimension Coordinate dimension (= number of coordinates.)
   * @param coord0[] Input coordinates.
   * @param scale[] Scale factors.
   * @param[out] coord1[] Output coordinates.
   *   - coord1[d] =  (\a scale[d] * \a coord0[d]).
   */
  template <typename DTYPE, typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void scale_coord
  (const DTYPE dimension, const STYPE scale[], const CTYPE0 coord0[],
   CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = scale[d]*coord0[d]; };
  }

  /*!
   * Reverse scaling of each \a coord0[d].
   * @param dimension Coordinate dimension (= number of coordinates.)
   * @param coord0[] Input coordinates.
   * @param scale[] Scale factors.
   * @param[out] coord1[] Output coordinates.
   *   - coord1[d] =  (\a coord0[d] / \a scale[d]).
   * @pre \a scale[d] is not zero for any \a d.
   */
  template <typename DTYPE, typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void reverse_coord_scaling
  (const DTYPE dimension, const STYPE scale[], const CTYPE0 coord0[],
   CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { coord1[d] = coord0[d]/scale[d]; };
  }

  /*!
   *  Add \a s * \a coord0[] to \a coord1[]).
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param s Scaling factor.
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param[out] coord2 Output coordinate.
   *    - coord2[] = \a s* \a coord[0] + \a coord1[].
   */
  template <typename DTYPE, typename STYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void add_scaled_coord
  (const DTYPE dimension, const STYPE s, const CTYPE0 coord0[],
   const CTYPE1 coord1[], CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { coord2[d] = s*coord0[d] + coord1[d]; };
  }

  /// Compute absolute value of the each coord.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  void abs_coord
  (const DTYPE dimension, const CTYPE0 coord0[], CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      coord1[d] = coord0[d];
      if (coord1[d] < 0) { coord1[d] = -coord1[d]; }
    }
  }


  /*!
   *  @brief Return index of coordinate with max absolute value.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *    @pre dimension > 0.
   */
  template <typename DTYPE, typename CTYPE0>
  DTYPE get_index_of_max_abs_coord
  (const DTYPE dimension, const CTYPE0 coord0[])
  {
    DTYPE index = 0;
    if (dimension < 2) 
      { return(index); }

    CTYPE0 c = std::abs(coord0[0]);

    for (DTYPE d = 1; d < dimension; d++) {
      const CTYPE0 c_d = std::abs(coord0[d]);

      if (c < c_d) {
        index = d;
        c = c_d;
      }
    }

    return(index);
  }


  /*!
   *  @brief Return index of coordinate with min absolute value.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *    @pre dimension > 0.
   */
  template <typename DTYPE, typename CTYPE0>
  DTYPE get_index_of_min_abs_coord
  (const DTYPE dimension, const CTYPE0 coord0[])
  {
    DTYPE index = 0;
    if (dimension < 2) 
      { return(index); }

    CTYPE0 c = std::abs(coord0[0]);

    for (DTYPE d = 1; d < dimension; d++) {
      const CTYPE0 c_d = std::abs(coord0[d]);
      if (c > c_d) {
        index = d;
        c = c_d;
      }
    }

    return(index);
  }


  /*!
   *  @brief Return index of coordinate with max value.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *    @pre dimension > 0.
   */
  template <typename DTYPE, typename CTYPE0>
  DTYPE get_index_of_max_coord
  (const DTYPE dimension, const CTYPE0 coord0[])
  {
    DTYPE index = 0;
    if (dimension < 2) 
      { return(index); }

    CTYPE0 c = coord0[0];

    for (DTYPE d = 1; d < dimension; d++) {
      const CTYPE0 c_d = coord0[d];

      if (c < c_d) {
        index = d;
        c = c_d;
      }
    }

    return(index);
  }


  /*!
   *  @brief Compute the midpoint of two coordinates.
   *  @param dimension = Coordinate dimension (= number of coordinates.)
   *  @param coord0 = Input coordinates.
   *  @param coord1 = Input coordinates.
   *  @param[out] midpoint_coord = Coordinates of midpoint.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void compute_midpoint
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   CTYPE2 midpoint_coord[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { midpoint_coord[d] = (coord0[d] + coord1[d])/2.0; };
  }


  /*!
   *  @brief Compute midpoint of line segment midpoints.
   *  - Line segments are (endpoint[0], endpoint[1]) and (endpoint[2], endpoint[3]).
   *  - Equivalent to computing the centroid of 4 points.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPE4>
  void compute_midpoint_of_midpoints
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const CTYPE2 coord2[], const CTYPE3 coord3[],
   CTYPE4 centroid_coord[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { centroid_coord[d] = (coord0[d] + coord1[d] + coord2[d] + coord3[d])/4.0; };
  }


  /*!
   *  @brief Compute sum of squares of coordinates.
   *  - Equivalent to computing the square of the magnitude.
   *  @param dimension = Coordinate dimension (= number of coordinates.)
   *  @param coord0 = Input coordinates.
   *  @param[out] sum = sum of squares
   */
  template <typename DTYPE, typename CTYPE0, typename STYPE>
  void compute_sum_of_squares(const DTYPE dimension, const CTYPE0 coord0[],
                              STYPE & sum)
  {
    sum = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { sum = sum + coord0[d]*coord0[d]; }
  }


  /*!
   *  @brief Compute magnitude of coordinate vector.
   *  - Equivalent to computing the square root of the sum 
   *    of squares of the coordinates.
   *  - Fast implementation, but returns 0 when the value of (magnitude^2) 
   *    is less than the machine epsilon of CTYPE0 or CTYPE1.
   *  @param dimension = Coordinate dimension (= number of coordinates.)
   *  @param coord0 = Input coordinates.
   *  @param[out] magnitude = magnitude
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  void compute_magnitude(const DTYPE dimension, const CTYPE0 coord0[],
                         CTYPE1 & magnitude)
  {
    compute_sum_of_squares(dimension, coord0, magnitude);
    magnitude = std::sqrt(magnitude);
  }

  //@}


  // **************************************************
  //! @name Vector operations
  // **************************************************

  //@{

  /*!
   *  @brief Compute inner product of \a coord0[] and \a coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param[out] product Inner product.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename STYPE>
  void compute_inner_product(const DTYPE dimension, const CTYPE0 coord0[],
                             const CTYPE1 coord1[], STYPE & product)
  {
    product = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { product = product + coord0[d]*coord1[d]; }
  }


  /*!
   *  @brief Project vector \a v onto direction \a dir[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param v[] Input vector coordinates.
   *  @param dir[] Direction.
   *  @param[out] v_proj[] Projection of \a v onto \a dir[].
   *     - v_proj[] = inner_product(v[],dir[]) * dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2>
  void project_vector(const DTYPE dimension, const VTYPE0 v[], 
                      const VTYPE1 dir[], VTYPE2 v_proj[])
  {
    VTYPE0 inner_product;
    compute_inner_product(dimension, v, dir, inner_product);
    for (DTYPE d = 0; d < dimension; d++) 
      { v_proj[d] = dir[d]*inner_product; }
  }


  /*!
   *  @brief Normalize vector.
   *  - If magnitude of \a coord0[] is less than or equal 
   *    to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input vector coordinates.
   *  @param max_small_magnitude If magnitude is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param[out] coord1[] Normalized coordinates or the zero vector.
   *  @param[out] magnitude Magnitude of \a coord0[].
   *  @param[out] flag_zero True if \a coord1[] is set to the zero vector.
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename MTYPE1>
  void normalize_vector
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const MTYPE0 max_small_magnitude, CTYPE1 coord1[],
   MTYPE1 & magnitude, bool & flag_zero)
  { 
    compute_magnitude(dimension, coord0, magnitude); 
    if (magnitude > max_small_magnitude) {
      multiply_coord(dimension, 1.0/magnitude, coord0, coord1);
      flag_zero = false;
    }
    else {
      set_coord(dimension, 0, coord1);
      flag_zero = true;
    }
  }


  /*!
   *  @brief Normalize vector.
   *  If magnitude of \a coord0[] is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set 
   *     to the zero vector. \n
   *  Version which does not return \a flag_zero.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input vector coordinates.
   *  @param max_small_magnitude If magnitude is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param[out] coord1[] Normalized coordinates or the zero vector.
   *  @param[out] magnitude Magnitude of \a coord0[].
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename MTYPE1>
  void normalize_vector
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const MTYPE0 max_small_magnitude, CTYPE1 coord1[],
   MTYPE1 & magnitude)
  { 
    compute_magnitude(dimension, coord0, magnitude); 
    if (magnitude > max_small_magnitude) {
      multiply_coord(dimension, 1.0/magnitude, coord0, coord1);
    }
    else {
      set_coord(dimension, 0, coord1);
    }
  }


  /*!
   *  @brief Normalize vector.
   *  -  If magnitude of \a coord0[] is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set 
   *     to the zero vector. \n
   *  -  Version which does not return \a magnitude or \a flag_zero.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input vector coordinates.
   *  @param max_small_magnitude If magnitude is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param[out] coord1[] Normalized coordinates or the zero vector.
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE, 
            typename CTYPE1>
  void normalize_vector
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const MTYPE max_small_magnitude, CTYPE1 coord1[])
  { 
    double magnitude;
    normalize_vector
      (dimension, coord0, max_small_magnitude, coord1, magnitude);
  }


  /*!
   *  @brief Normalize vector.
   *  -  Robust, more numerically stable, but more expensive, computation.
   *  If magnitude of \a coord0[] is less than or equal 
   *     to \a max_small_magnitude, then \a coord1[] is set to the zero vector.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input vector coordinates.
   *  @param max_small_magnitude If magnitude is less than or equal 
   *    to \a max_small_magnitude, set coord1[] to the zero vector.
   *  @param[out] coord1[] Normalized coordinates or the zero vector.
   *  @param[out] magnitude Magnitude of \a coord0[].
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename MTYPE1>
  void normalize_vector_robust
  (const DTYPE dimension, const CTYPE0 coord0[], 
   const MTYPE0 max_small_magnitude, CTYPE1 coord1[], MTYPE1 & magnitude)
  { 
    MTYPE1 max_abs_coord;  // Use MTYPE1 in case CTYPE0 is int.

    magnitude = 0;
    if (dimension <= 0) { return; }

    // compute maximum absolute value of any coordinate.
    max_abs_coord = coord0[0];
    if (max_abs_coord < 0) { max_abs_coord = -max_abs_coord; };
    for (DTYPE d = 1; d < dimension; d++) {
      MTYPE1 abs_coord = coord0[d];
      if (abs_coord < 0) { abs_coord = -abs_coord; };
      if (max_abs_coord < abs_coord)
        { max_abs_coord = abs_coord; }
    }

    if (max_abs_coord > 0) {
      for (DTYPE d = 0; d < dimension; d++) {
        MTYPE1 c = coord0[d]/max_abs_coord;
        magnitude += c*c;
      }
      magnitude = std::sqrt(magnitude);
      magnitude = max_abs_coord*magnitude;
    }

    if (magnitude > max_small_magnitude) {
      divide_coord(dimension, magnitude, coord0, coord1);
    }
    else {
      set_coord(dimension, 0, coord1);
    }
  }


  /*!
   *  @brief Compute unit vector from coord0[] to coord1[].
   *  -  If distance from \a coord0[] to \a coord1[]is less than or equal 
   *     to \a max_small_magnitude, then \a coord2[] is set to the zero vector.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param max_small_magnitude If distance from \a coord0[] to \a coord1[]
   *     is less than or equal to \a max_small_magnitude, then \a coord2[] 
   *      is set to the zero vector.
   *  @param[out] coord2[] Normalized coordinates or the zero vector.
   *  @param[out] magnitude Magnitude of (\a coord1[] - \a coord0[]).
   *  @param[out] flag_zero True if \a coord2[] is set to the zero vector.
   *  @pre max_small_magnitude >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename CTYPE2, typename MTYPE1>
  void compute_unit_vector
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const MTYPE0 max_small_magnitude, CTYPE2 coord2[],
   MTYPE1 & magnitude, bool & flag_zero)
  {
    subtract_coord(dimension, coord1, coord0, coord2);
    normalize_vector
      (dimension, coord2, max_small_magnitude, coord2, magnitude, flag_zero);
  }


  /*!
   *  @brief Compute component of vector orthogonal to given direction
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param v[] Input vector coordinates.
   *  @param dir[] Direction.
   *  @param[out] v_orth[] Component of \a v orthogonal to \a dir[].
   *     - v_orth[] = (\a v - v_proj[]) where v_proj[] is the projection
   *                  of \a v onto \a dir[].
   *  @pre dir[] is a unit vector or a zero vector.
   *     - If dir[] is 0 vector, then v_orth[] = v[].
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2>
  void compute_orthogonal_vector_component
  (const DTYPE dimension, const VTYPE0 v[], 
   const VTYPE1 dir[], VTYPE2 v_orth[])
  {
    VTYPE0 inner_product;
    compute_inner_product(dimension, v, dir, inner_product);
    for (DTYPE d = 0; d < dimension; d++) 
      { v_orth[d] = v[d]-dir[d]*inner_product; }
  }


  /*!
   *  @brief Compute component of vector orthogonal to given direction
   *  - DEPRECATED. Should be compute_orthogonal_vector_component.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param v[] Input vector coordinates.
   *  @param dir[] Direction.
   *  @param[out] v_orth[] Component of \a v orthogonal to \a dir[].
   *     - v_orth[] = (\a v - v_proj[]) where v_proj[] is the projection
   *                  of \a v onto \a dir[].
   *  @pre dir[] is a unit vector or a zero vector.
   *     - If dir[] is 0 vector, then v_orth[] = v[].
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2>
  void compute_orthogonal_vector
  (const DTYPE dimension, const VTYPE0 v[], 
   const VTYPE1 dir[], VTYPE2 v_orth[])
  {
    compute_orthogonal_vector_component(dimension, v, dir, v_orth);
  }


  /*!
   *  @brief Compute component of vector orthogonal to two given directions.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *    @pre dimension >= 3.
   *  @param v[] Input vector coordinates.
   *  @param dir0[] Direction 0. 
   *  @pre dir0[] is a unit vector or a zero vector.
   *  @param dir1[] Direction 2. 
   *    @pre dir1[] is a unit vector orthogonal to dir0[], or a zero vector.
   *  @param[out] v_orth[] Component of \a v orthogonal 
   *    to \a dir0[] and \a dir1[].
   */
  template <typename DTYPE, typename VTYPE0, typename VTYPE1,
            typename VTYPE2, typename VTYPE3>
  void compute_orthogonalII_vector_component
  (const DTYPE dimension, const VTYPE0 v[], 
   const VTYPE1 dir0[], VTYPE2 dir1[], VTYPE3 v_orth[])
  {
    compute_orthogonal_vector_component(dimension, v, dir0, v_orth);
    compute_orthogonal_vector_component(dimension, v_orth, dir1, v_orth);
  }


  /*!
   *  @brief Compute three orthogonal unit vectors.
   *  - Return false if v0, v1 or v2 are (near) zero 
   *    or v0 and v1 are collinear or v0, v1 and v2 are coplanar.
   *  - Note: max_small_magnitude is at or near zero can result
   *    in extremely incorrect results.
   *  @param[out] fail_type Type of failure.
   *  - If fail_type == 0, then v0 is zero vector.
   *  - If fail_type == 1, then v0 and v1 are collinear.
   *  - If fail_type == 2, then v0, v1 and v2 are coplanar.
   *  - If function returns true, then fail_type is undefined.
   */
  template <typename DTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE,
            typename WTYPE0, typename WTYPE1, typename WTYPE2,
            typename FTYPE>
  bool compute_orthogonal_unit_vectors
  (const DTYPE dimension, 
   const VTYPE0 v0[], const VTYPE1 v1[], const VTYPE2 v2[],
   const MTYPE max_small_magnitude,
   WTYPE0 w0[], WTYPE1 w1[], WTYPE2 w2[],
   FTYPE & fail_type)
  {
    double magnitude;
    double tolerance_for_check(0.001);
    double inner_product;

    set_coord(dimension, 0, w0);
    set_coord(dimension, 0, w1);
    set_coord(dimension, 0, w2);

    fail_type = 0;

    // Compute unit vector w0 = v0/|v0|.
    normalize_vector(dimension, v0, max_small_magnitude, w0, magnitude);
    if (magnitude <= max_small_magnitude) { return(false); }

    fail_type = 1;

    // Compute vector w1 orthogonal to w0.
    compute_orthogonal_vector(dimension, v1, w0, w1);
    normalize_vector(dimension, w1, max_small_magnitude, w1, magnitude);
    if (magnitude <= max_small_magnitude) { return(false); }

    // Check w1 is orthogonal to w0.
    compute_inner_product(dimension, w0, w1, inner_product);
    if (inner_product > tolerance_for_check || 
        inner_product < -tolerance_for_check) {
      set_coord(dimension, 0, w1);
      return(false);
    }

    fail_type = 2;

    // Compute vector w2 orthogonal to w0 and w1.
    compute_orthogonal_vector(dimension, v2, w0, w2);
    compute_orthogonal_vector(dimension, w2, w1, w2);
    normalize_vector(dimension, w2, max_small_magnitude, w2, magnitude);
    if (magnitude <= max_small_magnitude) { return(false); }

    // Check w2 is orthogonal to w1 and w0.
    compute_inner_product(dimension, w0, w2, inner_product);
    if (inner_product > tolerance_for_check || 
        inner_product < -tolerance_for_check) {
      set_coord(dimension, 0, w2);
      return(false);
    }
    compute_inner_product(dimension, w1, w2, inner_product);
    if (inner_product > tolerance_for_check || 
        inner_product < -tolerance_for_check) {
      set_coord(dimension, 0, w2);
      return(false);
    }

    return(true);
  }


  /*!
   *  @brief Compute three orthogonal unit vectors.
   *  - Does not set a fail type.
   */
  template <typename DTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename MTYPE,
            typename WTYPE0, typename WTYPE1, typename WTYPE2>
  bool compute_orthogonal_unit_vectors
  (const DTYPE dimension, 
   const VTYPE0 v0[], const VTYPE1 v1[], const VTYPE2 v2[],
   const MTYPE max_small_magnitude,
   WTYPE0 w0[], WTYPE1 w1[], WTYPE2 w2[])
  {
    int fail_type;

    bool flag =
      compute_orthogonal_unit_vectors
      (dimension, v0, v1, v2, max_small_magnitude, w0, w1, w2, fail_type);

    return(flag);
  }


  /*!
   *  Compute 3 vectors formed by 3 triangle edges.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param[out] u01[] Vector from coord0 to coord1.
   *  @param[out] u02[] Vector from coord0 to coord2.
   *  @param[out] u12[] Vector from coord1 to coord2.
   *  @param[out] mag01 Magnitude of u01.
   *  @param[out] mag02 Magnitude of u02.
   *  @param[out] mag12 Magnitude of u12.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, typename CTYPE5,
            typename MTYPE0, typename MTYPE1, typename MTYPE2>
  void compute_triangle_edge_vectors
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], const CTYPE2 coord2[],
   CTYPE3 u01[], CTYPE4 u02[], CTYPE5 u12[],
   MTYPE0 & mag01, MTYPE1 & mag02, MTYPE2 & mag12)
  {
    IJK::subtract_coord(dimension, coord1, coord0, u01);
    IJK::subtract_coord(dimension, coord2, coord0, u02);
    IJK::subtract_coord(dimension, coord2, coord1, u12);

    compute_magnitude(dimension, u01, mag01);
    compute_magnitude(dimension, u02, mag02);
    compute_magnitude(dimension, u12, mag12);
  }


  /*!
   *  Compute 3 vectors formed by 3 triangle edges.
   *  - Version using array mag[], instead of mag01, mag02 and mag12.
   *  @param[out] mag[] Vector magnitudes.
   *  - mag[0]: Magnitude of vector u01[].
   *  - mag[1]: Magnitude of vector u02[].
   *  - mag[2]: Magnitude of vector u12[].
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, typename CTYPE5,
            typename MTYPE>
  void compute_triangle_edge_vectors
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], const CTYPE2 coord2[],
   CTYPE3 u01[], CTYPE4 u02[], CTYPE5 u12[],
   MTYPE mag[3])
  {
    compute_triangle_edge_vectors
      (dimension, coord0, coord1, coord2, u01, u02, u12,
       mag[0], mag[1], mag[2]);
  }


  /*!
   *  Compute 3 vectors formed by 3 triangle edges.
   *  - Compute vectors in direction given by triangle orientation
   *    (coord0[], coord[1], coord2[]).
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param[out] u01[] Vector from coord0 to coord1.
   *  @param[out] u12[] Vector from coord1 to coord2.
   *  @param[out] u20[] Vector from coord2 to coord0.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPEU0, typename CTYPEU1, typename CTYPEU2>
  void compute_triangle_edge_vectors_oriented
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], const CTYPE2 coord2[],
   CTYPEU0 u01[], CTYPEU1 u12[], CTYPEU2 u20[])
  {
    IJK::subtract_coord(dimension, coord1, coord0, u01);
    IJK::subtract_coord(dimension, coord2, coord1, u12);
    IJK::subtract_coord(dimension, coord0, coord2, u20);
  }


  /*!
   *  Compute 3 vectors formed by 3 triangle edges.
   *  - Compute vectors in direction given by triangle orientation
   *    (coord0[], coord[1], coord2[]).
   *  - Version that returns magnitudes of each vector.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param[out] u01[] Vector from coord0 to coord1.
   *  @param[out] u12[] Vector from coord1 to coord2.
   *  @param[out] u20[] Vector from coord2 to coord0.
   *  @param[out] length[i] Length of i'th edge.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPEU0, typename CTYPEU1, typename CTYPEU2,
            typename LTYPE>
  void compute_triangle_edge_vectors_oriented
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], const CTYPE2 coord2[],
   CTYPEU0 u01[], CTYPEU1 u12[], CTYPEU2 u20[],
   LTYPE length[3])
  {
    compute_triangle_edge_vectors_oriented
      (dimension, coord0, coord1, coord2, u01, u12, u20);

    compute_magnitude(dimension, u01, length[0]);
    compute_magnitude(dimension, u12, length[1]);
    compute_magnitude(dimension, u20, length[2]);
  }

  //@}


  // **************************************************
  //! @name Compute angles
  // **************************************************

  //@{

  /*!
   *  @brief Compute cosine of the angle between two vectors.
   *  - If either vector has magnitude <= max_small_magnitude, 
   *    return zero and set flag_zero to true.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_angle Cosine of angle between coord0[] and coord1[].
   *  @param[out] flag_zero True, if coord0[] or coord1[]
   *    have zero magnitude. False, otherwise.
   *    - Note: Computed magnitude may be zero even if vectors are not zero.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_angle
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const MTYPE max_small_magnitude, COS_TYPE & cos_angle, bool & flag_zero)
  {
    CTYPE0 mag0;
    CTYPE1 mag1;
    compute_magnitude(dimension, coord0, mag0);
    compute_magnitude(dimension, coord1, mag1);
    if (mag0 > max_small_magnitude && mag1 > max_small_magnitude) {
      compute_inner_product(dimension, coord0, coord1, cos_angle);
      flag_zero = false;
      cos_angle = (cos_angle/mag0)/mag1;

      // Bound cos_angle to [-1,1].
      // Numerical error could cause cos_angle to be outside bounds.
      if (cos_angle < -1) { cos_angle = -1; }
      if (cos_angle > 1) { cos_angle = 1; }
    }
    else {
      cos_angle = 0;
      flag_zero = true;
    }
  }


  /*!
   *  @brief Compute cosine of the angle between two vectors.
   *  - Version without max_small_magnitude parameter.
   *  - If either vector has (near) zero magnitude, return zero
   *    and set flag_zero to true.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename COS_TYPE>
  void compute_cos_angle
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   COS_TYPE & cos_angle, bool & flag_zero)
  {
    compute_cos_angle(dimension, coord0, coord1, 0.0, cos_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of angle (coord0, coord1, coord2).
   *  - If either coord1 or coord2 have (nearly) the same coordinates
   *    as coord0, return 0 and set flag_identical to true.
   * @param[out] cos_angle Cosine of angle(coord0, coord1, coord2).
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename MTYPE>
  void compute_cos_triangle_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   const MTYPE max_small_magnitude, CTYPE3 & cos_angle, bool & flag_identical)
  {
    IJK::ARRAY<CTYPE0> u1(dimension);
    IJK::ARRAY<CTYPE1> u2(dimension);

    IJK::subtract_coord(dimension, coord0, coord1, u1.Ptr());
    IJK::subtract_coord(dimension, coord2, coord1, u2.Ptr());
    IJK::compute_cos_angle(dimension, u1.PtrConst(), u2.PtrConst(),
                           max_small_magnitude, cos_angle, flag_identical);
  }


  /*!
   *  Compute the cosine of angle (coord0, coord1, coord2).
   *  - Version without max_small_magnitude parameter.
   *  - If either coord1 or coord2 have (nearly) the same coordinates
   *    as coord0, return zero and set flag_identical to true.
   * @param[out] cos_angle Cosine of angle(coord0, coord1, coord2).
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3>
  void compute_cos_triangle_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   CTYPE3 & cos_angle, bool & flag_identical)
  {
    IJK::compute_cos_triangle_angle
      (dimension, coord0, coord1, coord2, 0.0, cos_angle, flag_identical);
  }


  /*!
   *  @brief Compute the cosine of angle (iv0, iv1, iv2).
   * - Vertex coordinates are coord[iv0*dimension], coord[iv1*dimension],
   *   and coord[iv2*dimension].
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename COS_TYPE, typename MTYPE>
  void compute_cos_triangle_angle_coord_list
  (const DTYPE dimension, const CTYPE * vertex_coord, 
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2, 
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_angle, bool & flag_identical)
  {
    const CTYPE * coord0 = vertex_coord+dimension*iv0;
    const CTYPE * coord1 = vertex_coord+dimension*iv1;
    const CTYPE * coord2 = vertex_coord+dimension*iv2;

    compute_cos_triangle_angle
      (dimension, coord0, coord1, coord2, max_small_magnitude, 
       cos_angle, flag_identical);
  }


  /*!
   *  @brief Compute the cosine of angle (iv0, iv1, iv2).
   *  - Version without max_small_magnitude parameter.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename COS_TYPE>
  void compute_cos_triangle_angle_coord_list
  (const DTYPE dimension, const CTYPE * vertex_coord, 
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2, 
   COS_TYPE & cos_angle, bool & flag_identical)
  {
    compute_cos_triangle_angle_coord_list
      (dimension, vertex_coord, iv0, iv1, iv2, 0.0, cos_angle, flag_identical);
  }


  /*!
   *  @brief Compute the cosine of the angle at the k'th polygon vertex.
   *  - If vertex vertex_list[k] has (nearly) the same coordinates
   *    as either of its neighbors, returns cosine equal to 0 
   *    and sets flag_identical to true.
   * @param vertex_list[] Array of polygon vertices.
   * @param num_poly_vertices Number of polygon vertices.
   * @param[out] cos_angle
   *    Cosine of angle(vertex_list[k-1],vertex_list[k]],vertex_list[k+1]).
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename NTYPE, typename ITYPE, typename MTYPE,
            typename COS_TYPE>
  void compute_cos_angle_at_polygon_vertex
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const VTYPE * vertex_list, const NTYPE num_poly_vertices,
   const ITYPE k1, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_angle, bool & flag_identical)
  {
    if (num_poly_vertices == 0) {
      cos_angle = 0.0;
      flag_identical = true;
      return;
    }

    const ITYPE k0 = (k1 + (num_poly_vertices-1))%num_poly_vertices;
    const ITYPE k2 = (k1 + 1)%num_poly_vertices;

    compute_cos_triangle_angle_coord_list
      (dimension, vertex_coord, vertex_list[k0], vertex_list[k1],
       vertex_list[k2], max_small_magnitude, cos_angle, flag_identical);
  }


  /*!
   *  @brief Compute the cosine of the angle at the k'th polygon vertex.
   *  - Version using list of pointers to vertex coordinates.
   *  - If vertex vertex_list[k] has (nearly) the same coordinates
   *    as either of its neighbors, returns cosine equal to 0 
   *    and sets flag_identical to true.
   * @param vertex_list[] Array of polygon vertices.
   * @param num_poly_vertices Number of polygon vertices.
   * @param[out] cos_angle
   *    Cosine of angle(vertex_list[k-1],vertex_list[k]],vertex_list[k+1]).
   */
  template <typename DTYPE, typename CTYPE,
            typename NTYPE, typename ITYPE, typename MTYPE,
            typename COS_TYPE>
  void compute_cos_angle_at_polygon_vertex
  (const DTYPE dimension, const NTYPE num_poly_vertices,
   const CTYPE * const vcoord_ptr[], const ITYPE k1, 
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_angle, bool & flag_identical)
  {
    if (num_poly_vertices == 0) {
      cos_angle = 0.0;
      flag_identical = true;
      return;
    }

    const ITYPE k0 = (k1 + (num_poly_vertices-1))%num_poly_vertices;
    const ITYPE k2 = (k1 + 1)%num_poly_vertices;

    compute_cos_triangle_angle
      (dimension, vcoord_ptr[k0], vcoord_ptr[k1], vcoord_ptr[k2], 
       max_small_magnitude, cos_angle, flag_identical);
  }


  /*!
   *  @brief Compute the angle in degrees at the k'th polygon vertex.
   *  - Should be used only for output of angle.
   *  - If vertex vertex_list[k] has (nearly) the same coordinates
   *    as either of its neighbors, returns angle equal to 0 
   *    and sets flag_identical to true.
   * @param vertex_list[] Array of polygon vertices.
   * @param num_poly_vertices Number of polygon vertices.
   * @param[out] angle
   *    angle(vertex_list[k-1],vertex_list[k]],vertex_list[k+1]).
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename NTYPE, typename ITYPE, typename MTYPE,
            typename ATYPE>
  void compute_angle_in_degrees_at_polygon_vertex
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const VTYPE * vertex_list, const NTYPE num_poly_vertices,
   const ITYPE k1, const MTYPE max_small_magnitude, 
   ATYPE & angle, bool & flag_identical)
  {
    CTYPE cos_angle;
    compute_cos_angle_at_polygon_vertex
      (dimension, vertex_coord, vertex_list, num_poly_vertices,
       k1, max_small_magnitude, cos_angle, flag_identical);

    angle = std::acos(cos_angle) * 180.0/std::numbers::pi;
  }


  /*!
   *  @brief Compute the cosine of smallest angle of triangle (coord0, coord1, coord2).
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if two or three triangles edges have
   *    length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    IJK::ARRAY<CTYPE0> u01(dimension);
    IJK::ARRAY<CTYPE1> u02(dimension);
    IJK::ARRAY<CTYPE2> u12(dimension);
    CTYPE0 mag[3];
    COS_TYPE cos_angleA, cos_angleB, cos_angleC;

    compute_triangle_edge_vectors
      (dimension, coord0, coord1, coord2, u01.Ptr(), u02.Ptr(), u12.Ptr(),
       mag);

    if (mag[0] > max_small_magnitude && mag[1] > max_small_magnitude &&
        mag[2] > max_small_magnitude) {

      flag_zero = false;

      compute_inner_product
        (dimension, u01.PtrConst(), u02.PtrConst(), cos_angleA);
      cos_angleA = (cos_angleA/mag[0])/mag[1];

      // Note: Angle b is formed by vectors -u01 and u12.
      compute_inner_product
        (dimension, u01.PtrConst(), u12.PtrConst(), cos_angleB);
      cos_angleB = -(cos_angleB/mag[0])/mag[2];

      compute_inner_product
        (dimension, u02.PtrConst(), u12.PtrConst(), cos_angleC);
      cos_angleC = (cos_angleC/mag[1])/mag[2];

      if (cos_angleA > cos_angleB) { cos_min_angle = cos_angleA; }
      else { cos_min_angle = cos_angleB; }
      if (cos_min_angle < cos_angleC) { cos_min_angle = cos_angleC; }
    }
    else if (mag[0] > max_small_magnitude && mag[1] > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u01.PtrConst(), u02.PtrConst(), cos_angleA);
      cos_min_angle = (cos_angleA/mag[0])/mag[1];
    }
    else if (mag[0] > max_small_magnitude && mag[2] > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u01.PtrConst(), u12.PtrConst(), cos_angleB);
      cos_min_angle = -(cos_angleB/mag[0])/mag[2];
    }
    else if (mag[1] > max_small_magnitude && mag[2] > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u02.PtrConst(), u12.PtrConst(), cos_angleC);
      cos_min_angle = (cos_angleC/mag[1])/mag[2];
    }
    else {
      cos_min_angle = 0.0;
      flag_zero = true;
    }

    // Bound cos_min_angle to [-1,1] in case
    //   numerical error causes cos_min_angle to be outside bounds.
    if (cos_min_angle < -1) { cos_min_angle = -1; }
    if (cos_min_angle > 1) { cos_min_angle = 1; }
  }

  /*!
   *  @brief Compute the cosine of smallest angle of triangle (coord0, coord1, coord2).
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord[i] Pointers to coordinate of i'th vertex.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] flag_zero True, if two or three triangles edges have
   *    length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPEV,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, const CTYPEV * coord[3],
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_triangle_angle
      (dimension, coord[0], coord[1], coord[2], max_small_magnitude,
       cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of triangle (iv0,iv1,iv2).
   *  - Version using array vcoord of vertex coordinates.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   */
  template <typename DTYPE, typename CTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2,  
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const CTYPE * vcoord0 = vertex_coord+iv0*dimension;
    const CTYPE * vcoord1 = vertex_coord+iv1*dimension;
    const CTYPE * vcoord2 = vertex_coord+iv2*dimension;

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoord2, max_small_magnitude, 
       cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of triangle (iv0,iv1,iv2).
   *  - Version using array of triangle vertices.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const VTYPE triangle_vert[3],
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_triangle_angle
      (dimension, vertex_coord,
       triangle_vert[0], triangle_vert[1], triangle_vert[2],
       max_small_magnitude, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of triangle (iv0,iv1,iv2).
   *  - Version using STL vector vcoord of vertex coordinates.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   */
  template <typename DTYPE, typename CTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2,  
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const CTYPE * vertex_coord_ptr = IJK::vector2pointer(vertex_coord);
    const CTYPE * vcoord0 = vertex_coord_ptr+iv0*dimension;
    const CTYPE * vcoord1 = vertex_coord_ptr+iv1*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+iv2*dimension;

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoord2, max_small_magnitude, 
       cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of triangle (iv0,iv1,iv2).
   *  - Version using STL vector vcoord of vertex coordinates.
   *  - Version using array of triangle vertices.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   * @param triangle_vert[] Triangle vertices.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const VTYPE triangle_vert[3],
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_min_triangle_angle
      (dimension, vertex_coord, 
       triangle_vert[0], triangle_vert[1], triangle_vert[2],
       max_small_magnitude, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of triangle (coord0, coord1, coord2).
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   *  - Version that returns index of vertex with smallest angle.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] iangle Index (0, 1 or 2) of vertex with smallest angle.
   *  @param[out] flag_zero True, if two or three triangles edges have
   *    length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename MTYPE, typename COS_TYPE, typename ITYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, 
   const CTYPE0 * coord0, const CTYPE1 * coord1, const CTYPE2 * coord2,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, 
   ITYPE & iangle, bool & flag_zero)
  {
    IJK::ARRAY<CTYPE0> u01(dimension);
    IJK::ARRAY<CTYPE1> u02(dimension);
    IJK::ARRAY<CTYPE2> u12(dimension);
    CTYPE0 mag[3];
    COS_TYPE cos_angleA, cos_angleB, cos_angleC;

    compute_triangle_edge_vectors
      (dimension, coord0, coord1, coord2, u01.Ptr(), u02.Ptr(), u12.Ptr(),
       mag);

    if (mag[0] > max_small_magnitude && mag[1] > max_small_magnitude &&
        mag[2] > max_small_magnitude) {

      flag_zero = false;

      compute_inner_product
        (dimension, u01.PtrConst(), u02.PtrConst(), cos_angleA);
      cos_angleA = (cos_angleA/mag[0])/mag[1];

      // Note: Angle b is formed by vectors -u01 and u12.
      compute_inner_product
        (dimension, u01.PtrConst(), u12.PtrConst(), cos_angleB);
      cos_angleB = -(cos_angleB/mag[0])/mag[2];

      compute_inner_product
        (dimension, u02.PtrConst(), u12.PtrConst(), cos_angleC);
      cos_angleC = (cos_angleC/mag[1])/mag[2];

      if (cos_angleA > cos_angleB) {
        if (cos_angleA > cos_angleC) {
          cos_min_angle = cos_angleA;
          iangle = 0;
        }
        else {
          cos_min_angle = cos_angleC;
          iangle = 2;
        }
      }
      else if (cos_angleB > cos_angleC) {
        cos_min_angle = cos_angleB;
        iangle = 1;
      }
      else {
        cos_min_angle = cos_angleC;
        iangle = 2;
      }
    }
    else if (mag[0] > max_small_magnitude && mag[1] > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u01.PtrConst(), u02.PtrConst(), cos_angleA);
      cos_min_angle = (cos_angleA/mag[0])/mag[1];
      iangle = 0;
    }
    else if (mag[0] > max_small_magnitude && mag[2] > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u01.PtrConst(), u12.PtrConst(), cos_angleB);
      cos_min_angle = -(cos_angleB/mag[0])/mag[2];
      iangle = 1;
    }
    else if (mag[1] > max_small_magnitude && mag[2] > max_small_magnitude) {
      flag_zero = false;
      compute_inner_product
        (dimension, u02.PtrConst(), u12.PtrConst(), cos_angleC);
      cos_min_angle = (cos_angleC/mag[1])/mag[2];
      iangle = 2;
    }
    else {
      cos_min_angle = 0.0;
      iangle = 0;
      flag_zero = true;
    }

    // Bound cos_min_angle to [-1,1] in case
    //   numerical error causes cos_min_angle to be outside bounds.
    if (cos_min_angle < -1) { cos_min_angle = -1; }
    if (cos_min_angle > 1) { cos_min_angle = 1; }
  }


  /*!
   *  @brief Compute the cosine of smallest angle of triangle (iv0,iv1,iv2).
   *  - Version that returns index of vertex with smallest angle.
   *  - Version using STL vector vcoord of vertex coordinates.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   */
  template <typename DTYPE, typename CTYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2,  
            typename MTYPE, typename COS_TYPE, typename ITYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2,
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, ITYPE & iangle, bool & flag_zero)
  {
    const CTYPE * vertex_coord_ptr = IJK::vector2pointer(vertex_coord);
    const CTYPE * vcoord0 = vertex_coord_ptr+iv0*dimension;
    const CTYPE * vcoord1 = vertex_coord_ptr+iv1*dimension;
    const CTYPE * vcoord2 = vertex_coord_ptr+iv2*dimension;

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, vcoord2, max_small_magnitude, 
       cos_min_angle, iangle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of triangle (iv0,iv1,iv2).
   *  - Version that returns index of vertex with smallest angle.
   *  - Version using STL vector vcoord of vertex coordinates.
   *  - Version using array of triangle vertices.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the three angles.
   * @param triangle_vert[] Triangle vertices.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename MTYPE, typename ITYPE, typename COS_TYPE>
  void compute_cos_min_triangle_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const VTYPE triangle_vert[3], const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, ITYPE & iangle, bool & flag_zero)
  {
    compute_cos_min_triangle_angle
      (dimension, vertex_coord, 
       triangle_vert[0], triangle_vert[1], triangle_vert[2],
       max_small_magnitude, cos_min_angle, iangle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of quadrilateral
   *    (coord0, coord1, coord2, coord3).
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the four angles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] iangle Index of vertex at apex of min angle.
   *    - If (iangle = 0), then min angle is (coord3,coord0,coord1).
   *    - If (iangle = 1), then min angle is (coord0,coord1,coord2).
   *    - If (iangle = 2), then min angle is (coord1,coord3,coord3).
   *    - If (iangle = 3), then min angle is (coord2,coord3,coord0).
   *  @param[out] flag_zero True, if some quadrilateral edges has
   *    length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename CTYPE2, typename CTYPE3, 
            typename ATYPE, typename ITYPE, typename MTYPE>
  void compute_cos_min_quad_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3,
   const MTYPE max_small_magnitude, ATYPE & cos_min_angle, 
   ITYPE & iangle, bool & flag_zero)
  {
    IJK::ARRAY<CTYPE0> u01(dimension);
    IJK::ARRAY<CTYPE1> u12(dimension);
    IJK::ARRAY<CTYPE2> u23(dimension);
    IJK::ARRAY<CTYPE2> u30(dimension);
    CTYPE0 mag01, mag12, mag23, mag30;
    ATYPE cos_angle0, cos_angle1, cos_angle2, cos_angle3;

    IJK::subtract_coord(dimension, coord1, coord0, u01.Ptr());
    IJK::subtract_coord(dimension, coord2, coord1, u12.Ptr());
    IJK::subtract_coord(dimension, coord3, coord2, u23.Ptr());
    IJK::subtract_coord(dimension, coord0, coord3, u30.Ptr());

    compute_magnitude(dimension, u01.PtrConst(), mag01);
    compute_magnitude(dimension, u12.PtrConst(), mag12);
    compute_magnitude(dimension, u23.PtrConst(), mag23);
    compute_magnitude(dimension, u30.PtrConst(), mag30);
    
    if (mag01 > max_small_magnitude && mag12 > max_small_magnitude &&
        mag23 > max_small_magnitude && mag30 > max_small_magnitude) {

      flag_zero = false;

      compute_inner_product
        (dimension, u01.PtrConst(), u30.PtrConst(), cos_angle0);
      cos_angle0 = -(cos_angle0/mag01)/mag30;

      compute_inner_product
        (dimension, u01.PtrConst(), u12.PtrConst(), cos_angle1);
      cos_angle1 = -(cos_angle1/mag01)/mag12;

      compute_inner_product
        (dimension, u12.PtrConst(), u23.PtrConst(), cos_angle2);
      cos_angle2 = -(cos_angle2/mag12)/mag23;

      compute_inner_product
        (dimension, u23.PtrConst(), u30.PtrConst(), cos_angle3);
      cos_angle3 = -(cos_angle3/mag23)/mag30;

      if (cos_angle0 > cos_angle1) { 
        cos_min_angle = cos_angle0; 
        iangle = 0;
      }
      else { 
        cos_min_angle = cos_angle1; 
        iangle = 1;
      }

      if (cos_min_angle < cos_angle2) { 
        cos_min_angle = cos_angle2; 
        iangle = 2;
      }

      if (cos_min_angle < cos_angle3) { 
        cos_min_angle = cos_angle3; 
        iangle = 3;
      }

      // Bound cos_min_angle to [-1,1] in case
      //   numerical error causes cos_min_angle to be outside bounds.
      if (cos_min_angle < -1) { cos_min_angle = -1; }
      if (cos_min_angle > 1) { cos_min_angle = 1; }
    }
    else {
      // Default value for iangle.
      iangle = 0;
      flag_zero = true;
    }
  }


  /*!
   *  @brief Compute the cosine of smallest angle of quadrilateral 
   *    (iv0,iv1,iv2,iv3).
   *  - Version using array of vcoord of vertex coordinates.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   */
  template <typename DTYPE, typename CTYPE,
            typename IVTYPE0, typename IVTYPE1, 
            typename IVTYPE2, typename IVTYPE3, 
            typename ATYPE, typename ITYPE, typename MTYPE>
  void compute_cos_min_quad_angle
  (const DTYPE dimension, const CTYPE vertex_coord[],
   const IVTYPE0 iv0, const IVTYPE1 iv1, 
   const IVTYPE2 iv2, const IVTYPE3 iv3,
   const MTYPE max_small_magnitude, ATYPE & cos_min_angle, 
   ITYPE & iangle, bool & flag_zero)
  {
    const CTYPE * vcoord0 = vertex_coord+iv0*dimension;
    const CTYPE * vcoord1 = vertex_coord+iv1*dimension;
    const CTYPE * vcoord2 = vertex_coord+iv2*dimension;
    const CTYPE * vcoord3 = vertex_coord+iv3*dimension;

    compute_cos_min_quad_angle
      (dimension, vcoord0, vcoord1, vcoord2, vcoord3,
       max_small_magnitude, cos_min_angle, iangle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of quadrilateral 
   *    (iv0,iv1,iv2,iv3).
   *  - Version using STL vector coord of vertex coordinates.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   */
  template <typename DTYPE, typename CTYPE,
            typename IVTYPE0, typename IVTYPE1, 
            typename IVTYPE2, typename IVTYPE3, 
            typename ATYPE, typename ITYPE, typename MTYPE>
  void compute_cos_min_quad_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const IVTYPE0 iv0, const IVTYPE1 iv1, 
   const IVTYPE2 iv2, const IVTYPE3 iv3,
   const MTYPE max_small_magnitude, ATYPE & cos_min_angle, 
   ITYPE & iangle, bool & flag_zero)
  {
    compute_cos_min_quad_angle
      (dimension, IJK::vector2pointer(vertex_coord), iv0, iv1, iv2, iv3,
       max_small_magnitude, cos_min_angle, iangle, flag_zero);
  }


  /*!
   *  @brief Compute the cosine of smallest angle of quadrilateral 
   *    (iv0,iv1,iv2,iv3).
   *  - Version using STL vector coord of vertex coordinates.
   *  - Version using array of vertices.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   */
  template <typename DTYPE, typename CTYPE, typename IVTYPE,
            typename ATYPE, typename ITYPE, typename MTYPE>
  void compute_cos_min_quad_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vertex_coord,
   const IVTYPE iv[],
   const MTYPE max_small_magnitude, ATYPE & cos_min_angle, 
   ITYPE & iangle, bool & flag_zero)
  {
    compute_cos_min_quad_angle
      (dimension, vertex_coord, iv[0], iv[1], iv[2], iv[3],
       max_small_magnitude, cos_min_angle, iangle, flag_zero);
  }


  /*!
   *  @brief Compute cosine of min/max of polygon angles of polygon.
   *  - If two adjacent polygon vertices have (near) identical coordinates,
   *    the angle at those vertices is skipped.
   *  @param[out] cos_min Cosine of minimum polygon angle.
   *     - Note: Min angle has maximimum cosine.
   *  @param[out] cos_max Cosine of maximum polygon angle.
   *     - Note: Max angle has minimum cosine.
   *  @param[out] iloc_min Location of vertex with min angle.
   *  @param[out] iloc_max Location of vertex with max angle.
   *  @param[out] num_angle Number of cos angles computed (not skipped).
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE, 
            typename NTYPE0, typename NTYPE1, typename MTYPE,
            typename ITYPE0, typename ITYPE1,
            typename COS_TYPE0, typename COS_TYPE1>
  void compute_cos_min_max_polygon_angles
  (const DTYPE dimension, const CTYPE vertex_coord[],
   const VTYPE poly_vert[], const NTYPE0 num_poly_vert,
   MTYPE max_small_magnitude,
   COS_TYPE0 & cos_min, COS_TYPE1 & cos_max,
   ITYPE0 & iloc_min, ITYPE1 & iloc_max, NTYPE1 & num_angle)
  {
    COS_TYPE0 cos_angle;

    cos_min = -1;
    cos_max = 1;
    iloc_min = 0;
    iloc_max = 0;
    num_angle = 0;

    if (num_poly_vert < 3) {
      // Ignore degenerate polygons with 1 or 2 vertices.
      return;
    }

    for (NTYPE0 i = 0; i < num_poly_vert; i++) {
      bool flag_identical;
      compute_cos_angle_at_polygon_vertex
        (dimension, vertex_coord, poly_vert, num_poly_vert,
         i, max_small_magnitude, cos_angle, flag_identical);

      if (!flag_identical) {
        num_angle++;
        if (cos_angle > cos_min) { 
          cos_min = cos_angle; 
          iloc_min = i;
        }
        if (cos_angle < cos_max) { 
          cos_max = cos_angle; 
          iloc_max = i;
        }
      }
    }
  }


  /*!
   *  @brief Compute cosine of min/max of polygon angles of polygon.
   *  - Version that returns only cos_min, cos_max and num_angle.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE, 
            typename NTYPE0, typename NTYPE1, typename MTYPE, 
            typename COS_TYPE0, typename COS_TYPE1>
  void compute_cos_min_max_polygon_angles
  (const DTYPE dimension, const CTYPE vertex_coord[],
   const VTYPE poly_vert[], const NTYPE0 num_poly_vert,
   MTYPE max_small_magnitude,
   COS_TYPE0 & cos_min, COS_TYPE1 & cos_max, NTYPE1 & num_angle)
  {
    NTYPE0 iloc_min, iloc_max;

    compute_cos_min_max_polygon_angles
      (dimension, vertex_coord, poly_vert, num_poly_vert,
       max_small_magnitude, cos_min, cos_max, 
       iloc_min, iloc_max, num_angle);
  }


  /*!
   *  @brief Compute cosine of min/max of polygon angles of polygon.
   *  - Version using list of pointers to vertex coordinates.
   *  - If two adjacent polygon vertices have (near) identical coordinates,
   *    the angle at those vertices is skipped.
   *  @param[out] cos_min Cosine of minimum polygon angle.
   *     - Note: Min angle has maximimum cosine.
   *  @param[out] cos_max Cosine of maximum polygon angle.
   *     - Note: Max angle has minimum cosine.
   *  @param[out] iloc_min Location of vertex with min angle.
   *  @param[out] iloc_max Location of vertex with max angle.
   *  @param[out] num_angle Number of cos angles computed (not skipped).
   */
  template <typename DTYPE, typename CTYPE,
            typename NTYPE0, typename NTYPE1, typename MTYPE,
            typename ITYPE0, typename ITYPE1,
            typename COS_TYPE0, typename COS_TYPE1>
  void compute_cos_min_max_polygon_angles
  (const DTYPE dimension, const NTYPE0 num_poly_vert,
   const CTYPE * const vcoord_ptr[], MTYPE max_small_magnitude,
   COS_TYPE0 & cos_min, COS_TYPE1 & cos_max,
   ITYPE0 & iloc_min, ITYPE1 & iloc_max, NTYPE1 & num_angle)
  {
    COS_TYPE0 cos_angle;

    cos_min = -1;
    cos_max = 1;
    iloc_min = 0;
    iloc_max = 0;
    num_angle = 0;

    if (num_poly_vert < 3) {
      // Ignore degenerate polygons with 1 or 2 vertices.
      return;
    }

    for (NTYPE0 i = 0; i < num_poly_vert; i++) {
      bool flag_identical;
      compute_cos_angle_at_polygon_vertex
        (dimension, num_poly_vert, vcoord_ptr,
         i, max_small_magnitude, cos_angle, flag_identical);

      if (!flag_identical) {
        num_angle++;
        if (cos_angle > cos_min) { 
          cos_min = cos_angle; 
          iloc_min = i;
        }
        if (cos_angle < cos_max) { 
          cos_max = cos_angle; 
          iloc_max = i;
        }
      }
    }
  }


  // *** SHOULD CHANGE TO min_and_max_ ...
  /*!
   *  @brief Compute min/max polygon angles (in degrees) of polgyon.
   *  - If two adjacent polygon vertices have (near) identical coordinates,
   *    the angle at those vertices is skipped.
   *  @param[out] num_angle Number of angles computed (not skipped).
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE, 
            typename NTYPE0, typename NTYPE1, typename MTYPE,
            typename ATYPE0, typename ATYPE1>
  void compute_min_max_polygon_angles
  (const DTYPE dimension, const CTYPE vertex_coord[],
   const VTYPE poly_vert[], const NTYPE0 num_poly_vert,
   const MTYPE max_small_magnitude,
   ATYPE0 & min_angle, ATYPE1 & max_angle, NTYPE1 & num_angle)
  {
    ATYPE0 cos_min, cos_max;

    compute_cos_min_max_polygon_angles
      (dimension, vertex_coord, poly_vert, num_poly_vert,
       max_small_magnitude, cos_min, cos_max, num_angle);
    
    min_angle = std::acos(cos_min) * 180.0/std::numbers::pi;
    max_angle = std::acos(cos_max) * 180.0/std::numbers::pi;
  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle
   *    in the triangulation of a hexahedron.
   *  - Hexadhedron is split into six triangles formed from coordX[]
   *    and each of the six hexahedron edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the six triangles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param coord4[] Input coordinates.
   *  @param coord5[] Input coordinates.
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
            typename CTYPE3, typename CTYPE4, typename CTYPE5,
            typename CTYPEX, typename COS_TYPE,
            typename MTYPE>
  void compute_cos_min_hexagon_tri6_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3, const CTYPE4 * coord4,
   const CTYPE5 * coord5, const CTYPEX * coordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_HEX = 6;
    COS_TYPE cos[NUM_VERT_PER_HEX];
    bool flag_tri_zero[NUM_VERT_PER_HEX];

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
      (dimension, coord3, coord4, coordX, max_small_magnitude, 
       cos[3], flag_tri_zero[3]);
    compute_cos_min_triangle_angle
      (dimension, coord4, coord5, coordX, max_small_magnitude, 
       cos[4], flag_tri_zero[4]);
    compute_cos_min_triangle_angle
      (dimension, coord5, coord0, coordX, max_small_magnitude, 
       cos[5], flag_tri_zero[5]);

    bool is_cos_min_angle_set = false;
    for (int i = 0; i < NUM_VERT_PER_HEX; i++) {

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
   *    in the triangulation of a hexahedron.
   *  - Hexadhedron is split into six triangles formed 
   *    from coordX0[] or coordX1[] or coordX2[] and each of 
   *    the six hexahedron edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the six triangles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param coord4[] Input coordinates.
   *  @param coord5[] Input coordinates.
   *  @param coordX0[] Input coordinates.
   *  @param coordX1[] Input coordinates.
   *  @param coordX2[] Input coordinates.
   *  @param max_small_magnitude Vectors with magnitude less than or
   *           equal to max_small_magnitude are set to 0.
   *  @pre max_small_magnitude >= 0.
   *  @param[out] cos_min_angle Cosine of min triange angle.
   *  @param[out] icoord Index (0, 1 or 2) of coordinate which maximizes
   *        min angle.
   *  @param[out] flag_zero True, if some triangle has two or three edges
   *    with length less than or equal to max_small_magnitude.
   */
  template <typename DTYPE, 
            typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, typename CTYPE5,
            typename CTYPEX0, typename CTYPEX1, typename CTYPEX2,
            typename COS_TYPE, typename ITYPE, typename MTYPE>
  void compute_cos_max_min_hex_tri6_angle_3coord
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3, const CTYPE4 * coord4,
   const CTYPE5 * coord5, 
   const CTYPEX0 * coordX0, const CTYPEX1 * coordX1, const CTYPEX2 * coordX2,
   const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, ITYPE & icoord, bool & flag_zero)
  {
    COS_TYPE cos_min_angle0, cos_min_angle1, cos_min_angle2;
    bool flag_zero0, flag_zero1, flag_zero2;

    compute_cos_min_hex_tri6_angle
      (dimension, coord0, coord1, coord2, coord3, coord4, coord5, 
       coordX0, max_small_magnitude, cos_min_angle0, flag_zero0);
    compute_cos_min_hex_tri6_angle
      (dimension, coord0, coord1, coord2, coord3, coord4, coord5, 
       coordX1, max_small_magnitude, cos_min_angle1, flag_zero1);
    compute_cos_min_hex_tri6_angle
      (dimension, coord0, coord1, coord2, coord3, coord4, coord5, 
       coordX2, max_small_magnitude, cos_min_angle2, flag_zero2);

    flag_zero = (flag_zero0 && flag_zero1 && flag_zero2);
    
    if (!flag_zero0) {
      cos_min_angle = cos_min_angle0;
      icoord = 0;

      if (!flag_zero1) {
        if (cos_min_angle1 <= cos_min_angle) {
          cos_min_angle = cos_min_angle1;
          icoord = 1;
        }
      }

      if (!flag_zero2) {
        if (cos_min_angle2 <= cos_min_angle) {
          cos_min_angle = cos_min_angle2;
          icoord = 2;
        }
      }
    }
    else if (!flag_zero1) {
      cos_min_angle = cos_min_angle1;
      icoord = 1;

      if (!flag_zero2) {
        if (cos_min_angle2 <= cos_min_angle) {
          cos_min_angle = cos_min_angle2;
          icoord = 2;
        }
      }
    }
    else if (!flag_zero2) {
      cos_min_angle = cos_min_angle2;
      icoord = 2;
    }
    else {
      // flag_zero0 && flag_zero1 && flag_zero2
      cos_min_angle = 0;
      icoord = 0;
    }

  }


  /*!
   *  @brief Compute the cosine of the smallest triangle angle
   *    in the triangulation of a septagon.
   *  - Septagon is split into seven triangles formed from coordX[]
   *    and each of the seven septagon edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the seven triangles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param coord4[] Input coordinates.
   *  @param coord5[] Input coordinates.
   *  @param coord6[] Input coordinates.
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
            typename CTYPE3, typename CTYPE4, typename CTYPE5,
            typename CTYPE6, typename CTYPEX, typename COS_TYPE,
            typename MTYPE>
  void compute_cos_min_septagon_tri7_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3, const CTYPE4 * coord4,
   const CTYPE5 * coord5, const CTYPE6 * coord6, const CTYPEX * coordX,
   const MTYPE max_small_magnitude, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_SEPTAGON(7);
    COS_TYPE cos[NUM_VERT_PER_SEPTAGON];
    bool flag_tri_zero[NUM_VERT_PER_SEPTAGON];

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
      (dimension, coord3, coord4, coordX, max_small_magnitude, 
       cos[3], flag_tri_zero[3]);
    compute_cos_min_triangle_angle
      (dimension, coord4, coord5, coordX, max_small_magnitude, 
       cos[4], flag_tri_zero[4]);
    compute_cos_min_triangle_angle
      (dimension, coord5, coord6, coordX, max_small_magnitude, 
       cos[5], flag_tri_zero[5]);
    compute_cos_min_triangle_angle
      (dimension, coord6, coord0, coordX, max_small_magnitude, 
       cos[6], flag_tri_zero[6]);

    bool is_cos_min_angle_set = false;
    for (int i = 0; i < NUM_VERT_PER_SEPTAGON; i++) {

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
   *    in the triangulation of an octagon.
   *  - Octagon is split into eight triangles formed from coordX[]
   *    and each of the eight octagon edges.
   *  - Note: Cosine of the smallest angle is the largest cosine 
   *    of the eight triangles.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param coord2[] Input coordinates.
   *  @param coord3[] Input coordinates.
   *  @param coord4[] Input coordinates.
   *  @param coord5[] Input coordinates.
   *  @param coord6[] Input coordinates.
   *  @param coord7[] Input coordinates.
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
            typename CTYPE3, typename CTYPE4, typename CTYPE5,
            typename CTYPE6, typename CTYPE7, typename CTYPEX, 
            typename COS_TYPE, typename MTYPE>
  void compute_cos_min_octagon_tri8_angle
  (const DTYPE dimension, const CTYPE0 * coord0, const CTYPE1 * coord1, 
   const CTYPE2 * coord2, const CTYPE3 * coord3, const CTYPE4 * coord4,
   const CTYPE5 * coord5, const CTYPE6 * coord6, const CTYPE7 * coord7,
   const CTYPEX * coordX, const MTYPE max_small_magnitude, 
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    const int NUM_VERT_PER_OCTAGON(8);
    COS_TYPE cos[NUM_VERT_PER_OCTAGON];
    bool flag_tri_zero[NUM_VERT_PER_OCTAGON];

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
      (dimension, coord3, coord4, coordX, max_small_magnitude, 
       cos[3], flag_tri_zero[3]);
    compute_cos_min_triangle_angle
      (dimension, coord4, coord5, coordX, max_small_magnitude, 
       cos[4], flag_tri_zero[4]);
    compute_cos_min_triangle_angle
      (dimension, coord5, coord6, coordX, max_small_magnitude, 
       cos[5], flag_tri_zero[5]);
    compute_cos_min_triangle_angle
      (dimension, coord6, coord7, coordX, max_small_magnitude, 
       cos[6], flag_tri_zero[6]);
    compute_cos_min_triangle_angle
      (dimension, coord7, coord0, coordX, max_small_magnitude, 
       cos[7], flag_tri_zero[0]);

    bool is_cos_min_angle_set = false;
    for (int i = 0; i < NUM_VERT_PER_OCTAGON; i++) {

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

  //@}


  // **************************************************
  //! @name Compute dihedral angles
  // **************************************************

  //@{

  /*!
   *  @brief Compute cosine of the dihedral angle between triangle (iv0, iv1, iv2)
   *    and (iv0, iv1, iv3).
   *  @param[out] flag_collinear True if (iv0, iv1, iv2) or (iv0, iv1, iv3)
   *    are (near) collinear.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3,
            typename COS_TYPE>
  void compute_cos_dihedral_angle
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const VTYPE0 iv0, const VTYPE1 iv1, const VTYPE2 iv2, const VTYPE3 iv3,
   COS_TYPE & cos_angle, bool & flag_fail)
  {
    const CTYPE * vertex0_coord = vertex_coord+dimension*iv0;
    const CTYPE * vertex1_coord = vertex_coord+dimension*iv1;
    const CTYPE * vertex2_coord = vertex_coord+dimension*iv2;
    const CTYPE * vertex3_coord = vertex_coord+dimension*iv3;
    CTYPE magnitude;
    IJK::ARRAY<CTYPE> w01(dimension);
    IJK::ARRAY<CTYPE> w02(dimension);
    IJK::ARRAY<CTYPE> w03(dimension);
    IJK::ARRAY<CTYPE> u02(dimension);
    IJK::ARRAY<CTYPE> u03(dimension);

    cos_angle = 1.0;

    compute_unit_vector
      (dimension, vertex0_coord, vertex1_coord, 0.0, 
       w01.Ptr(), magnitude, flag_fail);
    if (flag_fail) { return; }
    subtract_coord(dimension, vertex2_coord, vertex0_coord, w02.Ptr());
    subtract_coord(dimension, vertex3_coord, vertex0_coord, w03.Ptr());

    compute_orthogonal_vector(dimension, w02.Ptr(), w01.Ptr(), u02.Ptr());
    compute_orthogonal_vector(dimension, w03.Ptr(), w01.Ptr(), u03.Ptr());

    normalize_vector
      (dimension, u02.Ptr(), 0.0, u02.Ptr(), magnitude, flag_fail);
    if (flag_fail) { return; }
    normalize_vector
      (dimension, u03.Ptr(), 0.0, u03.Ptr(), magnitude, flag_fail);
    if (flag_fail) { return; }

    compute_inner_product(dimension, u02.Ptr(), u03.Ptr(), cos_angle);

    // Bound cos_angle to [-1,1].
    // Numerical error could cause cos_angle to be outside bounds.
    if (cos_angle < -1) { cos_angle = -1; }
    if (cos_angle > 1) { cos_angle = 1; }
  }


  /// Compute cosine of min/max dihedral angles of a tetrahedron.
  template <typename DTYPE, typename VTYPE, typename CTYPE,
            typename COS_TYPE0, typename COS_TYPE1, typename NUM_TYPE>
  void compute_cos_min_max_tetrahedron_dihedral_angles
  (const DTYPE dimension, const VTYPE tetrahedra_vert[], 
   const CTYPE * vertex_coord, COS_TYPE0 & cos_min, COS_TYPE1 & cos_max, 
   NUM_TYPE & num_angle)
  {
    const NUM_TYPE NUM_VERT_PER_TETRAHEDRON(4);

    cos_min = -1;
    cos_max = 1;
    num_angle = 0;

    for (int i0 = 0; i0 < NUM_VERT_PER_TETRAHEDRON; i0++) {
      const int iv0 = tetrahedra_vert[i0];

      for (int i1 = i0+1; i1 < NUM_VERT_PER_TETRAHEDRON; i1++) {
        const int iv1 = tetrahedra_vert[i1];

        int i2 = 0;
        while (i2 == i0 || i2 == i1) { i2++; }

        int i3 = 0;
        while (i3 == i0 || i3 == i1 || i3 == i2) { i3++; }

        const int iv2 = tetrahedra_vert[i2];
        const int iv3 = tetrahedra_vert[i3];

        // Compute dihedral angle between triangle (iv0, iv1, iv2) and
        //   triangle (iv0, iv1, iv3).

        double cos_angle;
        bool flag_collinear;
        compute_cos_dihedral_angle
          (dimension, vertex_coord, iv0, iv1, iv2, iv3, 
           cos_angle, flag_collinear);

        if (!flag_collinear) {
          num_angle++;
          if (cos_angle > cos_min) { cos_min = cos_angle; }
          if (cos_angle < cos_max) { cos_max = cos_angle; }
        }
      }
    }

  }

  //@}


  // **************************************************
  //! @name Distance/Length operations
  // **************************************************

  //@{

  /// Compute square of distance between two points.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename DIST_TYPE>
  void compute_distance_squared
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance_squared)
  {
    distance_squared = 0.0;
    for (DTYPE d = 0; d < dimension; d++) {
      CTYPE0 diff = coord0[d] - coord1[d];
      distance_squared += diff*diff;
    }
  }

  /// Compute distance between two points.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename DIST_TYPE>
  void compute_distance
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance)
  {
    compute_distance_squared(dimension, coord0, coord1, distance);
    distance = std::sqrt(distance);
  }


  /*!
   *  @brief Compute the distance between two iv0 and iv1.
   * - Vertex coordinates are coord[iv0*dimension] and coord[iv1*dimension].
   */
  template <typename DTYPE, typename CTYPE, 
            typename ITYPE0, typename ITYPE1, typename DIST_TYPE>
  void compute_distance_coord_list
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const ITYPE0 iv0, const ITYPE1 iv1, DIST_TYPE & distance)
  {
    const CTYPE * coord0 = vertex_coord+dimension*iv0;
    const CTYPE * coord1 = vertex_coord+dimension*iv1;

    compute_distance(dimension, coord0, coord1, distance);
  }


  /*!
   *  @brief Compute the length of the k'th polygon edge.
   * @param vertex_list[] Array of polygon vertices.
   * @param num_poly_vertices Number of polygon vertices.
   * @param[out] length Length of k'th polygon edge.
   */
  template <typename DTYPE, typename CTYPE, typename VTYPE,
            typename NTYPE, typename ITYPE, typename LTYPE>
  void compute_length_of_polygon_edge
  (const DTYPE dimension, const CTYPE * vertex_coord,
   const VTYPE * vertex_list, const NTYPE num_poly_vertices,
   const ITYPE k0, LTYPE & length)
  {
    if (num_poly_vertices == 0) {
      length = 0;
      return;
    }

    const ITYPE k1 = (k0+1)%num_poly_vertices;

    compute_distance_coord_list
      (dimension, vertex_coord, vertex_list[k0], vertex_list[k1], length);

    return;
  }


  /*!
   *  @brief Compute square of scaled distance between two points.
   *  Scale by coordinate-wise divide.
   *  @pre scale[d]>0 for all d in range [0,dimension-1].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename SCALE_TYPE, typename DIST_TYPE>
  void compute_scaled_distance_squared_divide
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   // *** SHOULD THIS BE "const SCALE_TYPE scale[]" ??? ***
   const SCALE_TYPE scale, DIST_TYPE & scaled_distance_squared)
  {
    scaled_distance_squared = 0.0;
    for (DTYPE d = 0; d < dimension; d++) {
      CTYPE0 diff = (coord0[d] - coord1[d])/scale[d];
      scaled_distance_squared += diff*diff;
    }
  }


  /*!
   *  @brief Compute scaled distance between two points.
   *  - Scale by coordinate-wise divide.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename SCALE_TYPE, typename DIST_TYPE>
  void compute_scaled_distance_divide
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const SCALE_TYPE scale, DIST_TYPE & scaled_distance)
  {
    compute_scaled_distance_squared_divide
      (dimension, coord0, coord1, scale, scaled_distance);
    scaled_distance = std::sqrt(scaled_distance);
  }

  /// Compute Linf distance between two points.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename DIST_TYPE, typename AXIS_TYPE>
  void compute_Linf_distance
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance, AXIS_TYPE & axis)
  {
    DIST_TYPE x;

    distance = 0;
    axis = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] < coord1[d])
        { x = coord1[d] - coord0[d]; }
      else
        { x = coord0[d] - coord1[d]; }

      if (x > distance) {
        distance = x; 
        axis = d;
      }
    }
  }

  /// Compute Linf distance between two points.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename DIST_TYPE>
  void compute_Linf_distance
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance)
  {
    DIST_TYPE x;

    distance = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] < coord1[d])
        { x = coord1[d] - coord0[d]; }
      else
        { x = coord0[d] - coord1[d]; }

      if (x > distance) { distance = x; }
    }
  }


  /*! 
   *  @brief Compute signed (L2) distance from coord0[] to hyperplane through coord1[].
   *  @tparam DIST_TYPE Distance type.  DIST_TYPE should be a signed type.
   *  @param dimension Coordinate and vector dimensions.
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Hyperplane passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to hyperplane.
   *  @param[out] distance Signed distance to hyperplane.
   *    - distance = inner product of (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_signed_distance_to_hyperplane
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    // Compute inner product of (coord0[]-coord1[]) and orth_dir[].
    distance = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { distance = distance + (coord0[d]-coord1[d])*orth_dir[d]; }
  }


  /*!
   *  @brief Compute (unsigned L2) distance from coord0[] 
   *     to hyperplane through coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Hyperplane passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to hyperplane.
   *  @param[out] distance Distance to hyperplane.
   *     - distance = the absolute value of inner product of 
   *                 (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_to_hyperplane
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    VTYPE signed_distance;

    compute_signed_distance_to_hyperplane
      (dimension, coord0, coord1, orth_dir, signed_distance);

    if (signed_distance >= 0) 
      { distance = signed_distance; }
    else
      { distance = -signed_distance; }
  }


  /*!
   *  @brief Compute (L2) distance squared from coord0[] to line through coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param dir[] Line direction.
   *  @param[out] distance_squared
   *    Distance squared from coord0[] to line through coord1[].
   *    - distance_squared = Magnitude squared of w where w is the component
   *       of (coord0[]-coord1[]) orthogonal to dir[].
   *    - w = (coord0[]-coord1[]) - x*dir[] where x is the inner product
   *       of (coord0[]-coord1[]) and dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_squared_to_line
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE dir[], DIST_TYPE & distance_squared)
  {
    // Compute inner product of (coord0[]-coord1[]) and dir[].
    DIST_TYPE x = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { x = x + (coord0[d]-coord1[d])*dir[d]; }

    // Compute magnitude squared of ((coord0[]-coord1[])-x*dir[]).
    distance_squared = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      DIST_TYPE y = (coord0[d]-coord1[d]) - x*dir[d];
      distance_squared = distance_squared + (y*y);
    }
  }


  /*!
   *  @brief Compute (unsigned, L2) distance from coord0[] to line through coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param dir[] Line direction.
   *  @param[out] distance Distance from coord0[] to line through coord1[].
   *    - distance = Magnitude of w where w is the component
   *       of (coord0[]-coord1[]) orthogonal to dir[].
   *    - w = (coord0[]-coord1[]) - x*dir[] where x is the inner product
   *       of (coord0[]-coord1[]) and dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_to_line
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE dir[], DIST_TYPE & distance)
  {
    compute_distance_squared_to_line
      (dimension, coord0, coord1, dir, distance);
    distance = std::sqrt(distance);
  }

  //@}


  // **************************************************
  //! @name Intersection operations
  // **************************************************  

  //@{

  /*!
   *  @brief Compute the intersection of a line and a hyperplane.
   *  - Line contains \a coord0[] and has direction \a line_dir[].
   *  - Hyperplane contains \a coord1[] 
   *    and has orthogonal direction \a orth_dir[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Coordinates of point contained by line. 
   *  @param line_dir[] Line direction.
   *  @param coord1[] Coordinates of point contained by hyperplane.
   *  @param orth_dir[] Direction orthgonal to \a line_dir[].
   *  @param min_abs_cos_angle_LO Minimum value of the absolute value
   *      of the cosine of the angle between \a line_dir[] and \a orth_dir[].
   *  @param[out] intersection_point[] 
   *                Intersection point of line and hyperplane.
   *  @param[out] flag_succeeded True if computation succeeded.
   *     Computations succeeds if:
   *       abs(inner_product(line_dir[],orth_dir[])) > min_abs_cos_angle_LO.
   *  @pre \a line_dir[] is a unit vector.
   *  @pre \a orth_dir[] is a unit vector.
   *  @pre \a min_abs_cos_angle_LO >= 0.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename VTYPE0, typename VTYPE1, typename MTYPE,
            typename CTYPE2>
  void intersect_line_hyperplane
  (const DTYPE dimension, const CTYPE0 coord0[], const VTYPE0 line_dir[],
   const CTYPE1 coord1[], const VTYPE1 orth_dir[],
   const MTYPE min_abs_cos_angle_LO,
   CTYPE2 intersection_point[], bool & flag_succeeded)
  {
    VTYPE0 cos_angle;
    compute_inner_product(dimension, line_dir, orth_dir, cos_angle);

    if (std::abs(cos_angle) > min_abs_cos_angle_LO) {
      VTYPE0 distance, t;
      flag_succeeded = true;

      compute_signed_distance_to_hyperplane
        (dimension, coord0, coord1, orth_dir, distance);
      t = -distance/cos_angle;
      add_scaled_coord(dimension, t, line_dir, coord0, intersection_point);
    }
    else {
      flag_succeeded = false;
      IJK::copy_coord(dimension, coord1, intersection_point);
    }
  }

  //@}


  // **************************************************
  //! @name Average operations
  // **************************************************  

  //@{

  /*!
   *  @brief Compute the weighted average of three coordinates.
   *  - Compute w0*coord0[]+w1*coord1[]+w2*coord2[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param w0 Weight of coord0[].
   *  @param coord0[] Input coordinates.
   *  @param w1 Weight of coord1[].
   *  @param coord1[] Input coordinates.
   *  @param w2 Weight of coord2[].
   *  @param coord2[] Input coordinates.
   *  @param[out] coord3 Output coordinate.
   *    - coord2[] = \a w0* \a coord[0] + \a w1*coord1[] + \a w2*coord2[].
   */
  template <typename DTYPE, 
            typename WTYPE0, typename WTYPE1, typename WTYPE2,
            typename CTYPE0, typename CTYPE1, typename CTYPE2, typename CTYPE3>
  void compute_weighted_average_3coord
  (const DTYPE dimension, 
   const WTYPE0 w0, const CTYPE0 coord0[],
   const WTYPE1 w1, const CTYPE1 coord1[], 
   const WTYPE2 w2, const CTYPE2 coord2[],
   CTYPE3 coord3[])
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { coord3[d] = w0*coord0[d] + w1*coord1[d] + w2*coord2[d]; };
  }
  
  //@}


  // **********************************************************
  //! @name Compute centroid.
  // **********************************************************

  //@{

  /*!
   *  @brief Compute centroid of a list of points.
   *  @param point_list[] List of points.  Indices into array coord[].
   *  @param coord[] List of point coordinates.
   *    i'th point has coordinates starting at coord[point_list[i]*dimension].
   *  @param[out] centroid_coord[] Coordinates of centroid.
   *  @pre Array centroid_coord[] is preallocated to length at least dimension.
   */
  template <typename DTYPE, typename PTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE>
  void compute_coord_centroid
  (const DTYPE dimension, const PTYPE point_list[], const NTYPE num_points,
   const CTYPE0 coord[], CTYPE1 centroid_coord[])
  {
    IJK::set_coord(dimension, 0, centroid_coord);
    
    for (NTYPE i = 0; i < num_points; i++) {
      PTYPE ipoint = point_list[i];
      IJK::add_coord
        (dimension, coord+ipoint*dimension, centroid_coord, centroid_coord);
    }

    if (num_points > 0) {
      IJK::divide_coord(dimension, num_points, centroid_coord, centroid_coord);
    }
  }


  /*!
   *  @brief Compute centroid of a list of points.
   *  Version using C++ vector of coordinates.
   */
  template <typename DTYPE, typename PTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE>
  void compute_coord_centroid
  (const DTYPE dimension, const PTYPE point_list[], const NTYPE num_points,
   const std::vector<CTYPE0> & coord, CTYPE1 centroid_coord[])
  {
    compute_coord_centroid
      (dimension, point_list, num_points, IJK::vector2pointer(coord),
       centroid_coord);
  }


  /// Compute centroid of quad.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4>
  void compute_quad_centroid
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const CTYPE2 coord2[], const CTYPE3 coord3[], 
   CTYPE4 centroid_coord[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      centroid_coord[d] = 
        (coord0[d] + coord1[d] + coord2[d] + coord3[d])/4.0; 
    }
  }


  /// Compute centroid of quad.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  void compute_quad_centroid
  (const DTYPE dimension, const CTYPE0 * quad_coord[4],
   CTYPE1 centroid_coord[])
  {
    compute_quad_centroid
      (dimension, quad_coord[0], quad_coord[1], 
       quad_coord[2], quad_coord[3], centroid_coord);
  }


  /// Compute centroid of quad.
  template <typename DTYPE, typename VTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1>
  void compute_quad_centroid
  (const DTYPE dimension, const VTYPE quad_vert[], 
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 centroid_coord[])
  {
    const int NUM_VERT_PER_QUAD(4);

    compute_coord_centroid
      (dimension, quad_vert, NUM_VERT_PER_QUAD, coord_array, centroid_coord);
  }


  /// Compute centroid of a pentagon.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2, 
            typename CTYPE3, typename CTYPE4, typename CTYPE5>
  void compute_pentagon_centroid
  (const DTYPE dimension, const CTYPE0 coord0[], const CTYPE1 coord1[],
   const CTYPE2 coord2[], const CTYPE3 coord3[], const CTYPE4 coord4[],
   CTYPE5 centroid_coord[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      centroid_coord[d] = 
        (coord0[d] + coord1[d] + coord2[d] + coord3[d] + coord4[d])/5.0; 
    };
  }


  /// Compute centroid of a pentagon.
  template <typename DTYPE, typename VTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1>
  void compute_pentagon_centroid
  (const DTYPE dimension, const VTYPE pentagon_vert[], 
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 centroid_coord[])
  {
    const int NUM_VERT_PER_PENTAGON(5);

    compute_coord_centroid
      (dimension, pentagon_vert, NUM_VERT_PER_PENTAGON, coord_array, 
       centroid_coord);
  }


  /// Compute centroid of hexahedron.
  template <typename DTYPE, typename VTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1>
  void compute_hexahedron_centroid
  (const DTYPE dimension, const VTYPE hexahedron_vert[], 
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 centroid_coord[])
  {
    const int NUM_VERT_PER_HEXAHEDRON(8);

    compute_coord_centroid
      (dimension, hexahedron_vert, NUM_VERT_PER_HEXAHEDRON, coord_array, 
       centroid_coord);
  }

  //@}


  // **************************************************
  //! @name Comparison operators
  // **************************************************

  //@{

  /*!
   *  @brief Return true if coordinates are the same.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool is_coord_equal(const DTYPE dimension, 
                      const CTYPE0 coord0[], const CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { if (coord0[d] != coord1[d]) { return(false); }; }

    return(true);
  }


  /*!
   *  @brief Return true if some coordinate is the same.
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool is_some_coord_equal(const DTYPE dimension, 
                           const CTYPE0 coord0[], const CTYPE1 coord1[])
  {
    for (DTYPE d = 0; d < dimension; d++)
      { if (coord0[d] == coord1[d]) { return(true); }; }

    return(false);
  }


  /*!
   *  @brief Return true if all elements of coord0[] are less than or equal
   *    to all elements of coord1[].
   */
  // *** SHOULD USE coord0[] and coord1[] IN PARAMETER LIST? ***
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool is_coord_less_than_or_equal
  (const DTYPE dimension, const CTYPE0 coord0, const CTYPE1 coord1)
  {
    for (DTYPE d = 0; d < dimension; d++)
      { if (coord0[d] > coord1[d]) { return(false); }; };

    return(true);
  }


  /*!
   *  @brief Return true if coord0[] is contained in the rectangular region
   *    with opposing corners coord1[] and coord2[].
   */
  template <typename DTYPE, typename CTYPE0, 
            typename CTYPE1, typename CTYPE2>
  bool is_coord_in_rectangular_region
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], const CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] < coord1[d]) {
        if (coord0[d] < coord2[d]) { return(false); }
      }
      else if (coord0[d] > coord1[d]) {
        if (coord0[d] > coord2[d]) { return(false); }
      }
    }

    return(true);
  }


  /*!
   *  @brief Return true if coord0[] is contained in the interior 
   *    of the rectangular region with opposing corners coord1[] and coord2[].
   */
  template <typename DTYPE, typename CTYPE0, 
            typename CTYPE1, typename CTYPE2>
  bool is_coord_in_interior_of_rectangular_region
  (const DTYPE dimension, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], const CTYPE2 coord2[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] < coord1[d]) {
        if (coord0[d] <= coord2[d]) { return(false); }
      }
      else if (coord0[d] > coord1[d]) {
        if (coord0[d] >= coord2[d]) { return(false); }
      }
      else {
        // coord0[d] == coord1[d]
        return(false); 
      }
    }

    return(true);
  }


  /*!
   *  @brief Return true if scalar s0 is contained in the range [s1,s2],
   *    i.e. Return true if (s1 <= s0 <= s2) or if (s1 >= s0 >= s2).
   */
  template <typename STYPE0, typename STYPE1, typename STYPE2>
  bool is_scalar_in_range
  (const STYPE0 s0, const STYPE1 s1, const STYPE2 s2)
  {
    if (s0 < s1) {
      if (s0 < s2) { return(false); }
    }
    else if (s0 > s1) {
      if (s0 > s2) { return(false); }
    }

    return(true);
  }


  /*!
   *  @brief Return true if coord0[d] is in range [coord1[d], coord2[d]].
   *  Note: coord1[d] is not necessarily less than coord2[d].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  bool is_dth_coord_in_range
  (const DTYPE d, CTYPE0 coord0[], const CTYPE1 coord1[], 
   const CTYPE2 coord2[])
  {
    return(is_scalar_in_range(coord0[d], coord1[d], coord2[d]));
  }


  /*
   *  @brief Compare coord0[d] and coord1[d].
   *  - Return -1 if coord0[d] < coord1[d].
   *  - Return 1 if coord0[d] > coord1[d].
   *  - Return 0 if coord0[d] == coord1[d].
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  inline int compare_dth_coord
  (const CTYPE0 coord0[], const CTYPE1 coord1[], const DTYPE d)
  {
    if (coord0[d] < coord1[d]) {
      return -1;
    }
    else if (coord0[d] > coord1[d]) {
      return 1;
    }
    else {
      // coord0[d] = coord0[d].
        return 0;
    }
  }

  //@}


  // **************************************************
  //! @name Rounding operators
  // **************************************************

  //@{

  /// Round x to nearest integer.  
  template <typename T>
  T round(const T x)
  {
    return x >= 0.0f ? floorf(x + 0.5f) : ceilf(x - 0.5f);
  }
  

  /// Round single coordinate
  template <int n, typename T>
  T round_coord(const T x)
  {
    T y = round(n*x);
    return(y/n);
  }

  /// Round single coordinate
  template <typename ITYPE, typename T>
  T round_coord(const ITYPE n, const T x)
  {
    T y = round(n*x);
    return(y/n);
  }

  /// Round coordinates.
  template <int n, typename DTYPE, typename T0, typename T1>
  void round_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round_coord<n>(v0[d]); }
  }

  /// Round coordinates.
  template <typename ITYPE, typename DTYPE, typename T0, typename T1>
  void round_coord
  (const ITYPE n, const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round_coord(n, v0[d]); }
  }

  // Specialized rounding operators.

  /// Round single coordinate to nearest 1/16'th.
  template <typename T>
  T round16_coord(const T x) { return(round_coord<16>(x)); }

  /// Round single coordinate to nearest 1/32'th.
  template <typename T>
  T round32_coord(const T x) { return(round_coord<32>(x)); }

  /// Round single coordinate to nearest 1/64'th.
  template <typename T>
  T round64_coord(const T x) { return(round_coord<64>(x)); }

  /// Round single coordinate to nearest 1/128'th.
  template <typename T>
  T round128_coord(const T x) { return(round_coord<128>(x)); }

  /// Round single coordinate to nearest 1/256'th.
  template <typename T>
  T round256_coord(const T x) { return(round_coord<256>(x)); }

  /// Round single coordinate to nearest 1/512'th.
  template <typename T>
  T round512_coord(const T x) { return(round_coord<512>(x)); }

  /// Round single coordinate to nearest 1/1024'th.
  template <typename T>
  T round1024_coord(const T x) { return(round_coord<1024>(x)); }

  /// Round coordinates to nearest 1/16'th.
  template <typename DTYPE, typename T0, typename T1>
  void round16_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round16_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/32'th.
  template <typename DTYPE, typename T0, typename T1>
  void round32_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round32_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/64'th.
  template <typename DTYPE, typename T0, typename T1>
  void round64_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round64_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/128'th.
  template <typename DTYPE, typename T0, typename T1>
  void round128_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round128_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/256'th.
  template <typename DTYPE, typename T0, typename T1>
  void round256_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round256_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/512'th.
  template <typename DTYPE, typename T0, typename T1>
  void round512_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round512_coord(v0[d]); }
  }

  /// Round coordinates to nearest 1/1024'th.
  template <typename DTYPE, typename T0, typename T1>
  void round1024_coord(const DTYPE dimension, const T0 * v0, T1 * v1)
  {
    for (DTYPE d = 0; d < dimension; d++) 
      { v1[d] = round1024_coord(v0[d]); }
  }

  //@}  

  // **************************************************
  //! @name Clamp operator
  // **************************************************

  //@{

  /*!
   * @brief Clamp x to interval [0,1].
   * - If (x < 0), then x = 0.
   * - If (x > 1), then x = 1.
   */
  template <typename T>
  T clamp01_coord(const T x)
  {
    if (x < 0) { return(0); }
    else if (x > 1) { return(1); }

    return(x);
  }

  /*!
   *  @brief Clamp x to range [a,b].
   *  @pre a <= b.
   */
  template <typename Tx, typename Ta, typename Tb>
  Tx clamp_coord_to_range(const Tx x, const Ta a, const Tb b)
  {
    if (x < a) { return(Tx(a)); }
    else if (x > b) { return(Tx(b)); }

    return(x);
  }

  //@}

  // **************************************************
  // @name OLD DEPRECATED Comparison operators
  // **************************************************

  //@{

  /// DEPRECATED: Replaced by is_coord_equal()
  /* OBSOLETE
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool compare_coord(const DTYPE dimension, 
                     const CTYPE0 coord0[], const CTYPE1 coord1[])
  {
    return(is_coord_equal(dimension, coord0, coord1));
  }

  /// DEPRECATED: Replaced by is_coord_less_than_or_equal()
  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  bool is_less_than_or_equal_to
  (const DTYPE dimension, const CTYPE0 coord0, const CTYPE1 coord1)
  {
    return(is_coord_less_than_or_equal(dimension, coord0, coord1));
  }
  */

  //@}

  // **************************************************
  //! @name Insert Coord in List
  // **************************************************

  //@{

  template <typename DTYPE, typename CTYPE0, typename CTYPE1>
  std::size_t insert_coord
  (const DTYPE dimension, const CTYPE0 coord[], 
   std::vector<CTYPE1> & coord_list)
  {
    const std::size_t list_size = coord_list.size();
    const std::size_t ic = list_size/dimension;

    // Allocate room for new coord.
    coord_list.resize(list_size + dimension);

    copy_coord(dimension, coord, &(coord_list[list_size]));

    return(ic);
  }

  //@}


  // **************************************************
  //! @name 2D versions
  // **************************************************

  //@{

  /*!
   *  @brief Subtract \a coord1[] from \a coord0[].
   *  @param coord0 = Input coordinates.
   *  @param coord1 = Input coordinates.
   *  @param[out] coord2 = Output coordinate equal to (\a coord0[] - \a coord1[]).
   */
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void subtract_coord_2D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], CTYPE2 coord2[])
  { 
    const int DIM2(2);
    subtract_coord(DIM2, coord0, coord1, coord2); 
  }

  //@}

  // **************************************************
  //! @name 2D vector operations
  // **************************************************

  //@{

  /// Compute 2D unit vector.
  template <typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename CTYPE2, typename MTYPE1>
  void compute_unit_vector_2D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const MTYPE0 max_small_magnitude, CTYPE2 coord2[],
   MTYPE1 & magnitude, bool & flag_zero)
  {
    const int DIM2(2);
    compute_unit_vector
      (DIM2, coord0, coord1, max_small_magnitude, coord2,
       magnitude, flag_zero);
  }

  
  //@}

  
  // **************************************************
  //! @name 2D distances
  // **************************************************

  //@{
  
  /*! 
   *  @brief Compute signed (L2) distance from coord0[] to line
   *    through coord1[].
   *  @tparam DIST_TYPE Distance type.  DIST_TYPE should be a signed type.
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to line.
   *  @param[out] distance Signed distance to line.
   *    - distance = inner product of (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_signed_distance_to_line_2D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    const int DIM2(2);
    compute_signed_distance_to_hyperplane
      (DIM2, coord0, coord1, orth_dir, distance);
  }

  //@}
  
  
  // **************************************************
  //! @name 3D versions
  // **************************************************

  //@{

  /*!
   *  @brief Set all vertex coordinates to \a c.
   *  @param c = Scalar constant.  Set vertex coordinates to \a c.
   *  @param[out] coord[] = Output coordinates.
   */
  template <typename STYPE, typename CTYPE>
  void set_coord_3D(const STYPE c, CTYPE coord[])
  { 
    const int DIM3(3);
    set_coord(DIM3, c, coord); 
  }


  /*!
   *  @brief Copy \a coord0[] to \a coord1[].
   *  @param coord0 = Input coordinates.
   *  @param[out] coord1 = Output coordinates.
   */
  template <typename CTYPE0, typename CTYPE1>
  void copy_coord_3D(const CTYPE0 coord0[], CTYPE1 coord1[])
  { 
    const int DIM3(3);
    copy_coord(DIM3, coord0, coord1); 
  }


  /*!
   *  @brief Add \a coord0[] to \a coord1[].
   *  @param coord0 = Input coordinates.
   *  @param coord1 = Input coordinates.
   *  @param[out] coord2 = Output coordinate equal to (\a coord0[] + \a coord1[]).
   */
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void add_coord_3D(const CTYPE0 coord0[],
                    const CTYPE1 coord1[], CTYPE2 coord2[])
  { 
    const int DIM3(3);
    add_coord(DIM3, coord0, coord1, coord2); 
  }


  /*!
   *  @brief Subtract \a coord1[] from \a coord0[].
   *  @param coord0 = Input coordinates.
   *  @param coord1 = Input coordinates.
   *  @param[out] coord2 = Output coordinate equal to (\a coord0[] - \a coord1[]).
   */
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void subtract_coord_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], CTYPE2 coord2[])
  { 
    const int DIM3(3);
    subtract_coord(DIM3, coord0, coord1, coord2); 
  }

  /// Multiply \a coord0[] by the scalar \a s.
  template <typename STYPE, typename CTYPE0, typename CTYPE1>
  inline void multiply_coord_3D
  (const STYPE s, const CTYPE0 coord0[], CTYPE1 coord1[])
  { 
    const int DIM3(3);
    multiply_coord(DIM3, s, coord0, coord1); 
  }

  /*!
   *  Add \a s * \a coord0[] to \a coord1[]).
   *  @param s Scaling factor.
   *  @param coord0[] Input coordinates.
   *  @param coord1[] Input coordinates.
   *  @param[out] coord2[] Output coordinate.
   *     - coord2[] = (\a s * \a coord0[] + \a coord1[]).
   */
  template <typename STYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  void add_scaled_coord_3D
  (const STYPE s, 
   const CTYPE0 coord0[], const CTYPE1 coord1[], CTYPE2 coord2[])
  {
    const int DIM3(3);
    add_scaled_coord(DIM3, s, coord0, coord1, coord2);
  }

  /// Compute absolute value of the each coord.
  template <typename CTYPE0, typename CTYPE1>
  void abs_coord_3D(const CTYPE0 coord0[], CTYPE1 coord1[])
  {
    const int DIM3(3);
    abs_coord(DIM3, coord0, coord1);
  }

  /// Compute the midpoint of two coordinates.
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2>
  inline void compute_midpoint_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], CTYPE2 midpoint_coord[])
  { 
    const int DIM3(3);
    compute_midpoint(DIM3, coord0, coord1, midpoint_coord); 
  }


  /*!
   *  @brief Compute inner product of \a coord0[] and \a coord1[].
   *  @param coord0 = Input coordinates.
   *  @param coord1 = Input coordinates.
   *  @param[out] product = Inner product.
   */
  template <typename CTYPE0, typename CTYPE1, typename STYPE>
  void compute_inner_product_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], STYPE & product)
  { 
    const int DIM3(3);
    compute_inner_product(DIM3, coord0, coord1, product); 
  }

  /// Compute sum of squares of coordinates.
  template <typename CTYPE0, typename STYPE>
  void compute_sum_of_squares_3D(const CTYPE0 coord0[], STYPE & sum)
  { 
    const int DIM3(3);
    compute_sum_of_squares(DIM3, coord0, sum); 
  }


  /*!
   *  @brief Compute magnitude of coordinate vector.
   *  @param coord0 = Input coordinates.
   *  @param[out] magnitude = magnitude
   */
  template <typename CTYPE0, typename STYPE>
  void compute_magnitude_3D(const CTYPE0 coord0[], STYPE & magnitude)
  { 
    const int DIM3(3);
    compute_magnitude(DIM3, coord0, magnitude); 
  }


  /*!
   *  @brief Return true if coordinates are the same.
   *  @param coord0 = Input coordinates.
   *  @param coord1 = Input coordinates.
   */
  template <typename CTYPE0, typename CTYPE1>
  bool is_coord_equal_3D(const CTYPE0 coord0[], const CTYPE1 coord1[])
  { 
    const int DIM3(3);
    return(is_coord_equal(DIM3, coord0, coord1)); 
  }


  /*!
   *  @brief Return true if all elements of coord0[] are less than 
   *    or equal to all elements of coord1[].
   */
  template <typename CTYPE0, typename CTYPE1>
  bool is_coord_less_than_or_equal_3D
  (const CTYPE0 coord0, const CTYPE1 coord1)
  { 
    const int DIM3(3);
    return(is_coord_less_than_or_equal(DIM3, coord0, coord1)); 
  }

  template <typename CTYPE0, typename CTYPE1>
  inline std::size_t insert_coord_3D
  (const CTYPE0 coord[], std::vector<CTYPE1> & coord_list)
  {
    const int DIM3(3);
    return(insert_coord(DIM3, coord, coord_list));
  }

  //@}


  // **************************************************
  //! @name 3D vector operations
  // **************************************************

  //@{

  /// Compute 3D unit vector.
  template <typename CTYPE0, typename MTYPE0, 
            typename CTYPE1, typename CTYPE2, typename MTYPE1>
  void compute_unit_vector_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const MTYPE0 max_small_magnitude, CTYPE2 coord2[],
   MTYPE1 & magnitude, bool & flag_zero)
  {
    const int DIM3(3);
    compute_unit_vector
      (DIM3, coord0, coord1, max_small_magnitude, coord2,
       magnitude, flag_zero);
  }

  //@}


  // **************************************************
  //! @name 3D distances
  // **************************************************

  //@{

  /// Compute distance between two 3D points.
  template <typename CTYPE0, typename CTYPE1, typename DIST_TYPE>
  void compute_distance_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   DIST_TYPE & distance)
  {
    const int DIM3(3);
    compute_distance(DIM3, coord0, coord1, distance);
  }


  /// Compute square of distance between two 3D points.
  template <typename CTYPE0, typename CTYPE1, typename DIST_TYPE>
  void compute_distance_squared_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[], DIST_TYPE & distance_squared)
  {
    const int DIM3(3);
    compute_distance_squared(DIM3, coord0, coord1, distance_squared);
  }


  /*! 
   *  @brief Compute signed (L2) distance from coord0[] to plane
   *    through coord1[].
   *  @tparam DIST_TYPE Distance type.  DIST_TYPE should be a signed type.
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Plane passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to plane.
   *  @param[out] distance Signed distance to plane.
   *    - distance = inner product of (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_signed_distance_to_plane_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    const int DIM3(3);
    compute_signed_distance_to_hyperplane
      (DIM3, coord0, coord1, orth_dir, distance);
  }


  /*!
   *  @brief Compute (unsigned, L2) distance from coord0[] to plane through coord1[].
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Plane passes through coord1[].
   *  @param orth_dir[] Direction orthogonal to hyperplane.
   *  @param[out] distance Distance to hyperplane.
   *     - distance = the absolute value of inner product of 
   *                 (coord0[]-coord1[]) and orth_dir[].
   *  @pre orth_dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_to_plane_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE orth_dir[], DIST_TYPE & distance)
  {
    const int DIM3(3);
    compute_distance_to_hyperplane
      (DIM3, coord0, coord1, orth_dir, distance);
  }


  /*!
   *  @brief Compute (L2) distance squared from coord0[] to line through coord1[].
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param dir[] Line direction.
   *  @param[out] distance_squared
   *    Distance squared from coord0[] to line through coord1[].
   *    - distance_squared = Magnitude squared of w where w is the component
   *       of (coord0[]-coord1[]) orthogonal to dir[].
   *    - w = (coord0[]-coord1[]) - x*dir[] where x is the inner product
   *       of (coord0[]-coord1[]) and dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_squared_to_line_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE dir[], DIST_TYPE & distance_squared)
  {
    const int DIM3(3);
    compute_distance_squared_to_line
      (DIM3, coord0, coord1, dir, distance_squared);
  }


  /*!
   *  @brief Compute (unsigned, L2) distance from coord0[] to line through coord1[].
   *  @param dimension Coordinate dimension (= number of coordinates.)
   *  @param coord0[] Compute distance from coord0[].
   *  @param coord1[] Line passes through coord1[].
   *  @param dir[] Line direction.
   *  @param[out] distance Distance from coord0[] to line through coord1[].
   *    - distance = Magnitude of w where w is the component
   *       of (coord0[]-coord1[]) orthogonal to dir[].
   *    - w = (coord0[]-coord1[]) - x*dir[] where x is the inner product
   *       of (coord0[]-coord1[]) and dir[].
   *  @pre dir[] is a unit vector.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE, typename DIST_TYPE>
  void compute_distance_to_line_3D
  (const CTYPE0 coord0[], const CTYPE1 coord1[],
   const VTYPE dir[], DIST_TYPE & distance)
  {
    const int DIM3(3);
    compute_distance_to_line(DIM3, coord0, coord1, dir, distance);
  }

  //@}


  // **************************************************
  //! @name 2D and 3D Matrix operations
  // **************************************************

  //@{

  /// Return determinant of a 2x2 matrix.
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void determinant_2x2
  (const COORD_TYPE a00, const COORD_TYPE a01,
   const COORD_TYPE a10, const COORD_TYPE a11,
   RESULT_TYPE & result)
  {
    result = a00*a11 - a01*a10;
  }


  /*!
   *  @overload
   *  @brief Return determinant of a 2x2 matrix. (C array.)
   *  - Version passing two C arrays.
   *  @param a0[] First matrix row.
   *  @param a1[] Second matrix row.
   *  @param[out] result 2x2 determinant.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename RESULT_TYPE>
  inline void determinant_2x2
  (const CTYPE0 a0[2], const CTYPE1 a1[2],
   RESULT_TYPE & result)
  {
    determinant_2x2(a0[0], a0[1], a1[0], a1[1], result);
  }


  /// Return determinant of a 3x3 matrix.
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void determinant_3x3
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], RESULT_TYPE & result)
  {
    RESULT_TYPE D;
    determinant_2x2(p1[1], p1[2], p2[1], p2[2], D);
    result = p0[0] * D;
    determinant_2x2(p2[1], p2[2], p0[1], p0[2], D);
    result += p1[0] * D;
    determinant_2x2(p0[1], p0[2], p1[1], p1[2], D);
    result += p2[0] * D;
  }


  /// Return the square of the Froebenius norm of a 3x3 matrix.
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void compute_froebenius_norm_squared_3x3
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], RESULT_TYPE & result)
  {
    RESULT_TYPE x;
    compute_sum_of_squares_3D(p0, result);
    compute_sum_of_squares_3D(p1, x);
    result += x;
    compute_sum_of_squares_3D(p2, x);
    result += x;
  }


  /// Return the Froebenius norm of a 3x3 matrix.
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void compute_froebenius_norm_3x3
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], RESULT_TYPE & result)
  {
    compute_froebenius_norm_squared_3x3(p0, p1, p2, result);
    result = std::sqrt(result);
  }


  /// Transpose a 3x3 matrix.
  template <typename COORD_TYPE>
  void transpose_3x3
  (COORD_TYPE p0[3], COORD_TYPE p1[2], COORD_TYPE p2[3])
  {
    std::swap(p0[1], p1[0]);
    std::swap(p0[2], p2[0]);
    std::swap(p1[2], p2[1]);
  }


  /// Return the first row of the co-factor matrix of a 3x3 matrix.
  template <typename COORD_TYPEA, typename COORD_TYPEB>
  void co_factor_3x3_row1
  (const COORD_TYPEA p0[3], const COORD_TYPEA p1[2], 
   const COORD_TYPEA p2[3], COORD_TYPEB q0[3])
  {
    const int DIM3(3);

    for (int j0 = 0; j0 < DIM3; j0++) {
      const int j1 = (j0+1)%DIM3;
      const int j2 = (j0+2)%DIM3;
      determinant_2x2(p1[j1], p1[j2], p2[j1], p2[j2], q0[j0]);
    }
  }


  /// Return the co-factor matrix of a 3x3 matrix.
  template <typename COORD_TYPEA, typename COORD_TYPEB>
  void co_factor_3x3
  (const COORD_TYPEA p0[3], const COORD_TYPEA p1[2], 
   const COORD_TYPEA p2[3],
   COORD_TYPEB q0[3], COORD_TYPEB q1[2], COORD_TYPEB q2[3])
  {
    co_factor_3x3_row1(p0, p1, p2, q0);
    co_factor_3x3_row1(p1, p2, p0, q1);
    co_factor_3x3_row1(p2, p0, p1, q2);
  }


  /// Return the transpose of the co-factor matrix of a 3x3 matrix.
  template <typename COORD_TYPEA, typename COORD_TYPEB>
  void transpose_co_factor_3x3
  (const COORD_TYPEA p0[3], const COORD_TYPEA p1[2], 
   const COORD_TYPEA p2[3],
   COORD_TYPEB q0[3], COORD_TYPEB q1[2], COORD_TYPEB q2[3])
  {
    co_factor_3x3(p0, p1, p2, q0, q1, q2);
    transpose_3x3(q0, q1, q2);
  }


  /*!
   *  @brief Compute inverse of a 3x3 matrix.
   *  If absolute values of the determinant is less than or equal to
   *     max_small_determinant, then result is set to the transpose
   *     of the co-factor matrix.
   *  @param[out] determinant Determinant.
   *  @param[out] flag_zero True if determinant is less
   *     than or equal to max_small_determinant.
   *  @pre max_small_determinant >= 0.
   */
  template <typename COORD_TYPEA, typename COORD_TYPEB,
            typename MTYPE, typename DET_TYPE>
  void compute_inverse_3x3
  (const COORD_TYPEA p0[3], const COORD_TYPEA p1[2], 
   const COORD_TYPEA p2[3],
   const MTYPE max_small_determinant, 
   COORD_TYPEB q0[3], COORD_TYPEB q1[2], COORD_TYPEB q2[3],
   DET_TYPE & determinant, bool & flag_zero)
  { 
    const int DIM3(3);
    DET_TYPE D;

    // Initialize flag_zero.
    flag_zero = false;

    transpose_co_factor_3x3(p0, p1, p2, q0, q1, q2);
    determinant_3x3(p0, p1, p2, determinant);
    D = determinant;
    if (D < 0) { D = -D; }

    if (D <= max_small_determinant) {
      flag_zero = true;
    }
    else {
      divide_coord(DIM3, determinant, q0, q0);
      divide_coord(DIM3, determinant, q1, q1);
      divide_coord(DIM3, determinant, q2, q2);
    }
  }


  /*!
   *  @brief Return determinant of three 2D points.
   *  - Return the determinant of the 2x2 matrix [(p1-p0),(p2-p0)].
   */
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void determinant_point_2D
  (const COORD_TYPE p0[2], const COORD_TYPE p1[2], const COORD_TYPE p2[2],
   RESULT_TYPE & result)
  {
    RESULT_TYPE v1[2], v2[2];

    subtract_coord_2D(p1, p0, v1);
    subtract_coord_2D(p2, p0, v2);

    determinant_2x2(v1[0], v1[1], v2[0], v2[1], result);
  }


  /*!
   *  @brief Return determinant of four 3D points.
   *  - Return the determinant of the 2x2 matrix [(p1-p0),(p2-p0)].
   *  - Return the determinant of the 3x3 matrix 
   *    [(p1-p0),(p2-p0),(p3-p0)].
   *  - Equivalent to -determinant((p0,1);(p1,1);(p2,1);(p3,1)).
   */
  template <typename COORD_TYPE, typename RESULT_TYPE>
  inline void determinant_point_3D
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], const COORD_TYPE p3[3],
   RESULT_TYPE & result)
  {
    RESULT_TYPE v1[3], v2[3], v3[3];

    subtract_coord_3D(p1, p0, v1);
    subtract_coord_3D(p2, p0, v2);
    subtract_coord_3D(p3, p0, v3);

    determinant_3x3(v1, v2, v3, result);
  }

  
  /*!
   *  @brief Return true if determinant of four 3D points is
   *    greater than or equal to zero.
   *  @tparam DET_TYPE Determinant type.
   *    - Type used for any internal variables in calculating
   *      the determinant.
   */
  template <typename DET_TYPE, typename COORD_TYPE>
  inline bool is_determinant_point_ge_zero_3D
  (const COORD_TYPE p0[3], const COORD_TYPE p1[3],
   const COORD_TYPE p2[3], const COORD_TYPE p3[3])
  {
    DET_TYPE D;
    determinant_point_3D(p0, p1, p2, p3, D);
    if (D >= 0) { return true; }
    else { return false; }
  }

  //@}


  // **************************************************
  //! @name Jacobian Determinant (3D)
  // **************************************************

  //@{

  /// @brief Compute the Jacobian matrix of a hexahedron at the hexahedron center.
  /// @param cube Cube with facet information.
  /// @pre     cube.Dimension() = 3. 
  /// @param[out] Jacobian 3x3 Jacobian matrix.
  template <typename VTYPE, typename CTYPE0, typename CTYPE1,
            typename CUBE_TYPE>
  void compute_Jacobian_at_hex_center_3D
  (const VTYPE hex_vert[],
   const CTYPE0 * vertex_coord,
   const CUBE_TYPE & cube,
   CTYPE1 Jacobian[3][3])
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE DIM3(3);
    CTYPE1 temp_coord[DIM3];

    for (DTYPE d = 0; d < DIM3; d++) {

      IJK::set_coord_3D(0, Jacobian[d]);

      for (NTYPE j = 0; j < cube.NumFacetVertices(); j++) {
        const NTYPE i0 = cube.FacetVertex(d, j);
        const NTYPE i1 = cube.VertexNeighbor(i0, d);
        const VTYPE iv0 = hex_vert[i0];
        const VTYPE iv1 = hex_vert[i1];
        const CTYPE0 * v0coord = vertex_coord + iv0*DIM3;
        const CTYPE0 * v1coord = vertex_coord + iv1*DIM3;

        IJK::subtract_coord_3D(v1coord, v0coord, temp_coord);
        IJK::add_coord_3D(temp_coord, Jacobian[d], Jacobian[d]);
      }

      // Reduce by factor of 4.
      IJK::divide_coord(DIM3, 4, Jacobian[d], Jacobian[d]);
    }
  }


  /// @brief Compute Jacobian matrix of a hexahedron at a given vertex
  /// @pre cube.Dimension() = 3. 
  /// @param icorner0 Cube corner index.  Possible values are 0,1,...,7.
  /// @param[out] v0coord Pointer to coordinates of vertex at icorner0.
  /// @param[out] w0coord Pointer to coordinates of vertex adjacent
  ///   to icorner0 in direction 0.
  /// @param[out] w1coord Pointer to coordinates of vertex adjacent
  ///   to icorner0 in direction 1.
  /// @param[out] w2coord Pointer to coordinates of vertex adjacent
  ///   to icorner0 in direction 2.
  template <typename VTYPE, typename CTYPE0, typename CTYPE1,
            typename CUBE_TYPE, typename CORNER_TYPE>
  void compute_Jacobian_at_hex_vertex_3D
  (const VTYPE hex_vert[],
   const CTYPE0 * vertex_coord,
   const CUBE_TYPE & cube,
   const CORNER_TYPE icorner0,
   CTYPE1 Jacobian[3][3])
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE DIM3(3);

    const NTYPE jneighbor0 = cube.VertexNeighbor(icorner0, 0);
    const NTYPE jneighbor1 = cube.VertexNeighbor(icorner0, 1);
    const NTYPE jneighbor2 = cube.VertexNeighbor(icorner0, 2);
    const VTYPE iv0 = hex_vert[icorner0];
    const VTYPE jw0 = hex_vert[jneighbor0];
    const VTYPE jw1 = hex_vert[jneighbor1];
    const VTYPE jw2 = hex_vert[jneighbor2];
    const CTYPE0 * v0coord = vertex_coord + iv0*DIM3;
    const CTYPE0 * w0coord = vertex_coord + jw0*DIM3;
    const CTYPE0 * w1coord = vertex_coord + jw1*DIM3;
    const CTYPE0 * w2coord = vertex_coord + jw2*DIM3;

    IJK::subtract_coord_3D(w0coord, v0coord, Jacobian[0]);
    IJK::subtract_coord_3D(w1coord, v0coord, Jacobian[1]);
    IJK::subtract_coord_3D(w2coord, v0coord, Jacobian[2]);
  }


  /// @brief Compute determinant of the Jacobian matrix of a hexahedron
  ///   at the hexahedron center.
  /// @param orientation Orientation of hexahedron. +1 or -1.
  /// @param cube Cube with facet information.
  /// @pre     cube.Dimension() = 3. 
  /// @param[out] Jacobian 3x3 Jacobian matrix.
  template <typename ORIENT_TYPE, typename VTYPE, 
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE>
  void compute_Jacobian_determinant_at_hex_center_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE & Jacobian_determinant,
   CTYPE Jacobian[3][3])
  {
    compute_Jacobian_at_hex_center_3D
      (hex_vert, vertex_coord, cube, Jacobian);

    IJK::determinant_3x3
      (Jacobian[0], Jacobian[1], Jacobian[2], Jacobian_determinant);

    if (orientation < 0) { Jacobian_determinant = -Jacobian_determinant; }
  }


  /// @brief Compute determinant of the Jacobian matrix of a hexahedron
  ///   at the hexahedron center.
  /// - Version which does not return Jacobian matrix.
  /// @param orientation Orientation of hexahedron. +1 or -1.
  /// @param cube Cube with facet information.
  /// @pre     cube.Dimension() = 3. 
  template <typename ORIENT_TYPE, typename VTYPE, 
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE>
  void compute_Jacobian_determinant_at_hex_center_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE & Jacobian_determinant)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;
    const DTYPE DIM3(3);
    CTYPE Jacobian[DIM3][DIM3];

    compute_Jacobian_determinant_at_hex_center_3D
      (hex_vert, orientation, vertex_coord, cube, Jacobian_determinant,
       Jacobian);
  }


  /// @brief Compute Jacobian matrix determinant of a hexahedron at a given corner.
  /// @pre cube.Dimension() = 3. 
  /// @param icorner0 Cube corner index.  Possible values are 0,1,...,7.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE0, typename CTYPE1,
            typename CUBE_TYPE, 
            typename CORNER_TYPE, typename DET_TYPE>
  void compute_Jacobian_determinant_at_hex_vertex_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE0 * vertex_coord,
   const CUBE_TYPE & cube,
   const CORNER_TYPE icorner0,
   DET_TYPE & Jacobian_determinant,
   CTYPE1 Jacobian[3][3])
  {
    // Multiple det at cube vertex i by orient_factor[i] 
    //   to get correct sign of Jacobian determinant.
    const static ORIENT_TYPE orient_factor[] = { 1, -1, -1, 1, -1, 1, 1, -1 };

    compute_Jacobian_at_hex_vertex_3D
      (hex_vert, vertex_coord, cube, icorner0, Jacobian);

    IJK::determinant_3x3
      (Jacobian[0], Jacobian[1], Jacobian[2], Jacobian_determinant);

    if (Jacobian_determinant == 0.0) {
      // Set to +0.0.
      Jacobian_determinant = 0.0;
    }
    else {
      Jacobian_determinant = Jacobian_determinant * orient_factor[icorner0];
      if (orientation < 0) { Jacobian_determinant = -Jacobian_determinant; }
    }
  }


  /// @brief Compute Jacobian matrix determinant of a hexahedron at a given corner.
  /// - Version which does not return Jacobian.
  /// @pre cube.Dimension() = 3. 
  /// @param icorner0 Cube corner index.  Possible values are 0,1,...,7.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename CORNER_TYPE, typename DET_TYPE>
  void compute_Jacobian_determinant_at_hex_vertex_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const CORNER_TYPE icorner0,
   DET_TYPE & Jacobian_determinant)
  {
    CTYPE Jacobian[3][3];

    compute_Jacobian_determinant_at_hex_vertex_3D
      (hex_vert, orientation, vertex_coord,  cube, icorner0,
       Jacobian_determinant, Jacobian);
  }


  /// @brief Compute Jacobian matrix determinant of a hexahedron at a given corner.
  /// - Version with C++ STL vector format for vertex_coord.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename CORNER_TYPE, typename DET_TYPE>
  void compute_Jacobian_determinant_at_hex_vertex_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const std::vector<CTYPE> & vertex_coord,
   const CUBE_TYPE & cube,
   const CORNER_TYPE icorner0,
   DET_TYPE & Jacobian_determinant)
  {
    const CTYPE * vcoord = IJK::vector2pointer(vertex_coord);

    compute_Jacobian_determinant_at_hex_vertex_3D
      (hex_vert, orientation, vcoord, cube, icorner0, Jacobian_determinant);
  }


  /// @brief Compute Jacobian matrix determinant of a hexahedron at a given corner.
  /// - Version with C++ STL vector format for vertex_coord[].
  /// - Version with C++ STL vector hex_vert[] containing 
  ///     an array of vertices of multiple hexahedra.
  template <typename VTYPE, typename NTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename CORNER_TYPE, typename DET_TYPE>
  void compute_Jacobian_determinant_at_hex_vertex_3D
  (const std::vector<VTYPE> & hex_vert,
   const NTYPE ihex,
   const ORIENT_TYPE orientation,
   const std::vector<CTYPE> & vertex_coord,
   const CUBE_TYPE & cube,
   const CORNER_TYPE icorner0,
   DET_TYPE & Jacobian_determinant)
  {
    const NTYPE NUM_VERT_PER_HEX(8);
    const CTYPE * vcoord = IJK::vector2pointer(vertex_coord);
    const VTYPE * hex_i_vert = &(hex_vert[ihex*NUM_VERT_PER_HEX]);

    compute_Jacobian_determinant_at_hex_vertex_3D
      (hex_i_vert, orientation, vcoord, cube, icorner0, Jacobian_determinant);
  }


  /// @brief Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
  /// @pre cube.Dimension() = 3. 
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE0, typename DET_TYPE1>
  void compute_min_max_hexahedron_Jacobian_determinant_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    min_Jacobian_determinant = 0;
    max_Jacobian_determinant = 0;

    compute_Jacobian_determinant_at_hex_center_3D
      (hex_vert, orientation, vertex_coord, cube,
       min_Jacobian_determinant);
    max_Jacobian_determinant = min_Jacobian_determinant;

    for (NTYPE i0 = 0; i0 < cube.NumVertices(); i0++) {
      DET_TYPE0 det;

      compute_Jacobian_determinant_at_hex_vertex_3D
        (hex_vert, orientation, vertex_coord, cube, i0, det);

      if (det < min_Jacobian_determinant) { min_Jacobian_determinant = det; }
      if (det > max_Jacobian_determinant) { max_Jacobian_determinant = det; }
    }
  }


  /// @brief Compute min/max of the eight Jacobian matrix determinants 
  ///   at the eight vertices of a hexahedron.
  /// @pre cube.Dimension() = 3. 
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE0, typename DET_TYPE1>
  void compute_min_max_hex_vert_Jacobian_determinant_3D
  (const VTYPE0 hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant,
   VTYPE1 & vert_with_min_Jacobian_determinant,
   VTYPE2 & vert_with_max_Jacobian_determinant)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    min_Jacobian_determinant = 0;
    max_Jacobian_determinant = 0;

    compute_Jacobian_determinant_at_hex_vertex_3D
      (hex_vert, orientation, vertex_coord, cube, 0, 
       min_Jacobian_determinant);
    max_Jacobian_determinant = min_Jacobian_determinant;
    vert_with_min_Jacobian_determinant = hex_vert[0];
    vert_with_max_Jacobian_determinant = hex_vert[0];

    for (NTYPE i0 = 1; i0 < cube.NumVertices(); i0++) {
      DET_TYPE0 det;
      compute_Jacobian_determinant_at_hex_vertex_3D
        (hex_vert, orientation, vertex_coord, cube, i0, det);

      if (det < min_Jacobian_determinant) { 
        min_Jacobian_determinant = det; 
        vert_with_min_Jacobian_determinant = hex_vert[i0];
      }

      if (det > max_Jacobian_determinant) { 
        max_Jacobian_determinant = det; 
        vert_with_max_Jacobian_determinant = hex_vert[i0];
      }
    }
  }


  /// @brief Compute min/max of the eight Jacobian matrix determinants 
  ///   at the eight vertices of a hexahedron.
  /// - Version which does not return vertices 
  ///     with min/max Jacobian determinants.
  /// @pre cube.Dimension() = 3. 
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE0, typename DET_TYPE1>
  void compute_min_max_hex_vert_Jacobian_determinant_3D
  (const VTYPE0 hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant)
  {
    VTYPE0 vert_with_min_Jacobian_determinant;
    VTYPE0 vert_with_max_Jacobian_determinant;

    compute_min_max_hex_vert_Jacobian_determinant_3D
      (hex_vert, orientation, vertex_coord, cube, 
       min_Jacobian_determinant, max_Jacobian_determinant,
       vert_with_min_Jacobian_determinant, 
       vert_with_max_Jacobian_determinant);
  }


  /// @brief Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
  /// - Version which returns all nine determinants.
  /// @pre cube.Dimension() = 3. 
  /// @param[out] Jacobian_determinant[] 
  ///   - If i < 8, Jacobian_determinant[i] is the Jacobian determinant
  ///     at corner i.
  ///   - Jacobian_determinant[8] is the Jacobian determinant
  ///    at the hexahedron center.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE0, typename DET_TYPE1, typename DET_TYPE2>
  void compute_min_max_hexahedron_Jacobian_determinant_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant,
   DET_TYPE2 Jacobian_determinant[9])
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    min_Jacobian_determinant = 0;
    max_Jacobian_determinant = 0;

    compute_Jacobian_determinant_at_hex_center_3D
      (hex_vert, orientation, vertex_coord, cube,
       Jacobian_determinant[8]);
    min_Jacobian_determinant = Jacobian_determinant[8];
    max_Jacobian_determinant = min_Jacobian_determinant;

    for (NTYPE i0 = 0; i0 < cube.NumVertices(); i0++) {
      DET_TYPE0 det;
      compute_Jacobian_determinant_at_hex_vertex_3D
        (hex_vert, orientation, vertex_coord, cube, i0, det);

      Jacobian_determinant[i0] = det;
      if (det < min_Jacobian_determinant) { min_Jacobian_determinant = det; }
      if (det > max_Jacobian_determinant) { max_Jacobian_determinant = det; }
    }

  }

  /// @brief Compute the eight Jacobian matrix determinants of the eight 
  ///   hexahedron vertices.
  /// @pre cube.Dimension() = 3. 
  /// @param[out] Jacobian_determinant[] 
  ///   - Jacobian_determinant[i] is the Jacobian determinant at corner i.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE>
  void compute_Jacobian_determinant_at_all_hex_vert_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE Jacobian_determinant[8])
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    for (NTYPE i0 = 0; i0 < cube.NumVertices(); i0++) {
      compute_Jacobian_determinant_at_hex_vertex_3D
        (hex_vert, orientation, vertex_coord, cube, i0, 
         Jacobian_determinant[i0]);
    }
  }


  /// @brief Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
  /// - Version with C++ STL vector format for vertex_coord.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE0, typename DET_TYPE1>
  void compute_min_max_hexahedron_Jacobian_determinant_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const std::vector<CTYPE> & vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant)
  {
    const CTYPE * vcoord = IJK::vector2pointer(vertex_coord);

    compute_min_max_hexahedron_Jacobian_determinant_3D
      (hex_vert, orientation, vcoord, cube,
       min_Jacobian_determinant, max_Jacobian_determinant);
  }


  /// @brief Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
  /// - Version which returns all nine determinants.
  /// - Version with C++ STL vector format for vertex_coord.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE0, typename DET_TYPE1, typename DET_TYPE2>
  void compute_min_max_hexahedron_Jacobian_determinant_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const std::vector<CTYPE> & vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant,
   DET_TYPE2 Jacobian_determinant[9])
  {
    const CTYPE * vcoord = IJK::vector2pointer(vertex_coord);

    compute_min_max_hexahedron_Jacobian_determinant_3D
      (hex_vert, orientation, vcoord, cube,
       min_Jacobian_determinant, max_Jacobian_determinant,
       Jacobian_determinant);
  }


  /// @brief Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
  /// - Version with C++ STL vector format for vertex_coord[].
  /// - Version with C++ STL vector hex_vert[] containing 
  ///     an array of vertices of multiple hexahedra.
  template <typename VTYPE, typename NTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE0, typename DET_TYPE1>
  void compute_min_max_hexahedron_Jacobian_determinant_3D
  (const std::vector<VTYPE> & hex_vert,
   const NTYPE ihex,
   const ORIENT_TYPE orientation,
   const std::vector<CTYPE> & vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant)
  {
    const NTYPE NUM_VERT_PER_HEX(8);
    const VTYPE * hex_i_vert = &(hex_vert[ihex*NUM_VERT_PER_HEX]);

    compute_min_max_hexahedron_Jacobian_determinant_3D
      (hex_i_vert, orientation, vertex_coord, cube,
       min_Jacobian_determinant, max_Jacobian_determinant);
  }


  /// @brief Compute min/max of the nine Jacobian matrix determinants of a hexahedron.
  /// - Version which returns all nine determinants.
  /// - Version with C++ STL vector format for vertex_coord.
  /// - Version with C++ STL vector hex_vert[] containing 
  ///     an array of vertices of multiple hexahedra.
  template <typename VTYPE, typename NTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename DET_TYPE0, typename DET_TYPE1, typename DET_TYPE2>
  void compute_min_max_hexahedron_Jacobian_determinant_3D
  (const std::vector<VTYPE> & hex_vert,
   const NTYPE ihex,
   const ORIENT_TYPE orientation,
   const std::vector<CTYPE> & vertex_coord,
   const CUBE_TYPE & cube,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant,
   DET_TYPE2 Jacobian_determinant[9])
  {
    const NTYPE NUM_VERT_PER_HEX(8);
    const VTYPE * hex_i_vert = &(hex_vert[ihex*NUM_VERT_PER_HEX]);

    compute_min_max_hexahedron_Jacobian_determinant_3D
      (hex_i_vert, orientation, vertex_coord, cube,
       min_Jacobian_determinant, max_Jacobian_determinant,
       Jacobian_determinant);
  }


  /// @brief Compute normalized determinant of the Jacobian matrix of a hexahedron
  ///   at the hexahedron center.
  /// @param orientation Orientation of hexahedron. +1 or -1.
  /// @param cube Cube with facet information.
  /// @pre     cube.Dimension() = 3. 
  template <typename ORIENT_TYPE, typename VTYPE, 
            typename CTYPE, typename CUBE_TYPE, typename MTYPE,
            typename DET_TYPE>
  void compute_normalized_Jacobian_determinant_at_hex_center_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   DET_TYPE & Jacobian_determinant,
   bool & flag_zero)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE DIM3(3);

    CTYPE Jacobian[DIM3][DIM3];
    DET_TYPE L0, L1, L2;

    compute_Jacobian_determinant_at_hex_center_3D
      (hex_vert, orientation, vertex_coord, cube, Jacobian_determinant,
       Jacobian);

    compute_magnitude_3D(Jacobian[0], L0);
    compute_magnitude_3D(Jacobian[1], L1);
    compute_magnitude_3D(Jacobian[2], L2);

    if (L0 > max_small_magnitude &&
        L1 > max_small_magnitude &&
        L2 > max_small_magnitude) {
      flag_zero = false;
      Jacobian_determinant = (Jacobian_determinant/(L0*L1*L2));
      return;
    }
    else {
      Jacobian_determinant = 0;
      flag_zero = true;
    }
  
  }


  /// @brief Compute normalized Jacobian matrix determinant 
  ///   of a hexahedron at a given corner.\br
  /// - Compute the determinant of the thee unit vectors in the directions
  ///   of the three edges incident on each hexahedron corner.
  /// @pre cube.Dimension() = 3. 
  /// @param icorner0 Cube corner index.  Possible values are 0,1,...,7.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, 
            typename CORNER_TYPE, typename MTYPE, typename DET_TYPE>
  void compute_normalized_Jacobian_determinant_at_hex_vertex_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const CORNER_TYPE icorner0,
   const MTYPE max_small_magnitude,
   DET_TYPE & Jacobian_determinant,
   bool & flag_zero)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE DIM3(3);
    DET_TYPE L0, L1, L2;
    CTYPE Jacobian[3][3];

    compute_Jacobian_determinant_at_hex_vertex_3D
      (hex_vert, orientation, vertex_coord,  cube, icorner0,
       Jacobian_determinant,  Jacobian);

    compute_magnitude_3D(Jacobian[0], L0);
    compute_magnitude_3D(Jacobian[1], L1);
    compute_magnitude_3D(Jacobian[2], L2);

    if (L0 > max_small_magnitude &&
        L1 > max_small_magnitude &&
        L2 > max_small_magnitude) {
      flag_zero = false;
      Jacobian_determinant = (Jacobian_determinant/(L0*L1*L2));
      return;
    }
    else {
      Jacobian_determinant = 0;
      flag_zero = true;
    }

  }


  /// @brief Compute min/max of the nine normalized Jacobian matrix determinants 
  ///   of a hexahedron.
  /// @pre cube.Dimension() = 3. 
  /// @param[out] num_determinants 
  ///   Number of determinants computed (not skipped).
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, typename MTYPE,
            typename DET_TYPE0, typename DET_TYPE1,
            typename NTYPE>
  void compute_min_max_hexahedron_normalized_Jacobian_determinant_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant,
   NTYPE & num_determinants)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE CUBE_NTYPE;

    bool flag_zero;
    DET_TYPE0 det;

    num_determinants = 0;
    min_Jacobian_determinant = 0;
    max_Jacobian_determinant = 0;

    compute_normalized_Jacobian_determinant_at_hex_center_3D
      (hex_vert, orientation, vertex_coord, cube, max_small_magnitude,
       det, flag_zero);

    if (!flag_zero) { 
      min_Jacobian_determinant = det;
      max_Jacobian_determinant = det;
      num_determinants++; 
    }

    for (CUBE_NTYPE i0 = 0; i0 < cube.NumVertices(); i0++) {

      compute_normalized_Jacobian_determinant_at_hex_vertex_3D
        (hex_vert, orientation, vertex_coord, cube, i0, 
         max_small_magnitude, det, flag_zero);

      if (!flag_zero) {

        if (num_determinants == 0) {
          min_Jacobian_determinant = det;
          max_Jacobian_determinant = det;
        }
        else {
          if (det < min_Jacobian_determinant) 
            { min_Jacobian_determinant = det; }
          if (det > max_Jacobian_determinant) 
            { max_Jacobian_determinant = det; }
        }

        num_determinants++;
      }
    }

  }

  /// @brief Compute min/max of the eight normalized Jacobian matrix determinants 
  ///   at the eight vertices of a hexahedron.
  /// @pre cube.Dimension() = 3. 
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, typename MTYPE,
            typename DET_TYPE0, typename DET_TYPE1,
            typename NTYPE>
  void compute_min_max_hex_vert_normalized_Jacobian_determinant_3D
  (const VTYPE0 hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   DET_TYPE0 & min_Jacobian_determinant,
   DET_TYPE1 & max_Jacobian_determinant,
   VTYPE1 & vert_with_min_Jacobian_determinant,
   VTYPE2 & vert_with_max_Jacobian_determinant,
   NTYPE & num_determinants)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE CUBE_NTYPE;

    DET_TYPE0 det;
    bool flag_zero;

    min_Jacobian_determinant = 0;
    max_Jacobian_determinant = 0;
    num_determinants = 0;

    compute_normalized_Jacobian_determinant_at_hex_vertex_3D
      (hex_vert, orientation, vertex_coord, cube, 0, max_small_magnitude,
       det, flag_zero);


    if (!flag_zero) {
      min_Jacobian_determinant = det;
      max_Jacobian_determinant = det;
      vert_with_min_Jacobian_determinant = hex_vert[0];
      vert_with_max_Jacobian_determinant = hex_vert[0];
      num_determinants++;
    }

    for (CUBE_NTYPE i0 = 1; i0 < cube.NumVertices(); i0++) {
      compute_normalized_Jacobian_determinant_at_hex_vertex_3D
        (hex_vert, orientation, vertex_coord, cube, i0, max_small_magnitude,
         det, flag_zero);

      if (!flag_zero) {
        
        if (num_determinants == 0) {
          min_Jacobian_determinant = det; 
          vert_with_min_Jacobian_determinant = hex_vert[i0];
          max_Jacobian_determinant = det; 
          vert_with_max_Jacobian_determinant = hex_vert[i0];
        }
        else {
          if (det < min_Jacobian_determinant) { 
            min_Jacobian_determinant = det; 
            vert_with_min_Jacobian_determinant = hex_vert[i0];
          }

          if (det > max_Jacobian_determinant) { 
            max_Jacobian_determinant = det; 
            vert_with_max_Jacobian_determinant = hex_vert[i0];
          }
        }
        
        num_determinants++;
      }
    }
  }

  /// @brief Compute the eight normalized Jacobian matrix determinants of the eight 
  ///   hexahedron vertices.
  /// @pre cube.Dimension() = 3. 
  /// @param[out] Jacobian_determinant[].
  ///   - Jacobian_determinant[i] is the normalized Jacobian determinant 
  ///       at corner i.
  /// @param[out] flag_zero[].
  ///   - flag_zero[i] is true if corner i is incident on a zero length edge.
  template <typename VTYPE, typename ORIENT_TYPE,
            typename CTYPE, typename CUBE_TYPE, typename MTYPE,
            typename DET_TYPE>
  void compute_normalized_Jacobian_determinant_at_all_hex_vert_3D
  (const VTYPE hex_vert[],
   const ORIENT_TYPE orientation,
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   DET_TYPE Jacobian_determinant[8],
   bool flag_zero[8])
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    for (NTYPE i0 = 0; i0 < cube.NumVertices(); i0++) {
      compute_normalized_Jacobian_determinant_at_hex_vertex_3D
        (hex_vert, orientation, vertex_coord, cube, i0, max_small_magnitude,
         Jacobian_determinant[i0], flag_zero[i0]);
    }
  }


  //@}


  // **************************************************
  //! @name Jacobian Determinant Shape (3D)
  // **************************************************

  //@{

  /// @brief Compute shape metric from 3x3 Jacobian matrix.
  /// @param orientation Orientation of hexahedron. +1 or -1.
  /// @param Jacobian 3x3 Jacobian matrix.
  template <typename CTYPE, typename MTYPE, 
            typename SHAPE_TYPE>
  void compute_shape_from_3x3_Jacobian
  (const CTYPE Jacobian[3][3],
   const MTYPE max_small_magnitude,
   SHAPE_TYPE & shape_value,
   bool & flag_zero)
  {
    const CTYPE ONE_THIRD(1.0/3.0);
    SHAPE_TYPE x0, x1, x2;

    IJK::determinant_3x3
      (Jacobian[0], Jacobian[1], Jacobian[2], shape_value);

    compute_sum_of_squares_3D(Jacobian[0], x0);
    compute_sum_of_squares_3D(Jacobian[1], x1);
    compute_sum_of_squares_3D(Jacobian[2], x2);

    SHAPE_TYPE y = x0+x1+x2;

    if (y > max_small_magnitude){
      flag_zero = false;
      shape_value = shape_value*shape_value;
      shape_value = std::pow(shape_value, ONE_THIRD);
      shape_value = 3.0*(shape_value/y);
      return;
    }
    else {
      shape_value = 0;
      flag_zero = true;
    }

  }


  /// @brief Compute hexahedron shape metric based on determinant 
  ///   of the Jacobian matrix of a hexahedron at the hexahedron center.
  /// @param cube Cube with facet information.
  /// @pre     cube.Dimension() = 3. 
  /// @param[out] Jacobian 3x3 Jacobian matrix.
  template <typename VTYPE, typename CTYPE, typename CUBE_TYPE, 
            typename MTYPE, typename SHAPE_TYPE>
  void compute_Jacobian_shape_at_hex_center_3D
  (const VTYPE hex_vert[],
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   SHAPE_TYPE & shape_value,
   bool & flag_zero,
   CTYPE Jacobian[3][3])
  {
    compute_Jacobian_at_hex_center_3D
      (hex_vert, vertex_coord, cube, Jacobian);

    compute_shape_from_3x3_Jacobian
      (Jacobian, max_small_magnitude, shape_value, flag_zero);
  }


  /*!
   *  @overload
   *  @brief Compute hexahedron shape metric based on determinant 
   *    of the Jacobian matrix of a hexahedron at the hexahedron center.
   *  - Version which does not return Jacobian matrix.
   *  @param cube Cube with facet information.
   *  @pre     cube.Dimension() = 3. 
   */
  template <typename VTYPE, typename CTYPE, typename CUBE_TYPE, 
            typename MTYPE, typename SHAPE_TYPE>
  void compute_Jacobian_shape_at_hex_center_3D
  (const VTYPE hex_vert[],
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   SHAPE_TYPE & shape_value,
   bool & flag_zero)
  {
    CTYPE Jacobian[3][3];

    compute_Jacobian_shape_at_hex_center_3D
      (hex_vert, vertex_coord, cube, max_small_magnitude, 
       shape_value, flag_zero, Jacobian);
  }


  /// @brief Compute hexahedron shape metric based on determinant 
  ///   of the Jacobian matrix of a hexahedron at a given corner.
  /// @pre cube.Dimension() = 3. 
  /// @param icorner0 Cube corner index.  Possible values are 0,1,...,7.
  template <typename VTYPE, typename CTYPE0, typename CTYPE1,
            typename CUBE_TYPE, typename CORNER_TYPE, 
            typename MTYPE, typename SHAPE_TYPE>
  void compute_Jacobian_shape_at_hex_vertex_3D
  (const VTYPE hex_vert[],
   const CTYPE0 * vertex_coord,
   const CUBE_TYPE & cube,
   const CORNER_TYPE icorner0,
   const MTYPE max_small_magnitude,
   SHAPE_TYPE & shape_value,
   bool & flag_zero,
   CTYPE1 Jacobian[3][3])
  {
    compute_Jacobian_at_hex_vertex_3D
      (hex_vert, vertex_coord, cube, icorner0, Jacobian);

    compute_shape_from_3x3_Jacobian
      (Jacobian, max_small_magnitude, shape_value, flag_zero);
  }


  /// @brief Compute hexahedron shape metric based on determinant 
  ///   of the Jacobian matrix of a hexahedron at a given corner.
  /// - Version which does not return Jacobian matrix.
  /// @pre cube.Dimension() = 3. 
  /// @param icorner0 Cube corner index.  Possible values are 0,1,...,7.
  template <typename VTYPE, typename CTYPE, 
            typename CUBE_TYPE, typename CORNER_TYPE, 
            typename MTYPE, typename SHAPE_TYPE>
  void compute_Jacobian_shape_at_hex_vertex_3D
  (const VTYPE hex_vert[],
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const CORNER_TYPE icorner0,
   const MTYPE max_small_magnitude,
   SHAPE_TYPE & shape_value,
   bool & flag_zero)
  {
    CTYPE Jacobian[3][3];

    compute_Jacobian_shape_at_hex_vertex_3D
      (hex_vert, vertex_coord, cube, icorner0, max_small_magnitude,
       shape_value, flag_zero, Jacobian);
  }


  /// @brief Compute min/max of the nine shape metrics based on the Jacobian
  ///   determinant of a hexahedron.
  /// @pre cube.Dimension() = 3. 
  template <typename VTYPE, typename CTYPE, 
            typename CUBE_TYPE, typename MTYPE,
            typename SHAPE_TYPE0, typename SHAPE_TYPE1,
            typename NTYPE>
  void compute_min_max_hexahedron_Jacobian_shape_3D
  (const VTYPE hex_vert[],
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   SHAPE_TYPE0 & min_Jacobian_shape,
   SHAPE_TYPE1 & max_Jacobian_shape,
   NTYPE & num_shape_metrics)
  {
    bool flag_zero;

    num_shape_metrics = 0;
    min_Jacobian_shape = 0;
    max_Jacobian_shape = 0;

    compute_Jacobian_shape_at_hex_center_3D
      (hex_vert, vertex_coord, cube, max_small_magnitude,
       min_Jacobian_shape, flag_zero);
    if (!flag_zero) {
      max_Jacobian_shape = min_Jacobian_shape;
      num_shape_metrics++;
    }

    for (NTYPE i0 = 0; i0 < cube.NumVertices(); i0++) {
      SHAPE_TYPE0 shape_value;

      compute_Jacobian_shape_at_hex_vertex_3D
        (hex_vert, vertex_coord, cube, i0, max_small_magnitude,
         shape_value, flag_zero);

      if (!flag_zero) {
        if ((num_shape_metrics == 0) || (shape_value < min_Jacobian_shape)) 
          { min_Jacobian_shape = shape_value; }
        if ((num_shape_metrics == 0) || (shape_value > max_Jacobian_shape))
          { max_Jacobian_shape = shape_value; }
        num_shape_metrics++;
      }
    }

  }


  /// @brief Compute min/max of the eight Jacobian shape metrics
  ///   at the eight vertices of a hexahedron.
  /// @pre cube.Dimension() = 3. 
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename CTYPE, typename CUBE_TYPE, typename MTYPE,
            typename SHAPE_TYPE0, typename SHAPE_TYPE1,
            typename NTYPE>
  void compute_min_max_hex_vert_Jacobian_shape_3D
  (const VTYPE0 hex_vert[],
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   SHAPE_TYPE0 & min_Jacobian_shape,
   SHAPE_TYPE1 & max_Jacobian_shape,
   VTYPE1 & vert_with_min_Jacobian_shape,
   VTYPE2 & vert_with_max_Jacobian_shape,
   NTYPE & num_shapes)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE CUBE_NTYPE;

    SHAPE_TYPE0 shape_metric;
    bool flag_zero;

    min_Jacobian_shape = 0;
    max_Jacobian_shape = 0;
    num_shapes = 0;

    compute_Jacobian_shape_at_hex_vertex_3D
      (hex_vert, vertex_coord, cube, 0, max_small_magnitude,
       shape_metric, flag_zero);


    if (!flag_zero) {
      min_Jacobian_shape = shape_metric;
      max_Jacobian_shape = shape_metric;
      vert_with_min_Jacobian_shape = hex_vert[0];
      vert_with_max_Jacobian_shape = hex_vert[0];
      num_shapes++;
    }

    for (CUBE_NTYPE i0 = 1; i0 < cube.NumVertices(); i0++) {
      compute_Jacobian_shape_at_hex_vertex_3D
        (hex_vert, vertex_coord, cube, i0, max_small_magnitude,
         shape_metric, flag_zero);

      if (!flag_zero) {
        
        if (num_shapes == 0) {
          min_Jacobian_shape = shape_metric; 
          vert_with_min_Jacobian_shape = hex_vert[i0];
          max_Jacobian_shape = shape_metric; 
          vert_with_max_Jacobian_shape = hex_vert[i0];
        }
        else {
          if (shape_metric < min_Jacobian_shape) { 
            min_Jacobian_shape = shape_metric; 
            vert_with_min_Jacobian_shape = hex_vert[i0];
          }

          if (shape_metric > max_Jacobian_shape) { 
            max_Jacobian_shape = shape_metric; 
            vert_with_max_Jacobian_shape = hex_vert[i0];
          }
        }
        
        num_shapes++;
      }
    }
  }

  /// @brief Compute the eight Jacobian shape metrics of the eight hex vertices.
  /// @pre cube.Dimension() = 3. 
  /// @param[out] shape_metric[].
  ///   - shape_metric[i] is the shape metric base on the Jacobian matrix
  ///       at corner i.
  /// @param[out] flag_zero[].
  ///   - flag_zero[i] is true if corner i is incident on a zero length edge.
  template <typename VTYPE, typename CTYPE, 
            typename CUBE_TYPE, typename MTYPE,
            typename SHAPE_TYPE>
  void compute_Jacobian_shape_at_all_hex_vert_3D
  (const VTYPE hex_vert[],
   const CTYPE * vertex_coord,
   const CUBE_TYPE & cube,
   const MTYPE max_small_magnitude,
   SHAPE_TYPE shape_metric[8],
   bool flag_zero[8])
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    for (NTYPE i0 = 0; i0 < cube.NumVertices(); i0++) {
      compute_Jacobian_shape_at_hex_vertex_3D
        (hex_vert, vertex_coord, cube, i0, max_small_magnitude,
         shape_metric[i0], flag_zero[i0]);
    }
  }


  //@}


  // **************************************************
  //! @name 3D intersections
  // **************************************************

  //@{

  /*!
   *  @brief Compute the intersection of a line and a plane.
   *  - Line contains \a coord0[] and has direction \a line_dir[].
   *  - Plane contains \a coord1[] and has orthogonal direction \a orth_dir[].
   *  @param coord0[] Coordinates of point contained by line. 
   *  @param line_dir[] Line direction.
   *  @param coord1[] Coordinates of point contained by hyperplane.
   *  @param orth_dir[] Direction orthgonal to \a line_dir[].
   *  @param min_abs_cos_angle_LO Minimum value of the absolute value
   *      of the cosine of the angle between \a line_dir[] and \a orth_dir[].
   *  @param[out] intersection_point[] 
   *                Intersection point of line and hyperplane.
   *  @param[out] flag_succeeded True if computation succeeded.
   *     Computations succeeds if:
   *       abs(inner_product(line_dir[],orth_dir[])) > min_abs_cos_angle_LO.
   *  @pre \a line_dir[] is a unit vector.
   *  @pre \a orth_dir[] is a unit vector.
   *  @pre \a min_abs_cos_angle_LO >= 0.
   */
  template <typename CTYPE0, typename CTYPE1, 
            typename VTYPE0, typename VTYPE1, typename MTYPE,
            typename CTYPE2>
  void intersect_line_plane_3D
  (const CTYPE0 coord0[], const VTYPE0 line_dir[],
   const CTYPE1 coord1[], const VTYPE1 orth_dir[],
   const MTYPE min_abs_cos_angle_LO,
   CTYPE2 intersection_point[], bool & flag_succeeded)
  {
    const int DIM3(3);
    intersect_line_hyperplane
      (DIM3, coord0, line_dir, coord1, orth_dir, min_abs_cos_angle_LO,
       intersection_point, flag_succeeded);
  }

  //@}


  // **************************************************
  //! @name 3D cross product
  // **************************************************

  //@{

  /// Compute the cross product of two 3D vectors.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void compute_cross_product_3D
  (const VTYPE0 v0[], const VTYPE1 v1[], VTYPE2 v2[])
  {
    determinant_2x2(v0[1], v0[2], v1[1], v1[2], v2[0]);
    determinant_2x2(v0[0], v0[2], v1[0], v1[2], v2[1]);
    v2[1] = -v2[1];
    determinant_2x2(v0[0], v0[1], v1[0], v1[1], v2[2]);
  }

  //@}

}

#endif
