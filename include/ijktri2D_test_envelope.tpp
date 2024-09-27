/*!
 *  \file ijktri2D_test_envelope.tpp
 *  @brief Test routines for testing envelopes formed
 *    by four quad edges and dual grid edge.
 *  - Version 0.4.0
 *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
 */



/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2022-2023 Rephael Wenger

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


#ifndef _IJKTRI2D_TEST_ENVELOPE_

#include "ijktri2D_envelope.tpp"
#include "ijkcheck_geom.tpp"

#include <cmath>
#include <functional>
#include <string>


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJKTRI2D_TEST_ENVELOPE {

  // **************************************************
  // ERROR MESSAGE ROUTINES
  // **************************************************

  /*!
   *  @brief Add quad vertex and edge endpoint coordinates 
   *    to error message.
   *  @param test_descriptionA First part of test_description.
   *    - May be NULL.
   *  @param test_descriptionB Second part of test_description.
   *    - May be NULL.
   *    - Ignored if test_descriptinA is NULL.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void add_quad_coord_to_error_message
  (const std::vector<std::string> & test_description,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   IJK::ERROR & error)
  {
    const int DIM3(3);

    for (int i = 0; i < test_description.size(); i++) {
      error.AddMessage("  ", test_description[i]);
    }

    error.AddArrayMessage("    w0coord: ", w0coord, DIM3, "");
    error.AddArrayMessage("    w1coord: ", w1coord, DIM3, "");
    error.AddArrayMessage("    w2coord: ", w2coord, DIM3, "");
    error.AddArrayMessage("    w3coord: ", w3coord, DIM3, "");
  }


  /*!
   *  @brief Add quad vertex and edge endpoint coordinates 
   *    to error message.
   *  @param test_descriptionA First part of test_description.
   *    - May be NULL.
   *  @param test_descriptionB Second part of test_description.
   *    - May be NULL.
   *    - Ignored if test_descriptinA is NULL.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1>
  void add_quad_and_endpoint_coord_to_error_message
  (const std::vector<std::string> & test_description,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   IJK::ERROR & error)
  {
    const int DIM3(3);

    add_quad_coord_to_error_message
      (test_description, w0coord, w1coord, w2coord, w3coord, error);
    error.AddArrayMessage("    edge_endpoint0_coord: ", v0coord, DIM3, "");
    error.AddArrayMessage("    edge_endpoint1_coord: ", v1coord, DIM3, "");
  }


  /*!
   *  @brief DEPRECATED: Add quad vertex and edge endpoint coordinates 
   *    to error message.
   *  @param test_descriptionA First part of test_description.
   *    - May be NULL.
   *  @param test_descriptionB Second part of test_description.
   *    - May be NULL.
   *    - Ignored if test_descriptinA is NULL.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void add_quad_coord_to_error_message
  (const char * test_descriptionA,
   const char * test_descriptionB,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   IJK::ERROR & error)
  {
    const int DIM3(3);
    
    if (test_descriptionA != NULL) {
      if (test_descriptionB != NULL) {
        error.AddMessage
          ("  ", test_descriptionA, test_descriptionB);
      }
      else {
        error.AddMessage
          ("  ", test_descriptionA);
      }
    }

    error.AddArrayMessage("    w0coord: ", w0coord, DIM3, "");
    error.AddArrayMessage("    w1coord: ", w1coord, DIM3, "");
    error.AddArrayMessage("    w2coord: ", w2coord, DIM3, "");
    error.AddArrayMessage("    w3coord: ", w3coord, DIM3, "");
  }

  
  /*!
   *  @brief Add quad vertex and edge endpoint coordinates 
   *    to error message.
   *  @param test_descriptionA First part of test_description.
   *    - May be NULL.
   *  @param test_descriptionB Second part of test_description.
   *    - May be NULL.
   *    - Ignored if test_descriptinA is NULL.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1>
  void add_quad_and_endpoint_coord_to_error_message
  (const char * test_descriptionA,
   const char * test_descriptionB,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   IJK::ERROR & error)
  {
    const int DIM3(3);

    add_quad_coord_to_error_message
      (test_descriptionA, test_descriptionB,
       w0coord, w1coord, w2coord, w3coord, error);
    error.AddArrayMessage("    edge_endpoint0_coord: ", v0coord, DIM3, "");
    error.AddArrayMessage("    edge_endpoint1_coord: ", v1coord, DIM3, "");
  }


  /*!
   *  @overload
   *  @brief Add quad vertex and edge endpoint coordinates 
   *    to error message.
   *  - Version with only test description in a single char string.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1>
  void add_quad_and_endpoint_coord_to_error_message
  (const char * test_description,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   IJK::ERROR & error)
  {
    add_quad_and_endpoint_coord_to_error_message
      (test_description, NULL,
       w0coord, w1coord, w2coord, w3coord,
       v0coord, v1coord,
       error);
  }


  /*!
   *  @brief Add quad orienation and quadrant containing w0
   *    to error message.
   */
  template <typename ITYPE>
  void add_quad_orientation_and_quadrant_to_error_message
  (const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   IJK::ERROR & error)
  {
    if (flag_quad_pos_orient)
      { error.AddMessage("    Quad orientation: +1."); }
    else
      { error.AddMessage("    Quad orientation: -1."); }

    error.AddMessage
      ("    Quadrant containing w0: ", iquadrant_w0, "");
  }

  
  // **************************************************
  // CHECK ENVELOPE
  // **************************************************

  /*!
   *  @brief Test are quad diagonals in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool check_are_quad_diagonals_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const std::vector<std::string> & test_description,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    bool is_diag02_in_envelope, is_diag13_in_envelope;

    if (!check_edge_direction
        (DIM3, v0coord, v1coord, edge_direction, error)) {    
      add_quad_and_endpoint_coord_to_error_message
        (test_description, w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord, error);
      return false;
    }

    if (!check_dual_quadrilateral_orientation_and_quadrants_3D
        (w0coord, w1coord, w2coord, w3coord,
         v0coord, edge_direction,
         flag_quad_pos_orient, iquadrant_w0, error)) {
      add_quad_and_endpoint_coord_to_error_message
        (test_description, w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord, error);
      return false;
    }

    const bool flag_both_diag_in_envelope =
      envelope_function.are_both_quad_diagonals_in_quad_edge_envelope
      (w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0, 
       is_diag02_in_envelope, is_diag13_in_envelope);

    if (flag_diag02_in_envelope != is_diag02_in_envelope) {
      error.AddMessage("Incorrect value for is_diag02_in_envelope.");
      error.AddMessage
        ("  is_diag02_in_envelope: ", int(is_diag02_in_envelope), "");
      error.AddMessage
        ("  Value should be: ", int(flag_diag02_in_envelope), "");
      add_quad_and_endpoint_coord_to_error_message
        (test_description, w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord, error);
      add_quad_orientation_and_quadrant_to_error_message
        (flag_quad_pos_orient, iquadrant_w0, error);
      return false;
    }

    if (flag_diag13_in_envelope != is_diag13_in_envelope) {
      error.AddMessage("Incorrect value for is_diag13_in_envelope.");
      error.AddMessage
        ("  is_diag13_in_envelope: ", int(is_diag13_in_envelope), "");
      error.AddMessage
        ("  Value should be: ", int(flag_diag13_in_envelope), "");
      add_quad_and_endpoint_coord_to_error_message
        (test_description, w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord, error);
      add_quad_orientation_and_quadrant_to_error_message
        (flag_quad_pos_orient, iquadrant_w0, error);
      return false;
    }
    
    if (flag_both_diag_in_envelope) {
      if (!is_diag02_in_envelope) {
        error.AddMessage("Inconsistency in returned flags.");
        error.AddMessage
          ("  Returned value indicates both diagonals are in envelope but");
        error.AddMessage("  is_diag02_in_envelope is false.");
        add_quad_and_endpoint_coord_to_error_message
          (test_description, w0coord, w1coord, w2coord, w3coord,
           v0coord, v1coord, error);
        add_quad_orientation_and_quadrant_to_error_message
          (flag_quad_pos_orient, iquadrant_w0, error);
        return false;
      }

      if (!is_diag13_in_envelope) {
        error.AddMessage("Inconsistency in returned flags.");
        error.AddMessage
          ("  Returned value indicates both diagonals are in envelope but");
        error.AddMessage("  is_diag13_in_envelope is false.");
        add_quad_and_endpoint_coord_to_error_message
          (test_description, w0coord, w1coord, w2coord, w3coord,
           v0coord, v1coord, error);
        add_quad_orientation_and_quadrant_to_error_message
          (flag_quad_pos_orient, iquadrant_w0, error);
        return false;
      }
    }
    else {
      if (is_diag02_in_envelope && is_diag13_in_envelope) {
        error.AddMessage("Inconsistency in returned flags.");
        error.AddMessage
          ("  Returned value indicates both diagonals are not in envelope but");
        error.AddMessage
          ("  Returned value indicates both diagonals are in envelope but");
        error.AddMessage
          ("  both is_diag02_in_envelope and is_diag13_in_envelope are true.");
        add_quad_and_endpoint_coord_to_error_message
          (test_description, w0coord, w1coord, w2coord, w3coord,
           v0coord, v1coord, error);
        add_quad_orientation_and_quadrant_to_error_message
          (flag_quad_pos_orient, iquadrant_w0, error);
        return false;
      }
    }

    return true;
  }

  
  /*!
   *  @overload
   *  @brief Test are quad diagonals in quad edge envelope.
   *  - Version with single string for test description.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool check_are_quad_diagonals_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const std::string & test_descriptionI,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   IJK::ERROR & error)
  {
    std::vector<std::string> test_description;

    test_description.push_back(test_descriptionI);
    return check_are_quad_diagonals_in_envelope
      (envelope_function, test_description,
       w0coord, w1coord, w2coord, w3coord,
       v0coord, v1coord, edge_direction,
       flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);
  }

  
  /*!
   *  @brief DEPRECATED: Test are quad diagonals in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool check_are_quad_diagonals_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_descriptionA,
   const char * test_descriptionB,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    bool is_diag02_in_envelope, is_diag13_in_envelope;

    if (!check_edge_direction(DIM3, v0coord, v1coord,
                              edge_direction, error)) {    
      add_quad_and_endpoint_coord_to_error_message
        (test_descriptionA, test_descriptionB,
         w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord, error);
      return false;
    }

    if (!check_dual_quadrilateral_orientation_and_quadrants_3D
        (w0coord, w1coord, w2coord, w3coord,
         v0coord, edge_direction,
         flag_quad_pos_orient, iquadrant_w0, error)) {
      add_quad_and_endpoint_coord_to_error_message
        (test_descriptionA, test_descriptionB,
         w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord, error);
      return false;
    }

    const bool flag_both_diag_in_envelope =
      envelope_function.are_both_quad_diagonals_in_quad_edge_envelope
      (w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0, 
       is_diag02_in_envelope, is_diag13_in_envelope);

    if (flag_diag02_in_envelope != is_diag02_in_envelope) {
      error.AddMessage("Incorrect value for is_diag02_in_envelope.");
      error.AddMessage
        ("  is_diag02_in_envelope: ", int(is_diag02_in_envelope), "");
      error.AddMessage
        ("  Value should be: ", int(flag_diag02_in_envelope), "");
      add_quad_and_endpoint_coord_to_error_message
        (test_descriptionA, test_descriptionB,
         w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord, error);
      add_quad_orientation_and_quadrant_to_error_message
        (flag_quad_pos_orient, iquadrant_w0, error);
      return false;
    }

    if (flag_diag13_in_envelope != is_diag13_in_envelope) {
      error.AddMessage("Incorrect value for is_diag13_in_envelope.");
      error.AddMessage
        ("  is_diag13_in_envelope: ", int(is_diag13_in_envelope), "");
      error.AddMessage
        ("  Value should be: ", int(flag_diag13_in_envelope), "");
      add_quad_and_endpoint_coord_to_error_message
        (test_descriptionA, test_descriptionB,
         w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord, error);
      add_quad_orientation_and_quadrant_to_error_message
        (flag_quad_pos_orient, iquadrant_w0, error);
      return false;
    }
    
    if (flag_both_diag_in_envelope) {
      if (!is_diag02_in_envelope) {
        error.AddMessage("Inconsistency in returned flags.");
        error.AddMessage
          ("  Returned value indicates both diagonals are in envelope but");
        error.AddMessage("  is_diag02_in_envelope is false.");
        add_quad_and_endpoint_coord_to_error_message
          (test_descriptionA, test_descriptionB,
           w0coord, w1coord, w2coord, w3coord,
           v0coord, v1coord, error);
        add_quad_orientation_and_quadrant_to_error_message
          (flag_quad_pos_orient, iquadrant_w0, error);
        return false;
      }

      if (!is_diag13_in_envelope) {
        error.AddMessage("Inconsistency in returned flags.");
        error.AddMessage
          ("  Returned value indicates both diagonals are in envelope but");
        error.AddMessage("  is_diag13_in_envelope is false.");
        add_quad_and_endpoint_coord_to_error_message
          (test_descriptionA, test_descriptionB,
           w0coord, w1coord, w2coord, w3coord,
           v0coord, v1coord, error);
        add_quad_orientation_and_quadrant_to_error_message
          (flag_quad_pos_orient, iquadrant_w0, error);
        return false;
      }
    }
    else {
      if (is_diag02_in_envelope && is_diag13_in_envelope) {
        error.AddMessage("Inconsistency in returned flags.");
        error.AddMessage
          ("  Returned value indicates both diagonals are not in envelope but");
        error.AddMessage
          ("  Returned value indicates both diagonals are in envelope but");
        error.AddMessage
          ("  both is_diag02_in_envelope and is_diag13_in_envelope are true.");
        add_quad_and_endpoint_coord_to_error_message
          (test_descriptionA, test_descriptionB,
           w0coord, w1coord, w2coord, w3coord,
           v0coord, v1coord, error);
        add_quad_orientation_and_quadrant_to_error_message
          (flag_quad_pos_orient, iquadrant_w0, error);
        return false;
      }
    }

    return true;
  }


  /*!
   *  @overload
   *  @brief DEPRECATED: Test are quad diagonals in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  - Version with std::string for test_description_suffix.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool check_are_quad_diagonals_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const std::string & test_description_suffix,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   IJK::ERROR & error)
  {
    if (test_description_suffix.empty()) {
      return check_are_quad_diagonals_in_envelope
        (envelope_function, test_description_prefix, NULL,
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, flag_diag13_in_envelope,
         error);
    }
    else {
      return check_are_quad_diagonals_in_envelope
        (envelope_function, test_description_prefix,
         test_description_suffix.c_str(),
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, flag_diag13_in_envelope,
         error);
    }
  }

  
  /*!
   *  @brief Multiple tests for are quad diagonals in quad edge envelope. 
   *    Rotate vertex arguments.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool check_are_quad_diagonals_in_envelope_rotate_isov
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const std::string & test_description_prefix,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const bool flag_clockwise,
   IJK::ERROR & error)
  {
    std::string test_description_strA;
    std::string test_description_suffix;
    std::string test_description;

    if (flag_clockwise)
      { test_description_strA = "CW"; }
    else
      { test_description_strA = "CCW"; }

    test_description_suffix = " (" + test_description_strA + "0123):";
    test_description = test_description_prefix + test_description_suffix;
    if (!check_are_quad_diagonals_in_envelope
        (envelope_function, test_description,
         w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, flag_diag13_in_envelope,
         error))
      { return false; }

    const ITYPE iquadrant_w1 =
      IJK::next_quadrant(flag_quad_pos_orient, iquadrant_w0);
    test_description_suffix = " (" + test_description_strA + "1230):";
    test_description = test_description_prefix + test_description_suffix;
    if (!check_are_quad_diagonals_in_envelope
        (envelope_function, test_description,
         w1coord, w2coord, w3coord, w0coord,
         v0coord, v1coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w1,
         flag_diag13_in_envelope, flag_diag02_in_envelope, 
         error))
      { return false; }

    const ITYPE iquadrant_w2 =
      IJK::next_quadrant(flag_quad_pos_orient, iquadrant_w1);
    test_description_suffix = " (" + test_description_strA + "2301):";
    test_description = test_description_prefix + test_description_suffix;
    if (!check_are_quad_diagonals_in_envelope
        (envelope_function, test_description,
         w2coord, w3coord, w0coord, w1coord,
         v0coord, v1coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w2,
         flag_diag02_in_envelope, flag_diag13_in_envelope,
         error))
      { return false; }

    const ITYPE iquadrant_w3 =
      IJK::next_quadrant(flag_quad_pos_orient, iquadrant_w2);
    test_description_suffix = " (" + test_description_strA + "3012):";
    test_description = test_description_prefix + test_description_suffix;
    if (!check_are_quad_diagonals_in_envelope
        (envelope_function, test_description,
         w3coord, w0coord, w1coord, w2coord,
         v0coord, v1coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w3,
         flag_diag13_in_envelope, flag_diag02_in_envelope, 
         error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Multiple tests for are quad diagonals in quad edge envelope. 
   *    Rotate and reverse vertex arguments.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool check_are_quad_diagonals_in_envelope_reorder_isov
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   IJK::ERROR & error)
  {
    if (!check_are_quad_diagonals_in_envelope_rotate_isov
        (envelope_function,
         test_description_prefix,
         w0coord, w1coord, w2coord, w3coord,
         v0coord, v1coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, flag_diag13_in_envelope,
         false, error))
      { return false; }

    const ITYPE iquadrant_w3 =
      IJK::prev_quadrant(flag_quad_pos_orient, iquadrant_w0);
    if (!check_are_quad_diagonals_in_envelope_rotate_isov
        (envelope_function,
         test_description_prefix,
         w3coord, w2coord, w1coord, w0coord,
         v0coord, v1coord,
         edge_direction, !flag_quad_pos_orient, iquadrant_w3,
         flag_diag13_in_envelope, flag_diag02_in_envelope, 
         true, error))
      { return false; }

    return true;
  }


  namespace {

    /// Replace (c[0], c[1], c[2]) with (c[1], c[2], c[0]).
    template <typename CTYPE>
    void rotate_coord_3D(CTYPE coord[3])
    {
      std::swap(coord[0], coord[1]);
      std::swap(coord[2], coord[1]);
    }
  }

  
  /*!
   *  @brief Multiple tests for are quad diagonals in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool multi_check_are_quad_diagonals_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    
    CTYPEW0 w0coordII[DIM3];
    CTYPEW1 w1coordII[DIM3];
    CTYPEW2 w2coordII[DIM3];
    CTYPEW3 w3coordII[DIM3];
    CTYPEV0 v0coordII[DIM3];
    CTYPEV1 v1coordII[DIM3];
         
    if (!check_are_quad_diagonals_in_envelope_reorder_isov
        (envelope_function, test_description_prefix,
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, flag_diag13_in_envelope, error))
      { return false; }

    if (edge_direction != 0) {
      // No grid rotation.
      return true;
    }

    std::copy(w0coord, w0coord+DIM3, w0coordII);
    std::copy(w1coord, w1coord+DIM3, w1coordII);
    std::copy(w2coord, w2coord+DIM3, w2coordII);
    std::copy(w3coord, w3coord+DIM3, w3coordII);
    std::copy(v0coord, v0coord+DIM3, v0coordII);
    std::copy(v1coord, v1coord+DIM3, v1coordII);

    rotate_coord_3D(w0coordII);
    rotate_coord_3D(w1coordII);
    rotate_coord_3D(w2coordII);
    rotate_coord_3D(w3coordII);
    rotate_coord_3D(v0coordII);
    rotate_coord_3D(v1coordII);

    // Edge direction is 1.
    if (!check_are_quad_diagonals_in_envelope_reorder_isov
        (envelope_function, test_description_prefix,
         w0coordII, w1coordII, w2coordII, w3coordII, v0coordII, v1coordII,
         1, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, flag_diag13_in_envelope, error))
      { return false; }

    // Rotate again
    rotate_coord_3D(w0coordII);
    rotate_coord_3D(w1coordII);
    rotate_coord_3D(w2coordII);
    rotate_coord_3D(w3coordII);
    rotate_coord_3D(v0coordII);
    rotate_coord_3D(v1coordII);

    // Edge direction is 2.
    if (!check_are_quad_diagonals_in_envelope_reorder_isov
        (envelope_function, test_description_prefix,
         w0coordII, w1coordII, w2coordII, w3coordII, v0coordII, v1coordII,
         2, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, flag_diag13_in_envelope, error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Multiple tests for are quad diagonals in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool multi_check_are_quad_diagonals_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   IJK::ERROR & error)
  {
    return multi_check_are_quad_diagonals_in_envelope
      (envelope_function, test_description_prefix,
       wcoord[0], wcoord[1], wcoord[2], wcoord[3], v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);
  }


  // **************************************************
  // SET flag_diag_in_envelope
  // **************************************************

  /*!
   *  @brief Set flag_is_diag_in_envelope to false if offset > 0.
   */ 
  template <typename CTYPE>
  void set_flag_diag_in_envelope
  (const CTYPE offset, bool & flag_is_diag_in_envelope)
  {
    if (offset > 0)
      { flag_is_diag_in_envelope = false; }
    else
      { flag_is_diag_in_envelope = true; }
  }
  
  
  // **************************************************
  // SET COORD
  // **************************************************

  /// Set coordinates of wcoord[].
  template <typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPEW>
  void set_coord_3D
  (const CTYPE0 c0, const CTYPE1 c1, const CTYPE2 c2,
   CTYPEW wcoord[3])
  {
    wcoord[0] = c0;
    wcoord[1] = c1;
    wcoord[2] = c2;
  }


  /// Add (c0,c1,c2) to coordinates of wcoord[].
  template <typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPEW>
  void add_to_coord_3D
  (const CTYPE0 c0, const CTYPE1 c1, const CTYPE2 c2,
   CTYPEW wcoord[3])
  {
    wcoord[0] += c0;
    wcoord[1] += c1;
    wcoord[2] += c2;
  }


  /// @brief Set quad vertices to form symmetric square
  ///    with z-coordinate cz.
  template <typename CTYPEZ,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_symmetric_square
  (const CTYPEZ cz,
   CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_coord_3D(0.5, 0.5, cz, w0coord);
    set_coord_3D(1.5, 0.5, cz, w1coord);
    set_coord_3D(1.5, 1.5, cz, w2coord);
    set_coord_3D(0.5, 1.5, cz, w3coord);
  }


  /// @overload
  /// @brief Set quad vertices to form symmetric square
  ///    with z-coordinate cz.
  template <typename CTYPEZ, typename CTYPEW>
  void set_quad_coord_symmetric_square
  (const CTYPEZ cz, CTYPEW * wcoord[4])
  {
    set_quad_coord_symmetric_square
      (cz, wcoord[0], wcoord[1], wcoord[2], wcoord[3]);
  }


  /*!
   *  @brief Set quad coordinates to rectangle.
   *  @param x0 Coordinates of w0coord[0] and w3coord[0].
   *  @param x1 Coordinates of w1coord[0] and w2coord[0].
   *  @param y0 Coordinates of w0coord[1] and w1coord[1].
   *  @param y1 Coordinates of w2coord[1] and w2coord[1].
   *  @param cz Coordinates of w*coord[2].
   */
  template <typename CTYPEX, typename CTYPEY, typename CTYPEZ,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_rectangle
  (const CTYPEX x0, const CTYPEX x1,
   const CTYPEY y0, const CTYPEX y1,
   const CTYPEZ cz,
   CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    // Set x-coordinates.
    w0coord[0] = w3coord[0] = x0;
    w1coord[0] = w2coord[0] = x1;

    // Set y-coordinates.
    w0coord[1] = w1coord[1] = y0;
    w2coord[1] = w3coord[1] = y1;

    // Set z-coordinates.
    w0coord[2] = w1coord[2] = w2coord[2] = w3coord[2] = cz;
  }

  /*!
   *  @brief Set quad coordinates to rectangle.
   *  @overload
   *  @param x0 Coordinates of w0coord[0] and w3coord[0].
   *  @param x1 Coordinates of w1coord[0] and w2coord[0].
   *  @param y0 Coordinates of w0coord[1] and w1coord[1].
   *  @param y1 Coordinates of w2coord[1] and w2coord[1].
   *  @param cz Coordinates of w*coord[2].
   */
  template <typename CTYPEX, typename CTYPEY, typename CTYPEZ,
            typename CTYPEW>
  void set_quad_coord_rectangle
  (const CTYPEX x0, const CTYPEX x1,
   const CTYPEY y0, const CTYPEX y1,
   const CTYPEZ cz,
   CTYPEW * wcoord[4])
  {
    set_quad_coord_rectangle
      (x0, x1, y0, y1, cz,
       wcoord[0], wcoord[1], wcoord[2], wcoord[3]);
  }

    
  /*!
   *  @brief Set quad coordinates to vertical rectangle.
   *  @param x0 Coordinates of w0coord[0] and w3coord[0].
   *  @param x1 Coordinates of w1coord[0] and w2coord[0].
   *  @param cz Coordinates of w*coord[2].
   */
  template <typename CTYPEX0, typename CTYPEX1, typename CTYPEZ,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_vertical_rectangle
  (const CTYPEX0 x0, const CTYPEX1 x1,
   const CTYPEZ cz,
   CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_quad_coord_symmetric_square
      (cz, w0coord, w1coord, w2coord, w3coord);

    // Set x-coordinates.
    w0coord[0] = w3coord[0] = x0;
    w1coord[0] = w2coord[0] = x1;
  }
  
    
  /*!
   *  @brief Set quad coordinates to horizontal rectangle.
   *  @param y0 Coordinates of w0coord[1] and w3coord[1].
   *  @param y1 Coordinates of w1coord[1] and w2coord[1].
   *  @param cz Coordinates of w*coord[2].
   */
  template <typename CTYPEX0, typename CTYPEX1,
            typename CTYPEZ,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_horizontal_rectangle
  (const CTYPEX0 y0, const CTYPEX1 y1,
   const CTYPEZ cz,
   CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_quad_coord_symmetric_square
      (cz, w0coord, w1coord, w2coord, w3coord);

    // Set y-coordinates.
    w0coord[1] = w1coord[1] = y0;
    w2coord[1] = w3coord[1] = y1;
  }

  
  template <typename CTYPEZ,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_left_vertical_rectangle
  (const CTYPEZ cz,
   CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_quad_coord_vertical_rectangle
      (0.5, 1.1, cz, w0coord, w1coord, w2coord, w3coord);
  }

    
  template <typename CTYPEZ,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_right_vertical_rectangle
  (const CTYPEZ cz, CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_quad_coord_vertical_rectangle
      (0.9, 1.5, cz, w0coord, w1coord, w2coord, w3coord);
  }

    
  template <typename CTYPEZ,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_low_horizontal_rectangle
  (const CTYPEZ cz,
   CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_quad_coord_horizontal_rectangle
      (0.5, 1.1, cz, w0coord, w1coord, w2coord, w3coord);
  }


  template <typename CTYPEZ,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_high_horizontal_rectangle
  (const CTYPEZ cz,
   CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_quad_coord_horizontal_rectangle
      (0.9, 1.5, cz, w0coord, w1coord, w2coord, w3coord);
  }

    
  /*!
   *  @brief Set quadrilateral coordinates to form
   *    horizontal trapezoid (parallel to y = 1 plane.)
   *  @param half_length0 Half length of trapezoid bottom.
   *  @param half_length1 Half length of trapezoid top.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename LTYPE>
  void set_quad_coord_horizontal_trapezoid
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3],
   const LTYPE half_length0, const LTYPE half_length1)
  {
    set_coord_3D(1-half_length0, 0.5, 0.5, w0coord);
    set_coord_3D(1+half_length0, 0.5, 0.5, w1coord);
    set_coord_3D(1+half_length1, 1.5, 0.5, w2coord);
    set_coord_3D(1-half_length1, 1.5, 0.5, w3coord);
  }


  /*!
   *  @brief Set quadrilateral coordinates to form
   *    vertical trapezoid (parallel to x = 1 plane.)
   *  @param half_length0 Half length of trapezoid left.
   *  @param half_length1 Half length of trapezoid right.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename LTYPE>
  void set_quad_coord_vertical_trapezoid
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3],
   const LTYPE half_length0, const LTYPE half_length1)
  {
    set_coord_3D(0.5, 1-half_length0, 0.5, w0coord);
    set_coord_3D(1.5, 1-half_length1, 0.5, w1coord);
    set_coord_3D(1.5, 1+half_length1, 0.5, w2coord);
    set_coord_3D(0.5, 1+half_length0, 0.5, w3coord);
  }

    
  /*!
   *  @brief Skew lower left corner down.
   *  - Also skew upper right corner down
   *    and lower right and upper left corners up.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void skew_LL_corner_down
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    w0coord[2] = 0.1;
    w1coord[2] = 0.9;
    w2coord[2] = 0.1;
    w3coord[2] = 0.9;
  }

    
  /*!
   *  @brief Skew lower left corner up.
   *  - Also skew upper right corner up
   *    and lower right and upper left corners down.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void skew_LL_corner_up
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    w0coord[2] = 0.9;
    w1coord[2] = 0.1;
    w2coord[2] = 0.9;
    w3coord[2] = 0.1;
  }

    
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_concave_v0
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_coord_3D(0.9, 0.9, 0.5, w0coord);
    set_coord_3D(1.1, 0.1, 0.5, w1coord);
    set_coord_3D(1.5, 1.5, 0.5, w2coord);
    set_coord_3D(0.1, 1.1, 0.5, w3coord);
  }

    
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_concave_v1
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_coord_3D(0.9, 0.1, 0.5, w0coord);
    set_coord_3D(1.1, 0.9, 0.5, w1coord);
    set_coord_3D(1.9, 1.1, 0.5, w2coord);
    set_coord_3D(0.5, 1.5, 0.5, w3coord);
  }

    
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_concave_v2
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_coord_3D(0.5, 0.5, 0.5, w0coord);
    set_coord_3D(1.9, 0.9, 0.5, w1coord);
    set_coord_3D(1.1, 1.1, 0.5, w2coord);
    set_coord_3D(0.9, 1.9, 0.5, w3coord);
  }

    
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_concave_v3
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_coord_3D(0.1, 0.9, 0.5, w0coord);
    set_coord_3D(1.5, 0.5, 0.5, w1coord);
    set_coord_3D(1.1, 1.9, 0.5, w2coord);
    set_coord_3D(0.9, 1.1, 0.5, w3coord);
  }

  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3>
  void set_quad_coord_on_axes
  (CTYPEW0 w0coord[3], CTYPEW1 w1coord[3],
   CTYPEW2 w2coord[3], CTYPEW3 w3coord[3])
  {
    set_coord_3D(0.5, 1, 0.5, w0coord);
    set_coord_3D(1, 0.5, 0.5, w1coord);
    set_coord_3D(1.5, 1, 0.5, w2coord);
    set_coord_3D(1, 1.5, 0.5, w3coord);
  }

  
  // **************************************************
  // SET COORD AND CHECK ENVELOPE
  // **************************************************

  /*!
   *  @brief Set rectangle and test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param zcoord[] Array of z-coordinates. Test each z-coordinate.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEX, typename CTYPEY, typename CTYPEZ>
  bool set_rectangle_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const CTYPEX x0,
   const CTYPEX x1,
   const CTYPEY y0,
   const CTYPEY y1,
   const std::vector<CTYPEZ> zcoord,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    
    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    const int edge_direction = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;

    set_quad_coord_rectangle(x0, x1, y0, y1, 0, wcoord);
    
    for (int i = 0; i < zcoord.size(); i++) {

      w0coord[2] = zcoord[i];
      w1coord[2] = zcoord[i];
      w2coord[2] = zcoord[i];
      w3coord[2] = zcoord[i];
      
      if (!multi_check_are_quad_diagonals_in_envelope
          (envelope_function,
           test_description_prefix,
           wcoord, v0coord, v1coord,
           edge_direction, flag_quad_pos_orient, iquadrant_w0,
           flag_diag02_in_envelope, flag_diag13_in_envelope,
           error))
        { return false; }
    }

    return true;
  }


  /*!
   *  @brief Set left vertical rectangle and 
   *    test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param zcoord[] Array of z-coordinates. Test each z-coordinate.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE, typename CTYPEZ>
  bool set_left_vertical_rectangle_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const std::vector<CTYPEZ> zcoord,
   IJK::ERROR & error)
  {
    return set_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
      (envelope_function, test_description_prefix,
       flag_diag02_in_envelope,
       flag_diag13_in_envelope,
       0.5, 1.1, 0.5, 1.5, zcoord, error);
  }


  /*!
   *  @brief Set right vertical rectangle and 
   *    test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param zcoord[] Array of z-coordinates. Test each z-coordinate.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE, typename CTYPEZ>
  bool set_right_vertical_rectangle_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const std::vector<CTYPEZ> zcoord,
   IJK::ERROR & error)
  {
    return set_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
      (envelope_function, test_description_prefix,
       flag_diag02_in_envelope,
       flag_diag13_in_envelope,
       0.9, 1.5, 0.5, 1.5, zcoord, error);
  }


  /*!
   *  @brief Set low horizontal rectangle and 
   *    test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param zcoord[] Array of z-coordinates. Test each z-coordinate.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE, typename CTYPEZ>
  bool set_low_horizontal_rectangle_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const std::vector<CTYPEZ> zcoord,
   IJK::ERROR & error)
  {
    return set_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
      (envelope_function, test_description_prefix,
       flag_diag02_in_envelope,
       flag_diag13_in_envelope,
       0.5, 1.5, 0.5, 1.1, zcoord, error);
  }


  /*!
   *  @brief Set high horizontal rectangle and 
   *    test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE, typename CTYPEZ>
  bool set_high_horizontal_rectangle_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const std::vector<CTYPEZ> zcoord,
   IJK::ERROR & error)
  {
    return set_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
      (envelope_function, test_description_prefix,
       flag_diag02_in_envelope,
       flag_diag13_in_envelope,
       0.5, 1.5, 0.9, 1.5, zcoord, error);
  }


  /*!
   *  @brief Set symmetric square and
   *    test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param zcoord[] Array of z-coordinates. Test each z-coordinate.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE, typename CTYPEZ>
  bool set_symmetric_square_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const std::vector<CTYPEZ> zcoord,
   IJK::ERROR & error)
  {
    return set_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
      (envelope_function, test_description_prefix,
       flag_diag02_in_envelope,
       flag_diag13_in_envelope,
       0.5, 1.5, 0.5, 1.5, zcoord, error);
  }


  /*!
   *  @brief Symmetrically offset quad vertices around edge 
   *    and test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param xoffset02,yoffset02,zoffset02 Offsets of w0 and w2.
   *    - wcoord[0][0] += xoffset02.
   *    - wcoord[0][1] += yoffset02.
   *    - wcoord[0][2] += zoffset02.
   *    - wcoord[2][0] += -xoffset02.
   *    - wcoord[2][1] += -yoffset02.
   *    - wcoord[2][2] += -zoffset02.
   *  @param xoffset13,yoffset13,zoffset13 Offsets of w1 and w3.
   *    - wcoord[1][0] += -xoffset13.
   *    - wcoord[1][1] += yoffset13.
   *    - wcoord[1][2] += zoffset13.
   *    - wcoord[3][0] += xoffset13.
   *    - wcoord[3][1] += -yoffset13.
   *    - wcoord[3][2] += -zoffset13.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE,
            typename CTYPEX, typename CTYPEY, typename CTYPEZ>
  bool offset_vertices_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const CTYPEX xoffset02,
   const CTYPEY yoffset02,
   const CTYPEZ zoffset02,
   const CTYPEX xoffset13,
   const CTYPEY yoffset13,
   const CTYPEZ zoffset13,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    CTYPEW new_wcoord_array[NUM_VERT_PER_QUADRILATERAL*DIM3];
    CTYPEW * new_wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { new_wcoord_array, new_wcoord_array+DIM3,
        new_wcoord_array+2*DIM3, new_wcoord_array+3*DIM3 };

    for (int iw = 0; iw < NUM_VERT_PER_QUADRILATERAL; iw++)
      { std::copy(wcoord[iw], wcoord[iw]+DIM3, new_wcoord[iw]); }

    add_to_coord_3D(xoffset02, yoffset02, zoffset02, new_wcoord[0]);
    add_to_coord_3D(-xoffset13, yoffset13, zoffset13, new_wcoord[1]);
    add_to_coord_3D(-xoffset02, -yoffset02, -zoffset02, new_wcoord[2]);
    add_to_coord_3D(xoffset13, -yoffset13, -zoffset13, new_wcoord[3]);
    
    return multi_check_are_quad_diagonals_in_envelope
      (envelope_function,
       test_description_prefix,
       new_wcoord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);
  }


  /*!
   *  @overload
   *  @brief Symmetrically offset quad vertices around edge 
   *    and test are quads in quad edge envelope.
   *  - Version where (xoffset13 == yoffset02) and (yoffset13 == xoffset02)
   *    and (zoffset13 == zoffset02 == 0).
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW,
            typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE,
            typename CTYPEX, typename CTYPEY>
  bool offset_vertices_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const CTYPEX xoffset,
   const CTYPEY yoffset,
   IJK::ERROR & error)
  {
    return
      offset_vertices_multi_check_in_envelope
      (envelope_function, test_description_prefix, wcoord,
       v0coord, v1coord, edge_direction,
       flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       xoffset, yoffset, 0, yoffset, xoffset, 0, error);
  }


  /*!
   *  @brief Set rhombus and test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param xoffset,yoffset Offsets of vertices.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEX, typename CTYPEY>
  bool set_rhombus_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const CTYPEX xoffset,
   const CTYPEY yoffset,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    
    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    const int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;

    set_quad_coord_symmetric_square(0.5, wcoord);

    return
      offset_vertices_multi_check_in_envelope
      (envelope_function, test_description_prefix,
       wcoord, v0coord, v1coord, edge_dir,
       flag_quad_pos_orient, iquadrant_w0, 
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       xoffset, yoffset, error);
  }


  /*!
   *  @brief Set quad whose projection in (x,y) is a parallelogram 
   *    and test are quads in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param xoffset02,yoffset02 Offsets of w0 and w2.
   *  @param xoffset13,yoffset13 Offsets of w1 and w3.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEX, typename CTYPEY, typename CTYPEZ>
  bool set_parallelogram_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const CTYPEX xoffset02,
   const CTYPEY yoffset02,
   const CTYPEZ zoffset02,
   const CTYPEX xoffset13,
   const CTYPEY yoffset13,
   const CTYPEZ zoffset13,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    
    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    const int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;

    set_quad_coord_symmetric_square(0.5, wcoord);

    return
      offset_vertices_multi_check_in_envelope
      (envelope_function, test_description_prefix,
       wcoord, v0coord, v1coord, edge_dir,
       flag_quad_pos_orient, iquadrant_w0, 
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       xoffset02, yoffset02, zoffset02,
       xoffset13, yoffset13, zoffset13, error);
  }

  
  /*!
   *  @brief Set iw on cube boundary and test are quad diagonals 
   *    in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param iw Index of vertex to move to cube boundary.
   *  @param ic Set wcoord[iw][ic] to 1.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW, typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool set_on_boundary_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const int iw,
   const int ic,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    
    CTYPEW wcoord_i[DIM3];
    CTYPEW * new_wcoord[NUM_VERT_PER_QUADRILATERAL];

    std::copy(wcoord[iw], wcoord[iw]+DIM3, wcoord_i);
    std::copy(wcoord, wcoord+NUM_VERT_PER_QUADRILATERAL,
              new_wcoord);
    
    // Set iw to be on cube boundary.
    wcoord_i[ic] = 1;
    new_wcoord[iw] = wcoord_i;

    return multi_check_are_quad_diagonals_in_envelope
      (envelope_function,
       test_description_prefix,
       new_wcoord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);

    return true;
  }


  /*!
   *  @brief Set two vertices on cube boundary and test are quad diagonals 
   *    in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param iw0 Index of first vertex to move to cube boundary.
   *  @param iw1 Index of second vertex to move to cube boundary.
   *  @param ic Set wcoord[iw0][ic] and wcoord[iw1][ic] to 1.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW, typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool setII_on_boundary_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const int iw0,
   const int iw1,
   const int ic,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    
    CTYPEW wcoord_i0[DIM3];
    CTYPEW wcoord_i1[DIM3];
    CTYPEW * new_wcoord[NUM_VERT_PER_QUADRILATERAL];

    std::copy(wcoord[iw0], wcoord[iw0]+DIM3, wcoord_i0);
    std::copy(wcoord[iw1], wcoord[iw1]+DIM3, wcoord_i1);
    std::copy(wcoord, wcoord+NUM_VERT_PER_QUADRILATERAL,
              new_wcoord);

    // Set iw to be on cube boundary.
    wcoord_i0[ic] = wcoord_i1[ic] = 1;
    new_wcoord[iw0] = wcoord_i0;
    new_wcoord[iw1] = wcoord_i1;

    return multi_check_are_quad_diagonals_in_envelope
      (envelope_function,
       test_description_prefix,
       new_wcoord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);

    return true;
  }


  /*!
   *  @brief Set three vertices on cube boundary and test are quad diagonals 
   *    in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param iw0 Index of first vertex to move to cube boundary.
   *  @param iw1 Index of second vertex to move to cube boundary.
   *  @param iw1 Index of third vertex to move to cube boundary.
   *  @param ic Set wcoord[iw0][ic] and wcoord[iw1][ic]
   *    and wcoord[iw2][ic] to 1.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW, typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool setIII_on_boundary_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const int iw0,
   const int iw1,
   const int iw2,
   const int ic,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    
    CTYPEW wcoord_i0[DIM3];
    CTYPEW wcoord_i1[DIM3];
    CTYPEW wcoord_i2[DIM3];
    CTYPEW * new_wcoord[NUM_VERT_PER_QUADRILATERAL];

    std::copy(wcoord[iw0], wcoord[iw0]+DIM3, wcoord_i0);
    std::copy(wcoord[iw1], wcoord[iw1]+DIM3, wcoord_i1);
    std::copy(wcoord[iw2], wcoord[iw2]+DIM3, wcoord_i2);
    std::copy(wcoord, wcoord+NUM_VERT_PER_QUADRILATERAL,
              new_wcoord);

    // Set iw to be on cube boundary.
    wcoord_i0[ic] = wcoord_i1[ic] = wcoord_i2[ic] = 1;
    new_wcoord[iw0] = wcoord_i0;
    new_wcoord[iw1] = wcoord_i1;
    new_wcoord[iw2] = wcoord_i2;

    return multi_check_are_quad_diagonals_in_envelope
      (envelope_function,
       test_description_prefix,
       new_wcoord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);

    return true;
  }



  /*!
   *  @brief Set x and z coordinates of two vertices 
   *    and test are quad diagonals in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param iw0 Index of first vertex to move to grid edge.
   *  @param iw1 Index of second vertex to move to grid edge.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEV0, typename CTYPEV1,
            typename CTYPEW, typename DIR_TYPE, typename ITYPE,
            typename CTYPEX, typename CTYPEZ>
  bool setIIxz_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const int iw0,
   const int iw1,
   const CTYPEX xcoord,
   const CTYPEZ zcoord,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    
    CTYPEW wcoord_i0[DIM3];
    CTYPEW wcoord_i1[DIM3];
    CTYPEW * new_wcoord[NUM_VERT_PER_QUADRILATERAL];

    std::copy(wcoord[iw0], wcoord[iw0]+DIM3, wcoord_i0);
    std::copy(wcoord[iw1], wcoord[iw1]+DIM3, wcoord_i1);
    std::copy(wcoord, wcoord+NUM_VERT_PER_QUADRILATERAL,
              new_wcoord);

    // Set iw to be on cube boundary.
    wcoord_i0[0] = wcoord_i1[0] = xcoord;
    wcoord_i0[2] = wcoord_i1[2] = zcoord;
    new_wcoord[iw0] = wcoord_i0;
    new_wcoord[iw1] = wcoord_i1;

    return multi_check_are_quad_diagonals_in_envelope
      (envelope_function,
       test_description_prefix,
       new_wcoord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);

    return true;
  }

  
  /*!
   *  @brief Set wcoord[iz][2] to zcoord.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param iw Index of vertex to move to axis.
   *  @param ic Set wcoord[iw][ic] to 1.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW, typename CTYPEV0, typename CTYPEV1,
            typename CTYPEZ, typename DIR_TYPE, typename ITYPE>
  bool set_z_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const int iw,
   const CTYPEZ zcoord,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    
    CTYPEW wcoord_i[DIM3];
    CTYPEW * new_wcoord[NUM_VERT_PER_QUADRILATERAL];

    std::copy(wcoord[iw], wcoord[iw]+DIM3, wcoord_i);
    std::copy(wcoord, wcoord+NUM_VERT_PER_QUADRILATERAL,
              new_wcoord);

    // Set iw to be on cube boundary.
    wcoord_i[2] = zcoord;
    new_wcoord[iw] = wcoord_i;

    return multi_check_are_quad_diagonals_in_envelope
      (envelope_function,
       test_description_prefix,
       new_wcoord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);

    return true;
  }


  /*!
   *  @brief Set iw on cube boundary and test are quad diagonals 
   *    in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param _are_both_qdiag_in_envelope() Function to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param iw Index of vertex to move to axis.
   */
  template <typename ENVELOPE_FUNCTION_TYPE,
            typename CTYPEW, typename CTYPEV0, typename CTYPEV1,
            typename DIR_TYPE, typename ITYPE>
  bool set_on_grid_edge_multi_check_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const char * test_description_prefix,
   CTYPEW * const wcoord[4],
   const CTYPEV0 v0coord[3],
   const CTYPEV1 v1coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const bool flag_diag02_in_envelope,
   const bool flag_diag13_in_envelope,
   const int iw,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;
    
    CTYPEW wcoord_i[DIM3];
    CTYPEW * new_wcoord[NUM_VERT_PER_QUADRILATERAL];

    std::copy(wcoord[iw], wcoord[iw]+DIM3, wcoord_i);
    std::copy(wcoord, wcoord+NUM_VERT_PER_QUADRILATERAL,
              new_wcoord);

    // Set iw to be on grid edge
    wcoord_i[d1] = 1;
    wcoord_i[d2] = 1;
    new_wcoord[iw] = wcoord_i;

    return multi_check_are_quad_diagonals_in_envelope
      (envelope_function,
       test_description_prefix,
       new_wcoord, v0coord, v1coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       flag_diag02_in_envelope, flag_diag13_in_envelope,
       error);

    return true;
  }


  // **************************************************
  // TEST ENVELOPE
  // **************************************************

  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Rectangle cases.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_rectangle
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    std::vector<CTYPEW> zcoord;

    zcoord.push_back(0.5);
    zcoord.push_back(0.4);
    zcoord.push_back(0.6);

      
    if (!set_left_vertical_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "Left vertical rectangle.",
         true, true, zcoord, error))
      { return false; }

    if (!set_right_vertical_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "Right vertical rectangle.",
         true, true, zcoord, error))
      { return false; }

    if (!set_low_horizontal_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "Low horizontal rectangle.",
         true, true, zcoord, error))
      { return false; }

    if (!set_high_horizontal_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "High horizontal rectangle.",
         true, true, zcoord, error))
      { return false; }

    if (!set_symmetric_square_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Symmetric square.",
         true, true, zcoord, error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Parallelogram cases.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param flag_only_projection_test If true, only use envelope
   *    projection test.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_parallelogram
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const bool flag_only_projection_test,
   IJK::ERROR & error)
  {
    bool flag_diag02_in_envelope, flag_diag13_in_envelope;
    
    if (!set_rhombus_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Rhombus (offset(0.1,0))",
         true, true, 0.1, 0.0, error))
      { return false; }

    if (!set_rhombus_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Rhombus (offset(0,0.1))", 
         true, true, 0.0, 0.1, error))
      { return false; }

    if (!set_rhombus_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Rhombus (on axes)",
         true, true, 0.5, 0.0, error))
      { return false; }

    if (!set_rhombus_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Rhombus (on axes)",
         true, true, 0.0, 0.5, error))
      { return false; }

    // Note: Projection tests on the following sometimes return false,
    //   even though the diagonal is in the envelope.
    flag_diag02_in_envelope = !flag_only_projection_test;
    flag_diag13_in_envelope = !flag_only_projection_test;
    
    if (!set_parallelogram_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Parallelogram (w0,w2) offset (0.1, 0.1, 0.0)",
         true, flag_diag13_in_envelope,
         0.1, 0.1, 0.0, 0.0, 0.0, 0.0, error))
      { return false; }

    if (!set_parallelogram_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Parallelogram (w1,w3) offset (0.1, 0.1, 0.0)",
         flag_diag02_in_envelope, true,
         0.0, 0.0, 0.0, 0.1, 0.1, 0.0, error))
      { return false; }

    if (!set_parallelogram_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Skew parallelogram (w0,w2) offset (0.1, 0.1, 0.2)",
         true, flag_diag13_in_envelope,
         0.1, 0.1, 0.0, 0.0, 0.0, 0.0, error))
      { return false; }

    if (!set_parallelogram_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function, "Skew parallelogram (w1,w3) offset (0.1, 0.1, 0.2)",
         flag_diag02_in_envelope, true,
         0.0, 0.0, 0.0, 0.1, 0.1, 0.2, error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Skew cases.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_skew
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    const int DIM3(3);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;
    
    set_quad_coord_left_vertical_rectangle
      (0.5, w0coord, w1coord, w2coord, w3coord);
    skew_LL_corner_down(w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Left vertical rectangle, skew LL down",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }

    skew_LL_corner_up(w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Left vertical rectangle, skew LL up",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }

    set_quad_coord_right_vertical_rectangle
      (0.5, w0coord, w1coord, w2coord, w3coord);
    skew_LL_corner_down(w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Right vertical rectangle, skew LL down",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }

    skew_LL_corner_up(w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Right vertical rectangle, skew LL up",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }


    set_quad_coord_low_horizontal_rectangle
      (0.5, w0coord, w1coord, w2coord, w3coord);
    skew_LL_corner_down(w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Low horizontal rectangle, skew LL down",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }


    skew_LL_corner_up(w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Low horizontal rectangle, skew LL up",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }

    set_quad_coord_high_horizontal_rectangle
      (0.5, w0coord, w1coord, w2coord, w3coord);
    skew_LL_corner_down(w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "High horizontal rectangle, skew LL down",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }

    skew_LL_corner_up(w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "High horizontal rectangle, skew LL up",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Concave projection cases.
   *  - Cases where projection of quadrilateral is concave.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_concave_projection
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    const int DIM3(3);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    const int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;
    
    set_quad_coord_concave_v0
      (w0coord, w1coord, w2coord, w3coord);

    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Projection concave at v0",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, error))
      { return false; }

    set_quad_coord_concave_v1
      (w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Projection concave at v1",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, error))
      { return false; }

    set_quad_coord_concave_v2
      (w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Projection concave at v2",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, error))
      { return false; }

    set_quad_coord_concave_v3
      (w0coord, w1coord, w2coord, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "Projection concave at v3",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Boundary cases, one vertex.
   *  - Cases where one vertex is on quadrant boundaries.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_boundaryI
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);
    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;


    // Vertical rectangle tests.

    set_quad_coord_left_vertical_rectangle
      (0.5, wcoord[0], wcoord[1], wcoord[2], wcoord[3]);

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Left vertical rectangle; w0 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Left vertical rectangle; w0 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 0, error))
      { return false; }

    if (!set_on_grid_edge_multi_check_in_envelope
        (envelope_function,
         "Left vertical rectangle; w0 on grid edge",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, 0, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Left vertical rectangle; w3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 3, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Left vertical rectangle; w3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 3, 0, error))
      { return false; }

    if (!set_on_grid_edge_multi_check_in_envelope
        (envelope_function,
         "Left vertical rectangle; w3 on grid edge",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, 3, error))
      { return false; }

    
    set_quad_coord_right_vertical_rectangle
      (0.5, wcoord[0], wcoord[1], wcoord[2], wcoord[3]);

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Right vertical rectangle; w1 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Right vertical rectangle; w1 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 0, error))
      { return false; }

    if (!set_on_grid_edge_multi_check_in_envelope
        (envelope_function,
         "Right vertical rectangle; w1 on grid edge",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Right vertical rectangle; w2 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Right vertical rectangle; w2 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 0, error))
      { return false; }

    if (!set_on_grid_edge_multi_check_in_envelope
        (envelope_function,
         "Right vertical rectangle; w2 on grid edge",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, 2, error))
      { return false; }

    
    // Horizontal rectangle tests.

    set_quad_coord_low_horizontal_rectangle
      (0.5, wcoord[0], wcoord[1], wcoord[2], wcoord[3]);

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Low horizontal rectangle; w0 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Low horizontal rectangle; w0 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 0, error))
      { return false; }

    if (!set_on_grid_edge_multi_check_in_envelope
        (envelope_function,
         "Low horizontal rectangle; w0 on grid edge",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, 0, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Low horizontal rectangle; w1 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Low horizontal rectangle; w1 on x=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 0, error))
      { return false; }

    if (!set_on_grid_edge_multi_check_in_envelope
        (envelope_function,
         "Low horizontal rectangle; w1 on grid edge",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, 1, error))
      { return false; }

    set_quad_coord_high_horizontal_rectangle
      (0.5, wcoord[0], wcoord[1], wcoord[2], wcoord[3]);

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "High horizontal rectangle; w2 on y=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "High horizontal rectangle; w2 on x=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 0, error))
      { return false; }

    if (!set_on_grid_edge_multi_check_in_envelope
        (envelope_function,
         "High horizontal rectangle; w2 on grid edge",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, 2, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "High horizontal rectangle; w3 on y=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 3, 1, error))
      { return false; }

    if (!set_on_boundary_multi_check_in_envelope
        (envelope_function,
         "High horizontal rectangle; w3 on x=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 3, 0, error))
      { return false; }

    if (!set_on_grid_edge_multi_check_in_envelope
        (envelope_function,
         "High horizontal rectangle; w3 on grid edge",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, 3, error))
      { return false; }

    
    // Change z-coord.

    set_quad_coord_vertical_rectangle
      (0.4, 1.5, 0.5, wcoord[0], wcoord[1], wcoord[2], wcoord[3]);
    
    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w0 on z=0 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 0, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w0 on z=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 1, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w1 on z=0 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 0, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w1 on z=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 1, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w2 on z=0 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 0, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w2 on z=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 1, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w3 on z=0 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 3, 0, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w3 on z=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 3, 1, error))
      { return false; }
    
    set_quad_coord_vertical_rectangle
      (0.5, 1.6, 0.5, wcoord[0], wcoord[1], wcoord[2], wcoord[3]);

        if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w0 on z=0 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 0, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w0 on z=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 1, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w1 on z=0 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 0, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w1 on z=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 1, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w2 on z=0 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 0, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w2 on z=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 1, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w3 on z=0 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 3, 0, error))
      { return false; }

    if (!set_z_multi_check_in_envelope
        (envelope_function,
         "w3 on z=1 plane",
         wcoord, v0coord, v1coord, 
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 3, 1, error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Boundary cases, two vertices.
   *  - Cases where two vertices are on quadrant boundaries.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_boundaryII
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;
    
    set_quad_coord_symmetric_square
      (0.5, wcoord[0], wcoord[1], wcoord[2], wcoord[3]);
    
    if (!setII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices w0 and w1 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 1, 1, error))
      { return false; }

    if (!setII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices w1 and w2 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 2, 0, error))
      { return false; }

    if (!setII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices w2 and w3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 2, 3, 1, error))
      { return false; }

    if (!setII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices w0 and w3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 3, 0, error))
      { return false; }

    if (!setII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices w0 and w2 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 2, 1, error))
      { return false; }

    if (!setII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices w0 and w2 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 0, 2, 0, error))
      { return false; }

    if (!setII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices w1 and w3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 3, 1, error))
      { return false; }

    if (!setII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices w1 and w3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, 1, 3, 0, error))
      { return false; }

    return true;
  }

  
  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Boundary cases, three vertices.
   *  - Cases where three vertices are on quadrant boundaries.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_boundaryIII
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_IN_QUADRILATERAL(4);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_IN_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;

    // Horizontal trapezoid. Short bottom.
    set_quad_coord_horizontal_trapezoid
      (w0coord, w1coord, w2coord, w3coord, 0.4, 0.6);

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v1 and v2 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, 0, 1, 2, 1, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v1, v2 and v3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, 1, 2, 3, 1, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v2 and v3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, 0, 2, 3, 1, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v1 and v3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, 0, 1, 3, 1, error))
      { return false; }


    // Horizontal trapezoid. Long bottom.
    set_quad_coord_horizontal_trapezoid
      (w0coord, w1coord, w2coord, w3coord, 0.6, 0.4);
    
    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v1 and v2 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, 0, 1, 2, 1, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v1, v2 and v3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, 1, 2, 3, 1, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v1 and v3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, 0, 1, 3, 1, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v2 and v3 on y=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, 0, 2, 3, 1, error))
      { return false; }


    // Vertical trapezoid. Short left.
    set_quad_coord_vertical_trapezoid
      (w0coord, w1coord, w2coord, w3coord, 0.4, 0.6);

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v1 and v2 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, 0, 1, 2, 0, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v1, v2 and v3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, 1, 2, 3, 0, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v1 and v3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, 0, 1, 3, 0, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v2 and v3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, 0, 2, 3, 0, error))
      { return false; }


    // Vertical trapezoid. Long left.
    set_quad_coord_vertical_trapezoid
      (w0coord, w1coord, w2coord, w3coord, 0.6, 0.4);

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v1 and v2 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, 0, 1, 2, 0, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v1, v2 and v3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, 1, 2, 3, 0, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v1 and v3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, 0, 1, 3, 0, error))
      { return false; }

    if (!setIII_on_boundary_multi_check_in_envelope
        (envelope_function,
         "Vertices v0, v2 and v3 on x=1 plane",
         wcoord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, 0, 2, 3, 0, error))
      { return false; }


    return true;
  }

  
  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Two vertices on cube edge.
   *  - Cases where two vertices are on cube edges.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_cube_edgeII
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    const int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;

    set_quad_coord_symmetric_square
      (0.5, w0coord, w1coord, w2coord, w3coord);
    set_coord_3D(1, 1, 0.25, w0coord);
    set_coord_3D(1, 1, 0.75, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "w0 and w3 on grid edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }

    set_quad_coord_symmetric_square
      (0.5, w0coord, w1coord, w2coord, w3coord);
    set_coord_3D(1, 1, 0.25, w0coord);
    set_coord_3D(1, 1, 0.75, w1coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "w0 and w1 on grid edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, false, error))
      { return false; }

    set_quad_coord_symmetric_square
      (0.5, w0coord, w1coord, w2coord, w3coord);
    set_coord_3D(1, 1, 0.25, w0coord);
    set_coord_3D(1, 1, 0.75, w2coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "w0 and w2 on grid edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, false, error))
      { return false; }

    set_quad_coord_symmetric_square
      (0.5, w0coord, w1coord, w2coord, w3coord);
    set_coord_3D(1, 1, 0.25, w1coord);
    set_coord_3D(1, 1, 0.75, w3coord);
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "w1 and w3 on grid edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         false, true, error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    if (!test_are_quad_diagonals_in_envelope_rectangle
        <CTYPEW,CTYPEV>
        (envelope_function, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_parallelogram
        <CTYPEW,CTYPEV>
        (envelope_function, false, error))
      { return false; }
    
    if (!test_are_quad_diagonals_in_envelope_skew
        <CTYPEW,CTYPEV>
        (envelope_function, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_concave_projection
        <CTYPEW,CTYPEV>
        (envelope_function, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_boundaryI
        <CTYPEW,CTYPEV>
        (envelope_function, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_boundaryII
        <CTYPEW,CTYPEV>
        (envelope_function, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_boundaryIII
        <CTYPEW,CTYPEV>
        (envelope_function, error))
      { return false; }

    // Add test_are_quad_diagonals_in_envelope_boundaryIV

    return true;
  }


  // **************************************************
  // QUAD ORTHOGONAL TO DUAL GRID EDGE TESTS
  // **************************************************

  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope - 
   *    Quad in plane orthogonal to dual grid edge.
   *  - Cases where quad is orthogonal to dual grid edge.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_quad_orth_to_dual_edge
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_IN_QUADRILATERAL(4);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_IN_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;
    std::vector<CTYPEW> zcoord;

    zcoord.push_back(0.5);
    zcoord.push_back(0);
    zcoord.push_back(1);
    zcoord.push_back(0.001);
    zcoord.push_back(0.999);

    if (!set_left_vertical_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "Quad orth to edge dir. Left vertical rectangle.",
         true, true, zcoord, error))
      { return false; }

    if (!set_right_vertical_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "Quad orth to edge dir. Right vertical rectangle.",
         true, true, zcoord, error))
      { return false; }


    if (!set_low_horizontal_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "Quad orth to edge dir. Low horizontal rectangle.",
         true, true, zcoord, error))
      { return false; }


    if (!set_high_horizontal_rectangle_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "Quad orth to edge dir. High horizontal rectangle.",
         true, true, zcoord, error))
      { return false; }

    if (!set_symmetric_square_multi_check_in_envelope<CTYPEW,CTYPEV>
        (envelope_function,
         "Quad orth to edge dir. Symmetric square.",
         true, true, zcoord, error))
      { return false; }

    // TO BE CONTINUED

    return true;
  }

  
  // **************************************************
  // NUMERICAL ENVELOPE TESTS
  // **************************************************

  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope -
   *    Diagonals through grid edge.
   *  - Cases where two vertices are very near cube edges.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_coplanar
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    const int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;

    set_quad_coord_symmetric_square
      (0.5, w0coord, w1coord, w2coord, w3coord);
    
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function, "Symmetric square",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, error))
      { return false; }

    w0coord[0] = w0coord[1] = 2.0/3.0;
    w2coord[0] = w2coord[1] = 4.0/3.0;

    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function, "Diagonals through edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, true, error))
      { return false; }

    return true;
  }


  /*!
   *  @brief Test are_both_diagonals_in_quad_edge_envelope -
   *    Two vertices very near cube edge.
   *  - Cases where two vertices are very near cube edges.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param flag_only_projection_test If true, only use envelope
   *    projection test.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_near_cube_edgeII
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const bool flag_only_projection_test,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 1, 1, 0 };
    const CTYPEV v1coord[DIM3] = { 1, 1, 1 };
    const int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;
    const CTYPEW offset = envelope_function.TranslationAddend();
    bool flag_diag02_in_envelope, flag_diag13_in_envelope;

    // xLT1 is the largest number less than 1.
    const CTYPEW xLT1 =
      std::nextafter(CTYPEW(1.0), CTYPEW(0.0));
    
    // xGT1 is the largest number greater than 1.
    const CTYPEW xGT1 =
      std::nextafter(CTYPEW(1.0), CTYPEW(2.0));

    error.SetPrecision((std::numeric_limits<CTYPEW>::digits10)+2);

    set_quad_coord_symmetric_square
      (0.5, w0coord, w1coord, w2coord, w3coord);
    set_coord_3D(xLT1, xLT1, 0.5, w0coord);
    set_coord_3D(xGT1, xGT1, 0.5, w2coord);

    set_flag_diag_in_envelope(offset, flag_diag13_in_envelope);
    if (flag_only_projection_test)
      { flag_diag13_in_envelope = false; }
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "w0 and w2 very near grid edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, flag_diag13_in_envelope, error))
      { return false; }

    set_quad_coord_symmetric_square
      (0.5, w0coord, w1coord, w2coord, w3coord);
    set_coord_3D(xGT1, xLT1, 0.5, w1coord);
    set_coord_3D(xLT1, xGT1, 0.5, w3coord);

    set_flag_diag_in_envelope(offset, flag_diag02_in_envelope);
    if (flag_only_projection_test)
      { flag_diag02_in_envelope = false; }
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "w1 and w3 very near grid edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, true, error))
      { return false; }

    error.UnsetPrecision();

    return true;
  }

  
  /*!
   *  @brief Testing routine for are diagonals in quad edge envelope -
   *    Two vertices very near cube edge.
   *  - Cases where two vertices are very near cube edges.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   *  @param flag_only_projection_test If true, only use envelope
   *    projection test.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_near_cube_edgeIIB
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   const bool flag_only_projection_test,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUADRILATERAL(4);

    CTYPEW w0coord[DIM3];
    CTYPEW w1coord[DIM3];
    CTYPEW w2coord[DIM3];
    CTYPEW w3coord[DIM3];
    CTYPEW * wcoord[NUM_VERT_PER_QUADRILATERAL] =
      { w0coord, w1coord, w2coord, w3coord };
    const CTYPEV v0coord[DIM3] = { 0, 0, 0 };
    const CTYPEV v1coord[DIM3] = { 0, 0, 1 };
    const int edge_dir = 2;
    const bool flag_quad_pos_orient = true;
    const int iquadrant_w0 = 0;
    const CTYPEW offset = envelope_function.TranslationAddend();
    bool flag_diag02_in_envelope, flag_diag13_in_envelope;

    // xLT0 is the largest number less than 0.
    const CTYPEW xLT0 =
      std::nextafter(CTYPEW(0.0), CTYPEW(-1.0));
    
    // xGT1 is the largest number greater than 0.
    const CTYPEW xGT0 =
      std::nextafter(CTYPEW(0.0), CTYPEW(1.0));

    error.SetPrecision((std::numeric_limits<CTYPEW>::digits10)+2);

    set_coord_3D(xLT0, xLT0, 0.5, w0coord);
    set_coord_3D(0.5, -0.5, 0.5, w1coord);
    set_coord_3D(xGT0, xGT0, 0.5, w2coord);
    set_coord_3D(-0.5, 0.5, 0.5, w3coord);

    set_flag_diag_in_envelope(offset, flag_diag13_in_envelope);
    if (flag_only_projection_test)
      { flag_diag13_in_envelope = false; }
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "w0 and w2 very near grid edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         true, flag_diag13_in_envelope, error))
      { return false; }

    set_coord_3D(-0.5, -0.5, 0.5, w0coord);
    set_coord_3D(xGT0, xLT0, 0.5, w1coord);
    set_coord_3D(0.5, 0.5, 0.5, w2coord);
    set_coord_3D(xLT0, xGT0, 0.5, w3coord);

    set_flag_diag_in_envelope(offset, flag_diag02_in_envelope);
    if (flag_only_projection_test)
      { flag_diag02_in_envelope = false; }
    if (!multi_check_are_quad_diagonals_in_envelope
        (envelope_function,
         "w1 and w3 very near grid edge",
         w0coord, w1coord, w2coord, w3coord, v0coord, v1coord,
         edge_dir, flag_quad_pos_orient, iquadrant_w0,
         flag_diag02_in_envelope, true, error))
      { return false; }

    error.UnsetPrecision();

    return true;
  }
  
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_are_quad_diagonals_in_envelope_numeric
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,
   IJK::ERROR & error)
  {
    if (!test_are_quad_diagonals_in_envelope_coplanar
        <CTYPEW, CTYPEV> (envelope_function, error))
      { return false; }
        
    if (!test_are_quad_diagonals_in_envelope_near_cube_edgeII
        <CTYPEW, CTYPEV> (envelope_function, false, error))
      { return false; }
    
    if (!test_are_quad_diagonals_in_envelope_near_cube_edgeIIB
        <CTYPEW, CTYPEV> (envelope_function, false, error))
      { return false; }

    return true;
  }


  // **************************************************
  // TEST ENVELOPE PROJECTION TESTS
  // **************************************************

  /*!
   *  @brief Testing routine for projection based are diagonals in quad edge envelope.
   *  - Return true if routine passes all tests.
   *  - Return false and set error message if routines fails a test.
   *  @param envelope_function Envelope function class to determine
   *    if both quad diagonals are in quad-edge envelope.
   */
  template <typename CTYPEW, typename CTYPEV,
            typename ENVELOPE_FUNCTION_TYPE>
  bool test_projection_are_quad_diagonals_in_envelope
  (const ENVELOPE_FUNCTION_TYPE & envelope_function,   
   IJK::ERROR & error)
  {

    if (!test_are_quad_diagonals_in_envelope_parallelogram
        <CTYPEW,CTYPEV>
        (envelope_function, true, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_cube_edgeII
        <CTYPEW,CTYPEV>
        (envelope_function, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_concave_projection
        <CTYPEW,CTYPEV>
        (envelope_function, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_near_cube_edgeII
        <CTYPEW, CTYPEV>
        (envelope_function, true, error))
      { return false; }

    if (!test_are_quad_diagonals_in_envelope_near_cube_edgeIIB
        <CTYPEW, CTYPEV>
        (envelope_function, true, error))
      { return false; }

    return true;
  }

  

}

#endif
