/*!
 *  @file ijkcheck_geom.tpp
 *  @brief ijk templates for checking geometric objects.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2023 Rephael Wenger

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

#ifndef _IJKCHECK_GEOM_
#define _IJKCHECK_GEOM_

#include "ijk.tpp"
#include "ijkquadrant.tpp"


namespace IJK {

  // *****************************************************************
  //! @name Add error message.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Add edge endpoint0 coord to error message.
   */
  template <typename CTYPEP, typename DIR_TYPE>
  void add_edge_endpoint0_coord_to_error_message
  (const CTYPEP edge_endpoint0_coord[],
   const DIR_TYPE d,
   IJK::ERROR & error)
  {
    error.AddArrayElementMessage
      ("  ", "edge_endpoint0_coord", edge_endpoint0_coord, d, ".");
  }

  
  /*!
   *  @brief Add edge endpoint1 coord to error message.
   */
  template <typename CTYPEP, typename DIR_TYPE>
  void add_edge_endpoint1_coord_to_error_message
  (const CTYPEP edge_endpoint1_coord[],
   const DIR_TYPE d,
   IJK::ERROR & error)
  {
    error.AddArrayElementMessage
      ("  ", "edge_endpoint1_coord", edge_endpoint1_coord, d, ".");
  }

    
  /*!
   *  @brief Add coord of both endpoints to error message.
   */
  template <typename CTYPEP0, typename CTYPEP1, typename DIR_TYPE>
  void add_edge_endpoint_coord_to_error_message
  (const CTYPEP0 edge_endpoint0_coord[],
   const CTYPEP1 edge_endpoint1_coord[],
   const DIR_TYPE d,
   IJK::ERROR & error)
  {
    add_edge_endpoint0_coord_to_error_message
      (edge_endpoint0_coord, d, error);
    add_edge_endpoint1_coord_to_error_message
      (edge_endpoint1_coord, d, error);
  }


  /*!
   *  @brief Add quadrant error message.
   */
  template <typename ITYPEW, typename ITYPEQ,
            typename CTYPEW, typename CTYPEP, typename DTYPE>
  void add_quadrant_error_message
  (const ITYPEW iw,
   const CTYPEW wcoord[],
   const CTYPEP edge_endpoint0_coord[],
   const DTYPE d,
   const ITYPEQ iquadrant_w,
   IJK::ERROR & error)
  {
    error.AddMessage
      ("Quadrilateral vertex w", iw, " is not in quadrant ",
       iquadrant_w, ".");
    error.AddMessage
      ("  wcoord", iw, "[", d, "] = ", wcoord[d], ".");
    add_edge_endpoint0_coord_to_error_message
      (edge_endpoint0_coord, d, error);
  }

  
  /*!
   *  @brief Add quadrant should be less than or equal to error message.
   */
  template <typename ITYPEW, typename ITYPEQ,
            typename CTYPEW, typename CTYPEP, typename DTYPE>
  void add_quadrant_should_be_LE_error_message
  (const ITYPEW iw,
   const CTYPEW wcoord[],
   const CTYPEP edge_endpoint0_coord[],
   const DTYPE d,
   const ITYPEQ iquadrant_w,
   IJK::ERROR & error)
  {
    add_quadrant_error_message
      (iw, wcoord, edge_endpoint0_coord, d, iquadrant_w, error);
    error.AddMessage
      ("  wcoord", iw, "[", d,
       "] should be less than or equal to edge_endpoint0_coord[",
       d, "].");
  }


  /*!
   *  @brief Add quadrant should be greater than or equal to error message.
   */
  template <typename ITYPEW, typename ITYPEQ,
            typename CTYPEW, typename CTYPEP, typename DTYPE>
  void add_quadrant_should_be_GE_error_message
  (const ITYPEW iw,
   const CTYPEW wcoord[],
   const CTYPEP edge_endpoint0_coord[],
   const DTYPE d,
   const ITYPEQ iquadrant_w,
   IJK::ERROR & error)
  {
    add_quadrant_error_message
      (iw, wcoord, edge_endpoint0_coord, d, iquadrant_w, error);
    error.AddMessage
      ("  wcoord", iw, "[", d,
       "] should be greater than or equal to edge_endpoint0_coord[",
       d, "].");
  }


  /*!
   *  @brief Add quadrilateral orientation error message.
   */
  template <typename ITYPEW, typename ITYPEQ,
            typename CTYPEW, typename CTYPEP, typename DTYPE>
  void add_quadrilateral_orientation_error_message
  (const ITYPEW iw,
   const CTYPEW wcoord[],
   const CTYPEP edge_endpoint0_coord[],
   const DTYPE d,
   const ITYPEQ iquadrant_w,
   const bool flag_quad_pos_orient,
   IJK::ERROR & error)
  {
    error.AddMessage("Quadrilateral orientation error.");
    if (flag_quad_pos_orient) {
      error.AddMessage
        ("  Based on given positive quad orientation, quadrilateral vertex w", iw);
    }
    else {
      error.AddMessage
        ("  Based on given negative quad orientation, quadrilateral vertex w", iw);
    }
    
    error.AddMessage
      ("  is expected but is not in quadrant ", iquadrant_w, ".");
    error.AddMessage
      ("  wcoord", iw, "[", d, "] = ", wcoord[d], ".");
    add_edge_endpoint0_coord_to_error_message
      (edge_endpoint0_coord, d, error);
  }

  
  /*!
   *  @brief Add quadrilateral orientation error message.
   *  - Add message that if iw was in iquadrant_w,
   *    then wcoord[d] would be less than or equal 
   *    to edge_endpoint0_coord[d].
   */
  template <typename ITYPEW, typename ITYPEQ,
            typename CTYPEW, typename CTYPEP, typename DTYPE>
  void add_quadrilateral_orientation_LE_error_message
  (const ITYPEW iw,
   const CTYPEW wcoord[],
   const CTYPEP edge_endpoint0_coord[],
   const DTYPE d,
   const ITYPEQ iquadrant_w,
   const bool flag_quad_pos_orient,
   IJK::ERROR & error)
  {
    add_quadrilateral_orientation_error_message
      (iw, wcoord, edge_endpoint0_coord, d,
       iquadrant_w, flag_quad_pos_orient, error);
    error.AddMessage
      ("  If vertex w", iw, " was in quadrant ", iquadrant_w,
       ", then wcoord", iw, "[", d, "]");
    error.AddMessage
      ("  would be less than or equal to edge_endpoint0_coord[",
       d, "].");
  }


  /*!
   *  @brief Add quadrilateral orientation error message.
   *  - Add message that if iw was in iquadrant_w,
   *    then wcoord[d] would be greater than or equal 
   *    to edge_endpoint0_coord[d].
   */
  template <typename ITYPEW, typename ITYPEQ,
            typename CTYPEW, typename CTYPEP, typename DTYPE>
  void add_quadrilateral_orientation_GE_error_message
  (const ITYPEW iw,
   const CTYPEW wcoord[],
   const CTYPEP edge_endpoint0_coord[],
   const DTYPE d,
   const ITYPEQ iquadrant_w,
   const bool flag_quad_pos_orient,
   IJK::ERROR & error)
  {
    add_quadrilateral_orientation_error_message
      (iw, wcoord, edge_endpoint0_coord, d,
       iquadrant_w, flag_quad_pos_orient, error);
    error.AddMessage
      ("  If vertex w", iw, " was in quadrant ", iquadrant_w,
       ", then wcoord", iw, "[", d, "]");
    error.AddMessage
      ("  would be less than or equal to edge_endpoint0_coord[",
       d, "].");
  }
  
  
  // *****************************************************************
  //! @name Check grid edge.
  // *****************************************************************
  
  /*!
   *  @brief Check edge direction.
   */
  template <typename DTYPE, typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE>
  bool check_edge_direction
  (const DTYPE dimension,
   const CTYPEP0 edge_endpoint0_coord[],
   const CTYPEP1 edge_endpoint1_coord[],
   const DIR_TYPE edge_direction,
   IJK::ERROR & error)
  {
    const CTYPEP0 e0d0 =
      edge_endpoint0_coord[edge_direction];
    const CTYPEP1 e1d0 =
      edge_endpoint1_coord[edge_direction];

    if (e1d0 != e0d0 + 1) {

      if (e0d0 == e1d0) {
        error.AddMessage
          ("Incorrect edge direction ", edge_direction, ".");
        error.AddMessage
          ("  edge_endpoint0_coord[", edge_direction,
           "] = edge_endpoint1_coord[", edge_direction, "] = ",
           e0d0, ".");
        return false;
      }
      else if (e0d0 > e1d0) {
        error.AddMessage
          ("Incorrect order of edge endpoints in direction ",
           edge_direction, ".");
        error.AddMessage
          ("  edge_endpoint0_coord[", edge_direction,
           "] > edge_endpoint1_coord[", edge_direction, "].");
        error.AddMessage
          ("  edge_endpoint0_coord[", edge_direction,
           "] should be less than edge_endpoint1_coord[",
           edge_direction, "].");
        add_edge_endpoint_coord_to_error_message
          (edge_endpoint0_coord, edge_endpoint1_coord,
           edge_direction, error);
        return false;
      }
      else {
        error.AddMessage("Incorrect upper/right edge endpoint.");
        add_edge_endpoint_coord_to_error_message
          (edge_endpoint0_coord, edge_endpoint1_coord,
           edge_direction, error);
        error.AddMessage
          ("  edge_endpoint1[", edge_direction,
           "] should be ", e0d0+1, ".");
        return false;
      }
    }

    for (int d = 0; d < dimension; d++) {

      if (d == edge_direction) { continue; }

      const CTYPEP0 e0d = edge_endpoint0_coord[d];
      const CTYPEP1 e1d = edge_endpoint1_coord[d];

      if (e0d != e1d) {
        error.AddMessage
          ("Incorrect edge endpoint coordinates.");
        error.AddMessage
          ("  Endpoint1 is not in direction ",
           edge_direction, " from endpoint0.");
        add_edge_endpoint_coord_to_error_message
          (edge_endpoint0_coord, edge_endpoint1_coord,
           d, error);
        error.AddMessage
          ("  edge_endpoint1[", d,
           "] should be ", e0d, ".");
        return false;
      }
    }

    return true;
  }


  // *****************************************************************
  //! @name Check quadrilateral.
  // *****************************************************************

  /*!
   *  @brief Check quadrant containing vertex is consistent
   *    with vertex and edge endpoint0 coordinates.
   *  @param iw Vertex index.
   *  @param flag_orientation_error If true, error is caused
   *     by incorrect quadrilateral orientation.
   */
  template <typename ITYPEW, typename ITYPEQ,
            typename CTYPEW, typename CTYPEP,
            typename DIR_TYPE>
  bool check_quadrant_containing_vertex_3D
  (const ITYPEW iw,
   const CTYPEW wcoord[3],
   const CTYPEP edge_endpoint0_coord[3],
   const DIR_TYPE edge_direction,
   const ITYPEQ iquadrant_w,
   const bool flag_quad_pos_orient,
   const bool flag_orientation_error,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;
    
    if (iquadrant_w < 2) {
      // iquadrant_w is 0 or 1.
      if (wcoord[d2] > edge_endpoint0_coord[d2]) {
        if (flag_orientation_error) {
          add_quadrilateral_orientation_LE_error_message
            (iw, wcoord, edge_endpoint0_coord, d2,
             iquadrant_w, flag_quad_pos_orient, error);
        }
        else {
          add_quadrant_should_be_LE_error_message
            (iw, wcoord, edge_endpoint0_coord, d2, iquadrant_w, error);
        }
        return false;
      }
    }
    else {
      // iquadrant_w is 2 or 3.
      if (wcoord[d2] < edge_endpoint0_coord[d2]) {
        if (flag_orientation_error) {
          add_quadrilateral_orientation_GE_error_message
            (iw, wcoord, edge_endpoint0_coord, d2,
             iquadrant_w, flag_quad_pos_orient, error);
        }
        else {
          add_quadrant_should_be_GE_error_message
            (iw, wcoord, edge_endpoint0_coord, d2, iquadrant_w, error);
        }
        return false;
      }
    }

    if ((iquadrant_w == 0) || (iquadrant_w == 3)) {
      if (wcoord[d1] > edge_endpoint0_coord[d1]) {
        if (flag_orientation_error) {
          add_quadrilateral_orientation_LE_error_message
            (iw, wcoord, edge_endpoint0_coord, d1,
             iquadrant_w, flag_quad_pos_orient, error);
        }
        else {
          add_quadrant_should_be_LE_error_message
            (iw, wcoord, edge_endpoint0_coord, d1, iquadrant_w, error);
        }
        return false;
      }
    }
    else {
      // iquadrant_w is 2 or 3.
      if (wcoord[d1] < edge_endpoint0_coord[d1]) {
        if (flag_orientation_error) {
          add_quadrilateral_orientation_GE_error_message
            (iw, wcoord, edge_endpoint0_coord, d1,
             iquadrant_w, flag_quad_pos_orient, error);
        }
        else {
          add_quadrant_should_be_GE_error_message
            (iw, wcoord, edge_endpoint0_coord, d1, iquadrant_w, error);
        }
        return false;
      }
    }

    return true;
  }

  
  /*!
   *  @brief Check that quadrilateral coordinates are consistent with given
   *    quad orientation and quadrant containing w0.
   *  @pre edge_endpoint1_coord[edge_direction] == 
   *    edge_endpoint0_coord[edge_direction] + 1.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP,
            typename DIR_TYPE, typename ITYPE>
  bool check_dual_quadrilateral_orientation_and_quadrants_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP edge_endpoint0_coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   IJK::ERROR & error)
  {
    const int DIM3(3);
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;
    const ITYPE iquadrant_w1 =
      IJK::next_quadrant(flag_quad_pos_orient, iquadrant_w0);
    const ITYPE iquadrant_w2 =
      IJK::next_quadrant(flag_quad_pos_orient, iquadrant_w1);
    const ITYPE iquadrant_w3 =
      IJK::next_quadrant(flag_quad_pos_orient, iquadrant_w2);

    if (!check_quadrant_containing_vertex_3D
        (0, w0coord, edge_endpoint0_coord, edge_direction,
         iquadrant_w0, flag_quad_pos_orient, false, error))
      { return false; }

    if (!check_quadrant_containing_vertex_3D
        (1, w1coord, edge_endpoint0_coord, edge_direction,
         iquadrant_w1, flag_quad_pos_orient, true, error))
      { return false; }

    if (!check_quadrant_containing_vertex_3D
        (2, w2coord, edge_endpoint0_coord, edge_direction,
         iquadrant_w2, flag_quad_pos_orient, true, error))
      { return false; }

    if (!check_quadrant_containing_vertex_3D
        (3, w3coord, edge_endpoint0_coord, edge_direction,
         iquadrant_w3, flag_quad_pos_orient, true, error))
      { return false; }

    return true;
  }

  ///@}

}

#endif
