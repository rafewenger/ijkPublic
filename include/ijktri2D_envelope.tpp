/*!
 *  \file ijktri2D_envelope.tpp
 *  @brief Determine if quad diagonals are within envelope formed 
 *    by four quad edges and dual grid edge.
 *  - Version 0.4.1
 *  \anchor envelopeDoc
 *  - Envelope is union of four tetrahedron where each tetrahedron
 *    is the convex hull of each quad edge and the grid edge dual 
 *    to the quad.
 *  - See "Intersection free contouring on an octree grid" by Ju and Udeshi
 *    Proceedings of Pacific Graphics, 2006.
 */



/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2023-2024 Rephael Wenger

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


#ifndef _IJKTRI2D_ENVELOPE_
#define _IJKTRI2D_ENVELOPE_

#include "ijkcoord.tpp"
#include "ijkquadrant.tpp"


// *** DEBUG ***
#include "ijkprint.tpp"
/*
bool flag_debug = false;
*/


namespace IJK {

  // *****************************************************************
  //! @name Determine projected triangle orientation.
  // *****************************************************************

  ///@{

  /*!
   *  @brief Return true if determinant of projection of three points is
   *    greater than or equal to zero.
   *  @tparam DET_TYPE Determinant type.
   *    - Type used for any internal variables in calculating
   *      the determinant.
   */
  template <typename DET_TYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename DTYPE0, typename DTYPE1>
  inline bool is_determinant_projected_point_ge_zero
  (const CTYPE0 p0[], const CTYPE1 p1[], const CTYPE2 p2[],
   const DTYPE0 d0, const DTYPE1 d1)
  {
    const int DIM2(2);
    const DET_TYPE v0coord[DIM2] = { p1[d0] - p0[d0], p1[d1] - p0[d1] };
    const DET_TYPE v1coord[DIM2] = { p2[d0] - p0[d0], p2[d1] - p0[d1] };
    double D;

    IJK::determinant_2x2(v0coord, v1coord, D);
    if (D >= 0) { return true; }
    else { return false; }
  }


  /*!
   *  @brief Return true if determinant of projection of three points is
   *    less than or equal to zero.
   *  @tparam DET_TYPE Determinant type.
   *    - Type used for any internal variables in calculating
   *      the determinant.
   */
  template <typename DET_TYPE,
            typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename DTYPE0, typename DTYPE1>
  inline bool is_determinant_projected_point_le_zero
  (const CTYPE0 p0[], const CTYPE1 p1[], const CTYPE2 p2[],
   const DTYPE0 d0, const DTYPE1 d1)
  {
    const int DIM2(2);
    const DET_TYPE v0coord[DIM2] = { p1[d0] - p0[d0], p1[d1] - p0[d1] };
    const DET_TYPE v1coord[DIM2] = { p2[d0] - p0[d0], p2[d1] - p0[d1] };
    double D;

    IJK::determinant_2x2(v0coord, v1coord, D);
    if (D <= 0) { return true; }
    else { return false; }
  }


  /*!
   *  @brief Return true if projection of point p equals 
   *    projection of point q.
   */
  template <typename CTYPEP, typename CTYPEQ,
            typename DIR_TYPE1, typename DIR_TYPE2>
  bool does_point_projection_equal
  (const CTYPEP pcoord[],
   const CTYPEQ qcoord[],
   const DIR_TYPE1 d1,
   const DIR_TYPE2 d2)
  {
    return ((pcoord[d1] == qcoord[d1]) &&
            (pcoord[d2] == qcoord[d2]));
  }
  
  ///@}


  // *****************************************************************
  //! @name Simple envelope routines
  // *****************************************************************

  ///@{
  
  /*!
   *  @brief Return true if quad diagonal (w0coord[], w2coord[]) 
   *    is in interior of envelope.
   *  - Simple version with no check for diagonal "deep"
   *    in envelope interior.
   *  - Quad has positive orientation around grid edge.
   *  - Expected sign of determinant sign is +1 if (w0,w2) is in envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad vertices in counterclockwise order around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1>
  inline bool is_pos_orient_quad_diagonal02_in_quad_edge_envelope_3D_simple
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3])
  {
    if (IJK::is_determinant_point_ge_zero_3D<double>
        (edge_endpoint1_coord, w3coord, w0coord, w2coord))
      { return false; }

    if (IJK::is_determinant_point_ge_zero_3D<double>
        (w3coord, edge_endpoint0_coord, w0coord, w2coord))
      { return false; }

    if (IJK::is_determinant_point_ge_zero_3D<double>
        (edge_endpoint0_coord, w1coord, w0coord, w2coord))
      { return false; }

    if (IJK::is_determinant_point_ge_zero_3D<double>
        (w1coord, edge_endpoint1_coord, w0coord, w2coord))
      { return false; }

    return true;
  }


  /*!
   *  @brief Return true if quad diagonal (w0coord[], w2coord[]) 
   *    is in interior of envelope.
   *  - Simple version with no check for diagonal "deep"
   *    in envelope interior.
   *  - Quad has positive orientation around grid edge.
   *  - Expected sign of determinant sign is -1 if (w0,w2) is in envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad vertices in clockwise order around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1>
  inline bool is_neg_orient_quad_diagonal02_in_quad_edge_envelope_3D_simple
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3])
  {
    // Switch w0 and w2 to reverse sign.
    return is_pos_orient_quad_diagonal02_in_quad_edge_envelope_3D_simple
      (w2coord, w1coord, w0coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord);
  }
  
  
  /*!
   *  @brief Return true if quad diagonal (w0coord[], w2coord[]) 
   *    is in interior of envelope.
   *  - Simple version with no check for diagonal "deep"
   *    in envelope interior.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad vertices in clockwise or counterclockwise
   *      order around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param orientation Orientation around grid edge. 
   *    - Expected determinant sign if (w0, w2) is in envelope.
   *    - Positive orientation: +1.
   *    - Negative orientation: -1.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1>
  inline bool is_quad_diagonal02_in_quad_edge_envelope_3D_simple
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const bool flag_quad_pos_orient)
  {
    if (flag_quad_pos_orient) {
      return is_pos_orient_quad_diagonal02_in_quad_edge_envelope_3D_simple
        (w0coord, w1coord, w2coord, w3coord,
         edge_endpoint0_coord, edge_endpoint1_coord);
    }
    else {
      return is_neg_orient_quad_diagonal02_in_quad_edge_envelope_3D_simple
        (w0coord, w1coord, w2coord, w3coord,
         edge_endpoint0_coord, edge_endpoint1_coord);
    }
  }

  
  /*!
   *  @brief Return true if quad diagonal (w1coord[], w3coord[]) 
   *    is in interior of envelope.
   *  - Simple version with no check for diagonal "deep"
   *    in envelope interior.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad vertices in clockwise or counterclockwise
   *      order around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param orientation Orientation around grid edge. 
   *    - Expected determinant sign if (w0, w2) is in envelope.
   *    - Positive orientation: +1.
   *    - Negative orientation: -1.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1>
  inline bool is_quad_diagonal13_in_quad_edge_envelope_3D_simple
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const bool flag_quad_pos_orient)
  {
    return is_quad_diagonal02_in_quad_edge_envelope_3D_simple
      (w1coord, w2coord, w3coord, w0coord,
       edge_endpoint0_coord, edge_endpoint1_coord, flag_quad_pos_orient);
  }

  
  /*!
   *  @brief Return true if both quad diagonals are in interior of envelope.
   *  - Simple version with no check for diagonal "deep"
   *    in envelope interior.
   *  - Quad diagonals are (w0coord[],w2coord[]) and (w1coord[],w3coord[]).
   *  - Version that returns is_diag02_in_envelope and is_diag13_in_envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad has positive or negative orientation around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *    @pre edge_endpoint_coord0[] and edge_endpoint_coord1[] differ
   *       only in the coordinate in edge_direction.
   *  @param[out] is_diag02_in_envelope True if diag(quad_vert[0],quad_vert[2])
   *    is in envelope.
   *  @param[out] is_diag13_in_envelope True if diag(quad_vert[1],quad_vert[3])
   *    is in envelope.
   *   - Function returns (is_diag02_in_envelope && is_diag13_in_envelope).
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1>
  bool are_both_quad_diagonals_in_quad_edge_envelope_3D_simple
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const bool flag_quad_pos_orient,
   bool & is_diag02_in_envelope,
   bool & is_diag13_in_envelope)
  {
    is_diag02_in_envelope =
      is_quad_diagonal02_in_quad_edge_envelope_3D_simple
      (w0coord, w1coord, w2coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       flag_quad_pos_orient);

    is_diag13_in_envelope =
      is_quad_diagonal13_in_quad_edge_envelope_3D_simple
      (w0coord, w1coord, w2coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       flag_quad_pos_orient);

    return (is_diag02_in_envelope && is_diag13_in_envelope);
  }
  
  ///@}


  // *****************************************************************
  //! @name Determine if diagonal is deep in envelope.
  // *****************************************************************

  namespace {

    template <typename CTYPEW, 
              typename DTYPE0, typename DTYPE1, typename DTYPE2,
              typename CTYPET0, typename CTYPET1, typename CTYPET2,
              typename CTYPEWII>
    void translate_point_3D
      (const CTYPEW wcoord[3],
       const DTYPE0 d0, const DTYPE1 d1, const DTYPE2 d2,
       const CTYPET0 trans0, const CTYPET1 trans1, const CTYPET2 trans2,
       CTYPEWII wcoordII[3])
    {
      wcoordII[d0] = wcoord[d0] + trans0;
      wcoordII[d1] = wcoord[d1] + trans1;
      wcoordII[d2] = wcoord[d2] + trans2;
    }
    
        
    template <typename CTYPEW0, typename CTYPEW2,
              typename DTYPE0, typename DTYPE1, typename DTYPE2,
              typename CTYPET0, typename CTYPET1, typename CTYPET2,
              typename CTYPEW0II, typename CTYPEW2II>
    void translate_two_points_3D
      (const CTYPEW0 w0coord[3], const CTYPEW2 w2coord[3],
       const DTYPE0 d0, const DTYPE1 d1, const DTYPE2 d2,
       const CTYPET0 trans0, const CTYPET1 trans1, const CTYPET2 trans2,
       CTYPEW0II w0coordII[3], CTYPEW2II w2coordII[3])
    {
      translate_point_3D
        (w0coord, d0, d1, d2, trans0, trans1, trans2, w0coordII);
      translate_point_3D
        (w2coord, d0, d1, d2, trans0, trans1, trans2, w2coordII);
    }


    /*!
     *  @brief Get translation addend from quadrantA and quadrantC
     *    toward quadrantB.
     *  - Quadrants are positively oriented (CCW) around center point
     *    in order quadA, quadB, quadC.
     *  @param iquadrantA Index (0,1,2 or 3) of quadrant corresponding to quadA.
     *  @param translation_addend Add/subtract translation_addend to x,y.
     *    @pre translation_addend >= 0.
     *  @param[out] transx Translation in x-direction.
     *  @param[out] transsy Translation in y-direction.
     */
    template <typename ITYPE, typename CTYPET,
              typename CTYPEX, typename CTYPEY>
    void get_pos_orient_translation_quadrantA_quadrantC_to_quadrantB
    (const ITYPE iquadrantA, const CTYPET translation_addend,
     CTYPEX & transx, CTYPEY & transy)
    {
      const int NUM_QUADRANTS(4);
      static const int translation_signx[NUM_QUADRANTS] =
        { 1, 1, -1, -1 };
      static const int translation_signy[NUM_QUADRANTS] =
        { -1, 1, 1, -1 };

      transx = translation_signx[iquadrantA]*translation_addend;
      transy = translation_signy[iquadrantA]*translation_addend;
    }
    
  }

  
  ///@{

  /*!
   *  @brief Return true if quad diagonal (w0, w2) is deep in envelope .
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  - Quad has positive orientation around grid edge.
   *  - Expected sign of determinant sign is -1 if (w0,w2) is in envelope.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param edge_direction Direction (0, 1 or 2) of edge.
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - See \ref quadrantDoc.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE, typename ITYPE,
            typename CTYPET>
  bool is_pos_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const DIR_TYPE edge_direction,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend)
  {
    const int DIM3(3);
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;
    CTYPEW0 w0coordII[DIM3];
    CTYPEW2 w2coordII[DIM3];

    if (does_point_projection_equal
        (w0coord, edge_endpoint0_coord, d1, d2) ||
        does_point_projection_equal
        (w2coord, edge_endpoint0_coord, d1, d2)) {
      // w0 or w2 is on the dual grid edge.
      
      if (does_point_projection_equal
          (w1coord, edge_endpoint0_coord, d1, d2) ||
          does_point_projection_equal
          (w3coord, edge_endpoint0_coord, d1, d2)) {
        // Two adjacent vertices are on the dual grid edge.
        // Diagonal (w0,w2) is on the envelope boundary.
        return false;
      }

      return true;
    }

    CTYPET trans1, trans2;
    get_pos_orient_translation_quadrantA_quadrantC_to_quadrantB
      (iquadrant_w0, translation_addend, trans1, trans2);

    // Check (w0,w2) against (endpoint1,w3).
    translate_two_points_3D
      (w0coord, w2coord, edge_direction, d1, d2,
       translation_addend, -trans1, -trans2,
       w0coordII, w2coordII);
    if (IJK::is_determinant_point_ge_zero_3D<double>
        (edge_endpoint1_coord, w3coord, w0coordII, w2coordII))
      { return false; }


    // Check (w0,w2) against (w3,endpoint0).
    w0coordII[edge_direction] = w0coord[edge_direction] - translation_addend;
    w2coordII[edge_direction] = w2coord[edge_direction] - translation_addend;

    if (IJK::is_determinant_point_ge_zero_3D<double>
        (w3coord, edge_endpoint0_coord, w0coordII, w2coordII))
      { return false; }

    // Check (w0,w2) against (endpoint0,w1).
    translate_two_points_3D
      (w0coord, w2coord, edge_direction, d1, d2,
       -translation_addend, trans1, trans2,
       w0coordII, w2coordII);
    if (IJK::is_determinant_point_ge_zero_3D<double>
        (edge_endpoint0_coord, w1coord, w0coordII, w2coordII))
      { return false; }

    // Check (w0,w2) against (w1,endpoint1).
    w0coordII[edge_direction] = w0coord[edge_direction] + translation_addend;
    w2coordII[edge_direction] = w2coord[edge_direction] + translation_addend;
    if (IJK::is_determinant_point_ge_zero_3D<double>
        (w1coord, edge_endpoint1_coord, w0coordII, w2coordII))
      { return false; }

    return true;
  }


  /*!
   *  @brief Return true if quad diagonal (w0, w2) is deep in envelope .
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  - Quad has negative orientation around grid edge.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param edge_direction Direction (0, 1 or 2) of edge.
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - See \ref quadrantDoc.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE, typename ITYPE, typename CTYPET>
  bool is_neg_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const DIR_TYPE edge_direction,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend)
  {
    const ITYPE iquadrant_w2 =
      diagonally_opposite_quadrant(iquadrant_w0);

    // Switch w0 and w2 to reverse sign.
    return is_pos_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
      (w2coord, w1coord, w0coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       edge_direction, iquadrant_w2, translation_addend);
  }


  /*!
   *  @brief Return true if quad diagonal (w0, w2) is deep in envelope .
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad vertices in clockwise or counterclockwise
   *      order around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param edge_direction Direction (0, 1 or 2) of edge.
   *  @param flag_quad_pos_orient
   *    If true, quadrilateral is positively oriented around dual edge.
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - See \ref quadrantDoc.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE, typename ITYPE,  typename CTYPET>
  bool is_quad_diagonal02_deep_in_quad_edge_envelope_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend)
  {
    if (flag_quad_pos_orient) {
      return is_pos_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
        (w0coord, w1coord, w2coord, w3coord,
         edge_endpoint0_coord, edge_endpoint1_coord,
         edge_direction, iquadrant_w0, translation_addend);
    }
    else {
      return is_neg_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
        (w0coord, w1coord, w2coord, w3coord,
         edge_endpoint0_coord, edge_endpoint1_coord,
         edge_direction, iquadrant_w0, translation_addend);
    }
  }


  /*!
   *  @brief Return true if projection of quad diagonal (w0, w2) 
   *    is deep in projection of quad.
   *  - Quad has positive orientation around grid edge.
   *    (Projection is CCW around center.)
   *  - Expected sign of determinant sign is +1 if (w0,w2) is in envelope.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param edge_direction Direction (0, 1 or 2) of edge.
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - See \ref quadrantDoc.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename DIR_TYPE, typename ITYPE,
            typename CTYPET>
  bool is_projected_pos_orient_quad_diagonal02_deep_in_quad_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const DIR_TYPE edge_direction,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend)
  {
    const int DIM3(3);
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;
    CTYPEW0 w0coordII[DIM3];
    CTYPEW2 w2coordII[DIM3];

    CTYPET trans1, trans2;
    get_pos_orient_translation_quadrantA_quadrantC_to_quadrantB
      (iquadrant_w0, translation_addend, trans1, trans2);

    // Check projection of (w0,w2) against projection of w1.
    translate_two_points_3D
      (w0coord, w2coord, edge_direction, d1, d2,
       0, trans1, trans2, w0coordII, w2coordII);
    if (is_determinant_projected_point_ge_zero<double>
        (w0coordII, w2coordII, w1coord, d1, d2))
      { return false; }
    

    // Check projection of (w0,w2) against projection of w3.
    translate_two_points_3D
      (w0coord, w2coord, edge_direction, d1, d2,
       0, -trans1, -trans2, w0coordII, w2coordII);
    if (is_determinant_projected_point_le_zero<double>
        (w0coordII, w2coordII, w3coord, d1, d2))
      { return false; }
    
    return true;
  }


  /*!
   *  @brief Return true if quad diagonal (w0, w2) is deep in envelope 
   *    or if (w0,w2) is deep in relative interior of a quad 
   *    that is in plane orthogonal to its dual grid edge.
   *  - Quad has positive orientation around grid edge.
   *  - Note: May return true for quad diagonal on envelope boundary,
   *    when quad is in plane orthogonal to its dual grid edge.
   *    the edge endpoints.
   *  - Note: This routine may have return different results
   *    than is_quad_diagonal02_in_quad_edge_envelope_3D().
   *  - In particular, it could allow diagonals on grid facets.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param edge_direction Direction (0, 1 or 2) of edge.
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - See \ref quadrantDoc.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   *  @param flag_process_quad_orth_to_dual_edge Special processing 
   *    if quadrilateral lies in a plane orthogonal to the dual grid edge.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE, typename ITYPE, typename CTYPET>
  bool is_pos_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const DIR_TYPE edge_direction,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge)
  {
    if (flag_process_quad_orth_to_dual_edge) {

      const CTYPEW0 w0cz = w0coord[edge_direction];
      if ((w1coord[edge_direction] == w0cz) &&
          (w2coord[edge_direction] == w0cz) &&
          (w3coord[edge_direction] == w0cz)) {

        // Isosurface quadrilateral is in a plane orthogonal
        //   to the dual grid edge.
        return is_projected_pos_orient_quad_diagonal02_deep_in_quad_3D
          (w0coord, w1coord, w2coord, w3coord,
           edge_direction, iquadrant_w0, translation_addend);
      }
    }

    // Either no special processing of quad orthogonal to dual grid edge,
    // or quad is not orthogonal to grid edge.
    return is_pos_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
      (w0coord, w1coord, w2coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       edge_direction, iquadrant_w0,
       translation_addend);
  }


  /*!
   *  @brief Return true if quad diagonal (w0, w2) is deep in envelope 
   *    or if (w0,w2) is deep in relative interior of a quad 
   *    that is in plane orthogonal to its dual grid edge.
   *  - Quad has negative orientation around grid edge.
   *  - Note: May return true for quad diagonal on envelope boundary,
   *    when quad is in plane orthogonal to its dual grid edge.
   *    the edge endpoints.
   *  - Note: This routine may have return different results
   *    than is_quad_diagonal02_in_quad_edge_envelope_3D().
   *  - In particular, it could allow diagonals on grid facets.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param edge_direction Direction (0, 1 or 2) of edge.
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - See \ref quadrantDoc.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   *  @param flag_process_quad_orth_to_dual_edge Special processing 
   *    if quadrilateral lies in a plane orthogonal to the dual grid edge.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE, typename ITYPE, typename CTYPET>
  bool is_neg_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const DIR_TYPE edge_direction,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge)
  {
    const ITYPE iquadrant_w2 =
      diagonally_opposite_quadrant(iquadrant_w0);

    // Switch w0 and w2 to reverse sign.
    return is_pos_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
      (w2coord, w1coord, w0coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       edge_direction, iquadrant_w2, translation_addend,
       flag_process_quad_orth_to_dual_edge);
  }

  
  /*!
   *  @brief Return true if quad diagonal (w0, w2) is deep in envelope 
   *    or if (w0,w2) is deep in relative interior of a quad 
   *    that is in plane orthogonal to its dual grid edge.
   *  - Note: May return true for quad diagonal on envelope boundary,
   *    when quad is in plane orthogonal to its dual grid edge.
   *    the edge endpoints.
   *  - Note: This routine may have return different results
   *    than is_quad_diagonal02_in_quad_edge_envelope_3D().
   *  - In particular, it could allow diagonals on grid facets.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad has positive or negative orientation around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param edge_direction Direction (0, 1 or 2) of edge.
   *  @param flag_quad_pos_orient
   *    If true, quadrilateral is positively oriented around dual edge.
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - See \ref quadrantDoc.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   *  @param flag_process_quad_orth_to_dual_edge Special processing 
   *    if quadrilateral lies in a plane orthogonal to the dual grid edge.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE, typename ITYPE, typename CTYPET>
  bool is_quad_diagonal02_deep_in_quad_edge_envelope_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge)
  {
    if (flag_quad_pos_orient) {
      return is_pos_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
      (w0coord, w1coord, w2coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       edge_direction, iquadrant_w0, translation_addend,
       flag_process_quad_orth_to_dual_edge);
    }
    else {
      return is_neg_orient_quad_diagonal02_deep_in_quad_edge_envelope_3D
      (w0coord, w1coord, w2coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       edge_direction, iquadrant_w0, translation_addend,
       flag_process_quad_orth_to_dual_edge);
    }

  }



  /*!
   *  @brief Return true if quad diagonal (w1, w3) is deep in envelope 
   *    or if (w1,w3) is deep in relative interior of a quad 
   *    that is in plane orthogonal to its dual grid edge.
   *  - Note: May return true for quad diagonal on envelope boundary,
   *    when quad is in plane orthogonal to its dual grid edge.
   *    the edge endpoints.
   *  - Note: This routine may have return different results
   *    than is_quad_diagonal02_in_quad_edge_envelope_3D().
   *  - In particular, it could allow diagonals on grid facets.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad has positive or negative orientation around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *  @param edge_direction Direction (0, 1 or 2) of edge.
   *  @param flag_quad_positive_orientation
   *    If true, quadrilateral is positively oriented around dual edge.
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - See \ref quadrantDoc.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   *  @param flag_process_quad_orth_to_dual_edge Special processing 
   *    if quadrilateral lies in a plane orthogonal to the dual grid edge.
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE, typename ITYPE,
            typename CTYPET>
  bool is_quad_diagonal13_deep_in_quad_edge_envelope_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge)
  {
    const ITYPE iquadrant_w1 =
      next_quadrant(flag_quad_pos_orient, iquadrant_w0);
    
    // Start with w1coord so that (w1coord, w3coord)
    //   becomes diagonal (0,2).
    return is_quad_diagonal02_deep_in_quad_edge_envelope_3D
      (w1coord, w2coord, w3coord, w0coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w1,
       translation_addend, flag_process_quad_orth_to_dual_edge);
  }


  /*!
   *  @brief Return true if both quad diagonals are deep in envelope 
   *    or if diagonals are deep in relative interior of a quad 
   *    that is in plane orthogonal to its dual grid edge.
   *  - Quad diagonals are (w0coord[],w2coord[]) and (w1coord[],w3coord[]).
   *  - Version that returns is_diag02_in_envelope and is_diag13_in_envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param w0coord[] Coordinates of first quad vertex.
   *  @param w1coord[] Coordinates of second quad vertex.
   *  @param w2coord[] Coordinates of third quad vertex.
   *  @param w3coord[] Coordinates of fourth quad vertex.
   *    - Quad has positive or negative orientation around grid edge.
   *  @param edge_endpoint_coord0[] Coordinates of lower/leftmost
   *         endpoint of grid edge dual to isosurface quadrilatral.
   *  @param edge_endpoint_coord1[] Coordinates of upper/rightmost
   *         endpoint of grid edge dual to isosurface quadrilateral
   *    @pre edge_endpoint_coord0[] and edge_endpoint_coord1[] differ
   *       only in the coordinate in edge_direction.
   *  @param edge_direction Direction of edge dual
   *    to isosurface quadrilateral.
   *  @param[out] is_diag02_in_envelope True if diag(quad_vert[0],quad_vert[2])
   *    is in envelope.
   *  @param[out] is_diag13_in_envelope True if diag(quad_vert[1],quad_vert[3])
   *    is in envelope.
   *   - Function returns (is_diag02_in_envelope && is_diag13_in_envelope).
   */
  template <typename CTYPEW0, typename CTYPEW1,
            typename CTYPEW2, typename CTYPEW3,
            typename CTYPEP0, typename CTYPEP1,
            typename DIR_TYPE, typename ITYPE,
            typename CTYPET>
  bool are_both_quad_diagonals_deep_in_quad_edge_envelope_3D
  (const CTYPEW0 w0coord[3],
   const CTYPEW1 w1coord[3],
   const CTYPEW2 w2coord[3],
   const CTYPEW3 w3coord[3],
   const CTYPEP0 edge_endpoint0_coord[3],
   const CTYPEP1 edge_endpoint1_coord[3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge,
   bool & is_diag02_in_envelope,
   bool & is_diag13_in_envelope)
  {
    is_diag02_in_envelope =
      is_quad_diagonal02_deep_in_quad_edge_envelope_3D
      (w0coord, w1coord, w2coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       translation_addend, flag_process_quad_orth_to_dual_edge);

    is_diag13_in_envelope =
      is_quad_diagonal13_deep_in_quad_edge_envelope_3D
      (w0coord, w1coord, w2coord, w3coord,
       edge_endpoint0_coord, edge_endpoint1_coord,
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       translation_addend, flag_process_quad_orth_to_dual_edge);

    return (is_diag02_in_envelope && is_diag13_in_envelope);
  }


  /*!
   *  @overload
   *  @brief Return true if both quad diagonals are deep in envelope 
   *    or if diagonals are deep in relative interior of a quad 
   *    that is in plane orthogonal to its dual grid edge.
   *  - Quad diagonals are (quad_vert[0],quad_vert[2]) and (quad_vert[1],quad_vert[3]).
   *  - Version that returns is_diag02_in_envelope and is_diag13_in_envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param quad_vert[] Four quad vertices in clockwise or counter-clockwise
   *    order around grid edge.
   */
  template <typename WTYPE, typename COORD_TYPEW, typename COORD_TYPEV,
            typename DIR_TYPE, typename ITYPE, typename CTYPET>
  bool are_both_quad_diagonals_deep_in_quad_edge_envelope_3D
  (const WTYPE quad_vert[],
   const COORD_TYPEV * const edge_endpoint_coord[2],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const COORD_TYPEW * coord_array,
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge,
   bool & is_diag02_in_envelope,
   bool & is_diag13_in_envelope)
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUAD(4);
    const COORD_TYPEW * const wcoord[NUM_VERT_PER_QUAD] = 
      { coord_array + quad_vert[0]*DIM3,
        coord_array + quad_vert[1]*DIM3,
        coord_array + quad_vert[2]*DIM3,
        coord_array + quad_vert[3]*DIM3 };

    return are_both_quad_diagonals_deep_in_quad_edge_envelope_3D
      (wcoord[0], wcoord[1], wcoord[2], wcoord[3],
       edge_endpoint_coord[0], edge_endpoint_coord[1],
       edge_direction, flag_quad_pos_orient, iquadrant_w0,
       translation_addend, flag_process_quad_orth_to_dual_edge,
       is_diag02_in_envelope, is_diag13_in_envelope);
  }


  /*!
   *  @overload
   *  @brief Return true if both quad diagonals are deep in envelope 
   *    or if diagonals are deep in relative interior of a quad 
   *    that is in plane orthogonal to its dual grid edge.
   *  - Version using argument edge_endpoint_coord[2][3].
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param edge_endpoint_coord[] Coordinates of grid edge endpoints
   *         to isosurface quadrilateral.
   *    - edge_endpoint_coord[0] is coordinates of lower/leftmost endpoint.
   *    @pre edge_endpoint_coord[0] and edge_endpoint_coord[1] differ
   *       only in the coordinate in edge_direction.
   */
  template <typename WTYPE, typename COORD_TYPEW, typename COORD_TYPEV,
            typename DIR_TYPE, typename ITYPE, typename CTYPET>
  bool are_both_quad_diagonals_deep_in_quad_edge_envelope_3D
  (const WTYPE quad_vert[],
   const COORD_TYPEV edge_endpoint_coord_3D[2][3],
   const DIR_TYPE edge_direction,
   const bool flag_quad_pos_orient,
   const ITYPE iquadrant_w0,
   const COORD_TYPEW coord_array[],
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge,
   bool & is_diag02_in_envelope,
   bool & is_diag13_in_envelope)
  {
    const COORD_TYPEV * edge_endpoint_coord[2] = 
      { edge_endpoint_coord_3D[0], edge_endpoint_coord_3D[1] };

    return (are_both_quad_diagonals_deep_in_quad_edge_envelope_3D
            (quad_vert, edge_endpoint_coord, edge_direction,
             flag_quad_pos_orient, iquadrant_w0, coord_array, 
             translation_addend, flag_process_quad_orth_to_dual_edge,
             is_diag02_in_envelope, is_diag13_in_envelope));
  }

  ///@}


  // *****************************************************************
  //! @name Compute quad diagonals deeply in envelope interior.
  // *****************************************************************

  ///@{
  
  /*!
   *  @brief Determine which isosurface quadrilaterals are deeply
   *    in envelope interior.
   *  - Stores flags in quad_info[] indicating diagonals deeply 
   *    inside envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre Vertices are in dimension 3.
   *  @param translation_addend Add/subtract translation_addend to x,y,z,
   *    to diagonal to determine if diagonal is deep in envelope.
   *    @pre translation_addend >= 0.
   *  @param flag_process_quad_orth_to_dual_edge Special processing 
   *    if quadrilateral lies in a plane orthogonal to the dual grid edge.
   *  @param[out] quad_info[i] Information about isosurface_quadrilateral i.
   *    - This routine sets quad_info[i].is_diagonal_in_envelope[0]
   *      and quad_info[i].is_diagonal_in_envelope[1].
   *  @param[out] num_quads_with_diag_outside_envelope Number of isosurface
   *    quadrilaterals with one or both diagonals outside envelope.
   */
  template <typename GRID_TYPE, typename VTYPEQ,
            typename DUAL_ISOVERT_TYPE, typename CTYPEV,
            typename ISOPOLY_INFO_TYPE,
            typename CTYPET, typename NTYPEQ, typename NTYPEQ2>
  void compute_quad_diagonals_deep_in_envelope
  (const GRID_TYPE & grid,
   const VTYPEQ quad_vert[], const NTYPEQ num_quad,
   const DUAL_ISOVERT_TYPE isov_list[],
   const CTYPEV coord_array[],
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge,
   ISOPOLY_INFO_TYPE quad_info[], 
   NTYPEQ2 & num_quads_with_diag_outside_envelope)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE GRID_VTYPE;

    const int NUM_VERT_PER_QUAD(4);
    const int DIM3(3);
    GRID_VTYPE edge_endpoint[2];
    int edge_direction;
    CTYPEV edge_endpoint_coord[2][DIM3];
    IJK::PROCEDURE_ERROR error
      ("compute_quad_diagonals_deep_in_envelope");

    // Initialize.
    num_quads_with_diag_outside_envelope = 0;

    if (num_quad == 0) {
      // Nothing to do.
      return;
    }

    if (!IJK::check_array_non_empty(quad_vert, "quad_vert", error))
      { throw error; }    
    if (!IJK::check_array_non_empty(isov_list, "isov_list", error))
      { throw error; }
    if (!IJK::check_array_non_empty(coord_array, "coord_array", error))
      { throw error; }
    if (!IJK::check_array_non_empty(quad_info, "quad_info", error))
      { throw error; }
    
    for (NTYPEQ iquad = 0; iquad < num_quad; iquad++) {
      const VTYPEQ * quad_vert_i = quad_vert + iquad*NUM_VERT_PER_QUAD;
      bool is_diag02_in_envelope, is_diag13_in_envelope;

      edge_endpoint[0] = quad_info[iquad].GridEdgeEndpoint0();
      edge_direction = quad_info[iquad].GridEdgeDirection();
      edge_endpoint[1] = 
        grid.NextVertex(edge_endpoint[0], edge_direction);

      grid.ComputeCoord(edge_endpoint[0], edge_endpoint_coord[0]);
      grid.ComputeCoord(edge_endpoint[1], edge_endpoint_coord[1]);

      bool flag_quad_pos_orient;
      int iquadrant_w0;
      get_quad_orientation_and_quadrant_containing_w0_3D
        (grid, quad_vert_i, isov_list, edge_direction,
         flag_quad_pos_orient, iquadrant_w0);
      
      if (!are_both_quad_diagonals_deep_in_quad_edge_envelope_3D
          (quad_vert_i, edge_endpoint_coord, edge_direction,
           flag_quad_pos_orient, iquadrant_w0, coord_array, 
           translation_addend, flag_process_quad_orth_to_dual_edge,
           is_diag02_in_envelope, is_diag13_in_envelope)) {

        num_quads_with_diag_outside_envelope++;
      }

      quad_info[iquad].is_diagonal_in_envelope[0] = is_diag02_in_envelope;
      quad_info[iquad].is_diagonal_in_envelope[1] = is_diag13_in_envelope;
    }
  }


  /*!
   *  @overload
   *  @brief Determine which isosurface quadrilaterals are deeply
   *    in envelope interior.
   *  - Stores flags in quad_info[] indicating diagonals deeply 
   *    inside envelope.
   *  - Version using C++ vector for quad_vert[], quad_info[] 
   *    and coord_array[].
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   */
  template <typename GRID_TYPE, typename VTYPEQ,
            typename DUAL_ISOVERT_TYPE, typename CTYPE,
            typename ISOPOLY_INFO_TYPE,
            typename CTYPET, typename NTYPEQ>
  void compute_quad_diagonals_deep_in_envelope
  (const GRID_TYPE & grid, 
   const std::vector<VTYPEQ> & quad_vert,
   const std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   const std::vector<CTYPE> & coord_array,
   const CTYPET translation_addend,
   const bool flag_process_quad_orth_to_dual_edge,
   std::vector<ISOPOLY_INFO_TYPE> & quad_info, 
   NTYPEQ & num_quads_with_diag_outside_envelope)
  {
    typedef typename std::vector<VTYPEQ>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;
    IJK::PROCEDURE_ERROR error("compute_quad_diagonals_outside_envelope");

    if (quad_info.size() != num_quad) {
      error.AddMessage("Programming error. Incorrect size of quad_info.");
      error.AddMessage
        ("  quad_info.size() = ", quad_info.size(), "");
      error.AddMessage
        ("  Should equal ", num_quad, ", number of quads.");
      
      throw error;
    }

    if (num_quad == 0) {
      // Nothing to do.
      return;
    }

    compute_quad_diagonals_deep_in_envelope
      (grid, IJK::vector2pointer(quad_vert), num_quad,
       IJK::vector2pointer(isov_list),
       IJK::vector2pointer(coord_array),
       translation_addend, flag_process_quad_orth_to_dual_edge,
       IJK::vector2pointer(quad_info),
       num_quads_with_diag_outside_envelope);
  }


  ///@}

  
  // *****************************************************************
  //! @name ENVELOPE FUNCTION OBJECT
  // *****************************************************************

  ///@{

  /*!
   *  @brief Base class for envelope function class.
   *  - Default routines returning default values.
   */
  template <typename CTYPE>
  class ENVELOPE_BASE {

  public:

    CTYPE TranslationAddend() const
    { return 0; }

    bool FlagProcessQuadOrthToDualEdge() const
    { return false; }
    
    template <typename CTYPEW0, typename CTYPEW1,
              typename CTYPEW2, typename CTYPEW3,
              typename CTYPEP0, typename CTYPEP1,
              typename DIR_TYPE, typename ITYPE>
    bool are_both_quad_diagonals_in_quad_edge_envelope
    (const CTYPEW0 w0coord[3],
     const CTYPEW1 w1coord[3],
     const CTYPEW2 w2coord[3],
     const CTYPEW3 w3coord[3],
     const CTYPEP0 edge_endpoint0_coord[3],
     const CTYPEP1 edge_endpoint1_coord[3],
     const DIR_TYPE edge_direction,
     const bool flag_quad_pos_orient,
     const ITYPE iquadrant_w0,
     bool & is_diag02_in_envelope,
     bool & is_diag13_in_envelope) const;
  };


  /*!
   *  @brief Class that calls simple envelope routines.
   */
  template <typename CTYPE>
  class SIMPLE_ENVELOPE:public ENVELOPE_BASE<CTYPE> {

  public:
    
    /// Constructor.
    SIMPLE_ENVELOPE() {}

    template <typename CTYPEW0, typename CTYPEW1,
              typename CTYPEW2, typename CTYPEW3,
              typename CTYPEP0, typename CTYPEP1,
              typename DIR_TYPE, typename ITYPE>
    bool are_both_quad_diagonals_in_quad_edge_envelope
    (const CTYPEW0 w0coord[3],
     const CTYPEW1 w1coord[3],
     const CTYPEW2 w2coord[3],
     const CTYPEW3 w3coord[3],
     const CTYPEP0 edge_endpoint0_coord[3],
     const CTYPEP1 edge_endpoint1_coord[3],
     const DIR_TYPE edge_direction,
     const bool flag_quad_pos_orient,
     const ITYPE iquadrant_w0,
     bool & is_diag02_in_envelope,
     bool & is_diag13_in_envelope) const
    {
      return are_both_quad_diagonals_in_quad_edge_envelope_3D_simple
        (w0coord, w1coord, w2coord, w3coord,
         edge_endpoint0_coord, edge_endpoint1_coord,
         flag_quad_pos_orient,
         is_diag02_in_envelope, is_diag13_in_envelope);
      
      return true;
    };

  };

  
  /*!
   *  @brief Class that calls deep envelope routines.
   */
  template <typename CTYPE>
  class DEEP_ENVELOPE:public ENVELOPE_BASE<CTYPE> {

  public:
    CTYPE translation_addend;
    bool flag_process_quad_orth_to_dual_edge;

    void Clear() {
      translation_addend = 0.0;
      flag_process_quad_orth_to_dual_edge = false;
    };
  
  public:
    
    /// Constructor.
    DEEP_ENVELOPE() { Clear(); }

    CTYPE TranslationAddend() const
    { return translation_addend; }

    bool FlagProcessQuadOrthToDualEdge() const
    { return flag_process_quad_orth_to_dual_edge; }

    template <typename CTYPEW0, typename CTYPEW1,
              typename CTYPEW2, typename CTYPEW3,
              typename CTYPEP0, typename CTYPEP1,
              typename DIR_TYPE, typename ITYPE>
    bool are_both_quad_diagonals_in_quad_edge_envelope
    (const CTYPEW0 w0coord[3],
     const CTYPEW1 w1coord[3],
     const CTYPEW2 w2coord[3],
     const CTYPEW3 w3coord[3],
     const CTYPEP0 edge_endpoint0_coord[3],
     const CTYPEP1 edge_endpoint1_coord[3],
     const DIR_TYPE edge_direction,
     const bool flag_quad_pos_orient,
     const ITYPE iquadrant_w0,
     bool & is_diag02_in_envelope,
     bool & is_diag13_in_envelope) const
    {
      return are_both_quad_diagonals_deep_in_quad_edge_envelope_3D
        (w0coord, w1coord, w2coord, w3coord,
         edge_endpoint0_coord, edge_endpoint1_coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w0,
         translation_addend, flag_process_quad_orth_to_dual_edge,
         is_diag02_in_envelope, is_diag13_in_envelope);
    };

  };

  
  /*!
   *  @brief Class encapsulating envelope function object.
   *  - Checks at class creation that ENVELOPE_TYPE
   *    has desired member function.
   */
  template <typename ENVELOPE_TYPE>
  class ENVELOPE_FUNCTION:public ENVELOPE_TYPE {

  public:
    
    template <typename CTYPEW0, typename CTYPEW1,
              typename CTYPEW2, typename CTYPEW3,
              typename CTYPEP0, typename CTYPEP1,
              typename DIR_TYPE, typename ITYPE>
    inline bool are_both_quad_diagonals_in_quad_edge_envelope
    (const CTYPEW0 w0coord[3],
     const CTYPEW1 w1coord[3],
     const CTYPEW2 w2coord[3],
     const CTYPEW3 w3coord[3],
     const CTYPEP0 edge_endpoint0_coord[3],
     const CTYPEP1 edge_endpoint1_coord[3],
     const DIR_TYPE edge_direction,
     const bool flag_quad_pos_orient,
     const ITYPE iquadrant_w0,
     bool & is_diag02_in_envelope,
     bool & is_diag13_in_envelope) const
    {
      return ENVELOPE_TYPE::are_both_quad_diagonals_in_quad_edge_envelope
        (w0coord, w1coord, w2coord, w3coord,
         edge_endpoint0_coord, edge_endpoint1_coord,
         edge_direction, flag_quad_pos_orient, iquadrant_w0,
         is_diag02_in_envelope, is_diag13_in_envelope);
    }

  };


  ///@}


}

#endif

