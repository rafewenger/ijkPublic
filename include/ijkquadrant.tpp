/*!
 *  \file ijkquadrant.tpp
 *  @brief Determine quadrant containing a point.
 *  - Version 0.4.1
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


#ifndef _IJKQUADRANT_
#define _IJKQUADRANT_

// Needed for routing IJK::get_isovert_grid_cube_index.
// *** Not a very logical dependency, so maybe should modify in the future.
#include "ijkisopoly.tpp"

namespace IJK {

  namespace {

    /*!
     *  @brief Get orientation of cubes around grid edge and
     *    quadrant containing cube0 if cube0d1 != cube1.d1.
     *  - Return true if cube0d1 != cube1d1.
     *  - Return false if cube0d1 == cube1d1.
     *  @param flag_quad_pos_orient
     *    If true, quad has positive orientation around dual grid edge.
     *  @param iquadrant_c0 Index (0,1,2 or 3) of quadrant
     *    containing cube0.
     */
    template <typename CTYPE01, typename CTYPE02,
              typename CTYPE11, typename CTYPE12,
              typename CTYPE21, typename CTYPE22,
              typename ITYPEQ>
    inline bool get_cube_orientation_and_quadrant_c0d1_neq_c1d1
    (const CTYPE01 cube0d1, const CTYPE02 cube0d2,
     const CTYPE11 cube1d1, const CTYPE12 cube1d2,
     const CTYPE21 cube2d1, const CTYPE22 cube2d2,
     bool & flag_quad_pos_orient,
     ITYPEQ & iquadrant_c0)
    {
      if (cube0d1 < cube1d1) {
        if (cube1d2 < cube2d2) {
          flag_quad_pos_orient = true;
          iquadrant_c0 = 0;
        }
        else {
          // cube1d2 > cube2d2.
          flag_quad_pos_orient = false;
          iquadrant_c0 = 3;
        }

        return true;
      }
      else if (cube0d1 > cube1d1) {
        if (cube1d2 < cube2d2) {
          flag_quad_pos_orient = false;
          iquadrant_c0 = 1;

        }
        else {
          // cube1d2 > cube2d2
          flag_quad_pos_orient = true;
          iquadrant_c0 = 2;
        }

        return true;
      }
      else {
        return false;
      }
    }

  }

  
  /*!
   *  @brief Get orientation of cubes around grid edge and
   *    quadrant containing cube0.
   *  @param edge_direction Direction (0, 1, or 2) of edge
   *    contained in all three cubes.
   *  @param icube0 First cube.
   *  @param icube1 Second cube.
   *    @pre icube1 shares a facet with icube0 and a facet with icube2.
   *  @param icube2 Third cube.
   *    @pre All three cubes share a grid edge.
   *  @param flag_quad_pos_orient
   *    If true, quad has positive orientation around dual grid edge.
   *  \anchor quadrantDoc
   *  @param iquadrant_c0 Index (0,1,2 or 3) of quadrant
   *    containing cube0.
   *    - Orthogonal axes are d1 and d2 where d1=(edge_direction+1)%3
   *      and d2=(edge_direction+2)%3.
   *    - Quadrant 0 is lower/leftmost quadrant
   *        (min coord[d1], min coord[d2]).
   *    - Quadrant 1 is lower/rightmost 
   *        (min coord[d1], max coord[d2]) quadrant.
   *    - Quadrant 2 is upper/rightmost quadrant.
   *    - Quadrant 3 is upper/leftmost quadrant.
   */
  template <typename GRID_TYPE,
            typename ITYPEC0, typename ITYPEC1, typename ITYPEC2,
            typename DIR_TYPE, typename ITYPEQ>
  inline void get_cube_orientation_and_quadrant_containing_cube0_3D
  (const GRID_TYPE & grid,
   const ITYPEC0 icube0,
   const ITYPEC1 icube1,
   const ITYPEC2 icube2,
   const DIR_TYPE edge_direction,
   bool & flag_quad_pos_orient,
   ITYPEQ & iquadrant_w0)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE GRID_COORD_TYPE;
    
    const int DIM3(3);
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;

    const GRID_COORD_TYPE cube0d1 = grid.CoordD(icube0, d1);
    const GRID_COORD_TYPE cube0d2 = grid.CoordD(icube0, d2);
    const GRID_COORD_TYPE cube1d1 = grid.CoordD(icube1, d1);
    const GRID_COORD_TYPE cube1d2 = grid.CoordD(icube1, d2);
    const GRID_COORD_TYPE cube2d1 = grid.CoordD(icube2, d1);
    const GRID_COORD_TYPE cube2d2 = grid.CoordD(icube2, d2);

    // Initialize
    flag_quad_pos_orient = true;
    iquadrant_w0 = 0;

    if (get_cube_orientation_and_quadrant_c0d1_neq_c1d1
        (cube0d1, cube0d2, cube1d1, cube1d2, cube2d1, cube2d2,
         flag_quad_pos_orient, iquadrant_w0)) {
      return;
    }
    else {
      // cube0d1 == cube1d1

      if (cube0d2 < cube1d2) {
        if (cube1d1 < cube2d1) {
          flag_quad_pos_orient = false;
          iquadrant_w0 = 0;
        }
        else {
          // cube1d1 > cube2d1
          flag_quad_pos_orient = true;
          iquadrant_w0 = 1;
        }
      }
      else {
        // cube0d2 > cube1d2
        if (cube1d1 < cube2d1) {
          flag_quad_pos_orient = true;
          iquadrant_w0 = 3;
        }
        else {
          // cube1d1 > cube2d1
          flag_quad_pos_orient = false;
          iquadrant_w0 = 2;
        }
      }
    }

  }

    
  /*!
   *  @brief Get orientation of quad around grid edge and
   *    quadrant containing projection of quad_vert[0].
   *  @param edge_direction Direction (0, 1, or 2) of edge dual to quadrilateral.
   *  @param[out] flag_quad_pos_orient
   *    If true, quad has positive (CCW) orientation around dual grid edge.
   *  \anchor quadrantDoc
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - Orthogonal axes are d1 and d2 where d1=(edge_direction+1)%3
   *      and d2=(edge_direction+2)%3.
   *    - Quadrant 0 is lower/leftmost quadrant
   *        (min coord[d1], min coord[d2]).
   *    - Quadrant 1 is lower/rightmost 
   *        (min coord[d1], max coord[d2]) quadrant.
   *    - Quadrant 2 is upper/rightmost quadrant.
   *    - Quadrant 3 is upper/leftmost quadrant.
   */
  template <typename GRID_TYPE, typename VTYPEQ,
            typename DUAL_ISOVERT_TYPE, typename DIR_TYPE,
            typename ITYPEQ>
  inline void get_quad_orientation_and_quadrant_containing_w0_3D
  (const GRID_TYPE & grid,
   const VTYPEQ quad_vert[],
   const DUAL_ISOVERT_TYPE isov_list[],
   const DIR_TYPE edge_direction,
   bool & flag_quad_pos_orient,
   ITYPEQ & iquadrant_w0)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE CUBE_INDEX_TYPE;

    const VTYPEQ isov0 = quad_vert[0];
    const VTYPEQ isov1 = quad_vert[1];
    const VTYPEQ isov2 = quad_vert[2];
    const CUBE_INDEX_TYPE icube0 =
      IJK::get_isovert_grid_cube_index(isov_list[isov0]);
    const CUBE_INDEX_TYPE icube1 = 
      IJK::get_isovert_grid_cube_index(isov_list[isov1]);
    const CUBE_INDEX_TYPE icube2 = 
      IJK::get_isovert_grid_cube_index(isov_list[isov2]);
    
    get_cube_orientation_and_quadrant_containing_cube0_3D
      (grid, icube0, icube1, icube2, edge_direction,
       flag_quad_pos_orient, iquadrant_w0);
  }


  /*!
   *  @brief Get orientation of cubes around grid edge and
   *    quadrant containing cube0.
   *  - Version B. Recursively call routine if 
   *      cube0.coord[d1] == cube1.coord[d1].
   *  @param edge_direction Direction (0, 1, or 2) of edge
   *    contained in all three cubes.
   *  @param icube0 First cube.
   *  @param icube1 Second cube.
   *    @pre icube1 shares a facet with icube0 and a facet with icube2.
   *  @param icube2 Third cube.
   *    @pre icube2 shares a facet with icube1 and a facet with icube3.
   *  @param icube3 Fourth cube.
   *    @pre icube3 shares a facet with icube2 and a facet with icube0.
   *    @pre All three cubes share a grid edge.
   *  @param flag_quad_pos_orient
   *    If true, quad has positive orientation around dual grid edge.
   *  \anchor quadrantDoc
   *  @param iquadrant_c0 Index (0,1,2 or 3) of quadrant
   *    containing cube0.
   *    - Orthogonal axes are d1 and d2 where d1=(edge_direction+1)%3
   *      and d2=(edge_direction+2)%3.
   *    - Quadrant 0 is lower/leftmost quadrant
   *        (min coord[d1], min coord[d2]).
   *    - Quadrant 1 is lower/rightmost 
   *        (min coord[d1], max coord[d2]) quadrant.
   *    - Quadrant 2 is upper/rightmost quadrant.
   *    - Quadrant 3 is upper/leftmost quadrant.
   */
  template <typename GRID_TYPE, typename ITYPEC0, typename ITYPEC1,
            typename ITYPEC2, typename ITYPEC3,
            typename DIR_TYPE, typename ITYPEQ>
  void get_cube_orientation_and_quadrant_containing_cube0_3D_B
  (const GRID_TYPE & grid,
   const ITYPEC0 icube0,
   const ITYPEC1 icube1,
   const ITYPEC2 icube2,
   const ITYPEC3 icube3,
   const DIR_TYPE edge_direction,
   bool & flag_quad_pos_orient,
   ITYPEQ & iquadrant_w0)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE GRID_COORD_TYPE;
    
    const int DIM3(3);
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;
    const int NUM_VERTICES_PER_QUADRILATERAL(4);

    const GRID_COORD_TYPE cube0d1 = grid.CoordD(icube0, d1);
    const GRID_COORD_TYPE cube0d2 = grid.CoordD(icube0, d2);
    const GRID_COORD_TYPE cube1d1 = grid.CoordD(icube1, d1);
    const GRID_COORD_TYPE cube1d2 = grid.CoordD(icube1, d2);
    const GRID_COORD_TYPE cube2d1 = grid.CoordD(icube2, d1);
    const GRID_COORD_TYPE cube2d2 = grid.CoordD(icube2, d2);

    // Initialize
    flag_quad_pos_orient = true;
    iquadrant_w0 = 0;

    if (get_cube_orientation_and_quadrant_c0d1_neq_c1d1
        (cube0d1, cube0d2, cube1d1, cube1d2, cube2d1, cube2d2,
         flag_quad_pos_orient, iquadrant_w0)) {
      return;
    }
    else {
      // cube0d1 == cube1d1

      const GRID_COORD_TYPE cube3d1 = grid.CoordD(icube3, d1);
      const GRID_COORD_TYPE cube3d2 = grid.CoordD(icube3, d2);

      get_cube_orientation_and_quadrant_c0d1_neq_c1d1
          (cube1d1, cube1d2, cube2d1, cube2d2, cube3d1, cube3d2,
           flag_quad_pos_orient, iquadrant_w0);

      if (flag_quad_pos_orient) {
        iquadrant_w0 =
          (iquadrant_w0 + 3)%NUM_VERTICES_PER_QUADRILATERAL;
      }
      else {
        iquadrant_w0 =
          (iquadrant_w0 + 1)%NUM_VERTICES_PER_QUADRILATERAL;
      }
    }

  }

  
  /*!
   *  @brief Get orientation of quad around grid edge and
   *    quadrant containing projection of quad_vert[0].
   *  - Version B. Recursively call routine if 
   *      cube0.coord[d1] == cube1.coord[d1].
   *  @param edge_direction Direction (0, 1, or 2) of edge dual to quadrilateral.
   *  @param flag_quad_pos_orient
   *    If true, quad has positive orientation around dual grid edge.
   *  \anchor quadrantDoc
   *  @param iquadrant_w0 Index (0,1,2 or 3) of quadrant
   *    containing projection of quad_vert[0] in direction edge_direction.
   *    - Orthogonal axes are d1 and d2 where d1=(edge_direction+1)%3
   *      and d2=(edge_direction+2)%3.
   *    - Quadrant 0 is lower/leftmost quadrant.
   *    - Quadrant 1 is next quadrant after quadrant 0
   *      in counter-clockwise order around projection of dual grid edge.
   *    - Quadrant 2 is next quadrant after quadrant 1
   *      in counter-clockwise order.
   *    - Quadrant 3 is next quadrant after quadrant 2
   *      in counter-clockwise order.
   */
  template <typename GRID_TYPE, typename VTYPEQ,
            typename DUAL_ISOVERT_TYPE, typename DIR_TYPE,
            typename ITYPEQ>
  void get_quad_orientation_and_quadrant_containing_w0_3D_B
  (const GRID_TYPE & grid,
   const VTYPEQ quad_vert[],
   const DUAL_ISOVERT_TYPE isov_list[],
   const DIR_TYPE edge_direction,
   bool & flag_quad_pos_orient,
   ITYPEQ & iquadrant_w0)
  {
    typedef typename DUAL_ISOVERT_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;

    const VTYPEQ isov0 = quad_vert[0];
    const VTYPEQ isov1 = quad_vert[1];
    const VTYPEQ isov2 = quad_vert[2];
    const VTYPEQ isov3 = quad_vert[3];
    const CUBE_INDEX_TYPE icube0 = isov_list[isov0].GridCubeIndex();
    const CUBE_INDEX_TYPE icube1 = isov_list[isov1].GridCubeIndex();
    const CUBE_INDEX_TYPE icube2 = isov_list[isov2].GridCubeIndex();
    const CUBE_INDEX_TYPE icube3 = isov_list[isov3].GridCubeIndex();

    get_cube_orientation_and_quadrant_containing_cube0_3D_B
      (grid, icube0, icube1, icube2, icube3, edge_direction,
       flag_quad_pos_orient, iquadrant_w0);
  }

    
  // *** CHANGE ORIENT TO CCW AND MAKE CCW A PARAM.
  /*!
   *  @brief Return next quadrant after iquadrant. (Positive quad orientation.)
   *  - Quadrant vertices have positive (counter-clockwise) orientation.
   */
  constexpr int next_quadrant_pos_orient(const int iquadrant)
  { return (iquadrant+1)%4; }


  // *** REPLACE WITH constexpr 
  /*!
   *  @brief Return next quadrant after iquadrant. (Negative quad orientation.)
   *  - Quadrant vertices have negative (clockwise) orientation.
   */
  template <typename ITYPE>
  inline ITYPE next_quadrant_neg_orient(const ITYPE iquadrant)
  {
    const int NUM_VERT_PER_QUADRILATERAL(4);
    const ITYPE quadrant_map[NUM_VERT_PER_QUADRILATERAL] =
      { 3, 0, 1, 2 };

    return quadrant_map[iquadrant];
  }

  
  // *** REPLACE WITH constexpr 
  /*!
   *  @brief Return next quadrant after iquadrant.
   */
  template <typename ITYPE>
  inline ITYPE next_quadrant
  (const bool flag_quad_pos_orient, const ITYPE iquadrant)
  {
    if (flag_quad_pos_orient)
      { return next_quadrant_pos_orient(iquadrant); }
    else
      { return next_quadrant_neg_orient(iquadrant); }
  }

  
  // *** REPLACE WITH constexpr 
  /*!
   *  @brief Return quadrant before iquadrant.
   */
  template <typename ITYPE>
  inline ITYPE prev_quadrant
  (const bool flag_quad_pos_orient, const ITYPE iquadrant)
  {
    if (flag_quad_pos_orient)
      { return next_quadrant_neg_orient(iquadrant); }
    else
      { return next_quadrant_pos_orient(iquadrant); }
  }
  

  // *** REPLACE WITH constexpr 
  /*!
   *  @brief Return diagonally opposite quadrant.
   *  - Note: Quad orientation does not matter.
   *  @param iquadrant Quadrant index. 0,1,2 or 3.
   *    @pre 0 <= iquadrant <= 3.
   */
  template <typename ITYPE>
  inline ITYPE diagonally_opposite_quadrant
  (const ITYPE iquadrant)
  {
    const int NUM_VERT_PER_QUADRILATERAL(4);
    const ITYPE quadrant_map[NUM_VERT_PER_QUADRILATERAL] =
      { 2, 3, 0, 1 };

    return quadrant_map[iquadrant];
  }

}

#endif

