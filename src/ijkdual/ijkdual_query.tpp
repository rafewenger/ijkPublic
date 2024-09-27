/*!
 *  @file ijkdual_query.tpp
 *  @brief Query routines for ijkdual.
 * - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2023 Rephael Wenger

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


#ifndef _IJKDUAL_QUERY_TPP_
#define _IJKDUAL_QUERY_TPP_

#include <algorithm>
#include <vector>

#include "ijk.tpp"


namespace IJKDUAL {

  // **************************************************
  /// @name GET INFORMATION ABOUT SHARED FACET
  // **************************************************

  ///@{

  /// Get orthogonal direction of facet shared by c0 and c1.
  /// @param[out] facet_orth_dir Orthogonal direction of shared facet.
  template <typename GRID_CUBE_TYPE0, typename GRID_CUBE_TYPE1,
            typename DTYPE>
  void get_shared_facet_orth_dir_3D
  (const GRID_CUBE_TYPE0 & c0, const GRID_CUBE_TYPE1 & c1,
   DTYPE & facet_orth_dir)
  {
    const DTYPE DIM3(3);

    for (DTYPE d = 0; d < DIM3; d++) {
      if (c0.CubeCoord(d) != c1.CubeCoord(d)) {
        facet_orth_dir = d;
        return;
      }
    }

    IJK::PROCEDURE_ERROR error("get_shared_facet_orth_dir_3D");
    error.AddMessage
      ("Programming error.  Cubes c0 and c1 have identical coordinates.");
    error.AddMessage
      ("  Coordinates: (", c0.CubeCoord(0), ",", c0.CubeCoord(1),
       ",", c0.CubeCoord(2), ").");
    throw error;
  }


  /// Get orthogonal direction of facet shared by cubes 
  ///   containing isov0 and isov1.
  template <typename GRID_CUBE_TYPE, typename DUAL_ISOV_TYPE,
            typename VTYPE0, typename VTYPE1,
            typename DTYPE>
  void get_shared_facet_orth_dir_of_cubes_containing_isov_3D
  (const std::vector<GRID_CUBE_TYPE> & cube_isov_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const VTYPE0 isov0, const VTYPE1 isov1,
   DTYPE & facet_orth_dir)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE index_to_c0 = isov_list[isov0].cube_list_index;
    const SIZE_TYPE index_to_c1 = isov_list[isov1].cube_list_index;
    
    get_shared_facet_orth_dir_3D
      (cube_isov_list[index_to_c0], cube_isov_list[index_to_c1],
       facet_orth_dir);
  }


  /// Get orthogonal direction and coordinate of facet shared by c0 and c1.
  /// @param[out] facet_orth_dir Orthogonal direction of shared facet.
  /// @param[out] orth_dir_coord Coordinate orth_dir_coord of facet.
  template <typename GRID_CUBE_TYPE0, typename GRID_CUBE_TYPE1,
            typename DTYPE, typename CTYPE>
  void get_shared_facet_orth_dir_coord_3D
  (const GRID_CUBE_TYPE0 & c0, const GRID_CUBE_TYPE1 & c1,
   DTYPE & facet_orth_dir, CTYPE & orth_dir_coord)
  {
    get_shared_facet_orth_dir_3D(c0, c1, facet_orth_dir);
    
    orth_dir_coord = std::max(c0.CubeCoord(facet_orth_dir), 
                              c1.CubeCoord(facet_orth_dir));
  }

  /// Get orthogonal direction and coordinate of facet shared 
  ///   by cubes containing isov0 and isov1.
  template <typename GRID_CUBE_TYPE, typename DUAL_ISOV_TYPE,
            typename VTYPE0, typename VTYPE1,
            typename DTYPE, typename CTYPE>
  void get_shared_facet_orth_dir_coord_of_cubes_containing_isov_3D
  (const std::vector<GRID_CUBE_TYPE> & cube_isov_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const VTYPE0 isov0, const VTYPE1 isov1,
   DTYPE & facet_orth_dir,
   CTYPE & orth_dir_coord)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    SIZE_TYPE index_to_c0 = isov_list[isov0].cube_list_index;
    SIZE_TYPE index_to_c1 = isov_list[isov1].cube_list_index;
    
    get_shared_facet_orth_dir_coord_3D
      (cube_isov_list[index_to_c0], cube_isov_list[index_to_c1],
       facet_orth_dir, orth_dir_coord);
  }

  /// Get orthogonal direction and coordinate of facet shared by c0 and c1.<br>
  /// Get also orientation from c0 to c1.
  /// @param[out] facet_orth_dir Orthogonal direction of shared facet.
  /// @param[out] orth_dir_coord Coordinate orth_dir_coord of facet.
  /// @param[out] orient Orientation of c0 and c1.
  ///   - If c0.CubeCoord(facet_orth_dir) < c1.CubeCoord(facet_orth_dir), 
  ///       then orient = 1.
  ///   - If c0.CubeCoord(facet_orth_dir) > c1.CubeCoord(facet_orth_dir), 
  ///       then orient = -1.
  /// @pre c0.CubeCoord(facet_orth_dir) != c1.CubeCoord(facet_orth_dir).
  template <typename GRID_CUBE_TYPE0, typename GRID_CUBE_TYPE1,
            typename DTYPE, typename CTYPE, typename ORIENT_TYPE>
  void get_shared_facet_coord_and_orient_3D
  (const GRID_CUBE_TYPE0 & c0, const GRID_CUBE_TYPE1 & c1,
   DTYPE & facet_orth_dir, CTYPE & orth_dir_coord,
   ORIENT_TYPE & orient)
  {
    const DTYPE DIM3(3);

    for (DTYPE d = 0; d < DIM3; d++) {
      if (c0.CubeCoord(d) < c1.CubeCoord(d)) {
        facet_orth_dir = d;
        orth_dir_coord = c1.CubeCoord(d);
        orient = 1;
        return;
      }
      else if (c0.CubeCoord(d) > c1.CubeCoord(d)) {
        facet_orth_dir = d;
        orth_dir_coord = c0.CubeCoord(d);
        orient = -1;
        return;
      }
    }

    IJK::PROCEDURE_ERROR error
      ("get_shared_facet_coord_and_orient_3D");
    error.AddMessage
      ("Programming error.  Cubes c0 and c1 have identical coordinates.");
    error.AddMessage
      ("  Coordinates: (", c0.CubeCoord(0), ",", c0.CubeCoord(1),
       ",", c0.CubeCoord(2), ").");
    throw error;
  }

  /// Get orthogonal direction and coordinate of facet shared 
  ///   by cubes containing isov0 and isov1.<br>
  /// Get also orientation from cube containing isov0 
  ///   to cube containing isov1.
  /// @pre Some facet separates cube containing isov0 from isov1.
  template <typename GRID_CUBE_TYPE, typename DUAL_ISOV_TYPE,
            typename VTYPE0, typename VTYPE1,
            typename DTYPE, typename CTYPE, typename ORIENT_TYPE>
  void get_shared_facet_coord_and_orient_of_cubes_containing_isov_3D
  (const std::vector<GRID_CUBE_TYPE> & cube_isov_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const VTYPE0 isov0, const VTYPE1 isov1,
   DTYPE & facet_orth_dir,
   CTYPE & orth_dir_coord,
   ORIENT_TYPE & orient)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE index_to_c0 = isov_list[isov0].cube_list_index;
    const SIZE_TYPE index_to_c1 = isov_list[isov1].cube_list_index;
    
    get_shared_facet_coord_and_orient_3D
      (cube_isov_list[index_to_c0], cube_isov_list[index_to_c1],
       facet_orth_dir, orth_dir_coord, orient);
  }


  /// Get cube c0 index of facet shared by c0 and c1
  /// @param[out] facet_index Index in c0 of facet shared by c0 and c1.
  ///             In range [0..5].
  template <typename GRID_CUBE_TYPE0, typename GRID_CUBE_TYPE1,
            typename ITYPE>
  void get_shared_facet_index_3D
  (const GRID_CUBE_TYPE0 & c0, const GRID_CUBE_TYPE1 & c1,
   ITYPE & c0_facet_index)
  {
    typedef typename GRID_CUBE_TYPE0::CUBE_COORD_TYPE CUBE_COORD_TYPE0;
    typedef typename GRID_CUBE_TYPE1::CUBE_COORD_TYPE CUBE_COORD_TYPE1;

    const ITYPE DIM3(3);

    for (ITYPE d = 0; d < DIM3; d++) {
      const CUBE_COORD_TYPE0 coord0 = c0.CubeCoord(d);
      const CUBE_COORD_TYPE1 coord1 = c1.CubeCoord(d);
      if (coord0 < coord1) {
        c0_facet_index = d+DIM3;
        return;
      }
      else if (coord0 > coord1) {
        c0_facet_index = d;
        return;
      }
    }

    IJK::PROCEDURE_ERROR error("get_shared_facet_index");
    error.AddMessage
      ("Programming error.  Cubes c0 and c1 have identical coordinates.");
    error.AddMessage
      ("  Coordinates: (", c0.CubeCoord(0), ",", c0.CubeCoord(1),
       ",", c0.CubeCoord(2), ").");
    throw error;
  }


  /// Get cube c0 index of facet shared by cubes containing isov0 and isov1.
  /// @param[out] facet_index Index in c0 of facet shared by c0 and c1.
  ///             In range [0..5].
  template <typename GRID_CUBE_TYPE, typename DUAL_ISOV_TYPE,
            typename VTYPE0, typename VTYPE1,
            typename ITYPE>
  void get_shared_facet_index_of_cubes_containing_isov_3D
  (const std::vector<GRID_CUBE_TYPE> & cube_isov_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const VTYPE0 isov0, const VTYPE1 isov1,
   ITYPE & c0_facet_index)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE index_to_c0 = isov_list[isov0].cube_list_index;
    const SIZE_TYPE index_to_c1 = isov_list[isov1].cube_list_index;
    
    get_shared_facet_index_3D
      (cube_isov_list[index_to_c0], cube_isov_list[index_to_c1],
       c0_facet_index);
  }

  ///@}


  // **************************************************
  /// @name QUERY LOCATION ROUTINES
  // **************************************************

  ///@{

  /// Get distance from isov to axis-orthogonal plane.
  /// @param orth_dir Direction orthogonal to plane.
  /// @param c Coordinate of plane in direction orth_dir.
  template <typename CTYPE0, typename CTYPE1,
            typename VTYPE, typename DTYPE>
  CTYPE0 get_distance_to_orth_plane_3D
  (const std::vector<CTYPE0> & vertex_coord,
   const VTYPE isov,
   const DTYPE orth_dir,
   const CTYPE1 c)
  {
    const DTYPE DIM3(3);
    CTYPE0 diff;

    const CTYPE0 isov_coord = vertex_coord[isov*DIM3 + orth_dir];

    if (isov_coord < c) { diff = c - isov_coord; }
    else { diff = isov_coord - c; }

    return(diff);
  }


  /// Return true if isosurface vertex is at most max_distance 
  ///   to given plane.
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename VTYPE, typename DTYPE>
  bool is_close_3D
  (const std::vector<CTYPE0> & vertex_coord,
   const VTYPE isov,
   const DTYPE orth_dir,
   const CTYPE1 c,
   const CTYPE2 max_distance)
  {
    const CTYPE0 dist = 
      get_distance_to_orth_plane_3D(vertex_coord, isov, orth_dir, c);

    if (dist <= max_distance) { return(true); }
    else { return(false); }

  }


  /// Return true if both isosurface vertices are at most max_distance 
  ///   to given plane.
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename VTYPE0, typename VTYPE1, typename DTYPE>
  bool are_both_close_3D
  (const std::vector<CTYPE0> & vertex_coord,
   const VTYPE0 isov0,
   const VTYPE1 isov1,
   const DTYPE orth_dir,
   const CTYPE1 c,
   const CTYPE2 max_distance)
  {
    const bool flag = 
      (is_close_3D(vertex_coord, isov0, orth_dir, c, max_distance) &&
       is_close_3D(vertex_coord, isov1, orth_dir, c, max_distance));

    return(flag);
  }

  /// Return true if either isosurface vertex is at most max_distance 
  ///   to given plane.
  template <typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename VTYPE0, typename VTYPE1, typename DTYPE>
  bool is_either_close_3D
  (const std::vector<CTYPE0> & vertex_coord,
   const VTYPE0 isov0,
   const VTYPE1 isov1,
   const DTYPE orth_dir,
   const CTYPE1 c,
   const CTYPE2 max_distance)
  {
    const bool flag = 
      (is_close_3D(vertex_coord, isov0, orth_dir, c, max_distance) ||
       is_close_3D(vertex_coord, isov1, orth_dir, c, max_distance));

    return(flag);
  }

  ///@}


  // **************************************************
  /// @name QUERY CUBE ROUTINES
  // **************************************************

  ///@{

  /// Return coordinates of cube containing isosurface vertex isov.
  template <typename GRID_CUBE_TYPE, typename DUAL_ISOV_TYPE,
            typename VTYPE>
  const typename GRID_CUBE_TYPE::CUBE_COORD_TYPE * 
  get_coord_of_cube_containing_isov
  (const std::vector<GRID_CUBE_TYPE> & cube_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const VTYPE isov)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE cube_list_index = isov_list[isov].cube_list_index;
    return(cube_list[cube_list_index].CubeCoord());
  }

  /// Get shared facet lower/left vertex.
  /// @pre icube0 and icube1 are indices of distinct cubes
  ///   and share a grid vertex.
  template <typename VTYPE0, typename VTYPE1>
  VTYPE0 get_shared_facet_LL_vertex
  (const VTYPE0 icube0, const VTYPE1 icube1)
  {
    if (icube0 < icube1) { return(icube1); }
    else { return(icube0); }
  }

  /// Get shared facet lower/left vertex.
  template <typename GRID_CUBE_TYPE>
  typename GRID_CUBE_TYPE::CUBE_INDEX_TYPE get_shared_facet_LL_vertex
  (const GRID_CUBE_TYPE & c0, const GRID_CUBE_TYPE & c1)
  {
    return(get_shared_facet_LL_vertex(c0.cube_index, c1.cube_index));
  }

  /// Get shared facet lower/left vertex of cubes containing isov0 and isov1.
  template <typename DUAL_ISOV_TYPE, 
            typename VTYPE0, typename VTYPE1>
  typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE 
  get_shared_facet_LL_vertex_of_cubes_containing_isov
  (const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const VTYPE0 isov0,
   const VTYPE1 isov1)
  {
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;

    const CUBE_INDEX_TYPE icube0 = isov_list[isov0].cube_index;
    const CUBE_INDEX_TYPE icube1 = isov_list[isov1].cube_index;

    return(get_shared_facet_LL_vertex(icube0, icube1));
  }

  /// Get lower-left coordinate of facet shared by cubes c0 and c1.
  /// @pre Cubes c0 and c1 are distinct and share a facet.
  template <typename GRID_CUBE_TYPE>
  const typename GRID_CUBE_TYPE::CUBE_COORD_TYPE * 
  get_shared_facet_LL_coord
  (const GRID_CUBE_TYPE & c0, const GRID_CUBE_TYPE & c1)
  {
    if (c0.cube_index < c1.cube_index) { return(c1.CubeCoord()); }
    else { return(c0.CubeCoord()); }
  }

  /// Get lower-left coordinate of facet shared by cubes 
  ///   containing isov0 and isov1.
  template <typename GRID_CUBE_TYPE, typename DUAL_ISOV_TYPE, 
            typename VTYPE0, typename VTYPE1>
  const typename GRID_CUBE_TYPE::CUBE_COORD_TYPE * 
  get_shared_facet_LL_coord_of_cubes_containing_isov
  (const std::vector<GRID_CUBE_TYPE> & cube_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const VTYPE0 isov0, const VTYPE1 isov1)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE index_to_c0 = isov_list[isov0].cube_list_index;
    const SIZE_TYPE index_to_c1 = isov_list[isov1].cube_list_index;
    
    return(get_shared_facet_LL_coord
           (cube_list[index_to_c0], cube_list[index_to_c1]));
  }


  /// Copy lower-left coordinate of facet shared by cubes 
  ///   containing isov0 and isov1 into facet_coord
  template <typename GRID_CUBE_TYPE, typename DUAL_ISOV_TYPE, 
            typename VTYPE0, typename VTYPE1, typename CTYPE>
  void copy_shared_facet_LL_coord_of_cubes_containing_isov_3D
  (const std::vector<GRID_CUBE_TYPE> & cube_list,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const VTYPE0 isov0, const VTYPE1 isov1,
   CTYPE facet_LL_coord[])
  {
    typedef typename GRID_CUBE_TYPE::CUBE_COORD_TYPE CUBE_COORD_TYPE;

    const int DIM3(3);
    const CUBE_COORD_TYPE * coord =
      get_shared_facet_LL_coord_of_cubes_containing_isov
      (cube_list, isov_list, isov0, isov1);

    std::copy(coord, coord+DIM3, facet_LL_coord);
  }


  ///@}

}

#endif
