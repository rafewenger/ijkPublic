/*!
 *  \file ijkdual_position.tpp
 *  @brief Position dual isosurface vertices.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2011-2024 Rephael Wenger

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

#ifndef IJKDUAL_POSITION_TPP_
#define IJKDUAL_POSITION_TPP_

#include <bitset>
#include <vector>

#ifdef INCLUDE_RANDOM
#include "ijkrandom.tpp"
#endif

#include "ijkdual_datastruct.h"

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkcube.tpp"
#include "ijkinterpolate.tpp"
#include "ijkisocoord.tpp"
#include "ijkisopoly.tpp"
#include "ijkscalar_grid.tpp"
#include "ijksort.tpp"


namespace IJKDUAL {

  // ******************************************************************
  //! @name Position single isosurface vertex in each grid cube.
  // ******************************************************************

  ///@{

  /*!
   *  @brief Position dual isosurface vertices in cube centers.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube active_cube_list[i].
   *    @pre Array coord[] is preallocated to size at least
   *    active_cube_list.size() * grid.Dimension().
   */
  template <typename GRID_TYPE, typename ITYPE,
            typename CTYPE>
  void position_all_dual_isovertices_cube_center
  (const GRID_TYPE & grid,
   const std::vector<ITYPE> & active_cube_list, CTYPE * coord)
  {
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;

    const int dimension = grid.Dimension();

    for (SIZE_TYPE i = 0; i < active_cube_list.size(); i++) {
      const ITYPE icube = active_cube_list[i];

      grid.ComputeCubeCenterCoord(icube, coord+i*dimension);
    }
  }


  /// @brief Position dual isosurface vertices in cube centers.
  /// - C++ STL vector format for array coord[].
  template <typename GRID_TYPE, typename ITYPE,
            typename CTYPE>
  void position_all_dual_isovertices_cube_center
  (const GRID_TYPE & grid,
   const std::vector<ITYPE> & active_cube_list,
   std::vector<CTYPE> & coord)
  {
    const int dimension = grid.Dimension();

    coord.resize(active_cube_list.size()*dimension);
    position_all_dual_isovertices_cube_center
      (grid, active_cube_list, IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertex in centroid
   *  of isosurface-edge intersections. (Simple version.)
   *  - Simple version that does not precompute intersections
   *    of isosurface and grid edges.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube icube
   *    @pre Array coord[] is preallocated to size at least
   *      grid.Dimension().
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename CTYPEV, typename CTYPET>
  void position_dual_isov_centroid_simple
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const ITYPE icube,
   CTYPEV * vcoord,
   CTYPET * temp_coord0, CTYPET * temp_coord1, CTYPET * temp_coord2)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const int dimension = scalar_grid.Dimension();

    NTYPE num_intersected_edges = 0;
    IJK::set_coord(dimension, 0.0, vcoord);

    for (int edge_dir = 0; edge_dir < dimension; edge_dir++)
      for (NTYPE k = 0; k < scalar_grid.NumCubeFacetVertices(); k++) {
        const VERTEX_INDEX_TYPE iend0 =
          scalar_grid.FacetVertex(icube, edge_dir, k);
        const VERTEX_INDEX_TYPE iend1 =
          scalar_grid.NextVertex(iend0, edge_dir);

        const STYPE s0 = scalar_grid.Scalar(iend0);
        bool is_end0_positive = true;
        if (s0 < isovalue)
          { is_end0_positive = false; };

        const STYPE s1 = scalar_grid.Scalar(iend1);
        bool is_end1_positive = true;
        if (s1 < isovalue)
          { is_end1_positive = false; };

        if (is_end0_positive != is_end1_positive) {

          scalar_grid.ComputeCoord(iend0, temp_coord0);
          scalar_grid.ComputeCoord(iend1, temp_coord1);

          IJK::linear_interpolate_coord
            (dimension, s0, temp_coord0, s1, temp_coord1,
             isovalue, temp_coord2);

          IJK::add_coord
            (dimension, vcoord, temp_coord2, vcoord);

          num_intersected_edges++;
        }

      }

    if (num_intersected_edges > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_intersected_edges, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /*!
   *  @brief Position all dual isosurface vertices in centroid
   *  of isosurface-edge intersections. (Simple version.)
   *  - Simple version that does not precompute intersections
   *    of isosurface and grid edges.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename CTYPE>
  void position_all_dual_isovertices_centroid_simple
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   CTYPE * coord)
  {
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> temp_coord0(dimension);
    IJK::ARRAY<CTYPE> temp_coord1(dimension);
    IJK::ARRAY<CTYPE> temp_coord2(dimension);

    for (SIZE_TYPE i = 0; i < active_cube_list.size(); i++) {
      const ITYPE icube = active_cube_list[i];

      position_dual_isov_centroid_simple
        (scalar_grid, isovalue, icube, coord+i*dimension,
         temp_coord0.Ptr(), temp_coord1.Ptr(), temp_coord2.Ptr());
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *    of isosurface-edge intersections. (Simple version. C++ vector.)
   *  - Simple version that does not precompute intersections
   *    of isosurface and grid edges.
   *  - C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename CTYPE>
  void position_all_dual_isovertices_centroid_simple
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(active_cube_list.size()*dimension);
    position_all_dual_isovertices_centroid_simple
      (scalar_grid, isovalue, active_cube_list,
       IJK::vector2pointerNC(coord));
  }


  /*!
   *  @brief Position dual isosurface vertices in centroid
   *  of isosurface-edge intersections.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube active_cube_list[i].
   *    @pre Array coord[] is preallocated to size at least
   *    active_cube_list.size() * grid.Dimension().
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename ISOV_INDEX_TYPE, typename ISOPOLY_INFO_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_centroid
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   CTYPE * isov_coord)
  {
    typedef typename std::vector<ITYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const int dimension = scalar_grid.Dimension();
    const int num_vertices_per_isopoly =
      scalar_grid.NumFacetVertices();
    const NTYPE num_isov = active_cube_list.size();
    const NTYPE num_isopoly = isopoly_info.size();
    IJK::ARRAY<CTYPE> temp_coord0(dimension);
    IJK::ARRAY<CTYPE> temp_coord1(dimension);
    IJK::ARRAY<CTYPE> temp_coord2(dimension);
    IJK::PROCEDURE_ERROR
      error("position_all_dual_isovertices_centroid");

    if (!IJK::check_array_allocated
        (isov_coord, "isov_coord", error)) { throw error; }

    if (!check_isopoly_vert_size
        (isopoly_vert, num_vertices_per_isopoly, isopoly_info.size(),
         error)) { throw error; }

    // num_incident_isopoly[i] =
    //   Number of isosurface polytopes incident on isov_list[i].
    std::vector<ISO_VERTEX_DEGREE>
      num_incident_isopoly(num_isov, 0);

    // Set all isov_coord entries to 0.0.
    for (NTYPE i = 0; i < num_isov*dimension; i++)
      { isov_coord[i] = 0.0; }

    SIZE_TYPE k = 0;
    for (NTYPE ipoly = 0; ipoly < num_isopoly; ipoly++) {

      compute_isosurface_grid_edge_intersection_linear
        (scalar_grid, isovalue, isopoly_info[ipoly],
         temp_coord0.Ptr(), temp_coord1.Ptr(), temp_coord2.Ptr());

      for (int j = 0; j < num_vertices_per_isopoly; j++) {
        const ISOV_INDEX_TYPE kw = isopoly_vert[k];
        CTYPE * isov_coord_ptr = isov_coord + dimension*kw;

        IJK::add_coord
          (dimension, isov_coord_ptr, temp_coord0.Ptr(), isov_coord_ptr);
        num_incident_isopoly[kw]++;
        k++;
      }
    }

    if (k != isopoly_vert.size()) {
      error.AddMessage
        ("Programming error. Processed incorrect number of isopoly vertices.");
      error.AddMessage
        ("  Isosurface polytopes are incident on ",
         isopoly_vert.size(), " vertices.");
      error.AddMessage("  Processed ", k, " vertices.");
      throw error;
    }

    for (NTYPE kw = 0; kw < num_isov; kw++) {
      const ITYPE icube = active_cube_list[kw];
      CTYPE * isov_coord_ptr = isov_coord + dimension*kw;

      if (scalar_grid.IsCubeOnGridBoundary(icube)) {
        // Isosurface vertex is in cube on grid boundary.
        // Some dual facets may be missing from isopoly_info.
        // Use slower method to compute isosurface vertex position.
        position_dual_isov_centroid_simple
          (scalar_grid, isovalue, icube, isov_coord_ptr,
           temp_coord0.Ptr(), temp_coord1.Ptr(), temp_coord2.Ptr());
      }
      else {
        IJK::multiply_coord
          (dimension, 1.0/num_incident_isopoly[kw],
           isov_coord_ptr, isov_coord_ptr);
      }
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices in centroid
   *  of isosurface-edge intersections. (C++ vector.)
   *  - C++ STL vector format for array coord[].
   *  @param[out] coord[] Isosurface vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate
   *      of isosurface vertex in cube active_cube_list[i].
   */
  template <typename GRID_TYPE, typename STYPE, typename ITYPE,
            typename ISOV_INDEX_TYPE, typename ISOPOLY_INFO_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_centroid
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(active_cube_list.size()*dimension);
    position_all_dual_isovertices_centroid
      (scalar_grid, isovalue, active_cube_list,
       isopoly_vert, isopoly_info, IJK::vector2pointerNC(coord));
  }


  ///@}


  // ******************************************************************
  //! @name Position multiple isosurface vertices in each grid cube.
  // ******************************************************************

  ///@{

  /*!
   *  @brief Position dual isosurface vertex at centroid.
   *    (Simple version.)
   *  - Simple version that does not precompute intersections
   *    of isosurface and grid edges.
   *  @pre scalar_grid.Dimension() == cube.Dimension().
   *  @param[out] vcoord[d] d'th coordinate of isosurface vertex.
   *    @pre Array vcoord[] is preallocated to size at least
   *      scalar_grid.Dimension().
   *  @param temp_coord0[] Temporary coordinate array.
   *    @pre Pre-allocated to size at least dimension.
   *  @param temp_coord1[] Temporary coordinate array.
   *    @pre Pre-allocated to size at least dimension.
   *  @param temp_coord2[] Temporary coordinate array.
   *    @pre Pre-allocated to size at least dimension.
   */
  template <typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE,
            typename STYPE,
            typename DUAL_ISOV_TYPE,
            typename CUBE_TYPE,
            typename CTYPE, typename CTYPE2>
  void position_dual_isov_centroid_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const DUAL_ISOV_TYPE & isov_info,
   const CUBE_TYPE & cube,
   CTYPE * vcoord,
   CTYPE2 * temp_coord0, CTYPE2 * temp_coord1, CTYPE2 * temp_coord2)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = scalar_grid.Dimension();
    const VTYPE icube = isov_info.cube_index;
    const NTYPE ipatch = isov_info.patch_index;
    const TABLE_INDEX table_index = isov_info.table_index;
    const NTYPE num_incident_isopoly =
      isodual_table.NumIncidentIsoPoly(table_index, ipatch);

    IJK::set_coord(dimension, 0.0, vcoord);

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
      if (isodual_table.IsBipolar(table_index, ie)) {
        if (isodual_table.IncidentIsoVertex(table_index, ie) == ipatch) {
          const NTYPE k0 = cube.EdgeEndpoint(ie, 0);
          const NTYPE k1 = cube.EdgeEndpoint(ie, 1);
          const VTYPE iend0 = scalar_grid.CubeVertex(icube, k0);
          const VTYPE iend1 = scalar_grid.CubeVertex(icube, k1);

          compute_isosurface_line_segment_intersection_linear
            (scalar_grid, isovalue, iend0, iend1,
             temp_coord0, temp_coord1, temp_coord2);

          IJK::add_coord(dimension, vcoord, temp_coord0, vcoord);
        }
      }
    }

    if (num_incident_isopoly > 0) {
      IJK::multiply_coord
        (dimension, 1.0/num_incident_isopoly, vcoord, vcoord);
    }
    else {
      scalar_grid.ComputeCubeCenterCoord(icube, vcoord);
    }
  }


  /*!
   *  @brief Position all dual isosurface vertices using centroids.
   *    (Simple version.)
   */
  template <typename GRID_TYPE,
            typename ISODUAL_TABLE_TYPE,
            typename STYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_centroid_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> coord0(dimension);
    IJK::ARRAY<CTYPE> coord1(dimension);
    IJK::ARRAY<CTYPE> coord2(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

    for (SIZE_TYPE i = 0; i < isov_list.size(); i++) {
      position_dual_isov_centroid_multi_simple
        (scalar_grid, isodual_table, isovalue, isov_list[i], cube,
         coord+i*dimension, coord0.Ptr(), coord1.Ptr(), coord2.Ptr());
    }
  }


  /*
   *  @overload
   *  @brief Position all dual isosurface vertices using centroids.
   *    (Simple version. C++ vector.)
   *  - Version using C++ STL vector for coord.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_centroid_multi_simple
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_dual_isovertices_centroid_multi_simple
      (scalar_grid, isodual_table,isovalue, isov_list,
       IJK::vector2pointerNC(coord));

  }


  /*!
   *  @brief Position all dual isosurface vertices using centroids.
   *  @param isov_list[] List of isosurface vertices
   *    including containing cube and associated isosurface patch.
   *  @param isopoly_vert[i*numv_per_isopoly+k]
   *    Index of k'th vertex of isosurface polytope i.
   *  @param isopoly_info[] Isososurface polygon information,
   *    including dual grid edge.
   *  @param[out] coord[] Array of vertex coordinates.
   *    - coord[i*dimension+d] = d'th coordinate of isosurface vertex i.
   *    @pre coord[] is pre-allocated to size at least
   *      isov_list.size()*scalar_grid.Dimension().
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_centroid_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;

    const int dimension = scalar_grid.Dimension();
    const int num_vertices_per_isopoly =
      scalar_grid.NumFacetVertices();
    IJK::ARRAY<CTYPE> temp_coord0(dimension);
    IJK::ARRAY<CTYPE> temp_coord1(dimension);
    IJK::ARRAY<CTYPE> temp_coord2(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    IJK::PROCEDURE_ERROR
      error("position_all_dual_isovertices_centroid_multi");

    if (!check_isopoly_vert_size
        (isopoly_vert, num_vertices_per_isopoly, isopoly_info.size(),
         error)) { throw error; }

    // num_incident_isopoly[i] =
    //   Number of isosurface polytopes incident on isov_list[i].
    std::vector<ISO_VERTEX_DEGREE>
      num_incident_isopoly(isov_list.size(), 0);

    // Set all coord entries to 0.0.
    for (SIZE_TYPE i = 0; i < isov_list.size()*dimension; i++)
      { coord[i] = 0.0; }

    SIZE_TYPE k = 0;
    for (SIZE_TYPE ipoly = 0; ipoly < isopoly_info.size(); ipoly++) {

      compute_isosurface_grid_edge_intersection_linear
        (scalar_grid, isovalue, isopoly_info[ipoly],
         temp_coord0.Ptr(), temp_coord1.Ptr(), temp_coord2.Ptr());

      for (int j = 0; j < num_vertices_per_isopoly; j++) {
        const ISOV_INDEX_TYPE kw = isopoly_vert[k];
        CTYPE * coord_ptr = coord + dimension*kw;

        IJK::add_coord
          (dimension, coord_ptr, temp_coord0.Ptr(), coord_ptr);
        num_incident_isopoly[kw]++;
        k++;
      }
    }

    if (k != isopoly_vert.size()) {
      error.AddMessage
        ("Programming error. Processed incorrect number of isopoly vertices.");
      error.AddMessage
        ("  Isosurface polytopes are incident on ",
         isopoly_vert.size(), " vertices.");
      error.AddMessage("  Processed ", k, " vertices.");
      throw error;
    }

    for (SIZE_TYPE kw = 0; kw < isov_list.size(); kw++) {
      const int ipatch = isov_list[kw].patch_index;
      const TABLE_INDEX table_index = isov_list[kw].table_index;
      const ISO_VERTEX_DEGREE num_incident_isopolyB =
        isodual_table.NumIncidentIsoPoly(table_index, ipatch);

      CTYPE * coord_ptr = coord + dimension*kw;

      if (num_incident_isopoly[kw] == num_incident_isopolyB) {
        IJK::multiply_coord
          (dimension, 1.0/num_incident_isopoly[kw], coord_ptr, coord_ptr);
      }
      else {
        // Isosurface vertex is in grid cube with some bipolar edges
        //   on the grid boundary.
        // Dual facets are missing from isopoly_info.
        // Use slower method to compute isosurface vertex position.
        position_dual_isov_centroid_multi_simple
          (scalar_grid, isodual_table, isovalue, isov_list[kw], cube,
           coord_ptr, temp_coord0.Ptr(), temp_coord1.Ptr(), temp_coord2.Ptr());
      }
    }
  }


  /*
   *  @overload
   *  @brief Position all dual isosurface vertices using centroids.
   *    (C++ vector.)
   *  - Version using C++ STL vector for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_centroid_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_dual_isovertices_centroid_multi
      (scalar_grid, isodual_table,isovalue,
       isov_list, isopoly_vert, isopoly_info,
       IJK::vector2pointerNC(coord));

  }


  /*!
   *  @brief Position dual isosurface vertices near cube centers.
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    near but not on cube center.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE0, typename CTYPE1>
  void position_all_dual_isovertices_near_cube_center_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE0 offset,
   CTYPE1 * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int dimension = scalar_grid.Dimension();

    IJK::ARRAY<CTYPE1> vcoord(dimension);

    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);
    IJK::UNIT_CUBE<int,int,int> unit_cube(dimension);
    bool intersects_facet[2*dimension];

    for (SIZE_TYPE i = 0; i < isov_list.size(); i++) {
      const VTYPE icube = isov_list[i].cube_index;
      const NTYPE ipatch = isov_list[i].patch_index;
      const TABLE_INDEX it = isov_list[i].table_index;

      if (isodual_table.NumIsoVertices(it) == 1) {
        scalar_grid.ComputeCubeCenterCoord(icube, coord+i*dimension);
      }
      else {

        scalar_grid.ComputeCubeCenterCoord(icube, vcoord.Ptr());

        // Set intersects facet to false.
        for (int d = 0; d < dimension; d++) {
          intersects_facet[2*d] = false;
          intersects_facet[2*d+1] = false;
        }

        for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
          if (isodual_table.IsBipolar(it, ie)) {
            if (isodual_table.IncidentIsoVertex(it, ie) == ipatch) {
              NTYPE k0 = cube.EdgeEndpoint(ie, 0);
              NTYPE k1 = cube.EdgeEndpoint(ie, 1);

              for (int d = 0; d < dimension; d++) {
                NTYPE c = unit_cube.VertexCoord(k0, d);
                if (c == unit_cube.VertexCoord(k1, d)) {
                  if (c > 0)
                    { intersects_facet[2*d+1] = true; }
                  else
                    { intersects_facet[2*d] = true; }
                }
              }
            }
          }
        }

        for (int d = 0; d < dimension; d++) {

          if (intersects_facet[2*d] != intersects_facet[2*d+1]) {
            if (intersects_facet[2*d])
              { vcoord[d] -= offset; }
            else
              { vcoord[d] += offset; }
          }
        }
        IJK::copy_coord(dimension, vcoord.PtrConst(), coord+i*dimension);
      }
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices near cube centers.
   *    (C++ vector)
   *  - More than one vertex can be in a cube.
   *  - C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE0, typename CTYPE1>
  void position_all_dual_isovertices_near_cube_center_multi
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE0 offset,
   std::vector<CTYPE1> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_dual_isovertices_near_cube_center_multi
      (scalar_grid, isodual_table, isov_list, offset,
       IJK::vector2pointerNC(coord));
  }

  ///@}


  // ******************************************************************
  //! @name Position lifted interval volume vertices.
  // ******************************************************************

  ///@{

  /*!
   *  @brief Position interval volume vertices which have been lifted
   *    to one higher dimension.
   *  @param isovalue0 Lower isovalue.
   *  @param isovalue1 Upper isovalue.
   *  @pre isovalue0 < isovalue1.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISOVAL0_TYPE, typename ISOVAL1_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_ivol_lifted
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISOVAL0_TYPE isovalue0,
   const ISOVAL1_TYPE isovalue1,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const int dimension = scalar_grid.Dimension();

    if (dimension < 1) { return; }

    const VTYPE numv_in_grid_facet_maxd =
      scalar_grid.AxisIncrement(dimension-1);
    IJK::ARRAY<CTYPE> coord0(dimension);
    IJK::ARRAY<CTYPE> coord1(dimension);
    IJK::ARRAY<CTYPE> coord2(dimension);
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {

      const VTYPE icube = isov_list[isov].cube_index;

      if (icube < numv_in_grid_facet_maxd) {
        // Apply isovalue1 to (*,*,*,0) grid cubes.
        position_dual_isov_centroid_multi_simple
          (scalar_grid, isodual_table, isovalue1, isov_list[isov], cube,
           coord+isov*dimension, coord0.Ptr(), coord1.Ptr(), coord2.Ptr());
      }
      else {
        // Apply isovalue0 to (*,*,*,1) grid cubes.
        position_dual_isov_centroid_multi_simple
          (scalar_grid, isodual_table, isovalue0, isov_list[isov], cube,
           coord+isov*dimension, coord0.Ptr(), coord1.Ptr(), coord2.Ptr());

      }
    }
  }


  /*!
   *  @overload
   *  @brief Position interval volume vertices which have been lifted
   *    to one higher dimension. (C++ vector.)
   *  - Version using C++ STL vector for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE0, typename STYPE1,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_ivol_lifted
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE0 isovalue0,
   const STYPE1 isovalue1,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   std::vector<CTYPE> & coord)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = scalar_grid.Dimension();

    coord.resize(isov_list.size()*dimension);
    position_all_dual_isovertices_ivol_lifted
      (scalar_grid, isodual_table, isovalue0, isovalue1, isov_list,
       IJK::vector2pointerNC(coord));
  }

  ///@}


  // ******************************************************************
  //! @name Reposition isosurface vertices away from cube facets.
  // ******************************************************************


  /*!
   *  @brief Reposition isosurface vertices away from cube facets.
   *  - Note: Forces isosurface vertex to lie in interior
   *    of associated cube.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE,
            typename ISOV_INDEX_TYPE, typename CUBE_INDEX_TYPE,
            typename OFFSET_TYPE, typename CTYPE,
            typename NTYPE>
  void reposition_dual_isovertex_away_from_cube_facets
  (const GRID_TYPE & grid,
   const ISOV_INDEX_TYPE isov,
   const CUBE_INDEX_TYPE icube,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPE & num_repositions)
  {
    const int dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    // Initialize
    num_repositions = 0;

    grid.ComputeCoord(icube, cube_coord.Ptr());
    CTYPE * first_isov_coord = coord + dimension*isov;

    for (int d = 0; d < dimension; d++) {
      const CTYPE minc = cube_coord[d] + offset;
      const CTYPE maxc = cube_coord[d] + 1 - offset;
      if (first_isov_coord[d] < minc) {
        first_isov_coord[d] = minc;
        num_repositions++;
      }
      else if (first_isov_coord[d] > maxc) {
        first_isov_coord[d] = maxc;
        num_repositions++;
      }
    }
  }


  /*!
   *  @brief Reposition all isosurface vertices away from cube facets.
   *  - Note: Forces each isosurface vertices to lie in associated cube.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename NTYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPER>
  void reposition_all_dual_isovert_away_from_cube_facets
  (const GRID_TYPE & grid,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPE num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPER & num_repositions)
  {
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;

    num_repositions = 0;

    for (NTYPE isov = 0; isov < num_isov; isov++) {
      const CUBE_INDEX_TYPE icube = isov_list[isov].cube_index;
      NTYPER num_repositions_i;
      reposition_dual_isovertex_away_from_cube_facets
        (grid, isov, icube, offset, coord, num_repositions_i);
      num_repositions += num_repositions_i;
    }
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube facets.
   *    (Specialization of isov_list[] to int.)
   *  - Specialization for isov_list[] being a list of cube indices.
   *  - Note: Forces each isosurface vertices to lie in associated cube.
   *  @param isov_list[isov] Index of cube containing isosurface vertex isov.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE, typename NTYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPER>
  void reposition_all_dual_isovert_away_from_cube_facets
  (const GRID_TYPE & grid,
   const int isov_list[],
   const NTYPE num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPER & num_repositions)
  {
    num_repositions = 0;

    for (NTYPE isov = 0; isov < num_isov; isov++) {
      NTYPER num_repositions_i;
      reposition_dual_isovertex_away_from_cube_facets
        (grid, isov, isov_list[isov], offset, coord, num_repositions_i);
      num_repositions += num_repositions_i;
    }
  }


  /*!
   *  @brief Reposition all isosurface vertices away from cube facets.
   *    (C++ STL vector for isov_list[].)
   *  - Version with C++ STL vector format for array isov_list[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void reposition_all_dual_isovert_away_from_cube_facets
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPE & num_repositions)
  {
    reposition_all_dual_isovert_away_from_cube_facets
      (grid, IJK::vector2pointer(isov_list), isov_list.size(),
       offset, coord, num_repositions);
  }


  /*!
   *  @brief Reposition all isosurface vertices away from cube facets.
   *    (C++ STL vector for isov_list[] and coord[].)
   *  - Version with C++ STL vector format for array isov_list[].
   *  - Version with C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void reposition_all_dual_isovert_away_from_cube_facets
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord,
   NTYPE & num_repositions)
  {
    reposition_all_dual_isovert_away_from_cube_facets
      (grid, isov_list, offset, IJK::vector2pointerNC(coord),
       num_repositions);
  }


  // ******************************************************************
  //! @name Reposition isosurface vertices away from cube ridges.
  // ******************************************************************

  /*!
   *  @brief Reposition isosurface vertices away from cube ridges.
   *  - Note: Forces isosurface vertex to lie in interior
   *    of associated cube.
   *  - In 3D, cube ridges are cube edges.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE,
            typename ISOV_INDEX_TYPE, typename CUBE_INDEX_TYPE,
            typename OFFSET_TYPE, typename CTYPE,
            typename NTYPE>
  void reposition_dual_isovertex_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const ISOV_INDEX_TYPE isov,
   const CUBE_INDEX_TYPE icube,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPE & num_repositions)
  {
    const int dimension = grid.Dimension();
    const int TWO(2);
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    // Initialize
    num_repositions = 0;

    grid.ComputeCoord(icube, cube_coord.Ptr());
    CTYPE * first_isov_coord = coord + dimension*isov;

    int num_close_facets = 0;
    
    // Count number of facets close to isosurface vertex.
    for (int d = 0; d < dimension; d++) {
      const CTYPE minc = cube_coord[d] + offset;
      const CTYPE maxc = cube_coord[d] + 1 - offset;
      if (first_isov_coord[d] < minc)
        { num_close_facets++; }
      else if (first_isov_coord[d] > maxc) {
        { num_close_facets++; }
      }
    }

    if (num_close_facets >= TWO) {
      reposition_dual_isovertex_away_from_cube_facets
        (grid, isov, icube, offset, coord, num_repositions);
    }
  }


  /*!
   *  @brief Reposition all isosurface vertices away from cube ridges.
   *  - Note: Forces each isosurface vertices to lie in associated cube.
   *  @param offset Offset. Must be in range [0,0.5].
   *  - In 3D, cube ridges are cube edges.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename NTYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPER>
  void reposition_all_dual_isovert_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPE num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPER & num_repositions)
  {
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX_TYPE;

    num_repositions = 0;

    for (NTYPE isov = 0; isov < num_isov; isov++) {
      const CUBE_INDEX_TYPE icube = isov_list[isov].cube_index;
      NTYPER num_repositions_i;
      reposition_dual_isovertex_away_from_cube_ridges
        (grid, isov, icube, offset, coord, num_repositions_i);
      num_repositions += num_repositions_i;
    }
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube ridges.
   *    (Specialization of isov_list[] to int.)
   *  - Specialization for isov_list[] being a list of cube indices.
   *  - Note: Forces each isosurface vertices to lie in associated cube.
   *  - In 3D, cube ridges are cube edges.
   *  @param isov_list[isov] Index of cube containing isosurface vertex isov.
   *  @param offset Offset. Must be in range [0,0.5].
   */
  template <typename GRID_TYPE, typename NTYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPER>
  void reposition_all_dual_isovert_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const int isov_list[],
   const NTYPE num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPER & num_repositions)
  {
    num_repositions = 0;

    for (NTYPE isov = 0; isov < num_isov; isov++) {
      NTYPER num_repositions_i;
      reposition_dual_isovertex_away_from_cube_ridges
        (grid, isov, isov_list[isov], offset, coord, num_repositions_i);
      num_repositions += num_repositions_i;
    }
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube ridges.
   *    (C++ STL vector for isov_list[].)
   *  - Version with C++ STL vector format for array isov_list[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void reposition_all_dual_isovert_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord,
   NTYPE & num_repositions)
  {
    reposition_all_dual_isovert_away_from_cube_ridges
      (grid, IJK::vector2pointer(isov_list), isov_list.size(),
       offset, coord, num_repositions);
  }


  /*!
   *  @overload
   *  @brief Reposition all isosurface vertices away from cube ridges.
   *    (C++ STL vector for isov_list[] and coord[].)
   *  - Version with C++ STL vector format for array isov_list[].
   *  - Version with C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE, typename CTYPE, typename NTYPE>
  void reposition_all_dual_isovert_away_from_cube_ridges
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord,
   NTYPE & num_repositions)
  {
    reposition_all_dual_isovert_away_from_cube_ridges
      (grid, isov_list, offset, IJK::vector2pointerNC(coord),
       num_repositions);
  }
  

  // ******************************************************************
  //! @name Separate isosurface vertices that are near grid facets.
  // ******************************************************************

  /*!
   *  @brief Store information about isosurface vertex, grid facet pair.
   *  - Facet is store implicitly.
   */
  template <typename ISOV_TYPE, typename DIST_TYPE>
  class FACET_ISOV_PAIR {
  public:

    /*!
     *  @brief isov Isosurface vertex.
     *  - Note: isov_list[isov] contains index of cube containing isov.
     *  - Axis orthogonal to facet is fixed so there are 
     *    only two possible facets orthogonal to given axis.
     *  - Boolean flag_below determines which facet in cube
     *    in in this facet-isov pair.
     */
    ISOV_TYPE isov;

    /// @brief Distance from isov to facet
    DIST_TYPE distance;

    /// @brief True if facet contains lowest/leftmost cube vertex
    ///   of cube containing isov.
    bool flag_below;

    /// True if isosurface vertex was repositioned.
    bool flag_repositioned;

    void Set
    (const ISOV_TYPE isov, const bool flag_below,
     const DIST_TYPE distance)
    {
      this->isov = isov;
      this->flag_below = flag_below;
      this->distance = distance;
    }

    /// Constructor.
    FACET_ISOV_PAIR
    (const ISOV_TYPE isov, const bool flag_below,
     const DIST_TYPE distance)
    {
      flag_repositioned = false;
      Set(isov, flag_below, distance);
    }

  };


  /*!
   *  @brief Store information about distance from facets
   *    incident on a grid vertex to closest isosurface vertices.
   */
  template <typename DIST_TYPE>
  class FACET2ISOV_DISTANCE {

  public:

    /*!
     *  @brief Distance from facets incident on a grid vertex
     *    to closest isosurface vertex.
     *  - Note: Grid vertex and direction orthogonal to facets
     *    are implicit.
     */
    DIST_TYPE distance;

    /// @brief True if facet contains lowest/leftmost cube vertex
    ///   of cube containing isov.
    bool flag_below;

    void Set(const bool flag_below, const DIST_TYPE distance)
    {
      this->flag_below = flag_below;
      this->distance = distance;
    }

  };



  namespace {

    /*!
     *  @brief Add facet-isov pair to facet_isov_list[].
     *  - Add facet vertices to facet_vertex[].
     *  - Add locations in facet_isov_list[] to facet_isov_loc[].
     *  @param icube Cube containing isosurface vertex.
     *  @param ifacet Cube facet. Index in range [0,(num_cube_facets-1)].
     *  @param isov Isosurface vertex in cube \a icube.
     *  @param distance Distance from isov to facet.
     *  @param[out] facet_isov_list[] List of pairs of facets
     *    and isosurface vertices where the isosurface vertex
     *    is near the facet.
     *  @param[out] facet_vertex[] List of vertices of each facet
     *    in facet_isov_list[].
     *  @param[out] facet_isov_loc[i] Location in facet_isov_list[]
     *    of facet containing facet_vertex[i].
     */
    template <typename GRID_TYPE, typename ITYPEC, typename ITYPEF,
              typename ISOV_TYPE, typename DIST_TYPE,
              typename PAIR_TYPE, typename VTYPE, typename LTYPE>
    void add_facet_isov_pair
    (const GRID_TYPE & grid,
     const ITYPEC icube,
     const ITYPEF ifacet,
     const ISOV_TYPE isov,
     const DIST_TYPE distance,
     std::vector<PAIR_TYPE> & facet_isov_list,
     std::vector<VTYPE> & facet_vertex,
     std::vector<LTYPE> & facet_isov_loc)
    {
      const bool flag_below = grid.IsLowerFacet(ifacet);
      const LTYPE iloc = facet_isov_list.size();
      const PAIR_TYPE facet_isov_pair(isov, flag_below, distance);

      facet_isov_list.push_back(facet_isov_pair);

      for (int j = 0; j < grid.NumCubeFacetVertices(); j++) {
        const ITYPEC jv = grid.FacetVertex(icube, ifacet, j);
        facet_vertex.push_back(jv);
        facet_isov_loc.push_back(iloc);
      }
    }

  }


  /*!
   *  @brief Separate isosurface vertices that are near shared grid cube facets.
   *  - Separate isosurface vertices w0, w1 that are from two different
   *    grid cubes that intersect at a facet.
   *  - Move vertex that is furthest from the facet.
   *  @param min_distance Move vertices if both cubes have
   *    isosurface vertices within distance (Linf) \a min_distance
   *    of facet f.
   *  @param offset Move vertices so that they are distance
   *    \a offset from the facet f.
   *  @param min_num_bins Minimum number of bins used for sorting.
   *    - If list range is less than min_num_bins, use couting_sort().
   *  @param[out] num_repositions Number of times isosurface vertices
   *    are repositioned away from some grid cube facet.
   *    - Note: An isosurface vertex that is repositioned away
   *      from multiple grid cube facets will contribute
   *      multiple times to this count.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE,
            typename NTYPEV, typename NTYPEB, typename NTYPER>
  void separate_isov_near_shared_grid_cube_facets
  (const GRID_TYPE & grid,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPEV num_isov,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   const NTYPEB min_num_bins,
   CTYPE * coord,
   NTYPER & num_repositions)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE CUBE_INDEX_TYPE;
    typedef FACET_ISOV_PAIR<NUMBER_TYPE,CTYPE> FACET_ISOV_PAIR_TYPE;
    typedef FACET2ISOV_DISTANCE<CTYPE> FACET2ISOV_DISTANCE_TYPE;
    typedef typename std::vector<FACET_ISOV_PAIR_TYPE>::size_type
      SIZE_TYPE;

    const int dimension = grid.Dimension();
    IJK::PROCEDURE_ERROR error
      ("separate_isov_near_grid_cube_facets");

    // Initialize
    num_repositions = 0;

    // List of pairs of facets and isosurface vertices.
    std::vector<FACET_ISOV_PAIR_TYPE> facet_isov_list;

    // List of lower/leftmost vertices of facets in facet_isov_list[].
    std::vector<VERTEX_INDEX> facet_vertex;

    for (int d = 0; d < dimension; d++ ) {
      facet_isov_list.clear();
      facet_vertex.clear();
      for (NTYPEV isov = 0; isov < num_isov; isov++) {
        const CUBE_INDEX_TYPE icube =
          get_isovert_grid_cube_index(isov_list[isov]);

        const CTYPE isov_coord_d = coord[dimension*isov+d];
        const CTYPE cube_coord_d = grid.CoordD(icube, d);
        const CTYPE cdiff0 = isov_coord_d - cube_coord_d;
        const CTYPE cdiff1 = (cube_coord_d+1) - isov_coord_d;

        if (cdiff0 <= min_distance) {
          const FACET_ISOV_PAIR_TYPE facet_isov_pair(isov, true, cdiff0);
          facet_isov_list.push_back(facet_isov_pair);

          // is_primary = lower/leftmost vertex of facet.
          const VERTEX_INDEX_TYPE iv_primary =
            get_isovert_grid_cube_index(isov_list[isov]);
          facet_vertex.push_back(iv_primary);
        }

        if (cdiff1 <= min_distance) {
          const FACET_ISOV_PAIR_TYPE facet_isov_pair(isov, false, cdiff1);
          facet_isov_list.push_back(facet_isov_pair);
          const VERTEX_INDEX_TYPE icube =
            get_isovert_grid_cube_index(isov_list[isov]);

          // is_primary = lower/leftmost vertex of facet.
          const VERTEX_INDEX_TYPE iv_primary = grid.NextVertex(icube, d);
          facet_vertex.push_back(iv_primary);
        }
        
      }

      // Merge identical vertices in facet_vertex[].
      // Since the orthogonal direction is fixed, each vertex
      //   represents a unique facet.
      std::vector<VERTEX_INDEX_TYPE> facet_vertex_nodup;
      std::vector<SIZE_TYPE> facet_vertex_loc;
      IJKSORT::merge_identical_radix_sort_iter2
       (facet_vertex, 0, grid.NumVertices(), min_num_bins,
        facet_vertex_nodup, facet_vertex_loc);

      // Determine closest distance to each facet.
      std::vector<FACET2ISOV_DISTANCE_TYPE>
        isov2facet_distance(facet_vertex_nodup.size());

      // Initialize each entry in isov2facet_distance[].
      for (SIZE_TYPE i = 0; i < facet_vertex.size(); i++) {
        // May overwrite isov2facet_distance[iloc] multiple times,
        //   but final result will be some correct distance.
        const SIZE_TYPE iloc = facet_vertex_loc[i];
        const bool flag_below = facet_isov_list[i].flag_below;
        const CTYPE distance = facet_isov_list[i].distance;
        isov2facet_distance[iloc].Set(flag_below, distance);
      }

      // Determine minimum distance to each facet.
      for (SIZE_TYPE i = 0; i < facet_vertex.size(); i++) {
        const VERTEX_INDEX_TYPE iloc = facet_vertex_loc[i];
        const bool flag_below = facet_isov_list[i].flag_below;
        const CTYPE distance = facet_isov_list[i].distance;
        if (distance < isov2facet_distance[iloc].distance)
          { isov2facet_distance[iloc].Set(flag_below, distance); }
      }

      // Move vertices away from facets.
      for (SIZE_TYPE i = 0; i < facet_vertex.size(); i++) {
        const VERTEX_INDEX_TYPE iloc = facet_vertex_loc[i];

        const bool flag_below = facet_isov_list[i].flag_below;

        if (flag_below != isov2facet_distance[iloc].flag_below) {
          const NUMBER_TYPE isov = facet_isov_list[i].isov;
          const NUMBER_TYPE iv = facet_vertex[i];
          const CTYPE facet_coord_d = grid.CoordD(iv, d);
          facet_isov_list[i].flag_repositioned = true;

          // Move isov away from facet.
          if (flag_below) {
            // Facet is below isov.
            coord[dimension*isov+d] = facet_coord_d + offset;
          }
          else {
            // Facet is above isov.
            coord[dimension*isov+d] = facet_coord_d - offset;
          }

          num_repositions++;
        }
      }
    }
  }


  /*!
   *  @overload
   *  @brief Separate isosurface vertices that are near shared grid cube facets.
   *    (C++ vector.)
   *  - Version using C++ STL vector for isov_list[] and coord[].
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename MIN_DIST_TYPE, typename OFFSET_TYPE,
            typename CTYPE, typename NTYPEB, typename NTYPER>
  void separate_isov_near_shared_grid_cube_facets
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const MIN_DIST_TYPE min_distance,
   const OFFSET_TYPE offset,
   const NTYPEB min_num_bins,
   std::vector<CTYPE> & coord,
   NTYPER & num_repositioned)
  {
    separate_isov_near_shared_grid_cube_facets
    (grid, IJK::vector2pointer(isov_list), isov_list.size(),
     min_distance, offset, min_num_bins,
     IJK::vector2pointerNC(coord), num_repositioned);
  }
  

  // ******************************************************************
  //! @name Reposition "doubly connected" isosurface vertices.
  // ******************************************************************

  ///@{

  /*!
   *  @brief Return true if isov0 is "doubly connected" to the cube containing isov1.
   *  @pre isov0 and isov1 are in cubes sharing a facet.
   *  @pre isov0 is connected to isov1.
   *  @param ishared_facet Index of facet of icube0 shared with icube1.
   */
  template <typename ISODUAL_TABLE_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename FACET_INDEX_TYPE, typename CUBE_TYPE,
            typename CTYPE>
  bool is_doubly_connected
  (const ISODUAL_TABLE_TYPE & isodual_table,
   const ISOV_INDEX_TYPE isov0,
   const ISOV_INDEX_TYPE isov1,
   const DUAL_ISOV_TYPE isov_list[],
   const FACET_INDEX_TYPE ishared_facet,
   const CUBE_TYPE & cube,
   CTYPE * coord)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int dimension = cube.Dimension();

    const TABLE_INDEX it0 = isov_list[isov0].table_index;
    const TABLE_INDEX it1 = isov_list[isov1].table_index;

    if (isodual_table.NumIsoVertices(it1) <= 1) {
      // Cube containing isov1 has only one vertex.
      // isov0 cannot be connected to two vertices
      //   in cube containing isov1.
      return false;
    }

    if (!isodual_table.IsFacetAmbiguous(it0, ishared_facet)) {
      // Shared facet is not ambiguous.
      // isov0 is only connected to isov1 in cube containing isov1.
      return false;
    }

    const int ipatch0 = isov_list[isov0].patch_index;
    const int ipatch1 = isov_list[isov1].patch_index;
    const int facet_orth_dir = cube.FacetOrthDir(ishared_facet);

    // For each edge in shared facet.
    for (int d = 0; d < dimension; d++) {
      if (d == facet_orth_dir) { continue; }

      for (NUMBER_TYPE j = 0; j < cube.NumFacetVertices(); j++) {

        const NUMBER_TYPE je0 = d*cube.NumFacetVertices() + j;
        const NUMBER_TYPE jend0 = cube.EdgeEndpoint(je0, 0);

        if (cube.IsFacetIncidentOnVertex(ishared_facet, jend0)) {
          // Vertex jend0 and edge je0 is in facet ishared_facet.

          if (isodual_table.IsBipolar(it0, je0)) {
            if (isodual_table.IncidentIsoVertex(it0, je0) == ipatch0) {

              const NUMBER_TYPE je1 = cube.OppositeEdge(je0, facet_orth_dir);

              if (isodual_table.IncidentIsoVertex(it1, je1) != ipatch1)
                { return true; }
            }
          }
        }
      }
    }

    // All isosurface edges dual to ishared_facet and incident on isov0
    //   are also incident on isov1.
    // isov0 is not doubly connected to isosurface vertices in cube 1.
    return false;
  }


  /*!
   *  @brief Reposition isosurface vertex if is "doubly connected"
   *    to the cube containing isov1.
   *  @pre isov0 and isov1 are in cubes sharing a facet.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename CUBE_TYPE, typename CTYPE>
  void reposition_doubly_connected_dual_isovertex_cube
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const ISOV_INDEX_TYPE isov0,
   const ISOV_INDEX_TYPE isov1,
   const DUAL_ISOV_TYPE isov_list[],
   const CUBE_TYPE & cube,
   CTYPE * coord)
  {
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX;

    const int dimension = grid.Dimension();
    const CUBE_INDEX icube0 = isov_list[isov0].cube_index;
    const CUBE_INDEX icube1 = isov_list[isov1].cube_index;

    const int ishared_facet = grid.SharedFacet(icube0, icube1);

    if (is_doubly_connected
        (isodual_table, isov0, isov1,
         isov_list, ishared_facet, cube, coord)) {

      const TABLE_INDEX it1 = isov_list[isov1].table_index;
      const int ipatch1 = isov_list[isov1].patch_index;
      const int facet_orth_dir = cube.FacetOrthDir(ishared_facet);
      CTYPE * first_isov0_coord = coord + dimension*isov0;
      CTYPE * first_isov1_coord = coord + dimension*isov1;

      // Note: This does not correctly handle all cases in 4D.
      for (int d = 0; d < dimension; d++) {

        if (d == facet_orth_dir) {
          // Skip coordinate corresponding to facet_orth_dir.
          continue;
        }

        // At most one of these two flags should be true.
        bool flag_lower_pos = false;
        bool flag_upper_pos = false;

        const bool flag_intersects_lower_facet =
          isodual_table_vinfo.VertexInfo(it1, ipatch1).IntersectsFacet(d);
        const bool flag_intersects_upper_facet =
          isodual_table_vinfo.VertexInfo(it1, ipatch1).IntersectsFacet(d+dimension);

        if (flag_intersects_lower_facet != flag_intersects_upper_facet) {
          flag_lower_pos = flag_intersects_lower_facet;
          flag_upper_pos = flag_intersects_upper_facet;
        }

        if (flag_lower_pos) {
          if (first_isov0_coord[d] < first_isov1_coord[d])
            { first_isov0_coord[d] = first_isov1_coord[d]; }
        }
        else if (flag_upper_pos) {
          if (first_isov0_coord[d] > first_isov1_coord[d])
            { first_isov0_coord[d] = first_isov1_coord[d]; }
        }

      }
    }
  }


  /*!
   *  @brief Reposition "doubly connected" isosurface vertices.
   *  - An isosurface vertex is doubly connected if it is connected
   *    to two isosurface vertices in the same cube.
   *  @param isopoly_vert Isosurface poly vertices are listed
   *    by increasing coordinate, not clockwise or counter-clockwise
   *    around polygon.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename NTYPE,
            typename DUAL_ISOV_TYPE,
            typename CUBE_TYPE, typename CTYPE>
  void reposition_all_doubly_connected_dual_isovertices_cube
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const ISOV_INDEX_TYPE isopoly_vert[],
   const NTYPE num_isopoly,
   const DUAL_ISOV_TYPE isov_list[],
   const CUBE_TYPE & cube,
   CTYPE * coord)
  {
    const int dimension = grid.Dimension();

    // *** REPLACE BY grid.NumFacetVertices();
    const NTYPE num_vert_per_isopoly = cube.NumFacetVertices();

    if (cube.Dimension() != dimension) {
      IJK::PROCEDURE_ERROR error
        ("reposition_all_doubly_connected_dual_isovertices_cube");

      error.AddMessage
        ("Programming error. Cube dimension does not equal grid dimension.");
      error.AddMessage("  Cube dimension: ", cube.Dimension(), "");
      error.AddMessage("  Grid dimension: ", dimension, "");
      throw error;
    }

    for (NTYPE ipoly = 0; ipoly < num_isopoly; ipoly++) {

      // first_vert: First vertex in isosurface polytope ipoly.
      const ISOV_INDEX_TYPE * first_isov =
        isopoly_vert + ipoly*num_vert_per_isopoly;

      // For each edge in facet 0.
      for (int d = 0; d+1 < dimension; d++) {

        for (NTYPE j = 0; j < num_vert_per_isopoly; j++) {

          const NTYPE je = (d*cube.NumFacetVertices()) + j;
          const NTYPE jend0 = cube.EdgeEndpoint(je, 0);

          if (jend0 < num_vert_per_isopoly) {
            // Vertex jend0 and edge je are on facet 0

            const NTYPE jend1 = cube.EdgeEndpoint(je, 1);
            const ISOV_INDEX_TYPE isov0 = first_isov[jend0];
            const ISOV_INDEX_TYPE isov1 = first_isov[jend1];

            reposition_doubly_connected_dual_isovertex_cube
              (grid, isodual_table, isodual_table_vinfo,
               isov0, isov1, isov_list, cube, coord);

            reposition_doubly_connected_dual_isovertex_cube
              (grid, isodual_table, isodual_table_vinfo,
               isov1, isov0, isov_list, cube, coord);
          }
        }
      }
    }

  }


  /*!
   *  @brief Reposition "doubly connected" isosurface vertices.
   *  - An isosurface vertex is doubly connected if it is connected
   *    to two isosurface vertices in the same cube.
   *  - C++ STL vector format for arrays isopoly_vert[] and isov_list[].
   *  @param isopoly_vert Isosurface poly vertices are listed
   *    by increasing coordinate, not clockwise or counter-clockwise
   *    around polygon.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
	    typename ISODUAL_TABLE_VINFO_TYPE,
	    typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
	    typename CUBE_TYPE, typename CTYPE>
  void reposition_all_doubly_connected_dual_isovertices_cube
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CUBE_TYPE & cube,
   CTYPE * coord)
  {
    typedef typename std::vector<ISOV_INDEX_TYPE>::size_type SIZE_TYPE;

    // *** REPLACE BY grid.NumFacetVertices().
    const SIZE_TYPE num_vert_per_isopoly = cube.NumFacetVertices();
    const SIZE_TYPE num_isopoly = isopoly_vert.size()/num_vert_per_isopoly;

    reposition_all_doubly_connected_dual_isovertices_cube
      (grid, isodual_table, isodual_table_vinfo,
       IJK::vector2pointer(isopoly_vert), num_isopoly,
       IJK::vector2pointer(isov_list),
       cube, coord);
  }

  //@}


  // ******************************************************************
  //! @name Separate isosurface vertices by cube centers.
  // ******************************************************************

  ///@{

  /*!
   *  @brief Separate isosurface vertex by cube center from other isosurface vertices.
   *  - Reposition isosurface vertices in cube containing multiple isosurface vertices.
   *  - If isosurface vertices isov0 and isov1 are in the same cube
   *    and each has an incident edge dual to the same cube facet f,
   *    then reposition isov0, if necessary.
   *  - Reposition isov0 so that isov_coord[d] <= cube_coord[d]+0.5-offset
   *    or isov_coord[d] >= cube_coord[d]+0.5+offset, as appropriate.
   *  @param offset Offset. Must be in range [0,0.25].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename OFFSET_TYPE,
            typename CTYPE>
  void separate_dual_isov_by_cube_center
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const ISOV_INDEX_TYPE isov0,
   const DUAL_ISOV_TYPE isov_list[],
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef typename DUAL_ISOV_TYPE::CUBE_INDEX_TYPE CUBE_INDEX;
    typedef typename DUAL_ISOV_TYPE::PATCH_INDEX_TYPE PATCH_INDEX;

    const int dimension = grid.Dimension();
    const CUBE_INDEX icube0 = isov_list[isov0].cube_index;
    const TABLE_INDEX it0 = isov_list[isov0].table_index;
    const PATCH_INDEX ipatch0 = isov_list[isov0].patch_index;
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    if (isodual_table.NumIsoVertices(it0) < 2) {
      // No repositioning necessary.
      return;
    }

    grid.ComputeCoord(icube0, cube_coord.Ptr());
    CTYPE * first_isov0_coord = coord + dimension*isov0;

    for (PATCH_INDEX ipatch1 = 0;
        ipatch1 < isodual_table.NumIsoVertices(it0); ipatch1++) {

      if (ipatch1 == ipatch0) { continue; }

      for (int d = 0; d < dimension; d++) {

        // At most one of these two flags should be true.
        bool flag_lower_pos = false;
        bool flag_upper_pos = false;

        const bool flag_intersects_lower_facet0 =
          isodual_table_vinfo.VertexInfo(it0, ipatch0).IntersectsFacet(d);
        const bool flag_intersects_upper_facet0 =
          isodual_table_vinfo.VertexInfo(it0, ipatch0).IntersectsFacet(d+dimension);
        const bool flag_intersects_lower_facet1 =
          isodual_table_vinfo.VertexInfo(it0, ipatch1).IntersectsFacet(d);
        const bool flag_intersects_upper_facet1 =
          isodual_table_vinfo.VertexInfo(it0, ipatch1).IntersectsFacet(d+dimension);

        if (flag_intersects_lower_facet0 &&
            flag_intersects_upper_facet1) {
          if (!flag_intersects_upper_facet0 ||
              !flag_intersects_lower_facet1) {
            flag_lower_pos = true;
          }
        }
        else if (flag_intersects_upper_facet0 &&
                 flag_intersects_lower_facet1) {
          if (!flag_intersects_lower_facet0 ||
              !flag_intersects_upper_facet1) {
            flag_upper_pos = true;
          }
        }

        if (flag_lower_pos) {
          const CTYPE cmax = cube_coord[d] + 0.5 - offset;
          if (first_isov0_coord[d] > cmax)
            { first_isov0_coord[d] = cmax; }
        }
        else if (flag_upper_pos) {
            const CTYPE cmax = cube_coord[d] + 0.5 + offset;
            if (first_isov0_coord[d] < cmax)
              { first_isov0_coord[d] = cmax; }
        }
      }

    }


  }


  /*!
   *  @brief Separate all isosurface vertices by cube centers from other isosurface vertices.
   *  - Reposition isosurface vertices in cubes containing multiple isosurface vertices.
   *  - If isosurface vertices isov0 and isov1 are in the same cube
   *    and each has an incident edge dual to the same cube facet f,
   *    then reposition isov0, if necessary.
   *  - Reposition isov0 so that isov_coord[d] <= cube_coord[d]+0.5-offset
   *    or isov_coord[d] >= cube_coord[d]+0.5+offset, as appropriate.
   *  - Calls reposition_all_doubly_connected_dual_isovertices_cube().
   *  @param offset Offset. Must be in range [0,0.25].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename NTYPEP,
            typename DUAL_ISOV_TYPE, typename NTYPEV,
            typename OFFSET_TYPE,
            typename CTYPE>
  void separate_all_dual_isov_by_cube_centers
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const ISOV_INDEX_TYPE isopoly_vert[],
   const NTYPEP num_isopoly,
   const DUAL_ISOV_TYPE isov_list[],
   const NTYPEV num_isov,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    for (NTYPEV isov = 0; isov < num_isov; isov++) {
      separate_dual_isov_by_cube_center
        (grid, isodual_table, isodual_table_vinfo,
         isov, isov_list, offset, coord);
    }

    IJK::CUBE_FACE_INFO<int,int,int> cube(grid.Dimension());
    reposition_all_doubly_connected_dual_isovertices_cube
      (grid, isodual_table, isodual_table_vinfo,
       isopoly_vert, num_isopoly, isov_list, cube, coord);
  }


  /*!
   *  @brief Separate all isosurface vertices by cube centers from other isosurface vertices.
   *  @brief Reposition all isosurface vertices in cubes containing multiple isosurface vertices.
   *  - Version with C++ STL vector format for array isov_list[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE,
            typename DUAL_ISOV_TYPE, typename OFFSET_TYPE,
            typename CTYPE>
  void separate_all_dual_isov_by_cube_centers
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   CTYPE * coord)
  {
    typedef typename std::vector<ISOV_INDEX_TYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE num_vert_per_isopoly = grid.NumCubeFacetVertices();
    const SIZE_TYPE num_isopoly = isopoly_vert.size()/num_vert_per_isopoly;

    separate_all_dual_isov_by_cube_centers
    (grid, isodual_table, isodual_table_vinfo,
     IJK::vector2pointer(isopoly_vert), num_isopoly,
     IJK::vector2pointer(isov_list), isov_list.size(),
     offset, coord);
  }


  /*!
   *  @brief Separate all isosurface vertices by cube centers from other isosurface vertices.
   *  - Version with C++ STL vector format for array isov_list[].
   *  - Version with C++ STL vector format for array coord[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE,
            typename DUAL_ISOV_TYPE, typename OFFSET_TYPE,
            typename CTYPE>
  void separate_all_dual_isov_by_cube_centers
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const OFFSET_TYPE offset,
   std::vector<CTYPE> & coord)
  {
    separate_all_dual_isov_by_cube_centers
    (grid, isodual_table, isodual_table_vinfo, isopoly_vert, isov_list,
     offset, IJK::vector2pointerNC(coord));
  }

  ///@}


  // ******************************************************************
  //! @name Compute location of isosurface vertex on grid edge.
  // ******************************************************************

  ///@{

  /*!
   *  @brief Compute position of isosurface vertex on grid edge.
   *  - Compute position using linear interpolation on scalar values
   *    at grid edge endpoints.
   *  - Store coefficient indicating grid edge position in isopoly_info[].
   *  @tparam PTYPE  Precision type to be used in calculations.
   *    - Should be float or double.
   *  @param max_small_difference Maximum small difference.
   *    - If abs(s1-s0) <= max_small_diffence, return 0.5.
   *  @pre max_small_difference >= 0.
   */
  template <typename PTYPE, typename GRID_TYPE, typename STYPE,
	    typename ISOPOLY_INFO_TYPE, typename NTYPE>
  void compute_grid_edge_isov_position_interpolate_scalar
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const ISOPOLY_INFO_TYPE isopoly_info[],
   const NTYPE num_poly,
   const PTYPE max_small_difference)
  {
    for (NTYPE ipoly = 0; ipoly < num_poly; ipoly++) {
      isopoly_info[ipoly].grid_edge_isov_position_coef =
        IJK::compute_linear_interpolation_coef_scalar_grid_edge
        (scalar_grid, isopoly_info[ipoly].DualGridEdge(),
         isovalue, max_small_difference);
    }
  }


  ///@}


#ifdef INCLUDE_RANDOM

  // **************************************************
  //! @name Position iso vertices at random locations.
  // **************************************************

  ///@{

   /*!
     *  @brief Position dual isosurface vertex at random location in cube.
     *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
     *  @param offset Position isosurface vertices at distance offset
     *      from cube facets.
     *  @param random_engine Pseudorandom number generator.
     *  @param random_pos Random number generator with a given distribution.
     *  @pre 0 <= offset <= 0.5.
     *  @cube_coord[] Coordinate of lower left vertex of cube.
     *  - Set by this routine.
     *  @pre cube_coord[] is allocated to size at list grid.Dimension().
     */
    template <typename GRID_TYPE, typename VTYPE,
              typename CUBE_INDEX_TYPE, typename CTYPE_OFFSET,
	      typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
	      typename CTYPE0, typename CTYPE1>
    void position_dual_isovert_random01
    (const GRID_TYPE & grid,
     const VTYPE isov,
     const CUBE_INDEX_TYPE icube,
     const CTYPE_OFFSET offset,
     RANDOM_ENGINE & random_engine,
     RANDOM_POS_TYPE & random_pos,
     CTYPE0 cube_coord[],
     CTYPE1 * coord)
    {
      const int dimension = grid.Dimension();
      CTYPE1 * first_isov_coord = coord + dimension*isov;

      grid.ComputeCoord(icube, cube_coord);

      // Put isosurface vertex in random position.
      for (int d = 0; d < dimension; d++) {
        first_isov_coord[d] =
	  cube_coord[d] + random_pos.Random01(offset, random_engine);
      }
    }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *  @pre 0 <= offset <= 0.5.
   *  @param random_pos Distribution is determined by class random_pos.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
	    typename CTYPE>
  void position_all_dual_isovertices_random
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    // @pre 0 <= offset <= 0.5, but clamp offset, just in case.
    const CTYPE_OFFSET offsetB =
      IJK::clamp_coord_to_range(offset, 0, 0.5);

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      const auto icube =
	IJK::get_isovert_grid_cube_index(isov_list[isov]);

      position_dual_isovert_random01
        (grid, isov, icube, offsetB, random_engine, random_pos,
	 cube_coord.Ptr(), coord);
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename CTYPE_OFFSET, typename RANDOM_ENGINE,
	    typename RANDOM_POS_TYPE, typename CTYPE>
  void position_all_dual_isovertices_random
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const CTYPE_OFFSET offset,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    const int dimension = grid.Dimension();

    coord.resize(vlist.size()*dimension);
    position_all_dual_isovertices_random
      (grid, vlist, offset, random_engine, random_pos,
       IJK::vector2pointerNC(coord));
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *  - Initialize random engine with seed.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE, typename CTYPE_OFFSET,
            typename RANDOM_SEED, typename RANDOM_POS_TYPE,
	    typename CTYPE>
  void position_all_dual_isovertices_random_init
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const CTYPE_OFFSET offset,
   const RANDOM_SEED random_seed,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    std::minstd_rand random_engine(random_seed);

    position_all_dual_isovertices_random
      (grid, vlist, offset, random_engine, random_pos, coord);
  }


  /*!
    *  @brief Position dual isosurface vertices at random location.
    *  - Use uniform distribution.
    *  - C++ STL vector format for array coord[].
    *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
    */
   template <typename GRID_TYPE, typename VTYPE,
             typename CTYPE_OFFSET, typename RANDOM_SEED,
	     typename CTYPE>
   void position_all_dual_isovertices_random_uniform
   (const GRID_TYPE & grid,
    const std::vector<VTYPE> & vlist,
    const CTYPE_OFFSET offset,
    const RANDOM_SEED random_seed,
    std::vector<CTYPE> & coord)
   {
     IJK::GENERATE_UNIFORM_RANDOM random_pos;

     position_all_dual_isovertices_random_init
       (grid, vlist, offset, random_seed, random_pos, coord);
   }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  - Use U-quadratic distribution.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename CTYPE_OFFSET, typename RANDOM_SEED,
	    typename CTYPE>
  void position_all_dual_isovertices_random_U_quadratic
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const CTYPE_OFFSET offset,
   const RANDOM_SEED random_seed,
   std::vector<CTYPE> & coord)
  {
    IJK::GENERATE_U_QUADRATIC_RANDOM random_pos;

    position_all_dual_isovertices_random_init
      (grid, vlist, offset, random_seed, random_pos, coord);
  }


  /*!
   *  @brief Position dual isosurface vertex at random location in cube.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param bwidth Boundary width. Controls boundary sampling.
   *    - Generate random number in range [-boundary_width,boundary_width]
   *      and then clamp to [0,1].
   *  @param random_engine Pseudorandom number generator.
   *  @param random_pos Random number generator with a given distribution.
   *  @param cube_coord[] Coordinate of lower left vertex of cube.
   *    - Set by this routine.
   *  @pre cube_coord[] is allocated to size at list grid.Dimension().
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename CUBE_INDEX_TYPE, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE0, typename CTYPE1>
   void position_dual_isovert_random01_sample_boundary
   (const GRID_TYPE & grid,
    const VTYPE isov,
    const CUBE_INDEX_TYPE icube,
    const BWIDTH boundary_width,
    RANDOM_ENGINE & random_engine,
    RANDOM_POS_TYPE & random_pos,
    CTYPE0 cube_coord[],
    CTYPE1 * coord)
   {
     const int dimension = grid.Dimension();
     CTYPE1 * first_isov_coord = coord + dimension*isov;

     grid.ComputeCoord(icube, cube_coord);

     // Put isosurface vertex in random position.
     for (int d = 0; d < dimension; d++) {
       const CTYPE1 rpos01 = random_pos.RandomAndClamp
         (0, 1, boundary_width, random_engine);

       first_isov_coord[d] = cube_coord[d] + rpos01;
     }
   }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param bwidth Boundary width. Controls boundary sampling.
   *     - Generate random number in range [-boundary_width,boundary_width]
   *       and then clamp to [0,1].
   *  @param random_pos Distribution is determined by class random_pos.
   */
  template <typename GRID_TYPE, typename DUAL_ISOV_TYPE, typename BWIDTH,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_random_sample_boundary
  (const GRID_TYPE & grid,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      const auto icube =
        IJK::get_isovert_grid_cube_index(isov_list[isov]);

      position_dual_isovert_random01_sample_boundary
      (grid, isov, icube, boundary_width, random_engine, random_pos,
       cube_coord.Ptr(), coord);
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename BWIDTH, typename RANDOM_ENGINE,
            typename RANDOM_POS_TYPE, typename CTYPE>
  void position_all_dual_isovertices_random_sample_boundary
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const BWIDTH boundary_width,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    const int dimension = grid.Dimension();

    coord.resize(vlist.size()*dimension);
    position_all_dual_isovertices_random_sample_boundary
    (grid, vlist, boundary_width, random_engine, random_pos,
     IJK::vector2pointerNC(coord));
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices at random location.
   *  - Initialize random engine with seed.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename VTYPE, typename BWIDTH,
            typename RANDOM_SEED, typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_random_sample_boundary_init
  (const GRID_TYPE & grid,
   const std::vector<VTYPE> & vlist,
   const BWIDTH boundary_width,
   const RANDOM_SEED random_seed,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
  {
    std::minstd_rand random_engine(random_seed);

    position_all_dual_isovertices_random_sample_boundary
    (grid, vlist, boundary_width, random_engine, random_pos, coord);
  }


  /*!
   *  @brief Position dual isosurface vertices at random location
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   *  @param isodual_table_vinfo Isodual table vertex information.
   *    @pre compute_dual_cube_isotable_vertex_connectivity()
   *       has been called on isodual_table_vinfo.
   *  @param offset Position isosurface vertices at distance offset
   *      from cube facets.
   *    @pre 0 <= offset <= 0.25.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET,
            typename CUBE_TYPE,
            typename RANDOM_ENGINE,
            typename RANDOM_POS_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_multi_random
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   const CUBE_TYPE & cube,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   CTYPE * coord)
  {
    typedef typename std::vector<DUAL_ISOV_TYPE>::size_type SIZE_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ISODUAL_TABLE_TYPE::TABLE_INDEX TABLE_INDEX;

    const int dimension = grid.Dimension();
    IJK::ARRAY<CTYPE> cube_coord(dimension);
    IJK::PROCEDURE_ERROR error
      ("position_all_dual_isovertices_multi_random");

    if (offset < 0.0 || offset > 0.25) {
      error.AddMessage("Programming error. Offset not in range [0,0.25].");
      throw error;
    }

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      const VTYPE icube = isov_list[isov].cube_index;
      const NTYPE ipatch = isov_list[isov].patch_index;
      const TABLE_INDEX it = isov_list[isov].table_index;

      CTYPE * first_isov_coord = coord + dimension*isov;

      grid.ComputeCoord(icube, cube_coord.Ptr());

      if (isodual_table.NumIsoVertices(it) == 1) {

        // Put isosurface vertex in random position.
        for (int d = 0; d < dimension; d++) {
          first_isov_coord[d] =
            cube_coord[d] + random_pos.Random01(offset, random_engine);
        }

      }
      else {

        NTYPE num_orth_dir_with_two_intersected_facets = 0;
        for (int d = 0; d < dimension; d++) {
          if (isodual_table_vinfo.VertexInfo(it,ipatch).IntersectsBothFacetsOrthogonalTo
              (d, dimension))
            { num_orth_dir_with_two_intersected_facets++; }
        }

        for (int d = 0; d < dimension; d++) {

          // At most one of these two flags should be true.
          bool flag_lower_pos = false;
          bool flag_upper_pos = false;

          // Note: This does not correctly handle all cases in 4D.
          if (num_orth_dir_with_two_intersected_facets >= dimension-1) {

            const NTYPE ipatch1 = (ipatch + 1)%(isodual_table.NumIsoVertices(it));

            // Use ipatch1 to determine position of vertex in ipatch0.
            // Note: This does not correctly handle all cases in 4D.
            // In particular, there may be 3 iso vertices,
            //   and the vertex should be positioned between two others.

            const bool flag_intersects_lower_facet =
              isodual_table_vinfo.VertexInfo(it, ipatch1).IntersectsFacet(d);
            const bool flag_intersects_upper_facet =
              isodual_table_vinfo.VertexInfo(it, ipatch1).IntersectsFacet(d+dimension);

            if (flag_intersects_lower_facet != flag_intersects_upper_facet) {
              flag_lower_pos = !flag_intersects_lower_facet;
              flag_upper_pos = !flag_intersects_upper_facet;
            }
          }
          else {
            const bool flag_intersects_lower_facet =
              isodual_table_vinfo.VertexInfo(it, ipatch).IntersectsFacet(d);
            const bool flag_intersects_upper_facet =
              isodual_table_vinfo.VertexInfo(it, ipatch).IntersectsFacet(d+dimension);

            if (flag_intersects_lower_facet != flag_intersects_upper_facet) {
              flag_lower_pos = flag_intersects_lower_facet;
              flag_upper_pos = flag_intersects_upper_facet;
            }

          }

          if (flag_lower_pos) {
            first_isov_coord[d] =
              cube_coord[d] +
              random_pos.Random(0.0, 0.5, offset, random_engine);
          }
          else if (flag_upper_pos) {
            first_isov_coord[d] =
              cube_coord[d] +
              random_pos.Random(0.5, 1.0, offset, random_engine);
          }
          else {
            first_isov_coord[d] =
              cube_coord[d] + random_pos.Random01(offset, random_engine);
          }

        }
      }
    }

    reposition_all_doubly_connected_dual_isovertices_cube
      (grid, isodual_table, isodual_table_vinfo,
       isopoly_vert, isov_list, cube, coord);
  }


  /*!
   *  @brief Position dual isosurface vertices at random location
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    in separate halves/quarters/eighths of cube.
   *  - C++ STL vector format for array coord[].
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
 	    typename ISODUAL_TABLE_VINFO_TYPE,
 	    typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
	    typename CTYPE_OFFSET, typename CUBE_TYPE,
	    typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
	    typename CTYPE>
  void position_all_dual_isovertices_multi_random
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const CTYPE_OFFSET offset,
   const CUBE_TYPE & cube,
   RANDOM_ENGINE & random_engine,
   RANDOM_POS_TYPE & random_pos,
   std::vector<CTYPE> & coord)
   {
     typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

     const DTYPE dimension = grid.Dimension();

     coord.resize(isov_list.size()*dimension);
     position_all_dual_isovertices_multi_random
       (grid, isodual_table, isodual_table_vinfo,
        isopoly_vert, isov_list, offset, cube,
	random_engine, random_pos,
        IJK::vector2pointerNC(coord));
   }


   /*!
    *  @brief Position dual isosurface vertices at random location
    *  - More than one vertex can be in a cube.
    *  - If cube contains multiple isosurface then vertices are positioned
    *    in separate halves/quarters/eighths of cube.
    *  - C++ STL vector format for array coord[].
    *  - Version that creates cube.
    *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
    */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE_OFFSET,
            typename RANDOM_ENGINE, typename RANDOM_POS_TYPE,
            typename CTYPE>
   void position_all_dual_isovertices_multi_random
   (const GRID_TYPE & grid,
    const ISODUAL_TABLE_TYPE & isodual_table,
    const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
    const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
    const std::vector<DUAL_ISOV_TYPE> & isov_list,
    const CTYPE_OFFSET offset,
    RANDOM_ENGINE & random_engine,
    RANDOM_POS_TYPE & random_pos,
    std::vector<CTYPE> & coord)
  {
    IJK::CUBE_FACE_INFO<int,int,int> cube(grid.Dimension());

    position_all_dual_isovertices_multi_random
    (grid, isodual_table, isodual_table_vinfo,
     isopoly_vert, isov_list, offset, cube,
     random_engine, random_pos, coord);
  }


  /*!
    *  @brief Position dual isosurface vertices at random location
    *  - More than one vertex can be in a cube.
    *  - If cube contains multiple isosurface then vertices are positioned
    *    in separate halves/quarters/eighths of cube.
    *  - C++ STL vector format for array coord[].
    *  - Version that creates cube.
    *  - Version that creates random_engine initialized with random_seed.
    *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
    */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
	    typename ISODUAL_TABLE_VINFO_TYPE,
  	    typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
 	    typename SEED_TYPE, typename CTYPE_OFFSET,
 	    typename RANDOM_POS_TYPE, typename CTYPE>
   void position_all_dual_isovertices_multi_random_init
   (const GRID_TYPE & grid,
    const ISODUAL_TABLE_TYPE & isodual_table,
    const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
    const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
    const std::vector<DUAL_ISOV_TYPE> & isov_list,
    const CTYPE_OFFSET offset,
    const SEED_TYPE random_seed,
    RANDOM_POS_TYPE & random_pos,
    std::vector<CTYPE> & coord)
  {
    std::minstd_rand random_engine(random_seed);

    position_all_dual_isovertices_multi_random
    (grid, isodual_table, isodual_table_vinfo,
     isopoly_vert, isov_list, offset,
     random_engine, random_pos, coord);
  }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  - Uniform distribution.
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    in separate halves/quarters/eighths of cube.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_multi_random_uniform
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const RANDOM_POS_PARAM & random_pos_param,
   std::vector<CTYPE> & coord)
  {
    const RANDOM_SEED_TYPE random_seed = random_pos_param.random_seed;
    IJK::GENERATE_UNIFORM_RANDOM random_pos;

    COORD_TYPE position_offset = 0.0;
    if (random_pos_param.flag_gen_isov_apply_offset)
      { position_offset = random_pos_param.position_offset; }

    if (random_pos_param.flag_gen_isov_on_cube_boundary) {
      position_all_dual_isovertices_random_sample_boundary_init
        (grid, isov_list, random_pos_param.boundary_width, random_seed,
         random_pos, coord);
    }
    else if (random_pos_param.flag_gen_isov_separated_by_cube_center) {
      position_all_dual_isovertices_multi_random_init
        (grid, isodual_table, isodual_table_vinfo,
         isopoly_vert, isov_list,
         position_offset, random_seed, random_pos, coord);
    }
    else {
      position_all_dual_isovertices_random_uniform
        (grid, isov_list, position_offset, random_seed, coord);
    }
  }


  /*!
    *  @brief Position dual isosurface vertices at random location
    *  - More than one vertex can be in a cube.
    *  - Use U quadratic distribution to position vertices.
    *  - If cube contains multiple isosurface then vertices are positioned ????
    *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
    *  @param isodual_table_vinfo Isodual table vertex information.
    *    @pre compute_dual_cube_isotable_vertex_connectivity()
    *       has been called on isodual_table_vinfo.
    */
   template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
             typename ISODUAL_TABLE_VINFO_TYPE,
             typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
             typename CTYPE>
   void position_all_dual_isovertices_multi_random_U_quadratic
   (const GRID_TYPE & grid,
    const ISODUAL_TABLE_TYPE & isodual_table,
    const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
    const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
    const std::vector<DUAL_ISOV_TYPE> & isov_list,
    const RANDOM_POS_PARAM & random_pos_param,
    std::vector<CTYPE> & coord)
   {
     const RANDOM_SEED_TYPE random_seed = random_pos_param.random_seed;
     IJK::GENERATE_U_QUADRATIC_RANDOM random_pos;

     COORD_TYPE position_offset = 0.0;
     if (random_pos_param.flag_gen_isov_apply_offset)
       { position_offset = random_pos_param.position_offset; }

     if (random_pos_param.flag_gen_isov_on_cube_boundary) {
       position_all_dual_isovertices_random_sample_boundary_init
         (grid, isov_list, random_pos_param.boundary_width, random_seed,
          random_pos, coord);
     }
     else {
       position_all_dual_isovertices_multi_random_init
         (grid, isodual_table, isodual_table_vinfo,
          isopoly_vert, isov_list,
          position_offset, random_seed, random_pos, coord);
     }
   }


  /*!
   *  @brief Position dual isosurface vertices at random location.
   *  - Uniform distribution.
   *  - More than one vertex can be in a cube.
   *  - If cube contains multiple isosurface then vertices are positioned
   *    in separate halves/quarters/eighths of cube.
   *  @pre Requires compilation with macro flag INCLUDE_RANDOM.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename ISODUAL_TABLE_VINFO_TYPE,
            typename ISOV_INDEX_TYPE, typename DUAL_ISOV_TYPE,
            typename CTYPE>
  void position_all_dual_isovertices_multi_random
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const ISODUAL_TABLE_VINFO_TYPE & isodual_table_vinfo,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const RANDOM_POS_PARAM & random_pos_param,
   std::vector<CTYPE> & coord)
  {
    const RANDOM_DISTRIBUTION random_pos_distribution =
      random_pos_param.distribution;
    IJK::PROCEDURE_ERROR
      error("position_all_dual_isovertices_multi_random");

    if (random_pos_distribution == UNIFORM_DISTRIBUTION) {
       position_all_dual_isovertices_multi_random_uniform
       (grid, isodual_table, isodual_table_vinfo,
        isopoly_vert, isov_list, random_pos_param,
        coord);
     }
     else if (random_pos_distribution == U_QUADRATIC_DISTRIBUTION) {
       position_all_dual_isovertices_multi_random_U_quadratic
       (grid, isodual_table, isodual_table_vinfo,
        isopoly_vert, isov_list, random_pos_param,
        coord);
      }
     else {
       error.AddMessage
         ("Programming error. Random position distribution not set.");
       throw error;
     }

  }


  ///@}

#endif

}

#endif
