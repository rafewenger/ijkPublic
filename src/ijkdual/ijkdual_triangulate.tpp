/*!
 *  @file ijkdual_triangulate.tpp
 *  @brief ijkdual triangulation routines.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2016-2024 Rephael Wenger

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


#ifndef _IJKDUAL_TRIANGULATE_TPP_
#define _IJKDUAL_TRIANGULATE_TPP_

#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"
#include "ijkisocoord.tpp"
#include "ijkmesh.tpp"
#include "ijktime.tpp"
#include "ijktri2D_max_min_angle.tpp"
#include "ijktri2D_envelope.tpp"
#include "ijktriangulate_poly2D.tpp"

// *** DEPRECATED ***
#include "ijktriangulate_geom.tpp"

#include "ijkdual_types.h"


namespace IJKDUAL {

  // ***********************************************************************
  //! @name COMPUTE INTERSECTION OF QUADRILATERAL AND GRID EDGE
  // ***********************************************************************

  //@{

  /*!
   *  @brief Compute intersection of bilinear surface patch and grid edge.
   *  - Bilinear surface patch is defined by four quad vertice.
   *  @pre Quad vertices "surround" grid edge, i.e. projection of grid edge
   *    in the grid edge direction is contained in the convex hull
   *    of the projection of the quad vertices.
   */
  template <typename CTYPE0, typename ISOV_INDEX_TYPE, 
            typename ISOPOLY_INFO_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1>
  void compute_intersection_of_bilinear_surface_and_grid_edge
  (const CTYPE0 endpoint0_coord[],
   const int edge_direction,
   const ISOV_INDEX_TYPE quad_vert[],
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 intersection_coord[])
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUAD(4);
    COORD_TYPE0 coord01[DIM3], coord23[DIM3];
    const COORD_TYPE0 * isov_coord[NUM_VERT_PER_QUAD];

    const int dir1 = (edge_direction+1)%DIM3;
    const int dir2 = (edge_direction+2)%DIM3;

    for (int i = 0; i < NUM_VERT_PER_QUAD; i++) {
      const ISOV_INDEX_TYPE isov = quad_vert[i];
      isov_coord[i] = &(coord_array[isov*DIM3]);
    }

    int dir01, dir12;
    if (IJK::is_dth_coord_in_range
        (dir1, endpoint0_coord, isov_coord[0], isov_coord[1]) &&
        IJK::is_dth_coord_in_range
        (dir1, endpoint0_coord, isov_coord[2], isov_coord[3])) {
      dir01 = dir1;
      dir12 = dir2;
    }
    else {
      dir01 = dir2;
      dir12 = dir1;
    }

    IJK::compute_line_segment_hyperplane_intersection
      (DIM3, coord_array, quad_vert[0], quad_vert[1],
       dir01, endpoint0_coord[dir01], coord01);
    IJK::compute_line_segment_hyperplane_intersection
      (DIM3, coord_array, quad_vert[2], quad_vert[3],
       dir01, endpoint0_coord[dir01], coord23);
    IJK::compute_line_segment_hyperplane_intersection
      (DIM3, coord01, coord23, dir12, endpoint0_coord[dir12], intersection_coord);
  }

  
  /*!
   *  @brief Compute intersection of bilinear surface patch and grid edge.
   *  - Bilinear surface patch is defined by four quad vertice.
   *  @pre Quad vertices "surround" grid edge, i.e. projection of grid edge
   *    in the grid edge direction is contained in the convex hull
   *    of the projection of the quad vertices.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename ISOPOLY_INFO_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1>
  void compute_intersection_of_bilinear_surface_and_grid_edge
  (const GRID_TYPE & grid, 
   const ISOV_INDEX_TYPE quad_vert[],
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const int DIM3(3);
    const VTYPE iv0 = isopoly_info.GridEdgeEndpoint0();
    const int edge_direction = isopoly_info.GridEdgeDirection();
    COORD_TYPE0 endpoint0_coord[DIM3];

    grid.ComputeCoord(iv0, endpoint0_coord);

    compute_intersection_of_bilinear_surface_and_grid_edge
      (endpoint0_coord, edge_direction, quad_vert, isopoly_info,
       coord_array, intersection_coord);
  }


  /*!
   *  @overload
   *  @brief Compute intersection of bilinear surface patch and grid edge.
   *  - Version with array of quad vertices.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename GRID_EDGE_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1, typename QUAD_INDEX_TYPE>
  void compute_intersection_of_bilinear_surface_and_grid_edge
  (const GRID_TYPE & grid, 
   const std::vector<ISOV_INDEX_TYPE> & quad_vert,
   const GRID_EDGE_TYPE & grid_edge,
   std::vector<COORD_TYPE0> & coord_array,
   const QUAD_INDEX_TYPE iquad,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE NUM_VERT_PER_QUAD(4);

    const ISOV_INDEX_TYPE * first_quad_vert = 
      &(quad_vert[iquad*NUM_VERT_PER_QUAD]);

    compute_intersection_of_bilinear_surface_and_grid_edge
      (grid, first_quad_vert, grid_edge, coord_array, intersection_coord);

  }


  /*!
   *  @brief Compute scaled intersection of bilinear surface patch and grid edge.
   *  - Scale vertex coordinates by grid.Spacing().
   *  @pre Coordinates in coord_array[] are already scaled by grid.Spacing().
   *  - Bilinear surface patch is defined by four quad vertice.
   *  @pre Quad vertices "surround" grid edge, i.e. projection of grid edge
   *    in the grid edge direction is contained in the convex hull
   *    of the projection of the quad vertices.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename ISOPOLY_INFO_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1>
  void compute_scaled_intersection_of_bilinear_surface_and_grid_edge
  (const GRID_TYPE & grid, 
   const ISOV_INDEX_TYPE quad_vert[],
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const int DIM3(3);
    const VTYPE iv0 = isopoly_info.GridEdgeEndpoint0();
    const int edge_direction = isopoly_info.GridEdgeDirection();
    COORD_TYPE0 endpoint0_coord[DIM3];

    grid.ComputeScaledCoord(iv0, endpoint0_coord);

    compute_intersection_of_bilinear_surface_and_grid_edge
      (endpoint0_coord, edge_direction, quad_vert, isopoly_info,
       coord_array, intersection_coord);
  }


  /*!
   *  @overload
   *  @brief Compute scaled intersection of bilinear surface patch and grid edge.
   *  - Scale vertex coordinates by grid.Spacing().
   *  - Version with array of quad vertices.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename GRID_EDGE_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1, typename QUAD_INDEX_TYPE>
  void compute_scaled_intersection_of_bilinear_surface_and_grid_edge
  (const GRID_TYPE & grid, 
   const std::vector<ISOV_INDEX_TYPE> & quad_vert,
   const GRID_EDGE_TYPE & grid_edge,
   std::vector<COORD_TYPE0> & coord_array,
   const QUAD_INDEX_TYPE iquad,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE NUM_VERT_PER_QUAD(4);

    const ISOV_INDEX_TYPE * first_quad_vert = 
      &(quad_vert[iquad*NUM_VERT_PER_QUAD]);

    compute_scaled_intersection_of_bilinear_surface_and_grid_edge
      (grid, first_quad_vert, grid_edge, coord_array, intersection_coord);

  }


  /*!
   *  @brief Project quad vertices on grid edge and average projections.
   *  - Use average location as approximation of intersection of
   *    quad and grid edge.
   *  @param dual_edge Edge dual to quadrilateral.
   *    - dual_edge.direction = Direction of dual edge.
   *  @param[out] intersection_coord[] Intersection coordinates.
   */
  template <typename CTYPE0, typename ISOV_INDEX_TYPE, 
            typename ISOPOLY_INFO_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1>
  void compute_intersection_of_quad_and_grid_edge_avg
  (const CTYPE0 endpoint0_coord[],
   const int edge_direction,
   const ISOV_INDEX_TYPE quad_vert[],
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 intersection_coord[])
  {
    const int DIM3(3);
    const int NUM_VERT_PER_QUAD(4);
    COORD_TYPE x;

    IJK::copy_coord_3D(endpoint0_coord, intersection_coord);

    x = 0;
    for (int i = 0; i < NUM_VERT_PER_QUAD; i++) {
      const ISOV_INDEX_TYPE isov = quad_vert[i];
      const COORD_TYPE c = coord_array[isov*DIM3+edge_direction];
      x += c;
    }

    intersection_coord[edge_direction] = x/float(NUM_VERT_PER_QUAD);
  }

  
  /*!
   *  @brief Project quad vertices on grid edge and average projections.
   *  - Use average location as approximation of intersection of
   *    quad and grid edge.
   *  @param dual_edge Edge dual to quadrilateral.
   *    - dual_edge.direction = Direction of dual edge.
   *  @param[out] intersection_coord[] Intersection coordinates.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename ISOPOLY_INFO_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1>
  void compute_intersection_of_quad_and_grid_edge_avg
  (const GRID_TYPE & grid, 
   const ISOV_INDEX_TYPE quad_vert[],
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const int DIM3(3);
    const VTYPE iv0 = isopoly_info.GridEdgeEndpoint0();
    const int edge_direction = isopoly_info.GridEdgeDirection();
    COORD_TYPE0 endpoint0_coord[DIM3];

    grid.ComputeCoord(iv0, endpoint0_coord);

    compute_intersection_of_quad_and_grid_edge_avg
      (endpoint0_coord, edge_direction, quad_vert, isopoly_info,
       coord_array, intersection_coord);
  }


  /*!
   *  @overload
   *  @brief Compute intersection of quad and grid edge.
   *  - Version with array of quad vertices.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename GRID_EDGE_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1, typename QUAD_INDEX_TYPE>
  void compute_intersection_of_quad_and_grid_edge_avg
  (const GRID_TYPE & grid, 
   const std::vector<ISOV_INDEX_TYPE> & quad_vert,
   const GRID_EDGE_TYPE & grid_edge,
   std::vector<COORD_TYPE0> & coord_array,
   const QUAD_INDEX_TYPE iquad,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE NUM_VERT_PER_QUAD(4);

    const ISOV_INDEX_TYPE * first_quad_vert = 
      &(quad_vert[iquad*NUM_VERT_PER_QUAD]);

    compute_intersection_of_quad_and_grid_edge_avg
      (grid, first_quad_vert, grid_edge, coord_array, intersection_coord);

  }


  /*!
   *  @brief Project scaled quad vertices on grid edge and average projections.
   *  - Scale vertex coordinates by grid.Spacing().
   *  @pre Coordinates in coord_array[] are already scaled by grid.Spacing().
   *  - Use average location as approximation of intersection of
   *    quad and grid edge.
   *  @param dual_edge Edge dual to quadrilateral.
   *    - dual_edge.direction = Direction of dual edge.
   *  @param[out] intersection_coord[] Intersection coordinates.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename ISOPOLY_INFO_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1>
  void compute_scaled_intersection_of_quad_and_grid_edge_avg
  (const GRID_TYPE & grid, 
   const ISOV_INDEX_TYPE quad_vert[],
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const std::vector<COORD_TYPE0> & coord_array,
   COORD_TYPE1 intersection_coord[])
  {
        typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const int DIM3(3);
    const VTYPE iv0 = isopoly_info.GridEdgeEndpoint0();
    const int edge_direction = isopoly_info.GridEdgeDirection();
    COORD_TYPE0 endpoint0_coord[DIM3];

    grid.ComputeScaledCoord(iv0, endpoint0_coord);

    compute_intersection_of_quad_and_grid_edge_avg
      (endpoint0_coord, edge_direction, quad_vert, isopoly_info,
       coord_array, intersection_coord);
  }


  /*!
   *  @overload
   *  @brief Compute scaled intersection of quad and grid edge.
   *  - Scale vertex coordinates by grid.Spacing().
   *  - Version with array of quad vertices.
   */
  template <typename GRID_TYPE, typename ISOV_INDEX_TYPE, 
            typename GRID_EDGE_TYPE, typename COORD_TYPE0, 
            typename COORD_TYPE1, typename QUAD_INDEX_TYPE>
  void compute_scaled_intersection_of_quad_and_grid_edge_avg
  (const GRID_TYPE & grid, 
   const std::vector<ISOV_INDEX_TYPE> & quad_vert,
   const GRID_EDGE_TYPE & grid_edge,
   std::vector<COORD_TYPE0> & coord_array,
   const QUAD_INDEX_TYPE iquad,
   COORD_TYPE1 intersection_coord[])
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE NUM_VERT_PER_QUAD(4);

    const ISOV_INDEX_TYPE * first_quad_vert = 
      &(quad_vert[iquad*NUM_VERT_PER_QUAD]);

    compute_scaled_intersection_of_quad_and_grid_edge_avg
      (grid, first_quad_vert, grid_edge, coord_array, intersection_coord);

  }

  //@}



  // ***********************************************************************
  //! @name COMPUTE ENVELOPES
  // ***********************************************************************

  //@{

  /*!
   *  @brief Determine which isosurface quadrilateral diagonals
   *    are outside envelope.
   *  @pre Quadrilateral vertices are listed in circular order.
   *  @param[out] isoquad_info[iquad] Information 
   *    about isosurface quadrilateral iquad.
   *    - Routine sets isoquad_info[iquad].is_diagonal_in_envelope[0]
   *      and isoquad_info[iquad].is_diagonal_in_envelope[1].
   */
  template <typename GRID_TYPE, typename ISOVERT_INDEX_TYPE,
            typename ISOVERT_TYPE, typename ISOQUAD_INFO_TYPE,
            typename CTYPE>
  void compute_isoquad_diagonals_deep_in_envelope
  (const GRID_TYPE & grid,
   const std::vector<ISOVERT_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOVERT_TYPE> & isov_list,
   const std::vector<CTYPE> & vertex_coord,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISOQUAD_INFO_TYPE> & isoquad_info,
   DUALISO_INFO & dualiso_info)
  {
    const COORD_TYPE envelope_offset =
        param.envelope_offset;
    const bool flag_use_relative_positions =
      param.envelope_flag_use_relative_positions;

    const std::clock_t t0 = clock();
    
    int num_quads_with_diag_outside_envelope(0);
      
    compute_quad_diagonals_deep_in_envelope
      (grid, isopoly_vert, isov_list, vertex_coord,
       envelope_offset, flag_use_relative_positions,
       isoquad_info, num_quads_with_diag_outside_envelope);

    const std::clock_t t1 = clock();
    
    dualiso_info.triangulation.num_iso_cubes_with_diag_outside_envelope =
      num_quads_with_diag_outside_envelope;
    IJK::clock2seconds(t1-t0, dualiso_info.time.envelope);
  }


  /*!
   *  @brief Determine which isosurface quadrilateral diagonals
   *    are outside envelope.
   *  @pre Quadrilateral vertices are listed in circular order.
   *  @param dual_isosurface Data structure storing dual isosurface.
   *     @pre Quadrilateral vertices are listed in circular order.
   *  @param[out] isoquad_info[iquad] Information 
   *    about isosurface quadrilateral iquad.
   *    - Routine sets isoquad_info[iquad].is_diagonal_in_envelope[0]
   *      and isoquad_info[iquad].is_diagonal_in_envelope[1].
   */
  template <typename GRID_TYPE, typename DUAL_ISOSURFACE_TYPE,
            typename ISOVERT_TYPE, typename ISOQUAD_INFO_TYPE>
  void compute_isoquad_diagonals_deep_in_envelope
  (const GRID_TYPE & grid,
   const DUAL_ISOSURFACE_TYPE & dual_isosurface,
   const std::vector<ISOVERT_TYPE> & isov_list,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISOQUAD_INFO_TYPE> & isoquad_info,
   DUALISO_INFO & dualiso_info)
  {
    if (dual_isosurface.CubeVertexOrder() != CIRCULAR_VERTEX_ORDER) {
      IJK::PROCEDURE_ERROR
        error("compute_isoquad_diagonals_deep_in_envelope");

      error.AddMessage
        ("Programming error. Quad vertices are NOT in circular order.");
      throw error;
    }
        
    compute_isoquad_diagonals_deep_in_envelope
      (grid, dual_isosurface.isopoly_vert, isov_list,
       dual_isosurface.vertex_coord, param,
       isoquad_info, dualiso_info);
  }
    
  //@}

  
  // ***********************************************************************
  //! @name ADD TRI4 INTERIOR VERTICES ON GRID EDGES
  // ***********************************************************************

  //@{
  
  /// @brief Add isosurface vertices on grid edges.
  template <typename DATA_TYPE, typename STYPE,
            typename ISOPOLY_INFO_TYPE, typename ISOSURFACE_TYPE>
  void add_isov_on_grid_edges
  (const DATA_TYPE & dualiso_data, const STYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface)
  {
    const IJK::POLY_EDGE_INTERSECTION_METHOD intersection_method =
      dualiso_data.interior_vertex_param.poly_edge_intersection_method;

    if (intersection_method == IJK::POLY_EDGE_MULTILINEAR_INTERPOLATION) {
      add_isov_on_grid_edges_bilinear_surface
        (dualiso_data, isopoly_info, dual_isosurface);
    }
    else if (intersection_method == IJK::POLY_EDGE_AVERAGE_PROJECTION) {
      add_isov_on_grid_edges_average_projection
        (dualiso_data, isopoly_info, dual_isosurface);
    }
    else {
      add_isov_on_grid_edges_interpolate_scalar
        (dualiso_data, isovalue, isopoly_info, dual_isosurface);
    }
  }


  /// @brief Add isosurface vertices on grid edges
  ///   using linear interpolation on scalar values.
  template <typename DATA_TYPE, typename STYPE,
            typename ISOPOLY_INFO_TYPE, typename ISOSURFACE_TYPE>
  void add_isov_on_grid_edges_interpolate_scalar
  (const DATA_TYPE & dualiso_data, const STYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface)
  {
    const int dimension = dual_isosurface.Dimension();
    const bool flag_move_all_isov_away_from_cube_facets =
      dualiso_data.iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets;
    const bool flag_move_all_isov_away_from_cube_ridges =
      dualiso_data.iso_vertex_position_param.flag_move_all_isov_away_from_cube_ridges;
    IJK::PROCEDURE_ERROR error("add_isov_on_grid_edges_interpolate_scalar");

    if (!dualiso_data.CheckIsoQuadDualToGridEdges
        ("add_isov_on_grid_edges_interpolate_scalar", error))
      { throw error; }

    if (!dual_isosurface.CheckIsoPolyInfo(isopoly_info,error))
      { throw error; }
    
    dual_isosurface.AllocateIsovDualToIsopoly();
    const ISO_VERTEX_INDEX first_isov_on_edge =
      dual_isosurface.FirstIsovDualToIsopoly();

    if (dual_isosurface.vertex_coord.size() > 0) {

      COORD_TYPE * first_isov_on_edge_coord =
        IJK::vector2pointerNC(dual_isosurface.vertex_coord)
        + dimension*first_isov_on_edge;

      if (flag_move_all_isov_away_from_cube_facets ||
          flag_move_all_isov_away_from_cube_ridges) {

        const COORD_TYPE offset =
          dualiso_data.iso_vertex_position_param.position_offset;
        
        IJK::compute_all_isov_coord_on_grid_edge_linear
          (dualiso_data.ScalarGrid(), isovalue, isopoly_info,
           offset, first_isov_on_edge_coord);        
      }
      else {
        IJK::compute_all_isov_coord_on_grid_edge_linear
          (dualiso_data.ScalarGrid(), isovalue, isopoly_info,
           first_isov_on_edge_coord);
      }
    }

  }


  /*!
   *  @brief Add isosurface vertices on grid edges 
   *    using bilinear interpolation.
   *  - Locate isosurface vertices at intersection of grid edges 
   *    and bilinear surface patches.
   */
  template <typename DATA_TYPE,
            typename ISOPOLY_INFO_TYPE, typename ISOSURFACE_TYPE>
  void add_isov_on_grid_edges_bilinear_surface
  (const DATA_TYPE & dualiso_data,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface)
  {
    const int dimension = dual_isosurface.Dimension();
    IJK::PROCEDURE_ERROR error("add_isov_on_grid_edges_bilinear_surface");

    if (!dualiso_data.CheckIsoQuadDualToGridEdges
        ("add_isov_on_grid_edges_bilinear_surface", error))
      { throw error; }

    if (!dual_isosurface.CheckIsoPolyInfo(isopoly_info,error))
      { throw error; }

    dual_isosurface.AllocateIsovDualToIsopoly();
    const ISO_VERTEX_INDEX first_isov_on_edge =
      dual_isosurface.FirstIsovDualToIsopoly();

    if (dual_isosurface.vertex_coord.size() > 0) {

      for (int ipoly = 0; ipoly < dual_isosurface.NumIsoPoly(); ipoly++) {

        const ISO_VERTEX_INDEX isov = first_isov_on_edge + ipoly;
        COORD_TYPE * isov_coord =
          IJK::vector2pointerNC(dual_isosurface.vertex_coord)
          + isov*dimension;

        compute_intersection_of_bilinear_surface_and_grid_edge
          (dualiso_data.ScalarGrid(), dual_isosurface.isopoly_vert,
           isopoly_info[ipoly], dual_isosurface.vertex_coord, 
           ipoly, isov_coord);        
      }
    }
  }


  /*!
   *  @brief Add isosurface vertices on grid edges 
   *    by averaging projected vertices.
   *  - Compute isosurface vertex location by projecting 
   *    quadrilateral vertices onto a grid edge and averaging the
   *    projection locations along the grid edge.
   */
  template <typename DATA_TYPE,
            typename ISOPOLY_INFO_TYPE, typename ISOSURFACE_TYPE>
  void add_isov_on_grid_edges_average_projection
  (const DATA_TYPE & dualiso_data,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface)
  {
    const int dimension = dual_isosurface.Dimension();
    IJK::PROCEDURE_ERROR error("add_isov_on_grid_edges_average_projection");

    if (!dualiso_data.CheckIsoQuadDualToGridEdges
        ("add_isov_on_grid_edges_interpolate_coord", error))
      { throw error; }
    
    if (!dual_isosurface.CheckIsoPolyInfo(isopoly_info,error))
      { throw error; }

    dual_isosurface.AllocateIsovDualToIsopoly();
    const ISO_VERTEX_INDEX first_isov_on_edge =
      dual_isosurface.FirstIsovDualToIsopoly();
    
    if (dual_isosurface.vertex_coord.size() > 0) {

      for (int ipoly = 0; ipoly < dual_isosurface.NumIsoPoly(); ipoly++) {

        const ISO_VERTEX_INDEX isov = first_isov_on_edge + ipoly;
        COORD_TYPE * isov_coord =
          IJK::vector2pointerNC(dual_isosurface.vertex_coord)
          + isov*dimension;

        compute_intersection_of_quad_and_grid_edge_avg
          (dualiso_data.ScalarGrid(), dual_isosurface.isopoly_vert,
           isopoly_info[ipoly], dual_isosurface.vertex_coord, 
           ipoly, isov_coord);        
      }
    }
  }

  //@}

  
  // ***********************************************************************
  //! @name ADD TRI4 INTERIOR VERTICES AT QUADRILATERAL CENTROIDS
  // ***********************************************************************


  //@{
  
  /// Add new isosurface vertices at centroid of isosurface poly vertices.
  template <typename DATA_TYPE, typename ISOSURFACE_TYPE>
  void add_isov_at_poly_centroids
  (const DATA_TYPE & dualiso_data, ISOSURFACE_TYPE & dual_isosurface)
  {
    const int dimension = dual_isosurface.Dimension();
    const int numv_per_iso_poly = dual_isosurface.NumVerticesPerIsoPoly();
    IJK::PROCEDURE_ERROR error("add_isov_at_poly_centroids");

    dual_isosurface.AllocateIsovDualToIsopoly();
    const ISO_VERTEX_INDEX first_new_isov =
      dual_isosurface.FirstIsovDualToIsopoly();
    
    if (dual_isosurface.vertex_coord.size() > 0) {

      for (int ipoly = 0; ipoly < dual_isosurface.NumIsoPoly(); ipoly++) {

        const ISO_VERTEX_INDEX isov = first_new_isov + ipoly;
        COORD_TYPE * isov_coord =
          &(dual_isosurface.vertex_coord[isov*dimension]);

        const ISO_VERTEX_INDEX * first_isopoly_vert = 
          &(dual_isosurface.isopoly_vert[ipoly*numv_per_iso_poly]);
        
        IJK::compute_quad_centroid
          (dimension, first_isopoly_vert, dual_isosurface.vertex_coord, 
           isov_coord);
      }
    }

  }


  /// Add new vertex coordinate at edge midpoint.
  template <typename VTYPE0, typename VTYPE1, typename CTYPE>
  VTYPE0 add_vertex_coord_at_edge_midpoint
  (const VTYPE0 edge_endpoint0, const VTYPE1 edge_endpoint1,
   std::vector<CTYPE> & vertex_coord)
  {
    const int DIM3(3);
    const COORD_TYPE * vcoord0 = &(vertex_coord[edge_endpoint0*DIM3]);
    const COORD_TYPE * vcoord1 = &(vertex_coord[edge_endpoint1*DIM3]);
    CTYPE mid_coord[DIM3];

    IJK::compute_midpoint_3D(vcoord0, vcoord1, mid_coord);
    VTYPE0 iv_new = IJK::insert_coord_3D(mid_coord, vertex_coord);

    return(iv_new);
  }


  /// @brief Add new vertex coordinate at midpoint of line segment 
  ///   whose endpoints are edge midpoints.
  template <typename VTYPEA0, typename VTYPEA1, 
            typename VTYPEB0, typename VTYPEB1, typename CTYPE>
  VTYPEA0 add_vertex_coord_at_midpoint_of_edge_midpoints
  (const VTYPEA0 edgeA_endpoint0, const VTYPEA1 edgeA_endpoint1,
   const VTYPEB0 edgeB_endpoint0, const VTYPEB1 edgeB_endpoint1,
   std::vector<CTYPE> & vertex_coord)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int DIM3(3);
    const int NUM_VERT_PER_QUAD(4);
    const VTYPEA0 quad_vert[NUM_VERT_PER_QUAD] =
      { edgeA_endpoint0, edgeA_endpoint1, edgeB_endpoint0, edgeB_endpoint1 };
    const SIZE_TYPE vertex_coord_size = vertex_coord.size();
    const VTYPEA0 iv_new = vertex_coord_size/DIM3;

    vertex_coord.resize(vertex_coord_size + DIM3);
    COORD_TYPE * new_coord = &(vertex_coord[vertex_coord_size]);

    IJK::compute_quad_centroid(DIM3, quad_vert, vertex_coord, new_coord);

    return(iv_new);
  }



  /// Add new vertex coordinates at centroid of quad vertices.
  template <typename VTYPE, typename CTYPE>
  VTYPE add_vertex_coord_at_quad_centroid
  (const VTYPE quad_vert[],
   std::vector<CTYPE> & vertex_coord)
  {
    const int DIM3(3);
    COORD_TYPE quad_centroid[DIM3];

    IJK::compute_quad_centroid(DIM3, quad_vert, vertex_coord, quad_centroid);
    VTYPE iv_new = IJK::insert_coord_3D(quad_centroid, vertex_coord);

    return(iv_new);
  }


  /// Add new vertex coordinates at centroid of pentagon vertices.
  template <typename VTYPE, typename NTYPE, typename CTYPE>
  void add_vertex_coord_at_pentagon_centroid
  (const VTYPE pentagon_vert[],
   const NTYPE num_pentagon,
   std::vector<CTYPE> & vertex_coord)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int DIM3(3);
    const int NUM_VERT_PER_PENTAGON(5);
    const SIZE_TYPE vertex_coord_size = vertex_coord.size();

    vertex_coord.resize(vertex_coord_size + DIM3*num_pentagon);

    for (int ipoly = 0; ipoly < num_pentagon; ipoly++) {

      COORD_TYPE * new_coord =
        &(vertex_coord[vertex_coord_size + ipoly*DIM3]);

      const ISO_VERTEX_INDEX * first_pentagon_vert = 
        &(pentagon_vert[ipoly*NUM_VERT_PER_PENTAGON]);
        
      IJK::compute_pentagon_centroid
        (DIM3, first_pentagon_vert, vertex_coord, new_coord);
    }

  }


  /// @brief Add new vertex coordinates at centroid of pentagon vertices.
  /// - Version using C++ STL vector for array pentagon_vert[].
  template <typename VTYPE, typename CTYPE>
  void add_vertex_coord_at_pentagon_centroid
  (const std::vector<VTYPE> & pentagon_vert,
   std::vector<CTYPE> & vertex_coord)
  {
    const int NUM_VERT_PER_PENTAGON(5);
    const int num_pentagon = pentagon_vert.size()/NUM_VERT_PER_PENTAGON;

    add_vertex_coord_at_pentagon_centroid
      (IJK::vector2pointer(pentagon_vert), num_pentagon, vertex_coord);
  }


  // *** NEED TO TEST THIS CODE ***
  /// Add new vertex coordinates at centroid of hexahedra vertices.
  template <typename DTYPE, typename VTYPE, 
            typename NTYPE, typename ITYPE, typename CTYPE>
  void add_vertex_coord_at_hexahedra_centroid
  (const DTYPE dimension,
   const VTYPE hexahedra_vert[],
   const NTYPE num_hexahedra,
   const ITYPE first_isov_dual_to_hex,
   std::vector<CTYPE> & vertex_coord)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int NUM_VERT_PER_HEXAHEDRA(8);
    const SIZE_TYPE vertex_coord_size = vertex_coord.size();

    const VERTEX_INDEX * first_coord =
      IJK::vector2pointerNC(vertex_coord) +
      dimension*first_isov_dual_to_hex;
    
    for (int ipoly = 0; ipoly < num_hexahedra; ipoly++) {

      COORD_TYPE * new_coord = first_coord + dimension*ipoly;

      const ISO_VERTEX_INDEX * first_hexahedron_vert = 
        hexahedra_vert + ipoly*NUM_VERT_PER_HEXAHEDRA;
        
      IJK::compute_hexahedron_centroid
        (dimension, first_hexahedron_vert, vertex_coord, new_coord);
    }

  }


  /// Add new vertex coordinates at centroid of hexahedra vertices.
  template <typename DATA_TYPE, typename ISOSURFACE_TYPE>
  void add_vertex_coord_at_hexahedra_centroid
  (const DATA_TYPE & dualiso_data, ISOSURFACE_TYPE & dual_isosurface)
  {
    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const ISO_VERTEX_INDEX num_hexahedra = dual_isosurface.NumIsoPoly();

    dual_isosurface.AllocateIsovDualToIsopoly();
    const ISO_VERTEX_INDEX first_isov_dual_to_hex =
      dual_isosurface.FirstIsovDualToIsopoly();
    
    const ISO_VERTEX_INDEX * first_hex_vert = 
      &(dual_isosurface.isopoly_vert.front());

    add_vertex_coord_at_hexahedra_centroid
      (dimension, first_hex_vert, num_hexahedra,
       first_isov_dual_to_hex,
       dual_isosurface.vertex_coord);
  }

  //@}


  // ***********************************************************************
  //! @name ADD TRI4 INTERIOR VERTICES INSIDE ENVELOPES
  // ***********************************************************************

  //@{
  
  namespace {

    template <typename VTYPE, typename DIR_TYPE, typename CTYPEV,
              typename CTYPEMIN, typename CTYPEMAX>
    void get_min_max_quad_vert_coord_3D
    (const VTYPE quad_vert[], const DIR_TYPE d,
     const CTYPEV * isov_coord,
     CTYPEMIN & cmin, CTYPEMAX & cmax)
    {
      const int DIM3(3);
      const int NUM_VERT_PER_QUADRILATERAL(4);
      
      const VTYPE iv0 = quad_vert[0];
      cmax = cmin = isov_coord[iv0*DIM3 + d];

      for (int i = 1; i < NUM_VERT_PER_QUADRILATERAL; i++) {
        const VTYPE iv = quad_vert[i];
        const CTYPEV c = isov_coord[iv*DIM3+d];
        if (cmin > c) { cmin = c; }
        if (cmax < c) { cmax = c; }
      }
    }


    /*!
     *  @brief Return min distance from quad vertices to grid facet.
     *  - Return min distance from (quad_vert[i0], quad_vert[i1])
     *    to grid facet.
     *  - Return min distance from (quad_vert[i2], quad_vert[i3])
     *    to grid facet.
     *  @pre quad_vert[0] is in quadrant left/below grid vertex.
     *  - Other quad vertices could be in any quadrant.
     *  @param facet_coord[] Coordinates of lowest/leftmost
     *    facet vertex.
     *  @param[out] distLB Minimum distance from quad vertices
     *    left/below facet to facet.
     *  @param[out] distRA Minimum distance  from quad vertices
     *    right/above facet to facet.
     */
    template <typename VTYPE, typename DIR_TYPE,
              typename CTYPEV, typename CTYPEF,
              typename DIST_TYPELB, typename DIST_TYPERA>
    void get_min_dist_quad_vert_to_facet_3D
    (const VTYPE quad_vert[], const DIR_TYPE d,
     const CTYPEV isov_coord[], const CTYPEF facet_coord[],
     DIST_TYPELB & min_distLB, DIST_TYPERA & min_distRA)
    {
      const int DIM3(3);
      const int NUM_VERT_PER_QUADRILATERAL(4);
      const int TWO(2);
      
      const CTYPEF c = facet_coord[d];
      const CTYPEV quad_vcoord[NUM_VERT_PER_QUADRILATERAL] =
        { isov_coord[quad_vert[0]*DIM3 + d],
          isov_coord[quad_vert[1]*DIM3 + d],
          isov_coord[quad_vert[2]*DIM3 + d],
          isov_coord[quad_vert[3]*DIM3 + d] };

      // Number of coordinates less than c.
      int num_lt = 0;

      // Number of coordinates greater than c.
      int num_gt = 0;

      for (int i = 0; i < NUM_VERT_PER_QUADRILATERAL; i++) {
        if (quad_vcoord[i] < c) {
          const CTYPEF dist = c - quad_vcoord[i];
          if (num_lt == 0)
            { min_distLB = dist; }
          else
            { min_distLB = std::min(min_distLB, dist); }

          num_lt++;
        }
        else if (quad_vcoord[i] > c) {
          const CTYPEF dist = quad_vcoord[i] - c;
          if (num_gt == 0)
            { min_distRA = dist; }
          else
            { min_distRA = std::min(min_distRA, dist); }
          num_gt++;
        }
      }
      
      if (num_lt < TWO) {
        if (num_gt < TWO) {
          // Facet contains isosurface vertices from quads
          // left&right c or below&above c.
          min_distLB = 0.0;
          min_distRA = 0.0;
        }
        else {
          // Facet contains isosurface vertex from quad
          // left/below c.
          min_distLB = 0.0;
        }
      }
      else if (num_gt < TWO) {
        // Facet contains isosurface vertex from quad
        // left/below c.
        min_distRA = 0.0;
      }
    }


    /*!
     *  @brief Position tri4 interior coord.
     *  @param d Coordinate index.
     *  @param min_distLB Min distance from quad vertices
     *    left/below facet to facet.
     *  @param min_distRA Min distance from quad vertices
     *    right/above facet to facet.
     *  @param tri4_interior_coord[] Coordinate of tri4
     *    interior vertex.
     */
    template <typename DIR_TYPE, typename CTYPEF, typename CTYPET,
              typename DIST_TYPE, typename RTYPE>
    void position_tri4_interior_coord_in_edge_envelope
    (const DIR_TYPE d, const CTYPEF facet_coord[],
     const DIST_TYPE min_distLB, const DIST_TYPE min_distRA,
     const RTYPE tri4_ratio_to_cube_facet,
     CTYPET tri4_interior_coord[])
    {
      if (min_distLB <= min_distRA) {
        // Place tri4 interior vert to right/above of facet_coord[].
        const COORD_TYPE displacement =
          tri4_ratio_to_cube_facet * (min_distRA-min_distLB);
        tri4_interior_coord[d] = facet_coord[d] + displacement;
      }
      else {
        // Place tri4 interior vert to left/below of facet_coord[].
        const COORD_TYPE displacement =
          tri4_ratio_to_cube_facet * (min_distLB-min_distRA);
        tri4_interior_coord[d] = facet_coord[d] - displacement;
      }
    }
      
  };

  
  /// Compute location of tri4 interior vertex in quad-edge envelope.
  template <typename GRID_TYPE, typename VTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename ISOPOLY_INFO_TYPE, typename RTYPE>
  void compute_tri4_coord_in_quad_edge_envelope
  (const GRID_TYPE & grid, 
   const VTYPE quad_vert[],
   const COORD_TYPE0 * isov_coord,
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const RTYPE tri4_ratio_to_cube_facet,
   COORD_TYPE1 tri4_interior_coord[])
  {
    const int DIM3(3);
    const int edge_direction = isopoly_info.GridEdgeDirection();

    const VTYPE iv0 = isopoly_info.GridEdgeEndpoint0();
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;
    COORD_TYPE0 facet_coord[DIM3];
    COORD_TYPE0 cmin, cmax;
    COORD_TYPE0 min_distLB, min_distRA;

    grid.ComputeCoord(iv0, facet_coord);

    get_min_max_quad_vert_coord_3D
      (quad_vert, edge_direction, isov_coord, cmin, cmax);
    tri4_interior_coord[edge_direction] = (cmin+cmax)/2.0;

    get_min_dist_quad_vert_to_facet_3D
      (quad_vert, d1, isov_coord, facet_coord,
       min_distLB, min_distRA);

    position_tri4_interior_coord_in_edge_envelope
      (d1, facet_coord, min_distLB, min_distRA,
       tri4_ratio_to_cube_facet, tri4_interior_coord);
    
    get_min_dist_quad_vert_to_facet_3D
      (quad_vert, d2, isov_coord, facet_coord,
       min_distLB, min_distRA);

    position_tri4_interior_coord_in_edge_envelope
      (d2, facet_coord, min_distLB, min_distRA,
       tri4_ratio_to_cube_facet, tri4_interior_coord);
  }


  /*!
   *  @overload
   *  @brief Compute location of tri4 interior vertex in quad-edge envelope.
   *  - Version using C++ STL vector for coord_array[].
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename ISOPOLY_INFO_TYPE, typename RTYPE>
  void compute_tri4_coord_in_quad_edge_envelope
  (const GRID_TYPE & grid, 
   const VTYPE quad_vert[], 
   const std::vector<COORD_TYPE0> & coord_array,
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const RTYPE tri4_ratio_to_cube_facet,
   COORD_TYPE1 tri4_interior_coord[])
  {
    compute_tri4_coord_in_quad_edge_envelope
      (grid, quad_vert, IJK::vector2pointer(coord_array),
       isopoly_info, tri4_ratio_to_cube_facet,
       tri4_interior_coord);
  }


  /*!
   *  @brief Compute scaled location of tri4 interior vertex in quad-edge envelope.
   *  - Quad coordinates scaled by grid.Spacing().
   *  @pre Coordinates in vertex_coord[] are already scaled by grid.Spacing().
   */
  template <typename GRID_TYPE, typename VTYPE, 
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename ISOPOLY_INFO_TYPE, typename RTYPE>
  void compute_scaled_tri4_coord_in_quad_edge_envelope
  (const GRID_TYPE & grid, 
   const VTYPE quad_vert[],
   const COORD_TYPE0 * isov_coord,
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const RTYPE tri4_ratio_to_cube_facet,
   COORD_TYPE1 tri4_interior_coord[])
  {
    const int DIM3(3);
    const int edge_direction = isopoly_info.GridEdgeDirection();

    const VTYPE iv0 = isopoly_info.GridEdgeEndpoint0();
    const int d1 = (edge_direction+1)%DIM3;
    const int d2 = (edge_direction+2)%DIM3;
    COORD_TYPE0 facet_coord[DIM3];
    COORD_TYPE0 cmin, cmax;
    COORD_TYPE0 min_distLB, min_distRA;

    grid.ComputeScaledCoord(iv0, facet_coord);

    get_min_max_quad_vert_coord_3D
      (quad_vert, edge_direction, isov_coord, cmin, cmax);
    tri4_interior_coord[edge_direction] = (cmin+cmax)/2.0;

    get_min_dist_quad_vert_to_facet_3D
      (quad_vert, d1, isov_coord, facet_coord,
       min_distLB, min_distRA);

    position_tri4_interior_coord_in_edge_envelope
      (d1, facet_coord, min_distLB, min_distRA,
       tri4_ratio_to_cube_facet, tri4_interior_coord);
    
    get_min_dist_quad_vert_to_facet_3D
      (quad_vert, d2, isov_coord, facet_coord,
       min_distLB, min_distRA);

    position_tri4_interior_coord_in_edge_envelope
      (d2, facet_coord, min_distLB, min_distRA,
       tri4_ratio_to_cube_facet, tri4_interior_coord);
  }


  /*!
   *  @overload
   *  @brief Compute scaled location of tri4 interior vertex in quad-edge envelope.
   *  - Version using C++ STL vector for coord_array[].
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename COORD_TYPE0, typename COORD_TYPE1,
            typename ISOPOLY_INFO_TYPE, typename RTYPE>
  void compute_scaled_tri4_coord_in_quad_edge_envelope
  (const GRID_TYPE & grid, 
   const VTYPE quad_vert[], 
   const std::vector<COORD_TYPE0> & coord_array,
   const ISOPOLY_INFO_TYPE & isopoly_info,
   const RTYPE tri4_ratio_to_cube_facet,
   COORD_TYPE1 tri4_interior_coord[])
  {
    compute_scaled_tri4_coord_in_quad_edge_envelope
      (grid, quad_vert, IJK::vector2pointer(coord_array),
       isopoly_info, tri4_ratio_to_cube_facet,
       tri4_interior_coord);
  }
  
  
  /// @brief Add new isosurface vertices inside quad-edge envelopes.
  template <typename DATA_TYPE, typename ISOPOLY_INFO_TYPE,
            typename ISOSURFACE_TYPE>
  void add_isov_in_quad_edge_envelopes
  (const DATA_TYPE & dualiso_data,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface)
  {
    const int dimension = dual_isosurface.Dimension();
    const int numv_per_iso_poly = dual_isosurface.NumVerticesPerIsoPoly();
    const COORD_TYPE tri4_ratio_to_cube_facet =
      dualiso_data.interior_vertex_param.interior_vertex_ratio_to_cube_facet;
    IJK::PROCEDURE_ERROR error("add_isov_in_quad_edge_envelopes");

    dual_isosurface.AllocateIsovDualToIsopoly();

    const ISO_VERTEX_INDEX first_new_isov =
      dual_isosurface.FirstIsovDualToIsopoly();

    if (dual_isosurface.vertex_coord.size() > 0) {

      for (int ipoly = 0; ipoly < dual_isosurface.NumIsoPoly(); ipoly++) {
        const ISO_VERTEX_INDEX isov = first_new_isov + ipoly;
        COORD_TYPE * isov_coord =
          &(dual_isosurface.vertex_coord[isov*dimension]);

        const ISO_VERTEX_INDEX * first_isopoly_vert = 
          &(dual_isosurface.isopoly_vert[ipoly*numv_per_iso_poly]);

        compute_tri4_coord_in_quad_edge_envelope
          (dualiso_data.ScalarGrid(), first_isopoly_vert,
           dual_isosurface.vertex_coord,
           isopoly_info[ipoly],
           tri4_ratio_to_cube_facet, isov_coord);
      }
    }

  }

  //@}

  
  // ***********************************************************************
  //! @name ADD TRI4 INTERIOR VERTICES
  // ***********************************************************************

  //@{

  /*!
   *  @brief Set scaled tri4 interior vertex on grid edge.
   *    - Scale vertex coordinates by grid.Spacing().
   *  @pre Coordinates in vertex_coord[] are already scaled by grid.Spacing().
   *  @param isoquad_vert[] Vertices of the isosurface quadrilateral.
   *  @param isoquad_info Isosurface quadrilateral information,
   *    including dual grid edge.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename ISOQ_INFO_TYPE,
            typename STYPE, typename CTYPE0, typename CTYPE1>
  void set_scaled_tri4_interior_vertex_on_grid_edge
  (const GRID_TYPE & scalar_grid,
   const VTYPE isoquad_vert[],
   const ISOQ_INFO_TYPE & isoquad_info,
   const STYPE isovalue,
   const INTERIOR_VERTEX_PARAM & interior_vertex_param,
   const std::vector<CTYPE0> & vertex_coord,
   CTYPE1 new_vertex_coord[])
  {
    const IJK::POLY_EDGE_INTERSECTION_METHOD intersection_method =
      interior_vertex_param.poly_edge_intersection_method;

    if (intersection_method == IJK::POLY_EDGE_MULTILINEAR_INTERPOLATION) {
      compute_scaled_intersection_of_bilinear_surface_and_grid_edge
        (scalar_grid, isoquad_vert, isoquad_info, vertex_coord,
         new_vertex_coord);
    }
    else if (intersection_method == IJK::POLY_EDGE_AVERAGE_PROJECTION) {
      compute_scaled_intersection_of_quad_and_grid_edge_avg
        (scalar_grid, isoquad_vert, isoquad_info,
         vertex_coord, new_vertex_coord);
    }
    else {
      IJK::compute_isov_coord_on_scaled_grid_edge_linear_I
        (scalar_grid, isovalue, isoquad_info, new_vertex_coord);
    }
  }

  
  /*!
   *  @brief Set scaled tri4 interior vertex coordinate. 
   *    - Scale vertex coordinates by grid.Spacing().
   *  @pre Coordinates in vertex_coord[] are already scaled by grid.Spacing().
   *  @param isoquad_vert[] Vertices of the isosurface quadrilateral.
   *  @param isoquad_info Isosurface quadrilateral information,
   *    including dual grid edge.
   *  @param[out] new_vertex_coord[] 
   *    Pointer to array to store coordinates of new vertex.
   *    @pre new_vertex_coord[] is preallocated to length 
   *       at least grid.Dimension().
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename ISOQ_INFO_TYPE,
            typename STYPE, typename CTYPE0, typename CTYPE1>
  void set_scaled_tri4_interior_vertex_coord
  (const GRID_TYPE & grid,
   const VTYPE isoquad_vert[],
   const ISOQ_INFO_TYPE & isoquad_info,
   const STYPE isovalue,
   const INTERIOR_VERTEX_PARAM & interior_vertex_param,
   const std::vector<CTYPE0> & vertex_coord,
   CTYPE1 new_vertex_coord[])
  {
    if (interior_vertex_param.PositionAtCentroid()) {
      // new_vertex_coord[] automatically scaled since
      // quad vertex coordinates are already scaled.
      IJK::compute_quad_centroid
        (grid.Dimension(), isoquad_vert, vertex_coord,
         new_vertex_coord);
    }
    else if (interior_vertex_param.PositionInEnvelope()) {
      compute_scaled_tri4_coord_in_quad_edge_envelope
        (grid, isoquad_vert, vertex_coord, isoquad_info,
         interior_vertex_param.interior_vertex_ratio_to_cube_facet,
         new_vertex_coord);
    }
    else {
      set_scaled_tri4_interior_vertex_on_grid_edge
        (grid, isoquad_vert, isoquad_info,
         isovalue, interior_vertex_param,
         vertex_coord, new_vertex_coord);
    }
  }


  /*!
   *  @brief Add one tri4 interior vertex with scaled coordinates. 
   *    - Scale vertex coordinates scaled by grid.Spacing().
   *  @pre Coordinates in vertex_coord[] are already scaled by grid.Spacing().
   *  @param isoquad_vert[] Vertices of the isosurface quadrilateral.
   *  @param isoquad_info Isosurface quadrilateral information,
   *    including dual grid edge.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename ISOP_INFO_TYPE,
            typename STYPE, typename CTYPE>
  VTYPE add_tri4_interior_vertex_scale
  (const GRID_TYPE & grid,
   const VTYPE isopoly_vert[],
   const ISOP_INFO_TYPE & isopoly_info,
   const STYPE isovalue,
   const INTERIOR_VERTEX_PARAM & interior_vertex_param,
   std::vector<CTYPE> & vertex_coord)
  {
    typedef std::vector<VERTEX_INDEX>::size_type SIZE_TYPE;

    const int dimension = grid.Dimension();
    const SIZE_TYPE vertex_coord_size = vertex_coord.size();
    const VTYPE new_isov = (vertex_coord_size/dimension);

    vertex_coord.resize(vertex_coord_size+dimension);
    
    CTYPE * new_vertex_coord =
      IJK::vector2pointerNC(vertex_coord) + vertex_coord_size;

    set_scaled_tri4_interior_vertex_coord
      (grid, isopoly_vert, isopoly_info, isovalue,
       interior_vertex_param, vertex_coord, new_vertex_coord);

    return new_isov;
  }


  /// @brief For each quadrilateral, add a vertex to split
  ///   the quadrilateral into four triangles.
  template <typename DATA_TYPE, typename STYPE,
            typename ISOPOLY_INFO_TYPE, typename ISOSURFACE_TYPE>
  void add_tri4_interior_vertices
  (const DATA_TYPE & dualiso_data, const STYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface)
  {
    if (dualiso_data.interior_vertex_param.PositionAtCentroid()) {
      add_isov_at_poly_centroids(dualiso_data, dual_isosurface); 
    }
    else if (dualiso_data.interior_vertex_param.PositionInEnvelope()) {
      add_isov_in_quad_edge_envelopes
        (dualiso_data, isopoly_info, dual_isosurface);
    }
    else {
      add_isov_on_grid_edges
        (dualiso_data, isovalue, isopoly_info, dual_isosurface); 
    }
  }

  //@}


  // ***********************************************************************
  //! @name CLASS FOR ADDING TRI4 INTERIOR VERTICES
  // ***********************************************************************

  template <typename GRID_TYPE, typename NTYPEV>
  class ADD_ISOV_DUAL_TO_QUAD_SCALE:public IJK::ADD_DUAL_ISOV_BASE {

  protected:

    const GRID_TYPE * grid_ptr;

  public:
    INTERIOR_VERTEX_PARAM interior_vertex_param;

    
  public:

    /// Constructor.
    ADD_ISOV_DUAL_TO_QUAD_SCALE(const GRID_TYPE * _grid_ptr):
      grid_ptr(_grid_ptr)
    {}

    /*!
     *  @brief Add isosurface vertex dual to polytope ipoly.
     *  @param isopoly_vert[] Vertices of the isosurface polytope.
     *  @param isopoly_info Information about isosurface polytope.
     */
    template <typename DTYPE, typename VTYPE,
              typename ISOP_INFO_TYPE, typename STYPE,
              typename CTYPE>
    VTYPE AddDualIsov
    (const DTYPE dimension, const VTYPE isopoly_vert[],
     const ISOP_INFO_TYPE & isopoly_info,
     const STYPE isovalue, std::vector<CTYPE> & vertex_coord) const;
  };

  
  // ***********************************************************************
  //! @name CONVERT QUADRILATERALS TO TRIANGLES
  // ***********************************************************************

  //@{

  /*!
   *  @brief Triangulate isosurface quadrilaterals.
   *  - Assumes input isosurface is quadrilaterals embedded in 3D.
   *  - Reorders isopoly_vert[] so that quad vertices are
   *    listed in circular order around the quadrilateral.
   *  @param flag_skip_collapsed_quads If true, skip any quads
   *    that have been collapsed to a single vertex.
   *    - (Still triangulates quads that are collapsed to an edge
   *       or two edges.)
   *    - Default: true.
   */
  template <typename DATA_TYPE, typename ISOVAL_TYPE,
            typename ISOPOLY_INFO_TYPE,
            typename ISOSURFACE_TYPE, typename INFO_TYPE>
  void triangulate_isosurface_quadrilaterals
  (const DATA_TYPE & dualiso_data, const ISOVAL_TYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface,
   INFO_TYPE & triangulation_info,
   bool flag_skip_collapsed_quads=true)
  {
    const int NUM_VERT_PER_QUAD = 4;
    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const COORD_TYPE max_small_magnitude = dualiso_data.MaxSmallMagnitude();
    const QUAD_TRI_METHOD quad_tri_method = 
      dualiso_data.QuadTriangulationMethod();
    MESH_SIZE_TYPE num_triangulated_quad;
    IJK::PROCEDURE_ERROR error("triangulate_isosurface_quadrilaterals");

    if (!dual_isosurface.CheckNumVerticesPerIsoPoly
        (NUM_VERT_PER_QUAD, error)) { throw error; }

    dual_isosurface.SetQuadVertexOrder(CIRCULAR_VERTEX_ORDER);

    if (quad_tri_method == MAX_MIN_ANGLE) {      
      IJK::triangulate_quad_tri2_max_min_angle
        (dimension, dual_isosurface.vertex_coord,
         dual_isosurface.isopoly_vert, 
         max_small_magnitude, dual_isosurface.simplex_vert,
         num_triangulated_quad, flag_skip_collapsed_quads);

      triangulation_info.num_iso_cubes_tri_no_add +=
        num_triangulated_quad;
      triangulation_info.num_iso_cubes_tri_total +=
        num_triangulated_quad;
    }
    else if (quad_tri_method == SPLIT_MAX_ANGLE) {
      IJK::triangulate_quad_split_max_angle
        (dimension, dual_isosurface.vertex_coord,
         dual_isosurface.isopoly_vert,
         max_small_magnitude, dual_isosurface.simplex_vert,
         num_triangulated_quad);

      triangulation_info.num_iso_cubes_tri_no_add +=
        num_triangulated_quad;
      triangulation_info.num_iso_cubes_tri_total +=
        num_triangulated_quad;
    }
    else if (dualiso_data.flag_tri4_quad) {
      const VERTEX_INDEX num_iso_vert = dual_isosurface.NumIsoVert();
      const VERTEX_INDEX num_iso_poly = dual_isosurface.NumIsoPoly();

      if (num_iso_vert < num_iso_poly+dual_isosurface.FirstIsovDualToIsopoly()) {
        error.AddMessage
          ("Programming error.  Incorrect number of isosurface vertices.");
        error.AddMessage
          ("  Check that add_tri4_interior_vertices() is called before convert_quad_to_tri().");
        throw error;
      }

      if (!dual_isosurface.FlagIsovDualToIsopoly()) {
        error.AddMessage("Programming error. Missing isosurface vertices.");
        error.AddMessage
          ("  Missing isosurface vertices dual to isosurface quadrilaterals.");
        error.AddMessage
          ("  Call add_tri4_interior_vertices() before calling this routine.");
        throw error;
      }
      
      if (dual_isosurface.vertex_coord.size() > 0) {

        const ISO_VERTEX_INDEX first_isov_dual_to_iso_poly =
          dual_isosurface.FirstIsovDualToIsopoly();

        if (quad_tri_method == TRI4_ALL_QUADS) {
          IJK::triangulate_quad_list_from_interior_vertices
            (dual_isosurface.isopoly_vert, first_isov_dual_to_iso_poly, 
             dual_isosurface.simplex_vert, num_triangulated_quad);
          
          triangulation_info.num_iso_cubes_tri_add_interior1 +=
            num_triangulated_quad;
          triangulation_info.num_iso_cubes_tri_total +=
            num_triangulated_quad;
        }
        else {
          MESH_SIZE_TYPE num_tri2, num_tri4;
          
          triangulate_quad_tri2_or_tri4_max_min_angle
            (dualiso_data.ScalarGrid(), dual_isosurface.vertex_coord, 
             dual_isosurface.isopoly_vert, first_isov_dual_to_iso_poly,
             max_small_magnitude, 
             dual_isosurface.simplex_vert, num_tri2, num_tri4);

          triangulation_info.num_iso_cubes_tri_no_add +=
            num_tri2;
          triangulation_info.num_iso_cubes_tri_add_interior1 +=
            num_tri4;
          triangulation_info.num_iso_cubes_tri_total +=
            num_tri2 + num_tri4;
        }
      }
    }
    else {
      // Apply uniform triangulation.
      IJK::triangulate_quad_list_using_diagonal02
        (dual_isosurface.isopoly_vert, dual_isosurface.simplex_vert,
         num_triangulated_quad);

      triangulation_info.num_iso_cubes_tri_no_add +=
        num_triangulated_quad;
      triangulation_info.num_iso_cubes_tri_total +=
        num_triangulated_quad;
    }
  }

  //@}


  // ***********************************************************************
  //! @nam ENVELOPE BASED TRIANGULATION
  // ***********************************************************************

  //@{
  
  /*!
   *  @brief Triangulate isosurface quadrilaterals with some diagonal outside the envelope.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Each quadrilateral is dual to a grid edge.
   *  - Only use diagonals that are inside the envelope formed
   *    by the quad edges and the dual grid edge.
   *  - If no diagonals are inside the envelope, triangulate into four triangles.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   */
  template <typename DATA_TYPE, typename ISOPOLY_INFO_TYPE,
            typename ISOSURFACE_TYPE,
            typename NTYPES2, typename NTYPES4>
  void tri2_or_tri4_quads_with_diagonals_outside_envelope
  (const DATA_TYPE & dualiso_data,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface,
   NTYPES2 & num_tri2_split, NTYPES4 & num_tri4_split)
  {
    const int NUM_VERT_PER_QUAD(4);
    IJK::PROCEDURE_ERROR error("tri2_or_tri4_quads_with_diagonals_outside_envelope");

    if (!dual_isosurface.CheckNumVerticesPerIsoPoly
        (NUM_VERT_PER_QUAD, error)) { throw error; }

    if (!dual_isosurface.FlagIsovDualToIsopoly()) {
      error.AddMessage("Programming error. Missing isosurface vertices.");
      error.AddMessage
        ("  Missing isosurface vertices dual to isosurface quadrilaterals.");
      error.AddMessage
        ("  Call add_tri4_interior_vertices() before calling this routine.");
      throw error;
    }
        
    dual_isosurface.SetQuadVertexOrder(CIRCULAR_VERTEX_ORDER);

    // *** REORDER SO FirstIsovDualToIsopoly() is later in list.
    tri2_or_tri4_dual_quad_list_with_diagonals_outside_envelope
      (dualiso_data.ScalarGrid(), dual_isosurface.vertex_coord,
       dual_isosurface.FirstIsovDualToIsopoly(),
       dualiso_data.max_small_magnitude,
       dual_isosurface.isopoly_vert, isopoly_info,
       dual_isosurface.simplex_vert, num_tri2_split, num_tri4_split);
  }


  /*!
   *  @brief Triangulate isosurface quadrilaterals with some diagonal outside the envelope.
   *  - Choose triangulation that maximizes the minimum triangle angle.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If exactly one diagonal is outside the envelope 
   *    formed by the quad edges and the dual grid edge, 
   *    split the quad into two triangles using the other diagonal.
   *  - If both diagonals are outside the envelope,
   *    split the quad into four triangles.
   *  - Ignore quads that have both diagonals inside the envelope.
   *  - See \ref envelopeDoc "detailed description (file ijktri2D_envelope.tpp)"
   *    for envelope definition and references.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   *  @pre Vertices are in dimension 3.
   */
  template <typename DATA_TYPE, typename STYPE,
            typename ISOPOLY_INFO_TYPE,
            typename ISOSURFACE_TYPE,
            typename NTYPES2, typename NTYPES4>
  void tri2_or_tri4_quads_with_diagonals_outside_envelope_prefer_tri2
  (const DATA_TYPE & dualiso_data, const STYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface,
   NTYPES2 & num_tri2_split, NTYPES4 & num_tri4_split)
  {
    const int NUM_VERT_PER_QUAD(4);
    IJK::PROCEDURE_ERROR error
      ("tri2_or_tri4_quads_with_diagonals_outside_envelope_prefer_tri2");

    if (!dual_isosurface.CheckNumVerticesPerIsoPoly
        (NUM_VERT_PER_QUAD, error)) { throw error; }

    dual_isosurface.SetQuadVertexOrder(CIRCULAR_VERTEX_ORDER);

    if (dual_isosurface.FlagIsovDualToIsopoly()) {
      tri2_or_tri4_dual_quad_list_DoutsideE_prefer_tri2
        (dualiso_data.ScalarGrid(), dual_isosurface.vertex_coord,
         dual_isosurface.isopoly_vert, isopoly_info,
         dual_isosurface.FirstIsovDualToIsopoly(),
         dual_isosurface.simplex_vert, num_tri2_split, num_tri4_split);
    }
    else {
      ADD_ISOV_DUAL_TO_QUAD_SCALE
        <DUALISO_SCALAR_GRID_BASE,MESH_SIZE_TYPE>
        add_isov_dual_to_quad_scale(&(dualiso_data.ScalarGrid()));

      add_isov_dual_to_quad_scale.interior_vertex_param =
        dualiso_data.interior_vertex_param;

      tri2_or_tri4_dual_quad_list_DoutsideE_prefer_tri2_addV
        (dualiso_data.ScalarGrid(), 
         dual_isosurface.isopoly_vert, isopoly_info,
         isovalue, dual_isosurface.vertex_coord,
         dual_isosurface.simplex_vert, num_tri2_split, num_tri4_split,
         add_isov_dual_to_quad_scale);
    }
  }


  /*!
   *  @brief Split isosurface quadrilaterals into 4 triangles 
   *     if either quad diagonal is outside the envelope.
   */
  template <typename DATA_TYPE, typename STYPE,
            typename ISOPOLY_INFO_TYPE, typename ISOSURFACE_TYPE,
            typename NTYPE>
  void tri4_quads_with_diagonals_outside_envelope
  (const DATA_TYPE & dualiso_data, const STYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface, NTYPE & num_split)
  {
    const int NUM_VERT_PER_QUAD(4);
    IJK::PROCEDURE_ERROR error("tri4_quads_with_diagonals_outside_envelope");

    if (!dual_isosurface.CheckNumVerticesPerIsoPoly
        (NUM_VERT_PER_QUAD, error)) { throw error; }

    dual_isosurface.SetQuadVertexOrder(CIRCULAR_VERTEX_ORDER);

    if (dual_isosurface.FlagIsovDualToIsopoly()) {
      tri4_dual_quad_list_with_diagonals_outside_envelope
        (dualiso_data.ScalarGrid(), dual_isosurface.vertex_coord,
         dual_isosurface.isopoly_vert, isopoly_info,
         dual_isosurface.FirstIsovDualToIsopoly(),
         dual_isosurface.simplex_vert, num_split);
    }
    else {
      ADD_ISOV_DUAL_TO_QUAD_SCALE
        <DUALISO_SCALAR_GRID_BASE,MESH_SIZE_TYPE>
        add_isov_dual_to_quad_scale(&(dualiso_data.ScalarGrid()));

      add_isov_dual_to_quad_scale.interior_vertex_param =
        dualiso_data.interior_vertex_param;
      
      tri4_dual_quad_list_with_diagonals_outside_envelope_add_vertices
        (dualiso_data.ScalarGrid(), 
         dual_isosurface.isopoly_vert, isopoly_info,
         isovalue,
         dual_isosurface.vertex_coord,
         dual_isosurface.simplex_vert, num_split,
         add_isov_dual_to_quad_scale);      
    }
  }


  /*!
   *  @brief Split isosurface quadrilaterals with exactly one diagonal 
   *    outside the envelope.
   *  - Each quadrilateral is dual to a grid edge.
   *  - If one diagonal is inside the envelope formed by the quad edges
   *    and the dual grid edge and one is outside the envelope,
   *    split using the diagonal inside the envelope.
   *  - Ignore quads that have both diagonals inside the envelope or
   *    both diagonals outside the envelope.
   *  @pre compute_quad_diagonals_outside_envelope() has already
   *    been called to set quad_info[*].is_diagonal_in_envelope[].
   */
  template <typename DATA_TYPE, typename ISOPOLY_INFO_TYPE,
            typename ISOSURFACE_TYPE, typename NTYPE>
  void tri2_dual_quad_list_with_exactly_one_diagonal_outside_envelope
  (const DATA_TYPE & dualiso_data,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   ISOSURFACE_TYPE & dual_isosurface, NTYPE & num_split)
  {
    const int NUM_VERT_PER_QUAD(4);
    IJK::PROCEDURE_ERROR error
      ("tri2_quads_with_exactly_one_diagonal_inside_envelope");

    if (!dual_isosurface.CheckNumVerticesPerIsoPoly
        (NUM_VERT_PER_QUAD, error)) { throw error; }

    dual_isosurface.SetQuadVertexOrder(CIRCULAR_VERTEX_ORDER);

    tri2_dual_quad_list_with_exactly_one_diagonal_outside_envelope
      (dualiso_data.ScalarGrid(), dual_isosurface.vertex_coord,
       dual_isosurface.isopoly_vert, isopoly_info,
       dual_isosurface.simplex_vert, num_split);
  }

  //@}


  // ***********************************************************************
  // CLASS ADD_ISOV_DUAL_TO_QUAD_SCALE MEMBER FUNCTIONS
  // ***********************************************************************

  /// Add isosurface vertex dual to polytope ipoly. Scale coordinates by grid.Spacing().
  template <typename GRID_TYPE, typename NTYPEV>
  template <typename DTYPE, typename VTYPE,
            typename ISOP_INFO_TYPE, typename STYPE,
            typename CTYPE>
  VTYPE ADD_ISOV_DUAL_TO_QUAD_SCALE<GRID_TYPE,NTYPEV>::AddDualIsov
  (const DTYPE dimension, const VTYPE isopoly_vert[],
   const ISOP_INFO_TYPE & isopoly_info,
   const STYPE isovalue,
   std::vector<CTYPE> & vertex_coord) const
  {
    const VTYPE new_isov =
      add_tri4_interior_vertex_scale
      (*grid_ptr, isopoly_vert, isopoly_info, isovalue,
       interior_vertex_param, vertex_coord);

    return new_isov;
  }
  
}

#endif
