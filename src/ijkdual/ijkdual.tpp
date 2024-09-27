/*!
 *  @file ijkdual.tpp
 *  @brief Templates for dual contouring isosurface generation in arbitrary dimensions.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2006-2024 Rephael Wenger

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

/*!
  \mainpage IJKDUAL: IJK DUAL CONTOURING ALGORITHM

  IJKDUAL is a program for generating isosurfaces using dual contouring.
  The program computes isosurfaces from two, three or four dimensional
  volumetric grid data.
  Options include triangulating the quadrilateral mesh (3D only),
  and applying envelopes to guarantee no self intersections.
  (See \ref envelopeDoc for definition and references on envelopes.)

  The use of envelopes to guarantee no self intersections
  in 3D may not be necessary, since there are no known examples
  of self intersections generated using centroid positioning.
  However, there is currently no proof that centroid positioning
  guarantees no self intersections.

  With no self intersections and default -multi_isov option,
  the surface is guaranteed to be a manifold.
*/

#ifndef _IJKDUAL_TPP_
#define _IJKDUAL_TPP_

#include <string>

#include "ijk.tpp"
#include "ijkfacet_intersection_table.tpp"
#include "ijkmesh.tpp"
#include "ijktime.tpp"

#include "ijkdual_extract.tpp"
#include "ijkdual_position.tpp"

#include "ijkdual_types.h"
#include "ijkdual_datastruct.h"


/// ijkdual template functions
namespace IJKDUAL {

  // **************************************************
  // FORWARD DECLARATIONS
  // **************************************************

  template <typename GRID_TYPE, typename STYPE,
            typename VTYPE, typename CTYPE>
  void position_all_dual_isovertices_single_isov
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<VTYPE> & isov_list,
   const ISO_VERTEX_POSITION_PARAM & iso_vertex_position_param,
   std::vector<CTYPE> & vertex_coord,
   DUALISO_INFO & dualiso_info);


  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_multi_isov
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const ISO_VERTEX_POSITION_PARAM & iso_vertex_position_param,
   std::vector<CTYPE> & vertex_coord,
   DUALISO_INFO & dualiso_info);


  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_multi_isov
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const ISO_VERTEX_POSITION_PARAM & iso_vertex_position_param,
   std::vector<CTYPE> & vertex_coord,
   DUALISO_INFO & dualiso_info);


  // **************************************************
  //! @name DUAL CONTOURING ALGORITHM (HYPERCUBES)
  // **************************************************

  ///@{

  /*!
   *  @brief Dual Contouring Algorithm - single isosurface vertex
   *    per grid cube.
   *  @param dualiso_data Data for constructing dual isosurface.
   *    - Includes scalar grid and algorithm parameters.
   *  @param isovalue Isovalue.
   *  @param dual_isosurface Mesh storing the dual isosurface.
   *    - Includes array of isosurface vertex coordinates 
   *      and array of isosurface polytope vertices.
   *  @param[out] active_cube_list[kw] 
   *      Index of cube containing isosurface vertex kw.
   *    - Required information for computing envelopes.
   *  @param[out] isopoly_info[ipoly] Information on isosurface
   *    polytope kw.
   *    - Required information for triangulating isosurfaces.
   */
  template <typename DUAL_ISOSURFACE_TYPE,
            typename CUBE_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE>
  void dual_contouring_single_isov
  (const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue,
   DUAL_ISOSURFACE_TYPE & dual_isosurface,
   std::vector<CUBE_INDEX_TYPE> & active_cube_list,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   DUALISO_INFO & dualiso_info)
  {
    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const AXIS_SIZE_TYPE * axis_size = dualiso_data.ScalarGrid().AxisSize();
    IJK::PROCEDURE_ERROR error("dual_contouring_single_isov");

    std::clock_t t_start = std::clock();

    if (!dualiso_data.Check(error)) { throw error; };

    dual_isosurface.Clear();
    dualiso_info.time.Clear();
    
    ISO_MERGE_DATA merge_data(dimension, axis_size);

    dual_contouring_single_isov
      (dualiso_data.ScalarGrid(), isovalue,
       dualiso_data.iso_vertex_position_param,
       dual_isosurface.isopoly_vert, isopoly_info, active_cube_list,
       dual_isosurface.vertex_coord, merge_data, dualiso_info);

    // store times
    std::clock_t t_end = std::clock();
    IJK::clock2seconds(t_end-t_start, dualiso_info.time.total);    
  }


  /*!
   *  @brief Dual Contouring Algorithm - grid cubes may contain
   *    multiple isosurface vertices.
   *  @param dualiso_data Data for constructing dual isosurface.
   *    - Includes scalar grid and algorithm parameters.
   *  @param isovalue Isovalue.
   *  @param dual_isosurface Mesh storing the dual isosurface.
   *    - Includes array of isosurface vertex coordinates 
   *      and array of isosurface polytope vertices.
   *  @param[out] isov_list[kw] Information on isosurface vertex kw.
   *    - Required information for computing envelopes.
   *    - Must contain cube index, lookup table index,
   *      and isosurface patch index, and is usually derived
   *      from class DUAL_ISOVERT in ijkisopoly.tpp.
   *  @param[out] isopoly_info[ipoly] Information on isosurface
   *    polytope kw.
   *    - Required information for triangulating isosurfaces.
   */
  template <typename DUAL_ISOSURFACE_TYPE,
            typename DUAL_ISOVERT_TYPE,
            typename ISOPOLY_INFO_TYPE>
  void dual_contouring_multi_isov
  (const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue,
   DUAL_ISOSURFACE_TYPE & dual_isosurface,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   DUALISO_INFO & dualiso_info)
  {
    const int dimension = dualiso_data.ScalarGrid().Dimension();
    const AXIS_SIZE_TYPE * axis_size = dualiso_data.ScalarGrid().AxisSize();
    IJK::PROCEDURE_ERROR error("dual_contouring_multi_isov");

    std::clock_t t_start = std::clock();

    if (!dualiso_data.Check(error)) { throw error; };

    dual_isosurface.Clear();
    dualiso_info.time.Clear();
    
    ISO_MERGE_DATA merge_data(dimension, axis_size);

    dual_contouring_multi_isov
      (dualiso_data.ScalarGrid(), isovalue, dualiso_data,
       isov_list, dual_isosurface.isopoly_vert, isopoly_info,
       dual_isosurface.vertex_coord, dualiso_info);

    // store times
    std::clock_t t_end = std::clock();
    IJK::clock2seconds(t_end-t_start, dualiso_info.time.total);    
  }
  
  
  /*!
   *  @brief Dual Contouring Algorithm.
   *  @param dualiso_data Data for constructing dual isosurface.
   *    - Includes scalar grid and algorithm parameters.
   *  @param isovalue Isovalue.
   *  @param dual_isosurface Mesh storing the dual isosurface.
   *    - Includes array of isosurface vertex coordinates 
   *      and array of isosurface polytope vertices.
   *  @param[out] isov_list[kw] Information on isosurface vertex kw.
   *    - Required information for computing envelopes.
   *    - Could be just the index of the grid cube containing
   *      isosurface vertex kw.
   *    - If multiple isosurface vertices can be in a grid cube,
   *      then must contain cube index, lookup table index,
   *      and isosurface patch index, and is usually derived
   *      from class DUAL_ISOVERT in ijkisopoly.tpp.
   *  @param[out] isopoly_info[ipoly] Information on isosurface
   *    polytope kw.
   *    - Required information for triangulating isosurfaces.
   */
  template <typename DUAL_ISOSURFACE_TYPE,
            typename DUAL_ISOVERT_TYPE,
            typename ISOPOLY_INFO_TYPE>
  void dual_contouring
  (const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue,
   DUAL_ISOSURFACE_TYPE & dual_isosurface,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   DUALISO_INFO & dualiso_info)
  {
    const bool allow_multiple_isov =
      dualiso_data.AllowMultipleIsoVertices();
    
    if (allow_multiple_isov) {
      // Grid cubes may contain multiple isosurface vertices.
      dual_contouring_multi_isov
        (dualiso_data, isovalue, dual_isosurface,
         isov_list, isopoly_info, dualiso_info);
    }
    else {
      // Only one isosurface vertex per grid cube.
      dual_contouring_single_isov
        (dualiso_data, isovalue, dual_isosurface,
         isov_list, isopoly_info, dualiso_info);      
    }
  }

    
  /*!
   *  @overload
   *  @brief Dual Contouring Algorithm.
   ** - Version that does not return isov_list[] or isopoly_info[].
   */
  template <typename DUAL_ISOSURFACE_TYPE>
  void dual_contouring
  (const DUALISO_DATA & dualiso_data, const SCALAR_TYPE isovalue,
   DUAL_ISOSURFACE_TYPE & dual_isosurface, DUALISO_INFO & dualiso_info)
  {
    typedef typename DUALISO_DATA::SCALAR_GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename DUALISO_DATA::SCALAR_GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename IJK::GRID_EDGE<VTYPE,DTYPE> ISOPOLY_INFO_TYPE;

    const bool allow_multiple_isov =
      dualiso_data.AllowMultipleIsoVertices();
    std::vector<ISOPOLY_INFO_TYPE> isopoly_info;

    if (allow_multiple_isov) {
      // Grid cubes may contain multiple isosurface vertices.
      // Use version of isov_list storing cube index, tableindex
      //   and patch index.
      std::vector<DUAL_ISOVERT> isov_list;
      
      dual_contouring_multi_isov
        (dualiso_data, isovalue, dual_isosurface,
         isov_list, isopoly_info, dualiso_info);
    }
    else {
      // Only one isosurface vertex per grid cube.
      // Use version of isov_list storing only cube index.
      std::vector<VTYPE> active_cube_list;

      dual_contouring_single_isov
        (dualiso_data, isovalue, dual_isosurface,
         active_cube_list, isopoly_info, dualiso_info);      
    }
  }


  ///@}


  // **************************************************
  //! @name CONSTRUCT DUAL CONTOURING ISOSURFACE MESH
  // **************************************************

  ///@{

  /*!
   *  @brief Create isosurface vertices.
   *  - Create isosurface polytope vertices.
   *  - Allow multiple isosurface vertices per grid cube.
   *  @param isopoly_cube[i*numv_per_isopoly+k] Index of grid cube
   *    containing k'th vertex of i'th isosurface polygon.
   *  @param isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param[out] isopoly_vert[i*numv_per_isopoly+k]
   *    Index of k'th vertex of isosurface polytope i.
   *  @param[out] active_cube_list[] List of active cubes,
   *    i.e., cubes intersected by the isosurface.
   *  @param[out] isov_list[] List of isosurface vertices.
   *  @param[out] isopoly_cube_loc[i] Location of the cube isopoly_cube[i]
   *    in active_cube_index[].
   *  @param[out] merge_data Auxiliary data structure for merging
   *    identical elements in a list.
   *    - Used for merging identical cubes in isopoly_cube[].
   *  @param[out] dualiso_info Information collected in constructing
   *      dual contouring isosurface.
   */
  template <typename GRID_CUBE_DATA_TYPE, typename DUAL_ISOVERT_TYPE>
  void create_isosurface_vertices_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const std::vector<CUBE_INDEX> &  isopoly_cube,
   const std::vector<CUBE_EDGE_INDEX> & isopoly_cube_edge,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   std::vector<ISO_VERTEX_INDEX> & isopoly_cube_loc,
   MERGE_DATA & merge_data,
   DUALISO_INFO & dualiso_info)
  {
    const int dimension = scalar_grid.Dimension();
    const bool flag_split_non_manifold = param.SplitNonManifoldFlag();
    const bool flag_split_non_manifold_ambig =
      param.SplitNonManifoldAmbigFlag();
    const bool flag_select_split = param.SelectSplitFlag();
    const bool flag_separate_neg = param.SeparateNegFlag();
    const bool flag_connect_ambiguous = param.ConnectAmbiguousFlag();
    std::clock_t t0, t1, t2;
    IJK::PROCEDURE_ERROR error("create_isosurface_vertices_multi_isov");

    if (scalar_grid.Dimension() != isodual_table.Dimension()) {
      error.AddMessage
        ("Programming error.  Incorrect isodual table dimension.");
      error.AddMessage
        ("  Isodual table dimension does not match scalar grid dimension.");
      error.AddMessage
        ("    Isodual table dimension: ", isodual_table.Dimension(), "");
      error.AddMessage
        ("    Scalar grid dimension: ", scalar_grid.Dimension(), "");
      throw error;
    }

    t0 = clock();

    // active_cube_index[] is a list of active grid cube indices.
    std::vector<CUBE_INDEX> active_cube_index;

    merge_identical(isopoly_cube, active_cube_index, isopoly_cube_loc,
                    merge_data);
    t1 = clock();

    set_grid_cube_indices(active_cube_index, active_cube_list);

    VERTEX_INDEX num_split;
    if (flag_split_non_manifold || flag_select_split ||
        flag_connect_ambiguous || flag_split_non_manifold_ambig) {

      int num_non_manifold_split = 0;
      int num_1_2_changed = 0;
      int num_connect_changed = 0;
      int num_ambig_ridge_cubes_changed = 0;
      int num_non_ambig_ridge_cubes_changed = 0;

      IJK::compute_cube_isotable_info
        (scalar_grid, isodual_table, isovalue, active_cube_list);

      if (flag_split_non_manifold) {
        IJK::FACET_INTERSECTION_TABLE
          <DIRECTION_TYPE,VERTEX_INDEX,BOUNDARY_BITS_TYPE>
          facet_intersection_table(dimension);
        facet_intersection_table.Create();

        // Recompute isotable indices for cubes on grid ridges
        //   to avoid non-manifold vertices.
        IJK::recompute_cube_isotable_info_on_grid_ridges
          (scalar_grid, isodual_table, isovalue, facet_intersection_table,
           flag_separate_neg, active_cube_list, num_ambig_ridge_cubes_changed,
           num_non_ambig_ridge_cubes_changed);
      }

      if (flag_split_non_manifold || flag_split_non_manifold_ambig) {
        IJK::split_non_manifold_isov_pairs
          (scalar_grid, isodual_table, active_cube_list,
           num_non_manifold_split);
      }

      if (flag_select_split) {
        IJK::select_split_1_2_ambig
          (scalar_grid, isodual_table, isovalue, active_cube_list,
           num_1_2_changed);
      }

      if (flag_connect_ambiguous) {
        IJK::select_ambig_to_connect_isosurface
          (scalar_grid, isodual_table, isovalue, active_cube_list,
           num_connect_changed);
      }

      IJK::split_dual_isovert
        (isodual_table, isopoly_cube_loc, isopoly_cube_edge, active_cube_list,
         isov_list, isopoly_vert, num_split);

      dualiso_info.multi_isov.num_non_manifold_split = num_non_manifold_split;
      dualiso_info.multi_isov.num_1_2_changed = num_1_2_changed;
      dualiso_info.multi_isov.num_connect_changed = num_connect_changed;
      dualiso_info.multi_isov.num_ambig_ridge_cubes_changed =
        num_ambig_ridge_cubes_changed;
      dualiso_info.multi_isov.num_non_ambig_ridge_cubes_changed =
        num_non_ambig_ridge_cubes_changed;
    }
    else {
      IJK::split_dual_isovert
        (scalar_grid, isodual_table, isovalue,
         isopoly_cube_loc, isopoly_cube_edge, active_cube_list,
         isov_list, isopoly_vert, num_split);
    }

    t2 = clock();

    dualiso_info.scalar.num_non_empty_cubes = active_cube_list.size();
    dualiso_info.multi_isov.num_cubes_multi_isov = num_split;
    dualiso_info.multi_isov.num_cubes_single_isov =
      active_cube_list.size() - num_split;

    // store time
    IJK::clock2seconds(t1-t0, dualiso_info.time.merge);
    IJK::clock2seconds(t2-t1, dualiso_info.time.split_isov);
  }


  /*!
   *  @overload
   *  @brief Create isosurface vertices.
   *  - Version that does not return isopoly_cube_loc[].
   *  - Create isosurface polytope vertices.
   *  - Allow multiple isosurface vertices per grid cube.
   *  @param isopoly_cube[i*numv_per_isopoly+k] Index of grid cube
   *    containing k'th vertex of i'th isosurface polygon.
   *  @param isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param[out] isopoly_vert[i*numv_per_isopoly+k]
   *    Index of k'th vertex of isosurface polytope i.
   *  @param[out] active_cube_list[] List of active cubes,
   *    i.e., cubes intersected by the isosurface.
   *  @param[out] isov_list[] List of isosurface vertices.
   *  @param[out] merge_data Auxiliary data structure for merging
   *    identical elements in a list.
   *    - Used for merging identical cubes in isopoly_cube[].
   *  @param[out] dualiso_info Information collected in constructing
   *      dual contouring isosurface.
   */
  template <typename GRID_CUBE_DATA_TYPE, typename DUAL_ISOVERT_TYPE>
  void create_isosurface_vertices_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const std::vector<CUBE_INDEX> &  isopoly_cube,
   const std::vector<CUBE_EDGE_INDEX> & isopoly_cube_edge,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   MERGE_DATA & merge_data,
   DUALISO_INFO & dualiso_info)
  {
    // isopoly_cube_loc[i] = Location of the cube isopoly_cube[i]
    //   in active_cube_index[].
    std::vector<ISO_VERTEX_INDEX> isopoly_cube_loc;

    create_isosurface_vertices_multi_isov
      (scalar_grid, isovalue, isodual_table,
       isopoly_cube, isopoly_cube_edge, param,
       isopoly_vert, active_cube_list, isov_list,
       isopoly_cube_loc, merge_data, dualiso_info);
  }


  /*!
   *  @brief Construct isosurface mesh using Dual Contouring algorithm.
   *  - Allow multiple isosurface vertices per grid cube.
   *  - Returns list of isosurface polytope vertices.
   *  @param[out] isopoly_vert[ip*num_vert_per_poly+k]
   *    Index of k'th vertex of ip'th isosurface polytope.
   *  @param[out] isov_cube_info[i] Information about cube
   *    containing isopoly_vert[i].
   *  @param[out] active_cube_list[] List of active grid cubes.
   *  @param[out] isov_list[i] Information about isopoly_vert[i].
   */
  template <typename GRID_CUBE_DATA_TYPE, typename DUAL_ISOVERT_TYPE>
  void construct_multi_isov_mesh
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   MERGE_DATA & merge_data,
   DUALISO_INFO & dualiso_info)
  {
    std::clock_t t0, t1;

    t0 = clock();

    isopoly_vert.clear();
    dualiso_info.time.Clear();

    IJK::EXTRACT_DUAL_ISOPOLY_CUBE_AND_DUAL_EDGE
      <CUBE_INDEX,CUBE_EDGE_INDEX> extract_result;
    extract_dual_isopoly
      (scalar_grid, isovalue, extract_result, dualiso_info);

    t1 = clock();

    create_isosurface_vertices_multi_isov
      (scalar_grid, isovalue, isodual_table,
       extract_result.isopoly_cube, extract_result.isopoly_cube_edge,
       param, isopoly_vert, active_cube_list, isov_list,
       merge_data, dualiso_info);    

    // store times
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }


  /*!
   *  @brief Construct isosurface mesh using Dual Contouring algorithm.
   *  - Allow multiple isosurface vertices per grid cube.
   *  - Returns list of isosurface polytope vertices.
   *  - Return isosurface polytope information.
   *  @param[out] isopoly_vert[ip*num_vert_per_poly+k]
   *    Index of k'th vertex of ip'th isosurface polytope.
   *  @param isopoly_info[i] Isososurface polygon information,
   *    including dual grid edge.
   *  @param[out] isov_cube_info[i] Information about cube
   *    containing isopoly_vert[i].
   *  @param[out] active_cube_list[] List of active grid cubes.
   *  @param[out] isov_list[i] Information about isopoly_vert[i].
   */
  template <typename ISOPOLY_INFO_TYPE, typename GRID_CUBE_DATA_TYPE,
            typename DUAL_ISOVERT_TYPE>
  void construct_multi_isov_mesh_I
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   MERGE_DATA & merge_data,
   DUALISO_INFO & dualiso_info)
  {
    std::clock_t t0, t1;

    t0 = clock();

    isopoly_vert.clear();
    dualiso_info.time.Clear();

    std::vector<CUBE_INDEX> isopoly_cube;
    std::vector<CUBE_EDGE_INDEX> isopoly_cube_edge;
    extract_dual_isopoly
      (scalar_grid, isovalue, isopoly_cube, isopoly_cube_edge,
       isopoly_info, dualiso_info);

    t1 = clock();

    create_isosurface_vertices_multi_isov
      (scalar_grid, isovalue, isodual_table,
       isopoly_cube, isopoly_cube_edge,
       param, isopoly_vert, active_cube_list, isov_list,
       merge_data, dualiso_info);    
    
    // store times
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }

  ///@}


  // *******************************************************************
  //! @name DUAL CONTOURING (HYPERCUBES) MULTI ISOSURFACE VERTICES
  // *******************************************************************

  ///@{

  /*!
   *  @brief Extract isosurface using Dual Contouring algorithm.
   *  - Allow multiple isosurface vertices per grid cube.
   *  - Returns list of isosurface polytope vertices
   *    and list of isosurface vertex coordinates.
   *  - Return isosurface polytope information.
   *  @tparam GRID_CUBE_DATA_TYPE Type containing grid cube information.
   *    - Should be derived from class GRID_CUBE_ISOVERT.
   *  @tparam DUAL_ISOVERT_TYPE Type containing dual isosurface vertex.
   *    - Should be derived from template class DUAL_ISOVERT.
   *  @param[out] isopoly_vert[] Array of isosurface polygon vertices.
   *    - isopoly_vert[numv_per_poly*ip+k] =
   *      k'th vertex of isosurface polygon ip.
   *  @param[out] active_cube_list[] Information about active grid cubes.
   *  @param[out] isov_list[i] Information about isosurface vertex i.
   *    - Includes grid index of cube containing isosurface vertex i,
   *      isosurface patch containing isosurface vertex i, and
   *      isosurface lookup table index.
   */
  template <typename ISOPOLY_INFO_TYPE, typename GRID_CUBE_DATA_TYPE,
            typename DUAL_ISOVERT_TYPE>
  void dual_contouring_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data,
   DUALISO_INFO & dualiso_info)
  {
    IJK::PROCEDURE_ERROR error("dual_contouring_multi_isov");

    const std::clock_t t0 = clock();

    construct_multi_isov_mesh_I
      (scalar_grid, isovalue, isodual_table, param, isopoly_vert,
       isopoly_info, active_cube_list, isov_list, merge_data, dualiso_info);

    const std::clock_t t1 = clock();

    position_all_dual_isovertices_multi_isov
      (scalar_grid, isodual_table, isovalue,
       isov_list, isopoly_vert, isopoly_info,
       param.iso_vertex_position_param, vertex_coord, dualiso_info);

    const std::clock_t t2 = clock();

    // store times
    IJK::clock2seconds(t2-t1, dualiso_info.time.position.total);
    IJK::clock2seconds(t2-t0, dualiso_info.time.total);
  }


  /*!
   *  @overload
   *  @brief Extract isosurface using Dual Contouring algorithm.
   *  - Allow multiple isosurface vertices per grid cube.
   *  - Return isosurface vertex and isosurface polytope information.
   *  - Version that creates isodual_table and merge_data.
   *  - Version that does not return active_cube_list.
   */
  template <typename ISOPOLY_INFO_TYPE, typename DUAL_ISOVERT_TYPE>
  void dual_contouring_multi_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const DUALISO_DATA_PARAM & param,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   COORD_ARRAY & vertex_coord,
   DUALISO_INFO & dualiso_info)
  {
    const int dimension = scalar_grid.Dimension();
    const AXIS_SIZE_TYPE * axis_size = scalar_grid.AxisSize();
    const bool flag_separate_neg = param.SeparateNegFlag();
    const bool flag_always_separate_opposite(true);
    ISODUAL_CUBE_TABLE_AMBIG
      isodual_table(dimension, flag_separate_neg,
                    flag_always_separate_opposite);
    ISO_MERGE_DATA merge_data(dimension, axis_size);
    std::vector<GRID_CUBE_DATA> active_cube_list;

    dual_contouring_multi_isov
      (scalar_grid, isovalue, isodual_table, param, isopoly_vert,
       isopoly_info, active_cube_list, isov_list, vertex_coord,
       merge_data, dualiso_info);
  }

  ///@}


  // *******************************************************************
  //! @name DUAL CONTOURING MANIFOLD
  // *******************************************************************

  ///@{

  /*!
   *  @brief Extract isosurface using Manifold Dual Contouring algorithm.
   *  - Isosurface must be 2D surface in 3D volume.
   *  - Allow multiple isosurface vertices per grid cube.
   *  - Returns list of isosurface polytope vertices
   *    and list of isosurface vertex coordinates.
   *  - Return isosurface polytope information.
   *  - Sets algorithm flags to ensure that isosurface mesh
   *    has a manifold representation.
   */
  template <typename ISOPOLY_INFO_TYPE, typename GRID_CUBE_DATA_TYPE,
            typename DUAL_ISOVERT_TYPE>
  void dual_contouring_manifold_I
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<GRID_CUBE_DATA_TYPE> & active_cube_list,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data,
   DUALISO_INFO & dualiso_info)
  {
    DUALISO_DATA_PARAM manifold_param = param;

    manifold_param.allow_multiple_iso_vertices = true;
    manifold_param.flag_split_non_manifold = true;

    dual_contouring_multi_isov
      (scalar_grid, isovalue, isodual_table, manifold_param,
       isopoly_vert, isopoly_info, active_cube_list, isov_list,
       vertex_coord, merge_data, dualiso_info);
  }


  /*!
   *  @brief Extract isosurface using Manifold Dual Contouring algorithm.
   *  - Isosurface must be 2D surface in 3D volume.
   *  - Allow multiple isosurface vertices per grid cube.
   *  - Returns list of isosurface polytope vertices
   *    and list of isosurface vertex coordinates.
   *  - Return isosurface polytope information.
   *  - Sets algorithm flags to ensure that isosurface mesh
   *    has a manifold representation.
   *  - Version without active_cube_list.
   */
  template <typename ISOPOLY_INFO_TYPE, typename DUAL_ISOVERT_TYPE>
  void dual_contouring_manifold_I
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISODUAL_CUBE_TABLE_AMBIG & isodual_table,
   const DUALISO_DATA_PARAM & param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<DUAL_ISOVERT_TYPE> & isov_list,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data,
   DUALISO_INFO & dualiso_info)
  {
    DUALISO_DATA_PARAM manifold_param = param;

    manifold_param.allow_multiple_iso_vertices = true;
    manifold_param.flag_split_non_manifold = true;

    dual_contouring_multi_isov
      (scalar_grid, isovalue, isodual_table, manifold_param,
       isopoly_vert, isopoly_info, isov_list,
       vertex_coord, merge_data, dualiso_info);
  }

  ///@}


  // *******************************************************************
  //! @name DUAL CONTOURING (HYPERCUBES) SINGLE ISOSURFACE VERTICES
  // *******************************************************************

  ///@{

  /*!
   *  @brief Extract isosurface using Dual Contouring algorithm.
   *  - Single isosurface vertex per grid cube.
   *  - Returns list of isosurface polytope vertices
   *   and list of isosurface vertex coordinates.
   *  - Return isosurface polytope information.
   *  @tparam GRID_CUBE_DATA_TYPE Type containing grid cube information.
   *    - Should be derived from class GRID_CUBE_ISOVERT.
   *  @tparam DUAL_ISOVERT_TYPE Type containing dual isosurface vertex.
   *    - Should be derived from template class DUAL_ISOVERT.
   *  @param[out] isopoly_vert[] Array of isosurface polygon vertices.
   *    - isopoly_vert[numv_per_poly*ip+k] =
   *      k'th vertex of isosurface polygon ip.
   *  @param[out] active_cube_list[] Information about active grid cubes.
   */
  template <typename ISOPOLY_INFO_TYPE, typename CUBE_INDEX_TYPE>
  void dual_contouring_single_isov
  (const DUALISO_SCALAR_GRID_BASE & scalar_grid,
   const SCALAR_TYPE isovalue,
   const ISO_VERTEX_POSITION_PARAM & iso_vertex_position_param,
   std::vector<ISO_VERTEX_INDEX> & isopoly_vert,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   std::vector<CUBE_INDEX_TYPE> & active_cube_list,
   COORD_ARRAY & vertex_coord,
   MERGE_DATA & merge_data,
   DUALISO_INFO & dualiso_info)
  {
    std::clock_t t0, t1, t2, t3;

    t0 = clock();

    isopoly_vert.clear();
    vertex_coord.clear();
    dualiso_info.time.Clear();

    std::vector<CUBE_INDEX> isopoly_cube;
    extract_dual_isopoly
      (scalar_grid, isovalue, isopoly_cube, isopoly_info,
       dualiso_info);
    
    t1 = clock();

    merge_identical(isopoly_cube, active_cube_list,
                    isopoly_vert, merge_data);    
    t2 = clock();

    position_all_dual_isovertices_single_isov
      (scalar_grid, isovalue, active_cube_list,
       isopoly_vert, isopoly_info,
       iso_vertex_position_param, vertex_coord, dualiso_info);

    t3 = clock();

    dualiso_info.scalar.num_non_empty_cubes = active_cube_list.size();

    // store times
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
    IJK::clock2seconds(t2-t1, dualiso_info.time.merge);
    IJK::clock2seconds(t3-t2, dualiso_info.time.position.total);
    IJK::clock2seconds(t3-t0, dualiso_info.time.total);
  }

  ///@}


  // ******************************************************************
  //! @name Position isosurface vertices.
  // ******************************************************************

  ///@{

  /*!
   *  @brief Position dual isosurface vertices.
   *  @pre Each grid cube contains at most one isosurface vertex.
   *  @param active_cube_list[kw] Cube containing isosurface vertex kw.
   *    - One isosurface vertex per cube.
   */
   template <typename GRID_TYPE, typename STYPE,
             typename ITYPE, typename CTYPE>
   void position_all_dual_isovertices_single_isov
   (const GRID_TYPE & scalar_grid,
    const STYPE isovalue,
    const std::vector<ITYPE> & active_cube_list,
    const ISO_VERTEX_POSITION_PARAM & iso_vertex_position_param,
    std::vector<CTYPE> & vertex_coord,
    DUALISO_INFO & dualiso_info)
   {
     const VERTEX_POSITION_METHOD vertex_position_method =
       iso_vertex_position_param.vertex_position_method;
     const bool flag_move_all_isov_away_from_cube_facets =
       iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets;     
     const COORD_TYPE position_offset =
       iso_vertex_position_param.position_offset;
     const RANDOM_SEED_TYPE random_seed =
       iso_vertex_position_param.random_pos_param.random_seed;
     const RANDOM_DISTRIBUTION random_pos_distribution =
       iso_vertex_position_param.random_pos_param.distribution;
     IJK::PROCEDURE_ERROR error("position_all_dual_isovertices_single_isov");

     std::clock_t t0, t1;

     t0 = clock();

     if (vertex_position_method == CUBE_CENTER) {
       position_all_dual_isovertices_cube_center
         (scalar_grid, active_cube_list, vertex_coord);
     }

 #ifdef INCLUDE_RANDOM

     else if (vertex_position_method == RANDOM_ISOV_POS) {
       COORD_TYPE random_offset = position_offset;
       if (!iso_vertex_position_param.random_pos_param.flag_gen_isov_apply_offset)
         { random_offset = 0; }

       if (random_pos_distribution == UNIFORM_DISTRIBUTION) {
         position_all_dual_isovertices_random_uniform
           (scalar_grid, active_cube_list, random_offset, random_seed,
            vertex_coord);
       }
       else if (random_pos_distribution == U_QUADRATIC_DISTRIBUTION) {
         position_all_dual_isovertices_random_U_quadratic
           (scalar_grid, active_cube_list, random_offset, random_seed,
            vertex_coord);
       }
       else {
         error.AddMessage
           ("Programming error. Random position distribution not set.");
         throw error;
       }
     }

 # endif

     else {
       // default
       position_all_dual_isovertices_centroid_simple
         (scalar_grid, isovalue, active_cube_list, vertex_coord);
     }

     t1 = clock();
     IJK::clock2seconds(t1-t0, dualiso_info.time.position.basic);

     if (flag_move_all_isov_away_from_cube_facets) {
       reposition_all_dual_isovert_away_from_cube_facets
         (scalar_grid, active_cube_list, position_offset, vertex_coord,
          dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets);

       std::clock_t t2 = clock();
       IJK::clock2seconds
         (t2-t1, dualiso_info.time.position.move_isov_away_from_cube_faces);
     }

     // *** IMPLEMENT separate_isov_near_shared_grid_cube_facets...
   }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices. 
   *    (Return isopoly_vert[] and isopoly_info[].)
   *  - Version returning isopoly_vert[] and isopoly_info[].
   *  @pre Each grid cube contains at most one isosurface vertex.
   */
  template <typename GRID_TYPE, typename STYPE,
            typename ITYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE>
  void position_all_dual_isovertices_single_isov
  (const GRID_TYPE & scalar_grid,
   const STYPE isovalue,
   const std::vector<ITYPE> & active_cube_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const ISO_VERTEX_POSITION_PARAM & iso_vertex_position_param,
   std::vector<CTYPE> & vertex_coord,
   DUALISO_INFO & dualiso_info)
  {
    const VERTEX_POSITION_METHOD vertex_position_method =
      iso_vertex_position_param.vertex_position_method;
    const bool flag_move_all_isov_away_from_cube_facets =
      iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets;         
    const COORD_TYPE position_offset =
      iso_vertex_position_param.position_offset;
    const RANDOM_SEED_TYPE random_seed =
      iso_vertex_position_param.random_pos_param.random_seed;
    const RANDOM_DISTRIBUTION random_pos_distribution =
      iso_vertex_position_param.random_pos_param.distribution;
    IJK::PROCEDURE_ERROR error("position_all_dual_isovertices_single_isov");

    std::clock_t t0, t1;

    t0 = clock();

    if (vertex_position_method == CUBE_CENTER) {
      position_all_dual_isovertices_cube_center
        (scalar_grid, active_cube_list, vertex_coord);
    }

#ifdef INCLUDE_RANDOM

    else if (vertex_position_method == RANDOM_ISOV_POS) {
      COORD_TYPE random_offset = position_offset;
      if (!iso_vertex_position_param.random_pos_param.flag_gen_isov_apply_offset)
        { random_offset = 0; }

      if (random_pos_distribution == UNIFORM_DISTRIBUTION) {
        position_all_dual_isovertices_random_uniform
          (scalar_grid, active_cube_list, random_offset, random_seed,
           vertex_coord);
      }
      else if (random_pos_distribution == U_QUADRATIC_DISTRIBUTION) {
        position_all_dual_isovertices_random_U_quadratic
          (scalar_grid, active_cube_list, random_offset, random_seed,
           vertex_coord);
      }
      else {
        error.AddMessage
          ("Programming error. Random position distribution not set.");
        throw error;
      }
    }

# endif

    else {
      // default
      position_all_dual_isovertices_centroid
        (scalar_grid, isovalue, active_cube_list,
         isopoly_vert, isopoly_info, vertex_coord);
    }

    t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.position.basic);

    if (flag_move_all_isov_away_from_cube_facets) {
      reposition_all_dual_isovert_away_from_cube_facets
        (scalar_grid, active_cube_list, position_offset, vertex_coord,
         dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets);

      std::clock_t t2 = clock();
      IJK::clock2seconds
        (t2-t1, dualiso_info.time.position.move_isov_away_from_cube_faces);
    }
  }


  /*!
   *  @brief Position dual isosurface vertices.
   *  - Allows more than one isosurface vertex per cube.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_multi_isov
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const ISO_VERTEX_POSITION_PARAM & iso_vertex_position_param,
   std::vector<CTYPE> & vertex_coord,
   DUALISO_INFO & dualiso_info)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const VERTEX_POSITION_METHOD vertex_position_method =
      iso_vertex_position_param.vertex_position_method;
    const bool flag_move_all_isov_away_from_cube_facets =
      iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets;
    const bool flag_move_all_isov_away_from_cube_ridges =
      iso_vertex_position_param.flag_move_all_isov_away_from_cube_ridges;    
    const bool flag_separate_isov_near_cube_facets =
      iso_vertex_position_param.flag_separate_isov_near_cube_facets;
    const bool flag_separate_isov_by_cube_center =
      iso_vertex_position_param.flag_separate_isov_by_cube_center;
    const COORD_TYPE min_distance =
       iso_vertex_position_param.position_offset;
    const COORD_TYPE position_offset =
      iso_vertex_position_param.position_offset;
    IJK::PROCEDURE_ERROR error("position_all_dual_isovertices_multi_isov");

    // min_num_bins for radix_sort().
    // *** WHAT SHOULD THIS BE? ***
    const NUMBER_TYPE min_num_bins = 100;

    std::clock_t t0, t1;

    t0 = clock();

    if (vertex_position_method == CUBE_CENTER) {
      position_all_dual_isovertices_near_cube_center_multi
      (scalar_grid, isodual_table, isov_list, position_offset,
       vertex_coord);
    }
    else if (vertex_position_method == IVOL_LIFTED02) {
      SCALAR_TYPE isovalue0 = isovalue/2;
      SCALAR_TYPE isovalue1 = (isovalue+2)/2;
      position_all_dual_isovertices_ivol_lifted
      (scalar_grid, isodual_table, isovalue0, isovalue1,
       isov_list, vertex_coord);
    }

 #ifdef INCLUDE_RANDOM

    else if (vertex_position_method == RANDOM_ISOV_POS) {
      DUAL_TABLE_VERTEX_CONNECTIVITY
        isodual_table_vertex_connectivity(isodual_table);

      compute_dual_cube_isotable_vertex_connectivity
        (isodual_table, isodual_table_vertex_connectivity);

      position_all_dual_isovertices_multi_random
        (scalar_grid, isodual_table,
         isodual_table_vertex_connectivity,
         isopoly_vert, isov_list,
         iso_vertex_position_param.random_pos_param,
         vertex_coord);
    }

#endif

    else {
      position_all_dual_isovertices_centroid_multi_simple
        (scalar_grid, isodual_table, isovalue, isov_list, vertex_coord);
    }

    t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.position.basic);

    if (flag_move_all_isov_away_from_cube_facets) {
      const std::clock_t t0 = clock();

      reposition_all_dual_isovert_away_from_cube_facets
        (scalar_grid, isov_list, position_offset, vertex_coord,
         dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets);

      const std::clock_t t1 = clock();
      IJK::clock2seconds
        (t1-t0, dualiso_info.time.position.move_isov_away_from_cube_faces);
    }
    else if (flag_separate_isov_near_cube_facets) {
      const std::clock_t t0 = clock();

      separate_isov_near_shared_grid_cube_facets
        (scalar_grid, isov_list, min_distance, position_offset,
         min_num_bins, vertex_coord,
         dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets);

      const std::clock_t t1 = clock();
      IJK::clock2seconds
        (t1-t0, dualiso_info.time.position.separate_isov_near_cube_facets);
    }


    if (flag_move_all_isov_away_from_cube_ridges) {
      const std::clock_t t0 = clock();

      int num_moved = 0;
      
      reposition_all_dual_isovert_away_from_cube_ridges
        (scalar_grid, isov_list, position_offset, vertex_coord, num_moved);
 

      dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets += num_moved;
      
      const std::clock_t t1 = clock();
      IJK::clock2seconds
        (t1-t0, dualiso_info.time.position.move_isov_away_from_cube_faces);      
    }
    
    if (flag_separate_isov_by_cube_center) {
      const std::clock_t t0 = clock();

      DUAL_TABLE_VERTEX_CONNECTIVITY
        isodual_table_vertex_connectivity(isodual_table);
      compute_dual_cube_isotable_vertex_connectivity
        (isodual_table, isodual_table_vertex_connectivity);

      separate_all_dual_isov_by_cube_centers
        (scalar_grid, isodual_table,
         isodual_table_vertex_connectivity,
         isopoly_vert, isov_list, position_offset, vertex_coord);

      const std::clock_t t1 = clock();
      IJK::clock2seconds
        (t1-t0, dualiso_info.time.position.separate_by_cube_center);
    }
  }


  /*!
   *  @overload
   *  @brief Position dual isosurface vertices. (Use isopoly_info[].)
   *  - Allows more than one isosurface vertex per cube.
   *  - Version that uses array isopoly_info[].
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE_TYPE,
            typename STYPE, typename ISOV_INDEX_TYPE,
            typename ISOPOLY_INFO_TYPE,
            typename DUAL_ISOV_TYPE, typename CTYPE>
  void position_all_dual_isovertices_multi_isov
  (const GRID_TYPE & scalar_grid,
   const ISODUAL_TABLE_TYPE & isodual_table,
   const STYPE isovalue,
   const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const std::vector<ISOV_INDEX_TYPE> & isopoly_vert,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const ISO_VERTEX_POSITION_PARAM & iso_vertex_position_param,
   std::vector<CTYPE> & vertex_coord,
   DUALISO_INFO & dualiso_info)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUMBER_TYPE;

    const VERTEX_POSITION_METHOD vertex_position_method =
      iso_vertex_position_param.vertex_position_method;
    const bool flag_move_all_isov_away_from_cube_facets =
      iso_vertex_position_param.flag_move_all_isov_away_from_cube_facets;    
    const bool flag_move_all_isov_away_from_cube_ridges =
      iso_vertex_position_param.flag_move_all_isov_away_from_cube_ridges;
    const bool flag_separate_isov_near_cube_facets =
      iso_vertex_position_param.flag_separate_isov_near_cube_facets;
    const bool flag_separate_isov_by_cube_center =
      iso_vertex_position_param.flag_separate_isov_by_cube_center;
    const COORD_TYPE min_distance =
       iso_vertex_position_param.position_offset;
    const COORD_TYPE position_offset =
      iso_vertex_position_param.position_offset;
    IJK::PROCEDURE_ERROR error("position_all_dual_isovertices_multi_isov");

    // min_num_bins for radix_sort().
    // *** WHAT SHOULD THIS BE? ***
    const NUMBER_TYPE min_num_bins = 100;

    std::clock_t t0, t1;

    t0 = clock();

    if (vertex_position_method == CUBE_CENTER) {
      position_all_dual_isovertices_near_cube_center_multi
      (scalar_grid, isodual_table, isov_list, position_offset,
       vertex_coord);
    }
    else if (vertex_position_method == IVOL_LIFTED02) {
      SCALAR_TYPE isovalue0 = isovalue/2;
      SCALAR_TYPE isovalue1 = (isovalue+2)/2;
      position_all_dual_isovertices_ivol_lifted
      (scalar_grid, isodual_table, isovalue0, isovalue1,
       isov_list, vertex_coord);
    }

 #ifdef INCLUDE_RANDOM

    else if (vertex_position_method == RANDOM_ISOV_POS) {
      DUAL_TABLE_VERTEX_CONNECTIVITY
        isodual_table_vertex_connectivity(isodual_table);

      compute_dual_cube_isotable_vertex_connectivity
        (isodual_table, isodual_table_vertex_connectivity);

      position_all_dual_isovertices_multi_random
        (scalar_grid, isodual_table,
         isodual_table_vertex_connectivity,
         isopoly_vert, isov_list,
         iso_vertex_position_param.random_pos_param,
         vertex_coord);
    }

#endif

    else {
      position_all_dual_isovertices_centroid_multi
        (scalar_grid, isodual_table, isovalue,
         isov_list, isopoly_vert, isopoly_info, vertex_coord);
    }

    t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.position.basic);

    if (flag_move_all_isov_away_from_cube_facets) {    
      const std::clock_t t0 = clock();

      reposition_all_dual_isovert_away_from_cube_facets
        (scalar_grid, isov_list, position_offset, vertex_coord,
         dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets);

      const std::clock_t t1 = clock();
      IJK::clock2seconds
        (t1-t0, dualiso_info.time.position.move_isov_away_from_cube_faces);
    }
    else if (flag_separate_isov_near_cube_facets) {
      const std::clock_t t0 = clock();

      separate_isov_near_shared_grid_cube_facets
        (scalar_grid, isov_list, min_distance, position_offset,
         min_num_bins, vertex_coord,
         dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets);

      const std::clock_t t1 = clock();
      IJK::clock2seconds
        (t1-t0, dualiso_info.time.position.separate_isov_near_cube_facets);
    }

    if (flag_move_all_isov_away_from_cube_ridges) {
      const std::clock_t t0 = clock();
      int num_moved = 0;

      reposition_all_dual_isovert_away_from_cube_ridges
        (scalar_grid, isov_list, position_offset, vertex_coord, num_moved);

      dualiso_info.isov.num_times_isov_moved_away_from_grid_cube_facets += num_moved;
      
      const std::clock_t t1 = clock();
      IJK::clock2seconds
        (t1-t0, dualiso_info.time.position.move_isov_away_from_cube_faces);      
    }
    
    if (flag_separate_isov_by_cube_center) {
      const std::clock_t t0 = clock();

      DUAL_TABLE_VERTEX_CONNECTIVITY
        isodual_table_vertex_connectivity(isodual_table);
      compute_dual_cube_isotable_vertex_connectivity
        (isodual_table, isodual_table_vertex_connectivity);

      separate_all_dual_isov_by_cube_centers
        (scalar_grid, isodual_table,
         isodual_table_vertex_connectivity,
         isopoly_vert, isov_list, position_offset, vertex_coord);

      const std::clock_t t1 = clock();
      IJK::clock2seconds
        (t1-t0, dualiso_info.time.position.separate_by_cube_center);
    }
  }

  ///@}

}

#endif
