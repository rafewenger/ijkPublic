/*!
 *  @file ijkdual_extract.tpp
 *  @brief Template functions for extracting dual isosurface mesh
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


#ifndef _IJKDUAL_EXTRACT_TPP_
#define _IJKDUAL_EXTRACT_TPP_

#include "ijkgrid_macros.h"
#include "ijkdual_types.h"

#include "ijkinterpolate.tpp"
#include "ijkisopoly.tpp"
#include "ijkmerge.tpp"
#include "ijktime.tpp"


namespace IJKDUAL {

  // ***************************************************
  // EXTRACT ROUTINES
  // ***************************************************

  /*!
   *  @brief Extract isosurface polytopes
   *  - Returns list representing isosurface polytopes.
   *  @param scalar_grid Scalar grid data.
   *  @param isovalue Isosurface scalar value.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polygon.
   */
  template <typename GTYPE, typename STYPE,
            typename EXTRACT_TYPE,
            typename INFO_TYPE>
  void extract_dual_isopoly
  (const GTYPE & scalar_grid, const STYPE isovalue,
   EXTRACT_TYPE & extract_result,
   INFO_TYPE & dualiso_info)
  {
    typedef typename GTYPE::VERTEX_INDEX_TYPE VTYPE;
    
    const int dimension = scalar_grid.Dimension();
    dualiso_info.time.extract = 0;

    if (dimension < 1) { return; }
    
    clock_t t0 = clock();

    // initialize output
    extract_result.Clear();

    if (scalar_grid.NumCubeVertices() < 1) { return; }
    IJK::FACET_INTERIOR_VERTEX_LIST<VERTEX_INDEX>
      facet_vertex_list(scalar_grid, 0, true);

    for (int edge_dir = 0; edge_dir < dimension; edge_dir++) {

      if (scalar_grid.AxisSize(edge_dir) < 1) {
        // No interior edges.
        return;
      }

      const int num_edges = scalar_grid.AxisSize(edge_dir)-1;

      facet_vertex_list.GetVertices(scalar_grid, edge_dir);

      if (edge_dir+1 < dimension) {
        // For each vertex in facet_vertex_list,
        //   process each edge in direction edge_dir.
        for (VTYPE i = 0; i < facet_vertex_list.NumVertices(); i++) {
          VTYPE iend0 = facet_vertex_list.VertexIndex(i);

          for (int j = 0; j < num_edges; j++) {
            const VTYPE iend1 = scalar_grid.NextVertex(iend0, edge_dir);
          
            if ((scalar_grid.Scalar(iend0) < isovalue) !=
                (scalar_grid.Scalar(iend1) < isovalue)) {
              extract_result.extract_dual_isopoly
                (scalar_grid, isovalue, iend0, edge_dir);
            }
            
            iend0 = iend1;
          }
        }
      }
      else {
        // For each edge on axis in direction edge_dir,
        //   process each vertex in facet_vertex_list.
        // Change order to processes vertices/edges
        //   that are close together in memory and significantly
        //   speed up this routine.

        // Add increment0 to each vertex in facet_vertex_list to get iend0.
        // Add increment to iend0 to get iend1.
        VTYPE increment0 = 0;
        VTYPE increment1 = scalar_grid.AxisIncrement(edge_dir);
        for (int j = 0; j < num_edges; j++) {
          for (VTYPE i = 0; i < facet_vertex_list.NumVertices(); i++) {
            const VTYPE iend0 = facet_vertex_list.VertexIndex(i) + increment0;
            const VTYPE iend1 = iend0 + increment1;
            if ((scalar_grid.Scalar(iend0) < isovalue) !=
                (scalar_grid.Scalar(iend1) < isovalue)) {
              extract_result.extract_dual_isopoly
                (scalar_grid, isovalue, iend0, edge_dir);
            }
          }
          increment0 = increment0 + increment1;
        }
      }
    }

    clock_t t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }


  /*!
   *  @brief Extract isosurface polytopes
   *  - Returns list representing isosurface polytopes.
   *  @param scalar_grid Scalar grid data.
   *  @param isovalue Isosurface scalar value.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polygon.
   *  @param[out] isopoly_info[i] Isososurface polygon information,
   *    including dual grid edge.
   */
  template <typename GTYPE, typename STYPE,
            typename ICUBE_TYPE, typename ISOPOLY_INFO_TYPE,
            typename INFO_TYPE>
  void extract_dual_isopoly
  (const GTYPE & scalar_grid, const STYPE isovalue,
   std::vector<ICUBE_TYPE> & isopoly_cube,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   INFO_TYPE & dualiso_info)
  {
    typedef typename GTYPE::VERTEX_INDEX_TYPE VTYPE;
    
    const int dimension = scalar_grid.Dimension();
    dualiso_info.time.extract = 0;

    if (dimension < 1) { return; }
    
    clock_t t0 = clock();

    // initialize output
    isopoly_cube.clear();
    isopoly_info.clear();

    if (scalar_grid.NumCubeVertices() < 1) { return; }
    IJK::FACET_INTERIOR_VERTEX_LIST<VERTEX_INDEX>
      facet_vertex_list(scalar_grid, 0, true);

    for (int edge_dir = 0; edge_dir < dimension; edge_dir++) {

      if (scalar_grid.AxisSize(edge_dir) < 1) {
        // No interior edges.
        return;
      }

      const int num_edges = scalar_grid.AxisSize(edge_dir)-1;

      facet_vertex_list.GetVertices(scalar_grid, edge_dir);

      if (edge_dir+1 < dimension) {
        // For each vertex in facet_vertex_list,
        //   process each edge in direction edge_dir.
        for (VTYPE i = 0; i < facet_vertex_list.NumVertices(); i++) {
          VTYPE iend0 = facet_vertex_list.VertexIndex(i);

          for (int j = 0; j < num_edges; j++) {
            extract_dual_isopoly_around_bipolar_edge_I
              (scalar_grid, isovalue, iend0, edge_dir,
               isopoly_cube, isopoly_info);
            iend0 = scalar_grid.NextVertex(iend0, edge_dir);
          }
        }
      }
      else {
        // For each edge on axis in direction edge_dir,
        //   process each vertex in facet_vertex_list.
        // Change order to processes vertices/edges
        //   that are close together in memory and significantly
        //   speed up this routine.

        // Add increment0 to each vertex in facet_vertex_list to get iend0.
        // Add increment to iend0 to get iend1.
        VTYPE increment0 = 0;
        const VTYPE increment1 = scalar_grid.AxisIncrement(edge_dir);
        for (int j = 0; j < num_edges; j++) {
          for (VTYPE i = 0; i < facet_vertex_list.NumVertices(); i++) {
            const VTYPE iend0 = facet_vertex_list.VertexIndex(i) + increment0;
            extract_dual_isopoly_around_bipolar_edge_I
              (scalar_grid, isovalue, iend0, edge_dir,
               isopoly_cube, isopoly_info);              
          }
          increment0 = increment0 + increment1;
        }
      }
    }

    clock_t t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }
  

  /*!
   *  @brief Extract isosurface polytopes
   *  - Returns list representing isosurface polytopes.
   *  @param scalar_grid Scalar grid data.
   *  @param isovalue Isosurface scalar value.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polygon.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param[out] isopoly_info[i] Isososurface polygon information,
   *    including dual grid edge.
   */
  template <typename GTYPE, typename STYPE,
            typename ICUBE_TYPE, typename IEDGE_TYPE,
            typename ISOPOLY_INFO_TYPE, typename INFO_TYPE>
  void extract_dual_isopoly
  (const GTYPE & scalar_grid, const STYPE isovalue,
   std::vector<ICUBE_TYPE> & isopoly_cube,
   std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   INFO_TYPE & dualiso_info)
  {
    typedef typename GTYPE::VERTEX_INDEX_TYPE VTYPE;
    
    const int dimension = scalar_grid.Dimension();
    dualiso_info.time.extract = 0;

    if (dimension < 1) { return; }
    
    clock_t t0 = clock();

    // initialize output
    isopoly_cube.clear();
    isopoly_cube_edge.clear();
    isopoly_info.clear();

    if (scalar_grid.NumCubeVertices() < 1) { return; }
    IJK::FACET_INTERIOR_VERTEX_LIST<VERTEX_INDEX>
      facet_vertex_list(scalar_grid, 0, true);

    for (int edge_dir = 0; edge_dir < dimension; edge_dir++) {

      if (scalar_grid.AxisSize(edge_dir) < 1) {
        // No interior edges.
        return;
      }

      const int num_edges = scalar_grid.AxisSize(edge_dir)-1;

      facet_vertex_list.GetVertices(scalar_grid, edge_dir);

      if (edge_dir+1 < dimension) {
        // For each vertex in facet_vertex_list,
        //   process each edge in direction edge_dir.
        for (VTYPE i = 0; i < facet_vertex_list.NumVertices(); i++) {
          VTYPE iend0 = facet_vertex_list.VertexIndex(i);

          for (int j = 0; j < num_edges; j++) {

            const VTYPE iend1 = scalar_grid.NextVertex(iend0, edge_dir);
            
            if (scalar_grid.Scalar(iend0) < isovalue) {
              if (scalar_grid.Scalar(iend1) >= isovalue) {
                extract_dual_isopoly_around_edge_multi
                  (scalar_grid, iend0, iend1, edge_dir,
                   isopoly_cube, isopoly_cube_edge, isopoly_info);
              }              
            }
            else {
              // scalar_grid.Scalar(iend0) >= isovalue) {
              if (scalar_grid.Scalar(iend1) < isovalue) {
                extract_dual_isopoly_around_edge_reverse_orient_multi
                  (scalar_grid, iend0, iend1, edge_dir,
                   isopoly_cube, isopoly_cube_edge, isopoly_info);
              }
            }
            
            iend0 = iend1;
          }
        }
      }
      else {
        // For each edge on axis in direction edge_dir,
        //   process each vertex in facet_vertex_list.
        // Change order to processes vertices/edges
        //   that are close together in memory and significantly
        //   speed up this routine.

        // Add increment0 to each vertex in facet_vertex_list to get iend0.
        // Add increment to iend0 to get iend1.
        VTYPE increment0 = 0;
        const VTYPE increment1 = scalar_grid.AxisIncrement(edge_dir);
        for (int j = 0; j < num_edges; j++) {          
          for (VTYPE i = 0; i < facet_vertex_list.NumVertices(); i++) {
            const VTYPE iend0 = facet_vertex_list.VertexIndex(i) + increment0;
            const VTYPE iend1 = scalar_grid.NextVertex(iend0, edge_dir);

            if (scalar_grid.Scalar(iend0) < isovalue) {
              if (scalar_grid.Scalar(iend1) >= isovalue) {
                extract_dual_isopoly_around_edge_multi
                  (scalar_grid, iend0, iend1, edge_dir,
                   isopoly_cube, isopoly_cube_edge, isopoly_info);
              }              
            }
            else {
              // scalar_grid.Scalar(iend0) >= isovalue) {
              if (scalar_grid.Scalar(iend1) < isovalue) {
                extract_dual_isopoly_around_edge_reverse_orient_multi
                  (scalar_grid, iend0, iend1, edge_dir,
                   isopoly_cube, isopoly_cube_edge, isopoly_info);
              }
            }
            
          }
          increment0 = increment0 + increment1;
        }
      }
    }

    clock_t t1 = clock();
    IJK::clock2seconds(t1-t0, dualiso_info.time.extract);
  }

};

#endif
