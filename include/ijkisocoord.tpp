/*!
 *  \file ijkisocoord.tpp
 *  @brief ijk templates for computing isosurface coordinates
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2015-2024 Rephael Wenger

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

#ifndef _IJKISOCOORD_TPP_
#define _IJKISOCOORD_TPP_

#include <vector>

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkinterpolate.tpp"


namespace IJK {

  // *******************************************************************
  //! @name Compute intersections of isosurface and grid edges.
  // *******************************************************************

  ///@{

  /*!
   *  @brief Compute intersection of isosurface and line segment (iv0, iv1)
   *    using linear interpolation.
   *  @param scalar_grid Scalar grid data.
   *  @param ivA First grid vertex.
   *  @param ivB Second grid vertex
   *  @param[out] isov_coord[] Array of isosurface vertex coordinates.
   *    - coord[d] = d'th coordinate of vertex i'th vertex.
   *    @pre Array isov_coord[] is preallocated to length dimension.
   *  @param temp_coord0[] Temporary array for storing coordinates.
   *    @pre Array temp_coord0[] is preallocated to length dimension.
   *  @param temp_coord1[] Temporary array for storing coordinates.
   *    @pre Array temp_coord1[] is preallocated to length dimension.
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename VTYPE, typename CTYPEV, typename CTYPET>
  inline void compute_isosurface_line_segment_intersection_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const VTYPE iv0, const VTYPE iv1, CTYPEV isov_coord[],
   CTYPET temp_coord0[], CTYPET temp_coord1[])
  {
    typedef typename SGRID_TYPE::SCALAR_TYPE STYPE;

    const int dimension = scalar_grid.Dimension();
    const STYPE * scalar = scalar_grid.ScalarPtrConst();
  
    const STYPE s0 = scalar[iv0];
    const STYPE s1 = scalar[iv1];

    scalar_grid.ComputeCoord(iv0, temp_coord0);
    scalar_grid.ComputeCoord(iv1, temp_coord1);

    linear_interpolate_coord
      (dimension, s0, temp_coord0, s1, temp_coord1, isovalue,
       isov_coord);
  }

    
  /*!
   *  @overload
   *  @brief Compute intersection of isosurface and line segment (iv0, iv1)
   *    using linear interpolation.
   *    (No temp_coord*.)
   *  @details
   *  - Version that does not have temp_coord* as arguments.
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename VTYPE, typename CTYPE>
  inline void compute_isosurface_line_segment_intersection_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const VTYPE iv0, const VTYPE iv1, CTYPE * coord)
  {
    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> temp_coord0(dimension);
    IJK::ARRAY<CTYPE> temp_coord1(dimension);

    compute_isosurface_line_segment_intersection_linear
      (scalar_grid, isovalue, iv0, iv1, coord,
       temp_coord0.Ptr(), temp_coord1.Ptr());
  }


  /*!
   *  @brief Compute coordinates of isosurface vertices 
   *    on list of line segments using linear interpolation.
   *  @param scalar_grid Scalar grid data.
   *  @param endpoint[] Array of line segment endpoints.
   *    Line segment i has endpoints (2*endpoint[i], 2*endpoint[i]+1).
   *  @param[out] coord[] Array of isosurface vertex coordinates.
   *   coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
   *  @pre Array coord[] is preallocated to length at least
   *     dimension*(endpoint.size()/2).
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename VTYPE, typename CTYPE>
  void compute_isov_coord_on_line_segment_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const std::vector<VTYPE> & endpoint, CTYPE * coord)
  {
    typedef typename SGRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SGRID_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    const NTYPE num_line_segments = endpoint.size()/2;
  
    for (NTYPE i = 0; i < num_line_segments; i++) {
      compute_isosurface_line_segment_intersection_linear
        (scalar_grid, isovalue, endpoint[2*i], endpoint[2*i+1], 
         coord+dimension*i);
    }
  }


  /*!
   *  @brief Compute intersection of isosurface and grid edge
   *         using linear interpolation.
   *  @param scalar_grid Scalar grid data.
   *  @param grid_edge Grid edge.
   *   - Edge i has lower/leftmost endpoint 
   *     grid_edge.GridEdgeEndpoint0()
   *     and direction grid_edge.GridEdgeDirection().
   *  @param[out] isov_coord[] Array of isosurface vertex coordinates.
   *    - isov_coord[d] = d'th coordinate of isosurface vertex.
   *  @pre Array isov_coord[] is preallocated to length at least dimension.
   *  @param temp_coord0[] Temporary array for storing coordinates.
   *    @pre Array temp_coord0[] is preallocated to length dimension.
   *  @param temp_coord1[] Temporary array for storing coordinates.
   *    @pre Array temp_coord1[] is preallocated to length dimension.
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename GRID_EDGE_TYPE,
            typename CTYPEV, typename CTYPET>
  inline void compute_isosurface_grid_edge_intersection_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const GRID_EDGE_TYPE & grid_edge, CTYPEV isov_coord[],
   CTYPET temp_coord0[], CTYPET temp_coord1[])
  {
    typedef typename SGRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const VTYPE iv0 = grid_edge.GridEdgeEndpoint0();
    const int edge_direction = grid_edge.GridEdgeDirection();
    const VTYPE iv1 = scalar_grid.NextVertex(iv0, edge_direction);

    compute_isosurface_line_segment_intersection_linear
      (scalar_grid, isovalue, iv0, iv1,
       isov_coord, temp_coord0, temp_coord1);
  }


  /*!
   *  @overload
   *  @brief Compute intersection of isosurface and grid edge
   *         using linear interpolation.
   *    (No temp_coord*.)
   *  @details
   *  - Version that does not have temp_coord* as arguments.
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename GRID_EDGE_TYPE, typename CTYPE>
  void compute_isosurface_grid_edge_intersection_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const GRID_EDGE_TYPE & grid_edge, CTYPE isov_coord[])
  {
    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPE> temp_coord0(dimension);
    IJK::ARRAY<CTYPE> temp_coord1(dimension);

    compute_isosurface_grid_edge_intersection_linear
      (scalar_grid, isovalue, grid_edge, isov_coord,
       temp_coord0.Ptr(), temp_coord1.Ptr());
  }


  /*!
   *  @overload
   *  @brief Compute intersection of isosurface and grid edge
   *         using linear interpolation. (Apply offset.)
   *  @details
   *  - Version that applies offset to keep vertex from end of edge.
   *  @param offset Minimum distance to edge endpoint.
   *    @pre Must be in range [0,0.5].
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename GRID_EDGE_TYPE,
            typename CTYPEO, typename CTYPEV, typename CTYPET>
  inline void compute_isosurface_grid_edge_intersection_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const GRID_EDGE_TYPE & grid_edge,
   const CTYPEO offset,
   CTYPEV isov_coord[],
   CTYPET temp_coord0[], CTYPET temp_coord1[])
  {
    typedef typename SGRID_TYPE::SCALAR_TYPE STYPE;
    typedef typename SGRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const int dimension = scalar_grid.Dimension();
    const STYPE * scalar = scalar_grid.ScalarPtrConst();
    const VTYPE iv0 = grid_edge.GridEdgeEndpoint0();
    const int edge_direction = grid_edge.GridEdgeDirection();
    const VTYPE iv1 = scalar_grid.NextVertex(iv0, edge_direction);
    const STYPE s0 = scalar[iv0];
    const STYPE s1 = scalar[iv1];

    scalar_grid.ComputeCoord(iv0, temp_coord0);
    scalar_grid.ComputeCoord(iv1, temp_coord1);

    linear_interpolate_coord
      (dimension, s0, temp_coord0, s1, temp_coord1, isovalue,
       isov_coord);
    
    compute_isosurface_line_segment_intersection_linear
      (scalar_grid, isovalue, iv0, iv1,
       isov_coord, temp_coord0, temp_coord1);

    const CTYPET minc = temp_coord0[edge_direction] + offset;
    const CTYPET maxc = temp_coord1[edge_direction] - offset;
    isov_coord[edge_direction] =
      std::max(isov_coord[edge_direction], minc);
    isov_coord[edge_direction] =
      std::min(isov_coord[edge_direction], maxc);    
  }
  

  /*!
   *  @overload
   *  @brief Compute intersection of isosurface and grid edge
   *         using linear interpolation.
   *    (Apply offset. No temp_coord*.)
   *  @details
   *  - Version that applies offset to keep vertex from end of edge.
   *  - Version that does not have temp_coord* as arguments.
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename GRID_EDGE_TYPE, typename CTYPEO, typename CTYPEV>
  void compute_isosurface_grid_edge_intersection_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const GRID_EDGE_TYPE & grid_edge, const CTYPEO offset,
   CTYPEV isov_coord[])
  {
    const int dimension = scalar_grid.Dimension();
    IJK::ARRAY<CTYPEV> temp_coord0(dimension);
    IJK::ARRAY<CTYPEV> temp_coord1(dimension);

    compute_isosurface_grid_edge_intersection_linear
      (scalar_grid, isovalue, grid_edge, offset, isov_coord,
       temp_coord0.Ptr(), temp_coord1.Ptr());
  }
  

  /*!
   *  Compute all coordinates of isosurface vertices on list of edges
   *    using linear interpolation.
   *  - Version using grid edge stored in isopoly_info.
   *  @param scalar_grid Scalar grid data.
   *  @param isopoly_info[i] Isosurface polytope information containing
   *   dual grid edge to polytope i.
   *   - Edge i has lower/leftmost endpoint 
   *     isopoly_info[i].GridEdgeEndpoint0()
   *     and direction isopoly_info[i].GridEdgeDirection().
   *  @param[out] isov_coord[] Array of isosurface vertex coordinates.
   *   - isov_coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
   *  @pre Array isov_coord[] is preallocated to length at least
   *     dimension*isopoly_info.size().
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE>
  void compute_all_isov_coord_on_grid_edge_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   CTYPE isov_coord[])
  {
    typedef typename SGRID_TYPE::NUMBER_TYPE NTYPE;

    const int dimension = scalar_grid.Dimension();    
    const NTYPE nume = isopoly_info.size();
    IJK::ARRAY<CTYPE> temp_coord0(dimension);
    IJK::ARRAY<CTYPE> temp_coord1(dimension);

    for (NTYPE i = 0; i < nume; i++) {
      CTYPE * isov_coord_i = isov_coord + (dimension*i);
      compute_isosurface_grid_edge_intersection_linear
        (scalar_grid, isovalue, isopoly_info[i],
         isov_coord_i, temp_coord0.Ptr(), temp_coord1.Ptr());
    }
  }


  /*!
   *  @overload
   *  Compute all coordinates of isosurface vertices on list of edges
   *    using linear interpolation. (Apply offset.)
   *  @details
   *  - Version using grid edge stored in isopoly_info.
   *  - Version that applies offset to keep vertex from end of edge.
   *  @param offset Minimum distance to edge endpoint.
   *    @pre Must be in range [0,0.5].
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPEO, typename CTYPEV>
  void compute_all_isov_coord_on_grid_edge_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info,
   const CTYPEO offset, CTYPEV isov_coord[])
  {
    typedef typename SGRID_TYPE::NUMBER_TYPE NTYPE;

    const int dimension = scalar_grid.Dimension();    
    const NTYPE nume = isopoly_info.size();
    IJK::ARRAY<CTYPEV> temp_coord0(dimension);
    IJK::ARRAY<CTYPEV> temp_coord1(dimension);

    for (NTYPE i = 0; i < nume; i++) {
      CTYPEV * isov_coord_i = isov_coord + (dimension*i);
      compute_isosurface_grid_edge_intersection_linear
        (scalar_grid, isovalue, isopoly_info[i], offset,
         isov_coord_i, temp_coord0.Ptr(), temp_coord1.Ptr());
    }
  }

    
  ///@}


  // *******************************************************************
  //! @name Compute isosurface vertex coordinates on scaled grid edges.
  // *******************************************************************

  ///@{

  /*!
   *  @brief Compute coordinates of isosurface vertices 
   *    on scaled line segment (iv0, iv1) using linear interpolation.
   *  @param scalar_grid Scalar grid data with scaling data.
   *  @param ivA First grid vertex.
   *  @param ivB Second grid vertex
   *  @param[out] coord[] Array of isosurface vertex coordinates.
   *   coord[d] = d'th coordinate of vertex i'th vertex.
   *  @pre Array coord[] is preallocated to length dimension.
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename VTYPE, typename CTYPE>
  void compute_isov_coord_on_scaled_line_segment_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const VTYPE iv0, const VTYPE iv1, CTYPE * coord)
  {
    typedef typename SGRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SGRID_TYPE::SCALAR_TYPE STYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    const STYPE * scalar = scalar_grid.ScalarPtrConst();
    IJK::ARRAY<CTYPE> coord0(dimension);
    IJK::ARRAY<CTYPE> coord1(dimension);
  
    STYPE s0 = scalar[iv0];
    STYPE s1 = scalar[iv1];

    scalar_grid.ComputeScaledCoord(iv0, coord0.Ptr());
    scalar_grid.ComputeScaledCoord(iv1, coord1.Ptr());

    linear_interpolate_coord
      (dimension, s0, coord0.PtrConst(), s1, coord1.PtrConst(), 
       isovalue, coord);
  }


  /*!
   *  @brief Compute coordinates of isosurface vertices 
   *    on list of scaled line segments using linear interpolation.
   *  @param scalar_grid Scalar grid data.
   *  @param endpoint[] Array of line segment endpoints.
   *    Line segment i has endpoints (2*endpoint[i], 2*endpoint[i]+1).
   *  @param[out] coord[] Array of isosurface vertex coordinates.
   *   coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
   *  @pre Array coord[] is preallocated to length at least
   *     dimension*(endpoint.size()/2).
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename VTYPE, typename CTYPE>
  void compute_isov_coord_on_scaled_line_segment_linear
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const std::vector<VTYPE> & endpoint, CTYPE * coord)
  {
    typedef typename SGRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SGRID_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    const NTYPE num_line_segments = endpoint.size()/2;
  
    for (NTYPE i = 0; i < num_line_segments; i++) {
      compute_isov_coord_on_scaled_line_segment_linear
        (scalar_grid, isovalue, endpoint[2*i], endpoint[2*i+1], 
         coord+dimension*i);
    }
  }


  /*!
   *  @brief Compute coordinates of isosurface vertices 
   *    on a single scaled grid edge using linear interpolation.
   *  - Version using grid edge stored in isopoly_info.
   *  @param scalar_grid Scalar grid data.
   *  @param isopoly_info Isosurface polytope information containing
   *   dual grid edge.
   *   - Edge i has lower/leftmost endpoint 
   *     isopoly_info.GridEdgeEndpoint0()
   *     and direction isopoly_info.GridEdgeDirection().
   *  @param[out] coord[] Array of isosurface vertex coordinates.
   *   coord[d] = d'th coordinate of isosurface vertex.
   *  @pre Array coord[] is preallocated to length at least dimension.
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE>
  void compute_isov_coord_on_scaled_grid_edge_linear_I
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const ISOPOLY_INFO_TYPE & isopoly_info, CTYPE * coord)
  {
    typedef typename SGRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SGRID_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const VTYPE iv0 = isopoly_info.GridEdgeEndpoint0();
    const DTYPE edge_direction = isopoly_info.GridEdgeDirection();
    const VTYPE iv1 = scalar_grid.NextVertex(iv0, edge_direction);

    compute_isov_coord_on_scaled_line_segment_linear
      (scalar_grid, isovalue, iv0, iv1, coord);
  }


  /*!
   *  Compute coordinates of isosurface vertices on list of scaled edges
   *    using linear interpolation.
   *  - Version using grid edge stored in isopoly_info.
   *  @param scalar_grid Scalar grid data.
   *  @param isopoly_info[i] Isosurface polytope information containing
   *   dual grid edge to polytope i.
   *   - Edge i has lower/leftmost endpoint 
   *     isopoly_info[i].GridEdgeEndpoint0()
   *     and direction isopoly_info[i].GridEdgeDirection().
   *  @param[out] coord[] Array of isosurface vertex coordinates.
   *   coord[i*dimension+j] = j'th coordinate of vertex i'th vertex
   *  @pre Array coord[] is preallocated to length at least
   *     dimension*isopoly_info.size().
   */
  template <typename SGRID_TYPE, typename ISOVALUE_TYPE,
            typename ISOPOLY_INFO_TYPE, typename CTYPE>
  void compute_isov_coord_on_scaled_grid_edge_linear_I
  (const SGRID_TYPE & scalar_grid, const ISOVALUE_TYPE isovalue,
   const std::vector<ISOPOLY_INFO_TYPE> & isopoly_info, CTYPE * coord)
  {
    typedef typename SGRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename SGRID_TYPE::NUMBER_TYPE NTYPE;

    const DTYPE dimension = scalar_grid.Dimension();
    const NTYPE nume = isopoly_info.size();

    for (NTYPE i = 0; i < nume; i++) {
      compute_isov_coord_on_scaled_grid_edge_linear_I
        (scalar_grid, isovalue, isopoly_info[i], coord+dimension*i);
    }
  }

}

#endif

