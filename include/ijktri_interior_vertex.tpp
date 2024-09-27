/*!
 *  @file ijktri_interior_vertex.tpp
 *  @brief ijk types, data structures and templates for determining
 *    coordinates of vertices in interior of isosurface polytopes.
 *  - Used in adding vertices to triangulate isosurface polytopes.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2023 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#ifndef _IJKTRI_INTERIOR_VERTEX_
#define _IJKTRI_INTERIOR_VERTEX_

#include "ijk.tpp"

namespace IJK {

  // ***************************************************************
  //! @name TYPES AND DATA STRUCTURES
  // ***************************************************************

  //@{

  /*!
   *  @brief Position method for additional interior vertices to isosurface polytopes.
   *  - INTERNAL_VERTEX_ON_GRID_EDGE: Add vertex on grid edge dual 
   *    to isosurface polytope.
   *  - INTERNAL_VERTEX_IN_ENVELOPE: Add vertex so that edges lie
   *    in (relative) interior of envelope.
   *  - INTERNAL_VERTEX_AT_CENTROID: Add vertex at centroid 
   *    of isosurface polytope.
   */
  typedef enum
    { INTERIOR_VERTEX_ON_GRID_EDGE, INTERIOR_VERTEX_IN_ENVELOPE,
      INTERIOR_VERTEX_AT_CENTROID }
    INTERIOR_VERTEX_POSITION_METHOD;

  /*!
   *  @brief Polytope-edge intersection method.
   *  - Method for determining intersection of isosurface polytope and grid edge.
   */
  typedef enum
    { POLY_EDGE_MULTILINEAR_INTERPOLATION, INTERPOLATE_EDGE_ENDPOINT_SCALARS,
      POLY_EDGE_AVERAGE_PROJECTION, POLY_EDGE_WEIGHTED_AVERAGE_PROJECTION }
    POLY_EDGE_INTERSECTION_METHOD;


  /// Return string (const char *) representing interior vertex position method.
  inline const char * get_interior_vertex_position_method_str
  (const INTERIOR_VERTEX_POSITION_METHOD position_method)
  {
    if (position_method == INTERIOR_VERTEX_ON_GRID_EDGE)
      { return "OnGridEdge"; }
    else if (position_method == INTERIOR_VERTEX_IN_ENVELOPE)
      { return "InEnvelope"; }
    else if (position_method == INTERIOR_VERTEX_AT_CENTROID)
      { return "AtCentroid"; }
    else
      { return "Unknown"; }
  }

  
  /// Return string (const char *) representing poly edge intersection method.
  inline const char * get_poly_edge_intersection_method_str
  (const POLY_EDGE_INTERSECTION_METHOD intersection_method)
  {
    if (intersection_method == POLY_EDGE_MULTILINEAR_INTERPOLATION)
      { return "MultilinearInterpolation"; }
    else if (intersection_method == INTERPOLATE_EDGE_ENDPOINT_SCALARS)
      { return "InterpolateEdgeEndpointScalars"; }
    else if (intersection_method == POLY_EDGE_AVERAGE_PROJECTION)
      { return "AverageProjection"; }
    else if (intersection_method == POLY_EDGE_WEIGHTED_AVERAGE_PROJECTION)
      { return "WeightedAverageProjection"; }
    else
      { return "Unknown"; }
  }
  

  /// Parameters for positioning vertex in interior of an isosurface polytope.
  template <typename RTYPE>
  class INTERIOR_VERTEX_PARAM_T {
    
  public:
    INTERIOR_VERTEX_POSITION_METHOD interior_vertex_position_method;
    POLY_EDGE_INTERSECTION_METHOD poly_edge_intersection_method;

    /*!
     *  @brief Ratio of distance of interior vertex to cube facet 
     *    over distance of isosurface polytope vertex to cube facet.
     *  - Must be less than 0.5.
     */
    RTYPE interior_vertex_ratio_to_cube_facet;

  public:
    
    /// Constructor.
    INTERIOR_VERTEX_PARAM_T();

    /// Return default interior vertex ratio to cube facet.
    RTYPE DefaultInteriorVertexRatioToCubeFacet() const
    { return 0.375; }

    /// @brief Return true if interior_vertex_position_method ==
    ///   INTERIOR_VERTEX_AT_CENTROID.
    bool PositionAtCentroid() const
    {
      return (interior_vertex_position_method ==
              INTERIOR_VERTEX_AT_CENTROID);
    }

    /// @brief Return true if interior_vertex_position_method ==
    ///   INTERIOR_VERTEX_IN_ENVELOPE
    bool PositionInEnvelope() const
    {
      return (interior_vertex_position_method ==
              INTERIOR_VERTEX_IN_ENVELOPE);
    }

    /// @brief Return true if interior_vertex_position_method ==
    ///   INTERIOR_VERTEX_ON_GRID_EDGE
    bool PositionOnGridEdge() const
    {
      return (interior_vertex_position_method ==
              INTERIOR_VERTEX_ON_GRID_EDGE);
    }

      
    // Print routines (mainly for debugging.)

    template <typename OSTREAM_TYPE, typename STYPE0>
    void Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

    template <typename OSTREAM_TYPE, typename STYPE0>
    void PrintInteriorVertexPositionMethod
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;

    template <typename OSTREAM_TYPE, typename STYPE0>
    void PrintPolyEdgeIntersectionMethod
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const;
  };
  
  //@}

  
  // ***************************************************************
  //! @name BASE CLASS FOR ADDING DUAL ISOSURFACE VERTICES.
  // ***************************************************************

  //@{

  /// Base class for adding dual isosurface vertices.
  class ADD_DUAL_ISOV_BASE {
      
  public:

    /*!
     *  @brief Function to add isosurface vertex dual 
     *    to isosurface polytope.
     *  - Implement this routine in derived classes.
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

  //@}

  
  // ***************************************************************
  //! @name CLASS INTERIOR_VERTEX_PARAM_T MEMBER FUNCTIONS
  // ***************************************************************

  //@{

  // Constructor
  template <typename RTYPE>
  INTERIOR_VERTEX_PARAM_T<RTYPE>::INTERIOR_VERTEX_PARAM_T()
  {
    interior_vertex_position_method = INTERIOR_VERTEX_ON_GRID_EDGE;
    poly_edge_intersection_method = POLY_EDGE_MULTILINEAR_INTERPOLATION;
    interior_vertex_ratio_to_cube_facet =
      DefaultInteriorVertexRatioToCubeFacet();
  }


  template <typename RTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void INTERIOR_VERTEX_PARAM_T<RTYPE>::
  Print(OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    PrintInteriorVertexPositionMethod(out, line_prefix);
    PrintPolyEdgeIntersectionMethod(out, line_prefix);
    out << line_prefix << "interior_vertex_ratio_to_cube_facet: "
        << interior_vertex_ratio_to_cube_facet << "\n";
  }
  

  template <typename RTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void INTERIOR_VERTEX_PARAM_T<RTYPE>::
  PrintInteriorVertexPositionMethod
  (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    const char * position_method_str =
      get_interior_vertex_position_method_str(interior_vertex_position_method);
    
    out << line_prefix
        << "interior_vertex_position_method: "
        << position_method_str << "\n";
  }


  template <typename RTYPE>
  template <typename OSTREAM_TYPE, typename STYPE0>
  void INTERIOR_VERTEX_PARAM_T<RTYPE>::
  PrintPolyEdgeIntersectionMethod
    (OSTREAM_TYPE & out, const STYPE0 & line_prefix) const
  {
    const char * intersection_method_str =
      get_poly_edge_intersection_method_str
      (poly_edge_intersection_method);
    
    out << line_prefix
        << "poly_edge_intersection_method: "
        << intersection_method_str << "\n";
  }

  
  //@}

}

#endif

