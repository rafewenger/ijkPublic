/*!
 *  @file ijkmesh2D_triangulate.tpp
 *  @brief ijk template classes for triangulating 2D half-edge mesh.
 *  - Version 0.4.0
 */

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2019-2022 Rephael Wenger

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

#ifndef _IJKMESH2D_TRIANGULATE_
#define _IJKMESH2D_TRIANGULATE_

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkmesh2D_datastruct.tpp"
#include "ijkmesh2D_splitII.tpp"
#include "ijkmesh2D_geom.tpp"
#include "ijkmesh2D_angle.tpp"
#include "ijktri2D_info.tpp"
#include "ijktri2D_angle.tpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // *****************************************************************
  // Class VERTEX_TRI_INFO_BASE
  // *****************************************************************

  /// Vertex base type.
  template <typename IV_TYPE, typename HTYPE, typename NTYPE>
  class VERTEX_TRI_INFO_BASE:public VERTEX_BASE<IV_TYPE,HTYPE,NTYPE> {

  public:

    /// Triangulation vertex type.
    typedef enum { TRIV_INPUT, TRIV_INTERIOR, TRIV_SPLITE } 
      TRIANGULATION_VERTEX_TYPE;

    TRIANGULATION_VERTEX_TYPE triangulation_vertex_type;

    /// Return true if triangulation type is TRIV_SPLITE.
    bool SplitsEdge() const
    { return((triangulation_vertex_type == TRIV_SPLITE)); }
  };


  // *****************************************************************
  // Class POLYGON_TRI_INFO_BASE
  // *****************************************************************

  /*!
   *  @brief Polygon with triangulation information
   *  @tparam DIMENSION Dimension of space containing mesh.
   *     Equivalently, dimension of the vertex coordinates.
   */
  template <int DIMENSION, int BIT_SET_SIZE,
            typename COORD_TYPE_X, typename COS_TYPE_X,
            typename ITYPE, typename NTYPE>
  class POLYGON_TRI_INFO_BASE:public POLYGON_BASE {

  protected:
    
  public:

    /// Dimension type.
    typedef int DIMENSION_TYPE;

    /// Coord type.
    typedef COORD_TYPE_X COORD_TYPE;

    /// Cosine type.
    typedef COS_TYPE_X COS_TYPE;

    /// Triangulation information.
    typedef POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE, COS_TYPE_X, NTYPE> 
    TRIANGULATION_INFO_TYPE;

    static const int BitSetSize() { return(BIT_SET_SIZE); };

  public:

    /// Number of vertices added to split edges.
    NTYPE num_split_edge_vertices;

  public:

    /// Coordinates of interior point for possible use in the polygon triangulation.
    COORD_TYPE_X interior_coord[DIMENSION];

    /// Polygon triangulation information.
    TRIANGULATION_INFO_TYPE triangulation_info;

    /// Redefine Init().
    void Init();

    /// Redefine Init().
    template <typename NTYPE2>
    void Init(const NTYPE2 num_edges) {
      Init();
      POLYGON_BASE::Init(num_edges);
    };

  public:

    static DIMENSION_TYPE Dimension()
    { return(DIMENSION); }


  public:
    POLYGON_TRI_INFO_BASE() { Init(); };

    const COORD_TYPE_X * InteriorCoord() const
    { return(interior_coord); }

    COS_TYPE_X CosMinTriangulationAngle() const
    { return(triangulation_info.cos_min_triangulation_angle); }

    ITYPE TriVertexIndex() const
    { return(triangulation_info.tri_vertex_index); }

    NTYPE NumInteriorTriVertices() const
    { return(triangulation_info.num_interior_tri_vertices); }

  };


  // *****************************************************************
  // Class MESH2D_TRIANGULATION_BASE
  // *****************************************************************

  /*!
   *  @brief Mesh of 2D polygons for computing triangulations.
   *  \tparam VTYPE Vertex type.
   *          Should be derived from VERTEX_TRI_INFO_BASE.
   *  \tparam HALFE_TYPE Half edge type.
   *          Should be derived from HALF_EDGE_BASE.
   *  \tparam PTYPE Polygon type.
   *          Should be derived from POLYGON_TRI_INFO_BASE.
   *  *** PROBLEM: BIT_SET_SIZE is set.  SHOULD BE DERIVED FROM PTYPE.
   */
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  class MESH2D_TRIANGULATION_BASE:
    public MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE> {

  public:

    typedef typename PTYPE::DIMENSION_TYPE DIMENSION_TYPE;
    typedef typename PTYPE::COS_TYPE COS_TYPE;
    typedef typename PTYPE::TRIANGULATION_INFO_TYPE 
    TRIANGULATION_INFO_TYPE;

    typedef MESH2D_SPLIT_II_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>
    MESH2D_SPLIT_II_BASE_TYPE;

    typedef typename MESH2D_SPLIT_II_BASE_TYPE::POLYGON_INDEX_TYPE 
    POLYGON_INDEX_TYPE;

    typedef typename MESH2D_SPLIT_II_BASE_TYPE::VERTEX_INDEX_TYPE 
    VERTEX_INDEX_TYPE;

    typedef typename MESH2D_SPLIT_II_BASE_TYPE::HALF_EDGE_INDEX_TYPE 
    HALF_EDGE_INDEX_TYPE;

    typedef typename MESH2D_SPLIT_II_BASE_TYPE::NUMBER_TYPE NUMBER_TYPE;

    typedef typename VTYPE::TRIANGULATION_VERTEX_TYPE
    TRIANGULATION_VERTEX_TYPE;

    static constexpr int BitSetSize() { return(BIT_SET_SIZE); }


  private:

    // Implementations that do not enforce BIT_SET_SIZE.
    // Called by routines that require BIT_SET_SIZE 
    //   to match PTYPE::BitSetSize().
    // Should never be called outside this class.


  public:

    /// @brief Return dimension.
    DIMENSION_TYPE Dimension() const
    { return(PTYPE::Dimension()); }

    /// @brief Return number of polygon vertices added to split edges.
    template <typename ITYPE2>
    NTYPE NumSplitEdgeVertices(const ITYPE2 ipoly) const
    { return(this->polygon_list[ipoly].num_split_edge_vertices); }

    /*!
     *  @brief Return pointer to coordinates of possible triangulation vertex in interior of polygon ipoly.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on polygon ipoly.
     */
    template <typename ITYPE2>
    const typename PTYPE::COORD_TYPE * PolygonInteriorCoord
    (const ITYPE2 ipoly) const
    { return(this->polygon_list[ipoly].InteriorCoord()); }

    /*!
     *  @brief Return cosine of the min angle of the triangulation specified for polygon ipoly.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on polygon ipoly.
     */
    template <typename ITYPE2>
    COS_TYPE PolygonCosMinTriangulationAngle(const ITYPE2 ipoly) const
    { return(this->polygon_list[ipoly].CosMinTriangulationAngle()); }

    /*!
     *  @brief Return cosine of the min angle of the triangulation specified for two polygons.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on polygon ipoly.
     */
    template <typename ITYPE2, typename ITYPE3>
    COS_TYPE PolygonIICosMinTriangulationAngle
    (const ITYPE2 ipoly0, const ITYPE3 ipoly1) const
    { return(std::max(PolygonCosMinTriangulationAngle(ipoly0),
                      PolygonCosMinTriangulationAngle(ipoly1))); }

    /*!
     *  @brief Return cosine of the min angle of the triangulation specified for three polygons.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on polygon ipoly.
     */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    COS_TYPE PolygonIIICosMinTriangulationAngle
    (const ITYPE2 ipoly0, const ITYPE3 ipoly1, const ITYPE4 ipoly2) const
    { return(std::max(PolygonIICosMinTriangulationAngle(ipoly0, ipoly1),
                      PolygonCosMinTriangulationAngle(ipoly2))); }

    /*!
     *  @brief Return cosine of the min angle of the triangulation specified for three polygons.
     *  - Version with array ipoly[].
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on polygon ipoly.
     */
    template <typename ITYPEP>
    COS_TYPE PolygonIIICosMinTriangulationAngle(const ITYPEP ipoly[3])
    { return(PolygonIIICosMinTriangulationAngle
             (ipoly[0], ipoly[1], ipoly[2])); }

    /*!
     *  @brief Return cosine of the min angle of the triangulation specified for four polygons.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on polygon ipoly.
     */
    template <typename ITYPE2, typename ITYPE3, 
              typename ITYPE4, typename ITYPE5>
    COS_TYPE PolygonIVCosMinTriangulationAngle
    (const ITYPE2 ipoly0, const ITYPE3 ipoly1, const ITYPE4 ipoly2,
     const ITYPE5 ipoly3) const
    { return(std::max(PolygonIICosMinTriangulationAngle(ipoly0, ipoly1),
                      PolygonIICosMinTriangulationAngle(ipoly2, ipoly3))); }

    /*!
     *  @brief Return cosine of the min angle of the triangulation specified for polygon ipoly and its neighbors.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on polygon ipoly.
     */
    template <typename ITYPE2>
    COS_TYPE PolygonAndNeighborsCosMinTriangulationAngle
    (const ITYPE2 ipoly) const;


    /// @brief Return index of polygon vertex such that triangulation from vertex tri_vertex_index maximizes the minimun triangulation angle.
    /// - tri_vertex_index is in range [0..(NumPolygonVert(ipoly)-1)].
    template <typename ITYPE2>
    ITYPE TriVertexIndex(const ITYPE2 ipoly) const
    { return(this->polygon_list[ipoly].TriVertexIndex()); }

    /// @brief Return the number of interior vertices in the triangulation that maximizes the minimum triangulation angle of polygon ipoly.
    template <typename ITYPE2>
    NTYPE NumInteriorTriVertices(const ITYPE2 ipoly) const
    { return(this->polygon_list[ipoly].NumInteriorTriVertices()); }

    /// Return vertex type.
    template <typename ITYPE2>
    TRIANGULATION_VERTEX_TYPE
    TriangulationVertexType(const ITYPE2 iv) const
    { return(this->vertex_list[iv].triangulation_vertex_type); }

    /// Return triangulation information.
    template <typename ITYPE2>
    const TRIANGULATION_INFO_TYPE &
    PolygonTriangulationInfo(const ITYPE2 ipoly) const
    { return(this->polygon_list[ipoly].triangulation_info); }

    /// Return true if vertex splits edge.
    template <typename ITYPE2>
    bool DoesVertexSplitEdge(const ITYPE2 iv) const
    { return(this->vertex_list[iv].SplitsEdge()); }

    /// Return true if from vertex splits edge.
    template <typename ITYPE2>
    bool DoesFromVertexSplitEdge(const ITYPE2 ihalf_edge) const
    { return(DoesVertexSplitEdge(this->FromVertexIndex(ihalf_edge))); }

    /// Return true if to vertex splits edge.
    template <typename ITYPE2>
    bool DoesToVertexSplitEdge(const ITYPE2 ihalf_edge) const
    { return(DoesVertexSplitEdge(this->ToVertexIndex(ihalf_edge))); }

    /// Returns true if some edge endpoint is a split vertex.
    template <typename ITYPE2>
    bool DoesSomeEndpointSplitEdge(const ITYPE2 ihalf_edge) const
    { return(DoesFromVertexSplitEdge(ihalf_edge) || 
             DoesToVertexSplitEdge(ihalf_edge)); }

    /// Return true if next edge has a split vertex endpoint.
    template <typename ITYPE2>
    bool IsNextEdgeSplit(const ITYPE2 ihalf_edge) const
    { return(DoesSomeEndpointSplitEdge
             (this->IndexOfNextHalfEdgeInPolygon(ihalf_edge))); }

    /// Return true if previous edge has a split vertex endpoint.
    template <typename ITYPE2>
    bool IsPrevEdgeSplit(const ITYPE2 ihalf_edge) const
    { return(DoesSomeEndpointSplitEdge
             (this->IndexOfPrevHalfEdgeInPolygon(ihalf_edge))); }

    /*!
     *  @brief Return true if two adjacent polygons are different.
     *  - Return true if NextHalfEdgeAroundEdge(ihalf_edgeA) is
     *  different from NextHalfEdgeAroundEdge(ihalf_edgeB).
     *  - Return true if ihalf_edgeA or ihalf_edgeB are boundary edges.
     */
    template <typename ITYPEH>
    bool AreAdjacentPolygonsDifferent
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB) const;


    // Add vertices.

    /// Set vertex type.
    template <typename ITYPE2, typename TRIV_TYPE>
    void SetVertexType(const ITYPE2 iv, const TRIV_TYPE triv_type)
    { this->vertex_list[iv].triangulation_vertex_type = triv_type; }

    /// Add interior vertex.
    template <typename CTYPE2, typename CTYPE3>
    VERTEX_INDEX_TYPE AddInteriorVertex
    (const CTYPE2 interior_vcoord[], std::vector<CTYPE3> & vertex_coord);

    /// @brief Add interior vertex.
    /// - Version using C++ STL vector interior_vcoord[].
    template <typename CTYPE2, typename CTYPE3>
    VERTEX_INDEX_TYPE AddInteriorVertex
    (const std::vector<CTYPE2> & interior_vcoord, 
     std::vector<CTYPE3> & vertex_coord);

    /// Add split vertex.
    template <typename CTYPE2, typename CTYPE3>
    VERTEX_INDEX_TYPE AddSplitVertex
    (const CTYPE2 split_vcoord[], std::vector<CTYPE3> & vertex_coord);

    /// @brief Add split vertex.
    /// - Version using C++ STL vector split_vcoord[].
    template <typename CTYPE2, typename CTYPE3>
    VERTEX_INDEX_TYPE AddSplitVertex
    (const std::vector<CTYPE2> & split_vcoord, 
     std::vector<CTYPE3> & vertex_coord);

    /// @brief Add split vertex and split polygon edge.
    /// - Set polygon triangulation information.
    /// - Allow polygon triangulations using interior vertices.
    template <typename CTYPE2, typename CTYPE3,
              typename CTYPE4, typename CTYPE5,
              typename ITYPE2, typename ITYPE3, 
              typename ITYPE4, typename ITYPE5, 
              typename COS_TYPE2, typename NTYPE2>
    void AddVertexSplitPolygonEdgeIIAllowIV
    (const std::vector<CTYPE2> & split_vcoord, 
     const std::vector<CTYPE3> & interior_coordA, 
     const std::vector<CTYPE4> & interior_coordB, 
     std::vector<CTYPE5> & vertex_coord,
     const ITYPE2 ihalf_edge0, 
     POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE, COS_TYPE2, NTYPE2> & polyA_tri_info,
     POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE, COS_TYPE2, NTYPE2> & polyB_tri_info,
     ITYPE3 & iv_split, ITYPE4 & ipolyA, ITYPE5 & ipolyB);

    /// @brief Add split vertex and split polygon edge.
    /// - Set polygon triangulation information.
    /// - Allow polygon triangulations using interior vertices.
    /// - Version that does not return iv_split or ipolyA or ipolyB.
    template <typename CTYPE2, typename CTYPE3,
              typename CTYPE4, typename CTYPE5, 
              typename ITYPE2, typename COS_TYPE2, typename NTYPE2>
    void AddVertexSplitPolygonEdgeIIAllowIV
    (const std::vector<CTYPE2> & split_vcoord, 
     const std::vector<CTYPE3> & interior_coordA, 
     const std::vector<CTYPE4> & interior_coordB, 
     std::vector<CTYPE5> & vertex_coord,
     const ITYPE2 ihalf_edge0, 
     POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE, COS_TYPE2, NTYPE2> & 
     polyA_tri_info,
     POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE, COS_TYPE2, NTYPE2> & 
     polyB_tri_info);

    /// @brief Add new triangle vertex and replace triangle with triangle edge.
    /// - Set polygon triangulation information.
    /// - Allow polygon triangulations using interior vertices.
    template <typename CTYPE2, typename CTYPE3, 
              typename CTYPE4, typename CTYPE5,
              typename ITYPE2, typename ITYPE3, typename ITYPE4,
              typename RESULT_TYPEA, typename RESULT_TYPEB,
              typename RESULT_TYPEC>
    void AddVertexReplaceTriangleWithTriangleEdge
    (const std::vector<CTYPE2> & new_triangle_vcoord,
     const std::vector<CTYPE3> & interior_vcoord0,
     const std::vector<CTYPE4> & interior_vcoord2,
     std::vector<CTYPE5> & vertex_coord,
     const ITYPE2 ihalf_edge0, 
     const RESULT_TYPEA & polyA_tri_info,
     const RESULT_TYPEB & triangleB_tri_info,
     const RESULT_TYPEC & polyC_tri_info,
     ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0);

    /// @brief Add new triangle vertex and replace triangle with triangle edge.
    /// - Set polygon triangulation information.
    /// - Allow polygon triangulations using interior vertices.
    /// - Version using array poly_tri_info[3].
    template <typename CTYPE2, typename CTYPE3, 
              typename CTYPE4, typename CTYPE5,
              typename ITYPE2, typename ITYPE3, typename ITYPE4,
              typename RESULT_TYPE>
    void AddVertexReplaceTriangleWithTriangleEdge
    (const std::vector<CTYPE2> & new_triangle_vcoord,
     const std::vector<CTYPE3> & interior_vcoord0,
     const std::vector<CTYPE4> & interior_vcoord2,
     std::vector<CTYPE5> & vertex_coord,
     const ITYPE2 ihalf_edge0, const RESULT_TYPE poly_tri_info[3],
     ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0);

    /// @brief Add new triangle vertex and replace triangle with triangle edge.
    /// - Allow polygon triangulations using interior vertices.
    /// - Set polygon triangulation information.
    /// - Version using array poly_tri_info[3].
    /// - Version that does not return iv_new or ihalf_edge_new0.
    template <typename CTYPE2, typename CTYPE3, 
              typename CTYPE4, typename CTYPE5,
              typename ITYPE2, typename RESULT_TYPE>
    void AddVertexReplaceTriangleWithTriangleEdge
    (const std::vector<CTYPE2> & new_triangle_vcoord,
     const std::vector<CTYPE3> & interior_vcoord0,
     const std::vector<CTYPE4> & interior_vcoord2,
     std::vector<CTYPE5> & vertex_coord,
     const ITYPE2 ihalf_edge0, const RESULT_TYPE poly_tri_info[3]);

    /// @brief Add new triangle vertex and replace split edge triangle
    ///   with triangle edge.
    /// - Allow polygon triangulations using interior vertices.
    /// - Set polygon triangulation information.
    /// - Version using array poly_tri_info[3].
    template <typename CTYPE2, typename CTYPE3, 
              typename CTYPE4, typename CTYPE5,
              typename ITYPE2, typename ITYPE3, typename ITYPE4,
              typename RESULT_TYPE>
    void AddVertexReplaceSplit2EdgeTriangle
    (const std::vector<CTYPE2> & new_triangle_vcoord,
     const std::vector<CTYPE3> & interior_vcoord0,
     const std::vector<CTYPE4> & interior_vcoord2,
     std::vector<CTYPE5> & vertex_coord,
     const ITYPE2 ihalf_edge0, const RESULT_TYPE poly_tri_info[3],
     ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0);

    /// @brief Add new triangle vertex and replace split edge triangle with triangle edge.
    /// - Allow polygon triangulations using interior vertices.
    /// - Set polygon triangulation information.
    /// - Version using array poly_tri_info[3].
    /// - Version that does not return iv_new or ihalf_edge_new0.
    template <typename CTYPE2, typename CTYPE3, 
              typename CTYPE4, typename CTYPE5,
              typename ITYPE2, typename RESULT_TYPE>
    void AddVertexReplaceSplit2EdgeTriangle
    (const std::vector<CTYPE2> & new_triangle_vcoord,
     const std::vector<CTYPE3> & interior_vcoord0,
     const std::vector<CTYPE4> & interior_vcoord2,
     std::vector<CTYPE5> & vertex_coord,
     const ITYPE2 ihalf_edge0, const RESULT_TYPE poly_tri_info[3]);


    // Compute angles.

    /// @brief Compute cosine of the min angle of triangle ipoly.
    /// - Store the min angle in triangle ipoly.
    /// @pre Polygon ipoly is a triangle.
    template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
    void ComputeCosMinTriangleAngle
      (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
       const MTYPE2 max_small_magnitude);

    /// @brief Compute cosine of the min angle of polygon ipoly.
    /// - Does NOT store information in data structure.
    /// @param[out] cos_min_angle Cosine of the min polygon angle.
    /// @param[out] iloc Location in polygon of vertex with min angle.
    template <typename CTYPE2, typename ITYPE2, typename ITYPE3,
              typename MTYPE2, typename COS_TYPE2>
    void ComputeCosMinPolygonAngle
      (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
       const MTYPE2 max_small_magnitude,
       COS_TYPE2 & cos_min_angle, ITYPE3 & iloc) const;
    

    // Compute midpoint/centroid.

    /// @brief Compute midpoint of line segment between two mesh vertices.
    /// @param[out] midpoint_coord[] Midpoint coordinates.
    /// @pre Array midpoint_coord[] is preallocated to length 
    ///      at least Dimension().
    template <typename CTYPE2, typename CTYPE3, 
              typename ITYPE2, typename ITYPE3>
    void ComputeLineSegmentMidpoint
    (const std::vector<CTYPE2> & vcoord, 
     const ITYPE2 iv0, const ITYPE3 iv1,
     std::vector<CTYPE3> & midpoint_coord) const;

    /// @brief Compute edge midpoint.
    /// @param[out] midpoint_coord[] Midpoint coordinates.
    /// @pre Array midpoint_coord[] is preallocated to length 
    ///      at least Dimension().
    template <typename CTYPE2, typename CTYPE3, typename ITYPE2>
    void ComputeEdgeMidpoint
    (const std::vector<CTYPE2> & vcoord, const ITYPE2 ihalf_edge, 
     std::vector<CTYPE3> & midpoint_coord) const;

    /// @brief Compute midpoint of edge midpoints.
    /// @param[out] midpoint_coord[] Midpoint coordinates.
    /// @pre Array midpoint_coord[] is preallocated to length 
    ///      at least Dimension().
    template <typename CTYPE2, typename CTYPE3, 
              typename ITYPE2, typename ITYPE3>
    void ComputeMidpointOfEdgeMidpoints
    (const std::vector<CTYPE2> & vcoord, 
     const ITYPE2 ihalf_edge2, const ITYPE3 ihalf_edge3,
     std::vector<CTYPE3> & midpoint_coord) const;

    /// @brief Compute polygon centroid.
    /// @param[out] centroid_coord[] Centroid coordinates.
    /// @pre Array centroid_coord[] is preallocated to length 
    ///      at least Dimension().
    template <typename CTYPE2, typename CTYPE3, typename ITYPE2>
    void ComputePolygonCentroid
    (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
     CTYPE3 centroid_coord[]) const;

    /// @brief Compute polygon centroid.
    /// - Version using C++ STL vector vcoord[].
    /// - Version using C++ STL vector centroid_coord[].
    /// @param[out] centroid_coord[] Centroid coordinates.
    /// @pre Array centroid_coord[] is preallocated to length 
    ///      at least Dimension().
    template <typename CTYPE2, typename CTYPE3, typename ITYPE2>
    void ComputePolygonCentroid
    (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly,
     std::vector<CTYPE3> & centroid_coord) const;


    // Compute triangulations.

    /// @brief Compute polygon triangulation that maximizes the minimum angle.
    /// - Store the triangulation information in polygon ipoly.
    template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
    void ComputeCosMaxMinAngle
    (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
     const MTYPE2 max_small_magnitude);

    /// @brief Compute the polygon triangulation that maximizes the minimum angle.
    /// - Store the triangulation information in polygon ipoly.
    /// - Version using C++ STL vector vcoord[].
    template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
    void ComputeCosMaxMinAngle
    (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly, 
     const MTYPE2 max_small_magnitude);

    /// @brief Compute the polygon triangulation that maximizes the minimum angle.
    /// - Store the triangulation information in polygon ipoly.
    /// - Version that uses POLYGON_TRIANGULATION_SETTINGS.
    ///   to specify allowable trianguations.
    template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
    void ComputeCosMaxMinAngle
    (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
     const MTYPE2 max_small_magnitude,
     const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings);

    /// @brief Compute the polygon triangulation that maximizes the minimum angle.
    /// - Store the triangulation information in polygon ipoly.
    /// - Version that uses POLYGON_TRIANGULATION_SETTINGS.
    ///   to specify allowable trianguations.
    /// - Version using C++ STL vector vcoord[].
    template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
    void ComputeCosMaxMinAngle
    (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly, 
     const MTYPE2 max_small_magnitude,
     const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings);

    /// @brief Compute the polygon triangulation that maximizes the minimum angle.
    /// - Consider triangulations from a single polygon vertex
    ///     or from the polygon centroid.
    /// - Store the triangulation information in polygon ipoly.
    template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
    void ComputeCosMaxMinAngleAllowCentroid
    (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
     const MTYPE2 max_small_magnitude);

    /// @brief Compute triangulations of each polygon, maximizing the minimum angle in each triangulation.
    /// - Store the triangulation information in polygon ipoly.
    template <typename CTYPE2, typename MTYPE2>
    void ComputeAllPolygonsCosMaxMinAngle
    (const CTYPE2 vcoord[], const MTYPE2 max_small_magnitude);

    /// @brief Compute triangulations of each polygon, maximizing the minimum angle in each triangulation.
    ///  - Version using C++ STL vector vertex_coord[].
    template <typename CTYPE2, typename MTYPE2>
    void ComputeAllPolygonsCosMaxMinAngle
    (const std::vector<CTYPE2> & vcoord, const MTYPE2 max_small_magnitude);

    /// @brief Compute triangulations of each polygon, maximizing the minimum angle in each triangulation.
    /// - Consider triangulations from a single polygon vertex
    ///     or from the polygon centroid.
    /// - Store the triangulation information in polygon ipoly.
    template <typename CTYPE2, typename MTYPE2>
    void ComputeAllPolygonsCosMaxMinAngleAllowCentroid
    (const CTYPE2 vcoord[], const MTYPE2 max_small_magnitude);

    /// @brief Compute triangulations of each polygon, maximizing the minimum angle in each triangulation.
    /// - Consider triangulations from a single polygon vertex
    ///     or from the polygon centroid.
    ///  - Version using C++ STL vector vertex_coord[].
    template <typename CTYPE2, typename MTYPE2>
    void ComputeAllPolygonsCosMaxMinAngleAllowCentroid
    (const std::vector<CTYPE2> & vcoord, const MTYPE2 max_small_magnitude);


    // Triangulate polygons.
    /*!
     *  @brief Triangulate polygon ipoly based on triangulation encoding.
     */
    template <typename ITYPE2, int BIT_SIZE_TYPE_X>
    void TriangulatePolygon
    (const ITYPE2 ipoly, 
     const POLYGON_TRIANGULATION_ENCODING<BIT_SIZE_TYPE_X> & 
     triangulation_encoding);

    /*!
     *  @brief Triangulate polygon ipoly.
     *  - Use polygon triangulation information to determine triangulation.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on ipoly to determine triangulation.
     */
    template <typename ITYPE2, typename CTYPE2>
    void TriangulatePolygon
    (const ITYPE2 ipoly, std::vector<CTYPE2> & vcoord);

    /*!
     *  @brief Triangulate all polygons.
     *  - Use polygon triangulation information to determine triangulation.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on all undeleted polygons with 4 or more vertices
     *     to determine triangulation.
     */
    template <typename CTYPE2>
    void TriangulateAllPolygons(std::vector<CTYPE2> & vcoord);


    // Compute edge/polygon characteristics.

    /// @brief Return true if edge is "long".
    /// - Return true if edge length squared is greater than or equal to L.
    template <typename CTYPE2, typename ITYPE2, typename LTYPE2>
    bool IsLongEdge(const std::vector<CTYPE2> & vcoord,
                    const ITYPE2 ihalf_edge, const LTYPE2 L) const;

    /// @brief Return true if quad is "long".
    /// @param[out] ilongest_half_edge Longest half edge of quad.
    template <typename CTYPE2, typename ITYPE2, 
              typename ITYPE3, typename RTYPE>
    bool IsLongQuad(const std::vector<CTYPE2> & vcoord,
                    const ITYPE2 iquad, const RTYPE Rsquared,
                    ITYPE3 & ilongest_half_edge) const;

    /// @brief Return true if shortest edge is between two long edges.
    /// @param[out] ishortest_half_edge Shortest half edge of quad.
    template <typename CTYPE2, typename ITYPE2, 
              typename ITYPE3, typename RTYPE>
    bool DoesQuadHaveLongShortLongEdges
    (const std::vector<CTYPE2> & vcoord,
     const ITYPE2 iquad, const RTYPE Rsquared,
     ITYPE3 & ishortest_half_edge) const;

    /// @brief Return true if quad has one edge much shorter than all other edges.
    /// @param[out] ishortest_half_edge Shortest half edge of quad.
    template <typename CTYPE2, typename ITYPE2, 
              typename ITYPE3, typename RTYPE>
    bool DoesQuadHaveVeryShortEdge
    (const std::vector<CTYPE2> & vcoord,
     const ITYPE2 ipolyg, const RTYPE Rsquared,
     ITYPE3 & ishortest_half_edge) const;

    /// @brief Return index of shortest half edge in polygon.
    /// @param[out] shortest_length_squared Length squared of shortest edge.
    template <typename CTYPE2, typename ITYPE2, typename LTYPE2>
    HALF_EDGE_INDEX_TYPE ComputeShortestHalfEdge
    (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly,
     LTYPE2 & shortest_length_squared) const;

    /// @brief Get index of shortest half edge in polygon.
    /// - Version that does not return shortest_length_squared.
    template <typename CTYPE2, typename ITYPE2>
    HALF_EDGE_INDEX_TYPE ComputeShortestHalfEdge
    (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly) const;

    /// @brief Get index of shortest non-split half edge in polygon.
    /// - Neither endpoint of the half edge is a split vertex.
    /// - Return false if no such edge is found.
    /// @param[out] shortest_length_squared Length squared of shortest edge.
    template <typename CTYPE2, typename ITYPE2, typename ITYPE3, 
              typename LTYPE2>
    bool ComputeShortestNonSplitHalfEdge
    (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly, 
     ITYPE3 & ishortest,
     LTYPE2 & shortest_length_squared) const;

    /// Return true if polygon is a triangle with one split edge.
    template <typename ITYPE2>
    bool IsTriangleWithOneSplitEdge(const ITYPE2 ipoly) const;

    /// @brief Return true if polygon is a triangle with two split edges.
    /// @param[out] If true, then ihalf_edge0 is non-split triangle edge.
    template <typename ITYPE2, typename ITYPE3>
    bool IsTriangleWithTwoSplitEdges
    (const ITYPE2 ipoly, ITYPE3 & ihalf_edge0) const;


    /// Return number of vertices of a given type in a polygon.
    template <typename ITYPE2>
    NUMBER_TYPE CountNumPolygonVerticesOfType
    (const ITYPE2 ipoly, const TRIANGULATION_VERTEX_TYPE vtype) const;

    /// Return number of input vertices in a polygon.
    template <typename ITYPE2>
    NUMBER_TYPE CountNumInputPolygonVertices(const ITYPE2 ipoly) const
    { return (CountNumPolygonVerticesOfType
              (ipoly, TRIANGULATION_VERTEX_TYPE::TRIV_INPUT)); }

    /// Return number of split edge vertices in a polygon.
    template <typename ITYPE2>
    NUMBER_TYPE CountNumSplitEdgePolygonVertices(const ITYPE2 ipoly) const
    { return (CountNumPolygonVerticesOfType
              (ipoly, TRIANGULATION_VERTEX_TYPE::TRIV_SPLITE)); }


    /// @brief Count and set number of split edge vertices in a polygon.
    /// - Stores number in data structure.
    template <typename ITYPE2>
    void SetNumSplitEdgePolygonVertices(const ITYPE2 ipoly);


    // Set triangulation information.

    /// @brief Set triangulation vertex index in triangulation_info.
    /// - Triangulate polygon by adding diagonals from tri_vertex_index.
    template <typename ITYPE2, typename ITYPE3>
    void SetTriVertexIndex
    (const ITYPE2 ipoly, const ITYPE3 tri_vertex_index)
    { this->polygon_list[ipoly].triangulation_info.tri_vertex_index =
        tri_vertex_index; }

    /// Set cos min triangulation angle.
    template <typename ITYPE2, typename COS_TYPE2>
    void SetCosMinTriangulationAngle
    (const ITYPE2 ipoly, const COS_TYPE2 cos_min_angle)
    { this->polygon_list[ipoly].triangulation_info.cos_min_triangulation_angle =
        cos_min_angle; }

    /// Set triangulation to include one interior vertex.
    template <typename ITYPE2, typename COS_TYPE2, typename CTYPE2>
    void SetInteriorTriVertex
    (const ITYPE2 ipoly, const COS_TYPE2 cos_min, 
     const CTYPE2 interior_coord[]);

    /// @brief Set triangulation to include one interior vertex.
    /// - Version using C++ STL vector interior_coord[].
    template <typename ITYPE2, typename COS_TYPE2, typename CTYPE2>
    void SetInteriorTriVertex
    (const ITYPE2 ipoly, const COS_TYPE2 cos_min, 
     const std::vector<CTYPE2> & interior_coord);

    /// @brief Set triangulation information with 0 interior vertices.
    /// - Set num_interior_tri_vertices to 0.
    /// - Set flag_zero to false.
    /// - Set number_triangles to NumPolygonVertices(ipoly)-2.
    template <typename ITYPE2, typename ITYPE3, typename COS_TYPE2>
    void SetTriInfoNoInterior
    (const ITYPE2 ipoly, const COS_TYPE2 cos_min, const ITYPE3 triv_index)
    {
      this->polygon_list[ipoly].triangulation_info.SetNoInterior
        (cos_min, this->NumPolygonVertices(ipoly)-2, triv_index);
    }

    /// @brief Copy triangulation information and interior_coord.
    /// - Copy triangulation from poly_tri_info and interior_coord
    ///   into triangulation info for polygon ipoly.
    template <typename INFO_TYPE, typename CTYPE2, typename ITYPE2>
    void CopyTriInfoAndInteriorCoord
    (const INFO_TYPE & poly_tri_info, 
     const CTYPE2 interior_coord[],
     const ITYPE2 ipoly);

    /// @brief Copy triangulation information and interior_coord.
    /// - Copyt triangulation from poly_tri_info and interior_coord
    ///   into triangulation info for polygon ipoly.
    /// - Version using C++ STL vector interior_coord[].
    template <typename INFO_TYPE, typename CTYPE2, typename ITYPE2>
    void CopyTriInfoAndInteriorCoord
    (const INFO_TYPE & poly_tri_info, 
     const std::vector<CTYPE2> & interior_coord,
     const ITYPE2 ipoly);

    /// @brief Copy triangulation information.
    ///   - Copy triangulation from poly_tri_info
    ///     into triangulation info for polygon ipoly.
    template <typename INFO_TYPE, typename ITYPE2>
    void CopyTriInfo(const INFO_TYPE & poly_tri_info, const ITYPE2 ipoly);

    /// @brief Set triangulation information of polygon with split edge.
    /// @param ihalf_edge0 Split edge in original polygon.
    /// @param poly0_tri_info Triangulation information of original polygon.
    /// @param ipoly1 New polygon containing split vertex.
    template <typename ITYPE2, typename ITYPE3, typename INFO_TYPE>
    void SetSplitEdgePolygonTriInfo
    (const ITYPE2 ihalf_edge0, const INFO_TYPE & poly0_tri_info,
     const ITYPE3 ipoly1);

    /// @brief Set triangulation information of polygon with split edge.
    /// - Allow triangulation with interior triangulation vertices.
    /// @param ihalf_edge0 Split edge in original polygon.
    /// @param poly0_tri_info Triangulation information of original polygon.
    /// @param interior_coord[] Coordinates of interior triangulation vertex.
    /// @param ipoly1 New polygon containing split vertex.
    template <typename ITYPE2, typename ITYPE3, typename INFO_TYPE,
              typename CTYPE2>
    void SetSplitEdgePolygonTriInfoAllowIV
    (const ITYPE2 ihalf_edge0, const INFO_TYPE & poly0_tri_info,
     const CTYPE2 interior_coord[], const ITYPE3 ipoly1);

    /// @brief Set triangulation information of polygon with split edge.
    /// - Allow triangulation with interior triangulation vertices.
    /// - Version using C++ STL vector interior_coord[].
    /// @param ihalf_edge0 Split edge in original polygon.
    /// @param poly0_tri_info Triangulation information of original polygon.
    /// @param interior_coord[] Coordinates of interior triangulation vertex.
    /// @param ipoly1 New polygon containing split vertex.
    template <typename ITYPE2, typename ITYPE3, typename INFO_TYPE,
              typename CTYPE2>
    void SetSplitEdgePolygonTriInfoAllowIV
    (const ITYPE2 ihalf_edge0, const INFO_TYPE & poly0_tri_info,
     const std::vector<CTYPE2> & interior_coord, const ITYPE3 ipoly1);


    // Delete edges.

    /// @brief Delete edge to improve max min triangulation.
    /// - Return true if edge is deleted.
    /// @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
    ///    has been called on ipoly to determine triangulation.
    template <typename ITYPE2, typename MTYPE2, typename CTYPE2>
    bool DeleteEdgeToMaxMinTriangulationAngle
    (const ITYPE2 ihalf_edge, const MTYPE2 max_small_magnitude,
     std::vector<CTYPE2> & vcoord);

    /// @brief Delete longest polygon edge to improve max min triangulation.
    /// - Return true if edge is split.
    /// @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
    ///    has been called on ipoly to determine triangulation.
    template <typename ITYPE2, typename MTYPE2, typename CTYPE2>
    bool DeleteLongestEdgeToMaxMinTriangulationAngle
    (const ITYPE2 ipoly, const MTYPE2 max_small_magnitude, 
     std::vector<CTYPE2> & vcoord);

    /// @brief Delete (long) triangle edges to improve max min triangulation.
    /// @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
    ///    has been called on all polygons to determine triangulation.
    /// @param cos_angle_threshold Attempt to delete longest triangle edge
    ///    if cos_min_triangulation_angle is above cos_angle_threshold.
    template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
    void DeleteTriangleEdgesToMaxMinTriangulationAngle
    (const COS_TYPE2 cos_angle_threshold, const MTYPE2 max_small_magnitude,
     std::vector<CTYPE2> & vcoord);

    /// @brief Delete (long) triangle edges to improve max min triangulation.
    /// - Version with array cos_angle_threshold[].
    template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
    void DeleteTriangleEdgesToMaxMinTriangulationAngle
    (const std::vector<COS_TYPE2> & cos_angle_threshold, 
     const MTYPE2 max_small_magnitude,
     std::vector<CTYPE2> & vcoord);

    /// @brief Delete edge and add triangle to improve max min triangulation.
    /// - Return true if edge is deleted.
    /// - Variation of delete edge that cuts off an ear of the quad
    ///     containing ihalf_edge.
    /// - Does nothing if iear is an endpoint of ihalf_edgeA.
    /// @pre Polygon containing ihalf_edgeA is a quadrilateral.
    /// @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
    ///    has been called on ipoly to determine triangulation.
    /// @param ihalf_edgeA Edge to be deleted.
    /// @param iear Index (location, 0,1,2 or 3) of vertex cut off by ear.
    ///    Added triangle is iear and vertices preceding and following iear.
    template <typename ITYPE2, typename ITYPE3,
              typename MTYPE2, typename CTYPE2>
    bool DeleteEdgeAddTriangleToMaxMinTriangulationAngle
    (const ITYPE2 ihalf_edgeA, const ITYPE3 iear, 
     const MTYPE2 max_small_magnitude, std::vector<CTYPE2> & vcoord);

    /// @brief Delete (long) quad edges to improve max min triangulation.
    /// - Calls DeleteEdgeAddTriangleToMaxMinTriangulationAngle().
    /// @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
    ///    has been called on all polygons to determine triangulation.
    /// @param cos_angle_threshold Attempt to delete longest quad edge (u,v)
    ///    if quad angle at u or at v is less than 
    ///    cos^{-1}(cos_angle_threshold).
    template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
    void DeleteQuadEdgesToMaxMinTriangulationAngle
    (const COS_TYPE2 cos_angle_threshold, const MTYPE2 max_small_magnitude,
     std::vector<CTYPE2> & vcoord);

    /// @brief Delete (long) quad edges to improve max min triangulation.
    /// - Calls DeleteEdgeAddTriangleToMaxMinTriangulationAngle().
    /// - Version with array cos_angle_threshold[].
    template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
    void DeleteQuadEdgesToMaxMinTriangulationAngle
    (const std::vector<COS_TYPE2> & cos_angle_threshold, 
     const MTYPE2 max_small_magnitude,
     std::vector<CTYPE2> & vcoord);


    // Replace triangle.

    /*!
     *  @brief Replace triangle with triangle-edge to max min triangulation angle.
     *  - Return true if triangle is replaced.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on ipoly to determine triangulation.
     *  @param ihalf_edge Triangle base.
     *    - Replace vertex opposite triangle base.
     *  @param flag_max_min_all If true, replace triangle if it
     *      improves the min triangulation of the triangle and
     *      its two neighbors.
     *    - Otherwise, replace triangle only if it increase the min
     *      triangulation angle of the triangle.
     */
    template <typename ITYPE2, typename MTYPE2, typename CTYPE2>
    bool ReplaceTriangleWithTriangleEdgeToMaxMinTriangulationAngle
    (const ITYPE2 ihalf_edge, const MTYPE2 max_small_magnitude,
     const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
     const bool flag_allow_centroid2,
     const bool flag_max_min_all, std::vector<CTYPE2> & vcoord);




    /*!
     *  @brief Replace triangle that has two split edges with triangle-edg to max min triangulation angle.
     *  - Return true if triangle is replaced.
     *  - Allow triangulations from polygon centroids.
     *  - Consider two possible centroids for each adjacent polygon.
     *    Either centroid of polygon or centroid of polygon and the new
     *    triangle vertex.
     *  @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
     *     has been called on ipoly to determine triangulation.
     *  @param ihalf_edge Triangle base.
     *    - Replace vertex opposite triangle base.
     *  @param flag_max_min_all If true, replace triangle if it
     *      improves the min triangulation of the triangle and
     *      its two neighbors.
     *    - Otherwise, replace triangle only if it increase the min
     *      triangulation angle of the triangle.
     */
    template <typename ITYPE2, typename MTYPE2, typename CTYPE2>
    bool 
    ReplaceSplitEdge2TriangleToMaxMinTriangulationAngle
    (const ITYPE2 ihalf_edge, const MTYPE2 max_small_magnitude,
     const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
     const bool flag_allow_centroidx2,
     const bool flag_max_min_all, std::vector<CTYPE2> & vcoord);

    /*!
     *  @brief Replace triangles that have two split edges to max min triangulation angle.
     *  @param flag_allow_centroid  If true, allow adding triangulation
     *      vertices at centroids of adjacent polygons.
     *  @param flag_triangle_has_min_angle  If true, replace only if
     *      min triangle angle is less than min angle in adjacent polygons.
     *  @param flag_max_min_all If true, replace triangle if it
     *      improves the min triangulation of the triangle and
     *      its two neighbors.
     *    - Otherwise, replace triangle only if it increase the min
     *      triangulation angle of the triangle.
     */
    template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
    void ReplaceSplitEdge2TrianglesToMaxMinTriangulationAngle
    (const COS_TYPE2 cos_angle_threshold, 
     const MTYPE2 max_small_magnitude,
     const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
     const bool flag_allow_centroidx2,
     const bool flag_triangle_has_min_angle,
     const bool flag_max_min_all, 
     std::vector<CTYPE2> & vcoord);

    /*!
     *  @brief Replace triangles that have two split edges to improve
     *    max min triangulation angle.
     *  - Version with array cos_angle_threshold[].
    */
    template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
    void ReplaceSplitEdge2TrianglesToMaxMinTriangulationAngle
    (const std::vector<COS_TYPE2> & cos_angle_threshold, 
     const MTYPE2 max_small_magnitude,
     const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
     const bool flag_allow_centroidx2,
     const bool flag_triangle_has_min_angle,
     const bool flag_max_min_all, 
     std::vector<CTYPE2> & vcoord);


    // Check routines.

    /*!
     *  @brief Check that two adjacent polygons are different.
     *  - Return true if two adjacent polygons are different.
     *  - Return true if ihalf_edgeA or ihalf_edgeB are boundary edges.
    */
    template <typename ITYPEH>
    bool CheckAreAdjacentPolygonsDifferent
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     IJK::ERROR & error) const;

    /*!
     *  @brief Check that three adjacent polygons are different.
     *  - Return true if three adjacent polygons are different.
     *  - Return true if two or three edges are boundary edges.
     *  - Return true if one edge is a boundary edge and the other
     *    two adjacent polygons are different.
    */
    template <typename ITYPEH>
    bool CheckAreAdjacentPolygonsDifferent
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     const ITYPEH ihalf_edgeC, IJK::ERROR & error) const;

  };


  // *****************************************************************
  // Class MESH2D_TRIANGULATION_A
  // *****************************************************************

  /*!
   *  @brief Mesh of 2D polygons embedded for computing triangulations.
   *  - Use VERTEX_TRI_INFO_BASE, HALF_EDGE_BASE, POLYGON_TRI_INFO_BASE
   *   to represent vertices, half edges and polygons.
   */
  template <int DIMENSION, int BIT_SET_SIZE,
            typename COORD_TYPE, typename COS_TYPE,
            typename ITYPE, typename NTYPE>
  class MESH2D_TRIANGULATION_A:
    public MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,
   VERTEX_TRI_INFO_BASE<ITYPE,ITYPE,NTYPE>, 
   HALF_EDGE_BASE<ITYPE,ITYPE,ITYPE>,
   POLYGON_TRI_INFO_BASE
   <DIMENSION, BIT_SET_SIZE, COORD_TYPE, COS_TYPE, ITYPE, NTYPE>,
   ITYPE, NTYPE>
  {
  public:
    MESH2D_TRIANGULATION_A() {};
  };


  // *****************************************************************
  // MESH2D triangulation routines
  // *****************************************************************

  /*!
   *  @brief Triangulate polygon ipoly from a single vertex, maximizing the minimum angle.
   *  - Return triangulation information.
   */
  template <typename DTYPE, typename CTYPE, typename MESH_TYPE,
            typename ITYPE0, typename ITYPE1, typename MTYPE,
            typename COS_TYPE>
  void triangulate_polygon_max_min_angle
  (const DTYPE dimension, const CTYPE vcoord[], const ITYPE0 ipoly, 
   const MTYPE max_small_magnitude,
   MESH_TYPE & mesh, 
   ITYPE1 & tri_vertex_index, COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_fan_triangulation_to_max_min_angle_M
      (dimension, vcoord, mesh, ipoly, max_small_magnitude,
       tri_vertex_index, cos_min_angle, flag_zero);

    mesh.TriangulatePolygonFromVertex(ipoly, tri_vertex_index);
  }


  /*!
   *  @brief Triangulate polygon ipoly from a single vertex, maximizing the minimum angle.
   *  - Version that does not return triangulation information.
   */
  template <typename DTYPE, typename CTYPE, typename MESH_TYPE,
            typename ITYPE0, typename MTYPE>
  void triangulate_polygon_max_min_angle
  (const DTYPE dimension, const CTYPE vcoord[], const ITYPE0 ipoly, 
   const MTYPE max_small_magnitude, MESH_TYPE & mesh)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    CTYPE cos_min_angle;
    NUMBER_TYPE tri_vertex_index;
    bool flag_zero;

    triangulate_polygon_max_min_angle
      (dimension, vcoord, ipoly, max_small_magnitude, mesh,
       tri_vertex_index, cos_min_angle, flag_zero);
  }


  /*!
   *  @brief Triangulate polygon ipoly, maximizing the minimum angle.
   *  - Consider triangulations from a single polygon vertex
   *      or from an interior vertex.
   *  - Return triangulation information.
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename MESH_TYPE,
            typename ITYPE0, typename ITYPE1, typename NTYPE,
            typename MTYPE, typename COS_TYPE>
  void triangulate_polygon_max_min_angle_allow_interior_vertex
  (const DTYPE dimension, const COORD_TYPE1 interior_vcoord[], 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   std::vector<COORD_TYPE0> & vcoord, 
   MESH_TYPE & mesh, COS_TYPE & cos_min_angle, 
   NTYPE & num_interior, ITYPE1 & tri_vertex_index, bool & flag_zero)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

    compute_cos_max_min_polygon_triangulation_angle_allow_interior_vertex
      (dimension, vcoord, interior_vcoord, mesh, ipoly, max_small_magnitude,
       cos_min_angle, num_interior, tri_vertex_index, flag_zero);

    if (num_interior == 0) {
      mesh.TriangulatePolygonFromVertex(ipoly, tri_vertex_index);
    }
    else {
      
      const VERTEX_INDEX_TYPE
        ivX = IJK::insert_coord(dimension, interior_vcoord, vcoord);
      mesh.AddVertex();

      mesh.TriangulatePolygonFromInteriorVertex(ipoly, ivX);
    }
  }


  /*!
   *  @brief Triangulate polygon ipoly, maximizing the minimum angle.
   *  - Consider triangulations from a single polygon vertex
   *      or from an interior vertex.
   *  - Version that does not return triangulation information.
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename ITYPE0,  typename MTYPE,
            typename MESH_TYPE>
  void triangulate_polygon_max_min_angle_allow_interior_vertex
  (const DTYPE dimension, const COORD_TYPE1 interior_vcoord[], 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   std::vector<COORD_TYPE0> & vcoord, MESH_TYPE & mesh)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    COORD_TYPE1 cos_min_angle;
    NUMBER_TYPE num_interior;
    NUMBER_TYPE tri_vertex_index;
    bool flag_zero;

    triangulate_polygon_max_min_angle_allow_interior_vertex
      (dimension, interior_vcoord, ipoly, max_small_magnitude,
       vcoord, mesh, cos_min_angle, num_interior, 
       tri_vertex_index, flag_zero);
  }


  /*!
   *  @brief Triangulate polygon ipoly, maximizing the minimum angle.
   *  - Consider triangulations from a single polygon vertex
   *      or from an interior vertex.
   *  - Version that does not return triangulation information.
   *  - Version using C++ STL vector interior_vcoord[].
   */
  template <typename DTYPE, typename COORD_TYPE0, typename COORD_TYPE1,
            typename ITYPE0,  typename MTYPE,
            typename MESH_TYPE>
  void triangulate_polygon_max_min_angle_allow_interior_vertex
  (const DTYPE dimension, const std::vector<COORD_TYPE1> & interior_vcoord, 
   const ITYPE0 ipoly, const MTYPE max_small_magnitude,
   std::vector<COORD_TYPE0> & vcoord, MESH_TYPE & mesh)
  {
    triangulate_polygon_max_min_angle_allow_interior_vertex
      (dimension, IJK::vector2pointer(interior_vcoord), ipoly,
       max_small_magnitude, vcoord, mesh);
  }


  /*!
   *  @brief Triangulate all polygons in mesh, maximizing the minimum angle.
   *  - Triangulate each polygon from a single vertex.
   */
  template <typename DTYPE, typename CTYPE, typename MTYPE,
            typename MESH_TYPE>
  void triangulate_mesh_max_min_angle
  (const DTYPE dimension, const CTYPE vcoord[], 
   const MTYPE max_small_magnitude, MESH_TYPE & mesh)
  {
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < mesh.NumPolygons();
         ipoly++) {

      if (mesh.IsTriangle(ipoly)) { continue; }
      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      triangulate_polygon_max_min_angle
        (dimension, vcoord, ipoly, max_small_magnitude, mesh);
    }
  }


  /*!
   *  @brief Triangulate all polygons in mesh, maximizing the minimum angle.
   *  - Triangulate each polygon from a single vertex.
   *  - Version using C++ STL vector vertex_coord[].
   */
  template <typename DTYPE, typename CTYPE, typename MTYPE,
            typename MESH_TYPE>
  void triangulate_mesh_max_min_angle
  (const DTYPE dimension, const std::vector<CTYPE> & vcoord, 
   const MTYPE max_small_magnitude, MESH_TYPE & mesh)
  {
    triangulate_mesh_max_min_angle
      (dimension, IJK::vector2pointer(vcoord), 
       max_small_magnitude, mesh);
  }


  /*!
   *  @brief Triangulate all polygons in mesh, maximizing the minimum angle.
   *  - Consider triangulations from a single polygon vertex
   *    or from the polygon centroid.
   */
  template <typename DTYPE, typename CTYPE, typename MTYPE,
            typename MESH_TYPE>
  void triangulate_mesh_max_min_angle_allow_centroid
  (const DTYPE dimension, const MTYPE max_small_magnitude, 
   std::vector<CTYPE> & vcoord, MESH_TYPE & mesh)
  {
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

    std::vector<CTYPE> centroid_coord(dimension);

    for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < mesh.NumPolygons();
         ipoly++) {

      if (mesh.IsTriangle(ipoly)) { continue; }
      if (mesh.IsPolygonDeleted(ipoly)) { continue; }

      ComputePolygonCentroid(vcoord, ipoly, centroid_coord);

      triangulate_polygon_max_min_angle_allow_interior_vertex
        (dimension, centroid_coord, ipoly, max_small_magnitude, 
         vcoord, mesh);
    }
  }


  // *****************************************************************
  // Class POLYGON_TRI_INFO_BASE member functions.
  // *****************************************************************

  template <int DIMENSION, int BIT_SET_SIZE,
            typename COORD_TYPE_X, typename COS_TYPE_X,
            typename ITYPE, typename NTYPE>
  void POLYGON_TRI_INFO_BASE
  <DIMENSION,BIT_SET_SIZE,COORD_TYPE_X,COS_TYPE_X,ITYPE,NTYPE>::
  Init()
  {
    this->num_split_edge_vertices = 0;
  }


  // *****************************************************************
  // MACROS FOR TEMPLATE MEMBER FUNCTIONS MESH2D_TRIANGULATION_BASE
  // *****************************************************************

  /// Template parameters for template MESH2D_TRIANGULATION_BASE
# define _MESH2D_TRIANGULATION_BASE_T_PARAM_                          \
  <int BIT_SET_SIZE,                                                  \
   typename VTYPE, typename HALFE_TYPE, typename PTYPE,               \
   typename ITYPE, typename NTYPE>

  /// Class name for template MESH2D_TRIANGULATION_BASE
# define _MESH2D_TRIANGULATION_BASE_T_NAME_                           \
  MESH2D_TRIANGULATION_BASE                                           \
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>



  // *****************************************************************
  // Class MESH2D_TRIANGULATION_BASE add vertices.
  // *****************************************************************

  // Return cosine of the min angle of the triangulation
  //   specified for polygon ipoly and its neighbors.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2>
  typename _MESH2D_TRIANGULATION_BASE_T_NAME_::COS_TYPE 
  _MESH2D_TRIANGULATION_BASE_T_NAME_::
  PolygonAndNeighborsCosMinTriangulationAngle(const ITYPE2 ipoly) const
  {
    COS_TYPE cos_poly_angle = 
      PolygonCosMinTriangulationAngle(ipoly);

    for (NUMBER_TYPE j = 0; j < this->NumPolygonEdges(ipoly); j++) {
      const HALF_EDGE_INDEX_TYPE jhalf_edge = 
        this->HalfEdgeIndex(ipoly, j);
      const POLYGON_INDEX_TYPE kpoly =
        this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(jhalf_edge);

      cos_poly_angle = 
        std::max(cos_poly_angle, PolygonCosMinTriangulationAngle(kpoly));
    }

    return(cos_poly_angle);
  }


  // Add interior vertex.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename CTYPE3>
  typename _MESH2D_TRIANGULATION_BASE_T_NAME_::VERTEX_INDEX_TYPE
  _MESH2D_TRIANGULATION_BASE_T_NAME_::AddInteriorVertex
  (const CTYPE2 interior_vcoord[], std::vector<CTYPE3> & vertex_coord)
  {
    const VERTEX_INDEX_TYPE iv_interior =
      IJK::insert_coord(Dimension(), interior_vcoord, vertex_coord);
    this->AddVertex();

    this->vertex_list[iv_interior].triangulation_vertex_type = 
      VTYPE::TRIV_INTERIOR;

    return(iv_interior);
  }


  // Add interior vertex.
  // - Version using C++ STL vector interior_vcoord[].
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename CTYPE3>
  typename _MESH2D_TRIANGULATION_BASE_T_NAME_::VERTEX_INDEX_TYPE
  _MESH2D_TRIANGULATION_BASE_T_NAME_::AddInteriorVertex
  (const std::vector<CTYPE2> & interior_vcoord, 
   std::vector<CTYPE3> & vertex_coord)
  {
    return(AddInteriorVertex(IJK::vector2pointer(interior_vcoord), 
                             vertex_coord));
  }


  // Add split vertex.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename CTYPE3>
  typename _MESH2D_TRIANGULATION_BASE_T_NAME_::VERTEX_INDEX_TYPE
  _MESH2D_TRIANGULATION_BASE_T_NAME_::AddSplitVertex
  (const CTYPE2 split_vcoord[], std::vector<CTYPE3> & vertex_coord)
  {
    const VERTEX_INDEX_TYPE iv_split =
      IJK::insert_coord(Dimension(), split_vcoord, vertex_coord);
    this->AddVertex();

    this->vertex_list[iv_split].triangulation_vertex_type = 
      VTYPE::TRIV_SPLITE;

    return(iv_split);
  }

  // Add split vertex.
  // - Version using C++ STL vector split_vcoord[].
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename CTYPE3>
  typename _MESH2D_TRIANGULATION_BASE_T_NAME_::VERTEX_INDEX_TYPE
  _MESH2D_TRIANGULATION_BASE_T_NAME_::AddSplitVertex
  (const std::vector<CTYPE2> & split_vcoord, 
   std::vector<CTYPE3> & vertex_coord)
  {
    return(AddSplitVertex(IJK::vector2pointer(split_vcoord), 
                          vertex_coord));
  }


  // *****************************************************************
  // Class MESH2D_TRIANGULATION_BASE split functions.
  // *****************************************************************

  // Add split vertex and split polygon edge.
  // - Set polygon triangulation information.
  // - Allow polygon triangulations using interior vertices.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename CTYPE3,
            typename CTYPE4, typename CTYPE5, 
            typename ITYPE2, typename ITYPE3, 
            typename ITYPE4, typename ITYPE5,
            typename COS_TYPE2, typename NTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  AddVertexSplitPolygonEdgeIIAllowIV
  (const std::vector<CTYPE2> & split_vcoord, 
   const std::vector<CTYPE3> & interior_coordA, 
   const std::vector<CTYPE4> & interior_coordB, 
   std::vector<CTYPE5> & vertex_coord,
   const ITYPE2 ihalf_edge0, 
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE2,NTYPE2> & polyA_tri_info,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE2,NTYPE2> & polyB_tri_info,
   ITYPE3 & iv_split, ITYPE4 & ipolyA, ITYPE5 & ipolyB)
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edge0X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge0);

    MESH2D_SPLIT_II_BASE_TYPE::AddVertexSplitPolygonEdgeII
      (Dimension(), split_vcoord, vertex_coord, ihalf_edge0, 
       iv_split, ipolyA, ipolyB);
    SetVertexType(iv_split, VTYPE::TRIV_SPLITE);

    SetNumSplitEdgePolygonVertices(ipolyA);
    SetNumSplitEdgePolygonVertices(ipolyB);

    SetSplitEdgePolygonTriInfoAllowIV
      (ihalf_edge0, polyA_tri_info, interior_coordA, ipolyA);
    SetSplitEdgePolygonTriInfoAllowIV
      (ihalf_edge0X, polyB_tri_info, interior_coordB, ipolyB);
  }


  // Add split vertex and split polygon edge.
  // - Set polygon triangulation information.
  // - Allow polygon triangulations using interior vertices.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3,
            typename CTYPE4, typename CTYPE5, 
            typename ITYPE2, typename COS_TYPE2, typename NTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  AddVertexSplitPolygonEdgeIIAllowIV
  (const std::vector<CTYPE2> & split_vcoord, 
   const std::vector<CTYPE3> & interior_coordA, 
   const std::vector<CTYPE4> & interior_coordB, 
   std::vector<CTYPE5> & vertex_coord,
   const ITYPE2 ihalf_edge0, 
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE2,NTYPE2> & polyA_tri_info,
   POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE2,NTYPE2> & polyB_tri_info)
  {
    VERTEX_INDEX_TYPE iv_split;
    POLYGON_INDEX_TYPE ipolyA, ipolyB;

    AddVertexSplitPolygonEdgeIIAllowIV
      (split_vcoord, interior_coordA, interior_coordB,
       vertex_coord, ihalf_edge0, polyA_tri_info, polyB_tri_info,
       iv_split, ipolyA, ipolyB);
  }


  // Add new triangle vertex and replace triangle with triangle edge.
  // - Set polygon triangulation information.
  // - Allow polygon triangulations using interior vertices.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, 
            typename CTYPE4, typename CTYPE5,
            typename ITYPE2, typename ITYPE3, typename ITYPE4,
            typename RESULT_TYPEA, typename RESULT_TYPEB,
            typename RESULT_TYPEC>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  AddVertexReplaceTriangleWithTriangleEdge
    (const std::vector<CTYPE2> & new_triangle_vcoord,
     const std::vector<CTYPE3> & interior_vcoord0,
     const std::vector<CTYPE4> & interior_vcoord2,
     std::vector<CTYPE5> & vertex_coord,
     const ITYPE2 ihalf_edge0, 
     const RESULT_TYPEA & polyA_tri_info,
     const RESULT_TYPEB & triangleB_tri_info,
     const RESULT_TYPEC & polyC_tri_info,
     ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0)
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edge1 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge1X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge2);

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices
      (cerr, "*** Replacing triangle: ",
       this->IndexOfPolygonContainingHalfEdge(ihalf_edge1), "\n");
    this->PrintPolygonIndexAndVertices
      (cerr, "  polyA: ",
       this->IndexOfPolygonContainingHalfEdge(ihalf_edge1X), "\n");
    this->PrintPolygonIndexAndVertices
      (cerr, "  polyC: ",
       this->IndexOfPolygonContainingHalfEdge(ihalf_edge2X), "\n");
    IJK::print_coord3D
      (cerr, "    interior_vcoord0: ", interior_vcoord0, "\n");
    IJK::print_coord3D
      (cerr, "    interior_vcoord0: ", interior_vcoord2, "\n");
    */
       

    MESH2D_SPLIT_II_BASE_TYPE::AddVertexReplaceTriangleWithTriangleEdge
      (Dimension(), new_triangle_vcoord, vertex_coord, 
       ihalf_edge0, iv_new, ihalf_edge_new0);
    SetVertexType(iv_new, VTYPE::TRIV_SPLITE);

    const HALF_EDGE_INDEX_TYPE ihalf_edge_new1 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge_new2 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new1);
    const POLYGON_INDEX_TYPE ipolyA =
      this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge_new1);
    const POLYGON_INDEX_TYPE itriangleB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge_new0);
    const POLYGON_INDEX_TYPE ipolyC =
      this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge_new2);

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices
      (cerr, "  New triangle: ", itriangleB, "\n");
    this->PrintPolygonIndexAndVertices
      (cerr, "  New polyA: ", ipolyA, "\n");
    polyA_tri_info.Print(cerr, "  ");
    polyA_tri_info.PrintEar(cerr, "  ear: ", "\n");
    this->PrintPolygonIndexAndVertices
      (cerr, "  New polyC: ", ipolyC, "\n");
    polyC_tri_info.Print(cerr, "  ");
    */

    SetNumSplitEdgePolygonVertices(ipolyA);    
    SetNumSplitEdgePolygonVertices(itriangleB);
    SetNumSplitEdgePolygonVertices(ipolyC);

    CopyTriInfoAndInteriorCoord
      (polyA_tri_info, interior_vcoord0, ipolyA);
    CopyTriInfoAndInteriorCoord
      (polyC_tri_info, interior_vcoord2, ipolyC);

    this->polygon_list[itriangleB].triangulation_info.Copy
      (triangleB_tri_info);
  }


  // Add new triangle vertex and replace triangle
  //   with triangle edge.
  // - Allow polygon triangulations using interior vertices.
  // - Set polygon triangulation information.
  // - Version using C++ STL vector for new_triangle_vcoord[].
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, 
            typename CTYPE4, typename CTYPE5,
            typename ITYPE2, typename ITYPE3, typename ITYPE4,
            typename RESULT_TYPE>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  AddVertexReplaceTriangleWithTriangleEdge
  (const std::vector<CTYPE2> & new_triangle_vcoord,
   const std::vector<CTYPE3> & interior_vcoord0,
   const std::vector<CTYPE4> & interior_vcoord2,
   std::vector<CTYPE5> & vertex_coord,
   const ITYPE2 ihalf_edge0, const RESULT_TYPE poly_tri_info[3],
   ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0)
  {
    AddVertexReplaceTriangleWithTriangleEdge
      (new_triangle_vcoord, interior_vcoord0, interior_vcoord2,
       vertex_coord, ihalf_edge0,
       poly_tri_info[0], poly_tri_info[1], poly_tri_info[2],
       iv_new, ihalf_edge_new0);
  }


  // Add new triangle vertex and replace triangle
  //   with triangle edge.
  // - Allow polygon triangulations using interior vertices.
  // - Set polygon triangulation information.
  // - Version using C++ STL vector for new_triangle_vcoord[]
  //     and for interior_vcoord0[] and for interior_vcoord2[].
  // - Version that does not return iv_new or ihalf_edge_new0.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, 
            typename CTYPE4, typename CTYPE5,
            typename ITYPE2, typename RESULT_TYPE>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  AddVertexReplaceTriangleWithTriangleEdge
  (const std::vector<CTYPE2> & new_triangle_vcoord,
   const std::vector<CTYPE3> & interior_vcoord0,
   const std::vector<CTYPE4> & interior_vcoord2,
   std::vector<CTYPE5> & vertex_coord,
   const ITYPE2 ihalf_edge0, const RESULT_TYPE poly_tri_info[3])
  {
    VERTEX_INDEX_TYPE iv_new;
    HALF_EDGE_INDEX_TYPE ihalf_edge_new0;

    AddVertexReplaceTriangleWithTriangleEdge
      (new_triangle_vcoord, interior_vcoord0, interior_vcoord2,
       vertex_coord, ihalf_edge0, poly_tri_info, 
       iv_new, ihalf_edge_new0);
  }

  // Add new triangle vertex and replace split edge triangle
  //   with triangle edge.
  // - Allow polygon triangulations using interior vertices.
  // - Set polygon triangulation information.
  // - Version using C++ STL vector for new_triangle_vcoord[]
  //     and for interior_vcoord0[] and for interior_vcoord2[].
  // - Version that does not return ihalf_edge_new0.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, 
            typename CTYPE4, typename CTYPE5,
            typename ITYPE2, typename ITYPE3, typename ITYPE4,
            typename RESULT_TYPE>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  AddVertexReplaceSplit2EdgeTriangle
  (const std::vector<CTYPE2> & new_triangle_vcoord,
   const std::vector<CTYPE3> & interior_vcoord0,
   const std::vector<CTYPE4> & interior_vcoord2,
   std::vector<CTYPE5> & vertex_coord,
   const ITYPE2 ihalf_edge0, const RESULT_TYPE poly_tri_info[3],
   ITYPE3 & iv_new, ITYPE4 & ihalf_edge_new0)
  {
    HALF_EDGE_INDEX_TYPE ihalf_edge1, ihalf_edge2, ihalf_edge3;

    this->GetIndicesOfNextThreeHalfEdgesInPolygon
      (ihalf_edge0, ihalf_edge1, ihalf_edge2, ihalf_edge3);

    const HALF_EDGE_INDEX_TYPE ihalf_edge1X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge3X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge3);

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "*** Replacing triangle ", ihalf_edge0, "\n");
    this->PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "  Adjacent polygon ", ihalf_edge1X, "\n");
    this->PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "  Adjacent polygon ", ihalf_edge3X, "\n");
    */

    MESH2D_SPLIT_II_BASE_TYPE::
      AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
      (Dimension(), new_triangle_vcoord, vertex_coord, 
       ihalf_edge0, iv_new, ihalf_edge_new0);
    SetVertexType(iv_new, VTYPE::TRIV_SPLITE);

    HALF_EDGE_INDEX_TYPE 
      ihalf_edge_new1, ihalf_edge_new2;

    this->GetIndicesOfNextTwoHalfEdgesInPolygon
      (ihalf_edge_new0, ihalf_edge_new1, ihalf_edge_new2);

    const POLYGON_INDEX_TYPE ipolyA =
      this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge_new1);
    const POLYGON_INDEX_TYPE itriangleB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge_new0);
    const POLYGON_INDEX_TYPE ipolyC =
      this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge_new2);

    SetNumSplitEdgePolygonVertices(ipolyA);
    SetNumSplitEdgePolygonVertices(itriangleB);
    SetNumSplitEdgePolygonVertices(ipolyC);

    // Note: Polygons ipolyA and ipolyC have not changed.
    CopyTriInfoAndInteriorCoord(poly_tri_info[0], interior_vcoord0, ipolyA);
    this->polygon_list[itriangleB].triangulation_info.Copy
      (poly_tri_info[1]);
    CopyTriInfoAndInteriorCoord(poly_tri_info[2], interior_vcoord2, ipolyC);
  }


  // Add new triangle vertex and replace split edge triangle
  //   with triangle edge.
  // - Allow polygon triangulations using interior vertices.
  // - Set polygon triangulation information.
  // - Version using C++ STL vector for new_triangle_vcoord[]
  //     and for interior_vcoord0[] and for interior_vcoord2[].
  // - Version that does not return iv_new or ihalf_edge_new0.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, 
            typename CTYPE4, typename CTYPE5,
            typename ITYPE2, typename RESULT_TYPE>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  AddVertexReplaceSplit2EdgeTriangle
  (const std::vector<CTYPE2> & new_triangle_vcoord,
   const std::vector<CTYPE3> & interior_vcoord0,
   const std::vector<CTYPE4> & interior_vcoord2,
   std::vector<CTYPE5> & vertex_coord,
   const ITYPE2 ihalf_edge0, const RESULT_TYPE poly_tri_info[3])
  {
    VERTEX_INDEX_TYPE iv_new;
    HALF_EDGE_INDEX_TYPE ihalf_edge_new0;

    AddVertexReplaceSplit2EdgeTriangle
      (new_triangle_vcoord, interior_vcoord0, interior_vcoord2,
       vertex_coord, ihalf_edge0, poly_tri_info, 
       iv_new, ihalf_edge_new0);
  }


  // *****************************************************************
  // Class MESH2D_TRIANGULATION_BASE member functions
  // *****************************************************************

  // Return true if two adjacent polygons are different.
  // - Return true if ihalf_edgeA or ihalf_edgeB are boundary edges.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPEH>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  AreAdjacentPolygonsDifferent
  (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB) const
  {
    if (this->IsBoundaryEdge(ihalf_edgeA)) { return(true); }
    if (this->IsBoundaryEdge(ihalf_edgeB)) { return(true); }

    if (this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeA) ==
        this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeB)) 
      { return(false); }

    return(true);
  }


  // Compute cosine of the min angle of triangle ipoly.
  // - Store the min angle in triangle ipoly.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeCosMinTriangleAngle
  (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
   const MTYPE2 max_small_magnitude)
  {
    IJK::compute_cos_min_triangle_angle
      (Dimension(), vcoord, this->PolygonVertexList(ipoly),
       max_small_magnitude,
       this->polygon_list[ipoly].triangulation_info.cos_min_triangulation_angle,
       this->polygon_list[ipoly].triangulation_info.flag_zero);
  }


  // Compute cosine of the min angle of polygon ipoly.
  // - Does NOT store information in data structure.
  // @param[out] cos_min_angle Cosine of the min polygon angle.
  // @param[out] iloc Location in polygon of vertex with min angle.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename ITYPE2, typename ITYPE3,
            typename MTYPE2, typename COS_TYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeCosMinPolygonAngle
  (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
   const MTYPE2 max_small_magnitude,
   COS_TYPE2 & cos_min_angle, ITYPE3 & iloc_min) const
  {
    const NTYPE num_poly_vert;
    const VERTEX_INDEX_TYPE * poly_vert = 
      this->PolygonVertexList(ipoly);
    COS_TYPE2 cos_max_angle;
    ITYPE iloc_max;
    NTYPE num_angle;

    compute_cos_min_max_polygon_angles
      (Dimension(), vcoord, poly_vert, num_poly_vert,
       max_small_magnitude, cos_min_angle, cos_max_angle,
       iloc_min, iloc_max, num_angle);
  }


  // Compute midpoint of line segment between two mesh vertices.
  // @param[out] midpoint_coord[] Midpoint coordinates.
  // @pre Array midpoint_coord[] is preallocated to length 
  //      at least Dimension().
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, 
            typename ITYPE2, typename ITYPE3>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeLineSegmentMidpoint
  (const std::vector<CTYPE2> & vcoord, 
   const ITYPE2 iv0, const ITYPE3 iv1,
   std::vector<CTYPE3> & midpoint_coord) const
  {
    compute_mesh2D_line_segment_midpoint
      (this->Dimension(), vcoord, *this, iv0, iv1, midpoint_coord);
  }


  // Compute edge midpoint.
  // @param[out] midpoint_coord[] Centroid coordinates.
  // @pre Array midpoint_coord[] is preallocated to length 
  //      at least Dimension().
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, typename ITYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeEdgeMidpoint
  (const std::vector<CTYPE2> & vcoord, const ITYPE2 ihalf_edge, 
   std::vector<CTYPE3> & midpoint_coord) const
  {
    compute_mesh2D_edge_midpoint
      (this->Dimension(), vcoord, *this, ihalf_edge, midpoint_coord);
  }


  // Compute midpoint of edge midpoints.
  // @param[out] midpoint_coord[] Midpoint coordinates.
  // @pre Array midpoint_coord[] is preallocated to length 
  //      at least Dimension().
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, 
            typename ITYPEA, typename ITYPEB>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeMidpointOfEdgeMidpoints
  (const std::vector<CTYPE2> & vcoord, 
   const ITYPEA ihalf_edgeA, const ITYPEB ihalf_edgeB,
   std::vector<CTYPE3> & midpoint_coord) const
  {
    compute_mesh2D_midpoint_of_edge_midpoints
      (this->Dimension(), vcoord, *this, ihalf_edgeA, ihalf_edgeB,
       midpoint_coord);
  }

  // Compute polygon centroid.
  // @param[out] centroid_coord Centroid coordinates.
  // @pre Array centroid_coord[] is preallocated to length 
  //      at least Dimension().
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, typename ITYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputePolygonCentroid
  (const CTYPE2 vcoord[], const ITYPE2 ipoly,
   CTYPE3 centroid_coord[]) const
  {
    compute_mesh2D_polygon_centroid
      (this->Dimension(), vcoord, *this, ipoly, centroid_coord);
  }


  // Compute polygon centroid.
  //  - Version using C++ STL vector centroid_coord.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename CTYPE3, typename ITYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputePolygonCentroid
  (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly,
   std::vector<CTYPE3> & centroid_coord) const
  {
    compute_mesh2D_polygon_centroid
      (this->Dimension(), vcoord, *this, ipoly, centroid_coord);
  }


  // Compute the triangulation of polygon ipoly, that maximizes
  //   the minimum angle.
  // - Store the triangulation information in polygon ipoly.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeCosMaxMinAngle
  (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
   const MTYPE2 max_small_magnitude)
  {
    compute_fan_triangulation_to_max_min_angle_M
      (Dimension(), vcoord, *this, ipoly, max_small_magnitude,
       this->polygon_list[ipoly].triangulation_info);
  }


  // Compute the triangulation of polygon ipoly, that maximizes
  //   the minimum angle.
  // - Store the triangulation information in polygon ipoly.
  // - Version using C++ STL vector vcoord[].
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeCosMaxMinAngle
  (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly, 
   const MTYPE2 max_small_magnitude)
  {
    ComputeCosMaxMinAngle
      (IJK::vector2pointer(vcoord), ipoly, max_small_magnitude);
  }


  // Compute the triangulation of polygon ipoly, that maximizes
  //   the minimum angle.
  // - Store the triangulation information in polygon ipoly.
  // - Version that uses POLYGON_TRIANGULATION_SETTINGS.
  //   to specify allowable trianguations.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeCosMaxMinAngle
  (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
   const MTYPE2 max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings)
  {
    ComputePolygonCentroid
      (vcoord, ipoly, this->polygon_list[ipoly].interior_coord);

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices(cerr, "Polygon: ", ipoly, "\n");
    */

    compute_triangulation_to_max_min_angle
      (Dimension(), vcoord, PolygonInteriorCoord(ipoly),
       this->NumPolygonVertices(ipoly),
       this->PolygonVertexList(ipoly), 
       max_small_magnitude, triangulation_settings,
       this->polygon_list[ipoly].triangulation_info);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "Triangulation info:" << endl;
    this->polygon_list[ipoly].triangulation_info.Print(cerr, "  ");
    */
  }


  // Compute the triangulation of polygon ipoly, that maximizes
  //   the minimum angle.
  // - Store the triangulation information in polygon ipoly.
  // - Version that uses POLYGON_TRIANGULATION_SETTINGS.
  //   to specify allowable trianguations.
  // - Version using C++ STL vector vcoord[].
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeCosMaxMinAngle
  (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly, 
   const MTYPE2 max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS & triangulation_settings)
  {
    ComputeCosMaxMinAngle
      (IJK::vector2pointer(vcoord), ipoly, max_small_magnitude,
       triangulation_settings);
  }

  // Compute the triangulation of polygon ipoly, that maximizes
  //   the minimum angle.
  // - Consider triangulations from a single polygon vertex
  //     or from the polygon centroid.
  // - Store the triangulation information in polygon ipoly.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename ITYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeCosMaxMinAngleAllowCentroid
  (const CTYPE2 vcoord[], const ITYPE2 ipoly, 
   const MTYPE2 max_small_magnitude)
  {
    ComputePolygonCentroid
      (vcoord, ipoly, this->polygon_list[ipoly].interior_coord);

    compute_fan_or_interior_vertex_triangulation_to_max_min_angle_M
      (Dimension(), vcoord, PolygonInteriorCoord(ipoly),
       *this, ipoly, max_small_magnitude,
       this->polygon_list[ipoly].triangulation_info);
  }


  // Compute triangulations of each polygon, maximizing
  //   the minimum angle in each triangulation.
  // - Store the triangulation information in polygon ipoly.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeAllPolygonsCosMaxMinAngle
  (const CTYPE2 vcoord[], const MTYPE2 max_small_magnitude)
  {
    bool flag_zero;

    for (POLYGON_INDEX_TYPE ipoly = 0; 
         ipoly < this->NumPolygons(); ipoly++) {

      if (this->IsPolygonDeleted(ipoly)) { continue; }

      if (this->IsTriangle(ipoly)) {
        // Compute min triangle angle.
        ComputeCosMinTriangleAngle(vcoord, ipoly, max_small_magnitude);
      }
      else {
        ComputeCosMaxMinAngle(vcoord, ipoly, max_small_magnitude);
      }
    }
  }


  // Compute triangulations of each polygon, maximizing
  //   the minimum angle in each triangulation.
  //  - Version using C++ STL vector vertex_coord[].
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeAllPolygonsCosMaxMinAngle
  (const std::vector<CTYPE2> & vcoord, const MTYPE2 max_small_magnitude)
  {
    ComputeAllPolygonsCosMaxMinAngle
      (IJK::vector2pointer(vcoord), max_small_magnitude);
  }


  // Compute triangulations of each polygon, maximizing
  //   the minimum angle in each triangulation.
  // - Consider triangulations from a single polygon vertex
  //     or from the polygon centroid.
  // - Store the triangulation information in polygon ipoly.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename MTYPE2>
  void MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  ComputeAllPolygonsCosMaxMinAngleAllowCentroid
  (const CTYPE2 vcoord[], const MTYPE2 max_small_magnitude)
  {
    bool flag_zero;

    for (POLYGON_INDEX_TYPE ipoly = 0; 
         ipoly < this->NumPolygons(); ipoly++) {

      if (this->IsPolygonDeleted(ipoly)) { continue; }

      if (this->IsTriangle(ipoly)) {
        // Compute min triangle angle.
        ComputeCosMinTriangleAngle
          (vcoord, ipoly, max_small_magnitude);
      }
      else {
        ComputeCosMaxMinAngleAllowCentroid
          (vcoord, ipoly, max_small_magnitude);
      }
    }
  }


  // Compute triangulations of each polygon, maximizing
  //   the minimum angle in each triangulation.
  // - Consider triangulations from a single polygon vertex
  //     or from the polygon centroid.
  //  - Version using C++ STL vector vertex_coord[].
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename MTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  ComputeAllPolygonsCosMaxMinAngleAllowCentroid
  (const std::vector<CTYPE2> & vcoord, const MTYPE2 max_small_magnitude)
  {
    ComputeAllPolygonsCosMaxMinAngleAllowCentroid
      (IJK::vector2pointer(vcoord), max_small_magnitude);
  }


  // Triangulate polygon ipoly based on triangulation encoding.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, int BIT_SIZE_TYPE_X>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::TriangulatePolygon
  (const ITYPE2 ipoly, 
   const POLYGON_TRIANGULATION_ENCODING<BIT_SIZE_TYPE_X> & tri_encoding)
  {
    const NTYPE NUM_VERTICES_PER_TRIANGLE(3);
    const NTYPE numv = this->NumPolygonVertices(ipoly);
    VERTEX_INDEX_TYPE triangle_vert[NUM_VERTICES_PER_TRIANGLE];
    std::vector<HALF_EDGE_INDEX_TYPE> half_edge_stack(numv);
    std::vector<VERTEX_INDEX_TYPE> tri_list;
    IJK::PROCEDURE_ERROR error("MESH2D_TRIANGULATION_BASE::TriangulatePolygon");

    gen_polygon_triangulation_list_based_on_encoding
      (numv, tri_encoding, tri_list);

    if (tri_list.size() == 0) {
      // No triangles
      return;
    }

    const NTYPE num_new_tri = tri_list.size()/3;

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices(cerr, "Polygon: ", ipoly, "\n");
    cerr << " tri_list:";
    for (int i = 0; i < num_new_tri; i++) {
      cerr << "  (" << tri_list[3*i] << "," << tri_list[3*i+1] << "," << tri_list[3*i+2] << ")";
    }
    cerr << endl;
    */

    triangle_vert[0] = this->IndexOfVertexInPolygon(ipoly, tri_list[0]);
    triangle_vert[1] = this->IndexOfVertexInPolygon(ipoly, tri_list[1]);
    triangle_vert[2] = this->IndexOfVertexInPolygon(ipoly, tri_list[2]);

    const POLYGON_INDEX_TYPE itriangle0 = this->_AddTriangle(triangle_vert);

    for (NTYPE i = 1; i < num_new_tri; i++) {
      triangle_vert[0] = this->IndexOfVertexInPolygon(ipoly, tri_list[3*i]);
      triangle_vert[1] = this->IndexOfVertexInPolygon(ipoly, tri_list[3*i+1]);
      triangle_vert[2] = this->IndexOfVertexInPolygon(ipoly, tri_list[3*i+2]);

      this->_AddTriangle(triangle_vert);
    }

    NTYPE nstack = 0;

    // Replace/link half edges.
    for (NTYPE i = 0; i < num_new_tri; i++) {
      const VERTEX_INDEX_TYPE jloc[NUM_VERTICES_PER_TRIANGLE] = 
        { tri_list[3*i], tri_list[3*i+1], tri_list[3*i+2] };

      const POLYGON_INDEX_TYPE itriangle = itriangle0+i;

      for (NTYPE k0 = 0; k0 < NUM_VERTICES_PER_TRIANGLE; k0++) {

        // *** DEBUG ***
        /*
        using namespace std;
        cerr << "  nstack: " << nstack
             << "  itriangle: " << itriangle 
             << "  k0: " << k0 << endl;
        if (nstack > 0) {
          this->PrintHalfEdgeIndexAndEndpoints
            (cerr, "  Top: ", half_edge_stack[nstack-1], "\n");
        }
        */

        const NTYPE k1 = (k0+1)%NUM_VERTICES_PER_TRIANGLE;
        if (((jloc[k0]+1)%numv) == jloc[k1]) {
          this->_ReplaceHalfEdge
            (this->HalfEdgeIndex(ipoly,jloc[k0]),
             this->HalfEdgeIndex(itriangle, k0));
        }
        else if (nstack == 0) {
          half_edge_stack[nstack] = this->HalfEdgeIndex(itriangle, k0);
          nstack++;
        }
        else {
          const HALF_EDGE_INDEX_TYPE ihalf_edgeA = half_edge_stack[nstack-1];
          const HALF_EDGE_INDEX_TYPE ihalf_edgeB = this->HalfEdgeIndex(itriangle, k0);

          if ((this->FromVertexIndex(ihalf_edgeA) ==
               this->ToVertexIndex(ihalf_edgeB)) &&
              (this->ToVertexIndex(ihalf_edgeA) ==
               this->FromVertexIndex(ihalf_edgeB))) {

            this->_LinkTwoHalfEdgesAroundEdge
              (ihalf_edgeA, ihalf_edgeB);
            nstack--;
          }
          else {
            half_edge_stack[nstack] = ihalf_edgeB;
            nstack++;
          }
        }
      }
    }

    this->_DeletePolygon(ipoly);

    if (nstack != 0) {
      const HALF_EDGE_INDEX_TYPE ihalf_edgeA = half_edge_stack[nstack-1];

      error.AddMessage
        ("Programming error. Not all triangle edges are linked.");
      error.AddMessage
        ("  Error in triangulating polygon ", ipoly, ".");
      error.AddMessage
        ("  Internal half edge (", this->FromVertexIndex(ihalf_edgeA), ",",
         this->ToVertexIndex(ihalf_edgeA), ") is not linked to another half edge.");
      throw error;
    }

  }


  // Triangulate polygon ipoly.
  // - Use polygon triangulation information to determine triangulation.
  // @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
  //    has been called on ipoly to determine triangulation.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  TriangulatePolygon(const ITYPE2 ipoly, std::vector<CTYPE2> & vcoord)
  {
    NUMBER_TYPE j_ear;
    IJK::PROCEDURE_ERROR error("MESH2D_TRIANGULATION_BASE::TriangulatePolygon");

    // *** DEBUG
    /*
    using namespace std;
    cerr << "In " << __func__ << " Poly: " << ipoly << endl;
    */

    if (this->polygon_list[ipoly].triangulation_info.GetFirstEar(j_ear)) {

      const HALF_EDGE_INDEX_TYPE jhalf_edge1 = 
        this->HalfEdgeIndex(ipoly, j_ear);

      // *** DEBUG ***
      /*
      using namespace std;
      this->PrintPolygonIndexAndVertices(cerr, "  Polygon: ", ipoly, "\n");
      cerr << "    Split vertex at loc: " << j_ear
           << "  vertex: " << this->PolygonVertex(ipoly, j_ear) << endl;
      cerr << "    Triangulation vertex location: " << TriVertexIndex(ipoly)
           << "  Triangulation vertex: " << this->PolygonVertex(ipoly, TriVertexIndex(ipoly))
           << endl;
      */

      HALF_EDGE_INDEX_TYPE jhalf_edge1_new, jhalf_edge2_new;
      this->CutPolygonEar
        (jhalf_edge1, jhalf_edge1_new, jhalf_edge2_new);

      const POLYGON_INDEX_TYPE ipoly2 = 
        this->IndexOfPolygonContainingHalfEdge(jhalf_edge2_new);

      if (NumInteriorTriVertices(ipoly) == 0) {

        if (j_ear == TriVertexIndex(ipoly)) {
          error.AddMessage
            ("Programming error.  Triangulation vertex is in ear cut from polygon.");
          error.AddMessage
            ("  Ear vertex: ", this->PolygonVertex(ipoly, j_ear), ".");
          error.AddMessage
            ("  Triangulation vertex: ", 
             this->PolygonVertex(ipoly, TriVertexIndex(ipoly)), ".");
          throw error;
        }

        const NUMBER_TYPE numv = this->NumPolygonVertices(ipoly);
        const NUMBER_TYPE jloc2 = (j_ear+1)%numv;

        const ITYPE tri_vertex_index =
          this->LocationRelativeToLocation(ipoly, TriVertexIndex(ipoly), jloc2);

        // *** DEBUG ***
        /*
        using namespace std;
        this->PrintPolygonIndexAndVertices(cerr, "  New polygon: ", ipoly2, "\n");
        cerr << "    Triangulation vertex index: " << tri_vertex_index
             << "  Triangulation vertex: " << this->PolygonVertex(ipoly2, tri_vertex_index)
             << endl;
        */
        
        this->TriangulatePolygonFromVertex(ipoly2, tri_vertex_index);
      }
      else {
        const VERTEX_INDEX_TYPE ivX = IJK::insert_coord
          (Dimension(), PolygonInteriorCoord(ipoly), vcoord);
        this->AddVertex();
        this->TriangulatePolygonFromInteriorVertex(ipoly2, ivX);
      }

    }
    else {

      if (NumInteriorTriVertices(ipoly) == 0) {

        if (PolygonTriangulationInfo(ipoly).triangulation_encoding.IsZero()) {

          // *** DEBUG ***
	  /*
          using namespace std;
          cerr << "  Calling TriangulatePolygonFromVertex on poly "
               << ipoly
               << "  Vertex loc: " << TriVertexIndex(ipoly) << endl;
	  */

          this->TriangulatePolygonFromVertex(ipoly, TriVertexIndex(ipoly));
        }
        else {
          // *** DEBUG ***
	  /*
          using namespace std;
          this->PrintPolygonIndexAndVertices
	    (cerr, "  Triangulating polygon: ", ipoly, " using encoding ");
          PolygonTriangulationInfo(ipoly).triangulation_encoding.PrintTriangulationRepresentation
            (cerr, this->NumPolygonVertices(ipoly), "", "\n");  
	  */

          // Use triangulation_encoding.
          this->TriangulatePolygon
	    (ipoly, PolygonTriangulationInfo(ipoly).triangulation_encoding);
        }
      }
      else {
 
        // *** DEBUG ***
        /*
        using namespace std;
        cerr << "  Calling TriangulatePolygonFromInteriorVertex "
             << " poly " << ipoly << endl;
        */

        const VERTEX_INDEX_TYPE ivX = IJK::insert_coord
          (Dimension(), PolygonInteriorCoord(ipoly), vcoord);
        this->AddVertex();
        this->TriangulatePolygonFromInteriorVertex(ipoly, ivX);
      }

    }

  }


  // Triangulate all polygons.
  // - Use polygon triangulation information to determine triangulation.
  // @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
  //    has been called on all undeleted polygons with 4 or more vertices
  //    to determine triangulation.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  TriangulateAllPolygons(std::vector<CTYPE2> & vcoord)
  {
    for (POLYGON_INDEX_TYPE ipoly = 0; 
         ipoly < this->NumPolygons(); ipoly++) {

      if (this->IsTriangle(ipoly)) { continue; }
      if (this->IsPolygonDeleted(ipoly)) { continue; }

      TriangulatePolygon(ipoly, vcoord);
    }

  }


  // *****************************************************************
  // MESH2D_TRIANGULATION_BASE compute polygon characteristics
  // *****************************************************************

  // Return true if edge is "long".
  // - Return true if edge length squared is greater than or equal to L.
  template <int BIT_SET_SIZE,
            typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename CTYPE2, typename ITYPE2, typename LTYPE2>
  bool MESH2D_TRIANGULATION_BASE
  <BIT_SET_SIZE,VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  IsLongEdge(const std::vector<CTYPE2> & vcoord,
             const ITYPE2 ihalf_edge, const LTYPE2 L) const
  {
    return(is_long_distance
           (Dimension(), vcoord, this->FromVertexIndex(ihalf_edge),
            this->ToVertexIndex(ihalf_edge), L));
  }


  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename ITYPE2, 
            typename ITYPE3, typename RTYPE>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  IsLongQuad(const std::vector<CTYPE2> & vcoord,
             const ITYPE2 iquad, const RTYPE Rsquared,
             ITYPE3 & ilongest_half_edge) const
  {
    ITYPE3 jlongest_loc;  // Location of longest edge in quad.
    bool flag_long;

    if (!this->IsQuadrilateral(iquad)) { return(false); }

    flag_long = is_long_quad
      (Dimension(), vcoord, this->PolygonVertexList(iquad), 
       Rsquared, jlongest_loc);

    ilongest_half_edge = this->HalfEdgeIndex(iquad, jlongest_loc);

    return(flag_long);
  }


  // Return true if shortest edge is between two long edges.
  // @param[out] ishortest_half_edge Shortest half edge of quad.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename ITYPE2, 
            typename ITYPE3, typename RTYPE>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  DoesQuadHaveLongShortLongEdges
  (const std::vector<CTYPE2> & vcoord,
   const ITYPE2 iquad, const RTYPE Rsquared,
   ITYPE3 & ishortest_half_edge) const
  {
    ITYPE3 jshortest_loc;  // Location of shortest edge in quad.

    if (!this->IsQuadrilateral(iquad)) { return(false); }

    const bool flag_long_short_long = 
      is_shortest_quad_edge_between_long_edges
      (Dimension(), vcoord, this->PolygonVertexList(iquad), 
       Rsquared, jshortest_loc);

    ishortest_half_edge = this->HalfEdgeIndex(iquad, jshortest_loc);

    return(flag_long_short_long);
  }

  // Return true if quad has one edge much shorter 
  //   than all other edges.
  // @param[out] ishortest_half_edge Shortest half edge of polygon.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename ITYPE2, 
            typename ITYPE3, typename RTYPE>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::DoesQuadHaveVeryShortEdge
  (const std::vector<CTYPE2> & vcoord,
   const ITYPE2 iquad, const RTYPE Rsquared,
   ITYPE3 & ishortest_half_edge) const
  {
    ITYPE3 jshortest_loc;  // Location of shortest edge in polygon
    bool flag_short;

    if (!this->IsQuadrilateral(iquad)) { return(false); }

    flag_short = does_quad_have_very_short_edge
      (Dimension(), vcoord, this->PolygonVertexList(iquad), 
       Rsquared, jshortest_loc);

    ishortest_half_edge = this->HalfEdgeIndex(iquad, jshortest_loc);

    return(flag_short);
  }


  // Return index of shortest half edge in polygon.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename ITYPE2, typename LTYPE2>
  typename _MESH2D_TRIANGULATION_BASE_T_NAME_::HALF_EDGE_INDEX_TYPE 
  _MESH2D_TRIANGULATION_BASE_T_NAME_::ComputeShortestHalfEdge
  (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly,
   LTYPE2 & shortest_length_squared) const
  {
    return(compute_shortest_polygon_half_edge
           (Dimension(), vcoord, *this, ipoly, shortest_length_squared));
  }


  // Return index of shortest half edge in polygon.
  // - Version that does not return shortest_length_squared.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename ITYPE2>
  typename _MESH2D_TRIANGULATION_BASE_T_NAME_::HALF_EDGE_INDEX_TYPE 
  _MESH2D_TRIANGULATION_BASE_T_NAME_::ComputeShortestHalfEdge
  (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly) const
  {
    return(compute_shortest_polygon_half_edge
           (Dimension(), vcoord, *this, ipoly));
  }


  // Get index of shortest non-split half edge in polygon.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename CTYPE2, typename ITYPE2, 
            typename ITYPE3, typename LTYPE2>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::ComputeShortestNonSplitHalfEdge
  (const std::vector<CTYPE2> & vcoord, const ITYPE2 ipoly,
   ITYPE3 & ishortest, LTYPE2 & shortest_length_squared) const
  {
    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices
      (cerr, "  Computing shortest non-split half edge for polygon ",
       ipoly, "\n");
    */

    return(compute_shortest_non_split_polygon_half_edge
           (Dimension(), vcoord, *this, ipoly, ishortest,
            shortest_length_squared));
  }


  // Return true if polygon is a triangle with one split edge.
  // @param[out] If true, then ihalf_edge0 is non-split triangle edge.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  IsTriangleWithOneSplitEdge(const ITYPE2 ipoly) const
  {
    const NTYPE ONE(1);

    if (this->NumSplitEdgeVertices(ipoly) != ONE)
      { return(false); }
    if (!this->IsQuadrilateral(ipoly)) { return(false); }

    // Polygon has 4 vertices, but two adjacent vertices are splitting vertices.
    for (ITYPE i = 0; i < this->NumPolygonEdges(ipoly); i++) {
      
      const HALF_EDGE_INDEX_TYPE ihalf_edge = this->HalfEdgeIndex(ipoly, i);
      const VERTEX_INDEX_TYPE iv0 = this->FromVertexIndex(ihalf_edge);
      const VERTEX_INDEX_TYPE iv1 = this->ToVertexIndex(ihalf_edge);

      if (this->TriangulationVertexType(iv0) != VTYPE::TRIV_INPUT &&
          this->TriangulationVertexType(iv1) != VTYPE::TRIV_INPUT) {
        // A triangle edge is split twice.
        return(false);
      }
    }

    return(true);
  }


  // Return true if polygon is a triangle with two split edges.
  // @param[out] If true, then ihalf_edge0 is non-split triangle edge.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename ITYPE3>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  IsTriangleWithTwoSplitEdges
  (const ITYPE2 ipoly, ITYPE3 & ihalf_edge0) const
  {
    const NTYPE TWO(2);

    // Initialize
    ihalf_edge0 = this->HalfEdgeIndex(ipoly, 0);

    if (this->NumSplitEdgeVertices(ipoly) != TWO) 
      { return(false); }
    if (!this->IsPentagon(ipoly)) { return(false); }

    // Polygon has 5 vertices, but two adjacent vertices are splitting vertices.
    for (ITYPE i = 0; i < this->NumPolygonEdges(ipoly); i++) {
      
      const HALF_EDGE_INDEX_TYPE ihalf_edge = this->HalfEdgeIndex(ipoly, i);
      const VERTEX_INDEX_TYPE iv0 = this->FromVertexIndex(ihalf_edge);
      const VERTEX_INDEX_TYPE iv1 = this->ToVertexIndex(ihalf_edge);

      if (this->TriangulationVertexType(iv0) != VTYPE::TRIV_INPUT &&
          this->TriangulationVertexType(iv1) != VTYPE::TRIV_INPUT) {
        // A triangle edge is split twice.
        return(false);
      }

      if (this->TriangulationVertexType(iv0) == VTYPE::TRIV_INPUT &&
          this->TriangulationVertexType(iv1) == VTYPE::TRIV_INPUT) {
        ihalf_edge0 = ihalf_edge;
      }
    }

    return(true);
  }


  // Return number of vertices of a given type in a polygon.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2>
  typename _MESH2D_TRIANGULATION_BASE_T_NAME_::NUMBER_TYPE
  _MESH2D_TRIANGULATION_BASE_T_NAME_::CountNumPolygonVerticesOfType
  (const ITYPE2 ipoly, 
   const typename 
   _MESH2D_TRIANGULATION_BASE_T_NAME_::TRIANGULATION_VERTEX_TYPE vtype) const
  {
    NUMBER_TYPE numv_of_type = 0;

    for (NTYPE iloc = 0; iloc < this->NumPolygonVertices(ipoly); iloc++) {
      const VERTEX_INDEX_TYPE iv = 
        this->IndexOfVertexInPolygon(ipoly, iloc);
      if (TriangulationVertexType(iv) == vtype) 
        { numv_of_type++; }
    }

    return(numv_of_type);
  }


  // Count and set number of split edge vertices in a polygon.
  // - Stores number in data structure.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  SetNumSplitEdgePolygonVertices(const ITYPE2 ipoly)
  {
    this->polygon_list[ipoly].num_split_edge_vertices =
      CountNumSplitEdgePolygonVertices(ipoly);
  }


  // *****************************************************************
  // MESH2D_TRIANGULATION_BASE Set triangulation information.
  // *****************************************************************


  // Set triangulation to include one interior vertex.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename COS_TYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  SetInteriorTriVertex
  (const ITYPE2 ipoly, const COS_TYPE2 cos_min, 
   const CTYPE2 interior_coord[])
  { 
    this->polygon_list[ipoly].triangulation_info.SetInterior
      (cos_min, ipoly, 1); 
    std::copy(interior_coord, interior_coord+Dimension(),
              this->polygon_list[ipoly].interior_coord);
  }


  // Set triangulation to include one interior vertex.
  // - Version using C++ STL vector interior_coord[].
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename COS_TYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  SetInteriorTriVertex
  (const ITYPE2 ipoly, const COS_TYPE2 cos_min, 
   const std::vector<CTYPE2> & interior_coord)
  { 
    SetInteriorTriVertex(ipoly, cos_min, IJK::vector2pointer(interior_coord));
  }


  // Copy triangulation information from poly_tri_info and interior_coord
  //   into triangulation info for polygon ipoly.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename INFO_TYPE, typename CTYPE2, typename ITYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  CopyTriInfoAndInteriorCoord
  (const INFO_TYPE & poly_tri_info, const CTYPE2 interior_coord[],
   const ITYPE2 ipoly)
  {
    CopyTriInfo(poly_tri_info, ipoly);

    if (poly_tri_info.num_interior_tri_vertices > 0) {
      std::copy(interior_coord, interior_coord+Dimension(),
                this->polygon_list[ipoly].interior_coord);
    }
  }

  // Copy triangulation information from poly_tri_info and interior_coord
  //   into triangulation info for polygon ipoly.
  // - Version using C++ STL vector interior_coord[].
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename INFO_TYPE, typename CTYPE2, typename ITYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  CopyTriInfoAndInteriorCoord
  (const INFO_TYPE & poly_tri_info, 
   const std::vector<CTYPE2> & interior_coord, const ITYPE2 ipoly)
  {
    CopyTriInfoAndInteriorCoord
      (poly_tri_info, IJK::vector2pointer(interior_coord), ipoly);
  }


  // Copy triangulation information from poly_tri_info
  //   into triangulation info for polygon ipoly.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename INFO_TYPE, typename ITYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  CopyTriInfo(const INFO_TYPE & poly_tri_info, const ITYPE2 ipoly)
  {
    this->polygon_list[ipoly].triangulation_info.Copy(poly_tri_info);
  }


  // Set triangulation information of polygon with split edge.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename ITYPE3, typename INFO_TYPE>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  SetSplitEdgePolygonTriInfo
  (const ITYPE2 ihalf_edge0, const INFO_TYPE & poly0_tri_info,
   const ITYPE3 ipoly1)
  {
    const POLYGON_INDEX_TYPE ipoly0 =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    const NTYPE numv0 = this->NumPolygonVertices(ipoly0);
    IJK::PROCEDURE_ERROR
      error("MESH2D_TRIANGULATION_BASE::SetSplitEdgePolygonTriInfo");

    CopyTriInfo(poly0_tri_info, ipoly1);

    if (poly0_tri_info.HasEar()) {

      // iloc0 is location in ipoly0 of first edge in ipoly1
      const ITYPE iloc0 =
        this->LocationOfNextHalfEdgeInPolygon(ihalf_edge0);

      const POLYGON_INDEX_TYPE ipoly0 =
        this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);

      this->polygon_list[ipoly1].triangulation_info.RotateRightEarFlags
        (iloc0, this->NumPolygonVertices(ipoly0));
    }


    if (poly0_tri_info.num_interior_tri_vertices == 0) {

      // Compute tri_vertex_index.

      if (poly0_tri_info.tri_vertex_index < numv0) {

        // iloc0 is location in ipoly0 of first edge in ipoly1
        const ITYPE iloc0 =
          this->LocationOfNextHalfEdgeInPolygon(ihalf_edge0);

        // Compute tri_vertex_index as location of tri_vertex_index
        //   relative to iloc0.
        const ITYPE tri_vertex_index = 
          this->LocationRelativeToLocation
          (ipoly0, poly0_tri_info.tri_vertex_index, iloc0);

        // *** DEBUG ***
        /*
        using namespace std;
        this->PrintHalfEdge(cerr, "  ihalf_edge0: ", ihalf_edge0, "\n");
        cerr << "  ihalf_edge0: " << ihalf_edge0
             << " iloc0: " << iloc0 << endl;
        cerr << "    poly0_tri_info.tri_vertex_index: "
             << poly0_tri_info.tri_vertex_index
             << " new tri_vertex_index: " << tri_vertex_index << endl;
        */

        SetTriVertexIndex(ipoly1, tri_vertex_index);
      }
      else {
        // Set tri_vertex_index to location of split vertex.
        // Split vertex has location numv0.
        SetTriVertexIndex(ipoly1, numv0);
      }
    }
    else {
      error.AddMessage
        ("Programming error. Wrong number of interior triangulation vertices.");
      error.AddMessage
        ("  Number of interior triangulation vertices: ",
         poly0_tri_info.NumInteriorTriVertices(), ".");
      error.AddMessage
        ("  Expected 0 interior triangulation vertices.");
      throw error;
    }

  }


  // Set triangulation information of polygon with split edge.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename ITYPE3, typename INFO_TYPE,
            typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  SetSplitEdgePolygonTriInfoAllowIV
  (const ITYPE2 ihalf_edge0, const INFO_TYPE & poly0_tri_info,
   const CTYPE2 interior_coord[], const ITYPE3 ipoly1)
  {
    if (poly0_tri_info.num_interior_tri_vertices == 0) {
      SetSplitEdgePolygonTriInfo(ihalf_edge0, poly0_tri_info, ipoly1);
    }
    else {
      INFO_TYPE poly1_tri_info = poly0_tri_info;

      if (poly0_tri_info.HasEar()) {

        // iloc0 is location in ipoly0 of first edge in ipoly1
        const ITYPE iloc0 =
          this->LocationOfNextHalfEdgeInPolygon(ihalf_edge0);

        const POLYGON_INDEX_TYPE ipoly0 =
          this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);

        poly1_tri_info.RotateRightEarFlags
          (iloc0, this->NumPolygonVertices(ipoly0));
      }

      CopyTriInfoAndInteriorCoord
        (poly1_tri_info, interior_coord, ipoly1);
    }

  }


  // Set triangulation information of polygon with split edge.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename ITYPE3, typename INFO_TYPE,
            typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  SetSplitEdgePolygonTriInfoAllowIV
  (const ITYPE2 ihalf_edge0, const INFO_TYPE & poly0_tri_info,
   const std::vector<CTYPE2> & interior_coord, const ITYPE3 ipoly1)
  {
    SetSplitEdgePolygonTriInfoAllowIV
      (ihalf_edge0, poly0_tri_info, IJK::vector2pointer(interior_coord),
       ipoly1);
  }


  // *****************************************************************
  // MESH2D_TRIANGULATION_BASE delete edge member functions
  // *****************************************************************

  // Delete edge to improve max min triangulation.
  // - Return true if edge is deleted.
  // @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
  //    has been called on ipoly to determine triangulation.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename MTYPE2, typename CTYPE2>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  DeleteEdgeToMaxMinTriangulationAngle
  (const ITYPE2 ihalf_edgeA, const MTYPE2 max_small_magnitude,
   std::vector<CTYPE2> & vcoord)
  {
    const POLYGON_INDEX_TYPE ipolyA =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const VERTEX_INDEX_TYPE iv0 = this->FromVertexIndex(ihalf_edgeA);
    const VERTEX_INDEX_TYPE iv1 = this->ToVertexIndex(ihalf_edgeA);
    const CTYPE2 * vcoord_ptr = IJK::vector2pointer(vcoord);
    const CTYPE2 * v0coord = vcoord_ptr + iv0*Dimension();
    const CTYPE2 * v1coord = vcoord_ptr + iv1*Dimension();
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);
    std::vector<CTYPE2> centroidA_coord(Dimension());
    std::vector<CTYPE2> centroidB_coord(Dimension());
    std::vector<CTYPE2> centroidC_coord(Dimension());
    COS_TYPE cos_min_angle, cosA, cosB, cosC;
    bool flag_zero, flagA_zero, flagB_zero, flagC_zero;
    IJK::PROCEDURE_ERROR error
      ("DeleteEdgeToMaxMinTriangulationAngle");
    
    if (this->IsBoundaryEdge(ihalf_edgeA)) {
      // Never delete a boundary edge.
      return(false);
    }

    if (this->DoPolygonsShareThreeVertices(ihalf_edgeA)) {
      // Don't delete ihalf_edgeA if some third vertex is shared
      //   by the two polygons incident on ihalf_edgeA.
      return(false);
    }

    const POLYGON_INDEX_TYPE ipolyB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeB);

    compute_mesh2D_polygon_centroid
      (Dimension(), vcoord, *this, ipolyA, centroidA_coord);
    compute_mesh2D_polygon_centroid
      (Dimension(), vcoord, *this, ipolyB, centroidB_coord);
    compute_mesh2D_two_polygon_centroid
      (Dimension(), vcoord, *this, ihalf_edgeA, centroidC_coord);

    compute_cos_min_merge_two_polygons_triangulation_angle
      (Dimension(), vcoord, centroidA_coord, *this,
       ihalf_edgeA, max_small_magnitude, cosA, flagA_zero);
    compute_cos_min_merge_two_polygons_triangulation_angle
      (Dimension(), vcoord, centroidB_coord, *this,
       ihalf_edgeA, max_small_magnitude, cosB, flagB_zero);
    compute_cos_min_merge_two_polygons_triangulation_angle
      (Dimension(), vcoord, centroidC_coord, *this,
       ihalf_edgeA, max_small_magnitude, cosC, flagC_zero);

    ITYPE index_selected;
    IJK::select_minIII
      (cosA, flagA_zero,
       cosB, flagB_zero,
       cosC, flagC_zero,
       cos_min_angle, index_selected, flag_zero);

    if (!flag_zero) {
      if (cos_min_angle < PolygonCosMinTriangulationAngle(ipolyA)) {
        POLYGON_INDEX_TYPE ipoly_new;
        this->DeleteEdge(ihalf_edgeA, ipoly_new);

        const NTYPE numv = this->NumPolygonVertices(ipoly_new);

        this->polygon_list[ipoly_new].triangulation_info.SetInterior
          (cos_min_angle, numv, 1);

        switch(index_selected) {
          
        case 0:
          IJK::copy_coord
            (Dimension(), IJK::vector2pointer(centroidA_coord), 
             this->polygon_list[ipoly_new].interior_coord);
          break;

        case 1:
          IJK::copy_coord
            (Dimension(), IJK::vector2pointer(centroidB_coord), 
             this->polygon_list[ipoly_new].interior_coord);
          break;

        case 2:
        default:
          IJK::copy_coord
            (Dimension(), IJK::vector2pointer(centroidC_coord), 
             this->polygon_list[ipoly_new].interior_coord);
          break;
        }

        return(true);
      }
    }


    return(false);
  }


  /// Delete longest polygon edge to improve max min triangulation.
  /// - Return true if edge is split.
  /// @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
  ///    has been called on ipoly to determine triangulation.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename MTYPE2, typename CTYPE2>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  DeleteLongestEdgeToMaxMinTriangulationAngle
  (const ITYPE2 ipoly, const MTYPE2 max_small_magnitude, 
   std::vector<CTYPE2> & vcoord)
  {
    const HALF_EDGE_INDEX_TYPE ilongest_half_edge = 
      compute_longest_polygon_half_edge
      (Dimension(), vcoord, *this, ipoly);

    return(DeleteEdgeToMaxMinTriangulationAngle
           (ilongest_half_edge, max_small_magnitude, vcoord));
  }


  // Delete (long) triangle edges to improve max min triangulation.
  // @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
  //    has been called on all polygons to determine triangulation.
  // @param cos_angle_threshold Attempt to delete longest triangle edge
  //    if cos_min_triangulation_angle is above cos_angle_threshold.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  DeleteTriangleEdgesToMaxMinTriangulationAngle 
  (const COS_TYPE2 cos_angle_threshold, const MTYPE2 max_small_magnitude,
   std::vector<CTYPE2> & vcoord)
  {
    const NUMBER_TYPE num_poly = this->NumPolygons();

    for (ITYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (this->IsPolygonDeleted(ipoly)) { continue; }

      if (!this->IsTriangle(ipoly)) { continue; }

      if (PolygonCosMinTriangulationAngle(ipoly) > cos_angle_threshold) {

        DeleteLongestEdgeToMaxMinTriangulationAngle
          (ipoly, max_small_magnitude, vcoord);
      }
    }
  }


  // Delete (long) triangle edges to improve max min triangulation.
  // - Version with array cos_angle_threshold[].
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  DeleteTriangleEdgesToMaxMinTriangulationAngle 
  (const std::vector<COS_TYPE2> & cos_angle_threshold, 
   const MTYPE2 max_small_magnitude,
   std::vector<CTYPE2> & vcoord)
  {
    for (NTYPE i = 0; i < cos_angle_threshold.size(); i++) {
      DeleteTriangleEdgesToMaxMinTriangulationAngle 
        (cos_angle_threshold[i], max_small_magnitude, vcoord);
    }
  }


  // Delete edge and add triangle to improve max min triangulation.
  // - Return true if edge is deleted.
  // - Variation of delete edge that cuts off an ear of the quad
  //     containing ihalf_edge.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename ITYPE3,
            typename MTYPE2, typename CTYPE2>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  DeleteEdgeAddTriangleToMaxMinTriangulationAngle
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iear, 
   const MTYPE2 max_small_magnitude, std::vector<CTYPE2> & vcoord)
  {
    const NTYPE NUM_VERTICES_PER_QUAD(4);
    const POLYGON_INDEX_TYPE iquadA =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const VERTEX_INDEX_TYPE jv_ear = this->PolygonVertex(iquadA, iear);
    const VERTEX_INDEX_TYPE iv0 = this->FromVertexIndex(ihalf_edgeA);
    const VERTEX_INDEX_TYPE iv1 = this->ToVertexIndex(ihalf_edgeA);
    const CTYPE2 * vcoord_ptr = IJK::vector2pointer(vcoord);
    const CTYPE2 * v0coord = vcoord_ptr + iv0*Dimension();
    const CTYPE2 * v1coord = vcoord_ptr + iv1*Dimension();
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);
    const POLYGON_INDEX_TYPE ipolyB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeB);
    std::vector<CTYPE2> centroidX_coord(Dimension());
    std::vector<CTYPE2> centroidY_coord(Dimension());
    COS_TYPE cos_min_angle, cos_min_triangle_angle, cosX, cosY;
    bool flag_zero, flagT_zero, flagX_zero, flagY_zero;
    POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE, COS_TYPE, NTYPE> quad_tri_result;
    POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE, COS_TYPE, NTYPE> poly_tri_result;
    IJK::PROCEDURE_ERROR error
      ("DeleteEdgeAddTriangleToMaxMinTriangulationAngle");

    if (!this->IsQuadrilateral(iquadA)) { return(false); }

    const ITYPE j0 = 
      (iear+(NUM_VERTICES_PER_QUAD-1))%NUM_VERTICES_PER_QUAD;
    const ITYPE j1 = (iear+1)%NUM_VERTICES_PER_QUAD;
    VERTEX_INDEX_TYPE jv0 = this->PolygonVertex(iquadA, j0);
    VERTEX_INDEX_TYPE jv1 = this->PolygonVertex(iquadA, j1);

    if (jv_ear == iv0 || jv_ear == iv1) {
      // Vertex jv_ear is an endpoint of proposed deleted edge.
      // Do nothing.
      return false;
    }

    if (this->AreTwoHalfEdgesAroundVertex(jv_ear)) {
      // Creating triangle (iv0,jv_ear,iv1) could create a duplicate triangle.
      return false;
    }

    if (jv1 == iv0 || jv1 == iv1) 
      { std::swap(jv0, jv1); }
    // jv0 now lies on ihalf_edgeA

    IJK::compute_cos_min_triangle_angle
      (Dimension(), vcoord, jv0, jv_ear, jv1, max_small_magnitude,
       cos_min_triangle_angle, flagT_zero);

    compute_mesh2D_polygon_centroid
      (Dimension(), vcoord, *this, ipolyB, centroidX_coord);
    compute_mesh2D_polygon_plus_vertex_centroid
      (Dimension(), vcoord, *this, ipolyB, jv1, centroidY_coord);

    // Merge triangle jv1 and endpoints of ihalf_edgeB
    //   and polygon containing ihalf_edgeB.
    compute_cos_min_merge_triangle_polygon_triangulation_angle
      (Dimension(), vcoord, centroidX_coord, *this,
       jv1, ihalf_edgeB, max_small_magnitude, cosX, flagX_zero);

    // Merge again, but use centroidC_coord.
    compute_cos_min_merge_triangle_polygon_triangulation_angle
      (Dimension(), vcoord, centroidY_coord, *this,
       jv1, ihalf_edgeB, max_small_magnitude, cosY, flagY_zero);

    ITYPE index_selected;
    IJK::select_min(cosX, flagX_zero,
                    cosY, flagY_zero,
                    cos_min_angle, index_selected, flag_zero);

    flag_zero = (flag_zero || flagT_zero);
    cos_min_angle = std::max(cos_min_angle, cos_min_triangle_angle);

    if (!flag_zero) {
      if (cos_min_angle < PolygonCosMinTriangulationAngle(iquadA)) {

        POLYGON_INDEX_TYPE itriangle_new, ipoly_new;

        this->DeleteEdgeAddTriangle
          (ihalf_edgeA, iear, itriangle_new, ipoly_new);
        SetCosMinTriangulationAngle(itriangle_new, cos_min_triangle_angle);

        if (index_selected == 0) 
          { SetInteriorTriVertex(ipoly_new, cosX, centroidX_coord); }
        else
          { SetInteriorTriVertex(ipoly_new, cosY, centroidY_coord); }

        return(true);
      }
    }

    return(false);
  }


  // Delete (long) quad edges to improve max min triangulation.
  // - Calls DeleteEdgeAddTriangleToMaxMinTriangulationAngle().
  // @pre ComputeCosMaxMinAngle() or ComputeCosMaxMinAngleAllowCentroid()
  //    has been called on all polygons to determine triangulation.
  // @param cos_angle_threshold Attempt to delete longest quad edge (u,v)
  //    if quad angle at u or at v is less than 
  //    cos^{-1}(cos_angle_threshold).
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  DeleteQuadEdgesToMaxMinTriangulationAngle
  (const COS_TYPE2 cos_angle_threshold, const MTYPE2 max_small_magnitude,
   std::vector<CTYPE2> & vcoord)
  {
    const NTYPE NUM_VERTICES_PER_QUAD(4);
    const NUMBER_TYPE num_poly = this->NumPolygons();
    COS_TYPE2 cos_min_angle;
    ITYPE iangle;
    bool flag_zero;

    for (ITYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (this->IsPolygonDeleted(ipoly)) { continue; }

      if (!this->IsQuadrilateral(ipoly)) { continue; }

      if (PolygonCosMinTriangulationAngle(ipoly) <= cos_angle_threshold)
        { continue; }

      const HALF_EDGE_INDEX_TYPE ilongest_half_edge = 
        compute_longest_polygon_half_edge
        (Dimension(), vcoord, *this, ipoly);

      IJK::compute_cos_min_quad_angle
        (Dimension(), vcoord, this->PolygonVertexList(ipoly), 
         max_small_magnitude, cos_min_angle, iangle, flag_zero);

      if (cos_min_angle <= cos_angle_threshold) { continue; }


      const ITYPE iear = (iangle+2)%NUM_VERTICES_PER_QUAD;
      DeleteEdgeAddTriangleToMaxMinTriangulationAngle
        (ilongest_half_edge, iear, max_small_magnitude, vcoord);
    }
  }


  // Delete (long) quad edges to improve max min triangulation.
  // - Calls DeleteEdgeAddTriangleToMaxMinTriangulationAngle().
  // - Version with array cos_angle_threshold[].
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  DeleteQuadEdgesToMaxMinTriangulationAngle
  (const std::vector<COS_TYPE2> & cos_angle_threshold, 
   const MTYPE2 max_small_magnitude,
   std::vector<CTYPE2> & vcoord)
  {
    for (NTYPE i = 0; i < cos_angle_threshold.size(); i++) {
      DeleteQuadEdgesToMaxMinTriangulationAngle 
        (cos_angle_threshold[i], max_small_magnitude, vcoord);
    }
  }


  // *****************************************************************
  // MESH2D_TRIANGULATION_BASE replace triangle member functions
  // *****************************************************************

  // Replace triangle with triangle-edge to improve 
  //   max min triangulation angle.
  // - Return true if triangle is replaced.
  // @param ihalf_edge Triangle base.
  //    - Replace vertex opposite triangle base.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename MTYPE2, typename CTYPE2>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  ReplaceTriangleWithTriangleEdgeToMaxMinTriangulationAngle
  (const ITYPE2 ihalf_edge0, const MTYPE2 max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroid2,
   const bool flag_max_min_all, std::vector<CTYPE2> & vcoord)
  {
    const NUMBER_TYPE THREE(3);

    const POLYGON_INDEX_TYPE itriangleB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1 = this->IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge2 = this->IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    const POLYGON_INDEX_TYPE ipolyA =
      this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge1);
    const POLYGON_INDEX_TYPE ipolyC =
      this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge2);
    POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE, COS_TYPE, NTYPE> poly_tri_result[3];
    COS_TYPE cos_min_angle;
    bool flag_zero;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    this->PrintPolygonIndexAndVertices(cerr, "  Triangle: ", itriangleB, "\n");
    */

    std::vector<CTYPE2> midpoint_coord(Dimension());
    std::vector<CTYPE2> centroidAI_coord(Dimension());
    std::vector<CTYPE2> centroidCI_coord(Dimension());

    if (this->IsBoundaryEdge(ihalf_edge1)) { return(false); }
    if (this->IsBoundaryEdge(ihalf_edge2)) { return(false); }
    if (this->AreThreeHalfEdgesAroundFromVertex(ihalf_edge2)) {
      // Polygons ipolyA and ipolyC already share an edge.
      // Replacing the triangle would cause them to share two edges.
      return(false);
    }

    const COS_TYPE cos_min_angle_poly3 = 
      PolygonIIICosMinTriangulationAngle(ipolyA, itriangleB, ipolyC);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  CosMinAngle(" << ipolyA << "): "
         << PolygonCosMinTriangulationAngle(ipolyA) << endl;
    cerr << "  CosMinAngle(" << itriangleB << "): "
         << PolygonCosMinTriangulationAngle(itriangleB) << endl;
    cerr << "  CosMinAngle(" << ipolyC << "): "
         << PolygonCosMinTriangulationAngle(ipolyC) << endl;
    PolygonTriangulationInfo(ipolyC).Print(cerr, "    ");
    */
    

    ComputeMidpointOfEdgeMidpoints
      (vcoord, ihalf_edge1, ihalf_edge2, midpoint_coord);
    ComputePolygonCentroid(vcoord, ipolyA, centroidAI_coord);
    ComputePolygonCentroid(vcoord, ipolyC, centroidCI_coord);

    if (flag_allow_centroid2) {

      std::vector<CTYPE2> centroidAII_coord(Dimension());
      std::vector<CTYPE2> centroidCII_coord(Dimension());
      std::vector<CTYPE2> selectedA_coord(Dimension());
      std::vector<CTYPE2> selectedC_coord(Dimension());

      compute_mesh2D_polygon_plus_point_centroid
        (Dimension(), vcoord, midpoint_coord, *this, ipolyA, centroidAII_coord);
      compute_mesh2D_polygon_plus_point_centroid
        (Dimension(), vcoord, midpoint_coord, *this, ipolyC, centroidCII_coord);

      // Replace triangle itriangleB with triangle and edge.
      compute_cos_max_min_replace_triangle_triangulation_angle_IVx2
        (Dimension(), vcoord, midpoint_coord, 
         centroidAI_coord, centroidAII_coord, 
         centroidCI_coord, centroidCII_coord, 
         *this, ihalf_edge0, max_small_magnitude, 
         triangulation_settings,
         poly_tri_result, cos_min_angle, flag_zero,
         selectedA_coord, selectedC_coord);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "  cos_min_angle: " << cos_min_angle << endl;
      cerr << "  CosMinAngle(" << itriangleB << "): "
           << PolygonCosMinTriangulationAngle(itriangleB) << endl;
      cerr << "  cos_min_angle_poly3: " << cos_min_angle_poly3 << endl;
      cerr << "  flag_max_min_all: " << int(flag_max_min_all) << endl;
      */

      if (!flag_zero) {
        if (cos_min_angle < PolygonCosMinTriangulationAngle(itriangleB) ||
            (flag_max_min_all &&
             cos_min_angle < cos_min_angle_poly3)) {

          AddVertexReplaceTriangleWithTriangleEdge
            (midpoint_coord, selectedA_coord, selectedC_coord,
             vcoord, ihalf_edge0, poly_tri_result);

          return(true);
        }
      }
    }
    else {

      // Replace triangle itriangleB with triangle and edge.
      compute_cos_max_min_replace_triangle_triangulation_angle
        (Dimension(), vcoord, midpoint_coord, 
         centroidAI_coord, centroidCI_coord, 
         *this, ihalf_edge0, max_small_magnitude, 
         triangulation_settings, flag_allow_centroid2,
         poly_tri_result, cos_min_angle, flag_zero);

      // *** DEBUG ***
      /*
      using namespace std;
      this->PrintPolygonIndexAndVertices
        (cerr, "Triangle: ", itriangleB, "\n");
      this->PrintPolygonIndexAndVertices
        (cerr, "ipolyA: ", ipolyA, "\n");
      poly_tri_result[0].Print(cerr, "  ");
      this->PrintPolygonIndexAndVertices
        (cerr, "ipolyC: ", ipolyC, "\n");
      poly_tri_result[2].Print(cerr, "  ");
      */


      if (!flag_zero) {
        if (cos_min_angle < PolygonCosMinTriangulationAngle(itriangleB) ||
            (flag_max_min_all &&
             cos_min_angle < cos_min_angle_poly3)) {

          AddVertexReplaceTriangleWithTriangleEdge
            (midpoint_coord, centroidAI_coord, centroidCI_coord,
             vcoord, ihalf_edge0, poly_tri_result);

          return(true);
        }
      }
    }

    return(false);
  }


  // *****************************************************************
  // MESH2D_TRIANGULATION_BASE replace triangle with split edges
  // *****************************************************************

  // Replace triangle that has two split edges with triangle-edge 
  //   to improve max min triangulation angle.
  // - Return true if triangle is replaced.
  // - Allow triangulations from polygon centroids.
  // - Consider two possible centroids for each adjacent polygon.
  //   Either centroid of polygon or centroid of polygon and the new
  //   triangle vertex.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPE2, typename MTYPE2, typename CTYPE2>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  ReplaceSplitEdge2TriangleToMaxMinTriangulationAngle
  (const ITYPE2 ihalf_edge0, const MTYPE2 max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_max_min_all, std::vector<CTYPE2> & vcoord)
  {
    const NUMBER_TYPE THREE(3);

    const POLYGON_INDEX_TYPE itriangleB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    HALF_EDGE_INDEX_TYPE ihalf_edge1, ihalf_edge2, ihalf_edge3;

    if (!this->IsPentagon(itriangleB)) { return(false); }

    this->GetIndicesOfNextThreeHalfEdgesInPolygon
      (ihalf_edge0, ihalf_edge1, ihalf_edge2, ihalf_edge3);

    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1X = this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge3X = this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge3);
    const POLYGON_INDEX_TYPE ipolyA =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge1X);
    const POLYGON_INDEX_TYPE ipolyC =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge3X);
    const HALF_EDGE_INDEX_TYPE iv_split0 =
      this->FromVertexIndex(ihalf_edge2);
    const HALF_EDGE_INDEX_TYPE iv_split1 =
      this->ToVertexIndex(ihalf_edge3);

    POLYGON_TRIANGULATION_RESULT_E<BIT_SET_SIZE,COS_TYPE,NTYPE> 
      poly_tri_result[3];
    COS_TYPE cos_min_angle, cosA, cosC;
    bool flag_zero, flagA_zero, flagC_zero;

    std::vector<CTYPE2> midpoint_coord(Dimension());
    std::vector<CTYPE2> centroidA_coord(Dimension());
    std::vector<CTYPE2> centroidC_coord(Dimension());
    

    if (this->IsBoundaryEdge(ihalf_edge1)) { return(false); }
    if (this->IsBoundaryEdge(ihalf_edge3)) { return(false); }
    if (this->AreThreeHalfEdgesAroundFromVertex(ihalf_edge3)) {
      // Polygons ipolyA and ipolyC already share an edge.
      // Replacing the triangle would cause them to share two edges.
      return(false);
    }

    ComputeLineSegmentMidpoint
      (vcoord, iv_split0, iv_split1, midpoint_coord);
    ComputePolygonCentroid(vcoord, ipolyA, centroidA_coord);
    ComputePolygonCentroid(vcoord, ipolyC, centroidC_coord);

    // Replace triangle itriangleB with triangle and edge.
    compute_cos_max_min_replace_splitE2T_triangulation_angle
      (Dimension(), vcoord, midpoint_coord, 
       centroidA_coord, centroidC_coord, 
       *this, ihalf_edge0, max_small_magnitude, 
       triangulation_settings, poly_tri_result, 
       cos_min_angle, flag_zero);

    COS_TYPE cos_min_angle_poly3 =
      PolygonIIICosMinTriangulationAngle(ipolyA, itriangleB, ipolyC);

    if (!flag_zero) {
      if (cos_min_angle < PolygonCosMinTriangulationAngle(itriangleB) ||
          (flag_max_min_all &&
           cos_min_angle < cos_min_angle_poly3)) {

        // *** DEBUG ***
        /*
        using namespace std;
        this->PrintPolygonIndexAndVertices
          (cerr, "*** Replace triangle ", itriangleB, "\n");
        this->PrintPolygonIndexAndVertices
          (cerr, "  ipolyA: ", ipolyA, "\n");
        this->PrintPolygonIndexAndVertices
          (cerr, "  ipolyC: ", ipolyC, "\n");
        */

        AddVertexReplaceSplit2EdgeTriangle
          (midpoint_coord, centroidA_coord, centroidC_coord,
           vcoord, ihalf_edge0, poly_tri_result);

        return(true);
      }
    }

    return(false);
  }


  // Replace triangles that two have split edges to improve 
  //   max min triangulation angle.
  // @param flag_allow_centroid  If true, allow adding triangulation
  //     vertices at centroids of adjacent polygons.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  ReplaceSplitEdge2TrianglesToMaxMinTriangulationAngle
  (const COS_TYPE2 cos_angle_threshold, 
   const MTYPE2 max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   const bool flag_triangle_has_min_angle,
   const bool flag_max_min_all, 
   std::vector<CTYPE2> & vcoord)
  {
    const NTYPE NUM_VERTICES_PER_TRIANGLE(3);
    const NTYPE TWO(2);
    const NTYPE num_poly = this->NumPolygons();
    HALF_EDGE_INDEX_TYPE ihalf_edge0, ihalf_edge1, ihalf_edge2, ihalf_edge3;
    bool flag_zero;

    for (ITYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (this->IsPolygonDeleted(ipoly)) { continue; }
      if (!IsTriangleWithTwoSplitEdges(ipoly, ihalf_edge0))
        { continue; }

      this->GetIndicesOfNextThreeHalfEdgesInPolygon
        (ihalf_edge0, ihalf_edge1, ihalf_edge2, ihalf_edge3);

      if (this->IsBoundaryEdge(ihalf_edge1) ||
          this->IsBoundaryEdge(ihalf_edge3)) {
        // Don't replace triangle if a split edge is on the boundary.
        continue;
      }

      const HALF_EDGE_INDEX_TYPE ihalf_edge4 =
        this->IndexOfNextHalfEdgeInPolygon(ihalf_edge3);
      if (!this->AreTwoHalfEdgesAroundFromVertex(ihalf_edge2) ||
          !this->AreTwoHalfEdgesAroundFromVertex(ihalf_edge4)) {
        // Only replace vertex if both split vertices are incident
        //  on exactly two polygons.
        continue;
      }

      const POLYGON_INDEX_TYPE ipolyA =
        this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge1);
      const POLYGON_INDEX_TYPE ipolyC =
        this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge3);

      const COS_TYPE2 cos_min_angle =
        this->PolygonCosMinTriangulationAngle(ipoly);

      if (flag_triangle_has_min_angle) {
        if (cos_min_angle <
            PolygonCosMinTriangulationAngle(ipolyA)) { continue; }

        if (cos_min_angle <
            PolygonCosMinTriangulationAngle(ipolyC)) { continue; }

        if (cos_min_angle <= cos_angle_threshold)
          { continue; } 
      }
      else {
        COS_TYPE2 cos_min_angle_poly3 = 
          std::max(PolygonCosMinTriangulationAngle(ipolyA),
                   cos_min_angle);
        cos_min_angle_poly3 = 
          std::max(cos_min_angle_poly3,
                   PolygonCosMinTriangulationAngle(ipolyC));

        if (cos_min_angle_poly3 <= cos_angle_threshold) 
          { continue; }
      }

      ReplaceSplitEdge2TriangleToMaxMinTriangulationAngle
        (ihalf_edge0, max_small_magnitude, triangulation_settings,
         flag_allow_centroidx2, flag_max_min_all, vcoord);
    }
  }


  // Replace triangles to improve max min triangulation angle.
  // - Version with array cos_angle_threshold[].
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename COS_TYPE2, typename MTYPE2, typename CTYPE2>
  void _MESH2D_TRIANGULATION_BASE_T_NAME_::
  ReplaceSplitEdge2TrianglesToMaxMinTriangulationAngle
  (const std::vector<COS_TYPE2> & cos_angle_threshold, 
   const MTYPE2 max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroid2,
   const bool flag_triangle_has_min_angle,
   const bool flag_max_min_all,
   std::vector<CTYPE2> & vcoord)
  {
    for (NTYPE i = 0; i < cos_angle_threshold.size(); i++) {
      ReplaceSplitEdge2TrianglesToMaxMinTriangulationAngle
        (cos_angle_threshold[i], max_small_magnitude, 
         triangulation_settings, flag_allow_centroid2,
         flag_triangle_has_min_angle, flag_max_min_all, vcoord);
    }
  }



  // *****************************************************************
  // MESH2D_TRIANGULATION_BASE check routines
  // *****************************************************************

  // Check that two adjacent polygons are different.
  // - Return true if two adjacent polygons are different.
  // - Return true if ihalf_edgeA or ihalf_edgeB are boundary edges.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPEH>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  CheckAreAdjacentPolygonsDifferent
  (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
   IJK::ERROR & error) const
  {
    if (AreAdjacentPolygonsDifferent(ihalf_edgeA, ihalf_edgeB)) 
      { return(true); }
    else {
      error.AddMessage
        ("Programming error. Two edges are incident on same polygons.");
      error.AddMessage
        ("  Polygons sharing two edges: ", 
         this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA), ", ",
         this->IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edgeA),
         "");
      return(false);
    }
  }


  // Check that three adjacent polygons are different.
  // - Return true if three adjacent polygons are different.
  template _MESH2D_TRIANGULATION_BASE_T_PARAM_
  template <typename ITYPEH>
  bool _MESH2D_TRIANGULATION_BASE_T_NAME_::
  CheckAreAdjacentPolygonsDifferent
  (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
   const ITYPEH ihalf_edgeC, IJK::ERROR & error) const
  {
    if (!CheckAreAdjacentPolygonsDifferent
        (ihalf_edgeA, ihalf_edgeB, error))
      { return(false); }

    if (!CheckAreAdjacentPolygonsDifferent
        (ihalf_edgeA, ihalf_edgeC, error))
      { return(false); }

    if (!CheckAreAdjacentPolygonsDifferent
        (ihalf_edgeB, ihalf_edgeC, error))
      { return(false); }

    return(true);
  }

}

#endif

