/*!
 *  @file ijkmesh2D_split.tpp
 *  @brief ijk template classes for splitting edges and polygons 
 *         in 2D polygonal mesh data structures. 
 *  @details
 *  - Uses half-edges.  Totally distinct from POLYMESH.
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

#ifndef _IJKMESH2D_SPLIT_POLY_
#define _IJKMESH2D_SPLIT_POLY_

#include "ijk.tpp"
#include "ijklist.tpp"
#include "ijkmesh2D_datastruct.tpp"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>


// *** DEBUG ***
#include "ijkprint.tpp"


namespace IJK {

  // *****************************************************************
  // Class MESH2D_SPLIT_BASE
  // *****************************************************************

  /**
     @brief Mesh of 2D polygons supporting polygon splits/replacement.
     @tparam VTYPE Vertex type.  Should be derived from VERTEX_BASE.
     @tparam HALFE_TYPE Half edge type.  
       - Should be derived from HALF_EDGE_BASE.
     @tparam PTYPE POLYGON type.  
       - Should be derived from POLYGON_BASE.
  */
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  class MESH2D_SPLIT_BASE:
    public MESH2D_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE> {

  public:

    typedef MESH2D_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>
    MESH2D_BASE_TYPE;

    typedef typename MESH2D_BASE_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;
    typedef typename MESH2D_BASE_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH2D_BASE_TYPE::HALF_EDGE_INDEX_TYPE 
    HALF_EDGE_INDEX_TYPE;
    typedef typename MESH2D_BASE_TYPE::INDEX_TYPE INDEX_TYPE;
    typedef typename MESH2D_BASE_TYPE::NUMBER_TYPE NUMBER_TYPE;


  protected:

    /// @brief Initialize.
    void Init();

    // Replace vertex.

    /// @brief Replace from vertex.
    template <typename ITYPE2, typename ITYPE3>
    void ReplaceFromVertex
    (const ITYPE2 ihalf_edge, const ITYPE3 iv_new);

    /// @brief Replace to vertex.
    template <typename ITYPE2, typename ITYPE3>
    void ReplaceToVertex
    (const ITYPE2 ihalf_edge, const ITYPE3 iv_new);

    /**
       @brief Split polygon into two polygons.
         - Split polygon with diagonal 
           (FromVertexIndex(ihalf_edge0), FromVertexIndex(ihalf_edge1)).
    */
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5>
    void SplitPolygonByDiagonal
    (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1,
     ITYPE4 & ihalf_edge_new0, ITYPE5 & ihalf_edge_new1);


    // Split edges or polygons.

    /**
      @brief Split half edge ihalf_edge of polygon, creating new polygon.
      @param[out] ihalf_edge_new Index of first new half edge that
                  replaces ihalf_edge0.
        - IndexOfNextHalfEdgeInPolygon(ihalf_edge_new) is index 
          of the second half edge that replaces ihalf_edge0.
    */
    // *** CHANGE TO RETURN new half edge instead of polygon ***
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    POLYGON_INDEX_TYPE SplitPolygonHalfEdge
    (const ITYPE2 ihalf_edge0, const ITYPE3 iv_split,
     ITYPE4 & ihalf_edge_new);

    /**
       @brief Split half edge ihalf_edge of polygon, creating new polygon,
              and link two new half edges to ihalf_edgeXC and ihalf_edgeXD.
       @param ihalf_edgeXC First half edge in opposite direction.
       @param ihalf_edgeXD Second half edge in opposite direction.
         - ToVertexIndex(ihalf_edgeXC) = iv_split.
         - FromVertexIndex(ihalf_edgeXD) = iv_split.
       @param[out] ihalf_edge_newA Index of first new half edge that
                   replaces ihalf_edgeA.
         - IndexOfNextHalfEdgeInPolygon(ihalf_edge_newA) is index 
           of the second half edge that replaces ihalf_edgeA.
    */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4,
              typename ITYPE5, typename ITYPE6>
    POLYGON_INDEX_TYPE SplitAndLinkPolygonHalfEdge
    (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split,
     const ITYPE4 ihalf_edgeXC, const ITYPE5 ihalf_edgeXD,
     ITYPE6 & ihalf_edge_newA);

    /** 
       @brief Split polygon half edges at given locations.
       - Return index of new polygon.
       - First vertex in new polygon equals first vertex
         in original polygon.
       - Order of vertices (clockwise or counter-clockwise)
         in new polygon matches order of vertices in original polygon.
       @param isplit_loc[] Location of split half edges.
         - Values are in range [0..(NumPolygonEdges(ipoly)-1)].
         @pre isplit_loc[] are sorted in increasing order.
       @param iv_split[] Array of split vertices.
         - iv_split[j] is vertex splitting half edge at isplit_loc[j].
       @param num_split Number of half edges to split.
         @pre Must be at least 1.
    */
    template <typename ITYPEP, typename ITYPEL, 
              typename ITYPEV, typename NTYPE2>
    POLYGON_INDEX_TYPE SplitPolygonHalfEdgesAtLoc
    (const ITYPEP ipoly, const ITYPEL isplit_loc[], 
     const ITYPEV iv_split[], const NTYPE2 num_split);

    /** 
       @brief Split polygon half edges at given locations.
       - Return index of new polygon.
       - Version using C++ STL vector isplit_loc[] and iv_split[].
       @param isplit_loc[] Location of split half edges.
       @param iv_split[] Array of split vertices.
         @pre iv_split.size() equals isplit_loc.size().
    */
    template <typename ITYPEP, typename ITYPEL, typename ITYPEV>
    POLYGON_INDEX_TYPE SplitPolygonHalfEdgesAtLoc
    (const ITYPEP ipoly, const std::vector<ITYPEL> & isplit_loc, 
     const std::vector<ITYPEV> & iv_split);

    /**
       @brief Split two polygon half edges.
       - Return index of new polygon.
    */
    template <typename ITYPEHA, typename ITYPEHB, 
              typename ITYPEVA, typename ITYPEVB,
              typename ITYPEHA2, typename ITYPEHB2>
    POLYGON_INDEX_TYPE SplitTwoHalfEdges
    (const ITYPEHA ihalf_edgeA, const ITYPEHB ihalf_edgeB,
     const ITYPEVA iv_splitA, const ITYPEVB iv_splitB,
     ITYPEHA2 & ihalf_edgeA_new, ITYPEHB2 & ihalf_edgeB_new);

    /**
       @brief Split edges of one polygon.
       - Return indices of new polygons created by splitting edges.
       @param ihalf_edge[] List of half edges to split.
         - All half edges must be in same polygon.
       @param iv_split[] List of split vertices.
       @param num_split Number of half edges to split.
         @pre Must be at least 1.
       @param[out] ipoly_new[] Indices of the new polygons.
       - Number of new polygons equals num_split+1.
       - ipoly_new[0] Polygon containing half edges in ihalf_edge[].
       - ipoly_new[j] Polygons with split edges adjacent to ipoly_new[0].
    */
    template <typename ITYPEH, typename ITYPEV, typename ITYPEP,
              typename NTYPES>
    void SplitEdgesOfPolygon
    (const ITYPEH ihalf_edge[], const ITYPEV iv_split[],
     const NTYPES num_split, ITYPEP ipoly_new[]);

    /**
       @brief Split edges of one polygon.
       - Return indices of new polygons created by splitting edges.
       @param ihalf_edge[] List of half edges to split.
         - All half edges must be in same polygon.
       @param iv_split[] List of split vertices.
       @param num_split Number of half edges to split.
         @pre Must be at least 1.
       @param[out] ihalf_edge_new[] List of indices of some
         of the new half edges.
       - Number of elements of ihalf_edge_new[] is num_split.
       - ihalf_edge_new[j] and NextHalfEdgeInPolygon(ihalf_edge_new[j])
         are the two new half edges replacing ihalf_edge[j].
    */
    template <typename ITYPEH, typename ITYPEV, typename NTYPES,
              typename ITYPEH2>
    void SplitEdgesOfPolygonNew
    (const ITYPEH ihalf_edge[], const ITYPEV iv_split[],
     const NTYPES num_split, ITYPEH2 ihalf_edge_new[]);

    /**
       \brief Add vertices and split half edges and polygons.
       @param ihalf_edge[] List of half edges to split.
         - All half edges must be in same polygon.
       @param num_split Number of half edges to split.
         @pre Must be at least 1.
       @param[out] ipoly_new[i] Index of i'th new polygon.
         @pre Must be preallocated to length at least num_split+1.
         - ipoly_new[0] is the new polygon replacing
           the original polygon containing half_edge[0].
    */
    template <typename ITYPEH, typename ITYPEP, typename NTYPES>
    void AddVerticesSplitEdgesOfPolygon
    (const ITYPEH ihalf_edge[], const NTYPES num_split,
     ITYPEP ipoly_new[]);

    /**
       \brief Split and link edge of polygon adjacent to split triangle.
       @param ivloc0 Location (0, 1 or 2) of to vertex of ihalf_edgeX.
       @param ihalf_edgeX Split edge of polygon.
       @param iv_split Index of edge split vertex.
       @param itriangle_new[] Four triangles forming split triangle.
         - itriangle_new[j] contains j'th triangle vertex.
    */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4, 
              typename ITYPE5, typename ITYPE6, typename ITYPE7>
    void _SplitAndLinkEdgeOfPolygonAdjacentToSplitTriangle
    (const ITYPE2 ivloc0, const ITYPE3 ihalf_edgeX, 
     const ITYPE4 iv_split, const ITYPE5 itriangle_new[4], 
     ITYPE6 & ipoly_new, ITYPE7 & ihalf_edgeX_new);


    // Cut ear from polygon.

    /**    
       @brief Cut ear from polygon.
         - Separate ear vertex FromVertexIndex(ihalf_edge1) from the
           rest of the polygon.
       @param[out] ihalf_edge1_new New edge replacing ihalf_edge1.
          - Note: Input half edge to SplitPolygonEar is different
            from input half edges to SplitPolygonByDiagonal().
       @param[out] ihalf_edge2_new New edge 
         IndexOfNextEdgeInPolygon(ihalf_edge1).
         - ihalf_edge2_new is in different polygon from ihalf_edge1_new.
    */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    void CutPolygonEar
    (const ITYPE2 ihalf_edge1, ITYPE3 & ihalf_edge1_new, 
     ITYPE4 & ihalf_edge2_new);


    /**
       @brief Replace triangle.
       @param ihalf_edge0 Base edge of triangle.
         - Replace vertex opposite ihalf_edge0 with iv_new.
       @param[out] Replace ihalf_edge0 with ihalf_edge_new0.
    */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    POLYGON_INDEX_TYPE ReplaceTriangle
    (const ITYPE2 ihalf_edge0, const ITYPE3 iv_new, 
     ITYPE4 & ihalf_edge_new0);

    /**
       @brief Replace pentagon with triangle.
         - Merge two pentagon vertices to form triangle.
       @param ihalf_edge0 Pentagon edge that becomes base edge of triangle.
       @param iv_new Replace merged vertices with iv_new.
    */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    POLYGON_INDEX_TYPE ReplacePentagonWithTriangle
    (const ITYPE2 ihalf_edge0, const ITYPE3 iv_new, 
     ITYPE4 & ihalf_edge_new0);


  public:
    /// Constructor.
    MESH2D_SPLIT_BASE() { Init(); };

    // Query functions.

    /**
       @brief Location of j'th vertex in polygon created by splitting half edge ihalf_edge0.
       - ToVertex(ihalf_edgeA) is 0'th vertex of new polygon.
    */
    template <typename ITYPE2, typename ITYPE3>
    INDEX_TYPE LocationInSplitEdgePolygon
    (const ITYPE2 ihalf_edgeA, const ITYPE3 j) const;


    /// Copy.
    template <typename MTYPE>
    void Copy(const MTYPE & meshA);

    // Split functions.

    /**
       @brief Split edge contained in two polygon.
       - If edge in only one polygon, split one polygon.
       @param[out] ihalf_edgeA_new Index of first new half edge that
                   replaces ihalf_edgeA.
         - IndexOfNextHalfEdgeInPolygon(ihalf_edgeA_new) is index 
           of the second half edge that replaces ihalf_edgeA.
    */
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5, typename ITYPE6>
    void SplitPolygonEdgeII
    (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split,
     ITYPE4 & ipoly0, ITYPE5 & ipoly1, ITYPE6 & ihalf_edgeA_new);

    /**
       @brief Split edge contained in two polygons. 
              (Does not return ihalf_edgeA_new).
       - Version that does not return ihalf_edgeA_new.
    */
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5>
    void SplitPolygonEdgeII
    (const ITYPE2 ihalf_edge0, const ITYPE3 iv_split,
     ITYPE4 & ipoly0, ITYPE5 & ipoly1);

    /**
       @brief Split two edges of a polygon.
         - Splits half edges in 3 polygons.
    */
    template <typename ITYPEH, typename ITYPEV, typename ITYPEH2>
    void SplitTwoEdgesOfPolygon
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     const ITYPEV ivA_split, const ITYPEV ivB_split,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new);

    /**
       @brief Split three edges of a polygon.
         - Splits half edges in 4 polygons.
    */
    template <typename ITYPEH, typename ITYPEV, typename ITYPEH2>
    void SplitThreeEdgesOfPolygon
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     const ITYPEH ihalf_edgeC,
     const ITYPEV ivA_split, const ITYPEV ivB_split,
     const ITYPEV ivC_split,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new,
     ITYPEH2 & ihalf_edgeC_new);

    /**
       @brief Split three edges of a  polygon.
         - Splits half edges in 4 polygons.
       @param ipolyABC Index of new polygon replacing polygon incident
         on ihalf_edgeA, ihalf_edgeB and ihalf_edgeC.
       @param ipolyA Index of new polygon replacing polygon incident
         on ihalf_edgeA.
       @param ipolyB Index of new polygon replacing polygon incident
         on ihalf_edgeB.
       @param ipolyC Index of new polygon replacing polygon incident
         on ihalf_edgeC.
    */
    template <typename ITYPEH, typename ITYPEV, typename ITYPEP>
    void SplitThreePolygonEdgesIV
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     const ITYPEH ihalf_edgeC,
     const ITYPEV ivA_split, const ITYPEV ivB_split,
     const ITYPEV ivC_split,
     ITYPEP & ipolyABC, ITYPEP & ipolyA, ITYPEP & ipolyB, ITYPEP & ipolyC);

    /**
       @brief Split three edges incident on the same vertex.
       @pre ihalf_edgeB1 is not a boundary edge.
          - Modifies four polygons.
    */
    template <typename ITYPE2, 
              typename ITYPEA, typename ITYPEB, typename ITYPEC,
              typename ITYPEP>
    void SplitThreeIncidentEdgesIV
    (const ITYPE2 ihalf_edgeB1, 
     const ITYPEA iv_splitA, const ITYPEB iv_splitB, 
     const ITYPEC iv_splitC,
     ITYPEP ipoly[4]);

    /**
       @brief Add vertex and split edge contained in two polygon.
       @param[out] ihalf_edgeA_new Index of first new half edge that
                   replaces ihalf_edgeA.
         - IndexOfNextHalfEdgeInPolygon(ihalf_edgeA_new) is index 
           of the second half edge that replaces ihalf_edgeA.
    */
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5, typename ITYPE6>
    void AddVertexSplitPolygonEdgeII
    (const ITYPE2 ihalf_edgeA, ITYPE3 & iv_split,
     ITYPE4 & ipoly0, ITYPE5 & ipoly1, ITYPE6 & ihalf_edgeA_new);

    /**
       @brief Add vertex and split edge contained in two polygon.
              (Does not return ihalf_edgeA_new.)
         - Version that does not return ihalf_edgeA_new.
    */
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5>
    void AddVertexSplitPolygonEdgeII
    (const ITYPE2 ihalf_edge0, ITYPE3 & iv_split,
     ITYPE4 & ipolyA, ITYPE5 & ipolyB);

    /**
       @brief Add two vertices and split two edges of a polygon.
         - Splits half edges in 3 polygons.
    */
    template <typename ITYPEH, typename ITYPEV, typename ITYPEH2>
    void AddTwoVerticesSplitTwoEdgesOfPolygon
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     ITYPEV & ivA_split, ITYPEV & ivB_split,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new);

    /**
       @brief Add three vertices and split three edges of a polygon.
         - Splits half edges in 4 polygons.
    */
    template <typename ITYPEH, typename ITYPEV, typename ITYPEH2>
    void AddThreeVerticesSplitThreeEdgesOfPolygon
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     const ITYPEH ihalf_edgeC,
     ITYPEV & ivA_split, ITYPEV & ivB_split,
     ITYPEV & ivC_split,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new,
     ITYPEH2 & ihalf_edgeC_new);

    /**
       @brief Add vertices and split three edges that are incident
              on the same vertex.
         - Modifies four polygons.
    */
    template <typename ITYPE2, 
              typename ITYPEA, typename ITYPEB, typename ITYPEC, 
              typename ITYPEP>
    void AddVerticesSplitThreeIncidentEdgesIV
    (const ITYPE2 ihalf_edgeB1, 
     ITYPEA & iv_splitA, ITYPEB & iv_splitB, ITYPEC & iv_splitC, 
     ITYPEP ipoly[4]);

    /**
       @brief Split triangle and its half edges.
         - Split triangle into four subtriangles.
    */
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4>
    void SplitTriangleAndHalfEdges
    (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split[3], 
     ITYPE4 itriangle_new[4]);

    /**
       @brief Split triangle and its three edges.
         - Splits four polygons.
    */
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5>
    void SplitTriangleAndEdgesIV
    (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split[3], 
     ITYPE4 itriangle_new[4], ITYPE5 ipoly_new[3]);

    /// @brief Add vertices and split triangle and its edges.
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5>
    void AddVerticesSplitTriangleAndEdgesIV
    (const ITYPE2 ihalf_edgeA, ITYPE3 iv_split[3], 
     ITYPE4 itriangle_new[4], ITYPE5 ipoly_new[3]);

    /**
       @brief Split pentagon half edge and triangulate resulting hexagon
              into 4 triangles.
       @pre Polygon containing ihalf_edgeA is a pentagon.
         - Hexagon triangulation has 3 ears and one "internal" triangle.
         - Vertex iv_split is incident on the "internal" triangle.
    */
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4>
    void SplitPentagonHalfEdgeAndPentagon
    (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split, 
     ITYPE4 itriangle_new[4]);
    

    // Delete edge functions.

    /**
       @brief Delete edge contained in two polygon.
       @pre Edge ihalf_edge0 is not a boundary edge.
    */
    template <typename ITYPE2, typename ITYPE3>
    void DeleteEdge
    (const ITYPE2 ihalf_edgeA, ITYPE3 & ipoly_new);

    /**
       @brief Delete edge and add triangle.
         @pre Edge ihalf_edgeA is not a boundary edge.
         @pre Quadrilateral contains ihalf_edgeA.
       @param iear Location of vertex in quad containing ihalf_edgeA.
       @param[out] itriangle_new Index of new triangle.
    */
    template <typename ITYPE2, typename ITYPE3, 
              typename ITYPE4, typename ITYPE5>
    void DeleteEdgeAddTriangle
    (const ITYPE2 ihalf_edgeA, const ITYPE3 iear,
     ITYPE4 & itriangle_new, ITYPE5 & ipoly_new);


    // Replace triangle function.

    /**
       @brief Replace triangle with triangle edge.
       @param ihalf_edge0 Base of triangle.
         - Replace vertex opposite ihalf_edge0 with iv_new.
       @param iv_new New triangle vertex.
       @param[out] Replace ihalf_edge0 with ihalf_edge_new0.
    */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    void ReplaceTriangleWithTriangleEdge
    (const ITYPE2 ihalf_edge0, const ITYPE3 iv_new,
     ITYPE4 & ihalf_edge_new0);

    /// \brief Add vertex and replace triangle with triangle edge.
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    void AddVertexReplaceTriangleWithTriangleEdge
    (const ITYPE2 ihalf_edge0, ITYPE3 & iv_new,
     ITYPE4 & ihalf_edge_new0);


    /**
       @brief Replace split edge triangle with triangle edge.
         - Triangle has two split edges.
       @param ihalf_edge0 Base of triangle.
       @param iv_new New triangle vertex.
    */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    void ReplaceSplit2EdgeTriangleWithTriangleEdge
    (const ITYPE2 ihalf_edge0, const ITYPE3 iv_new,
     ITYPE4 & ihalf_edge_new0);

    /**
       @brief Add vertex and replace split edge triangle with triangle edge.
         - Triangle has two split edges.
    */
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    void AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
    (const ITYPE2 ihalf_edge0, ITYPE3 & iv_new,
     ITYPE4 & ihalf_edge_new0);

  };


  // *****************************************************************
  // Class MESH2D_SPLIT_A
  // *****************************************************************

  /**
     @brief Mesh of 2D polygons. 
            (Uses defaults for vertex, half_edge and polygon types.)

     @details
     - Use VERTEX_BASE, HALF_EDGE_BASE, POLYGON_BASE
       to represent vertices, half edges and polygons.
  */
  template <typename ITYPE, typename NTYPE>
  class MESH2D_SPLIT_A:
    public MESH2D_SPLIT_BASE<VERTEX_BASE<ITYPE,ITYPE,NTYPE>, 
                             HALF_EDGE_BASE<ITYPE,ITYPE,ITYPE>,
                             POLYGON_BASE,ITYPE,NTYPE>
  {
  public:
    MESH2D_SPLIT_A() {};
  };


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE member functions
  // *****************************************************************

  // constructor
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  void MESH2D_SPLIT_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  Init()
  {
    // Nothing to do, so far.
  }

  
  // Copy.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename MTYPE>
  void MESH2D_SPLIT_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  Copy(const MTYPE & meshA)
  {
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::Copy(meshA);
  }


  // Location of j'th vertex in polygon created 
  //   by splitting half edge ihalf_edge0.
  // - ToVertex(ihalf_edgeA) is 0'th vertex of new polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3>
  typename MESH2D_SPLIT_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::INDEX_TYPE 
  MESH2D_SPLIT_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
  LocationInSplitEdgePolygon
  (const ITYPE2 ihalf_edgeA, const ITYPE3 j) const
  {
    const POLYGON_INDEX_TYPE ipoly = 
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NUMBER_TYPE numv = this->NumPolygonVertices(ipoly);

    if (j == numv) {
      // j is vertex splitting ihalf_edgeA.
      return(numv);
    }
    else {
      // jloc0 = Location of first vertex in polygon with split half edge.
      const INDEX_TYPE jloc0 = 
        this->LocationOfNextHalfEdgeInPolygon(ihalf_edgeA);

      // jloc_new = Location of j'th index in polygon with split half edge.
      const INDEX_TYPE jloc_new = 
        this->LocationRelativeToLocation(ipoly, j, jloc0);

      return(jloc_new);
    }
  }


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE Replace vertex member functions
  // *****************************************************************

  // Replace from vertex.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  ReplaceFromVertex
  (const ITYPE2 ihalf_edge, const ITYPE3 iv_new)
  {
    const VERTEX_INDEX_TYPE iv_old = this->FromVertexIndex(ihalf_edge);

    this->_RemoveHalfEdgeIncidentOnFromVertex(ihalf_edge);
    this->_SetHalfEdgeFromVertex(ihalf_edge, iv_new);
  }


  // Replace to vertex.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  ReplaceToVertex
  (const ITYPE2 ihalf_edge, const ITYPE3 iv_new)
  {
    ReplaceFromVertex
      (this->IndexOfNextHalfEdgeInPolygon(ihalf_edge), iv_new);
  }


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE Split member functions
  // *****************************************************************

  // Split polygon into two polygons.
  // - Split polygon with diagonal 
  //   (FromVertexIndex(ihalf_edge0), FromVertexIndex(ihalf_edge1)).
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitPolygonByDiagonal
  (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1,
   ITYPE4 & ihalf_edge_new0, ITYPE5 & ihalf_edge_new1)
  {
    const POLYGON_INDEX_TYPE ipolyA = 
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    const NTYPE numv0 = this->NumPolygonVertices(ipolyA);
    const VERTEX_INDEX_TYPE iv0 = this->FromVertexIndex(ihalf_edge0);
    const VERTEX_INDEX_TYPE iv1 = this->FromVertexIndex(ihalf_edge1);
    const ITYPE iloc0 = this->LocationOfHalfEdgeInPolygon(ihalf_edge0);
    const ITYPE iloc1 = this->LocationOfHalfEdgeInPolygon(ihalf_edge1);
    NTYPE numvB0, numvB1;
    IJK::PROCEDURE_ERROR error("MESH2D_SPLIT_BASE::SplitPolygonByDiagonal");

    if (!this->CheckAreHalfEdgesInSamePolygon
        (ihalf_edge0, ihalf_edge1, error)) { throw error; };

    if (!this->CheckAreHalfEdgesDifferent
        (ihalf_edge0, ihalf_edge1, error)) { throw error; };

    if (iloc1 > iloc0) {
      numvB0 = (iloc1-iloc0)+1;
      numvB1 = (numv0-iloc1) + iloc0 + 1;
    }
    else {
      numvB0 = (numv0-iloc0) + iloc1 + 1;
      numvB1 = (iloc0-iloc1)+1;
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "  Splitting at vertices: "
         << this->PolygonVertex(ipolyA, iloc0)
         << " and " << this->PolygonVertex(ipolyA, iloc1) << endl;
    this->PrintHalfEdge(cerr, "  ihalf_edge0: ", ihalf_edge0, "\n");
    this->PrintHalfEdge(cerr, "  ihalf_edge1: ", ihalf_edge1, "\n");
    this->PrintPolygonVertices
      (cerr, "  Polygon vert: ", ipolyA, "\n");
    cerr << "  iv0: " << iv0 << endl;
    cerr << "  iv1: " << iv1 << endl;
    */

    // Add two new polygons.
    POLYGON_INDEX_TYPE ipolyB0, ipolyB1;
    this->_AddTwoPolygonsX(numvB0, numvB1, ipolyB0, ipolyB1);

    ihalf_edge_new0 =
      this->polygon_half_edge_list.first_element[ipolyB0];
    this->_ReplaceHalfEdges(ihalf_edge0, ihalf_edge_new0, numvB0-1);

    ihalf_edge_new1 =
      this->polygon_half_edge_list.first_element[ipolyB1];
    this->_ReplaceHalfEdges(ihalf_edge1, ihalf_edge_new1, numvB1-1);

    const HALF_EDGE_INDEX_TYPE ihalf_edge_diagonal_B0 =
      ihalf_edge_new0 + numvB0-1;
    const HALF_EDGE_INDEX_TYPE ihalf_edge_diagonal_B1 =
      ihalf_edge_new1 + numvB1-1;

    this->_SetHalfEdgeFromVertex(ihalf_edge_diagonal_B0, iv1);
    this->_SetHalfEdgeFromVertex(ihalf_edge_diagonal_B1, iv0);

    this->_LinkTwoHalfEdgesAroundEdge
      (ihalf_edge_diagonal_B0, ihalf_edge_diagonal_B1);

    // *** DEBUG ***
    /*
    this->PrintPolygonVertices
      (cerr, "  Polygon B0 vert: ", ipolyB0, "\n");
    this->PrintPolygonVertices
      (cerr, "  Polygon B1 vert: ", ipolyB1, "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  ihalf_edge_diagonal_B0: ", ihalf_edge_diagonal_B0, "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  ihalf_edge_diagonal_B1: ", ihalf_edge_diagonal_B1, "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  ihalf_edge: ",
       this->polygon_half_edge_list.first_element[ipolyB0],
       "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  ihalf_edge: ",
       this->polygon_half_edge_list.first_element[ipolyB0]+1,
       "\n");
    */

    this->_DeletePolygon(ipolyA);
  }


  // Split half edge ihalf_edge of polygon, creating new polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  typename MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  POLYGON_INDEX_TYPE
  MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitPolygonHalfEdge
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split,
   ITYPE4 & ihalf_edge_new)
  {
    const POLYGON_INDEX_TYPE ipoly =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NTYPE nume = this->NumPolygonEdges(ipoly);
    const NTYPE nume_new = nume+1;

    const POLYGON_INDEX_TYPE ipoly_new = this->_AddPolygonX(nume_new);

    const HALF_EDGE_INDEX_TYPE ihalf_edge_newX =
      this->polygon_half_edge_list.first_element[ipoly_new];

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  Replacing half edges starting at: ", ihalf_edgeB, "\n");
    */

    const NUMBER_TYPE ilocA = 
      this->LocationOfHalfEdgeInPolygon(ihalf_edgeA);

    if (ilocA > 0) {
      this->_ReplaceHalfEdges
        (this->FirstHalfEdgeIndex(ipoly), 
         this->FirstHalfEdgeIndex(ipoly_new), ilocA);
    }

    if (ilocA+1 < nume) {
      this->_ReplaceHalfEdges
        (this->HalfEdgeIndex(ipoly, ilocA+1), 
         this->HalfEdgeIndex(ipoly_new, ilocA+2), (nume-ilocA-1));
    }

    ihalf_edge_new = this->HalfEdgeIndex(ipoly_new, ilocA);
    this->_ReplaceHalfEdgeWithTwoHalfEdges
      (ihalf_edgeA, iv_split, ihalf_edge_new);

    this->_DeletePolygon(ipoly);

    return(ipoly_new);
  }


  // Split half edge ihalf_edge of polygon, creating new polygon,
  //   and link two new half edges to ihalf_edgeXC and ihalf_edgeXD.
  // @param ihalf_edgeXC First half edge in opposite direction.
  // @param ihalf_edgeXD Second half edge in opposite direction.
  //   - ToVertexIndex(ihalf_edgeXC) = iv_split.
  //   - FromVertexIndex(ihalf_edgeXD) = iv_split.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4,
            typename ITYPE5, typename ITYPE6>
  typename MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  POLYGON_INDEX_TYPE
  MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitAndLinkPolygonHalfEdge
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split,
   const ITYPE4 ihalf_edgeXC, const ITYPE5 ihalf_edgeXD,
   ITYPE6 & ihalf_edge_newA)
  {
    const POLYGON_INDEX_TYPE ipoly_new =
      SplitPolygonHalfEdge(ihalf_edgeA, iv_split, ihalf_edge_newA);

    const HALF_EDGE_INDEX_TYPE ihalf_edge_newB =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_newA);

    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeA);

    this->_LinkTwoHalfEdgesAroundEdge(ihalf_edge_newA, ihalf_edgeXD);
    this->_LinkTwoHalfEdgesAroundEdge(ihalf_edge_newB, ihalf_edgeXC);

    return(ipoly_new);
  }


  // Split edge contained in two polygon.
  // - If edge in only one polygon, split one polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5,
            typename ITYPE6>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitPolygonEdgeII
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split,
   ITYPE4 & ipoly0, ITYPE5 & ipoly1,
   ITYPE6 & ihalf_edgeA_new)
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edgeAX =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);

    HALF_EDGE_INDEX_TYPE ihalf_edgeA_new_next;
    HALF_EDGE_INDEX_TYPE ihalf_edgeAX_new;

    if (this->IsBoundaryEdge(ihalf_edgeA)) {
      ipoly0 = SplitPolygonHalfEdge
        (ihalf_edgeA, iv_split, ihalf_edgeA_new);
      this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly0);
      ipoly1 = ipoly0;
    }
    else {
      ipoly0 = SplitPolygonHalfEdge
        (ihalf_edgeA, iv_split, ihalf_edgeA_new);
      const HALF_EDGE_INDEX_TYPE ihalf_edgeA_new_next = 
        this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeA_new);
      ipoly1 = SplitAndLinkPolygonHalfEdge
        (ihalf_edgeAX, iv_split, ihalf_edgeA_new, ihalf_edgeA_new_next, 
         ihalf_edgeAX_new);

      // Call UpdateHalfEdges*() after both old polygons are deleted
      //   in SplitPolygonHalfEdge().
      this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly0);
      this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly1);
    }
  }


  // Split edge contained in two polygon.
  // - Version that does not return new half edge.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitPolygonEdgeII
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split,
   ITYPE4 & ipoly0, ITYPE5 & ipoly1)
  {
    HALF_EDGE_INDEX_TYPE ihalf_edgeA_new;

    SplitPolygonEdgeII
      (ihalf_edgeA, iv_split, ipoly0, ipoly1, ihalf_edgeA_new);
  }


  // Add vertex and split edge contained in two polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5, typename ITYPE6>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexSplitPolygonEdgeII
  (const ITYPE2 ihalf_edgeA, ITYPE3 & iv_split,
   ITYPE4 & ipoly0, ITYPE5 & ipoly1, ITYPE6 & ihalf_edgeA_new)
  {
    iv_split = this->AddVertex();
    SplitPolygonEdgeII
      (ihalf_edgeA, iv_split, ipoly0, ipoly1, ihalf_edgeA_new);
  }


  // Add vertex and split edge contained in two polygon.
  // - Version that does not return ihalf_edgeA_new.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexSplitPolygonEdgeII
  (const ITYPE2 ihalf_edgeA, ITYPE3 & iv_split,
   ITYPE4 & ipoly0, ITYPE5 & ipoly1)
  {
    HALF_EDGE_INDEX_TYPE ihalf_edgeA_new;

    iv_split = this->AddVertex();
    SplitPolygonEdgeII
      (ihalf_edgeA, iv_split, ipoly0, ipoly1, ihalf_edgeA_new);
  }


  // Add two vertices and split two polygon edges.
  // - Splits half edges in 3 polygons.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEH, typename ITYPEV, typename ITYPEH2>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddTwoVerticesSplitTwoEdgesOfPolygon
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     ITYPEV & ivA_split, ITYPEV & ivB_split,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new)
  {
    ivA_split = this->AddVertex();
    ivB_split = this->AddVertex();
    SplitTwoEdgesOfPolygon
      (ihalf_edgeA, ihalf_edgeB, ivA_split, ivB_split,
       ihalf_edgeA_new, ihalf_edgeB_new);
  }


  // Add three vertices and split three polygon edges.
  // - Splits half edges in 4 polygons.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEH, typename ITYPEV, typename ITYPEH2>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddThreeVerticesSplitThreeEdgesOfPolygon
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     const ITYPEH ihalf_edgeC,
     ITYPEV & ivA_split, ITYPEV & ivB_split,
     ITYPEV & ivC_split,
     ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new,
     ITYPEH2 & ihalf_edgeC_new)
  {
    ivA_split = this->AddVertex();
    ivB_split = this->AddVertex();
    ivC_split = this->AddVertex();
    SplitThreeEdgesOfPolygon
      (ihalf_edgeA, ihalf_edgeB, ihalf_edgeC,
       ivA_split, ivB_split, ivC_split,
       ihalf_edgeA_new, ihalf_edgeB_new, ihalf_edgeC_new);
  }


  // Split two edges of a polygon.
  // - Splits half edges in three polygons.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEH, typename ITYPEV, typename ITYPEH2>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitTwoEdgesOfPolygon
  (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
   const ITYPEV ivA_split, const ITYPEV ivB_split,
   ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new)
  {
    const NUMBER_TYPE TWO(2);
    const VERTEX_INDEX_TYPE iv_split[TWO] = { ivA_split, ivB_split };
    const HALF_EDGE_INDEX_TYPE ihalf_edge[TWO] = { ihalf_edgeA, ihalf_edgeB };
    HALF_EDGE_INDEX_TYPE ihalf_edge_new[TWO];

    SplitEdgesOfPolygonNew(ihalf_edge, iv_split, TWO, ihalf_edge_new);

    ihalf_edgeA_new = ihalf_edge_new[0];
    ihalf_edgeB_new = ihalf_edge_new[1];
  }


  // Split three edges of a polygon.
  // - Splits half edges in three polygons.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEH, typename ITYPEV, typename ITYPEH2>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitThreeEdgesOfPolygon
  (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
   const ITYPEH ihalf_edgeC,
   const ITYPEV ivA_split, const ITYPEV ivB_split,
   const ITYPEV ivC_split,
   ITYPEH2 & ihalf_edgeA_new, ITYPEH2 & ihalf_edgeB_new,
   ITYPEH2 & ihalf_edgeC_new)
  {
    const NUMBER_TYPE THREE(3);
    const VERTEX_INDEX_TYPE iv_split[THREE] = 
      { ivA_split, ivB_split, ivC_split };
    const HALF_EDGE_INDEX_TYPE ihalf_edge[THREE] = 
      { ihalf_edgeA, ihalf_edgeB, ihalf_edgeC };
    HALF_EDGE_INDEX_TYPE ihalf_edge_new[THREE];

    SplitEdgesOfPolygonNew(ihalf_edge, iv_split, THREE, ihalf_edge_new);

    ihalf_edgeA_new = ihalf_edge_new[0];
    ihalf_edgeB_new = ihalf_edge_new[1];
    ihalf_edgeC_new = ihalf_edge_new[2];
  }


  // Split three edges of a polygon.
  // - Splits half edges in four polygons.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEH, typename ITYPEV, typename ITYPEP>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitThreePolygonEdgesIV
  (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
   const ITYPEH ihalf_edgeC,
   const ITYPEV ivA_split, const ITYPEV ivB_split,
   const ITYPEV ivC_split,
   ITYPEP & ipolyABC, ITYPEP & ipolyA, ITYPEP & ipolyB,
   ITYPEP & ipolyC)
  {
    const NUMBER_TYPE THREE(3);
    const NUMBER_TYPE FOUR(4);
    const VERTEX_INDEX_TYPE iv_split[THREE] = 
      { ivA_split, ivB_split, ivC_split };
    const HALF_EDGE_INDEX_TYPE ihalf_edge[THREE] = 
      { ihalf_edgeA, ihalf_edgeB, ihalf_edgeC };
    POLYGON_INDEX_TYPE ipoly_new[FOUR];

    // TO BE CONTINUED...

    /* MODIFY
    SplitEdgesOfPolygon(ihalf_edge, iv_split, THREE, ipoly_new);

    const NUMBER_TYPE ilocA = 
      this->LocationOfHalfEdgeInPolygon(ihalf_edgeA);
    const NUMBER_TYPE ilocB = 
      this->LocationOfHalfEdgeInPolygon(ihalf_edgeB);

    ipolyAB = ipoly_new[0];

    if (ilocA < ilocB) {
      ipolyA = ipoly_new[1];
      ipolyB = ipoly_new[2];
    }
    else {
      ipolyA = ipoly_new[2];
      ipolyB = ipoly_new[1];
    }
    */
  }


  // Split polygon half edges at given locations.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEP, typename ITYPEL,
            typename ITYPEV, typename NTYPE2>
  typename MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  POLYGON_INDEX_TYPE
  MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitPolygonHalfEdgesAtLoc
  (const ITYPEP ipoly0, const ITYPEL * isplit_loc, 
   const ITYPEV * iv_split, const NTYPE2 num_split)
  {
    const NUMBER_TYPE THREE(3);
    const NTYPE nume0 = this->NumPolygonEdges(ipoly0);
    const NTYPE nume1 = nume0+num_split;
    IJK::PROCEDURE_ERROR error("MESH2D_SPLIT_BASE::SplitPolygonHalfEdges");

    if (num_split < 1) {
      error.AddMessage("Programming error. No splits specified.");
      throw error;
    }

    if (nume1 < THREE) {
      error.AddMessage
        ("Programming error. Polygon has fewer than ", THREE, " edges.");
      error.AddMessage("  Polygon ", ipoly0, " has only ", nume0, " edges.");
      throw error;
    }

    const POLYGON_INDEX_TYPE ipoly1 = this->_AddPolygonX(nume1);

    NUMBER_TYPE j = 0;
    NUMBER_TYPE i1 = 0;
    for (NUMBER_TYPE i0 = 0; i0 < nume0; i0++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge0 = 
        this->HalfEdgeIndex(ipoly0, i0);
      const VERTEX_INDEX_TYPE iv0 = this->FromVertexIndex(ihalf_edge0);
      const HALF_EDGE_INDEX_TYPE ihalf_edge1 = 
        this->HalfEdgeIndex(ipoly1, i1);
      this->_SetHalfEdgeFromVertex(ihalf_edge1, iv0);

      if ((j < num_split) &&
          (isplit_loc[j] == i0)) {
        const HALF_EDGE_INDEX_TYPE ihalf_edge1B =
          this->IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
        this->_ReplaceHalfEdgeWithTwoHalfEdges
          (ihalf_edge0, iv_split[j], ihalf_edge1);
        this->_SetIndexOfHalfEdgeIncidentOnVertex
          (iv_split[j], ihalf_edge1B);
        j++;
        i1 += 2;
      }
      else {
        this->_ReplaceHalfEdge(ihalf_edge0, ihalf_edge1);
        i1++;
      }
    }

    if (j != num_split) {
      error.AddMessage("Programming error. Failed in creating new polygon.");
      error.AddMessage("  Only split ", j, " edges.");
      error.AddMessage("  Should split ", num_split, " edges.");
      throw error;
    }

    if (i1 != nume1) {
      error.AddMessage("Programming error. Failed in creating new polygon.");
      error.AddMessage
        ("  Only created ", i1, " edges in new polygon ", ipoly1, ".");
      error.AddMessage("  New polygon should have ", nume1, " edges.");
      throw error;
    }

    this->_DeletePolygon(ipoly0);

    return(ipoly1);
  }


  // Split polygon half edges at given locations.
  // - Version using C++ STL vector isplit_loc[] and iv_split[].
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEP, typename ITYPEL, typename ITYPEV>
  typename MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  POLYGON_INDEX_TYPE
  MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitPolygonHalfEdgesAtLoc
  (const ITYPEP ipoly0, const std::vector<ITYPEL> & isplit_loc,
   const std::vector<ITYPEV> & iv_split)
  {
    const POLYGON_INDEX_TYPE ipoly_new = SplitPolygonHalfEdgesAtLoc
      (ipoly0, IJK::vector2pointer(isplit_loc),
       IJK::vector2pointer(iv_split), isplit_loc.size());

    return(ipoly_new);
  }


  // Split two half edges of a polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEHA, typename ITYPEHB, 
            typename ITYPEVA, typename ITYPEVB,
            typename ITYPEHA2, typename ITYPEHB2>
  typename MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  POLYGON_INDEX_TYPE
  MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitTwoHalfEdges
  (const ITYPEHA ihalf_edgeA, const ITYPEHB ihalf_edgeB,
   const ITYPEVA iv_splitA, const ITYPEVB iv_splitB,
   ITYPEHA2 & ihalf_edgeA_new, ITYPEHB2 & ihalf_edgeB_new)
  {
    const NUMBER_TYPE TWO(2);
    const POLYGON_INDEX_TYPE ipoly = 
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const NTYPE nume = this->NumPolygonEdges(ipoly);

    NUMBER_TYPE isplit_loc[TWO];
    ITYPEVA iv_split[TWO];
    IJK::PROCEDURE_ERROR error("MESH2D_SPLIT_BASE::SplitTwoHalfEdges");

    // If true, (flocation of ihalf_edgeA) < (location of ihalf_edgeB)
    bool flag_orderAB;     

    if (!this->CheckAreHalfEdgesInSamePolygon
        (ihalf_edgeA, ihalf_edgeB, error)) 
      { throw error; }

    isplit_loc[0] = this->LocationOfHalfEdgeInPolygon(ihalf_edgeA);
    isplit_loc[1] = this->LocationOfHalfEdgeInPolygon(ihalf_edgeB);

    if (isplit_loc[0] < isplit_loc[1]) {
      flag_orderAB = true;
      iv_split[0] = iv_splitA;
      iv_split[1] = iv_splitB;
    }
    else {
      flag_orderAB = false;
      std::swap(isplit_loc[0], isplit_loc[1]);
      iv_split[0] = iv_splitB;
      iv_split[1] = iv_splitA;
    }

    const POLYGON_INDEX_TYPE ipoly_new =
      SplitPolygonHalfEdgesAtLoc(ipoly, isplit_loc, iv_split, TWO);

    if (flag_orderAB) {
      ihalf_edgeA_new = this->HalfEdgeIndex(ipoly_new, isplit_loc[0]);
      ihalf_edgeB_new = this->HalfEdgeIndex(ipoly_new, isplit_loc[1]+1);
    }
    else {
      ihalf_edgeB_new = this->HalfEdgeIndex(ipoly_new, isplit_loc[0]);
      ihalf_edgeA_new = this->HalfEdgeIndex(ipoly_new, isplit_loc[1]+1);
    }

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "Splitting half edges ", ihalf_edgeA, "");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, " and ", ihalf_edgeB, "");
    cerr << "  flag_orderAB: " << int(flag_orderAB) << endl;
    this->PrintPolygonIndexAndVertices(cerr, "  New poly: ", ipoly_new, "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "    ihalf_edgeA_new: ", ihalf_edgeA_new, "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "    ihalf_edgeB_new: ", ihalf_edgeB_new, "\n");
    */

    return(ipoly_new);
  }


  // Split edges of one polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEH, typename ITYPEV, typename ITYPEP,
            typename NTYPES>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitEdgesOfPolygon
  (const ITYPEH ihalf_edge[], const ITYPEV iv_split[],
   const NTYPES num_split, ITYPEP ipoly_new[])
  {
    IJK::PROCEDURE_ERROR error("MESH2D_SPLIT_BASE::SplitEdgesOfPolygon");

    if (num_split < 1) {
      error.AddMessage("Programming error. No splits specified.");
      throw error;
    }

    const POLYGON_INDEX_TYPE ipoly0 = 
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge[0]);
    const NTYPE nume0 = this->NumPolygonEdges(ipoly0);

    std::vector< std::pair<NUMBER_TYPE,VERTEX_INDEX_TYPE> >
      split_info(num_split);
    std::vector<NUMBER_TYPE> isplit_loc(num_split);
    std::vector<VERTEX_INDEX_TYPE> iv_splitII(num_split);
    std::vector<HALF_EDGE_INDEX_TYPE> ihalf_edge_newX(num_split);
    std::vector<bool> is_boundary_edge(num_split);

    for (NUMBER_TYPE j = 0; j < num_split; j++) {
      split_info[j].first = 
        this->LocationOfHalfEdgeInPolygon(ihalf_edge[j]);
      split_info[j].second = iv_split[j];
    }

    // Reorder split locations, so that splits are in increasing order
    //   by location.
    std::sort(split_info.begin(), split_info.end());

    for (NUMBER_TYPE j = 0; j < num_split; j++) {
      isplit_loc[j] = split_info[j].first;
      iv_splitII[j] = split_info[j].second;
    }

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices
      (cerr, "***Splitting polygon ", ipoly0, "\n");
    IJK::print_list(cerr, "  isplit_loc: ", isplit_loc, "\n");
    IJK::print_list(cerr, "  iv_splitII: ", iv_splitII, "\n");
    */

    for (NUMBER_TYPE j = 0; j < num_split; j++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge_j =
        this->HalfEdgeIndex(ipoly0, isplit_loc[j]);

      // Store flag indicating which split edges are boundary edges.
      is_boundary_edge[j] = this->IsBoundaryEdge(ihalf_edge_j);

      if (!is_boundary_edge[j]) {

        const HALF_EDGE_INDEX_TYPE ihalf_edgeX =
          this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge_j);

        // *** DEBUG ***
        /*
        using namespace std;
        this->PrintHalfEdgeIndexAndEndpoints
          (cerr, "  Splitting half edge: ", ihalf_edgeX, "\n");
        this->PrintPolygonIndexAndVertices
          (cerr, "  Splitting 1 edge of polygon ",
           this->IndexOfPolygonContainingHalfEdge(ihalf_edgeX), "\n");
        */

        ipoly_new[j+1] = SplitPolygonHalfEdge
          (ihalf_edgeX, iv_splitII[j], ihalf_edge_newX[j]);
      }
    }

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices
      (cerr, "  Splitting multiple edges of polygon ", ipoly0, "\n");
    */

    const POLYGON_INDEX_TYPE ipoly1 =
      SplitPolygonHalfEdgesAtLoc(ipoly0, isplit_loc, iv_splitII);

    ipoly_new[0] = ipoly1;

    NUMBER_TYPE i1 = 0;
    NUMBER_TYPE j = 0;
    for (NUMBER_TYPE i0 = 0; i0 < nume0; i0++) {
      if ((j < num_split) && (i0 == isplit_loc[j])) {

        if (is_boundary_edge[j]) {
          // No new polygon j+1.
          ipoly_new[j+1] = ipoly1;
        }
        else {
          // Link ipoly1 with ipoly_new[j+1].
          const HALF_EDGE_INDEX_TYPE ihalf_edge_new1A =
            this->HalfEdgeIndex(ipoly1, i1);
          const HALF_EDGE_INDEX_TYPE ihalf_edge_new1B =
            this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new1A);
          const HALF_EDGE_INDEX_TYPE ihalf_edge_newXA =
            ihalf_edge_newX[j];
          const HALF_EDGE_INDEX_TYPE ihalf_edge_newXB =
            this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_newXA);

          this->_LinkTwoHalfEdgesAroundEdge
            (ihalf_edge_newXA, ihalf_edge_new1B);
          this->_LinkTwoHalfEdgesAroundEdge
            (ihalf_edge_newXB, ihalf_edge_new1A);
        }
        j++;
        i1 += 2;
      }
      else {
        i1++;
      }
        
    }

    if (j != num_split) {
      error.AddMessage("Programming error. Failed linking split edges.");
      error.AddMessage("  Only processed ", j, " split edges.");
      throw error;
    }

    if (i1 != this->NumPolygonEdges(ipoly1)) {
      error.AddMessage("Programming error. Failed linking split edges.");
      error.AddMessage
        ("  Did not process all edges of new polygon ", ipoly1, ".");
      throw error;
    }
  }


  // Split edges of one polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEH, typename ITYPEV,
            typename NTYPES, typename ITYPEH2>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitEdgesOfPolygonNew
  (const ITYPEH ihalf_edge[], const ITYPEV iv_split[],
   const NTYPES num_split, ITYPEH2 ihalf_edge_new[])
  {
    IJK::PROCEDURE_ERROR error("MESH2D_SPLIT_BASE::SplitEdgesOfPolygon");

    if (num_split < 1) {
      error.AddMessage("Programming error. No splits specified.");
      throw error;
    }

    const POLYGON_INDEX_TYPE ipoly0 = 
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge[0]);
    const NTYPE nume0 = this->NumPolygonEdges(ipoly0);

    std::vector< std::pair<NUMBER_TYPE,VERTEX_INDEX_TYPE> >
      split_info(num_split);
    std::vector<NUMBER_TYPE> isplit_loc(num_split);
    std::vector<VERTEX_INDEX_TYPE> iv_splitII(num_split);
    std::vector<HALF_EDGE_INDEX_TYPE> ihalf_edge_newX(num_split);
    std::vector<bool> is_boundary_edge(num_split);

    for (NUMBER_TYPE j = 0; j < num_split; j++) {
      split_info[j].first = 
        this->LocationOfHalfEdgeInPolygon(ihalf_edge[j]);
      split_info[j].second = j;
    }

    // Reorder split locations, so that splits are in increasing order
    //   by location.
    std::sort(split_info.begin(), split_info.end());

    for (NUMBER_TYPE j = 0; j < num_split; j++) {
      isplit_loc[j] = split_info[j].first;
      iv_splitII[j] = iv_split[split_info[j].second];
    }

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices
      (cerr, "***Splitting polygon ", ipoly0, "\n");
    IJK::print_list(cerr, "  isplit_loc: ", isplit_loc, "\n");
    IJK::print_list(cerr, "  iv_splitII: ", iv_splitII, "\n");
    */

    for (NUMBER_TYPE j = 0; j < num_split; j++) {
      const HALF_EDGE_INDEX_TYPE ihalf_edge_j =
        this->HalfEdgeIndex(ipoly0, isplit_loc[j]);

      // Store flag indicating which split edges are boundary edges.
      is_boundary_edge[j] = this->IsBoundaryEdge(ihalf_edge_j);

      if (!is_boundary_edge[j]) {
        const HALF_EDGE_INDEX_TYPE ihalf_edgeX =
          this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge_j);

        // *** DEBUG ***
        /*
        using namespace std;
        this->PrintHalfEdgeIndexAndEndpoints
          (cerr, "  Splitting half edge: ", ihalf_edgeX, "\n");
        this->PrintPolygonIndexAndVertices
          (cerr, "  Splitting 1 edge of polygon ",
           this->IndexOfPolygonContainingHalfEdge(ihalf_edgeX), "\n");
        */

        SplitPolygonHalfEdge
          (ihalf_edgeX, iv_splitII[j], ihalf_edge_newX[j]);
      }
    }

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices
      (cerr, "  Splitting multiple edges of polygon ", ipoly0, "\n");
    */

    const POLYGON_INDEX_TYPE ipoly0_new =
      SplitPolygonHalfEdgesAtLoc(ipoly0, isplit_loc, iv_splitII);

    NUMBER_TYPE i1 = 0;
    NUMBER_TYPE j = 0;
    for (NUMBER_TYPE i0 = 0; i0 < nume0; i0++) {
      if ((j < num_split) && (i0 == isplit_loc[j])) {

        const HALF_EDGE_INDEX_TYPE ihalf_edge_new1A =
          this->HalfEdgeIndex(ipoly0_new, i1);

        const NUMBER_TYPE j2 = split_info[j].second;
        ihalf_edge_new[j2] = ihalf_edge_new1A;

        if (!is_boundary_edge[j]) {
          // Link ipoly1 with ipoly_new[j+1].
          const HALF_EDGE_INDEX_TYPE ihalf_edge_new1B =
            this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new1A);
          const HALF_EDGE_INDEX_TYPE ihalf_edge_newXA =
            ihalf_edge_newX[j];
          const HALF_EDGE_INDEX_TYPE ihalf_edge_newXB =
            this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_newXA);

          this->_LinkTwoHalfEdgesAroundEdge
            (ihalf_edge_newXA, ihalf_edge_new1B);
          this->_LinkTwoHalfEdgesAroundEdge
            (ihalf_edge_newXB, ihalf_edge_new1A);
        }

        j++;
        i1 += 2;
      }
      else {
        i1++;
      }
        
    }

    if (j != num_split) {
      error.AddMessage("Programming error. Failed linking split edges.");
      error.AddMessage("  Only processed ", j, " split edges.");
      throw error;
    }

    if (i1 != this->NumPolygonEdges(ipoly0_new)) {
      error.AddMessage("Programming error. Failed linking split edges.");
      error.AddMessage
        ("  Did not process all edges of new polygon ", ipoly0_new, ".");
      throw error;
    }
  }


  // Add vertices and split edges and polygons.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPEH, typename ITYPEP, typename NTYPES>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVerticesSplitEdgesOfPolygon
  (const ITYPEH ihalf_edge[], const NTYPES num_split,
   ITYPEP ipoly_new[])
  {
    IJK::PROCEDURE_ERROR error
      ("MESH2D_SPLIT_BASE::AddVerticesSplitEdgesOfPolygon");

    if (num_split < 1) {
      error.AddMessage("Programming error. No splits specified.");
      throw error;
    }

    std::vector<VERTEX_INDEX_TYPE> iv_split(num_split);

    for (NUMBER_TYPE j = 0; j < num_split; j++) 
      { iv_split[j] = this->AddVertex();  }

    SplitEdgesOfPolygon
      (ihalf_edge, IJK::vector2pointer(iv_split), num_split, ipoly_new);
  }


  // Split three edges incident on the same vertex.
  // *** SHOULD MODIFY TO HANDLE CASE WHERE POLY 0 == POLY 4
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, 
            typename ITYPEA, typename ITYPEB, typename ITYPEC,
            typename ITYPEP>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitThreeIncidentEdgesIV
  (const ITYPE2 ihalf_edgeB1, 
   const ITYPEA iv_splitA, const ITYPEB iv_splitB, const ITYPEC iv_splitC,
   ITYPEP ipoly[4])
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC1 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeB1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC0 =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeC1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB2 =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeB1);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeA2 =
      this->IndexOfPrevHalfEdgeInPolygon(ihalf_edgeB2);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeA3 =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA2);
    const bool is_boundary_edgeA = this->IsBoundaryEdge(ihalf_edgeA2);
    const bool is_boundary_edgeC = this->IsBoundaryEdge(ihalf_edgeC1);
    IJK::PROCEDURE_ERROR error("MESH2D_SPLIT_BASE::SplitThreeIncidentEdgesIV");

    HALF_EDGE_INDEX_TYPE ihalf_edgeB1_new, ihalf_edgeB1_new2;
    HALF_EDGE_INDEX_TYPE ihalf_edgeC1_new, ihalf_edgeC1_new2;
    HALF_EDGE_INDEX_TYPE ihalf_edgeA2_new, ihalf_edgeA2_new2;
    HALF_EDGE_INDEX_TYPE ihalf_edgeB2_new, ihalf_edgeB2_new2;

    if (!this->IsBoundaryEdge(ihalf_edgeC1) && 
        !this->IsBoundaryEdge(ihalf_edgeB2)) {
      if (this->IndexOfPolygonContainingHalfEdge(ihalf_edgeC0) ==
          this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA3)) {
        error.AddMessage("Programming error. Polygon 0 and 3 are identical.");
        error.AddMessage
          ("  Splitting three edges incident on vertex ",
           this->ToVertexIndex(ihalf_edgeB1), ".\n");
        throw error;
      }
    }

    if (!this->IsBoundaryEdge(ihalf_edgeC1)) {
      if (this->IndexOfPolygonContainingHalfEdge(ihalf_edgeC0) ==
          this->IndexOfPolygonContainingHalfEdge(ihalf_edgeB2)) {
        error.AddMessage("Programming error. Polygon 0 and 2 are identical.");
        error.AddMessage
          ("  Splitting three edges incident on vertex ",
           this->ToVertexIndex(ihalf_edgeB1), ".\n");
        throw error;
      }
    }

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  iv_splitA: " << iv_splitA << endl;
    cerr << "  iv_splitB: " << iv_splitB << endl;
    cerr << "  iv_splitC: " << iv_splitC << endl;
    */

    ipoly[1] = 
      SplitTwoHalfEdges
      (ihalf_edgeB1, ihalf_edgeC1, iv_splitB, iv_splitC, 
       ihalf_edgeB1_new, ihalf_edgeC1_new);

    ipoly[2] = 
      SplitTwoHalfEdges
      (ihalf_edgeA2, ihalf_edgeB2, iv_splitA, iv_splitB, 
       ihalf_edgeA2_new, ihalf_edgeB2_new);

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices(cerr, "New poly 1: ", ipoly[1], "\n");
    this->PrintPolygonIndexAndVertices(cerr, "New poly 2: ", ipoly[2], "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  ihalf_edgeB1_new: ", ihalf_edgeB1_new, "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  ihalf_edgeC1_new: ", ihalf_edgeC1_new, "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  ihalf_edgeA2_new: ", ihalf_edgeA2_new, "\n");
    this->PrintHalfEdgeIndexAndEndpoints
      (cerr, "  ihalf_edgeB2_new: ", ihalf_edgeB2_new, "\n");
    */

    ihalf_edgeB1_new2 = this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeB1_new);
    ihalf_edgeC1_new2 = this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeC1_new);
    ihalf_edgeA2_new2 = this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeA2_new);
    ihalf_edgeB2_new2 = this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeB2_new);

    // Link poly[1] and poly[2]
    this->_LinkTwoHalfEdgesAroundEdge(ihalf_edgeB1_new, ihalf_edgeB2_new2);
    this->_LinkTwoHalfEdgesAroundEdge(ihalf_edgeB1_new2, ihalf_edgeB2_new);

    if (is_boundary_edgeC) 
      { ipoly[0] = ipoly[1]; }
    else {
      HALF_EDGE_INDEX_TYPE ihalf_edgeC0_new;
      ipoly[0] = SplitPolygonHalfEdge
        (ihalf_edgeC0, iv_splitC, ihalf_edgeC0_new);

      const HALF_EDGE_INDEX_TYPE ihalf_edgeC0_new2 =
        this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeC0_new);

      this->_LinkTwoHalfEdgesAroundEdge(ihalf_edgeC1_new, ihalf_edgeC0_new2);
      this->_LinkTwoHalfEdgesAroundEdge(ihalf_edgeC1_new2, ihalf_edgeC0_new);
    }

    if (is_boundary_edgeA) 
      { ipoly[3] = ipoly[2]; }
    else {
      HALF_EDGE_INDEX_TYPE ihalf_edgeA3_new;
      ipoly[3] = SplitPolygonHalfEdge
        (ihalf_edgeA3, iv_splitA, ihalf_edgeA3_new);

      const HALF_EDGE_INDEX_TYPE ihalf_edgeA3_new2 =
        this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeA3_new);

      this->_LinkTwoHalfEdgesAroundEdge(ihalf_edgeA2_new, ihalf_edgeA3_new2);
      this->_LinkTwoHalfEdgesAroundEdge(ihalf_edgeA2_new2, ihalf_edgeA3_new);
    }

    // Call UpdateHalfEdges*() after all old polygons are deleted
    //   in SplitTwoAdjacentPolygonHalfEdges() or SplitPolygonHalfEdge().
    this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly[1]);
    this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly[2]);
    if (!is_boundary_edgeC) 
      { this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly[0]); }
    if (!is_boundary_edgeA) 
      { this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly[3]); }
  }


  // Add vertices and split three polygon edges that are incident
  //   on the same vertex.
  // - Modifies four polygons.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, 
            typename ITYPEA, typename ITYPEB, typename ITYPEC, 
            typename ITYPEP>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVerticesSplitThreeIncidentEdgesIV
  (const ITYPE2 ihalf_edgeB1, 
   ITYPEA & iv_splitA, ITYPEB & iv_splitB, ITYPEC & iv_splitC, 
   ITYPEP ipoly[4])
  {
    iv_splitA = this->AddVertex();
    iv_splitB = this->AddVertex();
    iv_splitC = this->AddVertex();
    SplitThreeIncidentEdgesIV
      (ihalf_edgeB1, iv_splitA, iv_splitB, iv_splitC, ipoly);
  }


  // Split triangle and its half edges.
  // - Split triangle into four subtriangles.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitTriangleAndHalfEdges
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split[3], 
   ITYPE4 itriangle_new[4])
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeB);
    const POLYGON_INDEX_TYPE itriangle =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    VERTEX_INDEX_TYPE iv[3];

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    this->PrintPolygonIndexAndVertices(cerr, "  triangle: ", itriangle, "\n");
    */

    iv[0] = this->FromVertexIndex(ihalf_edgeA);
    iv[1] = this->FromVertexIndex(ihalf_edgeB);
    iv[2] = this->FromVertexIndex(ihalf_edgeC);

    this->_AddFourTriangles_T024
      (iv_split[0], iv[1], iv_split[1], iv[2], iv_split[2], iv[0],
       itriangle_new);

    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeA);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeB);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeC);

    this->_DeletePolygon(itriangle);
  }


  // Split triangle and its edges.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitTriangleAndEdgesIV
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split[3], 
   ITYPE4 itriangle_new[4], ITYPE5 ipoly_new[3])
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeB);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeAX =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeBX =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeB);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeCX =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeC);
    const bool is_boundary_edgeA =
      this->IsBoundaryEdge(ihalf_edgeA);
    const bool is_boundary_edgeB =
      this->IsBoundaryEdge(ihalf_edgeB);
    const bool is_boundary_edgeC =
      this->IsBoundaryEdge(ihalf_edgeC);
    HALF_EDGE_INDEX_TYPE ihalf_edgeAX_new, ihalf_edgeBX_new, ihalf_edgeCX_new;

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << endl;
    cerr << "In " << __func__ << endl;
    this->PrintPolygonIndexAndVertices
      (cerr, "  triangle: ", this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA), "\n");
    */

    this->SplitTriangleAndHalfEdges(ihalf_edgeA, iv_split, itriangle_new);

    // Initialize, in cse some triangle edge is a boundary edge.
    ipoly_new[0] = itriangle_new[0];
    ipoly_new[1] = itriangle_new[1];
    ipoly_new[2] = itriangle_new[2];

    if (!is_boundary_edgeA) {
      _SplitAndLinkEdgeOfPolygonAdjacentToSplitTriangle
        (0, ihalf_edgeAX, iv_split[0], itriangle_new, 
         ipoly_new[0], ihalf_edgeAX_new);
    }

    
    if (!is_boundary_edgeB) {
      _SplitAndLinkEdgeOfPolygonAdjacentToSplitTriangle
        (1, ihalf_edgeBX, iv_split[1], itriangle_new, 
         ipoly_new[1], ihalf_edgeBX_new);
    }

    
    if (!is_boundary_edgeC) {
      _SplitAndLinkEdgeOfPolygonAdjacentToSplitTriangle
        (2, ihalf_edgeCX, iv_split[2], itriangle_new, 
         ipoly_new[2], ihalf_edgeCX_new);
    }

    // Call UpdateHalfEdges*() after all old polygons are deleted
    //   in SplitTriangleAndHalfEdges() or _SplitAndLinkEdge*.
    this->UpdateHalfEdgesIncidentOnPolygonVertices(itriangle_new[0]);
    this->UpdateHalfEdgesIncidentOnPolygonVertices(itriangle_new[1]);
    this->UpdateHalfEdgesIncidentOnPolygonVertices(itriangle_new[2]);
    this->UpdateHalfEdgesIncidentOnPolygonVertices(itriangle_new[3]);
    if (!is_boundary_edgeA) 
      { this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly_new[0]); }
    if (!is_boundary_edgeB)
      { this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly_new[1]); }
    if (!is_boundary_edgeC)
      { this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly_new[2]); }
  }


  // Split and link edge of polygon adjacent to split triangle.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4, 
            typename ITYPE5, typename ITYPE6, typename ITYPE7>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  _SplitAndLinkEdgeOfPolygonAdjacentToSplitTriangle
  (const ITYPE2 ivloc0, const ITYPE3 ihalf_edgeX, 
   const ITYPE4 iv_split, const ITYPE5 itriangle_new[4], 
   ITYPE6 & ipoly_new, ITYPE7 & ihalf_edgeX_new0)
  {
    const ITYPE2 NUM_VERTICES_IN_TRIANGLE(3);
    const ITYPE2 ivloc1 = (ivloc0+1) % NUM_VERTICES_IN_TRIANGLE;

    const HALF_EDGE_INDEX_TYPE ihalf_edge_newC =
      this->HalfEdgeIndex(itriangle_new[ivloc0],0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge_newD =
      this->HalfEdgeIndex(itriangle_new[ivloc1],2);

    ipoly_new = this->SplitAndLinkPolygonHalfEdge
      (ihalf_edgeX, iv_split, ihalf_edge_newC, ihalf_edge_newD,
       ihalf_edgeX_new0);
  }


  // Add vertices and split triangle and its edges.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3,
            typename ITYPE4, typename ITYPE5>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVerticesSplitTriangleAndEdgesIV
  (const ITYPE2 ihalf_edgeA, ITYPE3 iv_split[3], 
   ITYPE4 itriangle_new[4], ITYPE5 ipoly_new[3])
  {
    iv_split[0] = this->AddVertex();
    iv_split[1] = this->AddVertex();
    iv_split[2] = this->AddVertex();
    SplitTriangleAndEdgesIV(ihalf_edgeA, iv_split, itriangle_new, ipoly_new);
  }


  /// Split pentagon half edge and triangulate resulting hexagon
  ///   into 4 triangles.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3,
            typename ITYPE4>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  SplitPentagonHalfEdgeAndPentagon
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iv_split,
   ITYPE4 itriangle_new[4])
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeB);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeD =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeC);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeE =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeD);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeBX =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeB);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeCX =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeC);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeDX =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeD);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeEX =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeE);

    const POLYGON_INDEX_TYPE ipentagon =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    VERTEX_INDEX_TYPE iv[5];
    IJK::PROCEDURE_ERROR error("SplitPentagonHalfEdgeAndPentagon");

    if (!this->IsPentagon(ipentagon)) { 
      error.AddMessage
        ("  Polygon ", ipentagon, "containing half edge ", ihalf_edgeA,
         " is not a pentagon.");
      error.AddMessage
        ("  Polygon ", ipentagon, " has ", 
         this->NumPolygonVertices(ipentagon), " vertices.");
      throw error;
    }      

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    this->PrintPolygonIndexAndVertices(cerr, "  ipentagon: ", ipentagon, "\n");
    */

    iv[0] = this->FromVertexIndex(ihalf_edgeA);
    iv[1] = this->FromVertexIndex(ihalf_edgeB);
    iv[2] = this->FromVertexIndex(ihalf_edgeC);
    iv[3] = this->FromVertexIndex(ihalf_edgeD);
    iv[4] = this->FromVertexIndex(ihalf_edgeE);

    this->_AddFourTriangles_T024
      (iv_split, iv[1], iv[2], iv[3], iv[4], iv[0],
       itriangle_new);

    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeA);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeB);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeC);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeD);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edgeE);

    this->_LinkTwoHalfEdgesAroundEdge
      (ihalf_edgeBX, this->HalfEdgeIndex(itriangle_new[1],0));
    this->_LinkTwoHalfEdgesAroundEdge
      (ihalf_edgeCX, this->HalfEdgeIndex(itriangle_new[2],2));
    this->_LinkTwoHalfEdgesAroundEdge
      (ihalf_edgeDX, this->HalfEdgeIndex(itriangle_new[2],0));
    this->_LinkTwoHalfEdgesAroundEdge
      (ihalf_edgeEX, this->HalfEdgeIndex(itriangle_new[0],2));

    this->_DeletePolygon(ipentagon);
  }


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE cut ear member functions
  // *****************************************************************

  // Cut ear from polygon.
  // - Separate ear vertex FromVertexIndex(ihalf_edge1) from the
  //   rest of the polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  CutPolygonEar
  (const ITYPE2 ihalf_edge1, ITYPE3 & ihalf_edge1_new, 
     ITYPE4 & ihalf_edge2_new)
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edge0 =
      this->IndexOfPrevHalfEdgeInPolygon(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    HALF_EDGE_INDEX_TYPE ihalf_edge0_new;

    // Note: Input to SplitPolygonByDiagonal is two half edges
    //   whose from vertices form the diagonal.
    SplitPolygonByDiagonal
      (ihalf_edge0, ihalf_edge2, ihalf_edge0_new, ihalf_edge2_new);

    ihalf_edge1_new = this->IndexOfNextHalfEdgeInPolygon(ihalf_edge0_new);
  }


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE replace triangle member functions
  // *****************************************************************

  // Replace triangle.
  // @param ihalf_edge0 Base edge of triangle.
  //   - Replace vertex opposite ihalf_edge0 with iv_new.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  typename MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  POLYGON_INDEX_TYPE
  MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  ReplaceTriangle
  (const ITYPE2 ihalf_edge0, const ITYPE3 iv_new, 
   ITYPE4 & ihalf_edge_new0)
  {
    const NTYPE NUM_VERTICES_PER_TRIANGLE(3);
    const POLYGON_INDEX_TYPE itriangle =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge1 = 
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2 = 
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2X = 
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge2);
    const HALF_EDGE_INDEX_TYPE ihalf_edge3X = 
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge2X);
    const ITYPE iloc0 = this->LocationOfHalfEdgeInPolygon(ihalf_edge0);
    HALF_EDGE_INDEX_TYPE ihalf_edge_new1, ihalf_edge_new2;
    IJK::PROCEDURE_ERROR error("ReplaceTriangle");

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintPolygonIndexAndVertices(cerr, "Replacing triangle ", itriangle, "\n");
    this->PrintHalfEdgeIndexAndEndpoints(cerr, "  half_edge0: ", ihalf_edge0, "\n");
    const VERTEX_INDEX_TYPE iw0 = this->FromVertexIndex(ihalf_edge0);
    const VERTEX_INDEX_TYPE iw1 = this->FromVertexIndex(ihalf_edge1);
    const VERTEX_INDEX_TYPE iw2 = this->FromVertexIndex(ihalf_edge2);
    cerr << "  NumPolygonsIncidentOnVertex(" << iw0 << "): "
         << this->NumPolygonsIncidentOnVertex(iw0) << endl;
    cerr << "  NumPolygonsIncidentOnVertex(" << iw1 << "): "
         << this->NumPolygonsIncidentOnVertex(iw1) << endl;
    cerr << "  NumPolygonsIncidentOnVertex(" << iw2 << "): "
         << this->NumPolygonsIncidentOnVertex(iw2) << endl;
    */

    if (!this->IsTriangle(itriangle)) { 
      error.AddMessage
        ("  Polygon ", itriangle, "containing half edge ", ihalf_edge0,
         " is not a triangle.");
      error.AddMessage
        ("  Polygon ", itriangle, " has ", 
         this->NumPolygonVertices(itriangle), " vertices.");
      throw error;
    }      

    // Edge ihalf_edge2 will be deleted.
    // Replace half edge incident on vertex FromVertexIndex(ihalf_edge2),
    //   if necessary.
    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edge2, ihalf_edge3X);

    const POLYGON_INDEX_TYPE itriangle_new =
      this->_AddPolygonX(NUM_VERTICES_PER_TRIANGLE);

    ihalf_edge_new0 = this->HalfEdgeIndex(itriangle_new, iloc0);
    this->_SetHalfEdgeFromVertex
      (ihalf_edge_new0, this->FromVertexIndex(ihalf_edge0));

    this->_ReplaceHalfEdge(ihalf_edge0, ihalf_edge_new0);

    this->_UnlinkHalfEdgeAroundEdge(ihalf_edge1);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edge2);

    // Note: IndexOfNextHalfEdgeInPolygon(ihalf_edge_new1) will be
    //  incorrect, until ihalf_edge_new1 is initialized.
    ihalf_edge_new1 = this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new0);
    ihalf_edge_new2 = 
      this->IndexOfKthNextHalfEdgeInPolygon(ihalf_edge_new0, 2);
    this->_SetHalfEdgeFromVertex
      (ihalf_edge_new1, this->ToVertexIndex(ihalf_edge0));
    this->_SetHalfEdgeFromVertex(ihalf_edge_new2, iv_new);

    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edge1, ihalf_edge_new1);

    this->_DeletePolygon(itriangle);

    return(itriangle_new);
  }


  // Replace pentagon with triangle.
  // - Merge two pentagon vertices to form triangle.
  // @param ihalf_edge0 Pentagon edge that becomes base edge of triangle.
  // @param iv_new Replace merged vertices with iv_new.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  typename MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  POLYGON_INDEX_TYPE
  MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  ReplacePentagonWithTriangle
  (const ITYPE2 ihalf_edge0, const ITYPE3 iv_new, 
   ITYPE4 & ihalf_edge_new0)
  {
    const NUMBER_TYPE NUM_VERTICES_PER_TRIANGLE(3);
    const POLYGON_INDEX_TYPE ipentagon =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    HALF_EDGE_INDEX_TYPE ihalf_edge1, ihalf_edge2, ihalf_edge3, ihalf_edge4;
    HALF_EDGE_INDEX_TYPE ihalf_edge_new1, ihalf_edge_new2;
    IJK::PROCEDURE_ERROR error("ReplacePentagonWithTriangle");

    if (!this->IsPentagon(ipentagon)) { 
      error.AddMessage
        ("  Polygon ", ipentagon, "containing half edge ", ihalf_edge0,
         " is not a pentagon.");
      error.AddMessage
        ("  Polygon ", ipentagon, " has ", 
         this->NumPolygonVertices(ipentagon), " vertices.");
      throw error;
    }      

    this->GetIndicesOfNextFourHalfEdgesInPolygon
      (ihalf_edge0, ihalf_edge1, ihalf_edge2, ihalf_edge3, ihalf_edge4);

    const HALF_EDGE_INDEX_TYPE ihalf_edge3X = 
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge3);
    const HALF_EDGE_INDEX_TYPE ihalf_edge3XN = 
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge3X);

    // Edge ihalf_edge3 will be deleted.
    // Replace half edge incident on vertex FromVertexIndex(ihalf_edge3),
    //   if necessary.
    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edge3, ihalf_edge3XN);

    const POLYGON_INDEX_TYPE itriangle_new =
      this->_AddPolygonX(NUM_VERTICES_PER_TRIANGLE);

    // New half edges start from iv_new.
    ihalf_edge_new0 = this->HalfEdgeIndex(itriangle_new, 2);
    this->_SetHalfEdgeFromVertex
      (ihalf_edge_new0, this->FromVertexIndex(ihalf_edge0));
    this->_ReplaceHalfEdge(ihalf_edge0, ihalf_edge_new0);

    this->_UnlinkHalfEdgeAroundEdge(ihalf_edge1);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edge2);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edge3);
    this->_UnlinkHalfEdgeAroundEdge(ihalf_edge4);

    // Note: IndexOfNextHalfEdgeInPolygon(ihalf_edge_new1) will be
    //  incorrect, until ihalf_edge_new1 is initialized.
    ihalf_edge_new1 = this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new0);
    ihalf_edge_new2 = 
      this->IndexOfKthNextHalfEdgeInPolygon(ihalf_edge_new0, 2);
    this->_SetHalfEdgeFromVertex
      (ihalf_edge_new1, this->ToVertexIndex(ihalf_edge0));
    this->_SetHalfEdgeFromVertex(ihalf_edge_new2, iv_new);

    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edge0, ihalf_edge_new0);
    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edge1, ihalf_edge_new1);

    this->_DeletePolygon(ipentagon);

    return(itriangle_new);
  }


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE Delete member functions
  // *****************************************************************


  // Delete edge contained in two polygon.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  DeleteEdge(const ITYPE2 ihalf_edgeA, ITYPE3 & ipoly_new)
  {
    const POLYGON_INDEX_TYPE ipolyA =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);
    const POLYGON_INDEX_TYPE ipolyB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeB);
    const NUMBER_TYPE numvA = this->NumPolygonVertices(ipolyA);
    const NUMBER_TYPE numvB = this->NumPolygonVertices(ipolyB);

    const NUMBER_TYPE numv = numvA + numvB - 2;

    ipoly_new = this->_AddPolygonX(numv);

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintHalfEdge(cerr, "Deleting edge", ihalf_edgeA, "\n");
    */

    const HALF_EDGE_INDEX_TYPE ihalf_edgeC =
      this->polygon_half_edge_list.first_element[ipoly_new];
    const HALF_EDGE_INDEX_TYPE ihalf_edgeD = ihalf_edgeC+numvA-1;

    // Set half edges from polygon ipolyA.
    const HALF_EDGE_INDEX_TYPE ihalf_edgeAN =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeA);

    this->_ReplaceHalfEdges(ihalf_edgeAN, ihalf_edgeC, numvA-1);

    // Set half edges from polygon ipolyB.
    const HALF_EDGE_INDEX_TYPE ihalf_edgeBN =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeB);

    this->_ReplaceHalfEdges(ihalf_edgeBN, ihalf_edgeD, numvB-1);

    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edgeA, ihalf_edgeD);
    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edgeB, ihalf_edgeC);

    this->_DeletePolygon(ipolyA);
    this->_DeletePolygon(ipolyB);

    // Call UpdateHalfEdges*() after both old polygons are deleted.
    this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly_new);

    return;
  }


  // Delete edge and add triangle.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, 
            typename ITYPE4, typename ITYPE5>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  DeleteEdgeAddTriangle
  (const ITYPE2 ihalf_edgeA, const ITYPE3 iear,
   ITYPE4 & itriangle_new, ITYPE5 & ipoly_new)
  {
    const NTYPE NUM_VERTICES_PER_TRIANGLE(3);
    const NTYPE NUM_VERTICES_PER_QUAD(4);
    const POLYGON_INDEX_TYPE iquadA =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);
    const POLYGON_INDEX_TYPE ipolyB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edgeB);
    const NUMBER_TYPE numvB = this->NumPolygonVertices(ipolyB);
    const VERTEX_INDEX_TYPE jv_ear = this->PolygonVertex(iquadA, iear);
    const VERTEX_INDEX_TYPE iv0 = this->FromVertexIndex(ihalf_edgeA);
    const VERTEX_INDEX_TYPE iv1 = this->ToVertexIndex(ihalf_edgeA);
    const ITYPE j0 = 
      (iear+(NUM_VERTICES_PER_QUAD-1))%NUM_VERTICES_PER_QUAD;
    const ITYPE j1 = (iear+1)%NUM_VERTICES_PER_QUAD;
    VERTEX_INDEX_TYPE jv0 = this->PolygonVertex(iquadA,j0);
    VERTEX_INDEX_TYPE jv1 = this->PolygonVertex(iquadA,j1);
    IJK::PROCEDURE_ERROR error("DeleteEdgeAddTriangle");

    if (!this->IsQuadrilateral(iquadA)) { 
      error.AddMessage
        ("  Polygon ", iquadA, "containing half edge ", ihalf_edgeA,
         " is not a quadrilateral.");
      error.AddMessage
        ("  Polygon ", iquadA, " has ", 
         this->NumPolygonVertices(iquadA), " vertices.");
      throw error;
    }      

    // *** DEBUG ***
    /*
    using namespace std;
    this->PrintHalfEdge(cerr, "Deleting edge", ihalf_edgeA, "\n");
    cerr << "  Adding triangle (" << jv0 << "," 
	 << jv_ear << "," << jv1 << ")" << endl;
    */

    itriangle_new = this->_AddPolygonX(NUM_VERTICES_PER_TRIANGLE);
    ipoly_new = this->_AddPolygonX(numvB+1);


    // Store triangle half edges.

    const HALF_EDGE_INDEX_TYPE ihalf_edgeT =
      this->polygon_half_edge_list.first_element[itriangle_new];

    const HALF_EDGE_INDEX_TYPE ihalf_edgeAN =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeA);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeANN =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeAN);
    const HALF_EDGE_INDEX_TYPE ihalf_edgeAP =
      this->IndexOfPrevHalfEdgeInPolygon(ihalf_edgeA);

    if (jv0 == iv1) {
      // Replace two edges.
      this->_ReplaceHalfEdges(ihalf_edgeAN, ihalf_edgeT, 2);
    }
    else {
      // Replace two edges.
      this->_ReplaceHalfEdges(ihalf_edgeANN, ihalf_edgeT, 2);
    }

    // Set last half edge in itriangle_new.
    this->_SetHalfEdgeFromVertex(ihalf_edgeT+2, jv1);

    // Store new polygon half edges.

    const HALF_EDGE_INDEX_TYPE ihalf_edgeC =
      this->polygon_half_edge_list.first_element[ipoly_new];

    const HALF_EDGE_INDEX_TYPE ihalf_edgeBN =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edgeB);

    this->_ReplaceHalfEdges(ihalf_edgeBN, ihalf_edgeC, numvB-1);

    // Set last two half edges in inew polygon.
    const HALF_EDGE_INDEX_TYPE ihalf_edgeC_new = ihalf_edgeC+numvB-1;
    this->_SetHalfEdgeFromVertex(ihalf_edgeC_new, iv1);

    if (jv0 == iv1) {
      this->_SetHalfEdgeFromVertex(ihalf_edgeC_new+1, jv1);
      this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edgeAN, ihalf_edgeC_new);
      this->_LinkTwoHalfEdgesAroundEdge
        (this->HalfEdgeIndex(itriangle_new, 2), ihalf_edgeC_new);

      this->_ReplaceHalfEdge(ihalf_edgeAP, ihalf_edgeC_new+1);
    }
    else {
      this->_ReplaceHalfEdge(ihalf_edgeAN, ihalf_edgeC_new);
      this->_SetHalfEdgeFromVertex(ihalf_edgeC_new+1, jv0);
      this->_LinkTwoHalfEdgesAroundEdge
        (this->HalfEdgeIndex(itriangle_new, 2), ihalf_edgeC_new+1);
    }

    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edgeA, ihalf_edgeC);
    this->_ReplaceHalfEdgeIncidentOnVertex(ihalf_edgeB, ihalf_edgeC_new);

    this->_DeletePolygon(iquadA);
    this->_DeletePolygon(ipolyB);

    // Call UpdateHalfEdges*() after both old polygons are deleted.
    this->UpdateHalfEdgesIncidentOnPolygonVertices(ipoly_new);
    this->UpdateHalfEdgesIncidentOnPolygonVertices(itriangle_new);

    return;

  }


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE Replace triangle functions
  // *****************************************************************


  // Replace triangle with triangle edge.
  // @param ihalf_edge0 Base of triangle.
  // @param iv_new New triangle vertex.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  ReplaceTriangleWithTriangleEdge
  (const ITYPE2 ihalf_edge0, const ITYPE3 iv_new,
   ITYPE4 & ihalf_edge_new0)
  {
    const HALF_EDGE_INDEX_TYPE ihalf_edge1 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge1X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge2);
    const POLYGON_INDEX_TYPE ipolyA =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge1X);
    const POLYGON_INDEX_TYPE itriangleB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    const POLYGON_INDEX_TYPE ipolyC =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge2X);
    const NUMBER_TYPE numvA = this->NumPolygonVertices(ipolyA);
    const NUMBER_TYPE numvC = this->NumPolygonVertices(ipolyC);
    IJK::PROCEDURE_ERROR error("ReplaceTriangleWithTriangleEdge");

    if (!this->IsTriangle(itriangleB)) { 
      error.AddMessage
        ("  Polygon ", itriangleB, "containing half edge ", ihalf_edge0,
         " is not a triangle.");
      error.AddMessage
        ("  Polygon ", itriangleB, " has ", 
         this->NumPolygonVertices(itriangleB), " vertices.");
      throw error;
    }      

    // Replace triangle.
    const POLYGON_INDEX_TYPE itriangle_new =
      ReplaceTriangle(ihalf_edge0, iv_new, ihalf_edge_new0);

    const HALF_EDGE_INDEX_TYPE ihalf_edge_new1 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge_new2 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new1);

    // Split half edges of adjacent polygons, replacing adjacent polygons.

    HALF_EDGE_INDEX_TYPE ihalf_edge_new1X, ihalf_edge_new2X;

    const POLYGON_INDEX_TYPE ipolyA_new = 
      SplitPolygonHalfEdge(ihalf_edge1X, iv_new, ihalf_edge_new1X);
    const POLYGON_INDEX_TYPE ipolyB_new = 
      SplitPolygonHalfEdge(ihalf_edge2X, iv_new, ihalf_edge_new2X);

    // Link half edges of new triangle and polygons.
    const HALF_EDGE_INDEX_TYPE ihalf_edge_new1XN =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new1X);
    const HALF_EDGE_INDEX_TYPE ihalf_edge_new2XN =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new2X);

    this->_LinkTwoHalfEdgesAroundEdge
      (ihalf_edge_new1X, ihalf_edge_new2XN);
    this->_LinkTwoHalfEdgesAroundEdge
      (ihalf_edge_new1, ihalf_edge_new1XN);
    this->_LinkTwoHalfEdgesAroundEdge
      (ihalf_edge_new2, ihalf_edge_new2X);

    // Call UpdateHalfEdges*() after both old polygons are deleted.
    this->UpdateHalfEdgesIncidentOnPolygonVertices(itriangle_new);
    this->UpdateHalfEdgesIncidentOnPolygonVertices(ipolyA_new);
    this->UpdateHalfEdgesIncidentOnPolygonVertices(ipolyB_new);

  }


  // Add vertex and replace triangle with triangle edge.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexReplaceTriangleWithTriangleEdge
  (const ITYPE2 ihalf_edge0, ITYPE3 & iv_new,
     ITYPE4 & ihalf_edge_new0)
  {
    iv_new = this->AddVertex();
    ReplaceTriangleWithTriangleEdge
      (ihalf_edge0, iv_new, ihalf_edge_new0);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "New vertex: " << iv_new << endl;
    */
  }




  // Replace split edge triangle with triangle edge.
  // - Triangle has two split edges.
  // @param ihalf_edge0 Base of triangle.
  // @param iv_new New triangle vertex.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  ReplaceSplit2EdgeTriangleWithTriangleEdge
  (const ITYPE2 ihalf_edge0, const ITYPE3 iv_new,
   ITYPE4 & ihalf_edge_new0)
  {
    HALF_EDGE_INDEX_TYPE ihalf_edge1, ihalf_edge2, ihalf_edge3, ihalf_edge4;

    this->GetIndicesOfNextFourHalfEdgesInPolygon
      (ihalf_edge0, ihalf_edge1, ihalf_edge2, ihalf_edge3, ihalf_edge4);


    const HALF_EDGE_INDEX_TYPE ihalf_edge1X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE ihalf_edge2X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge2);
    const HALF_EDGE_INDEX_TYPE ihalf_edge3X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge3);
    const HALF_EDGE_INDEX_TYPE ihalf_edge4X =
      this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge4);
    const POLYGON_INDEX_TYPE itriangleB =
      this->IndexOfPolygonContainingHalfEdge(ihalf_edge0);
    IJK::PROCEDURE_ERROR error("ReplaceSplit2EdgeTriangleWithTriangleEdge");

    if (!this->IsPentagon(itriangleB)) { 
      error.AddMessage
        ("  Polygon ", itriangleB, "containing half edge ", ihalf_edge0,
         " is not a triangle with two split edges.");
      error.AddMessage
        ("  Polygon ", itriangleB, " has ", 
         this->NumPolygonVertices(itriangleB), " vertices.");
      error.AddMessage
        ("  Polygon ", itriangleB, " should have 5 vertices.");
      throw error;
    }

    if (this->AreFourOrMoreHalfEdgesAroundFromVertex(ihalf_edge2)) {
      error.AddMessage
        ("Programming error. Four or more half edges around vertex ",
         this->FromVertexIndex(ihalf_edge2), ".");
      throw error;
    } else if (this->AreFourOrMoreHalfEdgesAroundFromVertex(ihalf_edge4)) {
      error.AddMessage
        ("Programming error. Four or more half edges around vertex ",
         this->FromVertexIndex(ihalf_edge4), ".");
      throw error;
    }

    // Replace split edge triangle with triangle.
    const POLYGON_INDEX_TYPE itriangle_new = 
      ReplacePentagonWithTriangle(ihalf_edge0, iv_new, ihalf_edge_new0);

    const HALF_EDGE_INDEX_TYPE ihalf_edge_new1 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new0);
    const HALF_EDGE_INDEX_TYPE ihalf_edge_new2 =
      this->IndexOfNextHalfEdgeInPolygon(ihalf_edge_new1);

    // Replace vertices of polygons with iv_new.
    // Sets half edges incident on iv_new.
    // Changes half edges incident on FromVertex(ihalf_edge1X) 
    //   and FromVertex(ihalf_edge3X).
    ReplaceFromVertex(ihalf_edge1X, iv_new);
    ReplaceFromVertex(ihalf_edge3X, iv_new);

    // Replace vertices incident on ToVertex(ihalf_edge2X) and
    //   ToVertex(ihalf_edge4X).
    // Called in case different polygons contain ihalf_edge1X 
    //   and ihalf_edge2X or different polygons contain
    //   ihalf_edge3X and ihalf_edge4X.
    ReplaceToVertex(ihalf_edge2X, iv_new);
    ReplaceToVertex(ihalf_edge4X, iv_new);

    // Link half edges of polygons ipolyA and ipolyC.
    this->_LinkTwoHalfEdgesAroundEdge(ihalf_edge2X, ihalf_edge3X);

    // Link half edges of new triangle and polygons.
    this->_LinkTwoHalfEdgesAroundEdge(ihalf_edge1X, ihalf_edge_new1);
    this->_LinkTwoHalfEdgesAroundEdge(ihalf_edge4X, ihalf_edge_new2);

    // Call UpdateHalfEdges*() after all old polygons are deleted.
    this->UpdateHalfEdgesIncidentOnPolygonVertices(itriangle_new);

    return;
  }


  // Add vertex and replace split edge triangle with triangle edge.
  // - Triangle has two split edges.
  template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
            typename ITYPE, typename NTYPE>
  template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
  void MESH2D_SPLIT_BASE<VTYPE, HALFE_TYPE, PTYPE, ITYPE, NTYPE>::
  AddVertexReplaceSplit2EdgeTriangleWithTriangleEdge
  (const ITYPE2 ihalf_edge0, ITYPE3 & iv_new,
   ITYPE4 & ihalf_edge_new0)
  {
    iv_new = this->AddVertex();
    ReplaceSplit2EdgeTriangleWithTriangleEdge
      (ihalf_edge0, iv_new, ihalf_edge_new0);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    cerr << "New vertex: " << iv_new << endl;
    */
  }


  // *****************************************************************
  // Class MESH2D_SPLIT_BASE Check member functions
  // *****************************************************************

  // None yet.

}


#endif
