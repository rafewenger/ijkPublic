/*!
 *  @file ijkmesh2D_datastruct.tpp
 *  @brief ijk template classes for 2D polygonal mesh data structures.
 *  - Uses half-edges.  Totally distinct from POLYMESH.
 *  - Version 0.4.0
 */


  /*
    IJK: Isosurface Jeneration Kode
    Copyright (C) 2019-2024 Rephael Wenger

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


  #ifndef _IJKMESH2D_DATASTRUCT_
  #define _IJKMESH2D_DATASTRUCT_

  #include "ijk.tpp"
  #include "ijklist.tpp"
  #include "ijkmesh.tpp"

  #include <algorithm>
  #include <numeric>
  #include <tuple>
  #include <utility>
  #include <vector>


  // *** DEBUG ***
  #include "ijkprint.tpp"


  namespace IJK {

    // *****************************************************************
    // Class HALF_EDGE_BASE
    // *****************************************************************

    /// @brief Half edge base type.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE>
    class HALF_EDGE_BASE {

    public:

      /// @brief Half edge index type.
      typedef HALFE_TYPE HALF_EDGE_INDEX_TYPE;

      /// @brief Polygon index type.
      typedef PTYPE POLYGON_INDEX_TYPE;

    public:

      POLYGON_INDEX_TYPE polygon_index;

      /*!
       *   - If the index equals the index of this half edge,
       *     the half edge is a boundary edge.
      */
      HALF_EDGE_INDEX_TYPE index_of_next_half_edge_around_edge;

    public:

      /// @brief Return polygon index.
      PTYPE PolygonIndex() const
      { return(polygon_index); }

      /// @brief Return index of next half edge around edge.
      HALFE_TYPE IndexOfNextHalfEdgeAroundEdge() const
      { return(index_of_next_half_edge_around_edge); }
    };


    // *****************************************************************
    // Class VERTEX_BASE
    // *****************************************************************

    /// @brief Vertex base type.
    template <typename VI_TYPE, typename HTYPE, typename NTYPE>
    class VERTEX_BASE {

    public:

      /// Half edge index type.
      typedef HTYPE HALF_EDGE_INDEX_TYPE;

      /// Vertex index type.
      typedef VI_TYPE VERTEX_INDEX_TYPE;

    protected:

      /*!
       *  @brief Index of some half edge with *this as from vertex.
       *  - Index of boundary half edge if *this is a from vertex
       *    of some boundary half edge.
       */
      HTYPE incident_half_edge_index;

      /// @brief If true, incident_half_edge_index is set.
      bool is_set_incident_half_edge_index;


    public:

      /// @brief Constructor.
      VERTEX_BASE() {
        is_set_incident_half_edge_index = false;
      };

      /// @brief Return true if incident half edge index is set.
      bool IsSetIncidentHalfEdgeIndex() const
      { return(is_set_incident_half_edge_index); }

      /// @brief Return index of some half edge index incident on vertex.
      HTYPE IncidentHalfEdgeIndex() const
      { return(incident_half_edge_index); }

      /// @brief Set incident half edge index to ihalf_edge.
      template <typename ITYPE>
      void SetIncidentHalfEdgeIndex(const ITYPE ihalf_edge)
      {
        incident_half_edge_index = ihalf_edge;
        is_set_incident_half_edge_index = true;
      }

      /// @brief Unset incident half edge index.
      void ClearIncidentHalfEdgeIndex()
      { is_set_incident_half_edge_index = false; }

      /// @brief Clear incident half edge index, if it is set
      ///   and equals ihalf_edge.
      template <typename ITYPE>
      void ClearIncidentHalfEdgeIndexIfEquals(const ITYPE ihalf_edge)
      {
        if (IsSetIncidentHalfEdgeIndex()) {
          if (ihalf_edge == incident_half_edge_index)
            { ClearIncidentHalfEdgeIndex(); }
        }
      }

    };


    // *****************************************************************
    // Class POLYGON_BASE
    // *****************************************************************

    /// @brief Polygon base type.
    class POLYGON_BASE {

    public:
      /// @brief Deletion flag.
      bool is_deleted;

    public:
      /// @brief Constructor.
      POLYGON_BASE() { is_deleted = false; }

      /// @brief Return true if polygon is deleted.
      bool IsDeleted() const
      { return(is_deleted); }

      /*!
       * @brief Initialize polygon.
       * @param num_edges Argument num_edges is ignored by the default function.
       *   - Argument num_edges may be used in derived classes to store
       *     number of polygon  edges.
      */
      template <typename NTYPE2>
      void Init(const NTYPE2 num_edges) {
        // Argument num_edges not used.
        // May be reset by derived classes to store number of initial edges.
        is_deleted = false;
      }

    };


    // *****************************************************************
    // Class MESH2D_BASE
    // *****************************************************************

    /*!
     * @brief Mesh of 2D polygons.
     * @tparam VTYPE Vertex type.  Should be derived from VERTEX_BASE.
     * @tparam HALFE_TYPE Half edge type.
     *   - Should be derived from HALF_EDGE_BASE.
     * @tparam PTYPE POLYGON type.  Should be derived from POLYGON_BASE.
    */
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    class MESH2D_BASE {

    public:

      /// @brief Vertex type.
      typedef VTYPE VERTEX_TYPE;

      /// @brief Half edge type.
      typedef HALFE_TYPE HALF_EDGE_TYPE;

      /// @brief Vertex index type.
      typedef typename VERTEX_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

      /// @brief Half edge index type.
      typedef typename HALF_EDGE_TYPE::HALF_EDGE_INDEX_TYPE
      HALF_EDGE_INDEX_TYPE;

      /// @brief Polygon index type.
      typedef typename HALF_EDGE_TYPE::POLYGON_INDEX_TYPE POLYGON_INDEX_TYPE;

      /// @brief Polygon type.
      typedef PTYPE POLYGON_TYPE;

      /// @brief Index type.
      typedef ITYPE INDEX_TYPE;

      /// @brief Number type.
      typedef NTYPE NUMBER_TYPE;

    protected:

      /*!
       *  @brief List of polygon vertices.
       *  - Vertices are in clockwise or counter-clockwise order
       *    around each polygon.
       */
      LIST_OF_LISTS<VERTEX_INDEX_TYPE,NTYPE> polygon_vertex_list;

      /*!
       *  @brief List of polygon half edges.
       *  - Half edges are in clockwise or counter-clockwise order
       *    around each polygon.
       *  - From vertex of half edge i is polygon_vertex_list[i].
       */
      LIST_OF_LISTS<HALF_EDGE_TYPE,NTYPE> polygon_half_edge_list;

      /// List of vertices.
      std::vector<VERTEX_TYPE> vertex_list;

      /// List of polygons.
      std::vector<POLYGON_TYPE> polygon_list;

      /// If true, half edges are linked.
      bool are_half_edges_linked;

      /// If true, mesh is an oriented 2D-manifold
      bool is_oriented_manifold;

    protected:

      // Internal functions.

      /// Set index of next half edge around edge.
      template <typename ITYPE2, typename ITYPE3>
      void _SetIndexOfNextHalfEdgeAroundEdge
      (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1)
      {
        polygon_half_edge_list.element[ihalf_edge0].
          index_of_next_half_edge_around_edge = ihalf_edge1;
      }

      /// Set index of half edge incident on vertex.
      template <typename ITYPE2, typename ITYPE3>
      void _SetIndexOfHalfEdgeIncidentOnVertex
      (const ITYPE2 iv, const ITYPE3 ihalf_edge)
      { vertex_list[iv].SetIncidentHalfEdgeIndex(ihalf_edge); }

      /// Link two half edges around edge.
      template <typename ITYPE2, typename ITYPE3>
      void _LinkTwoHalfEdgesAroundEdge
      (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1);

      /// @brief Initialize half edge.
      /// - Does not set polygon vertices.
      template <typename ITYPE2, typename ITYPE3>
      void _InitHalfEdgeX
      (const ITYPE2 ihalf_edge, const ITYPE3 ipoly);

      /// Set half edge from vertex
      template <typename ITYPE2, typename ITYPE3>
      void _SetHalfEdgeFromVertex
      (const ITYPE2 ihalf_edge, const ITYPE3 iv)
      { polygon_vertex_list.element[ihalf_edge] = iv; }

      /// Set polygon vertices.
      template <typename ITYPE2, typename ITYPE3, typename NTYPE2>
      void _SetPolygonVertices
      (const ITYPE2 ipoly, const ITYPE3 poly_vert_list[], const NTYPE2 list_length);

      /// @brief Remove half edge ihalf_edge from half_edges incident
      ///   on FromVertex(ihalf_edge).
      /// - Decrement number of polygons incident on from vertex.
      template <typename ITYPE2>
      void _RemoveHalfEdgeIncidentOnFromVertex
      (const ITYPE2 ihalf_edge);

      /// @brief Allocate memory for k additional half edges.
      /// - No initialization or setting of half edge vertices.
      template <typename NTYPE2>
      void _AllocateAdditionalHalfEdges(const NTYPE2 k);

      /// @brief Add polygon with numv vertices.
      /// - Internal, protected version.
      /// - Does not set polygon vertices.
      template <typename NTYPE2>
      POLYGON_INDEX_TYPE _AddPolygonX(const NTYPE2 numv);

      /// @brief Add a polygon with list_length vertices.
      /// - Internal, protected version.
      /// - Does not modify flags are_half_edges_linked
      ///     or is_oriented_manifold.
      template <typename ITYPE2, typename NTYPE2>
      POLYGON_INDEX_TYPE _AddPolygon
      (const ITYPE2 poly_vert_list[], const NTYPE2 list_length);

      /// @brief Add two polygon with numv0 and numv1 vertices.
      /// - Internal, protected version.
      /// - Does not set polygon edges.
      template <typename NTYPE2, typename NTYPE3,
                typename ITYPE2, typename ITYPE3>
      void _AddTwoPolygonsX
      (const NTYPE2 numv0, const NTYPE3 numv1,
       ITYPE2 & ipoly0, ITYPE3 & ipoly1);

      /// @brief Add a triangle.
      /// - Does not modify flags are_half_edges_linked
      ///     or is_oriented_manifold.
      template <typename ITYPE2>
      POLYGON_INDEX_TYPE _AddTriangle(const ITYPE2 triangle_vert_list[])
      { return(_AddPolygon(triangle_vert_list, 3)); }

      /// @brief Add a triangle.
      /// - Does not modify flags are_half_edges_linked
      ///     or is_oriented_manifold.
      template <typename ITYPE2>
      POLYGON_INDEX_TYPE _AddTriangle
      (const ITYPE2 iv0, const ITYPE2 iv1, const ITYPE2 iv2);

      /// @brief Add a quadrilateral
      /// - Does not modify flags are_half_edges_linked
      ///     or is_oriented_manifold.
      template <typename ITYPE2>
      POLYGON_INDEX_TYPE _AddQuadrilateral
      (const ITYPE2 iv0, const ITYPE2 iv1, const ITYPE2 iv2, const ITYPE2 iv3);

      /*!
       * @brief Add four triangles forming the triangulation of a hexagon.
       * - Triangulation has 3 ears and one center triangle (iv0, iv2, iv4).
       * - Links edges in triangulation interior.
       * - Does not modify flags are_half_edges_linked
       *   or is_oriented_manifold.
      */
      /// @param itriangle_new[i] Index of i'th new triangle.
      template <typename ITYPE0, typename ITYPE1, typename ITYPE2,
                typename ITYPE3, typename ITYPE4, typename ITYPE5,
                typename ITYPE6>
      void _AddFourTriangles_T024
      (const ITYPE0 iv0, const ITYPE1 iv1, const ITYPE2 iv2,
       const ITYPE3 iv3, const ITYPE4 iv4, const ITYPE5 iv5,
       ITYPE6 itriangle_new[4]);

      /// Unlink a half edge from the list of half edges around an edge.
      template <typename ITYPE2>
      void _UnlinkHalfEdgeAroundEdge(const ITYPE2 ihalf_edge);

      /// Delete polygon.
      template <typename ITYPE2>
      void _DeletePolygon(const ITYPE2 ipoly);

      /// Replace half edge ihalf_edge0 with ihalf_edge1.
      template <typename ITYPE2, typename ITYPE3>
      void _ReplaceHalfEdge
      (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1);

      /// @brief Replace numR half edges.
      /// @param Start by replacing edge ihalf_edge0.
      /// @param Half edge ihalf_edge1 replaces ihalf_edge0.
      template <typename ITYPE2, typename ITYPE3, typename NTYPE2>
      void _ReplaceHalfEdges
      (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1,
       const NTYPE2 numR);

      /// Replace half edge by two half edges.
      template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
      void _ReplaceHalfEdgeWithTwoHalfEdges
      (const ITYPE2 ihalf_edge, const ITYPE3 iv_split,
       const ITYPE4 ihalf_edge_new0);

      /// Replace half edge in list of half edges around edge.
      template <typename ITYPE2, typename ITYPE3>
      void _ReplaceHalfEdgeAroundEdge
      (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1);

      /// @brief Replace half edge incident on vertex.
      /// - If incident_half_edge_index of FromVertexIndex(ihalf_edge0)
      ///   equals ihalf_edge0, replace by ihalf_edge1.
      /// @pre FromVertexIndex(ihalf_edge0) == FromVertexIndex(ihalf_edge1);
      template <typename ITYPE2, typename ITYPE3>
      void _ReplaceHalfEdgeIncidentOnVertex
      (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1);


    public:
      /// constructor
      MESH2D_BASE();

      // Query functions.

      /// Number of polygons (including deleted polygons.)
      NTYPE NumPolygons() const
      { return(polygon_list.size()); }

      /// Number of vertices (includes vertices that are in no polygons.)
      NTYPE NumVertices() const
      { return(vertex_list.size()); }

      /// Total number of half edge (includes edges of deleted polygons.)
      NTYPE NumHalfEdges() const
      { return(polygon_half_edge_list.element.size()); }

      /// @brief Return true if half-edges are linked.
      /// - Adding polygons sets this to false.
      bool AreHalfEdgesLinked() const
      { return(are_half_edges_linked); }

      /// @brief Return true if mesh is an oriented, 2D-manifold.
      /// - Adding polygons sets this to false.
      bool IsOrientedManifold() const
      { return(is_oriented_manifold); }

      /// Return reference to polygon_vertex_list.
      const LIST_OF_LISTS<VERTEX_INDEX_TYPE,NTYPE> &
      PolygonVertexList() const
      { return(polygon_vertex_list); }

      /// Return reference to vertex list of polygon ipoly.
      template <typename ITYPE2>
      const VERTEX_INDEX_TYPE *
      PolygonVertexList(const ITYPE2 ipoly) const
      { return(polygon_vertex_list.List(ipoly)); }

      // *** DEPRECATED ***
      /// @brief DEPRECATED: Return vertex j of polygon ipoly.
      template <typename ITYPE2, typename ITYPE3>
      VERTEX_INDEX_TYPE PolygonVertexList
      (const ITYPE2 ipoly, const ITYPE3 j) const
      { return(polygon_vertex_list.Element(ipoly,j)); }

      // *** DEPRECATED ***
      /// @brief DEPRECATED: Return vertex j of polygon ipoly.
      template <typename ITYPE2, typename ITYPE3>
      VERTEX_INDEX_TYPE PolygonVertex
      (const ITYPE2 ipoly, const ITYPE3 j) const
      { return(polygon_vertex_list.Element(ipoly,j)); }

      /// Return index of vertex j of polygon ipoly.
      template <typename ITYPE2, typename ITYPE3>
      VERTEX_INDEX_TYPE VertexIndex
      (const ITYPE2 ipoly, const ITYPE3 j) const
      { return(polygon_vertex_list.Element(ipoly,j)); }

      /// Return reference to polygon_half_edge_list.
      const LIST_OF_LISTS<HALF_EDGE_TYPE,NTYPE> & PolygonHalfEdgeList() const
      { return(polygon_half_edge_list); }

      /// Return reference to list of vertices.
      const std::vector<VERTEX_TYPE> & VertexList() const
      { return(vertex_list); }

      /// Return reference to list of polygons.
      const std::vector<POLYGON_TYPE> & PolygonList() const
      { return(polygon_list); }

      /// Number of (half) edges of polygon i.
      NTYPE NumPolygonEdges(const ITYPE ipoly) const
      { return(polygon_half_edge_list.ListLength(ipoly)); }

      /// Number of vertices of polygon i (= number of polygon edges.)
      NTYPE NumPolygonVertices(const ITYPE ipoly) const
      { return(NumPolygonEdges(ipoly)); }

      /// Return true if polygon is deleted.
      bool IsPolygonDeleted(const ITYPE ipoly) const
      { return(polygon_list[ipoly].is_deleted); }

      /// Return index of j'th half edge of polygon ipoly.
      const HALF_EDGE_INDEX_TYPE HalfEdgeIndex
      (const ITYPE ipoly, const ITYPE j) const
      { return(polygon_half_edge_list.ElementIndex(ipoly,j)); }

      /// Return reference to j'th half edge of polygon ipoly.
      const HALF_EDGE_TYPE & HalfEdge(const ITYPE ipoly, const ITYPE j) const
      { HALF_EDGE_INDEX_TYPE ihalf_edge = HalfEdgeIndex(ipoly,j);
        return(polygon_half_edge_list.element[ihalf_edge]); }

      /// Return index of first half edge of polygon ipoly.
      const HALF_EDGE_INDEX_TYPE FirstHalfEdgeIndex (const ITYPE ipoly) const
      { return(polygon_half_edge_list.FirstElement(ipoly)); }

      /// Return pointer to first half edge of polygon ipoly.
      const HALF_EDGE_TYPE * HalfEdgeList(const ITYPE ipoly) const
      { return(polygon_half_edge_list.List(ipoly)); }

      /// Return index of polygon containing half edge.
      POLYGON_INDEX_TYPE IndexOfPolygonContainingHalfEdge
      (const ITYPE ihalf_edge) const
      { return(polygon_half_edge_list.element[ihalf_edge].polygon_index); }

      /// Return index of polygon containing next half edge around edge.
      POLYGON_INDEX_TYPE IndexOfPolygonContainingNextHalfEdgeAroundEdge
      (const ITYPE ihalf_edge) const
      { const HALF_EDGE_INDEX_TYPE ihalf_edge2 =
          IndexOfNextHalfEdgeAroundEdge(ihalf_edge);
        return(IndexOfPolygonContainingHalfEdge(ihalf_edge2)); }

      /*!
       * @brief Location of half edge in list of half edges in containing polygon.
       * - Returns values in range [0..NumPolygonEdges(ipoly)-1]
       *   where polygon ipoly contains the half edge.
      */
      ITYPE LocationOfHalfEdgeInPolygon(const ITYPE ihalf_edge) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge);
        return(ihalf_edge - FirstHalfEdgeIndex(ipoly));
      }

      /*!
       * @brief Location of next half edge in list of half edges
       *       *in containing polygon.
       * - Returns values in range [0..NumPolygonEdges(ipoly)-1]
       *   where polygon ipoly contains the half edge.
      */
      ITYPE LocationOfNextHalfEdgeInPolygon(const ITYPE ihalf_edge) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge);
        ITYPE iloc = ihalf_edge - FirstHalfEdgeIndex(ipoly);
        iloc = iloc+1;
        if (iloc == NumPolygonEdges(ipoly)) { iloc = 0; }
        return(iloc);
      }

      /*!
       * @brief Location of from vertex in list of vertices in containing polygon.
       * - Returns values in range [0..NumPolygonEdges(ipoly)-1]
       *   where polygon ipoly contains the half edge.
      */
      ITYPE LocationOfFromVertexInPolygon(const ITYPE ihalf_edge) const
      { return(LocationOfHalfEdgeInPolygon(ihalf_edge)); }

      /*!
       * @brief Location of to vertex in list of vertices in containing polygon.
       * - Returns values in range [0..NumPolygonEdges(ipoly)-1]
       *   where polygon ipoly contains the half edge.
      */
      ITYPE LocationOfToVertexInPolygon(const ITYPE ihalf_edge) const
      { return(LocationOfNextHalfEdgeInPolygon(ihalf_edge)); }

      /// Return location in polygon ipoly of iloc relative to jloc.
      template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
      ITYPE LocationRelativeToLocation
      (const ITYPE2 ipoly, const ITYPE3 iloc, const ITYPE4 jloc) const
      {
        if (iloc >= jloc) { return(iloc-jloc); }
        else { return(iloc+NumPolygonVertices(ipoly)-jloc); }
      }

      /*!
       *  @brief Return index of some half edge with from vertex iv.
       *  - Returns boundary half edge if some boundary half edge
       *    has from vertex iv.
       */
      template <typename ITYPE2>
      HALF_EDGE_INDEX_TYPE IndexOfHalfEdgeIncidentOnVertex
      (const ITYPE2 iv) const
      { return(vertex_list[iv].IncidentHalfEdgeIndex()); }

      /// Return true if index of half edge incident on vertex is set.
      template <typename ITYPE2>
      bool IsSetIndexOfHalfEdgeIncidentOnVertex
      (const ITYPE2 iv) const
      { return(vertex_list[iv].IsSetIncidentHalfEdgeIndex()); }

      /// Return index of next half edge around edge.
      template <typename ITYPE2>
      HALF_EDGE_INDEX_TYPE
      IndexOfNextHalfEdgeAroundEdge(const ITYPE2 ihalf_edge) const
      { return((polygon_half_edge_list.element[ihalf_edge].
                IndexOfNextHalfEdgeAroundEdge())); }

      /*!
       *  Return index of previous half edge around edge.
       *  - This routine iterates through the next half edges around edge,
       *    until it finds the half edge preceding the current one.
       *  - Running time is proportional to the number of half edges
       *    around the edge.
       */
      template <typename ITYPE2>
      HALF_EDGE_INDEX_TYPE
      IndexOfPrevHalfEdgeAroundEdge(const ITYPE2 ihalf_edge) const;

      /*!
       *  @brief Return true if half_edge is boundary edge.
       *  @pre AreHalfEdgesLinked() is true.
       */
      template <typename ITYPE2>
      bool IsBoundaryEdge(const ITYPE2 ihalf_edge) const
      { return((ihalf_edge == IndexOfNextHalfEdgeAroundEdge(ihalf_edge))); }

      /*!
       * @brief Return true if polygon contains a boundary edge.
       * @pre AreHalfEdgesLinked() is true.
       */
      template <typename ITYPE2>
      bool IsBoundaryPolygon(const ITYPE2 ipoly) const;

      /*!
       *  @brief Return true if vertex is a from vertex of a boundary half edge.
       */
      template <typename ITYPE2>
      bool IsBoundaryVertex(const ITYPE2 iv) const
      {
        return IsBoundaryEdge(IndexOfHalfEdgeIncidentOnVertex(iv));
      }

      /// Return index of next half edge in polygon.
      HALF_EDGE_INDEX_TYPE IndexOfNextHalfEdgeInPolygon
      (const ITYPE ihalf_edge) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge);
        const ITYPE j =
          (LocationOfHalfEdgeInPolygon(ihalf_edge)+1) % NumPolygonEdges(ipoly);
        return(j + FirstHalfEdgeIndex(ipoly));
      }

      /// Return index of k'th next half edge in polygon.
      template <typename ITYPE2, typename ITYPE3>
      HALF_EDGE_INDEX_TYPE IndexOfKthNextHalfEdgeInPolygon
      (const ITYPE2 ihalf_edge, const ITYPE3 k) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge);
        const ITYPE j =
          (LocationOfHalfEdgeInPolygon(ihalf_edge)+k) % NumPolygonEdges(ipoly);
        return(j + FirstHalfEdgeIndex(ipoly));
      }


      /// Return true if ihalf_edgeB is next half edge of ihalf_edgeA
      template <typename ITYPEA, typename ITYPEB>
      bool IsNextHalfEdge
      (const ITYPEA ihalf_edgeA, const ITYPEB ihalf_edgeB) const
      {
        if (IndexOfNextHalfEdgeInPolygon(ihalf_edgeA) == ihalf_edgeB)
          { return(true); }
        else
          { return(false); }
      }


      /// Return from vertex of half edge.
       VERTEX_INDEX_TYPE FromVertexIndex(const ITYPE ihalf_edge) const
       { return(polygon_vertex_list.element[ihalf_edge]); }

       /// Return index of j'th vertex of polygon ipoly.
       /// - j'th vertex is from_vertex of j'th half edge.
       VERTEX_INDEX_TYPE IndexOfVertexInPolygon
       (const ITYPE ipoly, const ITYPE j) const
       { return(FromVertexIndex(HalfEdgeIndex(ipoly,j))); }

       /// Return to vertex of half edge.
       VERTEX_INDEX_TYPE ToVertexIndex(const ITYPE ihalf_edge) const
       { return(FromVertexIndex(IndexOfNextHalfEdgeInPolygon(ihalf_edge))); }

       /// @brief Return index of next half edge around vertex.
       /// @pre AreHalfEdgesLinked() is true.
       /// @pre IsOrientedManifold() is true.
       template <typename ITYPE2>
       HALF_EDGE_INDEX_TYPE IndexOfNextHalfEdgeAroundVertex
       (const ITYPE2 ihalf_edge0) const
       {
         HALF_EDGE_INDEX_TYPE ihalf_edge1 =
           IndexOfPrevHalfEdgeInPolygon(ihalf_edge0);
         ihalf_edge1 = IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
         return(ihalf_edge1);
       }

       /// @brief Return true if iv is endpoint of half edge.
       template <typename ITYPE2, typename ITYPE3>
       bool IsHalfEdgeEndpoint
       (const ITYPE2 iv, const ITYPE3 ihalf_edge) const
       {
         return ((FromVertexIndex(ihalf_edge) == iv) ||
           (ToVertexIndex(ihalf_edge) == iv));
       }

       /// Return true if polygon is a triangle, i.e. has 3 vertices.
       template <typename ITYPE2>
       bool IsTriangle(const ITYPE2 ipoly) const
       { return((NumPolygonEdges(ipoly) == 3)); }

       /// Return true if polygon is a quadrilateral, i.e. has 4 vertices.
       template <typename ITYPE2>
       bool IsQuadrilateral(const ITYPE2 ipoly) const
       { return((NumPolygonEdges(ipoly) == 4)); }

       /// Return true if polygon is a pentagon, i.e. has 5 vertices.
       template <typename ITYPE2>
       bool IsPentagon(const ITYPE2 ipoly) const
       { return((NumPolygonEdges(ipoly) == 5)); }

      /// Get index of adjacent half edge in polygon (next or previous.)
      template <typename ITYPE0>
      HALF_EDGE_INDEX_TYPE GetIndexOfAdjacentHalfEdgeInPolygon
      (const ITYPE0 ihalf_edge, const bool flag_next) const
      {
        if (flag_next)
          { return(IndexOfNextHalfEdgeInPolygon(ihalf_edge)); }
        else
          { return(IndexOfPrevHalfEdgeInPolygon(ihalf_edge)); }
      }


      /// Get index of endpoint of half edge (to or from.)
      template <typename ITYPE0>
      VERTEX_INDEX_TYPE GetIndexOfHalfEdgeEndpoint
      (const ITYPE0 ihalf_edge, const bool flag_to) const
      {
        if (flag_to)
          { return(ToVertexIndex(ihalf_edge)); }
        else
          { return(FromVertexIndex(ihalf_edge)); }
      }


      /// Get location of endpoint of half edge (to or from) in polygon.
      template <typename ITYPE0>
      VERTEX_INDEX_TYPE GetLocationOfHalfEdgeEndpointInPolygon
      (const ITYPE0 ihalf_edge, const bool flag_to) const
      {
        if (flag_to)
          { return(LocationOfHalfEdgeInPolygon
                   (IndexOfNextHalfEdgeInPolygon(ihalf_edge))); }
        else
          { return(LocationOfHalfEdgeInPolygon(ihalf_edge)); }
      }

      /// Get indices of next two half edges in polygon.
      template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
      void GetIndicesOfNextTwoHalfEdgesInPolygon
      (const ITYPE0 ihalf_edge0, ITYPE1 & ihalf_edge1,
       ITYPE2 & ihalf_edge2) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge0);
        const HALF_EDGE_INDEX_TYPE ifirst_half_edge = FirstHalfEdgeIndex(ipoly);
        const ITYPE iloc0 = LocationOfHalfEdgeInPolygon(ihalf_edge0);
        const NTYPE nume = NumPolygonEdges(ipoly);

        const ITYPE iloc1 = (iloc0 + 1) % nume;
        const ITYPE iloc2 = (iloc1 + 1) % nume;
        ihalf_edge1 = ifirst_half_edge + iloc1;
        ihalf_edge2 = ifirst_half_edge + iloc2;
      }


      /// Get indices of next three half edges in polygon.
      template <typename ITYPE0, typename ITYPE1, typename ITYPE2,
                typename ITYPE3>
      void GetIndicesOfNextThreeHalfEdgesInPolygon
      (const ITYPE0 ihalf_edge0, ITYPE1 & ihalf_edge1,
       ITYPE2 & ihalf_edge2, ITYPE3 & ihalf_edge3) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge0);
        const HALF_EDGE_INDEX_TYPE ifirst_half_edge = FirstHalfEdgeIndex(ipoly);
        const ITYPE iloc0 = LocationOfHalfEdgeInPolygon(ihalf_edge0);
        const NTYPE nume = NumPolygonEdges(ipoly);

        const ITYPE iloc1 = (iloc0 + 1) % nume;
        const ITYPE iloc2 = (iloc1 + 1) % nume;
        const ITYPE iloc3 = (iloc2 + 1) % nume;
        ihalf_edge1 = ifirst_half_edge + iloc1;
        ihalf_edge2 = ifirst_half_edge + iloc2;
        ihalf_edge3 = ifirst_half_edge + iloc3;
      }


      /// Get indices of next three half edges in polygon.
      template <typename ITYPE0, typename ITYPE1, typename ITYPE2,
                typename ITYPE3, typename ITYPE4>
      void GetIndicesOfNextFourHalfEdgesInPolygon
      (const ITYPE0 ihalf_edge0, ITYPE1 & ihalf_edge1,
       ITYPE2 & ihalf_edge2, ITYPE3 & ihalf_edge3,
       ITYPE4 & ihalf_edge4) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge0);
        const HALF_EDGE_INDEX_TYPE ifirst_half_edge = FirstHalfEdgeIndex(ipoly);
        const ITYPE iloc0 = LocationOfHalfEdgeInPolygon(ihalf_edge0);
        const NTYPE nume = NumPolygonEdges(ipoly);

        const ITYPE iloc1 = (iloc0 + 1) % nume;
        const ITYPE iloc2 = (iloc1 + 1) % nume;
        const ITYPE iloc3 = (iloc2 + 1) % nume;
        const ITYPE iloc4 = (iloc3 + 1) % nume;
        ihalf_edge1 = ifirst_half_edge + iloc1;
        ihalf_edge2 = ifirst_half_edge + iloc2;
        ihalf_edge3 = ifirst_half_edge + iloc3;
        ihalf_edge4 = ifirst_half_edge + iloc4;
      }


      /// Return index of previous half edge around polygon.
      HALF_EDGE_INDEX_TYPE IndexOfPrevHalfEdgeInPolygon
      (const ITYPE ihalf_edge) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge);
        const ITYPE j =
          (LocationOfHalfEdgeInPolygon(ihalf_edge)+(NumPolygonEdges(ipoly)-1)) %
          NumPolygonEdges(ipoly);
        return(j + FirstHalfEdgeIndex(ipoly));
      }


      /// Get indices of previous two half edges in polygon.
      template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
      void GetIndicesOfPrevTwoHalfEdgesInPolygon
      (const ITYPE0 ihalf_edge0, ITYPE1 & ihalf_edge1,
       ITYPE2 & ihalf_edge2) const
      {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge0);
        const HALF_EDGE_INDEX_TYPE ifirst_half_edge = FirstHalfEdgeIndex(ipoly);
        const ITYPE iloc0 = LocationOfHalfEdgeInPolygon(ihalf_edge0);
        const NTYPE nume = NumPolygonEdges(ipoly);

        const ITYPE iloc1 = (iloc0 + nume-1) % nume;
        const ITYPE iloc2 = (iloc1 + nume-1) % nume;
        ihalf_edge1 = ifirst_half_edge + iloc1;
        ihalf_edge2 = ifirst_half_edge + iloc2;
      }

      /// @brief Update half edges incident on polygon vertices.
      /// - Only set half edge incident on vertex iv, if the half edge incident
      ///   on vertex iv is not set.
      /// - When altering mesh, should be called on new polygons after all
      ///   old polygons have been deleted by call to _DeletePolygon().
      template <typename ITYPE2>
      void UpdateHalfEdgesIncidentOnPolygonVertices
      (const ITYPE2 ipoly);

      /// Return true if polygon ipoly contains vertex iv.
      template <typename ITYPE2, typename ITYPE3>
      bool DoesPolygonContainVertex(const ITYPE2 ipoly, const ITYPE3 iv) const
      { return(polygon_vertex_list.DoesListContain(ipoly, iv)); }

      /// @brief Return true if polygon ipoly contains edge (iv0,iv1).
      /// @pre Half edges are listed in cylic order around the polygon.
      template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
      bool DoesPolygonContainEdge
      (const ITYPE2 ipoly, const ITYPE3 iv0, const ITYPE4 iv1) const;

      /*!
       *  @brief Return true if polygon ipoly share two adjacent edges with some other polygon.
       *  @pre Polygon contains at least one half edge.
       *  @param ipoly Polygon index.
       *  @param[out] iloc Location of second shared half_edge.
       *    - HalfEdgeIndex(ipoly,iloc) and the previous half edge
       *      are shared with the other polygon.
       */
      template <typename ITYPE2, typename ITYPE3>
      bool DoesPolygonShareTwoAdjacentEdges
      (const ITYPE2 ipoly, ITYPE3 & iloc) const;

      /*!
       *  @brief Return true if two polygons, ipolyA and ipolyB, share three or more vertices.
       *  @pre Polygons ipolyA and ipolyB share an edge.
       *  @param ihalf_edgeA Half edge of polygon ipolyA corresponding
       *     to edge shared by ipolyA and ipolyB.
       */
      template <typename ITYPE2>
      bool DoPolygonsShareThreeVertices(const ITYPE2 ihalf_edgeA) const;

      /*!
       *  @brief Get triangle vertex and opposing half_edge.
       *  @pre Polygon itriangle is a triangle.
       *  @param itriangle Polygon index.
       *  @param iloc Location (0, 1 or 2) of vertex in triangle.
       *  @param[out] iv Vertex at location iloc.
       *  @param[out] ihalf_edge0 Half-edge opposite vertex iv.
       */
      template <typename ITYPE2, typename ITYPE3,
                typename ITYPE4, typename ITYPE5>
      void GetTriangleVertexAndOpposingEdge
      (const ITYPE2 itriangle, const ITYPE3 iloc,
       ITYPE4 & iv, ITYPE5 & ihalf_edge0);

      // *** CHECK THAT RETURNS false if ihalf_edgeA is a boundary half edge.
      /// @brief Return true if there are exactly two half edges in the cap containing half edge ihalf_edgeA.
      /// @pre Assumes half edges are linked.
      template <typename ITYPE2>
      bool AreTwoHalfEdgesAroundFromVertex
      (const ITYPE2 ihalf_edgeA) const;

      /// @brief Return true if there are exactly three half edges in the cap containing half edge ihalf_edgeA.
      /// @pre Assumes half edges are linked.
      template <typename ITYPE2>
      bool AreThreeHalfEdgesAroundFromVertex
      (const ITYPE2 ihalf_edgeA) const;

      /// @brief Return true if there are four or more half edges in the cap containing half edge ihalf_edgeA.
      /// @pre Assumes half edges are linked.
      template <typename ITYPE2>
      bool AreFourOrMoreHalfEdgesAroundFromVertex
      (const ITYPE2 ihalf_edgeA) const;

      /// @brief Return true if there are exactly two half edges in the cap containing vertex iv.
      /// - Note: Argument is a vertex index, not a half edge index.
      template <typename ITYPE2>
      bool AreTwoHalfEdgesAroundVertex(const ITYPE2 iv) const;

      /// @brief Return true if there are exactly three half edges in the cap containing vertex iv.
      /// - Note: Argument is a vertex index, not a half edge index.
      template <typename ITYPE2>
      bool AreThreeHalfEdgesAroundVertex(const ITYPE2 iv) const;

      /*!
       *  @brief Append triangle vertices and triangle indices to lists.
       *  - Append vertices of triangles in mesh to triangle_vertex_list[].
       *  - Append corresponding triangle indices to triangle_index[].
       *  - Ignores/skips any polygons that are not triangles.
       *  - Ignores/skips any deleted triangles.
       *  @param[out] triangle_vertex_list[] List of triangle vertices.
       *    - i'th triangle in list has vertices:
       *      (triangle_vertex_list[3*i],triangle_vertex_list[3*i+1],
       *       triangle_vertex_list[3*i+2]).
       *  @param[out] triangle_index_list[] List of triangle indices.
       *    - i'th triangle in list is triangle_index[i].
       */
      template <typename VTYPE2, typename PTYPE2>
      void AppendToListOfVerticesOfAllTriangles
      (std::vector<VTYPE2> & triangle_vertex_list,
       std::vector<PTYPE2> & triangle_index) const;

      /*!
       *  @brief Get triangle vertices and corresponding triangle indices.
       *  - Similar to AppendToListOfAllTriangles() but first clears
       *  - Ignores/skips any polygons that are not triangles.
       *  - Ignores/skips any deleted triangles.
       */
      template <typename VTYPE2, typename PTYPE2>
      void GetListOfVerticesOfAllTriangles
      (std::vector<VTYPE2> & triangle_vertex_list,
       std::vector<PTYPE2> & triangle_index) const;

      /*!
       *  @brief Sort triangle vertices.
       *  @param[out] triangle_vertex_list[] Array of triangle vertices.
       *    - Routine sorts vertices of each triangle.
       *    @pre triangle_vertex_list.size()%3 = 0.
       */
      template <typename VTYPE2>
      void SortTriangleVertices
      (std::vector<VTYPE2> & triangle_vertex_list) const;


      /*!
       *  @brief Get sorted list of sorted triangle vertices and corresponding triangle indices.
       *  - Linear time sorting (assuming bounded vertex degree) using LIST_OF_LISTS.
       *  @param[out] sorted_triangle_vertex_list[] Sorted list of sorted triangle vertices.
       *    - i'th triangle in list has vertices:
       *      (sorted_triangle_vertex_list[3*i],
       *       sorted_triangle_vertex_list[3*i+1],
       *       sorted_triangle_vertex_list[3*i+2]).
       *    - sorted_triangle_vertex_list[3*i] <= sorted_triangle_vertex_list[3*i+1].
       *    - sorted_triangle_vertex_list[3*i+1] <= sorted_triangle_vertex_list[3*i+2].
       *    - Vertex tuple of i'th triangle is lexicographically less than or equal to
       *      vertex tuple of triangle i+1.
       *  @param[out] triangle_index_list[] List of triangle indices
       *    corresponding to triangles in sorted_triangle_vertex_list[].
       *    - i'th triangle in list is triangle_index[i].
       */
      template <typename VTYPE2, typename PTYPE2>
      void GetSortedListOfVerticesOfAllTriangles
      (std::vector<VTYPE2> & sorted_triangle_vertex_list,
       std::vector<PTYPE2> & triangle_index) const;

      /*!
       *  @brief Get list of vertices of duplicate triangles.
       *  @param[out] duplicate_triangle_vertex[]
       *    Array containing vertices of duplicate triangles.
       *    - Duplicate triangle i has vertices:
       *      (duplicate_triangle_vertex[3*i],
       *       duplicate_triangle_vertex[3*i+1],
       *       duplicate_triangle_vertex[3*i+2])
       *    - Note: Vertex order is arbitrary.
       *  @param[out] indices_of_duplicate_triangle[i]
       *    List of indices of triangles that have vertices
       *    given by i'th triangle represented by duplicate_triangle_vertex[].
       */
      template <typename VTYPE2, typename NTYPE2>
      void GetListOfVerticesOfDuplicateTriangles
      (std::vector<VTYPE2> & duplicate_triangle_vertex,
       IJK::LIST_OF_LISTS<NTYPE2,NTYPE2> &
         indices_of_duplicate_triangles) const;

      // Add vertices/polygons functions.

      /// Add a single vertex.
      VERTEX_INDEX_TYPE AddVertex();

      /// Add k vertices
      template <typename NTYPE2>
      void AddVertices(const NTYPE2 k);

      /// @brief Add a polygon with list_length vertices.
      /// - Sets flag are_half_edges_linked to false.
      /// - Sets flag is_oriented_manifold to false.
      template <typename ITYPE2, typename NTYPE2>
      POLYGON_INDEX_TYPE AddPolygon
      (const ITYPE2 poly_vert_list[], const NTYPE2 list_length);

      /// @brief Add a polygon.
      /// - C++ STL vector type for array poly_vert_list[].
      /// @param poly_vert_list List of polygon vertices in C++ STL vector.
      template <typename ITYPE2>
      POLYGON_INDEX_TYPE AddPolygon(const std::vector<ITYPE2> & poly_vert_list)
      { return(AddPolygon(IJK::vector2pointer(poly_vert_list),
                          poly_vert_list.size())); }

      /// Add polygons.
      template <typename NTYPE2, typename NTYPE3,
                typename ITYPE2, typename ITYPE3>
      void AddPolygons(const NTYPE2 * num_poly_vert,
                       const NTYPE3 num_poly,
                       const ITYPE2 * poly_vert_list,
                       const ITYPE3 * first_poly_vert);

      /// Add polygons.
      template <typename NTYPE2, typename ITYPE2, typename ITYPE3>
      void AddPolygons(const std::vector<NTYPE2> & num_poly_vert,
                       const std::vector<ITYPE2> & poly_vert_list,
                       const std::vector<ITYPE3> & first_poly_vert)
      { AddPolygons(IJK::vector2pointer(num_poly_vert),
                    num_poly_vert.size(),
                    IJK::vector2pointer(poly_vert_list),
                    IJK::vector2pointer(first_poly_vert)); }

      /// @brief Add polygon from list where each polygon has num_vert_per_poly vertices.
      /// @param num_vert Total number of vertices in poly_vert_list[].
      /// @param num_vert_per_poly Number of vertices in each polygon.
      template <typename ITYPE2, typename NTYPE2, typename NTYPE3>
      void AddPolygons(const ITYPE2 poly_vert_list[],
                       const NTYPE2 num_vert,
                       const NTYPE3 num_vert_per_poly);

      /// @brief Add polygons from list where each polygon has num_vert_per_poly vertices.
      /// - C++ STL vector type for array poly_vert_list.
      /// @param poly_vert_list List of polygon vertices.
      /// @pre poly_vert_list.size() is not 0.
      /// @param num_vert_per_poly Number of vertices in each polygon.
      template <typename ITYPE2>
      void AddPolygons(const std::vector<ITYPE2> & poly_vert_list,
                       const NTYPE num_vert_per_poly)
      { AddPolygons(IJK::vector2pointer(poly_vert_list), poly_vert_list.size(),
                    num_vert_per_poly); }


      /// @brief Compute number of polygons incident on each vertex.
      template <typename NTYPE2>
      void ComputeNumberOfPolygonsIncidentOnEachVertex
      (std::vector<NTYPE2> & num_incident_polygons) const;

      /// @brief Return number of boundary edges.
      NTYPE CountNumberOfBoundaryEdges() const;
      
      /*!
       *  @brief Link all half edges in data structure.
       *  - May change half_edge_incident_index for each vertex.
      */
      void LinkHalfEdges();

      /// Copy.
      template <typename MTYPE>
      void Copy(const MTYPE & meshA);

      /*!
       * @brief Determine if mesh is an oriented manifold.
       * - Return true if mesh is an oriented manifold.
       * - Sets flag is_oriented_manifold.
       * @pre AreHalfEdgesLinked() is true.
       * - (Call LinkHalfEdges() before calling this function.)
      */
      bool DetermineOrientedManifold();

      /// @brief Determine if each mesh edge is in at most two polygons.
      /// - Return true if each mesh edge is in at most two polygons.
      template <typename ITYPE2>
      bool DetermineManifoldEdges(ITYPE2 & half_edge_index) const;

      /*!
       * @brief Determine if polygons have matching orientations.
       * - Return true if polygons have matching orientations.
       * @param[out] is_non_manifold If true, some non-manifold condition
       *             was found.
       * @param[out] jpoly,kpoly Polygons with mismatched orientations.
       *             Set only if false is returned and is_non_manifold is false;
      */
      template <typename ITYPE2, typename ITYPE3>
      bool DetermineMatchingPolygonOrientations
      (bool & is_non_manifold, ITYPE2 & jpoly, ITYPE3 & kpoly) const;

      /*!
       * @brief Determine if each mesh vertex is in a disk of polygons.
       * - Return true if each mesh vertex is in a disk of polygons.
       * @param[out] is_oriented If not true, mesh is not oriented.
       * @param[out] vertex_index Index of non-manifold vertex.
       *             Set only if false is returned and is_oriented is true;
      */
      template <typename ITYPE2>
      bool DetermineOrientedManifoldVertices
      (bool & is_oriented, ITYPE2 & vertex_index) const;


      // Triangulate functions.

      /// @brief Triangulate quadrilateral ipoly from vertex 0.
      /// @pre Polygon ipoly is a quadrilateral.
      template <typename ITYPE2>
      void TriangulateQuadFromVertex0(const ITYPE2 ipoly);

      /// Triangulate polygon ipoly from vertex 0.
      template <typename ITYPE2>
      void TriangulatePolygonFromVertex0(const ITYPE2 ipoly);

      /*!
       *  @brief Uniformly triangulate mesh polygons from vertex 0.
       *  - First calls TriangulatePolygonSharingTwoAdjacentEdges().
      */
      void TriangulateUniform();

      /*!
       *  @brief Triangulate any polygons sharing two adjacent edges.
       *  - Triangulates polygons from vertex between the two shared edges.
       */
      void TriangulatePolygonsSharingTwoAdjacentEdges();

      /// Triangulate polygon ipoly from j'th vertex.
      template <typename ITYPE2, typename ITYPE3>
      void TriangulatePolygonFromVertex(const ITYPE2 ipoly, ITYPE3 j);

      /// Triangulate polygon ipoly from interior vertex iv.
      template <typename ITYPE2, typename ITYPE3>
      void TriangulatePolygonFromInteriorVertex
      (const ITYPE2 ipoly, ITYPE3 iv);


      // Print routines (mainly for debugging.)

      /*!
       *  @brief Print half edge. (Extended version.)
       *  - Mainly for debugging.
       *  - Extended version.  Define left/right delimiters and separator.
       *  @param c0 Left delimiter.
       *  @param c1 Separator.
       *  @param c2 Right delimiter.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2>
      void PrintHalfEdgeX
      (OSTREAM_TYPE & out, const ITYPE2 ihalf_edge,
       const char c0, const char c1, const char c2) const;

      /*!
       *  @brief Print half edge. (Default delimiters and separator.)
       *  - Mainly for debugging.
       *  - Version using delimiters '(' and ')' and separator ','.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2>
      void PrintHalfEdge(OSTREAM_TYPE & out, const ITYPE2 ihalf_edge) const
      { PrintHalfEdgeX(out, ihalf_edge, '(', ',', ')'); }

      /*!
       *  @brief Print half edge. (Set prefix and suffix.)
       *  - Mainly for debugging.
       *  - Version adding prefix and suffix strings.
       *  @param s0 Prefix string.
       *  @param s1 Suffix string.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2,
                typename STYPE0, typename STYPE1>
      void PrintHalfEdge(OSTREAM_TYPE & out, const STYPE0 & s0,
                          const ITYPE2 ihalf_edge, const STYPE1 & s1) const;

      /*!
       *  @brief Print half edge index and endpoints. (Extended version.)
       *  - Mainly for debugging.
       *  - Extended version.  Define left/right delimiters and separator.
       *  @param c0 Left delimiter.
       *  @param c1 Separator.
       *  @param c2 Right delimiter.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2>
      void PrintHalfEdgeIndexAndEndpointsX
      (OSTREAM_TYPE & out, const ITYPE2 ihalf_edge,
       const char c0, const char c1, const char c2) const;

      /*!
       *  @brief Print half edge index and endpoints.
       *  (Default delimiters and separator.)
       *  - Mainly for debugging.
       *  - Version using delimiters '(' and ')' and separator ','.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2>
      void PrintHalfEdgeIndexAndEndpoints
      (OSTREAM_TYPE & out, const ITYPE2 ihalf_edge) const
      { PrintHalfEdgeIndexAndEndpointsX(out, ihalf_edge, '(', ',', ')'); }

      /*!
       *  @brief Print half edge index and endpoints. (Set prefix and suffix.)
       *  - Mainly for debugging.
       *  - Version adding prefix and suffix strings.
       *  @param s0 Prefix string.
       *  @param s1 Suffix string.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2,
                typename STYPE0, typename STYPE1>
      void PrintHalfEdgeIndexAndEndpoints
      (OSTREAM_TYPE & out, const STYPE0 & s0,
       const ITYPE2 ihalf_edge, const STYPE1 & s1) const;

      /*!
       *  @brief Print half edge information.
       *  - Include endpoints, containing polygons,
       *    next and previous half edges, and next half edge around edge.
       */
      template <typename OSTREAM_TYPE, typename STYPE0, typename ITYPE2>
      void PrintHalfEdgeInfo
      (OSTREAM_TYPE & out, const STYPE0 & prefix,
       const ITYPE2 ihalf_edge) const;

      /*!
       *  @brief Print information of all half edges.
       *  - Include endpoints, containing polygons,
       *    next and previous half edges, and next half edge around edge.
       */
      template <typename OSTREAM_TYPE>
      void PrintInfoOfAllHalfEdges(OSTREAM_TYPE & out) const;

      /*!
       *  @brief Print polygon vertices.
       *  - Extended version.  Define left/right delimiters and separator.
       *  - Mainly for debugging.
       *  @param c0 Left delimiter.
       *  @param c1 Separator.
       *  @param c2 Right delimiter.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2>
      void PrintPolygonVerticesX
      (OSTREAM_TYPE & out, const ITYPE2 ipoly,
       const char c0, const char c1, const char c2) const;

      /*!
       *  @brief Print polygon vertices.
       *    (Default delimiters and separator.)
       *  - Mainly for debugging.
       *  - Version using delimiters '(' and ')' and separator ','.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2>
      void PrintPolygonVertices
      (OSTREAM_TYPE & out, const ITYPE2 ipoly) const
      { PrintPolygonVerticesX(out, ipoly, '(', ',', ')'); }

      /*!
       *  @brief Print polygon vertices. (Set prefix and suffix.)
       *  - Mainly for debugging.
       *  - Version adding prefix and suffix strings.
       *  @param s0 Prefix string.
       *  @param s1 Suffix string.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2,
                typename STYPE0, typename STYPE1>
      void PrintPolygonVertices
      (OSTREAM_TYPE & out, const STYPE0 & s0,
       const ITYPE2 ipoly, const STYPE1 & s1) const;

      /*!
       *  @brief Print polygon index and vertices. (Extended version.)
       *  - Extended version.  Define left/right delimiters and separator.
       *  - Mainly for debugging.
       *  @param c0 Left delimiter.
       *  @param c1 Separator.
       *  @param c2 Right delimiter.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2>
      void PrintPolygonIndexAndVerticesX
      (OSTREAM_TYPE & out, const ITYPE2 ipoly,
       const char c0, const char c1, const char c2) const;

      /*!
       *  @brief Print polygon index and vertices.
       *         (Default delimiters and separator.)
       *  - Mainly for debugging.
       *  - Version using delimiters '(' and ')' and separator ','.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2>
      void PrintPolygonIndexAndVertices
      (OSTREAM_TYPE & out, const ITYPE2 ipoly) const
      { PrintPolygonIndexAndVerticesX(out, ipoly, '(', ',', ')'); }

      /*!
       *  @brief Print polygon index and vertices. (Set prefix and suffix.)
       *  - Mainly for debugging.
       *  - Version adding prefix and suffix strings.
       *  @param s0 Prefix string.
       *  @param s1 Suffix string.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2,
                typename STYPE0, typename STYPE1>
      void PrintPolygonIndexAndVertices
      (OSTREAM_TYPE & out, const STYPE0 & s0,
       const ITYPE2 ipoly, const STYPE1 & s1) const;

      /*!
       *  @brief Print index and vertices of polygon containing half edge. (Set prefix and suffix.)
       *  - Mainly for debugging
       */
      template <typename OSTREAM_TYPE, typename ITYPE2,
                typename STYPE0, typename STYPE1>
      void PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (OSTREAM_TYPE & out, const STYPE0 & s0,
       const ITYPE2 ihalf_edge, const STYPE1 & s1) const
      {
        PrintPolygonIndexAndVertices
          (out, s0, this->IndexOfPolygonContainingHalfEdge(ihalf_edge), s1);
      }

      /*!
       *  @brief Print polygon vertices and half edges.  
       *  - Mainly for debugging.
       *  @param prefix Prefix string.
       */
      template <typename OSTREAM_TYPE, typename ITYPE2,
                typename STYPE0>
      void PrintPolygonInfo
      (OSTREAM_TYPE & out, const STYPE0 & prefix,
       const ITYPE2 ipoly) const;

      /*!
       *  @brief Print information of all polygons.
       */
      template <typename OSTREAM_TYPE>
      void PrintInfoOfAllPolygons(OSTREAM_TYPE & out) const;


      // Check functions.

      /// @brief Check data structure.
      /// - Return true if data structure passes check.
      bool Check(IJK::ERROR & error) const;

      /// @brief Check number of polygons.
      /// - Return true if data structure passes check.
      bool CheckNumPolygons(IJK::ERROR & error) const;

      /// @brief Check that iv is a vertex index.
      /// - Return true if data structure passes check.
      template <typename ITYPE2>
      bool CheckVertex(const ITYPE2 iv, IJK::ERROR & error) const;

      /// @brief Check polygon_half_edge_list and polygon_vertex_list match.
      /// - Check same number of lists and corresponding lists have
      ///   the same length.
      /// - Return true if data structure passes check.
      bool CheckHalfEdges(IJK::ERROR & error) const;

      /// @brief Check that some half edge is incident on from vertex of ihalf_edge.
      /// - Return true if data structure passes check.
      template <typename ITYPE2>
      bool CheckHalfEdgeIncidentOnVertex
      (const ITYPE2 ihalf_edge, IJK::ERROR & error) const;

      /// @brief Check half edges incident on vertices.
      /// - Return true if data structure passes check.
      bool CheckHalfEdgesIncidentOnVertices(IJK::ERROR & error) const;

      /// Check polygons containing half edges.
      /// - Check that the polygon index is correct in each half edge.
      bool CheckPolygonsContainingHalfEdges(IJK::ERROR & error) const;

      /*!
       *  @brief Check that iloc is a valid location in polygon ipoly.
       *  - Return true if data structure passes check.
       *  @param ipoly Polygon index.
       *    @pre ipoly is a valid polygon index.
       *  @param iloc Polygon location.
       *  @param error Error class.
      */
      template <typename PTYPE2, typename ITYPE2>
      bool CheckPolygonLocation
      (const PTYPE2 ipoly, const ITYPE2 iloc, IJK::ERROR & error) const;

      /// @brief Check half edge links.
      /// - Return true if data structure passes check.
      bool CheckHalfEdgeLinks(IJK::ERROR & error) const;

      /// @brief Check that AreHalfEdgesLinked() is true.
      /// - Return true if AreHalfEdgesLinked() is true.
      bool CheckAreHalfEdgesLinked(IJK::ERROR & error) const;

      /// @brief Check that data structure is an oriented manifold.
      /// - Return true if data structure passes check.
      bool CheckOrientedManifold(IJK::ERROR & error) const;

      /// @brief Check that two half edges are in the same polygon.
      /// - Return true if two half edges are in the same polygon.
      template <typename ITYPEH>
      bool CheckAreHalfEdgesInSamePolygon
      (const ITYPEH ihalf_edge0, const ITYPEH ihalf_edge1,
       IJK::ERROR & error) const;

      /// @brief Check that three half edges are in the same polygon.
      /// - Return true if three half edges are in the same polygon.
      template <typename ITYPEH>
      bool CheckAreHalfEdgesInSamePolygon
      (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
       const ITYPEH ihalf_edgeC, IJK::ERROR & error) const;

      /// @brief Check that two half edges are different.
      /// - Return true if two half edges are different
      template <typename ITYPEH>
      bool CheckAreHalfEdgesDifferent
      (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
       IJK::ERROR & error) const;

      /// @brief Check that three half edges are different.
      /// - Return true if three half edges are different
      template <typename ITYPEH>
      bool CheckAreHalfEdgesDifferent
      (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
       const ITYPEH ihalf_edgeC, IJK::ERROR & error) const;

    };


    // *****************************************************************
    // Class MESH2D_A
    // *****************************************************************

    /*!
     *  @brief Mesh of 2D polygons.
     *         (Uses defaults for vertex, half_edge and polygon types.)
     *  @details
     *  - Use VERTEX_BASE, HALF_EDGE_BASE, POLYGON_BASE
     *    to represent vertices, half edges and polygons.
    */
    template <typename ITYPE, typename NTYPE>
    class MESH2D_A:
      public MESH2D_BASE<VERTEX_BASE<ITYPE,ITYPE,NTYPE>,
                         HALF_EDGE_BASE<ITYPE,ITYPE,ITYPE>,
                         POLYGON_BASE, ITYPE, NTYPE>
    {
    public:
      MESH2D_A() {};
    };


    // *****************************************************************
    // Class MESH2D_BASE member functions
    // *****************************************************************

    // constructor
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::MESH2D_BASE()
    {
      are_half_edges_linked = false;
      is_oriented_manifold = false;
    }


    // Return true if polygon ipoly contains edge (iv0,iv1).
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::DoesPolygonContainEdge
    (const ITYPE2 ipoly, const ITYPE3 iv0, const ITYPE4 iv1) const
    {
      if (NumPolygonVertices(ipoly) <= 1) { return(false); }

      for (ITYPE j0 = 0; j0 < NumPolygonVertices(ipoly); j0++) {
        VERTEX_INDEX_TYPE jv0 = IndexOfVertexInPolygon(ipoly, j0);
        ITYPE j1 = (j0+1)%NumPolygonVertices(ipoly);
        VERTEX_INDEX_TYPE jv1 = IndexOfVertexInPolygon(ipoly, j1);

        if (jv0 == iv0 && jv1 == iv1) { return(true); }
        if (jv1 == iv0 && jv0 == iv1) { return(true); }
      }

      return(false);
    }


    // Return true if polygon ipoly share two adjacent half edges
    //   with some other polygon.
    // @param[out] iloc Location of second shared half_edge.
    //   - HalfEdgeIndex(ipoly,iloc) and the previous half edge
    //     are shared with the other polygon.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    DoesPolygonShareTwoAdjacentEdges
    (const ITYPE2 ipoly, ITYPE3 & iloc) const
    {
      iloc = 0;
      for (NUMBER_TYPE i = 0; i < NumPolygonEdges(ipoly); i++) {
        const HALF_EDGE_INDEX_TYPE ihalf_edge =
          HalfEdgeIndex(ipoly, i);

        if (AreTwoHalfEdgesAroundFromVertex(ihalf_edge)) {
          iloc = i;
          return(true);
        }
      }

      return(false);
    }


    // Return true if two polygons, ipolyA and ipolyB, share three
    //   or more vertices.
    // @pre Polygons ipolyA and ipolyB share an edge.
    // @param ihalf_edgeA Half edge of polygon ipolyA corresponding
    //    to edge shared by ipolyA and ipolyB
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    DoPolygonsShareThreeVertices(const ITYPE2 ihalf_edgeA) const
    {
      if (IsBoundaryEdge(ihalf_edgeA)) { return(false); }

      const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
        IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA);
      const POLYGON_INDEX_TYPE ipolyA =
        IndexOfPolygonContainingHalfEdge(ihalf_edgeA);
      const POLYGON_INDEX_TYPE ipolyB =
        IndexOfPolygonContainingHalfEdge(ihalf_edgeB);
      const NTYPE numvA = NumPolygonVertices(ipolyA);
      const NTYPE numvB = NumPolygonVertices(ipolyB);

      ITYPE ilocA = LocationOfNextHalfEdgeInPolygon(ihalf_edgeA);
      for (ITYPE jA = 2; jA < numvA; jA++) {
        ilocA = (ilocA+1)%numvA;
        const VERTEX_INDEX_TYPE ivA = PolygonVertex(ipolyA, ilocA);

        ITYPE ilocB = LocationOfNextHalfEdgeInPolygon(ihalf_edgeB);
        for (ITYPE jB = 2; jB < numvB; jB++) {
          ilocB = (ilocB+1)%numvB;
          const VERTEX_INDEX_TYPE ivB = PolygonVertex(ipolyB, ilocB);

          if (ivA == ivB) { return(true); }
        }
      }

      return(false);
    }


    // Get triangle vertex and opposing half_edge.
    // @pre Polygon itriangle is a triangle.
    // @param itriangle Polygon index.
    // @param iloc Location (0, 1 or 2) of vertex in triangle.
    // @param[out] iv Vertex at location iloc.
    // @param[out] ihalf_edge0 Half-edge opposite vertex iv.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3,
              typename ITYPE4, typename ITYPE5>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    GetTriangleVertexAndOpposingEdge
    (const ITYPE2 itriangle, const ITYPE3 iloc,
     ITYPE4 & iv, ITYPE5 & ihalf_edge0)
    {
      const NTYPE NUM_VERTICES_PER_TRIANGLE(3);

      iv = IndexOfVertexInPolygon(itriangle, iloc);

      // Assumes precondition that itriagnle is a triangle.
      const ITYPE3 iloc1 = (iloc+1) % NUM_VERTICES_PER_TRIANGLE;

      ihalf_edge0 = HalfEdgeIndex(itriangle, iloc1);
    }


    // Return true if there are exactly two half edges
    //   in the cap containing half edge ihalf_edgeA.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    AreTwoHalfEdgesAroundFromVertex
    (const ITYPE2 ihalf_edgeA) const
    {
      const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
        IndexOfNextHalfEdgeAroundVertex(ihalf_edgeA);

      if (ihalf_edgeB == ihalf_edgeA) {
        // Half edge ihalf_edgeA is only half edge in cap.
        return(false);
      }

      if (IndexOfNextHalfEdgeAroundVertex(ihalf_edgeB)
          == ihalf_edgeA)
        { return(true); }

      return(false);
    }


    // Return true if there are exactly three half edges
    //   in the cap containing half edge ihalf_edgeA.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    AreThreeHalfEdgesAroundFromVertex
    (const ITYPE2 ihalf_edgeA) const
    {
      const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
        IndexOfNextHalfEdgeAroundVertex(ihalf_edgeA);

      if (ihalf_edgeB == ihalf_edgeA) {
        // Half edge ihalf_edgeA is only half edge in cap.
        return(false);
      }

      const HALF_EDGE_INDEX_TYPE ihalf_edgeC =
        IndexOfNextHalfEdgeAroundVertex(ihalf_edgeB);

      if (ihalf_edgeC == ihalf_edgeA) {
        // Only two half edges in cap.
        return(false);
      }

      if (IndexOfNextHalfEdgeAroundVertex(ihalf_edgeC)
          == ihalf_edgeA)
        { return(true); }

      return(false);
    }


    // Return true if there are four or more half edges
    //   in the cap containing half edge ihalf_edgeA.
    // @pre Assumes half edges are linked.
    //   in the cap containing half edge ihalf_edgeA.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    AreFourOrMoreHalfEdgesAroundFromVertex
    (const ITYPE2 ihalf_edgeA) const
    {
      const HALF_EDGE_INDEX_TYPE ihalf_edgeB =
        IndexOfNextHalfEdgeAroundVertex(ihalf_edgeA);

      if (ihalf_edgeB == ihalf_edgeA) {
        // Half edge ihalf_edgeA is only half edge in cap.
        return(false);
      }

      const HALF_EDGE_INDEX_TYPE ihalf_edgeC =
        IndexOfNextHalfEdgeAroundVertex(ihalf_edgeB);

      if (ihalf_edgeC == ihalf_edgeA) {
        // Exactly two half edges in cap.
        return(false);
      }

      if (IndexOfNextHalfEdgeAroundVertex(ihalf_edgeC)
          == ihalf_edgeA) {
        // Exactly three half edges in cap.
        return(false);
      }

      return(true);
    }


    // Return true if there are exactly two half edges
    //   in the cap containing vertex iv.
    // - Note: Argument is a vertex index, not a half edge index.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    AreTwoHalfEdgesAroundVertex(const ITYPE2 iv) const
    {
      if (!IsSetIndexOfHalfEdgeIncidentOnVertex(iv))
        { return false; }

      return AreTwoHalfEdgesAroundFromVertex
        (IndexOfHalfEdgeIncidentOnVertex(iv));
    }


    // Return true if there are exactly three half edges
    //   in the cap containing vertex iv.
    // - Note: Argument is a vertex index, not a half edge index.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    AreThreeHalfEdgesAroundVertex(const ITYPE2 iv) const
    {
      if (!IsSetIndexOfHalfEdgeIncidentOnVertex(iv))
        { return false; }

      return AreThreeHalfEdgesAroundFromVertex
        (IndexOfHalfEdgeIncidentOnVertex(iv));
    }


    // Return index of previous half edge around edge.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    typename MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    HALF_EDGE_INDEX_TYPE
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    IndexOfPrevHalfEdgeAroundEdge(const ITYPE2 ihalf_edge0) const
    {
      HALF_EDGE_INDEX_TYPE ihalf_edge =
        IndexOfNextHalfEdgeAroundEdge(ihalf_edge0);

      // Use k as a guard against an infinite loop,
      //   in case data structure is corrupted.
      NTYPE k = 0;
      while (IndexOfNextHalfEdgeAroundEdge(ihalf_edge) != ihalf_edge0 &&
             k < NumHalfEdges()) {
        ihalf_edge = IndexOfNextHalfEdgeAroundEdge(ihalf_edge);
        k++;
      };

      if (IndexOfNextHalfEdgeAroundEdge(ihalf_edge) != ihalf_edge0) {
        IJK::PROCEDURE_ERROR error("MESH2D::IndexOfPrevHalfEdgeAroundEdge");

        error.AddMessage
          ("Programming error. Inconsistency in data structure MESH2D_BASE.");
        error.AddMessage
          ("  Half edges around edge (", FromVertexIndex(ihalf_edge0), ",",
           ToVertexIndex(ihalf_edge0), ") do not form a cycle.");
        error.AddMessage("  Starting half edge ", ihalf_edge0,
                         " contained in polygon ",
                         IndexOfPolygonContainingHalfEdge(ihalf_edge0), ".");
      }

      return(ihalf_edge);
    }


    // Return true if polygon contains a boundary edge.
    // @pre AreHalfEdgesLinked() is true.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
      bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
      IsBoundaryPolygon(const ITYPE2 ipoly) const
    {
      for (NUMBER_TYPE j = 0; j < NumPolygonEdges(ipoly); j++) {
        const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);

        if (IsBoundaryEdge(jhalf_edge)) { return(true); }
      }

      return(false);
    }


    // Append triangle vertices and triangle indices to lists.
    // Append vertices of triangles in mesh to triangle_vertex_list[].
    // - Append corresponding triangle indices to triangle_index[].
    // - Ignores/skips any polygons that are not triangles.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
               typename ITYPE, typename NTYPE>
    template <typename VTYPE2, typename PTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    AppendToListOfVerticesOfAllTriangles
    (std::vector<VTYPE2> & triangle_vertex_list,
     std::vector<PTYPE2> & triangle_index) const
    {
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

	if (IsPolygonDeleted(ipoly)) { 
	  // Skip deleted triangles.
	  continue; 
	}

        if (IsTriangle(ipoly)) {
          const VERTEX_INDEX_TYPE iv0 = IndexOfVertexInPolygon(ipoly,0);
          const VERTEX_INDEX_TYPE iv1 = IndexOfVertexInPolygon(ipoly,1);
          const VERTEX_INDEX_TYPE iv2 = IndexOfVertexInPolygon(ipoly,2);
          add_triangle_vertices(iv0,iv1,iv2,triangle_vertex_list);
          triangle_index.push_back(ipoly);
        }
      }
    }


    // Get triangle vertices and corresponding triangle indices.
    // - Similar to AppendToListOfAllTriangles() but first clears
    // - Ignores/skips any polygons that are not triangles.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
               typename ITYPE, typename NTYPE>
    template <typename VTYPE2, typename PTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    GetListOfVerticesOfAllTriangles
    (std::vector<VTYPE2> & triangle_vertex_list,
     std::vector<PTYPE2> & triangle_index) const
    {
      triangle_vertex_list.clear();
      triangle_index.clear();
      AppendToListOfVerticesOfAllTriangles
        (triangle_vertex_list, triangle_index);
    }


    // Sort triangle vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
               typename ITYPE, typename NTYPE>
    template <typename VTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    SortTriangleVertices
    (std::vector<VTYPE2> & triangle_vertex_list) const
    {
      typedef typename std::vector<VTYPE2>::size_type SIZE_TYPE;

      const int NUM_VERTICES_PER_TRIANGLE(3);

      for (SIZE_TYPE i = 0; i < triangle_vertex_list.size();
           i+=NUM_VERTICES_PER_TRIANGLE) {
        if (triangle_vertex_list[i] > triangle_vertex_list[i+1])
          { std::swap(triangle_vertex_list[i], triangle_vertex_list[i+1]); }
        if (triangle_vertex_list[i] > triangle_vertex_list[i+2])
          { std::swap(triangle_vertex_list[i], triangle_vertex_list[i+2]); }
        if (triangle_vertex_list[i+1] > triangle_vertex_list[i+2])
          { std::swap(triangle_vertex_list[i+1], triangle_vertex_list[i+2]); }
      }
    }


    namespace {

      // Compare vertices of two triangles, given that first vertex
      //   of each is the same.
      template <typename VTYPE>
      class _TRIANGLE_LESS_THAN_COMPARE2 {

      protected:
        const VTYPE * triangle_vertex_array;

      public:
        _TRIANGLE_LESS_THAN_COMPARE2(const VTYPE * vertex_array):
          triangle_vertex_array(vertex_array) {};

        bool operator ()(const int iA, const int iB) const
        {
          const int TWO(2);
          const int NUM_VERTICES_PER_TRIANGLE(3);

          // First triangle vertex always matches.
          // Ignore first triangle vertex.
          const VTYPE * vA =
            triangle_vertex_array + iA*NUM_VERTICES_PER_TRIANGLE + 1;
          const VTYPE * vB =
             triangle_vertex_array + iB*NUM_VERTICES_PER_TRIANGLE + 1;

          return(std::lexicographical_compare(vA, vA+TWO,vB,vB+TWO));
        }
      };

    };


    // Get sorted list of sorted triangle vertices and corresponding triangle indices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
               typename ITYPE, typename NTYPE>
    template <typename VTYPE2, typename PTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    GetSortedListOfVerticesOfAllTriangles
    (std::vector<VTYPE2> & sorted_triangle_vertex_list,
     std::vector<PTYPE2> & triangle_index) const
    {
      const int NUM_VERTICES_PER_TRIANGLE(3);
      std::vector<VTYPE2> triangle_vertex_list;
      std::vector<PTYPE2> triangle_indexB;
      IJK::PROCEDURE_ERROR error
        ("MESH2D_BASE::GetSortedListOfVerticesOfAllTriangles");

      // triangle_list[iv] is list of triangles incident on vertex iv.
      IJK::LIST_OF_LISTS<VERTEX_INDEX_TYPE,NUMBER_TYPE> triangle_list;

      // Initialize
      sorted_triangle_vertex_list.clear();
      triangle_index.clear();

      GetListOfVerticesOfAllTriangles
      (triangle_vertex_list, triangle_indexB);
      SortTriangleVertices(triangle_vertex_list);

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "Unsorted list of vertices of all triangles:" << endl;
      for (int i = 0; i < triangle_indexB.size(); i++) {
        cerr << "Triangle " << triangle_indexB[i] << ": ";
        cerr << triangle_vertex_list[3*i] << ","
       << triangle_vertex_list[3*i+1] << ","
       << triangle_vertex_list[3*i+2] << endl;
      }
      cerr << endl;
      */

      // triangle_list has a list for each vertex.
      triangle_list.SetNumLists(NumVertices());

      // For each vertex iv, count number of triangles with first vertex iv.
      for (NUMBER_TYPE i = 0; i < triangle_vertex_list.size();
           i+=NUM_VERTICES_PER_TRIANGLE) {
        const VERTEX_INDEX_TYPE iv = triangle_vertex_list[i];
        triangle_list.list_length[iv]++;
      }

      // *** DEBUG ***
      /*
      using namespace std;
      cerr << "triangle_list.list_length[]:" << endl;
      for (int iv = 0; iv < triangle_list.NumLists(); iv++) {
        cout << "  Vertex " << iv << ", list length: "
       << triangle_list.ListLength(iv) << endl;
      }
      */

      triangle_list.SetFirstElement();
      triangle_list.AllocArrayElement();

      // current_triangle_loc[iv] =
      //   Location of current triangle in List(iv).
      std::vector<NUMBER_TYPE>
        current_triangle_list_loc(NumVertices(), 0);

      for (NUMBER_TYPE i = 0; i < triangle_vertex_list.size();
           i+=NUM_VERTICES_PER_TRIANGLE) {
        const VERTEX_INDEX_TYPE iv = triangle_vertex_list[i];
        const NUMBER_TYPE iloc = current_triangle_list_loc[iv];
        triangle_list.ElementRef(iv, iloc) = i/NUM_VERTICES_PER_TRIANGLE;
        current_triangle_list_loc[iv]++;
      }

      // Check that store correct number of triangles for each vertex.
      for (NUMBER_TYPE iv = 0; iv < triangle_list.NumLists(); iv++) {
        if (current_triangle_list_loc[iv] != triangle_list.ListLength(iv)) {
          error.AddMessage
            ("Programming error. Incorrect storage of triangles.");
          error.AddMessage
            ("  Number of triangles with lowest vertex ", iv, " = ",
             triangle_list.ListLength(iv), ".");
          error.AddMessage
            ("  Stored ", current_triangle_list_loc[iv],
             " triangles for vertex ", iv, ".");
          throw error;
        }
      }

      _TRIANGLE_LESS_THAN_COMPARE2<VERTEX_INDEX_TYPE> triangle_lt
      (IJK::vector2pointer(triangle_vertex_list));

      // Sort each list in triangle_list based on triangle vertices.
      for (NUMBER_TYPE iv = 0; iv < triangle_list.NumLists(); iv++) {
        const NUMBER_TYPE list_length = triangle_list.ListLength(iv);

        if (list_length > 1) {
	  const NUMBER_TYPE first_element =
	    triangle_list.FirstElement(iv);
          VERTEX_INDEX_TYPE * list =
            IJK::vector2pointer(triangle_list.element) + first_element;
          std::sort(list, list+list_length, triangle_lt);
        }
      }

      sorted_triangle_vertex_list.resize(triangle_vertex_list.size());
      triangle_index.resize(triangle_indexB.size());
  
      auto sorted_triangle_vertex_list_iterator =
        sorted_triangle_vertex_list.begin();
      NUMBER_TYPE k = 0;
      for (NUMBER_TYPE iv = 0; iv < triangle_list.NumLists(); iv++) {
        const NUMBER_TYPE list_length = triangle_list.ListLength(iv);
  
        for (NUMBER_TYPE jloc = 0; jloc < list_length; jloc++) {
          const NUMBER_TYPE j = triangle_list.Element(iv, jloc);
          const auto triangle_vertex_list_iterator =
            triangle_vertex_list.begin()+j*NUM_VERTICES_PER_TRIANGLE;
          triangle_index[k] = triangle_indexB[j];
          std::copy(triangle_vertex_list_iterator,
                    triangle_vertex_list_iterator+NUM_VERTICES_PER_TRIANGLE,
                    sorted_triangle_vertex_list_iterator);
          sorted_triangle_vertex_list_iterator += 3;
          k++;
        }
      }

      if (sorted_triangle_vertex_list_iterator !=
          sorted_triangle_vertex_list.end()) {
        const NUMBER_TYPE num_added =
          (sorted_triangle_vertex_list_iterator - sorted_triangle_vertex_list.begin())/
          NUM_VERTICES_PER_TRIANGLE;
        const NUMBER_TYPE num_expected =
          (sorted_triangle_vertex_list.end() - sorted_triangle_vertex_list.begin())/
          NUM_VERTICES_PER_TRIANGLE;
        error.AddMessage("Programming error. Added incorrect number of triangles.");
        error.AddMessage("  Added ", num_added, " triangles to sorted_triangle_vertex_list.");
        error.AddMessage("  Should have added ", num_expected, " triangles.");
        throw error;
      }

    }


    namespace {

      template <typename VTYPE>
      bool _triangle_equals(const VTYPE * vA, const VTYPE * vB)
      {
        if (vA[0] != vB[0]) { return false; }
        if (vA[1] != vB[1]) { return false; }
        if (vA[2] != vB[2]) { return false; }

        // Otherwise:
        return true;
      }
    }


    // Get list of vertices of duplicate triangles.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
               typename ITYPE, typename NTYPE>
    template <typename VTYPE2, typename NTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    GetListOfVerticesOfDuplicateTriangles
    (std::vector<VTYPE2> & duplicate_triangle_vertex,
     IJK::LIST_OF_LISTS<NTYPE2,NTYPE2> &
       indices_of_duplicate_triangles) const
    {
      const int NUM_VERTICES_PER_TRIANGLE(3);
      std::vector<VERTEX_INDEX_TYPE> sorted_triangle_vertex_list;
      std::vector<VERTEX_INDEX_TYPE> triangle_index;

      duplicate_triangle_vertex.clear();
      indices_of_duplicate_triangles.Clear();

      GetSortedListOfVerticesOfAllTriangles
      (sorted_triangle_vertex_list, triangle_index);

      const NUMBER_TYPE num_triangles = triangle_index.size();

      NUMBER_TYPE i0 = 0;
      const VERTEX_INDEX_TYPE * vptr =
        IJK::vector2pointer(sorted_triangle_vertex_list);
      const NUMBER_TYPE * triangle_index_ptr =
        IJK::vector2pointer(triangle_index);
      while (i0 < num_triangles) {
        NUMBER_TYPE i1 = i0+1;
        while ((i1 < num_triangles) &&
               _triangle_equals(vptr+NUM_VERTICES_PER_TRIANGLE*i0,
                                vptr+NUM_VERTICES_PER_TRIANGLE*i1)) {
          i1++;
        }

        if (i0+1 != i1) {
          const NUMBER_TYPE j0 = NUM_VERTICES_PER_TRIANGLE*i0;
          const VERTEX_INDEX_TYPE iv0 = sorted_triangle_vertex_list[j0];
          const VERTEX_INDEX_TYPE iv1 = sorted_triangle_vertex_list[j0+1];
          const VERTEX_INDEX_TYPE iv2 = sorted_triangle_vertex_list[j0+2];
          add_triangle_vertices(iv0,iv1,iv2,duplicate_triangle_vertex);
          indices_of_duplicate_triangles.AddList(triangle_index_ptr+i0, i1-i0);
        }

	i0 = i1;
      }
    }

    // Add a single vertex.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    typename VTYPE::VERTEX_INDEX_TYPE
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::AddVertex()
    {
      const NTYPE numv = vertex_list.size();
      vertex_list.resize(numv+1);

      return(numv);
    }


    // Add k vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename NTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    AddVertices(const NTYPE2 k)
    {
      IJK::PROCEDURE_ERROR error("MESH2D_BASE::AddVertices");

      const NTYPE numv = vertex_list.size();

      if (k < 1) {
        error.AddMessage("Programming error.  Must add at least one vertex.");
        error.AddMessage("  Num vertices added = ", k, ".");
        throw error;
      }

      vertex_list.resize(numv+k);
    }


    // Initialize half edge.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::_InitHalfEdgeX
    (const ITYPE2 ihalf_edge, const ITYPE3 ipoly)
    {
      polygon_half_edge_list.element[ihalf_edge].polygon_index = ipoly;

      //  Next half edge is initially the current half edge.
      //  Represents edge on border.
      polygon_half_edge_list.element[ihalf_edge].
        index_of_next_half_edge_around_edge = ihalf_edge;
    }


    // Set polygon vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3, typename NTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::_SetPolygonVertices
    (const ITYPE2 ipoly, const ITYPE3 poly_vert_list[], const NTYPE2 list_length)
    {
      IJK::PROCEDURE_ERROR error("MESH2D_BASE::_AddPolygon");

      for (NTYPE2 j = 0; j < list_length; j++) {

        const ITYPE2 jv = poly_vert_list[j];

        if (!CheckVertex(jv, error)) { throw error; }

        const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);
        _SetHalfEdgeFromVertex(jhalf_edge, jv);
      }

    }


    // Allocate memory for k additional half edges.
    // - No initialization or setting of half edge vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename NTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _AllocateAdditionalHalfEdges(const NTYPE2 k)
    {
      polygon_vertex_list.AddList(k);
      polygon_half_edge_list.AddList(k);
    }


    // Remove half edge ihalf_edge from half_edges incident
    //   on FromVertex(ihalf_edge).
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _RemoveHalfEdgeIncidentOnFromVertex
    (const ITYPE2 ihalf_edge)
    {
      const VERTEX_INDEX_TYPE iv = FromVertexIndex(ihalf_edge);

      vertex_list[iv].ClearIncidentHalfEdgeIndexIfEquals(ihalf_edge);
    }

    // Add polygon with numv vertices.
    // - Does not set polygon vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename NTYPE2>
    typename MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::POLYGON_INDEX_TYPE
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::_AddPolygonX
    (const NTYPE2 numv)
    {
      const POLYGON_INDEX_TYPE ipoly = polygon_list.size();
      IJK::PROCEDURE_ERROR error("MESH2D_BASE::_AddPolygonX");

      if (!CheckNumPolygons(error)) { throw error; }

      // Increase size of polygon list.
      polygon_list.resize(ipoly+1);
      polygon_list[ipoly].Init(numv);

      // Add list of polygon vertices.
      _AllocateAdditionalHalfEdges(numv);

      // Store vertices.  Set half edges.
      for (NTYPE2 j = 0; j < numv; j++) {
        const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);
        _InitHalfEdgeX(jhalf_edge, ipoly);
      }

      return(ipoly);
    }


    // Add a polygon with list_length vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename NTYPE2>
    typename MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::POLYGON_INDEX_TYPE
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::_AddPolygon
    (const ITYPE2 poly_vert_list[], const NTYPE2 list_length)
    {
      const POLYGON_INDEX_TYPE ipoly = _AddPolygonX(list_length);
      _SetPolygonVertices(ipoly, poly_vert_list, list_length);
      polygon_list[ipoly].is_deleted = false;

      return(ipoly);
    }


    // Add two polygon with numv0 and numv1 vertices.
    // - Internal, protected version.
    // - Does not set polygon edges.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename NTYPE2, typename NTYPE3,
              typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _AddTwoPolygonsX
    (const NTYPE2 numv0, const NTYPE3 numv1,
     ITYPE2 & ipoly0, ITYPE3 & ipoly1)
    {
      ipoly0 = _AddPolygonX(numv0);
      ipoly1 = _AddPolygonX(numv1);
    }


    // Add a triangle.
    // - Does not modify flags are_half_edges_linked
    //     or is_oriented_manifold.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    typename MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::POLYGON_INDEX_TYPE
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::_AddTriangle
    (const ITYPE2 iv0, const ITYPE2 iv1, const ITYPE2 iv2)
    {
      const NTYPE NUM_VERTICES_PER_TRIANGLE(3);

      const VERTEX_INDEX_TYPE vertex_list[NUM_VERTICES_PER_TRIANGLE] =
        { iv0, iv1, iv2 };

      const POLYGON_INDEX_TYPE ipoly_new =
        _AddPolygon(vertex_list, NUM_VERTICES_PER_TRIANGLE);

      return(ipoly_new);
    }


    // Add a quadrilateral.
    // - Does not modify flags are_half_edges_linked
    //     or is_oriented_manifold.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    typename MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::POLYGON_INDEX_TYPE
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::_AddQuadrilateral
    (const ITYPE2 iv0, const ITYPE2 iv1, const ITYPE2 iv2, const ITYPE2 iv3)
    {
      const NTYPE NUM_VERTICES_PER_QUAD(4);

      const VERTEX_INDEX_TYPE vertex_list[NUM_VERTICES_PER_QUAD] =
        { iv0, iv1, iv2, iv3 };

      // *** SHOULD BE _AddPolygon()?
      const POLYGON_INDEX_TYPE ipoly_new =
        AddPolygon(vertex_list, NUM_VERTICES_PER_QUAD);

      return(ipoly_new);
    }


    // Add four triangles forming the triangulation of a hexagon.
    // - Triangulation has 3 ears and one center triangle.
    // - Links edges in triangulation interior.
    // @param itriangle[i] Index of i'th triangle.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE0, typename ITYPE1, typename ITYPE2,
              typename ITYPE3, typename ITYPE4, typename ITYPE5,
              typename ITYPE6>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _AddFourTriangles_T024
    (const ITYPE0 iv0, const ITYPE1 iv1, const ITYPE2 iv2,
     const ITYPE3 iv3, const ITYPE4 iv4, const ITYPE5 iv5,
     ITYPE6 itriangle_new[4])
    {
      itriangle_new[0] = _AddTriangle(iv5, iv0, iv4);
      itriangle_new[1] = _AddTriangle(iv1, iv2, iv0);
      itriangle_new[2] = _AddTriangle(iv3, iv4, iv2);
      itriangle_new[3] = _AddTriangle(iv0, iv2, iv4);

      const HALF_EDGE_INDEX_TYPE ihalf_edgeD0 =
        this->HalfEdgeIndex(itriangle_new[3], 0);
      const HALF_EDGE_INDEX_TYPE ihalf_edgeD1 =
        this->HalfEdgeIndex(itriangle_new[3], 1);
      const HALF_EDGE_INDEX_TYPE ihalf_edgeD2 =
        this->HalfEdgeIndex(itriangle_new[3], 2);

      this->_LinkTwoHalfEdgesAroundEdge
        (ihalf_edgeD0, this->HalfEdgeIndex(itriangle_new[1], 1));
      this->_LinkTwoHalfEdgesAroundEdge
        (ihalf_edgeD1, this->HalfEdgeIndex(itriangle_new[2], 1));
      this->_LinkTwoHalfEdgesAroundEdge
        (ihalf_edgeD2, this->HalfEdgeIndex(itriangle_new[0], 1));
    }


    // Add a polygon with list_length vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename NTYPE2>
    typename MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::POLYGON_INDEX_TYPE
    MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::AddPolygon
    (const ITYPE2 poly_vert_list[], const NTYPE2 list_length)
    {
      const POLYGON_INDEX_TYPE ipoly =
        _AddPolygon(poly_vert_list, list_length);

      // Half edges of polygon ipoly are not linked.
      are_half_edges_linked = false;

      // Manifold may still be oriented, but need to verify.
      is_oriented_manifold = false;

      return(ipoly);
    }


    // Add polygons.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename NTYPE2, typename NTYPE3, typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::AddPolygons
    (const NTYPE2 * num_poly_vert,
     const NTYPE3 num_poly,
     const ITYPE2 * poly_vert_list,
     const ITYPE3 * first_poly_vert)
    {
      for (NTYPE2 ipoly = 0; ipoly < num_poly; ipoly++) {
        AddPolygon(poly_vert_list+first_poly_vert[ipoly], num_poly_vert[ipoly]);
      }
    }


    // Compute number of polygons incident on each vertex.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename NTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    ComputeNumberOfPolygonsIncidentOnEachVertex
    (std::vector<NTYPE2> & num_incident_polygons) const
    {
      num_incident_polygons.clear();
      num_incident_polygons.resize(NumVertices(), 0);

      for (POLYGON_INDEX_TYPE jpoly = 0; jpoly < NumPolygons(); jpoly++) {

        if (IsPolygonDeleted(jpoly)) { continue; }

        for (ITYPE k = 0; k < NumPolygonEdges(jpoly); k++) {
          const VERTEX_INDEX_TYPE kv = IndexOfVertexInPolygon(jpoly, k);
          num_incident_polygons[kv]++;
        }
      }
    }


    // Return number of boundary edges.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    NTYPE MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CountNumberOfBoundaryEdges() const
    {
      NTYPE num_boundary_edges = 0;
      
      for (POLYGON_INDEX_TYPE jpoly = 0; jpoly < NumPolygons(); jpoly++) {

        if (IsPolygonDeleted(jpoly)) { continue; }

        for (ITYPE k = 0; k < NumPolygonEdges(jpoly); k++) {
          const HALF_EDGE_INDEX_TYPE ihalf_edge = HalfEdgeIndex(jpoly, k);
          
          if (IsBoundaryEdge(ihalf_edge))
            { num_boundary_edges++; }
        }
      }

      return num_boundary_edges;
    }
    
    
    // Update half edges incident on polygon vertices.
    // - Only set half edge incident on vertex iv, if the half edge incident
    //   on vertex iv is not set.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    UpdateHalfEdgesIncidentOnPolygonVertices
      (const ITYPE2 ipoly)
    {
      for (NUMBER_TYPE j = 0; j < this->NumPolygonVertices(ipoly); j++) {
        const VERTEX_INDEX_TYPE jv = IndexOfVertexInPolygon(ipoly, j);
        const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);

        if (IsSetIndexOfHalfEdgeIncidentOnVertex(jv)) {
          if (IsBoundaryEdge(jhalf_edge)) {
            const HALF_EDGE_INDEX_TYPE incident_half_edge =
              IndexOfHalfEdgeIncidentOnVertex(jv);

            if (!IsBoundaryEdge(incident_half_edge)) {
              _SetIndexOfHalfEdgeIncidentOnVertex(jv, jhalf_edge);
            }
          }
        }
        else {
          _SetIndexOfHalfEdgeIncidentOnVertex(jv, jhalf_edge);
        }
      }
    }


    // Link two half edges around edge.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _LinkTwoHalfEdgesAroundEdge
    (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1)
    {
      const HALF_EDGE_INDEX_TYPE ihalf_edge2 =
        this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge0);
      const HALF_EDGE_INDEX_TYPE ihalf_edge3 =
        this->IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);

      this->polygon_half_edge_list.element[ihalf_edge0].
        index_of_next_half_edge_around_edge = ihalf_edge3;
      this->polygon_half_edge_list.element[ihalf_edge1].
        index_of_next_half_edge_around_edge = ihalf_edge2;
    }


    // Link half edges.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::LinkHalfEdges()
    {
      // For each vertex, list of half edges containing the vertex.
      IJK::LIST_OF_LISTS<ITYPE,NTYPE> vertex_half_edge_incidence;
      std::vector<NTYPE> current_element(NumVertices());
      IJK::PROCEDURE_ERROR error("MESH2D_BASE::LinkHalfEdges");


      // Default. Initialize are_half_edges_linked to false.
      are_half_edges_linked = false;

      // Initialize all half edges to boundary.
      for (HALF_EDGE_INDEX_TYPE ihalf_edge = 0;
           ihalf_edge < polygon_half_edge_list.element.size(); ihalf_edge++)
        { _SetIndexOfNextHalfEdgeAroundEdge(ihalf_edge, ihalf_edge); }

      ComputeNumberOfPolygonsIncidentOnEachVertex
        (vertex_half_edge_incidence.list_length);

      vertex_half_edge_incidence.SetFirstElement();
      vertex_half_edge_incidence.AllocArrayElement();

      // Initialize current.
      for (VERTEX_INDEX_TYPE iv = 0; iv < NumVertices(); iv++)
        { current_element[iv] = vertex_half_edge_incidence.FirstElement(iv); }

      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

        if (IsPolygonDeleted(ipoly)) { continue; }

        for (NTYPE j = 0; j < NumPolygonVertices(ipoly); j++) {

          const NTYPE jv = IndexOfVertexInPolygon(ipoly, j);

          vertex_half_edge_incidence.element[current_element[jv]] =
            HalfEdgeIndex(ipoly, j);
          current_element[jv]++;
        }
      }

      // Check that all elements are set.
      for (VERTEX_INDEX_TYPE iv = 0; iv < NumVertices(); iv++) {
        if (current_element[iv] !=
            vertex_half_edge_incidence.FirstElement(iv) +
            vertex_half_edge_incidence.ListLength(iv)) {
          error.AddMessage
            ("Programming error.  Problem finding half_edges from vertex ",
             iv, ".");
          error.AddMessage
            ("  Expected ", vertex_half_edge_incidence.ListLength(iv),
             " half edges");
          error.AddMessage
            ("  but found ",
             current_element[iv]-vertex_half_edge_incidence.FirstElement(iv),
             " half edges.");
          throw error;
        }
      }

      // Link half edges.
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

        if (IsPolygonDeleted(ipoly)) { continue; }

        for (NTYPE j0 = 0; j0 < NumPolygonVertices(ipoly); j0++) {

          const ITYPE jhalf_edge0 = HalfEdgeIndex(ipoly, j0);

          if (IndexOfNextHalfEdgeAroundEdge(jhalf_edge0) != jhalf_edge0) {
            // Half edge already linked.
            continue;
          }

          const NTYPE jv0 = IndexOfVertexInPolygon(ipoly, j0);
          const ITYPE j1 = (j0+1)%NumPolygonEdges(ipoly);
          const NTYPE jv1 = IndexOfVertexInPolygon(ipoly, j1);

          // Link jhalf_edge0 with edges from jv1.
          for (NTYPE k = 0;
               k < vertex_half_edge_incidence.ListLength(jv1); k++) {

            const ITYPE jhalf_edge1 =
              vertex_half_edge_incidence.Element(jv1, k);

            if (ToVertexIndex(jhalf_edge1) == jv0) {
              _LinkTwoHalfEdgesAroundEdge(jhalf_edge0, jhalf_edge1);
            }
          }

          // Link jhalf_edge0 with edges from jv0.
          // (Links polygons whose orientations do not match.)
          for (ITYPE k = 0;
               k < vertex_half_edge_incidence.ListLength(jv0); k++) {

            const HALF_EDGE_INDEX_TYPE jhalf_edge1 =
              vertex_half_edge_incidence.Element(jv0, k);

            if (this->IndexOfPolygonContainingHalfEdge(jhalf_edge1) != ipoly &&
                ToVertexIndex(jhalf_edge1) == jv1) {
              _LinkTwoHalfEdgesAroundEdge(jhalf_edge0, jhalf_edge1);
            }
          }
        }
      }

      // Reset vertex_list[iv].incident_half_edge_index.
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

        if (IsPolygonDeleted(ipoly)) { continue; }

        UpdateHalfEdgesIncidentOnPolygonVertices(ipoly);
      }

      are_half_edges_linked = true;
    }


    // Determine if each mesh edge is in at most two polygons.
    // - Return true if each mesh edge is in at most two polygons.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    DetermineManifoldEdges(ITYPE2 & half_edge_index) const
    {
      IJK::PROCEDURE_ERROR error("MESH2D_BASE::DetermineManifoldEdges");

      // Initialize.
      half_edge_index = 0;

      if (!CheckAreHalfEdgesLinked(error)) { throw(error); }

      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

        if (IsPolygonDeleted(ipoly)) { continue; }

        for (ITYPE j = 0; j < NumPolygonEdges(ipoly); j++) {

          const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);
          const HALF_EDGE_INDEX_TYPE khalf_edge =
            IndexOfNextHalfEdgeAroundEdge(jhalf_edge);

          if (khalf_edge == jhalf_edge) {
            // jhalf_edge is boundary edge.
            continue;
          }

          if (IndexOfNextHalfEdgeAroundEdge(khalf_edge) != jhalf_edge) {

            // Three or more polygons incident on the same edge.
            half_edge_index = jhalf_edge;
            return(false);
          }
        }
      }

      return(true);
    }


    // Determine if polygons have matching orientations.
    // - Return true if polygons have matching orientations.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    DetermineMatchingPolygonOrientations
    (bool & is_non_manifold, ITYPE2 & jpoly, ITYPE3 & kpoly) const
    {
      IJK::PROCEDURE_ERROR error("MESH2D_BASE::DetermineOrientedPolygons");

      // Initialize.
      is_non_manifold = false;
      jpoly = 0;
      kpoly = 0;

      if (!CheckAreHalfEdgesLinked(error)) { throw(error); }

      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

        if (IsPolygonDeleted(ipoly)) { continue; }

        for (ITYPE j = 0; j < NumPolygonEdges(ipoly); j++) {

          const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);
          const HALF_EDGE_INDEX_TYPE khalf_edge =
            IndexOfNextHalfEdgeAroundEdge(jhalf_edge);

          if (khalf_edge == jhalf_edge) {
            // jhalf_edge is boundary edge.
            continue;
          }

          if (IndexOfNextHalfEdgeAroundEdge(khalf_edge) != jhalf_edge) {

            // Three or more polygons incident on the same edge.
            is_non_manifold = true;
            return(false);
          }

          const VERTEX_INDEX_TYPE jv0 = FromVertexIndex(jhalf_edge);
          const VERTEX_INDEX_TYPE jv1 = ToVertexIndex(jhalf_edge);
          const VERTEX_INDEX_TYPE kv0 = FromVertexIndex(khalf_edge);
          const VERTEX_INDEX_TYPE kv1 = ToVertexIndex(khalf_edge);

          if ((jv0 == kv0) && (jv1 == kv1)) {
            // Orientation mismatch.
            jpoly = IndexOfPolygonContainingHalfEdge(jhalf_edge);
            kpoly = IndexOfPolygonContainingHalfEdge(khalf_edge);
            return(false);
          }
          else if (!((jv0 == kv1) && (jv1 == kv0))) {
            // jhalf_edge and khalf_edge do not represent the same edge.
            error.AddMessage
              ("Programming error. Inconsistency in data structure MESH2D_BASE.");
            error.AddMessage
              ("  Half edges ", jhalf_edge, " and ", khalf_edge,
               " have different endpoints.");
            error.AddMessage
              ("  Half edge ", jhalf_edge, ".  Endpoints: (", jv0, ",", jv1, ".");
            error.AddMessage
              ("  Half edge ", khalf_edge, ".  Endpoints: (", kv0, ",", kv1, ".");
            throw error;
          }


        }
      }

      return(true);
    }


    // Determine if each mesh vertex is in a disk of polygons.
    // - Return true if each mesh vertex is in a disk of polygons.
    // @pre Passed DetermineOrientedManifoldEdges() test before
    //   calling this routine.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    DetermineOrientedManifoldVertices
    (bool & is_oriented, ITYPE2 & vertex_index) const
    {
      std::vector<NUMBER_TYPE> num_incident_polygons;
      IJK::PROCEDURE_ERROR error
        ("MESH2D_BASE::DetermineOrientedManifoldVertices");

      // Initialize.
      vertex_index = 0;

      // Assume polygons orientations match.
      is_oriented = true;

      if (!CheckAreHalfEdgesLinked(error)) { throw(error); }

      ComputeNumberOfPolygonsIncidentOnEachVertex
        (num_incident_polygons);

      for (VERTEX_INDEX_TYPE iv = 0; iv < NumVertices(); iv++) {

        if (!IsSetIndexOfHalfEdgeIncidentOnVertex(iv))
          { continue; }

        const HALF_EDGE_INDEX_TYPE ihalf_edge0 =
          IndexOfHalfEdgeIncidentOnVertex(iv);

        NTYPE k = 0;
        HALF_EDGE_INDEX_TYPE ihalf_edge1 = ihalf_edge0;

        do {
          ihalf_edge1 = IndexOfNextHalfEdgeAroundVertex(ihalf_edge1);

          if (FromVertexIndex(ihalf_edge1) != iv) {

            // Check that iv is a vertex of ihalf_edge1.
            if (ToVertexIndex(ihalf_edge1) != iv) {
              error.AddMessage
                ("Programming error. Inconsistency in data structure MESH2D_BASE.");
              error.AddMessage("  Incorrect next half edge around vertex ", iv, ".");
              error.AddMessage("  Start half edge: ", ihalf_edge0, ".");
              error.AddMessage("  Current half edge: ", ihalf_edge1, ".");
              error.AddMessage("  FromVertexIndex(", ihalf_edge1, ": ",
                               FromVertexIndex(ihalf_edge1), ".");
              error.AddMessage("  ToVertexIndex(", ihalf_edge1, ": ",
                               ToVertexIndex(ihalf_edge1), ".");
              throw error;
            }

            if (!IsBoundaryEdge(ihalf_edge1)) {
              is_oriented = false;
              return(false);
            }
          }

          k++;
        } while (ihalf_edge1 != ihalf_edge0 &&
                 !IsBoundaryEdge(ihalf_edge1) &&
                 k < NumHalfEdges());

        if (k != num_incident_polygons[iv]) {
          // Non-manifold condition at vertex iv.
          // Vertex iv is contained in two sets of polygons
          //   that do not share any edges.
          vertex_index = iv;
          return(false);
        }
      }

      return(true);
    }


    // Determine if mesh is an oriented manifold.
    // - Return true if mesh is an oriented manifold.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    DetermineOrientedManifold()
    {
      bool is_non_manifold;
      bool is_oriented;
      HALF_EDGE_INDEX_TYPE ihalf_edge;
      POLYGON_INDEX_TYPE ipoly0, ipoly1;
      VERTEX_INDEX_TYPE iv;
      IJK::PROCEDURE_ERROR error("MESH2D_BASE::DetermineIsOrientedManifold");

      // Default.  Set is_oriented_manifold to false.
      is_oriented_manifold = false;

      if (!CheckAreHalfEdgesLinked(error)) { throw(error); }

      if (!DetermineManifoldEdges(ihalf_edge))
        { return(false); }

      if (!DetermineMatchingPolygonOrientations
          (is_non_manifold, ipoly0, ipoly1))
        { return(false); }

      if (!DetermineOrientedManifoldVertices(is_oriented, iv))
        { return(false); }

      is_oriented_manifold = true;

      return(true);
    }


    // Copy.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename MTYPE>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    Copy(const MTYPE & meshA)
    {
      polygon_vertex_list.Copy(meshA.PolygonVertexList());
      polygon_half_edge_list.Copy(meshA.PolygonHalfEdgeList());
      vertex_list = meshA.VertexList();
      polygon_list = meshA.PolygonList();
      are_half_edges_linked = meshA.AreHalfEdgesLinked();
      is_oriented_manifold = meshA.IsOrientedManifold();
    }


    // *****************************************************************
    // Class MESH2D_BASE Triangulate member functions
    // *****************************************************************

    // Delete polygon.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _DeletePolygon(const ITYPE2 ipoly)
    {
      polygon_list[ipoly].is_deleted = true;

      for (ITYPE j = 0; j < NumPolygonVertices(ipoly); j++) {

        const VERTEX_INDEX_TYPE jv = IndexOfVertexInPolygon(ipoly, j);
        const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);

        vertex_list[jv].ClearIncidentHalfEdgeIndexIfEquals(jhalf_edge);
      }
    }


    // Replace half edge ihalf_edge0 with ihalf_edge1.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _ReplaceHalfEdge
    (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1)
    {
      _ReplaceHalfEdgeAroundEdge(ihalf_edge0, ihalf_edge1);
      _ReplaceHalfEdgeIncidentOnVertex(ihalf_edge0, ihalf_edge1);
    }


    // Replace numR half edges.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3, typename NTYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _ReplaceHalfEdges
    (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1,
     const NTYPE2 numR)
    {
      const POLYGON_INDEX_TYPE ipoly0 =
        IndexOfPolygonContainingHalfEdge(ihalf_edge0);

      HALF_EDGE_INDEX_TYPE ihalf_edge0X = ihalf_edge0;
      HALF_EDGE_INDEX_TYPE ihalf_edge1X = ihalf_edge1;

      for (ITYPE i = 0; i < numR; i++) {
        const VERTEX_INDEX_TYPE iv = FromVertexIndex(ihalf_edge0X);
        _SetHalfEdgeFromVertex(ihalf_edge1X, iv);

        _ReplaceHalfEdge(ihalf_edge0X, ihalf_edge1X);

        ihalf_edge0X = IndexOfNextHalfEdgeInPolygon(ihalf_edge0X);
        ihalf_edge1X = IndexOfNextHalfEdgeInPolygon(ihalf_edge1X);
      }
    }


    // Replace half edge with two half edges.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3, typename ITYPE4>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _ReplaceHalfEdgeWithTwoHalfEdges
    (const ITYPE2 ihalf_edge, const ITYPE3 iv_split,
     const ITYPE4 ihalf_edge_new0)
    {
      const HALF_EDGE_INDEX_TYPE ihalf_edge_new1 =
        IndexOfNextHalfEdgeInPolygon(ihalf_edge_new0);
      const VERTEX_INDEX_TYPE iv0 = FromVertexIndex(ihalf_edge);
      const VERTEX_INDEX_TYPE iv1 = FromVertexIndex(ihalf_edge);

      _UnlinkHalfEdgeAroundEdge(ihalf_edge);

      _SetHalfEdgeFromVertex(ihalf_edge_new0, iv0);
      _SetHalfEdgeFromVertex(ihalf_edge_new1, iv_split);

      _ReplaceHalfEdgeIncidentOnVertex(ihalf_edge, ihalf_edge_new0);
    }


    // Replace half edge in list of half edges around edge.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _ReplaceHalfEdgeAroundEdge
    (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1)
    {
      if (IsBoundaryEdge(ihalf_edge0)) {
        // Nothing to replace.
        return;
      }
      else {
        const HALF_EDGE_INDEX_TYPE iprev_half_edge =
          IndexOfPrevHalfEdgeAroundEdge(ihalf_edge0);

        _SetIndexOfNextHalfEdgeAroundEdge(iprev_half_edge, ihalf_edge1);
        _SetIndexOfNextHalfEdgeAroundEdge
          (ihalf_edge1, IndexOfNextHalfEdgeAroundEdge(ihalf_edge0));
        _SetIndexOfNextHalfEdgeAroundEdge(ihalf_edge0, ihalf_edge0);
      }
    }


    // Replace half edge incident on vertex.
    // - If incident_half_edge_index of FromVertex(ihalf_edge0)
    //   equals ihalf_edge0, replace by ihalf_edge1.
    // @pre FromVertex(ihalf_edge0) = FromVertex(ihalf_edge1);
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _ReplaceHalfEdgeIncidentOnVertex
    (const ITYPE2 ihalf_edge0, const ITYPE3 ihalf_edge1)
    {
      const VERTEX_INDEX_TYPE iv = FromVertexIndex(ihalf_edge0);

      if (iv != FromVertexIndex(ihalf_edge1)) {
        IJK::PROCEDURE_ERROR error
          ("MESH2D::_ReplaceHalfEdgeIncidentOnVertex");
        error.AddMessage
          ("Programming error.  Half edges have different from vertices.");
        error.AddMessage
          ("  FromVertexIndex(", ihalf_edge0, "): ", iv, ".");
        error.AddMessage
          ("  FromVertexIndex(", ihalf_edge1, "): ",
           FromVertexIndex(ihalf_edge1), ".");
        throw error;
      }

      if (IndexOfHalfEdgeIncidentOnVertex(iv) == ihalf_edge0)
        { _SetIndexOfHalfEdgeIncidentOnVertex(iv, ihalf_edge1); }

    }


    // Triangulate quadrilateral ipoly from vertex 0.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    TriangulateQuadFromVertex0(const ITYPE2 ipoly)
    {
      const NTYPE NUM_VERTICES_PER_TRIANGLE(3);
      const NTYPE NUM_VERTICES_PER_QUAD(4);
      VERTEX_INDEX_TYPE triangle_vert[NUM_VERTICES_PER_TRIANGLE];
      IJK::PROCEDURE_ERROR error("TriangulateQuadFromVertex0");

      if (NumPolygonEdges(ipoly) != NUM_VERTICES_PER_QUAD) {
        error.AddMessage("Programming error. Polygon ", ipoly,
                         " is not a quadrilateral.");
        throw error;
      }

      triangle_vert[0] = IndexOfVertexInPolygon(ipoly, 0);
      triangle_vert[1] = IndexOfVertexInPolygon(ipoly, 1);
      triangle_vert[2] = IndexOfVertexInPolygon(ipoly, 2);

      const POLYGON_INDEX_TYPE itriangleA = _AddTriangle(triangle_vert);

      triangle_vert[0] = IndexOfVertexInPolygon(ipoly, 0);
      triangle_vert[1] = IndexOfVertexInPolygon(ipoly, 2);
      triangle_vert[2] = IndexOfVertexInPolygon(ipoly, 3);

      const POLYGON_INDEX_TYPE itriangleB = _AddTriangle(triangle_vert);

      _ReplaceHalfEdge
        (HalfEdgeIndex(ipoly, 0), HalfEdgeIndex(itriangleA, 0));
      _ReplaceHalfEdge
        (HalfEdgeIndex(ipoly, 1), HalfEdgeIndex(itriangleA, 1));
      _ReplaceHalfEdge
        (HalfEdgeIndex(ipoly, 2), HalfEdgeIndex(itriangleB, 1));
      _ReplaceHalfEdge
        (HalfEdgeIndex(ipoly, 3), HalfEdgeIndex(itriangleB, 2));

      _LinkTwoHalfEdgesAroundEdge
        (HalfEdgeIndex(itriangleA, 2), HalfEdgeIndex(itriangleB, 0));

      _DeletePolygon(ipoly);
    }


    // Unlink a half edge from the list of half edges around an edge.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    _UnlinkHalfEdgeAroundEdge(const ITYPE2 ihalf_edge)
    {
      const HALF_EDGE_INDEX_TYPE iprev_half_edge =
        IndexOfPrevHalfEdgeAroundEdge(ihalf_edge);
      const HALF_EDGE_INDEX_TYPE inext_half_edge =
        IndexOfNextHalfEdgeAroundEdge(ihalf_edge);

      _SetIndexOfNextHalfEdgeAroundEdge(iprev_half_edge, inext_half_edge);
      _SetIndexOfNextHalfEdgeAroundEdge(ihalf_edge, ihalf_edge);
    }


    // Triangulate polygon ipoly from vertex 0.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    TriangulatePolygonFromVertex0(const ITYPE2 ipoly)
    {
      const NTYPE NUM_VERTICES_PER_TRIANGLE(3);
      const NTYPE numv = NumPolygonVertices(ipoly);
      VERTEX_INDEX_TYPE triangle_vert[NUM_VERTICES_PER_TRIANGLE];
      IJK::PROCEDURE_ERROR error("TriangulatePolygonFromVertex0");

      if (numv <= NUM_VERTICES_PER_TRIANGLE) {
        // Skip.  No triangulation.
        return;
      }

      for (NTYPE i = 1; i+1 < numv; i++) {
        triangle_vert[0] = this->IndexOfVertexInPolygon(ipoly, 0);
        triangle_vert[1] = this->IndexOfVertexInPolygon(ipoly, i);
        triangle_vert[2] = this->IndexOfVertexInPolygon(ipoly, i+1);

        const POLYGON_INDEX_TYPE itriangle = _AddTriangle(triangle_vert);

        if (i == 1) {
          _ReplaceHalfEdge
            (HalfEdgeIndex(ipoly, 0), HalfEdgeIndex(itriangle, 0));
        }
        else if (i+2 == numv) {
          _ReplaceHalfEdge
            (HalfEdgeIndex(ipoly, i+1), HalfEdgeIndex(itriangle, 2));
        }

        _ReplaceHalfEdge
          (HalfEdgeIndex(ipoly, i), HalfEdgeIndex(itriangle, 1));

        if (i > 1) {
          // Link previous triangle with itriangle.
          _LinkTwoHalfEdgesAroundEdge
            (HalfEdgeIndex(itriangle-1, 2), HalfEdgeIndex(itriangle, 0));
        }

      }

      _DeletePolygon(ipoly);
    }


    // Uniformly triangulate mesh polygons.
    // - First calls TriangulatePolygonsSharingTwoAdjacentEdges().
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    TriangulateUniform()
    {
      const NTYPE nump = NumPolygons();

      TriangulatePolygonsSharingTwoAdjacentEdges();

      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < nump; ipoly++) {

        if (IsPolygonDeleted(ipoly)) { continue; }

        TriangulatePolygonFromVertex0(ipoly);
      }
    }


    // Triangulate any polygons sharing two adjacent edges.
    // - Triangulates polygons from vertex between the two shared edges.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    TriangulatePolygonsSharingTwoAdjacentEdges()
    {
      const NTYPE nump = NumPolygons();

      for (POLYGON_INDEX_TYPE ipoly0 = 0; ipoly0 < nump; ipoly0++) {

        if (IsPolygonDeleted(ipoly0)) { continue; }

        INDEX_TYPE ishared_loc0;
        if (DoesPolygonShareTwoAdjacentEdges(ipoly0, ishared_loc0)) {

          const HALF_EDGE_INDEX_TYPE ihalf_edgeA0 =
            HalfEdgeIndex(ipoly0, ishared_loc0);
          const HALF_EDGE_INDEX_TYPE ihalf_edgeA1 =
            IndexOfNextHalfEdgeAroundEdge(ihalf_edgeA0);
          const HALF_EDGE_INDEX_TYPE ihalf_edgeB1 =
            IndexOfNextHalfEdgeInPolygon(ihalf_edgeA1);
          const POLYGON_INDEX_TYPE ipoly1 =
            IndexOfPolygonContainingHalfEdge(ihalf_edgeB1);
          const INDEX_TYPE ishared_loc1 =
            LocationOfHalfEdgeInPolygon(ihalf_edgeB1);

          TriangulatePolygonFromVertex(ipoly0, ishared_loc0);
          TriangulatePolygonFromVertex(ipoly1, ishared_loc1);
        }
      }
    }


    // Triangulate polygon ipoly from j'th vertex.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    TriangulatePolygonFromVertex(const ITYPE2 ipoly, ITYPE3 j)
    {
      const NTYPE NUM_VERTICES_PER_TRIANGLE(3);
      const NTYPE numv = NumPolygonVertices(ipoly);
      VERTEX_INDEX_TYPE triangle_vert[NUM_VERTICES_PER_TRIANGLE];
      IJK::PROCEDURE_ERROR error("TriangulatePolygonFromVertex");

      if (IsTriangle(ipoly)) { return; }

      if (j >= numv) {
        error.AddMessage
          ("Programming error.  Polygon ", ipoly, " has ", numv, " vertices.");
        error.AddMessage
          ("  Cannot triangulate from vertex at location ", j, ".");
        throw error;
      }

      const NTYPE j0 = j;

      for (NTYPE i = 1; i+1 < numv; i++) {

        const NTYPE j1 = (j+i)%NumPolygonEdges(ipoly);
        const NTYPE j2 = (j+i+1)%NumPolygonEdges(ipoly);

        triangle_vert[0] = IndexOfVertexInPolygon(ipoly, j0);
        triangle_vert[1] = IndexOfVertexInPolygon(ipoly, j1);
        triangle_vert[2] = IndexOfVertexInPolygon(ipoly, j2);

        const POLYGON_INDEX_TYPE itriangle = _AddTriangle(triangle_vert);

        if (i == 1) {
          _ReplaceHalfEdge
            (HalfEdgeIndex(ipoly, j0), HalfEdgeIndex(itriangle, 0));
        }
        else if (i+2 == numv) {
          _ReplaceHalfEdge
            (HalfEdgeIndex(ipoly, j2), HalfEdgeIndex(itriangle, 2));
        }

        _ReplaceHalfEdge
          (HalfEdgeIndex(ipoly, j1), HalfEdgeIndex(itriangle, 1));

        if (i > 1) {
          // Link previous triangle with itriangle.
          _LinkTwoHalfEdgesAroundEdge
            (HalfEdgeIndex(itriangle-1, 2), HalfEdgeIndex(itriangle, 0));
        }

      }

      _DeletePolygon(ipoly);
    }

    // Triangulate polygon ipoly from interior vertex iv.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2, typename ITYPE3>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    TriangulatePolygonFromInteriorVertex(const ITYPE2 ipoly, ITYPE3 ivX)
    {
      const NTYPE NUM_VERTICES_PER_TRIANGLE(3);
      const NTYPE numv = NumPolygonVertices(ipoly);
      VERTEX_INDEX_TYPE triangle_vert[NUM_VERTICES_PER_TRIANGLE];
      POLYGON_INDEX_TYPE itriangle0;
      IJK::PROCEDURE_ERROR error("TriangulatePolygonFromInteriorVertex");

      for (NTYPE i0 = 0; i0 < numv; i0++) {

        const NTYPE i1 = (i0+1)%numv;

        triangle_vert[0] = ivX;
        triangle_vert[1] = IndexOfVertexInPolygon(ipoly, i0);
        triangle_vert[2] = IndexOfVertexInPolygon(ipoly, i1);

        const POLYGON_INDEX_TYPE itriangle = _AddTriangle(triangle_vert);

        _ReplaceHalfEdge
          (HalfEdgeIndex(ipoly, i0), HalfEdgeIndex(itriangle, 1));

        if (i0 == 0) {
          itriangle0 = itriangle;
        }
        else {
          // Link previous triangle with itriangle.
          _LinkTwoHalfEdgesAroundEdge
            (HalfEdgeIndex(itriangle-1, 2), HalfEdgeIndex(itriangle, 0));
        }

        if (i0+1 == numv) {
          // Link last triangle with first triangle.
          _LinkTwoHalfEdgesAroundEdge
            (HalfEdgeIndex(itriangle, 2), HalfEdgeIndex(itriangle0, 0));
        }
      }

      _DeletePolygon(ipoly);
    }


    // *****************************************************************
    // Class MESH2D_BASE Print member functions
    // *****************************************************************

    // Print half edge.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::PrintHalfEdgeX
      (OSTREAM_TYPE & out, const ITYPE2 ihalf_edge,
       const char c0, const char c1, const char c2) const
    {
      out << c0 << FromVertexIndex(ihalf_edge)
          << c1 << ToVertexIndex(ihalf_edge) << c2;
    }


    // Print half edge.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::PrintHalfEdge
      (OSTREAM_TYPE & out, const STYPE0 & s0,
       const ITYPE2 ihalf_edge, const STYPE1 & s1) const
    {
      out << s0;
      PrintHalfEdge(out, ihalf_edge);
      out << s1;
    }


    // Print half edge index and endpoints.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintHalfEdgeIndexAndEndpointsX
      (OSTREAM_TYPE & out, const ITYPE2 ihalf_edge,
       const char c0, const char c1, const char c2) const
    {
      out << ihalf_edge << " ";
      PrintHalfEdgeX(out, ihalf_edge, c0, c1, c2);
    }


    // Print half edge index and endpoints.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintHalfEdgeIndexAndEndpoints
      (OSTREAM_TYPE & out, const STYPE0 & s0,
       const ITYPE2 ihalf_edge, const STYPE1 & s1) const
    {
      out << s0;
      PrintHalfEdgeIndexAndEndpoints(out, ihalf_edge);
      out << s1;
    }


    // @brief Print half edge information.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename STYPE0, typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintHalfEdgeInfo
    (OSTREAM_TYPE & out, const STYPE0 & prefix,
     const ITYPE2 ihalf_edge) const
    {
      const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge);
      const HALF_EDGE_INDEX_TYPE iprev_in_poly =
          IndexOfPrevHalfEdgeInPolygon(ihalf_edge);
      const HALF_EDGE_INDEX_TYPE inext_in_poly =
          IndexOfNextHalfEdgeInPolygon(ihalf_edge);
      const HALF_EDGE_INDEX_TYPE inext_half_edge_around_edge =
          IndexOfNextHalfEdgeAroundEdge(ihalf_edge);

      out << prefix;
      PrintHalfEdgeIndexAndEndpoints(out, "Half edge ", ihalf_edge, "\n");
      out << prefix;
      PrintPolygonIndexAndVertices
        (out, "  Contained in polygon ", ipoly, "\n");
      out << prefix;
      PrintHalfEdgeIndexAndEndpoints
         (out, "  Previous half edge in polygon: ", iprev_in_poly, "\n");
      out << prefix;
      PrintHalfEdgeIndexAndEndpoints
        (out, "  Next half edge in polygon: ", inext_in_poly, "\n");
      out << prefix;
      PrintHalfEdgeIndexAndEndpoints
        (out, "  Next half edge around edge: ", inext_half_edge_around_edge, "\n");
    }



    // Print information of all half edges.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
               typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintInfoOfAllHalfEdges(OSTREAM_TYPE & out) const
    {
      for (HALF_EDGE_INDEX_TYPE ihalf_edge = 0;
           ihalf_edge < NumHalfEdges(); ihalf_edge++) {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge);

        if (IsPolygonDeleted(ipoly)) {
          out << "  Half edge " << ihalf_edge
              << " is in a deleted polygon.\n";
        }
        else {
          PrintHalfEdgeInfo(out, "  ", ihalf_edge);
        }
      }
    }


    // Print polygon vertices.
    // - Extended version.  Define left/right delimiters and separator.
    // @param c0 Left delimiter.
    // @param c1 Separator.
    // @param c2 Right delimiter.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintPolygonVerticesX
    (OSTREAM_TYPE & out, const ITYPE2 ipoly,
     const char c0, const char c1, const char c2) const
    {
      out << c0;

      if (NumPolygonVertices(ipoly) > 0) {
        out << IndexOfVertexInPolygon(ipoly, 0);
        for (NTYPE i = 1; i < NumPolygonVertices(ipoly); i++) {
          out << c1 << IndexOfVertexInPolygon(ipoly, i);
        }
      }

      out << c2;
    }


    // Print polygon vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintPolygonVertices
      (OSTREAM_TYPE & out, const STYPE0 & s0,
       const ITYPE2 ipoly, const STYPE1 & s1) const
    {
      out << s0;
      PrintPolygonVertices(out, ipoly);
      out << s1;
    }


    // Print polygon index and vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintPolygonIndexAndVerticesX
      (OSTREAM_TYPE & out, const ITYPE2 ipoly,
       const char c0, const char c1, const char c2) const
    {
      out << ipoly << " ";
      PrintPolygonVerticesX(out, ipoly, c0, c1, c2);
    }



    // Print polygon index and vertices.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintPolygonIndexAndVertices
      (OSTREAM_TYPE & out, const STYPE0 & s0,
       const ITYPE2 ipoly, const STYPE1 & s1) const
    {
      out << s0;
      PrintPolygonIndexAndVertices(out, ipoly);
      out << s1;
    }


    // Print polygon vertices and half edges.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE, typename ITYPE2,
	      typename STYPE0>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintPolygonInfo
    (OSTREAM_TYPE & out, const STYPE0 & prefix,
     const ITYPE2 ipoly) const
    {
      out << prefix;
      PrintPolygonIndexAndVertices(out, "Polygon ", ipoly, "\n");
      out << prefix;
      out << "  Half edges: ";
      for (NUMBER_TYPE j = 0; j < NumPolygonVertices(ipoly); j++) {
	const HALF_EDGE_INDEX_TYPE jhalf_edge =
	  HalfEdgeIndex(ipoly, j);
	PrintHalfEdgeIndexAndEndpoints(out, "  ", jhalf_edge, "");
      }
      out << "\n";
    }


    // Print information of all polygons.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename OSTREAM_TYPE>
    void MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    PrintInfoOfAllPolygons(OSTREAM_TYPE & out) const
    {
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {
	if (IsPolygonDeleted(ipoly)) {
	  out << "  Polygon " << ipoly << " is deleted." << "\n";
	}
	else {
	  PrintPolygonInfo(out, "  ", ipoly);
	}
      }
    }


    // *****************************************************************
    // Class MESH2D_BASE Check member functions
    // *****************************************************************

    // Check data structure.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    Check(IJK::ERROR & error) const
    {
      if (!CheckNumPolygons(error)) { return(false); }
      if (!CheckHalfEdgesIncidentOnVertices(error)) { return(false); }
      if (!CheckPolygonsContainingHalfEdges(error)) { return(false); }

      if (AreHalfEdgesLinked()) {
        if (!CheckHalfEdgeLinks(error)) { return(false); }
      }

      if (IsOrientedManifold()) {
        if (!CheckOrientedManifold(error)) { return(false); }
      }

      return(true);
    }


    // Check number of polygons.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckNumPolygons(IJK::ERROR & error) const
    {
      if (polygon_half_edge_list.NumLists() != polygon_list.size()) {
        error.AddMessage
          ("Programming error. Inconsistency in data structure MESH2D_BASE.");
        error.AddMessage
          ("  Number of polygon half edge lists does not match number of polygons.");
        error.AddMessage("  Number of polygons: ", polygon_list.size(), "");
        error.AddMessage("  Number of polygon half edge lists: ",
                         polygon_half_edge_list.NumLists(), "");
        return(false);
      }

      return(true);
    }

    // Check that iv is a vertex index.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckVertex(const ITYPE2 iv, IJK::ERROR & error) const
    {
      if (iv < 0) {
        error.AddMessage
          ("Error.  Illegal vertex index value ", iv, ".");
        error.AddMessage("  Vertex indices must be non-negative.");
        return(false);
      }

      if (iv >= NumVertices()) {
        if (NumVertices() <= 0) {
          error.AddMessage
            ("Programming error.  Add vertices before adding any polygons.");
          return(false);
        }
        else {
          error.AddMessage
            ("Error.  Illegal vertex index value ", iv, ".");
          error.AddMessage
            ("  Vertex index must be in bounds [",
             0, ",", NumVertices()-1, "].");
          error.AddMessage("  Vertex must be added to data structure");
          error.AddMessage("    before it is included in any polygon.");
          return(false);
        }
      }

      return(true);
    }


    // Check that some half edge is incident on vertex iv.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckHalfEdgeIncidentOnVertex
    (const ITYPE2 ihalf_edge, IJK::ERROR & error) const
    {
      const VERTEX_INDEX_TYPE iv = FromVertexIndex(ihalf_edge);

      if (!IsSetIndexOfHalfEdgeIncidentOnVertex(iv)) {
        const POLYGON_INDEX_TYPE ipoly =
          IndexOfPolygonContainingHalfEdge(ihalf_edge);

        error.AddMessage
          ("Programming error. Inconsistency in data structure MESH2D_BASE.");
        error.AddMessage
          ("  Half edge ", ihalf_edge, " in polygon ", ipoly,
           " is incident on vertex ", iv, ".");
        error.AddMessage
          ("  Flag indicates that half edge incident on vertex ", iv,
           " is not set.");

        return(false);
      }

      return(true);
    }


    // Check polygon_half_edge_list and polygon_vertex_list match.
    // - Check same number of lists and corresponding lists have
    //   the same length.
    // - Return true if data structure passes check.

    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckHalfEdges(IJK::ERROR & error) const
    {
      if (polygon_vertex_list.NumLists() != NumPolygons()) {
        error.AddMessage
          ("Programming error. Inconsistency in data structure MESH2D_BASE.");
        error.AddMessage
          ("  Number of polygon vertex lists does not equal number of polygons.");
        error.AddMessage("  Number of polygons: ", NumPolygons(), ".");
        error.AddMessage("  Number of polygon vertex lists: ",
                         polygon_vertex_list.NumLists(), ".");

        return(false);
      }

      if (polygon_half_edge_list.NumLists() != NumPolygons()) {
        error.AddMessage
          ("Programming error. Inconsistency in data structure MESH2D_BASE.");
        error.AddMessage
          ("  Number of polygon half edge lists does not equal number of polygons.");
        error.AddMessage("  Number of polygons: ", NumPolygons(), ".");
        error.AddMessage("  Number of polygon half edge lists: ",
                         polygon_half_edge_list.NumLists(), ".");

        return(false);
      }

      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {
        if (polygon_vertex_list.ListLength(ipoly) !=
            polygon_half_edge_list.ListLength(ipoly)) {

          error.AddMessage
            ("Programming error. Inconsistency in data structure MESH2D_BASE.");
          error.AddMessage
            ("  Polygon ", ipoly, " has ",
             polygon_vertex_list.ListLength(ipoly), " vertices");
          error.AddMessage
            ("    but ", polygon_half_edge_list.ListLength(ipoly),
             " half edges.");
          return(false);
        }
      }

      return(true);
    }


    // Check vertices.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckHalfEdgesIncidentOnVertices(IJK::ERROR & error) const
    {
      for (VERTEX_INDEX_TYPE iv = 0; iv < NumVertices(); iv++) {

        if (IsSetIndexOfHalfEdgeIncidentOnVertex(iv)) {

          const HALF_EDGE_INDEX_TYPE ihalf_edge =
            IndexOfHalfEdgeIncidentOnVertex(iv);

          if (FromVertexIndex(ihalf_edge) != iv) {
            error.AddMessage
              ("Programming error. Inconsistency in data structure MESH2D_BASE.");
            error.AddMessage
              ("  vertex[", iv, "].incident_half_edge_index: ",
               ihalf_edge, ".");
            error.AddMessage
              ("  Half edge ", ihalf_edge, " endpoints: (",
               FromVertexIndex(ihalf_edge), ",",
               ToVertexIndex(ihalf_edge), ").");
            error.AddMessage
              ("  Vertex ", iv, " should be first vertex of half edge ",
               ihalf_edge, ".");

            return(false);
          }

          const POLYGON_INDEX_TYPE ipoly =
            IndexOfPolygonContainingHalfEdge(ihalf_edge);
          if (IsPolygonDeleted(ipoly)) {
            error.AddMessage
              ("Programming error. Inconsistency in data structure MESH2D_BASE.");
            error.AddMessage
              ("  vertex[", iv, "].incident_half_edge_index: ",
               ihalf_edge, ".");
            error.AddMessage
              ("  Half edge ", ihalf_edge, " (", FromVertexIndex(ihalf_edge),
               ",", ToVertexIndex(ihalf_edge), ") is contained in deleted polygon ",
               ipoly, ".");
            error.AddMessage
              ("  vertex[", iv,
               "].incident_half_edge_index should be contained");
            error.AddMessage
              ("    in an undeleted polygon.");
            return(false);
          }
        }
      }


      if (AreHalfEdgesLinked()) {

        // Check that vertex_list[jv].incident_half_edge is a boundary edge
        //   if jv is on the boundary.

        for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

          if (IsPolygonDeleted(ipoly)) { continue; }

          for (ITYPE j = 0; j < NumPolygonEdges(ipoly); j++) {

            const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);

            if (IsBoundaryEdge(jhalf_edge)) {
              const VERTEX_INDEX_TYPE jv = FromVertexIndex(jhalf_edge);

              if (!CheckHalfEdgeIncidentOnVertex(jhalf_edge, error))
                { return(false); }

              const HALF_EDGE_INDEX_TYPE jhalf_edge2 =
                IndexOfHalfEdgeIncidentOnVertex(jv);

              if (!IsBoundaryEdge(jhalf_edge2)) {
                error.AddMessage
                  ("Programming error. Inconsistency in data structure MESH2D_BASE.");
                error.AddMessage
                  ("  Vertex ", jv, " is on boundary half edge ", jhalf_edge, "");
                error.AddMessage
                  ("  but vertex_list[", jv,
                   "].incident_half_edge_index is non-boundary half edge ", jhalf_edge2, ".");
                error.AddMessage
                  ("  jhalf_edge ", jhalf_edge,
                   ".  Endpoints: (", FromVertexIndex(jhalf_edge), ",",
                   ToVertexIndex(jhalf_edge), ").  In polygon ",
                   IndexOfPolygonContainingHalfEdge(jhalf_edge), ".");
                error.AddMessage
                  ("  jhalf_edge2 ", jhalf_edge2,
                   ".  Endpoints: (", FromVertexIndex(jhalf_edge2), ",",
                   ToVertexIndex(jhalf_edge2), ").  In polygon ",
                   IndexOfPolygonContainingHalfEdge(jhalf_edge2), ".");

                return(false);
              }
            }
          }
        }
      }

      return(true);
    }


    // Check polygons containing half edges.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckPolygonsContainingHalfEdges(IJK::ERROR & error) const
    {
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

        if (IsPolygonDeleted(ipoly)) { continue; }

        for (ITYPE j = 0; j < NumPolygonEdges(ipoly); j++) {

          const HALF_EDGE_INDEX_TYPE khalf_edge = HalfEdgeIndex(ipoly, j);

          if (IndexOfPolygonContainingHalfEdge(khalf_edge) != ipoly) {

            error.AddMessage
              ("Programming error. Inconsistency in data structure MESH2D_BASE.");
            error.AddMessage
              ("  Polygon ", ipoly, " contains half edge ", khalf_edge, ".");
            error.AddMessage
              ("  half_edge[", khalf_edge, "].polygon_index: ",
               IndexOfPolygonContainingHalfEdge(khalf_edge), ".");
            return(false);
          }
        }
      }

      return(true);
    }


    // Check that iloc is a valid location in polygon ipoly.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename PTYPE2, typename ITYPE2>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckPolygonLocation
    (const PTYPE2 ipoly, const ITYPE2 iloc, IJK::ERROR & error) const
    {
      if (iloc < 0) {
        error.AddMessage
          ("Error.  Illegal location ", iloc, " in polygon ",
           ipoly, ".");
        error.AddMessage("  Polygon location must be non-negative.");
        return(false);
      }

      if (iloc >= NumPolygonVertices(ipoly)) {
        error.AddMessage
          ("Error.  Illegal polygon location ", iloc, " in polygon ",
           ipoly, ".");
        error.AddMessage
          ("  Polygon only has ", NumPolygonVertices(ipoly), " vertices.");
        error.AddMessage("  Polygon location must be less than ",
                         NumPolygonVertices(ipoly), ".");
        return(false);
      }

      return(true);
    }


    // Check half edge links.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckHalfEdgeLinks(IJK::ERROR & error) const
    {
      if (!AreHalfEdgesLinked()) {
        error.AddMessage("Programming error.  Half edges are not linked.");
        error.AddMessage
          ("  Call LinkHalfEdges() after adding polygons.");
        return(false);
      }

      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < NumPolygons(); ipoly++) {

        if (IsPolygonDeleted(ipoly)) { continue; }

        for (ITYPE j = 0; j < NumPolygonEdges(ipoly); j++) {

          const HALF_EDGE_INDEX_TYPE jhalf_edge = HalfEdgeIndex(ipoly, j);
          const VERTEX_INDEX_TYPE jv0 = FromVertexIndex(jhalf_edge);
          const VERTEX_INDEX_TYPE jv1 = ToVertexIndex(jhalf_edge);
          HALF_EDGE_INDEX_TYPE khalf_edge =
            IndexOfNextHalfEdgeAroundEdge(jhalf_edge);

          NTYPE k = 0;
          do {
            const VERTEX_INDEX_TYPE kv0 = FromVertexIndex(khalf_edge);
            const VERTEX_INDEX_TYPE kv1 = ToVertexIndex(khalf_edge);
            const POLYGON_INDEX_TYPE kpoly =
              IndexOfPolygonContainingHalfEdge(khalf_edge);

            if (!((jv0 == kv0) && (jv1 == kv1)) &&
                !((jv0 == kv1) && (jv1 == kv0))) {
              // jhalf_edge and khalf_edge do not represent the same edge.
              error.AddMessage
                ("Programming error. Inconsistency in data structure MESH2D_BASE.");
              error.AddMessage
                ("  Half edges ", jhalf_edge, " and ", khalf_edge,
                 " represent the same edge but have different endpoints.");
              error.AddMessage
                ("  Half edge ", jhalf_edge,
                 ".  Endpoints: (", jv0, ",", jv1, ").");
              error.AddMessage
                ("  Half edge ", khalf_edge,
                 ".  Endpoints: (", kv0, ",", kv1, ").");
              return(false);
            }

            if (IsPolygonDeleted(kpoly)) {
              error.AddMessage
                ("Programming error. Inconsistency in data structure MESH2D_BASE.");
              error.AddMessage
                ("  Half edges ", khalf_edge, " and ", jhalf_edge,
                 " are linked but polygon ", kpoly, " is deleted.");
              error.AddMessage
                ("  Half edge ", khalf_edge,
                 ".  Endpoints: (", kv0, ",", kv1, ").  In deleted polygon ", kpoly, ".");
              error.AddMessage
                ("  Half edge ", jhalf_edge,
                 ".  Endpoints: (", jv0, ",", jv1, ").  In active polygon ", ipoly, ".");
              return(false);
            }

            khalf_edge = IndexOfNextHalfEdgeAroundEdge(khalf_edge);
            k++;
          } while (khalf_edge != jhalf_edge && k < NumHalfEdges());

          if (jhalf_edge != khalf_edge) {
            error.AddMessage
              ("Programming error. Inconsistency in data structure MESH2D_BASE.");
            error.AddMessage
              ("  Half edges around edge (", jv0, ",", jv1,
               ") do not form a cycle.");
            return(false);
          }
        }
      }

      return(true);
    }


    // Check that AreHalfEdgesLinked() is true.
    // - Return true if AreHalfEdgesLinked() is true.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckAreHalfEdgesLinked(IJK::ERROR & error) const
    {
      if (AreHalfEdgesLinked()) {
        return(true);
      }
      else {
        error.AddMessage("Programming error.  Half edges are not linked.");
        error.AddMessage
          ("  Call LinkHalfEdges() after adding polygons.");
        return(false);
      }
    }

    // Check that data structure is an oriented manifold.
    // - Return true if data structure passes check.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckOrientedManifold(IJK::ERROR & error) const
    {
      bool is_non_manifold;
      bool is_oriented;
      HALF_EDGE_INDEX_TYPE ihalf_edge;
      POLYGON_INDEX_TYPE ipoly0, ipoly1;
      VERTEX_INDEX_TYPE iv;

      if (!CheckAreHalfEdgesLinked(error)) { return(false); }

      if (!DetermineManifoldEdges(ihalf_edge)) {
        error.AddMessage
          ("Programming error.  Mesh is not a manifold.");
        error.AddMessage
          ("  Edge (", FromVertexIndex(ihalf_edge), ",",
           ToVertexIndex(ihalf_edge),
           ") is in three or more polygons.");
        return(false);
      }

      if (!DetermineMatchingPolygonOrientations
          (is_non_manifold, ipoly0, ipoly1)) {
        error.AddMessage
          ("Programming error.  Mesh is not oriented.");
        error.AddMessage
          ("  Mismatched orientations of polygons ",
           IndexOfPolygonContainingHalfEdge(ipoly0), " and ",
           IndexOfPolygonContainingHalfEdge(ipoly1), ".");

        return(false);
      }


      if (!DetermineOrientedManifoldVertices(is_oriented, iv)) {

        if (is_oriented) {
          error.AddMessage("  Non-manifold vertex ", iv, ".");
        }
        else {
          error.AddMessage
            ("  Non-manifold condition detected in DetermineOrientedManifoldVertices.");
        }

        return(false);
      }

      return(true);
    }


    // Check that two half edges are in the same polygon.
    // - Return true if two half edges are in the same polygon.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPEH>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckAreHalfEdgesInSamePolygon
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     IJK::ERROR & error) const
    {
      if (IndexOfPolygonContainingHalfEdge(ihalf_edgeA) !=
          IndexOfPolygonContainingHalfEdge(ihalf_edgeB)) {

        error.AddMessage
          ("Programming error.  Half edges ", ihalf_edgeA, " and ",
           ihalf_edgeB, " are not in the same polygon.");

        return(false);
      }

      return(true);
    }


    // Check that three half edges are in the same polygon.
    // - Return true if two half edges are in the same polygon.
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPEH>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckAreHalfEdgesInSamePolygon
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     const ITYPEH ihalf_edgeC, IJK::ERROR & error) const
    {
      POLYGON_INDEX_TYPE ipoly =
        IndexOfPolygonContainingHalfEdge(ihalf_edgeA);

      if (IndexOfPolygonContainingHalfEdge(ihalf_edgeB) != ipoly) {
        error.AddMessage
          ("Programming error.  Half edges ", ihalf_edgeA, " and ",
           ihalf_edgeB, " are not in the same polygon.");
        return(false);
      }
      else if (IndexOfPolygonContainingHalfEdge(ihalf_edgeC) != ipoly) {
        error.AddMessage
          ("Programming error.  Half edges ", ihalf_edgeA, " and ",
           ihalf_edgeC, " are not in the same polygon.");
        return(false);
      }

      return(true);
    }


    // Check that two half edges are different.
    // - Return true if two half edges are different
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPEH>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckAreHalfEdgesDifferent
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     IJK::ERROR & error) const
    {
      if (ihalf_edgeA == ihalf_edgeB) {
        error.AddMessage("Programming error.  Half edges must be different.");
        error.AddMessage("  ihalf_edgeA: ", ihalf_edgeA, "");
        error.AddMessage("  ihalf_edgeB: ", ihalf_edgeB, "");

        return(false);
      }

      return(true);
    }


    // Check that three half edges are all different.
    // - Return true if three half edges are all different
    template <typename VTYPE, typename HALFE_TYPE, typename PTYPE,
              typename ITYPE, typename NTYPE>
    template <typename ITYPEH>
    bool MESH2D_BASE<VTYPE,HALFE_TYPE,PTYPE,ITYPE,NTYPE>::
    CheckAreHalfEdgesDifferent
    (const ITYPEH ihalf_edgeA, const ITYPEH ihalf_edgeB,
     const ITYPEH ihalf_edgeC, IJK::ERROR & error) const
    {
      if (!CheckAreHalfEdgesDifferent(ihalf_edgeA, ihalf_edgeB, error))
        { return(false); }

      if (!CheckAreHalfEdgesDifferent(ihalf_edgeA, ihalf_edgeC, error))
        { return(false); }

      if (!CheckAreHalfEdgesDifferent(ihalf_edgeB, ihalf_edgeC, error))
        { return(false); }

      return(true);
    }


    // *****************************************************************
    // Functions on MESH2D
    // *****************************************************************

    /*!
     *  @brief Return true if meshA and meshB have same number of polygons.
     */
    template <typename MESHA_TYPE, typename MESHB_TYPE>
    bool check_are_num_mesh2D_polygons_equal
    (const MESHA_TYPE & meshA, const MESHB_TYPE & meshB,
     IJK::ERROR & error)
    {
      if (meshA.NumPolygons() != meshB.NumPolygons()) {
        error.AddMessage("Meshes have different number of polygons.");
        error.AddMessage
          ("  Num triangles in meshA: ", meshA.NumPolygons(), ".");
        error.AddMessage
          ("  Num triangles in meshB: ", meshB.NumPolygons(), ".");
        return false;
      }

      return true;
    }


    /*!
     *  @brief Return true if meshA and meshB have same number 
     *    of polygons, (including deleted polygons,) 
     *    and corresponding polygons have same vertices.
     *  - Check that same polygons are deleted.
     *  - Don't check vertices of deleted polygons
     */
    template <typename MESHA_TYPE, typename MESHB_TYPE>
    bool check_are_mesh2D_polygons_equal
    (const MESHA_TYPE & meshA, const MESHB_TYPE & meshB,
     IJK::ERROR & error)
    {
      typedef typename MESHA_TYPE::VERTEX_INDEX_TYPE
        VERTEX_INDEX_TYPEA;
      typedef typename MESHB_TYPE::VERTEX_INDEX_TYPE
        VERTEX_INDEX_TYPEB;
      typedef typename MESHA_TYPE::POLYGON_INDEX_TYPE
	POLYGON_INDEX_TYPE;
      typedef typename MESHA_TYPE::NUMBER_TYPE
	NUMBER_TYPE;

      // Check meshes have same number of polygons.
      if (!check_are_num_mesh2D_polygons_equal(meshA, meshB, error)) 
	{ return false; }

      // Check that same polygons are deleted in each mesh.
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < meshA.NumPolygons(); 
	   ipoly++) {

        if (meshA.IsPolygonDeleted(ipoly) != meshB.IsPolygonDeleted(ipoly)) {
          error.AddMessage("Meshes have different deleted polygons.");
          if (meshA.IsPolygonDeleted(ipoly)) {
            error.AddMessage
              ("  Polygon ", ipoly, " is deleted in meshA but not in meshB.");
          }
          else {
            error.AddMessage
              ("  Polygon ", ipoly, " is deleted in meshB but not in meshA.");
          }
          return false;
        }
      }


      // Check that same polygons have same number of vertices.
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < meshA.NumPolygons(); 
	   ipoly++) {

	if (meshA.IsPolygonDeleted(ipoly) &&
	    meshB.IsPolygonDeleted(ipoly)) {
	  // Skip deleted polygons.
	  continue;
	}

        if (meshA.NumPolygonVertices(ipoly) != meshB.NumPolygonVertices(ipoly)) {
          error.AddMessage("Mesh polygons have different number of vertices.");
          error.AddMessage
            ("  In meshA, polygon ", ipoly, " has ", 
	     meshA.NumPolygonVertices(ipoly), " vertices.");
          error.AddMessage
            ("  In meshB, polygon ", ipoly, " has ", 
	     meshB.NumPolygonVertices(ipoly), " vertices.");
          return false;
        }
      }

      // Check that same polygons have same vertices.
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < meshA.NumPolygons(); 
	   ipoly++) {

	if (meshA.IsPolygonDeleted(ipoly) &&
	    meshB.IsPolygonDeleted(ipoly)) {
	  // Skip deleted polygons.
	  continue;
	}

	for (NUMBER_TYPE jloc = 0; jloc < meshA.NumPolygonVertices(ipoly); 
	     jloc++) {

	  const VERTEX_INDEX_TYPEA jvA = meshA.VertexIndex(ipoly, jloc);
	  const VERTEX_INDEX_TYPEB jvB = meshB.VertexIndex(ipoly, jloc);

	  if (jvA != jvB) {
	    error.AddMessage("Mesh polygons have different vertices.");
	    error.AddMessage
	      ("  Mesh A, polygon ", ipoly, ", location ", jloc, 
	       ", has vertex ", jvA, ".");
	    error.AddMessage
	      ("  Mesh B, polygon ", ipoly, ", location ", jloc, 
	       ", has vertex ", jvB, ".");

	    return false;
	  }
	}
      }

      return true;
    }


    /*!
     *  @brief Return true if from and to vertices of half edges are equal.
     *  - @param ihalf_edgeA Index of half edge in mesh A.
     *  - @param ihalf_edgeB Index of half edge in mesh B.
     */
    template <typename MESHA_TYPE, typename MESHB_TYPE,
	      typename ITYPEA, typename ITYPEB>
    bool are_mesh2D_half_edges_equal
    (const MESHA_TYPE & meshA, const ITYPEA ihalf_edgeA,
     const MESHB_TYPE & meshB, const ITYPEB ihalf_edgeB)
    {
      if (meshA.FromVertexIndex(ihalf_edgeA) !=
	  meshB.FromVertexIndex(ihalf_edgeB)) 
	{ return false; }

      if (meshA.ToVertexIndex(ihalf_edgeA) !=
	  meshB.ToVertexIndex(ihalf_edgeB)) 
	{ return false; }

      return true;
    }

    
    /*!
     *  @brief Return true if all polygons have same number of
     *    half edges, and corresponding polygons have same half_edges,
     *    and corresponding half edges have same previous, next
     *    half edges in polygon and same next half edge around polygon.
     *  - Don't check half edges of deleted polygons
     */
    template <typename MESHA_TYPE, typename MESHB_TYPE>
    bool check_are_mesh2D_polygon_half_edges_equal
    (const MESHA_TYPE & meshA, const MESHB_TYPE & meshB,
     IJK::ERROR & error)
    {
      typedef typename MESHA_TYPE::VERTEX_INDEX_TYPE
        VERTEX_INDEX_TYPEA;
      typedef typename MESHB_TYPE::VERTEX_INDEX_TYPE
        VERTEX_INDEX_TYPEB;
      typedef typename MESHA_TYPE::HALF_EDGE_INDEX_TYPE
        HALF_EDGE_INDEX_TYPEA;
      typedef typename MESHB_TYPE::HALF_EDGE_INDEX_TYPE
        HALF_EDGE_INDEX_TYPEB;
      typedef typename MESHA_TYPE::POLYGON_INDEX_TYPE
	POLYGON_INDEX_TYPE;
      typedef typename MESHA_TYPE::NUMBER_TYPE
	NUMBER_TYPE;

      // Check that meshes have same number of polygons.
      if (!check_are_num_mesh2D_polygons_equal(meshA, meshB, error)) 
	{ return false; }

      // Check that same polygons have same number of edges.
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < meshA.NumPolygons(); 
	   ipoly++) {

	if (meshA.IsPolygonDeleted(ipoly) &&
	    meshB.IsPolygonDeleted(ipoly)) {
	  // Skip deleted polygons.
	  continue;
	}

        if (meshA.NumPolygonEdges(ipoly) != meshB.NumPolygonEdges(ipoly)) {
          error.AddMessage("Mesh polygons have different number of edges.");
          error.AddMessage
            ("  In meshA, polygon ", ipoly, " has ", 
	     meshA.NumPolygonEdges(ipoly), " edges.");
          error.AddMessage
            ("  In meshB, polygon ", ipoly, " has ", 
	     meshB.NumPolygonEdges(ipoly), " edges.");
          return false;
        }
      }

      // Check that same polygons have same half edges.
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < meshA.NumPolygons(); 
	   ipoly++) {

	if (meshA.IsPolygonDeleted(ipoly) &&
	    meshB.IsPolygonDeleted(ipoly)) {
	  // Skip deleted polygons.
	  continue;
	}

	for (NUMBER_TYPE jloc = 0; jloc < meshA.NumPolygonVertices(ipoly); 
	     jloc++) {

	  const HALF_EDGE_INDEX_TYPEA jhalf_edgeA =
	    meshA.HalfEdgeIndex(ipoly, jloc);
	  const HALF_EDGE_INDEX_TYPEB jhalf_edgeB =
	    meshB.HalfEdgeIndex(ipoly, jloc);

	  if (!are_mesh2D_half_edges_equal
	      (meshA, jhalf_edgeA, meshB, jhalf_edgeB)) {

	    const VERTEX_INDEX_TYPEA jfromA =
	      meshA.FromVertexIndex(jhalf_edgeA);
	    const VERTEX_INDEX_TYPEA jfromB =
	      meshB.FromVertexIndex(jhalf_edgeB);
	    const VERTEX_INDEX_TYPEA jtoA =
	      meshA.ToVertexIndex(jhalf_edgeA);
	    const VERTEX_INDEX_TYPEA jtoB =
	      meshB.ToVertexIndex(jhalf_edgeB);

	    error.AddMessage("Mesh polygons have different half edges.");
	    error.AddMessage
	      ("  Mesh A, polygon ", ipoly, ", location ", jloc, 
	       ", has half edge ", jhalf_edgeA, 
	       " (", jfromA, ",", jtoA, ").");
	    error.AddMessage
	      ("  Mesh B, polygon ", ipoly, ", location ", jloc, 
	       ", has half edge ", jhalf_edgeB, 
	       " (", jfromB, ",", jtoB, ").");
	    return false;
	  }
	}
      }


      // Check that links (next half edge around edge) between polygons
      //   are the same.
      for (POLYGON_INDEX_TYPE ipoly = 0; ipoly < meshA.NumPolygons(); 
	   ipoly++) {

	if (meshA.IsPolygonDeleted(ipoly) &&
	    meshB.IsPolygonDeleted(ipoly)) {
	  // Skip deleted polygons.
	  continue;
	}

	for (NUMBER_TYPE jloc = 0; jloc < meshA.NumPolygonVertices(ipoly); 
	     jloc++) {

	  const HALF_EDGE_INDEX_TYPEA jhalf_edgeA =
	    meshA.HalfEdgeIndex(ipoly, jloc);
	  const HALF_EDGE_INDEX_TYPEB jhalf_edgeB =
	    meshB.HalfEdgeIndex(ipoly, jloc);
	  const HALF_EDGE_INDEX_TYPEA jnext_around_edgeA =
	    meshA.IndexOfNextHalfEdgeAroundEdge(jhalf_edgeA);
	  const HALF_EDGE_INDEX_TYPEA jnext_around_edgeB =
	    meshA.IndexOfNextHalfEdgeAroundEdge(jhalf_edgeB);

	  if (!are_mesh2D_half_edges_equal
	      (meshA, jnext_around_edgeA, 
	       meshB, jnext_around_edgeB)) {

	    const VERTEX_INDEX_TYPEA jfromA0 =
	      meshA.FromVertexIndex(jhalf_edgeA);
	    const VERTEX_INDEX_TYPEA jfromB0 =
	      meshB.FromVertexIndex(jhalf_edgeB);
	    const VERTEX_INDEX_TYPEA jtoA0 =
	      meshA.ToVertexIndex(jhalf_edgeA);
	    const VERTEX_INDEX_TYPEA jtoB0 =
	      meshB.ToVertexIndex(jhalf_edgeB);
	    const VERTEX_INDEX_TYPEA jfromA1 =
	      meshA.FromVertexIndex(jnext_around_edgeA);
	    const VERTEX_INDEX_TYPEA jfromB1 =
	      meshB.FromVertexIndex(jnext_around_edgeB);
	    const VERTEX_INDEX_TYPEA jtoA1 =
	      meshA.ToVertexIndex(jnext_around_edgeA);
	    const VERTEX_INDEX_TYPEA jtoB1 =
	      meshB.ToVertexIndex(jnext_around_edgeB);

	    error.AddMessage
	      ("Links between mesh polygons (next edge around edge) are different.");
	    error.AddMessage
	      ("  Mesh A, half edge ", jhalf_edgeA,
	       " (", jfromA0, ",", jtoA0, "),"
	       " next half edge around edge: ", jnext_around_edgeA,
	       " (", jfromA1, ",", jtoA1, ").");
	    error.AddMessage
	      ("  Mesh B, half edge ", jhalf_edgeB,
	       " (", jfromB0, ",", jtoB0, "),"
	       " next half edge around edge: ", jnext_around_edgeB,
	       " (", jfromB1, ",", jtoB1, ").");
	    return(false); 
	  }
	}
      }

      return true;
    }


    /*!
     *  @brief Return true if all mesh and meshB have same number 
     *    of vertices, edges and polygons (including deleted polygons) 
     *    and all polygon vertices and half edges are the same.
     */
    template <typename MESHA_TYPE, typename MESHB_TYPE>
    bool check_are_mesh2D_meshes_equal
    (const MESHA_TYPE & meshA, const MESHB_TYPE & meshB,
     IJK::ERROR & error)
    {
      typedef typename MESHA_TYPE::HALF_EDGE_INDEX_TYPE
        HALF_EDGE_INDEX_TYPEA;
      typedef typename MESHB_TYPE::HALF_EDGE_INDEX_TYPE
        HALF_EDGE_INDEX_TYPEB;

      if (!check_are_mesh2D_polygons_equal(meshA, meshB, error))
	{ return false; }

      if (!check_are_mesh2D_polygon_half_edges_equal(meshA, meshB, error))
	{ return false; }

      return true;
    }

  }


  #endif
