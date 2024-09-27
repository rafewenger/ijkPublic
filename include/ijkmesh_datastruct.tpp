/*!
 *  @file ijkmesh_datastruct.tpp
 *  @brief ijk template classes for polyhedral mesh data structures.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2010-2024 Rephael Wenger

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

#ifndef _IJKMESH_DATASTRUCT_
#define _IJKMESH_DATASTRUCT_

#include "ijk.tpp"
#include "ijklist.tpp"
#include "ijkmesh.tpp"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>


namespace IJK {

  // *****************************************************************
  // Class POLYMESH 
  // *****************************************************************

  /*!
   *  @brief Mesh of polytopes.
   *  @details
   *  - Mesh is represented as lists of polytopes.
   *  - Polytopes are represented as lists of vertices.
   *  - Number of vertices and maximum vertex index are NOT represented.
   */
  template <typename VTYPE, typename NTYPE>
  class POLYMESH:public LIST_OF_LISTS<VTYPE,NTYPE> {

  public:

    /// Vertex index type.
    typedef VTYPE VERTEX_INDEX_TYPE;

    /// NUMBER type.
    typedef NTYPE NUMBER_TYPE;

  public:
    /// constructor
    POLYMESH(){};

    /// Number of polytopes.
    NTYPE NumPoly() const               
    { return(this->NumLists()); }

    /// Number of vertices of polytope i.
    NTYPE NumPolyVert(const NTYPE ipoly) const
    { return(this->ListLength(ipoly)); }

    /// Return j'th vertex of polytope ipoly.
    VTYPE Vertex(const NTYPE ipoly, const NTYPE j) const
    { return(this->Element(ipoly,j)); }

    /// Return pointer to first vertex of polytope ipoly.
    const VTYPE * VertexList(const NTYPE ipoly) const
    { return(this->List(ipoly)); }

    /// Return true if polytope ipoly contains vertex iv.
    template <typename VTYPE2>
    bool DoesPolyContainVertex(const NTYPE ipoly, const VTYPE2 iv) const
    { return(this->DoesListContain(ipoly, iv)); }

    /// @brief Return true if polytope ipoly contains vertex iv.
    /// - Version which returns location of vertex iv in list.
    template <typename VTYPE2, typename ITYPE>
    bool DoesPolyContainVertex
    (const NTYPE ipoly, const VTYPE2 iv, ITYPE & iloc) const
    { return(this->DoesListContain(ipoly, iv, iloc)); }

    /// @brief Return true if 2D polygon ipoly contains edge (iv0,iv1).
    /// @pre Polygon vertices are listed in cylic order around the polygon.
    template <typename VTYPE0, typename VTYPE1>
    bool DoesPoly2DContainEdge
    (const NTYPE ipoly, const VTYPE0 iv0, const VTYPE1 iv1) const;

    /// Return true if polytopes ipoly0 and ipoly1 share one or more vertices.
    bool DoPolytopesShareVertex(const NTYPE ipoly0, const NTYPE ipoly1) const;

    /// @brief Return true if ipolyA equals ipolyB.
    /// - Polytopes are only equal if vertices are listed 
    ///   in the same order in each polytope.
    bool ArePolytopesEqual(const NTYPE ipolyA, const NTYPE ipolyB) const
    { return(this->AreListsEqual(ipolyA, ipolyB)); }

    /// Count number of polytopes with given number of vertices.
    template <typename NTYPE2>
    NTYPE CountNumPolytopesOfSize(const NTYPE2 n) const;

    /// Count number of polytopes of size greater than or equal to n.
    template <typename NTYPE2>    
    NTYPE CountNumPolytopesOfSizeGE(const NTYPE2 n) const;

    /// Count number of triangles.
    NTYPE CountNumOfTriangles() const
    { return CountNumPolytopesOfSize(3); }

    /// Count number of quadrilaterals.
    NTYPE CountNumOfQuads() const
    { return CountNumPolytopesOfSize(4); }

    /// Count number of pentagons.
    NTYPE CountNumOfPentagons() const
    { return CountNumPolytopesOfSize(5); }

    /// Return string of all vertices of polytope ipoly.
    std::string GetVerticesStr
    (const NTYPE ipoly, const std::string & separator) const
    { return this->GetListStr(ipoly, separator); }

    /// @brief Return string of all vertices of polytope ipoly
    ///   bracketed by parentheses.
    std::string GetParenthesizedVerticesStr
    (const NTYPE ipoly, const std::string & separator) const;

    /// @brief Add a polytope with list_length vertices.
    template <typename VTYPE2, typename NTYPE2>
    void AddPolytope(const VTYPE2 poly_vert_list[],
                     const NTYPE2 list_length)
    { this->AddList(poly_vert_list, list_length); }

    /*!
     *  @overload
     *  @brief Add a polytope. (C++ STL vector)
     *  - C++ STL vector type for array poly_vert_list[].
     *  @param poly_vert_list List of polytope vertices 
     *    in C++ STL vector.
     */
    template <typename VTYPE2>
    void AddPolytope(const std::vector<VTYPE2> & poly_vert_list)
    { AddPolytope(&poly_vert_list.front(), poly_vert_list.size()); }

    /*!
     *  @brief Add polytopes from list where each polytope has 
     *     num_vert_per_poly vertices.
     *  @param num_vert Total number of vertices in poly_vert_list[].
     *  @param num_vert_per_poly Number of vertices in each polytope.
     */
    template <typename VTYPE2, typename NTYPE1, typename NTYPE2>
    void AddPolytopes(const VTYPE2 poly_vert_list[],
                      const NTYPE1 num_vert,
                      const NTYPE2 num_vert_per_poly)
    { this->AddLists(poly_vert_list, num_vert, num_vert_per_poly); }

    /*!
     *  @brief Add polytopes from list where each polytope has 
     *    num_vert_per_poly vertices.
     * @param num_vert_per_poly Number of vertices in each polytope.
     */
    template <typename VTYPE2>
    void AddPolytopes(const std::vector<VTYPE2> & poly_vert_list, 
                      const NTYPE num_vert_per_poly)
    { this->AddLists(poly_vert_list, num_vert_per_poly); }

    /*!
     *  @brief Sort vertex coordinates in lexicographic order
     *    and relabel polytope vertices to correspond to new order.
     *    @pre All mesh vertices must be in the range [0..(num_vert-1)].
     *  @param dimension Dimension. Number of coordinates per vertex.
     *  @param vertex_coord[] Array of vertex coordinates.
     *    vertex_coord[i*dimension+j] is j'th coordinate of vertex i.
     *    - Modified by this routine.
     *  @param num_vert Number of vertices.
     *  @param[out] loc_map[i] After sorting, coordinates of vertex i 
     *    are at vertex_coord[loc_map[i]*dimension,loc_map[i]*dimension+1,...].
     *    @pre Array loc_map[] is preallocated to size at least num_points.
     */
    template <typename CTYPE, typename NTYPE2, typename LTYPE>
    void SortVertexCoordinates
    (const int dimension, CTYPE * vertex_coord, const NTYPE2 num_vert,
     LTYPE * loc_map);

    /*!
     *  @brief Sort vertex coordinates in lexicographic order
     *    and relabel polytope vertices to correspond to new order. 
     *    (C++ STL vector.)
     *  - Version using C++ STL vector for loc_map[].
     */
    template <typename CTYPE, typename NTYPE2, typename LTYPE>
    void SortVertexCoordinates
    (const int dimension, CTYPE * vertex_coord, const NTYPE2 num_vert,
     std::vector<LTYPE> & loc_map);

    /*!
     *  @overload
     *  @brief Sort vertex coordinates in lexicographic order
     *    and relabel polytope vertices to correspond to new order. 
     *  - Version that does not return loc_map[].
     */
    template <typename CTYPE, typename NTYPE2>
    void SortVertexCoordinates
    (const int dimension, CTYPE * vertex_coord, const NTYPE2 num_vert);

    /*!
     *  @overload
     *  @brief Sort vertex coordinates in lexicographic order
     *    and relabel polytope vertices to correspond to new order. 
     *  - Version that does not return loc_map[].
     *  - Version using C++ STL vector for vertex_coord[].
     */
    template <typename CTYPE>
    void SortVertexCoordinates
    (const int dimension, std::vector<CTYPE> & vertex_coord);

    /*!
     *  @brief Sort the vertices of each polytope.
     *  @details 
     *  - Note: Sorting may change simplex orientation.
     */
    void SortPolyVert()
    { this->SortEachList(); }

    /*!
     *  @brief Return list of polytope indices in sorted order.
     *  @param sorted_poly[i] = Index of i'th poly in sorted order.
     *  @pre SortPolyVert() should be called before GetSortedPolytopeIndices.
     */
    template <typename NTYPE2>
    void GetSortedPolytopeIndices(std::vector<NTYPE2> & sorted_poly) const
    { this->GetSortedListIndices(sorted_poly); }

    /*!
     *  @brief Append triangle vertices and triangle indices to lists.
     *  - Append vertices of triangles in mesh to triangle_vertex_list[].
     *  - Append corresponding triangle indices to triangle_index[].
     *  - Ignores/skips any polygons that are not triangles.
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
     *  @brief Get sorted list of sorted triangle vertices 
     *    and corresponding triangle indices.
     *  @details
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


    // Print routines (mainly for debugging.)

    /*!
       @brief Print polygon vertices.
       - Extended version.  Define left/right delimiters and separator.
       - Mainly for debugging.
       @param c0 Left delimiter.
       @param c1 Separator.
       @param c2 Right delimiter.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void PrintPolygonVerticesX
    (OSTREAM_TYPE & out, const ITYPE2 ipoly,
     const char c0, const char c1, const char c2) const;

    /*!
       @brief Print polygon vertices.
              (Default delimiters and separator.)
       - Mainly for debugging.
       - Version using delimiters '(' and ')' and separator ','.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void PrintPolygonVertices
    (OSTREAM_TYPE & out, const ITYPE2 ipoly) const
    { PrintPolygonVerticesX(out, ipoly, '(', ',', ')'); }

    /*!
       @brief Print polygon vertices. (Set prefix and suffix.)
       - Mainly for debugging.
       - Version adding prefix and suffix strings.
       @param s0 Prefix string.
       @param s1 Suffix string.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void PrintPolygonVertices
    (OSTREAM_TYPE & out, const STYPE0 & s0, 
     const ITYPE2 ipoly, const STYPE1 & s1) const;

    /*!
       @brief Print polygon index and vertices. (Extended version.)
       - Extended version.  Define left/right delimiters and separator.
       - Mainly for debugging.
       @param c0 Left delimiter.
       @param c1 Separator.
       @param c2 Right delimiter.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void PrintPolygonIndexAndVerticesX
    (OSTREAM_TYPE & out, const ITYPE2 ipoly,
     const char c0, const char c1, const char c2) const;

    /*!
       @brief Print polygon index and vertices.
              (Default delimiters and separator.)
       - Mainly for debugging.
       - Version using delimiters '(' and ')' and separator ','.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void PrintPolygonIndexAndVertices
    (OSTREAM_TYPE & out, const ITYPE2 ipoly) const
    { PrintPolygonIndexAndVerticesX(out, ipoly, '(', ',', ')'); }

    /*!
       @brief Print polygon index and vertices. (Set prefix and suffix.)
       - Mainly for debugging.
       - Version adding prefix and suffix strings.
       @param s0 Prefix string.
       @param s1 Suffix string.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void PrintPolygonIndexAndVertices
    (OSTREAM_TYPE & out, const STYPE0 & s0, 
     const ITYPE2 ipoly, const STYPE1 & s1) const;

  };


  // *****************************************************************
  // Class POLYMESH_DATA
  // *****************************************************************

  /// Mesh of polytopes with data for each polytope.
  template <typename VTYPE, typename NTYPE, typename POLYDATA_TYPE>
  class POLYMESH_DATA:public POLYMESH<VTYPE,NTYPE> {

  public:
    typedef typename POLYMESH<VTYPE,NTYPE>::NUMBER_TYPE NUMBER_TYPE;

  protected:
    void _ResizePolyData(const NTYPE k)
    { const NTYPE num_poly = this->NumPoly(); poly_data.resize(num_poly+k); }

  public:
    std::vector<POLYDATA_TYPE> poly_data;

  public:

    /// constructor
    POLYMESH_DATA(){};

    /// Add a polytope with list_length vertices.
    template <typename VTYPE2, typename NTYPE2>
    void AddPolytope(const VTYPE2 poly_vert_list[],
                     const NTYPE2 list_length)
    { _ResizePolyData(1); 
      POLYMESH<VTYPE,NTYPE>::AddPolytope(poly_vert_list, list_length); }

    /// Add a polytope.
    /// - C++ STL vector type for array poly_vert_list[].
    /// @param poly_vert_list List of polytope vertices in C++ STL vector.
    template <typename VTYPE2>
    void AddPolytope(const std::vector<VTYPE2> & poly_vert_list)
    { _ResizePolyData(1); 
      POLYMESH<VTYPE,NTYPE>::AddPolytope(poly_vert_list); }

    /// Add polytopes from list where each polytope has 
    ///   num_vert_per_poly vertices.
    /// @param num_vert Total number of vertices in poly_vert_list[].
    /// @param num_vert_per_poly Number of vertices in each polytope.
    template <typename VTYPE2, typename NTYPE1, typename NTYPE2>
    void AddPolytopes(const VTYPE2 poly_vert_list[],
                      const NTYPE1 num_vert,
                      const NTYPE2 num_vert_per_poly);

    /// Add polytopes from list where each polytope has 
    ///   num_vert_per_poly vertices.
    /// @param num_vert_per_poly Number of vertices in each polytope.
    template <typename VTYPE2>
    void AddPolytopes(const std::vector<VTYPE2> & poly_vert_list, 
                      const NTYPE num_vert_per_poly);

    /// @brief Resize poly data to this->NumPoly().
    void ResizePolyData()
    { poly_data.resize(this->NumPoly()); }

    /// Set polymesh data from polymesh.
    template <typename VTYPE2, typename NTYPE2>
    void Set(const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    // Check functions.

    /// Check data structure.
    /// - Return true if data structure passes check.
    bool Check(IJK::ERROR & error);
  };


  // *****************************************************************
  // Class VERTEX_POLY_INCIDENCE_ELEMENT
  // *****************************************************************

  /// Base element of VERTEX_POLY_INCIDENCE_BASE
  template <typename PTYPE>
  class VERTEX_POLY_INCIDENCE_ELEMENT {

  protected:
    PTYPE poly_index;

  public:
    typedef PTYPE POLY_INDEX_TYPE;

  public:
    PTYPE PolyIndex() const { return(poly_index); };
    void SetPolyIndex(const PTYPE ipoly) { poly_index = ipoly; }
  };


  // *****************************************************************
  // Class VERTEX_POLY_INCIDENCE_ELEMENT_WITH_VLOC
  // *****************************************************************

  /// VERTEX_POLY_INCIDENCE_ELEMENT with vertex location
  ///   in list of vertices incident on polytope poly_index.
  template <typename PTYPE, typename LTYPE>
  class VERTEX_POLY_INCIDENCE_ELEMENT_WITH_VLOC:
    public VERTEX_POLY_INCIDENCE_ELEMENT<PTYPE> {

  protected:

    /// Location of vertex in polytope ipoly list of vertices.
    LTYPE vertex_location;

  public:
    typedef LTYPE VERTEX_LOCATION_TYPE;

  public:
    LTYPE VertexLocation() const { return(vertex_location); };
    void SetVertexLocation(const LTYPE iloc)
    { vertex_location = iloc; }
  };


  // *****************************************************************
  // Class VERTEX_POLY_INCIDENCE_BASE and VERTEX_POLY_INCIDENCE.
  // *****************************************************************

  /// \brief List of polytopes incident on each vertex. Base class.
  /// @tparam ETYPE List element type.
  ///    Usually derived from VERTEX_POLY_INCIDENCE_ELEMENT.
  /// @tparam NTYPE Number type.
  template <typename ETYPE, typename NTYPE>
  class VERTEX_POLY_INCIDENCE_BASE:protected LIST_OF_LISTS<ETYPE,NTYPE> {

  public:

    typedef typename ETYPE::POLY_INDEX_TYPE POLY_INDEX_TYPE;

  protected:

    /// Number of vertices.
    NTYPE num_vertices;

    /// Number of polytopes.
    NTYPE num_poly;

  protected:

    typedef typename ETYPE::POLY_INDEX_TYPE PTYPE;

    /// Compute number of polytopes incident on each vertex.
    /// Allocate arrays list_length[], first_element[] and element[].
    /// Set number of polytopes and number of vertices.
    /// Set arrays list_length[] and first_element[].
    template <typename VTYPE2, typename NTYPE2>
    void AllocateLists(const POLYMESH<VTYPE2,NTYPE2> & polymesh);

  public:

    // Constructor
    VERTEX_POLY_INCIDENCE_BASE()
    { Clear(); };

    template <typename VTYPE, typename NTYPE2>
    VERTEX_POLY_INCIDENCE_BASE
    (const POLYMESH<VTYPE,NTYPE2> & polymesh)
    { Set(polymesh); }

    NTYPE NumPoly() const
    { return(num_poly); }

    NTYPE NumVertices() const
    { return(num_vertices); }

    template <typename VTYPE>
    NTYPE NumIncidentPoly(const VTYPE iv) const
    { return(this->ListLength(iv)); }

    /// Return list element containing j'th poly incident on vertex iv.
    template <typename VTYPE, typename NTYPE2>
    const ETYPE & Element(const VTYPE iv, const NTYPE2 j) const
    { return(this->element[this->ElementIndex(iv,j)]); }

    /// Return j'th poly containing vertex iv.
    template <typename VTYPE, typename NTYPE2>
    PTYPE IncidentPoly(const VTYPE iv, const NTYPE2 j) const
    { return(Element(iv,j).PolyIndex()); }

    /// Return true if vertex iv is incident on polytope ipoly.
    /// @param[out] iloc Location (index) of ipoly in list
    ///   list of poly incident on iv.
    /// iloc is undefined if vertex iv is not incident on ipoly.
    template <typename VTYPE2, typename PTYPE2, typename NTYPE2>
    bool IsVertexIncidentOnPoly
    (const VTYPE2 iv, const PTYPE2 ipoly, NTYPE2 & iloc) const;

    /// Return true if vertex iv is incident on polytope ipoly.
    /// Version which does not return iloc.
    template <typename VTYPE2, typename PTYPE2>
    bool IsVertexIncidentOnPoly(const VTYPE2 iv, const PTYPE2 ipoly) const
    { NTYPE iloc; 
      return(IsVertexIncidentOnPoly(iv, ipoly, iloc)); }

    /// Get polytopes containing vertices iv0 and iv1.
    /// @param poly_list[] List of polytopes containing vertices iv0 and iv1.
    template <typename VTYPE0, typename VTYPE1, typename PTYPE2>
    void GetPolyContaining
    (const VTYPE0 iv0, const VTYPE1 iv1, 
     std::vector<PTYPE2> & poly_list) const;

    /// Set vertex poly incidence from polymesh.
    template <typename VTYPE2, typename NTYPE2>
    void Set(const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    /// Clear all lists.
    void Clear();
  };

  /// \brief List of polytopes incident on each vertex.
  /// @tparam PTYPE Polytope index type.
  /// @tparam NTYPE Number type.
  template <typename PTYPE, typename NTYPE>
  class VERTEX_POLY_INCIDENCE:
    public VERTEX_POLY_INCIDENCE_BASE
    <VERTEX_POLY_INCIDENCE_ELEMENT<PTYPE>,NTYPE> {

  public:
    VERTEX_POLY_INCIDENCE() {};

    template <typename VTYPE, typename NTYPE2>
    VERTEX_POLY_INCIDENCE
    (const POLYMESH<VTYPE,NTYPE2> & polymesh)
    { this->Set(polymesh); }
  };


  // *****************************************************************
  // Class VERTEX_POLY_INCIDENCE_WITH_VLOC_BASE 
  //   and VERTEX_POLY_INCIDENCE_WITH_VLOC
  // *****************************************************************

  /// \brief List of polytopes incident on each vertex. Base class.
  /// \brief - Includes location of each vertex in each polytope vertex list.
  /// @tparam ETYPE List element type.
  ///    Usually derived from VERTEX_POLY_INCIDENCE_WITH_VLOC_ELEMENT.
  ///    ETYPE must have member function SetVertexLocation()
  /// @tparam NTYPE Number type.
  template <typename ETYPE, typename NTYPE>
  class VERTEX_POLY_INCIDENCE_WITH_VLOC_BASE:
    public VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE> {

  protected:

    typedef typename ETYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename ETYPE::VERTEX_LOCATION_TYPE LTYPE;

  public:
    /// Return location of vertex iv in vertex list of poly PolyIndex(iv,j).
    template <typename VTYPE, typename NTYPE2>
    LTYPE VertexLocInPolyVertexList(const VTYPE iv, const NTYPE2 j) const
    { return(this->Element(iv,j).VertexLocation()); }

    /// Set vertex poly incidence from polymesh.
    /// - Redefine to also set location of each vertex
    ///   in each polytope vertex list.
    template <typename VTYPE2, typename NTYPE2>
    void Set(const POLYMESH<VTYPE2,NTYPE2> & polymesh);
  };


  /// \brief List of polytopes incident on each vertex.
  /// \brief - Includes location of each vertex in each polytope vertex list.
  /// @tparam PTYPE Polytope index type.
  /// @tparam NTYPE Number type.
  template <typename PTYPE, typename LTYPE, typename NTYPE>
  class VERTEX_POLY_INCIDENCE_WITH_VLOC:
    public VERTEX_POLY_INCIDENCE_WITH_VLOC_BASE
  <VERTEX_POLY_INCIDENCE_ELEMENT_WITH_VLOC<PTYPE,LTYPE>,NTYPE> {

  public:
    VERTEX_POLY_INCIDENCE_WITH_VLOC() {};

    template <typename VTYPE, typename NTYPE2>
    VERTEX_POLY_INCIDENCE_WITH_VLOC
    (const POLYMESH<VTYPE,NTYPE2> & polymesh)
    { this->Set(polymesh); }
  };


  // *****************************************************************
  // Class VERTEX_ADJACENCY_LIST_ELEMENT
  // *****************************************************************

  /// Base element of VERTEX_ADJACENCY_LIST_BASE
  template <typename VTYPE>
  class VERTEX_ADJACENCY_LIST_ELEMENT {

  protected:
    VTYPE vertex;

  public:
    typedef VTYPE VERTEX_INDEX_TYPE;

  public:
    VTYPE Vertex() const { return(vertex); };
    void SetVertex(const VTYPE iv) { vertex = iv; }
  };


  // *****************************************************************
  // Classes VERTEX_ADJACENCY_LIST_BASE and VERTEX_ADJACENCY_LIST.
  // *****************************************************************

  /// \brief List of vertices adjacent to each vertex. Base class.
  /// @tparam ETYPE List element type.
  ///   Usually derived from VERTEX_ADJACENCY_LIST_ELEMENT.
  /// @tparam NTYPE Number type.
  template <typename ETYPE, typename NTYPE>
  class VERTEX_ADJACENCY_LIST_BASE:protected LIST_OF_LISTS<ETYPE,NTYPE> {

  protected:

    /// Number of vertices.
    NTYPE num_vertices;

  protected:
    
    typedef typename ETYPE::VERTEX_INDEX_TYPE VTYPE;

    /// Allocate arrays list_length[], first_element[] and element[].
    /// Set number of vertices.
    /// Set arrays list_length[] and first_element[].
    template <typename NTYPE2>
    void AllocateLists(const std::vector<NTYPE2> & num_adjacent);


  public:

    // Constructor
    VERTEX_ADJACENCY_LIST_BASE() 
    { Clear(); };

    NTYPE NumVertices() const
    { return(num_vertices); }

    template <typename VTYPE2>
    NTYPE NumAdjacent(const VTYPE2 iv) const
    { return(this->ListLength(iv)); }

    /// Return list element containing j'th vertex adjacent to vertex iv.
    template <typename VTYPE, typename NTYPE2>
    const ETYPE & Element(const VTYPE iv, const NTYPE2 j) const
    { return(this->element[this->ElementIndex(iv,j)]); }

    /// Return j'th vertex adjacent to vertex iv.
    template <typename VTYPE2, typename NTYPE2>
    VTYPE AdjacentVertex(const VTYPE2 iv, const NTYPE2 j) const
    { return(Element(iv,j).Vertex()); }

    /// Return true if iv1 is adjacent to iv0.
    template <typename VTYPE0, typename VTYPE1>
    bool IsAdjacent(const VTYPE0 iv0, const VTYPE1 iv1) const;

    template <typename NTYPE2>
    void SetNumVertices(const NTYPE2 num_vertices);

    /// Set adjacency list from mesh of 2D polygons.
    /// @param polymesh Mesh of 2D polygons.
    ///   - Polygon vertices are listed in clockwise or counter-clockwise
    ///     order around the polygon.
    template <typename VTYPE2, typename NTYPE2>
    void SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    /// Set adjacency list from mesh of (hyper) cubes.
    /// @param polymesh Mesh of (hyper) cubes.
    ///   - Cube vertices are listed in order:
    ///     (0,0,0), (1,0,0), (0,1,0), (1,1,0), 
    ///     (0,0,1), (1,0,1), (0,1,1), (1,1,1).
    template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
    void SetFromMeshOfCubes
    (const POLYMESH<VTYPE2,NTYPE2> & polymesh,
     const CUBE_TYPE & cube);

    /// Clear all lists.
    void Clear();
  };

  /// \brief List of vertices adjacent to each vertex.
  /// @tparam VTYPE Vertex type.
  /// @tparam NTYPE Number type.
  template <typename VTYPE, typename NTYPE>
  class VERTEX_ADJACENCY_LIST:
    public VERTEX_ADJACENCY_LIST_BASE
    <VERTEX_ADJACENCY_LIST_ELEMENT<VTYPE>,NTYPE> {

  public:
    VERTEX_ADJACENCY_LIST() {};
  };


  // *****************************************************************
  // Class POLYMESH_ADJACENCY
  // *****************************************************************

  template <typename PTYPE>
  class ADJACENT_POLY {

  public:
    /// Polytope index type.
    typedef PTYPE POLY_INDEX_TYPE;

  public:
    /// Index of adjacent polytope.
    PTYPE poly;

    /// True if poly is set, i.e., poly is adjacent across corresponding facet.
    bool is_adjacent;

  public:
    ADJACENT_POLY()
    { is_adjacent = false; }

    template <typename PTYPE2>
    void SetPoly(const PTYPE2 poly)
    {
      this->poly = poly;
      is_adjacent = true;
    }
  };


  template <typename PTYPE, typename FTYPE>
  class ADJACENT_POLY_FACET:public ADJACENT_POLY<PTYPE> {

  public:
    /// Polytope facet index type.
    typedef FTYPE FACET_INDEX_TYPE;

  public:
    /// Index of facet of adjacent polytope.
    FTYPE facet;

    template <typename PTYPE2, typename FTYPE2>
    void SetPolyFacet(const PTYPE2 poly, const FTYPE2 facet)
    {
      this->poly = poly;
      this->facet = facet;
      this->is_adjacent = true;
    }

    // Undefine SetPoly
    template <typename PTYPE2>
    void SetPoly(const PTYPE2 poly);
  };

  /// Polytopes adjacent to each polytope facet.
  template <const int MAX_NUM_FACETS, typename ADJACENT_POLY_FACET_TYPE,
            typename NTYPE>
  class POLYMESH_ADJACENCY_BASE:
    public LIST_OF_LISTS<ADJACENT_POLY_FACET_TYPE,NTYPE> {

  public:

    /// Polytope index type.
    typedef typename ADJACENT_POLY_FACET_TYPE::POLY_INDEX_TYPE POLY_INDEX_TYPE;

    /// NUMBER type.
    typedef NTYPE NUMBER_TYPE;

  public:
    /// constructor
    POLYMESH_ADJACENCY_BASE(){};

    /// Number of polytopes.
    NTYPE NumPoly() const               
    { return(this->NumLists()); }

    /// Max number of facets.
    NTYPE MaxNumFacets() const
    { return(MAX_NUM_FACETS); }

    /// @brief Return true if some polytope shares facet jfacet with poly ipoly.
    template <typename PTYPE2, typename FTYPE>
    bool IsAdjacent(const PTYPE2 ipoly, const FTYPE jfacet) const
    { return(this->element[this->ElementIndex(ipoly, jfacet)].is_adjacent); }

    /// @brief Return j'th poly adjacent to polytope ipoly.
    /// @details
    /// - Undefined if IsAdjacent(ipoly, jfacet) is false.
    template <typename PTYPE2, typename FTYPE>
    POLY_INDEX_TYPE AdjacentPoly(const PTYPE2 ipoly, const FTYPE jfacet) const
    { return(this->element[this->ElementIndex(ipoly, jfacet)].poly); }

    /// @brief Return facet of j'th poly adjacent to polytope ipoly.
    /// @details
    /// - Undefined if IsAdjacent(ipoly, jfacet) is false.
    template <typename PTYPE2, typename FTYPE>
    POLY_INDEX_TYPE AdjacentPolyFacet
    (const PTYPE2 ipoly, const FTYPE jfacet) const
    { return(this->element[this->ElementIndex(ipoly, jfacet)].facet); }

    /*!
     *  @brief Set adjacency poly from mesh of (hyper) cubes.
     *  @pre Each facet is contained in at most two (hyper) cubes.
     *  @pre Each (hyper) cube facets has 2d distinct facets 
     *    (where d is cube dimension).
     *  @param polymesh Mesh of (hyper) cubes.
     *    - Cube vertices are listed in order:
     *      (0,0,0), (1,0,0), (0,1,0), (1,1,0), 
     *      (0,0,1), (1,0,1), (0,1,1), (1,1,1).
     */
    template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
    void SetFromMeshOfCubes
    (const POLYMESH<VTYPE2,NTYPE2> & polymesh,
     const CUBE_TYPE & cube);
  };

  /// Polytopes adjacent to each polytope facet.
  template <typename PTYPE, typename FTYPE, typename NTYPE>
  class HEXMESH_ADJACENCY:
    public POLYMESH_ADJACENCY_BASE<6, ADJACENT_POLY_FACET<PTYPE,FTYPE>, 
                                   NTYPE> 
  {
  };


  // *****************************************************************
  // Class POLYMESH_ADJACENCY_LISTS
  // *****************************************************************

  /*!
   *  @brief Lists of adjacent polytopes for each polytope in mesh.
   *  @details
   *  - Note: This is not derived from or related to POLYMESH_ADJACENCY.
   *  - This is a totally different class which stores polymesh adjacency
   *    as unordered lists with no facet information.
   */
  template <typename PTYPE, typename NTYPE>
  class POLYMESH_ADJACENCY_LISTS:public LIST_OF_LISTS<PTYPE,NTYPE> {

  public:

    /// Polytope index type.
    typedef PTYPE POLY_INDEX_TYPE;

    /// NUMBER type.
    typedef NTYPE NUMBER_TYPE;

  public:
    /// constructor
    POLYMESH_ADJACENCY_LISTS(){};

    /// Number of polytopes.
    NTYPE NumPoly() const               
    { return(this->NumLists()); }

    /// Number of polytopes adjacent to polytope i.
    NTYPE NumAdjacentPoly(const NTYPE ipoly) const
    { return(this->ListLength(ipoly)); }

    /// Return j'th poly adjacent to polytope ipoly.
    template <typename PTYPE2, typename NTYPE2>
    PTYPE AdjacentPoly(const PTYPE2 ipoly, const NTYPE2 j) const
    { return(this->Element(ipoly,j)); }

    /// Set adjacency poly from mesh of (hyper) cubes.
    /// @param polymesh Mesh of (hyper) cubes.
    ///   - Cube vertices are listed in order:
    ///     (0,0,0), (1,0,0), (0,1,0), (1,1,0), 
    ///     (0,0,1), (1,0,1), (0,1,1), (1,1,1).
    template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
    void SetFromMeshOfCubes
    (const POLYMESH<VTYPE2,NTYPE2> & polymesh,
     const CUBE_TYPE & cube);
  };


  // *****************************************************************
  // Class VERTEX_POLY_EDGE_INCIDENCE
  // *****************************************************************

  /// \brief List of edges incident on given polytope and vertex.
  /// 
  /// - List of polytopes incident on each vertex.
  /// - List of vertices adjacent to each vertex.
  /// - For each poly ipoly incident on vertex iv, list of edges
  ///     incident on ipoly and iv.
  /// @tparam VP_ELEMENT_TYPE Element type for vertex poly incidence list.
  ///   Class VP_ELEMENT_TYPE must include member functions SetPolyIndex(),
  ///     first_edge and num_edges and type POLY_INDEX_TYPE.
  /// @tparam VA_ELEMENT_TYPE Element type for vertex adjacent list.
  ///   Class VA_ELEMENT_TYPE must include member functions SetVertex()
  ///     and type VERTEX_INDEX_TYPE.
  template <typename VP_ELEMENT_TYPE, typename VA_ELEMENT_TYPE, typename NTYPE>
  class VERTEX_POLY_EDGE_INCIDENCE {

  public:
    typedef typename LIST_OF_LISTS<NTYPE,NTYPE>::SIZE_TYPE SIZE_TYPE;

  protected:
    LIST_OF_LISTS<VP_ELEMENT_TYPE,NTYPE> vertex_poly_incidence;
    LIST_OF_LISTS<VA_ELEMENT_TYPE,NTYPE> vertex_adjacency_list;

    /*!
     * @brief For each vertex iv and incident poly ipoly, list of edges incident on iv and ipoly.
     *  - incident_edge_list.Element(ia,ib)=t is a reference to the t'th element
     *    of some list in vertex_adjacency_list.
     *  - If (iv0,iv1) is incident on ipoly and
     *    vertex_poly_incidence.Element(iv0,j) = ipoly and
     *    vertex_adjacency_list.Element(iv0,k)=iv1,
     *    then edge_sublist.Element(IncidentEdgeListIndex(iv0,j),i)=k
     *    for some i.
     */
    LIST_OF_LISTS<NTYPE,NTYPE> incident_edge_list;

  protected:
    typedef typename VP_ELEMENT_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename VA_ELEMENT_TYPE::VERTEX_INDEX_TYPE VTYPE;

    /// Increment length of list of edges incident on vertex iv and
    ///   poly vertex_poly_incidence(iv,j).
    template <typename VTYPE2, typename NTYPE2>
    void IncrementIncidentEdgeListLength(const VTYPE2 iv, const NTYPE2 j)
    { incident_edge_list.list_length[IncidentEdgeListIndex(iv,j)]++; }

    /// Set k'th element of list of edges incident on vertex iv and
    ///   poly vertex_poly_incidence(iv,j) to x.
    template <typename VTYPE2, typename NTYPE2, typename NTYPE3,
              typename XTYPE>
    void SetIncidentEdgeListElement
    (const VTYPE2 iv, const NTYPE2 j, const NTYPE3 k, const XTYPE x)
    { 
      NTYPE list_index = IncidentEdgeListIndex(iv,j);
      incident_edge_list.ElementRef(list_index,k) = x;
    }


  public:
    VERTEX_POLY_EDGE_INCIDENCE() {};

    NTYPE NumVertices() const
    { return(vertex_poly_incidence.NumLists()); }

    template <typename VTYPE2, typename NTYPE2>
    void SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    template <typename VTYPE>
    NTYPE NumIncidentPoly(const VTYPE iv) const
    { return(vertex_poly_incidence.ListLength(iv)); }

    /// Return j'th poly containing vertex iv.
    template <typename VTYPE, typename NTYPE2>
    PTYPE IncidentPoly(const VTYPE iv, const NTYPE2 j) const
    { return(vertex_poly_incidence.ElementRefConst(iv,j).PolyIndex()); }

    template <typename VTYPE2>
    NTYPE NumAdjacent(const VTYPE2 iv) const
    { return(vertex_adjacency_list.ListLength(iv)); }

    /// Return j'th vertex adjacent to vertex iv.
    template <typename VTYPE2, typename NTYPE2>
    VTYPE AdjacentVertex(const VTYPE2 iv, const NTYPE2 j) const
    { return(vertex_adjacency_list.ElementRefConst(iv,j).Vertex()); }

    /// Return number of edges incident on vertex iv and 
    ///   poly IncidentPoly(iv,j).
    template <typename VTYPE2, typename NTYPE2>
    NTYPE NumIncidentEdges(const VTYPE2 iv, const NTYPE2 j) const
    { return(incident_edge_list.list_length[IncidentEdgeListIndex(iv,j)]); }; 

    /// Return k'th element in list of edges incident on vertex iv and
    ///   poly IncidentPoly(iv,j).
    template <typename VTYPE2, typename NTYPE2, typename NTYPE3>
    NTYPE IncidentEdge(const VTYPE2 iv, const NTYPE2 j, const NTYPE3 k) const
    { return(incident_edge_list.Element(IncidentEdgeListIndex(iv,j),k)); }; 

    /// Return index of list of edges incident on vertex iv and
    ///   poly vertex_poly_incidence(iv,j).
    template <typename VTYPE2, typename NTYPE2>
    NTYPE IncidentEdgeListIndex(const VTYPE2 iv, const NTYPE2 j) const
    { return(vertex_poly_incidence.ElementIndex(iv, j)); }

    /// Return index of k'th edge incident on vertex iv and iv
    ///   poly vertex_poly_incidence(iv,j).
    template <typename VTYPE2, typename NTYPE2, typename NTYPE3>
    NTYPE IncidentEdgeIndex
    (const VTYPE2 iv, const NTYPE2 j, const NTYPE3 k) const
    { return(incident_edge_list.ElementIndex
             (IncidentEdgeListIndex(iv, j), k)); }

    /// Return true if list of edges incident on vertex iv0 and 
    ///  poly vertex_poly_incidence(iv0,j) contains edge (iv0, iv1).
    template <typename VTYPE0, typename VTYPE1, typename NTYPE2>
    bool DoesIncidenceListContainEdge
    (const VTYPE0 iv0, const VTYPE1 iv1, const NTYPE2 j) const;

    void Clear();
  };


  // *****************************************************************
  // Class FACE_INFO_BASE
  // *****************************************************************

  /*!
   *  @brief Face information.
   *  @tparam PINDEX_TYPE Polytope index type.
   *  @tparam FINDEX_TYPE Face index type.
   */
  template <typename PINDEX_TYPE, typename FINDEX_TYPE>
  class FACE_INFO_BASE {

  public:
    typedef PINDEX_TYPE POLY_INDEX_TYPE;
    typedef FINDEX_TYPE FACE_INDEX_TYPE;

  public:
    PINDEX_TYPE poly_containing_face;
    FINDEX_TYPE face_index;

  public:
    FACE_INFO_BASE() {};
    template <typename PTYPE2, typename FTYPE2>
    FACE_INFO_BASE(const PTYPE2 ipoly, const FTYPE2 jface)
    { poly_containing_face = ipoly; face_index = jface; }
  };


  // *****************************************************************
  // Class FACE_LIST_BASE
  // *****************************************************************

  /*!
   *  @brief List of face information.
   *  - Note: A face which appears in k polytopes will usually appear
   *    k times in the face list.
   */
  template <typename VTYPE, typename NTYPE, typename FACE_INFO_TYPE>
  class FACE_LIST_BASE:public POLYMESH_DATA<VTYPE,NTYPE,FACE_INFO_TYPE>
  {

  public:
    typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE POLY_INDEX_TYPE;
    typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FACE_INDEX_TYPE;
    typedef typename POLYMESH_DATA<VTYPE,NTYPE,FACE_INFO_TYPE>::NUMBER_TYPE 
    NUMBER_TYPE;

  public:
    FACE_LIST_BASE(){};

    // Get functions.
    NUMBER_TYPE NumFaces() const
    { return(POLYMESH_DATA<VTYPE,NTYPE,FACE_INFO_TYPE>::NumPoly()); }

    NUMBER_TYPE NumFaceVert(const NUMBER_TYPE jf) const
    { return(POLYMESH_DATA<VTYPE,NTYPE,FACE_INFO_TYPE>::NumPolyVert(jf)); }

    template <typename FTYPE>
    POLY_INDEX_TYPE PolyContainingFace(const FTYPE jf) const
    { return(this->poly_data[jf].poly_containing_face); }

    template <typename FTYPE>
    POLY_INDEX_TYPE FaceIndex(const FTYPE jf) const
    { return(this->poly_data[jf].face_index); }


    // Disable.
    NUMBER_TYPE NumPoly() const;
    NUMBER_TYPE NumPolyVert() const;
  };


  // *****************************************************************
  // Class FACET_LIST_BASE (NOTE: FACET NOT FACE)
  // *****************************************************************

  /*!
   *  @brief List of facet information.
   *  @details
   *  - Note: Facet, not face, as in FACE_LIST_BASE.
   *  - Note: A facet which appears in k polytopes will usually appear
   *    k times in the facet list.
   */
  template <typename VTYPE, typename NTYPE, typename FACET_INFO_TYPE>
  class FACET_LIST_BASE:public FACE_LIST_BASE<VTYPE,NTYPE,FACET_INFO_TYPE>
  {

  public:
    typedef typename FACET_INFO_TYPE::POLY_INDEX_TYPE POLY_INDEX_TYPE;
    typedef typename FACET_INFO_TYPE::FACE_INDEX_TYPE FACET_INDEX_TYPE;
    typedef typename FACE_LIST_BASE<VTYPE,NTYPE,FACET_INFO_TYPE>::NUMBER_TYPE 
    NUMBER_TYPE;

  public:
    FACET_LIST_BASE() {};

    // Get functions.

    /// Return number of facets.
    NUMBER_TYPE NumFacets() const
    { return(this->NumFaces()); }

    /// Return number of facet vertices.
    NUMBER_TYPE NumFacetVert(const NUMBER_TYPE jf) const
    { return(this->NumFaceVert(jf)); }

    /// Return index of polytope containing facet.
    template <typename FTYPE>
    POLY_INDEX_TYPE PolyContainingFacet(const FTYPE jf) const
    { return(this->PolyContainingFace(jf)); }

    /// Return index of facet in polytope.
    template <typename FTYPE>
    FACET_INDEX_TYPE FacetIndex(const FTYPE jf) const
    { return(this->FaceIndex(jf)); }

    /// @brief Set facets from 2D mesh.
    /// @param polymesh Mesh of 2D polygons.
    template <typename VTYPE2, typename NTYPE2>
    void SetFrom2DMesh
    (const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    /*!
     *  @overload
     *  @brief Set facets from 2D mesh. (Skip some polygons.)
     *  @details
     *  - Version which skips some polygons.
     *  @param polymesh Mesh of 2D polygons.
     *  @param flag_skip_poly[] If flag_skip_poly[ipoly] is true,
     *    then skip polygon ipoly.
     */
    template <typename VTYPE2, typename NTYPE2>
    void SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh,
                       const std::vector<bool> & flag_skip_poly);

    /// @brief Set facets from mesh of simplices.
    /// @param polymesh Mesh of simplices.
    template <typename VTYPE2, typename NTYPE2>
    void SetFromMeshOfSimplices
    (const POLYMESH<VTYPE2,NTYPE2> & polymesh);

    /*!
     *  @brief Set facets from mesh of (hyper) cubes.
     *  @param polymesh Mesh of (hyper) cubes.
     *    - Cube vertices are listed in order:
     *      (0,0,0), (1,0,0), (0,1,0), (1,1,0), 
     *      (0,0,1), (1,0,1), (0,1,1), (1,1,1).
     *  @pre polymesh is mesh of cubes of same dimension as cube.
     */
    template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
    void SetFromMeshOfCubes
    (const POLYMESH<VTYPE2,NTYPE2> & polymesh, const CUBE_TYPE & cube);
  };


  // *****************************************************************
  // POLYMESH compare
  // *****************************************************************

  /// \brief Function class for comparing polytopes in POLYMESH.
  /// @pre Vertices for each poly are sorted.
  template <typename MTYPE> 
  class POLYMESH_LESS_THAN:public LIST_LESS_THAN<MTYPE> {

  public:
    POLYMESH_LESS_THAN(const MTYPE * pmesh):
      LIST_LESS_THAN<MTYPE>(pmesh)
    {};
  };


  // *****************************************************************
  // Class POLYMESH member functions
  // *****************************************************************

  template <typename VTYPE, typename NTYPE>
  std::string POLYMESH<VTYPE,NTYPE>::GetParenthesizedVerticesStr
    (const NTYPE ipoly, const std::string & separator) const
  {
    IJK::PROCEDURE_ERROR
      error("POLYMESH::GetParenthesizedVerticesStr()");
    
    if (ipoly >= NumPoly()) {
      if (NumPoly() == 0) {
        error.AddMessage("Programming error. No polytopes in mesh.");
        throw error;
      }
      else {
        error.AddMessage
          ("Programming error. No polytope ", ipoly, ".");
        error.AddMessage
          ("  All polytope indices in range: [0,", (NumPoly()-1), ".");
        throw error;
      }
    }

    return this->GetParenthesizedListStr(ipoly, separator);
  }

  
  // Return true if 2D polygon ipoly contains edge (iv0,iv1).
  template <typename VTYPE, typename NTYPE>
  template <typename VTYPE0, typename VTYPE1>
  bool POLYMESH<VTYPE,NTYPE>::DoesPoly2DContainEdge
  (const NTYPE ipoly, const VTYPE0 iv0, const VTYPE1 iv1) const
  {
    if (NumPolyVert(ipoly) <= 1) { return(false); }

    for (NTYPE j0 = 0; j0 < NumPolyVert(ipoly); j0++) {
      VTYPE jv0 = Vertex(ipoly, j0);
      NTYPE j1 = (j0+1)%NumPolyVert(ipoly);
      VTYPE jv1 = Vertex(ipoly, j1);

      if (jv0 == iv0 && jv1 == iv1) { return(true); }
      if (jv1 == iv0 && jv0 == iv1) { return(true); }
    }

    return(false);
  }


  // Return true if polytopes ipoly0 and ipoly1 share one or more vertices.
  template <typename VTYPE, typename NTYPE>
  bool POLYMESH<VTYPE,NTYPE>::
  DoPolytopesShareVertex(const NTYPE ipoly0, const NTYPE ipoly1) const
  {
    NTYPE ilocA, ilocB;

    return(do_lists_intersect
           (VertexList(ipoly0), NumPolyVert(ipoly0),
            VertexList(ipoly1), NumPolyVert(ipoly1),
            ilocA, ilocB));
  }


  // Count number of polytopes with given number of vertices.
  template <typename VTYPE, typename NTYPE>
  template <typename NTYPE2>
  NTYPE POLYMESH<VTYPE,NTYPE>::
  CountNumPolytopesOfSize(const NTYPE2 n) const
  {
    NTYPE kount = 0;
    for (NTYPE ipoly = 0; ipoly < NumPoly(); ipoly++) {
      if (NumPolyVert(ipoly) == n) { kount++; }
    }

    return kount;
  }


  // Count number of polytopes of size greater than or equal to n.
  template <typename VTYPE, typename NTYPE>
  template <typename NTYPE2>    
  NTYPE POLYMESH<VTYPE,NTYPE>::
  CountNumPolytopesOfSizeGE(const NTYPE2 n) const
  {
    NTYPE kount = 0;
    for (NTYPE ipoly = 0; ipoly < NumPoly(); ipoly++) {
      if (NumPolyVert(ipoly) >= n) { kount++; }
    }

    return kount;
  }


  // Sort vertex coordinates.
  template <typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename NTYPE2, typename LTYPE>
  void POLYMESH<VTYPE,NTYPE>::
  SortVertexCoordinates
  (const int dimension, CTYPE * vertex_coord, const NTYPE2 num_vert,
   LTYPE * loc_map)
  {
    IJK::sort_point_coordinates
      (dimension, vertex_coord, num_vert, loc_map);

    this->ReplaceElements(loc_map);
  }

  
  // Sort vertex coordinates. (C++ STL vector.)
  // - Version using C++ STL vector for loc_map[].
  template <typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename NTYPE2, typename LTYPE>
  void POLYMESH<VTYPE,NTYPE>::
  SortVertexCoordinates
  (const int dimension, CTYPE * vertex_coord, const NTYPE2 num_vert,
   std::vector<LTYPE> & loc_map)
  {
    SortVertexCoordinates
      (dimension, vertex_coord, num_vert, IJK::vector2pointer(loc_map));
  }
  
  
  // Sort vertex coordinates.
  // - Version that does not return loc_map[].
  template <typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename NTYPE2>
  void POLYMESH<VTYPE,NTYPE>::
  SortVertexCoordinates
  (const int dimension, CTYPE * vertex_coord, const NTYPE2 num_vert)
  {
    std::vector<NTYPE2> loc_map(num_vert);

    SortVertexCoordinates
      (dimension, vertex_coord, num_vert, loc_map);
  }

  
  // Sort vertex coordinates. (No loc_map[]. C++ STL vector.)
  // - Version using C++ STL vector for vertex_coord[].
  // - Version that does not return loc_map[].
  template <typename VTYPE, typename NTYPE>
  template <typename CTYPE>
  void POLYMESH<VTYPE,NTYPE>::
  SortVertexCoordinates
  (const int dimension, std::vector<CTYPE> & vertex_coord)
  {
    SortVertexCoordinates
      (dimension, IJK::vector2pointer(vertex_coord),
       vertex_coord.size()/dimension);
  }
  
  
  // Append triangle vertices and triangle indices to lists.
  // Append vertices of triangles in mesh to triangle_vertex_list[].
  // - Append corresponding triangle indices to triangle_index[].
  // - Ignores/skips any polygons that are not triangles.
  template <typename VTYPE, typename NTYPE>
  template <typename VTYPE2, typename PTYPE2>
  void POLYMESH<VTYPE,NTYPE>::
  AppendToListOfVerticesOfAllTriangles
  (std::vector<VTYPE2> & triangle_vertex_list,
   std::vector<PTYPE2> & triangle_index) const
  {
    const int NUM_VERTICES_PER_TRIANGLE(3);

    for (NUMBER_TYPE ipoly = 0; ipoly < NumPoly(); ipoly++) {

      if (NumVertices(ipoly) == NUM_VERTICES_PER_TRIANGLE) {
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
  template <typename VTYPE, typename NTYPE>
  template <typename VTYPE2, typename PTYPE2>
  void POLYMESH<VTYPE,NTYPE>::GetListOfVerticesOfAllTriangles
  (std::vector<VTYPE2> & triangle_vertex_list,
   std::vector<PTYPE2> & triangle_index) const
  {
    triangle_vertex_list.clear();
    triangle_index.clear();
    AppendToListOfVerticesOfAllTriangles
      (triangle_vertex_list, triangle_index);
  }


  // Sort triangle vertices.
  template <typename VTYPE, typename NTYPE>
  template <typename VTYPE2>
  void POLYMESH<VTYPE,NTYPE>::SortTriangleVertices
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


  // Print polygon vertices.
  // - Extended version.  Define left/right delimiters and separator.
  // @param c0 Left delimiter.
  // @param c1 Separator.
  // @param c2 Right delimiter.
  template <typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, typename ITYPE2>
  void POLYMESH<VTYPE,NTYPE>::PrintPolygonVerticesX
  (OSTREAM_TYPE & out, const ITYPE2 ipoly,
   const char c0, const char c1, const char c2) const
  {
    out << c0;

    if (NumPolyVert(ipoly) > 0) {
      out << Vertex(ipoly, 0);
      for (NTYPE i = 1; i < NumPolyVert(ipoly); i++) {
        out << c1 << Vertex(ipoly, i);
      }
    }

    out << c2;
  }


  // Print polygon vertices.
  template <typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, typename ITYPE2,
            typename STYPE0, typename STYPE1>
  void POLYMESH<VTYPE,NTYPE>::PrintPolygonVertices
    (OSTREAM_TYPE & out, const STYPE0 & s0, 
     const ITYPE2 ipoly, const STYPE1 & s1) const
  {
    out << s0;
    PrintPolygonVertices(out, ipoly);
    out << s1;
  }


  // Print polygon index and vertices.
  template <typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, typename ITYPE2>
  void POLYMESH<VTYPE,NTYPE>::PrintPolygonIndexAndVerticesX
    (OSTREAM_TYPE & out, const ITYPE2 ipoly,
     const char c0, const char c1, const char c2) const
  {
    out << ipoly << " ";
    PrintPolygonVerticesX(out, ipoly, c0, c1, c2);
  }


  // Print polygon index and vertices.
  template <typename VTYPE, typename NTYPE>

  template <typename OSTREAM_TYPE, typename ITYPE2,
            typename STYPE0, typename STYPE1>
  void POLYMESH<VTYPE,NTYPE>::PrintPolygonIndexAndVertices
    (OSTREAM_TYPE & out, const STYPE0 & s0, 
     const ITYPE2 ipoly, const STYPE1 & s1) const
  {
    out << s0;
    PrintPolygonIndexAndVertices(out, ipoly);
    out << s1;
  }


  // *****************************************************************
  // Class POLYMESH_DATA member functions
  // *****************************************************************

  template <typename VTYPE, typename NTYPE, typename POLY_DATA>
  template <typename VTYPE2, typename NTYPE1, typename NTYPE2>
  void POLYMESH_DATA<VTYPE,NTYPE,POLY_DATA>::
  AddPolytopes(const VTYPE2 poly_vert_list[],
               const NTYPE1 num_vert,
               const NTYPE2 num_vert_per_poly)
  {
    const NTYPE1 num_added_poly = num_vert/num_vert_per_poly;

    _ResizePolyData(num_added_poly);
    POLYMESH<VTYPE,NTYPE>::AddPolytopes
      (poly_vert_list, num_vert, num_vert_per_poly);
  }


  template <typename VTYPE, typename NTYPE, typename POLY_DATA>
  template <typename VTYPE2>
  void POLYMESH_DATA<VTYPE,NTYPE,POLY_DATA>::
  AddPolytopes(const std::vector<VTYPE2> & poly_vert_list, 
                    const NTYPE num_vert_per_poly)
  {
    const NTYPE num_added_poly = poly_vert_list.size()/num_vert_per_poly;

    _ResizePolyData(num_added_poly);
    POLYMESH<VTYPE,NTYPE>::AddPolytopes(poly_vert_list, num_vert_per_poly);
  }

  template <typename VTYPE, typename NTYPE, typename POLY_DATA>
  template <typename VTYPE2, typename NTYPE2>
  void POLYMESH_DATA<VTYPE,NTYPE,POLY_DATA>::
  Set(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      AddPolytope(polymesh.VertexList(ipoly), polymesh.NumPolyVert(ipoly));
    }
  }

  template <typename VTYPE, typename NTYPE, typename POLY_DATA>
  bool POLYMESH_DATA<VTYPE,NTYPE,POLY_DATA>::
  Check(IJK::ERROR & error)
  {
    if (poly_data.size() != this->NumPoly()) {
      error.AddMessage
        ("Programming error.  Incorrect size of vector POLY_DATA::poly_data.");
      error.AddMessage
        ("  POLY_DATA::poly_data.size() = ", poly_data.size(), "");
      error.AddMessage
        ("  Number of polytopes = ", this->NumPoly(), "");
      error.AddMessage
        ("  POLY_DATA::poly_data.size() should equal number of polytopes.");
      return(false);
    }

    return(true);
  }


  // *****************************************************************
  // Class VERTEX_POLY_INCIDENCE member functions
  // *****************************************************************

  // Return true if vertex iv is incident on polytope ipoly.
  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename PTYPE2, typename NTYPE2>
  bool VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::
  IsVertexIncidentOnPoly
  (const VTYPE2 iv, const PTYPE2 ipoly, NTYPE2 & iloc) const
  {
    const ETYPE * list = this->List(iv);
    for (iloc = 0; iloc < this->ListLength(iv); iloc++) {
      if (list[iloc].PolyIndex() == ipoly) { return(true); }
    }

    return(false);
  }

  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE0, typename VTYPE1, typename PTYPE2>
  void VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::
  GetPolyContaining
  (const VTYPE0 iv0, const VTYPE1 iv1, std::vector<PTYPE2> & poly_list) const
  {
    poly_list.clear();
    for (NTYPE j = 0; j < NumIncidentPoly(iv0);j++) {
      PTYPE ipoly = IncidentPoly(iv0, j);
      if (IsVertexIncidentOnPoly(iv1, ipoly)) 
        { poly_list.push_back(ipoly); }
    }
  }

  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::
  AllocateLists(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    this->num_poly = polymesh.NumPoly();
    this->num_vertices = 0;
    if (this->num_poly > 0) { this->num_vertices = get_max_vert(polymesh)+1; }
    this->SetNumLists(NumVertices());
    count_num_incident_poly(polymesh, this->list_length);
    this->SetFirstElement();
    this->AllocArrayElement();
  }


  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::
  Set(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    std::vector<PTYPE> current_poly_index;
    std::vector<PTYPE> last_poly;
    IJK::PROCEDURE_ERROR error("VERTEX_POLY_INCIDENCE_BASE::Set");

    Clear();
    AllocateLists(polymesh);

    last_poly.resize(NumVertices());
    current_poly_index.resize(NumVertices(), 0);
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      for (NTYPE i = 0; i < polymesh.NumPolyVert(ipoly); i++) {
        const VTYPE2 iv = polymesh.Vertex(ipoly, i);

        if (current_poly_index[iv] != 0) {
          if (last_poly[iv] == ipoly) {
            // Poly ipoly already inserted in list of poly incident on iv.
            continue;
          }
        }

        this->ElementRef(iv, current_poly_index[iv]).SetPolyIndex(ipoly);
        current_poly_index[iv]++;
        last_poly[iv] = ipoly;
      }
    }

    // Check stored correct number of poly for each vertex.
    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (current_poly_index[iv] != this->ListLength(iv)) {
        error.AddMessage
          ("Programming error.  Incorrect storage of poly vertices.");
        error.AddMessage
          ("  num_incident_poly[", iv, "] = ", this->ListLength(iv), ".");
        error.AddMessage
          ("  Stored ", current_poly_index[iv], " incident poly.");
        throw error;
      }
    }


  }

  template <typename ETYPE, typename NTYPE>
  void VERTEX_POLY_INCIDENCE_BASE<ETYPE,NTYPE>::Clear()
  {
    num_vertices = 0;
    num_poly = 0;
    LIST_OF_LISTS<ETYPE,NTYPE>::Clear();
  }


  // *****************************************************************
  // Class VERTEX_POLY_INCIDENCE_WITH_VLOC member functions
  // *****************************************************************

  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_POLY_INCIDENCE_WITH_VLOC_BASE<ETYPE,NTYPE>::
  Set(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    std::vector<PTYPE> current_poly_index;
    std::vector<PTYPE> last_poly;
    IJK::PROCEDURE_ERROR error("VERTEX_POLY_INCIDENCE_WITH_VLOC_BASE::Set");

    this->Clear();
    this->AllocateLists(polymesh);

    last_poly.resize(this->NumVertices());
    current_poly_index.resize(this->NumVertices(), 0);
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      for (NTYPE i = 0; i < polymesh.NumPolyVert(ipoly); i++) {
        const VTYPE2 iv = polymesh.Vertex(ipoly, i);

        if (current_poly_index[iv] != 0) {
          if (last_poly[iv] == ipoly) {
            // Poly ipoly already inserted in list of poly incident on iv.
            continue;
          }
        }

        this->ElementRef(iv, current_poly_index[iv]).SetPolyIndex(ipoly);
        this->ElementRef(iv, current_poly_index[iv]).SetVertexLocation(i);
        current_poly_index[iv]++;
        last_poly[iv] = ipoly;
      }
    }

    // Check stored correct number of poly for each vertex.
    for (NTYPE iv = 0; iv < this->NumVertices(); iv++) {
      if (current_poly_index[iv] != this->ListLength(iv)) {
        error.AddMessage
          ("Programming error.  Incorrect storage of poly vertices.");
        error.AddMessage
          ("  num_incident_poly[", iv, "] = ", this->ListLength(iv), ".");
        error.AddMessage
          ("  Stored ", current_poly_index[iv], " incident poly.");
        throw error;
      }
    }

  }


  // *****************************************************************
  // Class VERTEX_ADJACENCY_LIST_BASE member functions
  // *****************************************************************

  // Return true if iv1 is adjacent to iv0.
  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE0, typename VTYPE1>
  bool VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  IsAdjacent(const VTYPE0 iv0, const VTYPE1 iv1) const
  {
    const ETYPE * list = this->List(iv0);
    for (NTYPE i = 0; i < this->ListLength(iv0); i++) {
      if (list[i].Vertex() == iv1) { return(true); }
    }

    return(false);
  }


  template <typename ETYPE, typename NTYPE>
  template <typename NTYPE2>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  SetNumVertices(const NTYPE2 num_vertices)
  {
    this->num_vertices = num_vertices;
    this->SetNumLists(num_vertices);
  }


  // Allocate arrays list_length[], first_element[] and element[].
  // Set number of vertices.
  // Set arrays list_length[] and first_element[].
  template <typename ETYPE, typename NTYPE>
  template <typename NTYPE2>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  AllocateLists(const std::vector<NTYPE2> & num_adjacent)
  {
    num_vertices = num_adjacent.size();

    for (NTYPE i = 0; i < NumVertices(); i++) 
      { this->list_length[i] = num_adjacent[i]; }

    this->SetFirstElement();
    this->AllocArrayElement();
  }


  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    LIST_OF_LISTS<NTYPE,NTYPE> adjacent;
    std::vector<NTYPE> current_element;
    std::vector<NTYPE> num_adjacent;
    IJK::PROCEDURE_ERROR error("VERTEX_ADJACENCY_LIST::SetFrom2DMesh");
 
    Clear();

    if (polymesh.NumPoly() > 0) { 
      const int numv = get_max_vert(polymesh)+1; 
      SetNumVertices(numv);
    }

    if (NumVertices() == 0) { return; };

    adjacent.SetNumLists(NumVertices());
    count_num_poly2D_edges_incident_on_vertices
      (polymesh, adjacent.list_length);

    adjacent.SetFirstElement();
    adjacent.AllocArrayElement();

    current_element.resize(NumVertices());
    for (NTYPE iv = 0; iv < NumVertices(); iv++) 
      { current_element[iv] = adjacent.first_element[iv]; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);
      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          VTYPE2 iv0 = polymesh.Vertex(ipoly, i0);
          NTYPE i1 = (i0+1)%num_poly_vert;
          VTYPE2 iv1 = polymesh.Vertex(ipoly, i1);

          if (iv0 != iv1) {
            adjacent.element[current_element[iv0]] = iv1;
            current_element[iv0]++;
            adjacent.element[current_element[iv1]] = iv0;
            current_element[iv1]++;
          }
        }
      }
    }

    // Check stored correct number of adjacent vertices for each vertex.
    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (current_element[iv] != 
          adjacent.first_element[iv]+adjacent.list_length[iv]) {
        error.AddMessage
          ("Programming error.  Problem computing vertices adjacent to vertex ",
           iv, ".");
        error.AddMessage
          ("  Expected ", adjacent.list_length[iv], 
           " vertices (including duplicates)");
        error.AddMessage
          ("  but computed ", current_element[iv]-adjacent.first_element[iv],
           " adjacent vertices.");
        throw error;
      }
    }

    // Sort adjacent vertices.
    adjacent.SortEachList();

    // Set num_adjacent[i] to number of distinct vertices in adjacent list i.
    num_adjacent.resize(NumVertices());
    adjacent.CountNumDistinct(num_adjacent);

    AllocateLists(num_adjacent);

    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (this->list_length[iv] > 0) {
        NTYPE n = 0;
        NTYPE k = adjacent.first_element[iv];
        this->ElementRef(iv, n).SetVertex(adjacent.element[k]);
        n++;
        for (NTYPE j = 1; j < adjacent.list_length[iv]; j++) {
          if (adjacent.element[k+j] != adjacent.element[k+j-1]) {
            this->ElementRef(iv, n).SetVertex(adjacent.element[k+j]);
            n++;
          }
        }
        if (n != NumAdjacent(iv)) {
          error.AddMessage
            ("Programming error.  Problem computing vertices adjacent to vertex ",
             iv, ".");
          error.AddMessage
            ("  Expected ", NumAdjacent(iv), " vertices (no duplicates)");
          error.AddMessage
            ("  but computed ", n, " adjacent vertices.");
          throw error;
        }
      }
    }

  }


  template <typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::
  SetFromMeshOfCubes(const POLYMESH<VTYPE2,NTYPE2> & polymesh,
                     const CUBE_TYPE & cube)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = cube.Dimension();
    LIST_OF_LISTS<NTYPE,NTYPE> adjacent;
    std::vector<NTYPE> current_element;
    std::vector<NTYPE> num_adjacent;
    IJK::PROCEDURE_ERROR error("VERTEX_ADJACENCY_LIST::SetFromMeshOfCubes");
 
    Clear();

    if (polymesh.NumPoly() > 0) { 
      const int numv = get_max_vert(polymesh)+1; 
      SetNumVertices(numv);
    }

    if (NumVertices() == 0) { return; };

    adjacent.SetNumLists(NumVertices());
    count_num_cube_edges_incident_on_vertices
      (polymesh, cube, adjacent.list_length);

    adjacent.SetFirstElement();
    adjacent.AllocArrayElement();

    current_element.resize(NumVertices());
    for (NTYPE iv = 0; iv < NumVertices(); iv++) 
      { current_element[iv] = adjacent.first_element[iv]; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);
      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          const VTYPE iv0 = polymesh.Vertex(ipoly, i0);
          for (DTYPE d = 0; d < dimension; d++) {
            const NTYPE i1 = cube.VertexNeighbor(i0,d);
            const VTYPE iv1 = polymesh.Vertex(ipoly, i1);

            if (iv0 < iv1) {
              adjacent.element[current_element[iv0]] = iv1;
              current_element[iv0]++;
              adjacent.element[current_element[iv1]] = iv0;
              current_element[iv1]++;
            }
          }
        }
      }
    }

    // Check stored correct number of adjacent vertices for each vertex.
    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (current_element[iv] != 
          adjacent.first_element[iv]+adjacent.list_length[iv]) {
        error.AddMessage
          ("Programming error.  Problem computing vertices adjacent to vertex ",
           iv, ".");
        error.AddMessage
          ("  Expected ", adjacent.list_length[iv], 
           " vertices (including duplicates)");
        error.AddMessage
          ("  but computed ", current_element[iv]-adjacent.first_element[iv],
           " adjacent vertices.");
        throw error;
      }
    }

    // Sort adjacent vertices.
    adjacent.SortEachList();

    // Set num_adjacent[i] to number of distinct vertices in adjacent list i.
    num_adjacent.resize(NumVertices());
    adjacent.CountNumDistinct(num_adjacent);

    AllocateLists(num_adjacent);

    for (NTYPE iv = 0; iv < NumVertices(); iv++) {
      if (this->list_length[iv] > 0) {
        NTYPE n = 0;
        NTYPE k = adjacent.first_element[iv];
        this->ElementRef(iv, n).SetVertex(adjacent.element[k]);
        n++;
        for (NTYPE j = 1; j < adjacent.list_length[iv]; j++) {
          if (adjacent.element[k+j] != adjacent.element[k+j-1]) {
            this->ElementRef(iv, n).SetVertex(adjacent.element[k+j]);
            n++;
          }
        }
        if (n != NumAdjacent(iv)) {
          error.AddMessage
            ("Programming error.  Problem computing vertices adjacent to vertex ",
             iv, ".");
          error.AddMessage
            ("  Expected ", NumAdjacent(iv), " vertices (no duplicates)");
          error.AddMessage
            ("  but computed ", n, " adjacent vertices.");
          throw error;
        }
      }
    }

  }


  template <typename ETYPE, typename NTYPE>
  void VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE>::Clear()
  {
    num_vertices = 0;
    LIST_OF_LISTS<ETYPE,NTYPE>::Clear();
  }


  // *****************************************************************
  // Class POLYMESH_ADJACENCY_BASE member functions
  // *****************************************************************

  template <const int MAX_NUM_FACETS, typename ADJACENT_POLY_TYPE, 
            typename NTYPE>
  template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
  void POLYMESH_ADJACENCY_BASE<MAX_NUM_FACETS,ADJACENT_POLY_TYPE,NTYPE>::
  SetFromMeshOfCubes(const POLYMESH<VTYPE2,NTYPE2> & polymesh,
                     const CUBE_TYPE & cube)
  {
    typedef typename POLYMESH<VTYPE2,NTYPE2>::VERTEX_INDEX_TYPE VTYPE;
    typedef typename ADJACENT_POLY_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef FACE_INFO_BASE<PTYPE,NTYPE> FACET_INFO;

    const NTYPE num_poly = polymesh.NumPoly();
    const NTYPE num_cube_facets = cube.NumFacets();
    FACET_LIST_BASE<VTYPE,NTYPE,FACET_INFO> facets;
    std::vector<VTYPE> poly_facet(cube.NumFacetVertices());
    std::vector<PTYPE> sorted_facet;
    std::vector<PTYPE> adj_poly(num_poly*num_cube_facets);
    IJK::PROCEDURE_ERROR error("POLYMESH_ADJACENCY::SetFromMeshOfCubes");

    this->CreateUniformLengthLists(MAX_NUM_FACETS,num_poly);

    facets.SetFromMeshOfCubes(polymesh, cube);
    facets.SortPolyVert();
    facets.GetSortedPolytopeIndices(sorted_facet);

    NTYPE j0 = 0;
    while (j0+1 < facets.NumFacets()) {
      const NTYPE j1 = j0+1;

      // Note: This code assumes that at most two polytopes share a facet.
      if (facets.ArePolytopesEqual(sorted_facet[j0],sorted_facet[j1])) {

        const PTYPE jpoly0 = 
          facets.poly_data[sorted_facet[j0]].poly_containing_face;
        const NTYPE facet0 = facets.poly_data[sorted_facet[j0]].face_index;
        const PTYPE jpoly1 = 
          facets.poly_data[sorted_facet[j1]].poly_containing_face;
        const NTYPE facet1 = facets.poly_data[sorted_facet[j1]].face_index;

        if (jpoly0 != jpoly1) {
          this->element[this->ElementIndex(jpoly0,facet0)].
            SetPolyFacet(jpoly1, facet1);
          this->element[this->ElementIndex(jpoly1,facet1)].
            SetPolyFacet(jpoly0, facet0);
        }

        // Note: This code assumes that at most two polytopes share a facet.
        j0 += 2;
      } 
    else 
      { j0++; }
    }

  }


  // *****************************************************************
  // Class POLYMESH_ADJACENCY_LISTS member functions
  // *****************************************************************

  template <typename PTYPE, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
  void POLYMESH_ADJACENCY_LISTS<PTYPE,NTYPE>::
  SetFromMeshOfCubes(const POLYMESH<VTYPE2,NTYPE2> & polymesh,
                     const CUBE_TYPE & cube)
  {
    typedef typename POLYMESH<VTYPE2,NTYPE2>::VERTEX_INDEX_TYPE VTYPE;
    typedef FACE_INFO_BASE<PTYPE,NTYPE> FACET_INFO;

    const NTYPE num_poly = polymesh.NumPoly();
    const NTYPE num_cube_facets = cube.NumFacets();
    FACET_LIST_BASE<VTYPE,NTYPE,FACET_INFO> facets;
    std::vector<VTYPE> poly_facet(cube.NumFacetVertices());
    std::vector<PTYPE> sorted_facet;
    std::vector<PTYPE> adj_poly(num_poly*num_cube_facets);
    std::vector<std::pair<PTYPE,PTYPE>> overflow_adj_poly_list;
    std::vector<PTYPE> num_adj(num_poly, 0);
    IJK::PROCEDURE_ERROR error("POLYMESH_ADJACENCY_LISTS::SetFromMeshOfCubes");

    facets.SetFromMeshOfCubes(polymesh, cube);
    facets.SortPolyVert();
    facets.GetSortedPolytopeIndices(sorted_facet);

    NTYPE j0 = 0;
    while (j0 < facets.NumFacets()) {
      NTYPE j1 = j0+1;
      while (j1 < facets.NumFacets() && 
             facets.ArePolytopesEqual(sorted_facet[j0],sorted_facet[j1])) 
        { j1++; }

      const PTYPE jpoly0 = facets.PolyContainingFacet(sorted_facet[j0]);

      if (j0+2 == j1) {
        const PTYPE jpoly1 = facets.PolyContainingFacet(sorted_facet[j0+1]);

        if (jpoly0 != jpoly1) {
          IJK::add_to_list_store_overflow
            (jpoly0, jpoly1, &(adj_poly[jpoly0*num_cube_facets]),
             num_adj[jpoly0], num_cube_facets, overflow_adj_poly_list);
          IJK::add_to_list_store_overflow
            (jpoly1, jpoly0, &(adj_poly[jpoly1*num_cube_facets]),
             num_adj[jpoly1], num_cube_facets, overflow_adj_poly_list);
        }
      }
      else if (j0+2 < j1) {

        for (NTYPE k0 = j0; k0+1 < j1; k0++) {
          for (NTYPE k1 = k0+1; k1 < j1; k1++) {
            const PTYPE kpoly0 = facets.PolyContainingFacet(sorted_facet[k0]);
            const PTYPE kpoly1 = facets.PolyContainingFacet(sorted_facet[k1]);

            if (kpoly0 != kpoly1) {
              IJK::add_to_list_store_overflow
                (kpoly0, kpoly1, &(adj_poly[kpoly0*num_cube_facets]),
                 num_adj[kpoly0], num_cube_facets, overflow_adj_poly_list);
              IJK::add_to_list_store_overflow
                (kpoly1, kpoly0, &(adj_poly[kpoly1*num_cube_facets]),
                 num_adj[kpoly1], num_cube_facets, overflow_adj_poly_list);
            }
          }
        }
      }

      j0 = j1;
    }

    // Store adjacent polytopes.
    if (overflow_adj_poly_list.size() == 0) {
      for (PTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
        this->AddList(&(adj_poly[ipoly*num_cube_facets]), num_adj[ipoly]);
      }
    }
    else {
      std::vector<PTYPE> adj_poly_i(num_poly);

      std::sort(overflow_adj_poly_list.begin(), overflow_adj_poly_list.end());

      NTYPE k = 0;
      for (PTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
        NTYPE num_adj_i = num_adj[ipoly];
        std::copy(&(adj_poly[ipoly*num_cube_facets]),
                  &(adj_poly[ipoly*num_cube_facets+num_adj[ipoly]]),
                  &(adj_poly_i.front()));
        while (k < overflow_adj_poly_list.size() &&
               overflow_adj_poly_list[k].first == ipoly) {
          if (num_adj_i < adj_poly_i.size()) {
            adj_poly_i[num_adj_i] = overflow_adj_poly_list[k].second;
          }
          else {
            adj_poly_i.push_back(overflow_adj_poly_list[k].second);
          }
          num_adj_i++;
          k++;
        }

        this->AddList(&(adj_poly_i.front()), num_adj_i);
      }
    }
  }


  // *****************************************************************
  // Class VERTEX_POLY_EDGE_INCIDENCE member functions
  // *****************************************************************

  template <typename ETYPE0, typename ETYPE1, typename NTYPE>
  void VERTEX_POLY_EDGE_INCIDENCE<ETYPE0,ETYPE1,NTYPE>::Clear()
  {
    vertex_poly_incidence.Clear();
    vertex_adjacency_list.Clear();
    incident_edge_list.Clear();
  }


  template <typename ETYPE0, typename ETYPE1, typename NTYPE>
  template <typename VTYPE2, typename NTYPE2>
  void VERTEX_POLY_EDGE_INCIDENCE<ETYPE0,ETYPE1,NTYPE>::
  SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    // Vertex, incident poly, index of poly in vertex incident list.
    typedef std::tuple<VTYPE,PTYPE,NTYPE> VERTEX_POLY_TUPLE;

    LIST_OF_LISTS<VERTEX_POLY_TUPLE,NTYPE> adjacent;
    std::vector<NTYPE> current_element;
    std::vector<PTYPE> last_poly;
    std::vector<PTYPE> current_poly_index;
    std::vector<NTYPE> current_incident_edge_index;
    NTYPE numv, nump;
    IJK::PROCEDURE_ERROR error("VERTEX_POLY_EDGE_INCIDENCE::SetFrom2DMesh");

    Clear();

    if (polymesh.NumPoly() == 0) { return; }

    nump = polymesh.NumPoly();
    numv = 0;
    if (nump > 0) { numv = get_max_vert(polymesh)+1; }
    vertex_poly_incidence.SetNumLists(numv);
    count_num_incident_poly(polymesh, vertex_poly_incidence.list_length);
    vertex_poly_incidence.SetFirstElement();
    vertex_poly_incidence.AllocArrayElement();

    adjacent.SetNumLists(numv);
    count_num_poly2D_edges_incident_on_vertices(polymesh, adjacent.list_length);
    adjacent.SetFirstElement();
    adjacent.AllocArrayElement();

    current_poly_index.resize(numv, 0);
    current_element.resize(numv);

    for (NTYPE iv = 0; iv < numv; iv++) 
      { current_element[iv] = adjacent.FirstElement(iv); }


    last_poly.resize(numv);
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);

      // Insert ipoly in incident lists for its vertices.
      for (NTYPE i  = 0; i  < num_poly_vert; i ++) {
        const VTYPE2 iv = polymesh.Vertex(ipoly, i);

        if (current_poly_index[iv] != 0) {
          if (last_poly[iv] == ipoly) {
            // Poly ipoly already inserted in list of poly incident on iv.
            continue;
          }
        }

        const NTYPE k = 
          vertex_poly_incidence.ElementIndex(iv, current_poly_index[iv]);
        vertex_poly_incidence.element[k].SetPolyIndex(ipoly);
        current_poly_index[iv]++;
        last_poly[iv] = ipoly;
      }

      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          VTYPE2 iv0 = polymesh.Vertex(ipoly, i0);
          NTYPE i1 = (i0+1)%num_poly_vert;
          VTYPE2 iv1 = polymesh.Vertex(ipoly, i1);

          if (iv0 != iv1) {
            NTYPE k0 = current_poly_index[iv0]-1;
            NTYPE k1 = current_poly_index[iv1]-1;
            adjacent.element[current_element[iv0]] = 
              VERTEX_POLY_TUPLE(iv1,ipoly,k0);
            current_element[iv0]++;
            adjacent.element[current_element[iv1]] =
              VERTEX_POLY_TUPLE(iv0,ipoly,k1);
            current_element[iv1]++;
          }
        }
      }
    }

    // Check stored correct number of adjacent vertices for each vertex.
    for (NTYPE iv = 0; iv < numv; iv++) {
      if (current_element[iv] != 
          adjacent.FirstElement(iv) + adjacent.ListLength(iv)) {
        error.AddMessage
          ("Programming error.  Problem computing vertices adjacent to vertex ",
           iv, ".");
        error.AddMessage
          ("  Expected ", adjacent.ListLength(iv), 
           " vertices (including duplicates)");
        error.AddMessage
          ("  but computed ", current_element[iv]-adjacent.FirstElement(iv), 
           " adjacent vertices.");
        throw error;
      }
    }

    // Sort adjacent vertices.
    adjacent.SortEachList();

    vertex_adjacency_list.SetNumLists(numv);
    count_num_distinct_tuple0(adjacent, vertex_adjacency_list.list_length);
    vertex_adjacency_list.SetFirstElement();
    vertex_adjacency_list.AllocArrayElement();

    incident_edge_list.SetNumLists(vertex_poly_incidence.element.size());

    // Set arrays vertex_adjacency_list.element[] and 
    //   incident_edge_list.list_length[].
    for (NTYPE iv0 = 0; iv0 < NumVertices(); iv0++) {
      if (vertex_adjacency_list.ListLength(iv0) > 0) {
        NTYPE n = 0;
        NTYPE ifirst = adjacent.first_element[iv0];
        VTYPE2 iv1 = std::get<0>(adjacent.element[ifirst]);
        vertex_adjacency_list.ElementRef(iv0,n).SetVertex(iv1);
        NTYPE j = std::get<2>(adjacent.element[ifirst]);
        IncrementIncidentEdgeListLength(iv0, j);
        n++;
        for (NTYPE i = 1; i < adjacent.list_length[iv0]; i++) {
          VTYPE2 iv1 = std::get<0>(adjacent.element[ifirst+i]);
          NTYPE j1 = std::get<2>(adjacent.element[ifirst+i]);
          VTYPE2 iv2 = std::get<0>(adjacent.element[ifirst+i-1]);
          NTYPE j2 = std::get<2>(adjacent.element[ifirst+i-1]);
          if (iv1 != iv2) {
            vertex_adjacency_list.ElementRef(iv0,n).SetVertex(iv1);
            n++;
          }

          if (iv1 != iv2 || j1 != j2) 
            { IncrementIncidentEdgeListLength(iv0, j1); }
        }
        if (n != vertex_adjacency_list.ListLength(iv0)) {
          error.AddMessage
            ("Programming error.  Problem computing vertices adjacent to vertex ",
             iv0, ".");
          error.AddMessage
            ("  Expected ", vertex_adjacency_list.ListLength(iv0), 
             " vertices (no duplicates)");
          error.AddMessage
            ("    but computed ", n, " adjacent vertices.");
          throw error;
        }
      }
    }

    incident_edge_list.SetFirstElement();
    incident_edge_list.AllocArrayElement();

    current_incident_edge_index.resize(incident_edge_list.NumLists(), 0);

    // Set array incident_edge_list.element[].
    for (NTYPE iv0 = 0; iv0 < NumVertices(); iv0++) {
      if (vertex_adjacency_list.ListLength(iv0) > 0) {
        NTYPE n = 0;
        NTYPE ifirst = adjacent.first_element[iv0];
        NTYPE j = std::get<2>(adjacent.element[ifirst]);
        NTYPE k = IncidentEdgeListIndex(iv0,j);
        SetIncidentEdgeListElement(iv0, j, current_incident_edge_index[k], n);
        current_incident_edge_index[k]++;
        for (NTYPE i = 1; i < adjacent.list_length[iv0]; i++) {
          VTYPE2 iv1 = std::get<0>(adjacent.element[ifirst+i]);
          NTYPE j1 = std::get<2>(adjacent.element[ifirst+i]);
          VTYPE2 iv2 = std::get<0>(adjacent.element[ifirst+i-1]);
          NTYPE j2 = std::get<2>(adjacent.element[ifirst+i-1]);
          if (iv1 != iv2) { n++; }
          if (iv1 != iv2 || j1 != j2) {
            NTYPE k = IncidentEdgeListIndex(iv0,j1);
            SetIncidentEdgeListElement
              (iv0, j1, current_incident_edge_index[k], n);
            current_incident_edge_index[k]++;
          }
        }
        if (n+1 != vertex_adjacency_list.ListLength(iv0)) {
          error.AddMessage
            ("Programming error.  Problem setting incident edges for vertex ",
             iv0, ".");
          error.AddMessage
            ("  Expected ", vertex_adjacency_list.ListLength(iv0), 
             " vertices (no duplicates)");
          error.AddMessage
            ("    but computed ", n+1, " adjacent vertices.");
          throw error;
        }
      }
    }

    // Check stored correct number of edges for each incident edge list.
    for (SIZE_TYPE k = 0; k < incident_edge_list.NumLists(); k++) {
      if (current_incident_edge_index[k] != incident_edge_list.ListLength(k)) {
        error.AddMessage
          ("Programming error.  Incorrect storage of incident edges.");
        error.AddMessage
          ("  num_incident_edges[", k, "] = ", 
           incident_edge_list.ListLength(k), ".");
        error.AddMessage
          ("  Stored ", current_poly_index[k], " incident edges.");
        throw error;
      }
    }

  }

  // Return true if list of edges incident on vertex iv0 and 
  //  poly vertex_poly_incidence(iv0,j) contains edge (iv0, iv1).
  template <typename ETYPE0, typename ETYPE1, typename NTYPE>
  template <typename VTYPE0, typename VTYPE1, typename NTYPE2>
  bool VERTEX_POLY_EDGE_INCIDENCE<ETYPE0,ETYPE1,NTYPE>::
  DoesIncidenceListContainEdge
  (const VTYPE0 iv0, const VTYPE1 iv1, const NTYPE2 j) const
  {
    const NTYPE2 k = IncidentEdgeListIndex(iv0, j);
    for (NTYPE2 i = 0; i < incident_edge_list.ListLength(k); i++) {
      NTYPE2 n = incident_edge_list.Element(k, i);
      if (vertex_adjacency_list.Element(iv0, n).Vertex() == iv1) 
        { return(true); }
    }

    return(false);
  }


  // *****************************************************************
  // Class FACET_LIST_BASE member functions
  // *****************************************************************

  template <typename VTYPE, typename NTYPE, typename FACET_INFO_TYPE>
  template <typename VTYPE2, typename NTYPE2>
  void FACET_LIST_BASE<VTYPE,NTYPE,FACET_INFO_TYPE>::
  SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    const NTYPE NUM_VERT_PER_EDGE(2);
    NTYPE edge_endpoint[NUM_VERT_PER_EDGE];

    if (polymesh.NumPoly() == 0) { return; }
 
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      for (NTYPE j0 = 0; j0 < polymesh.NumPolyVert(ipoly); j0++) {
        const NTYPE j1 = (j0+1)%polymesh.NumPolyVert(ipoly);
        edge_endpoint[0] = polymesh.Vertex(ipoly, j0);
        edge_endpoint[1] = polymesh.Vertex(ipoly, j1);
        
        const NTYPE iedge = this->NumFacets();
        this->AddPolytope(edge_endpoint, NUM_VERT_PER_EDGE);
        this->poly_data[iedge].poly_containing_face = ipoly;
        this->poly_data[iedge].face_index = j0;
      }
    }
  }


  template <typename VTYPE, typename NTYPE, typename FACET_INFO_TYPE>
  template <typename VTYPE2, typename NTYPE2>
  void FACET_LIST_BASE<VTYPE,NTYPE,FACET_INFO_TYPE>::
  SetFrom2DMesh(const POLYMESH<VTYPE2,NTYPE2> & polymesh,
                const std::vector<bool> & flag_skip_poly)
  {
    const NTYPE NUM_VERT_PER_EDGE(2);
    NTYPE edge_endpoint[NUM_VERT_PER_EDGE];

    if (polymesh.NumPoly() == 0) { return; }
 
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {

      if (flag_skip_poly[ipoly]) { continue; }

      for (NTYPE j0 = 0; j0 < polymesh.NumPolyVert(ipoly); j0++) {
        const NTYPE j1 = (j0+1)%polymesh.NumPolyVert(ipoly);
        edge_endpoint[0] = polymesh.Vertex(ipoly, j0);
        edge_endpoint[1] = polymesh.Vertex(ipoly, j1);
        
        const NTYPE iedge = this->NumFacets();
        this->AddPolytope(edge_endpoint, NUM_VERT_PER_EDGE);
        this->poly_data[iedge].poly_containing_face = ipoly;
        this->poly_data[iedge].face_index = j0;
      }
    }
  }


  template <typename VTYPE, typename NTYPE, typename FACET_INFO_TYPE>
  template <typename VTYPE2, typename NTYPE2>
  void FACET_LIST_BASE<VTYPE,NTYPE,FACET_INFO_TYPE>::
  SetFromMeshOfSimplices(const POLYMESH<VTYPE2,NTYPE2> & polymesh)
  {
    if (polymesh.NumPoly() == 0) { return; }
 
    const NTYPE NUM_VERT_PER_SIMPLEX(polymesh.NumPolyVert(0));
    const NTYPE NUM_FACETS_PER_SIMPLEX(NUM_VERT_PER_SIMPLEX);
    const NTYPE NUM_VERT_PER_FACET(NUM_VERT_PER_SIMPLEX-1);
    std::vector<VTYPE> facet_vertices(NUM_VERT_PER_FACET);

    if (NUM_VERT_PER_SIMPLEX <= 1) { return; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      for (NTYPE jf = 0; jf < NUM_FACETS_PER_SIMPLEX; jf++) {
        for (NTYPE k = 0; k < NUM_VERT_PER_FACET; k++) {
          const VTYPE kv = (jf+1+k)%NUM_VERT_PER_SIMPLEX;
          facet_vertices[k] = polymesh.Vertex(ipoly, kv);
        }
        
        if (NUM_VERT_PER_FACET > 1) {
          if (jf%2 == 1) {
            // Correct orientation.
            std::swap(facet_vertices[NUM_VERT_PER_FACET-2],
                      facet_vertices[NUM_VERT_PER_FACET-1]);
          }
        }
        const NTYPE ifacet = this->NumFacets();
        this->AddPolytope(facet_vertices);
        this->poly_data[ifacet].poly_containing_face = ipoly;
        this->poly_data[ifacet].face_index = jf;
      }
    }
  }


  template <typename VTYPE, typename NTYPE, typename FACET_INFO_TYPE>
  template <typename VTYPE2, typename NTYPE2, typename CUBE_TYPE>
  void FACET_LIST_BASE<VTYPE,NTYPE,FACET_INFO_TYPE>::
  SetFromMeshOfCubes(const POLYMESH<VTYPE2,NTYPE2> & polymesh,
                     const CUBE_TYPE & cube)
  {
    std::vector<VTYPE> facet_vertices(cube.NumFacetVertices());

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      for (NTYPE jf = 0; jf < cube.NumFacets(); jf++) {
        for (NTYPE k = 0; k < cube.NumFacetVertices(); k++) {
          const VTYPE kv = cube.FacetVertex(jf, k);
          facet_vertices[k] = polymesh.Vertex(ipoly, kv);
        }
        const NTYPE ifacet = this->NumFacets();
        this->AddPolytope(facet_vertices);
        this->poly_data[ifacet].poly_containing_face = ipoly;
        this->poly_data[ifacet].face_index = jf;
      }
    }
  }


  // *****************************************************************
  // Functions on POLYMESH
  // *****************************************************************

  /// Return maximum poly vertex index.
  /// Return 0 if there are no polytope vertices.
  template <typename VTYPE, typename NTYPE>
  const VTYPE get_max_vert(const POLYMESH<VTYPE,NTYPE> & polymesh)
  {
    if (polymesh.element.size() > 0) {

      const VTYPE vmax = 
        *(std::max_element
          (polymesh.element.begin(), polymesh.element.end()));

      return(vmax);
    }
    else {
      return(0);
    }
  }

  /// Return sum of number of vertices in all polytopes.
  template <typename VTYPE, typename NTYPE>
  const NTYPE sum_num_poly_vert
  (const POLYMESH<VTYPE,NTYPE> & polymesh)
  {
    NTYPE sum = 0;
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) 
      { sum += polymesh.NumPolyVert(ipoly); }

    return(sum);
  }

  /// Count number of polytopes incident on each vertex.
  /// @param[out] num_incident_poly[] Array of number of incident polytopes.
  /// @pre num_incident_poly is preallocated to size = number of vertices.
  template <typename MESH_TYPE, typename NTYPE>
  void count_num_incident_poly
  (const MESH_TYPE & polymesh, std::vector<NTYPE> & num_incident_poly)
  {
    const NTYPE numv = num_incident_poly.size();
    std::vector<NTYPE> last_poly(numv);

    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < numv; iv++) 
      { num_incident_poly[iv] = 0; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      for (NTYPE i = 0; i < polymesh.NumPolyVert(ipoly); i++) {
        VTYPE iv = polymesh.Vertex(ipoly, i);

        if (num_incident_poly[iv] > 0) {
          if (last_poly[iv] == ipoly) {
            // Poly ipoly already counted as incident on iv.
            continue;
          }
        }

        num_incident_poly[iv]++;
        last_poly[iv] = ipoly;
      }
    }
  }


  /// Count number of edges incident on each vertex.
  /// If same edge is in k polygons, counts edge k times.
  /// @param polymesh Mesh of 2D polygons.  Vertices are listed in cyclic
  ///   order around each polygon.
  /// @param[out] num_incident_edges[] Array of number of incident edges.
  /// @pre num_incident_edges is preallocated to size = number of vertices.
  template <typename MESH_TYPE, typename NTYPE>
  void count_num_poly2D_edges_incident_on_vertices
  (const MESH_TYPE & polymesh, std::vector<NTYPE> & num_incident_edges)
  {
    const NTYPE numv = num_incident_edges.size();

    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (VTYPE iv = 0; iv < numv; iv++) 
      { num_incident_edges[iv] = 0; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);
      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          VTYPE iv0 = polymesh.Vertex(ipoly, i0);
          NTYPE i1 = (i0+1)%num_poly_vert;
          VTYPE iv1 = polymesh.Vertex(ipoly, i1);
          if (iv0 != iv1) {
            num_incident_edges[iv0]++;
            num_incident_edges[iv1]++;
          }
        }
      }
    }
  }

  /// Count number of (hyper)cube edges incident on each vertex.
  /// If same edge is in k polygons, counts edge k times.
  /// @param polymesh Mesh of (hyper) cubes.  
  ///   - Vertices are listed in order: 
  ///  (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1).
  /// @param[out] num_incident_edges[] Array of number of incident edges.
  /// @pre num_incident_edges is preallocated to size = number of vertices.
  template <typename MESH_TYPE, typename CUBE_TYPE, typename NTYPE>
  void count_num_cube_edges_incident_on_vertices
  (const MESH_TYPE & polymesh, const CUBE_TYPE & cube,
   std::vector<NTYPE> & num_incident_edges)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const DTYPE dimension = cube.Dimension();
    const NTYPE numv = num_incident_edges.size();

    for (VTYPE iv = 0; iv < numv; iv++) 
      { num_incident_edges[iv] = 0; }

    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      const NTYPE num_poly_vert = polymesh.NumPolyVert(ipoly);
      if (num_poly_vert > 1) {
        for (NTYPE i0 = 0; i0 < num_poly_vert; i0++) {
          const VTYPE iv0 = polymesh.Vertex(ipoly, i0);
          for (NTYPE d = 0; d < dimension; d++) {
            const NTYPE i1 = cube.VertexNeighbor(i0,d);
            const VTYPE iv1 = polymesh.Vertex(ipoly, i1);

            if (iv0 < iv1) {
              num_incident_edges[iv0]++;
              num_incident_edges[iv1]++;
            }
          }
        }
      }
    }
  }


  /// Compute vertex link in hexahedral mesh.
  /// @param cube Cube information.
  /// @pre cube has same dimension as cubes in mesh.
  /// @param[out] link_mesh Polygonal mesh of links.
  ///   Note that quadrilaterals are stored lower-left, lower-right,
  ///   upper-left, upper_right, NOT in clockwise or counter-clockwise order.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename VTYPE, typename CUBE_TYPE, typename POLYMESH2_TYPE>
  void compute_vertex_link_in_hex_mesh
  (const POLYMESH_TYPE & polymesh, 
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
   const VTYPE iv,
   const CUBE_TYPE & cube,
   POLYMESH2_TYPE & link_mesh)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE num_vertices_per_facet = cube.NumFacetVertices();
    std::vector<VTYPE> facet_vertex(num_vertices_per_facet);

    link_mesh.Clear();

    for (NTYPE j = 0; j < vertex_info.NumIncidentPoly(iv); j++) {
      const NTYPE jpoly = vertex_info.IncidentPoly(iv,j);

      for (NTYPE jfacet = 0; jfacet < cube.NumFacets(); jfacet++) {

        for (NTYPE k = 0; k < num_vertices_per_facet; k++) {
          const NTYPE k2 = cube.FacetVertex(jfacet, k);
          const VTYPE kv = polymesh.Vertex(jpoly, k2);
          facet_vertex[k] = kv;
        }

        if (!IJK::does_list_contain(facet_vertex, iv)) {
          link_mesh.AddPolytope(facet_vertex);
        }
      }
    }
  }


  /// Compute vertex link in tetrahedral mesh.
  /// @param[out] link_mesh Polygonal mesh of links.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename VTYPE, typename POLYMESH2_TYPE>
  void compute_vertex_link_in_tet_mesh
  (const POLYMESH_TYPE & polymesh, 
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
   const VTYPE iv,
   POLYMESH2_TYPE & link_mesh)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE num_vertices_per_tet(4);
    const NTYPE num_vertices_per_facet(num_vertices_per_tet-1);
    std::vector<VTYPE> facet_vertex(num_vertices_per_facet);

    link_mesh.Clear();

    for (NTYPE j = 0; j < vertex_info.NumIncidentPoly(iv); j++) {
      const NTYPE jpoly = vertex_info.IncidentPoly(iv,j);

      NTYPE k = 0;
      for (NTYPE k2 = 0; k2 < num_vertices_per_tet; k2++) {

        const VTYPE kv = polymesh.Vertex(jpoly, k2);
        if (kv != iv) {
          facet_vertex[k] = kv;
          k++;
        }
      }

      link_mesh.AddPolytope(facet_vertex);
    }
  }

}

#endif
