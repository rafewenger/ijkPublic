/*!
 *  @file ijkMCtable.h
 *  @brief Class containing a table of Marching Cubes isosurface patches 
 *    in a given polyhedron.
 *  - All 2^numv +/- patterns are stored in the table 
 *    where numv = # polyhedron vertices.
 *  - Version 0.5.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2001-2024 Rephael Wenger

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

#ifndef _IJKTABLE_
#define _IJKTABLE_

#include "ijk.tpp"
#include "ijkenum.tpp"

#include "ijkMCtable_poly.h"
#include "ijkMCtable_properties.h"

/// Classes and routines for storing and manipulating isosurface lookup table.
namespace IJKMCUBE_TABLE {

  typedef int TABLE_INDEX;    ///< Index of entry in isosurface lookup table.

  // ***************************************************************
  // ISOSURFACE VERTEX
  // ***************************************************************

  /// Isosurface vertex class.
  class ISOSURFACE_VERTEX {

  public:
    typedef enum {VERTEX, EDGE, FACET, POINT} ISOSURFACE_VERTEX_TYPE;
    typedef float COORD_TYPE;

  protected:
    ISOSURFACE_VERTEX_TYPE vtype;
    int face;
    int num_coord;
    COORD_TYPE * coord;
    std::string label;
    bool is_label_set;

  public:
    ISOSURFACE_VERTEX();        // constructor
    ~ISOSURFACE_VERTEX();       // destructor

    // Get Functions

    /// Return isosurface vertex type.
    ISOSURFACE_VERTEX_TYPE Type() const { return(vtype); };

    /// Return index of face (vertex, edge, facet) containing isosurface.
    /// Valid only if vertex type is VERTEX, EDGE or FACET.
    int Face() const { return(face); };

    /// Return coordinate of isosurface vertex.
    /// Valid only if vertex type is POINT.
    COORD_TYPE Coord(const int d) const { return(coord[d]); };

    /// Return number of isosurface vertex coordinates.
    /// Valid only if vertex type is POINT.
    int NumCoord() const { return(num_coord); };

    /// Return label of isosurface vertex.
    /// Used for extending vertex types.
    std::string Label() const { return(label); };

    /// Return true if label is set.
    bool IsLabelSet() const { return(is_label_set); };

    // Set Functions
    void SetType(const ISOSURFACE_VERTEX_TYPE t) { vtype = t; };
    void SetFace(const int index) { face = index; };
    void SetNumCoord(const int numc);
    void SetCoord(const int ic, const COORD_TYPE c) 
    { coord[ic] = c; };
    void SetLabel(const std::string & s) 
    { label = s; is_label_set = true; };
  };


  // ***************************************************************
  // ISOSURFACE TABLE
  // ***************************************************************

  /*!
   *  @brief Isosurface lookup table.
   *  - Stores isosurface patches for each configuration
   *    of +/- labels at polyhedron vertices.
   */
  class ISOSURFACE_TABLE {

  protected:

    /// Entry in the isosurface lookup table.
    class ISOSURFACE_TABLE_ENTRY {

    public:
      int num_simplices;
      ISOSURFACE_VERTEX_INDEX * simplex_vertex_list;
      ISOSURFACE_TABLE_ENTRY();                  // constructor
      ~ISOSURFACE_TABLE_ENTRY();                 // destructor

      bool Check(IJK::ERROR & error_msg) const;
      void FreeAll();                            // free all memory
    };


  public:

    /// @brief Index of entry in isosurface lookup table.
    /// - Define within ISOSURFACE_TABLE for use in templates.
    typedef IJKMCUBE_TABLE::TABLE_INDEX TABLE_INDEX;    


  protected:

    /// @brief Isosurface table properties.
    ISOSURFACE_TABLE_PROPERTIES table_properties;

    ISOTABLE_POLY polytope;      ///< Mesh polytope.
    int simplex_dimension;       ///< Simplex dimension.

    /// Array of isosurface vertex descriptors.
    ISOSURFACE_VERTEX * isosurface_vertex;

    /// Number of vertices in array isosurface_vertex[].
    int num_isosurface_vertices;

    /// Array of isosurface table entries.
    ISOSURFACE_TABLE_ENTRY * entry;

    /// Number of entries in table.
    long num_table_entries;         

    /// Maximum number of vertices allowed for table polyhedron.
    int max_num_vertices; 

    /// True, if array num_table_entries[] is allocated.
    bool is_table_allocated;  

    /// @brief Check if isosurface vertices are allocated.
    /// - Throw error if not enough isosurface vertices are allocated.
    void CheckIsoVerticesAlloc
    (const char * procname, const int vstart, const int numv);

    /// Initialization routine.
    void Init(const int dimension, const int simplex_dimension);


  public:

    /// @name Constructors and destructors.
    ///@{
    ISOSURFACE_TABLE();
    ISOSURFACE_TABLE(const int d);
    ISOSURFACE_TABLE(const int dimension, const int simplex_dimension);

    ~ISOSURFACE_TABLE();                ///< Destructor

    ///@}

    
    /// @name Get Functions
    ///@{

    /// @brief Return table encoding.
    ENCODING Encoding() const 
    { return table_properties.encoding; };

    /// @brief Return string for table encoding.
    std::string EncodingName() const 
      { return table_properties.EncodingString(Encoding()); };

    /// @overload
    /// @brief Return string for encoding _encoding.
    std::string EncodingName(const ENCODING _encoding) const
      { return table_properties.EncodingString(_encoding); };

    /// Return isosurface table properties.
    const ISOSURFACE_TABLE_PROPERTIES & Properties() const
    { return table_properties; }

    /// Return polytope dimension.
    int Dimension() const { return polytope.Dimension() ; };

    /// Return isosurface simplex dimension.
    int SimplexDimension() const { return simplex_dimension; };

    /// Return true if lookup table represents an interval volume.
    bool IsIntervalVolume() const
    { return (Dimension() == SimplexDimension()); }
    
    /// Return number of vertices in each isosurface simplex.
    int NumVerticesPerSimplex() const { return (SimplexDimension()+1); };

    /// @brief Return number of isosurface vertices
    ///   in array isosurface_vertex[].
    int NumIsosurfaceVertices() const
    { return num_isosurface_vertices; };

    /// Return number of lookup table entries.
    int NumTableEntries() const { return num_table_entries; };

    /// Access isosurface table polytope.
    const ISOTABLE_POLY_BASE & Polytope() const
    { return polytope; };

    /// Return polytope shape.
    POLYTOPE_SHAPE PolyShape() const
    { return Polytope().Shape(); }

    /// Access i'th isosurface vertex.
    const ISOSURFACE_VERTEX & IsosurfaceVertex(const int i) const
    { return isosurface_vertex[i]; }; 

    /// @brief Return number of simplices in isosurface patch
    ///   for table entry \a it.
    int NumSimplices(const TABLE_INDEX it) const
    { return entry[it].num_simplices; }; 

    /*!
     *  @brief Return \a k'th vertex of isosurface simplex \a is, 
     *    table entry \a it.
     *  @param it Index of table entry.
     *  @param is Simplex \a is of table entry \a it.
     *  @param k Return \a k'th vertex of simplex \a is.
     */
    ISOSURFACE_VERTEX_INDEX SimplexVertex
    (const TABLE_INDEX it, const int is, const int k) const
    { return(entry[it].simplex_vertex_list[is*NumVerticesPerSimplex()+k]); };

    /*!
     *  @brief Return pointer to vertices of all simplices
     *    in table entry table_index.
     */
    const ISOSURFACE_VERTEX_INDEX *
    SimplexVertices(const TABLE_INDEX table_index) const
    { return entry[table_index].simplex_vertex_list; }
      
    /*!
     *  @brief Return pointer to vertices of simplex isimplex
     *    in table entry table_index.
     */
    const ISOSURFACE_VERTEX_INDEX *
    SimplexVertices(const TABLE_INDEX table_index,
                    const int isimplex) const
    { return (entry[table_index].simplex_vertex_list +
              NumVerticesPerSimplex()*isimplex); }

    /*!
     *  @brief Return maximum number of polytope vertices permitted 
     *    in any table.
     *  - Note: Even tables for polyhedra of this size are probably impossible
     *    to compute/store.
     */
    int MaxNumVertices() const { return max_num_vertices; };

    /// Return true if table memory is allocated.
    bool IsTableAllocated() const
    { return(is_table_allocated); };

    /// @brief Return opposite separation type.
    ISOSURFACE_SEPARATION_TYPE OppositeSeparationType() const
    { return Properties().OppositeSeparationType(); }

    /// @brief Return orientation of isosurface polytopes.
    ISO_POLY_ORIENTATION IsoPolyOrientation() const
    { return Properties().IsoPolyOrientation(); }

    /// @brief Return opposite orientation of isosurface polytopes.
    ISO_POLY_ORIENTATION OppositeIsoPolyOrientation() const
    { return Properties().OppositeIsoPolyOrientation(); }

    /// Return base of isotable encoding.
    int Base() const
    {
      if (Encoding() == BASE3) { return 3; }
      else { return 2; }
    }
      
    /*!
     *  @brief Value representing negative "-" label in isosurface table
     *    with BINARY encoding.
     *  - IJK::convert2base converts table index to array of digits,
     *    each one representing label of a vertex.
     *  - BinaryNegative() is the digit indicating a negative label.
     */
    constexpr int BinaryNegative() const { return 0; }

    /// @brief Value representing positive "+" label in isosurface table
    ///   with BINARY encoding.
    constexpr int BinaryPositive() const { return 1; }

    /*!
     *  @brief Value representing negative "-" label in isosurface table
     *    with BASE3 encoding.
     *  - IJK::convert2base converts table index to array of digits,
     *    each one representing label of a vertex.
     *  - Base3Negative() is the digit indicating a negative label.
     */
    constexpr int Base3Negative() const { return 0; }

    /// @brief Value representing positive "+" label in isosurface table
    ///   with BASE3 encoding.
    constexpr int Base3Positive() const { return 2; }
    
    /// @brief Value representing positive "=" label in isosurface table
    ///   with BASE3 encoding.
    /// - For interval volumes, "=" represents vertex in volume interior.
    constexpr int Base3Equals() const { return 1; }

    /*!
     *  @brief Return value representing negative label.
     *  - IJK::convert2base converts table index to array of digits,
     *    each one representing label of a vertex.
     *  - NegativeLabelValue() returns digit indicate a negative label.
     */
    int NegativeLabelValue() const
    {
      if (Encoding() == BASE3) { return Base3Negative(); }
      else { return BinaryNegative(); }
    }

    /*!
     *  @brief Return value representing positive label.
     *  - IJK::convert2base converts table index to array of digits,
     *    each one representing label of a vertex.
     *  - PositiveLabelValue() returns digit indicating a positive label.
     */
    int PositiveLabelValue() const
    {
      if (Encoding() == BASE3) { return Base3Positive(); }
      else { return BinaryPositive(); }
    }

    /*!
     *  @brief Return true if facet vertex labels are identical
     *    in table entries table_indexA and table_indexB.
     */
    bool AreAllFacetVertexLabelsIdentical
    (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     const int ifacet) const;

    ///@}

    
    /// @name Set Polytope Functions
    ///@{
    
    void SetDimension(const int d) { polytope.SetDimension(d); };
    void SetPolyShape(const POLYTOPE_SHAPE & shape)
    { polytope.SetShape(shape); }
    void SetNumPolyVertices(const int numv) 
    { polytope.SetNumVertices(numv); };
    void SetNumPolyEdges(const int nume) { polytope.SetNumEdges(nume); };
    void SetNumPolyFacets(const int numf) { polytope.SetNumFacets(numf); };
    void SetPolySize(const int numv, const int nume, const int numf)
    { SetNumPolyVertices(numv); SetNumPolyEdges(nume); 
      SetNumPolyFacets(numf); };
    void SetPolyVertexCoord(const int iv, const int ic, const int coord)
    { polytope.SetVertexCoord(iv, ic, coord); };
    // Note: SetNumPolyVertices or SetPolySize must be called before 
    //   SetPolyVertexCoord
    void SetPolyEdge(const int ie, const int iv0, const int iv1)
    { polytope.SetEdge(ie, iv0, iv1); };
    // Note: SetNumPolyEdges or SetPolySize must be called before SetPolyEdge
    void SetPolyNumFacetVertices(const int jf, const int numv)
    { polytope.SetNumFacetVertices(jf, numv); }
    void SetPolyFacetVertex(const int jf, const int k, const int iv)
    { polytope.SetFacetVertex(jf, k, iv); };
    // Note: SetPolyNumFacetVertices must be called before SetPolyFacetVertex

    /// @brief Set polytope to _poly.
    void Set(const ISOTABLE_POLY_BASE & _poly)
    { this->polytope = _poly; };

    ///@}

    /// @name Set Isosurface Vertices Functions
    ///@{
    void SetNumIsosurfaceVertices(const int num_vertices);
    void SetIsoVertexType(const int i, 
                          const ISOSURFACE_VERTEX::ISOSURFACE_VERTEX_TYPE t) 
    { isosurface_vertex[i].SetType(t); };
    void SetIsoVertexFace(const int i, const int index) 
    { isosurface_vertex[i].SetFace(index); };
    void SetIsoVertexNumCoord(const int i, const int numc)
    { isosurface_vertex[i].SetNumCoord(numc); };
    void SetIsoVertexCoord(const int i, 
                           const int ic, const ISOSURFACE_VERTEX::COORD_TYPE c) 
    { isosurface_vertex[i].SetCoord(ic, c); };
    void SetIsoVertexLabel(const int i, const std::string & s) 
    { isosurface_vertex[i].SetLabel(s); };

    /// @brief Set isosurface vertex.
    /// @pre SetNumIsosurfaceVertices() has already been called.
    void SetIsosurfaceVertex
    (const int iv, const ISOSURFACE_VERTEX & isosurface_vertex);
    
    /// @brief Copy isosurface vertices from isotable.
    void CopyIsosurfaceVertices(const ISOSURFACE_TABLE & isotable);
                                
    // store polytope vertices, edges or faces as isosurface vertices
    void StorePolyVerticesAsIsoVertices(const int vstart);
    void StorePolyEdgesAsIsoVertices(const int vstart);
    void StorePolyFacetsAsIsoVertices(const int vstart);
    ///@}

    /// @name Set Isosurface Table Functions
    ///@{
    void SetSimplexDimension(const int d) { this->simplex_dimension = d; };
    void SetEncoding(const ENCODING encoding);
    void SetBinaryEncoding() { SetEncoding(BINARY); };
    void SetBase3Encoding() { SetEncoding(BASE3); };
    void SetEncoding(const std::string & encoding_str) 
    { table_properties.SetEncoding(encoding_str); }

    virtual void SetNumTableEntries(const int num_table_entries);
    void SetNumSimplices(const TABLE_INDEX it, const int nums);

    /// @brief Set vertex iv of simplex is in table entry it to isov.
    void SetSimplexVertex(const TABLE_INDEX it, const int is, 
                          const int iv, const ISOSURFACE_VERTEX_INDEX isov);

    /// @brief Set simplex vertices.
    void SetSimplexVertices
    (const TABLE_INDEX it, const ISOSURFACE_VERTEX_INDEX simplex_vertices[],
     const int num_simplices);

    void SetTableType(const LOOKUP_TABLE_TYPE lookup_table_type);
    void SetGridVertexLabelType
      (const GRID_VERTEX_LABEL_TYPE grid_vertex_label_type);
    void SetGridVertexLabelType(const std::string s)
    { table_properties.SetGridVertexLabelType(s); }
    void SetSeparationType
      (const ISOSURFACE_SEPARATION_TYPE separation_type);
    void SetSeparationType(const std::string s)
    { table_properties.SetSeparationType(s); }
    void SetTriangulationType
      (const ISOSURFACE_TRIANGULATION_TYPE triangulation_type);
    void SetTriangulationType(const std::string s)
    { table_properties.SetTriangulationType(s); }
    void SetIsoPolyOrientation
      (const ISO_POLY_ORIENTATION iso_poly_orientation);
    void SetIsoPolyOrientation(const std::string s)
    { table_properties.SetIsoPolyOrientation(s); }
    void SetSeparateOpposite
      (const SEPARATE_OPPOSITE_TYPE separate_opposite);
    void SetSeparateOpposite(const std::string s)
    { table_properties.SetSeparateOpposite(s); }
    void SetSeparateOpposite(const bool flag);

    ///@}

    /// @name Copy Isosurface Table Functions
    ///@{

    /// @brief Copy polytope from isotable.
    void CopyPolytope(const ISOSURFACE_TABLE & isotable)
    { Set(isotable.Polytope()); }

    /// @brief Copy isosurface table properties from isotable.
    void CopyProperties(const ISOSURFACE_TABLE & isotable)
    { table_properties.Copy(isotable.Properties()); }

    ///@}

    /// @name Generate Polytope Functions
    ///@{
    void GenCube(const int cube_dimension) 
    { polytope.GenCube(cube_dimension); };
    void GenCubeOrderA(const int cube_dimension) 
    { polytope.GenCubeOrderA(cube_dimension); };
    // Note: Cubes of dimension > 4 will have too many vertices
    void GenSimplex(const int simplex_dimension) 
    { polytope.GenSimplex(simplex_dimension); };
    void GenPyramid(const int pyramid_dimension) 
    { polytope.GenPyramid(pyramid_dimension); };
    ///@}

    /// @name Process Simplex Orientations
    ///@{

    /*!
     *  @brief Sort vertices of simplex isimplex in increasing order.
     *  - Used to orient initial simplex for setting table orientation.
     */
    void SortSimplexVertices
    (const TABLE_INDEX it, const int isimplex);
    
    /*!
     *  @brief Flip isosurface polytope orientations from +1 to -1 
     *    or from -1 to +1.
     *  - Changes simplx orientation by swapping last two vertices
     *    in simplex.
     *  @param it Table entry index.  In range [0..NumTableEntries()-1].
     *  @param ipoly Polytope index.
     */
    void FlipIsoPolyOrientation(const TABLE_INDEX it, const int ipoly);

    /*!
     *  @brief Flip all isosurface polytope orientations 
     *    at table entry table_index.
     *  - Changes simplex orientation by swapping last two vertices
     *    in simplex.
     *  @param it Table entry index.  In range [0..NumTableEntries()-1].
     */
    void FlipAllIsoPolyOrientations(const TABLE_INDEX table_index);

    /*!
     *  @brief Flip all isosurface polytope orientations 
     *    in isosurface lookup table.
     *  - Changes simplex orientation by swapping last two vertices
     *    in simplex.
     *  - Reverses iso_table_orientation property, if that property
     *    is set.
     */
    void FlipAllIsoPolyOrientations();

    /*!
     *  @brief Orient simplices in table entry.
     *  - Note: Only simplices in same connected component as
     *    simplex istart are oriented.
     *  @param table_index Index of table entry.
     *  @param istart Index of starting simplex.
     *    @pre istart is in range [0..(NumSimplices(it)-1].
     */
    void OrientSimplicesInTableEntry
    (const TABLE_INDEX table_index, const TABLE_INDEX istart);

    /*!
     *  @overload
     *  @brief Orient simplices in table entry. (Argument is_oriented[].)
     *  - Version that uses and sets array is_orientedA[].
     *  @param[out] is_orientedA[isimplexA] 
     *    True if simplex isimplexA is oriented.
     *    - Note: OrientIsoTableWithSimplexFacet() ignores any simplices
     *      where is_orientedA[isimplexA] is true.
     */
    void OrientSimplicesInTableEntry
    (const TABLE_INDEX table_index, const TABLE_INDEX istart,
     std::vector<bool> & is_oriented);

    /*!
     *  @brief Orient all simplices in table entry.
     *  - Orientation of each connected component is arbitrary.
     *  @param[out] num_components Number of connected components
     *    in entry[table_index].simplex_vertex_list[].
     */
    void OrientAllSimplicesInTableEntry
    (const TABLE_INDEX table_index, int & num_components);
    
    /*!
     *  @brief Return false if simplices in table entry have 
     *     inconsistent orientations.
     *  @param table_index Index of table entry.
     */
    bool AreSimplicesConsistentlyOriented
    (const TABLE_INDEX table_indexA,
     int & isimplexA, int & isimplexB) const;

    ///@}

    
    /// @name Memory Management Functions
    
    ///@{
    virtual void FreeAll();                     /// Free all memory.
    ///@}


    /// @name Check Functions
    
    ///@{
    bool CheckDimension(const int d) const;
    bool CheckDimension() const
    { return(CheckDimension(Dimension())); };
    bool CheckTable(IJK::ERROR & error_msg) const;
    bool Check(IJK::ERROR & error_msg) const;

    /// @brief Return false and set error message if table_index
    ///   is not in range [0..NumTableEntries()-1].
    bool CheckTableIndex
    (const int table_index, IJK::ERROR & error) const;
    
    ///@}

  };

  typedef ISOSURFACE_TABLE * ISOSURFACE_TABLE_PTR;


  // ***************************************************************
  // ISOSURFACE EDGE TABLE
  // ***************************************************************

  /// @brief Isosurface edge table.
  /// - Store list of edges containing isosurface vertices.
  class ISOSURFACE_EDGE_TABLE:public ISOSURFACE_TABLE {

  protected:

    /// Entry in isosurface edge table.
    class ISOSURFACE_EDGE_TABLE_ENTRY {

    public:
      int num_edges;
      EDGE_INDEX * edge_endpoint_list;

      ISOSURFACE_EDGE_TABLE_ENTRY();                  // constructor
      ~ISOSURFACE_EDGE_TABLE_ENTRY();                 // destructor

      bool Check(IJK::ERROR & error_msg) const;
      void FreeAll();                            // free all memory
    };

  protected:
    ISOSURFACE_EDGE_TABLE_ENTRY * edge_entry;

    // initialization routine
    void Init(const int d);

  public:
    ISOSURFACE_EDGE_TABLE(const int d);      // constructor
    ~ISOSURFACE_EDGE_TABLE();                // destructor

    // get functions

    int NumEdges(const TABLE_INDEX it) const 
    { return(edge_entry[it].num_edges); };

    // it = table entry index. ie = edge index. iend = endpoint index (0 or 1)
    EDGE_INDEX EdgeEndpoint
    (const TABLE_INDEX it, const int ie, const int iend) const
    { return(edge_entry[it].edge_endpoint_list[2*ie+iend]); };

    // set isosurface table functions
    virtual void SetNumTableEntries(const int num_table_entries);

    // generate edge lists
    void GenEdgeLists();

    // check functions
    bool CheckTable(IJK::ERROR & error_msg) const;
    bool Check(IJK::ERROR & error_msg) const;

    // free memory
    virtual void FreeAll();                     // free all memory
  };

  typedef ISOSURFACE_EDGE_TABLE * ISOSURFACE_EDGE_TABLE_PTR;


  // ***************************************************************
  // TABLE MANIPULATION ROUTINES
  // ***************************************************************

  /*!
   *  @brief Invert marching cubes lookup table.
   *  - Swap isosurface patches in table[i] and table[numEntry-1-i]
   *    where numEntry is the number of isosurface table entries.
   *  - Note: Swapping patches swaps SEPARATE_NEG and SEPARATE_POS
   *    and flips isosurface polytope orientations.
   *  @param isotableA Input isosurface table.
   *  @param[out] Inverted table.
   *    - isotableB[i] has isosurface patch isotableA[numEntry-1-i]
   *      where numEntry is the number of isosurface table entries.
   *    - isotableB.SeparationType() is the opposite of
   *      isotableA.SeparationType().
   *    - isotableB.IsoPolyOrientation() is the opposite of
   *      isotableA.IsoPolyOrientation().
   */
  void invert_mcube_isotable
    (const ISOSURFACE_TABLE & isotableA, ISOSURFACE_TABLE & isotableB);

  // ***************************************************************
  // UTILITY FUNCTIONS
  // ***************************************************************

  // calculate number of entries required in ISOSURFACE_TABLE
  unsigned long calculate_num_entries(const int num_vert, const int num_colors);

}

#endif
