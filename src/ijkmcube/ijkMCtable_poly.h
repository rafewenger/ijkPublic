/*!
 *  @file ijkMCtable_poly.h
 *  @brief Class containing isosurface table polyhedron.
 *  - Version 0.5.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2024 Rephael Wenger

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


#ifndef _IJKTABLE_POLY_
#define _IJKTABLE_POLY_

#include <cctype>
#include <string>
#include <vector>

#include "ijk.tpp"
#include "ijkenum.tpp"

/// Classes and routines for storing and manipulating isosurface lookup table.
namespace IJKMCUBE_TABLE {

  typedef unsigned char 
  ISOSURFACE_VERTEX_INDEX;  ///< Index of isosurface vertex.
  typedef unsigned char EDGE_INDEX;    ///< Index of edge.
  typedef unsigned char FACET_INDEX;   ///< Index of facet.
  typedef int FACET;          ///< Bits representing vertices in facet.
  typedef int FACET_SET;      ///< Bits representing set of facets.

  const int NO_VERTEX = -1;

  // *****************************************************************
  // Enum type POLYTOPE_SHAPE
  // *****************************************************************

  /// @brief Isosurface table polytope shape.
  typedef enum { CUBE, SIMPLEX, PYRAMID, SIMPLEX_PRISM, UNDEFINED_SHAPE }
  POLYTOPE_SHAPE;

  namespace {
    static std::vector<IJK::ENUM_STR<POLYTOPE_SHAPE> > _poly_shape_list =
      { {CUBE,"Cube"}, {SIMPLEX,"Simplex"}, {PYRAMID,"Pyramid"},
	{SIMPLEX_PRISM,"SimplexPrism"}, {UNDEFINED_SHAPE,"UndefinedShape"} };
  };
  
  // *****************************************************************
  // MARCHING CUBES ISOSURFACE TABLE POLYTOPE BASE
  // *****************************************************************

  /// @brief Base class for Marching Cubes isosurface table polytope.
  /// - All set/generate functions are protected.
  class ISOTABLE_POLY_BASE {

  protected:

    /// @brief Strings corresponding to enum POLYTOPE_SHAPE.
    /// - Set shape_list from _poly_shape_list[] in constructor.
    IJK::ENUM_LIST<POLYTOPE_SHAPE> shape_list;

  protected:
    POLYTOPE_SHAPE shape;        ///< Polytope shape.
    int dimension;               ///< Polytope dimension.
    int num_vertices;            ///< Number of polytope vertices.
    int num_edges;               ///< Number of polytope edges.
    int num_facets;              ///< Number of polytope facets.
    int * vertex_coord;          ///< Polytope vertex coordinates.
    int * edge_endpoint;         ///< Polytope edge endpoints.
    int * num_facet_vertices;    ///< Number of vertices of each facet.
    int ** facet_vertex_list;    ///< List of vertices in each facet.

    // *** OBSOLETE? NEVER USED? ***
    int * num_vertex_neighbors;  ///< Number of neighbors of each vertex.
    int ** vertex_neighbor_list; ///< List of neighors of each vertex.

    int * num_incident_edges;    ///< Number of edges incident on vertex.
    int ** incident_edge_list;   ///< Lists of edges incident on each vertex.
    
    FACET * facet;               ///< Polytope facets.
    void Init();                 ///< Initialize.
    void FreeFacets();           ///< Free all facet arrays.
    void FreeIncidentEdges();    ///< Free all incident edge arrays.
    
    // *** OBSOLETE? NEVER USED? ***
    void FreeVertexNeighbors();  ///< Free all vertex neighbors arrays.

  protected:
    
    /// @name Set/Compute Functions
    //@{

    /// Set polytope shape.
    void SetShape(const POLYTOPE_SHAPE & shape);

    /// Set polytope dimension.
    void SetDimension(const int d);

    /// @brief Set number of polytope vertices.
    /// - Must be called before setting polytope vertices.
    void SetNumVertices(const int numv); 

    /// @brief Set number of polytope edges.
    /// - Must be called before setting polytope edges.
    void SetNumEdges(const int nume);

    /// @brief Set number of polytope facets.
    /// - Must be called before setting polytope facets.
    void SetNumFacets(const int numf);

    /// Set number of polytope vertices, edges and facets.
    void SetSize(const int numv, const int nume, const int numf)
    { SetNumVertices(numv); SetNumEdges(nume); SetNumFacets(numf); };

    /// @brief Set \a ic'th coordinate of vertex \a iv.
    /// @pre SetNumVertices() or SetSize() must be called
    ///   before SetVertexCoord().
    void SetVertexCoord 
    (const int iv, const int ic, const int coord);

    /*!
     *  @brief Set endpoints of edge \a ie.
     *  @pre SetNumEdges or SetSize must be called before SetEdge.
     *  @param ie Edge index.  In range [0..NumEdges()-1].
     *  @param iv0 Endpoint 0 index.  In range [0..NumVertices()-1].
     *  @param iv1 Endpoint 1 index.  In range [0..NumVertices()-1].
     */
    void SetEdge(const EDGE_INDEX ie, const int iv0, const int iv1);

    /// @brief Set number of vertices in facet \a jf.
    /// @pre SetNumFacets() or SetSize() must be called
    ///   before SetNumFacetVertices().
    void SetNumFacetVertices 
    (const FACET_INDEX jf, const int numv);

    /// @brief Set \a k'th facet vertex of facet \a jf to vertex \a iv.
    /// @pre SetNumFacetVertices(jf) must be called before SetFacetVertex().
    void SetFacetVertex(const FACET_INDEX jf, const int k, const int iv);

    /// @brief Compute edges incident on each vertex.
    /// - Note: All edges should be set before calling this routine.
    void ComputeIncidentEdges();
    
    // *** OBSOLETE? NEVER USED? ***
    /// @brief Compute vertex neighbors.
    /// - Note: All edges should be set before calling this routine.
    void ComputeVertexNeighbors();
    
    //@}

    
    /// @name Memory Management Functions
    //@{
    void FreeAll();             ///< Free all memory.
    //@}

    
    /// @name Generate polytope.
    //@{

    /// Generate a square, cube or hypercube.
    void GenCube(const int cube_dimension); 

    /// @brief Generate a square, cube or hypercube.
    /// - Use order of facets and edges given by template class CUBE_FACE_INFO.
    void GenCubeOrderA(const int cube_dimension);

    /// Generate a triangle, tetrahedron or simplex.      
    void GenSimplex(const int simplex_dimension);

    /// Generate a pyramid over a square, cube or hypercube base.
    void GenPyramid(const int pyramid_dimension);

    //@}

    
  public:
    ISOTABLE_POLY_BASE(const int d);           ///< Constructor
    ~ISOTABLE_POLY_BASE();                   ///< Destructor
    ISOTABLE_POLY_BASE                       ///< Copy constructor.
      (const ISOTABLE_POLY_BASE & init);
    
    const ISOTABLE_POLY_BASE & operator =    ///< Assignment. 
      (const ISOTABLE_POLY_BASE & right);

    /// @name Get Functions
    //@{
    int Dimension() const   /// Polytope dimension.
    { return dimension; };
    int NumVertices() const /// Number of polytope vertices.
    { return num_vertices; };
    int NumEdges() const    /// Number of polytope edges.
    { return num_edges; };
    int NumFacets() const   /// Number of polytope facets.
    { return num_facets; };
    int NumFacetVertices    /// Number of facet vertices of facet \a jf.
    (const FACET_INDEX jf) const
    { return num_facet_vertices[jf]; };
    int VertexCoord         /// \a ic'th vertex coordinate of vertex \a iv.
    (const int iv, const int ic) const
    { return vertex_coord[iv*dimension + ic]; };
    int EdgeEndpoint        /// \a j'th endpoint of edge \a ie. \a j = 0 or 1.
    (const EDGE_INDEX ie, const int j) const
    { return edge_endpoint[int(ie)*2 + j]; };
    POLYTOPE_SHAPE Shape() const   /// Return polytope shape.
    { return shape; }

    /// @brief Return string for poly shape shape
    std::string ShapeString() const
      { return shape_list.String(shape); }

    /// @brief Return enum shape named by shape_str.
    POLYTOPE_SHAPE Shape(const std::string & shape_str) const
    { return shape_list.EnumValue(shape_str); }

    /// @brief Return true if shape is undefined.
    bool IsShapeUndefined() const
    { return shape_list.IsUndefined(shape); }

    // *** OBSOLETE? ***
    int NumNeighbors        /// Number of neighbors of vertex \a iv.
    (const int iv) const
    { return num_vertex_neighbors[iv]; }
    
    // *** OBSOLETE? ***
    int Neighbor            /// Neighbor \a j of vertex \a iv.
    (const int iv, const int j)
    { return vertex_neighbor_list[iv][j]; }

    
    /// @brief Return j (0 or 1) where Endpoint(ie,j) == iv.
    /// @pre iv is an endpoint of edge ie.
    inline int EdgeEndpointIndex
    (const EDGE_INDEX ie, const int iv) const
    {
      if (iv == EdgeEndpoint(ie,0)) { return 0; }
      else { return 1; }
    }

    int NumIncidentEdges    /// Number of edges incident on vertex \a iv.
    (const int iv) const
    { return num_incident_edges[iv]; }

    int IncidentEdge        /// \a j'th edge incident on vertex \a iv.
    (const int iv, const int j) const
    { return incident_edge_list[iv][j]; }
    
    /*!
     *  @brief Return ic'th coordinate of midpoint of edge ie.
     *  - Note: Vertex coordinates must all be even so midpoint coordinate
     *  is an integer.
     */
    int MidpointCoord
    (const EDGE_INDEX ie, const int ic) const;
    
    FACET Facet             /// Bits representing vertices in facet \a jf.
    (const FACET_INDEX jf) const
    { return(facet[jf]); };
    bool IsVertexInFacet    /// Return true if vertex \a iv is in facet \a jf.
    (const FACET_INDEX jf, const int iv) const
    { return((facet[jf] & ((1L) << iv)) != 0); };
    int FacetVertex         /// Return \a k'th vertex in facet \a jf.
    (const FACET_INDEX jf, const int k) const
    { return(facet_vertex_list[jf][k]); };
    
    //@}


    /// @name Check Functions
    //@{

    /// Return true if dimension >= 1.
    bool CheckDimension() const;

    /*!
     *  @brief Check polytope. Return true if polytope passes check.
     *  - Return false if polytope has 0 vertices or edges,
     *    or arrays are not allocated or edge endpoints have
     *    illegal vertex indices, or...
     */
    bool Check(IJK::ERROR & error_msg) const;
    
    //@}

    
    /// @name Print Routines (mainly for debugging).
    //@{
    
    /// Print vertex coordinates (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _VTYPE>
    void PrintVertexCoord(OSTREAM_TYPE & out, const _VTYPE iv) const;

    /// Print vertex index and coordinates (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _VTYPE>
    void PrintVertexIndexAndCoord
    (OSTREAM_TYPE & out, const _VTYPE iv) const;

    /// @overload
    /// @brief Print vertex index and coordinates preceded
    ///   and followed by character strings (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _VTYPE>
    void PrintVertexIndexAndCoord
    (OSTREAM_TYPE & out, const char * prefix, 
     const _VTYPE iv, const char * suffix) const;

    /// @brief Print all vertex coordinates.
    template <typename OSTREAM_TYPE>
    void PrintAllVertexCoord
    (OSTREAM_TYPE & out, const char * line_prefix) const;

    /// Print edge endpoints (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _ETYPE>
    void PrintEdgeEndpoints(OSTREAM_TYPE & out, const _ETYPE ie) const;

    /// Print edge index and endpoints (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _ETYPE>
    void PrintEdgeIndexAndEndpoints
    (OSTREAM_TYPE & out, const _ETYPE ie) const;

    /// @overload
    /// @brief Print edge index and endpoints preceded 
    ///   and followed by character strings (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _ETYPE>
    void PrintEdgeIndexAndEndpoints
    (OSTREAM_TYPE & out, const char * prefix, 
     const _ETYPE ie, const char * suffix) const;

    /// @brief Print all edge endpoints.
    template <typename OSTREAM_TYPE>
    void PrintAllEdgeEndpoints
    (OSTREAM_TYPE & out, const char * line_prefix) const;
    
    /// Print facet vertices (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _FTYPE>
    void PrintFacetVertices
    (OSTREAM_TYPE & out, const _FTYPE jfacet) const;

    /// Print facet index and vertices (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _FTYPE>
    void PrintFacetIndexAndVertices
    (OSTREAM_TYPE & out, const _FTYPE jfacet) const;

    /// @overload
    /// @brief Print facet index and vertices preceded
    ///   and followed by character strings (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _FTYPE>
    void PrintFacetIndexAndVertices
    (OSTREAM_TYPE & out, const char * prefix, 
     const _FTYPE jfacet, const char * suffix) const;

    /// @brief Print vertices of every facet.
    template <typename OSTREAM_TYPE>
    void PrintAllFacetVertices
    (OSTREAM_TYPE & out, const char * line_prefix) const;

    /// @brief Print incident edges (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _VTYPE>
    void PrintIncidentEdges
    (OSTREAM_TYPE & out, const _VTYPE iv) const;

    /// @brief Print vertex index and incident edges
    ///   (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _VTYPE>
    void PrintVertexIndexAndIncidentEdges
    (OSTREAM_TYPE & out, const _VTYPE iv) const;
    
    /// @brief Print vertex index and incident edges preceded
    ///   and followed by character strings (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _VTYPE>
    void PrintVertexIndexAndIncidentEdges
    (OSTREAM_TYPE & out, const char * prefix, 
     const _VTYPE iv, const char * suffix) const;

    /// @brief Print all incident edges.
    template <typename OSTREAM_TYPE>
    void PrintAllIncidentEdges
    (OSTREAM_TYPE & out, const char * line_prefix) const;
    
    //@}

  };

  typedef ISOTABLE_POLY_BASE * ISOTABLE_POLY_BASE_PTR;


  // *****************************************************************
  // MARCHING CUBES ISOSURFACE TABLE POLYTOPE
  // *****************************************************************
  
  /// @brief Marching Cubes isosurface table polytope.
  /// - Public set/generate functions.
  class ISOTABLE_POLY:public ISOTABLE_POLY_BASE {

  public:

    /// Constructor.
    ISOTABLE_POLY(const int _dimension):
      ISOTABLE_POLY_BASE(_dimension) {}

    ISOTABLE_POLY(const ISOTABLE_POLY_BASE & init):
      ISOTABLE_POLY_BASE(init) {};  ///< Copy constructor.

    /// @name Set Functions
    //@{

    /// Set polytope shape.
    void SetShape(const POLYTOPE_SHAPE & shape)
    { ISOTABLE_POLY_BASE::SetShape(shape); };

    /// Set polyhedron dimension.
    void SetDimension(const int _dimension)
    { ISOTABLE_POLY_BASE::SetDimension(_dimension); }

    /// @brief Set number of polyhedron vertices.
    /// - Must be called before setting polyhedron vertices.
    void SetNumVertices(const int numv)
    { ISOTABLE_POLY_BASE::SetNumVertices(numv); }

    /// @brief Set number of polyhedron edges.
    /// - Must be called before setting polyhedron edges.
    void SetNumEdges(const int nume)
    { ISOTABLE_POLY_BASE::SetNumEdges(nume); }

    /// @brief Set number of polyhedron facets.
    /// - Must be called before setting polyhedron facets.
    void SetNumFacets(const int numf)
    { ISOTABLE_POLY_BASE::SetNumFacets(numf); }

    /// Set number of polyhedron vertices, edges and facets.
    void SetSize(const int numv, const int nume, const int numf)
    { SetNumVertices(numv); SetNumEdges(nume); SetNumFacets(numf); };

    /// @brief Set \a ic'th coordinate of vertex \a iv.
    /// @pre SetNumVertices() or SetSize() must be called
    ///   before SetVertexCoord().
    void SetVertexCoord 
    (const int iv, const int ic, const int coord)
    { ISOTABLE_POLY_BASE::SetVertexCoord(iv, ic, coord); }

    /*!
     *  @brief Set endpoints of edge \a ie.
     *  @pre SetNumEdges or SetSize must be called before SetEdge.
     *  @param ie Edge index.  In range [0..NumEdges()-1].
     *  @param iv0 Endpoint 0 index.  In range [0..NumVertices()-1].
     *  @param iv1 Endpoint 1 index.  In range [0..NumVertices()-1].
     */
    void SetEdge(const EDGE_INDEX ie, const int iv0, const int iv1)
    { ISOTABLE_POLY_BASE::SetEdge(ie, iv0, iv1); }

    /// @brief Set number of vertices in facet \a jf.
    /// @pre SetNumFacets() or SetSize() must be called
    ///   before SetNumFacetVertices().
    void SetNumFacetVertices 
    (const FACET_INDEX jf, const int numv)
    { ISOTABLE_POLY_BASE::SetNumFacetVertices(jf, numv); }

    /// @brief Set \a k'th facet vertex of facet \a jf to vertex \a iv.
    /// @pre SetNumFacetVertices(jf) must be called before SetFacetVertex().
    void SetFacetVertex(const FACET_INDEX jf, const int k, const int iv)
    { ISOTABLE_POLY_BASE::SetFacetVertex(jf, k, iv); }

    //@}

    
    /// @name Generate Polyhedron
    //@{

    /// Generate a square, cube or hypercube.
    void GenCube(const int cube_dimension)
    { ISOTABLE_POLY_BASE::GenCube(cube_dimension); }

    /// @brief Generate a square, cube or hypercube.
    /// - Use order of facets and edges given by template class CUBE_FACE_INFO.
    void GenCubeOrderA(const int cube_dimension)
    { ISOTABLE_POLY_BASE::GenCubeOrderA(cube_dimension); }

    /// Generate a triangle, tetrahedron or simplex.      
    void GenSimplex(const int simplex_dimension)
    { ISOTABLE_POLY_BASE::GenSimplex(simplex_dimension); }

    /// Generate a pyramid over a square, cube or hypercube base.
    void GenPyramid(const int pyramid_dimension)
    { ISOTABLE_POLY_BASE::GenPyramid(pyramid_dimension); }

    //@}

    
    /// @name Memory Management Functions
    //@{
    
    /// Free all memory.
    void FreeAll()
    { ISOTABLE_POLY_BASE::FreeAll(); }
    
    //@}
    
  };

  typedef ISOTABLE_POLY * ISOTABLE_POLY_PTR;


  // *****************************************************************
  // MARCHING CUBES ISOTABLE HALF EDGE 3D POLYTOPE
  // *****************************************************************

  /// @brief Marching cubes table 3D polyhedron with half edge support.
  /// - Only for polytope dimension 3.
  class ISOTABLE_HALF_EDGE_POLY_3D:public ISOTABLE_POLY_BASE {

  public:
    
    /// Half edge class.
    class HALF_EDGE {
    public:

      /*!
       *  @brief Half edge index.
       *  - For edge ie directed from EdgeEndpoint(ie,0) 
       *    to EdgeEndpoint(ie,1), ihalf_edge = 2*ie
       *  - For edge ie directed from EdgeEndpoint(ie,1) 
       *    to EdgeEndpoint(ie,0), ihalf_edge = 2*ie+1
       */
      EDGE_INDEX ihalf_edge;

      /// Return half edge index.
      EDGE_INDEX Index() const
      { return ihalf_edge; }
      
      /// Return edge index.
      EDGE_INDEX EdgeIndex() const
      { return (int(ihalf_edge)/2); }

      /// Return half edge from endpoint.
      EDGE_INDEX HalfEdgeFrom() const
      { return (int(ihalf_edge)%2); }

      /// Return half edge to endpoint.
      EDGE_INDEX HalfEdgeTo() const
      { return ((int(ihalf_edge)%2) +1); }
    };
    

  protected:

    /// @brief Index of next half edge in facet.
    HALF_EDGE * next_half_edge_in_facet;

    /// Index of previous half edge in facet.
    HALF_EDGE * prev_half_edge_in_facet;
    
  protected:

    /// Initialize.
    void Init();

    /// @brief Allocate memory.
    /// @pre Number of polyhedron edges already set.
    void Allocate();

    /// Free only memory allocated in this class, not in parent classes.
    void FreeLocal();
    
    /// Free all memory.
    void FreeAll();

  public:

    /// Constructor.
    ISOTABLE_HALF_EDGE_POLY_3D();

    /// Destructor.
    ~ISOTABLE_HALF_EDGE_POLY_3D();

    /// Return next half edge in facet.
    HALF_EDGE NextHalfEdgeInFacet(const HALF_EDGE half_edge)
    { return next_half_edge_in_facet[half_edge.Index()]; }

    /// Return previous half edge in facet.
    HALF_EDGE PrevHalfEdgeInFacet(const HALF_EDGE half_edge)
    { return prev_half_edge_in_facet[half_edge.Index()]; }
        
    /// @brief Check half edge data structures.
    bool CheckHalfEdge(IJK::ERROR & error) const;

    /// Check ISOTABLE_HALF_EDGE_POLY.
    bool Check(IJK::ERROR & error) const;
  };
    

  
  // *****************************************************************
  // MARCHING CUBES ISOSURFACE TABLE 3D CUBE
  // *****************************************************************

  /// Isosurface table 3D cube.
  class ISOTABLE_CUBE_3D:
    public ISOTABLE_POLY_BASE {

  public:

    /// Constructor.
    ISOTABLE_CUBE_3D():ISOTABLE_POLY_BASE(3) {};

    /// Generate 3D cube.
    void GenCube3D()
    { GenCube(3); }
  };


  // *****************************************************************
  // ROUTINES FOR GENERATING POLYHEDRA
  // *****************************************************************

  /*!
   *  @brief Generate a prism with base base_polyhedron
   *  - First numv vertices have last coordinate = 0.
   *  - Last numv vertices have last coordinate = 2.
   *  - First nume edges connect first numv vertices.
   *  - Second nume edges connect second numv vertices.
   *  - Third numv edges connect first to second set of vertices.
   *  - (numv = # vertices in base_polyhedron; nume = # edges in base_polyhedron)
   */
  void generate_prism(const ISOTABLE_POLY_BASE & base_polyhedron,
                      ISOTABLE_POLY & prism);

  
  // *****************************************************************
  // ISOTABLE_POLY_BASE TEMPLATE MEMBER FUNCTIONS
  // *****************************************************************

  // Print vertex coordinates (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _VTYPE>
  void ISOTABLE_POLY_BASE::
  PrintVertexCoord(OSTREAM_TYPE & out, const _VTYPE iv) const
  {
    out << "(";
    for (int d = 0; d < Dimension(); d++) {
      out << VertexCoord(iv, d);
      if (d+1 < Dimension()) { out << ","; }
    }
    out << ")";
  }

  
  // Print vertex index and coordinates (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _VTYPE>
  void ISOTABLE_POLY_BASE::
  PrintVertexIndexAndCoord(OSTREAM_TYPE & out, const _VTYPE iv) const
  {
    out << int(iv) << " ";
    PrintVertexCoord(out, iv);
  }

  
  // Print vertex index and coordinates preceded
  //   and followed by character strings (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _VTYPE>
  void ISOTABLE_POLY_BASE::
  PrintVertexIndexAndCoord
  (OSTREAM_TYPE & out, const char * prefix, 
   const _VTYPE iv, const char * suffix) const
  {
    out << prefix;
    PrintVertexIndexAndCoord(out, iv);
    out << suffix;
  }


  // Print all vertex coordinates.
  template <typename OSTREAM_TYPE>
  void ISOTABLE_POLY_BASE::
  PrintAllVertexCoord
  (OSTREAM_TYPE & out, const char * line_prefix) const
  {
    for (int iv = 0; iv < NumVertices(); iv++) {
      PrintVertexIndexAndCoord(out, line_prefix, iv, "\n");
    }
  }

  
  // Print edge endpoints (mainly for debugging).
    template <typename OSTREAM_TYPE, typename _ETYPE>
  void ISOTABLE_POLY_BASE::
  PrintEdgeEndpoints(OSTREAM_TYPE & out, const _ETYPE ie) const
  {
    out << "(";
    out << int(EdgeEndpoint(ie,0)) << "," << int(EdgeEndpoint(ie,1));
    out << ")";
  }

  
  // Print edge index and endpoints (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _ETYPE>
  void ISOTABLE_POLY_BASE::
  PrintEdgeIndexAndEndpoints
  (OSTREAM_TYPE & out, const _ETYPE ie) const
  {
    out << int(ie) << " ";
    PrintEdgeEndpoints(out, ie);
  }


  //  Print edge index and endpoints preceded 
  //   and followed by character strings (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _ETYPE>
  void ISOTABLE_POLY_BASE::
  PrintEdgeIndexAndEndpoints
  (OSTREAM_TYPE & out, const char * prefix, 
   const _ETYPE ie, const char * suffix) const
  {
    out << prefix;
    PrintEdgeIndexAndEndpoints(out, ie);
    out << suffix;
  }

  
  // @brief Print all edge endpoints.
  template <typename OSTREAM_TYPE>
  void ISOTABLE_POLY_BASE::
  PrintAllEdgeEndpoints
  (OSTREAM_TYPE & out, const char * line_prefix) const
  {
    for (EDGE_INDEX ie = 0; ie < NumEdges(); ie++) {
      PrintEdgeIndexAndEndpoints(out, line_prefix, ie, "\n");
    }
  }

  
  // Print facet vertices (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _FTYPE>
  void ISOTABLE_POLY_BASE::
  PrintFacetVertices
  (OSTREAM_TYPE & out, const _FTYPE jfacet) const
  {
    out << "(";
    for (int iv = 0; iv < NumFacetVertices(jfacet); iv++) {
      out << int(FacetVertex(jfacet, iv));
      if (iv+1 < NumFacetVertices(jfacet)) { out << ","; }
    }
    out << ")";
  }
  

  // Print facet index and vertices (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _FTYPE>
  void ISOTABLE_POLY_BASE::
  PrintFacetIndexAndVertices
  (OSTREAM_TYPE & out, const _FTYPE jfacet) const
  {
    out << int(jfacet) << " ";
    PrintFacetVertices(out, jfacet);
  }


  // Print facet index and vertices preceded
  //   and followed by character strings (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _FTYPE>
  void ISOTABLE_POLY_BASE::
  PrintFacetIndexAndVertices
  (OSTREAM_TYPE & out, const char * prefix, 
   const _FTYPE jfacet, const char * suffix) const
  {
    out << prefix;
    PrintFacetIndexAndVertices(out, jfacet);
    out << suffix;
  }


  // @brief Print all facet vertices.
  template <typename OSTREAM_TYPE>
  void ISOTABLE_POLY_BASE::
  PrintAllFacetVertices
  (OSTREAM_TYPE & out, const char * line_prefix) const
  {
    for (int jfacet = 0; jfacet < NumFacets(); jfacet++) {
      PrintFacetIndexAndVertices
        (out, line_prefix, jfacet, "\n");
    }
  }


  // Print incident edges (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _VTYPE>
  void ISOTABLE_POLY_BASE::
  PrintIncidentEdges
  (OSTREAM_TYPE & out, const _VTYPE iv) const
  {
    out << "(";
    for (int j = 0; j < NumIncidentEdges(iv); j++) {
      out << int(IncidentEdge(iv,j));
      if (j+1 < NumIncidentEdges(iv)) { out << ","; }
    }
    out << ")";
  }

    
  // Print vertex index and incident edges (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _VTYPE>
  void ISOTABLE_POLY_BASE::
    PrintVertexIndexAndIncidentEdges
  (OSTREAM_TYPE & out, const _VTYPE iv) const
  {
    out << int(iv) << " ";
    PrintIncidentEdges(out, iv);
  }

  
  // Print vertex index and incident edges preceded
  //   and followed by character strings (mainly for debugging).
  template <typename OSTREAM_TYPE, typename _VTYPE>
  void ISOTABLE_POLY_BASE::
  PrintVertexIndexAndIncidentEdges
    (OSTREAM_TYPE & out, const char * prefix, 
     const _VTYPE iv, const char * suffix) const
  {
    out << prefix;
    PrintVertexIndexAndIncidentEdges(out, iv);
    out << suffix;
  }

  
  // Print all incident edges.
  template <typename OSTREAM_TYPE>
  void ISOTABLE_POLY_BASE::
  PrintAllIncidentEdges
    (OSTREAM_TYPE & out, const char * line_prefix) const
  {
    for (int iv = 0; iv < NumVertices(); iv++) {
      PrintVertexIndexAndIncidentEdges
        (out, line_prefix, iv, "\n");
    }
  }

}

#endif
