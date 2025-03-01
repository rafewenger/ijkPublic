/*!
 *  @file ijkcube.tpp
 *  @brief ijk templates defining cube (hypercube) classes and functions.
 *  - Version 0.4.1
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

#ifndef _IJKCUBE_
#define _IJKCUBE_

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "ijk.tpp"

namespace IJK {

  // **************************************************
  // TEMPLATE CUBE FACET FUNCTIONS
  // **************************************************

  /// Return orthogonal direction to facet.
  template <typename DTYPE, typename FTYPE>
  inline FTYPE cube_facet_orth_dir(const DTYPE dimension, const FTYPE kf)
  {
    return (kf%dimension);
  }


  /*!
   *  @brief Return side containing facet.
   *  - Return 0 if facet contains vertex 0.
   *  - Return 1 if facet does not contain vertex 0.
   */
  template <typename DTYPE, typename FTYPE>
  inline FTYPE cube_facet_side(const DTYPE dimension, const FTYPE kf)
  {
    return (kf/dimension);
  }
  
  /// Return true if facet jf contains vertex i.
  template <typename DTYPE, typename FTYPE, typename VTYPE>
  bool cube_facet_contains
  (const DTYPE dimension, const FTYPE jf, const VTYPE iv)
  {
    FTYPE side = jf/dimension;
    FTYPE orth_dir = jf%dimension;
    FTYPE mask = (FTYPE(1) << orth_dir);

    bool flag_contains = ((side << orth_dir) == (mask & iv));
    return flag_contains;
  }


  // **************************************************
  // TEMPLATE CLASS CUBE_INFO
  // **************************************************

  /// Base cube class
  template <typename DTYPE, typename NTYPE> 
  class CUBE_INFO {

  protected:
    DTYPE dimension;           ///< Cube dimension
    NTYPE num_vertices;        ///< Number of cube vertices
    NTYPE num_edges;           ///< Number of cube edges.
    NTYPE num_facets;          ///< Number of cube facets.
    NTYPE num_facet_vertices;  ///< Number of cube facet vertices.
    NTYPE num_facet_edges;     ///< Number of cube facet edges.
    NTYPE num_ridge_vertices;  ///< Number of cube ridge vertices.

    void Init                  /// Initialize cube.
    (const DTYPE dimension);
    void FreeAll();            ///< Free all allocated memory.


  public:
    typedef DTYPE DIMENSION_TYPE;         ///< Dimension type.
    typedef NTYPE NUMBER_TYPE;            ///< Number type.

  public:
    // Constructors, destructor.
    CUBE_INFO(const DTYPE dimension);     ///< Constructor.
    CUBE_INFO();                          ///< Constructor.
    ~CUBE_INFO();                         ///< Destructor.

    // set functions
    template <typename DTYPE2>
    void SetDimension                     /// Set cube dimension
    (const DTYPE2 dimension);

    // get functions
    int Dimension() const                 /// Dimension.
    { return dimension; }
    NTYPE NumVertices() const             /// Number of cube vertices. 
    { return num_vertices; }
    NTYPE NumEdges() const                /// Number of cube edges.
    { return num_edges; }
    NTYPE NumFacets() const               /// Number of cube facets. 
    { return num_facets; }
    NTYPE NumFacetVertices() const        /// Number of cube facet vertices. 
    { return num_facet_vertices; }
    NTYPE NumFacetEdges() const           /// Number of cube facet edges. 
    { return num_facet_edges; }
    NTYPE NumRidgeVertices() const        /// Number of cube ridge vertices. 
    { return num_ridge_vertices; }
    NTYPE NumDiagonals() const            /// Number of cube diagonals.
    { return NumFacetVertices(); }

    /// Number of edges incident on a vertex.
    NTYPE NumIncidentEdges() const
    { return Dimension(); }

    /// Maximum cube vertex index. Undefined if dimension < 1.
    NTYPE MaxVertexIndex() const
    { return (NumVertices()-1); }

    /// Index of vertex diagonally opposite iv in cube.
    NTYPE OppositeVertex
    (const NTYPE iv) const
    { return (MaxVertexIndex()-iv); }

    /// Index of facet parallel to ifacet.
    NTYPE OppositeFacet
    (const NTYPE ifacet) const
    { return ((ifacet+Dimension())%NumFacets()); }

    /*!
     *  @brief Return neighbor of vertex \a iv0.
     *  @param dir Edge direction.
     *  @param iv0 Vertex.
     */
    NTYPE VertexNeighbor
    (const NTYPE iv0, const DTYPE dir) const
    {
      int mask = (1L << dir);
      return ((iv0^mask));
    };

    /// @brief Return direction orthogonal to facet.
    int FacetOrthDir(const int ifacet) const
    { return (ifacet%Dimension()); }
    
    /*!
     *  @brief Return side containing facet.
     *  - Return 0 if facet contains vertex 0.
     *  - Return 1 if facet does not contain vertex 0.
     */
    int FacetSide(const int ifacet) const
    { return cube_facet_side(Dimension(), ifacet); }

  };

  // **************************************************
  // TEMPLATE CLASS CUBE_FACE_INFO
  // **************************************************

  /// Information about cube faces.
  template <typename DTYPE, typename NTYPE, typename VTYPE>
  class CUBE_FACE_INFO:public CUBE_INFO<DTYPE,NTYPE> {

  protected:
    VTYPE * facet_vertex;                 ///< Facet vertices.
    VTYPE * facet_edge;                   ///< Facet edges.
    VTYPE * edge_endpoint;                ///< Edge endoints.


    /*!
     *  @brief Edges incident on a vertex.
     *  - incident_edge[iv*NumIncidentEdge()+d] is edge in direction d
     *    incident on vertex iv.
     */
    VTYPE * incident_edge;                

    void Init();
    void Init(const int dimension);
    void FreeAll();

    /// Set facet vertices, edge endpoints and incident edges.
    void SetFaces();

    /// Set facet vertices.
    void SetFacetVertices();

    /*!
     *  @brief Set facet edges.
     *  - for each facet f do:
     *    -# for each direction d that is not orthogonal to f do:
     *      -# f' = lower facet orthogonal to d.
     *      -# for each vertex v of f' do:
     *      -# if v is in f, then
     *        -# ie = edge incident on v in direction d.
     *        -# Add ie to list of edge in facet f.
     *  @pre Arrays facet_vertex[] and incident_edge[] should be set 
     *    before facet_edge[].
     */
    void SetFacetEdges();
    
    /// @brief Set edge endoints.
    /// @pre Array facet_vertex[] should be set before edge_endpoint[].
    void SetEdgeEndpoints();

    /// @brief Set edges incident on each vertex.
    /// @pre Array edge_endpoint should be set before incident_edge[].
    void SetIncidentEdges();

  public:
    typedef VTYPE VERTEX_INDEX;           ///< Vertex index type.

  public:
    CUBE_FACE_INFO();
    CUBE_FACE_INFO(const DTYPE dimension);
    ~CUBE_FACE_INFO() { FreeAll(); };

    // set functions
    template <typename DTYPE2>
    void SetDimension                     ///< Set cube dimension
    (const DTYPE2 dimension);

    /// Return pointer to edge endpoints
    const VTYPE * EdgeEndpoint() const
    { return(edge_endpoint); }

    /*!
     *  @brief Return edge endpoint.
     *  @param ie Edge index.
     *  @param iend Endpoint (0 or 1).
     */
    VTYPE EdgeEndpoint(const VTYPE ie, const NTYPE iend) const
    { return(edge_endpoint[2*ie+iend]); };

    /// @brief Return pointer to facet vertices.
    const VTYPE * FacetVertex() const
    { return facet_vertex; }

    /*!
     *  @brief Return k'th facet vertex.
     *  @param ifacet Facet index. 
     *    - Facet is orthogonal to axis (ifacet\%dimension).
     *    - If ifacet < dimension, facet contains vertex 0.
     *    - If ifacet >= dimension, facet contains largest vertex.
     */
    VTYPE FacetVertex(const NTYPE ifacet, const NTYPE k) const
    { return facet_vertex[ifacet*this->NumFacetVertices()+k]; };

    /// @brief Return pointer to facet edges.
    const VTYPE * FacetEdge() const
    { return facet_edge; }

    /*!
     *  @brief Return k'th facet vertex.
     *  @param ifacet Facet index. 
     *    - Facet is orthogonal to axis (ifacet\%dimension).
     *    - If ifacet < dimension, facet contains vertex 0.
     *    - If ifacet >= dimension, facet contains largest vertex.
     */
    VTYPE FacetEdge(const NTYPE ifacet, const NTYPE k) const
    { return facet_edge[ifacet*this->NumFacetEdges()+k]; };

    /// Index of vertex diagonally opposite vertex k in facet.
    NTYPE OppositeFacetVertex
    (const NTYPE ifacet, const NTYPE k) const
    { return (this->NumFacetVertices()-k-1); }

    /// Return pointer to array of incident edges.
    const VTYPE * IncidentEdge() const
    { return incident_edge; }

    /// Return incident edge with edge direction d.
    VTYPE IncidentEdge(const VTYPE iv, const NTYPE d) const
    { return incident_edge[iv*this->NumIncidentEdges()+d]; };

    /// Return direction orthogonal to facet ifacet.
    const int FacetOrthDir(const VTYPE ifacet) const
    { return (ifacet%this->Dimension()); }

    /// Return direction of edge ie.
    const int EdgeDir(const VTYPE ie) const
    { return (int(ie)/int(this->NumFacetVertices())); }

    /// Return index of facet containing vertex facet_vertex_index.
    VTYPE FacetIndex(const VTYPE facet_vertex_index) const
    { return (int(facet_vertex_index)/int(this->NumFacetVertices())); }

    /*!
     *  @brief Return index of edge "opposite" to ieA in direction dir.
     *  - Edge opposite ieA is translate of ieA in direction dir.
     *  - Edge opposite ieA is contained in facet opposite to facet
     *    containing ieA with orthogonal direction dir.
     *  @pre Direction of edge ieA does not equal dir.
     */
    int OppositeEdge(const VTYPE ieA, const DTYPE dir) const;

    /// Return true if edge is incident on vertex.
    bool IsEdgeIncidentOnVertex(const VTYPE ie, const VTYPE iv) const;

    /// Return true if facet is incident on vertex.
    bool IsFacetIncidentOnVertex(const VTYPE kf, const VTYPE iv) const;


    // Print routines (mainly for debugging.)

    /*!
       @brief Print edge endpoints.
       - Extended version.  Define left/right delimiters and separator.
       - Mainly for debugging.
       @param c0 Left delimiter.
       @param c1 Separator.
       @param c2 Right delimiter.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void PrintEdgeEndpointsX
    (OSTREAM_TYPE & out, const ITYPE2 iedge,
     const char c0, const char c1, const char c2) const;

    /*!
       @brief Print edge endpoints.
              (Default delimiters and separator.)
       - Mainly for debugging.
       - Version using delimiters '(' and ')' and separator ','.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2>
    void PrintEdgeEndpoints
    (OSTREAM_TYPE & out, const ITYPE2 iedge) const
    { PrintEdgeEndpointsX(out, iedge, '(', ',', ')'); }

    /*!
       @brief Print edge endpoints. (Set prefix and suffix.)
       - Mainly for debugging.
       - Version adding prefix and suffix strings.
       @param s0 Prefix string.
       @param s1 Suffix string.
    */
    template <typename OSTREAM_TYPE, typename ITYPE2,
              typename STYPE0, typename STYPE1>
    void PrintEdgeEndpoints
    (OSTREAM_TYPE & out, const STYPE0 & s0, 
     const ITYPE2 ipoly, const STYPE1 & s1) const;

  };


  // **************************************************
  // TEMPLATE CLASS CUBEV
  // **************************************************

  /// Cube vertices.
  template <typename DTYPE, typename NTYPE, typename CTYPE>
  class CUBEV:public CUBE_INFO<DTYPE,NTYPE> {

  protected:
    /// vertex_coord[dimension*k+j] = j'th coordinate of k'th vertex of cube
    CTYPE * vertex_coord;

    /// max_vertex_index.  Stored for faster processing.
    NTYPE max_vertex_index;

    void ZeroLocal();
    void InitLocal();
    void FreeLocal();
    void AllocateLocal();           ///< Allocate data structures in CUBE.
    void SetToUnitCube();           ///< Set coordinates to unit cube.
    void SetMaxVertexIndex();       ///< Set max vertex index.

  public:
    typedef CTYPE COORD_TYPE;       ///< Coordinate type.

  public:
    // Constructors and destructors.
    template <typename CTYPE2>
    CUBEV(const DTYPE dimension, const CTYPE2 * vertex0_coord);
    CUBEV(const DTYPE dimension);
    CUBEV();
    ~CUBEV();

    // set functions
    template <typename DTYPE2>
    void SetDimension               ///< Set cube dimension.
    (const DTYPE2 dimension);

    // *** get functions ***

    NTYPE MaxVertexIndex() const    ///< Maximum cube vertex index.
    { return(max_vertex_index); }

    NTYPE OppositeVertex            ///< Index of vertex opposite iv
    (const NTYPE iv) const
    { return(MaxVertexIndex()-iv); }

    /// Return pointer to vertex coordinates
    const CTYPE * VertexCoord() const
    { return(vertex_coord); }

    /// Return pointer to coordinates of k'th cube vertex
    const CTYPE * VertexCoord(const NTYPE k) const
    { return(vertex_coord+this->Dimension()*k); }

    /// Return j'th coordinate of k'th vertex
    const CTYPE VertexCoord         
    (const NTYPE k, const NTYPE j) const
    { return(vertex_coord[this->Dimension()*k+j]); }


    // *** get diagonal functions ***

    /// Return pointer to coordinates of endpoint 0 of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE * DiagonalEnd0Coord(const NTYPE k) const
    { return(vertex_coord+this->Dimension()*k); }

    /// Return pointer to coordinates of endpoint 1 of diagonal k.
    /// Corner opposite cube corner k is endpoint 1 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE * DiagonalEnd1Coord(const NTYPE k) const
    { return(vertex_coord+this->Dimension()*OppositeVertex(k)); }

    /// Return d'th coordinate of endpoint 0 of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE DiagonalEnd0Coord
    (const NTYPE k, const DTYPE d) const
    { return(vertex_coord[this->Dimension()*k+d]); }

    /// Return d'th coordinate of endpoint 1 of diagonal k.
    /// Cube corner k is endpoint 1 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE DiagonalEnd1Coord
    (const NTYPE k, const DTYPE d) const
    { return(vertex_coord[this->Dimension()*OppositeVertex(k)+d]); }

    /// Return pointer to coordinates of endpoint iend of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE * DiagonalCoord
    (const NTYPE k, const NTYPE iend) const
    { if (iend == 0) { return(DiagonalEnd0Coord(k)); }
      else { return(DiagonalEnd1Coord(k)); }
    }

    /// Return d'th coordinate of endpoint iend of diagonal k.
    /// Cube corner k is endpoint 0 of diagonal k.
    /// Note: Each diagonal is listed twice pointing in opposite directions.
    const CTYPE DiagonalCoord
    (const NTYPE k, const NTYPE iend, const DTYPE d) const
    { if (iend == 0) { return(DiagonalEnd0Coord(k, d)); }
      else { return(DiagonalEnd1Coord(k, d)); }
    }

  };


  // **************************************************
  // TEMPLATE CLASS CUBE
  // **************************************************

  /// Template class cube.  Stores cube vertices and cube edge length.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  class CUBE:public CUBEV<DTYPE,NTYPE,CTYPE> {

  protected:

    /// cube edge length
    LTYPE edge_length;

    void Init();
    template <typename CTYPE2, typename LTYPE2>
    void Init
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);    
    void SetToUnitCube();           ///< Set coordinates to unit cube.

  public:
    typedef CTYPE COORD_TYPE;       ///< Coordinate type.

  public:
    // Constructors and destructors.
    template <typename CTYPE2>
    CUBE(const DTYPE dimension, const CTYPE2 * vertex0_coord);
    template <typename CTYPE2, typename LTYPE2>
    CUBE(const DTYPE dimension, 
         const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);
    CUBE(const DTYPE dimension);
    CUBE();
    ~CUBE();

    // set functions
    template <typename DTYPE2>
    void SetDimension               ///< Set cube dimension.
    (const DTYPE2 dimension);
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set cube vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);
    template <typename DTYPE2, typename CTYPE2, typename LTYPE2>
    void Set                        ///< Set dimension and vertex coord.
    (const DTYPE2 dimension, 
     const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);

    // *** get functions ***

    LTYPE EdgeLength() const
    { return(edge_length); }        ///< Cube edge length.


    /*!
     *  @brief Return true if cube contains point.
     *  @param point_coord[] Point coordinates.
     *  @pre point_coord[] contains Dimension() point coordinates.
     */
    template <typename CTYPE2>
    bool Contains(const CTYPE2 * point_coord) const;

    /*!
     *  @brief Return true if interior of cube contains point.
     *  @param point_coord[] Point coordinates.
     *  @pre point_coord[] contains Dimension() point coordinates.
     */
    template <typename CTYPE2>
    bool InteriorContains(const CTYPE2 * point_coord) const;

  };


  // **************************************************
  // TEMPLATE CLASS HRECT
  // **************************************************

  /// Hyper-rectangle.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  class HRECT:public CUBEV<DTYPE,NTYPE,CTYPE> {

  protected:

    /// hyper-rectangle edge lengths
    LTYPE * edge_length;

    void ZeroLocal();
    void Init();
    void FreeLocal();
    void AllocateLocal();           ///< Allocate data structures in HRECT.
    void SetToUnitCube();           ///< Set coordinates to unit cube.

  public:
    typedef CTYPE COORD_TYPE;       ///< Coordinate type.

  public:
    // Constructors and destructors.
    HRECT(const DTYPE dimension);
    HRECT();
    ~HRECT();

    // set functions
    template <typename DTYPE2>
    void SetDimension               ///< Set hyper-rectangle dimension.
    (const DTYPE2 dimension);
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set hyper-rectangle vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 * edge_length);
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set hyper-rectangle vertex coordinates.
    (const CTYPE2 * vertex0_coord, const std::vector<LTYPE2> & edge_length);


    // *** get functions ***

    template <typename DTYPE2>
    LTYPE EdgeLength(const DTYPE2 d) const
    { return(edge_length[d]); }      ///< Cube edge length d.


    /// Return true if hyper-rectangle contains point.
    /// @param point_coord[] Point coordinates.
    /// @pre point_coord[] contains Dimension() point coordinates.
    template <typename CTYPE2>
    bool Contains(const CTYPE2 * point_coord) const;

    /// Return true if interior of hyper-rectangle contains point.
    /// @param point_coord[] Point coordinates.
    /// @pre point_coord[] contains Dimension() point coordinates.
    template <typename CTYPE2>
    bool InteriorContains(const CTYPE2 * point_coord) const;

  };


  // **************************************************
  // TEMPLATE CLASS UNIT_CUBE
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  class UNIT_CUBE:public CUBE<DTYPE,NTYPE,CTYPE,CTYPE> {

  public:
    // Constructors and destructors.
    UNIT_CUBE(const DTYPE dimension);
    UNIT_CUBE();
    ~UNIT_CUBE() {};

    // get functions
    CTYPE EdgeLength() const
    { return(1); }                  ///< Unit cube edge length.

    // Undefine set functions which do not apply to UNIT_CUBE.

    /// Undefine SetVertexCoord().
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set cube vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);

    /// Undefine Set().
    template <typename DTYPE2, typename CTYPE2, typename LTYPE2>
    void Set                        ///< Set dimension and vertex coord.
    (const DTYPE2 dimension, 
     const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);

  };

  // **************************************************
  // TEMPLATE CLASS CUBE3D
  // **************************************************

  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  class CUBE3D:public CUBE<DTYPE,NTYPE,CTYPE,LTYPE> {

  public:
    CUBE3D():CUBE<DTYPE,NTYPE,CTYPE,LTYPE>(3) {};
    template <typename CTYPE2>
    CUBE3D(const CTYPE2 * vertex0_coord):
      CUBE<DTYPE,NTYPE,CTYPE,LTYPE> (3, vertex0_coord) {};
    template <typename CTYPE2, typename LTYPE2>
    CUBE3D(const CTYPE2 * vertex0_coord, const LTYPE2 edge_length):
      CUBE<DTYPE,NTYPE,CTYPE,LTYPE> (vertex0_coord, edge_length) {};

    // Undefine set functions which do not apply to CUBE3D.

    /// Undefine SetDimension().
    template <typename DTYPE2>
    void SetDimension(const DTYPE2 dimension);

    /// Undefine Set().
    template <typename DTYPE2, typename CTYPE2, typename LTYPE2>
    void Set                        ///< Set dimension and vertex coord.
    (const DTYPE2 dimension, 
     const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);

  };

  // **************************************************
  // TEMPLATE CLASS UNIT_CUBE3D
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  class UNIT_CUBE3D:public CUBE3D<DTYPE,NTYPE,CTYPE,CTYPE> {

  public:
    // Constructors and destructors.
    UNIT_CUBE3D() {};

    // Undefine set functions which do not apply to UNIT_CUBE3D.

    /// Undefine SetVertexCoord().
    template <typename CTYPE2, typename LTYPE2>
    void SetVertexCoord             ///< Set cube vertex coordinates.
    (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length);
  };


  // **************************************************
  // TEMPLATE CLASS CUBE_CENTER_TO_VERTEX_DIRECTIONS
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  class CUBE_CENTER_TO_VERTEX_DIRECTIONS:public CUBEV<DTYPE,NTYPE,CTYPE>
  {
  protected:
    void Init(const DTYPE dimension);

  public:
    CUBE_CENTER_TO_VERTEX_DIRECTIONS();
    CUBE_CENTER_TO_VERTEX_DIRECTIONS(const DTYPE dimension);

    void SetDimension(const DTYPE dimension);

    /// Return pointer to directions.
    const CTYPE * Direction() const
    { return(this->VertexCoord()); }

    /// Return pointer to direction to k'th cube vertex.
    const CTYPE * Direction(const NTYPE k) const
    { return(this->VertexCoord(k)); }

    /// Return j'th coordinate of direction to k'th vertex.
    const CTYPE Direction         
    (const NTYPE k, const NTYPE j) const
    { return(this->VertexCoord(k,j)); }

  };


  // ****************************************************************
  // TEMPLATE CLASS CUBE_CENTER_TO_VERTEX_DIRECTIONS_3D
  // ****************************************************************

  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  class CUBE_CENTER_TO_VERTEX_DIRECTIONS_3D:
    public CUBE_CENTER_TO_VERTEX_DIRECTIONS<DTYPE,NTYPE,CTYPE>
  {
  public:
    CUBE_CENTER_TO_VERTEX_DIRECTIONS_3D():
      CUBE_CENTER_TO_VERTEX_DIRECTIONS<DTYPE,NTYPE,CTYPE>(3) {};

    // Unset SetDimension().
    void SetDimension(const DTYPE dimension);
  };


  // **************************************************
  // TEMPLATE FUNCTIONS: COUNTING
  // **************************************************

  /// Return number of cube vertices
  template <typename DTYPE> 
  long compute_num_cube_vertices(const DTYPE dimension)
  { return(1L << dimension); }

  /// @brief Return number of cube facet vertices. (Fast version.)
  /// @pre dimension > 0.
  template <typename DTYPE>
  long compute_num_cube_facet_verticesF(const DTYPE dimension)
  { return(1L << (dimension-1)); }

  /// @brief Return number of cube ridge vertices. (Fast version.)
  /// @pre dimension > 1.
  template <typename DTYPE>
  long compute_num_cube_ridge_verticesF(const DTYPE dimension)
  { return(1L << (dimension-2)); }

  /// @brief Return number of cube facet vertices.
  template <typename DTYPE> 
  long compute_num_cube_facet_vertices(const DTYPE dimension)
  {
    if (dimension < 1) { return 0; }
    return compute_num_cube_facet_verticesF(dimension);
  }

  /// @brief Return number of cube ridge vertices.
  template <typename DTYPE> 
  long compute_num_cube_ridge_vertices(const DTYPE dimension)
  {
    if (dimension < 2) { return 0; }
    return compute_num_cube_ridge_verticesF(dimension);
  }

  /// Return number of cube facets
  template <typename DTYPE> 
  long compute_num_cube_facets(const DTYPE dimension)
  { return(2*dimension); }

  /// Return number of cube edges
  template <typename DTYPE> 
  long compute_num_cube_edges(const DTYPE dimension)
  { 
    long numv = compute_num_cube_vertices(dimension);
    return((numv*dimension)/2);
  }


  // **************************************************
  // TEMPLATE FUNCTIONS: COMPUTE CUBE DIMENSION
  // **************************************************

  /*!
   *  @brief Compute dimension of cube from number of cube vertices.
   *  @param numv Number of cube vertices.
   *  @pre numv is a power of 2.
   */
  template <typename NTYPE>
  NTYPE compute_cube_dimension_from_num_vertices(const NTYPE numv)
  {
    NTYPE n = numv;
    NTYPE dimension = 0;
    while (n != 0) {
      n = (n >> 1);
      dimension++;
    }

    return(dimension-1);
  }


  // **************************************************
  // TEMPLATE FUNCTIONS: COMPUTE COORDINATES
  // **************************************************

  /*!
   *  @brief Compute cube vertex coordinates.
   *  @param dimension  Dimension of grid.
   *  @param vertex0_coord Coordinates of lowest/leftmost cube vertex.
   *  @param cube_edge_length Cube edge length.
   *  @param[out] coord[] Cube vertex coordinates.
   *  @pre Array coord[] is allocated with size at least 
   *       (number of cube vertices)*dimension
   */
  template <typename DTYPE, typename CTYPE0, typename LTYPE0, typename CTYPE1>
  void compute_cube_vertex_coord
  (const DTYPE dimension, const CTYPE0 * vertex0_coord,
   const LTYPE0 cube_edge_length, CTYPE1 * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_cube_vertex_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(vertex0_coord, "vertex0_coord", error)) 
      { throw error; }
    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    
    const long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long j = 0; j < num_cube_vertices; j++) {
      long j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        int iend = j0 % 2;
        coord[j*dimension+d] = vertex0_coord[d] + iend*cube_edge_length;
        j0 = j0/2;
      }
    }
  }


  /*!
   *  @brief Compute unit cube vertex coordinates (0,0,...,0) to (1,1,...,1)
   *  @param dimension  Dimension of grid.
   *  @param[out] coord[] = Unit cube vertex coordinates.
   *  @pre Array coord[] is allocated with size at least 
   *       (number of cube vertices)*dimension
   */
  template <typename DTYPE, typename CTYPE>
  void compute_unit_cube_vertex_coord(const DTYPE dimension, CTYPE * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_unit_cube_vertex_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    
    const long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long j = 0; j < num_cube_vertices; j++) {
      long j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        coord[j*dimension+d] = j0 % 2;
        j0 = j0/2;
      }
    }
  }


  /*!
   *  @brief Compute hyper-rectangle vertex coordinates.
   *  @param dimension  Dimension of grid.
   *  @param vertex0_coord Coordinates of lowest/leftmost hyper-rectangle vertex.
   *  @param hyper_rectangle[d] Length of d'th edge of hyper-rectangle.
   *  @param[out] coord[] Hyper-rectangle vertex coordinates.
   *  @pre Array coord[] is allocated with size at least 
   *       (number of hyper_rectangle vertices)*dimension
   */
  template <typename DTYPE, typename CTYPE0, typename LTYPE0, typename CTYPE1>
  void compute_hrect_vertex_coord
  (const DTYPE dimension, const CTYPE0 * vertex0_coord,
   const LTYPE0 * edge_length, CTYPE1 * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_hrect_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(vertex0_coord, "vertex0_coord", error)) 
      { throw error; }
    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    if (!check_array_allocated(edge_length, "edge_length", error)) 
      { throw error; }
    
    const long num_hrect_vertices = compute_num_cube_vertices(dimension);

    for (long j = 0; j < num_hrect_vertices; j++) {
      long j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        int iend = j0 % 2;
        coord[j*dimension+d] = vertex0_coord[d] + iend*edge_length[d];
        j0 = j0/2;
      }
    }
  }


  /*!
   *  @brief Compute direction of each vertex from cube center.
   *  @param dimension  Dimension of grid.
   *  @param[out] coord[] = Unit cube vertex coordinates.
   *  @pre Array coord[] is allocated with size at least 
   *       (number of cube vertices)*dimension
   */
  template <typename DTYPE, typename CTYPE>
  void compute_cube_center_to_vertex_directions
  (const DTYPE dimension, CTYPE * direction)
  {
    IJK::PROCEDURE_ERROR error("compute_cube_directions");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(direction, "direction", error)) 
      { throw error; }
    
    const long num_cube_vertices = compute_num_cube_vertices(dimension);
    const CTYPE pos_coord = 1.0/std::sqrt(CTYPE(dimension));

    for (long j = 0; j < num_cube_vertices; j++) {
      long j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        direction[j*dimension+d] = (2*(j0%2)-1)*pos_coord;
        j0 = j0/2;
      }
    }
  }


  /*!
   *  @brief Compute coordinates of cube diagonal endpoints.
   *  @param dimension  Dimension of grid.
   *  @param vertex0_coord Coordinates of lowest/leftmost cube vertex.
   *  @param cube_edge_length Cube edge length.
   *  @param[out] coord[] = Unit cube diagonal coordinates.
   *    coord[(2*k+i)*dimension+j] = 
   *      j'th coordinate of endpoint i of diagonal k.
   *  @pre Array coord[] is allocated with size at least 
   *       2*(number of cube vertices)*dimension
   */
  template <typename DTYPE, typename CTYPE0, 
            typename LTYPE0, typename CTYPE1>
  void compute_cube_diagonal_coord
  (const DTYPE dimension, const CTYPE0 * vertex0_coord,
   const LTYPE0 cube_edge_length, CTYPE1 * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_cube_diagonal_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(vertex0_coord, "vertex0_coord", error)) 
      { throw error; }
    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    
    const long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long k = 0; k < num_cube_vertices; k++) {
      long k0 = k;
      for (DTYPE d = 0; d < dimension; d++) {
        int iend = k0 % 2;
        // Diagonal k is pointing away from cube vertex k.
        coord[(2*k)*dimension+d] = vertex0_coord[d] + iend*cube_edge_length;
        // Diagonal k is pointing toward the opposite vertex.
        coord[(2*k+1)*dimension+d] = 
          vertex0_coord[d] + (1-iend)*cube_edge_length;
        k0 = k0/2;
      }
    }

  }


  /*!
   *  @brief Compute coordinates of unit cube diagonal endpoints.
   *  @param dimension  Dimension of grid.
   *  @param[out] coord[] = Cube diagonal coordinates.
   *    coord[(2*k+i)*dimension+j] = 
   *      j'th coordinate of endpoint i of diagonal k.
   *  @pre Array coord[] is allocated with size at least 
   *       2*(number of cube vertices)*dimension
   */
  template <typename DTYPE, typename CTYPE>
  void compute_unit_cube_diagonal_coord
  (const DTYPE dimension, CTYPE * coord)
  {
    IJK::PROCEDURE_ERROR error("compute_unit_cube_diagonal_coord");

    if (dimension <= 0) { return; };

    if (!check_array_allocated(coord, "coord", error)) { throw error; }
    
    const long num_cube_vertices = compute_num_cube_vertices(dimension);

    for (long k = 0; k < num_cube_vertices; k++) {
      long k0 = k;
      for (DTYPE d = 0; d < dimension; d++) {
        CTYPE c = k0 % 2;
        // Diagonal k is pointing away from cube vertex k.
        coord[(2*k)*dimension+d] = c;
        // Diagonal k is pointing toward the opposite vertex.
        coord[(2*k+1)*dimension+d] = 1-c;
        k0 = k0/2;
      }
    }

  }


  /*!
   *  @brief Compute Coxeter-Freudenthal-Kuhn triangulation of the cube.
   *  @pre dimension >= 2.
   *  @param dimension Cube dimension.
   *  @param cube_vertices[] List of cube vertices in coordinate vertex order,
   *    increasing x, then increasing y, then ...
   *    - Note: Order of x,y,z and increasing or decreasing is arbitrary.
   *      Could also be increasing y, then increasing x, then ...
   *      or increasing y, then decreasing z, then decreasing x, then ...
   *    @pre cube_vertices[] contains exactly 2^(dimension) vertices.
   *  @param[out] Add simplex vertices to simplex_vertex_list[].
   *    - Note: Does not clear simplex_vertex_list[].
   */
  template <typename DTYPE, typename VTYPE, typename VTYPES>
  void compute_cube_CFK_triangulation
  (const DTYPE dimension, const VTYPE cube_vertex[],
   std::vector<VTYPES> & simplex_vertex_list)
  {
    typedef typename std::vector<VTYPES>::size_type SIZE_TYPE;

    const int DIM2(2);
    const int num_vert_per_simplex = dimension+1;

    if (dimension == DIM2) {
      simplex_vertex_list.push_back(cube_vertex[0]);
      simplex_vertex_list.push_back(cube_vertex[1]);
      simplex_vertex_list.push_back(cube_vertex[3]);
      simplex_vertex_list.push_back(cube_vertex[0]);
      simplex_vertex_list.push_back(cube_vertex[2]);
      simplex_vertex_list.push_back(cube_vertex[3]);
    }
    else if (dimension > DIM2) {
      const CUBE_FACE_INFO<int,int,VTYPE> cube(dimension);
      const int num_facet_vertices = cube.NumFacetVertices();
      
      for (DTYPE d = 0; d < dimension; d++) {
        const VTYPE * facet_vertex_list =
          cube.FacetVertex() + d*num_facet_vertices;
        std::vector<VTYPES> triangulated_facet_vertex_list;

        compute_cube_CFK_triangulation
          (dimension-1, facet_vertex_list,
           triangulated_facet_vertex_list);

        // Add simplices to cube.
        const int num_vert_per_facet_simplex = num_vert_per_simplex-1;
        const SIZE_TYPE num_simplices_in_facet =
          triangulated_facet_vertex_list.size()/num_vert_per_facet_simplex;
        for (SIZE_TYPE jsimplex = 0; jsimplex < num_simplices_in_facet;
             jsimplex++) {
          const VTYPES * simplexA_vert =
            triangulated_facet_vertex_list.data() +
            jsimplex*num_vert_per_facet_simplex;

          for (int j = 0; j < num_vert_per_facet_simplex; j++) {
            const VTYPE jv = cube_vertex[simplexA_vert[j]];
            simplex_vertex_list.push_back(jv);
          }
          simplex_vertex_list.push_back(cube_vertex[cube.NumVertices()-1]);
        }
      }
    }
    else {
      // dimension < 2.
      IJK::PROCEDURE_ERROR error("compute_cube_CFK_triangulation");

      error.AddMessage
        ("Programming error. Illegal dimension ", dimension, ".");
      error.AddMessage
        ("  Dimension must be at least two.");
      throw error;
    }
  }


  /*!
   *  @overload
   *  @brief Compute Coxeter-Freudenthal-Kuhn triangulation of the cube.
   *  - Version computing CFK triangulation on the unit cube.
   *  @param[out] Add simplex vertices to simplex_vertex_list[].
   *    - Note: Does not clear simplex_vertex_list[].
   */
  template <typename DTYPE, typename VTYPE>
  void compute_cube_CFK_triangulation
  (const DTYPE dimension, std::vector<VTYPE> & simplex_vertex_list)
  {
    const VTYPE num_cube_vertices =
      compute_num_cube_vertices(dimension);
    std::vector<VTYPE> cube_vertex(num_cube_vertices);

    for (int i = 0; i < num_cube_vertices; i++)
      { cube_vertex[i] = i; }

    compute_cube_CFK_triangulation
      (dimension, cube_vertex.data(), simplex_vertex_list);
  }
    
  
  // **************************************************
  // TEMPLATE FUNCTION: VERTEX NEIGHBOR
  // **************************************************

  /// Compute the neighbor of cube vertex iv in given direction.
  template <typename VTYPE1, typename DTYPE, typename VTYPE2>
  void compute_cube_vertex_neighbor
  (const VTYPE1 iv1, const DTYPE d, VTYPE2 & iv2)
    {
      VTYPE1 mask = (VTYPE1(1) << d);
      iv2 = (iv1^mask);
    }

  template <typename NTYPE, typename FTYPE1, typename FTYPE2>
  void compute_opposite_cube_facet
  (const NTYPE num_facets, const FTYPE1 kf, FTYPE2 & kf2)
  { 
    kf2 = (kf+(num_facets/2))%num_facets;
  }


  // **************************************************
  // TEMPLATE BOUNDARY BIT FUNCTIONS
  // **************************************************


  /// Return index of boundary bit.
  template <typename DTYPE, typename STYPE, typename ITYPE>
  inline void compute_boundary_bit_index
  (const DTYPE orth_dir, const STYPE side, ITYPE & bit_index)
  {
    bit_index = 2*orth_dir + DTYPE(side);
  }


  // **************************************************
  // TEMPLATE FUNCTIONS: COMPUTE CORNER IN DIRECTION
  // **************************************************

  /// Compute cube corner nearest to given direction.
  template <typename DTYPE, typename CTYPE, typename NTYPE>
  void compute_corner_nearest_direction
  (const DTYPE dimension, const CTYPE * dir, NTYPE & icorner)
  {
    icorner = 0;
    NTYPE bit_flag = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (dir[d] > 0) 
        { icorner = (icorner | bit_flag); }
      bit_flag = (bit_flag << 1);
    }
  }

  // **************************************************
  // TEMPLATES TO CHECK VALUES
  // **************************************************

  /// @brief Check dimension.
  /// - Return true if dimension is non-negative.
  template <typename DTYPE>
  bool check_cube_dimension(const DTYPE dimension, IJK::ERROR & error)
  {
    if (dimension < 0) {
      error.AddMessage("Illegal dimension ", dimension, ".");
      error.AddMessage("Dimension must be non-negative.");
      return(false);
    }

    return(true);
  }

  // **************************************************
  // TEMPLATE CLASS CUBE_INFO MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE> 
  CUBE_INFO<DTYPE,NTYPE>::CUBE_INFO(const DTYPE dimension)
  {
    Init(dimension);
  }

  /// Constructor
  template <typename DTYPE, typename NTYPE> 
  CUBE_INFO<DTYPE,NTYPE>::CUBE_INFO()
  {
    Init(0);
  }

  /// Destructor
  template <typename DTYPE, typename NTYPE> 
  CUBE_INFO<DTYPE,NTYPE>::~CUBE_INFO()
  {
    dimension = 0;
    num_vertices = 0;
    num_edges = 0;
    num_facets = 0;
    num_facet_vertices = 0;
    num_facet_edges = 0;
    num_ridge_vertices = 0;
  }

  /// Initialize
  template <typename DTYPE, typename NTYPE> 
  void CUBE_INFO<DTYPE,NTYPE>::Init(const DTYPE dimension)
  {
    SetDimension(dimension);
  }

  /// Set cube dimension
  template <typename DTYPE, typename NTYPE> 
  template <typename DTYPE2>
  void CUBE_INFO<DTYPE,NTYPE>::SetDimension(const DTYPE2 dimension)
  {
    this->dimension = dimension;
    this->num_vertices = compute_num_cube_vertices(dimension);
    this->num_facets = compute_num_cube_facets(dimension);
    this->num_edges = compute_num_cube_edges(dimension);
    this->num_facet_vertices = num_vertices/2;
    this->num_ridge_vertices = num_facet_vertices/2;

    if (dimension > 0) {
      this->num_facet_edges =
        (this->num_edges - this->num_facet_vertices)/2;
    }
    else
      { this->num_facet_edges = 0; }
  }

  // **************************************************
  // TEMPLATE CLASS CUBE_FACE_INFO MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  CUBE_FACE_INFO(const DTYPE dimension):CUBE_INFO<DTYPE,NTYPE>(dimension)
  {
    Init(dimension);
  }

  /// Constructor
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::CUBE_FACE_INFO()
  {
    Init();
  }

  // Initialize
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::Init()
  {
    facet_vertex = NULL;
    facet_edge = NULL;
    edge_endpoint = NULL;
    incident_edge = NULL;
  }

  // Initialize
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::Init(const int dimension)
  {
    facet_vertex = NULL;
    facet_edge = NULL;
    edge_endpoint = NULL;
    incident_edge = NULL;

    SetDimension(dimension);
  }

  // Free all arrays
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::FreeAll()
  {
    if (facet_vertex != NULL) {
      delete [] facet_vertex;
      facet_vertex = NULL;
    }

    if (facet_edge != NULL) {
      delete [] facet_edge;
      facet_edge = NULL;
    }    

    if (edge_endpoint != NULL) {
      delete [] edge_endpoint;
      edge_endpoint = NULL;
    }

    if (incident_edge != NULL) {
      delete [] incident_edge;
      incident_edge = NULL;
    }
  }


  // Set dimension.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  template <typename DTYPE2>
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetDimension
  (const DTYPE2 dimension)
  {
    FreeAll();
    CUBE_INFO<DTYPE,NTYPE>::SetDimension(dimension);
    facet_vertex = new VTYPE[this->NumFacetVertices()*this->NumFacets()];
    facet_edge = new VTYPE[this->NumFacetEdges()*this->NumFacets()];
    edge_endpoint = new VTYPE[this->NumEdges()*2];
    incident_edge = new VTYPE[this->NumVertices()*this->NumIncidentEdges()];

    SetFaces();
  }

  // Set facet vertices, edge endpoints and incident edges.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetFaces()
  {    
    SetFacetVertices();

    // Must call SetFacetVertices() before SetEdgeEndpoints().
    SetEdgeEndpoints();

    // Must call SetEdgeEndpoints() before SetFacetVertices().
    SetIncidentEdges();

    // Must call SetFacetVertices() and SetIncidentEdges()
    //   before SetFacetEdges().
    SetFacetEdges();
  }

  
  // Set facet vertices.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetFacetVertices()
  {
    const DTYPE dimension = this->Dimension();
    const NTYPE num_cube_vertices = this->NumVertices();
    const NTYPE num_facet_vertices = this->NumFacetVertices();
    IJK::PROCEDURE_ERROR error("SetFacetVertices");

    if (dimension < 1) { return; };

    if (facet_vertex == NULL) {
      error.AddMessage("Programming error.  Array facet_vertex[] not allocated.");
      throw error;
    }

    for (NTYPE ifacet = 0; ifacet < this->NumFacets(); ifacet++) {

      /* OBSOLETE
      DTYPE orth_dir = ifacet%dimension;
      */
      const int orth_dir = FacetOrthDir(ifacet);

      // Multiply by s mod (num_cube_vertices-1) to compute vertex index.
      NTYPE s = 2;
      for (int d = 0; d < orth_dir; d++)
        { s = s*2; };
      s = s%(num_cube_vertices-1);

      NTYPE i0 = ifacet*num_facet_vertices;

      if (ifacet < dimension) {

        for (NTYPE i = 0; i < num_facet_vertices; i++) 
          { facet_vertex[i0+i] = (i*s)%(num_cube_vertices-1); }
      }
      else {

        // Translate by t mod num_cube_vertices to compute vertex index.
        NTYPE t = 1;
        for (DTYPE d = 0; d < orth_dir; d++)
          { t = t*2; };

        for (NTYPE i = 0; i < num_facet_vertices; i++) 
          { facet_vertex[i0+i] = ((i*s)%(num_cube_vertices-1))+t; }

        // Swap subfacets to get consistent facet orientation.
        if (num_facet_vertices > 1) {
          NTYPE num_subfacet_vertices = num_facet_vertices/2;
          for (NTYPE i = 0; i < num_subfacet_vertices; i++) {
            std::swap(facet_vertex[i0+i], 
                      facet_vertex[i0+i+num_subfacet_vertices]);
          }
        }
      }
    }

  }


  // Set facet edges
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetFacetEdges()
  {
    const DTYPE dimension = this->Dimension();
    const NTYPE num_facet_vertices = this->NumFacetVertices();
    const NTYPE num_facet_edges = this->NumFacetEdges();
    IJK::PROCEDURE_ERROR error("SetFacetEdges");

    if (dimension < 1) { return; };

    if (facet_edge == NULL) {
      error.AddMessage("Programming error.  Array facet_edge[] not allocated.");
      throw error;
    }

    for (NTYPE ifacet0 = 0; ifacet0 < this->NumFacets(); ifacet0++) {

      const int orth_dir = FacetOrthDir(ifacet0);

      int ikount = 0;
      VTYPE * facet0_edge = facet_edge + this->NumFacetEdges()*ifacet0;
      for (int edge_dir = 0; edge_dir < this->Dimension(); edge_dir++) {

        if (edge_dir == orth_dir) {
          // Edge with direction edge_dir is not in facet.
          continue;
        }

        // ifacet1 is lower facet orthogonal to edge_dir.
        const int ifacet1 = edge_dir;
        for (int j = 0; j < num_facet_vertices; j++) {
          const VTYPE jv = FacetVertex(ifacet1, j);
          if (this->IsFacetIncidentOnVertex(ifacet0, jv)) {
            const VTYPE iedge = IncidentEdge(jv, edge_dir);
            facet0_edge[ikount] = iedge;
            ikount++;
          }
        }
      }

      if (ikount != num_facet_edges) {
        error.AddMessage
          ("Programming error in determining edges for facet ",
           ifacet0, ".");
        error.AddMessage("  Dimension: ", dimension);
        error.AddMessage("  Added ", ikount, " edges.");
        error.AddMessage
          ("  Expected ", this->NumFacetEdges(), " edges.");
        throw error;
      }
    }
  }

      
  // Set edge endpoints.
  // @pre Array facet_vertex[] should be set before edge_endpoints[].
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetEdgeEndpoints()
  {
    const DTYPE dimension = this->Dimension();
    IJK::PROCEDURE_ERROR error("SetEdgeEndpoints");

    if (dimension < 1) { return; };

    if (edge_endpoint == NULL) {
      error.AddMessage
        ("Programming error.  Array edge_endpoint[] not allocated.");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      for (NTYPE i = 0; i < this->NumFacetVertices(); i++) {
        const NTYPE je = d*this->NumFacetVertices()+i;
        const VTYPE iv0 = FacetVertex(d, i);
        edge_endpoint[2*je] = iv0;
        const VTYPE iv1 = this->VertexNeighbor(iv0, d);
        edge_endpoint[2*je+1] = iv1;
      }
    }

  }

  // Set edges incident on vertices.
  // @pre Array facet_vertex[] should be set before edge_endpoints[].
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::SetIncidentEdges()
  {
    const DTYPE dimension = this->Dimension();
    IJK::PROCEDURE_ERROR error("SetIncidentEdges");

    if (dimension < 1) { return; };

    if (incident_edge == NULL) {
      error.AddMessage
        ("Programming error.  Array incident_edge[] not allocated.");
      throw error;
    }

    for (VTYPE ie = 0; ie < this->NumEdges(); ie++) {
      VTYPE iend0 = EdgeEndpoint(ie, 0);
      VTYPE iend1 = EdgeEndpoint(ie, 1);
      if (iend0 > iend1) { std::swap(iend0,iend1); }

      DTYPE edge_dir = EdgeDir(ie);
      incident_edge[iend0*this->NumIncidentEdges()+edge_dir] = ie;
      incident_edge[iend1*this->NumIncidentEdges()+edge_dir] = ie;
    }
  }

  
  // Return index of edge "opposite" to ieA in direction dir.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  int CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  OppositeEdge (const VTYPE ieA, const DTYPE dir) const
  {
    const DTYPE edge_dir = EdgeDir(ieA);
    const VTYPE ivA0 = EdgeEndpoint(ieA, 0);
    const VTYPE ivB0 = this->VertexNeighbor(ivA0, dir);
    const VTYPE ieB = IncidentEdge(ivB0, edge_dir);

    return ieB;
  }


  // Return true if edge is incident on vertex.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  bool CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  IsEdgeIncidentOnVertex(const VTYPE ie, const VTYPE iv) const
  {
    if (iv == EdgeEndpoint(ie, 0)) { return(true); }
    if (iv == EdgeEndpoint(ie, 1)) { return(true); }
    return(false);
  }


  // Return true if facet is incident on vertex
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  bool CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  IsFacetIncidentOnVertex(const VTYPE kf, const VTYPE iv) const
  {
    return (cube_facet_contains(this->Dimension(), kf, iv));
  }


  // Print edge endpoints.
  // - Extended version.  Define left/right delimiters and separator.
  // @param c0 Left delimiter.
  // @param c1 Separator.
  // @param c2 Right delimiter.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  template <typename OSTREAM_TYPE, typename ITYPE2>
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  PrintEdgeEndpointsX
  (OSTREAM_TYPE & out, const ITYPE2 iedge,
   const char c0, const char c1, const char c2) const
  {
    out << c0 << EdgeEndpoint(iedge, 0)
	<< c1 << EdgeEndpoint(iedge, 1) << c2;
  }


  // Print edge endpoints.
  template <typename DTYPE, typename NTYPE, typename VTYPE> 
  template <typename OSTREAM_TYPE, typename ITYPE2,
            typename STYPE0, typename STYPE1>
  void CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  PrintEdgeEndpoints
    (OSTREAM_TYPE & out, const STYPE0 & s0, 
     const ITYPE2 iedge, const STYPE1 & s1) const
  {
    out << s0;
    PrintEdgeEndpoints(out, iedge);
    out << s1;
  }


  // **************************************************
  // TEMPLATE CLASS CUBEV MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE, typename CTYPE>
  CUBEV<DTYPE,NTYPE,CTYPE>::CUBEV(const DTYPE dimension):
    CUBE_INFO<DTYPE,NTYPE>(dimension)
  {
    InitLocal();
  }

  /// Constructor.
  template <typename DTYPE, typename NTYPE, typename CTYPE>
  CUBEV<DTYPE,NTYPE,CTYPE>::CUBEV():CUBE_INFO<DTYPE,NTYPE>()
  {
    InitLocal();
  }

  /// Destructor
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  CUBEV<DTYPE,NTYPE,CTYPE>::~CUBEV()
  {
    FreeLocal();
  }

  /// Initialize.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  void CUBEV<DTYPE,NTYPE,CTYPE>::InitLocal()
  {
    ZeroLocal();
    AllocateLocal();
    SetMaxVertexIndex();
    SetToUnitCube();
  }

  /// Allocate local data structures in CUBEV
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  void CUBEV<DTYPE,NTYPE,CTYPE>::AllocateLocal()
  {
    IJK::PROCEDURE_ERROR error("CUBEV::AllocateLocal");

    if (!check_cube_dimension(this->Dimension(), error)) 
      { throw error; }
    if (!IJK::check_is_NULL(vertex_coord, "vertex_coord", error))
      { throw error; }

    NTYPE numv = this->NumVertices();
    this->vertex_coord = new CTYPE[numv*this->Dimension()];
  }

  // Set max vertex index.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  void CUBEV<DTYPE,NTYPE,CTYPE>::SetMaxVertexIndex()
  {
    if (this->NumVertices() > 0) 
      { max_vertex_index = this->NumVertices()-1; }
    else
      { max_vertex_index = 0; }
  }

  // Set vertex coordinates to unit cube vertex coordinates.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  void CUBEV<DTYPE,NTYPE,CTYPE>::SetToUnitCube()
  {
    IJK::PROCEDURE_ERROR error("CUBE::SetToUnitCube");

    if (!check_cube_dimension(this->Dimension(), error)) 
      { throw error; }
    if (!check_array_allocated(vertex_coord, "vertex_coord", error)) 
      { throw error; }

    compute_unit_cube_vertex_coord
      (this->Dimension(), this->vertex_coord);
  }

  /// Set pointers defined in CUBE to NULL.
  template <typename DTYPE, typename NTYPE, typename CTYPE>
  void CUBEV<DTYPE,NTYPE,CTYPE>::ZeroLocal()
  {
    max_vertex_index = 0;
    this->vertex_coord = NULL;
  }


  /// Free memory allocated in CUBEV.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  void CUBEV<DTYPE,NTYPE,CTYPE>::FreeLocal()
  {
    if (vertex_coord != NULL) { delete [] vertex_coord; };
    ZeroLocal();
  }

  /// Set cube dimension.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  template <typename DTYPE2>
  void CUBEV<DTYPE,NTYPE,CTYPE>::SetDimension
  (const DTYPE2 dimension)
  {
    FreeLocal();

    CUBE_INFO<DTYPE,NTYPE>::SetDimension(dimension);
    AllocateLocal();
    SetMaxVertexIndex();
    SetToUnitCube();
  }

  // **************************************************
  // TEMPLATE CLASS CUBE MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::CUBE(const DTYPE dimension):
    CUBEV<DTYPE,NTYPE,CTYPE>(dimension)
  {
    Init();
  }

  /// Constructor.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::CUBE():CUBEV<DTYPE,NTYPE,CTYPE>()
  {
    Init();
  }

  /// Destructor
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::~CUBE()
  {}

  /// Initialize.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::Init()
  {
    edge_length = 1;
    CUBEV<DTYPE,NTYPE,CTYPE>::InitLocal();
  }

  // Return true if cube contains coordinate.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2>
  bool CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::
  Contains(const CTYPE2 * point_coord) const
  {
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      COORD_TYPE c = this->VertexCoord(0,d);
      if (point_coord[d] < c)
        { return(false); }

      if (point_coord[d] > c + EdgeLength())
        { return(false); }
    }

    return(true);
  }

  // Return true if interior of cube contains coordinate.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2>
  bool CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::
  InteriorContains(const CTYPE2 * point_coord) const
  {
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      COORD_TYPE c = this->VertexCoord(0,d);
      if (point_coord[d] <= c)
        { return(false); }

      if (point_coord[d] >= c + EdgeLength())
        { return(false); }
    }

    return(true);
  }

  // Set vertex coordinates to unit cube vertex coordinates.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::SetToUnitCube()
  {
    edge_length = 1;
    CUBEV<DTYPE,NTYPE,CTYPE>::SetToUnitCube();
  }

  /// Set cube dimension.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename DTYPE2>
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::SetDimension
  (const DTYPE2 dimension)
  {
    CUBEV<DTYPE,NTYPE,CTYPE>::SetDimension(dimension);
    edge_length = 1;
  }

  /// Set cube vertex coordinates.
  /// @param vertex0_coord Coordinates of lowest/leftmost vertex.
  /// @param edge_length Length of cube edge.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2, typename LTYPE2>
  void CUBE<DTYPE,NTYPE,CTYPE,LTYPE>::
  SetVertexCoord             ///< Set cube vertex coordinates.
  (const CTYPE2 * vertex0_coord, const LTYPE2 edge_length)
  {
    this->edge_length = edge_length;
    compute_cube_vertex_coord
      (this->Dimension(), vertex0_coord, edge_length, this->vertex_coord);
  }


  // **************************************************
  // TEMPLATE CLASS HRECT MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::HRECT(const DTYPE dimension):
    CUBEV<DTYPE,NTYPE,CTYPE>(dimension)
  {
    Init();
  }

  /// Constructor.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::HRECT():CUBEV<DTYPE,NTYPE,CTYPE>()
  {
    Init();
  }

  /// Destructor
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::~HRECT()
  {
    FreeLocal();
  }

  /// Initialize.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::Init()
  {
    ZeroLocal();
    CUBEV<DTYPE,NTYPE,CTYPE>::InitLocal();
    AllocateLocal();
    SetToUnitCube();
  }

  /// Set pointers defined in HRECT to NULL.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE>
  void HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::ZeroLocal()
  {
    this->edge_length = NULL;
  }

  /// Allocate local data structures in HRECT.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::AllocateLocal()
  {
    IJK::PROCEDURE_ERROR error("HRECT::AllocateLocal");

    if (!check_cube_dimension(this->Dimension(), error)) 
      { throw error; }
    if (!IJK::check_is_NULL(edge_length, "edge_length", error))
      { throw error; }

    this->edge_length = new CTYPE[this->Dimension()];
  }

  /// Free memory allocated in HRECT
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::FreeLocal()
  {
    if (edge_length != NULL) { delete [] edge_length; };
    ZeroLocal();
  }

  // Return true if cube contains coordinate.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2>
  bool HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::
  Contains(const CTYPE2 * point_coord) const
  {
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      COORD_TYPE c = VertexCoord(0,d);
      if (point_coord[d] < c)
        { return(false); }

      if (point_coord[d] > c + EdgeLength(d))
        { return(false); }
    }

    return(true);
  }

  // Return true if interior of cube contains coordinate.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2>
  bool HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::
  InteriorContains(const CTYPE2 * point_coord) const
  {
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      COORD_TYPE c = VertexCoord(0,d);
      if (point_coord[d] <= c)
        { return(false); }

      if (point_coord[d] >= c + EdgeLength(d))
        { return(false); }
    }

    return(true);
  }

  // Set vertex coordinates to unit cube vertex coordinates.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  void HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::SetToUnitCube()
  {
    for (DTYPE d = 0; d < this->Dimension(); d++)
      { edge_length[d] = 1; }
    CUBEV<DTYPE,NTYPE,CTYPE>::SetToUnitCube();
  }

  /// Set cube dimension.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename DTYPE2>
  void HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::SetDimension
  (const DTYPE2 dimension)
  {
    CUBEV<DTYPE,NTYPE,CTYPE>::SetDimension(dimension);
    for (DTYPE d = 0; d < this->Dimension(); d++)
      { edge_length[d] = 1; }
  }

  /// Set cube vertex coordinates.
  /// @param vertex0_coord Coordinates of lowest/leftmost vertex.
  /// @param edge_length[d] Length of d'th hyper-rectangle edge.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2, typename LTYPE2>
  void HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::
  SetVertexCoord             ///< Set cube vertex coordinates.
  (const CTYPE2 * vertex0_coord, const LTYPE2 * edge_length)
  {
    for (DTYPE d = 0; d < this->Dimension(); d++)
      { this->edge_length[d] = edge_length[d]; }

    compute_hrect_vertex_coord
      (this->Dimension(), vertex0_coord, this->edge_length, this->vertex_coord);
  }

  /// Set cube vertex coordinates.
  template <typename DTYPE, typename NTYPE, 
            typename CTYPE, typename LTYPE> 
  template <typename CTYPE2, typename LTYPE2>
  void HRECT<DTYPE,NTYPE,CTYPE,LTYPE>::
  SetVertexCoord             ///< Set cube vertex coordinates.
  (const CTYPE2 * vertex0_coord, const std::vector<LTYPE2> & edge_length)
  {
    SetVertexCoord(vertex0_coord, IJK::vector2pointer(edge_length));
  }


  // **************************************************
  // TEMPLATE CLASS UNIT_CUBE MEMBER FUNCTIONS
  // **************************************************

  /// Constructor.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  UNIT_CUBE<DTYPE,NTYPE,CTYPE>::UNIT_CUBE(const DTYPE dimension):
    CUBE<DTYPE,NTYPE,CTYPE,CTYPE>(dimension)
  {}

  /// Constructor.
  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  UNIT_CUBE<DTYPE,NTYPE,CTYPE>::UNIT_CUBE():
    CUBE<DTYPE,NTYPE,CTYPE,CTYPE>()
  {}


  // **************************************************
  // TEMPLATE CLASS CUBE_CENTER_TO_VERTEX_DIRECTIONS MEMBER FUNCTIONS
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  CUBE_CENTER_TO_VERTEX_DIRECTIONS<DTYPE,NTYPE,CTYPE>::
  CUBE_CENTER_TO_VERTEX_DIRECTIONS()
  {
    Init(0);
  }


  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  CUBE_CENTER_TO_VERTEX_DIRECTIONS<DTYPE,NTYPE,CTYPE>::
  CUBE_CENTER_TO_VERTEX_DIRECTIONS(const DTYPE dimension)
  {
    Init(dimension);
  }


  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  void CUBE_CENTER_TO_VERTEX_DIRECTIONS<DTYPE,NTYPE,CTYPE>::
  Init(const DTYPE dimension)
  {
    SetDimension(dimension);
  }


  template <typename DTYPE, typename NTYPE, typename CTYPE> 
  void CUBE_CENTER_TO_VERTEX_DIRECTIONS<DTYPE,NTYPE,CTYPE>::
  SetDimension(const DTYPE dimension)
  {
    CUBEV<DTYPE,NTYPE,CTYPE>::SetDimension(dimension);

    compute_cube_center_to_vertex_directions(dimension, this->vertex_coord);
  }

}

#endif
