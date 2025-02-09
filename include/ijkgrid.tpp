/*!
 *  @file ijkgrid.tpp
 *  @brief ijk templates defining regular grid classes and functions.
 *  - Version 0.4.1
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2008-2024 Rephael Wenger

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

#ifndef _IJKGRID_TPP_
#define _IJKGRID_TPP_

#include <algorithm>
#include <limits>
#include <vector>

#include "ijk.tpp"
#include "ijkgrid_func.tpp"
#include "ijkcube.tpp"


namespace IJK {

  // *****************************************************************
  // TYPE DEFINITIONS
  // *****************************************************************

  /// Default type for grid size
  typedef long GRID_SIZE_TYPE;

  // *****************************************************************
  // TEMPLATE CLASS GRID
  // *****************************************************************

  /// Base grid class.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  class GRID {

  protected:
    DTYPE dimension;     ///< grid dimension
    ATYPE * axis_size;   ///< axis_size[i] = # grid points along axis i
    NTYPE num_vertices;  ///< number of grid vertices

    void Init            /// Initialize grid.
    (const DTYPE dimension, const ATYPE * axis_size);
    void FreeAll();      ///< Free all allocated memory.


  public:
    typedef DTYPE DIMENSION_TYPE;         ///< Dimension type.
    typedef ATYPE AXIS_SIZE_TYPE;         ///< Axis size type.
    typedef VTYPE VERTEX_INDEX_TYPE;      ///< Vertex index type.
    typedef NTYPE NUMBER_TYPE;            ///< Number type.

  public:
    // Constructors, destructors, assignment.
    GRID(const DTYPE dimension, const ATYPE * axis_size);
    GRID();
    ~GRID();
    GRID(const GRID & grid);

    // Copy function.
    template <typename GTYPE2>
    const GRID & Copy(const GTYPE2 & right);
    template <typename GTYPE2>
    const GRID & operator = (const GTYPE2 & right)
    { return(Copy(right)); }

    // set functions
    template <typename DTYPE2, typename ATYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // get functions
    inline DTYPE Dimension() const         /// Dimension.
    { return dimension; }
    inline const ATYPE * AxisSize() const  /// Array axis size.
    { return axis_size; }
    inline ATYPE AxisSize(const DTYPE i) const /// Axis size[i].
    { return axis_size[i]; }
    inline NTYPE NumVertices() const       /// Number of grid vertices. 
    { return num_vertices; }

    /*!
     *  @brief Return true if cube facet ifacet
     *    is the index of a lower/left cube facet.
     */
    template <typename ITYPEF>
    bool IsLowerFacet(const ITYPEF ifacet) const
    { return (ifacet < Dimension()); }

    /*!
     *  @brief Return cube facet index.
     *  - Cube facet index is in range [0,(num_cube_facets)-1].
     *  @param orth_dir Direction orthogonal to facet.
     *  @param flag_below True if facet is below cube.
     */
    template <typename _DIR_TYPE>
    _DIR_TYPE CubeFacetIndex
    (const _DIR_TYPE orth_dir, const bool flag_below) const
    { return (orth_dir + (1-int(flag_below))*Dimension()); }

    // compute functions
    NTYPE ComputeNumCubes() const;
    NTYPE ComputeNumEdges() const;
    NTYPE ComputeNumInteriorVertices() const;
    NTYPE ComputeNumInteriorCubes() const;
    NTYPE ComputeNumBoundaryCubes() const;
    template <typename DTYPE2>
    NTYPE ComputeNumVerticesInFacet(const DTYPE2 orth_dir) const;
    template <typename DTYPE2, typename WTYPE>
    NTYPE ComputeNumVerticesInFacet
    (const DTYPE2 orth_dir, const WTYPE boundary_width) const;
    template <typename DTYPE2>
    NTYPE ComputeNumCubesInFacet(const DTYPE2 orth_dir) const;
    template <typename DTYPE2>
    NTYPE ComputeNumFacetsInGridFacet(const DTYPE2 orth_dir) const;
    template <typename GTYPE>
    VTYPE ComputeVertexIndex(const GTYPE * coord) const;
    template <typename GTYPE>
    VTYPE ComputeVertexIndex(const std::vector<GTYPE> coord) const;
    template <typename GTYPE>
    void ComputeCoord(const VTYPE iv, GTYPE * coord) const;
    template <typename GTYPE>
    void ComputeCubeCenterCoord(const VTYPE iv, GTYPE * coord) const;
    template <typename BTYPE>
    void ComputeBoundaryBits(const VTYPE iv, BTYPE & boundary_bits) const;
    template <typename BTYPE>
    void ComputeBoundaryCubeBits
    (const VTYPE icube, BTYPE & boundary_bits) const;
    template <typename CTYPE, typename DIST_TYPE>
    void ComputeCubeDistanceToGridBoundary
    (const CTYPE cube_coord[], DIST_TYPE & distance) const;
    template <typename PTYPE, typename ATYPE2>
    void ComputeSubsampledAxisSizes
    (const PTYPE subsample_period, ATYPE2 subsampled_axis_size[]) const;

    // compare
    template <typename DTYPE2, typename ATYPE2>
    bool CompareSize(const DTYPE2 dimension, const ATYPE2 * axis_size) const;
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    bool CompareSize(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid) const;

    /// Return true if grid contains specified point.
    template <typename CTYPE>
    bool ContainsPoint(const CTYPE * coord) const;
    template <typename CTYPE>
    bool ContainsPoint(const std::vector<CTYPE> & coord) const;

    /// Return true if cube contains specified point.
    template <typename CTYPE, typename VTYPE2>
    bool CubeContainsPoint
    (const VTYPE2 icube, const CTYPE * coord) const;
    template <typename CTYPE, typename VTYPE2>
    bool CubeContainsPoint
    (const VTYPE2 icube, const std::vector<CTYPE> & coord) const;

    /// Return true if grid contains specified region.
    /// @pre Box.MinCoord(d) <= Box.MaxCoord(d) for every d < Box.Dimension().
    template <typename BOX_TYPE>
    bool ContainsRegion(const BOX_TYPE & region) const;

    /// @brief Return true if cube is on grid boundary.
    inline bool IsCubeOnGridBoundary(const VTYPE icube) const
    { return is_cube_on_grid_boundary
        (icube, Dimension(), AxisSize()); }

    /// Return true if cube facet is on grid boundary.
    template <typename CTYPE, typename DTYPE2>
    bool IsCubeFacetOnGridBoundary
    (const CTYPE cube_index, const DTYPE2 facet_orth_dir, 
     const bool facet_side) const;


    // Print routines (mainly for debugging.)

    /// Print vertex coordinates (mainly for debugging).
    template <typename OSTREAM_TYPE, typename VTYPE0>
    void PrintCoord(OSTREAM_TYPE & out, const VTYPE0 iv) const;

    /// Print vertex coordinates (mainly for debugging).
    template <typename OSTREAM_TYPE, typename VTYPE0>
    void PrintCoord(OSTREAM_TYPE & out, const char * prefix, 
                    const VTYPE0 iv, const char * suffix) const;

    /// Print vertex index and coordinates (mainly for debugging).
    template <typename OSTREAM_TYPE, typename VTYPE0>
    void PrintIndexAndCoord(OSTREAM_TYPE & out, const VTYPE0 iv) const;

    /// Print vertex index and coordinates (mainly for debugging).
    template <typename OSTREAM_TYPE, typename VTYPE0>
    void PrintIndexAndCoord(OSTREAM_TYPE & out, const char * prefix,
                            const VTYPE0 iv, const char * suffix) const;

    /// Print two vertex indices and coordinates (mainly for debugging).
    template <typename OSTREAM_TYPE, typename VTYPE0, typename VTYPE1>
    void PrintIndexAndCoord
    (OSTREAM_TYPE & out, const char * text0,
     const VTYPE0 iv0, const char * text1, 
     const VTYPE1 iv1, const char * text2) const;

    /// Print three vertex indices and coordinates (mainly for debugging).
    template <typename OSTREAM_TYPE, 
              typename VTYPE0, typename VTYPE1, typename VTYPE2>
    void PrintIndexAndCoord
    (OSTREAM_TYPE & out, const char * text0,
     const VTYPE0 iv0, const char * text1, 
     const VTYPE1 iv1, const char * text2,
     const VTYPE2 iv2, const char * text3) const;

    // Check functions

    /// Return true if grid *this has same dimension as grid. 
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    bool CheckDimension
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid,
     const char * grid_label1, const char * grid_label2,
     IJK::ERROR & error) const;

    /// Return true if grid *this has given dimension and axis_size[].
    bool Check(const DTYPE dimension, const ATYPE * axis_size,
               IJK::ERROR & error) const;

    /// Return true if grid *this has same dimension and axis_size[] as grid. 
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    bool Check(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid,
               IJK::ERROR & error) const;

    /// Return true if current grid size (grid1) matches grid2 size.
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    bool Check(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2,
               const char * grid_label1, const char * grid_label2,
               IJK::ERROR & error) const;

    /// Return true if coord is the coordinate of a grid vertex.
    template <typename GTYPE>
    bool CheckCoord(const GTYPE * coord, IJK::ERROR & error) const;
    template <typename GTYPE>
    bool CheckCoord(const std::vector<GTYPE> & coord, 
                    IJK::ERROR & error) const;

    /// Return true if coord is the coordinate of a grid cube.
    template <typename GTYPE>
    bool CheckCubeCoord(const GTYPE * coord, IJK::ERROR & error) const;
    template <typename GTYPE>
    bool CheckCubeCoord(const std::vector<GTYPE> & coord, 
                    IJK::ERROR & error) const;

    /*!
     *  @brief Return true if vertex_index is the index of a grid vertex.
     *  - Set error and return false if vertex_index is not the index 
     *    of a grid vertex.
     */
    template <typename ITYPE>
    bool CheckVertexIndex(const ITYPE vertex_index,
                          IJK::ERROR & error) const;

    /*!
     *  @brief Return true if cube_index is the index of a grid cube.
     *  - Set error and return false if cube_index is not the index 
     *    of a grid cube.
     */
    template <typename ITYPE>
    bool CheckCubeIndex(const ITYPE cube_index,
                        IJK::ERROR & error) const;

    /// Check that grid contains specified region.
    template <typename VTYPE2, typename ATYPE2>
    bool CheckContainsRegion
    (const VTYPE2 region_v0, const ATYPE2 * region_axis_size,
     IJK::ERROR & error) const;
  };


  // **************************************************
  // TEMPLATE CLASS GRID_PLUS
  // **************************************************

  /*!
   *  @brief Inherits class GRID and adds other indexes and operators 
   *    for fast accessing of grid vertices.
   *  @tparam DTYPE  Dimension data type.
   *  @tparam ATYPE  Axis size type.
   *  @tparam VTYPE  Vertex index type.
   *  @tparam NTYPE  Number type.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  class GRID_PLUS:public GRID<DTYPE,ATYPE,VTYPE,NTYPE> {

    /// @nosubgrouping
    
  protected:

    /// iv+axis_increment[d] is vertex next to iv
    VTYPE * axis_increment;    

    /// iv0+cube_vertex_increment[k] k'th vertex of cube with primary vertex iv0
    VTYPE * cube_vertex_increment;

    /*!
     *  \brief Increment for computing facet vertices.
     *  iv0+facet_vertex_increment[k+num_facet_vertices*d] = 
     *    k'th vertex of facet orthogonal to d with primary vertex iv0
     */
    VTYPE * facet_vertex_increment;

    /*!
     *  \brief Orientation of facet 0, (facet in hyperplane x_{dim-1} = 0.
     *  - All other facets are oriented to match facet0 orientation.
     *  - If grid dimension is 3, then true is counter-clockwise orientation.
     */
    bool facet0_orientation;

    /*!
     *  \brief Increment for computing ridge vertices.
     *  iv0+ridge_vertex_increment[k+nv*d0+nv*nv*d1] = 
     *    k'th vertex of ridge orthogonal to d0 and d1 with primary vertex iv0
     *    nv = number of ridge vertices.
     */
    VTYPE * ridge_vertex_increment;

    /// unit_cube_coord[dimension*k+j] = j'th coordinate of k'th vertex of unit cube
    NTYPE * unit_cube_coord;

    NTYPE num_cube_vertices;        ///< Number of cube vertices.

    // *** SHOULD RENAME AS num_cube_facet_vertices. ***
    NTYPE num_facet_vertices;       ///< Number of cube facet vertices.
    
    NTYPE num_cube_facets;          ///< Number of cube facets.
    NTYPE num_cube_edges;           ///< Number of cube edges.

    /// Number of cube ridge vertices.
    NTYPE num_cube_ridge_vertices;

    void ZeroLocal();
    void InitLocal();
    void FreeLocal();
    void Create();            ///< Allocate and set data in GRID_PLUS.

    // Not defined.
    template <typename GRID2>
    void Copy(const GRID2 & right);

  public:
    /// Constructors.
    GRID_PLUS(const DTYPE dimension, const ATYPE * axis_size);
    GRID_PLUS();
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    GRID_PLUS(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);
    GRID_PLUS(const GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> & grid2);

    ~GRID_PLUS();                       ///< Destructor

    // set functions
    template <typename DTYPE2, typename ATYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // get functions
    inline NTYPE NumCubeVertices() const       /// Return number of cube vertices.
    { return num_cube_vertices; };
    inline NTYPE NumCubeFacets() const         /// Return number of cube facets.
    { return num_cube_facets; };
    inline NTYPE NumCubeEdges() const          /// Return number of cube edges.
    { return num_cube_edges; };
    inline NTYPE NumDiagonals() const          /// Return number of cube diagonals.
    { return num_facet_vertices; };
    const VTYPE * AxisIncrement() const /// Return array axis_increment[]
    { return(axis_increment); }
    inline const VTYPE AxisIncrement(const DTYPE d) const /// Return axis_increment[d]
    { return axis_increment[d]; }
    inline const VTYPE * CubeVertexIncrement() const /// Return cube_vertex_increment[]
    { return cube_vertex_increment; }
    inline const VTYPE CubeVertexIncrement     /// Return cube_vertex_increment[k]
    (const VTYPE k) const              
    { return cube_vertex_increment[k]; }

    /// @brief Return the direction of a grid edge.
    /// @pre endpoint0 and endpoint1 are the endpoints of a grid edge.
    template <typename VTYPE2, typename VTYPE3>
    DTYPE EdgeDirection(const VTYPE2 endpoint0, const VTYPE3 endpoint1) const;

    /// @brief Return number of cube facet vertices. *** DEPRECATED ***
    inline NTYPE NumFacetVertices() const
    { return(num_facet_vertices); };

    /// @brief Return number of cube facet vertices.
    inline NTYPE NumCubeFacetVertices() const
    { return(num_facet_vertices); };

    /// @brief Return number of cube ridge vertices.
    inline NTYPE NumCubeRidgeVertices() const
    { return(num_cube_ridge_vertices); };

    /*!
     *  @brief Return vertex increment of k'th vertex of facet ifacet.
     *  @param ifacet Facet index. In range [0..(Nf-1)]
     *         where Nf is the number of cube facets.
     *         - If \a ifacet < (Nf/2), then facet \a ifacet is the
     *           lower facet orthogonal to direction \a ifacet.
     *         - If \a ifacet >= (Nf/2), then facet \a ifacet is the
     *           upper facet orthogonal to direction (\a ifacet/2).
     */
    inline const VTYPE FacetVertexIncrement
    (const DTYPE ifacet, const VTYPE k) const              
    { return facet_vertex_increment[k+ifacet*num_facet_vertices]; }

    /*!
     *  @brief Return pointer to array of vertex increments of facet ifacet.
     *  - Returns facet_vertex_increment[] where facet_vertex_increment[k]
     *    is increment for k'th vertex of facet ifacet.
     *  @param ifacet Facet index. In range [0..(Nf-1)]
     *         where Nf is the number of cube facets.
     *         - If \a ifacet < (Nf/2), then facet \a ifacet is the
     *           lower facet orthogonal to direction \a ifacet.
     *         - If \a ifacet >= (Nf/2), then facet \a ifacet is the
     *           upper facet orthogonal to direction (\a ifacet/2).
     */
    inline const VTYPE * FacetVertexIncrement(const DTYPE ifacet) const    
    { return (facet_vertex_increment+ifacet*num_facet_vertices); }

    /// @brief Return facet 0 orientation.
    inline const bool FacetZeroOrientation() const
    { return facet0_orientation; }

    /// @brief Return vertex increment of k'th vertex of ridge
    ///   orthogonal to orth_dir0 and orth_dir1.
    const VTYPE RidgeVertexIncrement    /// Return ridge_vertex_increment[k]
    (const DTYPE orth_dir0, const DTYPE orth_dir1, const VTYPE k) const
    { 
      const DTYPE dimension = this->Dimension();
      NTYPE j = 
        k+this->num_cube_ridge_vertices*(orth_dir0+dimension*orth_dir1);
      return(ridge_vertex_increment[j]);
    }

    // *** DEPRECATED. REPLACE BY class UNIT_CUBE. ***
    /// Return pointer to coordinates of k'th cube vertex
    const NTYPE * UnitCubeCoord(const NTYPE k) const
    { return(unit_cube_coord+this->Dimension()*k); }

    // *** DEPRECATED. REPLACE BY class UNIT_CUBE. ***
    /// Return j'th coordinate of k'th vertex
    const NTYPE UnitCubeCoord         
    (const NTYPE k, const NTYPE j) const
    { return(unit_cube_coord[this->Dimension()*k+j]); }

    /// \brief Return next vertex in direction d.
    /// @pre iv is not the last vertex in direction d.
    VTYPE NextVertex(const VTYPE iv, const DTYPE d) const  
    { return(iv+axis_increment[d]); }

    /// \brief Return previous vertex in direction d.
    /// @pre iv is not the first vertex in direction d.
    VTYPE PrevVertex(const VTYPE iv, const DTYPE d) const  
    { return(iv-axis_increment[d]); }

    /// \brief Return adjacent vertex in direction d.
    /// @param side = 0 (or false) or 1 (or true).
    template <typename STYPE>
    VTYPE AdjacentVertex(const VTYPE iv, const DTYPE d, const STYPE side) const
    {
      if (DTYPE(side) == 0) { return(PrevVertex(iv, d)); }
      else { return(NextVertex(iv, d)); };
    }

    /*!
     *  \brief Return k'th cube vertex.
     *  - k'th cube vertex is diagonally opposite vertex (num_cube_vertices-1-k).
     *  @param iv0 is a primary cube vertex.
     *  @param k k'th cube vertex.
     *  @pre k is less than the number of unit cube vertices.
     */
    VTYPE CubeVertex(const VTYPE iv0, const int k) const  
    { return(iv0+cube_vertex_increment[k]); }

    /*!
     *  \brief Return k'th facet vertex.
     *  @param iv0 is a primary cube vertex.
     *  @param ifacet Facet index. In range [0..(Nf-1)]
     *         where Nf is the number of cube facets.
     *         - If \a ifacet < (Nf/2), then facet \a ifacet is the
     *           lower facet orthogonal to direction \a ifacet.
     *         - If \a ifacet >= (Nf/2), then facet \a ifacet is the
     *           upper facet orthogonal to direction (\a ifacet%Dimension()).
     *  @param k \a k'th facet vertex.
     *  @pre \a ifacet is less than the number of cube facets.
     *  @pre \a k is less than the number of facet vertices.
     */
    VTYPE FacetVertex(const VTYPE iv0, const DTYPE ifacet, const int k) const
    { return(iv0+facet_vertex_increment[k+ifacet*num_facet_vertices]); }

    /// \brief Return k'th ridge vertex.
    VTYPE RidgeVertex(const VTYPE iv0, 
                      const DTYPE orth_dir0, const DTYPE orth_dir1, 
                      const int k) const
    { return(iv0+RidgeVertexIncrement(orth_dir0, orth_dir1, k)); }

    /// \brief Compute vertex neighbors.
    template <typename VTYPE0, typename VTYPE1, typename DIST_TYPE>
    void GetVertexNeighbors
    (const VTYPE0 iv0, const DIST_TYPE distance, std::vector<VTYPE1> & vlist)
      const;

    /*!
     *  @brief Return index of facet shared by icube0 and icube1.
     *  @pre icube0 and icube1 share some facet.
     *  - Index is some integer in range [0..(2*dimension-1)]
     *    representing cube facet.
     */
    template <typename ITYPE>
    NTYPE SharedFacet(const ITYPE icube0, const ITYPE icube1) const;

    /*!
     *  @brief Return direction orthogonal to facet shared by icube0 and icube1.
     *  @pre icube0 and icube1 share some facet.
     */
    template <typename ITYPE>
    NTYPE SharedFacetOrthDir(const ITYPE icube0, const ITYPE icube1) const;    

    /*!
     *  @brief Return d'th coordinate of vertex \a iv.
     *  @param iv Vertex index.
     *  @param d Axis index.
     */
    VTYPE CoordD(const VTYPE iv, const DTYPE d) const;
  };

  
  // **************************************************
  // TEMPLATE CLASS GRID_NEIGHBORS
  // **************************************************

  /*!
   *  @brief Inherits class GRID_PLUS and adds other indexes and operators
   *    for fast accessing neighbors of grid vertices.
   *  @tparam DTYPE  Dimension data type.
   *  @tparam ATYPE  Axis size type.
   *  @tparam VTYPE  Vertex index type.
   *  @tparam DIFFTYPE  Index difference type.  Must be signed.
   *  @tparam NTYPE  Number type.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename DIFFTYPE, 
            typename NTYPE> 
  class GRID_NEIGHBORS:public GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> {

  protected:

    /// Number of neighbors of a vertex in cubes containing the vertex.
    NTYPE num_vertex_neighborsC; 

    /// iv + vertex_neighborC[k] = k'th neighbor of vertex iv in cubes containing iv.
    DIFFTYPE * vertex_neighborC;

    /// Number of vertices which share an edge with a vertex.
    NTYPE num_vertex_neighborsE;

    /// iv + vertex_neighborE[k] = k'th neighbor sharing an edge with vertex iv.
    DIFFTYPE * vertex_neighborE;

    /// cube_index + cube_neighborE[k] = 
    ///   k'th neighbor sharing an edgex with cube cube_index.
    DIFFTYPE * cube_neighborE;

    /// cube_index + cube_neighborV[k] = 
    ///   k'th neighbor sharing a vertex with cube cube_index.
    DIFFTYPE * cube_neighborV;

    /// \brief Number of vertices in cubes containing a facet, 
    ///   not including facet vertices.
    NTYPE num_facet_neighborsC;

    /// \brief iv + facet_neighborC[d*num_facet_neighborsC+k] = 
    ///   k'th neighbor of facet d, primary vertex iv
    DIFFTYPE * facet_neighborC;

    /// \brief Number of vertices in 2-faces containing an edge,
    ///   not including edge vertices.
    NTYPE num_edge_neighborsF2;

    /// \brief iv + edge_neighborC[d*num_edge_neighborsF2 + K] =
    ///   k'th neighbor of edge d, primary vertex iv
    DIFFTYPE * edge_neighborF2;

    void ZeroLocal();
    void InitLocal();
    void FreeLocal();
    void CreateLocal();       ///< Allocate and set data in GRID_NEIGHBORS.

    // Not defined.
    template <typename GRID2>
    void Copy(const GRID2 & right);

  public:
    /// Constructors.
    GRID_NEIGHBORS(const DTYPE dimension, const ATYPE * axis_size);
    GRID_NEIGHBORS();
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    GRID_NEIGHBORS(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);
    GRID_NEIGHBORS
    (const GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE, NTYPE> & grid2);

    ~GRID_NEIGHBORS();        ///< Desctructor.

    // set functions

    /// Set dimensions and axis size.
    template <typename DTYPE2, typename ATYPE2>
    void SetSize       
    (const DTYPE2 dimension, const ATYPE2 * axis_size);

    /// Set dimensions and axis size.
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize       
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    // get functions

    /// \brief Return number of neighbors of a vertex.
    /// \details All vertices in a cube containing the vertex 
    ///    are counted as neighbors.
    /// The vertex itself is not counted in this number.
    NTYPE NumVertexNeighborsC() const   
    { return(num_vertex_neighborsC); };

    /// \brief Return number of vertices which share an edge with a vertex.
    /// \details The vertex itself is not counted in this number.
    NTYPE NumVertexNeighborsE() const   
    { return(num_vertex_neighborsE); };

    /// \brief Return number of cubes which share a facet with a cube.
    NTYPE NumCubeNeighborsF() const   
    { return(this->NumCubeFacets()); };

    /// \brief Return number of cubes which share an edge with a cube.
    NTYPE NumCubeNeighborsE() const   
    { return(this->NumCubeEdges()); };

    /// \brief Return number of cubes which share a vertex with a cube.
    NTYPE NumCubeNeighborsV() const   
    { return(this->NumCubeVertices()); };

    /// \brief Return number of neighbors of a facet.
    /// \details All vertices in a cube containing the facet are counted as neighbors,
    ///          not including vertices lying on the facet.
    NTYPE NumFacetNeighborsC() const   
    { return(num_facet_neighborsC); };

    /// \brief Return number of neighbors of an edge.
    /// \details All vertices in a 2-face containing the fedge are counted as neighbors,
    ///          not including the edge endpoints.
    NTYPE NumEdgeNeighborsF2() const   
    { return(num_edge_neighborsF2); };

    /// \brief Return k'th neighbor of vertex \a iv in cubes containing \a iv.
    /// \details A grid vertex is not a neighbor of itself.
    /// @pre Vertex \a iv must be an internal grid vertex, 
    ///    i.e., not on the grid boundary.
    VTYPE VertexNeighborC               
    (const VTYPE iv, const NTYPE k) const
    { return(iv+vertex_neighborC[k]); };

    /// \brief Return k'th vertex which shares an edge with the vertex.
    /// \details A grid vertex is not a neighbor of itself.
    /// @pre Vertex \a iv must be an internal grid vertex, 
    ///    i.e., not on the grid boundary.
    VTYPE VertexNeighborE
    (const VTYPE iv, const NTYPE k) const
    { return(iv+vertex_neighborE[k]); };

    /// \brief Return k'th cube which shares a facet with cube_index.
    /// @pre Cube \a cube_index must be an internal grid cube, 
    ///    i.e., not on the grid boundary.
    VTYPE CubeNeighborF
    (const VTYPE cube_index, const NTYPE k) const
    { return(this->VertexNeighborE(cube_index, k)); };

    /// \brief Return k'th cube which shares an edge with cube_index.
    /// @pre Cube \a cube_index must be an internal grid cube, 
    ///    i.e., not on the grid boundary.
    VTYPE CubeNeighborE
    (const VTYPE cube_index, const NTYPE k) const
    { return(cube_index+cube_neighborE[k]); };

    /// \brief Return k'th cube which shares a vertex with cube_index.
    /// @pre Cube \a cube_index must be an internal grid cube, 
    ///    i.e., not on the grid boundary.
    VTYPE CubeNeighborV
    (const VTYPE cube_index, const NTYPE k) const
    { return(cube_index+ cube_neighborV[k]); };

    /// \brief Return k'th neighbor of facet.
    /// @param iv  Primary facet vertex (facet vertex with lowest coordinates.)
    /// @param orth_dir  Direction orthogonal to facet.
    /// @param k   Return \a k'th neighbor.
    /// @pre Facet must NOT be contained in the grid boundary.
    VTYPE FacetNeighborC
    (const VTYPE iv, const DTYPE orth_dir, const NTYPE k) const
    { return(iv+facet_neighborC[orth_dir*num_facet_neighborsC+k]); };

    /// \brief Return k'th neighbor of edge.
    /// @param iv  Primary edge vertex (edge endpoint with lowest coordinates.)
    /// @param dir  Direction of edge.
    /// @param k   Return \a k'th neighbor.
    /// @pre Edge must NOT be contained in the grid boundary.
    VTYPE EdgeNeighborF2
    (const VTYPE iv, const DTYPE dir, const NTYPE k) const
    { return(iv+edge_neighborF2[dir*num_edge_neighborsF2+k]); };

  };


  // **************************************************
  // TEMPLATE CLASS GRID_SPACING
  // **************************************************

  /*!
   *  @brief Template to add grid spacing to grid.
   *  @tparam SP_TYPE Spacing type. (Should be float/double.)
   */
  template <typename SP_TYPE, typename GRID_TYPE>
  class GRID_SPACING:public GRID_TYPE
  {
  protected:
    SP_TYPE * spacing;   ///< spacing[i] = grid spacing along axis i.

    void InitLocal();                    ///< Initialize GRID_SPACING.
    void CreateLocal();                  ///< Create local data structure.
    void FreeLocal();                    ///< Free memory.

  public:
    typedef SP_TYPE SPACING_TYPE;


  public:
    GRID_SPACING();                      ///< Constructor.
    template <typename DTYPE2, typename ATYPE2>
    GRID_SPACING(const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    GRID_SPACING(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);
    GRID_SPACING(const GRID_SPACING<SP_TYPE, GRID_TYPE> & grid2);
    ~GRID_SPACING();                     ///< Destructor.

    // Copy function.
    template <typename GTYPE2>
    const GRID_SPACING & Copy(const GTYPE2 & right);
    template <typename GTYPE2>
    const GRID_SPACING & operator = (const GTYPE2 & right)
    { return(Copy(right)); }

    template <typename DTYPE2, typename SP_TYPE2>    ///< Set spacing[d] to \a c.
    void SetSpacing(const DTYPE2 d, const SP_TYPE2 c)
    { spacing[d] = c; };

    template <typename SP_TYPE2>           ///< Set all grid spacing to \a c.
    void SetAllSpacing(const SP_TYPE2 c);

    template <typename SP_TYPE2>           ///< Set spacing[d] to spacing2[d].
    void SetSpacing(const SP_TYPE2 * spacing2);

    /// @brief Set spacing[d] to spacing2[d].
    /// - Version using C++ STL vector spacing2[].
    template <typename SP_TYPE2>
    void SetSpacing(const std::vector<SP_TYPE2> & spacing2)
    { SetSpacing(IJK::vector2pointer(spacing2)); }

    /// Set spacing[d] to (c*spacing2[d]).
    template <typename SCALE_TYPE, typename SP_TYPE2>
    void SetSpacing(const SCALE_TYPE c, const SP_TYPE2 * spacing2);

    /// Set spacing[d] to grid2.Spacing(d).
    template <typename SP_TYPE2, typename GRID_TYPE2>
    void SetSpacing(const GRID_SPACING<SP_TYPE2, GRID_TYPE2> & grid2);

    template <typename DTYPE2, typename ATYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize       /// Set dimensions and axis size.
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    /// Get spacing.
    template <typename DTYPE2>
    SP_TYPE Spacing(const DTYPE2 d) const
    { return(spacing[d]); };

    /// Return true if (Spacing(d) == 1) for all d.
    bool IsUnitSpacing() const;

    /// Get const pointer to spacing.
    const SP_TYPE * SpacingPtrConst() const
    { return(spacing); }

    /// Compute scaled coord.
    template <typename VTYPE2, typename CTYPE>
    void ComputeScaledCoord
    (const VTYPE2 iv, CTYPE * coord) const;

    /// Compute scaled coordinates of cube center.
    template <typename VTYPE2, typename CTYPE>
    void ComputeCubeCenterScaledCoord
    (const VTYPE2 iv, CTYPE * coord) const;
  };


  // **************************************************
  // TEMPLATE CLASS GEOM_GRID
  // **************************************************

  /// @brief Template to add origin and grid spacing to grid.
  /// - GEOM_GRID has both origina and grid spacing.
  template <typename CTYPE, typename GRID_TYPE>
  class GEOM_GRID:public GRID_SPACING<CTYPE,GRID_TYPE>
  {
  protected:

    /// Coordinates of grid origin.
    CTYPE * origin;      

    void InitLocal();                    ///< Initialize GEOM_GRID.
    void CreateLocal();                  ///< Create local data structure.
    void FreeLocal();                    ///< Free memory.

  public:
    GEOM_GRID();                         ///< Constructor.
    template <typename DTYPE2, typename ATYPE2>
    GEOM_GRID(const DTYPE2 dimension, const ATYPE2 * axis_size);
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    GEOM_GRID(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);
    GEOM_GRID(const GEOM_GRID<CTYPE, GRID_TYPE> & grid2);
    ~GEOM_GRID();                     ///< Destructor.

    // Copy function.
    template <typename GTYPE2>
    const GEOM_GRID & Copy(const GTYPE2 & right);
    template <typename GTYPE2>
    const GEOM_GRID & operator = (const GTYPE2 & right)
    { return(Copy(right)); }

    template <typename DTYPE2, typename CTYPE2>    ///< Set origin[d] to \a c.
    void SetOrigin(const DTYPE2 d, const CTYPE2 c)
    { origin[d] = c; };

    template <typename CTYPE2>           ///< Set all origin[d] to \a c.
    void SetAllOrigin(const CTYPE2 c);

    template <typename CTYPE2>           ///< Set origin[d] to origin2[d].
    void SetOrigin(const CTYPE2 * origin2);

    /// @brief Set origin[d] to origin2[d].
    /// - Version using C++ STL vector origin2[].
    template <typename CTYPE2>
    void SetOrigin(const std::vector<CTYPE2> & origin2)
    { SetOrigin(IJK::vector2pointer(origin2)); }

    /// Set origin[d] to grid2.Origin(d).
    template <typename CTYPE2, typename GRID_TYPE2>
    void SetOrigin(const GEOM_GRID<CTYPE2, GRID_TYPE2> & grid2);

    /// Set dimensions and axis size.
    template <typename DTYPE2, typename ATYPE2>
    void SetSize
    (const DTYPE2 dimension, const ATYPE2 * axis_size);

    /// @brief Set dimensions and axis size.
    /// - Version using C++ STL vector axis_size[].
    template <typename DTYPE2, typename ATYPE2>
    void SetSize       
    (const DTYPE2 dimension, const std::vector<ATYPE2> & axis_size)
    { SetSize(dimension, IJK::vector2pointer(axis_size)); }

    /// @brief Set dimensions and axis size from grid2
    template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
              typename NTYPE2>
    void SetSize
    (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2);

    /// Get origin.
    template <typename DTYPE2>
    CTYPE Origin(const DTYPE2 d) const
    { return(origin[d]); };

    /// Get const pointer to spacing.
    const CTYPE * OriginPtrConst() const
    { return(origin); }
  };


  // *****************************************************************
  // FUNCTIONS WITH CLASS GRID ARGUMENT.
  // *****************************************************************

  // *** SHOULD BE MODIFIED TO MEMBER FUNCTIONS.
  
  /// DEPRECATED
  /*!
   *  @brief Compute increment to add to index of current vertex to get
   *    next vertex along each axis
   *  @pre Array increment[] is pre-allocated to size at least \a dimension.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE, typename ITYPE>
  void compute_increment
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, ITYPE * increment)
  {
    compute_increment(grid.Dimension(), grid.AxisSize(), increment);
  }

  /*!
   *  @brief Compute increment to add to vertex 0 to compute vertex i of hypercube.
   *  @param grid Grid.
   *  @param[out] cube_vertex_increment[] Cube vertex increment.
   *    - iv0+cube_vertex_increment[i] is the i'th vertex of the hypercube
   *      with primary vertex iv0.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE, typename ITYPE>
  void compute_cube_vertex_increment
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, ITYPE * cube_vertex_increment)
  {
    const DTYPE dimension = grid.Dimension();
    VTYPE axis_increment[dimension];

    compute_increment(grid, axis_increment);
    compute_cube_vertex_increment(dimension, axis_increment, cube_vertex_increment);
  }


  /// Subsample vertices in subgrid.
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename VTYPE, typename NTYPE>
  void subsample_subgrid_vertices
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid,
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   const PTYPE subsample_period, VTYPE * vlist)
  {
    subsample_subgrid_vertices
      (grid.Dimension(), grid.AxisSize(), subgrid_origin, subgrid_axis_size,
       subsample_period, vlist);
  }


  /// @brief Return number of vertices in specified grid facet
  /// - Specify grid facet by the direction orthogonal to the facet.
  template <typename DTYPE, typename DTYPE2,
            typename ATYPE, typename VTYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, const DTYPE2 orth_dir,
   NTYPE & num_vertices)
  {
    compute_num_vertices_in_grid_facet
      (grid.Dimension(), grid.AxisSize(), orth_dir, num_vertices);
  }


  /// Get vertices in grid facet 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  void get_vertices_in_grid_facet0
  (const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid, VTYPE * vlist)
  {
    get_vertices_in_grid_facet0
      (grid.Dimension(), grid.AxisSize(), vlist);
  }

  
  // *****************************************************************
  // TEMPLATE CLASS GRID MEMBER FUNCTIONS
  // *****************************************************************

  /// Constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::GRID
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    Init(dimension, axis_size);
  }


  /// Default constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::GRID()
  {
    Init(0, NULL);
  }


  /// Destructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::~GRID()
  {
    FreeAll();
  }


  /// @brief Initialize grid.
  /// @param dimension  Dimension of grid.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::Init
  (const DTYPE dimension, const ATYPE * axis_size)
  {
    this->axis_size = NULL;
    this->dimension = 0;
    this->num_vertices = 1;
    if (dimension > 0) 
      { SetSize(dimension, axis_size); };
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    FreeAll();

    if (dimension < 0) {
      IJK::PROCEDURE_ERROR error("Grid::SetSize");
      error.AddMessage("Programming error.  Illegal dimension ",
                       dimension, ".");
      error.AddMessage("Dimension should be non-negative.");
      throw error;
    }

    this->dimension = dimension;
    this->axis_size = NULL;
    if (dimension > 0)
      { this->axis_size = new ATYPE[dimension]; }
    else {
      // allocate axis_size even if dimension equals 0
      this->axis_size = new ATYPE[1];
      this->axis_size[0] = 0;
    }
      
    for (DTYPE d = 0; d < dimension; d++)
      { this->axis_size[d] = axis_size[d]; }

    compute_num_grid_vertices(dimension, axis_size, this->num_vertices);
  }


  /// @brief Set size of \a grid to size of \a grid2.
  /// @param grid2  Grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }


  /// Free all allocated memory in typename GRID.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::FreeAll()
  {
    if (axis_size != NULL) { delete [] axis_size; };
    axis_size = NULL;
    dimension = 0;
    num_vertices = 0;
  }


  /// Copy constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  GRID(const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & grid)
  {
    Init(grid.Dimension(), grid.AxisSize());
  }


  /// Copy.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename GTYPE2>
  const GRID<DTYPE,ATYPE,VTYPE,NTYPE> & 
  GRID<DTYPE,ATYPE,VTYPE,NTYPE>::Copy(const GTYPE2 & right)
  {
    if (&right != this) {         // avoid self-assignment
      SetSize(right.Dimension(), right.AxisSize());
    }
    return *this;
  }


  /// Compute and return number of grid cubes.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumCubes() const
  {
    NTYPE num_grid_cubes;
    compute_num_grid_cubes(Dimension(), AxisSize(), num_grid_cubes);
    return(num_grid_cubes);
  }


  /// Compute and return number of grid edges
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumEdges() const
  {
    NTYPE num_grid_edges;
    compute_num_grid_edges(Dimension(), AxisSize(), num_grid_edges);
    return(num_grid_edges);
  }


  /// Compute and return number of vertices in grid interior.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumInteriorVertices() const
  {
    NTYPE num_interior_vertices;
    compute_num_interior_grid_vertices
      (Dimension(), AxisSize(), num_interior_vertices);
    return(num_interior_vertices);
  }


  /// Compute and return number of cubes in grid interior.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumInteriorCubes() const
  {
    NTYPE num_interior_cubes;
    compute_num_interior_grid_cubes
      (Dimension(), AxisSize(), num_interior_cubes);
    return(num_interior_cubes);
  }


  /// Compute and return number of cubes in grid boundary.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::ComputeNumBoundaryCubes() const
  {
    NTYPE num_boundary_cubes;
    compute_num_boundary_grid_cubes
      (Dimension(), AxisSize(), num_boundary_cubes);
    return(num_boundary_cubes);
  }


  /// Compute and return number of vertices in facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeNumVerticesInFacet(const DTYPE2 orth_dir) const
  {
    NTYPE num_vertices_in_facet;
    compute_num_vertices_in_grid_facet
      (Dimension(), AxisSize(), orth_dir, num_vertices_in_facet);
    return(num_vertices_in_facet);
  }


  /// Compute and return number of vertices in facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename WTYPE>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeNumVerticesInFacet
  (const DTYPE2 orth_dir, const WTYPE boundary_width) const
  {
    NTYPE num_vertices_in_facet;
    compute_num_vertices_in_grid_facet
      (Dimension(), AxisSize(), orth_dir, boundary_width,
       num_vertices_in_facet);
    return(num_vertices_in_facet);
  }


  /// Compute and return number of facets in grid facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeNumFacetsInGridFacet(const DTYPE2 orth_dir) const
  {
    NTYPE num_facets;
    compute_num_facets_in_grid_facet
      (Dimension(), AxisSize(), orth_dir, num_facets);
    return(num_facets);
  }


  /// Compute and return number of cubes in facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2>
  NTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeNumCubesInFacet(const DTYPE2 orth_dir) const
  {
    NTYPE num_cubes_in_facet;
    compute_num_cubes_in_grid_facet
      (Dimension(), AxisSize(), orth_dir, num_cubes_in_facet);
    return(num_cubes_in_facet);
  }


  /// @brief Compute index of vertex with given coordinates.
  /// @param coord  Array: <em>coord[d]</em> = d'th coordinate of vertex.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  VTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeVertexIndex(const GTYPE * coord) const
  {
    return(compute_vertex_index<VTYPE>(coord, Dimension(), AxisSize()));
  }


  /// @brief Compute index of vertex with given coordinates.
  /// @param coord  Array: <em>coord[d]</em> = d'th coordinate of vertex.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  VTYPE GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeVertexIndex(const std::vector<GTYPE> coord) const
  {
    return(compute_vertex_index<VTYPE>(&(coord[0]), Dimension(), AxisSize()));
  }


  /*!
   *  @brief Compute coordinates of given vertex.
   *  @param iv  Vertex index.
   *  @param[out] coord  Array: <em>coord[d]</em> = d'th coordinate of vertex.
   *  @pre Array <em>coord[]</em> is preallocated to length
   *       at least \a Dimension().
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeCoord(const VTYPE iv, GTYPE * coord) const
  {
    compute_coord(iv, Dimension(), AxisSize(), coord);
  }


  /*!
   *  @brief Compute coordinates of given cube center.
   *  - Cube center coord is (vertex coord) + (0.5,0.5,...,0.5).
   *  @param iv  Index of primary vertex of cube.
   *  @param[out] coord  Array: <em>coord[d]</em> = d'th coordinate of vertex.
   *  @pre Array <em>coord[]</em> is preallocated to length
   *       at least \a Dimension().
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeCubeCenterCoord(const VTYPE iv, GTYPE * coord) const
  {
    ComputeCoord(iv, coord);
    for (DTYPE d = 0; d < Dimension(); d++)
      { coord[d] += 0.5; }
  }


  /*!
   *  @brief Compute bits identifying which boundary contains vertex \a iv.
   *  @param iv  Vertex index.
   *  @param [out] boundary_bits Bits flagging boundaries containing
   *          vertex \a iv.
   *        If bit \a 2d is true, then <em>d</em>'th coordinate of
   *          vertex \a iv is zero.
   *        If bit <em>(2d+1)</em> is true, then <em>d</em>'th coordinate
   *          of vertex \a iv equals <em>axis_size[d]-1</em>.
   *  @pre \li Variable \a boundary_bits has at least
   *           <em>(2*dimension)</em> bits.
   *  @pre \li AxisSize(d) > 0 for all \a d = 0,..., \a dimension-1.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename BTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeBoundaryBits(const VTYPE iv, BTYPE & boundary_bits) const
  {
    compute_boundary_bits(iv, Dimension(), AxisSize(), boundary_bits);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename PTYPE, typename ATYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeSubsampledAxisSizes
  (const PTYPE subsample_period, ATYPE2 subsampled_axis_size[]) const
  {
    compute_subsample_axis_sizes
      (Dimension(), AxisSize(), subsample_period, subsampled_axis_size);
  }


  /*!
   *  @brief Compute bits identifying which boundary contains cube \a icube.
   *  @param icube  Cube index.
   *  @param [out] boundary_bits Bits flagging boundaries containing cube \a icube.
   *        If bit \a 2d is true, then <em>d</em>'th coordinate
   *               of cube \a icube is zero.
   *        If bit <em>(2d+1)</em> is true, then <em>d</em>'th coordinate
   *               of cube \a iv equals <em>axis_size[d]-2</em>.
   *  @pre \li Variable \a boundary_bits has at least <em>(2*dimension)</em> bits.
   *  @pre \li AxisSize(d) > 0 for all \a d = 0,..., \a dimension-1.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename BTYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeBoundaryCubeBits(const VTYPE icube, BTYPE & boundary_bits) const
  {
    compute_boundary_cube_bits(icube, Dimension(), AxisSize(), boundary_bits);
  }


  /// @brief Compute distance of cube to grid boundary.
  /// @pre Cube is contained in grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename DIST_TYPE>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ComputeCubeDistanceToGridBoundary
  (const CTYPE * cube_coord, DIST_TYPE & distance) const
  {
    if (Dimension() < 1) {
      distance = 0;
      return;
    }

    distance = AxisSize(0);

    for (int d = 0; d < Dimension(); d++) {
      if (cube_coord[d] < distance) { distance = cube_coord[d]; }

      DIST_TYPE d2 = AxisSize(d)-2-cube_coord[d];
      if (d2 < distance) { distance = d2; }
    }
  }

  /// @brief Return true if grid dimension and axis size match parameters.
  /// @param dimension  Dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CompareSize(const DTYPE2 dimension, const ATYPE2 * axis_size) const
  {
    if (dimension != this->Dimension()) { return(false); };
    for (int d = 0; d < dimension; d++) {
      if (axis_size[d] != this->AxisSize(d)) { return(false); };
    }

    return(true);
  }

  /// @brief Return true if dimensions and axis size match dimensions and axis size of \a grid2.
  /// @param grid2  Grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CompareSize(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2) const
  {
    return(this->CompareSize(grid2.Dimension(), grid2.AxisSize()));
  }


  /// Return true if grid contains specified point.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ContainsPoint(const CTYPE * coord) const
  {
    const DTYPE dimension = this->Dimension();

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] < 0 || coord[d]+1 > this->AxisSize(d))
        { return(false); }
    }

    return(true);
  }


  /// Return true if grid contains specified point.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ContainsPoint(const std::vector<CTYPE> & coord) const
  {
    return(this->ContainsPoint(&(coord[0])));
  }


  /// Return true if cube contains specified point.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename VTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CubeContainsPoint(const VTYPE2 icube, const CTYPE * coord) const
  {
    const DTYPE dimension = this->Dimension();
    IJK::ARRAY<VTYPE> cube_coord(dimension);

    this->ComputeCoord(icube, cube_coord.Ptr());

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] < cube_coord[d] || coord[d] > cube_coord[d]+1)
        { return(false); }
    }

    return(true);
  }


  /// Return true if cube contains specified point.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename VTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CubeContainsPoint
  (const VTYPE2 icube, const std::vector<CTYPE> & coord) const
  {
    return(this->CubeContainsPoint(icube, &(coord[0])));
  }


  /// Return true if grid contains specified region.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename BOX_TYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  ContainsRegion(const BOX_TYPE & region) const
  {
    DTYPE dimension = this->Dimension();
    if (region.Dimension() < dimension)
      { dimension = region.Dimension(); };

    for (DTYPE d = 0; d < dimension; d++) {
      if (region.MinCoord(d) < 0 ||
          this->AxisSize(d) <= region.MaxCoord(d))
        { return(false); }
    }

    return(true);
  }

  
  // Return true if cube facet is on grid boundary.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename CTYPE, typename DTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  IsCubeFacetOnGridBoundary
  (const CTYPE cube_index, const DTYPE2 facet_orth_dir, 
   const bool facet_side) const
  {
    return(is_cube_facet_on_grid_boundary
           (Dimension(), AxisSize(), cube_index, facet_orth_dir, facet_side));
  }


  // Print vertex coordinates (mainly for debugging).
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, typename VTYPE0>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  PrintCoord(OSTREAM_TYPE & out, const VTYPE0 iv) const
  {
    const DTYPE dimension = this->Dimension();
    IJK::ARRAY<ATYPE> coord(dimension);

    this->ComputeCoord(iv, coord.Ptr());
    
    out << "(";
    for (DTYPE d = 0; d < dimension; d++) {
      out << coord[d];
      if (d+1 < dimension) { out << ","; }
    }
    out << ")";
  }


  // Print vertex coordinates (mainly for debugging).
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, typename VTYPE0>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  PrintCoord(OSTREAM_TYPE & out, const char * prefix,
             const VTYPE0 iv, const char * suffix) const
  {
    out << prefix;
    PrintCoord(out, iv);
    out << suffix;
  }


  // Print vertex index and coordinates (mainly for debugging).
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, typename VTYPE0>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  PrintIndexAndCoord(OSTREAM_TYPE & out, const VTYPE0 iv) const
  {
    out << iv << " ";
    PrintCoord(out, iv);
  }


  // Print vertex index and coordinates (mainly for debugging).
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, typename VTYPE0>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  PrintIndexAndCoord(OSTREAM_TYPE & out, const char * prefix,
                     const VTYPE0 iv, const char * suffix) const
  {
    out << prefix;
    PrintIndexAndCoord(out, iv);
    out << suffix;
  }


  // Print two vertex indices and coordinates (mainly for debugging).
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, typename VTYPE0, typename VTYPE1>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  PrintIndexAndCoord(OSTREAM_TYPE & out, const char * text0,
                     const VTYPE0 iv0, const char * text1, 
                     const VTYPE1 iv1, const char * text2) const
  {
    out << text0;
    PrintIndexAndCoord(out, iv0);
    out << text1;
    PrintIndexAndCoord(out, iv1);
    out << text2;
  }


  /// Print three vertex indices and coordinates (mainly for debugging).
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename OSTREAM_TYPE, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  PrintIndexAndCoord(OSTREAM_TYPE & out, const char * text0,
                     const VTYPE0 iv0, const char * text1, 
                     const VTYPE1 iv1, const char * text2,
                     const VTYPE2 iv2, const char * text3) const
  {
    out << text0;
    PrintIndexAndCoord(out, iv0);
    out << text1;
    PrintIndexAndCoord(out, iv1);
    out << text2;
    PrintIndexAndCoord(out, iv2);
    out << text3;
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckDimension(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2,
                 const char * grid_label1, const char * grid_label2,
                 IJK::ERROR & error) const
  {
    if (this->dimension != grid2.Dimension()) {
      error.AddMessage("Mismatch of volume dimensions.");
      error.AddMessage("  ", grid_label1, " has dimension ",
                       this->dimension, ".");
      error.AddMessage("  ", grid_label2, " has dimension ",
                       grid2.Dimension(), ".");
      return(false);
    }

    return(true);
  }


  /// @brief Return true if grid dimension and axis size match parameters.
  /// @param dimension  Dimension.
  /// @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
  /// @param[out] error Error message if grid dimension or axis_size do not match.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::Check
  (const DTYPE dimension, const ATYPE * axis_size, IJK::ERROR & error) const
  {
    if (dimension != this->dimension) {
      error.AddMessage("Incorrect grid dimension ", this->dimension, ".");
      error.AddMessage("  Dimension should be ", dimension, ".");
      return(false);
    }

    for (int d = 0; d < dimension; d++) {
      if (axis_size[d] != this->axis_size[d]) {
        error.AddMessage("Illegal axis size[", d, "] = ", 
                         this->axis_size[d], ".");
        error.AddMessage("  Axis size[", d, "] should be ", axis_size[d], ".");
        return(false);
      }
    }

    NTYPE num_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_vertices);

    if (num_vertices != this->num_vertices) {
      error.AddMessage("Incorrect number of grid vertices ", 
                       this->num_vertices, ".");
      error.AddMessage("  Number of grid vertices should be ", 
                       num_vertices, ".");
      return(false);
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  Check(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2,
        IJK::ERROR & error) const
  {
    return(Check(grid2.Dimension(), grid2.AxisSize(), error));
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  Check(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2,
        const char * grid_label1, const char * grid_label2,
        IJK::ERROR & error) const
  {
    if (!CheckDimension(grid2, grid_label1, grid_label2, error))
      { return(false); }

    for (int d = 0; d < dimension; d++) {
      if (this->AxisSize(d) != grid2.AxisSize(d)) {
        error.AddMessage("Mismatch of axis size ", d, ".");
        error.AddMessage("  ", grid_label1, " axis_size[", d,
                         "] = ", this->AxisSize(d), ".");
        error.AddMessage("  ", grid_label2, " axis_size[", d,
                         "] = ", grid2.AxisSize(d), ".");
        return(false);
      }
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCoord(const GTYPE * coord, IJK::ERROR & error) const
  {
    for (DTYPE d = 0; d < Dimension(); d++) {
      if (coord[d] < 0) {
        error.AddMessage("Coordinate should be non-negative.");
        error.AddMessage("  coord[", d, "] = ", coord[d], ".");
        return(false);
      }
      else if (AxisSize(d) < 1) {
        error.AddMessage("Coordinate ", d, " is out of bounds.");
        error.AddMessage("  coord[", d, "] = ", coord[d], ".");
        error.AddMessage
          ("  axis size[", d, "] = ", AxisSize(d), 
           " so all coordinates on axis ", d, " are out of bounds.");
        return(false);

      }
      else if (coord[d] >= AxisSize(d)) {
        error.AddMessage("Coordinate ", d, " is out of bounds.");
        error.AddMessage("  coord[", d, "] = ", coord[d], 
                         ".  axis size[", d, "] = ", AxisSize(d), ".");
        error.AddMessage("  coord[", d, 
                         "] should be less than axis size[", d, "].");
        return(false);
      }
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCoord(const std::vector<GTYPE> & coord, IJK::ERROR & error) const
  {
    return(this->CheckCoord(&(coord[0]), error));
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCubeCoord(const GTYPE * coord, IJK::ERROR & error) const
  {
    for (DTYPE d = 0; d < Dimension(); d++) {
      if (coord[d] < 0) {
        error.AddMessage("Cube coordinates should be non-negative.");
        error.AddMessage("  coord[", d, "] = ", coord[d], ".");
        return(false);
      }
      else if (AxisSize(d) < 2) {
        error.AddMessage("Cube coordinate ", d, " is out of bounds.");
        error.AddMessage("  coord[", d, "] = ", coord[d], ".");
        error.AddMessage
          ("  axis size[", d, "] = ", AxisSize(d), 
           " so all cube coordinates on axis ", d, " are out of bounds.");
        return(false);

      }
      else if (coord[d]+1 >= AxisSize(d)) {
        error.AddMessage("Cube coordinate ", d, " is out of bounds.");
        error.AddMessage("  coord[", d, "] = ", coord[d], 
                         ".  axis size[", d, "] = ", AxisSize(d), ".");
        error.AddMessage("  coord[", d, 
                         "] should be less than (axis size[", d, "]-1).");
        return(false);
      }
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename GTYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCubeCoord(const std::vector<GTYPE> & coord, IJK::ERROR & error) const
  {
    return(this->CheckCubeCoord(&(coord[0]), error));
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename ITYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckVertexIndex(const ITYPE vertex_index, IJK::ERROR & error) const
  {
    if (vertex_index < 0) {
      error.AddMessage("Illegal vertex index ", vertex_index, ".");
      error.AddMessage("  Vertex index should be non-negative.");
      return(false);
    }

    if (vertex_index >= NumVertices()) {
      error.AddMessage("Illegal vertex index ", vertex_index, ".");
      error.AddMessage
        ("  Vertex index should be less than number of grid vertices.");
      error.AddMessage("  Number of grid vertices = ", 
                       NumVertices(), ".");
      return(false);
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename ITYPE>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::
  CheckCubeIndex(const ITYPE cube_index, IJK::ERROR & error) const
  {
    if (cube_index < 0) {
      error.AddMessage("Illegal cube index ", cube_index, ".");
      error.AddMessage("  Cube index should be non-negative.");
      return(false);
    }

    for (DTYPE d = 0; d < Dimension(); d++) {
      if (AxisSize(d) <= 0) {
        // Never reached since NumVertices() == 0, but just in case.
        error.AddMessage("Illegal cube index ", cube_index, ".");
        if (AxisSize(d) <= 0) {
          error.AddMessage("  Grid axis size[", d, "] = ",
                           AxisSize(d), " so grid has no cubes.");
          return(false);
        }
      }
    }

    if (cube_index >= NumVertices()) {
      error.AddMessage("Illegal cube index ", cube_index, ".");
      error.AddMessage
        ("  Cube index should be less than number of grid vertices.");
      error.AddMessage("  Number of grid vertices = ", 
                       NumVertices(), ".");
      return(false);
    }

    VTYPE iv = cube_index;
    for (DTYPE d = 0; d < Dimension(); d++) {
      // Note: Already checked that AxisSize(d) > 0.
      VTYPE c = iv%AxisSize(d);
      iv = iv/AxisSize(d);
      if (c+1 >= AxisSize(d)) {
        error.AddMessage("  Vertices on right/top of cube with index ",
                         cube_index, " would have coordinate[", d,
                         "] = ", c+1, ".");
        error.AddMessage
          ("  Maximum coordinate[", d, "] in grid is ", AxisSize(d)-1, ".");
        return(false);
      }
    }

    return(true);
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  template <typename VTYPE2, typename ATYPE2>
  bool GRID<DTYPE,ATYPE,VTYPE,NTYPE>::CheckContainsRegion
  (const VTYPE2 region_v0, const ATYPE2 * region_axis_size,
   IJK::ERROR & error) const
  {
    const DTYPE dimension = this->Dimension();

    if (region_v0 < 0 || region_v0 >= this->NumVertices()) {
      error.AddMessage("Illegal lower/leftmost region vertex ",
                       region_v0, ".");
      return(false);
    }

    IJK::ARRAY<VTYPE> coord0(dimension);
    ComputeCoord(region_v0, coord0.Ptr());

    for (DTYPE d = 0; d < dimension; d++) {
      if (region_axis_size[d] < 0) {
        error.AddMessage("Illegal region_axis_size[", d, 
                         "] = ", region_axis_size[d], ".");
        return(false); 
      }

      if (coord0[d] + region_axis_size[d] > this->AxisSize(d)) {
        error.AddMessage("Error.  Region extends beyond grid.");
        error.AddMessage("  lower/leftmost coord[", d, 
                         "] = ", coord0[d], ".");
        error.AddMessage("  region_axis_size[", d, "] = ",
                         region_axis_size[d], ".");
        error.AddMessage("  grid axis_size[", d, "] = ",
                         this->AxisSize(d), ".");
        return(false); 
      }
    }

    return(true);
  }


  // *****************************************************************
  // TEMPLATE CLASS GRID_PLUS MEMBER FUNCTIONS
  // *****************************************************************

  /// Constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::GRID_PLUS
  (const DTYPE dimension, const ATYPE * axis_size):
    GRID<DTYPE,ATYPE,VTYPE,NTYPE> (dimension,axis_size)
  {
    InitLocal();
  }

  /// Default constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::GRID_PLUS()
  {
    InitLocal();
  }

  /// Constructor from another GRID_PLUS grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::
  GRID_PLUS(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2):
    GRID<DTYPE, ATYPE, VTYPE, NTYPE>(grid2)
  {
    InitLocal();
  }

  /// Constructor from another grid with same type.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::
  GRID_PLUS(const GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> & grid2):
    GRID<DTYPE,ATYPE,VTYPE,NTYPE>(grid2)
  {
    InitLocal();
  }

  /// \brief Set all local (not inherited) arrays to NULL.
  /// Set all local variables to 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::ZeroLocal()
  {
    this->axis_increment = NULL;
    this->cube_vertex_increment = NULL;
    this->facet_vertex_increment = NULL;
    this->unit_cube_coord = NULL;
    this->num_cube_vertices = 0;
    this->num_facet_vertices = 0;
    this->num_cube_ridge_vertices = 0;
    this->num_cube_facets = 0;
    this->num_cube_edges = 0;
    facet0_orientation = true;
  }

  /// \brief Initialize data structures in GRID_PLUS.
  /// @pre \a dimension and \a axis_size are already set.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::InitLocal()
  {
    ZeroLocal();
    if (this->Dimension() > 0) { Create(); };
  }

  /// Destructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::~GRID_PLUS()
  {
    FreeLocal();
  }


  /// @brief Allocate arrays and compute data in GRID_PLUS.
  /// @pre \a dimension and \a axis_size[] are already set.
  /// @pre All other arrays are set to NULL.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::Create()
  {
    const DTYPE dimension = this->Dimension();
    IJK::PROCEDURE_ERROR error("GRID_PLUS::Create");

    if (!check_dimension(dimension, error)) { throw error; };
    if (!IJK::check_is_NULL(axis_increment, "axis_increment", error))
      { throw error; }
    if (!IJK::check_is_NULL
        (cube_vertex_increment, "cube_vertex_increment", error))
      { throw error; }
    if (!IJK::check_is_NULL
        (facet_vertex_increment, "facet_vertex_increment", error))
      { throw error; }
    if (!IJK::check_is_NULL
        (unit_cube_coord, "unit_cube_coord", error))
      { throw error; }

    this->num_cube_vertices = compute_num_cube_vertices(dimension);
    this->num_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);
    this->num_cube_ridge_vertices = 
      compute_num_cube_ridge_vertices(dimension);
    this->num_cube_facets = compute_num_cube_facets(dimension);
    this->num_cube_edges = compute_num_cube_edges(dimension);

    // Initialize.
    this->axis_increment = NULL;
    this->cube_vertex_increment = NULL;
    this->ridge_vertex_increment = NULL;
    this->facet_vertex_increment = NULL;
    this->unit_cube_coord = NULL;

    if (dimension > 0) {
      this->axis_increment = new VTYPE[dimension];
      compute_increment
         (dimension, this->AxisSize(), this->axis_increment);
    }

    if (num_cube_vertices > 0) {
      // Always true, but...
      this->cube_vertex_increment = new VTYPE[num_cube_vertices];
      compute_cube_vertex_increment
         (dimension, this->AxisIncrement(), this->cube_vertex_increment);
    }

    if (this->NumCubeFacets() > 0) {
      this->facet_vertex_increment =
        new VTYPE[num_facet_vertices*this->NumCubeFacets()];

      for (DTYPE ifacet = 0; ifacet < this->NumCubeFacets(); ifacet++) {
        compute_facet_vertex_increment
          (dimension, ifacet, this->cube_vertex_increment,
           this->facet_vertex_increment+this->num_facet_vertices*ifacet);
      }
    }

    if (num_cube_ridge_vertices > 0) {
      this->ridge_vertex_increment =
        new VTYPE[num_cube_ridge_vertices*dimension*dimension];

      // Initialize ridge_vertex_increment to zero.
      for (NTYPE i = 0; i < this->num_cube_ridge_vertices*dimension*dimension;
           i++)
        { this->ridge_vertex_increment[i] = 0; }

      // Compute ridge_vertex_increment.
      for (DTYPE orth_dir1 = 0; orth_dir1 < dimension; orth_dir1++)
        for (DTYPE orth_dir0 = 0; orth_dir0 < dimension; orth_dir0++)
          if (orth_dir0 != orth_dir1) {
            compute_ridge_vertex_increment
              (dimension, orth_dir0, orth_dir1,
               this->cube_vertex_increment,
               this->ridge_vertex_increment+
               this->num_cube_ridge_vertices*(orth_dir0+dimension*orth_dir1));
          }
    }

    if (num_cube_vertices > 0) {
      // Always true, but..l
      this->unit_cube_coord = new NTYPE[num_cube_vertices*dimension];
      compute_unit_cube_vertex_coord
        (dimension, this->unit_cube_coord);
    }

    // Set facet0 orientation.
    // - If grid dimension is 3, then true is counter-clockwise orientation.
    facet0_orientation = true;
  }


  /// Free memory in the derived typename GRID_PLUS.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::FreeLocal()
  {
    if (axis_increment != NULL) { delete [] axis_increment; }
    if (cube_vertex_increment != NULL) 
      { delete [] cube_vertex_increment; };
    if (facet_vertex_increment != NULL)
      { delete [] facet_vertex_increment; };
    if (unit_cube_coord != NULL) { delete [] unit_cube_coord; };
    ZeroLocal();
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2>
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    IJK::PROCEDURE_ERROR error("GRID_PLUS::SetSize");

    FreeLocal();

    GRID<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize(dimension, axis_size);
    Create();
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  /// @brief Return the direction of a grid edge.
  /// @pre endpoint0 and endpoint1 are the endpoints of a grid edge.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename VTYPE2, typename VTYPE3>
  DTYPE GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::EdgeDirection
  (const VTYPE2 endpoint0, const VTYPE3 endpoint1) const
  {
    VTYPE2 diff = endpoint0-endpoint1;
    if (diff < 0) { diff = -diff; }
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      if (diff == this->AxisIncrement(d)) { return(d); }
    }

    IJK::PROCEDURE_ERROR error("GRID_PLUS::EdgeDirection");
    if (diff == 0) {
      error.AddMessage("Programming error. Edge endpoints are identical.");
      throw error;
    }
    else {
      error.AddMessage("Programming error. Vertices are not edge endpoints.");
      throw error;
    }
  }


  // Compute vertex neighbors.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename VTYPE0, typename VTYPE1, typename DIST_TYPE>
  void GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::GetVertexNeighbors
  (const VTYPE0 iv0, const DIST_TYPE distance, std::vector<VTYPE1> & vlist) 
    const
  {
    get_grid_vertices_in_neighborhood
      (this->Dimension(), this->AxisSize(), this->AxisIncrement(), 
       iv0, distance, vlist);
  }


  // Return facet of icube0 shared by icube0 and icube1.
  // @pre icube0 and icube1 share some facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename ITYPE>
  NTYPE GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::
  SharedFacet(const ITYPE icube0, const ITYPE icube1) const
  {
    const DTYPE dimension = this->Dimension();

    if (icube0 < icube1) {
      const ITYPE cube_diff = icube1 - icube0;

      for (DTYPE d = 0; d < dimension; d++) {
        if (AxisIncrement(d) == cube_diff)
          { return (d+dimension); }
      }
    }
    else {
      const ITYPE cube_diff = icube0 - icube1;
      for (DTYPE d = 0; d < dimension; d++) {
        if (AxisIncrement(d) == cube_diff)
          { return d; }
      }
    }

    IJK::PROCEDURE_ERROR error("GRID_PLUS::SharedFacet");
    error.AddMessage
      ("Programming error. Cubes do not share a facet.");
    error.AddMessage
      ("  Cube ", icube0, " does not share a facet with cube ",
       icube1, ".");
    throw error;
  }


  // Return direction orthogonal to facet shared by icube0 and icube1.
  // @pre icube0 and icube1 share some facet.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE> 
  template <typename ITYPE>
  NTYPE GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::
  SharedFacetOrthDir(const ITYPE icube0, const ITYPE icube1) const
  {
    const DTYPE dimension = this->Dimension();

    if (icube0 < icube1) {
      const ITYPE cube_diff = icube1 - icube0;

      for (DTYPE d = 0; d < dimension; d++) {
        if (AxisIncrement(d) == cube_diff)
          { return d; }
      }
    }
    else {
      const ITYPE cube_diff = icube0 - icube1;
      for (DTYPE d = 0; d < dimension; d++) {
        if (AxisIncrement(d) == cube_diff)
          { return d; }
      }
    }

    IJK::PROCEDURE_ERROR error("GRID_PLUS::SharedFacet");
    error.AddMessage
      ("Programming error. Cubes do not share a facet.");
    error.AddMessage
      ("  Cube ", icube0, " does not share a facet with cube ",
       icube1, ".");
    throw error;
  }
    

  // Return d'th coordinate of vertex \a iv.
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  VTYPE GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::
  CoordD(const VTYPE iv, const DTYPE d) const
  {
    return compute_coord_d
        (iv, d, this->Dimension(), this->AxisIncrement());
  }


  // *****************************************************************
  // TEMPLATE CLASS GRID_NEIGHBORS MEMBER FUNCTIONS
  // *****************************************************************

  /// Constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::GRID_NEIGHBORS
  (const DTYPE dimension, const ATYPE * axis_size):
    GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> (dimension,axis_size)
  {
    InitLocal();
  }

  /// Default constructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::GRID_NEIGHBORS()
  {
    InitLocal();
  }

  /// Constructor from another GRID_NEIGHBORS grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::
  GRID_NEIGHBORS(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2):
    GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> (grid2)
  {
    InitLocal();
  }

  /// Constructor from another grid with same type.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::
  GRID_NEIGHBORS
  (const GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE, NTYPE> & grid2):
    GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE> (grid2)
  {
    InitLocal();
  }


  /// @brief Set all local (not inherited) arrays to NULL.
  /// - Set all local variables to 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::ZeroLocal()
  {
    this->vertex_neighborC = NULL;
    this->vertex_neighborE = NULL;
    this->cube_neighborE = NULL;
    this->cube_neighborV = NULL;
    this->facet_neighborC = NULL;
    this->edge_neighborF2 = NULL;
    this->num_vertex_neighborsC = 0;
    this->num_vertex_neighborsE = 0;
    this->num_facet_neighborsC = 0;
    this->num_edge_neighborsF2 = 0;
  }

  /// \brief Initialize data structures in GRID_NEIGHBORS
  /// @pre \a dimension and \a axis_size are already set
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::InitLocal()
  {
    ZeroLocal();
    if (this->Dimension() > 0) { CreateLocal(); };
  }

  /// Destructor.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::~GRID_NEIGHBORS()
  {
    FreeLocal();
  }


  /// @brief Allocate arrays and compute data in GRID_NEIGHBORS.
  /// @pre \a dimension and \a axis_size[] are already set.
  /// @pre All other arrays are set to NULL.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::CreateLocal()
  {
    const DTYPE dimension = this->Dimension();
    const ATYPE * axis_size = this->AxisSize();
    IJK::PROCEDURE_ERROR error("GRID_NEIGHBORS::CreateLocal");

    FreeLocal();

    vertex_neighborC = NULL;
    compute_num_vertex_neighborsC(dimension, num_vertex_neighborsC);
    if (num_vertex_neighborsC > 0) {
      vertex_neighborC = new DIFFTYPE[num_vertex_neighborsC];
      compute_vertex_neighborC(dimension, axis_size, vertex_neighborC);
    }

    vertex_neighborE = NULL;
    compute_num_vertex_neighborsE(dimension, num_vertex_neighborsE);
    if (num_vertex_neighborsE > 0) {
      vertex_neighborE = new DIFFTYPE[num_vertex_neighborsE];
      compute_vertex_neighborE(dimension, axis_size, vertex_neighborE);
    }

    cube_neighborV = NULL;
    if (this->NumCubeVertices() > 0) {
      // Always true, but...
      cube_neighborV = new DIFFTYPE[this->NumCubeVertices()];
      compute_cube_neighborV(dimension, axis_size, cube_neighborV);
    }

    cube_neighborE = NULL;
    if (this->NumCubeEdges() > 0) {
      cube_neighborE = new DIFFTYPE[this->NumCubeEdges()];
      compute_cube_neighborE(dimension, axis_size, cube_neighborE);
    }

    facet_neighborC = NULL;
    compute_num_facet_neighborsC(dimension, num_facet_neighborsC);
    if (num_facet_neighborsC > 0) {
      facet_neighborC =
        new DIFFTYPE[num_facet_neighborsC*(dimension)];
      compute_facet_neighborC(dimension, axis_size, facet_neighborC);
    }

    edge_neighborF2 = NULL;
    compute_num_edge_neighborsF2
      (dimension, num_edge_neighborsF2);
    if (num_edge_neighborsF2 > 0) {
      this->edge_neighborF2 =
        new DIFFTYPE[num_edge_neighborsF2*(dimension)];
      compute_edge_neighborF2(dimension, axis_size, edge_neighborF2);
    }
  }


  // \brief Free all local (not inherited) arrays.
  // - Set all local arrays to NULL and variables to 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::FreeLocal()
  {
    if (vertex_neighborC != NULL) { delete [] vertex_neighborC; };
    if (vertex_neighborE != NULL) { delete [] vertex_neighborE; };
    if (cube_neighborE != NULL) { delete [] cube_neighborE; };
    if (cube_neighborV != NULL) { delete [] cube_neighborV; };
    if (facet_neighborC != NULL) { delete [] facet_neighborC; };
    if (edge_neighborF2 != NULL) { delete [] edge_neighborF2; };
    ZeroLocal();
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2>
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    IJK::PROCEDURE_ERROR error("GRID_NEIGHBORS::SetSize");

    FreeLocal();

    GRID_PLUS<DTYPE,ATYPE,VTYPE,NTYPE>::SetSize(dimension, axis_size);
    CreateLocal();
  }

  template <typename DTYPE, typename ATYPE, typename VTYPE, 
            typename DIFFTYPE, typename NTYPE> 
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, typename NTYPE2>
  void GRID_NEIGHBORS<DTYPE,ATYPE,VTYPE,DIFFTYPE,NTYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }


  // *****************************************************************
  // TEMPLATE CLASS GRID_SPACING MEMBER FUNCTIONS
  // *****************************************************************

  /// Constructor.
  template <typename SP_TYPE, typename GRID_TYPE>
  GRID_SPACING<SP_TYPE, GRID_TYPE>::GRID_SPACING():GRID_TYPE()
  {
    InitLocal();
  }

  /// Constructor.
  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2>
  GRID_SPACING<SP_TYPE, GRID_TYPE>::
  GRID_SPACING(const DTYPE2 dimension, const ATYPE2 * axis_size):
    GRID_TYPE(dimension, axis_size)
  {
    InitLocal();
  }

  /// Constructor from another grid.
  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
            typename NTYPE2>
  GRID_SPACING<SP_TYPE, GRID_TYPE>::
  GRID_SPACING(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2):
    GRID_TYPE(grid2)
  {
    InitLocal();
  }

  /// Constructor from another grid with same type.
  template <typename SP_TYPE, typename GRID_TYPE>
  GRID_SPACING<SP_TYPE, GRID_TYPE>::
  GRID_SPACING(const GRID_SPACING<SP_TYPE, GRID_TYPE> & grid2):
    GRID_TYPE(grid2)
  {
    InitLocal();
    SetSpacing(grid2.SpacingPtrConst());
  }

  /// Copy function.
  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename GTYPE2>
  const GRID_SPACING<SP_TYPE, GRID_TYPE> &
  GRID_SPACING<SP_TYPE, GRID_TYPE>::Copy(const GTYPE2 & right)
  {
    if (&right != this) {         // avoid self-assignment
      FreeLocal();
      GTYPE2::Copy(right);
      InitLocal();
      SetSpacing(right.SpacingPtrConst());
    }
    return(*this);
  }

  /// Initialize.
  template <typename SP_TYPE, typename GRID_TYPE>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::InitLocal()
  {
    spacing = NULL;
    CreateLocal();
    SetAllSpacing(1);
  }

  /// Create data structure.
  template <typename SP_TYPE, typename GRID_TYPE>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::CreateLocal()
  {
    FreeLocal();
    if (this->Dimension() != 0) {
      spacing = new SP_TYPE[this->dimension];
    }
    else
    { spacing = NULL; }
  }

  /// Free memory.
  template <typename SP_TYPE, typename GRID_TYPE>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::FreeLocal()
  {
    if (spacing != NULL) {
      delete [] spacing;
      spacing = NULL;
    };
  }

  /// Set spacing along all axes to c.
  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename SP_TYPE2>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::SetAllSpacing(const SP_TYPE2 c)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) 
      { SetSpacing(d, c); }
  }

  /// Set spacing[d] to spacing2[d].
  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename SP_TYPE2>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::SetSpacing(const SP_TYPE2 * spacing2)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) 
      { SetSpacing(d, spacing2[d]); }
  }

  /// Set spacing[d] to (c*spacing2[d]).
  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename SCALE_TYPE, typename SP_TYPE2>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::
  SetSpacing(const SCALE_TYPE c, const SP_TYPE2 * spacing2)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) 
      { SetSpacing(d, c*spacing2[d]); }
  }

  /// Set spacing[d] to grid2.Spacing(d).
  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename SP_TYPE2, typename GRID_TYPE2>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::
  SetSpacing(const GRID_SPACING<SP_TYPE2, GRID_TYPE2> & grid2)
  {
    SetSpacing(grid2.SpacingPtrConst());
  }

  /// Return true if (Spacing(d) == 1) for all d.
  template <typename SP_TYPE, typename GRID_TYPE>
  bool GRID_SPACING<SP_TYPE, GRID_TYPE>::IsUnitSpacing() const
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) {
      if (Spacing(d) != 1) { return(false); }
    }

    return(true);
  }

  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    FreeLocal();
    GRID_TYPE::SetSize(dimension, axis_size);
    CreateLocal();
    SetAllSpacing(1);
  }

  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
            typename NTYPE2>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }

  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename VTYPE2, typename CTYPE>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::
  ComputeScaledCoord(const VTYPE2 iv, CTYPE * coord) const
  {
    compute_scaled_coord
      (iv, this->Dimension(), this->AxisSize(), SpacingPtrConst(), coord);
  }

  template <typename SP_TYPE, typename GRID_TYPE>
  template <typename VTYPE2, typename CTYPE>
  void GRID_SPACING<SP_TYPE, GRID_TYPE>::
  ComputeCubeCenterScaledCoord(const VTYPE2 iv, CTYPE * coord) const
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    ComputeScaledCoord(iv, coord);
    for (DTYPE d = 0; d < this->Dimension(); d++) 
      { coord[d] += Spacing(d)/2.0; }
  }

  /// Destructor.
  template <typename SP_TYPE, typename GRID_TYPE>
  GRID_SPACING<SP_TYPE, GRID_TYPE>::~GRID_SPACING()
  {
    FreeLocal();
  }


  // *****************************************************************
  // TEMPLATE CLASS GEOM_GRID MEMBER FUNCTIONS
  // *****************************************************************

  // Constructor.
  template <typename CTYPE, typename GRID_TYPE>
  GEOM_GRID<CTYPE, GRID_TYPE>::GEOM_GRID():GRID_SPACING<CTYPE,GRID_TYPE>()
  {
    InitLocal();
  }


  // Constructor.
  template <typename CTYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2>
  GEOM_GRID<CTYPE, GRID_TYPE>::
  GEOM_GRID(const DTYPE2 dimension, const ATYPE2 * axis_size):
    GRID_SPACING<CTYPE,GRID_TYPE>(dimension, axis_size)
  {
    InitLocal();
  }


  // Constructor from another grid.
  template <typename CTYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
            typename NTYPE2>
  GEOM_GRID<CTYPE, GRID_TYPE>::
  GEOM_GRID(const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2):
    GRID_SPACING<CTYPE,GRID_TYPE>(grid2)
  {
    InitLocal();
  }


  // Constructor from another grid with same type.
  template <typename CTYPE, typename GRID_TYPE>
  GEOM_GRID<CTYPE, GRID_TYPE>::
  GEOM_GRID(const GEOM_GRID<CTYPE, GRID_TYPE> & grid2):
    GRID_SPACING<CTYPE,GRID_TYPE>(grid2)
  {
    InitLocal();
    SetOrigin(grid2.OriginPtrConst());
  }


  // Copy function.
  template <typename CTYPE, typename GRID_TYPE>
  template <typename GTYPE2>
  const GEOM_GRID<CTYPE, GRID_TYPE> &
  GEOM_GRID<CTYPE, GRID_TYPE>::Copy(const GTYPE2 & right)
  {
    if (&right != this) {         // avoid self-assignment
      FreeLocal();
      InitLocal();
      GRID_SPACING<CTYPE,GRID_TYPE>::Copy(right);
      SetOrigin(right.OriginPtrConst());
    }
    return(*this);
  }


  // Initialize.
  template <typename CTYPE, typename GRID_TYPE>
  void GEOM_GRID<CTYPE, GRID_TYPE>::InitLocal()
  {
    origin = NULL;
    CreateLocal();
    SetAllOrigin(0);
  }


  // Create data structure.
  template <typename CTYPE, typename GRID_TYPE>
  void GEOM_GRID<CTYPE, GRID_TYPE>::CreateLocal()
  {
    FreeLocal();
    if (this->dimension > 0) {
       origin = new CTYPE[this->dimension];
    }
  }


  // Free memory.
  template <typename CTYPE, typename GRID_TYPE>
  void GEOM_GRID<CTYPE, GRID_TYPE>::FreeLocal()
  {
    if (origin != NULL) {
      delete [] origin;
      origin = NULL;
    };
  }

  
  // Set all origin coordinates to c.
  template <typename CTYPE, typename GRID_TYPE>
  template <typename CTYPE2>
  void GEOM_GRID<CTYPE, GRID_TYPE>::SetAllOrigin(const CTYPE2 c)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) 
      { SetOrigin(d, c); }
  }


  // Set origin[d] to origin2[d].
  template <typename CTYPE, typename GRID_TYPE>
  template <typename CTYPE2>
  void GEOM_GRID<CTYPE, GRID_TYPE>::SetOrigin(const CTYPE2 * origin2)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    for (DTYPE d = 0; d < this->dimension; d++) 
      { SetOrigin(d, origin2[d]); }
  }
  

  // Set origin[d] to grid2.Origin(d).
  template <typename CTYPE, typename GRID_TYPE>
  template <typename CTYPE2, typename GRID_TYPE2>
  void GEOM_GRID<CTYPE, GRID_TYPE>::
  SetOrigin(const GEOM_GRID<CTYPE2, GRID_TYPE2> & grid2)
  {
    SetOrigin(grid2.OriginPtrConst());
  }


  template <typename CTYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2>
  void GEOM_GRID<CTYPE, GRID_TYPE>::SetSize
  (const DTYPE2 dimension, const ATYPE2 * axis_size)
  {
    FreeLocal();
    GRID_TYPE::SetSize(dimension, axis_size);
    CreateLocal();
    SetAllOrigin(0);
  }

  template <typename CTYPE, typename GRID_TYPE>
  template <typename DTYPE2, typename ATYPE2, typename VTYPE2, 
            typename NTYPE2>
  void GEOM_GRID<CTYPE, GRID_TYPE>::SetSize
  (const GRID<DTYPE2,ATYPE2,VTYPE2,NTYPE2> & grid2)
  {
    SetSize(grid2.Dimension(), grid2.AxisSize());
  }


  // Destructor.
  template <typename CTYPE, typename GRID_TYPE>
  GEOM_GRID<CTYPE, GRID_TYPE>::~GEOM_GRID()
  {
    FreeLocal();
  }

}


#endif
