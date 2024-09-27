/*!
 *  @file ijkgrid_func.tpp
 *  @brief ijk templates defining regular grid functions.
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

#ifndef _IJKGRID_FUNC_TPP
#define _IJKGRID_FUNC_TPP

#include <algorithm>
#include <limits>
#include <vector>

#include "ijk.tpp"
#include "ijkcube.tpp"


namespace IJK {

  // **************************************************
  // TEMPLATE CLASS FACET_LIST
  // **************************************************

  /// Base class for list of facets.
  template <typename DTYPE>
  class FACET_LIST {

  protected:
    FACET_LIST(){};        ///< Constructor.

  public:
    inline bool Contains(const DTYPE d) const;
  };

  /// list containing zero facets
  template <typename DTYPE>
  class FACET_LIST0:public FACET_LIST<DTYPE> {

  public:
    FACET_LIST0() {};      ///< Constructor.

    /// @brief Return true if facet list contains facet \a d.
    /// @param d Facet index (= direction orthogonal to facet.)
    inline bool Contains(const DTYPE d) const { return(false); };
  };

  /// List containing one facet.
  template <typename DTYPE>
  class FACET_LIST1:public FACET_LIST<DTYPE> {

  protected:
    DTYPE facet0;       ///< Direction orthogonal to facet0.

  public:

    /// Set facet1.
    inline void Set(const DTYPE facet0) { this->facet0 = facet0; };

    /// Constructor.
    FACET_LIST1(const DTYPE facet0) { Set(facet0); };

    /// @brief Return true if facet list contains facet \a d.
    /// @param d Facet index (= direction orthogonal to facet.)
    inline bool Contains(const DTYPE d) const { 
      if (d == facet0) { return(true); };
      return(false);
    };
  };

  /// List containing two facets.
  template <typename DTYPE>
  class FACET_LIST2:public FACET_LIST<DTYPE> {

  protected:
    DTYPE facet0;       ///< Direction orthogonal to facet0.
    DTYPE facet1;       ///< Direction orthogonal to facet1.

  public:

    /// Set facet0 and facet1.
    inline void Set(const DTYPE facet0, const DTYPE facet1) 
    { this->facet0 = facet0; this->facet1 = facet1; };

    /// Constructor.
    FACET_LIST2(const DTYPE facet0, const DTYPE facet1) 
    { Set(facet0, facet1); };

    /// @brief Return true if facet list contains facet \a d.
    /// @param d Facet index (= direction orthogonal to facet.)
    inline bool Contains(const DTYPE d) const { 
      if (d == facet0 || d == facet1) { return(true); };
      return(false);
    };
  };

  // **************************************************
  // TEMPLATE CLASS GRID_VERTEX_LIST
  // **************************************************


  /// Class containing list of grid vertices.
  template <typename VTYPE>
  class GRID_VERTEX_LIST {

  protected:
    VTYPE * vertex_list;         ///< Vertex list.

    /// @brief Allocated length of array vertex_list[].
    /// @invariant list_length >= num_vertices.
    /// - Note: list_length could be greater than num_vertices.
    VTYPE list_length;

    /// @brief Number of vertices in vertex_list[].
    /// @invariant num_vertices <= list_length.
    /// - Note: num_vertices could be less than list_length.
    VTYPE num_vertices;

    void Init();                 ///< Initialize class.
    void FreeAll();              ///< Free all allocated memory.
    void AllocateList(const VTYPE n); ///< Allocate vertex list.

  public:
    GRID_VERTEX_LIST() 
    { Init(); };
    ~GRID_VERTEX_LIST()
    { FreeAll(); }

    // get functions
    VTYPE VertexIndex(const VTYPE i) const
    { return(vertex_list[i]); }

    /// @brief Return number of vertices in list.
    /// - Note: NumVertices() could be less than ListLength().
    VTYPE NumVertices() const
    { return(num_vertices); }

    /// @brief Return allocated length of list.
    /// - Note: ListLength() could be greater than NumVertices().
    VTYPE ListLength() const
    { return(list_length); }
  };

  /// Class containing list of grid boundary vertices.
  template <typename VTYPE>
  class GRID_BOUNDARY_VERTEX_LIST:public GRID_VERTEX_LIST<VTYPE> {

  public:
    /// GRID_BOUNDARY_VERTEX_LIST constructor.
    /// @param grid Grid.
    template<typename GCLASS>
    GRID_BOUNDARY_VERTEX_LIST(const GCLASS & grid);
  };


  /// Class containing list of grid cubes.
  template <typename VTYPE>
  class GRID_CUBE_LIST:protected GRID_VERTEX_LIST<VTYPE> {

  public:
    GRID_CUBE_LIST() {};
    ~GRID_CUBE_LIST() {};

    // get functions
    VTYPE CubeIndex(const VTYPE i) const
    { return(this->VertexIndex(i)); }

    VTYPE NumCubes() const
    { return(this->NumVertices()); }

    VTYPE ListLength() const
    { return(GRID_VERTEX_LIST<VTYPE>::ListLength()); }
  };


  /// Class containing list of facet 0 cubes.
  template <typename VTYPE>
  class FACET0_CUBE_LIST:public GRID_CUBE_LIST<VTYPE> {

  public:
    template <typename GCLASS>
    FACET0_CUBE_LIST(const GCLASS & grid)
    { GetCubes(grid); }

    /// @brief Get cubes from grid and store in list.
    /// - Reallocates list if (current list length < num cubes in facet0)
    template <typename GCLASS>
    void GetCubes(const GCLASS & grid);
  };


  /// Class containing list of grid facets orthogonal to a given direction.
  template <typename DTYPE, typename VTYPE>
  class GRID_FACET_LIST:protected GRID_VERTEX_LIST<VTYPE> {

  protected:
    DTYPE orth_dir;

  public:
    /// Constructor.
    GRID_FACET_LIST() {};

    /// Constructor. 
    /// Set list to facets in lower/leftmost grid boundary facet
    /// orthogonal to orth_dir.
    template <typename GCLASS, typename DTYPE0>
    GRID_FACET_LIST(const GCLASS & grid, const DTYPE0 orth_dir)
    { GetFacetsInGridFacet(grid,orth_dir); }

    /// Constructor. 
    /// Set list to facets in grid boundary facet orthogonal to orth_dir.
    /// @param side Side of grid boundary.
    ///   - If false, facets lie on lower/leftmost grid boundary facet.
    ///   - If true, facets lie on upper/rightmost grid boundary facet.
    template <typename GCLASS, typename DTYPE0>
    GRID_FACET_LIST
    (const GCLASS & grid, const DTYPE0 orth_dir, const bool side)
    { GetFacetsInGridFacet(grid,orth_dir,side); }

    ~GRID_FACET_LIST() {};

    // get functions
    VTYPE FacetIndex(const VTYPE i) const
    { return(this->VertexIndex(i)); }

    VTYPE NumFacets() const
    { return(this->NumVertices()); }

    DTYPE OrthDir() const
    { return(orth_dir); }

    VTYPE ListLength() const
    { return(GRID_VERTEX_LIST<VTYPE>::ListLength()); }

    /// @brief Get facets from grid and store in list.
    /// - Reallocates list if (current list length < num facets)
    template <typename GCLASS, typename DTYPE0>
    void GetFacetsInGridFacet
    (const GCLASS & grid, const DTYPE0 orth_dir, 
     const bool side=false);
  };


  /// Class containing list of grid facet vertices.
  template <typename VTYPE>
  class FACET_VERTEX_LIST:public GRID_VERTEX_LIST<VTYPE> {

  public:

    /*!
     *  @brief FACET_VERTEX_LIST constructor.
     *  Get and store vertices in facet orthogonal to \a orth_dir.
     *  @param grid Grid.
     *  @param orth_dir Directional orthogonal to facet.
     *  @param allocate_max If true, allocate the list length
     *             to be max number of interior vertices over all facets.
     */
    template <typename GCLASS>
    FACET_VERTEX_LIST
    (const GCLASS & grid, const VTYPE orth_dir,
     const bool allocate_max=true);

    /// @brief Get vertices from grid and store in list.
    /// - Reallocates list if (current list length < num vertices in facet)
    template <typename GCLASS>
    void GetVertices(const GCLASS & grid, const VTYPE orth_dir);
  };

  /// Class containing list of grid facet interior vertices.
  template <typename VTYPE>
  class FACET_INTERIOR_VERTEX_LIST:public GRID_VERTEX_LIST<VTYPE> {

  public:
    /*!
     *  @brief FACET_INTERIOR_VERTEX_LIST constructor.
     *  - Get and store interior vertices in facet orthogonal to \a orth_dir.
     *  @param grid Grid.
     *  @param orth_dir Directional orthogonal to facet.
     *  @param allocate_max If true, allocate the list length
     *             to be max number of interior vertices over all facets.
     *  @param flag_dim1_facet_vertex If true, list the facet vertex
     *             if grid has dimension 1.
     */
    template <typename GCLASS>
    FACET_INTERIOR_VERTEX_LIST
    (const GCLASS & grid, const VTYPE orth_dir,
     const bool allocate_max=true,
     const bool flag_dim1_facet_vertex=false);

    /*!
     *  @brief Get vertices from grid and store in list.
     *  - Reallocates list if (current list length < num vertices in facet)
     *  @param flag_dim1_facet_vertex If true, list the facet vertex
     *             if grid has dimension 1.
     */
    template <typename GCLASS>
    void GetVertices(const GCLASS & grid, const VTYPE orth_dir,
                     const bool flag_dim1_facet_vertex=false);
  };


  // *****************************************************************
  //! @name inline UTILITY FUNCTIONS
  // *****************************************************************

  ///@{

  /// Integer divide.
  template <typename ATYPE, typename BTYPE>
  inline long integer_divide(const ATYPE a, const BTYPE b)
  { return(long(a)/long(b)); }

  /// Integer divide.
  inline long integer_divide(const long a, const long b)
  { return(a/b); }

  /// Integer divide.
  inline int integer_divide(const int a, const int b)
  { return(a/b); }
  
  /// Integer divide.
  inline unsigned long 
  integer_divide(const unsigned long a, const unsigned long b)
  { return(a/b); }

  /// Integer divide.
  inline unsigned int
  integer_divide(const unsigned int a, const unsigned int b)
  { return(a/b); }

  ///@}


  // *****************************************************************
  //! @name THROW ERROR ROUTINES
  // *****************************************************************

  ///@{

  /// Throw subsample period error.
  template <typename STRING_TYPE, typename PTYPE>
  void throw_subsample_period_error
  (const STRING_TYPE proc_name, const PTYPE subsample_period)
  {
    IJK::PROCEDURE_ERROR error(proc_name);
    error.AddMessage("Subsample period must be a positive integer.");
    error.AddMessage("  Subsampling period = ", subsample_period, ".");
    throw error;
  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: COUNTING AND INDEXING
  // *****************************************************************

  ///@{

  // Forward declaration.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir,
   NTYPE & num_vertices);


  /*!
   *  @brief Return coordinates of vertex \a iv.
   *  @param iv Vertex index.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param[out] coord Array. <em>coord[k]</em> = \a k'th coordinate of vertex \a iv.
   *  @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
   *  @pre \li Array coord[] is pre-allocated to size at least \a dimension.
   */
  template <typename VTYPE, typename DTYPE, typename ATYPE, typename GTYPE>
  void compute_coord(const VTYPE iv, const DTYPE dimension,
                     const ATYPE * axis_size, GTYPE * coord)
  {
    VTYPE k = iv;
    for (DTYPE d = 0; d < dimension; d++) {
      coord[d] = GTYPE(k % axis_size[d]);
      k = k / axis_size[d];
    };
  }


  /*!
   *  @brief Return d'th coordinate of vertex \a iv.
   *  @param iv Vertex index.
   *  @param d Axis index.
   *    @pre 0 <= \a d < \a dimension.
   *  @param dimension Dimension of grid.
   *  @param axis_increment[] Axis increment.
   *    - iv+axis_increment[jd] is the next vertex after vertex iv
   *      along axis jd.
   */
  template <typename VTYPE, typename DTYPE, typename ITYPE>
  VTYPE compute_coord_d(const VTYPE iv, const DTYPE d,
                        const DTYPE dimension, const ITYPE * axis_increment)
  {
    VTYPE coord_d = iv;

    if ((d+1) < dimension)
      { coord_d = coord_d % axis_increment[d+1]; }

    coord_d = coord_d / axis_increment[d];

    return coord_d;
  }


  /*!
   *  @brief Return coordinate of vertex \a iv.
   *  @param iv Vertex index.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param scale[]
   *     <em>scale[d]</em> = Scale along axis \a d.
   *  @param[out] scaled_coord[]
   *     <em>scaled_coord[k]</em> = \a k'th coordinate of vertex \a iv.
   *  @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
   *  @pre \li scale[d] > 0 for all \a d = 0,..., \a dimension-1.
   *  @pre \li Array coord[] is pre-allocated to size at least \a dimension.
   */
  template <typename VTYPE, typename DTYPE, typename ATYPE, 
            typename STYPE, typename CTYPE>
  void compute_scaled_coord
  (const VTYPE iv, const DTYPE dimension, const ATYPE * axis_size, 
   const STYPE * scale, CTYPE * scaled_coord)
  {
    compute_coord(iv, dimension, axis_size, scaled_coord);
    for (DTYPE d = 0; d < dimension; d++)
      { scaled_coord[d] *= scale[d]; }
  }


  /// Return index of vertex with specified coord.
  template <typename VTYPE, typename GTYPE, typename DTYPE, typename ATYPE>
  VTYPE compute_vertex_index
  (const GTYPE * coord, const DTYPE dimension, const ATYPE * axis_size)
  {
    VTYPE iv = 0;
    VTYPE inc = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      iv += inc*coord[d];
      inc = inc*axis_size[d];
    }

    return(iv);
  }


  /*!
   *  @brief Return number of vertices in grid or subgrid.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param facet_list List of facets determining subgrid.
   *    - If empty, return number of vertices in entire grid.
   *    - Otherwise, return number of vertices contained in the intersection of all facets in facet_list.
   *  @param[out] num_vertices Number of vertices.
   */
  template <typename DTYPE, typename ATYPE, typename FTYPE, typename NTYPE>
  void compute_num_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const FTYPE & facet_list,
   NTYPE & num_vertices)
  {
    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) 
        { num_vertices = num_vertices * axis_size[d]; }
    }
  }


  /*!
   *  @brief Return number of subsampled vertices along axis
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param subsample_period Subsample period. 
   *    - Only count every k'th vertex along axis 
   *      where k = subsample_period.
   *    @pre subsample_period is a positive integer.
   */
  template <typename ATYPE, typename PTYPE>
  inline ATYPE compute_num_subsampled_vertices_along_axis
  (const ATYPE axis_size, const PTYPE subsample_period)
  { 
    return(integer_divide(axis_size+subsample_period-1, subsample_period)); 
  }

  
  /*!
   *  @brief Return number of vertices in subsampled grid or grid subspace.
   *  @param dimension Grid dimension.
   *  axis_size[d] Number of vertices along grid axis d.
   *  subsample_period[d] Only count every k'th vertex along axis d
   *    where k = subsample_period[d]
   *  facet_list Ignore facets in facet_list
   *  num_vertices Number of vertices.
   */
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename FTYPE, 
            typename NTYPE>
  void compute_num_subsampled_vertices
  (const DTYPE dimension, const ATYPE * axis_size,
   const PTYPE subsample_period, const FTYPE & facet_list,
   NTYPE & num_vertices)

  {
    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {

        if (subsample_period[d] < 1) 
          { throw_subsample_period_error
              ("compute_num_subsampled_vertices", subsample_period[d]); }

        ATYPE num_subsampled_vertices_along_axis =
          compute_num_subsampled_vertices_along_axis
          (axis_size[d], subsample_period[d]);
        num_vertices = num_vertices * num_subsampled_vertices_along_axis;
      }
    }
  }


  /// @brief Return number of grid vertices
  /// - NOTE: REPLACE BY FASTER VERSION WHICH DOES NOT USE FACET_LIST0.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
  { 
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_vertices(dimension, axis_size, facet_list0, num_vertices);
  }


  /*!
   *  @brief Return number of grid vertices between two vertices.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param iv0 Vertex index.
   *  @param iv1 Vertex index.
   *  @param[out] num_vertices = Number of grid vertices between \a iv0 and \a iv1.
   *  @pre 0 <= iv0 <= iv1 < (total number of grid vertices)
   */
  template <typename VTYPE, typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const VTYPE iv1, NTYPE & num_vertices)
  {
    VTYPE coord0[dimension];
    VTYPE coord1[dimension];
    IJK::PROCEDURE_ERROR error("compute_num_grid_vertices");

    if (!check_range(dimension, axis_size, iv0, iv1, error)) { throw error; };

    compute_coord(iv0, dimension, axis_size, coord0);
    compute_coord(iv1, dimension, axis_size, coord1);

    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (coord0[d] > coord1[d]) {
        error.AddMessage("Programming error in calculating ", d,
                         "'th coordinate.");
        error.AddMessage("  coord0 (", coord0[d], 
                         ") > coord1 (", coord1[d], ").");
        throw error;
      }

      num_vertices = num_vertices * (coord1[d]-coord0[d]+1);
    }
  }


  /// Return number of vertices in grid interior.
  template <typename DTYPE, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_num_interior_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const WTYPE boundary_width,
   NTYPE & num_vertices)
  {
    IJK::ARRAY<ATYPE> interior_axis_size(dimension);
    num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 2*boundary_width) { return; };
      interior_axis_size[d] = axis_size[d]-2*boundary_width;
    }
    compute_num_grid_vertices(dimension, interior_axis_size.Ptr(), num_vertices);
  }


  /// Return number of vertices in grid interior for boundary width 1.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_interior_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
  {
    compute_num_interior_grid_vertices
      (dimension, axis_size, 1, num_vertices);
  }


  /// Return number of vertices in grid boundary.
  template <typename DTYPE, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_num_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const WTYPE boundary_width, NTYPE & num_boundary_vertices)
  {
    NTYPE num_grid_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
    NTYPE num_interior_vertices;
    compute_num_interior_grid_vertices
      (dimension, axis_size, boundary_width, num_interior_vertices);
    num_boundary_vertices = num_grid_vertices - num_interior_vertices;
  }


  /// Return number of vertices in grid boundary.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_boundary_vertices)
  {
    compute_num_boundary_grid_vertices
      (dimension, axis_size, 1, num_boundary_vertices);
  }


  /// Return number of cubes in grid or grid subspace
  template <typename DTYPE, typename ATYPE, typename FTYPE, typename NTYPE>
  void compute_num_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   const FTYPE & facet_list, NTYPE & num_cubes)
    // facet_list = ignore facets in facet_list
  {
    if (dimension < 1) {
      num_cubes = 0;
      return;
    }

    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {
        if (axis_size[d] < 2) { 
          num_cubes = 0;
          return; 
        };
        num_cubes = num_cubes * (axis_size[d]-1);
      }
    }
  }


  /// Return number of grid cubes
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_cubes)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_cubes(dimension, axis_size, facet_list0, num_cubes);
  }


  /// Return number of grid edges
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_edges
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_edges)
  {
    num_edges = 0;
    for (DTYPE orth_dir = 0; orth_dir < dimension; orth_dir++) {

      if (axis_size[orth_dir] >= 2) {
        NTYPE num_vertices_in_facet;
        compute_num_vertices_in_grid_facet
          (dimension, axis_size, orth_dir, num_vertices_in_facet);

        num_edges += (num_vertices_in_facet*(axis_size[orth_dir]-1));
      }
    }
  }


  /// Return number of cubes in grid interior
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_interior_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_interior_cubes)
  {
    if (dimension <= 0) {
      num_interior_cubes = 0;
      return;
    }

    num_interior_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 3) { 
        num_interior_cubes = 0;
        return; 
      };
      num_interior_cubes = num_interior_cubes*(axis_size[d]-3);
    }
  }


  /// Return number of cubes in grid boundary
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_boundary_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_boundary_cubes)
  {
    NTYPE num_cubes, num_interior_cubes;

    compute_num_grid_cubes(dimension, axis_size, num_cubes);
    compute_num_interior_grid_cubes
      (dimension, axis_size, num_interior_cubes);

    num_boundary_cubes = num_cubes - num_interior_cubes;
  }

  /*!
   *  @brief Return number of inner vertices in grid or grid subspace.
   *  - Inner vertices are vertices which are not on 
   *    an outer facet of the grid.
   *  - Outer facets are grid facets which do not contain the origin.
   *  @param dimension Grid dimension
   *  @param axis_size[d] Number of vertices along grid axis d.
   */
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_inner_vertices
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_vertices)
  {
    if (dimension < 1) { 
      num_vertices = 0;
      return; 
    };

    num_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE asize = axis_size[d];
      if (asize < 2) { 
        num_vertices = 0;
        return; 
      }
      else { num_vertices = num_vertices * (asize-1);  }
    }
  }


  /*!
   *  @brief Return number of outer vertices in grid or grid subspace.
   *  - Outer vertices are vertices which are are on an outer facet of the grid.
   *  - Outer facets are grid facets which do not contain the origin.
   *  @param dimension Grid dimension
   *  @param axis_size[d] Number of vertices along grid axis d.
   */
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_outer_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   NTYPE & num_outer_vertices)

  {
    NTYPE num_grid_vertices;
    NTYPE num_inner_vertices;
    compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
    compute_num_inner_vertices(dimension, axis_size, num_inner_vertices);
    num_outer_vertices = num_grid_vertices-num_inner_vertices;
  }


  /*!
   *  @brief Return number of grid vertices in a region 
   *     whose edges all have the same length.
   *  @param dimension  Dimension of grid.
   *  @param region_edge_length  Number of grid edges contained in each region edge.
   *  @param[out] num_region_vertices  Number of grid vertices in a region.
   */
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE region_edge_length,
   NTYPE & num_region_vertices)
  {
    num_region_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      num_region_vertices = num_region_vertices * (region_edge_length+1);
    }
  }


  /*!
   *  Return number of grid cubes in a region whose edges all have the same length.
   *  @param dimension  Dimension of grid.
   *  @param region_edge_length  Number of grid edges contained in each region edge.
   *  @param[out] num_region_cubes  Number of grid cubes in a region.
   */
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_grid_cubes_in_region
  (const DTYPE dimension, const ATYPE region_edge_length,
   NTYPE & num_region_cubes)
  {
    if (region_edge_length < 1) { 
      num_region_cubes = 0; 
      return;
    };

    num_region_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      num_region_cubes = num_region_cubes * region_edge_length;
    }
  }


  /*!
   *  \brief Return number of vertices in a region.
   *  - Regions completely contained in the interior of the grid all have the
   *    same number of vertices.
   *  - Regions bounded by the grid boundary will have fewer vertices.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param iv0  Primary vertex (lowest x,y,z,... coordinates) in region.
   *  @param max_region_edge_length  Maximum number of grid edges contained in each region edge.
   *  @param[out] num_region_vertices  Number of vertices in region.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  void compute_num_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const ATYPE max_region_edge_length,
   NTYPE & num_region_vertices)
  {
    ATYPE coord[dimension];

    compute_coord(iv0, dimension, axis_size, coord);

    num_region_vertices = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE numv_along_axis = max_region_edge_length + 1;
      if (coord[d] + max_region_edge_length >= axis_size[d]) {
        if (coord[d] < axis_size[d]) 
          { numv_along_axis = axis_size[d] - coord[d]; }
        else
          { numv_along_axis = 0; };
      }
      num_region_vertices = num_region_vertices * numv_along_axis;
    }
  }


  /*!
   *  @brief Return number of regions along a single axis.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param region_edge_length Number of grid edges contained in each region  * e.
   *  @pre \a region_edge_length > 0.
   */
  template <typename ATYPE>
  ATYPE compute_num_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    ATYPE num_regions_along_axis = 
      long(axis_size+region_edge_length-2)/long(region_edge_length);
    return(num_regions_along_axis);
  }


  /*!
   *  @brief Return total number of regions in grid or subgrid.
   *  @param dimension  Dimension of grid.
   *  @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
   *  @param region_edge_length  Number of grid edges contained in each region edge.
   *  @param facet_list List of facets determining subgrid.
   *    - If empty, return number of regions in entire grid.
   *    - Otherwise, return number of regions in subgrid 
   *      formed by intersection of all facets in facet_list.
   *  @param[out] num_regions Number of regions.
   *  @pre \a region_edge_length > 0.
   */
  template <typename DTYPE, typename ATYPE, typename FTYPE, typename NTYPE>
  void compute_num_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_regions)
  {
    IJK::PROCEDURE_ERROR error("compute_num_regions");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    num_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {

      if (!facet_list.Contains(d)) {
        if (axis_size[d] <= 1) { 
          num_regions = 0;
          return; 
        };

        ATYPE num_regions_along_axis = 
          compute_num_regions_along_axis(axis_size[d], region_edge_length);
        num_regions = num_regions * num_regions_along_axis; 
      }
    }
  }


  /*!
   *  @brief Return total number of regions.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d]  Number of grid vertices along axis \a d.
   *  @param region_edge_length  Number of grid edges 
   *    contained in each region edge.
   *  @param[out] num_regions Number of regions.
   *  @pre \a region_edge_length > 0.
   */
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, NTYPE & num_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_regions
      (dimension, axis_size, region_edge_length, facet_list0, num_regions);
  }


  /*!
   *  @brief Return number of full regions along a single axis.
   *  @param axis_size[d]  Number of grid vertices along axis \a d.
   *  @param region_edge_length  Number of grid edges 
   *    contained in each region edge.
   *  @pre \a region_edge_length > 0.
   */
  template <typename ATYPE>
  ATYPE compute_num_full_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    if (axis_size < 1) { return(0); };
    ATYPE num_full_regions_along_axis = 
      long(axis_size-1)/long(region_edge_length);
    return(num_full_regions_along_axis);
  }


  /*!
   *  @brief Return number of full regions in grid or subgrid.
   *  @param dimension Dimension of grid.
   *  @param axis_size[d]  Number of grid vertices along axis \a d.
   *  @param region_edge_length Number of grid edges contained in each region edge.
   *  @param facet_list List of facets determining subgrid.
   *    - If empty, return number of vertices in entire grid.
   *    - Otherwise, return number of vertices in subgrid 
   *      formed by intersection of all facets in facet_list.
   *  @param[out] num_full_regions Number of full regions in grid or subgrid.
   *  @pre \a region_edge_length > 0.
   */
  template <typename DTYPE, typename ATYPE, typename FTYPE, typename NTYPE>
  void compute_num_full_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_full_regions)
  {
    IJK::PROCEDURE_ERROR error("compute_num_full_regions");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    num_full_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (!facet_list.Contains(d)) {
        NTYPE k = compute_num_full_regions_along_axis
          (axis_size[d], region_edge_length);
        num_full_regions = num_full_regions*k;
      }
    }
  }


  /// Return number of full regions in grid.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_full_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, NTYPE & num_full_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list0, num_full_regions);
  }


  /*!
   *  @brief Return number of partial regions along a single axis (0 or 1.)
   *  @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
   *  @param region_edge_length = Number of grid edges contained in each region edge.
   *  @pre \a region_edge_length > 0.
   */
  template <typename ATYPE>
  ATYPE compute_num_partial_regions_along_axis
  (const ATYPE axis_size, const ATYPE region_edge_length)
  {
    if (axis_size < 1) { return(0); };
    if (long(axis_size-1)%long(region_edge_length) == 0) { return(0); };
    return(1);
  }


  /// Return number of partial regions in grid or subgrid.
  template <typename DTYPE, typename ATYPE, 
            typename FTYPE, typename NTYPE>
  void compute_num_partial_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const FTYPE & facet_list,
   NTYPE & num_partial_regions)
  {
    NTYPE num_regions;
    compute_num_regions
      (dimension, axis_size, region_edge_length, 
       facet_list, num_regions);
    NTYPE num_full_regions;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list, num_full_regions);
    num_partial_regions = num_regions-num_full_regions;
  }


  /// Return number of partial regions in grid.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_partial_regions
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, 
   NTYPE & num_partial_regions)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, 
       facet_list0, num_partial_regions);
  }

  /*!
   *  @brief Return total number of regions in grid facet.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param region_edge_length = Number of grid edges contained in each region edge.
   *  @param orth_dir = Direction orthogonal to facet.
   */
  template <typename NTYPE, typename DTYPE, typename ATYPE>
  NTYPE compute_num_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE * region_edge_length, const DTYPE orth_dir)
  {
    NTYPE num_regions = 1;
    for (DTYPE d = 0; d < dimension; d++) {

      if (d != orth_dir) {
        if (axis_size[d] <= 1) { return(0); };

        ATYPE num_regions_along_axis = 
          compute_num_regions_along_axis(axis_size[d], region_edge_length[d]);
        num_regions = num_regions * num_regions_along_axis; 
      }
    }
    return(num_regions);
  }


  /*!
   *  @brief Return total number of full regions in grid facet.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param region_edge_length  Number of grid edges contained in each region edge.
   *  @param orth_dir  Direction orthogonal to facet.
   *  @param[out] num_full_regions  Number of full regions in grid or subgrid.
   */
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_full_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const DTYPE orth_dir,
   NTYPE & num_full_regions)
  {
    FACET_LIST1<DTYPE> facet_list1(orth_dir);
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, 
       facet_list1, num_full_regions);
  }


  /// Return total number of partial regions in grid facet.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_partial_regions_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const DTYPE orth_dir,
   NTYPE & num_partial_regions)
  {
    FACET_LIST1<DTYPE> facet_list1(orth_dir);
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, 
       facet_list1, num_partial_regions);
  }


  /*!
   *  @brief Return number of vertices along subsampled axis.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param subsample_period Only count every k'th vertex along axis 
   *    where k = \a subsample period.
   *  @pre \a subsample_period > 0.
   */
  template <typename ATYPE, typename PTYPE>
  ATYPE compute_subsample_size
  (const ATYPE axis_size, const PTYPE subsample_period)
  {
    ATYPE subsample_size =
      integer_divide(axis_size+subsample_period-1, subsample_period);
    return(subsample_size);
  }


  /*!
   *  @brief Compute axis size for each axis in subsampled grid.
   *  @param dimension Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param subsample_period Subsample period, i.e., how often grid
   *     is subsampled along each axis.<br>
   *     For instance:
   *    - If subsample_period = 2, every other vertex is subsampled
   *      along each axis.
   *    - If subsample_period = 3, every third vertex is subsampled
   *      along each axis.
   *  @param[out] subsampled_axis_size Array:
   *      <em>subsampled_axis_size[d]</em> =
   *         Number of vertices along subsampled axis \a d.
   *  @pre Array subsampled_axis_size[] is preallocated
   *         to size at least \a dimension.
   */
  template <typename DTYPE, typename ATYPE, typename PTYPE, 
            typename ATYPE2>
  void compute_subsample_axis_sizes
  (const DTYPE dimension, const ATYPE axis_size[], 
   const PTYPE subsample_period, ATYPE2 subsampled_axis_size[])
  {
    for (DTYPE d = 0; d < dimension; d++) {
      subsampled_axis_size[d] = 
        compute_subsample_size(axis_size[d], subsample_period); 
    }
  }


  /*!
   *  @brief Return number of vertices in subsampled grid.
   *  @param dimension  Dimension of grid.
   *  @param axis_size  
   *    Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
   *  @param subsample_period
   *    Array: <em>subsample_period[d]</em> = subsample period along axis \a d.
   *  @param[out] num_vertices Number of vertices.
   */
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename NTYPE>
  void compute_subsample_size
  (const DTYPE dimension, const ATYPE * axis_size, 
   const PTYPE subsample_period, NTYPE & num_vertices)
  {
    FACET_LIST0<DTYPE> facet_list0;
    compute_num_subsampled_vertices
      (dimension, axis_size, subsample_period, facet_list0, num_vertices);
  }


  /*!
   *  @brief Return number of vertices along supersampled axis.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param supersample_period = Only count every k'th vertex 
   *    along axis where k = \a supersample period.
   *  @pre \a supersample_period > 0.
   */
  template <typename ATYPE, typename PTYPE>
  ATYPE compute_supersample_size
  (const ATYPE axis_size, const PTYPE supersample_period)
  {
    if (axis_size < 1) { return(0); };
    ATYPE supersample_size = (axis_size-1)*supersample_period+1;
    return(supersample_size);
  }

  ///@}


  // *****************************************************************
  //! @name COMPUTE BOUNDARY BITS
  // *****************************************************************

  ///@{

  /*!
   *  @brief Compute boundary bits for vertex \a iv.
   *  @param iv Vertex index.
   *  @param dimension Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param [out] boundary_bits Bits flagging boundaries containing vertex \a iv.
   *    - If bit \a 2d is true, then <em>d</em>'th coordinate 
   *      of vertex \a iv is zero.
   *    - If bit <em>(2d+1)</em> is true, then <em>d</em>'th coordinate 
   *      of vertex \a iv equals <em>axis_size[d]-1</em>.
   *  @pre \li Variable \a boundary_bits has at least 
   *     <em>(2*dimension)</em> bits.
   *  @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
   */
  template <typename VTYPE, typename DTYPE, typename ATYPE, typename BTYPE>
  void compute_boundary_bits
  (const VTYPE iv, const DTYPE dimension,
   const ATYPE * axis_size, BTYPE & boundary_bits)
  {
    VTYPE k = iv;
    BTYPE flag = 1;
    boundary_bits = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE c = k % axis_size[d];
      k = k / axis_size[d];

      if (c == 0) { boundary_bits = boundary_bits | flag; };
      flag = (flag << 1);
      if (c+1 >= axis_size[d]) { boundary_bits = boundary_bits | flag; };
      flag = (flag << 1);
    };
  }


  /*!
   *  @brief Compute boundary bits for cube \a iv.
   *  @param iv Cube index.
   *  @param dimension Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param [out] boundary_bits 
   *    Bits flagging boundaries containing cube \a iv.
   *    - If bit \a 2d is true, then <em>d</em>'th coordinate
   *      of cube \a icube is zero.
   *    - If bit <em>(2d+1)</em> is true, then <em>d</em>'th coordinate
   *      of cube \a icube equals <em>axis_size[d]-2</em>.
   *  @pre \li Variable \a boundary_bits has at least <em>(2*dimension)</em> bits.
   *  @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
   */
  template <typename VTYPE, typename DTYPE, typename ATYPE, typename BTYPE>
  void compute_boundary_cube_bits
  (const VTYPE iv, const DTYPE dimension,
   const ATYPE * axis_size, BTYPE & boundary_bits)
  {
    VTYPE k = iv;
    BTYPE flag = 1;
    boundary_bits = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE c = k % axis_size[d];
      k = k / axis_size[d];

      if (c == 0) { boundary_bits = boundary_bits | flag; };
      flag = (flag << 1);
      if (c+2 >= axis_size[d]) { boundary_bits = boundary_bits | flag; };
      flag = (flag << 1);
    };
  }


  /*!
   *  @brief Return true if cube is on grid boundary.
   *  @param icube Cube index.
   *  @param dimension Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @pre \li axis_size[d] > 0 for all \a d = 0,..., \a dimension-1.
   */
  template <typename VTYPE, typename DTYPE, typename ATYPE>
  inline bool is_cube_on_grid_boundary
  (const VTYPE icube, const DTYPE dimension, const ATYPE * axis_size)
  {
    VTYPE k = icube;
    for (DTYPE d = 0; d < dimension; d++) {
      const ATYPE c = k % axis_size[d];
      k = k / axis_size[d];

      if (c == 0) { return true; }
      if (c+2 >= axis_size[d]) { return true; }
    };

    return false;
  }

  
  /// Return true if edge is on grid boundary.
  template <typename VTYPE, typename DIR_TYPE,
            typename DTYPE, typename ATYPE>
  bool is_edge_on_grid_boundary
  (const VTYPE iend0, const DIR_TYPE edge_dir, 
   const DTYPE dimension, const ATYPE * axis_size)
  {
    VTYPE k = iend0;
    for (DTYPE d = 0; d < dimension; d++) {
      ATYPE c = k % axis_size[d];
      k = k / axis_size[d];

      if (c == 0) {
        if (edge_dir != d)
          { return(true); }
      }
      if (c+1 >= axis_size[d]) { 
        if (edge_dir != d)
          { return(true); }
      }
    }

    return(false);
  }

  ///@}


  // *****************************************************************
  //! @name COMPUTE REGION
  // *****************************************************************

  ///@{

  /// Compute region within boundary around given cube.
  template <typename VTYPE0, typename DTYPE, typename ATYPE0,
            typename DIST_TYPE, typename VTYPE1, typename ATYPE1>
  void compute_region_around_cube
  (const VTYPE0 icube, const DTYPE dimension, const ATYPE0 * axis_size,
   const DIST_TYPE dist2cube, VTYPE1 & region_iv0, ATYPE1 * region_axis_size)
  {
    IJK::ARRAY<ATYPE0> cube_coord(dimension);
    IJK::ARRAY<ATYPE0> region_iv0_coord(dimension);

    compute_coord(icube, dimension, axis_size, cube_coord.Ptr());

    for (DTYPE d = 0; d < dimension; d++) {
      if (cube_coord[d] > dist2cube) 
        { region_iv0_coord[d] = cube_coord[d]-dist2cube; }
      else
        { region_iv0_coord[d] = 0; }

      ATYPE0 c = cube_coord[d]+dist2cube+2;
      if (c <= axis_size[d]) 
        { region_axis_size[d] = c-region_iv0_coord[d]; }
      else 
        { region_axis_size[d] = axis_size[d]-region_iv0_coord[d]; }
    }

    region_iv0 = 
      compute_vertex_index<VTYPE1>
      (region_iv0_coord.PtrConst(), dimension, axis_size);
  }

  /// Compute region within boundary around given vertex coord
  template <typename CTYPE0, typename DTYPE, typename ATYPE0, 
            typename DIST_TYPE, typename CTYPE1, typename ATYPE1>
  void compute_region_around_vertex_coord
  (const CTYPE0 * vertex_coord, 
   const DTYPE dimension, const ATYPE0 * axis_size, 
   const DIST_TYPE dist2vertex, 
   CTYPE1 * region_iv0_coord, ATYPE1 * region_axis_size)
  {
    for (DTYPE d = 0; d < dimension; d++) {
      if (vertex_coord[d] > dist2vertex) 
        { region_iv0_coord[d] = vertex_coord[d]-dist2vertex; }
      else
        { region_iv0_coord[d] = 0; }

      ATYPE0 c = vertex_coord[d]+dist2vertex;
      if (c < axis_size[d]) 
        { region_axis_size[d] = c-region_iv0_coord[d]+1; }
      else 
        { region_axis_size[d] = axis_size[d]-region_iv0_coord[d]; }
    }

  }

  ///@}


  // **************************************************
  //! @name TEMPLATES TO CHECK VALUES AND ARRAYS.
  // **************************************************

  ///@{

  /// @brief Return true if dimension is non-negative.
  /// - Return false and set error message if dimension is negative.
  template <typename DTYPE>
  bool check_dimension(const DTYPE dimension, IJK::ERROR & error)
  {
    if (dimension < 0) {
      error.AddMessage("Illegal dimension ", dimension, ".");
      error.AddMessage("Dimension must be non-negative.");
      return false;
    }

    return true;
  }


  /// @brief Return true if 0 <= iv0 <= iv1 < total_num_grid_vertices.
  /// - Return false and set error message if conditions are violated.
  template <typename VTYPE, typename DTYPE, typename ATYPE>
  bool check_range
  (const DTYPE dimension, const ATYPE * axis_size,
   const VTYPE iv0, const VTYPE iv1, IJK::ERROR & error)
  {
    VTYPE total_num_grid_vertices;
    compute_num_grid_vertices(dimension, axis_size, total_num_grid_vertices);

    if (iv0 > iv1) {
      error.AddMessage("Illegal vertex range. Vertex index ", iv0, 
                       " is greater than vertex index ", iv1, ".");
      return(false);
    }

    if (0 > iv0 || iv1 >= total_num_grid_vertices) {
      error.AddMessage("Illegal vertex indices: ", iv0, ", ", iv1, ".");
      error.AddMessage("Vertex indices should be in range: [0,",
                       total_num_grid_vertices, "].");
      return(false);
    }

    return(true);
  }


  /// @brief Return true if (coord0[d] <= coord1[d]) for all d.
  /// - Return false and set error message if (coord0[d] > coord1[d]
  ///   for some d.
  template <typename DTYPE, typename VTYPE, typename CTYPE>
  bool check_region_coordinates
  (const DTYPE dimension, const VTYPE iv0, const CTYPE * coord0, 
   const VTYPE iv1, const CTYPE * coord1, IJK::ERROR & error)

  {
    for (DTYPE d = 0; d < dimension; ++d) {
      if (coord0[d] > coord1[d]) {
        error.AddMessage("Illegal coordinates.  Coordinates of vertex 0 should be less than or equal to coordinates of vertex 1.");
        error.AddMessage(" Vertex 0 = ", iv0, 
                         ".  Coordinate ", d , " = ", coord0[d], ".");
        error.AddMessage(" Vertex 1 = ", iv1, 
                         ".  Coordinate ", d, " = ", coord1[d], ".");
        return(false);
      }
    }

    return(true);
  }

  /// Check that array vertex_list[] is not NULL.
  template <typename VTYPE>
  bool check_vertex_list(const VTYPE * vertex_list, IJK::ERROR & error)
  {
    if (vertex_list == NULL) {
      error.AddMessage("Vertex list is NULL");
      return(false);
    }
    else { return(true); }
  }


  /// Check that region edge length is positive.
  template <typename LTYPE>
  bool check_region_edge_length(const LTYPE length, IJK::ERROR & error)
  {
    if (length <= 0) {
      error.AddMessage("Region edge length must be positive.");
      return(false);
    }
    return(true);
  }

  /// Check number of vertices added to vertex list.
  template <typename NTYPE0, typename NTYPE1>
  bool check_num_vertices_added
  (const NTYPE0 num_added, const NTYPE1 numv, IJK::ERROR & error)
  {
    if (num_added != numv) {
      error.AddMessage("Added ", num_added, " vertices to vertex list.");
      error.AddMessage("Number of vertices in list should be ", numv, ".");
      return(false);
    }
    return(true);
  }

  /// Check number of vertices equals number of grid vertices
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  bool check_num_vertices
  (const DTYPE dimension, const ATYPE * axis_size, const NTYPE numv,
   IJK::ERROR & error)
  {
    NTYPE numv2;
    compute_num_grid_vertices(dimension, axis_size, numv2);
    if (numv != numv2) {
      error.AddMessage("Programming error. Incorrect number of vertices.");
      return(false);
    }
    return(true);
  }

  /// Check that DIFFTYPE is signed and has range [-n:n]
  template <typename DIFFTYPE, typename NTYPE>
  bool check_difftype(const NTYPE n, IJK::ERROR & error)
  {
    if (std::numeric_limits<DIFFTYPE>::min() >= 0) {
      error.AddMessage("Programming error. Template parameter DIFFTYPE must be signed.");
      return(false);
    }
    else if (n > std::numeric_limits<DIFFTYPE>::max()) {
      error.AddMessage("Error. Template parameter DIFFTYPE is not large enough to store integer ",
                       n, ".");
      return(false);
    }
    else if (-n < std::numeric_limits<DIFFTYPE>::min()) {
      error.AddMessage("Error. Template parameter DIFFTYPE is not large enough to store integer ",
                       -n, ".");
      return(false);
    }

    return(true);
  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: COMPUTING INCREMENTS
  // *****************************************************************

  ///@{

  // Forward declaration.
  template <typename DTYPE, typename ATYPE, typename ITYPE, typename VTYPE,
            typename NTYPE>
  void get_subgrid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ITYPE * axis_increment,
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   VTYPE * vlist, NTYPE & num_subgrid_vertices);

  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_subgrid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   VTYPE * vlist);


  /*!
   *  @brief Compute increment to add to index of a vertex to get next vertices
   *    along the axes.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param[out] increment[] = Axis increment.
   *     iv+increment[d] is the index of the vertex after iv along axis \a d.
   *  @pre Array increment[] is pre-allocated to size at least \a dimension.
   */
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_axis_increment
  (const DTYPE dimension, const ATYPE * axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_increment");

    if (dimension <= 0) { return; };

    if (axis_size == NULL || increment == NULL) {
      error.AddMessage
        ("Programming error. axis_size == NULL or increment == NULL.");
      throw error;
    }

    increment[0] = 1;
    for (DTYPE d = 1; d < dimension; d++)
      { increment[d] = increment[d-1]*axis_size[d-1]; }
  }

  /// DEPRECATED
  /// OLD function name.
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_increment
  (const DTYPE dimension, const ATYPE * axis_size, ITYPE * increment)
  {
    compute_axis_increment(dimension, axis_size, increment);
  }


  /*!
   *  @brief Compute increment to add to index of vertex 
   *    to get next subsampled vertex along each axis.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param subsample_period  
   *    Array: <em>subsample_period[d]</em> = Subsample period along axis \a d.
   *  @param[out] increment  Array: <em>increment[d]</em> = Increment to add to index of a vertex to get next subsampled vertex along axis \a d.
   *  @pre Array increment[] is pre-allocated to size at least \a dimension.
   */
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename ITYPE>
  void compute_subsample_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const PTYPE subsample_period, ITYPE * increment)
  {
    compute_increment(dimension, axis_size, increment);

    for (DTYPE d = 0; d < dimension; ++d)
      { increment[d] *= subsample_period[d]; };
  }


  /*!
   *  @brief Compute increment to add to index of primary vertex to get
   *    k'th corner vertex of region.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param region_axis_size  Array: <em>region_axis_size[d]</em> = Number of vertices along region axis \a d.
   *  @param[out] increment  Array: Region corner increment. iv0+increment[k] is k'th corner vertex of region.
   *  @pre Array increment[] is allocated with size at least number of corner regions.
   */
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_region_corner_increment
  (const DTYPE dimension, const ATYPE * grid_axis_size, 
   const ATYPE * region_axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_region_corner_increment");
    IJK::ARRAY<ITYPE> axis_increment(dimension);

    if (dimension <= 0) { return; };

    if (grid_axis_size == NULL) {
      error.AddMessage("Programming error. grid_axis_size == NULL.");
      throw error;
    }
    if (region_axis_size == NULL) {
      error.AddMessage("Programming error. region_axis_size == NULL.");
      throw error;
    }    
    if (increment == NULL) {
      error.AddMessage("Programming error. increment == NULL.");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      if (region_axis_size[d] < 1) {
        error.AddMessage
          ("Programming error.  Region axis size must be at least 1.");
        throw error;
      }
    }

    compute_increment(dimension, grid_axis_size, axis_increment.Ptr());
    ITYPE num_region_corners = compute_num_cube_vertices(dimension);

    for (ITYPE j = 0; j < num_region_corners; j++) {
      increment[j] = 0;
      ITYPE j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        if ((j0 % 2) == 1) {
          increment[j] = 
            increment[j] + (region_axis_size[d]-1)*axis_increment[d];
        };
        j0 = j0/2;
      };
    }
  }


  /*!
   *  @brief Compute increment to add to index of primary vertex to get
   *    k'th corner vertex of cubic region (all axes have same size.)
   *  @param dimension  Dimension of grid.
   *  @param grid_axis_size[d] Number of grid vertices along axis d.
   *  @param region_axis_size[d] Number of region vertices along any axis d.
   *  @param[out] increment Region corner increment. 
   *    - iv0+increment[k] is the k'th corner vertex of the region.
   *    @pre Array increment[] is allocated with size at least number 
   *       of corner regions.
   */
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_cubic_region_corner_increment
  (const DTYPE dimension, const ATYPE * grid_axis_size, 
   const ATYPE region_axis_size, ITYPE * increment)
  {
    IJK::PROCEDURE_ERROR error("compute_cubic_region_corner_increment");
    IJK::ARRAY<ITYPE> axis_increment(dimension);

    if (dimension <= 0) { return; };

    if (grid_axis_size == NULL) {
      error.AddMessage("Programming error. grid_axis_size == NULL.");
      throw error;
    }
    if (increment == NULL) {
      error.AddMessage("Programming error. increment == NULL.");
      throw error;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      if (region_axis_size < 1) {
        error.AddMessage("Programming error.  Region axis size must be at least 1.");
        throw error;
      }
    }

    compute_increment(dimension, grid_axis_size, axis_increment.Ptr());
    ITYPE num_region_corners = compute_num_cube_vertices(dimension);

    for (ITYPE j = 0; j < num_region_corners; j++) {
      increment[j] = 0;
      ITYPE j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        if ((j0 % 2) == 1) {
          increment[j] = increment[j] + (region_axis_size-1)*axis_increment[d];
        };
        j0 = j0/2;
      };
    }
  }


  /// @brief Compute increment to add to index of primary vertex
  ///   to get k'th vertex in region.
  /// @pre Array increment is allocated with size at least
  ///   number of region vertices.
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_region_vertex_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, const ATYPE scale, ITYPE * increment)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(scale);

    for (DTYPE d = 0; d < dimension; d++) {
      subgrid_axis_size[d] = region_edge_length*scale+1; 

      if (subgrid_axis_size[d] > axis_size[d]) {
        IJK::PROCEDURE_ERROR error("compute_region_vertex_increment");
        error.AddMessage("Region is larger than grid.");
        error.AddMessage("  Grid axis size[", d, "] = ", axis_size[d], ".");
        error.AddMessage("  Region size[", d, "] = ", subgrid_axis_size[d],
                         ".");
        throw error;
      }
    }

    subsample_subgrid_vertices
      (dimension, axis_size, ITYPE(0), subgrid_axis_size.PtrConst(), 
       subsample_period, increment);
  }


  /*!
   *  @brief Compute increment to add to index of primary vertex of region to get
   *    primary vertex of k'th grid cube in cubic region.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d]  Number of vertices along axis \a d.
   *  @param region_edge_length[d] 
   *    Number of grid edges along edge \a d of the cubic region.
   *  @param increment Region increment. 
   *       - iv0+increment[k] is the k'th cube of the region.
   *  @pre Array increment[] is allocated with size at least number 
   *       of grid cubes contained in the region.
   */
  template <typename DTYPE, typename ATYPE, typename ITYPE>
  void compute_region_cube_increment
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, ITYPE * increment)
  {
    IJK::ARRAY<ATYPE> subgrid_size(dimension);

    if (region_edge_length < 1)  // No cubes.
      { return; }

    for (DTYPE d = 0; d < dimension; d++) 
      { subgrid_size[d] = region_edge_length; }

    get_subgrid_vertices
      (dimension, axis_size, 0, subgrid_size.PtrConst(), increment);
  }


  /*!
   *  @brief Compute increment to add to vertex 0 to compute vertex i 
   *    of hypercube.
   *  @param dimension  Dimension of grid.
   *  @param axis_increment[] Axis increment.
   *    - iv+axis_increment[d] is the next vertex after vertex iv along axis d.
   *  @param[out] cube_vertex_increment[] Cube vertex increment.
   *    - iv0+cube_vertex_increment[i] is the i'th vertex of the hypercube
   *    with primary vertex iv0.
   *  @pre Array cube_vertex_increment[] is allocated with size 
   *    at least number of cube vertices.
   */
  template <typename DTYPE, typename ITYPE1, typename ITYPE2>
  void compute_cube_vertex_increment
  (const DTYPE dimension, const ITYPE1 * axis_increment, 
   ITYPE2 * cube_vertex_increment)
  {
    const ITYPE2 num_cube_vertices = compute_num_cube_vertices(dimension);
    IJK::PROCEDURE_ERROR error("compute_cube_vertex_increment");

    // Initialize cube_vertex_increment[] to 0.
    for (ITYPE2 j = 0; j < num_cube_vertices; j++)
      { cube_vertex_increment[j] = 0; }

    if (dimension < 1) { return; };

    if (axis_increment == NULL) {
      error.AddMessage("Programming error. Array axis_increment[] must be allocated and set before calling compute_cube_vertex_increment.");
      throw error;
    }

    if (cube_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array cube_vertex_increment[] must be allocated before calling compute_cube_vertex_increment.");
      throw error;
    }

    for (ITYPE2 j = 0; j < num_cube_vertices; j++) {
      cube_vertex_increment[j] = 0;
      ITYPE2 j0 = j;
      for (DTYPE d = 0; d < dimension; d++) {
        if ((j0 % 2) == 1) {
          cube_vertex_increment[j] += axis_increment[d];
        };
        j0 = j0/2;
      }
    }
  }


  /*!
   *  @brief Compute increment to add to vertex 0 to compute vertex i of facet.
   *  @param dimension  Dimension of grid.
   *  @param ifacet Facet index.
   *  @param cube_vertex_increment[] Cube vertex increment.
   *    - iv0+cube_vertex_increment[i] is the i'th vertex 
   *      of the hypercube with primary vertex iv0.
   *  @param[out] facet_vertex_increment[] Facet vertex increment.
   *    - iv0+facet_vertex_increment[i] is the i'th vertex 
   *      of the facet with primary vertex iv0.
   *  @pre Array facet_vertex_increment[] is allocated with size 
   *      at least number of facet vertices.
   */
  template <typename DTYPE, typename ITYPE>
  void compute_facet_vertex_increment
  (const DTYPE dimension, const DTYPE ifacet,
   const ITYPE * cube_vertex_increment, ITYPE * facet_vertex_increment)
  {
    IJK::PROCEDURE_ERROR error("compute_facet_vertex_increment");

    if (dimension <= 0) { return; };

    if (cube_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array cube_vertex_increment[] must be allocated and set before calling compute_facet_vertex_increment.");
      throw error;
    }

    if (facet_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array facet_vertex_increment[] must be allocated before calling compute_facet_vertex_increment.");
      throw error;
    }
    
    DTYPE orth_dir = ifacet%dimension;
    ITYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    ITYPE num_facet_vertices = compute_num_cube_facet_vertices(dimension);

    // Multiply by s mod (num_cube_vertices-1) to compute vertex index.
    ITYPE s = 2;
    for (ITYPE d = 0; d < orth_dir; d++)
      { s = s*2; };
    s = s%(num_cube_vertices-1);

    if (ifacet < dimension) {

      for (ITYPE i = 0; i < num_facet_vertices; i++) {
        ITYPE k = (i*s)%(num_cube_vertices-1);
        facet_vertex_increment[i] = cube_vertex_increment[k];
      }
    }
    else {

      // Translate by t mod num_cube_vertices to compute vertex index.
      ITYPE t = 1;
      for (ITYPE d = 0; d < orth_dir; d++)
        { t = t*2; };

      for (ITYPE i = 0; i < num_facet_vertices; i++) {
        ITYPE k = (i*s)%(num_cube_vertices-1);
        facet_vertex_increment[i] = cube_vertex_increment[k+t];
      }

      // Swap subfacets to get consistent facet orientation.
      if (num_facet_vertices > 1) {
        ITYPE num_subfacet_vertices = num_facet_vertices/2;
        for (ITYPE i = 0; i < num_subfacet_vertices; i++) {
          std::swap(facet_vertex_increment[i], 
                    facet_vertex_increment[i+num_subfacet_vertices]);
        }
      }
    }
  }


  /*!
   *  @brief Compute increment to add to vertex 0 to compute vertex i of ridge.
   *  @param dimension Dimension of grid.
   *  @param orth_dir0 Direction orthogonal to ridge.
   *  @param orth_dir1 Second direction orthogonal to ridge.
   *  @param cube_vertex_increment[] Cube vertex increment.
   *    - iv0+cube_vertex_increment[i] is the i'th vertex of the hypercube 
   *      with primary vertex iv0.
   *  @param[out] ridge_vertex_increment[]
   *    - Ridge vertex increment. iv0+ridge_vertex_increment[i] 
   *      is the i'th vertex of the ridge with primary vertex iv0.
   *  @pre Array ridge_vertex_increment[] is allocated with size 
   *      at least number of ridge vertices.
   */
  template <typename DTYPE, typename ITYPE>
  void compute_ridge_vertex_increment
  (const DTYPE dimension, const DTYPE orth_dir0, const DTYPE orth_dir1,
   const ITYPE * cube_vertex_increment, ITYPE * ridge_vertex_increment)
  {
    IJK::PROCEDURE_ERROR error("compute_ridge_vertex_increment");

    if (dimension <= 0) { return; };

    if (cube_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array cube_vertex_increment[] must be allocated and set before calling compute_ridge_vertex_increment.");
      throw error;
    }

    if (ridge_vertex_increment == NULL) {
      error.AddMessage("Programming error. Array ridge_vertex_increment[] must be allocated before calling compute_ridge_vertex_increment.");
      throw error;
    }
    
    ITYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    ITYPE num_cube_ridge_vertices = compute_num_cube_ridge_vertices(dimension);
    ITYPE mask0, mask1;

    mask0 = (ITYPE(1) << orth_dir0);
    mask1 = (ITYPE(1) << orth_dir1);

    ITYPE i = 0;
    for (ITYPE j = 0; j < num_cube_vertices; j++) {

      if (((j & mask0) == 0) && ((j & mask1) == 0)) {
        ridge_vertex_increment[i] = cube_vertex_increment[j];
        i++;
      }
    }

    if (i != num_cube_ridge_vertices) {
      error.AddMessage("Programming error.  Added ", i,
                       " values to ridge_vertex_increment[].");
      error.AddMessage("  Should have added ",
                       num_cube_ridge_vertices, " values.");
      throw error;
    }
  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: GET VERTICES
  // *****************************************************************

  ///@{

  /*!
   *  @brief Subsample vertices in subgrid.
   *  @param dimension  Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param  subgrid_origin  Subgrid origin.
   *  @param subgrid_axis_size 
   *    Array: <em>subgrid_axis_size[d]</em> = Number of vertices 
   *    along subgrid axis d.
   *  @param subsample_period 
   *    Array: <em>subsample_period[d]</em> = Only report every k'th vertex
   *    along subgrid axis \a d where k = \a subsample_period[d].
   *  @param[out] vlist[]  List of vertices.
   *  @pre \li Subgrid is contained in grid, i.e. 
   *    ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
   *  @pre \li \a subsample_period[d] is a positive integer.
   *  @pre \li Array vlist[] is preallocated to length at least number
   *    of vertices in grid or subgrid.
   */
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename VTYPE>
  void subsample_subgrid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   const PTYPE subsample_period, VTYPE * vlist)
  {
    IJK::PROCEDURE_ERROR error("subsample_subgrid_vertices");

    VTYPE num_vertices;
    compute_subsample_size
      (dimension, subgrid_axis_size, subsample_period, num_vertices);

    if (num_vertices < 1) { return; };
    // Note: subgrid_axis_size[d] >= 1 for all d

    if (!check_vertex_list(vlist, error)) { throw error; };

    IJK::ARRAY<VTYPE> subsample_increment(dimension);
    compute_subsample_increment
      (dimension, axis_size, subsample_period, subsample_increment.Ptr());

    // initialize prev_num_vertices
    vlist[0] = subgrid_origin;
    VTYPE prev_num_vertices = 1;

    // Process axes 0,1,2,..., dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 0; d < dimension; d++) {

      ATYPE num_subsampled_along_axis =
        compute_num_subsampled_vertices_along_axis
        (subgrid_axis_size[d], subsample_period[d]);

      VTYPE iv0 = subsample_increment[d];

      for (ATYPE i = 1; i < num_subsampled_along_axis; i++) {
        for (VTYPE * vprev_ptr = vlist; 
             vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
          *(vcur_ptr) = iv0 + *(vprev_ptr);
          vcur_ptr++;
        };
        iv0 = iv0 + subsample_increment[d];
      }

      prev_num_vertices = prev_num_vertices*num_subsampled_along_axis;
    }

    if (!check_num_vertices_added(prev_num_vertices, num_vertices, error))
      throw error;

  }


  /*
   *  @brief Subsample vertices in grid.
   *  @param dimension  Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param subsample_period 
   *    Array <em>subsample_period[d]</em> = Only report every k'th vertex 
   *    along subgrid axis \a d where k = \a subsample_period[d].
   *  @param[out] vlist[] = List of vertices.
   *  @pre \li \a subsample_period is a positive integer.
   *  @pre \li Array vlist[] is preallocated to length 
   *    at least number of vertices in grid or subgrid.
   */
  template <typename DTYPE, typename ATYPE, typename PTYPE, typename VTYPE>
  void subsample_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const PTYPE subsample_period, VTYPE * vlist)
  {
    subsample_subgrid_vertices(dimension, axis_size, VTYPE(0), axis_size, 
                               subsample_period, vlist);
  }


  /*!
   *  @brief Get vertices in subgrid.
   *  - Returns vertex indices sorted in increasing order.
   *  @param dimension  Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param axis_increment[] = Axis increment. iv+axis_increment[i]
   *         is the next vertex after vertex iv along axis i.
   *  @param  subgrid_origin  Subgrid origin.
   *  @param subgrid_axis_size
   *         - Array: <em>subgrid_axis_size[d]</em> =
   *           Number of vertices along subgrid axis d.
   *  @param[out] vlist[]  List of vertices.
   *  @pre \li Subgrid is contained in grid,
   *     i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
   *  @pre \li Array vlist[] is preallocated to length
   *     at least number of vertices in grid or subgrid.
   *  @pre dimension >= 1.
   */
  template <typename DTYPE, typename ATYPE, typename ITYPE, typename VTYPE,
            typename NTYPE>
  void get_subgrid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ITYPE * axis_increment,
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   VTYPE * vlist, NTYPE & num_subgrid_vertices)
  {
    IJK::PROCEDURE_ERROR error("get_subgrid_vertices");

    compute_num_grid_vertices
      (dimension, subgrid_axis_size, num_subgrid_vertices);
    if (num_subgrid_vertices < 1) { return; };
    // Note: subgrid_axis_size[d] >= 1 for all d

    if (!check_vertex_list(vlist, error)) { throw error; };

    // add vertices along dimension 0 to vlist.
    vlist[0] = subgrid_origin;
    for (ATYPE i = 1; i < subgrid_axis_size[0]; i++)
      { vlist[i] = subgrid_origin + i; }

    // initialize prev_num_vertices
    VTYPE prev_num_vertices = subgrid_axis_size[0];

    // Process axes 1,2,...,dimension-1
    VTYPE * vcur_ptr = vlist + prev_num_vertices;
    for (DTYPE d = 1; d < dimension; d++) {

      VTYPE iv0 = axis_increment[d];

      for (ATYPE i = 1; i < subgrid_axis_size[d]; i++) {
        for (VTYPE * vprev_ptr = vlist;
             vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
          *(vcur_ptr) = iv0 + *(vprev_ptr);
          vcur_ptr++;
        }
        iv0 = iv0+axis_increment[d];
      }

      prev_num_vertices = prev_num_vertices*subgrid_axis_size[d];
    }

    if (!check_num_vertices_added
        (prev_num_vertices, num_subgrid_vertices, error))
      throw error;
  }


  /*!
   *  @brief Get vertices in subgrid.
   *  - Returns vertex indices sorted in increasing order.
   *  @param dimension  Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param  subgrid_origin  Subgrid origin.
   *  @param subgrid_axis_size
   *         Array: <em>subgrid_axis_size[d]</em> =
   *         Number of vertices along subgrid axis d.
   *  @param[out] vlist[]  List of vertices.
   *  @pre \li Subgrid is contained in grid,
   *      i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
   *  @pre \li Array vlist[] is preallocated to length at least
   *           number of vertices in grid or subgrid.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_subgrid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   VTYPE * vlist)
  {
    VTYPE num_subgrid_vertices;

    IJK::ARRAY<VTYPE> axis_increment(dimension);
    compute_axis_increment(dimension, axis_size, axis_increment.Ptr());

    get_subgrid_vertices(dimension, axis_size, axis_increment.Ptr(),
                         subgrid_origin, subgrid_axis_size, 
                         vlist, num_subgrid_vertices);
  }


  /*!
   *  @brief Get cubes in subgrid.
   *  @param dimension Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param  subgrid_origin  Subgrid origin.
   *  @param subgrid_axis_size
   *         Array: <em>subgrid_axis_size[d]</em> = Number of vertices along subgrid axis d.
   *  @param[out] vlist[]  List of vertices.
   *  @pre \li Subgrid is contained in grid,
   *     i.e. ( \a d'th coord of \a subgrid_origin ) + \a subgrid_axis_size[d] < \a axis_size[d].
   *  @pre \li Array cube_list[] is preallocated to length
   *     at least number of cubes in grid or subgrid.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_subgrid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size, 
   VTYPE * cube_list)
  {
    IJK::ARRAY<ATYPE> subgrid2_axis_size(dimension);

    for (DTYPE d = 0; d < dimension; d++) {
      if (subgrid_axis_size[d] < 2) { return; }
      subgrid2_axis_size[d] = subgrid_axis_size[d]-1;
    }

    get_subgrid_vertices
      (dimension, axis_size, subgrid_origin, subgrid2_axis_size.PtrConst(), 
       cube_list);
  }


  /*!
   *  @brief Get vertices in grid.
   *  @param dimension Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param[out] vlist[] List of vertices.
   *  @pre \li Array vlist[] is preallocated to length
   *    at least number of vertices in grid or subgrid.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   VTYPE * vlist)
  {
    subsample_grid_vertices(dimension, axis_size, 1, vlist);
  }


  /*!
   *  @brief Get grid vertices in region between two grid vertices (inclusive).
   *  @param dimension = Grid dimension.
   *  @param axis_size  Array: <em>axis_size[d]</em> = Number of vertices along axis \a d.
   *  @param iv0 = Lower grid vertex.
   *  @param iv1 = Upper grid vertex.
   *  @param[out] vlist[] = List of vertices between \a iv0 and \a iv1.
   *  @pre \li 0 <= iv0 <= iv1 < total_num_grid_vertices.
   *  @pre \li Array vlist[] is preallocated to size at least number of region vertices.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_grid_vertices_between
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const VTYPE iv1, VTYPE * vlist)
  {
    ATYPE coord0[dimension];
    ATYPE coord1[dimension];
    ATYPE region_size[dimension];
    IJK::PROCEDURE_ERROR error("get_grid_vertices_between");

    if (dimension < 0) { return; };
    if (!check_range(dimension, axis_size, iv0, iv1, error)) { throw error; };

    compute_coord(iv0, dimension, axis_size, coord0);
    compute_coord(iv1, dimension, axis_size, coord1);

    for (DTYPE d = 0; d < dimension; ++d) {
      if (!check_region_coordinates
          (dimension, iv0, coord0, iv1, coord1, error)) { throw error; };

      region_size[d] = coord1[d]-coord0[d]+1;
    }

    get_subgrid_vertices(dimension, axis_size, iv0, region_size, vlist);
  }


  /*!
   *  @brief Get grid vertices in region.
   *  @param dimension = Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param iv0 = Primary vertex of region.
   *  @param max_region_edge_length  Maximum number of grid edges contained in each region edge.
   *  @param[out] vlist = List of vertices in region.
   *  @pre \li 0 <= iv0 < total_num_grid_vertices.
   *  @pre \li Array vlist[] is preallocated to size at least number of region vertices.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_grid_vertices_in_region
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const ATYPE max_region_edge_length, VTYPE * vlist)
  {
    ATYPE coord[dimension];
    ATYPE region_size[dimension];

    compute_coord(iv0, dimension, axis_size, coord);

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] + max_region_edge_length < axis_size[d])
        { region_size[d] = max_region_edge_length + 1; }
      else if (coord[d] < axis_size[d])
        { region_size[d] = axis_size[d] - coord[d]; }
      else {
        // Vertex iv0 does not lie inside grid.
        return;
      }
    }

    get_subgrid_vertices(dimension, axis_size, iv0, region_size, vlist);
  }


  /*!
   *  @brief Get grid cubes in region.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param iv0  Primary vertex of region.
   *  @param max_region_edge_length  Maximum number of grid edges contained in each region edge.
   *  @param[out] vlist[] List of primary vertices of cubes in region.
   *  @param num_cubes Number of cubes in region (= number of vertices in vlist[].)
   *  @pre \li 0 <= iv0 < total_num_grid_vertices.
   *  @pre \li Array vlist[] is preallocated to size at least number of region vertices.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename NTYPE>
  void get_grid_cubes_in_region
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE iv0, const ATYPE max_region_edge_length, VTYPE * vlist,
   NTYPE & num_cubes)
  {
    IJK::ARRAY<ATYPE> coord(dimension);
    IJK::ARRAY<ATYPE> subgrid_size(dimension);

    num_cubes = 0;
    if (max_region_edge_length < 1) { return; };

    compute_coord(iv0, dimension, axis_size, coord.Ptr());

    for (DTYPE d = 0; d < dimension; d++) {
      if (coord[d] + max_region_edge_length < axis_size[d])
        { subgrid_size[d] = max_region_edge_length; }
      else if (coord[d]+1 < axis_size[d])
        { subgrid_size[d] = axis_size[d] - coord[d] - 1; }
      else {
        // Region contains no cubes
        return;
      }
    }

    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) 
      { num_cubes *= subgrid_size[d]; }

    get_subgrid_vertices
      (dimension, axis_size, iv0, subgrid_size.PtrConst(), vlist);
  }


  /*!
   *  @brief Get grid vertices in neighborhood around vertex \a iv.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param axis_increment[] = Axis increment. iv+axis_increment[i]
   *         is the next vertex after vertex iv along axis i.
   *  @param iv  Neighborhood around vertex \a iv.
   *  @param distance Distance to \a iv.
   *  @param[out] vlist[] List of vertices in neighborhood around \a iv.
   *  @pre dimension > 0.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE0, typename VTYPE1,
            typename DIST_TYPE>
  void get_grid_vertices_in_neighborhood
  (const DTYPE dimension, const ATYPE * axis_size, const ATYPE * axis_increment,
   const VTYPE0 iv, const DIST_TYPE distance,
   std::vector<VTYPE1> & vlist)
  {
    VTYPE0 region_iv0;
    IJK::ARRAY<ATYPE> region_iv0_coord(dimension);
    IJK::ARRAY<ATYPE> region_axis_size(dimension);
    IJK::ARRAY<ATYPE> vertex_coord(dimension);
    VTYPE0 num_subgrid_vertices;

    compute_coord(iv, dimension, axis_size, vertex_coord.Ptr());

    compute_region_around_vertex_coord
      (vertex_coord.PtrConst(), dimension, axis_size, 
       distance, region_iv0_coord.Ptr(), region_axis_size.Ptr());
    region_iv0 =
      compute_vertex_index<VTYPE0>
      (region_iv0_coord.PtrConst(), dimension, axis_size);

    compute_num_grid_vertices
      (dimension, region_axis_size.PtrConst(), num_subgrid_vertices);

    if (num_subgrid_vertices == 0) {
      vlist.clear();
      return;
    }

    // get all vertices including iv, then remove iv.
    vlist.resize(num_subgrid_vertices);
    get_subgrid_vertices(dimension, axis_size, axis_increment,
                         region_iv0, region_axis_size.PtrConst(),
                         &(vlist[0]), num_subgrid_vertices);

    // remove iv from list
    std::remove(vlist.begin(), vlist.end(), iv);
    vlist.pop_back();
  }

  
  /// Get primary cube vertices in subgrid.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_primary_cube_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const VTYPE subgrid_origin, const ATYPE * subgrid_axis_size,
   VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size2(dimension);

    for (DTYPE d = 0; d < dimension; d++) {
      if (subgrid_axis_size[d] < 2) 
        { return; }                      // zero cubes
      subgrid_axis_size2[d] = subgrid_axis_size[d]-1;
    }
    get_subgrid_vertices
      (dimension, axis_size, subgrid_origin, subgrid_axis_size2.PtrConst(), 
       vlist);
  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: FACET VERTICES, CUBES AND REGIONS
  // *****************************************************************

  ///@{

  /*!
   *  @brief Return number of vertices in specified grid facet
   *  - Specify grid facet by the direction orthogonal to the facet.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param orth_dir  Direction orthogonal to the facet.
   *  @param boundary_width  Width of boundary, (Number of vertices. Must be non-negative.)
   *  @param[out] num_vertices Number of vertices.
   */
  template <typename DTYPE, typename DTYPE2, typename ATYPE, 
            typename WTYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE2 orth_dir, const WTYPE boundary_width,
   NTYPE & num_vertices)
  {
    if (dimension < 1) { num_vertices = 0; }
    else if (dimension == 1 && boundary_width > 0)
      { num_vertices = 0; }
    else { num_vertices = 1; }

    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {
        if (axis_size[d] > ATYPE(2*boundary_width))
          { num_vertices *= (axis_size[d]-2*boundary_width); }
        else
          { num_vertices = 0; };
      }
      else if (axis_size[d] < 1)
        { num_vertices = 0; }
    }
  }


  /// @brief Return number of vertices in specified grid facet
  /// - Specify grid facet by the direction orthogonal to the facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir,
   NTYPE & num_vertices)
  {
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, orth_dir, 0, num_vertices);
  }


  /// @brief Return number of vertices in grid facet interior.
  /// - Same as compute_num_vertices_in_grid_facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_num_vertices_in_grid_facet_interior
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE2 orth_dir, const WTYPE boundary_width,
   NTYPE & num_vertices)
  {
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, orth_dir, boundary_width, num_vertices);
  }


  /// @brief Return number of vertices in specified grid ridge
  /// - Specify grid facet by the directions orthogonal to ridge.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_vertices_in_grid_ridge
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir0, const DTYPE orth_dir1,
   NTYPE & num_vertices)
  {
    FACET_LIST2<DTYPE> facet_list2(orth_dir0, orth_dir1);
    compute_num_vertices(dimension, axis_size, facet_list2, num_vertices);
  }


  /// Return maximum number of vertices over all grid facets.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_max_num_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   NTYPE & max_num_vertices)
  {
    max_num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      NTYPE num_face_vertices; 
      compute_num_vertices_in_grid_facet
        (dimension, axis_size, d, num_face_vertices);
      if (num_face_vertices > max_num_vertices)
        max_num_vertices = num_face_vertices;
    };
  }


  /// Return maximum number of vertices over all grid facet interiors.
  template <typename DTYPE, typename ATYPE, typename WTYPE, typename NTYPE>
  void compute_max_num_vertices_in_grid_facet_interior
  (const DTYPE dimension, const ATYPE * axis_size, const WTYPE boundary_width,
   NTYPE & max_num_vertices)
  {
    max_num_vertices = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      NTYPE num_face_vertices; 
      compute_num_vertices_in_grid_facet_interior
        (dimension, axis_size, d, boundary_width, num_face_vertices);
      if (num_face_vertices > max_num_vertices)
        max_num_vertices = num_face_vertices;
    };
  }


  /// @brief Return number of facets in specified grid facet.
  /// - Facet is lower facet orthogonal to the specified direction
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_facets_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE orth_dir,
   NTYPE & num_facets)
  {
    num_facets = 1;

    if (axis_size[orth_dir] <= 0) {
      num_facets = 0;
      return;
    }

    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {
        if (axis_size[d] <= 1) { 
          num_facets = 0;
          return; 
        };
        num_facets = num_facets*(axis_size[d]-1);
      };
    }
  }


  /// Return number of facets in grid facet orthogonal to axis 0.
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_facets_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_facets)
  {
    if (dimension < 1) 
      { num_facets = 0; }
    else {
      compute_num_facets_in_grid_facet
        (dimension, axis_size, DTYPE(0), num_facets);
    }
  }


  /// @brief Return number of cubes in specified grid facet.
  /// - Facet is lower facet orthogonal to the specified direction
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE orth_dir,
   NTYPE & num_cubes)
  {
    num_cubes = 1;
    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= 1) { 
        num_cubes = 0;
        return; 
      };
      if (d != orth_dir) {
        num_cubes = num_cubes*(axis_size[d]-1);
      };
    }
  }


  /// Return number of cubes in grid facet orthogonal to axis 0
  template <typename DTYPE, typename ATYPE, typename NTYPE>
  void compute_num_cubes_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, NTYPE & num_cubes)
  {
    if (dimension < 1) 
      { num_cubes = 0; }
    else {
      compute_num_cubes_in_grid_facet
        (dimension, axis_size, DTYPE(0), num_cubes);
    }
  }


  /*!
   *  @brief Get vertices in specified grid facet.
   *  - Does not return any vertices if \a axis_size[orth_dir] == 0.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param orth_dir = Direction orthogonal to facet.
   *  @param side = Side of grid containing facet.  If false, facet contains (0,0,...,0) origin.
   *  @param boundary_width = Width of boundary, (Number of vertices. Must be non-negative.)
   *  @param[out] vlist[] = List of primary vertices.
   *  @pre Array vlist[] must be pre-allocated to size at least number of vertices in facet.
   */
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename BTYPE, typename VTYPE>
  void get_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir, 
   const bool side, const BTYPE boundary_width, VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);
    IJK::ARRAY<ATYPE> coord(dimension);

    if (dimension < 1) { return; };
    if (axis_size[orth_dir] < 1) { return; };
    if (dimension == 1 && boundary_width > 0) { return; };
    VTYPE subgrid_origin = 0;

    if (boundary_width == 0) {
      std::copy(axis_size, axis_size+dimension, 
                subgrid_axis_size.Ptr());
    
      if (side) {
        for (DTYPE d = 0; d < dimension; d++) { coord[d] = 0; };
        coord[orth_dir] = axis_size[orth_dir]-1;
        subgrid_origin = 
          compute_vertex_index<VTYPE>(coord.PtrConst(), dimension, axis_size);
      }
    }
    else {
      for (DTYPE d = 0; d < dimension; d++) {
        if (d != orth_dir) {
          if (axis_size[d] < ATYPE(2*boundary_width)) { return; };

          coord[d] = boundary_width; 
          subgrid_axis_size[d] = axis_size[d]-2*boundary_width;
        }
      }

      if (side) { coord[orth_dir] = axis_size[orth_dir]-1; }
      else { coord[orth_dir] = 0; }

      subgrid_origin = 
        compute_vertex_index<VTYPE>(coord.PtrConst(), dimension, axis_size);
    }
    subgrid_axis_size[orth_dir] = 1;

    get_subgrid_vertices(dimension, axis_size, subgrid_origin, 
                         subgrid_axis_size.PtrConst(), vlist);
  }

  /// Get vertices in specified grid facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename VTYPE>
  void get_vertices_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir, 
   const bool side, VTYPE * vlist)
  {
    get_vertices_in_grid_facet
      (dimension, axis_size, orth_dir, side, 0, vlist);
  }

  /// Get vertices in grid facet 0.
  template <typename DTYPE, typename ATYPE, typename BTYPE, typename VTYPE>
  void get_vertices_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, 
   const BTYPE boundary_width, VTYPE * vlist)
  {
    get_vertices_in_grid_facet
      (dimension, axis_size, DTYPE(0), false, boundary_width, vlist);
  }

  /// Get vertices in grid facet 0.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_vertices_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_vertices_in_grid_facet0(dimension, axis_size, 0, vlist);
  }


  /// @brief Get vertices in grid facet interior.
  /// - Same as get_vertices_in_grid_facet.
  template <typename DTYPE, typename DTYPE2, typename ATYPE, typename BTYPE, typename VTYPE>
  void get_vertices_in_grid_facet_interior
  (const DTYPE dimension, const ATYPE * axis_size, const DTYPE2 orth_dir, 
   const bool side, const BTYPE boundary_width, VTYPE * vlist)
  {
    get_vertices_in_grid_facet(dimension, axis_size, orth_dir, side,
                               boundary_width, vlist);
  }


  /*!
   *  @brief Get vertices in specified grid ridge
   *  - Ridge is lower ridge orthogonal to the specified directions
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param orth_dir0 is first direction orthogonal to the ridge.
   *  @param orth_dir1 is second direction orthogonal to the ridge.
   *  @param[out] vlist[] = List of vertices in grid ridge.
   *  @pre Array vlist[] is preallocated to length
   *    at least number of vertices in grid ridge.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_vertices_in_grid_ridge
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir0, const DTYPE orth_dir1, VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);

    if (axis_size[orth_dir0] < 1 || axis_size[orth_dir1] < 1) { return; }
    std::copy(axis_size, axis_size+dimension, subgrid_axis_size.Ptr());
    subgrid_axis_size[orth_dir0] = 1;
    subgrid_axis_size[orth_dir1] = 1;
    
    get_subgrid_vertices
      (dimension, axis_size, 0, subgrid_axis_size.PtrConst(), vlist);
  }


  /*!
   *  @brief Get primary vertices of facets in specified grid facet.
   *  - Does not return any vertices if \a axis_size[orth_dir] == 0.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param orth_dir = Direction orthogonal to facet.
   *  @param side = Side of grid containing facet.  If false, facet contains (0,0,...,0) origin.
   *  @param[out] vlist[] = List of primary vertices.
   *  @pre Array vlist[] must be pre-allocated to size
   *    at least number of primary vertices in facet.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_facets_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, const bool side, VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);

    if (axis_size[orth_dir] < 1) { return; };
    std::copy(axis_size, axis_size+dimension, subgrid_axis_size.Ptr());
    subgrid_axis_size[orth_dir] = 2;

    VTYPE subgrid_origin = 0;

    if (side) {
      IJK::ARRAY<VTYPE> axis_increment(dimension);
      compute_increment(dimension, axis_size, axis_increment.Ptr());

      subgrid_origin = axis_increment[orth_dir]*(axis_size[orth_dir]-1);
    }

    get_primary_cube_vertices
      (dimension, axis_size, subgrid_origin, subgrid_axis_size.PtrConst(), 
       vlist);
  }

  /// @brief Get primary vertices of facets in specified grid facet.
  /// - Version without side parameter.  (side=false is lower/leftmost side.)
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_facets_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, VTYPE * vlist)
  {
    get_facets_in_grid_facet(dimension, axis_size, orth_dir, false, vlist);
  }


  /*!
   *  @brief Get primary vertices of (d-1)-dimensional cubes
   *    in specified grid facet where d = \a dimension.
   *  - Does not return any vertices if \a axis_size[orth_dir] <= 0.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param orth_dir = Direction orthogonal to facet.
   *  @param side = Side of grid containing facet.  If false, facet contains (0,0,...,0) origin.
   *  @param[out] vlist[] = List of primary vertices.
   *  @pre Array vlist[] must be pre-allocated to size at least number of primary vertices in facet.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, const bool side, VTYPE * vlist)
  {
    IJK::ARRAY<ATYPE> subgrid_axis_size(dimension);

    if (axis_size[orth_dir] < 2) { return; };
    std::copy(axis_size, axis_size+dimension, subgrid_axis_size.Ptr());
    subgrid_axis_size[orth_dir] = 2;

    VTYPE subgrid_origin = 0;

    if (side) {
      IJK::ARRAY<VTYPE> axis_increment(dimension);
      compute_increment(dimension, axis_size, axis_increment.Ptr());

      subgrid_origin = axis_increment[orth_dir]*(axis_size[orth_dir]-2);
    }

    get_primary_cube_vertices
      (dimension, axis_size, subgrid_origin, subgrid_axis_size.PtrConst(), 
       vlist);
  }


  /*!
   *  @brief Get primary vertices of (d-1)-dimensional cubes
   *     in specified grid facet where d = \a dimension.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param orth_dir  Direction orthogonal to facet.  Facet contains (0,0,...,0).
   *  @param[out] vlist[] = List of primary vertices.
   *  @pre Array vlist[] must be pre-allocated to size at least number of primary vertices in facet.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_cubes_in_grid_facet
  (const DTYPE dimension, const ATYPE * axis_size,
   const DTYPE orth_dir, VTYPE * vlist)
  {
    get_cubes_in_grid_facet(dimension, axis_size, orth_dir, false, vlist);
  }


  /// @brief Get primary vertices of (d-1)-dimensional cubes in grid facet 0 where d = \a dimension.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_cubes_in_grid_facet0
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_cubes_in_grid_facet(dimension, axis_size, DTYPE(0), vlist);
  }


  /// Get outer grid vertices
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_outer_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
    // Precondition: vlist[] is preallocated to size 
    //   at least num_outer_vertices
  {
    VTYPE axis_increment[dimension];

    if (dimension < 1) { return; };

    DTYPE d_last = dimension - 1;

    if (dimension == 1) {
      if (axis_size[0] < 2) { return; }
      else {
        vlist[0] = axis_size[0]-1;
        return;
      }
    }
    else {
      if (axis_size[d_last] < 2) { return; }

      get_outer_grid_vertices(d_last, axis_size, vlist);
      VTYPE num_vertices; 
      compute_num_outer_vertices(d_last, axis_size, num_vertices);

      compute_increment(dimension, axis_size, axis_increment);

      for (VTYPE i = 1; i < axis_size[d_last]-1; i++) {
        VTYPE k = i*num_vertices;
        VTYPE k_increment = i*axis_increment[d_last];
        for (VTYPE j = 0; j < num_vertices; j++) {
          vlist[k+j] = vlist[j]+k_increment;
        } 
      }

      num_vertices = num_vertices * (axis_size[d_last]-1);

      VTYPE * vlist2 = vlist + num_vertices;
      FACET_LIST1<DTYPE> facet_list1(d_last);

      get_vertices_in_grid_facet(dimension, axis_size, d_last, true, vlist2);
    }
  }


  /*!
   *  @brief Get primary vertices of regions in grid or subgrid.
   *  @param dimension Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param region_edge_length = Number of grid edges contained in each region edge.
   *  @param[out] vlist[] Array of primary vertices of full regions.
   *  @pre Array vlist[] must be pre-allocated to size at least number of full regions.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);
    
    IJK::PROCEDURE_ERROR error("get_region_primary_vertices");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] < 2) { return; };

      subgrid_axis_size[d] = axis_size[d]-1;
    }

    subsample_subgrid_vertices
      (dimension, axis_size, VTYPE(0), subgrid_axis_size, subsample_period,
       vlist);
  }


  /*!
   *  @brief Get primary vertices of regions in grid or subgrid.
   *  @param dimension  Grid dimension.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param region_edge_length Number of grid edges contained in each region edge.
   *  @param[out] vlist[] Array of primary vertices of full regions.
   *  @param[out] is_full Array: <em>is_full[k]</em> = True
   *     if region \a k is a full \a LxLx...xL region where \a L = \a region_edge_length.
   *  @pre Array vlist[] must be pre-allocated to size at least number of regions.
   *  @pre Array is_full[] must be pre-allocated to size at least number of regions.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist, bool * is_full)
  {
    IJK::PROCEDURE_ERROR error("get_region_primary_vertices");
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    VTYPE num_regions;
    compute_num_regions
      (dimension, axis_size, region_edge_length, num_regions);

    if (num_regions < 1) { return; }
    // Note: axis_size[d] >= 2 for all d

    VTYPE num_full_regions;
    compute_num_full_regions
      (dimension, axis_size, region_edge_length, num_full_regions);

    if (!check_vertex_list(vlist, error)) { throw error; };

    IJK::ARRAY<VTYPE> subsample_increment(dimension);
    compute_subsample_increment
      (dimension, axis_size, subsample_period, subsample_increment.Ptr());

    // set vlist[0], is_full[0] and initialize prev_num_vertices
    vlist[0] = 0;
    if (num_full_regions > 0) { is_full[0] = true; }
    else { is_full[0] = false; };
    VTYPE prev_num_regions = 1;

    // Process axes 0,1,2,..., dimension-1
    for (DTYPE d = 0; d < dimension; d++) {

      ATYPE num_regions_along_axis =
        compute_num_regions_along_axis(axis_size[d], region_edge_length);

      ATYPE num_full_along_axis =
        compute_num_full_regions_along_axis(axis_size[d], region_edge_length);

      VTYPE iv0 = subsample_increment[d];
      for (VTYPE i = 1; i < num_full_along_axis; i++) {
        VTYPE i2 = i * prev_num_regions;
        for (VTYPE j = 0; j < prev_num_regions; j++) {
          VTYPE k = j + i2;
          vlist[k] = vlist[j] + iv0;
          is_full[k] = is_full[j];
        }
        iv0 = iv0 + subsample_increment[d];
      }

      if (num_regions_along_axis != num_full_along_axis &&
          num_regions_along_axis > 1) {
        VTYPE i2 = num_full_along_axis * prev_num_regions;
        for (VTYPE j = 0; j < prev_num_regions; j++) {
          VTYPE k = j + i2;
          vlist[k] = vlist[j] + iv0;
          is_full[k] = false;
        }
      }

      prev_num_regions = prev_num_regions*num_regions_along_axis;
    }

    if (!check_num_vertices_added(prev_num_regions, num_regions, error))
      throw error;

  }


  /*!
   *  @brief Get primary vertices of full regions in grid.
   *  @param dimension Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param region_edge_length = Number of grid edges contained in each region edge.
   *  @param[out] vlist[] Array of primary vertices of full regions.
   *  @pre Array vlist[] must be pre-allocated to size at least number of full regions.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_full_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    ATYPE subgrid_axis_size[dimension];
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    IJK::PROCEDURE_ERROR error("get_full_region_primary_vertices");
    if (!check_region_edge_length(region_edge_length, error)) { throw error; };

    for (DTYPE d = 0; d < dimension; d++) {
      if (axis_size[d] <= region_edge_length) { return; };

      subgrid_axis_size[d] = axis_size[d]-region_edge_length;
    }

    subsample_subgrid_vertices
      (dimension, axis_size, VTYPE(0), subgrid_axis_size, subsample_period,
       vlist);
  }


  /*!
   *  @brief Get primary vertices of partial regions in grid.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param region_edge_length = Number of grid edges contained in each region edge.
   *  @param[out] vlist[] = Array of primary vertices of partial regions.
   *  @pre Array vlist[] must be pre-allocated to size at least number of partial regions.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_partial_region_primary_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const ATYPE region_edge_length, VTYPE * vlist)
  {
    IJK::PROCEDURE_ERROR error("get_partial_region_primary_vertices");
    IJK::CONSTANT<ATYPE,ATYPE> subsample_period(region_edge_length);

    VTYPE num_vertices;
    compute_num_partial_regions
      (dimension, axis_size, region_edge_length, num_vertices);

    if (num_vertices < 1) { return; };
    // Note: axis_size[d] >= 1 for all d not in facet_list

    if (!check_vertex_list(vlist, error)) { throw error; };

    VTYPE region_axis_increment[dimension];

    compute_subsample_increment(dimension, axis_size, subsample_period,
                                region_axis_increment);

    // initialize prev_num_vertices
    VTYPE prev_num_vertices = 0;

    // Process axes 0,1,2,..., dimension-1
    for (DTYPE d = 0; d < dimension; d++) {
      VTYPE * vcur_ptr = vlist + prev_num_vertices;

      ATYPE num_full_regions_along_axis = 
        compute_num_full_regions_along_axis
        (axis_size[d], region_edge_length);

      VTYPE iv0 = region_axis_increment[d];
      for (VTYPE i = 1; i < num_full_regions_along_axis; i++) {
        for (VTYPE * vprev_ptr = vlist; 
             vprev_ptr != vlist+prev_num_vertices; vprev_ptr++) {
          *(vcur_ptr) = iv0 + *(vprev_ptr);
          vcur_ptr++;
        };
        iv0 = iv0 + region_axis_increment[d];
      }

      prev_num_vertices = prev_num_vertices*num_full_regions_along_axis;

      ATYPE num_partial_regions_along_axis =
        compute_num_partial_regions_along_axis
        (axis_size[d], region_edge_length);

      if (num_partial_regions_along_axis > 0) {

        VTYPE inc = 
          region_axis_increment[d]*num_full_regions_along_axis;

        if (d == 0) {
          vlist[prev_num_vertices] = inc;
          prev_num_vertices++;
        }
        else {

          get_region_primary_vertices(d, axis_size, region_edge_length, 
                                      vlist + prev_num_vertices);

          ATYPE k;
          compute_num_regions(d, axis_size, region_edge_length, k);
          for (VTYPE i = prev_num_vertices; i < prev_num_vertices+k; i++)
            { vlist[i] += inc; }

          prev_num_vertices += k;
        }
      }
    }

    if (prev_num_vertices != num_vertices) {
      error.AddMessage("Programming error.  Added ", prev_num_vertices, 
                       " vertices to vertex list.");
      error.AddMessage("Number of vertices in list should be ", 
                       num_vertices, ".");
      throw error;
    }

  }


  /*!
   *  @brief Get vertices on boundary of grid.
   *  - Allows boundary_width to be greater than 1.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param boundary_width Width of boundary, in vertices.
   *  @param[out] vlist[] List of boundary vertices.
   *  @pre Array vlist[] must be pre-allocated to size at least number
   *    of boundary vertices.
   */
  template <typename DTYPE, typename ATYPE, typename VTYPE, typename WTYPE>
  void get_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, 
   const WTYPE boundary_width, VTYPE * vlist)
  {
    if (dimension < 1) { return; }
    if (boundary_width < 1) { return; }

    const DTYPE d_last = dimension - 1;

    if (axis_size[d_last] <= 2*boundary_width) {
      // all vertices are on the boundary
      VTYPE num_grid_vertices;
      compute_num_grid_vertices(dimension, axis_size, num_grid_vertices);
      for (VTYPE j = 0; j < num_grid_vertices; j++) 
        { vlist[j] = j; }
      return;
    }

    if (dimension == 1) {
      for (VTYPE j = 0; j < boundary_width; j++) 
        { vlist[j] = j; }

      for (VTYPE j = 0; j < boundary_width; j++) {
        VTYPE iv = axis_size[0]-boundary_width + j;
        vlist[j+boundary_width] = iv;
      };
      return;
    }

    IJK::ARRAY<ATYPE> axis_increment(dimension);
    compute_increment(dimension, axis_size, axis_increment.Ptr());

    // get vertices in lower facet
    VTYPE num_vertices_in_grid_facet;
    compute_num_vertices_in_grid_facet
      (dimension, axis_size, d_last, num_vertices_in_grid_facet);
    get_vertices_in_grid_facet(dimension, axis_size, d_last, false, vlist);

    // get remaining vertices in lower boundary
    for (VTYPE i = 1; i < boundary_width; i++) {
      VTYPE inc = i*axis_increment[d_last];
      for (VTYPE j = 0; j < num_vertices_in_grid_facet; j++) {
        vlist[i*num_vertices_in_grid_facet + j] = vlist[j] + inc;
      }
    }

    VTYPE * vlist2 = vlist+boundary_width*num_vertices_in_grid_facet;
    get_boundary_grid_vertices(dimension-1, axis_size, 
                               boundary_width, vlist2);

    VTYPE num_boundary_grid_vertices;
    compute_num_boundary_grid_vertices
      (dimension-1, axis_size, boundary_width, num_boundary_grid_vertices);
    for (VTYPE * vcur_ptr = vlist2; 
         vcur_ptr != vlist2+num_boundary_grid_vertices; vcur_ptr++)
      { *vcur_ptr += boundary_width*axis_increment[d_last]; }

    VTYPE * vlist3 = vlist2+num_boundary_grid_vertices;
    for (ATYPE j = boundary_width+1; j+boundary_width < axis_size[d_last]; 
         j++) {
      VTYPE inc = axis_increment[d_last]*(j-boundary_width);

      for (VTYPE i = 0; i < num_boundary_grid_vertices; i++)  
        { vlist3[i] = vlist2[i] + inc; }

      vlist3 += num_boundary_grid_vertices;
    }

    VTYPE * vlist4 = vlist3 + (boundary_width-1)*num_vertices_in_grid_facet;

    // get vertices in upper facet
    get_vertices_in_grid_facet(dimension, axis_size, d_last, true, vlist4);

    // get remaining vertices in upper boundary
    for (VTYPE i = 0; i+1 < boundary_width; i++) {
      VTYPE inc = (boundary_width-i-1)*axis_increment[d_last];
      for (VTYPE j = 0; j < num_vertices_in_grid_facet; j++) {
        vlist3[i*num_vertices_in_grid_facet + j] = vlist4[j] - inc;
      }
    }
  
  }


  /// Get vertices on boundary of grid.
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_boundary_grid_vertices
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * vlist)
  {
    get_boundary_grid_vertices(dimension, axis_size, 1, vlist);
  }


  /// Get boundary grid cubes
  template <typename DTYPE, typename ATYPE, typename VTYPE>
  void get_boundary_grid_cubes
  (const DTYPE dimension, const ATYPE * axis_size, VTYPE * cube_list)
  {
    if (dimension < 1) { return; }
    if (dimension == 1) {
      if (axis_size[0] > 0) { cube_list[0] = 0; }
      if (axis_size[0] > 2) { cube_list[1] = axis_size[0]-2; };
      return;
    }

    DTYPE d_last = dimension - 1;
    if (axis_size[d_last] < 2) { return; };

    // get vertices in lower facet
    VTYPE num_cubes_in_grid_facet;
    compute_num_cubes_in_grid_facet(dimension, axis_size, d_last,
                                    num_cubes_in_grid_facet);
    get_cubes_in_grid_facet(dimension, axis_size, d_last, false, cube_list);

    VTYPE * cube_list2 = cube_list+num_cubes_in_grid_facet;
    VTYPE * cube_list3 = cube_list2;
    if (axis_size[d_last] > 3) {
      ATYPE axis_increment[dimension];

      compute_increment(dimension, axis_size, axis_increment);
      get_boundary_grid_cubes(dimension-1, axis_size, cube_list2);

      VTYPE n;
      compute_num_boundary_grid_cubes(dimension-1, axis_size, n);
      for (VTYPE * vcur_ptr = cube_list2; vcur_ptr != cube_list2+n; vcur_ptr++)
        { *vcur_ptr += axis_increment[d_last]; }

      cube_list3 = cube_list2 + n;
      for (ATYPE j = 2; j+2 < axis_size[d_last]; j++) {
        VTYPE inc = axis_increment[d_last]*(j-1);

        for (VTYPE i = 0; i < n; i++)  { cube_list3[i] = cube_list2[i] + inc; }

        cube_list3 += n;
      }
    }

    get_cubes_in_grid_facet(dimension, axis_size, d_last, true, cube_list3);
  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: GRID BOUNDARIES
  // *****************************************************************

  ///@{

  /// Return true if cube_facet is on grid boundary.
  template <typename DTYPE, typename ATYPE, 
            typename DTYPE2, typename CTYPE>
  bool is_cube_facet_on_grid_boundary
  (const DTYPE dimension, const ATYPE * axis_size,
   const CTYPE cube_index, const DTYPE2 facet_orth_dir,
   const bool facet_side)
  {
    long boundary_bits;
    compute_boundary_cube_bits
      (cube_index, dimension, axis_size, boundary_bits);

    // Convert kf to index into boundary bits.
    long bit_index;
    IJK::compute_boundary_bit_index(facet_orth_dir, facet_side, bit_index);
    long mask = (long(1) << bit_index);
    if ((boundary_bits & mask) == 0) 
      { return(false); }
    else 
      { return(true); }
  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: COMPUTING NEIGHBORS
  // *****************************************************************

  ///@{

  /*!
   *  \brief Compute number of neighbors of a vertex in all cubes
   *    containing the vertex.
   *  - Does not count the vertex itself.
   */
  template <typename DTYPE, typename NTYPE> 
  void compute_num_vertex_neighborsC
  (const DTYPE dimension, NTYPE & num_neighbors)
  { 
    num_neighbors = 1;
    for (DTYPE d = 0; d < dimension; d++) 
      { num_neighbors = num_neighbors*3; }
    num_neighbors = num_neighbors-1;
  }

  
  /// \brief Compute number of vertices which share an edge with a vertex.
  /// - Does not count the vertex itself.
  template <typename DTYPE, typename NTYPE> 
  void compute_num_vertex_neighborsE
  (const DTYPE dimension, NTYPE & num_neighbors)
  { 
    num_neighbors = 2*dimension;
  }

  
  /// \brief Compute number of vertices in cubes containing a facet,
  ///        not including facet vertices.
  template <typename DTYPE, typename NTYPE> 
  void compute_num_facet_neighborsC
  (const DTYPE dimension, NTYPE & num_neighbors)
  {
    NTYPE num_facet_vertices = compute_num_cube_facet_vertices(dimension);
    num_neighbors = 2*num_facet_vertices;
  }


  /*!
   *  \brief Compute number of vertices in 2-faces containing an edge,
   *         not including edge vertices.
   *  @param dimension  Dimension of grid.
   *  @param[out] num_neighbors Number of vertices in 2-faces 
   *     containing an edge, not including the edge vertices.
   */
  template <typename DTYPE, typename NTYPE> 
  void compute_num_edge_neighborsF2
  (const DTYPE dimension, NTYPE & num_neighbors)
  {
    num_neighbors = 0;
    if (dimension > 0) 
      { num_neighbors = 4*(dimension-1); }
  }


  /*!
   *  \brief Compute integer to add to vertex index to compute vertex neighbors.
   *  - Use only for vertex neighbors of internal vertices.
   *  @param dimension  Dimension of grid.
   *  @param axis_size[d] Number of vertices along grid axis d.
   *  @param[out] vertex_neighborC Array. 
   *     - iv + vertex_neighborC[k] = index of k'th vertex neighbor of iv
   *       where iv is the index of an internal grid vertex.
   */
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_vertex_neighborC
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * vertex_neighborC)
  {
    IJK::PROCEDURE_ERROR error("compute_vertex_neighborC");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    if (dimension == 0) { return; }

    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    compute_increment(dimension, axis_size, axis_increment.Ptr());

    // iv0 = index of vertex (1,1,...,1).
    DIFFTYPE iv0 = 0;
    for (DTYPE d = 0; d < dimension; d++) 
      { iv0 += axis_increment[d]; };

    vertex_neighborC[0] = -iv0;
    DIFFTYPE k = 1;

    // Process axes 0,1,2,..., dimension-1
    DIFFTYPE * vcur_ptr = vertex_neighborC + k;
    DIFFTYPE * vcenter_ptr = vertex_neighborC;
    for (DTYPE d = 0; d < dimension; d++) {

      for (ATYPE j = 1; j < 3; j++) {
        for (DIFFTYPE * vprev_ptr = vertex_neighborC; 
             vprev_ptr != vertex_neighborC+k; vprev_ptr++) {
          DIFFTYPE iv = j*axis_increment[d] + *(vprev_ptr);
          if ((d+1 < dimension) || (j != 1) || 
              (vprev_ptr != vcenter_ptr)) { 
            *(vcur_ptr) = iv; 
            vcur_ptr++;
          }
        }
      }
      vcenter_ptr += k;

      k = vcur_ptr - vertex_neighborC;
    }

    DIFFTYPE num_neighbors; 
    compute_num_vertex_neighborsC(dimension, num_neighbors);

    if (!check_num_vertices_added(k, num_neighbors, error))
      throw error;
  }

  /// @brief Compute integer to add to vertex index to compute vertex neighbors
  ///   across edges.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_vertex_neighborE
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * vertex_neighborE)
  {
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::PROCEDURE_ERROR error("compute_vertex_neighborE");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    DIFFTYPE num_neighbors; 
    compute_num_vertex_neighborsE(dimension, num_neighbors);

    compute_increment(dimension, axis_size, axis_increment.Ptr());

    for (DIFFTYPE d = 0; d < dimension; d++) 
      { vertex_neighborE[d] = -axis_increment[d]; };

    for (DIFFTYPE d = 0; d < dimension; d++) 
      { vertex_neighborE[d + dimension] = axis_increment[d]; };
  }

  /// @brief Compute integer to add to cube_index to compute cubes sharing
  ///   an edge with cube cube_index.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_cube_neighborE
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * cube_neighborE)
  {
    if (dimension < 1) {
      // No edges or facet_vertices.
      return;
    }

    const DIFFTYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    const DIFFTYPE num_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);
    const DIFFTYPE num_cube_edges = compute_num_cube_edges(dimension);
    const DIFFTYPE facet_vlast = num_facet_vertices-1;
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::ARRAY<DIFFTYPE> cube_vertex_increment(num_cube_vertices);
    IJK::ARRAY<DIFFTYPE> facet_vertex_increment(num_facet_vertices);
    IJK::PROCEDURE_ERROR error("compute_cube_neighborE");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    if (num_cube_vertices <= 0) { return; }
    if (num_cube_edges <= 0) { return; }

    compute_increment(dimension, axis_size, axis_increment.Ptr());
    compute_cube_vertex_increment
      (dimension, axis_increment.PtrConst(), cube_vertex_increment.Ptr());

    facet_vertex_increment[0] = -cube_vertex_increment[0];

    DIFFTYPE k = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      compute_facet_vertex_increment
        (dimension, d, cube_vertex_increment.PtrConst(), 
         facet_vertex_increment.Ptr());

      for (int i = 0; i < num_facet_vertices; i++) {
        cube_neighborE[k] = 
          2*facet_vertex_increment[i] - facet_vertex_increment[facet_vlast];
        k++;
      }
    }

    if (k != num_cube_edges) {
      error.AddMessage
        ("Programming error.  Wrong number of cubes added to cube_neighborE.");
      throw error;
    }

  }

  /// @brief Compute integer to add to cube_index to compute cubes sharing
  ///   a vertex with cube cube_index.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_cube_neighborV
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * cube_neighborV)
  {
    const DIFFTYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::ARRAY<DIFFTYPE> cube_vertex_increment(num_cube_vertices);
    IJK::PROCEDURE_ERROR error("compute_cube_neighborV");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    if (num_cube_vertices <= 0) { return; }

    compute_increment(dimension, axis_size, axis_increment.Ptr());
    compute_cube_vertex_increment
      (dimension, axis_increment.PtrConst(), cube_vertex_increment.Ptr());

    cube_neighborV[0] = 0;
    for (DTYPE d = 0; d < dimension; d++)
      { cube_neighborV[0] -= axis_increment[d]; }

    for (DIFFTYPE k = 0; k < num_cube_vertices; k++)
      { cube_neighborV[k] = cube_neighborV[0] + 2*cube_vertex_increment[k]; }
  }

  /// Compute integer to add to vertex index to compute facet neighbors.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_facet_neighborC
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * facet_neighborC)
  {
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::PROCEDURE_ERROR error("compute_facet_neighborC");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    compute_increment(dimension, axis_size, axis_increment.Ptr());

    DIFFTYPE num_neighbors;
    compute_num_facet_neighborsC(dimension, num_neighbors);

    const DIFFTYPE num_cube_vertices = 
      compute_num_cube_vertices(dimension);
    IJK::ARRAY<DIFFTYPE> cube_vertex_increment(num_cube_vertices);

    compute_cube_vertex_increment
      (dimension, axis_increment.PtrConst(), cube_vertex_increment.Ptr());

    const DIFFTYPE num_facet_vertices = 
      compute_num_cube_facet_vertices(dimension);
    IJK::ARRAY<DIFFTYPE> facet_vertex_increment(dimension*num_facet_vertices);
    
    // for each facet defined by a different orthogonal direction
    for (DTYPE orth_dir = 0; orth_dir < dimension; orth_dir++) {

      DIFFTYPE * facet_ptr = 
        facet_vertex_increment.Ptr()+num_facet_vertices*orth_dir;

      compute_facet_vertex_increment
        (dimension, orth_dir, cube_vertex_increment.PtrConst(), facet_ptr);
      
      // for each vertex in the prvious facet
      for (DIFFTYPE k = 0; k < num_facet_vertices; k++) {
        facet_neighborC[k+num_neighbors*orth_dir] =
          facet_ptr[k] - axis_increment[orth_dir];
      }

      // for each vertex in the next facet
      for (DIFFTYPE k = 0; k < num_facet_vertices; k++) {
        facet_neighborC[k+num_facet_vertices+num_neighbors*orth_dir] =
          facet_ptr[k] + axis_increment[orth_dir];
      }
    }
  }

  /// Compute integer to add to vertex index to compute edge neighbors.
  template <typename DTYPE, typename ATYPE, typename DIFFTYPE>
  void compute_edge_neighborF2
  (const DTYPE dimension, const ATYPE * axis_size,
   DIFFTYPE * edge_neighborF2)
  {
    IJK::ARRAY<DIFFTYPE> axis_increment(dimension);
    IJK::PROCEDURE_ERROR error("compute_edge_neighborF2");

    if (!check_difftype<DIFFTYPE>(1, error)) 
      { throw error; }

    compute_increment(dimension, axis_size, axis_increment.Ptr());

    DIFFTYPE num_neighbors;
    compute_num_edge_neighborsF2(dimension, num_neighbors);

    DIFFTYPE k = 0;
    // for each edge determined by a direction
    for (DTYPE dir = 0; dir < dimension; dir++) {

      // for each edge incident on the lower vertex
      for (DTYPE d = 0; d < dimension; d++) {
        if (d != dir) {
          edge_neighborF2[k] = -axis_increment[d];
          k++;
        }
      }

      for (DTYPE d = 0; d < dimension; d++) {
        if (d != dir) {
          edge_neighborF2[k] = axis_increment[d];
          k++;
        }
      }

      // for each edge incident on the upper vertex
      for (DTYPE d = 0; d < dimension; d++) {
        if (d != dir) {
          edge_neighborF2[k] = axis_increment[dir]-axis_increment[d];
          k++;
        }
      }

      for (DTYPE d = 0; d < dimension; d++) {
        if (d != dir) {
          edge_neighborF2[k] = axis_increment[dir]+axis_increment[d];
          k++;
        }
      }
    }

    if (k != num_neighbors*dimension) {
      error.AddMessage
        ("Programming error.  Wrong number of edge neighbors added to edge_neighborF2[].");
      throw error;
    }

  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: COMPUTING DISTANCE
  // *****************************************************************

  ///@{

  /// Compute L-infinity distance between two grid vertices.
  template <typename GTYPE, typename VTYPE0, typename VTYPE1,
            typename DIST_TYPE, typename DIR_TYPE>
  void compute_Linf_distance_between_grid_vertices
  (const GTYPE & grid, const VTYPE0 iv0, const VTYPE1 iv1,
   DIST_TYPE & Linf_distance, DIR_TYPE & axis)
  {
    typedef typename GTYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<DIST_TYPE> coord0(dimension), coord1(dimension);
    DIST_TYPE diff;

    grid.ComputeCoord(iv0, coord0.Ptr());
    grid.ComputeCoord(iv1, coord1.Ptr());

    Linf_distance = 0;
    axis = 0;
    for (DTYPE d = 0; d < dimension; d++) {

      if (coord0[d] < coord1[d])
        { diff = coord1[d] - coord0[d]; }
      else
        { diff = coord0[d] - coord1[d]; }

      if (diff > Linf_distance) {
        Linf_distance = diff;
        axis = d;
      }
    }
  }

  /// Compute L-infinity distance between two grid vertices.
  template <typename GTYPE, typename VTYPE0, typename VTYPE1,
            typename DIST_TYPE>
  void compute_Linf_distance_between_grid_vertices
  (const GTYPE & grid, const VTYPE0 iv0, const VTYPE1 iv1,
   DIST_TYPE & Linf_distance)
  {
    typedef typename GTYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<DIST_TYPE> coord0(dimension), coord1(dimension);
    DIST_TYPE diff;

    grid.ComputeCoord(iv0, coord0.Ptr());
    grid.ComputeCoord(iv1, coord1.Ptr());

    Linf_distance = 0;
    for (DTYPE d = 0; d < dimension; d++) {

      if (coord0[d] < coord1[d])
        { diff = coord1[d] - coord0[d]; }
      else
        { diff = coord0[d] - coord1[d]; }

      if (diff > Linf_distance) {
        Linf_distance = diff;
      }
    }
  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: QUERY
  // *****************************************************************

  ///@{

  /// @brief Return true if edge (endpointA, endpointB) contains vertex iv,
  ///   i.e., if (endpointA == iv) or (endpointB == iv).
  template <typename VTYPEA, typename VTYPEB, typename VTYPEC>
  bool does_edge_contain_vertex
  (const VTYPEA endpointA, const VTYPEB endpointB, const VTYPEC iv)
  {
    if (endpointA == iv || endpointB == iv) { return(true); }
    else { return(false); }
  }

  ///@}


  // *****************************************************************
  //! @name TEMPLATE FUNCTIONS: SET GRID COORDINATES
  // *****************************************************************

  ///@{
  
  /*!
   *  @brief Set grid coord in array.
   *  @tparam ELEMENT_TYPE Type of elements.
   *    ELEMENT_TYPE must have member function SetCoord().
   *  @param grid Grid.
   *  @param[out] list[] C++ vector.
   */
  template <typename GRID_TYPE, typename ELEMENT_TYPE>
  void set_grid_coord
  (const GRID_TYPE & grid, std::vector<ELEMENT_TYPE> & list)
  {
    typedef typename std::vector<ELEMENT_TYPE>::size_type SIZE_TYPE;

    for (SIZE_TYPE i = 0; i < list.size(); i++) 
      { list[i].SetCoord(grid); }
  }

  ///@}



  // *****************************************************************
  // TEMPLATE CLASS GRID_VERTEX_LIST MEMBER FUNCTIONS
  // *****************************************************************
  
  template <typename VTYPE>
  void GRID_VERTEX_LIST<VTYPE>::Init()
  {
    this->vertex_list = NULL;
    this->list_length = 0;
    this->num_vertices = 0;
  }

  template <typename VTYPE>
  void GRID_VERTEX_LIST<VTYPE>::FreeAll()
  {
    if (vertex_list != NULL) { delete [] vertex_list; };
    vertex_list = NULL;
    this->list_length = 0;
    num_vertices = 0;
  }

  template <typename VTYPE>
  void GRID_VERTEX_LIST<VTYPE>::AllocateList
  (const VTYPE list_length)
  {
    FreeAll();

    if (list_length > 0) {
      vertex_list = new VTYPE[list_length];
      this->list_length = list_length;
      num_vertices = 0;
    }
  }

  template <typename VTYPE>
  template<typename GCLASS>
  void FACET0_CUBE_LIST<VTYPE>::GetCubes(const GCLASS & grid)
  {
    VTYPE num_cubes;
    compute_num_cubes_in_grid_facet0
      (grid.Dimension(), grid.AxisSize(), num_cubes);

    if (num_cubes > this->ListLength()) 
      { this->AllocateList(num_cubes); }

    if (num_cubes > 0) {
      get_cubes_in_grid_facet0
        (grid.Dimension(), grid.AxisSize(), this->vertex_list);
    }

    this->num_vertices = num_cubes;
  }


  /// Get facets from grid and store in list.
  template <typename DTYPE, typename VTYPE>
  template <typename GCLASS, typename DTYPE0>
  void GRID_FACET_LIST<DTYPE,VTYPE>::GetFacetsInGridFacet
  (const GCLASS & grid, const DTYPE0 orth_dir, const bool side)
  {
    VTYPE num_facets;
    compute_num_facets_in_grid_facet
      (grid.Dimension(), grid.AxisSize(), orth_dir, num_facets);

    if (num_facets > this->ListLength()) 
      { this->AllocateList(num_facets); }

    if (num_facets > 0) {
      get_facets_in_grid_facet
        (grid.Dimension(), grid.AxisSize(), orth_dir, side,
         this->vertex_list);
    }

    this->num_vertices = num_facets;
  }


  /// GRID_BOUNDARY_VERTEX_LIST constructor.
  template <typename VTYPE>
  template<typename GCLASS>
  GRID_BOUNDARY_VERTEX_LIST<VTYPE>::GRID_BOUNDARY_VERTEX_LIST
  (const GCLASS & grid)
  {
    VTYPE numv = 0;

    compute_num_boundary_grid_vertices
      (grid.Dimension(), grid.AxisSize(), numv);

    this->AllocateList(numv);

    if (numv > 0) {
      get_boundary_grid_vertices
        (grid.Dimension(), grid.AxisSize(), this->vertex_list);
    }

    this->num_vertices = numv;
  }

  /// FACET_VERTEX_LIST constructor.
  template <typename VTYPE>
  template<typename GCLASS>
  FACET_VERTEX_LIST<VTYPE>::FACET_VERTEX_LIST
  (const GCLASS & grid, const VTYPE orth_dir,
   const bool allocate_max)
  {
    VTYPE numv = 0;

    if (allocate_max) {
      compute_max_num_vertices_in_grid_facet
        (grid.Dimension(), grid.AxisSize(), numv);
    }
    else {
      compute_num_vertices_in_grid_facet
        (grid.Dimension(), grid.AxisSize(), orth_dir, numv);
    }

    this->AllocateList(numv);
    GetVertices(grid, orth_dir);
  }

  /// Get vertices in grid facet
  template <typename VTYPE>
  template<typename GCLASS>
  void FACET_VERTEX_LIST<VTYPE>::GetVertices
  (const GCLASS & grid, const VTYPE orth_dir)
  {
    VTYPE numv;
    const bool side = false;

    compute_num_vertices_in_grid_facet
      (grid.Dimension(), grid.AxisSize(), orth_dir, numv);

    if (numv > this->ListLength()) 
      { this->AllocateList(numv); }

    if (numv > 0) {
      get_vertices_in_grid_facet
        (grid.Dimension(), grid.AxisSize(), orth_dir, side,
         this->vertex_list);
    }

    this->num_vertices = numv;
  }

  /// FACET_INTERIOR_VERTEX_LIST constructor.
  template <typename VTYPE>
  template<typename GCLASS>
  FACET_INTERIOR_VERTEX_LIST<VTYPE>::FACET_INTERIOR_VERTEX_LIST
  (const GCLASS & grid, const VTYPE orth_dir,
   const bool allocate_max, const bool flag_dim1_facet_vertex)
  {
    const VTYPE boundary_width = 1;
    VTYPE numv = 0;

    if (allocate_max) {
      compute_max_num_vertices_in_grid_facet_interior
        (grid.Dimension(), grid.AxisSize(), boundary_width, numv);
    }
    else {
      compute_num_vertices_in_grid_facet_interior
        (grid.Dimension(), grid.AxisSize(), orth_dir, boundary_width,
         numv);
    }

    if (grid.Dimension() == 1 && flag_dim1_facet_vertex)
      { numv = std::max(numv, 1); };

    this->AllocateList(numv);
    GetVertices(grid, orth_dir, flag_dim1_facet_vertex);
  }

  /// Get vertices in grid facet interior
  template <typename VTYPE>
  template<typename GCLASS>
  void FACET_INTERIOR_VERTEX_LIST<VTYPE>::GetVertices
  (const GCLASS & grid, const VTYPE orth_dir, 
   const bool flag_dim1_facet_vertex)
  {
    VTYPE numv;
    const VTYPE boundary_width = 1;
    const bool side = false;

    compute_num_vertices_in_grid_facet_interior
      (grid.Dimension(), grid.AxisSize(), orth_dir, boundary_width,
       numv);

    if (grid.Dimension() == 1 && flag_dim1_facet_vertex)
      { numv = std::max(numv, 1); };

    if (numv > this->ListLength()) 
      { this->AllocateList(numv); }

    if (numv > 0) {
      if (grid.Dimension() > 1) {
        get_vertices_in_grid_facet_interior
          (grid.Dimension(), grid.AxisSize(), orth_dir, side,
           boundary_width, this->vertex_list);
      }
      else if (grid.Dimension() == 1 && flag_dim1_facet_vertex) {
        this->vertex_list[0] = 0;
      }
    }

    this->num_vertices = numv;
  }


  // *****************************************************************
  //! @name TEMPLATE OUTPUT FUNCTIONS (deprecated)
  // *****************************************************************

  ///@{

  /// DEPRECATED. Use GRID::PrintCoord().
  /// Output coord (for debugging purposes)
  template <typename DTYPE, typename CTYPE>
  void ijkgrid_output_coord
  (std::ostream & out,
   const DTYPE dimension, const CTYPE * coord)
  {
    out << "(";
    for (DTYPE d = 0; d < dimension; d++) {
      out << coord[d];
      if (d+1 < dimension) 
        { out << ","; }
    }
    out << ")";
  }

  /// DEPRECATED. Use GRID::PrintCoord().
  /// Output vertex coord (for debugging purposes)
  template <typename GTYPE, typename VTYPE>
  void ijkgrid_output_vertex_coord
  (std::ostream & out, 
   const GTYPE & grid, const VTYPE iv)
  {
    typedef typename GTYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = grid.Dimension();
    IJK::ARRAY<VTYPE> coord(dimension);

    grid.ComputeCoord(iv, coord.Ptr());
    ijkgrid_output_coord(out, dimension, coord.PtrConst());
  }

  ///@}

}


#endif
