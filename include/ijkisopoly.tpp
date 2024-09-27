/*!
 *  @file ijkisopoly.tpp
 *  @brief ijk templates for extracting an isosurface patch 
 *    from a polytope.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009-2024 Rephael Wenger

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

#ifndef _IJKISOPOLY_
#define _IJKISOPOLY_

#include <numeric>
#include <vector>

#include "ijk.tpp"
#include "ijkbits.tpp"
#include "ijkcube.tpp"

// *** DEBUG ***
#include <iostream>

namespace IJK {

  // ***************************************************************
  //! @name SUBROUTINES FOR EXTRACTING ISOSURFACE PATCHES
  // ***************************************************************

  ///@{

  /*!
   *  @brief Return true if isosurface intersects polytope.
   *  - Version where iv0 + increment[i] is the i'th vertex
   *    of the polytope.
   *  @pre num_poly_vertices > 0.
   *  @param scalar[iv] Scalar value of vertex iv.
   *  @param iv0 Index of base for computing polytope vertex indices.
   *  @param increment[i] Increment for the (i+1)'st vertex of the polytope.
   *    - The (i+1)'st polytope vertex has index (iv0+increment[i]).
   */
  template <typename STYPE, typename STYPE2, typename VTYPE, typename ITYPE, 
            typename NTYPE>
  inline bool intersects_poly
  (const STYPE * scalar, const STYPE2 isovalue, 
   const VTYPE iv0, const ITYPE * increment, const NTYPE num_poly_vertices)

  {
    VTYPE iv1 = iv0 + increment[0];
    if (scalar[iv1] < isovalue) {
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
        iv1 = iv0 + increment[j];
        if (scalar[iv1] >= isovalue) { return(true); };
      }
    }
    else {
      // scalar[iv1] >= isovalue
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
        iv1 = iv0 + increment[j];
        if (scalar[iv1] < isovalue) { return(true); };
      }
    }

    return(false);
  }


  /*!
   *  @brief Return true if isosurface intersects polytope.
   *  - Version where scalar[i] is the scalar value at polytope vertex i.
   *  @pre num_poly_vertices > 0.
   *  @param scalar[i] Scalar value at polytope vertex i.
   */
  template <typename STYPE, typename STYPE2, typename NTYPE>
  inline bool intersects_poly
  (const STYPE * scalar, const STYPE2 isovalue, 
   const NTYPE num_poly_vertices)
  {
    if (scalar[0] < isovalue) {
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
        if (scalar[j] >= isovalue) { return(true); };
      }
    }
    else {
      // scalar[0] >= isovalue
      for (NTYPE j = 1; j < num_poly_vertices; j++) {
        if (scalar[j] < isovalue) { return(true); };
      }
    }

    return(false);
  }

  /// Compute isosurface table index of polytope with primary vertex iv0.
  template <typename STYPE, typename STYPE2, typename VTYPE, typename INC_TYPE, 
            typename NTYPE, typename ITYPE>
  inline void compute_isotable_index
  (const STYPE * scalar, const STYPE2 isovalue,
   const VTYPE iv0, const INC_TYPE * increment,
   const NTYPE num_poly_vertices, ITYPE & it)
  {
    it = 0;
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
      VTYPE iv1 = iv0 + increment[j];
      if (scalar[iv1] >= isovalue) {
        it = (it | (1L << j));
      };
    };
  }

  /// Compute isosurface table index of polytope.
  template <typename STYPE, typename STYPE2, typename NTYPE, typename ITYPE>
  void compute_isotable_index
  (const STYPE * scalar, const STYPE2 isovalue,
   const NTYPE num_poly_vertices, ITYPE & it)
  {
    it = 0;
    for (NTYPE j = 0; j < num_poly_vertices; j++) {
      if (scalar[j] >= isovalue) {
        it = (it | (1L << j));
      };
    };
  }


  /// Compute isosurface table index of grid cube icube0.
  template <typename SGRID_TYPE, typename STYPE, typename ICUBE_TYPE, 
            typename TABLE_INDEX_TYPE>
  inline void compute_isotable_index_of_grid_cube
  (const SGRID_TYPE & scalar_grid, const STYPE isovalue,
   const ICUBE_TYPE icube0, TABLE_INDEX_TYPE & table_index)
  {
    compute_isotable_index
      (scalar_grid.ScalarPtrConst(), isovalue, icube0, 
       scalar_grid.CubeVertexIncrement(), scalar_grid.NumCubeVertices(), 
       table_index);
  }


  /*!
   *  @brief For each cube, compute isosurface table index 
   *    and the number of isosurface vertices in the cube.
   *  @param[out] cube_list[] Array of cube information.
   *    - Each entry contains cube_index, table_index and
   *      number of isosurface vertices for some cube.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename SCALAR_TYPE, typename GRID_CUBE_TYPE>
  void compute_cube_isotable_info
  (const GRID_TYPE & scalar_grid, const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue, 
   std::vector<GRID_CUBE_TYPE> & cube_list)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;
 
    for (SIZE_TYPE i = 0; i < cube_list.size(); i++) {
      compute_isotable_index_of_grid_cube
        (scalar_grid, isovalue, cube_list[i].cube_index, 
         cube_list[i].table_index);

      cube_list[i].num_isov = 
        isodual_table.NumIsoVertices(cube_list[i].table_index);
    }
  }


  /*!
   *  @brief Recompute isosurface table index for each cube incident
   *    on a grid boundary ridge.
   *  - For each cube incident on a grid boundary ridge,
   *    recompute isosurface table index and number of isosurface vertices,
   *    to avoid non-manifold vertices.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename SCALAR_TYPE, typename FACET_INTERSECTION_TABLE_TYPE,
            typename GRID_CUBE_TYPE, typename NUM_TYPE>
  void recompute_cube_isotable_info_on_grid_ridges
  (const GRID_TYPE & scalar_grid, const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const FACET_INTERSECTION_TABLE_TYPE & facet_intersection_table,
   const bool flag_separate_neg,
   std::vector<GRID_CUBE_TYPE> & cube_list,
   NUM_TYPE & num_ambig_changed,
   NUM_TYPE & num_non_ambig_changed)
  {
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX;
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;
    typedef typename FACET_INTERSECTION_TABLE_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename FACET_INTERSECTION_TABLE_TYPE::VERTEX_BITS_TYPE
      VERTEX_BITS_TYPE;

    const int ONE(1);
    const int TWO(2);
    const int dimension = scalar_grid.Dimension();
    const NUMBER_TYPE num_cube_facets = scalar_grid.NumCubeFacets();
    const NUMBER_TYPE num_cube_vertices = scalar_grid.NumCubeVertices();
    IJK::PROCEDURE_ERROR error
      ("recompute_cube_isotable_info_on_grid_ridges_NEW");

    // Initialize.
    num_ambig_changed = 0;
    num_non_ambig_changed = 0;

    if (!facet_intersection_table.IsTableAllocated()) {
      error.AddMessage
        ("Programming error.  Facet intersection table is not allocated.");
      error.AddMessage
        ("  Call FACET_INTERSECTION_TABLE::Create() before passing it as an argument");
      error.AddMessage
        ("  to recompute_cube_isotable_on_grid_ridges.");
      throw error;
    }

    NUMBER_TYPE boundary_bits;
    for (NUMBER_TYPE i = 0; i < cube_list.size(); i++) {
      const VERTEX_INDEX cube_index = cube_list[i].cube_index;
      const NUMBER_TYPE num_isov = cube_list[i].num_isov;
      const TABLE_INDEX table_index = cube_list[i].table_index;
      const NUMBER_TYPE num_active_facets =
        isodual_table.NumActiveFacets(table_index);

      if ((num_isov != ONE) ||
          (num_active_facets != num_cube_facets)) {
        continue;
      }

      scalar_grid.ComputeBoundaryCubeBits(cube_index, boundary_bits);

      if (boundary_bits == 0) {
        // Cube is not on grid boundary.
        continue;
      }

      const VERTEX_BITS_TYPE vertex_ridge_bits =
        facet_intersection_table.InRidgeBits(boundary_bits);

      if (vertex_ridge_bits == 0) {
        // Cube is not on grid ridge.
        continue;
      }

      NUMBER_TYPE num_vertices_on_ridges;
      IJK::count_one_bits
        (vertex_ridge_bits, num_cube_vertices, num_vertices_on_ridges);

      if (num_vertices_on_ridges == TWO) {

        const NUMBER_TYPE iv0 =
          get_first_one_bit(vertex_ridge_bits, num_cube_vertices);
        const NUMBER_TYPE iv1 =
          get_last_one_bit(vertex_ridge_bits, num_cube_vertices);
        const bool flag_iv0 = is_bit_one(table_index, iv0);
        const bool flag_iv1 = is_bit_one(table_index, iv1);

        if (flag_iv0 != flag_iv1) { continue; }

        const TABLE_INDEX table_index2 =
          (table_index ^ vertex_ridge_bits);

        const NUMBER_TYPE num_ambig_facets =
          isodual_table.NumAmbiguousFacets(table_index);

        if (num_ambig_facets == ONE) {
          // Note: num_isov == 1; num_ambig_facets == 1;
          //       every facet is active.
          // If flag_sep_neg: Case of 3 positive vertices with one ambig facet.
          // If flag_sep_pos: Case of 3 negative vertices with one ambig facet.

          int ambig_facet_orth_dir, ambig_facet_side;
          isodual_table.GetAmbiguousFacet
            (table_index, ambig_facet_orth_dir, ambig_facet_side);

          if (facet_intersection_table.IsSomeFacetOrthogonalTo
              (boundary_bits, ambig_facet_orth_dir)) {
            cube_list[i].table_index = table_index2;
            cube_list[i].num_isov =
              isodual_table.NumIsoVertices(table_index2);
            num_ambig_changed++;
          }
        }
        else if (num_ambig_facets == 0) {
          const NUMBER_TYPE num_incident_bipolar_edges0 =
              isodual_table.ComputeNumIncidentBipolarEdges
              (table_index, iv0);
          const NUMBER_TYPE num_incident_bipolar_edges1 =
               isodual_table.ComputeNumIncidentBipolarEdges
               (table_index, iv1);

          if ((num_incident_bipolar_edges0 == ONE) &&
              (num_incident_bipolar_edges1 == ONE)) {
            // Case of four positive and four negative vertices.
            // Pos vertices form a path. Neg vertices form a path.
            // Both ridge vertices have the same sign.
            cube_list[i].table_index = table_index2;
            cube_list[i].num_isov =
              isodual_table.NumIsoVertices(table_index2);
            num_non_ambig_changed++;
          }
        }

      }
      else if (num_vertices_on_ridges == NUMBER_TYPE(dimension+1)) {

        if (isodual_table.NumAmbiguousFacets(table_index) != ONE) {
          // No change.
          continue;
        }

        const VERTEX_BITS_TYPE facet_intersection_bits =
          facet_intersection_table.InFacetIntersectionBits
          (boundary_bits);

        const NUMBER_TYPE icorner =
          get_first_one_bit(facet_intersection_bits, num_cube_vertices);

        if (icorner >= num_cube_vertices) {
          error.AddMessage
            ("Programming error. Can't find cube vertex at grid corner.");
          throw error;
        }

        // iv_interior is cube vertex opposite icorner.
        const NUMBER_TYPE iv_interior =
          (num_cube_vertices-1 - icorner);

        const bool flag_icorner = is_bit_one(table_index, icorner);
        const bool flag_iv_interior = is_bit_one(table_index, iv_interior);

        if (flag_icorner != flag_iv_interior) {
          const TABLE_INDEX new_table_index = (table_index ^ facet_intersection_bits);
          cube_list[i].table_index = new_table_index;
          cube_list[i].num_isov =
          isodual_table.NumIsoVertices(new_table_index);
          num_ambig_changed++;
        }
      }
    }

  }


  /// Add isosurface simplex vertices.
  template <typename ISOTABLE_TYPE, typename INDEX_TYPE, 
            typename VTYPE, typename ITYPE, typename STYPE, 
            typename VTYPE2, typename VTYPE3>
  inline void add_iso_simplex_vertices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv_primary, const STYPE is, 
   const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint)
  {
    for (VTYPE j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
      VTYPE jv = isotable.SimplexVertex(it, is, j);
      VTYPE je = isotable.IsosurfaceVertex(jv).Face();
      VTYPE poly_v0 = isotable.Polytope().EdgeEndpoint(je, 0);
      VTYPE poly_v1 = isotable.Polytope().EdgeEndpoint(je, 1);

      VTYPE iv0 = iv_primary + increment[poly_v0];
      VTYPE iv1 = iv_primary + increment[poly_v1];
      if (iv0 > iv1) { std::swap(iv0, iv1); };

      iso_simplices.push_back(endpoint.size()/2);
      endpoint.push_back(iv0);
      endpoint.push_back(iv1);
    };

  }

  
  /// Add isosurface simplices.
  template <typename ISOTABLE_TYPE, typename INDEX_TYPE, 
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void add_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    num_simplices = isotable.NumSimplices(it);
    for (NTYPE is = 0; is < num_simplices; is++) {
      add_iso_simplex_vertices(isotable, it, iv0, is, 
                               increment, iso_simplices, endpoint);
    };
  }

  
  /// Add isosurface simplices and reverse orientation.
  template <typename ISOTABLE_TYPE, typename INDEX_TYPE, 
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void add_reverse_orient_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    num_simplices = isotable.NumSimplices(it);
    for (NTYPE is = 0; is < num_simplices; is++) {
      add_iso_simplex_vertices(isotable, it, iv0, is, 
                               increment, iso_simplices, endpoint);

      // reverse orientation by swapping last two vertices
      NTYPE ilast = iso_simplices.size()-1;
      std::swap(iso_simplices[ilast], iso_simplices[ilast-1]);
    };
  }

  
  /// Add isosurface simplices.
  template <typename ISOTABLE_TYPE, typename INDEX_TYPE, 
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void add_iso_simplices
  (const ISOTABLE_TYPE & isotable, const INDEX_TYPE it,
   const VTYPE iv0, const ITYPE * increment,
   const bool orientation,
   std::vector<VTYPE2> & iso_simplices,
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    if (orientation) {
      add_iso_simplices(isotable, it, iv0, increment, iso_simplices,
                        endpoint, num_simplices);
    }
    else {
      add_reverse_orient_iso_simplices
        (isotable, it, iv0, increment, iso_simplices,
         endpoint, num_simplices);
    }
  }

  //@}


  // ***************************************************************
  //! @name EXTRACT ISOSURFACE PATCH FROM A POLYTOPE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Extract isosurface simplices from polytope.
   *  - Returns list of isosurface simplex vertices 
   *    and list of endpoints of grid edges containing simplex vertices.
   *  @param mesh_scalar Array of scalar values of all mesh vertices.
   *  @param isotable Isosurface lookup table.
   *  @param isovalue Isovalue.
   *  @param iv0 Index of first vertex in polytope.
   *  @param increment k'th polytope vertex has index iv0+increment[k].
   *  @param iso_simplices List of vertices of isosurface simplices.
   *  @param endpoint List of endpoints of edges containing isosurface vertices.
   *       (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
   *  @param num_simplices Number of simplices added to isosurface.
   */
  template <typename STYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void extract_isopatch_from_mesh_poly
  (const STYPE * mesh_scalar, const ISOTABLE_TYPE & isotable,
   const STYPE2 isovalue, const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polytope().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(mesh_scalar, isovalue, iv0, increment, 
                        num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(mesh_scalar, isovalue, iv0, increment, 
                             num_poly_vertices, it);

      add_iso_simplices(isotable, it, iv0, increment, 
                        iso_simplices, endpoint, num_simplices);
    }
  }


  /*!
   *  @brief Extract reverse orientation isosurface simplices
   *    from polytope.
   *  - Returns list of isosurface simplex vertices 
   *    and list of endpoints of grid edges containing simplex vertices.
   *  @param mesh_scalar Array of scalar values of all mesh vertices.
   *  @param isotable Isosurface lookup table.
   *  @param isovalue Isovalue.
   *  @param iv0 Index of first vertex in polytope.
   *  @param increment k'th polytope vertex has index iv0+increment[k].
   *  @param iso_simplices List of vertices of isosurface simplices.
   *  @param endpoint List of endpoints of edges containing isosurface vertices.
   *       (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
   *  @param num_simplices Number of simplices added to isosurface.
   */
  template <typename STYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE, typename ITYPE, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  inline void extract_isopatch_reverse_orient_from_mesh_poly
  (const STYPE * mesh_scalar, const ISOTABLE_TYPE & isotable,
   const STYPE2 isovalue, const VTYPE iv0, const ITYPE * increment,
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polytope().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(mesh_scalar, isovalue, iv0, increment, 
                        num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(mesh_scalar, isovalue, iv0, increment, 
                             num_poly_vertices, it);

      add_reverse_orient_iso_simplices
        (isotable, it, iv0, increment, 
         iso_simplices, endpoint, num_simplices);
    }
  }


  /*!
   *  @brief Extract isosurface simplices from polytope.
   *  - Returns list of isosurface simplex vertices 
   *    and list of endpoints of grid edges containing simplex vertices.
   *  @param scalar Array of polytope scalar values.
   *  @param poly_vertex Array of polytope vertex indices.
   *  @param isotable Isosurface lookup table.
   *  @param isovalue Isovalue.
   *  @param iso_simplices List of vertices of isosurface simplices.
   *  @param endpoint List of endpoints of edges containing isosurface vertices.
   *       (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
   *  @param num_simplices Number of simplices added to isosurface.
   */
  template <typename STYPE, typename VTYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE2, typename VTYPE3, typename NTYPE>
  inline void extract_isopatch_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices =
      isotable.Polytope().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, num_poly_vertices, it);

      add_iso_simplices(isotable, it, 0, poly_vertex,
                        iso_simplices, endpoint, num_simplices);
    }
  }


  /*!
   *  @brief Extract reverse orientation isosurface simplices 
   *    from polytope.
   *  - Returns list of isosurface simplex vertices 
   *    and list of endpoints of grid edges containing simplex vertices.
   *  @param scalar Array of polytope scalar values.
   *  @param poly_vertex Array of polytope vertex indices.
   *  @param isovalue Isovalue.
   *  @param iso_simplices List of vertices of isosurface simplices.
   *  @param endpoint List of endpoints of edges containing isosurface vertices.
   *       (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
   *  @param num_simplices Number of simplices added to isosurface.
   */
  template <typename STYPE, typename VTYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE2, typename VTYPE3, typename NTYPE>
  inline void extract_isopatch_reverse_orient_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    const NTYPE num_poly_vertices = isotable.Polytope().NumVertices();

    num_simplices = 0;

    typedef typename ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    // check whether cube intersects isosurface
    if (intersects_poly(scalar, isovalue, num_poly_vertices)) {

      TABLE_INDEX it;
      compute_isotable_index(scalar, isovalue, num_poly_vertices, it);

      add_reverse_orient_iso_simplices
        (isotable, it, 0, poly_vertex, iso_simplices, endpoint, num_simplices);
    }
  }


  /*!
   *  @brief Extract isosurface simplices from polytope.
   *  - Returns list of isosurface simplex vertices 
   *    and list of endpoints of grid edges containing simplex vertices.
   *  @param scalar Array of polytope scalar values.
   *  @param poly_vertex Array of polytope vertex indices.
   *  @param isotable Isosurface lookup table.
   *  @param isovalue Isovalue.
   *  @param orientation Orientation of isosurface simplices.
   *  @param iso_simplices List of vertices of isosurface simplices.
   *  @param endpoint List of endpoints of edges containing isosurface vertices.
   *       (endpoint[2*i],endpoint[2*i+1]) = endpoints of i'th edge.
   *  @param num_simplices Number of simplices added to isosurface.
   */
  template <typename STYPE, typename VTYPE, typename ISOTABLE_TYPE, typename STYPE2,
            typename VTYPE2, typename VTYPE3, typename NTYPE>
  inline void extract_isopatch_from_poly
  (const STYPE * scalar, const VTYPE * poly_vertex,
   const ISOTABLE_TYPE & isotable, const STYPE2 isovalue, 
   const bool orientation, std::vector<VTYPE2> & iso_simplices, 
   std::vector<VTYPE3> & endpoint, NTYPE & num_simplices)
  {
    if (orientation) {
      extract_isopatch_from_poly
        (scalar, poly_vertex, isotable, isovalue, 
         iso_simplices, endpoint, num_simplices);
    }
    else {
      extract_isopatch_reverse_orient_from_poly
        (scalar, poly_vertex, isotable, isovalue, 
         iso_simplices, endpoint, num_simplices);
    }
  }

  ///@}


  // ************************************************************
  // GRID_EDGE: DATA STRUCTURE
  // ************************************************************

  /*!
   *  @brief GRID_EDGE Data structure representing a grid edge.
   *  - Contains lower/leftmost edge endpoint and edge direction.
   */
  template <typename VTYPE, typename DTYPE>
  class GRID_EDGE {

  public:
    typedef VTYPE VERTEX_INDEX_TYPE;
    typedef DTYPE DIRECTION_TYPE;

  public:
    VTYPE endpoint0;
    DTYPE direction;

    /// @brief Constructor.
    GRID_EDGE() {};

    /// @brief Constructor with initialization.
    GRID_EDGE(const VTYPE end0, const VTYPE end1, const DTYPE direction)
    {
      this->endpoint0 = end0;
      // Ignore end1.
      this->direction = direction;
    }
    
    /// Return lower/leftmost endpoint.
    VTYPE Endpoint0() const { return(endpoint0); }

    /// Return edge direction.
    DTYPE Direction() const { return(direction); }

    /*!
     *  @brief Return lower/leftmost endpoint.
     *  - This member function is an alias for the Endpoint0() member function.
     *  - Functions may use GridEdgeEndpoint0() instead of Endpoint0()
     *    because the grid edge may be embedded in a larger class 
     *    with many other fields and member functions.
     */
    VTYPE GridEdgeEndpoint0() const
    { return(Endpoint0()); }

    /*!
     *  @brief Return edge direction.
     *  - This member function is an alias for the Direction() member function.
     *  - Functions may use GridEdgeDirection() instead of Direction()
     *    because the grid edge may be embedded in a larger class 
     *    with many other fields and member functions.
     */
    DTYPE GridEdgeDirection() const
    { return(Direction()); }

    /*!
     *  @brief Set grid edge endpoints and direction.
     *  @param edge_endpoint0 Lower/leftmost edge endpoint.
     *  @param edge_endpoint1 Upper/rightmost edge endpoint.
     *  - Note: edge_endpoint1 is ignored in this default set function.
     *  - It is included in the parameter list for multi-resolution cases
     *    where the dual edge may span more than one grid edge.
     *  @param edge_direction Edge direction.
     *  - Edge direction is axis parallel.
     *  - Edge direction is an integer indicating the axis parallel 
     *    to the edge.
     */
    template <typename VTYPE2, typename VTYPE3, typename DTYPE2>
    inline void Set(const VTYPE2 edge_endpoint0, const VTYPE3 edge_endpoint1,
                    const DTYPE2 edge_direction) {
      endpoint0 = edge_endpoint0;
      direction = edge_direction;
      // Note: Default GRID_EDGE set function ignores edge_endpoint1.
    }

    /*! @brief Set grid edge endpoints and direction.
     *  - This member function is an alias for the Set() member function.
     *  - Function extract_isopoly_around_dual_grid_edge_I
     *    uses SetGridEdge() instead of Set() because the grid edge
     *    may be embedded in a larger class with it's own set function.
     */
    template <typename VTYPE2, typename VTYPE3, typename DTYPE2>
    inline void SetGridEdge
    (const VTYPE2 edge_endpoint0, const VTYPE3 edge_endpoint1,
     const DTYPE2 edge_direction) 
    { Set(edge_endpoint0, edge_endpoint1, edge_direction); }

    /*!
     *  @brief Print fields in GRID_EDGE.
     *  - Mainly for debugging.
     */
    template <typename OSTREAM_TYPE, typename STYPE>
    void Print(OSTREAM_TYPE & out, const STYPE & line_prefix) const
    {
      out << line_prefix << "endpoint0: " << int(Endpoint0()) << "\n";
      out << line_prefix << "direction: " << int(Direction()) << "\n";
    }
      
  };


  // ************************************************************
  // ISODUAL_POLY_INFO: DATA STRUCTURE
  // ************************************************************

  /*!
   *  @brief ISODUAL_POLY_INFO Data structure template for
   *    dual isosurface polytope information.
   *  - Contains member dual grid edge.
   */
  template <typename GRID_EDGE_TYPE>
  class ISODUAL_POLY_INFO {

  public:
    typedef typename GRID_EDGE_TYPE::VERTEX_INDEX_TYPE 
    VERTEX_INDEX_TYPE;
    typedef typename GRID_EDGE_TYPE::DIRECTION_TYPE 
    DIRECTION_TYPE;

  public:

    /// Grid edge dual to isosurface polytope.
    GRID_EDGE_TYPE dual_grid_edge;

  public:

    /// @brief Constructor.
    ISODUAL_POLY_INFO() {}

    /// @brief Constructor.
    template <typename VTYPE, typename DTYPE>
    ISODUAL_POLY_INFO
    (const VTYPE iend0, const VTYPE iend1, const DTYPE dir):
      dual_grid_edge(iend0, iend1, dir) {};
    
    /// Return dual grid edge.
    const GRID_EDGE_TYPE & DualGridEdge() const
    { return dual_grid_edge; }


    /// @brief Return lower/leftmost endpoint of grid edge.
    VERTEX_INDEX_TYPE GridEdgeEndpoint0() const
    { return(dual_grid_edge.Endpoint0()); }

    /// @brief Return grid edge direction.
    DIRECTION_TYPE GridEdgeDirection() const
    { return(dual_grid_edge.Direction()); }

    template <typename VTYPE2, typename VTYPE3, typename DTYPE2>
    inline void SetGridEdge
    (const VTYPE2 edge_endpoint0, const VTYPE3 edge_endpoint1,
     const DTYPE2 edge_direction)
    { 
      dual_grid_edge.SetGridEdge
        (edge_endpoint0, edge_endpoint1, edge_direction);
    }

    /*!
     *  @brief Print fields in ISODUAL_POLY_INFO.
     *  - Mainly for debugging.
     */
    template <typename OSTREAM_TYPE, typename STYPE>
    void Print(OSTREAM_TYPE & out, const STYPE & line_prefix) const
    {
      const std::string line_prefix2 = 
        std::string(line_prefix) + "dual grid edge ";
      dual_grid_edge.Print(out, line_prefix2);
    };

  };


  // ***************************************************************
  //! @name EXTRACT DUAL ISOSURFACE POLYTOPE AROUND EDGE
  // ***************************************************************

  ///@{

  /*!
   *  @brief Extract dual isosurface polytope around an edge.
   *    (Single isosurface vertex per grid cube.)
   *  @details
   *  - Single isosurface vertex per grid cube.
   *  - Does not return index of dual cube edge.
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ICUBE_TYPE>
  void extract_dual_isopoly_around_edge_single
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir, std::vector<ICUBE_TYPE> & isopoly_cube)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    
    const NUM_TYPE num_facet_vertices = grid.NumCubeFacetVertices();

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    for (NUM_TYPE k = 0; k < num_facet_vertices; k++) 
      { isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k)); }
  }


  /*!
   *  @brief Extract dual isosurface polytope around an edge. 
   *    (Multiple isosurface vertices per grid cube.)
   *  @details
   *  - Each grid cube may contain multiple isosurface vertices.
   *  - Returns index of dual cube edge.
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ICUBE_TYPE, typename IEDGE_TYPE>
  void extract_dual_isopoly_around_edge_multi
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir, 
   std::vector<ICUBE_TYPE> & isopoly_cube, 
   std::vector<IEDGE_TYPE> & isopoly_cube_edge)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const NUM_TYPE num_facet_vertices = grid.NumCubeFacetVertices();

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    for (NUM_TYPE k = 0; k < num_facet_vertices; k++) {
      const NUM_TYPE k2 = num_facet_vertices-1-k;
      isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k));
      isopoly_cube_edge.push_back(edge_dir*num_facet_vertices+k2);
    }
  }


  /*!
   *  @overload
   *  @brief Extract dual isosurface polytope around an edge. 
   *    (Multiple isosurface vertices per grid cube. Return dual grid edge.)
   *  - Version that returns dual grid edge.
   *  @details
   *  - Each grid cube may contain multiple isosurface vertices.
   *  - Returns index of dual cube edge.
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param[out] dual_edge Edge dual to i'th isosurface polytope.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ICUBE_TYPE, typename IEDGE_TYPE,
            typename GRID_EDGE_TYPE>
  void extract_dual_isopoly_around_edge_multi
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir, 
   std::vector<ICUBE_TYPE> & isopoly_cube, 
   std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<GRID_EDGE_TYPE> & dual_edge)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const NUM_TYPE num_facet_vertices = grid.NumCubeFacetVertices();

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    for (NUM_TYPE k = 0; k < num_facet_vertices; k++) {
      const NUM_TYPE k2 = num_facet_vertices-1-k;
      isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k));
      isopoly_cube_edge.push_back(edge_dir*num_facet_vertices+k2);
    }

    const GRID_EDGE_TYPE dual_edge_j(iend0, iend1, edge_dir);
    dual_edge.push_back(dual_edge_j);
  }


  /*!
   *  @brief Extract dual isosurface polytope around an edge, 
   *    reverse orientation. (Single isosurface vertex per grid cube.)
   *  @details
   *  - Reverse isosurface polytope orientation.
   *  - Single isosurface vertex per grid cube.
   *  - Return index of dual cube edge.
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ICUBE_TYPE>
  void extract_dual_isopoly_around_edge_reverse_orient_single
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir, std::vector<ICUBE_TYPE> & isopoly_cube)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const NUM_TYPE num_facet_vertices = grid.NumCubeFacetVertices();
    const NUM_TYPE half_num_facet_vertices = num_facet_vertices/2;
    const VTYPE facet_vertex_increment =
      grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);
    
    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = iend0 - facet_vertex_increment;

    // Add last half_num_facet_vertices vertices first.
    for (NUM_TYPE k = half_num_facet_vertices; k < num_facet_vertices; k++) 
      { isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k)); }

    for (NUM_TYPE k = 0; k < half_num_facet_vertices; k++) 
      { isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k)); }
  }


  /*!
   *  @brief Extract dual isosurface polytope around an edge, 
   *    reverse orientation.
   *    (Multiple isosurface vertices per grid cube.)
   *  @details
   *  - Each grid cube may contain multiple isosurface vertices.
   *  - Returns index of dual cube edge.
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   */
  template <typename GRID_TYPE, typename VTYPE,
            typename DTYPE, typename ICUBE_TYPE, typename IEDGE_TYPE>
  void extract_dual_isopoly_around_edge_reverse_orient_multi
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube, 
   std::vector<IEDGE_TYPE> & isopoly_cube_edge)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const NUM_TYPE num_facet_vertices = grid.NumCubeFacetVertices();
    const NUM_TYPE half_num_facet_vertices = num_facet_vertices/2;

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    // Add last half_num_facet_vertices vertices first.
    for (NUM_TYPE k = half_num_facet_vertices; 
         k < num_facet_vertices; k++) {
      const NUM_TYPE k2 = num_facet_vertices-1-k;
      isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k));
      isopoly_cube_edge.push_back(edge_dir*num_facet_vertices+k2);      
    }

    for (NUM_TYPE k = 0; k < half_num_facet_vertices; k++) {
      const NUM_TYPE k2 = num_facet_vertices-1-k;
      isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k));
      isopoly_cube_edge.push_back(edge_dir*num_facet_vertices+k2);      
    }
  }


  /*!
   *  @overload
   *  @brief Extract dual isosurface polytope around an edge, 
   *    reverse orientation.
   *    (Multiple isosurface vertices per grid cube. Return dual grid edge.)
   *  - Version that returns dual grid edge.
   *  @details
   *  - Each grid cube may contain multiple isosurface vertices.
   *  - Returns index of dual cube edge.
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @pre Edge (iend0, iend1) is an interior edge with direction \a edge_dir.
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param[out] dual_edge Edge dual to i'th isosurface polytope.
   */
  template <typename GRID_TYPE, typename VTYPE, typename DTYPE,
            typename ICUBE_TYPE, typename IEDGE_TYPE,
            typename GRID_EDGE_TYPE>
  void extract_dual_isopoly_around_edge_reverse_orient_multi
  (const GRID_TYPE & grid, const VTYPE iend0, const VTYPE iend1, 
   const DTYPE edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube, 
   std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<GRID_EDGE_TYPE> & dual_edge)
  {
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const NUM_TYPE num_facet_vertices = grid.NumCubeFacetVertices();
    const NUM_TYPE half_num_facet_vertices = num_facet_vertices/2;

    if (num_facet_vertices == 0) { return; };

    VTYPE iv0 = 
      iend0 - grid.FacetVertexIncrement(edge_dir,num_facet_vertices-1);

    // Add last half_num_facet_vertices vertices first.
    for (NUM_TYPE k = half_num_facet_vertices; 
         k < num_facet_vertices; k++) {
      const NUM_TYPE k2 = num_facet_vertices-1-k;
      isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k));
      isopoly_cube_edge.push_back(edge_dir*num_facet_vertices+k2);
    }

    for (NUM_TYPE k = 0; k < half_num_facet_vertices; k++) {
      const NUM_TYPE k2 = num_facet_vertices-1-k;
      isopoly_cube.push_back(grid.FacetVertex(iv0, edge_dir, k));
      isopoly_cube_edge.push_back(edge_dir*num_facet_vertices+k2);
    }

    const GRID_EDGE_TYPE dual_edge_j(iend0, iend1, edge_dir);
    dual_edge.push_back(dual_edge_j);
  }

    
  /*!
   *  @brief Extract dual isosurface polytope around bipolar edge.
   *  - Extract polytope and return true if edge is bipolar.
   *  - Returns list of isosurface polytope vertices.
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @param scalar_grid scalar grid data
   *  @param isovalue isosurface scalar value
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   */
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE, typename DTYPE, typename ICUBE_TYPE>
  bool extract_dual_isopoly_around_bipolar_edge
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE iend0, const DTYPE edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube)
  {
    const VTYPE iend1 = scalar_grid.NextVertex(iend0, edge_dir);

    if (scalar_grid.Scalar(iend0) < isovalue) {
      if (scalar_grid.Scalar(iend1) >= isovalue) {
        extract_dual_isopoly_around_edge_single
          (scalar_grid, iend0, iend1, edge_dir, isopoly_cube);
        return true;
      }
    }
    else if (scalar_grid.Scalar(iend1) < isovalue) {
      if (scalar_grid.Scalar(iend0) >= isovalue) {
        extract_dual_isopoly_around_edge_reverse_orient_single
          (scalar_grid, iend0, iend1, edge_dir, isopoly_cube);
        return true;
      }
    }

    return false;
  }


  /*!
   *  @brief Extract dual isosurface polytope around bipolar edge.
   *  - Extract polytope and return true if edge is bipolar.
   *  - Returns list of isosurface polytope vertices.
   *  - Returns list of dual edges (endpoints and directions).
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @param scalar_grid scalar grid data
   *  @param isovalue isosurface scalar value
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] dual_edge Edge dual to i'th isosurface polytope.
   */
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE0, typename DIR_TYPE0, typename ICUBE_TYPE,
            typename GRID_EDGE_TYPE>
  bool extract_dual_isopoly_around_bipolar_edge_I
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE0 iend0, const DIR_TYPE0 edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube,
   GRID_EDGE_TYPE & dual_edge)
  {
    const VTYPE0 iend1 = scalar_grid.NextVertex(iend0, edge_dir);

    bool is_end0_positive = true;
    if (scalar_grid.Scalar(iend0) < isovalue) 
      { is_end0_positive = false; };

    bool is_end1_positive = true;
    if (scalar_grid.Scalar(iend1) < isovalue) 
      { is_end1_positive = false; };

    if (!is_end0_positive && is_end1_positive) {
      extract_dual_isopoly_around_edge_single
        (scalar_grid, iend0, iend1, edge_dir, isopoly_cube);
      dual_edge.SetGridEdge(iend0, iend1, edge_dir);
      return(true);
    }
    else if (is_end0_positive && !is_end1_positive) {
      extract_dual_isopoly_around_edge_reverse_orient_single
        (scalar_grid, iend0, iend1, edge_dir, isopoly_cube);
      dual_edge.SetGridEdge(iend0, iend1, edge_dir);
      return(true);
    }

    return(false);
  }


  /*!
   *  @brief Extract dual isosurface polytope around bipolar edge.
   *  - Extract polytope and return true if edge is bipolar.
   *  - Returns list of isosurface polytope vertices.
   *  - Returns list of dual edges (endpoints and directions).
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @param scalar_grid scalar grid data
   *  @param isovalue isosurface scalar value
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_info[i] Isosurface polytope i information,
   *    including grid edge dual to polytope i.
   */
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE0, typename DIR_TYPE, typename ICUBE_TYPE,
            typename ISOPOLY_INFO_TYPE>
  bool extract_dual_isopoly_around_bipolar_edge_I
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE0 iend0, const DIR_TYPE edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info)
  {
    ISOPOLY_INFO_TYPE isopoly_info_i;

    if (extract_dual_isopoly_around_bipolar_edge_I
        (scalar_grid, isovalue, iend0, edge_dir, 
         isopoly_cube, isopoly_info_i)) {
      isopoly_info.push_back(isopoly_info_i);
      return(true);
    }
    else
      { return(false); }
  }


  /*!
   *  @brief Extract dual isosurface polytope around bipolar edge.
   *  - Extract polytope and return true if edge is bipolar.
   *  - Returns list of isosurface polytope vertices.
   *  - Returns location of isosurface vertex on isosurface polytope.
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @param scalar_grid scalar grid data
   *  @param isovalue isosurface scalar value
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   */
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE, typename DTYPE,
            typename ICUBE_TYPE, typename IEDGE_TYPE>
  bool extract_dual_isopoly_around_bipolar_edge
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE iend0, const DTYPE edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube,
   std::vector<IEDGE_TYPE> & isopoly_cube_edge)
  {
    VTYPE iend1 = scalar_grid.NextVertex(iend0, edge_dir);

    bool is_end0_positive = true;
    if (scalar_grid.Scalar(iend0) < isovalue) 
      { is_end0_positive = false; };

    bool is_end1_positive = true;
    if (scalar_grid.Scalar(iend1) < isovalue) 
      { is_end1_positive = false; };

    if (is_end1_positive) {
      if (!is_end0_positive) {
        extract_dual_isopoly_around_edge_multi
          (scalar_grid, iend0, iend1, edge_dir,
           isopoly_cube, isopoly_cube_edge);
        return(true);
      }
    }
    else if (is_end0_positive) {
      if (!is_end1_positive) {
        extract_dual_isopoly_around_edge_reverse_orient_multi
          (scalar_grid, iend0, iend1, edge_dir,
           isopoly_cube, isopoly_cube_edge);
        return(true);
      }
    }

    return(false);
  }


  /*!
   *  @brief Extract dual isosurface polytope around bipolar edge.
   *  - Extract polytope and return true if edge is bipolar.
   *  - Returns list of isosurface polytope vertices.
   *  - Returns location of isosurface vertex on isosurface polytope.
   *  - Returns list of dual edges (endpoints and directions).
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @param scalar_grid scalar grid data
   *  @param isovalue isosurface scalar value
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param[out] dual_edge Edge dual to i'th isosurface polytope.
   */
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE0, typename DIR_TYPE0, 
            typename ICUBE_TYPE, typename IEDGE_TYPE,
            typename GRID_EDGE_TYPE>
  inline bool extract_dual_isopoly_around_bipolar_edge_I
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE0 iend0, const DIR_TYPE0 edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube,
   std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   GRID_EDGE_TYPE & dual_edge)
  {
    const VTYPE0 iend1 = scalar_grid.NextVertex(iend0, edge_dir);

    bool is_end0_positive = true;
    if (scalar_grid.Scalar(iend0) < isovalue) 
      { is_end0_positive = false; };

    bool is_end1_positive = true;
    if (scalar_grid.Scalar(iend1) < isovalue) 
      { is_end1_positive = false; };

    if (!is_end0_positive && is_end1_positive) {
      extract_dual_isopoly_around_edge_multi
        (scalar_grid, iend0, iend1, edge_dir,
         isopoly_cube, isopoly_cube_edge);
      dual_edge.SetGridEdge(iend0, iend1, edge_dir);
      return(true);
    }
    else if (is_end0_positive && !is_end1_positive) {
      extract_dual_isopoly_around_edge_reverse_orient_multi
        (scalar_grid, iend0, iend1, edge_dir,
         isopoly_cube, isopoly_cube_edge);
      dual_edge.SetGridEdge(iend0, iend1, edge_dir);
      return(true);
    }

    return(false);
  }


  /*!
   *  @brief Extract dual isosurface polytope around bipolar edge.
   *  - Extract polytope and return true if edge is bipolar.
   *  - Returns list of isosurface polytope vertices.
   *  - Returns list of dual edges (endpoints and directions).
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @param scalar_grid scalar grid data
   *  @param isovalue isosurface scalar value
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param[out] isopoly_info[i] Isosurface polytope i information,
   *    including grid edge dual to polytope i.
   */
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE0, typename DIR_TYPE, 
            typename ICUBE_TYPE, typename IEDGE_TYPE,
            typename ISOPOLY_INFO_TYPE>
  inline bool extract_dual_isopoly_around_bipolar_edge_I
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE0 iend0, const DIR_TYPE edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube,
   std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<ISOPOLY_INFO_TYPE> & isopoly_info)
  {
    ISOPOLY_INFO_TYPE isopoly_info_i;

    if (extract_dual_isopoly_around_bipolar_edge_I
        (scalar_grid, isovalue, iend0, edge_dir, 
         isopoly_cube, isopoly_cube_edge, isopoly_info_i)) {
      isopoly_info.push_back(isopoly_info_i);
      return(true);
    }
    else 
      { return(false); }
  }


  /*!
   *  @brief Extract dual isosurface polytope around bipolar edge.
   *  - Extract polytope and return true if edge is bipolar.
   *  - Returns list of isosurface polytope vertices.
   *  - Returns location of isosurface vertex on isosurface polytope.
   *  - Returns list of dual edges (endpoints and directions).
   *  - Note: Isosurface polytope (quad/cube) vertices are listed
   *    by increasing x, then increasing y, (then increasing z), coordinates.
   *  - Quad vertices are NOT listed in clockwise or counterclockwise
   *    around the quad. (Listing by increasing coordinates generalizes
   *    easily to higher dimensions.)
   *  @param scalar_grid scalar grid data
   *  @param isovalue isosurface scalar value
   *  @param[out] isopoly_cube[i*numv_per_isopoly+k] Index of grid cube 
   *    containing k'th vertex of i'th isosurface polytope.
   *  @param[out] isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param[out] dual_edge Edge dual to i'th isosurface polytope.
   */
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE0, typename DIR_TYPE0, 
            typename ICUBE_TYPE, typename IEDGE_TYPE,
            typename GRID_EDGE_TYPE>
  inline bool extract_dual_isopoly_around_bipolar_edge_II
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE isovalue, 
   const VTYPE0 iend0, const DIR_TYPE0 edge_dir,
   std::vector<ICUBE_TYPE> & isopoly_cube,
   std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<GRID_EDGE_TYPE> & dual_edge)
  {
    const VTYPE0 iend1 = scalar_grid.NextVertex(iend0, edge_dir);

    bool is_end0_positive = true;
    if (scalar_grid.Scalar(iend0) < isovalue) 
      { is_end0_positive = false; };

    bool is_end1_positive = true;
    if (scalar_grid.Scalar(iend1) < isovalue) 
      { is_end1_positive = false; };

    if (!is_end0_positive && is_end1_positive) {
      extract_dual_isopoly_around_edge_multi
        (scalar_grid, iend0, iend1, edge_dir,
         isopoly_cube, isopoly_cube_edge);
      GRID_EDGE_TYPE dual_edge_j;
      dual_edge_j.SetGridEdge(iend0, iend1, edge_dir);
      dual_edge.push_back(dual_edge_j);
      return(true);
    }
    else if (is_end0_positive && !is_end1_positive) {
      extract_dual_isopoly_around_edge_reverse_orient_multi
        (scalar_grid, iend0, iend1, edge_dir,
         isopoly_cube, isopoly_cube_edge);
      GRID_EDGE_TYPE dual_edge_j;
      dual_edge_j.SetGridEdge(iend0, iend1, edge_dir);
      dual_edge.push_back(dual_edge_j);
      return(true);
    }

    return(false);
  }

  
  ///@}

  // ***************************************************************
  // C++ EXTRACT_DUAL_ISOPOLY classes
  // ***************************************************************

  /// Base class for various versions of EXTRACT_DUAL_ISOPOLY.
  template <typename ITYPE0>
  class EXTRACT_DUAL_ISOPOLY_BASE {

  public:
    std::vector<ITYPE0> isopoly_cube;

  public:
    void Clear() { isopoly_cube.clear(); }
    
  protected:
    /// Function must be defined in derived classes.
    template <typename GTYPE, typename STYPE,
              typename VTYPE, typename DIR_TYPE>
    void extract_dual_isopoly
    (const GTYPE & scalar_grid, const STYPE isovalue, 
     const VTYPE iend0, const DIR_TYPE edge_dir);
  };

  
  /// Extract only cubes containing isopoly vertex.
  template <typename ITYPE0>
  class EXTRACT_DUAL_ISOPOLY_CUBE:
    public EXTRACT_DUAL_ISOPOLY_BASE<ITYPE0> {

  public:
    EXTRACT_DUAL_ISOPOLY_CUBE() {};

    /// @brief Extract cubes containing isopoly vertices.
    /// - Isopoly is dual to edge with endpoint iend0 and direction edge_dir.
    template <typename GTYPE, typename STYPE,
              typename VTYPE, typename DIR_TYPE>
    void extract_dual_isopoly
    (const GTYPE & scalar_grid, const STYPE isovalue, 
     const VTYPE iend0, const DIR_TYPE edge_dir)
    {
      extract_dual_isopoly_around_bipolar_edge
        (scalar_grid, isovalue, iend0, edge_dir, this->isopoly_cube);
    }
  };


  /*!
   *  @brief Extract cube containing isopoly vertex and index of cube edge
   *    dual to corresponding isosurface polytope.
   */
  template <typename ITYPE0, typename ITYPE1>
  class EXTRACT_DUAL_ISOPOLY_CUBE_AND_DUAL_EDGE:
    public EXTRACT_DUAL_ISOPOLY_BASE<ITYPE0> {

  public:

    /*!
     *  @brief Index of edge of cube isopoly_cube[j] dual to isosurface polytope
     *       containing isosurface vertex j.
     *    - Value is in range [0..(num_cube_edges-1)].
     */
    std::vector<ITYPE1> isopoly_cube_edge;

  public:
    EXTRACT_DUAL_ISOPOLY_CUBE_AND_DUAL_EDGE() {};

  public:
    void Clear() {
      EXTRACT_DUAL_ISOPOLY_BASE<ITYPE0>::Clear();
      isopoly_cube_edge.clear();
    }

    /*!
     *  @brief Extract cubes containing isopoly vertices and 
     *    indices of cube edges dual to corresponding 
     *    isosurface polytopes.
     *  @detailed
     *  - Isopoly is dual to edge with endpoint iend0 and direction edge_dir.
     */
    template <typename GTYPE, typename STYPE,
              typename VTYPE, typename DIR_TYPE>
    void extract_dual_isopoly
    (const GTYPE & scalar_grid, const STYPE isovalue, 
     const VTYPE iend0, const DIR_TYPE edge_dir)
    {
      extract_dual_isopoly_around_bipolar_edge
        (scalar_grid, isovalue, iend0, edge_dir,
         this->isopoly_cube, isopoly_cube_edge);
    }
  };


  // ***************************************************************
  //! @name COUNT SUBROUTINES
  // ***************************************************************

  //@{

  /*!
   *  @brief Count number of vertices with scalar value above/below isovalue.
   *  - Count vertices that are in given facet plane and 
   *    adjacent to the facet vertices.
   *  @param iv0 Index of leftmost/lowest vertex in facet.
   *  @param orth_dir Direction orthogonal to facet.
   *  @param[out] num_neg Number of vertices with scalar value below isovalue.
   *  @param[out] num_pos Number of vertices with scalar value equal to
   *    or above isovalue.
   */
  template <typename GRID_TYPE, typename SCALAR_TYPE, 
            typename VTYPE, typename ORTH_DIR_TYPE, typename NUM_TYPE>
  void count_adjacent_vertices_in_facet_plane
  (const GRID_TYPE & scalar_grid, const SCALAR_TYPE & isovalue,
   const VTYPE iv0, const ORTH_DIR_TYPE orth_dir,
   NUM_TYPE & num_neg, NUM_TYPE & num_pos)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;

    NUM_TYPE num_vertices_in_cube_ridge;
    unsigned long bit0, bit1;      // boundary bits
    unsigned long mask;

    const DTYPE dimension = scalar_grid.Dimension();
    VTYPE iv1;

    // Initialize.
    num_neg = 0;
    num_pos = 0;

    // Compute rightmost/highest vertex in facet.
    iv1 = iv0;
    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir)
        { iv1 = scalar_grid.NextVertex(iv0, d); }
    }

    num_vertices_in_cube_ridge = scalar_grid.NumCubeRidgeVertices();

    scalar_grid.ComputeBoundaryBits(iv0, bit0);
    scalar_grid.ComputeBoundaryBits(iv1, bit1);

    for (DTYPE d = 0; d < dimension; d++) {
      if (d != orth_dir) {

        mask = (1L << d);
        if ((bit0 & mask) == 0) {
          for (NUM_TYPE k = 0; k < num_vertices_in_cube_ridge; k++) {
            VTYPE iv = scalar_grid.RidgeVertex(iv0, orth_dir, d, k);
            iv = scalar_grid.PrevVertex(iv, d);

            if (scalar_grid.Scalar(iv) < isovalue)
              { num_neg++; }
            else
              { num_pos++; }
          }
        }

        mask = (1L << (d+dimension));
        if ((bit1 & mask) == 0) {
          for (NUM_TYPE k = 0; k < num_vertices_in_cube_ridge; k++) {
            VTYPE iv = scalar_grid.RidgeVertex(iv0, orth_dir, d, k);
            iv += (2*scalar_grid.AxisIncrement(d));

            if (scalar_grid.Scalar(iv) < isovalue)
              { num_neg++; }
            else
              { num_pos++; }
          }
        }
      }
    }

  }

  //@}


  // **************************************************
  //! @name SPLIT SUBROUTINES
  // **************************************************

  //@{

  /*!
   *  @brief Construct list of dual isosurface vertices from cube_list and cube_data.
   *  @param isodual_table Dual isosurface lookup table.
   *  @param cube_list[] List of cubes.
   *     of configuration for cube cube_list[i].
   */
  template <typename ISODUAL_TABLE, typename GRID_CUBE_TYPE,
            typename ISOV_TYPE>
  void construct_dual_isovert_list
  (const ISODUAL_TABLE & isodual_table,
   std::vector<GRID_CUBE_TYPE> & cube_list,
   std::vector<ISOV_TYPE> & isov_list)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    const SIZE_TYPE num_cubes = cube_list.size();

    // Count total number of isosurface vertices.
    SIZE_TYPE total_num_isov = 0;
    for (SIZE_TYPE i = 0; i < num_cubes; i++) {
      TABLE_INDEX it = cube_list[i].table_index;
      cube_list[i].num_isov = isodual_table.NumIsoVertices(it);
      total_num_isov += cube_list[i].num_isov;
    }

    isov_list.resize(total_num_isov);

    // Set isosurface vertices.
    SIZE_TYPE k = 0;
    for (SIZE_TYPE i = 0; i < num_cubes; i++) {
      cube_list[i].first_isov = k;
      TABLE_INDEX it = cube_list[i].table_index;
      SIZE_TYPE num_isov = cube_list[i].num_isov;

      for (SIZE_TYPE j = 0; j < num_isov; j++) {
        isov_list[k+j].cube_index = cube_list[i].cube_index;

        // *** TO DO: TRY SET TO 0, INSTEAD OF j? ***
        isov_list[k+j].patch_index = j;
        
        isov_list[k+j].table_index = it;
        isov_list[k+j].cube_list_index = i;
      }
      k = k+num_isov;
    }
  }


  /// Set isosurface polytope vertices.
  template <typename ISODUAL_TABLE, typename GRID_CUBE_INDEX_TYPE, 
            typename CINDEX_TYPE, typename IEDGE_TYPE, 
            typename ISOV_INDEX_TYPE>
  void set_dual_isopoly_vertices
  (const ISODUAL_TABLE & isodual_table,
   const std::vector<GRID_CUBE_INDEX_TYPE> & cube_list,
   const std::vector<CINDEX_TYPE> & isopoly_cube, 
   const std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<ISOV_INDEX_TYPE> & isopoly)
  {
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;
    typedef typename std::vector<CINDEX_TYPE>::size_type SIZE_TYPE;

    const int dimension = isodual_table.Dimension();
    IJK::CUBE_FACE_INFO<int,int,int> cube(dimension);

    isopoly.resize(isopoly_cube.size());

    for (SIZE_TYPE i = 0; i < isopoly_cube.size(); i++) {
      const ISOV_INDEX_TYPE k = isopoly_cube[i];
      const TABLE_INDEX it = cube_list[k].table_index;

      const int isopoly_cube_edge_i = isopoly_cube_edge[i];
      isopoly[i] = cube_list[k].first_isov + 
        isodual_table.IncidentIsoVertex(it, isopoly_cube_edge_i);      
    }
  }


  /// Complement a single cube table index.
  template <typename ISODUAL_TABLE, typename GRID_CUBE_TYPE,
            typename INDEX_TYPE>
  void complement_cube_table_index
  (const ISODUAL_TABLE & isodual_table, const INDEX_TYPE it0,
   GRID_CUBE_TYPE & cube_info)
  {
    cube_info.table_index = isodual_table.Complement(it0);
    cube_info.num_isov =
      isodual_table.NumIsoVertices(cube_info.table_index);
  }


  /*!
   *  @overload
   *  @brief Complement a single cube table index. (Complement cube_list[i].)
   *  - Version based on cube list.
   *  - Version that sets num_changed.
   */
  template <typename ISODUAL_TABLE, typename GRID_CUBE_TYPE,
            typename ITYPE, typename INDEX_TYPE, typename NUM_TYPE>
  void complement_cube_table_index
  (const ISODUAL_TABLE & isodual_table, const ITYPE i0, const INDEX_TYPE it0,
   std::vector<GRID_CUBE_TYPE> & cube_list, NUM_TYPE & num_changed)
  {
    complement_cube_table_index(isodual_table, it0, cube_list[i0]);
    num_changed += 1;
  }


  /// Complement two cube table index.
  template <typename ISODUAL_TABLE, typename GRID_CUBE_TYPE,
            typename ITYPE, typename INDEX0_TYPE, typename INDEX1_TYPE,
            typename NUM_TYPE>
  void complement_cube_table_indices
  (const ISODUAL_TABLE & isodual_table, 
   const ITYPE i0, const ITYPE i1, 
   const INDEX0_TYPE it0, const INDEX1_TYPE it1,
   std::vector<GRID_CUBE_TYPE> & cube_list, NUM_TYPE & num_changed)
  {
    complement_cube_table_index(isodual_table, it0, cube_list[i0]);
    complement_cube_table_index(isodual_table, it1, cube_list[i1]);
    num_changed += 2;
  }

  template <typename GRID_CUBE_TYPE, typename ILIST_TYPE>
  void set_index_to_cube_list
  (const std::vector<GRID_CUBE_TYPE> & cube_list,
   ILIST_TYPE & index_to_cube_list)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    for (SIZE_TYPE i = 0; i < cube_list.size(); i++) {
      SIZE_TYPE cube_index = cube_list[i].cube_index;
      index_to_cube_list[cube_index] = i;
    }
  }


  /*!
   *  @brief Return true if ambiguous cube facet is on grid boundary.
   *  @pre Cube has one ambiguous facet.
   *  @param[out] orth_dir Direction orthogonal to ambiguous facet.
   *  @param[out] side  Side (0 or 1) of cube containing ambiguous facet.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename CUBE_INDEX_TYPE, typename TABLE_INDEX_TYPE,
            typename DIR_TYPE, typename SIDE_TYPE>
  bool is_ambiguous_cube_facet_on_grid_boundary
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const CUBE_INDEX_TYPE cube_index,
   const TABLE_INDEX_TYPE table_index,
   DIR_TYPE & orth_dir,
   SIDE_TYPE & side)
  {
    isodual_table.GetAmbiguousFacet(table_index, orth_dir, side);
    return(grid.IsCubeFacetOnGridBoundary(cube_index, orth_dir, side));
  }


  /*!
   *  @brief Split isosurface vertex pairs which create non-manifold edges.
   *  - Cubes containing vertices have only one ambiguous facet.
   *  - Also splits isosurface vertices in cube with one ambiguous facet
   *    that lies on the grid boundary.
   *  @param isodual_table Dual isosurface lookup table.
   *    - Should include ambiguity information and functions.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename GRID_CUBE_TYPE, typename INDEX_TYPE,
            typename NTYPE>
  void split_non_manifold_isov_pairs_ambig1
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const INDEX_TYPE index_to_cube_list[],
   std::vector<GRID_CUBE_TYPE> & cube_list, 
   NTYPE & num_split)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    IJK::PROCEDURE_ERROR error("split_non_manifold_isov_pairs_ambig1");

    if (grid.Dimension() != isodual_table.Dimension()) {
      error.AddMessage
        ("Programming error.  Grid dimension does not match isotable dimension.");
      error.AddMessage("  Grid dimension: ", grid.Dimension());
      error.AddMessage("  Isotable dimension: ", isodual_table.Dimension());
      throw error;
    }

    num_split = 0;

    for (SIZE_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

      TABLE_INDEX it0 = cube_list[i0].table_index;
      NUM_TYPE num_isov0 = cube_list[i0].num_isov;

      if (num_isov0 == 1) {
        if (isodual_table.NumAmbiguousFacets(it0) == 1) {
          DTYPE orth_dir;
          NUM_TYPE side;
          const VTYPE cube_index0 = cube_list[i0].cube_index;
          if (!is_ambiguous_cube_facet_on_grid_boundary
              (grid, isodual_table, cube_index0, it0, orth_dir, side)) {

            const VTYPE cube_index1 =
              grid.AdjacentVertex(cube_index0, orth_dir, side);
            const SIZE_TYPE i1 = index_to_cube_list[cube_index1];
            TABLE_INDEX it1 = cube_list[i1].table_index;
            const NUM_TYPE num_isov1 = cube_list[i1].num_isov;
            if (num_isov1 == 1) {
              if (isodual_table.NumAmbiguousFacets(it1) == 1) {
                complement_cube_table_indices
                  (isodual_table, i0, i1, it0, it1, cube_list, num_split);
              }
            }
          }
          else {
            // Split isosurface vertices in cube_index0.
            complement_cube_table_index
              (isodual_table, i0, it0, cube_list, num_split);
          }
        }
      }
    }
  }


  /*!
   *  @brief Get interior ambiguous facet.
   *  @param[out] num_interior_ambiguous_facets
   *  @param[out] interior_ambig_facet Index of an interior ambiguous facet,
   *    if one exists.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename CTYPE, typename IT_TYPE, 
            typename NTYPE, typename FTYPE>
  void get_interior_ambiguous_facet
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const CTYPE cube_index,
   const IT_TYPE table_index,
   NTYPE & num_interior_ambiguous_facets,
   FTYPE & interior_ambig_facet)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;

    const DTYPE dimension = grid.Dimension();
    long boundary_bits;

    grid.ComputeBoundaryCubeBits(cube_index, boundary_bits);

    // Initialize
    num_interior_ambiguous_facets = 0;
    interior_ambig_facet = 0;

    for (int jf = 0; jf < grid.NumCubeFacets(); jf++) {

      if (isodual_table.IsFacetAmbiguous(table_index, jf)) {

        DTYPE orth_dir = IJK::cube_facet_orth_dir(dimension, jf);
        NUM_TYPE iside = IJK::cube_facet_side(dimension, jf);
        long mask = ((1L) << (2*orth_dir+iside));

        if ((boundary_bits&mask) == 0) {
          num_interior_ambiguous_facets++;
          interior_ambig_facet = jf;
        }
      }
    }
  }


  /*!
   *  @brief Split isosurface vertex pairs which create non-manifold edges.
   *  - Cubes containing vertices have more than one ambiguous facet,
   *    but have only ambiguous facet in the grid interior.
   *  @param isodual_table Dual isosurface lookup table.
   *    - Should include ambiguity information and functions.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename GRID_CUBE_TYPE, typename INDEX_TYPE,
            typename NTYPE>
  void split_non_manifold_isov_pairs_ambig2
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const INDEX_TYPE index_to_cube_list[],
   std::vector<GRID_CUBE_TYPE> & cube_list, 
   NTYPE & num_split)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = isodual_table.Dimension();
    IJK::PROCEDURE_ERROR error("split_non_manifold_isov_pairs_ambig2");

    if (grid.Dimension() != isodual_table.Dimension()) {
      error.AddMessage
        ("Programming error.  Grid dimension does not match isotable dimension.");
      error.AddMessage("  Grid dimension: ", grid.Dimension());
      error.AddMessage("  Isotable dimension: ", isodual_table.Dimension());
      throw error;
    }

    NUM_TYPE num_ambig_facets0, num_ambig_facets1, jfacet;

    num_split = 0;

    for (SIZE_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

      TABLE_INDEX it0 = cube_list[i0].table_index;
      NUM_TYPE num_isov0 = cube_list[i0].num_isov;

      if (num_isov0 == 1) {
        if (isodual_table.NumAmbiguousFacets(it0) > 1) {
          VTYPE cube_index0 = cube_list[i0].cube_index;

          get_interior_ambiguous_facet
            (grid, isodual_table, cube_index0, it0, num_ambig_facets0, jfacet);

          if (num_ambig_facets0 == 1) {

            DTYPE orth_dir = IJK::cube_facet_orth_dir(dimension, jfacet);
            NUM_TYPE side = IJK::cube_facet_side(dimension, jfacet);
            VTYPE cube_index1 =
              grid.AdjacentVertex(cube_index0, orth_dir, side);
            SIZE_TYPE i1 = index_to_cube_list[cube_index1];
            TABLE_INDEX it1 = cube_list[i1].table_index;
            NUM_TYPE num_isov1 = cube_list[i1].num_isov;
            if (num_isov1 == 1) {

              get_interior_ambiguous_facet
                (grid, isodual_table, cube_index1, it1, num_ambig_facets1, 
                 jfacet);

              if (num_ambig_facets1 == 1) {
                complement_cube_table_indices
                  (isodual_table, i0, i1, it0, it1, cube_list, num_split);
              }
            }
          }
        }
      }
    }

  }


  /*!
   *  @brief Split isosurface vertex pairs which create non-manifold edges.
   *  @param isodual_table Dual isosurface lookup table.
   *    - Should include ambiguity information and functions.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE,
            typename GRID_CUBE_TYPE,  typename NTYPE>
  void split_non_manifold_isov_pairs
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   std::vector<GRID_CUBE_TYPE> & cube_list, 
   NTYPE & num_split)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    const DTYPE dimension = grid.Dimension();
    const NUM_TYPE num_vertices = grid.NumVertices();
    NUM_TYPE num_split1, num_split2;

    IJK::ARRAY<SIZE_TYPE> index_to_cube_list(num_vertices);
    set_index_to_cube_list(cube_list, index_to_cube_list);

    // Initialize
    num_split1 = 0;
    num_split2 = 0;

    split_non_manifold_isov_pairs_ambig1
      (grid, isodual_table, index_to_cube_list.PtrConst(),
       cube_list, num_split1);

    if (dimension > 3) {
      split_non_manifold_isov_pairs_ambig2
        (grid, isodual_table, index_to_cube_list.PtrConst(),
         cube_list, num_split2);
    }

    num_split = num_split1+num_split2;
  }


  /*!
   *  @brief Select which cube has configuration of split isosurface vertices.
   *  - Adjacent cubes share an ambiguous facet and one will have
   *  one isosurface vertex while the other has two isosurface vertices.
   *  - Choosing the configuration improves reconstruction of sharp edges.
   *  @param isodual_table Dual isosurface lookup table.
   *  - Should include ambiguity information and functions.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename SCALAR_TYPE, 
            typename GRID_CUBE_TYPE, typename NTYPE>
  void select_split_1_2_ambig
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   std::vector<GRID_CUBE_TYPE> & cube_list, 
   NTYPE & num_changed)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = grid.Dimension();
    const NUM_TYPE num_cube_vertices = compute_num_cube_vertices(dimension);
    const NUM_TYPE num_vertices = grid.NumVertices();
    IJK::ARRAY<SIZE_TYPE> index_to_cube_list(num_vertices);

    num_changed = 0;

    // Set up index_to_cube_list.
    set_index_to_cube_list(cube_list, index_to_cube_list);

    for (SIZE_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

      TABLE_INDEX it0 = cube_list[i0].table_index;
      NUM_TYPE num_isov0 = isodual_table.NumIsoVertices(it0);

      if (isodual_table.NumAmbiguousFacets(it0) == 1) {
        DTYPE orth_dir;
        NUM_TYPE side;
        const VTYPE cube_index0 = cube_list[i0].cube_index;

        // Only apply to interior facets.
        if (!is_ambiguous_cube_facet_on_grid_boundary
            (grid, isodual_table, cube_index0, it0, orth_dir, side)) {

          // Only process facets on top/right to avoid processing a facet twice.
          if (side) {
            const VTYPE cube_index1 =
              grid.AdjacentVertex(cube_index0, orth_dir, side);
            const SIZE_TYPE i1 = index_to_cube_list[cube_index1];
            TABLE_INDEX it1 = cube_list[i1].table_index;
            const NUM_TYPE num_isov1 = isodual_table.NumIsoVertices(it1);

            // if (num_isov0 == 1 and num_isov1 == 2) or 
            //    (num_isov1 == 1 and num_isov0 == 2) or 
            if (num_isov0*num_isov1 == 2) {

              if (isodual_table.NumAmbiguousFacets(it1) == 1) {

                NUM_TYPE num_neg_in_plane, num_pos_in_plane;
                const VTYPE iv0 = 
                  grid.AdjacentVertex(cube_index0, orth_dir, side);

                count_adjacent_vertices_in_facet_plane
                  (grid, isovalue, iv0, orth_dir, 
                   num_neg_in_plane, num_pos_in_plane);

                NUM_TYPE num_neg_cube0, num_pos_cube0;
                NUM_TYPE num_neg_cube1, num_pos_cube1;
                count_bits
                  (it0, num_cube_vertices, num_neg_cube0, num_pos_cube0);
                count_bits
                  (it1, num_cube_vertices, num_neg_cube1, num_pos_cube1);

                if (num_neg_cube0 > num_neg_cube1) {
                  if (num_pos_in_plane > num_neg_in_plane) {
                    if (num_isov0 == 1) {
                      complement_cube_table_indices
                        (isodual_table, i0, i1, it0, it1, 
                         cube_list, num_changed);
                    }
                  }
                  else if (num_pos_in_plane < num_neg_in_plane) {
                    if (num_isov1 == 1) {
                      complement_cube_table_indices
                        (isodual_table, i0, i1, it0, it1, 
                         cube_list, num_changed);
                    }
                  }
                }
                else if (num_neg_cube0 < num_neg_cube1) {
                  if (num_pos_in_plane > num_neg_in_plane) {
                    if (num_isov0 == 1) {
                      complement_cube_table_indices
                        (isodual_table, i0, i1, it0, it1, 
                         cube_list, num_changed);
                    }
                  }
                  else if (num_pos_in_plane < num_neg_in_plane) {
                    complement_cube_table_indices
                      (isodual_table, i0, i1, it0, it1, 
                       cube_list, num_changed);
                  }
                }
              }
            }
          }
        }
      }
    }

  }


  /*!
   *  @brief Select configuration of ambiguous cubes to increase connectivity of the isosurface.
   *  @param isodual_table Dual isosurface lookup table.
   *    - Should include ambiguity information and functions.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename SCALAR_TYPE, 
            typename GRID_CUBE_TYPE, typename NTYPE>
  void select_ambig_to_connect_isosurface
  (const GRID_TYPE & grid,
   const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   std::vector<GRID_CUBE_TYPE> & cube_list, 
   NTYPE & num_changed)
  {
    typedef typename GRID_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename GRID_TYPE::NUMBER_TYPE NUM_TYPE;
    typedef typename GRID_TYPE::VERTEX_INDEX_TYPE VTYPE;
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    const DTYPE dimension = grid.Dimension();
    const NUM_TYPE num_vertices = grid.NumVertices();
    IJK::ARRAY<SIZE_TYPE> index_to_cube_list(num_vertices);

    num_changed = 0;

    // Set up index_to_cube_list.
    set_index_to_cube_list(cube_list, index_to_cube_list);

    for (SIZE_TYPE i0 = 0; i0 < cube_list.size(); i0++) {

      TABLE_INDEX it0 = cube_list[i0].table_index;
      NUM_TYPE num_isov0 = isodual_table.NumIsoVertices(it0);

      if (isodual_table.NumAmbiguousFacets(it0) > 1 &&
          num_isov0 > 1) {

        VTYPE cube_index0 = cube_list[i0].cube_index;

        // Check if possible to complement i0 and adjacent cubes.
        bool flag_complement = true;
        for (int kf = 0; kf < grid.NumCubeFacets(); kf++) {

          if (isodual_table.IsFacetAmbiguous(it0, kf)) {

            DTYPE orth_dir = IJK::cube_facet_orth_dir(dimension, kf);
            NUM_TYPE side = IJK::cube_facet_side(dimension, kf);

            if (!grid.IsCubeFacetOnGridBoundary(cube_index0, orth_dir, side)) {

              VTYPE cube_index1 = 
                grid.AdjacentVertex(cube_index0, orth_dir, side);
              SIZE_TYPE i1 = index_to_cube_list[cube_index1];
              TABLE_INDEX it1 = cube_list[i1].table_index;
              NUM_TYPE num_isov1 = isodual_table.NumIsoVertices(it1);

              if (isodual_table.NumAmbiguousFacets(it1) > 1) {
                flag_complement = false;
                break;
              }
              else if (isodual_table.NumAmbiguousFacets(it1) == 1 &&
                       num_isov1 != 2) {
                flag_complement = false;
                break;
              }
            }
          }
        }

        if (flag_complement) {

          // Complement i0 and adjacent cubes.
          complement_cube_table_index
            (isodual_table, i0, it0, cube_list, num_changed);

          for (int kf = 0; kf < grid.NumCubeFacets(); kf++) {

            if (isodual_table.IsFacetAmbiguous(it0, kf)) {

              DTYPE orth_dir = IJK::cube_facet_orth_dir(dimension, kf);
              NUM_TYPE side = IJK::cube_facet_side(dimension, kf);

              if (!grid.IsCubeFacetOnGridBoundary
                  (cube_index0, orth_dir, side)) {

                VTYPE cube_index1 = 
                  grid.AdjacentVertex(cube_index0, orth_dir, side);
                SIZE_TYPE i1 = index_to_cube_list[cube_index1];
                TABLE_INDEX it1 = cube_list[i1].table_index;
                NUM_TYPE num_isov1 = isodual_table.NumIsoVertices(it1);

                if (isodual_table.NumAmbiguousFacets(it1) == 1 &&
                    num_isov1 == 2) {
                  complement_cube_table_index
                    (isodual_table, i1, it1, cube_list, num_changed);
                }
              }
            }
          }

        }
      }
    }

  }


  /*!
   *  @brief Compute number of cubes which have more than one isosurface vertex.
   *  - *** DOES NOT USE isodual_table. ***
   *  - *** SHOULD BE NAMED compute_num_split ***
   */
  template <typename ISODUAL_TABLE, typename GRID_CUBE_TYPE,
            typename NTYPE>
  void compute_num_splitB
  (const ISODUAL_TABLE & isodual_table,
   const std::vector<GRID_CUBE_TYPE> & cube_list,
   NTYPE & num_split)
  {
    typedef typename std::vector<GRID_CUBE_TYPE>::size_type SIZE_TYPE;

    num_split = 0;
    for (SIZE_TYPE i = 0; i < cube_list.size(); i++) {
      const NTYPE num_isov = cube_list[i].num_isov;
      if (num_isov > 1) { num_split++; }
    }
  }                              

  //@}


  // ***************************************************************
  // SPLIT DUAL ISOSURFACE VERTICES: DATA STRUCTURES
  // ***************************************************************

  /*!
   *  @brief Dual isosurface vertex information.
   *  - Grid cube index of cube containing the isosurface, 
   *    index to cube in list of active cubes containing isosurface,
   *    isosurface patch containing the isosurface and 
   *    index in isosurface lookup table corresponding 
   *    to '+' and '-' configuration in grid cube.
   */
  template <typename CI_TYPE, typename PI_TYPE, typename TI_TYPE,
            typename CUBE_LIST_INDEX_TYPE>
  class DUAL_ISOVERT {

  public:
    typedef CI_TYPE CUBE_INDEX_TYPE;
    typedef PI_TYPE PATCH_INDEX_TYPE;
    typedef TI_TYPE TABLE_INDEX_TYPE;

  public:
    /// Index of grid cube containing isosurface vertex.
    CI_TYPE cube_index;

    /// Index of isosurface vertex (patch) within cube.
    PI_TYPE patch_index;

    /// Index of the configuration of the '+' and '-' cube vertices in the isosurface lookup table.
    TI_TYPE table_index;

    /*!
     *  @brief Index in active_cube_list[] of cube containing isosurface vertex.
     *  - active_cube_list[] is a list (array) of active cubes
     *    containing isosurface vertices.
     */
    CUBE_LIST_INDEX_TYPE cube_list_index;

    /// Return index of grid cube containing isosurface vertex.
    CI_TYPE GridCubeIndex() const
    { return cube_index; }

    /// Return index in active_cube_list[] of cube containing isosurface vertex.
    CUBE_LIST_INDEX_TYPE CubeListIndex() const
    { return cube_list_index; }
  };


  /*!
   *  @brief Return grid cube index of isovert.
   *  - Version for class T with member function GridCubeIndex().
   *  - Class T is usually, but not necessarily, derived from
   *    class DUAL_ISOVERT.
   */
  template <typename T>
  typename T::CUBE_INDEX_TYPE get_isovert_grid_cube_index
    (const T & isov,
     typename std::enable_if<std::is_class<T>::value, bool>::type = true)
  { return isov.GridCubeIndex(); }


  /*!
   *  @brief Return grid cube index of isovert.
   *  - Version where isov is the actual grid cube index,
   *  - Type T is integral (unsigned/signed short, int, long, etc.)
   */
  template <typename T>
  T get_isovert_grid_cube_index
    (const T & isov,
     typename std::enable_if<std::is_integral<T>::value, bool>::type = true)
  { return isov; }


  /// Dual isosurface vertex with boundary flag.
  template <typename CI_TYPE, typename PI_TYPE, typename TI_TYPE,
            typename CUBE_LIST_INDEX_TYPE>
  class DUAL_ISOVERT_BFLAG:
    public DUAL_ISOVERT<CI_TYPE,PI_TYPE,TI_TYPE,CUBE_LIST_INDEX_TYPE> {

  public:
    /// If true, then isosurface vertex is on isosurface boundary.
    bool flag_surface_boundary;
  };


  // ***************************************************************
  // INFORMATION ABOUT ISOSURFACE VERTICES IN A GRID CUBE 
  //    (DUAL CONTOURING)
  // ***************************************************************

  /// @brief Information about dual contouring isosurface vertices in a grid cube.
  template <typename CI_TYPE, typename TI_TYPE, 
            typename NTYPE, typename FTYPE>
  class GRID_CUBE_ISOVERT {

  public:
    typedef CI_TYPE CUBE_INDEX_TYPE;
    typedef TI_TYPE TABLE_INDEX_TYPE;
    typedef NTYPE NUMBER_TYPE;

  public:
    /*!
     *  @brief Cube index in grid.
     *  - Usually equal to the index of the lower, leftmost vertex
     *    in the cube.
     */
    CI_TYPE cube_index;

    /*!
     *  @brief Index into isosurface lookup table.
     *  - Index of the configuration of the '+' and '-' cube vertices
     *    in the isosurface lookup table.
     */
    TI_TYPE table_index;

    /// Number of isosurface vertices contained in the cube.
    NTYPE num_isov;

    /*!
     *  @brief Index of first isosurface vertex in cube.
     *  - Grid cube isosurface vertices are index first_isov
     *    to (first_isov+num_isov-1).
     */
    FTYPE first_isov;

    /// Set the cube index to cube_index.
    template <typename CI_TYPE2>
    void SetCubeIndex(const CI_TYPE2 cube_index)
    { this->cube_index = cube_index; }

    /// Return cube_index.
    CI_TYPE CubeIndex() const
    { return cube_index; }

  };


  /*!
   *  @brief Set cube indices in cube_listB.
   *  - Resize cube_listB[] to size of cube_listA[].
   *  @param cube_listA[] Array of cube indices.
   *  @param cube_listB[] Array of cubes with member function SetCubeIndex().
   */
  template <typename CTYPE, typename GRID_CUBE_TYPE>
  void set_grid_cube_indices
  (const std::vector<CTYPE> & cube_listA,
   std::vector<GRID_CUBE_TYPE> & cube_listB)
  {
    typedef typename std::vector<CTYPE>::size_type SIZE_TYPE;

    cube_listB.resize(cube_listA.size());

    for (SIZE_TYPE i = 0; i < cube_listA.size(); i++) {
      cube_listB[i].SetCubeIndex(cube_listA[i]);
    }
  };


  /// GRID_CUBE_ISOVERT and cube_coord.
  template <const int GRID_DIMENSION, typename CI_TYPE, typename TI_TYPE, 
            typename NTYPE, typename FTYPE, typename COORD_TYPE>
  class GRID_CUBE_ISOVERT_AND_CUBE_COORD:public
  GRID_CUBE_ISOVERT<CI_TYPE,TI_TYPE,NTYPE,FTYPE> {

  public:
    typedef COORD_TYPE CUBE_COORD_TYPE;

  protected:
    COORD_TYPE cube_coord[GRID_DIMENSION];

  public:

    /// @brief Set cube coord.
    /// @pre cube_index is already set.
    template <typename GRID_TYPE>
    void SetCoord(const GRID_TYPE & grid)
    { grid.ComputeCoord(this->cube_index, cube_coord); };

    /// Return cube coordinates.
    const COORD_TYPE * CubeCoord() const
    { return(cube_coord); }

    /// Return d'th cube coordinate.
    COORD_TYPE CubeCoord(const int d) const
    { return(cube_coord[d]); }

    /// Return index of facet shared by this->cube and adjacent_cube.
    template <typename GRID_CUBE_TYPE>
    CI_TYPE GetSharedFacetIndex
    (const GRID_CUBE_TYPE & adjacent_cube) const;
  };


  template <const int GRID_DIMENSION, typename CI_TYPE, typename TI_TYPE, 
            typename NTYPE, typename FTYPE, typename COORD_TYPE>
  template <typename GRID_CUBE_TYPE>
  CI_TYPE GRID_CUBE_ISOVERT_AND_CUBE_COORD
  <GRID_DIMENSION, CI_TYPE, TI_TYPE, NTYPE, FTYPE, COORD_TYPE>::
  GetSharedFacetIndex(const GRID_CUBE_TYPE & adjacent_cube) const
  {
    for (NTYPE d = 0; d < GRID_DIMENSION; d++) {
      if (this->CubeCoord(d)+1 == adjacent_cube.CubeCoord(d)) {
        return(d+GRID_DIMENSION);
      }
      else if (this->CubeCoord(d) == adjacent_cube.CubeCoord(d)+1) {
        return(d);
      }
    }

    // Error.  Cubes are not adjacent or are identical.
    // Throw error message.

    IJK::PROCEDURE_ERROR error
      ("GRID_CUBE_ISOVERT_AND_CUBE_COORD::GetSharedFacetIndex");
    for (NTYPE d = 0; d < GRID_DIMENSION; d++) {
      if (this->CubeCoord(d) < adjacent_cube.CubeCoord(d)) {
        error.AddMessage("Programming error.  Cubes are not adjacent.");
        error.AddMessage("  cube coord[", d, "] = ", this->CubeCoord(d), ".");
        error.AddMessage
          ("  adjacent cube coord[", d, "] = ", 
           adjacent_cube.CubeCoord(d), ".");
        throw error;
      }
    }

    error.AddMessage("Programming error.  Cubes have identical coordinates.");
    throw error;
  }


  // ***************************************************************
  //! @name SPLIT DUAL ISOSURFACE VERTICES: ROUTINES
  // ***************************************************************

  //@{

  /*!
   *  @brief Split dual isosurface vertices.
   *  @param isopoly_cube_loc[i*num_polyv+k] 
   *     Location in active_cube_list[] of cube containing 
   *     the k'th vertex of isosurface polytope i.
   *  @param isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param active_cube_list[] List of cubes containing isosurface vertices.
   *  @pre active_cube_list[i].table_index equals table index 
   *    of cube active_cube_ist[i].cube_index.
   *  @param isov_list[] List of isosurface vertices.
   *    - Isosurface patch in isov_list[i] is isosurface patch
   *      in cube isov_list[i].cube_index containing isosurface vertex i.
   */
  template <typename ISODUAL_TABLE, typename CINDEX_TYPE, typename IEDGE_TYPE,
            typename GRID_CUBE_TYPE, 
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE,  typename NTYPE>
  void split_dual_isovert
  (const ISODUAL_TABLE & isodual_table,
   const std::vector<CINDEX_TYPE> & isopoly_cube_loc, 
   const std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<GRID_CUBE_TYPE> & active_cube_list,
   std::vector<ISOV_TYPE> & isov_list, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split)
  {
    construct_dual_isovert_list
      (isodual_table, active_cube_list, isov_list);

    set_dual_isopoly_vertices
      (isodual_table, active_cube_list, isopoly_cube_loc, 
       isopoly_cube_edge, isopoly);

    compute_num_splitB(isodual_table, active_cube_list, num_split);
  }


  /*!
   *  @brief Split dual isosurface vertices.
   *  - Version that computes active_cube_list[] from scalar_grid and isovalue.
   *  - active_cube_list[] contains a list of cubes and associated information.
   *  @param isopoly_cube_loc[] isopoly_cube_loc[j*num_polyv+k] is the cube containing
   *    the k'th vertex of isosurface polytope j
   *  @param isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param active_cube_list[] List of cubes containing isosurface vertices.
   *    - Each entry contains a cube index, table index and number
   *      of vertices per cube.
   *  @param isov_list[] List of isosurface vertices.
   *     isosurface patch in isov_list_cube[i] containing vertex i.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename SCALAR_TYPE,
            typename CINDEX_TYPE, typename IEDGE_TYPE,
            typename GRID_CUBE_TYPE, 
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE,  typename NTYPE>
  void split_dual_isovert
  (const GRID_TYPE & scalar_grid, const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE> & isopoly_cube_loc, 
   const std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<GRID_CUBE_TYPE> & active_cube_list,
   std::vector<ISOV_TYPE> & isov_list, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split)
  {
    compute_cube_isotable_info
      (scalar_grid, isodual_table, isovalue, active_cube_list);

    split_dual_isovert
      (isodual_table, isopoly_cube_loc, isopoly_cube_edge, active_cube_list, 
       isov_list, isopoly, num_split);
  }


  /*!
   *  @brief Split dual isosurface vertices.
   *  - Split non-manifold isosurface vertex pairs.
   *  @param isodual_table Dual isosurface lookup table.
   *         Should include ambiguity information and functions.
   *  @param isopoly_cube_loc[] isopoly_cube_loc[j*num_polyv+k] is the cube containing
   *     the k'th vertex of isosurface polytope j
   *  @param isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param active_cube_list[] List of cubes containing isosurface vertices.
   *  @param isov_list[] List of isosurface vertices.
   *     isosurface patch in isov_list_cube[i] containing vertex i.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename SCALAR_TYPE,
            typename CINDEX_TYPE, typename IEDGE_TYPE,
            typename GRID_CUBE_TYPE, 
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_manifold
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE> & isopoly_cube_loc, 
   const std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<GRID_CUBE_TYPE> & active_cube_list, 
   std::vector<ISOV_TYPE> & isov_list, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_non_manifold_split)
  {
    compute_cube_isotable_info
      (scalar_grid, isodual_table, isovalue, active_cube_list);

    split_non_manifold_isov_pairs
      (scalar_grid, isodual_table, active_cube_list, num_non_manifold_split);

    split_dual_isovert
      (isodual_table, isopoly_cube_loc, isopoly_cube_edge, active_cube_list,
       isov_list, isopoly, num_split);
  }


  /*!
   *  @brief Split dual isosurface vertices.
   *  - Select which cube has configuration of split isosurface vertices
   *    where adjacent cubes share an ambiguous facet and one will have
   *    one isosurface vertex while the other has two isosurface vertices.
   *  @param isodual_table Dual isosurface lookup table.
   *         Should include ambiguity information and functions.
   *  @param active_cube_list[] List of cubes containing isosurface vertices.
   *  @param isopoly_cube_loc[] isopoly_cube_loc[j*num_polyv+k] is the cube containing
   *     the k'th vertex of isosurface polytope j
   *  @param isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param isov_list[] List of isosurface vertices.
   *     isosurface patch in isov_list_cube[i] containing vertex i.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename SCALAR_TYPE,
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename IEDGE_TYPE,
            typename TABLE_INDEX_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_select_split
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & active_cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube_loc, 
   const std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<TABLE_INDEX_TYPE> & table_index,
   std::vector<ISOV_TYPE> & isov_list, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_1_2_change)
  {
    table_index.resize(active_cube_list.size());

    compute_cube_isotable_index
      (scalar_grid, isodual_table, isovalue, active_cube_list, table_index);

    select_split_1_2_ambig
      (scalar_grid, isodual_table, isovalue, active_cube_list,
       table_index, num_1_2_change);

    split_dual_isovert
      (isodual_table, active_cube_list, table_index, isopoly_cube_loc, isopoly_cube_edge, 
       isov_list, isopoly, num_split);
  }


  /*!
   *  @brief Split dual isosurface vertices.
   *  - Select which cube has configuration of split isosurface vertices
   *    where adjacent cubes share an ambiguous facet and one will have
   *    one isosurface vertex while the other has two isosurface vertices.
   *  @param isodual_table Dual isosurface lookup table.
   *         Should include ambiguity information and functions.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename SCALAR_TYPE, 
            typename CINDEX_TYPE0, typename CINDEX_TYPE1, typename IEDGE_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_select_split
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE0> & active_cube_list, 
   const std::vector<CINDEX_TYPE1> & isopoly_cube_loc, 
   const std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<ISOV_TYPE> & isov_list, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_1_2_change)
  {
    typedef typename ISODUAL_TABLE::TABLE_INDEX TABLE_INDEX;

    std::vector<TABLE_INDEX> table_index;

    split_dual_isovert_select_split
      (scalar_grid, isodual_table, isovalue, active_cube_list, 
       isopoly_cube_loc, isopoly_cube_edge, table_index, isov_list, isopoly,
       num_split, num_1_2_change);
  }


  /*!
   *  @brief Split dual isosurface vertices with special processing of ambiguous facets.
   *  - Calls both split_non_manifold_isov_pairs and select_split_1_2_ambig.
   *  @param isodual_table Dual isosurface lookup table.
   *         Should include ambiguity information and functions.
   *  @param active_cube_list[] List of cubes containing isosurface vertices.
   *  @param isopoly_cube_loc[] isopoly_cube_loc[j*num_polyv+k] is the cube containing
   *     the k'th vertex of isosurface polytope j
   *  @param isopoly_cube_edge[j]
   *    Index of edge of cube isopoly_cube[j] dual to isosurface polytope
   *       containing isosurface vertex j.
   *    - Value is in range [0..(num_cube_edges-1)].
   *  @param isov_list[] List of isosurface vertices.
   *     isosurface patch in isov_list_cube[i] containing vertex i.
   */
  template <typename GRID_TYPE, typename ISODUAL_TABLE, typename SCALAR_TYPE, 
            typename CINDEX_TYPE, typename IEDGE_TYPE,
            typename GRID_CUBE_TYPE,
            typename ISOV_TYPE, typename ISOV_INDEX_TYPE, 
            typename NTYPE>
  void split_dual_isovert_ambig
  (const GRID_TYPE & scalar_grid, 
   const ISODUAL_TABLE & isodual_table,
   const SCALAR_TYPE isovalue,
   const std::vector<CINDEX_TYPE> & isopoly_cube_loc, 
   const std::vector<IEDGE_TYPE> & isopoly_cube_edge,
   std::vector<GRID_CUBE_TYPE> & active_cube_list, 
   std::vector<ISOV_TYPE> & isov_list, 
   std::vector<ISOV_INDEX_TYPE> & isopoly,
   NTYPE & num_split,
   NTYPE & num_non_manifold_split,
   NTYPE & num_1_2_change)
  {
    compute_cube_isotable_info
      (scalar_grid, isodual_table, isovalue, active_cube_list);

    split_non_manifold_isov_pairs
      (scalar_grid, isodual_table, active_cube_list, num_non_manifold_split);

    select_split_1_2_ambig
      (scalar_grid, isodual_table, isovalue, active_cube_list, num_1_2_change);

    split_dual_isovert
      (isodual_table, isopoly_cube_loc, isopoly_cube_edge, active_cube_list, 
       isov_list, isopoly, num_split);
  }

  //@}


  // ***************************************************************
  //! @name CHECK FUNCTIONS
  // ***************************************************************

  //@{

  /*!
   *  @brief Return true if number of vertices in array isopoly_vert[]
   *    equals (num_vertices_per_isopoly)*(num_isopoly).
   *  - Otherwise, set error message and return false.
   */
  template <typename VTYPE, typename NTYPEV, typename NTYPEP>
  bool check_isopoly_vert_size
  (const std::vector<VTYPE> & isopoly_vert, const NTYPEV num_vertices_per_isopoly,
   const NTYPEP num_isopoly, IJK::ERROR & error)
  {
    if (isopoly_vert.size() != num_isopoly*num_vertices_per_isopoly) {
      error.AddMessage
        ("Programming error. Incorrect size of array isopoly_vert[].");
      error.AddMessage("  isopoly_vert.size(): ", isopoly_vert.size());
      error.AddMessage("  Number of isosurface polytopes: ", num_isopoly);
      error.AddMessage
        ("  isopoly_vert.size() should be ", num_vertices_per_isopoly,
         "*(# isosurface polytopes) = ",
         num_vertices_per_isopoly*num_isopoly, ".");
      return false;
    }

    return true;
  }
  
  //@}
  

  // ***************************************************************
  //! @name DETERMINE BOUNDARY DUAL ISOVERT
  // ***************************************************************

  //@{

  /*!
   *  @brief Determine vertices on dual contouring isosurface boundary.
   *  @pre: Vertices of each isosurface polytope are distinct.
   */
  template <typename ISOTABLE_TYPE, typename ISOV_INDEX_TYPE,
            typename NTYPE, typename DUAL_ISOV_TYPE>
  void determine_dual_contouring_boundary_isovert
  (const ISOTABLE_TYPE & dual_table,
   const std::vector<ISOV_INDEX_TYPE> & poly_vert,
   const NTYPE num_poly_vertices,
   std::vector<DUAL_ISOV_TYPE> & isov_list)
  {
    typedef typename std::vector<ISOV_INDEX_TYPE>::size_type SIZE_TYPE;

    /// num_poly_incident_on_isov[j] = Number of polytopes incident
    ///   on isovert j.
    std::vector<NTYPE> num_poly_incident_on_isov(isov_list.size(), 0);

    for (SIZE_TYPE i = 0; i < poly_vert.size(); i++) {
      ISOV_INDEX_TYPE isov = poly_vert[i];
      num_poly_incident_on_isov[isov]++;
    }

    for (SIZE_TYPE isov = 0; isov < isov_list.size(); isov++) {
      const NTYPE degree = 
        dual_table.Degree(isov_list[isov].table_index,
                          isov_list[isov].patch_index);

      if (degree == num_poly_incident_on_isov[isov]) {
        isov_list[isov].flag_surface_boundary = false;
      }
      else {
        isov_list[isov].flag_surface_boundary = true;
      }
    }
  }

  //@}


  // ***************************************************************
  //! @name QUERY FUNCTIONS
  // ***************************************************************

  //@{

  /// Return true if some isosurface vertex in vlist[] is on surface boundary.
  template <typename DUAL_ISOV_TYPE, typename ISOV_INDEX_TYPE,
            typename NTYPE>
  bool is_some_isovert_on_surface_boundary
  (const std::vector<DUAL_ISOV_TYPE> & isov_list,
   const ISOV_INDEX_TYPE vlist[], const NTYPE list_length) 
  {
    for (NTYPE i = 0; i < list_length; i++) {
      if (isov_list[vlist[i]].flag_surface_boundary)
        { return(true); }
    }
    
    return(false);
  }

  //@}

}


#endif
