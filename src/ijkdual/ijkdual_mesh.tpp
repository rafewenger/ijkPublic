/*!
 *  @file ijkdual_mesh.tpp
 *  @brief ijkdual template classes for polyhedral mesh data structures.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2023 Rephael Wenger

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


#ifndef _IJKDUAL_MESH_TPP_
#define _IJKDUAL_MESH_TPP_

#include "ijk.tpp"
#include "ijkmesh_datastruct.tpp"

#include <vector>

namespace IJKDUAL {

  // *****************************************************************
  // Class ADJACENT_VERTEX_AND_DUAL_FACET
  // *****************************************************************

  /// Class storing adjacent vertex and facet dual to edge
  ///   between vertex and adjacent vertex.
  template <typename VTYPE, typename FTYPE>
  class ADJACENT_VERTEX_AND_DUAL_FACET:
    public IJK::VERTEX_ADJACENCY_LIST_ELEMENT<VTYPE> {

  public:

    typedef FTYPE FACET_INDEX_TYPE;

    /// Index of facet dual to edge connecting vertex and adjacent vertex.
    /// - In range [0..(2*dimension-1)].
    FTYPE dual_facet;

  public:

    /// Constructors.
    ADJACENT_VERTEX_AND_DUAL_FACET() {};

    template <typename FTYPE2>
    ADJACENT_VERTEX_AND_DUAL_FACET(const FTYPE2 dual_facet) 
    { this->dual_facet = dual_facet; };

    /// Set dual facet.
    template <typename FTYPE2>
    void SetDualFacet(const FTYPE2 dual_facet)
    { this->dual_facet = dual_facet; };
  };


  // *****************************************************************
  // Class ADJACENT_VERTEX_AND_DUAL_FACET_WITH_FLAG
  // *****************************************************************

  template <typename VTYPE, typename FTYPE>
  class ADJACENT_VERTEX_AND_DUAL_FACET_WITH_FLAG:
    public ADJACENT_VERTEX_AND_DUAL_FACET<VTYPE,FTYPE> {

  public:

    /// True if edge is dual to a grid facet.
    bool is_dual_to_facet;

  public:
    ADJACENT_VERTEX_AND_DUAL_FACET_WITH_FLAG() 
    { is_dual_to_facet = false; };

    /// Return true if edge is dual to a grid facet.
    bool IsDualToFacet() const
    { return(is_dual_to_facet); }

    /// Set dual facet.
    template <typename FTYPE2>
    void SetDualFacet(const FTYPE2 dual_facet) {
      this->dual_facet = dual_facet; 
      is_dual_to_facet = true;
    };

  };


  // *****************************************************************
  // Classes VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE and
  //   VERTEX_ADJACENCY_AND_DUAL_FACET_LIST.
  // *****************************************************************

  /// List of vertices adjacent to each vertex and the facets dual
  ///   to the edge connecting each vertex to an adjacent vertex.
  /// - Base class.
  /// @tparam ETYPE List element type.
  ///   Usually derived from ADJACENT_VERTEX_AND_DUAL_FACET.
  /// @tparam NTYPE Number type.
  template <typename DTYPE, typename ETYPE, typename NTYPE>
  class VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE:
    public IJK::VERTEX_ADJACENCY_LIST_BASE<ETYPE,NTYPE> {

  protected:

    /// Grid dimension.
    /// - Used in calculating facet.
    DTYPE dimension;

    void Init()
    { dimension = 3; }

  public:
    typedef typename ETYPE::FACET_INDEX_TYPE FACET_INDEX_TYPE;

  public:
    /// Constructors.
    VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE() { Init(); };

    template <typename DTYPE2>
    VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE(const DTYPE2 dimension)
    { SetDimension(dimension); }

    /// Return dimension.
    DTYPE Dimension() const
    { return(dimension); }

    /// Return index of facet dual to edge connecting iv and Adjacent(iv,j).
    /// Index is in range [0..(2*dimension-1)].
    template <typename VTYPE2, typename NTYPE2>
    FACET_INDEX_TYPE DualFacet(const VTYPE2 iv, const NTYPE2 j) const
    { return(this->element[this->ElementIndex(iv,j)].dual_facet); }

    /// Return facet.
    /// @param orientation +1 or -1.
    template <typename NTYPE2>
    FACET_INDEX_TYPE FacetIndex
    (const DTYPE direction, const NTYPE2 orientation) const
    {
      return(direction+((orientation+1)/2)*Dimension());
    }

    /// Set dimension.
    template <typename DTYPE2>
    void SetDimension(const DTYPE2 dimension)
    { this->dimension = dimension; }

    /// Set index of facet dual to edge connecting iv and Adjacent(iv,j).
    template <typename VTYPE2, typename FTYPE2, typename NTYPE2>
    void SetDualFacet(const VTYPE2 iv, const NTYPE2 j, 
                      const FTYPE2 dual_facet)
    { this->element[this->ElementIndex(iv,j)].SetDualFacet(dual_facet); }

    /// Set dual facets from a list of vertices and a list of cubes.
    /// @pre Endpoints of every edge are in different, adjacent cubes.
    /// @tparam VTYPE Must have a field .cube_list_index storing the location
    ///   in cube_list of the cube containing the vertex.
    /// @tparam GRID_CUBE_TYPE Must have a member function
    ///   GetSharedFacetIndex(cubeB) which returns the facet shared
    ///   with the cube.
    template <typename VTYPE, typename GRID_CUBE_TYPE>
    void SetAllDualFacetsFromVertexAndCubeLists
    (const std::vector<VTYPE> & vertex_list,
     const std::vector<GRID_CUBE_TYPE> & cube_list);
  };


  /// List of vertices adjacent to each vertex and the facets dual
  ///   to the edge connecting each vertex to an adjacent vertex.
  /// @tparam VTYPE Vertex index type.
  /// @tparam FTYPE Facet index type.
  /// @tparam NTYPE Number type.
  template <typename DTYPE, typename VTYPE, typename FTYPE, typename NTYPE>
  class VERTEX_ADJACENCY_AND_DUAL_FACET_LIST:
    public VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE
  <DTYPE,ADJACENT_VERTEX_AND_DUAL_FACET<VTYPE,FTYPE>,NTYPE> {

  public:
    VERTEX_ADJACENCY_AND_DUAL_FACET_LIST() {};

    template <typename DTYPE2>
    VERTEX_ADJACENCY_AND_DUAL_FACET_LIST(const DTYPE2 dimension):
      VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE
      <DTYPE,ADJACENT_VERTEX_AND_DUAL_FACET<VTYPE,FTYPE>,NTYPE>(dimension)
    {};

  };


  // *****************************************************************
  // Classes VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE and
  //   VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST.
  // *****************************************************************

  /// List of vertices adjacent to each vertex and the facets dual
  ///   to the edge connecting each vertex to an adjacent vertex.
  /// - Edges may not have dual facets.
  /// - Base class.
  /// @tparam ETYPE List element type.
  ///   Usually derived from ADJACENT_VERTEX_AND_DUAL_FACET_WITH_FLAG.
  /// @tparam NTYPE Number type.
  template <typename DTYPE, typename ETYPE, typename NTYPE>
  class VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE:
    public VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE<DTYPE,ETYPE,NTYPE> {

  public:
    typedef typename ETYPE::FACET_INDEX_TYPE FACET_INDEX_TYPE;
    typedef typename ETYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;

  public:
    /// Constructors.
    VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE(){};

    template <typename DTYPE2>
    VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE
    (const DTYPE2 dimension):
      VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE<DTYPE,ETYPE,NTYPE>(dimension)
    {};

    /// Return true if edge from iv to Adjacent(iv,j) 
    ///   is dual to some grid facet.
    template <typename VTYPE2, typename NTYPE2>
    bool IsDualToFacet(const VTYPE2 iv, const NTYPE2 j) const
    { return(this->element[this->ElementIndex(iv,j)].IsDualToFacet()); }

    /// Get adjacent vertex in oriented direction.
    template <typename VTYPE2, typename VTYPE3,
              typename DIR_TYPE, typename NTYPE2>
    bool GetAdjacentVertexInOrientedDirection
    (const VTYPE2 iv, const DIR_TYPE dir, const NTYPE2 orientation, 
     VTYPE3 & adj_vertex) const;

    /// Redefine SetAllDualFacets() to allow adjacent vertices
    ///   to be in the same cube.
    /// - If adjacent vertices are in the same cube, then they
    ///     have no dual facet.
    /// @pre Endpoints of every edge are either in the same cube or
    ///   are in adjacent cubes.
    /// @tparam VTYPE Must have a field .cube_list_index storing the location
    ///   in cube_list of the cube containing the vertex.
    /// @tparam GRID_CUBE_TYPE Must have a member function
    ///   GetSharedFacetIndex(cubeB) which returns the facet shared
    ///   with the cube.
    template <typename VTYPE, typename GRID_CUBE_TYPE>
    void SetAllDualFacetsFromVertexAndCubeLists
    (const std::vector<VTYPE> & vertex_list,
     const std::vector<GRID_CUBE_TYPE> & cube_list);
  };


  /// List of vertices adjacent to each vertex and the facets dual
  ///   to the edge connecting each vertex to an adjacent vertex.
  /// - Edges may not have dual facets.
  template <typename DTYPE,typename VTYPE, typename FTYPE, typename NTYPE>
  class VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST:
    public VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE
  <DTYPE,ADJACENT_VERTEX_AND_DUAL_FACET_WITH_FLAG<VTYPE,FTYPE>,NTYPE> {

  public:
    VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST(){};

    template <typename DTYPE2>
    VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST(const DTYPE2 dimension):
      VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE
      <DTYPE,ADJACENT_VERTEX_AND_DUAL_FACET_WITH_FLAG<VTYPE,FTYPE>,NTYPE>
      (dimension)
    {};
  };


  // *****************************************************************
  // Classes VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE 
  //   member functions.
  // *****************************************************************

  template <typename DTYPE, typename ETYPE, typename NTYPE>
  template <typename VTYPE, typename GRID_CUBE_TYPE>
  void VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE<DTYPE,ETYPE,NTYPE>::
  SetAllDualFacetsFromVertexAndCubeLists
  (const std::vector<VTYPE> & vertex_list,
   const std::vector<GRID_CUBE_TYPE> & cube_list)
  {
    IJK::PROCEDURE_ERROR error
      ("VERTEX_ADJACENCY_AND_DUAL_FACET_LIST_BASE::SetAllDualFacetsFromVertexAndCubeLists");

    if (vertex_list.size() == 0) { return; }

    if (cube_list.size() == 0) {
      error.AddMessage("Programming error. cube_list not set.");
      throw error;
    }

    for (NTYPE iv0 = 0; iv0 < this->NumVertices(); iv0++) {
      for (NTYPE j = 0; j < this->NumAdjacent(iv0); j++) {
        const NTYPE iv1 = this->AdjacentVertex(iv0, j);
        const NTYPE index_to_c0 = vertex_list[iv0].cube_list_index;
        const NTYPE index_to_c1 = vertex_list[iv1].cube_list_index;
        const NTYPE ifacet =
          cube_list[index_to_c0].GetSharedFacetIndex(cube_list[index_to_c1]);

        this->SetDualFacet(iv0, j, ifacet);
      }
    }

  }    


  // *****************************************************************
  // Classes VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE 
  //   member functions.
  // *****************************************************************

  template <typename DTYPE, typename ETYPE, typename NTYPE>
  template <typename VTYPE2, typename VTYPE3,
            typename DIR_TYPE, typename NTYPE2>
  bool VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE<DTYPE,ETYPE,NTYPE>::
  GetAdjacentVertexInOrientedDirection
  (const VTYPE2 iv, const DIR_TYPE dir, const NTYPE2 orientation, 
   VTYPE3 & adj_vertex) const
  {
    const NTYPE jfacet = this->FacetIndex(dir, orientation);

    // Initialize.
    adj_vertex = 0;

    for (NTYPE k = 0; k < this->NumAdjacent(iv); k++) {
      if (this->IsDualToFacet(iv,k)) {
        if (this->DualFacet(iv,k) == jfacet) {
          adj_vertex = this->AdjacentVertex(iv,k);
          return(true);
        }
      }
    }

    return(false);
  }


  template <typename DTYPE, typename ETYPE, typename NTYPE>
  template <typename VTYPE, typename GRID_CUBE_TYPE>
  void VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE<DTYPE,ETYPE,NTYPE>::
  SetAllDualFacetsFromVertexAndCubeLists
  (const std::vector<VTYPE> & vertex_list,
   const std::vector<GRID_CUBE_TYPE> & cube_list)
  {
    IJK::PROCEDURE_ERROR error
      ("VERTEX_ADJACENCY_AND_DUAL_FACET_WITH_FLAG_LIST_BASE::SetAllDualFacetsFromVertexAndCubeLists");

    if (vertex_list.size() == 0) { return; }

    if (cube_list.size() == 0) {
      error.AddMessage("Programming error. cube_list not set.");
      throw error;
    }

    for (NTYPE iv0 = 0; iv0 < this->NumVertices(); iv0++) {
      for (NTYPE j = 0; j < this->NumAdjacent(iv0); j++) {
        const NTYPE iv1 = this->AdjacentVertex(iv0, j);
        const NTYPE index_to_c0 = vertex_list[iv0].cube_list_index;
        const NTYPE index_to_c1 = vertex_list[iv1].cube_list_index;

        if (index_to_c0 != index_to_c1) {
          const NTYPE ifacet =
            cube_list[index_to_c0].GetSharedFacetIndex(cube_list[index_to_c1]);

          this->SetDualFacet(iv0, j, ifacet);
        }
      }
    }

  }    

}

#endif
