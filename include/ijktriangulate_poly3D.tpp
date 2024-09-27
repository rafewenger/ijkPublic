/*!
 *  @file ijktriangulate_poly3D.tpp
 *  @brief ijk templates for triangulating 3D polytopes
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2010-2019 Rephael Wenger

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

#ifndef _IJKTRIANGULATE_POLY3D_
#define _IJKTRIANGULATE_POLY3D_

#include "ijkcube.tpp"
#include "ijkmesh.tpp"
#include "ijkmesh_datastruct.tpp"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>


namespace IJK {


  // **************************************************
  // @name CLASS HEX_TRIANGULATION_INFO
  // **************************************************

  ///@{

  /// Information on hex triangulation.
  template <typename ANCHOR_TYPE, typename FLAGS_TYPE>
  class HEX_TRIANGULATION_INFO {
  public:

    /// Hex vertex (0,1,...,7) which is the triangulation anchor.
    /// All triangulation tetrahedra are incident on the triangulation_anchor.
    ANCHOR_TYPE triangulation_anchor;

    /// Flags indicating triangulation of hex facet diagonals.
    FLAGS_TYPE triangulation_flags;

    HEX_TRIANGULATION_INFO() 
    {
      triangulation_flags = 0;
    }

    template <typename FTYPE>
    bool TriangulationFlag(const FTYPE jfacet)
    {
      const FLAGS_TYPE mask = (FLAGS_TYPE(1) << jfacet);
      return(!((mask & triangulation_flags) == 0));
    }

    template <typename FTYPE>
    void FlipTriangulationFlag(const FTYPE jfacet)
    {
      const FLAGS_TYPE mask = (FLAGS_TYPE(1) << jfacet);
      triangulation_flags = (triangulation_flags ^ mask);
    }
  };

  ///@}



  // **************************************************
  //! @name TRIANGULATE SQUARE PYRAMIDS
  // **************************************************

  ///@{

  /// Triangulate square pyramid adding diagonal (v0,v2).
  /// - Vertices of pyramid base are listed in clockwise/counter-clockwise
  ///   order around the square base.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_square_pyramid_diagonal02
  (const VTYPE0 apex, const VTYPE1 v0, const VTYPE1 v1,
   const VTYPE1 v2, const VTYPE1 v3,
   std::vector<VTYPE2> & tet_vert_list)
  {
    add_tetrahedron_vertices(apex, v0, v2, v3, tet_vert_list);
    add_tetrahedron_vertices(apex, v2, v0, v1, tet_vert_list);
  }


  /// Triangulate square pyramid adding diagonal (v1,v3).
  /// - Vertices of pyramid base are listed in clockwise/counter-clockwise
  ///   order around the square base.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2>
  void triangulate_square_pyramid_diagonal13
  (const VTYPE0 apex, const VTYPE1 v0, const VTYPE1 v1,
   const VTYPE1 v2, const VTYPE1 v3,
   std::vector<VTYPE2> & tet_vert_list)
  {
    add_tetrahedron_vertices(apex, v1, v3, v0, tet_vert_list);
    add_tetrahedron_vertices(apex, v3, v1, v2, tet_vert_list);
  }


  /// Triangulate square pyramid.  Use facet_triangulation_bits
  ///   to determine triangulation.
  /// - Vertices of pyramid base are listed in clockwise/counter-clockwise
  ///   order around the square base.
  /// - If ((facet_triangulation_bits & mask) == 0), then
  ///   triangulate facet with diagonal (v0,v2).
  /// - Otherwise, triangulate facet with diagonal (v1,v2).
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename BITS_TYPE>
  void triangulate_square_pyramid_diagonal_bits
  (const VTYPE0 apex, const VTYPE1 v0, const VTYPE1 v1,
   const VTYPE1 v2, const VTYPE1 v3,
   const BITS_TYPE facet_triangulation_bits,
   const BITS_TYPE mask,
   std::vector<VTYPE2> & tet_vert_list)
  {
    if ((facet_triangulation_bits & mask) == 0) {
      triangulate_square_pyramid_diagonal02
        (apex, v0, v1, v2, v3, tet_vert_list);
    }
    else {
      triangulate_square_pyramid_diagonal13
        (apex, v0, v1, v2, v3, tet_vert_list);
    }
  }

  ///@}


  // **************************************************
  //! @name TRIANGULATE HEXAHEDRA
  // **************************************************

  ///@{

  /// Triangulate a hexahedron with tetrahedra all incident
  ///   on the diagonal (hex_vert[0], hex_vert[7])
  template <typename VTYPE0, typename VTYPE1>
  void triangulate_hexahedron_diagonal07
  (const VTYPE0 hex_vert[], std::vector<VTYPE1> & tet_vert_list)
  {
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[3], hex_vert[2], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[2], hex_vert[6], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[6], hex_vert[4], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[4], hex_vert[5], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[5], hex_vert[1], tet_vert_list);
    add_tetrahedron_vertices
      (hex_vert[0], hex_vert[7], hex_vert[1], hex_vert[3], tet_vert_list);
  }


  /// Triangulate a hexahedron using a triangulation whose tetrahedra
  ///   all contain vertex v0.  Triangulation of the facets not
  ///   containing v0 is determined by the facet_triangulation_bits.
  /// @param hex_vert[] Array of 8 vertices of the hexahedron.
  /// @param v0 Vertex index between 0 and 7, inclusive.
  ///    All tetrahedra are incident on v0.
  /// @param facet_triangulation_bits Bits indicating the triangulation
  ///   of facets not incident on v0.
  ///   - If bit i is 0, then the triangulation diagonal of facet i
  ///     is either incident on vertex 0 or on vertex 7.
  ///   - If bit i is 1, then the triangulation diagonal of facet i
  ///     is not incident on vertex 0 or on vertex 7.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename BITS_TYPE>
  void triangulate_hexahedron_anchored
  (const VTYPE0 hex_vert[], const VTYPE1 v0,
   const BITS_TYPE facet_triangulation_bits,
   std::vector<VTYPE2> & tet_vert_list)
  {
    const int NUM_CUBE_FACETS(6);
    const BITS_TYPE mask[NUM_CUBE_FACETS] = 
      { BITS_TYPE(1), BITS_TYPE(2), BITS_TYPE(4),
        BITS_TYPE(8), BITS_TYPE(16), BITS_TYPE(32) };

    if (((v0 & mask[0]) == 0)) {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[1], hex_vert[3], hex_vert[7], hex_vert[5],
         facet_triangulation_bits, mask[3], tet_vert_list);
    }
    else {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[0], hex_vert[4], hex_vert[6], hex_vert[2],
         facet_triangulation_bits, mask[0], tet_vert_list);
    }

    if (((v0 & mask[1]) == 0)) {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[2], hex_vert[6], hex_vert[7], hex_vert[3],
         facet_triangulation_bits, mask[4], tet_vert_list);
    }
    else {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[0], hex_vert[1], hex_vert[5], hex_vert[4],
         facet_triangulation_bits, mask[1], tet_vert_list);
    }

    if (((v0 & mask[2]) == 0)) {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[4], hex_vert[5], hex_vert[7], hex_vert[6],
         facet_triangulation_bits, mask[5], tet_vert_list);
    }
    else {
      triangulate_square_pyramid_diagonal_bits
        (hex_vert[v0], hex_vert[0], hex_vert[2], hex_vert[3], hex_vert[1],
         facet_triangulation_bits, mask[2], tet_vert_list);
    }

  }

  /// Triangulate all hexahedron in a hex mesh using triangulation anchors
  ///   and triangulation flags.<br>
  /// Return list of tetrahedra vertices.
  /// @param hex_mesh Hexahedral mesh. 
  ///    Contains member fields poly_data[].triangulation_anchor and
  ///    poly_data[].triangulation_flags.
  template <typename HEX_MESH_TYPE, typename VTYPE>
  void triangulate_anchored_hex_mesh
  (const HEX_MESH_TYPE & hex_mesh,
   std::vector<VTYPE> & tet_vert_list)
  {
    typedef typename HEX_MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;

    for (NUMBER_TYPE ipoly = 0; ipoly < hex_mesh.NumPoly(); ipoly++) {
      triangulate_hexahedron_anchored
        (hex_mesh.VertexList(ipoly), 
         hex_mesh.poly_data[ipoly].triangulation_anchor,
         hex_mesh.poly_data[ipoly].triangulation_flags,
         tet_vert_list);
    }
  }


  /// Set the triangulation anchors for a hexahedral mesh.
  /// - No facet contains two distinct triangulation anchors 
  ///        (associated with different hexahedra.)
  template <typename VERTEX_POLY_INCIDENCE_TYPE,
            typename HEX_MESH_TYPE>
  void set_hex_triangulation_anchors
  (const VERTEX_POLY_INCIDENCE_TYPE & vertex_poly_incidence,
   HEX_MESH_TYPE & hex_mesh)
  {
    typedef typename HEX_MESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename HEX_MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    std::vector<bool> flag_anchor_set(hex_mesh.NumPoly(), false);

    for (NTYPE ipoly = 0; ipoly < hex_mesh.NumPoly(); ipoly++) {
      if (!flag_anchor_set[ipoly]) {
        hex_mesh.poly_data[ipoly].triangulation_anchor = 0;
        flag_anchor_set[ipoly] = true;
        const VTYPE iv0 = hex_mesh.Vertex(ipoly, 0);

        for (NTYPE j = 0; j < vertex_poly_incidence.NumIncidentPoly(iv0); j++) {
          const NTYPE jpoly = vertex_poly_incidence.IncidentPoly(iv0, j);

          if (flag_anchor_set[jpoly]) { continue; }

          NTYPE kloc;
          if (does_list_contain
              (hex_mesh.VertexList(jpoly), hex_mesh.NumPolyVert(jpoly), 
               iv0, kloc)) {
            hex_mesh.poly_data[jpoly].triangulation_anchor = kloc;
            flag_anchor_set[jpoly] = true;
          }
        }
      }
    }

  }


  // Set triangulation flags for hexahedra in hex_mesh
  // @pre Triangulation anchors are already set.
  // @pre No facet contains two distinct triangulation anchors 
  //        (associated with different hexahedra.)
  template <typename HEX_MESH_TYPE>
  void set_hex_triangulation_flags(HEX_MESH_TYPE & hex_mesh)
  {
    typedef typename HEX_MESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename HEX_MESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const NTYPE DIM3(3);
    const IJK::CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube_face_info(DIM3);
    const NTYPE NUM_VERT_PER_TETRAHEDRON(4);
    const NTYPE NUM_CUBE_VERT(8);
    const NTYPE NUM_CUBE_FACETS(6);
    const char DEFAULT_TRIANGULATION_FLAG[NUM_CUBE_VERT]
      = { char(0), char(06), char(05), char(030),
          char(03), char(050), char(060), char(0) };
    const bool VERTEX_FACET_FLAGS[NUM_CUBE_VERT][NUM_CUBE_FACETS]
      = { { true, true, true, false, false, false},
          { false, true, true, true, false, false},
          { true, false, true, false, true, false},
          { false, false, true, true, true, false},
          { true, true, false, false, false, true},
          { false, true, false, true, false, true},
          { true, false, false, false, true, true},
          { false, false, false, true, true, true} };

    for (NTYPE ipoly = 0; ipoly < hex_mesh.NumPoly(); ipoly++) {
      const VTYPE iv = hex_mesh.poly_data[ipoly].triangulation_anchor;
      hex_mesh.poly_data[ipoly].triangulation_flags =
        DEFAULT_TRIANGULATION_FLAG[iv];
    }
  
    IJK::HEXMESH_ADJACENCY<NTYPE,NTYPE,NTYPE> hexmesh_adjacency;
    VTYPE diagonal0_endpoint[2];
    VTYPE diagonal1_endpoint[2];

    hexmesh_adjacency.SetFromMeshOfCubes(hex_mesh, cube_face_info);

    for (NTYPE ipoly0 = 0; ipoly0 < hex_mesh.NumPoly(); ipoly0++) {

      for (NTYPE jfacet0 = 0; jfacet0 < NUM_CUBE_FACETS; jfacet0++) {
        if (hexmesh_adjacency.IsAdjacent(ipoly0, jfacet0)) {
          const NTYPE ipoly1 = hexmesh_adjacency.AdjacentPoly(ipoly0, jfacet0);
          const NTYPE jfacet1 = 
            hexmesh_adjacency.AdjacentPolyFacet(ipoly0, jfacet0);
          const bool tri0_flag = 
            hex_mesh.poly_data[ipoly0].TriangulationFlag(jfacet0);
          const bool tri1_flag = 
            hex_mesh.poly_data[ipoly1].TriangulationFlag(jfacet1);

          get_hexahedron_facet_diagonal
            (hex_mesh.VertexList(ipoly0), jfacet0, tri0_flag, 
             diagonal0_endpoint);
          get_hexahedron_facet_diagonal
            (hex_mesh.VertexList(ipoly1), jfacet1, tri1_flag, 
             diagonal1_endpoint);

          if (diagonal0_endpoint[0] != diagonal1_endpoint[0]) {
            const VTYPE iv1 = hex_mesh.poly_data[ipoly1].triangulation_anchor;
            if (VERTEX_FACET_FLAGS[iv1][jfacet1]) {
              hex_mesh.poly_data[ipoly0].FlipTriangulationFlag(jfacet0);
            }
            else {
              hex_mesh.poly_data[ipoly1].FlipTriangulationFlag(jfacet1);
            }
          }
        }
      }
    }

  }

  /// Triangulate a hexahedral mesh.<br>
  /// Return list of tetrahedra vertices.
  /// - Version which takes a hex mesh with member fields storing
  ///   triangulation anchors and triangulation flags.
  template <typename HEX_MESH_WITH_TRI_DATA, typename VTYPE>
  void triangulate_hex_meshX
  (HEX_MESH_WITH_TRI_DATA & hex_mesh, std::vector<VTYPE> & tet_vert_list)
  {
    typedef typename HEX_MESH_WITH_TRI_DATA::NUMBER_TYPE NTYPE;

    IJK::VERTEX_POLY_INCIDENCE<NTYPE,NTYPE> vertex_poly_incidence;

    vertex_poly_incidence.Set(hex_mesh);

    set_hex_triangulation_anchors(vertex_poly_incidence, hex_mesh);
    set_hex_triangulation_flags(hex_mesh);

    triangulate_anchored_hex_mesh(hex_mesh, tet_vert_list);
  }


  /// Triangulate a hexahedral mesh.<br>
  /// Return list of tetrahedra vertices.
  template <typename HEX_MESH_TYPE, typename VTYPE>
  void triangulate_hex_mesh
  (const HEX_MESH_TYPE & hex_mesh, std::vector<VTYPE> & tet_vert_list)
  {
    typedef typename HEX_MESH_TYPE::NUMBER_TYPE NTYPE;

    IJK::VERTEX_POLY_INCIDENCE<NTYPE,NTYPE> vertex_poly_incidence;
    IJK::POLYMESH_DATA<VTYPE,NTYPE,
      IJK::HEX_TRIANGULATION_INFO<char,char>> hex_data;

    hex_data.Set(hex_mesh);
    triangulate_hex_meshX(hex_data, tet_vert_list);
  }

  ///@}

}

#endif
