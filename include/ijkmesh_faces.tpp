/*!
 *  @file ijkmesh_faces.tpp
 *  @brief ijk templates for processing polyhedral meshes faces.
 * - Version 0.4.0
 */

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2017-2024 Rephael Wenger

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

#ifndef _IJKMESH_FACES_
#define _IJKMESH_FACES_

#include "ijk.tpp"
#include "ijkhash.tpp"
#include "ijklist.tpp"
#include "ijkmesh.tpp"
#include "ijkmesh_datastruct.tpp"

#include <vector>

namespace IJK {


  // **************************************************
  // ORIENTATION FUNCTIONS
  // **************************************************

  /// Return true if orientation of polygon edgeA matches orientation
  ///   of polygon edgeB.
  /// @pre polymesh has at least one polygon.
  /// @pre edgeA and edgeB have same set of endpoints.
  template <typename POLYMESH_TYPE, 
            typename FACET_INFO_TYPEA, typename FACET_INFO_TYPEB>
  bool do_polygon_orientations_match
  (const POLYMESH_TYPE & polymesh, 
   const FACET_INFO_TYPEA & edgeA, const FACET_INFO_TYPEB & edgeB)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
 
    const NTYPE ipolyA = edgeA.poly_containing_face;
    const NTYPE ipolyB = edgeB.poly_containing_face;
    const NTYPE iedgeA = edgeA.face_index;
    const NTYPE iedgeB = edgeB.face_index;

    if (polymesh.Vertex(ipolyA, iedgeA) != polymesh.Vertex(ipolyB, iedgeB))
      { return(true); }
    else
      { return(false); }
  }


  /// Return true if orientation of tetrahedron facetA matches orientation
  ///   of tetrahedron facetB.
  /// @pre facetA and facetB have same set of vertices.
  /// @pre polymesh has at least one simplex.
  template <typename POLYMESH_TYPE, 
            typename FACET_INFO_TYPEA, typename FACET_INFO_TYPEB>
  bool do_tetrahedra_orientations_match
  (const POLYMESH_TYPE & polymesh, 
   const FACET_INFO_TYPEA & facetA, const FACET_INFO_TYPEB & facetB)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename POLYMESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const NTYPE NUM_VERT_PER_TET(4);
    const NTYPE NUM_VERT_PER_FACET(NUM_VERT_PER_TET-1);
    VTYPE facetA_list[NUM_VERT_PER_FACET];
    VTYPE facetB_list[NUM_VERT_PER_FACET];

    const NTYPE ipolyA = facetA.poly_containing_face;
    const NTYPE ipolyB = facetB.poly_containing_face;
    const NTYPE ifacetA = facetA.face_index;
    const NTYPE ifacetB = facetB.face_index;
    for (NTYPE i = 0; i < NUM_VERT_PER_FACET; i++) {
      const NTYPE jA = (ifacetA+i+1)%NUM_VERT_PER_TET;
      const NTYPE jB = (ifacetB+i+1)%NUM_VERT_PER_TET;
      facetA_list[i] = polymesh.Vertex(ipolyA, jA);
      facetB_list[i] = polymesh.Vertex(ipolyB, jB);
    }

    bool flag_even = true;
    if ((ifacetA+ifacetB)%2 == 1) { flag_even = false; }

    // @pre facetA_list[] and facetB_list[] contain the same set of vertices.
    if (facetA_list[0] == facetB_list[0]) {
      if (facetA_list[1] == facetB_list[1]) { return(!flag_even); }
      else { return(flag_even); }
    }
    else if (facetA_list[0] == facetB_list[1]) {
      if (facetA_list[1] == facetB_list[2]) { return(!flag_even); }
      else { return(flag_even); }
    }
    else {
      // facetA_list[0] == facetB_list[2]

      if (facetA_list[1] == facetB_list[0]) { return(!flag_even); }
      else { return(flag_even); }
    }

  }


  /// Return true if parity of simplex facetA is even.
  /// @param facetA Contains simplex and facet of simplex
  ///               Simplex facet j does not contain j'th vertex of simplex.
  template <typename POLYMESH_TYPE, 
            typename FACET_INFO_TYPE>
  bool is_simplex_facet_parity_even
  (const POLYMESH_TYPE & polymesh, const FACET_INFO_TYPE & facetA)
  {
    typedef typename FACET_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename FACET_INFO_TYPE::FACE_INDEX_TYPE FTYPE;
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;

    const PTYPE ipolyA = facetA.poly_containing_face;
    const FTYPE ifacetA = facetA.face_index;
    const NTYPE num_polyA_vert = polymesh.NumPolyVert(ipolyA);
    bool flag_even_parity;

    flag_even_parity = true;
    for (NTYPE j0 = 0; j0 < num_polyA_vert; j0++) {
      if (j0 == ifacetA) { continue; }
      for (NTYPE j1 = j0+1; j1 < num_polyA_vert; j1++) {
        if (j1 == ifacetA) { continue; }
        if (polymesh.Vertex(ipolyA, j0) > polymesh.Vertex(ipolyA, j1))
          { flag_even_parity = !flag_even_parity; }
      }
    }

    if (ifacetA%2 == 1) { flag_even_parity = !flag_even_parity; }

    return(flag_even_parity);
  }


  template <typename POLYMESH_TYPE, 
            typename FACET_INFO_TYPEA, typename FACET_INFO_TYPEB>
  bool do_simplex_orientations_match
  (const POLYMESH_TYPE & polymesh, 
   const FACET_INFO_TYPEA & facetA, const FACET_INFO_TYPEB & facetB)
  {
    const bool is_parityA_even = 
      is_simplex_facet_parity_even(polymesh, facetA);
    const bool is_parityB_even = 
      is_simplex_facet_parity_even(polymesh, facetB);
    return((is_parityA_even != is_parityB_even));
  }


  /// Return true if orientation of hexahedron facetA matches orientation
  ///   of hexahedron facetB.
  template <typename POLYMESH_TYPE, 
            typename FACET_INFO_TYPEA, typename FACET_INFO_TYPEB,
            typename CUBE_TYPE>
  bool do_hexahedra_orientations_match
  (const POLYMESH_TYPE & polymesh, 
   const FACET_INFO_TYPEA & facetA, const FACET_INFO_TYPEB & facetB,
   const CUBE_TYPE & cube)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename POLYMESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const NTYPE NUM_VERT_PER_CUBE_FACET(cube.NumFacetVertices());
    VTYPE facetA_list[NUM_VERT_PER_CUBE_FACET];
    VTYPE facetB_list[NUM_VERT_PER_CUBE_FACET];

    const NTYPE ipolyA = facetA.poly_containing_face;
    const NTYPE ipolyB = facetB.poly_containing_face;
    const NTYPE ifacetA = facetA.face_index;
    const NTYPE ifacetB = facetB.face_index;
    for (NTYPE i = 0; i < NUM_VERT_PER_CUBE_FACET; i++) {
      const NTYPE jA = cube.FacetVertex(ifacetA, i);
      const NTYPE jB = cube.FacetVertex(ifacetB, i);
      facetA_list[i] = polymesh.Vertex(ipolyA, jA);
      facetB_list[i] = polymesh.Vertex(ipolyB, jB);
    }

    std::swap(facetA_list[2], facetA_list[3]);
    std::swap(facetB_list[2], facetB_list[3]);


    // Location of facetB_list[0] in facetA_list.
    NTYPE jloc0;
    if (!does_list_contain
        (facetA_list, NUM_VERT_PER_CUBE_FACET, facetB_list[0], jloc0)) {
      // facetA and facetB are not the same.
      // Cannot compare orientation.
      return(true);
    }

    const NTYPE jloc1 = (jloc0+1)%NUM_VERT_PER_CUBE_FACET;

    if (facetA_list[jloc1] == facetB_list[1]) {
      // Polytopes have opposite orientation.
      return(false); 
    }
    else
      { return(true); }
  }


  // **************************************************
  // POLYTOPE ADJACENCY FUNCTIONS
  // **************************************************

  /// Return true if 2D polygons kpoly and jpoly are adjacent.
  template <typename POLYMESH_TYPE, typename PTYPE0, typename PTYPE1>
  bool are_poly2D_adjacent
  (const POLYMESH_TYPE & polymesh, const PTYPE0 kpoly, const PTYPE1 jpoly)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename POLYMESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    for (NTYPE k0 = 0; k0 < polymesh.NumPolyVert(kpoly); k0++) {
      VTYPE kv0 = polymesh.Vertex(kpoly, k0);
      NTYPE k1 = (k0+1)%polymesh.NumPolyVert(kpoly);
      VTYPE kv1 = polymesh.Vertex(kpoly, k1);

      for (NTYPE j0 = 0; j0 < polymesh.NumPolyVert(jpoly); j0++) {
        VTYPE jv0 = polymesh.Vertex(jpoly, j0);
        NTYPE j1 = (j0+1)%polymesh.NumPolyVert(jpoly);
        VTYPE jv1 = polymesh.Vertex(jpoly, j1);
        if (kv0 == jv0 && kv1 == jv1) { return(true); }
        if (kv0 == jv1 && kv1 == jv0) { return(true); }
      }
    }

    return(false);
  }


  /// Return true if simplices ks and js share facet.
  template <typename POLYMESH_TYPE, typename STYPE0, typename STYPE1>
  bool are_simplices_adjacent
  (const POLYMESH_TYPE & polymesh, const STYPE0 ks, const STYPE1 js)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename POLYMESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const NTYPE numv_per_simplex = polymesh.NumPolyVert(js);
    const NTYPE numv_per_facet = numv_per_simplex-1;
    VTYPE jfacet_vert[numv_per_simplex];
    VTYPE kfacet_vert[numv_per_simplex];

    for (NTYPE i = 0; i < numv_per_simplex; i++) {
      jfacet_vert[i] = polymesh.Vertex(js,i);
      kfacet_vert[i] = polymesh.Vertex(ks,i);
    }
			
    std::sort(jfacet_vert, jfacet_vert+numv_per_simplex);
    std::sort(kfacet_vert, kfacet_vert+numv_per_simplex);
  
    NTYPE j = 0;
    NTYPE k = 0;
    NTYPE num_match = 0;
    while (j < numv_per_simplex && k < numv_per_simplex) {
      if (jfacet_vert[j] == kfacet_vert[k]) {
        j++;
        k++;
        num_match++;
      }
      else if (jfacet_vert[j] > kfacet_vert[k]) {
        k++;
      }
      else {
        j++;
      }
    }

    if (num_match >= numv_per_facet)
      return(true);
    else
      return(false);
  }


  /// Return true if hexahedra ihexA and ihexB share facet.
  template <typename POLYMESH_TYPE, typename HTYPEA, typename HTYPEB,
            typename CUBE_TYPE>
  bool are_hexahedra_adjacent
  (const POLYMESH_TYPE & polymesh, const HTYPEA ihexA, const HTYPEB ihexB,
   const CUBE_TYPE & cube)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename POLYMESH_TYPE::VERTEX_INDEX_TYPE VTYPE;

    const NTYPE DIM3(3);
    const NTYPE NUMV_PER_CUBE = 8;
    const NTYPE NUMV_PER_CUBE_FACET = NUMV_PER_CUBE/2;
    const NTYPE NUM_CUBE_FACETS = 6;
    VTYPE facetA_vert[NUMV_PER_CUBE_FACET];
    VTYPE facetB_vert[NUMV_PER_CUBE_FACET];

    for (NTYPE jfA = 0; jfA < NUM_CUBE_FACETS; jfA++) {
    
      for (NTYPE i = 0; i < NUMV_PER_CUBE_FACET; i++) {
        NTYPE iA = cube.FacetVertex(jfA, i);
        facetA_vert[i] = polymesh.Vertex(ihexA, iA);
      }
      std::sort(facetA_vert, facetA_vert+NUMV_PER_CUBE_FACET);

      for (NTYPE jfB = 0; jfB < NUM_CUBE_FACETS; jfB++) {

        for (NTYPE i = 0; i < NUMV_PER_CUBE_FACET; i++) {
          NTYPE iB = cube.FacetVertex(jfB, i);
          facetB_vert[i] = polymesh.Vertex(ihexB, iB);
        }
        std::sort(facetB_vert, facetB_vert+NUMV_PER_CUBE_FACET);

        if (std::equal(facetA_vert, facetA_vert+NUMV_PER_CUBE_FACET,
                       facetB_vert)) 
          { return(true);  }
      }
    }

    return(false);
  }

  template <typename VTYPE, typename ETYPE>
  bool are_edges_adjacent
  (const std::vector<VTYPE> & edge_endpoint,
   const ETYPE je, const ETYPE ke)
  {
    const ETYPE j = 2*je;
    const ETYPE k = 2*ke;

    if (edge_endpoint[j] == edge_endpoint[k])
      { return(true); }
    if (edge_endpoint[j+1] == edge_endpoint[k])
      { return(true); }
    if (edge_endpoint[j] == edge_endpoint[k+1])
      { return(true); }
    if (edge_endpoint[j+1] == edge_endpoint[k+1])
      { return(true); }

    return(false);
  }


  template <typename VTYPE>
  bool are_edges_connected
  (const std::vector<VTYPE> & edge_endpoint)
  {
    typedef typename std::vector<VTYPE>::size_type SIZE_TYPE;

    const SIZE_TYPE num_edges = edge_endpoint.size()/2;

    if (num_edges == 0) { return(true); }

    std::vector<bool> is_reachable(num_edges,false);
    std::vector<SIZE_TYPE> reachable;
    reachable.push_back(0);

    while (reachable.size() > 0) {
      int je = reachable.back();
      reachable.pop_back();
      is_reachable[je] = true;

      for (int ke = 0; ke < num_edges; ke++) {
        if (!is_reachable[ke]) {
          if (are_edges_adjacent(edge_endpoint, je, ke)) 
            { reachable.push_back(ke); }
        }
      }
    }

    for (SIZE_TYPE ke = 0; ke < num_edges; ke++) {
      if (!is_reachable[ke]) 
        { return(false); }
    }

    return(true);
  }


  // **************************************************
  // SUBROUTINES FOR PROCESSING FACES
  // **************************************************

  namespace MESH_FACES {

    /// Store vertices of facet kf in C++ vector facet_vert.
    /// Store facet kf in C++ vector facet_list.
    template <typename FACET_LIST_TYPE, typename FACET_INDEX_TYPE, 
              typename VTYPE, typename FACE_INFO_TYPE>
    void store_facet
    (const FACET_LIST_TYPE & facets, const FACET_INDEX_TYPE kf,
     std::vector<VTYPE> & facet_vert, 
     std::vector<FACE_INFO_TYPE> & facet_list)
    {
      typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
      typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;

      // Store facet vertices.
      add_list(facets.VertexList(kf), facets.NumFacetVert(kf), facet_vert);

      // Store facet.
      const PTYPE jpoly = facets.PolyContainingFacet(kf);
      const FTYPE jf = facets.FacetIndex(kf);
      facet_list.push_back(FACE_INFO_TYPE(jpoly, jf));
    }


    /// Store vertices of facet kf0 in C++ vector facet_vert.
    /// Store facets kf0 and kf1 in C++ vector facet_list.
    template <typename FACET_LIST_TYPE, typename FACET_INDEX_TYPE, 
              typename VTYPE, typename FACE_INFO_TYPE>
    void store_two_facets
    (const FACET_LIST_TYPE & facets, 
     const FACET_INDEX_TYPE kf0, const FACET_INDEX_TYPE kf1,
     std::vector<VTYPE> & facet_vert, 
     std::vector<FACE_INFO_TYPE> & facet_list)
    {
      typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
      typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;

      // Store vertices of facet kf0.
      add_list(facets.VertexList(kf0), facets.NumFacetVert(kf0), facet_vert);

      // Store facets.
      const PTYPE jpoly0 = facets.PolyContainingFacet(kf0);
      const FTYPE jf0 = facets.FacetIndex(kf0);
      facet_list.push_back(FACE_INFO_TYPE(jpoly0, jf0));
      const PTYPE jpoly1 = facets.PolyContainingFacet(kf1);
      const FTYPE jf1 = facets.FacetIndex(kf1);
      facet_list.push_back(FACE_INFO_TYPE(jpoly1, jf1));
    }


    /// Store vertices of facet sorted_facet[k0] in C++ vector facet_vert.
    /// Store facets sorted_facet[k0], sorted_facet[k0+1], ...
    ///   sorted_facet[k1-1] in C++ vector facet_list.
    /// @pre k0 < k1.
    template <typename FACET_LIST_TYPE, typename FACET_INDEX_TYPE,
              typename NTYPE0, typename NTYPE1,
              typename VTYPE, typename FACE_INFO_TYPE>
    void store_facets_in_sorted_range
    (const FACET_LIST_TYPE & facets,
     const std::vector<FACET_INDEX_TYPE> & sorted_facet,
     const NTYPE0 k0, const NTYPE1 k1,
     std::vector<VTYPE> & facet_vert, 
     std::vector<FACE_INFO_TYPE> & facet_list)
    {
      typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
      typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;

      const NTYPE0 kf0 = sorted_facet[k0];

      // Store vertices of facet kf0.
      add_list(facets.VertexList(kf0), facets.NumFacetVert(kf0), facet_vert);

      // Store facets.
      for (int m = k0; m < k1; m++) {
        const FACET_INDEX_TYPE mf = sorted_facet[m];
        const PTYPE jpoly = facets.PolyContainingFacet(mf);
        const FTYPE jf = facets.FacetIndex(mf);
        facet_list.push_back(FACE_INFO_TYPE(jpoly, jf));
      }
    }


    template <typename POLYMESH_TYPE, typename PTYPE, typename NTYPE,
              typename VTYPE>
    void store_polygon_edges
    (const POLYMESH_TYPE & polymesh, const PTYPE ipoly,
     NTYPE & num_edges, std::vector<VTYPE> & facet_vert_list)
    {
      NTYPE k = 2*num_edges;

      for (NTYPE j0 = 0; j0 < polymesh.NumPolyVert(ipoly); j0++) {
        const NTYPE j1 = (j0+1)%polymesh.NumPolyVert(ipoly);
        facet_vert_list[k] = polymesh.Vertex(ipoly, j0);
        facet_vert_list[k+1] = polymesh.Vertex(ipoly, j1);
        if (facet_vert_list[k] > facet_vert_list[k+1])
          { std::swap(facet_vert_list[k], facet_vert_list[k+1]); }
        num_edges++;
        k += 2;
      }
    }

  }


  // **************************************************
  // GET LINK ROUTINES
  // **************************************************

  /// Get boundary edges of mesh of 2D polygons.
  template <typename POLYMESH_TYPE, typename VTYPE>
  void get_boundary_edges_in_2D_mesh
  (const POLYMESH_TYPE & polymesh, 
   std::vector<VTYPE> & boundary_edge_endpoint)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;

    using namespace IJK::MESH_FACES;

    const NTYPE num_edges = sum_num_poly_vert(polymesh);

    if (num_edges == 0) { return; }

    std::vector<VTYPE> edge_endpoint(2*num_edges);

    NTYPE num_edges2 = 0;
    for (NTYPE ipoly = 0; ipoly < polymesh.NumPoly(); ipoly++) {
      store_polygon_edges(polymesh, ipoly, num_edges2, edge_endpoint);
    }
    std::vector<NTYPE> index_sorted(num_edges);

    for (NTYPE i = 0; i < num_edges; i++) 
      { index_sorted[i] = i; }

    TUPLE_LESS_THAN<VTYPE,VTYPE> edge_less_than(2, &(edge_endpoint.front()));
    sort(index_sorted.begin(), index_sorted.end(), edge_less_than);

    NTYPE j = 0;
    const VTYPE * endpoint_front = &(edge_endpoint.front());
    while (j+1 < num_edges) {
      NTYPE k = j+1;
      NTYPE jf = index_sorted[j];
      NTYPE kf = index_sorted[k];
      while (k < num_edges && 
             are_lists_equal(endpoint_front+2*jf, endpoint_front+2*kf, 2)) {
        k++;
        kf = index_sorted[k];
      };

      NTYPE num_duplicate = k - j;

      if (num_duplicate == 1) {

        // store boundary edge
        for (NTYPE iv = jf*2; iv < (jf+1)*2; iv++) {
          NTYPE iv2 = edge_endpoint[iv];
          boundary_edge_endpoint.push_back(iv2);
        }
      }

      j = k;
    }
  }


  /// Get boundary edges of link of vertex iv in hexahedral mesh.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename VTYPE0, typename VTYPE1, typename CUBE_TYPE>
  void get_link_boundary_edges_in_hex_mesh
  (const POLYMESH_TYPE & polymesh, 
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
   const VTYPE0 iv,
   const CUBE_TYPE & cube,
   std::vector<VTYPE1> & boundary_edge_endpoint)
  {
    POLYMESH_TYPE link_mesh;

    boundary_edge_endpoint.clear();

    compute_vertex_link_in_hex_mesh
      (polymesh, vertex_info, iv, cube, link_mesh);

    // Convert quad vertices to cylic order around quadrilaterals.
    reorder_quad_vertices(link_mesh.element);

    get_boundary_edges_in_2D_mesh(link_mesh, boundary_edge_endpoint);
  }


  /// Get boundary edges of link of vertex iv in tetrahedral mesh.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename VTYPE0, typename VTYPE1>
  void get_link_boundary_edges_in_tet_mesh
  (const POLYMESH_TYPE & polymesh, 
   const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
   const VTYPE0 iv,
   std::vector<VTYPE1> & boundary_edge_endpoint)
  {
    POLYMESH_TYPE link_mesh;

    boundary_edge_endpoint.clear();

    compute_vertex_link_in_tet_mesh
      (polymesh, vertex_info, iv, link_mesh);

    get_boundary_edges_in_2D_mesh(link_mesh, boundary_edge_endpoint);
  }


  // **************************************************
  // GET NON-MANIFOLD AND BOUNDARY FACETS
  // **************************************************

  /// Get non-manifold and boundary edges of 2D mesh.
  template <typename POLYMESH_TYPE,
            typename VTYPE, typename FACE_INFO_TYPE>
  void get_non_manifold_and_boundary_edges_of_mesh2D
    (const POLYMESH_TYPE & polymesh,
     const std::vector<bool> & flag_skip_poly,
     std::vector<VTYPE> & non_manifold_facet_vert,
     std::vector<FACE_INFO_TYPE> & non_manifold_facets,
     std::vector<VTYPE> & boundary_facet_vert,
     std::vector<FACE_INFO_TYPE> & boundary_facets,
     std::vector<VTYPE> & orientation_mismatch_facet_vert,
     std::vector<FACE_INFO_TYPE> & orientation_mismatch_facets)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;
    typedef FACE_INFO_BASE<PTYPE,NTYPE> FINFO_TYPE;

    FACET_LIST_BASE<VTYPE,NTYPE,FINFO_TYPE> facets;
    std::vector<NTYPE> sorted_facet;

    using namespace IJK::MESH_FACES;

    facets.SetFrom2DMesh(polymesh, flag_skip_poly);
    facets.SortPolyVert();
    facets.GetSortedPolytopeIndices(sorted_facet);

    NTYPE k0 = 0;
    while (k0 < facets.NumFacets()) {
      NTYPE kf0 = sorted_facet[k0];
      NTYPE k1 = k0+1;
      while (k1 < facets.NumFacets()) {
        const NTYPE kf1 = sorted_facet[k1];
        if (!facets.ArePolytopesEqual(kf0, kf1))
          { break; }
        k1++;
      }

      NTYPE num_duplicate = k1 - k0;
      if (num_duplicate > 2) {
        store_facets_in_sorted_range
          (facets, sorted_facet, k0, k1, non_manifold_facet_vert,
           non_manifold_facets);
      }
      else if (num_duplicate == 1) {
        store_facet(facets, kf0, boundary_facet_vert, boundary_facets);
      }
      else if (num_duplicate == 2) {

        const NTYPE kf0 = sorted_facet[k0];
        const NTYPE kf1 = sorted_facet[k0+1];
        if (!do_polygon_orientations_match
            (polymesh, facets.poly_data[kf0], facets.poly_data[kf1])) {

          store_two_facets(facets, kf0, kf1, orientation_mismatch_facet_vert,
                           orientation_mismatch_facets);
        }
      }

      k0 = k1;
    }

  }


  /// Get non-manifold and boundary facets of tetrahedral mesh.
  template <typename POLYMESH_TYPE,
            typename VTYPE, typename FACE_INFO_TYPE>
  void get_non_manifold_and_boundary_facets_of_tet_mesh
    (const POLYMESH_TYPE & polymesh,
     std::vector<VTYPE> & non_manifold_facet_vert,
     std::vector<FACE_INFO_TYPE> & non_manifold_facets,
     std::vector<VTYPE> & boundary_facet_vert,
     std::vector<FACE_INFO_TYPE> & boundary_facets,
     std::vector<VTYPE> & orientation_mismatch_facet_vert,
     std::vector<FACE_INFO_TYPE> & orientation_mismatch_facets)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;
    typedef FACE_INFO_BASE<PTYPE,NTYPE> FINFO_TYPE;

    FACET_LIST_BASE<VTYPE,NTYPE,FINFO_TYPE> facets;
    std::vector<NTYPE> sorted_facet;

    using namespace IJK::MESH_FACES;

    facets.SetFromMeshOfSimplices(polymesh);
    facets.SortPolyVert();
    facets.GetSortedPolytopeIndices(sorted_facet);

    NTYPE k0 = 0;
    while (k0 < facets.NumFacets()) {
      NTYPE kf0 = sorted_facet[k0];
      NTYPE k1 = k0+1;
      while (k1 < facets.NumFacets()) {
        const NTYPE kf1 = sorted_facet[k1];
        if (!facets.ArePolytopesEqual(kf0, kf1))
          { break; }
        k1++;
      }

      NTYPE num_duplicate = k1 - k0;
      if (num_duplicate > 2) {
        store_facets_in_sorted_range
          (facets, sorted_facet, k0, k1, non_manifold_facet_vert,
           non_manifold_facets);
      }
      else if (num_duplicate == 1) {
        store_facet(facets, kf0, boundary_facet_vert, boundary_facets);
      }
      else if (num_duplicate == 2) {

        const NTYPE kf0 = sorted_facet[k0];
        const NTYPE kf1 = sorted_facet[k0+1];
        if (!do_tetrahedra_orientations_match
            (polymesh, facets.poly_data[kf0], facets.poly_data[kf1])) {

          store_two_facets(facets, kf0, kf1, orientation_mismatch_facet_vert,
                           orientation_mismatch_facets);
        }
      }

      k0 = k1;
    }

  }


  /// Get non-manifold and boundary facets of simplicial mesh.
  template <typename POLYMESH_TYPE,
            typename VTYPE, typename FACE_INFO_TYPE>
  void get_non_manifold_and_boundary_facets_of_simplicial_mesh
    (const POLYMESH_TYPE & polymesh,
     std::vector<VTYPE> & non_manifold_facet_vert,
     std::vector<FACE_INFO_TYPE> & non_manifold_facets,
     std::vector<VTYPE> & boundary_facet_vert,
     std::vector<FACE_INFO_TYPE> & boundary_facets,
     std::vector<VTYPE> & orientation_mismatch_facet_vert,
     std::vector<FACE_INFO_TYPE> & orientation_mismatch_facets)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;
    typedef FACE_INFO_BASE<PTYPE,NTYPE> FINFO_TYPE;

    FACET_LIST_BASE<VTYPE,NTYPE,FINFO_TYPE> facets;
    std::vector<NTYPE> sorted_facet;

    using namespace IJK::MESH_FACES;

    facets.SetFromMeshOfSimplices(polymesh);
    facets.SortPolyVert();
    facets.GetSortedPolytopeIndices(sorted_facet);

    NTYPE k0 = 0;
    while (k0 < facets.NumFacets()) {
      NTYPE kf0 = sorted_facet[k0];
      NTYPE k1 = k0+1;
      while (k1 < facets.NumFacets()) {
        const NTYPE kf1 = sorted_facet[k1];
        if (!facets.ArePolytopesEqual(kf0, kf1))
          { break; }
        k1++;
      }

      NTYPE num_duplicate = k1 - k0;
      if (num_duplicate > 2) {
        store_facets_in_sorted_range
          (facets, sorted_facet, k0, k1, non_manifold_facet_vert,
           non_manifold_facets);
      }
      else if (num_duplicate == 1) {
        store_facet(facets, kf0, boundary_facet_vert, boundary_facets);
      }
      else if (num_duplicate == 2) {

        const NTYPE kf0 = sorted_facet[k0];
        const NTYPE kf1 = sorted_facet[k0+1];
        if (!do_simplex_orientations_match
            (polymesh, facets.poly_data[kf0], facets.poly_data[kf1])) {

          store_two_facets(facets, kf0, kf1, orientation_mismatch_facet_vert,
                           orientation_mismatch_facets);
        }
      }

      k0 = k1;
    }

  }


  /// Get non-manifold and boundary facets of hex mesh.
  template <typename POLYMESH_TYPE, typename CUBE_TYPE,
            typename VTYPE, typename FACE_INFO_TYPE>
  void get_non_manifold_and_boundary_facets_of_hex_mesh
    (const POLYMESH_TYPE & polymesh, const CUBE_TYPE & cube,
     std::vector<VTYPE> & non_manifold_facet_vert,
     std::vector<FACE_INFO_TYPE> & non_manifold_facets,
     std::vector<VTYPE> & boundary_facet_vert,
     std::vector<FACE_INFO_TYPE> & boundary_facets,
     std::vector<VTYPE> & orientation_mismatch_facet_vert,
     std::vector<FACE_INFO_TYPE> & orientation_mismatch_facets)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename FACE_INFO_TYPE::POLY_INDEX_TYPE PTYPE;
    typedef typename FACE_INFO_TYPE::FACE_INDEX_TYPE FTYPE;
    typedef FACE_INFO_BASE<PTYPE,NTYPE> FINFO_TYPE;

    FACET_LIST_BASE<VTYPE,NTYPE,FINFO_TYPE> facets;
    std::vector<NTYPE> sorted_facet;

    using namespace IJK::MESH_FACES;

    facets.SetFromMeshOfCubes(polymesh, cube);
    facets.SortPolyVert();
    facets.GetSortedPolytopeIndices(sorted_facet);

    NTYPE k0 = 0;
    while (k0 < facets.NumFacets()) {
      NTYPE kf0 = sorted_facet[k0];
      NTYPE k1 = k0+1;
      while (k1 < facets.NumFacets()) {
        const NTYPE kf1 = sorted_facet[k1];
        if (!facets.ArePolytopesEqual(kf0, kf1))
          { break; }
        k1++;
      }

      NTYPE num_duplicate = k1 - k0;
      if (num_duplicate > 2) {
        store_facets_in_sorted_range
          (facets, sorted_facet, k0, k1, non_manifold_facet_vert,
           non_manifold_facets);
      }
      else if (num_duplicate == 1) {
        store_facet(facets, kf0, boundary_facet_vert, boundary_facets);
      }
      else if (num_duplicate == 2) {

        const NTYPE kf0 = sorted_facet[k0];
        const NTYPE kf1 = sorted_facet[k0+1];
        if (!do_hexahedra_orientations_match
            (polymesh, facets.poly_data[kf0], facets.poly_data[kf1], cube)) {

          store_two_facets(facets, kf0, kf1, orientation_mismatch_facet_vert,
                           orientation_mismatch_facets);
        }
      }

      k0 = k1;
    }

  }




  // **************************************************
  // GET NON-MANIFOLD VERTICES
  // **************************************************

  /// Get non-manifold vertices of 2D mesh.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename VTYPE>
  void get_non_manifold_vertices_of_mesh2D
    (const POLYMESH_TYPE & polymesh,
     const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
     const std::vector<bool> & flag_skip_vertex,
     std::vector<VTYPE> & non_manifold_vert_list)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename VERTEX_POLY_INCIDENCE_TYPE::POLY_INDEX_TYPE PTYPE;

    const NTYPE num_vertices = vertex_info.NumVertices();
    non_manifold_vert_list.clear();

    for (VTYPE iv = 0; iv < vertex_info.NumVertices(); iv++) {

      if (flag_skip_vertex[iv]) { continue; }

      const NTYPE num_incident = vertex_info.NumIncidentPoly(iv);
      if (num_incident > 0) {

        std::vector<bool> is_reachable(num_incident, false);
        std::vector<int> reachable;
        reachable.push_back(0);

        while (reachable.size() > 0) {
          NTYPE n = reachable.size()-1;
          NTYPE j = reachable[n];
          reachable.pop_back();
          is_reachable[j] = true;
          PTYPE jpoly = vertex_info.IncidentPoly(iv, j);

          for (NTYPE k = 0; k < num_incident; k++) {
            if (!is_reachable[k]) {
              PTYPE kpoly = vertex_info.IncidentPoly(iv, k);
              if (are_poly2D_adjacent(polymesh, kpoly, jpoly)) {
                reachable.push_back(k);
              }
            }
          }
        }

        for (int k = 0; k < num_incident; k++) {
          if (!is_reachable[k]) {
            non_manifold_vert_list.push_back(iv);
            break;
          }

        }
      }
    }

  }


  /// Get non-manifold vertices of simplicial mesh.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename VTYPE>
  void get_non_manifold_vertices_of_simplicial_mesh
    (const POLYMESH_TYPE & polymesh,
     const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
     const std::vector<bool> & flag_skip_vertex,
     std::vector<VTYPE> & non_manifold_vert_list)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename VERTEX_POLY_INCIDENCE_TYPE::POLY_INDEX_TYPE PTYPE;

    const NTYPE num_vertices = vertex_info.NumVertices();
    non_manifold_vert_list.clear();

    for (VTYPE iv = 0; iv < vertex_info.NumVertices(); iv++) {

      if (flag_skip_vertex[iv]) { continue; }

      const NTYPE num_incident = vertex_info.NumIncidentPoly(iv);
      if (num_incident > 0) {

        std::vector<bool> is_reachable(num_incident, false);
        std::vector<int> reachable;
        reachable.push_back(0);

        while (reachable.size() > 0) {
          NTYPE n = reachable.size()-1;
          NTYPE j = reachable[n];
          reachable.pop_back();
          is_reachable[j] = true;
          PTYPE jpoly = vertex_info.IncidentPoly(iv, j);

          for (NTYPE k = 0; k < num_incident; k++) {
            if (!is_reachable[k]) {
              PTYPE kpoly = vertex_info.IncidentPoly(iv, k);
              if (are_simplices_adjacent(polymesh, kpoly, jpoly)) {
                reachable.push_back(k);
              }
            }
          }
        }

        for (int k = 0; k < num_incident; k++) {
          if (!is_reachable[k]) {
            non_manifold_vert_list.push_back(iv);
            break;
          }

        }
      }
    }

  }


  /// Get non-manifold vertices of tetrahedral mesh.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename VTYPE>
  void get_non_manifold_vertices_of_tet_mesh
    (const POLYMESH_TYPE & polymesh,
     const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
     const std::vector<bool> & flag_skip_vertex,
     std::vector<bool> & non_manifold_vert,
     std::vector<VTYPE> & non_manifold_vert_list)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename VERTEX_POLY_INCIDENCE_TYPE::POLY_INDEX_TYPE PTYPE;

    const NTYPE num_vertices = vertex_info.NumVertices();
    non_manifold_vert_list.clear();

    for (VTYPE iv = 0; iv < vertex_info.NumVertices(); iv++) {

      if (flag_skip_vertex[iv]) { continue; }

      const NTYPE num_incident = vertex_info.NumIncidentPoly(iv);
      if (num_incident > 0) {

        std::vector<bool> is_reachable(num_incident, false);
        std::vector<int> reachable;
        reachable.push_back(0);

        while (reachable.size() > 0) {
          NTYPE n = reachable.size()-1;
          NTYPE j = reachable[n];
          reachable.pop_back();
          is_reachable[j] = true;
          PTYPE jpoly = vertex_info.IncidentPoly(iv, j);

          for (NTYPE k = 0; k < num_incident; k++) {
            if (!is_reachable[k]) {
              PTYPE kpoly = vertex_info.IncidentPoly(iv, k);
              if (are_simplices_adjacent(polymesh, kpoly, jpoly)) {
                reachable.push_back(k);
              }
            }
          }
        }

        for (int k = 0; k < num_incident; k++) {
          if (!is_reachable[k]) {
            non_manifold_vert[iv] = true;
            non_manifold_vert_list.push_back(iv);
            break;
          }

        }
      }
    }

    // Check that the boundary of each link is connected.
    for (VTYPE iv = 0; iv < vertex_info.NumVertices(); iv++) {

      if (non_manifold_vert[iv] || flag_skip_vertex[iv])
        { continue; }

      std::vector<VTYPE> boundary_edge_endpoint;

      get_link_boundary_edges_in_tet_mesh
        (polymesh, vertex_info, iv, boundary_edge_endpoint);

      if (!are_edges_connected(boundary_edge_endpoint)) {
        non_manifold_vert_list.push_back(iv);
        non_manifold_vert[iv] = true;
      }
    }

  }


  /// Get non-manifold vertices of hexahedral mesh.
  template <typename POLYMESH_TYPE, typename VERTEX_POLY_INCIDENCE_TYPE,
            typename CUBE_TYPE, typename VTYPE>
  void get_non_manifold_vertices_of_hex_mesh
    (const POLYMESH_TYPE & polymesh,
     const VERTEX_POLY_INCIDENCE_TYPE & vertex_info,
     const std::vector<bool> & flag_skip_vertex,
     CUBE_TYPE & cube,
     std::vector<bool> & non_manifold_vert,
     std::vector<VTYPE> & non_manifold_vert_list)
  {
    typedef typename POLYMESH_TYPE::NUMBER_TYPE NTYPE;
    typedef typename VERTEX_POLY_INCIDENCE_TYPE::POLY_INDEX_TYPE PTYPE;

    const NTYPE num_vertices = vertex_info.NumVertices();
    non_manifold_vert_list.clear();

    for (VTYPE iv = 0; iv < vertex_info.NumVertices(); iv++) {

      if (flag_skip_vertex[iv]) { continue; }

      const NTYPE num_incident = vertex_info.NumIncidentPoly(iv);
      if (num_incident > 0) {

        std::vector<bool> is_reachable(num_incident, false);
        std::vector<int> reachable;
        reachable.push_back(0);

        while (reachable.size() > 0) {
          NTYPE n = reachable.size()-1;
          NTYPE j = reachable[n];
          reachable.pop_back();
          is_reachable[j] = true;
          PTYPE jpoly = vertex_info.IncidentPoly(iv, j);

          for (NTYPE k = 0; k < num_incident; k++) {
            if (!is_reachable[k]) {
              PTYPE kpoly = vertex_info.IncidentPoly(iv, k);
              if (are_hexahedra_adjacent(polymesh, kpoly, jpoly, cube)) {
                reachable.push_back(k);
              }
            }
          }
        }

        for (int k = 0; k < num_incident; k++) {
          if (!is_reachable[k]) {
            non_manifold_vert[iv] = true;
            non_manifold_vert_list.push_back(iv);
            break;
          }

        }
      }
    }

    // Check that the boundary of each link is connected.
    for (VTYPE iv = 0; iv < vertex_info.NumVertices(); iv++) {

      if (non_manifold_vert[iv] || flag_skip_vertex[iv])
        { continue; }

      std::vector<VTYPE> boundary_edge_endpoint;

      get_link_boundary_edges_in_hex_mesh
        (polymesh, vertex_info, iv, cube, boundary_edge_endpoint);

      if (!are_edges_connected(boundary_edge_endpoint)) {
        non_manifold_vert_list.push_back(iv);
        non_manifold_vert[iv] = true;
      }
    }

  }

  // **************************************************
  // EDGE HASH TABLE
  // **************************************************

  template <typename VTYPE, typename EDGE_INDEX_TYPE>
  class UNORIENTED_EDGE_HASH_TABLE:
    public unordered_map_key_pair<VTYPE,VTYPE,EDGE_INDEX_TYPE> {

  protected:

    typedef unordered_map_key_pair<VTYPE,VTYPE,EDGE_INDEX_TYPE> 
    BASE_CLASS;

  public:
    typedef typename BASE_CLASS::iterator iterator;
    typedef typename BASE_CLASS::const_iterator const_iterator;

    // Undefine insert, find, findData;
    std::pair<const_iterator,bool> 
    insert (const VTYPE iv0, const VTYPE iv1, const EDGE_INDEX_TYPE ie);
    iterator find(const VTYPE iv0, const VTYPE iv1);
    const_iterator find(const VTYPE iv0, const VTYPE iv1) const;


  public:
    UNORIENTED_EDGE_HASH_TABLE(){};

    std::pair<const_iterator,bool>
    Insert(const VTYPE iv0, const VTYPE iv1, const EDGE_INDEX_TYPE ie)
    {
      if (iv0 <= iv1) { return(BASE_CLASS::insert(iv0, iv1, ie)); }
      else { return(BASE_CLASS::insert(iv1, iv0, ie)); }
    }

    iterator Find(const VTYPE iv0, const VTYPE iv1)
    {
      if (iv0 <= iv1) { return(BASE_CLASS::find(iv0, iv1)); }
      else { return(BASE_CLASS::find(iv1, iv0)); };
    }

    const_iterator Find(const VTYPE iv0, const VTYPE iv1) const
    {
      if (iv0 <= iv1) { return(BASE_CLASS::find(iv0, iv1)); }
      else { return(BASE_CLASS::find(iv1, iv0)); };
    }

    /// Return true if hash table contains edge (iv0,iv1).
    /// @param[out] it Iterator pointing to hash table element
    ///   or to end of hash table.
    bool Contains(const VTYPE iv0, const VTYPE iv1, 
                  const_iterator & it) const
    {
      it = Find(iv0,iv1);
      return(!(it == this->end())); 
    }

    /// Return true if hash table contains edge (iv0,iv1).
    /// - Version which does not return iterator.
    bool Contains(const VTYPE iv0, const VTYPE iv1) const
    {
      const_iterator it;
      return(Contains(iv0,iv1,it));
    }

    /// Find edge index of edge (iv0,iv1).
    /// Return false if edge (iv0,iv1) is not in the hash table.
    bool FindEdgeIndex(const VTYPE iv0, const VTYPE iv1,
                       EDGE_INDEX_TYPE & ie) const
    {
      const_iterator it = Find(iv0,iv1);
      if (it == this->end()) {
        ie = 0;
        return(false);
      }
      else {
        ie = it->second;
        return(true);
      }
    }

  };

}

#endif

