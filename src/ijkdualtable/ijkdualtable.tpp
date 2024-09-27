/*!
 *  @file ijkdualtable.tpp
 *  @brief Class containing dual lookup table of isosurface vertices.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2012-2024 Rephael Wenger

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

#ifndef _IJKDUALTABLE_TPP_
#define _IJKDUALTABLE_TPP_

#include <climits>
#include <limits>
#include <vector>

#include "ijk.tpp"
#include "ijkbits.tpp"
#include "ijkcube.tpp"

/// Classes and routines for storing and manipulating dual isosurface lookup table.
namespace IJKDUALTABLE {

  // **************************************************
  //! @name COMPUTE FUNCTIONS
  // **************************************************

  /// Compute complement index.
  template <typename ITYPE, typename NTYPE>
  inline ITYPE compute_complement
  (const ITYPE ival, const NTYPE num_table_entries)
  { return(num_table_entries-1-ival); }
    

  // **************************************************
  //! @name ISODUAL CUBE TABLE ROUTINES
  // **************************************************

  /// @brief Create isodual cube table entry ientry.
  /// @tparam TI_TYPE Table index type.
  template <typename TI_TYPE, typename FIND_TYPE,
            typename CUBE_TYPE, typename ENTRY_TYPE>
  void create_isodual_cube_table_entry
  (const TI_TYPE ientry,
   const bool flag_separate_neg, const bool flag_separate_opposite,
   const CUBE_TYPE & cube,
   FIND_TYPE & find_component, 
   ENTRY_TYPE & table_entry)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NUMBER_TYPE;
    
    const NUMBER_TYPE num_cube_vertices = cube.NumVertices();
    const bool flag_separate_pos = (!flag_separate_neg);

    find_component.ClearAll();
    find_component.SetVertexFlags(ientry);

    if (flag_separate_opposite) {
      if (flag_separate_neg) {
        if (IJK::is_two_opposite_ones(ientry, num_cube_vertices))
          { find_component.NegateVertexFlags(); };
      }
      else {
        if (IJK::is_two_opposite_zeros(ientry, num_cube_vertices))
          { find_component.NegateVertexFlags(); };
      }
    }


    NUMBER_TYPE num_components(0);
    for (NUMBER_TYPE i = 0; i < num_cube_vertices; i++) {
      if (find_component.Component(i) == 0) {

        if (find_component.VertexFlag(i) == flag_separate_pos) {
          num_components++;
          find_component.Search(i, num_components);
        }
      }
    }

    NUMBER_TYPE num_zeros, num_ones;
    IJK::count_bits(ientry, num_cube_vertices, num_zeros, num_ones);
    if (num_zeros == 0 || num_ones == 0) {
      // All vertices are negative or all vertices are positive.
      table_entry.num_vertices = 0;
    }
    else {
      table_entry.num_vertices = num_components;
    }

    for (NUMBER_TYPE ie = 0; ie < cube.NumEdges(); ie++) {
      const NUMBER_TYPE iv0 = cube.EdgeEndpoint(ie, 0);
      const NUMBER_TYPE iv1 = cube.EdgeEndpoint(ie, 1);

      if (find_component.VertexFlag(iv0) == find_component.VertexFlag(iv1)) {
        table_entry.is_bipolar[ie] = false;
        table_entry.incident_isovertex[ie] = 0;
      }
      else {
        table_entry.is_bipolar[ie] = true;
        int icomp = find_component.Component(iv0);
        if (find_component.VertexFlag(iv1) == flag_separate_pos) {
          // Vertex iv1 is negative.
          icomp = find_component.Component(iv1);
        }
        table_entry.incident_isovertex[ie] = icomp-1;
      }
    }

    table_entry.SetIncidentIsoPoly(cube.NumEdges());
  }


  // **************************************************
  //! @name AMBIGUITY ROUTINE
  // **************************************************

  /*!
   *  @brief  Return true if cube with given configuration is ambiguous.
   *  @param ientry Index of +/- configuration.
   *  @param[out] ambiguous_facet_bits Integer representing set 
   *    of ambiguous facets.
   *    - k'th bit of facet_set is 1 if facet k is ambiguous.
   *    - Facet k is orthogonal to axis (k%dimension).
   *    - If k < dimension, facet k contains cube vertex 0.
   *    - If k >= dimension, facet k contains cube vertex with largest index.
   *  @param[out] num_ambiguous_facets Number of ambiguous facets.
   */
  template <typename TI_TYPE, typename CUBE_TYPE, typename FBITS_TYPE,
            typename NTYPE, typename FIND_TYPE>
  bool is_cube_ambiguous
  (const TI_TYPE ientry, const CUBE_TYPE & cube,
   FBITS_TYPE & ambiguous_facet_bits,  NTYPE & num_ambiguous_facets, 
   FIND_TYPE & find_component)
  {
    const NTYPE num_cube_facets = cube.NumFacets();

    compute_ambiguous_cube_facets
      (ientry, num_cube_facets, ambiguous_facet_bits,
       num_ambiguous_facets, find_component);

    if (num_ambiguous_facets > 0) { return(true); }

    // Cube may still be ambiguous, even though it has no ambiguous facets.

    // Check if cube is ambiguous.
    NTYPE num_pos_components = 
      find_component.ComputeNumComponents(ientry, true);
    NTYPE num_neg_components = 
      find_component.ComputeNumComponents(ientry, false);

    if (num_pos_components > 1 || num_neg_components > 1)
      { return(true); }
    else
      { return(false); }
  }


  /*!
   *  @brief Return true if cube facet kf is ambiguous.
   *    - Facet kf is orthogonal to axis (kf%dimension).
   *    - If kf < dimension, facet kf contains cube vertex 0.
   *    - If kf >= dimension, facet kf contains cube vertex with largest index.
   */
  template <typename TI_TYPE, typename FTYPE, typename FIND_TYPE>
  bool is_cube_facet_ambiguous
  (const TI_TYPE ientry, const FTYPE kf, 
   FIND_TYPE & find_component)
  {
    typedef typename FIND_TYPE::NUMBER_TYPE NTYPE;

    NTYPE num_pos_components = 
      find_component.ComputeNumComponentsInFacet(ientry, kf, true);
    NTYPE num_neg_components = 
      find_component.ComputeNumComponentsInFacet(ientry, kf, false);

    if (num_pos_components > 1 || num_neg_components > 1)
      { return(true); }
    else
      { return(false); }
  }


  /*!
   *  @brief Compute ambiguous cube facets.
   *  @param[out] facet_bits Integer representing set of ambiguous facets.
   *    - k'th bit of facet_set is 1 if facet k is ambiguous.
   *    - Facet k is orthogonal to axis (k%dimension).
   *    - If k < dimension, facet k contains cube vertex 0.
   *    - If k >= dimension, facet k contains cube vertex with largest index.
   *  @param[out] num_ambiguous_facets Number of ambiguous facets.
   */
  template <typename TI_TYPE, typename NTYPE, typename NTYPE2,
            typename FBITS_TYPE, typename FIND_TYPE>
  void compute_ambiguous_cube_facets
  (const TI_TYPE ientry, 
   const NTYPE num_facets,
   FBITS_TYPE & facet_bits,  
   NTYPE2 & num_ambiguous_facets,
   FIND_TYPE & find_component)
  {
    facet_bits = 0;
    num_ambiguous_facets = 0;
    FBITS_TYPE mask = FBITS_TYPE(1);

    for (NTYPE kf = 0; kf < num_facets; kf++) {
      if (is_cube_facet_ambiguous(ientry, kf, find_component)) {
        num_ambiguous_facets++;
        facet_bits = (facet_bits | mask);
      }
      mask = (mask << FBITS_TYPE(1));
    }
  }


  // **************************************************
  //! @name ISODUAL CUBE TABLE AMBIG ROUTINES
  // **************************************************

  /// @brief Create isodual cube table ambig entry ientry.
  /// @tparam TI_TYPE Table index type.
  template <typename TI_TYPE, typename FIND_TYPE,
            typename CUBE_TYPE, typename ENTRY_TYPE>
  void create_isodual_cube_table_ambig_entry
  (const TI_TYPE ientry,
   const bool flag_separate_neg, const bool flag_separate_opposite,
   const CUBE_TYPE & cube,
   FIND_TYPE & find_component, 
   ENTRY_TYPE & table_entry)
  {
    create_isodual_cube_table_entry
      (ientry, flag_separate_neg, flag_separate_opposite,
       cube, find_component, table_entry);

    table_entry.is_ambiguous = 
      is_cube_ambiguous(ientry, cube, table_entry.ambiguous_facet_bits,
                        table_entry.num_ambiguous_facets, find_component);

    table_entry.num_active_facets =
      cube.ComputeNumActiveCubeFacets(ientry);
  }


  // **************************************************
  //! @name IVOLDUAL CUBE TABLE ROUTINES
  // **************************************************

  /// @brief Set position of polytope vertices relative to interval volume.
  template <typename TI_TYPE, typename NTYPE, typename ENTRY_TYPE>
  void set_position_relative_to_interval_volume4
  (const TI_TYPE ientry, const NTYPE num_poly_vertices, 
   ENTRY_TYPE & table_entry)
  {
    const int NUM_VERTEX_TYPES = 4;

    TI_TYPE x = ientry;
    for (NTYPE iv = 0; iv < num_poly_vertices; iv++) {
      TI_TYPE y = x%NUM_VERTEX_TYPES;

      if (y == 0) {
        table_entry.poly_vertex_info[iv].SetBelowIntervalVolume();
      }
      else if (y == 3) {
        table_entry.poly_vertex_info[iv].SetAboveIntervalVolume();
      }
      else {
        table_entry.poly_vertex_info[iv].SetInIntervalVolume();
      }
      x = x/NUM_VERTEX_TYPES;
    }
  }


  /// Set vertex type.
  template <typename TI_TYPE, typename NTYPE, typename ENTRY_TYPE>
  void set_vertex_type4
  (const TI_TYPE ientry, const NTYPE num_poly_vertices, 
   ENTRY_TYPE & table_entry)
  {
    const int NUM_VERTEX_TYPES = 4;

    TI_TYPE x = ientry;
    for (NTYPE iv = 0; iv < num_poly_vertices; iv++) {
      TI_TYPE y = x%NUM_VERTEX_TYPES;

      table_entry.poly_vertex_info[iv].vertex_type = y;

      x = x/NUM_VERTEX_TYPES;
    }
  }


  /// @brief Compute lifted ivol indices.
  /// @param num_cube_vertices Number of vertices of cube (not lifted cube).
  template <typename TI_TYPE, typename NTYPE, 
            typename TI_TYPE_L, typename TI_TYPE_U>
  void compute_lifted_ivol_indices4
  (const TI_TYPE ientry, const NTYPE num_cube_vertices, 
   TI_TYPE_L & ilower, TI_TYPE_U & iupper)
  {
    const int NUM_CUBE_VERTEX_TYPES = 4;
    TI_TYPE factorA, factorB;

    TI_TYPE x = ientry;
    factorA = 1;
    factorB = (TI_TYPE(1) << num_cube_vertices);
    ilower = 0;
    iupper = 0;

    while (x != 0) {
      TI_TYPE vertex_type = x%NUM_CUBE_VERTEX_TYPES;
      if (vertex_type == 1) {
        ilower += factorA;
      }
      else if (vertex_type == 2) {
        ilower += factorA;
        ilower += factorB;
        iupper += factorA;
      }
      else if (vertex_type == 3) {
        iupper += factorA;
        iupper += factorB;
        ilower += factorA;
        ilower += factorB;
      }

      factorA = (factorA << 1);
      factorB = (factorB << 1);

      x = x/NUM_CUBE_VERTEX_TYPES;
    }
  }


  /// @brief Compute lower and upper isosurface table indices from interval volume table index.
  /// @param num_cube_vertices Number of vertices of cube.
  template <typename TI_TYPE, typename NTYPE, 
            typename TI_TYPE_L, typename TI_TYPE_U>
  void compute_isov_indices_from_ivol_index
  (const TI_TYPE ivol_table_index, const NTYPE num_cube_vertices, 
   TI_TYPE_L & lower_isosurface_table_index, 
   TI_TYPE_U & upper_isosurface_table_index)
  {
    const int NUM_CUBE_VERTEX_TYPES = 4;
    TI_TYPE factor;

    TI_TYPE x = ivol_table_index;
    factor = 1;
    lower_isosurface_table_index = 0;
    upper_isosurface_table_index = 0;

    while (x != 0) {
      TI_TYPE vertex_type = x%NUM_CUBE_VERTEX_TYPES;
      if (vertex_type == 3) {
        lower_isosurface_table_index += factor;
        upper_isosurface_table_index += factor;
      }
      else if (vertex_type > 0) {
        lower_isosurface_table_index += factor;
      }

      factor = (factor << 1);

      x = x/NUM_CUBE_VERTEX_TYPES;
    }
  }


  /// Determine if interval volume vertices are on lower or upper isosurfaces.
  template <typename TI_TYPE, typename CUBE_TYPE, 
            typename ISODUAL_ENTRY_TYPE, typename ENTRY_TYPE>
  void determine_isosurface_containing_ivol_vertices4
  (const TI_TYPE ientry, 
   const CUBE_TYPE & cube, const CUBE_TYPE & lifted_cube,
   ISODUAL_ENTRY_TYPE & lower_isodual, ISODUAL_ENTRY_TYPE & upper_isodual,
   ENTRY_TYPE & table_entry)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {

      const DTYPE edge_dir = cube.EdgeDir(ie);
      NTYPE iend0 = cube.EdgeEndpoint(ie,0);
      NTYPE iend1 = cube.EdgeEndpoint(ie,1);

      if (!table_entry.IsInIntervalVolume(iend0))
        { std::swap(iend0, iend1); }

      if (!table_entry.IsInIntervalVolume(iend0)) {
        // Neither edge endpoint is contained in the interval volume.
        continue;
      }

      if (table_entry.IsInIntervalVolume(iend1)) {
        // Both edge endpoints are contained in the interval volume.
        continue;
      }

      const NTYPE ivolv0 = 
        table_entry.poly_vertex_info[iend0].incident_ivol_vertex;

      const NTYPE je_lower = lifted_cube.IncidentEdge(iend0, edge_dir);
      const NTYPE kend0 = iend0 + cube.NumVertices();
      const NTYPE je_upper = lifted_cube.IncidentEdge(kend0, edge_dir);

      if (table_entry.ivolv_info[ivolv0].flag_in_lower_lifted_cube) {
        if (lower_isodual.IsBipolar(je_lower)) 
          { table_entry.ivolv_info[ivolv0].SetToLowerIsosurface(); }
      }
      else {
        if (upper_isodual.IsBipolar(je_upper))
          { table_entry.ivolv_info[ivolv0].SetToUpperIsosurface(); }
      }
    }
  }


  /*!
   *  Modify the ivoldual table entries of the lower and upper lifted cubes
   *    if at least one of the cubes is ambiguous, no facet shared
   *    with a third cube is ambiguous and each cube contains only
   *    one isosurface vertex.
   *  Replace the table entries by their complements.
   */
  template <typename NTYPE, typename CUBE_TYPE,
            typename TI_TYPE_L, typename TI_TYPE_U, 
            typename ISODUAL_ENTRY_TYPE>
  void modify_ambiguous_lifted_ivoldual_entries
  (const NTYPE num_table_entries,
   const CUBE_TYPE & lifted_cube,
   const TI_TYPE_L ilower, const ISODUAL_ENTRY_TYPE & lower_isodual,
   const TI_TYPE_U iupper, const ISODUAL_ENTRY_TYPE & upper_isodual,
   TI_TYPE_L & new_ilower,
   TI_TYPE_U & new_iupper)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE dimension = lifted_cube.Dimension();
    const NTYPE num_facets = lifted_cube.NumFacets();
    NTYPE ilower_facet, iupper_facet;

    new_ilower = ilower;
    new_iupper = iupper;

    if (!lower_isodual.IsAmbiguous() &&
        !upper_isodual.IsAmbiguous()) { return; }

    if (lower_isodual.NumVertices() > 1 ||
        upper_isodual.NumVertices() > 1) { return; }

    if (dimension < 1) { return; }

    ilower_facet = dimension-1;
    iupper_facet = ilower_facet+dimension;

    for (NTYPE jf = 0; jf < num_facets; jf++) {
      if (jf == ilower_facet || jf == iupper_facet) { 
        // Skip the facet between the upper and lower lifted cubes
        //   and the facet on the lifted boundary.
        continue; 
      }

      if (lower_isodual.IsFacetAmbiguous(jf) ||
          upper_isodual.IsFacetAmbiguous(jf)) {
        // Lower or upper cube has some ambiguous facet adjacent
        //   to a third cube.
        // Don't change table entry.
        return;
      }
    }

    // Flip ivoldual table entries.
    new_ilower = compute_complement(ilower, num_table_entries);
    new_iupper = compute_complement(iupper, num_table_entries);
  }


  /// @brief Compute ivoldual cube ambiguity information
  /// @tparam TI_TYPE Table index type.
  template <typename CUBE_TYPE, typename ISODUAL_ENTRY_TYPE,
            typename ENTRY_TYPE>
  void compute_ivoldual_cube_ambiguity_information
  (const CUBE_TYPE & cube,
   const ISODUAL_ENTRY_TYPE & lower_isodual,
   const ISODUAL_ENTRY_TYPE & upper_isodual,
   ENTRY_TYPE & table_entry)
  {
    typedef typename CUBE_TYPE::NUMBER_TYPE NTYPE;
    typedef typename ENTRY_TYPE::FACET_BITS_TYPE FBITS_TYPE;

    table_entry.is_ambiguous = false;

    FBITS_TYPE facet_bits = 0;
    FBITS_TYPE lower_facet_bits = 0;
    FBITS_TYPE upper_facet_bits = 0;
    FBITS_TYPE mask = FBITS_TYPE(1);
    for (NTYPE jf = 0; jf < cube.NumFacets(); jf++) {

      NTYPE jf_lifted = jf;

      if (jf >= cube.Dimension()) { 
        // Skip facet orthogonal to lifting direction.
        jf_lifted++;  
      }

      if (lower_isodual.IsFacetAmbiguous(jf_lifted) ||
          upper_isodual.IsFacetAmbiguous(jf_lifted)) {
        facet_bits = (facet_bits | mask);
        table_entry.is_ambiguous = true;
        table_entry.num_ambiguous_facets++;

        if (lower_isodual.IsFacetAmbiguous(jf_lifted)) 
          { lower_facet_bits = (lower_facet_bits | mask); }

        if (upper_isodual.IsFacetAmbiguous(jf_lifted)) 
          { upper_facet_bits = (upper_facet_bits | mask); }
      }

      mask = (mask << FBITS_TYPE(1));
    }

    table_entry.ambiguous_facet_bits = facet_bits;
    table_entry.ambiguous_facet_bits_in_lower_lifted = lower_facet_bits;
    table_entry.ambiguous_facet_bits_in_upper_lifted = upper_facet_bits;
  }


  /// @brief Compute ivoldual cube non manifold information
  /// @tparam TI_TYPE Table index type.
  template <typename ISODUAL_ENTRY_TYPE, typename ENTRY_TYPE>
  void compute_ivoldual_cube_non_manifold_information
  (const ISODUAL_ENTRY_TYPE & lower_isodual,
   const ISODUAL_ENTRY_TYPE & upper_isodual,
   ENTRY_TYPE & table_entry)
  {
    table_entry.is_non_manifold = false;

    if (table_entry.NumVertices() <= 2) {
      if (lower_isodual.NumAmbiguousFacets() > 0 &&
          lower_isodual.NumVertices() < 2) 
        { table_entry.is_non_manifold = true; }


      if (upper_isodual.NumAmbiguousFacets() > 0 &&
          upper_isodual.NumVertices() < 2) 
        { table_entry.is_non_manifold = true; }

    }
    else if (table_entry.NumAmbiguousFacets() == 0) {
      if (lower_isodual.NumAmbiguousFacets() == 2 &&
          upper_isodual.NumAmbiguousFacets() == 2) {
        if (lower_isodual.NumVertices() < 2 ||
            upper_isodual.NumVertices() < 2)
          { table_entry.is_non_manifold = true; }
      }
    }


  }


  /// @brief Create ivoldual cube table entry ientry.
  /// @tparam TI_TYPE Table index type.
  template <typename TI_TYPE, typename NTYPE, typename FIND_TYPE,
            typename CUBE_TYPE, typename ISODUAL_ENTRY_TYPE,
            typename ENTRY_TYPE>
  void create_ivoldual_cube_table_entry
  (const TI_TYPE ientry,
   const NTYPE num_table_entries,
   const bool flag_separate_neg,
   const CUBE_TYPE & cube,
   const CUBE_TYPE & lifted_cube,
   FIND_TYPE & find_component,
   ISODUAL_ENTRY_TYPE & lower_isodual,
   ISODUAL_ENTRY_TYPE & upper_isodual,
   ENTRY_TYPE & table_entry)
  {
    typedef typename CUBE_TYPE::DIMENSION_TYPE DTYPE;

    const DTYPE lifted_dimension = lifted_cube.Dimension();
    const NTYPE num_cube_vertices = cube.NumVertices();
    const NTYPE num_lifted_cube_vertices = lifted_cube.NumVertices();
    TI_TYPE ilower, iupper;
    TI_TYPE lower_isosurface_table_index, upper_isosurface_table_index;

    set_position_relative_to_interval_volume4
      (ientry, num_cube_vertices, table_entry);
    set_vertex_type4(ientry, num_cube_vertices, table_entry);

    const bool flag_separate_opposite = true;
    compute_lifted_ivol_indices4
      (ientry, num_cube_vertices, ilower, iupper);

    create_isodual_cube_table_ambig_entry
      (ilower, flag_separate_neg, flag_separate_opposite,
       lifted_cube, find_component, lower_isodual);
    create_isodual_cube_table_ambig_entry
      (iupper, flag_separate_neg, flag_separate_opposite,
       lifted_cube, find_component, upper_isodual);

    TI_TYPE new_ilower, new_iupper;
    modify_ambiguous_lifted_ivoldual_entries
      (num_table_entries, lifted_cube, ilower, lower_isodual,
       iupper, upper_isodual, new_ilower, new_iupper);

    if (new_ilower != ilower || new_iupper != iupper) {
      ilower = new_ilower;
      iupper = new_iupper;
      create_isodual_cube_table_ambig_entry
        (ilower, flag_separate_neg, flag_separate_opposite,
         lifted_cube, find_component, lower_isodual);
      create_isodual_cube_table_ambig_entry
        (iupper, flag_separate_neg, flag_separate_opposite,
         lifted_cube, find_component, upper_isodual);
    }

    table_entry.CreateIntervalVolumeVertices
      (lower_isodual.NumVertices(), upper_isodual.NumVertices());

    for (NTYPE ie = 0; ie < lifted_cube.NumEdges(); ie++) {
      const DTYPE edge_dir = lifted_cube.EdgeDir(ie);
      NTYPE iend0 = lifted_cube.EdgeEndpoint(ie,0);
      NTYPE iend1 = lifted_cube.EdgeEndpoint(ie,1);
      if (iend0 > iend1) { std::swap(iend0, iend1); }

      if (edge_dir+1 == lifted_dimension) {
        if (lower_isodual.IsBipolar(ie)) {
          const NTYPE isov = lower_isodual.IncidentIsoVertex(ie);
          table_entry.poly_vertex_info[iend0].SetIncident(isov);
        }
        else if (upper_isodual.IsBipolar(ie)) {
          const NTYPE isov = upper_isodual.IncidentIsoVertex(ie) 
            + lower_isodual.NumVertices();
          table_entry.poly_vertex_info[iend0].SetIncident(isov);
        }
      }
      else {

        if (iend1 < cube.NumVertices()) {
          if (upper_isodual.IsBipolar(ie)) {
            const NTYPE je = cube.IncidentEdge(iend0, edge_dir);
            const NTYPE isov = upper_isodual.IncidentIsoVertex(ie)
              + lower_isodual.NumVertices();

            // Vertices on the upper isodual are on the upper isosurface.
            table_entry.poly_edge_info[je].SetUpperIncident(isov);

            if (table_entry.IsAboveIntervalVolume(iend0) ||
                table_entry.IsAboveIntervalVolume(iend1))
              { table_entry.ivolv_info[isov].SetToUpperIsosurface(); }
          }
        }
        else {
          if (lower_isodual.IsBipolar(ie)) {
            const NTYPE jend0 = iend0 - cube.NumVertices();
            const NTYPE jend1 = iend1 - cube.NumVertices();
            const NTYPE je = cube.IncidentEdge(jend0, edge_dir);
            const NTYPE isov = lower_isodual.IncidentIsoVertex(ie);

            // Vertices on the lower isodual are on the lower isosurface.
            table_entry.poly_edge_info[je].SetLowerIncident(isov);

            if (table_entry.IsBelowIntervalVolume(jend0) ||
                table_entry.IsBelowIntervalVolume(jend1))
              { table_entry.ivolv_info[isov].SetToLowerIsosurface(); }

          }
        }
      }
    }

    compute_ivoldual_cube_ambiguity_information
      (cube, lower_isodual, upper_isodual, table_entry);

    compute_ivoldual_cube_non_manifold_information
      (lower_isodual, upper_isodual, table_entry);

    determine_isosurface_containing_ivol_vertices4
      (ientry, cube, lifted_cube, lower_isodual, upper_isodual, table_entry);

    compute_isov_indices_from_ivol_index
      (ientry, num_cube_vertices, 
       lower_isosurface_table_index, upper_isosurface_table_index);
    table_entry.lower_isosurface_table_index = lower_isosurface_table_index;
    table_entry.upper_isosurface_table_index = upper_isosurface_table_index;
  }


  // **************************************************
  //! @name IVOLDUAL QUERY ROUTINES
  // **************************************************

  /// Return true if ivol vertex is incident on an isosurface polytope dual to edge ie.
  template <typename DUAL_TABLE_TYPE, typename TI_TYPE, typename VTYPE,
            typename ETYPE>
  bool is_ivol_vertex_on_iso_poly_dual_to_edge
  (const DUAL_TABLE_TYPE & table, const TI_TYPE it,
   const VTYPE ivolv, const ETYPE ie)
  {
    typename DUAL_TABLE_TYPE::CUBE_TYPE::NUMBER_TYPE iend[2];

    iend[0] = table.Cube().EdgeEndpoint(ie, 0);
    iend[1] = table.Cube().EdgeEndpoint(ie, 1);

    if (table.EdgeHasDualIVolPoly(it, ie)) {

      if (ivolv == table.LowerIncident(it, ie)) {
        if (table.IsBelowIntervalVolume(it, iend[0]) ||
            table.IsBelowIntervalVolume(it, iend[1])) {

          if (table.OnLowerIsosurface(it, ivolv)) 
            { return(true); }
        }
      }
      else if (ivolv == table.UpperIncident(it, ie)) {
        if (table.IsAboveIntervalVolume(it, iend[0]) ||
            table.IsAboveIntervalVolume(it, iend[1])) {

          if (table.OnUpperIsosurface(it, ivolv)) 
            { return(true); }
        }
      }

    }
    else if (table.IsInIntervalVolume(it, iend[0]) ||
          table.IsInIntervalVolume(it, iend[1])) {

      if (table.IsInIntervalVolume(it, iend[0]))
        { std::swap(iend[0], iend[1]); }

      if (table.IsInIntervalVolume(it, iend[0])) 
        { return(false); }

      if (ivolv == table.IncidentIVolVertex(it, iend[1])) {

        if (table.IsBelowIntervalVolume(it, iend[0])) {
          if (table.OnLowerIsosurface(it, ivolv))
            { return(true); }
        }
        else if (table.IsAboveIntervalVolume(it, iend[0])) {
          if (table.OnUpperIsosurface(it, ivolv))
            { return(true); }
        }
      }
    }

    return(false);
  }


  /// Return true if lower isosurface intersects edge ie.
  template <typename DUAL_TABLE_TYPE, typename TI_TYPE, typename ETYPE>
  bool does_lower_isosurface_intersect_edge
  (const DUAL_TABLE_TYPE & table, const TI_TYPE it, const ETYPE ie)
  {
    typedef typename DUAL_TABLE_TYPE::CUBE_TYPE::NUMBER_TYPE NTYPE;

    NTYPE iend[2];

    iend[0] = table.Cube().EdgeEndpoint(ie, 0);
    iend[1] = table.Cube().EdgeEndpoint(ie, 1);

    if (table.IsBelowIntervalVolume(it, iend[0])) {
      if (table.IsBelowIntervalVolume(it, iend[1])) { return(false); }
      return(true);
    }
    else if (table.IsBelowIntervalVolume(it, iend[1])) 
      { return(true); }

    return(false);
  }


  /*!
   *  @brief Return true if lower isosurface intersects edge ie.
   *  - If true, return interval volume vertex on lower isosurface 
   *    quadrilateral which is dual to edge ie.
   */
  template <typename DUAL_TABLE_TYPE, typename TI_TYPE, typename ETYPE,
            typename VTYPE>
  bool does_lower_isosurface_intersect_edge
  (const DUAL_TABLE_TYPE & table, const TI_TYPE it,
   const ETYPE ie, VTYPE & ivolv)
  {
    typedef typename DUAL_TABLE_TYPE::CUBE_TYPE::NUMBER_TYPE NTYPE;

    NTYPE iend[2];

    iend[0] = table.Cube().EdgeEndpoint(ie, 0);
    iend[1] = table.Cube().EdgeEndpoint(ie, 1);

    if (table.IsBelowIntervalVolume(it, iend[0])) {
      if (table.IsBelowIntervalVolume(it, iend[1])) 
        { return(false); }
    }
    else if (table.IsBelowIntervalVolume(it, iend[1])) 
      { std::swap(iend[0], iend[1]); }
    else 
      { return(false); }

    // iend[0] is below interval volume.
    // iend[1] is in or above interval volume.

    if (table.EdgeHasDualIVolPoly(it, ie)) 
      { ivolv = table.LowerIncident(it, ie); }
      else {
      // Vertex iend[1] must be dual to an interval volume polytope.
      ivolv = table.IncidentIVolVertex(it, iend[1]); 
    }

    return(true);
  }


  /// Return true if upper isosurface intersects edge ie.
  template <typename DUAL_TABLE_TYPE, typename TI_TYPE, typename ETYPE>
  bool does_upper_isosurface_intersect_edge
  (const DUAL_TABLE_TYPE & table, const TI_TYPE it, const ETYPE ie)
  {
    typedef typename DUAL_TABLE_TYPE::CUBE_TYPE::NUMBER_TYPE NTYPE;

    NTYPE iend[2];

    iend[0] = table.Cube().EdgeEndpoint(ie, 0);
    iend[1] = table.Cube().EdgeEndpoint(ie, 1);

    if (table.IsAboveIntervalVolume(it, iend[0])) {
      if (table.IsAboveIntervalVolume(it, iend[1])) { return(false); }
      return(true);
    }
    else if (table.IsAboveIntervalVolume(it, iend[1])) 
      { return(true); }

    return(false);
  }


  /*!
   *  @brief Return true if upper isosurface intersects edge ie.
   *  - If true, return interval volume vertex on upper isosurface 
   *    quadrilateral which is dual to edge ie.
   */
  template <typename DUAL_TABLE_TYPE, typename TI_TYPE, typename ETYPE,
            typename VTYPE>
  bool does_upper_isosurface_intersect_edge
  (const DUAL_TABLE_TYPE & table, const TI_TYPE it,
   const ETYPE ie, VTYPE & ivolv)
  {
    typedef typename DUAL_TABLE_TYPE::CUBE_TYPE::NUMBER_TYPE NTYPE;

    NTYPE iend[2];

    iend[0] = table.Cube().EdgeEndpoint(ie, 0);
    iend[1] = table.Cube().EdgeEndpoint(ie, 1);

    if (table.IsAboveIntervalVolume(it, iend[0])) {
      if (table.IsAboveIntervalVolume(it, iend[1])) 
        { return(false); }
    }
    else if (table.IsAboveIntervalVolume(it, iend[1])) 
      { std::swap(iend[0], iend[1]); }
    else 
      { return(false); }

    // iend[0] is above interval volume.
    // iend[1] is in or below interval volume.

    if (table.EdgeHasDualIVolPoly(it, ie)) 
      { ivolv = table.UpperIncident(it, ie); }
    else {
      // Vertex iend[1] must be dual to an interval volume polytope.
      ivolv = table.IncidentIVolVertex(it, iend[1]); 
    }

    return(true);
  }


  /// @brief Return true if ivol vertex is doubly connected across facet jfacet.
  /// @pre table dimension is 3.
  template <typename DUAL_TABLE_TYPE, typename TI_TYPE, typename VTYPE,
            typename FTYPE>
  bool is_ivol3D_vertex_doubly_connected_across_facet
  (const DUAL_TABLE_TYPE & table, const TI_TYPE it,
   const VTYPE ivolv, const FTYPE jfacet)
  {
    typedef typename DUAL_TABLE_TYPE::CUBE_TYPE::NUMBER_TYPE NTYPE;
    
    const NTYPE DIM3(3);
    const FTYPE facet_orth_dir = table.Cube().FacetOrthDir(jfacet);
    const NTYPE NUM_FACET_EDGES(4);
    NTYPE iv0 = 0;
    NTYPE iedge[NUM_FACET_EDGES];

    if (jfacet >= DIM3)
      { iv0 = table.Cube().VertexNeighbor(0, facet_orth_dir); }

    const FTYPE dir1 = (facet_orth_dir+1)%DIM3;
    const FTYPE dir2 = (facet_orth_dir+2)%DIM3;
    const NTYPE iv1 = table.Cube().VertexNeighbor(iv0, dir1);
    const NTYPE iv2 = table.Cube().VertexNeighbor(iv0, dir2);
    iedge[0] = table.Cube().IncidentEdge(iv0, dir1);
    iedge[1] = table.Cube().IncidentEdge(iv2, dir1);
    iedge[2] = table.Cube().IncidentEdge(iv0, dir2);
    iedge[3] = table.Cube().IncidentEdge(iv1, dir2);

    for (NTYPE i = 0; i < NUM_FACET_EDGES; i++) {

      if (!is_ivol_vertex_on_iso_poly_dual_to_edge
          (table, it, ivolv, iedge[i])) 
        { return(false); }
    }

    return(true);
  }


  /*!
   *  Return true if ivol vertex is doubly connected across some facet.
   *  @pre table dimension is 3.
   *  @param[out] dc_facet Ivol vertex is doubly connected across dc_facet.
   */
  template <typename DUAL_TABLE_TYPE, typename TI_TYPE, typename VTYPE,
            typename FTYPE>
  bool is_ivol3D_vertex_doubly_connected
  (const DUAL_TABLE_TYPE & table, const TI_TYPE it, const VTYPE ivolv,
   FTYPE & dc_facet)
  {
    typedef typename DUAL_TABLE_TYPE::CUBE_TYPE::NUMBER_TYPE NTYPE;

    const NTYPE NUM_CUBE_FACETS = table.Cube().NumFacets();

    // Initialize
    dc_facet = 0;

    if (table.NumAmbiguousFacets(it) == 0) { return(false); }

    if (!table.OnLowerIsosurface(it, ivolv) &&
        !table.OnUpperIsosurface(it, ivolv)) 
      { return(false); }

    for (dc_facet = 0; dc_facet < NUM_CUBE_FACETS; dc_facet++) {
      if (is_ivol3D_vertex_doubly_connected_across_facet
          (table, it, ivolv, dc_facet))
        { return(true); }
    }

    dc_facet = 0;
    return(false);
  }


  /*!
   *  @brief Return true if ivol vertex is doubly connected across some facet.
   *  - Version which does not return dc_facet.
   *  @pre table dimension is 3.
   */
  template <typename DUAL_TABLE_TYPE, typename TI_TYPE, typename VTYPE>
  bool is_ivol3D_vertex_doubly_connected
  (const DUAL_TABLE_TYPE & table, const TI_TYPE it, const VTYPE ivolv)
  {
    typedef typename DUAL_TABLE_TYPE::CUBE_TYPE::NUMBER_TYPE NTYPE;
    NTYPE dc_facet;

    return(is_ivol3D_vertex_doubly_connected(table, it, ivolv, dc_facet));
  }


  // **************************************************
  //! @name UTILITY FUNCTIONS
  // **************************************************

  /// Calculate number of entries required in ISOSURFACE_TABLE.
  template <typename NUM_ENTRIES_TYPE,
            typename NTYPEV, typename NTYPEC>
  NUM_ENTRIES_TYPE calculate_num_entries
  (const NTYPEV num_vert, const NTYPEC num_colors)
  {
    const char * procname = "calculate_num_entries";

    NUM_ENTRIES_TYPE num_table_entries = 0;

    if (num_colors < 1)
      throw IJK::PROCEDURE_ERROR
        (procname, "Number of colors must be positive.");

    const NUM_ENTRIES_TYPE max2 = 
      std::numeric_limits<NUM_ENTRIES_TYPE>::max()/num_colors;

    num_table_entries = 1;
    for (NTYPEV iv = 0; iv < num_vert; iv++) {
      if (num_table_entries > max2)
        throw IJK::PROCEDURE_ERROR
          (procname, "Number of entries is too large.");

      num_table_entries = num_table_entries * num_colors;
    };

    return(num_table_entries);
  }


  /*!
   *  @brief Convert integer to boolean flags.
   *  @param[out] flag[] Array of boolean flags.
   *  @pre flag[] is pre-allocated to size at least num_flags.
   */
  template <typename TI_TYPE, typename NTYPE>
  void convert2bool
  (const TI_TYPE ival, bool flag[], const NTYPE num_flags)
  {
    const TI_TYPE base = 2;

    TI_TYPE jval = ival;
    for (NTYPE i = 0; i < num_flags; i++) {
      flag[i] = bool(jval % base);
      jval = jval/base;
    };

    if (jval != 0) {
      IJK::PROCEDURE_ERROR error("convert2bool");
      error.AddMessage("Error converting ", ival, " to base 2.");
      error.AddMessage("Output has more than ", num_flags, " digits.");

      throw error;
    };
  }


  // **************************************************
  // ISODUAL TABLE ENTRY
  // **************************************************

  /// Entry in the dual isosurface lookup table.
  template <typename NTYPE, typename ISOV_TYPE,
            typename ISOP_TYPE>
  class ISODUAL_TABLE_ENTRY {

  public:
    typedef ISOV_TYPE ISO_VERTEX_INDEX_TYPE;
    typedef ISOP_TYPE ISO_POLY_INDEX_TYPE;

    NTYPE num_vertices;         ///< Number of dualiso vertices in cube.

    /// @brief incident_isovertex[kf] = Isosurface vertex incident on face kf.
    /// - Isosurface face kf is dual to polytope edge kf.
    ISOV_TYPE * incident_isovertex;

    /// @brief  is_bipolar[ke] = True if polytope edge ke is bipolar.
    /// - Cube edge ke is dual to isosurface face ke.
    bool * is_bipolar;

    /// @brief List of isosurface polytopes incident on isosurface vertices.
    /// - Isosurface face kf is dual to polytope edge kf.
    ISOP_TYPE * incident_isopoly;

    /// @brief first_incident_isopoly[j] points to first polytope
    ///   incident on isosurface vertex j.
    ISOP_TYPE ** first_incident_isopoly;

    /// @brief num_incident_isopoly[j] number of polytopes
    ///   incident on isosurface vertex j.
    NTYPE * num_incident_isopoly;


  public:
    ISODUAL_TABLE_ENTRY();      ///< constructor
    ~ISODUAL_TABLE_ENTRY();     ///< destructor
    
    /// Return number of dual isosurface vertices in cube.
    NTYPE NumVertices() const
    { return(num_vertices); }

    /// Return true if edge ke is bipolar.
    template <typename NTYPE2>
    bool IsBipolar(const NTYPE2 ke) const
    { return(is_bipolar[ke]); }

    /// Return incident isosurface vertex.
    template <typename NTYPE2>
    inline ISOV_TYPE IncidentIsoVertex(const NTYPE2 ke) const
    { return incident_isovertex[ke]; }

    /// @brief Return pointer to index of first isosurface polytope
    ///   incident on isosurface vertex k.
    template <typename NTYPE2>
    ISOP_TYPE * FirstIncidentIsoPoly(const NTYPE2 k) const
    { return first_incident_isopoly[k]; }

    /// @brief Return number of isosurface polytopes
    ///   incident on isosurface vertex k.
    template <typename NTYPE2>
    NTYPE NumIncidentIsoPoly(const NTYPE2 k) const
    { return num_incident_isopoly[k]; }

    /*!
     *  @brief Allocate incident_isovert[] and is_bipolar[].
     *  - num_poly_vertices is not used in this Allocate()
     *    but other table entries may use it.
     */
    template <typename NUMV_TYPE, typename NUME_TYPE>
    void Allocate(const NUMV_TYPE num_poly_vertices,
                  const NUME_TYPE num_poly_edges);

    /// @brief Set incident_isopoly and first_incident_isopoly.
    void SetIncidentIsoPoly(const NTYPE num_mesh_poly_edges);
    
    /// @brief Free arrays incident_isopoly[], first_incident_isopoly[]
    ///    and num_incident_isopoly[].
    void FreeIncidentIsoPoly();

    /// @brief Free all memory.
    void FreeAll();

    /// @brief Check data structure for inconsistencies.
    bool Check(const NTYPE num_mesh_poly_edges, IJK::ERROR & error_msg) const;
  };


  // **************************************************
  // ISODUAL AMBIG TABLE ENTRY
  // **************************************************

  /// Entry in the dual isosurface lookup table with ambig polytope/facet information.
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FBITS_TYPE>
  class ISODUAL_AMBIG_TABLE_ENTRY:
    public ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE> {

  public:

    typedef FBITS_TYPE FACET_BITS_TYPE;
    typedef ISOP_TYPE ISO_POLY_INDEX_TYPE;

    /// True if polytope is ambiguous.
    bool is_ambiguous;
    
    /// Number of ambiguous facets.
    NTYPE num_ambiguous_facets;

    /// Number of active facets (with both positive and negative vertices.)
    NTYPE num_active_facets;

    /*!
     *  @brief Integer representing set of ambiguous facets.
     *  - k'th bit of ambiguous_facet is 1 if facet k is ambiguous.
     *  - Facet k is orthogonal to axis (k%dimension).
     *  - If k < dimension, facet k contains cube vertex 0.
     *  - If k >= dimension, facet k contains cube vertex with largest index.
     */
    FBITS_TYPE ambiguous_facet_bits;    

    ISODUAL_AMBIG_TABLE_ENTRY();      ///< constructor
    ~ISODUAL_AMBIG_TABLE_ENTRY();     ///< destructor

    // get functions

    /// Return true if cube is ambiguous.
    bool IsAmbiguous() const
    { return(is_ambiguous); };

    /*!
     *  @brief Return true if facet jf is ambiguous.
     *  - Facet jf is orthogonal to axis (jf%dimension).
     *  - If jf < dimension, facet jf contains cube vertex 0.
     *  - If jf >= dimension, facet jf contains cube vertex with largest index.
     */
    template <typename FTYPE>
    bool IsFacetAmbiguous(const FTYPE jf) const
    { return((ambiguous_facet_bits & (FBITS_TYPE(1) << jf)) != 0); };

    /// Return ambiguous facet bits.
    FBITS_TYPE AmbiguousFacetBits() const
    { return(ambiguous_facet_bits); };

    /// Return number of ambiguous facets.
    NTYPE NumAmbiguousFacets() const
    { return(num_ambiguous_facets); };

    /// Return number of active facets.
    NTYPE NumActiveFacets() const
    { return(num_active_facets); };

    /// Clear ambiguity information.
    void ClearAmbig();
  };


  // **************************************************
  // IVOLDUAL TABLE ENTRY
  // **************************************************

  /// Information about the interval volume vertices in an isosurface cube dual to a grid edge.
  template <typename ISOV_TYPE>
  class IVOLDUAL_POLY_EDGE_INFO {

  protected:
    void Init();

  public:
    /// Lower isosurface vertex incident on interval volume polytope dual to edge.
    ISOV_TYPE lower_incident_isovertex;

    /// Upper isosurface vertex incident on interval volume polytope dual to edge.
    ISOV_TYPE upper_incident_isovertex;

    /*!
     *  @brief If flag_dual_ivolpoly is true, then some isosurface polytope is dual to the edge.
     *  - If true, then both lower_incident_isovertex
     *    and upper_incident_isovertex are set.
     */
    bool flag_dual_ivolpoly;

    IVOLDUAL_POLY_EDGE_INFO() 
    { Init(); }

    /// Set lower_incident_isovertex.
    template <typename ISOV_TYPE2>
    void SetLowerIncident(const ISOV_TYPE2 ilower);

    /// Set upper_incident_isovertex.
    template <typename ISOV_TYPE2>
    void SetUpperIncident(const ISOV_TYPE2 iupper);
  };


  /// Information about the interval volume vertices in an isosurface cube dual to a grid vertex.
  template <typename ISOV_TYPE, typename VERTEX_TYPE>
  class IVOLDUAL_POLY_VERTEX_INFO {

  protected:
    void Init();

  public:

    typedef enum 
      { BELOW_INTERVAL_VOLUME, IN_INTERVAL_VOLUME, ABOVE_INTERVAL_VOLUME}
      RELATIVE_POSITION;

    /// Position of polytope vertex relative to interval volume.
    RELATIVE_POSITION relative_position;

    /// Type of vertex.
    VERTEX_TYPE vertex_type;

    /*!
     *  @brief Interval volume vertex incident on isosurface polytope dual to edge.
     *  - incident_ivol_vertex is defined if and only if relative_position
     *      is IN_INTERVAL_VOLUME.
     */
    ISOV_TYPE incident_ivol_vertex;

    /// True, if grid vertex is below interval volume.
    bool IsBelowIntervalVolume() const
    { return(relative_position == BELOW_INTERVAL_VOLUME); }

    /// True, if grid vertex is in interval volume.
    bool IsInIntervalVolume() const
    { return(relative_position == IN_INTERVAL_VOLUME); }

    /// True, if grid vertex is above interval volume.
    bool IsAboveIntervalVolume() const
    { return(relative_position == ABOVE_INTERVAL_VOLUME); }

    /// @brief Return position of vertex with respect to interval volume.
    /// - Could be below, in, or above interval volume.
    RELATIVE_POSITION RelativePosition() const
    { return(relative_position); }

    /*!
     *  Return vertex_type.
     *  - When there are 4 vertex types, vertices with type 0 are below
     *    the interval volume, vertices with type 3 are above the interval
     *    volume, and vertices with type 1 or 2 are in the interval volume.
     */
    VERTEX_TYPE VertexType() const
    { return(vertex_type); }

    /// Set incident_ivol_vertex.
    template <typename IVOLV_TYPE2>
    void SetIncident(const IVOLV_TYPE2 ivolv);

    /// Set relative position to below interval volume.
    void SetBelowIntervalVolume()
    { relative_position = BELOW_INTERVAL_VOLUME; }

    /// Set relative position to in interval volume.
    void SetInIntervalVolume()
    { relative_position = IN_INTERVAL_VOLUME; }

    /// Set relative position to above interval volume.
    void SetAboveIntervalVolume()
    { relative_position = ABOVE_INTERVAL_VOLUME; }

    IVOLDUAL_POLY_VERTEX_INFO() 
    { Init(); }
  };


  /// Information about interval volume vertices.
  class IVOLDUAL_IVOLV_INFO {

  public:

    /// True if ivol vertex is on lower isosurface.
    bool flag_lower_isosurface;

    /*!
     *  @brief True if ivol vertex is on upper isosurface.
     *  - An ivol vertex cannot be on both the lower and upper isosurface
     *    although an internal ivol vertex could be on neither.
     */
    bool flag_upper_isosurface;

    /*!
     *  @brief True if ivol vertex was generated from lower lifted cube.
     *  - False if ivol vertex was generated from upper lifted cube.
     */
    bool flag_in_lower_lifted_cube;

    IVOLDUAL_IVOLV_INFO()
    {
      flag_lower_isosurface = false;
      flag_upper_isosurface = false;
      flag_in_lower_lifted_cube = true;
    }

    /// Set vertex to be on lower isosurface.
    void SetToLowerIsosurface()
    { flag_lower_isosurface = true; }

    /// Set vertex to be on upper isosurface.
    void SetToUpperIsosurface()
    { flag_upper_isosurface = true; }

  };

  /// Entry in the dual interval volume lookup table.
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FBITS_TYPE, typename TI_TYPE>
  class IVOLDUAL_TABLE_ENTRY {

  protected:

    /// Number of interval volume vertices in lower lifted cube.
    NTYPE num_vertices_in_lower_lifted;

    /// Number of interval volume vertices in upper lifted cube.
    NTYPE num_vertices_in_upper_lifted;

  public:

    typedef ISOV_TYPE IVOL_VERTEX_INDEX_TYPE;
    typedef ISOP_TYPE IVOL_POLY_INDEX_TYPE;
    typedef ISOP_TYPE ISO_POLY_INDEX_TYPE;
    typedef FBITS_TYPE FACET_BITS_TYPE;
    typedef typename 
    IVOLDUAL_POLY_VERTEX_INFO<ISOV_TYPE,NTYPE>::RELATIVE_POSITION
    RELATIVE_POSITION;

    /// Lower isosurface table index.
    TI_TYPE lower_isosurface_table_index;

    /// Upper isosurface table index.
    TI_TYPE upper_isosurface_table_index;

    IVOLDUAL_POLY_VERTEX_INFO<ISOV_TYPE,NTYPE> * poly_vertex_info;
    IVOLDUAL_POLY_EDGE_INFO<ISOV_TYPE> * poly_edge_info;
    IVOLDUAL_IVOLV_INFO * ivolv_info;

    /// True if interval volume is non-manifold.
    bool is_non_manifold;


    // Ambiguous/active facet information

    /// True if configuration is ambiguous.
    bool is_ambiguous;

    /// Number of ambiguous facets.
    NTYPE num_ambiguous_facets;

    /// Number of active facets (dual to some interval volume edge.)
    NTYPE num_active_facets;

    /*!
     *  @brief Integer representing set of ambiguous facets.
     *  - k'th bit of ambiguous_facet is 1 if facet k is ambiguous.
     */
    FBITS_TYPE ambiguous_facet_bits;    

    /// @brief Integer representing set of ambiguous facets in lower lifted cube.
    /// - k'th bit of ambiguous_facet is 1 if facet k is ambiguous.
    FBITS_TYPE ambiguous_facet_bits_in_lower_lifted;    

    /// @brief Integer representing set of ambiguous facets in upper lifted cube.
    /// - k'th bit of ambiguous_facet is 1 if facet k is ambiguous.
    FBITS_TYPE ambiguous_facet_bits_in_upper_lifted;


  public:

    IVOLDUAL_TABLE_ENTRY();      ///< constructor
    ~IVOLDUAL_TABLE_ENTRY();     ///< destructor

    // get functions

    NTYPE NumVerticesInLowerLifted() const
    { return(num_vertices_in_lower_lifted); }

    NTYPE NumVerticesInUpperLifted() const
    { return(num_vertices_in_upper_lifted); }

    /// Return number of dual interval volume vertices in cube.
    NTYPE NumVertices() const
    { return(NumVerticesInLowerLifted() + NumVerticesInUpperLifted()); }

    /// @brief True, if grid vertex is below interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename ITYPE>
    bool IsBelowIntervalVolume(const ITYPE iv) const
    { return(poly_vertex_info[iv].IsBelowIntervalVolume()); }

    /// @brief True, if grid vertex is in interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename ITYPE>
    bool IsInIntervalVolume(const ITYPE iv) const
    { return(poly_vertex_info[iv].IsInIntervalVolume()); }

    /// @brief True, if grid vertex is above interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename ITYPE>
    bool IsAboveIntervalVolume(const ITYPE iv) const
    { return(poly_vertex_info[iv].IsAboveIntervalVolume()); }

    /// @brief Return the relative position.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename ITYPE>
    RELATIVE_POSITION RelativePosition(const ITYPE iv) const
    { return(poly_vertex_info[iv].RelativePosition()); }

    /// @brief Return vertex type of grid vertex.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename ITYPE>
    NTYPE VertexType(const ITYPE iv) const
    { return(poly_vertex_info[iv].VertexType()); }

    /// Return lower isosurface table index.
    TI_TYPE LowerIsosurfaceTableIndex() const
    { return(lower_isosurface_table_index); }

    /// Return upper isosurface table index.
    TI_TYPE UpperIsosurfaceTableIndex() const
    { return(upper_isosurface_table_index); }

    /// True, if interval volume is non-manifold.
    bool IsNonManifold() const
    { return(is_non_manifold); }

    /// True, if configuration is ambiguous.
    bool IsAmbiguous() const
    { return(is_ambiguous); }

    /// Return number of ambiguous facets.
    NTYPE NumAmbiguousFacets() const
    { return(num_ambiguous_facets); };

    /// Return ambiguous facet bits.
    FBITS_TYPE AmbiguousFacetBits() const
    { return(ambiguous_facet_bits); };

    /// @brief Return ambiguous facet bits in lower lifted cube.
    /// - Note: Ignores facets orthogonal to lifting direction.
    FBITS_TYPE AmbiguousFacetBitsInLowerLifted() const
    { return(ambiguous_facet_bits_in_lower_lifted); };

    /// @brief Return ambiguous facet bits in upper lifted cube.
    /// - Note: Ignores facets orthogonal to lifting direction.
    FBITS_TYPE AmbiguousFacetBitsInUpperLifted() const
    { return(ambiguous_facet_bits_in_upper_lifted); };

    /// True, if facet jf is ambiguous.
    template <typename FTYPE>
    bool IsFacetAmbiguous(const FTYPE jf) const
    { return((ambiguous_facet_bits & (FBITS_TYPE(1) << jf)) != 0); };

    /// @brief True, if facet jf in lower lifted cube is ambiguous.
    /// - Note: Ignores facets orthogonal to lifting direction.
    template <typename FTYPE>
    bool IsFacetInLowerLiftedAmbiguous(const FTYPE jf) const
    { return((ambiguous_facet_bits_in_lower_lifted & 
              (FBITS_TYPE(1) << jf)) != 0); };

    /// @brief True, if facet jf in upper lifted cube is ambiguous.
    /// - Note: Ignores facets orthogonal to lifting direction.
    template <typename FTYPE>
    bool IsFacetInUpperLiftedAmbiguous(const FTYPE jf) const
    { return((ambiguous_facet_bits_in_upper_lifted & 
              (FBITS_TYPE(1) << jf)) != 0); };

    /// Create interval volume vertices.
    template <typename NUMV0_TYPE, typename NUMV1_TYPE>
    void CreateIntervalVolumeVertices
    (const NUMV0_TYPE numv_in_lower_lifted, 
     const NUMV1_TYPE numv_in_upper_lifted);

    /// Allocate poly_vertex_info[] and poly_edge_info[].
    template <typename NUMV_TYPE, typename NUME_TYPE>
    void Allocate(const NUMV_TYPE num_poly_vertices,
                  const NUME_TYPE num_poly_edges);


    /// @brief Free all memory.
    void FreeAll();                            // free all memory

    /// @brief Check data structure.
    /// @param num_mesh_poly_edges Number of edges of mesh polytope.
    /// - Currently ignored, but included for compatibility
    ///   with ISODUAL_TABLE_ENTRY.
    bool Check(const NTYPE num_mesh_poly_edges, IJK::ERROR & error_msg) const;
    
  };


  // **************************************************
  // CLASS ISODUAL_CUBE_FACE_INFO
  // *************************************************

  /// Class with routines to query number of positive or negative vertices in a cube facet and number of active facets.
  template <typename DTYPE, typename NTYPE, typename VTYPE>
  class ISODUAL_CUBE_FACE_INFO:public IJK::CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE> {

  public:
    ISODUAL_CUBE_FACE_INFO() {};
    ISODUAL_CUBE_FACE_INFO(const DTYPE dimension):
      IJK::CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>(dimension) {};

    /*!
     *  @brief Compute number of '-' vertices and number of '+' vertices in cube facet.
     *  @param[out] num_zeros Number of '-' vertices.
     *  @param[out] num_ones Number of '+' vertices.
     */
    template <typename TI_TYPE>
    void ComputeNumCubeFacetBits
    (const TI_TYPE ientry, const NTYPE ifacet,
     NTYPE & num_zeros, NTYPE & num_ones) const;

    /// Return true if cube facet has at least one '-' vertex and at least one '+' vertex.
    template <typename TI_TYPE>
    bool IsCubeFacetActive(const TI_TYPE ientry, const NTYPE ifacet) const;

    /// Return number of cube facets with at least one '-' vertex and at least one '+' vertex.
    template <typename TI_TYPE>
    NTYPE ComputeNumActiveCubeFacets(const TI_TYPE ientry) const;
  };


  // **************************************************
  // ISODUAL TABLE BASE
  // **************************************************

  /*!
   *  @brief Dual isosurface lookup table.
   *  - Stores isosurface vertices and incident faces for each configuration 
   *    of +/- labels at polytope vertex.
   */
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class ISODUAL_TABLE_BASE {

  public:

    /// @brief Index of entry in isosurface lookup table.
    /// - Define within ISODUAL_TABLE_BASE for use in templates.
    typedef TI_TYPE TABLE_INDEX;    

    typedef DTYPE DIMENSION_TYPE;
    typedef NTYPE NUMBER_TYPE;
    typedef typename ENTRY_TYPE::ISO_POLY_INDEX_TYPE
    ISO_POLY_INDEX_TYPE;
    
  protected:

    DTYPE dimension;               ///< Dimension.
    NTYPE num_poly_vertices;       ///< Number of polytope vertices.
    NTYPE num_poly_edges;          ///< Number of polytope edges;
    NTYPE num_poly_facets;          ///< Number of polytope facets;
    ENTRY_TYPE * entry;            ///< Array of dual isosurface table entries.
    TI_TYPE num_table_entries;     ///< Number of entries in table.

    /// Maximum number of vertices allowed for cube.
    NTYPE max_num_vertices; 

    /// True, if array entry[] is allocated.
    bool is_table_allocated;  

    /// Initialization routine.
    template <typename DTYPE2>
    void Init(const DTYPE2 dimension);


  public:
    ISODUAL_TABLE_BASE();
    template <typename DTYPE2>
    ISODUAL_TABLE_BASE(const DTYPE2 d);
    ~ISODUAL_TABLE_BASE();             ///< Destructor

    // Get functions.
    DTYPE Dimension() const            /// Return dimension.
    { return(dimension); };
    NTYPE NumPolyVertices() const      /// Return number of polytope vertices.
    { return(num_poly_vertices); };
    NTYPE NumPolyEdges() const         /// Return number of polytope edges.
    { return(num_poly_edges); };
    NTYPE NumPolyFacets() const        /// Return number of polytope facets.
    { return(num_poly_facets); };

    /// Return number of lookup table entries.
    TI_TYPE NumTableEntries() const { return(num_table_entries); };

    /// Return complement of table index it
    template <typename TI_TYPE2>
    TI_TYPE2 Complement(const TI_TYPE2 it) const
    { return(compute_complement(it, num_table_entries)); }

    /// Return number of vertices in isosurface patch for table entry \a it.
    template <typename TI_TYPE2>
    NTYPE NumIsoVertices(const TI_TYPE2 it) const
    { return(entry[it].NumVertices()); }; 

    /*!
     *  @brief Return true if edge ke is bipolar.
     *  @param it Index of table entry.
     *  @param ke Polytope edge ke, dual to isosurface face ke.
     */
    template <typename TI_TYPE2>
    bool IsBipolar(const TI_TYPE2 it, const NTYPE ke) const
    { return(this->entry[it].is_bipolar[ke]); };

    /*!
     *  @brief Return index of isosurface vertex incident on face kf.
     *  - Undefined if polytope edge k is not bipolar.
     *  @param it Index of table entry.
     *  @param kf Isosurface face kf, dual to polytope edge kf.
     */
    template <typename TI_TYPE2>
    inline NTYPE IncidentIsoVertex
    (const TI_TYPE2 it, const NTYPE kf) const
    { return entry[it].incident_isovertex[kf]; };

    /*!
     *  @brief Return pointer to index of first isosurface polytope
     *    incident on isosurface vertex kv.
     *  @param it Index of table entry.
     *  @param kv Isosurface vertex index.
     */
    template <typename TI_TYPE2>
    ISO_POLY_INDEX_TYPE const * FirstIncidentIsoPoly
    (const TI_TYPE2 it, const NTYPE kv) const
    { return(entry[it].first_incident_isopoly[kv]); };

    /*!
     *  @brief Return number of isosurface polytopes
     *    incident on isosurface vertex kv.
     *  @param it Index of table entry.
     *  @param kv Isosurface vertex index.
     */
    template <typename TI_TYPE2>
    NTYPE NumIncidentIsoPoly(const TI_TYPE2 it, const NTYPE kv) const
    { return(entry[it].num_incident_isopoly[kv]); };
    
    /// @brief Return true if vertex iv is positive.
    /// @param iv Vertex index.
    template <typename TI_TYPE2>
    bool IsPositive(const TI_TYPE2 it, const NTYPE iv) const;

    /*!
     *  @brief Return maximum number of polytope vertices permitted in any table.
     *  - Note: Even tables for polytopes of this size are probably impossible 
     *    to compute/store.
     */
    NTYPE MaxNumVertices() const { return(max_num_vertices); };

    /// Return true if table memory is allocated.
    bool IsTableAllocated() const
    { return(is_table_allocated); };

    /// Return reference to entry it.
    template <typename TI_TYPE2>
    const ENTRY_TYPE & Entry(const TI_TYPE2 it) const
    { return(this->entry[it]); }

    // Set functions.
    template <typename DTYPE2>
    void SetDimension(const DTYPE2 d);
    void SetNumPolyVertices(const NTYPE num_vertices);
    void SetNumPolyEdges(const NTYPE num_edges);
    void SetNumPolyFacets(const NTYPE num_facets);

    /// Set polytope to be a cube.
    template <typename DTYPE2>
    void SetToCube(const DTYPE2 d);

    /// Allocate table
    template <typename NTYPE2>
    void SetNumTableEntries(const NTYPE2 num_table_entries);

    // Check functions
    template <typename DTYPE2>
    bool CheckDimension(const DTYPE2 d) const;
    bool CheckDimension() const
    { return(CheckDimension(Dimension())); };
    bool CheckTable(IJK::ERROR & error_msg) const;
    bool Check(IJK::ERROR & error_msg) const;

    virtual void FreeAll();                     /// Free all memory.
  };


  // **************************************************
  // ISODUAL CUBE TABLE
  // **************************************************

  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class ISODUAL_CUBE_TABLE:
    public ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {

  protected:
    /// If true, separate negative vertices.
    bool flag_separate_neg;

    /// If true, always separate two diagonally opposite positive or negative vertices.
    bool flag_always_separate_opposite;     

    /*!
     *  @brief Create table entries.
     *  @param flag_separate_opposite If true, always separate two
     *         diagonally opposite positive or negative vertices.
     */
    void CreateTableEntries
      (const bool flag_separate_neg, const bool flag_separate_opposite);

    /// @brief Set incident isosurface polytopes for each table entry.
    void SetIncidentIsoPoly();

  public:
    ISODUAL_CUBE_TABLE() {};

    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE(const DTYPE2 dimension) :
      ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> (dimension)
    { Create(dimension); }
      
    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE
    (const DTYPE2 dimension, const bool flag_separate_opposite) :
      ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Create(dimension, flag_separate_opposite); }

    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE
      (const DTYPE2 dimension, const bool flag_separate_neg,
       const bool flag_separate_opposite) :
        ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Create(dimension, flag_separate_neg, flag_separate_opposite); }


    // Get functions.

    template <typename DTYPE2>
    void Create(const DTYPE2 dimension);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_opposite);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_neg, 
                const bool flag_separate_opposite);
  };


  // **************************************************
  // ISODUAL TABLE AMBIG BASE
  // **************************************************

  //! Isodual table plus ambiguity information.
  //! @tparam TI_TYPE Table index type.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class ISODUAL_TABLE_AMBIG_BASE:
    public ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {

  public:
    typedef typename ENTRY_TYPE::FACET_BITS_TYPE FACET_BITS_TYPE;

  public:
    
    // constructors
    ISODUAL_TABLE_AMBIG_BASE() {};

    template <typename DTYPE2>
    ISODUAL_TABLE_AMBIG_BASE(const DTYPE2 dimension) :
      ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension) {};


    // get functions

    /// Return true if cube in configuration it is ambiguous.
    bool IsAmbiguous(const TI_TYPE it) const
    { return(this->entry[it].IsAmbiguous()); }

    /// Return true if cube facet jf in configuration it is ambiguous.
    template <typename FTYPE>
    bool IsFacetAmbiguous
    (const TI_TYPE it, const FTYPE jf) const
    { return(this->entry[it].IsFacetAmbiguous(jf)); }

    FACET_BITS_TYPE AmbiguousFacetBits(const TI_TYPE it) const
    { return(this->entry[it].AmbiguousFacetBits()); }

    /// Return number of ambiguous facets.
    NTYPE NumAmbiguousFacets(const TI_TYPE it) const
    { return(this->entry[it].NumAmbiguousFacets()); }

    /// Return number of active facets.
    NTYPE NumActiveFacets(const TI_TYPE it) const
    { return(this->entry[it].NumActiveFacets()); }

  };


  // **************************************************
  // ISODUAL CUBE TABLE AMBIG
  // **************************************************

  /// @brief Isodual cube table plus ambiguity information.
  /// @tparam TI_TYPE Table index type.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class ISODUAL_CUBE_TABLE_AMBIG:
    public ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {

  public:
    typedef typename ENTRY_TYPE::FACET_BITS_TYPE FACET_BITS_TYPE;
    typedef NTYPE NUMBER_TYPE;
    typedef DTYPE DIMENSION_TYPE;

  protected:
    /// If true, separate negative vertices.
    bool flag_separate_neg;

    /// @brief If true, always separate two diagonally opposite positive or negative vertices.
    bool flag_always_separate_opposite;     

    /*!
     *  @brief Create table entries.
     *  @param flag_separate_opposite If true, always separate two
     *         diagonally opposite positive or negative vertices.
     */
    void CreateTableEntries
    (const bool flag_separate_neg, const bool flag_separate_opposite);


  public:

    ISODUAL_CUBE_TABLE_AMBIG() {};

    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE_AMBIG(const DTYPE2 dimension) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
      (dimension)
    { Create(dimension); }
      
    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE_AMBIG
    (const DTYPE2 dimension, const bool flag_separate_opposite) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Create(dimension, flag_separate_opposite); }

    template <typename DTYPE2>
    ISODUAL_CUBE_TABLE_AMBIG
    (const DTYPE2 dimension, const bool flag_separate_neg,
     const bool flag_separate_opposite) :
      ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension)
    { Create(dimension, flag_separate_neg, flag_separate_opposite); }

    // Get functions.

    /*!
     *  @brief Return true if edge ke is bipolar.
     *  @param it Index of table entry.
     *  @param ke Polytope edge ke, dual to isosurface face ke.
     */
    template <typename TI_TYPE2>
    bool IsBipolar(const TI_TYPE2 it, const NTYPE ke) const
    { return(this->entry[it].is_bipolar[ke]); };

    /// Return true if separate negative vertices.
    bool FlagSeparateNeg() const
    { return(flag_separate_neg); }

    /*!
     *  @brief Return number of isosurface polytopes incident on isov, and containing edge dual to facet kf.
     *  @tparam CUBE_TYPE Either IJK::CUBE_FACE_INFO, or a class
     *    derived from IJK::CUBE_FACE_INFO, or a class with similar
     *    vertices, edges, and facets.
     *  @param cube Cube used to determine isosurface polytopes
     *    dual to facet edges.
     */
    template <typename IT_TYPE2, typename ITYPE3, typename ITYPE4,
	      typename CUBE_TYPE>
    NTYPE CountNumIncidentIsoPolyDualToFacetEdges
    (const IT_TYPE2 it, const ITYPE3 isov, const ITYPE4 kf,
     const CUBE_TYPE & cube) const;

    /*!
     *  @brief Return number of isosurface edges incident on isov, and dual to facet kf.
     *  - Equals (CountNumIncidentIsoPolyDualToEdges()/2).
     *  @tparam CUBE_TYPE Either IJK::CUBE_FACE_INFO, or a class
     *    derived from IJK::CUBE_FACE_INFO, or a class with similar
     *    vertices, edges, and facets.
     *  @param cube Cube used to determine isosurface polytopes
     *    dual to facet edges.
     */
    template <typename IT_TYPE2, typename ITYPE3, typename ITYPE4,
	      typename CUBE_TYPE>
    NTYPE CountNumIncidentIsoEdgesDualToFacet
    (const IT_TYPE2 it, const ITYPE3 isov, const ITYPE4 kf,
     const CUBE_TYPE & cube) const
    { 
      return
	(CountNumIncidentIsoPolyDualToFacetEdges(it,isov,kf,cube)/2); 
    }

    /*!
     *  @brief Get orthogonal direction and side of (first) ambiguous facet.
     *  - Returns orthogonal direction and side of first ambiguous facet,
     *    if cube has more than one ambiguous facet.
     *  - Undefined, if cube has no ambiguous facets.
     *  @pre Cube has at least one ambiguous facet.
     *  @param[out] orth_dir Direction orthgonal to ambiguous facet.
     *  @param[out] side Side of cube. 0 or 1.
     */
    template <typename IT_TYPE2, typename DTYPE2,
              typename SIDE_TYPE>
    void GetAmbiguousFacet
    (const IT_TYPE2 it, DTYPE2 & orth_dir, SIDE_TYPE & side) const;

    /// Return number of bipolar edges incident on vertex.
    template <typename TI_TYPE2, typename VTYPE>
    NTYPE ComputeNumIncidentBipolarEdges
    (const TI_TYPE2 it, const VTYPE iv) const;

    template <typename DTYPE2>
    void Create(const DTYPE2 dimension);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_opposite);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_neg, 
                const bool flag_separate_opposite);
  };


  // **************************************************
  // IVOLDUAL CUBE TABLE BASE
  // **************************************************


  //! Base class for interval volume cube table.<br>
  //! Does not include any routines to create table entries.
  //! @tparam TI_TYPE Table index type.
  template <const int NUM_VERTEX_TYPES,
            typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class IVOLDUAL_CUBE_TABLE_BASE:
    public ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {

  public:
    typedef ISODUAL_CUBE_FACE_INFO<DTYPE,NTYPE,NTYPE> CUBE_TYPE;
    typedef typename ENTRY_TYPE::IVOL_VERTEX_INDEX_TYPE 
    IVOL_VERTEX_INDEX_TYPE;

  protected:
    CUBE_TYPE cube;

  protected:
    void Init();

    //! If true, separate negative vertices.
    bool flag_separate_neg;


  public:
    typedef typename ENTRY_TYPE::FACET_BITS_TYPE FACET_BITS_TYPE;
    typedef typename ENTRY_TYPE::RELATIVE_POSITION RELATIVE_POSITION;


  public:

    /// Constructor
    IVOLDUAL_CUBE_TABLE_BASE() { Init(); };

    template <typename DTYPE2>
    IVOLDUAL_CUBE_TABLE_BASE(const DTYPE2 dimension);

    template <typename DTYPE2>
    IVOLDUAL_CUBE_TABLE_BASE
    (const DTYPE2 dimension, const bool flag_separate_neg);

    // Get functions.

    /// Return number of vertex types.
    NTYPE NumVertexTypes() const
    {  return(NUM_VERTEX_TYPES); }

    /// Return flag_separate_neg.
    bool FlagSeparateNeg() const
    { return(flag_separate_neg); }

    /// Return true if interval volume is non-manifold.
    template <typename TI_TYPE2>
    bool IsNonManifold(const TI_TYPE2 ientry) const
    { return(this->entry[ientry].IsNonManifold()); }

    /// Return true if edge ie has dual interval volume polytope.
    template <typename TI_TYPE2, typename NTYPE2>
    bool EdgeHasDualIVolPoly(const TI_TYPE2 ientry, const NTYPE2 ie) const
    { return(this->entry[ientry].poly_edge_info[ie].flag_dual_ivolpoly); }

    /// Return lower incident isovertex on polytope dual to edge ie.
    template <typename TI_TYPE2, typename NTYPE2>
    NTYPE LowerIncident(const TI_TYPE2 ientry, const NTYPE2 ie) const
    { return(this->entry[ientry].poly_edge_info[ie].lower_incident_isovertex); }

    /// Return upper incident isovertex on polytope dual to edge ie.
    template <typename TI_TYPE2, typename NTYPE2>
    NTYPE UpperIncident(const TI_TYPE2 ientry, const NTYPE2 ie) const
    { return(this->entry[ientry].poly_edge_info[ie].upper_incident_isovertex); }

    /// Return incident interval volume vertex on polytope dual to cube vertex iv.
    template <typename TI_TYPE2, typename NTYPE2>
    NTYPE IncidentIVolVertex(const TI_TYPE2 ientry, const NTYPE2 iv) const
    { return(this->entry[ientry].poly_vertex_info[iv].incident_ivol_vertex); }

    /// Return true if isov is lower vertex.
    template <typename TI_TYPE2, typename ISOV2_TYPE>
    bool IsLowerVertex(const TI_TYPE2 ientry, const ISOV2_TYPE isov) const
    { return(this->entry[ientry].IsLowerVertex(isov)); }

    /// Return number of interval volume vertices.
    template <typename TI_TYPE2>
    NTYPE NumIVolVertices(const TI_TYPE2 ientry) const
    { return(this->entry[ientry].NumVertices()); }

    /// Return number of interval volume vertices in lower lifted cube.
    template <typename TI_TYPE2>
    NTYPE NumVerticesInLowerLifted(const TI_TYPE2 ientry) const
    { return(this->entry[ientry].NumVerticesInLowerLifted()); }

    /// Return number of interval volume vertices in upper lifted cube.
    template <typename TI_TYPE2>
    NTYPE NumVerticesInUpperLifted(const TI_TYPE2 ientry) const
    { return(this->entry[ientry].NumVerticesInUpperLifted()); }

    /// True, if interval volume vertex was generated in lower lifted cube.
    template <typename TI_TYPE2, typename VTYPE>
    bool IsInLowerLiftedCube
    (const TI_TYPE2 ientry, const VTYPE ivolv) const
    { return(this->entry[ientry].ivolv_info[ivolv].flag_in_lower_lifted_cube); }

    /// @brief True, if grid vertex is below interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    bool IsBelowIntervalVolume(const TI_TYPE2 ientry, const VTYPE iv) const
    { return(this->entry[ientry].IsBelowIntervalVolume(iv)); }

    /// @brief True, if grid vertex is in interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    bool IsInIntervalVolume(const TI_TYPE2 ientry, const VTYPE iv) const
    { return(this->entry[ientry].IsInIntervalVolume(iv)); }

    /// @brief True, if grid vertex is above interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    bool IsAboveIntervalVolume(const TI_TYPE2 ientry, const VTYPE iv) const
    { return(this->entry[ientry].IsAboveIntervalVolume(iv)); }

    /// @brief Return the relative position.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    RELATIVE_POSITION RelativePosition
    (const TI_TYPE2 ientry, const VTYPE iv) const
    { return(this->entry[ientry].RelativePosition(iv)); }

    /// @brief True, if some neighbor of grid vertex is below interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    bool IsSomeNeighborBelowIntervalVolume
    (const TI_TYPE2 ientry, const VTYPE iv) const;

    /// @brief True, if some neighbor of grid vertex is above interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    bool IsSomeNeighborAboveIntervalVolume
    (const TI_TYPE2 ientry, const VTYPE iv) const;


    /// @brief True, if some neighbor of grid vertex is below interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    NTYPE NumNeighborsBelowIntervalVolume
    (const TI_TYPE2 ientry, const VTYPE iv) const;

    /// @brief True, if some neighbor of grid vertex is above interval volume.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    NTYPE NumNeighborsAboveIntervalVolume
    (const TI_TYPE2 ientry, const VTYPE iv) const;

    /// @brief Return vertex type of grid vertex.
    /// @param iv Index of the vertex in the grid polytope.
    template <typename TI_TYPE2, typename VTYPE>
    NTYPE VertexType(const TI_TYPE2 ientry, const VTYPE iv) const
    { return(this->entry[ientry].VertexType(iv)); }

    /// Return lower isosurface table index.
    template <typename TI_TYPE2>
    TI_TYPE LowerIsosurfaceTableIndex(const TI_TYPE2 ientry) const
    { return(this->entry[ientry].LowerIsosurfaceTableIndex()); }

    /// Return upper isosurface table index.
    template <typename TI_TYPE2>
    TI_TYPE UpperIsosurfaceTableIndex(const TI_TYPE2 ientry) const
    { return(this->entry[ientry].UpperIsosurfaceTableIndex()); }

    /// @brief Return true if interval volume vertex ivolv is on lower isosurface.
    /// @param ivolv Index of the interval volume vertex.
    template <typename TI_TYPE2, typename VTYPE>
    bool OnLowerIsosurface(const TI_TYPE2 ientry, const VTYPE ivolv) const
    { return(this->entry[ientry].ivolv_info[ivolv].flag_lower_isosurface); }

    /// @brief Return true if interval volume vertex ivolv is on upper isosurface.
    /// @param ivolv Index of the interval volume vertex.
    template <typename TI_TYPE2, typename VTYPE>
    bool OnUpperIsosurface(const TI_TYPE2 ientry, const VTYPE ivolv) const
    { return(this->entry[ientry].ivolv_info[ivolv].flag_upper_isosurface); }

    /// Return true if lower isosurface intersects edge ie.
    template <typename TI_TYPE2, typename NTYPE2>
    bool DoesLowerIsosurfaceIntersectEdge
    (const TI_TYPE2 ientry, const NTYPE2 ie) const
    { return(does_lower_isosurface_intersect_edge(*this, ientry, ie)); }

    /*!
     *  @brief Return true if lower isosurface intersects edge ie.
     *  - If true, return interval volume vertex on lower isosurface 
     *    quadrilateral which is dual to edge ie.
     */
    template <typename TI_TYPE2, typename NTYPE2, typename IVOLV_TYPE2>
    bool DoesLowerIsosurfaceIntersectEdge
    (const TI_TYPE2 ientry, const NTYPE2 ie, IVOLV_TYPE2 & ivolv) const
    { return(does_lower_isosurface_intersect_edge(*this, ientry, ie, ivolv)); }

    /// @brief Return true if upper isosurface intersects edge ie.
    template <typename TI_TYPE2, typename NTYPE2, typename IVOLV_TYPE2>
    bool DoesUpperIsosurfaceIntersectEdge
    (const TI_TYPE2 ientry, const NTYPE2 ie, IVOLV_TYPE2 & ivolv) const
    { return(does_upper_isosurface_intersect_edge(*this, ientry, ie, ivolv)); }

    /*!
     *  @brief Return true if upper isosurface intersects edge ie.
     *  - If true, return interval volume vertex on upper isosurface 
     *    quadrilateral which is dual to edge ie.
     */
    template <typename TI_TYPE2, typename NTYPE2>
    bool DoesUpperIsosurfaceIntersectEdge
    (const TI_TYPE2 ientry, const NTYPE2 ie) const
    { return(does_upper_isosurface_intersect_edge(*this, ientry, ie)); }

    /// Return ambiguous facet bits in lower lifted cube.
    FACET_BITS_TYPE AmbiguousFacetBitsInLowerLifted(const TI_TYPE it) const
    { return(this->entry[it].AmbiguousFacetBitsInLowerLifted()); }

    /// @brief Return ambiguous facet bits in upper lifted cube.
    FACET_BITS_TYPE AmbiguousFacetBitsInUpperLifted(const TI_TYPE it) const
    { return(this->entry[it].AmbiguousFacetBitsInUpperLifted()); }

    /// Return true if facet jf in lower lifted cube is ambiguous.
    template <typename FTYPE>
    bool IsFacetInLowerLiftedAmbiguous
    (const TI_TYPE it, const FTYPE jf) const
    { return(this->entry[it].IsFacetInLowerLiftedAmbiguous(jf)); }

    /// Return true if facet jf in upper lifted cube is ambiguous.
    template <typename FTYPE>
    bool IsFacetInUpperLiftedAmbiguous
    (const TI_TYPE it, const FTYPE jf) const
    { return(this->entry[it].IsFacetInUpperLiftedAmbiguous(jf)); }

    /// Return reference to cube.
    const ISODUAL_CUBE_FACE_INFO<DTYPE,NTYPE,NTYPE> & Cube() const
    { return(cube); }

    /// Return true if configuration contains a vertex with label vtype.
    template <const int vtype>
    bool ContainsVertexType(const TI_TYPE it) const;

    /// Return true if configuration contains a vertex with label -.
    bool ContainsNegVertex(const TI_TYPE it) const
    { return(ContainsVertexType<0>(it)); }

    /// Return true if configuration contains a vertex with label +.
    bool ContainsPosVertex(const TI_TYPE it) const
    { return(ContainsVertexType<NUM_VERTEX_TYPES-1>(it)); }

    /// Return true if configuration contains a vertex with label I1.
    bool ContainsI1(const TI_TYPE it) const
    { return(ContainsVertexType<1>(it)); }

    /// Return true if configuration contains a vertex with label I2.
    bool ContainsI2(const TI_TYPE it) const
    { return(ContainsVertexType<2>(it)); }

    /// Redefine SetDimension to set cube dimension.
    template <typename DTYPE2>
    void SetDimension(const DTYPE2 d);

    /// @brief Check number of vertex types.
    /// - Should be 4.
    bool CheckNumVertexTypes(IJK::ERROR & error);

    /// Check that cube dimension matches table dimension.
    bool CheckCubeDimension(IJK::ERROR & error);

    /// Undefine IncidentIsoVertex
    template <typename TI_TYPE2>
    NTYPE IncidentIsoVertex
    (const TI_TYPE2 it, const NTYPE kf) const;
  };


  // **************************************************
  // IVOLDUAL CUBE TABLE
  // **************************************************


  /// @brief Interval volume cube table.
  /// @tparam TI_TYPE Table index type.
  template <const int NUM_VERTEX_TYPES,
            typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class IVOLDUAL_CUBE_TABLE:
    public IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES,DTYPE,NTYPE,
                                    TI_TYPE,ENTRY_TYPE> {

  protected:

    void Init();

    /// @brief Create table entries.
    /// @param flag_separate_neg Set flag_separate_neg to given value.
    void CreateTableEntries(const bool flag_separate_neg);

  public:

    /// Constructor
    IVOLDUAL_CUBE_TABLE() { Init(); };

    template <typename DTYPE2>
    IVOLDUAL_CUBE_TABLE(const DTYPE2 dimension);

    template <typename DTYPE2>
    IVOLDUAL_CUBE_TABLE
    (const DTYPE2 dimension, const bool flag_separate_neg);

    template <typename DTYPE2>
    void Create(const DTYPE2 dimension);
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_neg);
  };


  // **************************************************
  // IVOLDUAL CUBE DOUBLE TABLE
  // **************************************************

  /*!
   *  @brief Interval volume cube table storing both seprate negative and separate positive.
   *  - Interval volume cube table storing storing both interval volume patches
   *    which separate negative cube vertices and interval volume patches
   *    which separate positive cube vertices.  
   *    Equivalent, to two interval volume tables, one separating negative
   *    cube vertices and one separating positive cube vertices.
   *  - This table has double the entries of IVOLDUAL_CUBE_TABLE.
   *  - Half the entries have interval volume patches which separate 
   *      negative cube vertices and half have interval volume patches
   *      which separate positive cube vertices.
   *  - If flag_separate_neg is true, then the first half of the table
   *      stores interval volume patches which separate negative cube vertices 
   *      and the second half stores interval volume patches which separate
   *      positive cube vertices.
   *  - If flag_separate_neg is false, then the first half of the table
   *      stores interval volume patches which separate positive cube vertices 
   *      and the second half stores interval volume patches which separate
   *      negative cube vertices.
   *  @tparam TI_TYPE Table index type.
   */
  template <const int NUM_VERTEX_TYPES,
            typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  class IVOLDUAL_CUBE_DOUBLE_TABLE:
    public IVOLDUAL_CUBE_TABLE_BASE
  <NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE> {

  protected:

    /*!
     *  @brief Number of table entries with isosurface patches 
     *    separating negative cube vertices.
     *  - Also equals number of table entries with isosurface patches
     *    separating positive cube vertices.
     */
    NTYPE num_negative_table_entries;

    /*!
     *  @brief Create table entries.
     *  @param flag_separate_neg If true, first half of table stores
     *         interval volume patches separating negative cube vertices.
     */
    void CreateTableEntries(const bool flag_separate_neg);

    /// Initialize.
    void Init();

  public:

    /// Constructor
    IVOLDUAL_CUBE_DOUBLE_TABLE() { Init(); };

    template <typename DTYPE2>
    IVOLDUAL_CUBE_DOUBLE_TABLE(const DTYPE2 dimension);

    template <typename DTYPE2>
    IVOLDUAL_CUBE_DOUBLE_TABLE
    (const DTYPE2 dimension, const bool flag_separate_neg);

    /// @brief Return number of lookup table entries
    ///   whose interval volume patch separates negative vertices.
    NTYPE NumNegativeTableEntries() const 
    { return(num_negative_table_entries); };

    /// @brief Return number of lookup table entries
    ///   whose interval volume patch separates positive vertices.
    /// - Note: This number equals the number of negative table entries.
    NTYPE NumPositiveTableEntries() const 
    { return(num_negative_table_entries); };

    /// @brief Replace function SetNumTableEntries
    ///   to also set num_negative_table_entries.
    template <typename NTYPE2>
    void SetNumTableEntries(const NTYPE2 num_table_entries);

    /*!
     *  @brief Return opposite table index.
     *  - Opposite table index has same configuration of vertex labels
     *    as table_index.
     *  - Opposite table index separates positive vertices if table_index
     *    separates negative vertices.
     *  - Opposite table index separates negative vertices if table_index
     *    separates negative vertices.
     */
    template <typename TI_TYPE2>
    TI_TYPE OppositeTableIndex(const TI_TYPE2 table_index) const
    { return((table_index +
              NumNegativeTableEntries())%(this->NumTableEntries())); };

    /// @brief Create interval volume lookup table.
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension);

    /// @overload
    /// @brief Create interval volume lookup table. (Set flag_separate_neg.)
    /// - Version that sets flag_separate_neg.
    template <typename DTYPE2>
    void Create(const DTYPE2 dimension, const bool flag_separate_neg);
  };


  // **************************************************
  // CLASS FIND_COMPONENT
  // **************************************************

  /// Find connected component among cube vertices.
  template <typename DTYPE, typename NTYPE>
  class FIND_COMPONENT {

  protected:
    DTYPE dimension;
    NTYPE num_cube_vertices;
    bool * vertex_flag;
    NTYPE * component;

  public:
    typedef DTYPE DIMENSION_TYPE;         ///< Dimension type.
    typedef NTYPE NUMBER_TYPE;            ///< Number type.

  public:
    template <typename DTYPE2>
    FIND_COMPONENT(const DTYPE2 dimension);
    ~FIND_COMPONENT();

    // set functions
    template <typename TI_TYPE2>
    void SetVertexFlags(const TI_TYPE2 ival);
    void NegateVertexFlags();
    void ClearAll();

    // get functions
    DTYPE Dimension() const
    { return(dimension); }
    bool VertexFlag(const int i) const
    { return(vertex_flag[i]); }
    template <typename ITYPE>
    NTYPE Component(const ITYPE i) const
    { return(component[i]); }
    NTYPE NumCubeVertices() const
    { return(num_cube_vertices); }

    /// @brief Search starting at vertex i.
    /// @pre icomp is not zero.
    template <typename ITYPE0, typename ITYPE1>
    void Search(const ITYPE0 i, const ITYPE1 icomp);

    /// @brief Search facet starting at vertex i.
    /// @pre Facet kf contains vertex i.
    /// @pre icomp is not zero.
    template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
    void SearchFacet(const ITYPE0 kf, const ITYPE1 i, const ITYPE2 icomp);

    /// @brief Compute number of components.
    /// @param flag_positive If true, compute components of positive vertices.
    ///                      If false, compute components of negative vertices.
    template <typename TI_TYPE2>
    NTYPE ComputeNumComponents
    (const TI_TYPE2 ientry, const bool flag_positive);

    /// @brief Compute number of components in facet.
    /// @param flag_positive If true, compute components of positive vertices.
    ///                      If false, compute components of negative vertices.
    template <typename TI_TYPE2, typename ITYPE>
    NTYPE ComputeNumComponentsInFacet
    (const TI_TYPE2 ientry, const ITYPE kf, const bool flag_positive);
  };


  // **************************************************
  // ISODUAL TABLE ENTRY MEMBER FUNCTIONS
  // **************************************************

  // constructor
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE>
  ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE>::ISODUAL_TABLE_ENTRY()
  {
    num_vertices = 0;
    incident_isovertex = NULL;
    is_bipolar = NULL;
    incident_isopoly = NULL;
    first_incident_isopoly = NULL;
    num_incident_isopoly = NULL;
  }

  
  // Set arrays incident_isopoly and first_incident_isopoly.
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE>
  void ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE>::
  SetIncidentIsoPoly(const NTYPE num_mesh_poly_edges)
  {
    typedef ISOP_TYPE * ISOP_PTR_TYPE;

    FreeIncidentIsoPoly();
    
    if (NumVertices() <= 0) {
      // Nothing to set.
      return;
    }

    // Allocate incident arrays[].
    incident_isopoly = new ISOP_TYPE[num_mesh_poly_edges];
    first_incident_isopoly = new ISOP_PTR_TYPE[NumVertices()];
    num_incident_isopoly = new NTYPE[NumVertices()];

    ISOP_PTR_TYPE current_isopoly_ptr = incident_isopoly;
    for (NTYPE k = 0; k < NumVertices(); k++) {
      num_incident_isopoly[k] = 0;
      first_incident_isopoly[k] = current_isopoly_ptr;
      for (NTYPE ie = 0; ie < num_mesh_poly_edges; ie++) {
        if (IsBipolar(ie)) {
          if (IncidentIsoVertex(ie) == k) {
            // Add ie to incident_isopoly[].
            *current_isopoly_ptr = ie;
            current_isopoly_ptr++;
            num_incident_isopoly[k]++;
          }
        }

      }
    }

    if (current_isopoly_ptr > incident_isopoly + num_mesh_poly_edges) {
      IJK::PROCEDURE_ERROR error("ISODUAL_TABLE_ENTRY::SetIncidentIsoPoly");
      error.AddMessage
        ("Programming error. Added too many incident isosurface polytopes.");
      error.AddMessage("  Added ", current_isopoly_ptr-incident_isopoly,
                       " incident isosurface polytopes.");
      error.AddMessage("  At most ", num_mesh_poly_edges,
                       " are dual to edges of a single mesh element.");
      throw error;
    }
  }
  
  
  // destructor
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE>
  ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE>::~ISODUAL_TABLE_ENTRY()
  {
    FreeAll();
  }

  
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE>
  template <typename NUMV_TYPE, typename NUME_TYPE>
  void ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE>::
  Allocate(const NUMV_TYPE num_poly_vertices,
           const NUME_TYPE num_poly_edges)
  {
    FreeAll();

    incident_isovertex = new ISOV_TYPE[num_poly_edges];
    is_bipolar = new bool[num_poly_edges];
  }


  // @brief Free arrays incident_isopoly[] and first_incident_isopoly[].
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE>
  void ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE>::
  FreeIncidentIsoPoly()
  {
    if (incident_isopoly != NULL) {
      delete [] incident_isopoly;
      incident_isopoly = NULL;
    }

    if (first_incident_isopoly != NULL) {
      delete [] first_incident_isopoly;
      first_incident_isopoly = NULL;
    }

    if (num_incident_isopoly != NULL) {
      delete [] num_incident_isopoly;
      num_incident_isopoly = NULL;
    }
  }
  
  
  // Free all memory.
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE>
  void ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE>::FreeAll()
  {
    if (incident_isovertex != NULL) {
      delete [] incident_isovertex;
      incident_isovertex = NULL;
    }

    if (is_bipolar != NULL) {
      delete [] is_bipolar;
      is_bipolar = NULL;
    }

    FreeIncidentIsoPoly();

    num_vertices = 0;
  }


  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE>
  bool ISODUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE>::
  Check(const NTYPE num_mesh_poly_edges, IJK::ERROR & error_msg) const
  {
    if (num_vertices < 0) {
      error_msg.AddMessage
        ("Error.  Dual isosurface table entry contains negative number of isosurface vertices.");
      return false;
    }

    if (incident_isovertex == NULL) {
      error_msg.AddMessage
        ("Memory for incident isosurface vertices not allocated.");
      return false;
    }

    if (is_bipolar == NULL) {
      error_msg.AddMessage("Memory for bipolar edge flags not allocated.");
      return false;
    }


    if (NumVertices() > 0) {

      if (incident_isopoly == NULL) {
        error_msg.AddMessage
          ("Memory for incident isosurface polytopes not allocated.");
        return false;
      }

      if (num_incident_isopoly == NULL) {
        error_msg.AddMessage
          ("Memory for number of incident isosurface polytopes not allocated.");
        return false;
      }

      if (first_incident_isopoly == NULL) {
        error_msg.AddMessage
          ("Memory for first incident isosurface polytopes not allocated.");
        return false;
      }

      for (NTYPE k = 0; k < NumVertices(); k++) {
        ISOP_TYPE * current_ptr = FirstIncidentIsoPoly(k);
        for (NTYPE j = 0; j < NumIncidentIsoPoly(k); j++) {
          const NTYPE ie = *current_ptr;

          if (ie < 0) {
            error_msg.AddMessage
              ("Negative index ", ie, " in incident isosurface polytope.");
            return false;
          }
            
          
          if (ie >= num_mesh_poly_edges) {
            error_msg.AddMessage
              ("Incorrect index ", ie, " in incident isosurface polytope.");
            error_msg.AddMessage
              ("  Index must be less than ", num_mesh_poly_edges, ".");
            return false;
          }
          
          if (!IsBipolar(ie)) {
            error_msg.AddMessage("Edge ", ie, " is NOT bipolar.");
            error_msg.AddMessage
              ("  Isosurface mesh polytopes are only dual to bipolar edges.");
            return false;
          }

          if (IncidentIsoVertex(ie) != k) {
            error_msg.AddMessage
              ("Isosurface mesh polytope ", ie, " is incident on vertex ", k);
            error_msg.AddMessage
              ("  but IncidentIsoVertex(", ie, ") = ", IncidentIsoVertex(ie), ".");
            return false;
          }
        }
      }
    }
    
    return true;
  }

  
  // **************************************************
  // ISODUAL AMBIG TABLE ENTRY MEMBER FUNCTIONS
  // **************************************************

  // constructor
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FACET_BITS_TYPE>
  ISODUAL_AMBIG_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,FACET_BITS_TYPE>::
  ISODUAL_AMBIG_TABLE_ENTRY()
  {
    ClearAmbig();
  }

  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FACET_BITS_TYPE>
  void ISODUAL_AMBIG_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,
                                 FACET_BITS_TYPE>::
  ClearAmbig()
  {
    is_ambiguous = false;
    num_ambiguous_facets = 0;
    num_active_facets = 0;
    ambiguous_facet_bits = 0;
  }

  // destructor
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FACET_BITS_TYPE>
  ISODUAL_AMBIG_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,
                            FACET_BITS_TYPE>::
  ~ISODUAL_AMBIG_TABLE_ENTRY()
  {
    ClearAmbig();
  }


  // **************************************************
  // IVOLDUAL TABLE ENTRY MEMBER FUNCTIONS
  // **************************************************
 
  // IVOLDUAL_POLY_EDGE_INFO initalization routine
  template <typename ISOV_TYPE>
  void IVOLDUAL_POLY_EDGE_INFO<ISOV_TYPE>::Init()
  {
    flag_dual_ivolpoly = false;
  }

  template <typename ISOV_TYPE>
  template <typename ISOV_TYPE2>
  void IVOLDUAL_POLY_EDGE_INFO<ISOV_TYPE>::
  SetLowerIncident(const ISOV_TYPE2 ilower)
  {
    lower_incident_isovertex = ilower;
    flag_dual_ivolpoly = true;
  }

  template <typename ISOV_TYPE>
  template <typename ISOV_TYPE2>
  void IVOLDUAL_POLY_EDGE_INFO<ISOV_TYPE>::
  SetUpperIncident(const ISOV_TYPE2 iupper)
  {
    upper_incident_isovertex = iupper;
    flag_dual_ivolpoly = true;
  }

  // IVOLDUAL_POLY_VERTEX_INFO initalization routine
  template <typename IVOLV_TYPE, typename VERTEX_TYPE>
  void IVOLDUAL_POLY_VERTEX_INFO<IVOLV_TYPE,VERTEX_TYPE>::Init()
  {
    vertex_type = 0;
  }

  // Set incident_ivol_vertex.
  template <typename IVOLV_TYPE, typename VERTEX_TYPE>
  template <typename IVOLV_TYPE2>
  void IVOLDUAL_POLY_VERTEX_INFO<IVOLV_TYPE,VERTEX_TYPE>::
  SetIncident(const IVOLV_TYPE2 ivolv)
  {
    incident_ivol_vertex = ivolv;
  }

  // constructor
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FBITS_TYPE, typename TI_TYPE>
  IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,FBITS_TYPE,TI_TYPE>::
  IVOLDUAL_TABLE_ENTRY()
  {
    num_vertices_in_lower_lifted = 0;
    num_vertices_in_upper_lifted = 0;
    lower_isosurface_table_index = 0;
    upper_isosurface_table_index = 0;
    poly_vertex_info = NULL;
    poly_edge_info = NULL;
    ivolv_info = NULL;
    is_non_manifold = false;

    // Ambiguity or active facet information.
    is_ambiguous = false;
    num_ambiguous_facets = 0;
    num_active_facets = 0;
    ambiguous_facet_bits = 0;
    ambiguous_facet_bits_in_lower_lifted = 0;
    ambiguous_facet_bits_in_upper_lifted = 0;
  }

  // destructor
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FBITS_TYPE, typename TI_TYPE>
  IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,
                       FBITS_TYPE,TI_TYPE>::
  ~IVOLDUAL_TABLE_ENTRY()
  {
    FreeAll();
  }

  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FBITS_TYPE, typename TI_TYPE>
  template <typename NUMV_TYPE, typename NUME_TYPE>
  void IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,
                            FBITS_TYPE,TI_TYPE>::
  Allocate(const NUMV_TYPE num_poly_vertices,
           const NUME_TYPE num_poly_edges)
  {
    FreeAll();

    poly_vertex_info = 
      new IVOLDUAL_POLY_VERTEX_INFO<ISOV_TYPE,NTYPE>[num_poly_vertices];
    poly_edge_info = 
      new IVOLDUAL_POLY_EDGE_INFO<ISOV_TYPE>[num_poly_edges];
  }


  /// Create interval volume vertices.
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FBITS_TYPE, typename TI_TYPE>
  template <typename NUMV0_TYPE, typename NUMV1_TYPE>
  void IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,
                            FBITS_TYPE,TI_TYPE>::
  CreateIntervalVolumeVertices
  (const NUMV0_TYPE numv_in_lower_lifted, 
   const NUMV1_TYPE numv_in_upper_lifted)
  {
    const NTYPE numv = numv_in_lower_lifted + numv_in_upper_lifted;

    num_vertices_in_lower_lifted = numv_in_lower_lifted;
    num_vertices_in_upper_lifted = numv_in_upper_lifted;


    if (ivolv_info != NULL) {
      delete [] ivolv_info;
      ivolv_info = NULL;
    }

    if (numv > 0) 
      { ivolv_info = new IVOLDUAL_IVOLV_INFO[numv]; }

    for (NTYPE i = 0; i < numv; i++) {
      if (i < numv_in_lower_lifted) 
        { ivolv_info[i].flag_in_lower_lifted_cube = true; }
      else
        { ivolv_info[i].flag_in_lower_lifted_cube = false; }
    }
  }


  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FBITS_TYPE, typename TI_TYPE>
  bool IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,
                            FBITS_TYPE,TI_TYPE>::
  Check(const NTYPE num_mesh_poly_edges, IJK::ERROR & error_msg) const
  {
    if (num_vertices_in_lower_lifted < 0) {
      error_msg.AddMessage
        ("Programming error.  Dual interval volume table entry contains negative number");
      error_msg.AddMessage("  of interval volume vertices in lower lifted cube.");
      return(false);
    }

    if (num_vertices_in_upper_lifted < 0) {
      error_msg.AddMessage
        ("Programming error.  Dual interval volume table entry contains negative number");
      error_msg.AddMessage("  of interval volume vertices in upper lifted cube.");
      return(false);
    }

    if (poly_vertex_info == NULL) {
      error_msg.AddMessage
        ("Programming error.  Memory for isosurface elements dual to polytope vertices not allocated.");
      return(false);
    }

    if (poly_edge_info == NULL) {
      error_msg.AddMessage
        ("Programming errror.  Memory for isosurface elements dual to polytope edges not allocated.");
      return(false);
    }

    return(true);
  }

  // Free all memory.
  template <typename NTYPE, typename ISOV_TYPE, typename ISOP_TYPE,
            typename FBITS_TYPE, typename TI_TYPE>
  void IVOLDUAL_TABLE_ENTRY<NTYPE,ISOV_TYPE,ISOP_TYPE,
                            FBITS_TYPE,TI_TYPE>::FreeAll()
  {
    if (poly_vertex_info != NULL) {
      delete [] poly_vertex_info;
      poly_vertex_info = NULL;
    }

    if (poly_edge_info != NULL) {
      delete [] poly_edge_info;
      poly_edge_info = NULL;
    }

    if (ivolv_info != NULL) {
      delete [] ivolv_info;
      ivolv_info = NULL;
    }

    num_vertices_in_lower_lifted = 0;
    num_vertices_in_upper_lifted = 0;
    lower_isosurface_table_index = 0;
    upper_isosurface_table_index = 0;
  }


  // **************************************************
  // ISODUAL TABLE BASE MEMBER FUNCTIONS
  // **************************************************

  // default constructor. dimension = 3
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::ISODUAL_TABLE_BASE()
  {
    const DTYPE DIM3(3);

    Init(DIM3);
  }

  // Constructor.
  // d = dimension of space containing isosurface.  Should be 2, 3 or 4.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  ISODUAL_TABLE_BASE(const DTYPE2 d)
  {
    Init(d);
  }

  // Initialize
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Init(const DTYPE2 dimension)
  {
    max_num_vertices = 20;
    // Note: Even tables for polytopes of this size are probably impossible 
    //   to compute/store

    num_table_entries = 0;
    entry = NULL;
    is_table_allocated = false;

    SetDimension(dimension);
  }

  // destructor
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::~ISODUAL_TABLE_BASE()
  {
    FreeAll();
  }

  // Set dimension.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetDimension(const DTYPE2 dimension)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetDimension";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    this->dimension = dimension;

    if (!CheckDimension())
      throw IJK::PROCEDURE_ERROR(procname, "Illegal polytope dimension.");
  }

  // Set number of polytope vertices.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumPolyVertices(const NTYPE num_vertices)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetNumVertices";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    num_poly_vertices = num_vertices;
  }

  // Set number of polytope edges
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumPolyEdges(const NTYPE num_edges)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetNumEdges";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    num_poly_edges = num_edges;
  }

  // Set number of polytope edges
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumPolyFacets(const NTYPE num_facets)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetNumFacets";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    num_poly_facets = num_facets;
  }

  // Set polytope to be a cube.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetToCube(const DTYPE2 dimension)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetToCube";

    if (IsTableAllocated()) {
      throw IJK::PROCEDURE_ERROR
        (procname, "Dual isosurface table already allocated.");
    }

    this->SetDimension(dimension);

    const NTYPE num_vertices = IJK::compute_num_cube_vertices(dimension);
    this->SetNumPolyVertices(num_vertices);

    const NTYPE num_edges = IJK::compute_num_cube_edges(dimension);
    this->SetNumPolyEdges(num_edges);

    const NTYPE num_facets = IJK::compute_num_cube_facets(dimension);
    this->SetNumPolyFacets(num_facets);
  }


  // Allocate table
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename NTYPE2>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumTableEntries(const NTYPE2 num_table_entries)
  {
    const char * procname = "ISODUAL_TABLE_BASE::SetNumTableEntries";

    if (entry != NULL) delete [] entry;
    entry = NULL;
    this->num_table_entries = 0;
    is_table_allocated = false;

    entry = new ENTRY_TYPE[num_table_entries];
    if (entry == NULL && num_table_entries > 0)
      throw IJK::PROCEDURE_ERROR
        (procname, "Unable to allocate memory for dual isosurface table.");

    for (NTYPE2 i = 0; i < num_table_entries; i++) {
      entry[i].Allocate(NumPolyVertices(), NumPolyEdges());
    }

    this->num_table_entries = num_table_entries;
    is_table_allocated = true;
  }

  // Return true if vertex iv is positive.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename TI_TYPE2>
  bool ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  IsPositive(const TI_TYPE2 it, const NTYPE iv) const
  {
    const TI_TYPE2 mask = (TI_TYPE2(1) << iv);

    if ((it & mask) == 0) { return(false); }
    else { return(true); }
  }

  // check dimension
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  bool ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CheckDimension(const DTYPE2 d) const
  {
    if (d < 1)
      return(false);
    else
      return(true);
  }


  // check table
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  bool ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CheckTable(IJK::ERROR & error_msg) const
  {
    if (NumPolyVertices() >= CHAR_BIT*sizeof(TI_TYPE)) {
      error_msg.AddMessage("Too many polytope vertices");
      return(false);
    }

    if (NumPolyVertices() > MaxNumVertices()) {
      error_msg.AddMessage("Too many polytope vertices");
      return(false);
    }

    if (NumPolyVertices() < 1) {
      error_msg.AddMessage("Polytope must have at least one vertex.");
      return(false);
    }

    if (entry == NULL) {
      error_msg.AddMessage("Memory for dual isosurface table not allocated.");
      return(false);
    }

    for (TI_TYPE it = 0; it < NumTableEntries(); it++)
      if (!entry[it].Check(NumPolyEdges(), error_msg)) {
        error_msg.AddMessage
          ("Error detected at isosurface table entry ", it, ".");
        return(false);
      }

    return(true);
  }


  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  bool ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Check(IJK::ERROR & error_msg) const
  {
    if (!CheckTable(error_msg)) return(false);
    return(true);
  }


  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_TABLE_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::FreeAll()
  // free all memory
  {
    if (entry != NULL) {
      for (TI_TYPE i = 0; i < num_table_entries; i++)
        { entry[i].FreeAll(); }
      delete [] entry;
      entry = NULL;
    };
    num_table_entries = 0;
    is_table_allocated = false;
  }


  // **************************************************
  // ISODUAL CUBE TABLE MEMBER FUNCTIONS
  // **************************************************

  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_neg,
         const bool flag_separate_opposite)
  {
    this->SetToCube(dimension);

    TI_TYPE n = calculate_num_entries<TI_TYPE>(this->NumPolyVertices(), 2);
    this->SetNumTableEntries(n);
    CreateTableEntries(flag_separate_neg, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_opposite)
  {
    Create(dimension, true, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension)
  {
    Create(dimension, true);
  }


  // Create table entries.
  // @param flag_separate_neg  If true, separate negative vertices.
  // @param flag_separate_opposite If true, always separate two diagonally 
  //        opposite negative or positive vertices 
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CreateTableEntries
  (const bool flag_separate_neg, const bool flag_separate_opposite)
  {
    this->flag_separate_neg = flag_separate_neg;
    this->flag_always_separate_opposite = flag_separate_opposite;

    const DTYPE dimension = this->Dimension();

    FIND_COMPONENT<DTYPE,NTYPE> find_component(dimension);
    IJK::CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube(dimension);

    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); ientry++) {
      create_isodual_cube_table_entry
        (ientry, flag_separate_neg, flag_separate_opposite, cube,
         find_component, this->entry[ientry]);
    }
  }


  // Set incident isosurface polytopes for each table entry.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_CUBE_TABLE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetIncidentIsoPoly()
  {
    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); ientry++) {
      this->entry[ientry].SetIncidentIsoPoly(this->NumPolyVertices());
    }
  }
  
  
  // **************************************************
  // ISODUAL CUBE TABLE AMBIG MEMBER FUNCTIONS
  // **************************************************

  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_neg,
         const bool flag_separate_opposite)
  {
    this->SetToCube(dimension);
    
    TI_TYPE n = calculate_num_entries<TI_TYPE>(this->NumPolyVertices(), 2);
    this->SetNumTableEntries(n);

    CreateTableEntries(flag_separate_neg, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_opposite)
  {
    Create(dimension, true, flag_separate_opposite);
  }


  // Create table.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename DTYPE2>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension)
  {
    Create(dimension, true);
  }


  // Create table entries.
  // @param flag_separate_neg  If true, separate negative vertices.
  // @param flag_separate_opposite If true, always separate two diagonally 
  //        opposite negative or positive vertices 
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CreateTableEntries
  (const bool flag_separate_neg, const bool flag_separate_opposite)
  {
    this->flag_separate_neg = flag_separate_neg;
    this->flag_always_separate_opposite = flag_separate_opposite;

    const DTYPE dimension = this->Dimension();

    FIND_COMPONENT<DTYPE,NTYPE> find_component(dimension);
    ISODUAL_CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube(dimension);

    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); ientry++) {
      create_isodual_cube_table_ambig_entry
        (ientry, flag_separate_neg, flag_separate_opposite, cube,
         find_component, this->entry[ientry]);      
    }

  }


  // Return number of isosurface polytopes incident on isov,
  //   and containing edge dual to facet kf.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename IT_TYPE2, typename ITYPE3, typename ITYPE4,
	    typename CUBE_TYPE>
  NTYPE ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CountNumIncidentIsoPolyDualToFacetEdges
  (const IT_TYPE2 it, const ITYPE3 isov, const ITYPE4 kf,
   const CUBE_TYPE & cube) const
  {
    NTYPE num_incident_isopoly = 0;

    for (NTYPE ie = 0; ie < cube.NumEdges(); ie++) {
      if (!this->IsBipolar(it, ie)) { continue; }
      if (this->IncidentIsoVertex(it, ie) != isov) { continue; }

      const NTYPE edge_dir = cube.EdgeDir(ie);
      const NTYPE facet_orth_dir = cube.FacetOrthDir(kf);
      if (edge_dir == facet_orth_dir) {
	// Edge ke is not contained in facet.
	continue;
      }

      const NTYPE iv0 = cube.EdgeEndpoint(ie, 0);

      if (cube.IsFacetIncidentOnVertex(kf, iv0))
	{ num_incident_isopoly++; }
    }

    return num_incident_isopoly;
  }


  //@ Get orthogonal direction and side of (first) ambiguous facet.
  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename IT_TYPE2, typename DTYPE2,
            typename SIDE_TYPE>
  void ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  GetAmbiguousFacet
  (const IT_TYPE2 it, DTYPE2 & orth_dir, SIDE_TYPE & side) const
  {
    const DIMENSION_TYPE dimension = this->Dimension();
    const NUMBER_TYPE num_cube_facets = this->NumPolyFacets();

    const FACET_BITS_TYPE facet_set = this->AmbiguousFacetBits(it);
    const NUMBER_TYPE kf = IJK::get_first_one_bit(facet_set, num_cube_facets);
    orth_dir = IJK::cube_facet_orth_dir(dimension, kf);
    side = IJK::cube_facet_side(dimension, kf);
  }


  template <typename DTYPE, typename NTYPE, typename TI_TYPE,
            typename ENTRY_TYPE>
  template <typename TI_TYPE2, typename VTYPE>
  NTYPE ISODUAL_CUBE_TABLE_AMBIG<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  ComputeNumIncidentBipolarEdges
  (const TI_TYPE2 it, const VTYPE iv) const
  {
    const bool flag_pos = this->IsPositive(it, iv);

    NTYPE num_bipolar_edges = 0;

    for (int d = 0; d < this->Dimension(); d++) {
      VTYPE jv;
      IJK::compute_cube_vertex_neighbor(iv, d, jv);
      if (this->IsPositive(it, jv) != flag_pos)
        { num_bipolar_edges++; }
    }

    return num_bipolar_edges;
  }


  // **************************************************
  // IVOLDUAL CUBE TABLE BASE MEMBER FUNCTIONS
  // **************************************************

  // Constructor.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  IVOLDUAL_CUBE_TABLE_BASE(const DTYPE2 dimension) :
    ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension),
    cube(dimension)
  {
    Init();
  }


  // Constructor.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  IVOLDUAL_CUBE_TABLE_BASE
  (const DTYPE2 dimension, const bool flag_separate_neg) :
    ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>(dimension),
    cube(dimension)
  {
    Init();
  }

  // Initialize.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  void IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  Init()
  {
    IJK::PROCEDURE_ERROR error("IVOLDUAL_CUBE_TABLE_BASE::Init");

    if (!CheckNumVertexTypes(error)) { throw error; }
  }


  // True, if some neighbor of grid vertex is below interval volume.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename TI_TYPE2, typename VTYPE>
  bool IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  IsSomeNeighborBelowIntervalVolume
  (const TI_TYPE2 ientry, const VTYPE iv) const
  {
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      const VTYPE iv2 = cube.VertexNeighbor(iv, d);
      if (IsBelowIntervalVolume(ientry, iv2)) 
        { return(true); }
    }
    return(false);
  }


  // True, if some neighbor of grid vertex is above interval volume.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename TI_TYPE2, typename VTYPE>
  bool IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  IsSomeNeighborAboveIntervalVolume
  (const TI_TYPE2 ientry, const VTYPE iv) const
  {
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      const VTYPE iv2 = cube.VertexNeighbor(iv, d);
      if (IsAboveIntervalVolume(ientry, iv2)) 
        { return(true); }
    }
    return(false);
  }


  // Return number of neighbors of grid vertex which are below interval volume.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename TI_TYPE2, typename VTYPE>
  NTYPE IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  NumNeighborsBelowIntervalVolume
  (const TI_TYPE2 ientry, const VTYPE iv) const
  {
    NTYPE kount = 0;
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      const VTYPE iv2 = cube.VertexNeighbor(iv, d);
      if (IsBelowIntervalVolume(ientry, iv2)) 
        { kount++; }
    }
    return(kount);
  }


  // Return number of neighbors of grid vertex which are above interval volume.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename TI_TYPE2, typename VTYPE>
  NTYPE IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  NumNeighborsAboveIntervalVolume
  (const TI_TYPE2 ientry, const VTYPE iv) const
  {
    NTYPE kount = 0;
    for (DTYPE d = 0; d < this->Dimension(); d++) {
      const VTYPE iv2 = cube.VertexNeighbor(iv, d);
      if (IsAboveIntervalVolume(ientry, iv2)) 
        { kount++; }
    }
    return(kount);
  }


  // Return true if configuration contains a vertex with label vtype.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <const int vtype>
  bool IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  ContainsVertexType(const TI_TYPE it) const
  {
    for (NTYPE iv = 0; iv < this->NumPolyVertices(); iv++) {

      if (VertexType(it,iv) == vtype) 
        { return(true); }
    }

    return(false);
  }

  // Redefine SetDimension to set cube dimension.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  void IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  SetDimension(const DTYPE2 d)
  {
    cube.SetDimension(d);
    ISODUAL_TABLE_AMBIG_BASE<DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::SetDimension(d);
  }

  // Check number of vertex types.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  bool IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  CheckNumVertexTypes(IJK::ERROR & error)
  {
    if (NUM_VERTEX_TYPES != 4) {
      error.AddMessage
        ("Programming error.  Illegal number of vertex types.");
      error.AddMessage
        ("  Number of vertex types: ", NUM_VERTEX_TYPES, "");
      error.AddMessage
        ("  Number of vertex types should be 4.");
      return(false);
    }

    return(true);
  }

  // Check that cube dimension matches table dimension.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  bool IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES, DTYPE,NTYPE,
                                TI_TYPE,ENTRY_TYPE>::
  CheckCubeDimension(IJK::ERROR & error)
  {
    if (cube.Dimension() != this->Dimension()) {
      error.AddMessage
        ("Programming error.  Cube dimension does not match table dimension.");
      error.AddMessage("  Cube dimension: ", cube.Dimension(), ".");
      error.AddMessage("  Table dimension: ", this->Dimension(), ".");
      return(false);
    }

    return(true);
  }


  // **************************************************
  // IVOLDUAL CUBE TABLE MEMBER FUNCTIONS
  // **************************************************

  // Constructor.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  IVOLDUAL_CUBE_TABLE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  IVOLDUAL_CUBE_TABLE(const DTYPE2 dimension) :
    IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES,DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
    (dimension) 
  {
    Init();
    Create(dimension); 
  }

  // Constructor.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  IVOLDUAL_CUBE_TABLE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  IVOLDUAL_CUBE_TABLE
  (const DTYPE2 dimension, const bool flag_separate_neg) :
    IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES,DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
    (dimension)
  {
    Init();
    Create(dimension, flag_separate_neg); 
  }


  // Initialize.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  void IVOLDUAL_CUBE_TABLE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Init()
  {
    IJK::PROCEDURE_ERROR error("IVOLDUAL_CUBE_TABLE::Init");

    if (!this->CheckNumVertexTypes(error)) { throw error; }
  }


  // Create table.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  void IVOLDUAL_CUBE_TABLE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_neg)
  {
    this->SetToCube(dimension);
    TI_TYPE n = calculate_num_entries<TI_TYPE>
      (this->NumPolyVertices(), this->NumVertexTypes());
    this->SetNumTableEntries(n);
    CreateTableEntries(flag_separate_neg);
  }

  // Create table.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  void IVOLDUAL_CUBE_TABLE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension)
  {
    Create(dimension, true);
  }


  // Create table entries.
  // @param flag_separate_neg  If true, separate negative vertices.
  // @param flag_separate_opposite If true, always separate two diagonally 
  //        opposite negative or positive vertices 
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  void IVOLDUAL_CUBE_TABLE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CreateTableEntries(const bool flag_separate_neg)
  {
    typedef typename ENTRY_TYPE::IVOL_VERTEX_INDEX_TYPE IVOLV_TYPE;
    typedef typename ENTRY_TYPE::IVOL_POLY_INDEX_TYPE IVOLP_TYPE;
    typedef typename ENTRY_TYPE::FACET_BITS_TYPE FACET_BITS_TYPE;

    IJK::PROCEDURE_ERROR error("IVOLDUAL_CUBE_TABLE::CreateTableEntries");

    // Make sure that number of vertex types is correct (4).
    if (!this->CheckNumVertexTypes(error)) { throw error; }

    if (!this->CheckCubeDimension(error)) { throw error; }

    this->flag_separate_neg = flag_separate_neg;

    const DTYPE dimension = this->Dimension();

    FIND_COMPONENT<DTYPE,NTYPE> find_component(dimension+1);
    ISODUAL_CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> lifted_cube(dimension+1);

    ISODUAL_AMBIG_TABLE_ENTRY<NTYPE,IVOLV_TYPE,
                              IVOLP_TYPE,FACET_BITS_TYPE> 
      lower_isodual, upper_isodual;

    lower_isodual.Allocate
      (lifted_cube.NumVertices(), lifted_cube.NumEdges());
    upper_isodual.Allocate
      (lifted_cube.NumVertices(), lifted_cube.NumEdges());

    for (TI_TYPE ientry = 0; ientry < this->NumTableEntries(); 
         ientry++) {
      // Create ivoldual cube table entry
      create_ivoldual_cube_table_entry
        (ientry, this->NumTableEntries(), 
         flag_separate_neg, this->cube, lifted_cube,
         find_component, lower_isodual, upper_isodual, this->entry[ientry]);
    }

  }


  // **************************************************
  // IVOLDUAL CUBE DOUBLE TABLE MEMBER FUNCTIONS
  // **************************************************

  // Constructor.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  IVOLDUAL_CUBE_DOUBLE_TABLE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  IVOLDUAL_CUBE_DOUBLE_TABLE(const DTYPE2 dimension) :
    IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES,DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
    (dimension) 
  {
    Init();
    Create(dimension); 
  }


  // Constructor.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  IVOLDUAL_CUBE_DOUBLE_TABLE<NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  IVOLDUAL_CUBE_DOUBLE_TABLE
  (const DTYPE2 dimension, const bool flag_separate_neg) :
    IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES,DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>
    (dimension, flag_separate_neg)
  {
    Init();
    Create(dimension, flag_separate_neg);
  }


  // Initialize.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  void IVOLDUAL_CUBE_DOUBLE_TABLE
  <NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Init()
  {
    IJK::PROCEDURE_ERROR error("IVOLDUAL_CUBE_DOUBLE_TABLE::Init");

    if (!this->CheckNumVertexTypes(error)) { throw error; }
  }

  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename NTYPE2>
  void IVOLDUAL_CUBE_DOUBLE_TABLE
  <NUM_VERTEX_TYPES,DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  SetNumTableEntries(const NTYPE2 num_table_entries)
  {
    num_negative_table_entries = num_table_entries/2;
    IVOLDUAL_CUBE_TABLE_BASE<NUM_VERTEX_TYPES,DTYPE,NTYPE,TI_TYPE, ENTRY_TYPE>::
      SetNumTableEntries(num_table_entries);
  }

  // Create table.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  void IVOLDUAL_CUBE_DOUBLE_TABLE
  <NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension, const bool flag_separate_neg)
  {
    this->SetToCube(dimension);
    TI_TYPE n = calculate_num_entries<TI_TYPE>
      (this->NumPolyVertices(), this->NumVertexTypes());

    // Double number of table entries for double table.
    n = 2*n;

    this->SetNumTableEntries(n);
    CreateTableEntries(flag_separate_neg);
  }

  // Create table.
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  template <typename DTYPE2>
  void IVOLDUAL_CUBE_DOUBLE_TABLE
  <NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  Create(const DTYPE2 dimension)
  {
    Create(dimension, true);
  }


  // Create table entries.
  // @param flag_separate_neg  If true, separate negative vertices.
  // @param flag_separate_opposite If true, always separate two diagonally 
  //        opposite negative or positive vertices 
  template <const int NUM_VERTEX_TYPES, typename DTYPE, typename NTYPE, 
            typename TI_TYPE, typename ENTRY_TYPE>
  void IVOLDUAL_CUBE_DOUBLE_TABLE
  <NUM_VERTEX_TYPES, DTYPE,NTYPE,TI_TYPE,ENTRY_TYPE>::
  CreateTableEntries(const bool flag_separate_neg)
  {
    typedef typename ENTRY_TYPE::IVOL_VERTEX_INDEX_TYPE IVOLV_TYPE;
    typedef typename ENTRY_TYPE::IVOL_POLY_INDEX_TYPE IVOLP_TYPE;
    typedef typename ENTRY_TYPE::FACET_BITS_TYPE FACET_BITS_TYPE;

    IJK::PROCEDURE_ERROR error
      ("IVOLDUAL_CUBE_DOUBLE_TABLE::CreateTableEntries");

    // Make sure that number of vertex types is correct (4).
    if (!this->CheckNumVertexTypes(error)) { throw error; }

    this->flag_separate_neg = flag_separate_neg;

    const DTYPE dimension = this->Dimension();

    FIND_COMPONENT<DTYPE,NTYPE> find_component(dimension+1);
    ISODUAL_CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> cube(dimension);
    ISODUAL_CUBE_FACE_INFO<NTYPE,NTYPE,NTYPE> lifted_cube(dimension+1);

    ISODUAL_AMBIG_TABLE_ENTRY
      <NTYPE,IVOLV_TYPE,IVOLP_TYPE,FACET_BITS_TYPE> 
      lower_isodual, upper_isodual;

    lower_isodual.Allocate
      (lifted_cube.NumVertices(), lifted_cube.NumEdges());
    upper_isodual.Allocate
      (lifted_cube.NumVertices(), lifted_cube.NumEdges());

    // Create first half of table.
    for (TI_TYPE ientry = 0; ientry < this->NumNegativeTableEntries(); 
         ientry++) {
      // Create ivoldual cube table entry
      create_ivoldual_cube_table_entry
        (ientry, this->NumNegativeTableEntries(), 
         flag_separate_neg, cube, lifted_cube,
         find_component, lower_isodual, upper_isodual, this->entry[ientry]);
    }

    // Create second half of table.
    for (TI_TYPE ientry = 0; ientry < this->NumNegativeTableEntries(); 
         ientry++) {
      // Actual index of entry in lookup table.
      const NTYPE jentry = ientry+NumNegativeTableEntries();

      // Create ivoldual cube table entry
      create_ivoldual_cube_table_entry
        (ientry, this->NumNegativeTableEntries(), 
         !flag_separate_neg, cube, lifted_cube,
         find_component, lower_isodual, upper_isodual, this->entry[jentry]);
    }

  }


  // **************************************************
  // CLASS FIND_COMPONENT MEMBER FUNCTIONS
  // **************************************************

  template <typename DTYPE, typename NTYPE>
  template <typename DTYPE2>
  FIND_COMPONENT<DTYPE,NTYPE>::FIND_COMPONENT(const DTYPE2 dimension)
  {
    this->dimension = dimension;
    num_cube_vertices = NTYPE(1) << Dimension();

    vertex_flag = new bool[num_cube_vertices];
    component = new NTYPE[num_cube_vertices];
  }


  template <typename DTYPE, typename NTYPE>
  FIND_COMPONENT<DTYPE,NTYPE>::~FIND_COMPONENT()
  {
    if (component != NULL) { delete [] component; }
    component = NULL;

    if (vertex_flag != NULL) { delete [] vertex_flag; }
    vertex_flag = NULL;
  }


  template <typename DTYPE, typename NTYPE>
  void FIND_COMPONENT<DTYPE,NTYPE>::ClearAll()
  {
    for (NTYPE i = 0; i < NumCubeVertices(); i++) {
      vertex_flag[i] = false;
      component[i] = 0;
    }
  }


  template <typename DTYPE, typename NTYPE>
  template <typename TI_TYPE>
  void FIND_COMPONENT<DTYPE,NTYPE>::SetVertexFlags(const TI_TYPE ival)
  {
    convert2bool(ival, vertex_flag, num_cube_vertices);
  }


  template <typename DTYPE, typename NTYPE>
  void FIND_COMPONENT<DTYPE,NTYPE>::NegateVertexFlags()
  {
    for (NTYPE i = 0; i < num_cube_vertices; i++) 
      { vertex_flag[i] = (!vertex_flag[i]); }
  }


  template <typename DTYPE, typename NTYPE>
  template <typename ITYPE0, typename ITYPE1>
  void FIND_COMPONENT<DTYPE,NTYPE>::
  Search(const ITYPE0 i, const ITYPE1 icomp)
  {
    std::vector<NTYPE> vlist;

    vlist.reserve(num_cube_vertices);

    if (icomp == 0) {
      IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
      error.AddMessage("Programming error. Component number cannot be zero.");
      throw error;
    }

    bool flag = vertex_flag[i];
    vlist.push_back(i);
    component[i] = icomp;
    while(vlist.size() > 0) {
      NTYPE j = vlist.back();
      vlist.pop_back();
      for (int d = 0; d < dimension; d++) {
        int j2;
        IJK::compute_cube_vertex_neighbor(j, d, j2);
        if (vertex_flag[j2] == flag && component[j2] == 0) {
          vlist.push_back(j2);
          component[j2] = icomp;
        }
      }
    }
  }


  template <typename DTYPE, typename NTYPE>
  template <typename ITYPE0, typename ITYPE1, typename ITYPE2>
  void FIND_COMPONENT<DTYPE,NTYPE>::
  SearchFacet(const ITYPE0 kf, const ITYPE1 i, const ITYPE2 icomp)
  {
    std::vector<NTYPE> vlist;

    vlist.reserve(num_cube_vertices);

    if (icomp == 0) {
      IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
      error.AddMessage("Programming error. Component number cannot be zero.");
      throw error;
    }

    if (!IJK::cube_facet_contains(dimension, kf, i)) {
      IJK::PROCEDURE_ERROR error("FIND_COMPONENT::Search");
      error.AddMessage("Programming error. Facet ", kf,
                       " does not contain vertex ", i, ".");
      throw error;
    }

    bool flag = vertex_flag[i];
    vlist.push_back(i);
    component[i] = icomp;
    while(vlist.size() > 0) {
      NTYPE j = vlist.back();
      vlist.pop_back();
      for (int d = 0; d < dimension; d++) {
        int j2;
        IJK::compute_cube_vertex_neighbor(j, d, j2);
        if (vertex_flag[j2] == flag && component[j2] == 0 &&
            IJK::cube_facet_contains(dimension, kf, j2)) {
          vlist.push_back(j2);
          component[j2] = icomp;
        }
      }
    }
  }


  // Compute number of components of vertices.
  template <typename DTYPE, typename NTYPE>
  template <typename TI_TYPE>
  NTYPE FIND_COMPONENT<DTYPE,NTYPE>::ComputeNumComponents
  (const TI_TYPE ientry, const bool flag_positive)
  {
    ClearAll();
    SetVertexFlags(ientry);
    if (!flag_positive) { NegateVertexFlags(); }

    NTYPE num_components(0);
    for (NTYPE i = 0; i < NumCubeVertices(); i++) {
      if (Component(i) == 0) {

        if (VertexFlag(i)) {
          num_components++;
          Search(i, num_components);
        }
      }
    }

    return(num_components);
  }


  // Compute number of components of vertices in facet.
  template <typename DTYPE, typename NTYPE>
  template <typename TI_TYPE, typename ITYPE>
  NTYPE FIND_COMPONENT<DTYPE,NTYPE>::ComputeNumComponentsInFacet
  (const TI_TYPE ientry, const ITYPE kf, const bool flag_positive)
  {
    ClearAll();
    SetVertexFlags(ientry);
    if (!flag_positive) { NegateVertexFlags(); }

    NTYPE num_components(0);
    for (NTYPE i = 0; i < NumCubeVertices(); i++) {
      if (Component(i) == 0 &&
          IJK::cube_facet_contains(dimension, kf, i)) {

        if (VertexFlag(i)) {
          num_components++;
          SearchFacet(kf, i, num_components);
        }
      }
    }

    return(num_components);
  }


  // **************************************************
  // CLASS ISODUAL_CUBE_FACE_INFO MEMBER FUNCTIONS
  // **************************************************


  template <typename DTYPE, typename NTYPE, typename VTYPE>
  template <typename TI_TYPE>
  void ISODUAL_CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  ComputeNumCubeFacetBits(const TI_TYPE ientry, const NTYPE ifacet,
                          NTYPE & num_zeros, NTYPE & num_ones) const
  {
    num_ones = 0;
    num_zeros = 0;

    for (NTYPE j = 0; j < this->NumFacetVertices(); j++) {
      const NTYPE jv = this->FacetVertex(ifacet, j);
      // *** Replace by IsPositive(jv)
      const TI_TYPE mask = (TI_TYPE(1) << jv);
      if ((mask & ientry) == 0) {
        num_zeros++;
      }
      else {
        num_ones++;
      }
    }
  }


  template <typename DTYPE, typename NTYPE, typename VTYPE>
  template <typename TI_TYPE>
  bool ISODUAL_CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  IsCubeFacetActive(const TI_TYPE ientry, const NTYPE ifacet) const
  {
    NTYPE num_zeros, num_ones;

    ComputeNumCubeFacetBits(ientry, ifacet, num_zeros, num_ones);
    if (num_zeros > 0 && num_ones > 0) 
      { return(true); }
    else 
      { return(false); }
  }


  template <typename DTYPE, typename NTYPE, typename VTYPE>
  template <typename TI_TYPE>
  NTYPE ISODUAL_CUBE_FACE_INFO<DTYPE,NTYPE,VTYPE>::
  ComputeNumActiveCubeFacets(const TI_TYPE ientry) const
  {
    NTYPE num_active_facets = 0;
    for (NTYPE ifacet = 0; ifacet < this->NumFacets(); ifacet++) {
      if (IsCubeFacetActive(ientry, ifacet)) 
        { num_active_facets++; }
    }

    return(num_active_facets);
  }

}

#endif
