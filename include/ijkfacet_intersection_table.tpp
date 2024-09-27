/*!
 *  @file ijkfacet_intersection_table.tpp
 *  @brief Table for determining vertices on facet intersections.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2016-2023 Rephael Wenger

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


#ifndef _IJKFACET_INTERSECTION_TABLE_
#define _IJKFACET_INTERSECTION_TABLE_

#include "ijk.tpp"
#include "ijkcube.tpp"

namespace IJK {

  // *********************************************************************
  // TEMPLATE CLASS FACET_INTERSECTION_TABLE
  // *********************************************************************

  /// Class to report vertices in facet intersections.
  /// @tparam VBITS_TYPE Type to store vertex bits.
  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  class FACET_INTERSECTION_TABLE {

  protected:
    DTYPE dimension;
    NUM_TYPE num_cube_vertices;
    NUM_TYPE num_table_entries;
    VBITS_TYPE * in_ridge_bits;
    VBITS_TYPE * in_facet_intersection_bits;
    VBITS_TYPE * in_facet_bits;

    /// True if array entry[] is allocated.
    bool is_table_allocated;
    
    /// Initialization routine.
    template <typename DTYPE2>
    void Init(const DTYPE2 dimension);

    void Clear();                           ///< Clear routine.
    void ComputeInRidgeBits();              ///< Compute in_ridge_bits[].

    /// Compute in_facet_intersection_bits[].
    void ComputeInFacetIntersectionBits();

    /// Compute in_facet_bits[].
    void ComputeInFacetBits();


  public:
    typedef NUM_TYPE NUMBER_TYPE;
    typedef VBITS_TYPE VERTEX_BITS_TYPE;

  public:
    FACET_INTERSECTION_TABLE() { Init(3); };
    template <typename DTYPE2>
    FACET_INTERSECTION_TABLE(const DTYPE2 dimension)
    { Init(dimension); }
    ~FACET_INTERSECTION_TABLE() { Clear(); }


    // get functions
    DTYPE Dimension() const { return(dimension); }
    NUM_TYPE NumCubeVertices() const { return(num_cube_vertices); }
    NUM_TYPE NumTableEntries() const { return(num_table_entries); }
    NUM_TYPE NumCubeFacets() const { return 2*Dimension(); }
    DTYPE FacetNormalDir(const int jf) const
    { return int(jf/2); }

    /// @brief Return 0 if facet is incident on origin.
    /// - Return 1 if facet is not incident on origin.
    DTYPE FacetSide(const int jf) const
    { return (jf%2); }

    /// Return true if some facet in ientry is orthogonal to direction orth_dir.
    bool IsSomeFacetOrthogonalTo
    (const NUM_TYPE ientry, const int orth_dir) const;

    /*!
     *  @brief Return array of bits representing which cube vertices are in ridge.
     *  @param ientry Array of bits representing
     *    which cube facets are on grid boundary.
     *    - If bit 2d is true, then lower cube facet orthogonal
     *      to direction d is on grid boundary.
     *    - If bit 2d+1 is true, then upper cube facet orthogonal
     *      to direction d is on grid boundary.
     *    - Unfortunately, this facet indexing is different than
     *      the facet indexing for cube facets in ijkcube.tpp.
     */
    VBITS_TYPE InRidgeBits(const NUM_TYPE ientry) const
    { return in_ridge_bits[ientry]; }

    /*!
     *  @brief Return array of bits respresenting which cube vertices are in intersection of facets.
     *  - Bit iv is true if cube vertex iv is in all facets
     *    indicated by ientry.
     *  @param ientry Array of bits representing cube facets.
     *    - Bit 2d represents lower cube facet orthogonal
     *      to direction d.
     *    - Bit 2d+1 represents upper cube facet orthogonal
     *      to direction d.
     *    - Unfortunately, this facet indexing is different than
     *      the facet indexing for cube facets in ijkcube.tpp.
     */
    VBITS_TYPE InFacetIntersectionBits(const NUM_TYPE ientry) const
    { return in_facet_intersection_bits[ientry]; }

    /*!
     *  @brief Return array of bits representing which cube vertices are in facets.
     *  - Bit iv is true if cube vertex iv is in any facet
     *    indicated by ientry.
     *  @param ientry Array of bits representing cube facets.
     *    - Bit 2d represents lower cube facet orthogonal
     *      to direction d.
     *    - Bit 2d+1 represents upper cube facet orthogonal
     *      to direction d.
     *    - Unfortunately, this facet indexing is different than
     *      the facet indexing for cube facets in ijkcube.tpp.
     */
    VBITS_TYPE InFacetBits(const NUM_TYPE ientry) const
    { return in_facet_bits[ientry]; }

    /// Return true if table is allocated.
    bool IsTableAllocated() const { return(is_table_allocated); }

    template <typename DTYPE2>
    void SetDimension(const DTYPE2 dimension);

    void Create();

    // Check routines.
    template <typename DTYPE2>
    void CheckVertexBitsType(const DTYPE2 dimension) const;
    template <typename DTYPE2>
    void CheckNumType(const DTYPE2 dimension) const;

    /// Check arrays in_ridge_bits[], in_facet_intersection_bits, and in_facet_bits[].
    /// - Return true if no errors found (passes check).
    /// - Return false and set error if error found.
    bool CheckBits(IJK::ERROR & error) const;
  };


  // *********************************************************************
  // TEMPLATE CLASS FACET_INTERSECTION_TABLE MEMBER FUNCTIONS
  // *********************************************************************

  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  template  <typename DTYPE2>
  void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
  Init(const DTYPE2 dimension)
  {
    in_ridge_bits = NULL;
    in_facet_intersection_bits = NULL;
    in_facet_bits = NULL;
    is_table_allocated = false;
    num_table_entries = 0;
    SetDimension(dimension);
  }

  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  template  <typename DTYPE2>
  void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
  SetDimension(const DTYPE2 dimension)
  {
    Clear();
    this->dimension = dimension;
    num_cube_vertices = IJK::compute_num_cube_vertices(dimension);
    CheckVertexBitsType(dimension);
    CheckNumType(dimension);
  }


  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::Clear()
  {
    if (in_ridge_bits != NULL) { delete [] in_ridge_bits; }
    in_ridge_bits = NULL;
    if (in_facet_intersection_bits != NULL) 
      { delete [] in_facet_intersection_bits; }
    in_facet_intersection_bits = NULL;
    if (in_facet_bits != NULL)
      { delete [] in_facet_bits; }
    num_table_entries = 0;
    dimension = 0;
    num_cube_vertices = 1;
    is_table_allocated = false;
  }


  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
  ComputeInRidgeBits()
  {
    IJK::UNIT_CUBE<DTYPE,NUM_TYPE,NUM_TYPE> unit_cube(Dimension());

    for (NUM_TYPE ientry = 0; ientry < NumTableEntries(); ientry++) {
      in_ridge_bits[ientry] = 0;
      for (NUM_TYPE jv = 0; jv < NumCubeVertices(); jv++) {
        NUM_TYPE num_selected_facet_directions = 0;
        for (DTYPE orth_dir = 0; orth_dir < Dimension(); orth_dir++) {
          // Bit 2*orth_dir is 1, if lower/leftmost facet 
          //   with orth direction orth_dir is on grid boundary.
          NUM_TYPE lower_facet_mask = (NUM_TYPE(1) << (2*orth_dir));
          // Bit 2d is 1, if lower/leftmost facet 
          //   with orth direction orth_dir is on grid boundary.
          NUM_TYPE upper_facet_mask = (NUM_TYPE(1) << (2*orth_dir+1));
          NUM_TYPE c = unit_cube.VertexCoord(jv,orth_dir);
          if (c == 0) {
            if ((lower_facet_mask & ientry) != 0)
              { num_selected_facet_directions++; }
          }
          else {
            if ((upper_facet_mask & ientry) != 0)
              { num_selected_facet_directions++; }
          }
        }

        if (num_selected_facet_directions > 1) {
          // Vertex jv is on ridge.
          VBITS_TYPE vertex_mask = (VBITS_TYPE(1) << jv);
          in_ridge_bits[ientry] = (in_ridge_bits[ientry] | vertex_mask);
        }
      }
    }

  }


  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
  ComputeInFacetIntersectionBits()
  {
    IJK::UNIT_CUBE<DTYPE,NUM_TYPE,NUM_TYPE> unit_cube(Dimension());

    for (NUM_TYPE ientry = 0; ientry < NumTableEntries(); ientry++) {
      in_facet_intersection_bits[ientry] = 0;

      NUM_TYPE num_selected_facets = 0;
      NUM_TYPE ival = ientry;
      for (NUM_TYPE jf = 0; jf < NumCubeFacets(); jf++) {
        if ((ival %2) == 1) { num_selected_facets++; }
        ival = (ival >> 1);
      }

      if (num_selected_facets > 0) {

        for (NUM_TYPE jv = 0; jv < NumCubeVertices(); jv++) {
          NUM_TYPE num_facets_containing_vertex = 0;

          for (NUM_TYPE jf = 0; jf < NumCubeFacets(); jf++) {
            NUM_TYPE facet_mask = (NUM_TYPE(1) << jf);

            DTYPE orth_dir = NUM_TYPE(jf/2);
            if ((facet_mask & ientry) != 0) {
              NUM_TYPE c = unit_cube.VertexCoord(jv,orth_dir);

              if (c == 0 && (jf%2 == 0))
                { num_facets_containing_vertex++; }
              else if (c == 1 && (jf%2 == 1))
                { num_facets_containing_vertex++; }
            }
          }

          if (num_selected_facets == num_facets_containing_vertex) {
            // Vertex jv is in all selected facets.
            VBITS_TYPE vertex_mask = (VBITS_TYPE(1) << jv);
            in_facet_intersection_bits[ientry] = 
              (in_facet_intersection_bits[ientry] | vertex_mask);
          }
        }
      }
    }

  }


  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
   void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
   ComputeInFacetBits()
  {
    IJK::UNIT_CUBE<DTYPE,NUM_TYPE,NUM_TYPE> unit_cube(Dimension());

    for (NUM_TYPE ientry = 0; ientry < NumTableEntries(); ientry++) {
      in_facet_bits[ientry] = 0;

      for (NUM_TYPE jv = 0; jv < NumCubeVertices(); jv++) {
        const VBITS_TYPE vertex_mask = (VBITS_TYPE(1) << jv);

        for (NUM_TYPE jf = 0; jf < NumCubeFacets(); jf++) {
          const NUM_TYPE facet_mask = (NUM_TYPE(1) << jf);

          const DTYPE facet_normal_dir = FacetNormalDir(jf);
          if ((facet_mask & ientry) != 0) {
            const NUM_TYPE c =
              unit_cube.VertexCoord(jv,facet_normal_dir);

            if (c == (jf%2)) {
              in_facet_bits[ientry] =
                (in_facet_bits[ientry] | vertex_mask);
              break;
            }
          }
        }
      }
    }
  }


  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::Create()
  {
    IJK::PROCEDURE_ERROR error("FACET_INTERSECTION_TABLE::Create");

    if (is_table_allocated || in_ridge_bits != NULL ||
        in_facet_intersection_bits != NULL) {
      error.AddMessage
        ("Programming error.  Facet intersection table already created.");
      throw error;
    }

    CheckVertexBitsType(Dimension());
    CheckNumType(Dimension());

    NUM_TYPE num_cube_facets = compute_num_cube_facets(Dimension());

    num_table_entries = (NUM_TYPE(1) << num_cube_facets);
    in_ridge_bits = new VBITS_TYPE[num_table_entries];
    in_facet_intersection_bits = new VBITS_TYPE[num_table_entries];
    in_facet_bits = new VBITS_TYPE[num_table_entries];
    is_table_allocated = true;

    ComputeInRidgeBits();
    ComputeInFacetIntersectionBits();
    ComputeInFacetBits();
  }


  /// Return true if some facet in ientry is orthogonal to direction orth_dir.
  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  bool FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
  IsSomeFacetOrthogonalTo(const NUM_TYPE ientry, const int orth_dir) const
  {
    const NUM_TYPE num_cube_facets = compute_num_cube_facets(Dimension());

    for (NUM_TYPE jfacet = 0; jfacet < num_cube_facets; jfacet++) {
      const NUMBER_TYPE facet_mask = (NUMBER_TYPE(1) << jfacet);
      if ((facet_mask & ientry) != 0) {
        // Facet is in ientry.

        const int facet_normal_dir = FacetNormalDir(jfacet);
        if (facet_normal_dir == orth_dir)
          { return true; }
      }
    }

    return false;
  }


  /// Check VBITS_TYPE.
  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  template  <typename DTYPE2>
  void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
  CheckVertexBitsType(const DTYPE2 dimension) const
  {
    NUM_TYPE num_vert = IJK::compute_num_cube_vertices(dimension);

    if (!check_number_of_bits<VBITS_TYPE>(num_vert)) {
      IJK::ERROR error;

      error.AddMessage
        ("Programming error using template class FACET_INTERSECTION_TABLE.");
      error.AddMessage
        ("  Fewer than ", num_vert, 
         " bits in type template parameter VBITS_TYPE.");
      error.AddMessage
        ("  Replace type template parameter with type with more bits.");
      throw error;
    }
  }


  /// Check NUM_TYPE.
  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  template  <typename DTYPE2>
  void FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
  CheckNumType(const DTYPE2 dimension) const
  {
    const NUM_TYPE num_facets = IJK::compute_num_cube_facets(dimension);

    if (!check_number_of_bits<VBITS_TYPE>(num_facets)) {
      IJK::ERROR error;

      error.AddMessage
        ("Programming error using template class FACET_INTERSECTION_TABLE.");
      error.AddMessage
        ("  Fewer than ", num_facets,
         " bits in type template parameter NUM_TYPE.");
      error.AddMessage
        ("  Replace type template parameter with type with more bits.");

      throw error;
    }
  }
  

  // Check arrays in_ridge_bits[], in_facet_intersection_bits, and in_facet_bits[].
  // - Return true if no errors found (passes check).
  // - Return false and set error if error found.
  template <typename DTYPE, typename VBITS_TYPE, typename NUM_TYPE>
  bool FACET_INTERSECTION_TABLE<DTYPE,VBITS_TYPE,NUM_TYPE>::
  CheckBits(IJK::ERROR & error) const
  {
    const int TWO(2);
    IJK::UNIT_CUBE<DTYPE,NUM_TYPE,NUM_TYPE> unit_cube(Dimension());

    for (NUM_TYPE ientry = 1; ientry < NumTableEntries(); ientry++) {
      for (NUM_TYPE iv = 0; iv < NumCubeVertices(); iv++) {
        // Vertex mask
        const VBITS_TYPE vmask = (VBITS_TYPE(1) << iv);

        bool flag_in_all_facets = true;
        NUM_TYPE num_facets_containing_vertex = 0;
        for (NUM_TYPE jfacet = 0; jfacet < NumCubeFacets(); jfacet++) {

          const NUM_TYPE facet_mask = (NUM_TYPE(1) << jfacet);

          // Facet direction and side.
          const DTYPE facet_normal_dir = FacetNormalDir(jfacet);
          const DTYPE facet_side = FacetSide(jfacet);

          const NUM_TYPE c =
            unit_cube.VertexCoord(iv,facet_normal_dir);

          if ((facet_mask & ientry) != 0) {
            // Facet jfacet is selected.

            if (c == FacetSide(jfacet)) {
              // Vertex iv is in facet jfacet.

              num_facets_containing_vertex++;
              if ((vmask & in_facet_bits[ientry]) == 0) {
                // Error. Vertex iv is in some facet. Bit should be one.
                error.AddMessage
                  ("Programming error. Incorrect in_facet_bits[",
                   ientry, "].");
                error.AddMessage
                  (" Vertex ", iv, " is in facet ", jfacet, ".");
                error.AddMessage
                  ("  Bit for vertex ", iv, " should be 1.");
                error.AddMessage
                  ("  in_facet_bits[", ientry, "]: ",
                   in_facet_bits[ientry], "");
                return false;
              }
            }
            else {
              // Vertex iv is not in facet jfacet.

              flag_in_all_facets = false;
              if ((vmask & in_facet_intersection_bits[ientry]) != 0) {
                // Error. If vertex iv is not in facet,
                //   then vertex iv is not in facet intersection.
                // Bit should be zero.
                error.AddMessage
                  ("Programming error. Incorrect in_facet_intersection_bits[",
                   ientry, "].");
                error.AddMessage
                  ("  Vertex ", iv, " is not in facet ", jfacet, ".");
                error.AddMessage
                  ("  Bit for vertex ", iv, " should be 0.");
                error.AddMessage
                   ("  in_facet_intersection_bits[", ientry, "]: ",
                    in_facet_intersection_bits[ientry], "");

                return false;
              }
            }
          }
        }

        if (flag_in_all_facets && (num_facets_containing_vertex > 0)) {
          if ((vmask & in_facet_intersection_bits[ientry]) == 0) {
            // Error. Vertex iv is in all facets. Bit should be one.
            error.AddMessage
              ("Programming error. Incorrect in_facet_intersection_bits[",
               ientry, "].");
            error.AddMessage("  Vertex ", iv, " is in all facets.");
            error.AddMessage
              ("  Bit for vertex ", iv, " should be 1.");
            return false;
          }
        }

        if (num_facets_containing_vertex == 0) {
          if ((vmask & in_facet_bits[ientry]) != 0) {
            // Error. Vertex iv is in not in any facet. Bit should be zero.
            error.AddMessage
              ("Programming error. Incorrect in_facet_bits[",
               ientry, "].");
            error.AddMessage("  Vertex ", iv, " is not in any facet.");
            error.AddMessage
              ("  Bit for vertex ", iv, " should be 0.");
            return false;
          }
        }

        if (num_facets_containing_vertex >= TWO) {
          if ((vmask & in_ridge_bits[ientry]) == 0) {
            // Error. Vertex iv is in ridge. Bit should be one.
            error.AddMessage
              ("Programming error. Incorrect in_ridge_bits[",
               ientry, "].");
            error.AddMessage
              ("  Vertex ", iv, " is in exactly two facets.");
            error.AddMessage
              ("  Bit for vertex ", iv, " should be 1.");
            return false;
          }
        }
        else {
          if ((vmask & in_ridge_bits[ientry]) != 0) {
            // Error. Vertex iv is in ridge. Bit should be one.
            error.AddMessage
              ("Programming error. Incorrect in_ridge_bits[",
               ientry, "].");
            error.AddMessage
              ("  Vertex ", iv, " is in ",
               num_facets_containing_vertex, " facets.");
            error.AddMessage
              ("  Bit for vertex ", iv, " should be 0.");
            error.AddMessage
               ("  in_ridge_bits[", ientry, "]: ",
                in_ridge_bits[ientry], "");
            return false;
          }
        }
      }
    }

    return true;
  }

}

#endif

