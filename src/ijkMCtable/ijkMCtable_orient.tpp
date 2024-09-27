/*!
 *  @file ijkMCtable_orient.tpp
 *  @brief Template functions for orienting Marching Cubes lookup table.
 *  - Version 0.5.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2023-2024 Rephael Wenger

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

#ifndef _IJKTABLE_ORIENT_TPP_
#define _IJKTABLE_ORIENT_TPP_

#include <algorithm>
#include <bitset>
#include <stack>
#include <vector>

#include "ijk.tpp"
#include "ijksimplex.tpp"
#include "ijkMCtable.h"


/// Classes and routines for orienting Marching Cubes lookup table.
namespace IJKMCUBE_TABLE {

  /*!
   *  @brief Size of bitset used to represent isosurface vertices
   *    in isotable entries.
   *  - Note: Number of isosurface vertices can be 
   *    at most ISO_VERTEX_BITSET_SIZE.
   */
  const int ISO_VERTEX_BITSET_SIZE = 64;

  typedef std::bitset<ISO_VERTEX_BITSET_SIZE> ISO_VERTEX_BITSET;

  /*!
   *  @brief Size of bitset used to represent connected components
   *    in isotable entries.
   *  - Note: Number of isosurface vertices can be 
   *    at most ISO_CONNECTED_COMPONENT_BITSET_SIZE.
   */
  const int ISO_CONNECTED_COMPONENT_BITSET_SIZE = 32;

  typedef std::bitset<ISO_CONNECTED_COMPONENT_BITSET_SIZE>
  ISO_CONNECTED_COMPONENT_BITSET;

  /*!
   *  @brief Simplex orientation information.
   */
  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  class SIMPLEX_ORIENT_INFO_BASE
  {
    
  public:
    
    typedef _ISO_VERTEX_BITSET_TYPE ISO_VERTEX_BITSET_TYPE;

  public:
    
    /// @brief Indicates isosurface vertices in simplex;
    /// - in_simplex[iw] = 1, if isosurface vertex iw
    ///   is in simplex js.
    _ISO_VERTEX_BITSET_TYPE in_simplex;

    /// @brief Indicates boundary facets.
    /// - is_boundary_facet[iw] = 1, if facet formed by removing
    ///   vertex iw is a boundary facet.
    /// - Undefined if iw is not a vertex of the simplex.
    _ISO_VERTEX_BITSET_TYPE is_boundary_facet;

    /// @brief Indicates facet swap parity.
    /// - facet_swap_parity[iw] is parity of swaps that move iw
    ///   to last vertex and sort remaining vertices
    ///   in increasing order.
    /// - Undefined if iw is not a vertex of the simplex.
    _ISO_VERTEX_BITSET_TYPE facet_swap_parity;

    /// @brief Index of connected component containing simplex.
    _ITYPE index_of_connected_component;
    
  public:

    /// Constructor.
    SIMPLEX_ORIENT_INFO_BASE() {}

    // Set functions.

    /// @brief Set index of connected component containing simplex.
    void SetConnectedComponent(const _ITYPE icomponent)
    { index_of_connected_component = icomponent; }

    
    // Get functions.

    /// @brief Return true if iw is a simplex vertex.
    bool InSimplex(const int iw) const
    { return bool(in_simplex[iw]); }

    /// @brief Return true if formed by removing
    ///   vertex iw is a boundary facet.
    bool IsBoundaryFacet(const int iw) const
    { return bool(is_boundary_facet[iw]); }
    
    /// @brief Return facet swap parity (0 or 1) of vertex iw.
    /// - facet_swap_parity[iw] is parity of swaps that move iw
    ///   to last vertex and sort remaining vertices
    ///   in increasing order.
    int FacetSwapParity(const int iw) const
    { return facet_swap_parity[iw]; }

    /// @brief Return index of connected component.
    int IndexOfConnectedComponent() const
    { return index_of_connected_component; }

    /// @brief Return size of isosurface vertex bitset.
    /// - Note: Not all bits are used.
    int IsoVertexBitsetSize() const
    { return ISO_VERTEX_BITSET_TYPE::size(); }

    /// @brief Return bitset in_simplex.
    ISO_VERTEX_BITSET_TYPE InSimplex() const
    { return in_simplex; }

    /// @brief Return bitset is_boundary_facet.
    ISO_VERTEX_BITSET_TYPE IsBoundaryFacet() const
    { return is_boundary_facet; }

    /// @brief Return bitset facet_swap_parity.
    ISO_VERTEX_BITSET_TYPE FacetSwapParity() const
    { return facet_swap_parity; }

    /// @brief Return true if simplex has a boundary facet.
    bool HasBoundaryFacet() const
    { return is_boundary_facet.any(); }

  };

  
  /*!
   *  @brief Flag isosurface vertices in isotable polytope facets.
   */
  template <typename _ISO_VERTEX_BITSET_TYPE>
  class FACET_ISO_VERTEX {

  public:
    typedef _ISO_VERTEX_BITSET_TYPE ISO_VERTEX_BITSET_TYPE;
    
  protected:

    /// @brief Number of isotable polytope facets.
    int num_facets;

    /// @brief Number of isosurface vertices.
    int num_isosurface_vertices;

    /// @brief in_facet[jf][iw] = true, if isosurface vertex iw
    ///   is in facet jf.
    std::vector<_ISO_VERTEX_BITSET_TYPE> in_facet;

  protected:

    /// @brief Initialize.
    void Init(const ISOSURFACE_TABLE &  isotable);

  public:

    /// @brief Constructor.
    FACET_ISO_VERTEX(const ISOSURFACE_TABLE & isotable)
    { Init(isotable); }

    /// @brief Return bitset indicating isosurface vertices
    ///    in facet ifacet.
    _ISO_VERTEX_BITSET_TYPE InFacet(const int ifacet) const
    { return in_facet[ifacet]; }
                      
    /// @brief Return number of facets.
    int NumFacets() const
    { return num_facets; }

    /// @brief Return number of isosurface vertices.
    int NumIsosurfaceVertices() const
    { return num_isosurface_vertices; }

    /*!
     *  @brief Return true if isosurface vertices indicated by isov_bitset
     *    are contained in facet ifacet.
     */
    bool AreVerticesInFacet
    (const ISO_VERTEX_BITSET_TYPE & isov_bitset, const int ifacet) const
    {
      const ISO_VERTEX_BITSET_TYPE shared_vert =
        isov_bitset & InFacet(ifacet);

      const ISO_VERTEX_BITSET_TYPE not_in_facet =
        isov_bitset ^ shared_vert;

      return not_in_facet.none();
    }
    
    
    // *** Return string representation. Mainly for debugging. ***

    /// @brief Return string representing bitset in_facet[ifacet].
    /// - Only returns num_isosurface_vertices bits,
    ///   since all other bits are zero.
    std::string InFacetStr(const int ifacet) const;
  };


  /// @brief Table containing orientation information.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  class MCUBE_ISOTABLE_ORIENT_INFO {

  public:
    typedef _ISO_VERTEX_BITSET_TYPE ISO_VERTEX_BITSET_TYPE;
    typedef _ISO_CONNECTED_COMPONENT_BITSET_TYPE
    ISO_CONNECTED_COMPONENT_BITSET_TYPE;
    typedef SIMPLEX_ORIENT_INFO_BASE<_ISO_VERTEX_BITSET_TYPE,_ITYPE>
    SIMPLEX_ORIENT_INFO_TYPE;
    typedef FACET_ISO_VERTEX<_ISO_VERTEX_BITSET_TYPE>
    FACET_INFO_TYPE;
    
  protected:

    class ORIENT_INFO_ENTRY {

    public:
      
      /// @brief simplex_info[js] Orientation information
      ///   for simplex js.
      std::vector<SIMPLEX_ORIENT_INFO_TYPE> simplex_info;

      /// @brief Number of connected components.
      _ITYPE num_connected_components;

      /// @brief num_oriented_connected_components Number of oriented
      ///   connected components.
      _ITYPE num_oriented_connected_components;

      /// @brief is_connected_component_oriented[icomponent]
      ///    True if connected component icomponent is oriented.
      ISO_CONNECTED_COMPONENT_BITSET_TYPE
      is_connected_component_oriented;

    public:

      /// Constructor.      
      ORIENT_INFO_ENTRY()
      {
        num_connected_components = 0;
        num_oriented_connected_components = 0;
      }
        
      /// @brief Set index of connected component containing isimplex.
      void SetConnectedComponent
      (const int isimplex, const _ITYPE icomponent)
      { simplex_info[isimplex].SetConnectedComponent(icomponent); }

      /// @brief Return number of simplices in entry.
      int NumSimplices() const
      { return simplex_info.size(); }

      /// @brief Return number of oriented connected components.
      int NumOrientedConnectedComponents() const
      { return num_oriented_connected_components; }
      
      /// @brief Return index of connected component containing 
      ///   simplex isimplex.
      int IndexOfConnectedComponent(const int isimplex) const
      { return simplex_info[isimplex].IndexOfConnectedComponent(); }
      
      /// @brief Return bitset is_connected_component_oriented.
      ISO_CONNECTED_COMPONENT_BITSET_TYPE
      IsConnectedComponentOriented() const
      { return is_connected_component_oriented; }

      /// @brief Return is_connected_component_oriented[icomponent]
      bool IsConnectedComponentOriented(const int icomponent) const
      { return bool(is_connected_component_oriented[icomponent]); }

      /// @brief Return true if simplex isimplex is oriented.
      bool IsSimplexOriented(const int isimplex) const
      { return IsConnectedComponentOriented
          (IndexOfConnectedComponent(isimplex)); }

      /// @brief Return true if all simplices are oriented.
      /// - Equivalently, return true if all connected components are oriented.
      bool AreAllSimplicesOriented() const
      { return (num_connected_components ==
                num_oriented_connected_components); }

    };
    
    
  protected:

    /// @brief Number of vertices per simplex.
    int num_vertices_per_simplex;

    /// @brief Number of isosurface vertices.
    int num_isosurface_vertices;

    /// @brief Table entries.
    std::vector<ORIENT_INFO_ENTRY> entry;

    /// @brief Facet information.
    FACET_INFO_TYPE facet_info;
      
    /// @brief Initialize.
    void Init(const ISOSURFACE_TABLE & isotable);

    /// @brief Set connected component of simplex isimplex
    ///   in entry[table_index].
    void _SetConnectedComponent
    (const TABLE_INDEX table_index, const int isimplex,
     const _ITYPE icomponent)
    {
      entry[table_index].SetConnectedComponent(isimplex,icomponent);
    }

    /*!
     *  @brief Flag vertices in each simplex.
     */
    void _FlagVerticesInEachSimplex
    (const TABLE_INDEX table_index, 
     const ISOSURFACE_VERTEX_INDEX * simplex_vertex_list,
     const int num_simplices);

    /*!
     *  @brief Set facet swap parity in all facets of simplex isimplex
     *    in table entry table_index.
     *  @param table_index Table index.
     *  @param isimplex Simplex index.
     *  @param simplex_vertex[j] j'th vertex of simplex isimplex.
     *  @param temp_simplex_vert[] Temporary array to store simplex vertices.
     *    @pre temp_simplex_vert[] has size at least NumVertPerSimplex().
     */
    void _SetSimplexFacetSwapParity
    (const TABLE_INDEX table_index,
     const int isimplex,
     const ISOSURFACE_VERTEX_INDEX * simplex_vertex,
     ISOSURFACE_VERTEX_INDEX temp_simplex_vert[]);

    /*!
     *  @overload
     *  @brief Set facet swap parity in all facets of simplex isimplex
     *    in table entry table_index. (No temporary array.)
     *  - Version whose argument list does not include 
     *    temporary array temp_simplex_vert[].
     */
    void _SetSimplexFacetSwapParity
    (const TABLE_INDEX table_index,
     const int isimplex,
     const ISOSURFACE_VERTEX_INDEX * simplex_vertex);

    /*!
     *  @brief Set parity of swaps that move iw to last vertex
     *    and sort remaining vertices in increasing order.
     *  - Set boundary facets, i.e. facets thar are
     *    in only one simplex in simplex_vertex_list[].
     *  @pre Call _FlagVerticesInEachSimplex before
     *    calling this function.
     *  @param simplex_vertex_list[] List of simplices.
     *    - Each simplex has NumVerticesPerSimplex() vertices.
     */
    void _SetFacetSwapParity
    (const TABLE_INDEX table_index, 
     const ISOSURFACE_VERTEX_INDEX * simplex_vertex_list,
     const int num_simplices);

    /*!
     *  @brief Flag boundary facets, i.e. facets thar are
     *    in only one simplex in simplex_vertex_list[].
     *  @pre Call _FlagVerticesInEachSimplex before
     *    calling this function.
     *  @param simplex_vertex_list[] List of simplices.
     *    - Each simplex has NumVerticesPerSimplex() vertices.
     */
    void _FlagBoundaryFacets(const TABLE_INDEX table_index);

    /// @brief Set connected component containing each simplex
    ///    for entry[table_index].
    void _SetConnectedComponentIndices
    (const TABLE_INDEX table_index,
     const ISOSURFACE_VERTEX_INDEX * simplex_vertex_list,
     const int num_simplices);

    /// @brief Convert ISO_VERTEX_BITSET_TYPE to string.
    /// - Only returns num_isosurface_vertices bits,
    ///   since all other bits are zero.
    std::string _ConvertBitsetToString
    (const ISO_VERTEX_BITSET_TYPE & bitset) const;
    
    
  public:

    /// @brief Constructor.
    MCUBE_ISOTABLE_ORIENT_INFO(const ISOSURFACE_TABLE & isotable):
      facet_info(isotable)
    { Init(isotable); }

    
    // *** Get routines ***

    /// @brief Return number of table entries.
    TABLE_INDEX NumTableEntries() const
      { return entry.size(); }
    
    /// @brief Return number of vertices per simplex.
    int NumVerticesPerSimplex() const
    { return num_vertices_per_simplex; }

    /// @brief Return number of isosurface vertices.
    int NumIsosurfaceVertices() const
    { return num_isosurface_vertices; }

    /// Return number of simplices for entry table_index.
    int NumSimplices(const TABLE_INDEX table_index) const
    { return entry[table_index].NumSimplices(); }

    /// @brief Return orientation information for simplex isimplex
    ///   in table entry table_index.
    const SIMPLEX_ORIENT_INFO_TYPE & SimplexInfo
    (const TABLE_INDEX table_index, const int isimplex) const
    { return entry[table_index].simplex_info[isimplex]; }

    /// @brief Return facet information.
    const FACET_INFO_TYPE & FacetInfo() const
    { return facet_info; }

    /// @brief Return index of connected component containing simplex isimplex
    ///   in table entry table_index.
    int IndexOfConnectedComponent
    (const TABLE_INDEX table_index, const int isimplex) const
    { return SimplexInfo(table_index, isimplex).IndexOfConnectedComponent(); }
    
    /// @brief Return true if connected component icomponent is oriented.
    bool IsConnectedComponentOriented
    (const TABLE_INDEX table_index,
     const int icomponent) const
    { return entry[table_index].IsConnectedComponentOriented[icomponent]; };

    /// @brief Return true if simplex isimplex in table entry table_index
    ///   is oriented.
    bool IsSimplexOriented
    (const TABLE_INDEX table_index,
     const int isimplex) const
    { return entry[table_index].IsSimplexOriented(isimplex); }

    /// @brief Return flag are_all_simplices_oriented
    ///   for table entry table_index.
    bool AreAllSimplicesOriented(const TABLE_INDEX table_index) const
    { return entry[table_index].AreAllSimplicesOriented(); }

    /// @brief Return true if all table entries are marked as oriented.
    /// @param[out] table_index Index of table entry that has not
    ///    been completely oriented.
    bool AreAllTableEntriesOriented(TABLE_INDEX & table_index) const;

    /// @brief Return number of connected components
    ///   in table entry table_index.
    int NumConnectedComponents
    (const TABLE_INDEX table_index) const
    { return entry[table_index].num_connected_components; }

    /// @brief Return number of oriented connected components
    ///   in table entry table_index.
    int NumOrientedConnectedComponents
    (const TABLE_INDEX table_index) const
    { return entry[table_index].NumOrientedConnectedComponents(); }
    
    /// @brief Return number of facets per simplex.
    int NumFacetsPerSimplex() const
    { return NumVerticesPerSimplex(); }

    /// @brief Return number of vertices per simplex facet.
    int NumVerticesPerSimplexFacet() const
    { return NumVerticesPerSimplex()-1; }

    /// @brief Return true if some simplex in connected component
    ///    icomponent has a boundary facet.
    bool ComponentHasBoundaryFacet
    (const TABLE_INDEX table_index, const int icomponent) const;

    /// @brief Return true if some simplex in table entry has 
    ///    a boundary facet.
    bool SomeSimplexHasBoundaryFacet(const TABLE_INDEX table_index) const;

    /*!
     *  @brief Find index of table entry with a single
     *    connected component that has at least one
     *    boundary simplex facet.
     *  @param[out] table_index Index of table entry found.
     *    - Undefined if table entry not found.
     *  @param[out] flag_found True if table entry found.
     */
    template <typename _TABLE_INDEX>
    void FindEntryWithSingleConnectedComponentWithBoundary
    (_TABLE_INDEX & table_index, bool & flag_found) const;


    // *** Orient routines ***

    /*!
     *  @brief Set flag is_connected_component_oriented 
     *    for connected component icomponent.
     *  - Also modifies entry[table_index].num_oriented_connected_components.
     */
    void SetIsConnectedComponentOriented
    (const TABLE_INDEX table_index, const int icomponent,
     const bool flag);

    /*!
     *  @brief Update flag entry[table_index].is_connected_component_oriented. 
     *  - Update by "or" with is connected_component_oriented.
     *  - Also modifies entry[table_index].num_oriented_connected_components.
     */
    void UpdateIsConnectedComponentOriented
    (const TABLE_INDEX table_index, 
     const ISO_CONNECTED_COMPONENT_BITSET & is_component_oriented);

    /*!
     *  @brief Return true if two simplices share a facet.
     *  - Also returns bitset of shared vertices, number of shared vertices,
     *    comparison of facet swap parities, and bitsets indicating 
     *    unshared vertices.
     *  @param[out] shared_vert Bits indicating vertices shared 
     *    by two simplices.
     *  @param[out] num_shared_vertices Number of vertices shared 
     *    by two simplices.
     *  @param[out] are_shared_facet_swap_parities_equal
     *    - True if simplices share a facet and
     *      shared facet has the same swap parity in each simplex.
     *    - False if simplices share a facet and
     *      shared facet has different swap parities in each simplex.
     *    - Undefined if simplices do not share a facet.
     *  @param[out] not_in_facetA Bit indicating vertex of simplexA
     *    not shared with simplexB.
     *    - Undefined if simplices A and B do not share a facet
     *      or are identical.
     *  @param[out] not_in_facetB Bit indicating vertex of simplexB
     *    not shared with simplexA.
     *    - Undefined if simplices A and B do not share a facet
     *      or are identical.
     */
    bool DoSimplicesShareFacet
    (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     const int isimplexA, const int isimplexB,
     ISO_VERTEX_BITSET_TYPE & shared_vert,
     int & num_shared_vertices,
     bool & are_shared_facet_swap_parities_equal,
     ISO_VERTEX_BITSET_TYPE & not_in_facetA,
     ISO_VERTEX_BITSET_TYPE & not_in_facetB) const;

    /*!
     *  @overload
     *  @brief Return true if two simplices share a facet.
     *    (Return only comparison of swap parities.)
     *  - Version that returns only comparison of swap parities.
     */
    bool DoSimplicesShareFacet
    (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     const int isimplexA, const int isimplexB,
     bool & are_shared_facet_swap_parities_equal) const;
    
    /*!
     *  @brief Return true if two simplices share a boundary facet.
     *  - Also returns bitset of shared vertices, number of shared vertices,
     *    comparison of facet swap parities, and bitsets indicating 
     *    unshared vertices.
     *  @param[out] shared_vert Bits indicating vertices shared 
     *    by two simplices.
     *  @param[out] num_shared_vertices Number of vertices shared 
     *    by two simplices.
     *  @param[out] are_shared_facet_swap_parities_equal
     *    - True if simplices share a facet and
     *      shared facet has the same swap parity in each simplex.
     *    - False if simplices share a facet and
     *      shared facet has different swap parities in each simplex.
     *    - Undefined if simplices do not share a facet.
     *  @param[out] not_in_facetA Bit indicating vertex of simplexA
     *    not shared with simplexB.
     *    - Undefined if simplices A and B do not share a facet
     *      or are identical.
     *  @param[out] not_in_facetB Bit indicating vertex of simplexB
     *    not shared with simplexA.
     *    - Undefined if simplices A and B do not share a facet
     *      or are identical.
     */
    bool DoSimplicesShareBoundaryFacet
    (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     const int isimplexA, const int isimplexB,
     ISO_VERTEX_BITSET_TYPE & shared_vert,
     int & num_shared_vertices,
     bool & are_shared_facet_swap_parities_equal,
     ISO_VERTEX_BITSET_TYPE & not_in_facetA,
     ISO_VERTEX_BITSET_TYPE & not_in_facetB) const;


    /*!
     *  @overload
     *  @brief Return true if two simplices share a boundary facet.
     *  - Version that does not return bitsets indicating 
     *    unshared vertices.
     */
    bool DoSimplicesShareBoundaryFacet
    (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     const int isimplexA, const int isimplexB,
     ISO_VERTEX_BITSET_TYPE & shared_vert,
     int & num_shared_vertices,
     bool & are_shared_facet_swap_parities_equal) const;

    /*!
     *  @brief Orient all simplices in connected component containing isimplexA
     *    consistently with isimplexA.
     *  - Modifies both isotable and entry[table_index].simplex_info[].
     */
    void OrientConnectedComponent
    (const TABLE_INDEX table_index, const int isimplexA,
     ISOSURFACE_TABLE & isotable);
    
    /*!
     *  @brief Consistently orient all simplices in each
     *    connected component in table entry table_index.
     *  - Orientation of each connected component is arbitrary.
     *  - Modifies both isotable and entry[table_index].simplex_info[].
     */
    void OrientAllSimplicesInTableEntry
    (const TABLE_INDEX table_index, ISOSURFACE_TABLE & isotable);


    /*! 
     *  @brief Flip all simplices in connected component icomponent
     *    in table entry table_index.
     *  - Also recomputes facet_swap_parity for each modified simplex.
     */
    void FlipSimplicesInConnectedComponent
    (const TABLE_INDEX table_index, const int icomponent,
     ISOSURFACE_TABLE & isotable);
    
    /*!
     *  @brief Orient simplices in table entry table_indexB
     *    from orientation of some boundary facet of table_indexA.
     *    - Does nothing if some simplex in table_indexB has
     *      orientation matching orientation of boundary facet 
     *      from table_indexA.
     *  @pre OrientAllSimplicesInTableEntry() has been called
     *    on table_indexA and table_indexB.
     *  @param table_indexA Index of first table entry.
     *  @param table_indexB Index of second table entry.
     *    - Orient table_indexB from table_indexA.
     *    @pre table_indexB != table_indexA.
     */
    void OrientTwoTableEntries
    (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     ISOSURFACE_TABLE & isotable);

    
    // *** Check routines ***

    /// @brief Check against isosurface lookup table.
    /// - Return true if passes check.
    template <typename _ISOTABLE_TYPE>
    bool Check
    (const _ISOTABLE_TYPE & isotable, IJK::ERROR & error) const;
    
    /// @brief Check that bitset entry[table_index].simplex_info[*].in_simplex 
    ///   is set.
    /// - Return true if all bitsets are set.
    bool CheckInSimplexIsSet
    (const TABLE_INDEX table_index, IJK::ERROR & error) const;

    /*!
     *  @brief Check that orientations of all pairs of simplices 
     *    in entry[table_index] are consistent.
     *  - Return true if all orientations in entry[table_index] are consistent.
     *  @param[out] isimplexA Simplex in entry[table_index]
     *    whose orientation does not match isimplexB.
     *    - Defined only if check returns false.
     *  @param[out] isimplexB Simplex in entry[table_index]
     *    whose orientation does not match isimplexA.
     *    - Defined only if check returns false.
     *  - Note: isimplexA and isimplexB share a facet.
     */
    bool CheckOrientationsInTableEntry
    (const TABLE_INDEX table_index, int & isimplexA, int & isimplexB,
    IJK::ERROR & error) const;

    /*!
     *  @brief Check that orientations of all pairs of simplices 
     *    in entry[table_index] are consistent. (Doesn't return simplices.)
     *  - Version that does not return mismatched simplices.
     */
    bool CheckOrientationsInTableEntry
    (const TABLE_INDEX table_index, IJK::ERROR & error) const;    
    
    /*!
     *  @brief Check that orientations in every table entry
     *     are consistent.
     *  - Does NOT compare orientations of simplices in
     *    DIFFERENT table entries.
     *  @param[out] table_index Index of table entry
     *    with simplices with mismatched orientations.
     *    - Defined only if check returns false.
     *  @param[out] isimplexA Simplex in entry[table_index]
     *    whose orientation does not match isimplexB.
     *    - Defined only if check returns false.
     *  @param[out] isimplexB Simplex in entry[table_index]
     *    whose orientation does not match isimplexA.
     *    - Defined only if check returns false.
     *  - Note: isimplexA and isimplexB share a facet.

     */
    bool CheckOrientationsInEveryTableEntry
    (TABLE_INDEX & table_index, int & isimplexA, int & isimplexB,
    IJK::ERROR & error) const;

    /*!
     *  @brief Check that orientations in every table entry
     *     are consistent. (Doesn't return table index or simplices.)
     *  - Does NOT compare orientations of simplices in
     *    DIFFERENT table entries.
     *  - Version that does not return index of table entry or 
     *    mismatched simplices.
     */
    bool CheckOrientationsInEveryTableEntry(IJK::ERROR & error) const;


    /*!
     *  @brief  Return true if orientations of simplices 
     *    in entry[table_indexB] are consistent with orientations 
     *    of simplices in entry[table_indexA].
     *  @pre CheckOrientationsInTableEntry(table_indexA) returns true
     *    and CheckOrientationsInTableEntry(table_indexB) returns true.
     *  @param[out] isimplexA Simplex in entry[table_indexA]
     *    whose orientation does not match isimplexB.
     *    - Defined only if check returns false.
     *  @param[out] isimplexB Simplex in entry[table_indexB]
     *    whose orientation does not match isimplexB.
     *    - Defined only if check returns false.
     *    - Note: isimplexA and isimplexB share a facet.
     *  @param[out] componentB_checked Bitset indicating
     *    which components were checked.
     *    - componentB_checked[i] = true if some boundary facet
     *      of some simplex in component i of table_indexB
     *      matched a boundary facet of some simplex of table_indexA.
     */
    template <typename _ISOTABLE_TYPE>
    bool CheckOrientationsOfTwoTableEntries
    (const _ISOTABLE_TYPE & isotable,
     const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     int & isimplexA, int & isimplexB,
     ISO_CONNECTED_COMPONENT_BITSET_TYPE & componentB_checked,
     IJK::ERROR & error) const;

    /// @brief Return true if orientations of simplices in entry[table_indexB] 
    ///   are consistent with orientations of simplices in entry[table_indexA].
    /// - Version that does not return bitset componentB_checked.
    template <typename _ISOTABLE_TYPE>
    bool CheckOrientationsOfTwoTableEntries
    (const _ISOTABLE_TYPE & isotable,
     const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     int & isimplexA, int & isimplexB,
     IJK::ERROR & error) const;

    /*!
     *  @brief Check orientation of table entry table_indexA against
     *    all other table entries.
     *  @param[out] isimplexA Simplex in entry[table_indexA]
     *    whose orientation does not match isimplexB.
     *    - Defined only if check returns false.
     *  @param[out] table_indexB Table index whose orientation
     *      does not match orientation in table_indexA.
     *    - Defined only if check returns false.
     *  @param[out] isimplexB Simplex in entry[table_indexB]
     *    whose orientation does not match isimplexB.
     *    - Defined only if check returns false.
     *    - Note: isimplexA and isimplexB share a facet.
     */
    template <typename _ISOTABLE_TYPE>
    bool CheckOrientationOfTableEntryAgainstAllOthers
    (const _ISOTABLE_TYPE & isotable,
     const TABLE_INDEX table_indexA, 
     int & isimplexA, TABLE_INDEX & table_indexB, int & isimplexB,
     IJK::ERROR & error) const;

    /*!
     *  @brief Check orientation of table entry table_indexA against
     *    all other table entries. (Doesn't return table index or simplices.)
     *  - Version that does not return table index or simplices
     *    that have inconsistent orientations.
     */
    template <typename _ISOTABLE_TYPE>
    bool CheckOrientationOfTableEntryAgainstAllOthers
    (const _ISOTABLE_TYPE & isotable,
     const TABLE_INDEX table_indexA, 
     IJK::ERROR & error) const;


    // *** Return string representation. Mainly for debugging. ***

    /// @brief Return string representing bitset SimplexInfo().in_simplex.
    /// - Only returns num_isosurface_vertices bits,
    ///   since all other bits are zero.
    std::string InSimplexStr
    (const TABLE_INDEX table_index, const int isimplex) const
    {
      return _ConvertBitsetToString
        (SimplexInfo(table_index,isimplex).in_simplex);
    }

    /// @brief Return string representing bitset SimplexInfo().is_boudary_facet.
    /// - Only returns num_isosurface_vertices bits,
    ///   since all other bits are zero.
    std::string IsBoundaryFacetStr
    (const TABLE_INDEX table_index, const int isimplex) const
    {
      return _ConvertBitsetToString
        (SimplexInfo(table_index,isimplex).is_boundary_facet);
    }

    /// @brief Return string representing bitset SimplexInfo().facet_swap_parity.
    /// - Only returns num_isosurface_vertices bits,
    ///   since all other bits are zero.
    std::string FacetSwapParityStr
    (const TABLE_INDEX table_index, const int isimplex) const
    {
      return _ConvertBitsetToString
        (SimplexInfo(table_index,isimplex).facet_swap_parity);
    }
    
  };


  namespace {

    /*!
     *  @brief Return true passes check.
     *  - Table entry istart must have at least one simplex and
     *    exactly one connected component.
     */
    template <typename _ISOTABLE_TYPE, typename _ORIENT_INFO_TABLE>
    bool check_orient_simplices_starting_table_entry
    (const _ISOTABLE_TYPE & isotable,
     const _ORIENT_INFO_TABLE & orient_info,
     const TABLE_INDEX istart, IJK::ERROR & error)
    {
      if (isotable.NumSimplices(istart) < 1) {
        error.AddMessage
        ("Programming error. Marching Cubes lookup table entry ",
         istart, " has no simplices.");
        error.AddMessage
          ("  Table index must have at least one simplex to orient MC table.");
        return false;
        
      }

      const int numc = orient_info.NumConnectedComponents(istart);

      if (numc < 1) {
        error.AddMessage
          ("Programming error. Marching Cubes table entry ", istart,
           " has zero connected component.");
        error.AddMessage
          ("  Routine must start from table entry with exactly one connected component.");
        return false;
      }
      else if (numc != 1) {
        error.AddMessage
          ("Programming error. Marching Cubes table entry ", istart,
           " has ", numc, " connected components.");
        error.AddMessage
          ("  Routine must start from table entry with exactly one connected component.");
        return false;
      }

      return true;
    }

    
    /// @brief Swap a[i] with last element in a[].
    /// - Does nothing if a[i] is last element.
    /// @pre i < a.size().
    template <typename ETYPE>
    void _swap_with_last(const int i, std::vector<ETYPE> & a)
    {
      if (i+1 < a.size()) {
        const int ilast = a.size() - 1;
        std::swap(a[i], a[ilast]);
      }
    }


    /*!
     *  @brief Orient simplex lists in Marching Cubes lookup table. (Local.)
     *  @param isotable Marching Cubes lookup table.
     *  @param istart Starting table index.
     *    - Orient consistent with first simplex in isotable.entry[istart].
     *    @pre isotable.NumSimplices(table_index) > 0.
     *    @pre Simplices in isotable.entry[table_index] are
     *    in one connected component.
     *  @param flag_verbose Print messages if flag_verbose is true.
     *  @param output_trigger Print messages whenever more than output_trigger
     *    table entries are completed.
     */
    template <typename _OSTREAM_TYPE,
              typename _ISOTABLE_TYPE, typename _ITYPE,
              typename _NTYPE>
    void orient_mcube_table_local
    (_OSTREAM_TYPE & out,
     _ISOTABLE_TYPE & isotable, const _ITYPE istart,
     const bool flag_verbose, const _NTYPE output_trigger)
    {
      typedef MCUBE_ISOTABLE_ORIENT_INFO
        <ISO_VERTEX_BITSET, ISO_CONNECTED_COMPONENT_BITSET, unsigned char>
        ORIENT_INFO_TABLE_TYPE;
    
      ORIENT_INFO_TABLE_TYPE orient_info(isotable);
      IJK::PROCEDURE_ERROR error("orient_mcube_table_local");

      if (orient_info.NumVerticesPerSimplex() < 2) {
        // Nothing to orient.
        return;
      }

      if (!isotable.CheckTableIndex(istart, error))
        { throw error; }

      for (TABLE_INDEX table_index = 0;
           table_index < isotable.NumTableEntries(); table_index++) {
        orient_info.OrientAllSimplicesInTableEntry(table_index, isotable);
      }

      if (!check_orient_simplices_starting_table_entry
          (isotable, orient_info, istart, error))
        { throw error; }


      int num_completed = 1;

      // Mark bounded components as oriented.
      for (TABLE_INDEX table_index = 0; 
           table_index < orient_info.NumTableEntries(); table_index++) {
        for (int icomponent = 0; 
             icomponent < orient_info.NumConnectedComponents(table_index);
             icomponent++) {

          if (!orient_info.ComponentHasBoundaryFacet
              (table_index, icomponent)) {
            // Isosurface patch is closed?!? with no boundary.
            // Component was already oriented.
            // Nothing to check against other table entries.
            orient_info.SetIsConnectedComponentOriented
              (table_index, icomponent, true);

            if (orient_info.AreAllSimplicesOriented(table_index))
              { num_completed++; }
          }
        }
      }

      // Array of unoriented table entries.
      std::vector<TABLE_INDEX> unoriented_entry;
      
      for (TABLE_INDEX table_index = 0;
           table_index < orient_info.NumTableEntries();
           table_index++) {

        if (table_index == istart) {
          // Orient all other table entries from istart.
          continue;
        }
        
        if (orient_info.AreAllSimplicesOriented(table_index)) {
          // Either no simplices or no boundary simplex facets.
          continue;
        }

        unoriented_entry.push_back(table_index);
      }


      // Stack containing indices of table entries with 1 component that have been oriented.
      std::stack<TABLE_INDEX> stackI;

      // Stack containing indices of table entries with multiple components
      //   that have been oriented.
      std::stack<TABLE_INDEX> stack_multi;
      
      const int icomponent =
        orient_info.IndexOfConnectedComponent(istart, 0);
      orient_info.SetIsConnectedComponentOriented
        (istart, icomponent, true);

      stackI.push(istart);
      while (stackI.size() > 0) {
        const int table_indexA = stackI.top();
        stackI.pop();

        int j = 0;
        while (j < unoriented_entry.size()) {
          const TABLE_INDEX table_indexB = unoriented_entry[j];

          if ((orient_info.AreAllSimplicesOriented(table_indexB)) ||
              (table_indexA == table_indexB) ||
              (isotable.NumSimplices(table_indexB) == 0)) {
            // None of these should happen, but just in case.

            // Remove table_indexB from unoriented_entries.
            _swap_with_last(j, unoriented_entry);
            unoriented_entry.pop_back();

            // Process new element in unoriented_entry[j].
            continue;
          }

          orient_info.OrientTwoTableEntries
            (table_indexA, table_indexB, isotable);


          if (orient_info.AreAllSimplicesOriented(table_indexB)) {
            if (orient_info.NumConnectedComponents(table_indexB) == 1)
              { stackI.push(table_indexB); }
            else
              { stack_multi.push(table_indexB); }

            // Remove table_indexB from unoriented_entries.
            _swap_with_last(j, unoriented_entry);
            unoriented_entry.pop_back();

            num_completed++;
            if (flag_verbose && num_completed%output_trigger == 0) {
              out << "  Completed orientation of " << num_completed
                << " isosurface table entries.\n";
              out.flush();
            }

            // Process new element in unoriented_entry[j].
            continue;
          }

          // Go to next element in unoriented_entry[j].
          j = j+1;
        }
      }

      if (unoriented_entry.size() != 0) {
        // Try using stack of multi components to orient.
        while (stack_multi.size() > 0) {
          const int table_indexA = stack_multi.top();
          stack_multi.pop();

          int j = 0;
          while (j < unoriented_entry.size()) {
            const TABLE_INDEX table_indexB = unoriented_entry[j];

            if ((orient_info.AreAllSimplicesOriented(table_indexB)) ||
                (table_indexA == table_indexB) ||
                (isotable.NumSimplices(table_indexB) == 0)) {
              // None of these should happen, but just in case.

              // Remove table_indexB from unoriented_entries.
              _swap_with_last(j, unoriented_entry);
              unoriented_entry.pop_back();

              // Process new element in unoriented_entry[j].
              continue;
            };

            orient_info.OrientTwoTableEntries
              (table_indexA, table_indexB, isotable);

            if (orient_info.AreAllSimplicesOriented(table_indexB)) {
              { stack_multi.push(table_indexB); }

              // Remove table_indexB from unoriented_entries.
              _swap_with_last(j, unoriented_entry);
              unoriented_entry.pop_back();

              num_completed++;
              if (flag_verbose && num_completed%output_trigger == 0) {
                out << "  Completed orientation of " << num_completed
                    << " isosurface table entries.\n";
                out.flush();
              }

              // Process new element in unoriented_entry[j].
              continue;
            }

            // Go to next element in unoriented_entry[j].
            j = j+1;
          }
        }
      }


      if (flag_verbose) {
        TABLE_INDEX table_indexQ;
        if (orient_info.AreAllTableEntriesOriented(table_indexQ)) {
          if (num_completed > output_trigger) {
            out << "  Completed orientation of all isosurface table entries.\n";
            out.flush();
          }
        }
        else {
          out << "*** Warning: Unable to determine orientation for table index: "
              << table_indexQ << ".\n";
          out.flush();
        }
      }
    }


    
    // *****************************************************************
    // NO_OUTPUT_MESSAGE class
    // *****************************************************************

    /// @brief Class that does not output any message.
    class NO_OUTPUT_MESSAGE {
      
    public:
      /* Empty class */

      void flush()
      { /* Do nothing. */ }

    };


    /// @brief Overload << on NO_OUTPUT_MESSAGE.
    template <typename _T>
    NO_OUTPUT_MESSAGE & operator <<
    (NO_OUTPUT_MESSAGE & no_output_message, const _T x)
    {
      // Do nothing.
      return no_output_message;
    }

  }


  // *****************************************************************
  // Orient mcube table routines
  // *****************************************************************


  /*!
   *  @brief Orient simplex lists in Marching Cubes lookup table.
   *  @param isotable Marching Cubes lookup table.
   *  @param istart Starting table index.
   *    - Orient consistent with first simplex in isotable.entry[istart].
   *    @pre isotable.NumSimplices(table_index) > 0.
   *    @pre Simplices in isotable.entry[table_index] are
   *    in one connected component.
   *  @param flag_verbose Print messages if flag_verbose is true.
   *  @param output_trigger Print messages whenever more than output_trigger
   *    table entries are completed.
   */
  template <typename OSTREAM_TYPE,
            typename _ISOTABLE_TYPE, typename _ITYPE,
            typename _NTYPE>
  void orient_mcube_table
  (OSTREAM_TYPE & out,
   _ISOTABLE_TYPE & isotable, const _ITYPE istart,
   const bool flag_verbose, const _NTYPE output_trigger)
  {
    orient_mcube_table_local
      (out, isotable, istart, flag_verbose, output_trigger);
  }


  /*!
   *  @overload
   *  @brief Orient simplex lists in Marching Cubes lookup table. 
   *    (No output messages.)
   *  - Version that does not print output messages.
   *  @param isotable Marching Cubes lookup table.
   *  @param istart Starting table index.
   *    - Orient consistent with first simplex in isotable.entry[istart].
   *    @pre isotable.NumSimplices(table_index) > 0.
   *    @pre Simplices in isotable.entry[table_index] are
   *    in one connected component.
   */
  template <typename _ISOTABLE_TYPE, typename _ITYPE>
  void orient_mcube_table
  (_ISOTABLE_TYPE & isotable, const _ITYPE istart)
  {
    NO_OUTPUT_MESSAGE no_output_message;

    // Dummy parameters.
    bool flag_verbose = false;
    int output_trigger(1);

    orient_mcube_table_local
      (no_output_message, isotable, istart, flag_verbose, output_trigger);
  }


  // *****************************************************************
  // Check mcube table orientation
  // *****************************************************************

  /*!
   *  @brief Return true if all simplex lists in Marching Cubes lookup table
   *    are consistently oriented. 
   *  - Fast algorithm: Skips entries that have already matched
   *    some table entry.
   *  @param isotable Marching Cubes lookup table.
   */
  template <typename OSTREAM_TYPE, typename _ISOTABLE_TYPE,
	    typename _NTYPE>
  bool check_mcube_table_orientation
  (OSTREAM_TYPE & out, const _ISOTABLE_TYPE & isotable,
   const bool flag_verbose, const _NTYPE output_trigger,
   IJK::ERROR & error)
  {
    typedef typename _ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef MCUBE_ISOTABLE_ORIENT_INFO
      <ISO_VERTEX_BITSET, ISO_CONNECTED_COMPONENT_BITSET, unsigned char>
      ORIENT_INFO_TABLE_TYPE;

    const int ONE(1);
    const int num_table_entries = isotable.NumTableEntries();
    const int num_vert_per_simplex = isotable.NumVerticesPerSimplex();
    int isimplexA, isimplexB;

    int num_table_entries_checked = 0;

    if (num_vert_per_simplex < 2) {
      // Nothing to check
      num_table_entries_checked = num_table_entries;
      return true;
    }

    ORIENT_INFO_TABLE_TYPE orient_info(isotable);

    if (!orient_info.CheckOrientationsInEveryTableEntry(error))
      { return false; }

    if (flag_verbose && 
        orient_info.NumTableEntries() > output_trigger) {
      out << "  All table entries have internal orientation consistency.\n";
    }

    bool flag_found;
    TABLE_INDEX istart;
    orient_info.FindEntryWithSingleConnectedComponentWithBoundary
      (istart, flag_found);

    if (!flag_found) {
      // No table entries with one connected component
      //   where connected component has a boundary.
      if (flag_verbose) {
        out << "***  No table entries with single connected component\n"
            << "     where connected component has a boundary.\n";
      }
      return false;
    }

    // Stack containing indices of table entries that have been checked.
    std::stack<TABLE_INDEX> stack;

    num_table_entries_checked = 1;

    const int icomponent =
      orient_info.IndexOfConnectedComponent(istart, 0);
    orient_info.SetIsConnectedComponentOriented
      (istart, icomponent, true);

    stack.push(istart);
    while (stack.size() > 0) {
      const int table_indexA = stack.top();
      stack.pop();

      for (int table_indexB = 0;
	   table_indexB < isotable.NumTableEntries(); table_indexB++) {

	if (orient_info.AreAllSimplicesOriented(table_indexB)) {
	  // table_indexB has already been checked.
	  continue;
	}

	ISO_CONNECTED_COMPONENT_BITSET componentB_checked;
	
	if (!orient_info.CheckOrientationsOfTwoTableEntries
	    (isotable, table_indexA, table_indexB,
	     isimplexA, isimplexB, componentB_checked, error))
	  { return false; }

	orient_info.UpdateIsConnectedComponentOriented
	  (table_indexB, componentB_checked);

	if (orient_info.AreAllSimplicesOriented(table_indexB)) {

	  if (orient_info.NumConnectedComponents(table_indexB) == ONE) 
	    { stack.push(table_indexB); }

	  num_table_entries_checked++;


	  if (flag_verbose && (num_table_entries_checked > 0) &&
	      num_table_entries_checked%output_trigger == 0) {
	    out << "  Checked " << num_table_entries_checked
		<< " out of " << isotable.NumTableEntries()
		<< " isosurface table entry orientations.\n";
	    out.flush();
	  }
	}
      }
    }

    // Verify that all table entries have been checked.
    for (TABLE_INDEX table_index = 0;
	 table_index < orient_info.NumTableEntries();
	 table_index++) {

      if (orient_info.AreAllSimplicesOriented(table_index))
	{ continue; }

      if (!orient_info.SomeSimplexHasBoundaryFacet(table_index)) {
	// Isosurface patch is closed?!? with no boundary.
	// Consistency of table entry was already checked.
	// Nothing to check against other table entries.
	continue;
      }

      if (!orient_info.CheckOrientationOfTableEntryAgainstAllOthers
	  (isotable, table_index, error)) {
	return false;
      }

      num_table_entries_checked++;

      if (flag_verbose && (num_table_entries_checked > 0) &&
	  num_table_entries_checked%output_trigger == 0) {
	out << "  Checked " << num_table_entries_checked
	    << " out of " << isotable.NumTableEntries()
	    << " isosurface table entry orientations.\n";
	out.flush();
      }
    }

    if (flag_verbose &&
	isotable.NumTableEntries() > output_trigger) {
      out << "  Checked orientations on all "
	  << isotable.NumTableEntries()
	  << " isosurface table entries.\n";
      out.flush();
    }

    return true;
  }


  /*!
   *  @brief Return true if all simplex lists in Marching Cubes lookup table
   *    are consistently oriented.
   *  - Check every table entry against every other table entry. (Very slow).
   *  @param isotable Marching Cubes lookup table.
   */
  template <typename OSTREAM_TYPE, typename _ISOTABLE_TYPE,
	    typename _NTYPE>
  bool check_mcube_table_orientation_all_pairs
  (OSTREAM_TYPE & out, const _ISOTABLE_TYPE & isotable, 
   const bool flag_verbose, const _NTYPE output_trigger,
   IJK::ERROR & error)
  {
    typedef typename _ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
    typedef MCUBE_ISOTABLE_ORIENT_INFO
      <ISO_VERTEX_BITSET, ISO_CONNECTED_COMPONENT_BITSET, unsigned char>
      ORIENT_INFO_TABLE_TYPE;
      
    const int num_table_entries = isotable.NumTableEntries();
    const int num_vert_per_simplex = isotable.NumVerticesPerSimplex();
    int isimplexA, isimplexB;

    if (num_vert_per_simplex < 2) {
      // Nothing to check
      return true;
    }

    ORIENT_INFO_TABLE_TYPE orient_info(isotable);

    if (!orient_info.CheckOrientationsInEveryTableEntry(error))
      { return false; }

    if (flag_verbose && 
        orient_info.NumTableEntries() > output_trigger) {
      out << "  All table entries have internal orientation consistency.\n";
    }

    for (TABLE_INDEX table_indexA = 0; 
	 table_indexA < isotable.NumTableEntries(); table_indexA++) {

      for (TABLE_INDEX table_indexB = table_indexA+1;
	   table_indexB < isotable.NumTableEntries(); table_indexB++) {

	if (!orient_info.CheckOrientationsOfTwoTableEntries
	    (isotable, table_indexA, table_indexB,
	     isimplexA, isimplexB, error))
	  { return false; }
      }

      if (flag_verbose && (table_indexA > 0) &&
          table_indexA%output_trigger == 0) {
        out << "  Checked " << table_indexA
            << " out of " << isotable.NumTableEntries()
            << " isosurface table entry orientations.\n";
        out.flush();
      }
    }

    if (flag_verbose &&
        isotable.NumTableEntries() > output_trigger) {
      out << "  Checked orientations on all "
          << isotable.NumTableEntries()
          << " isosurface table entries.\n";
      out.flush();
    }

    return true;
  }

      
  /*!
   *  @brief Return true if all simplex lists in Marching Cubes lookup table
   *    are consistently oriented. (Local)
   *  @param isotable Marching Cubes lookup table.
   */
  template <typename OSTREAM_TYPE, typename _ISOTABLE_TYPE,
            typename _NTYPE>
  bool check_mcube_table_orientation
  (OSTREAM_TYPE & out, const _ISOTABLE_TYPE & isotable, 
   const bool flag_verbose,
   const _NTYPE output_trigger,
   const bool flag_check_all_pairs,
   IJK::ERROR & error)
  {
    if (flag_check_all_pairs) {
      return check_mcube_table_orientation_all_pairs
        (out, isotable, flag_verbose, output_trigger, error);
    }
    else {
      return check_mcube_table_orientation
        (out, isotable, flag_verbose, output_trigger, error);
    }
  }


  /*!
   *  @overload
   *  @brief Return true if all simplex lists in Marching Cubes lookup table
   *    are consistently oriented. (No output message.)
   *  - Version that does not print output messages.
   *  @param isotable Marching Cubes lookup table.
   */
  template <typename _ISOTABLE_TYPE>
  bool check_mcube_table_orientation
  (const _ISOTABLE_TYPE & isotable, const bool flag_check_all_pairs,
   IJK::ERROR & error)
  {
    NO_OUTPUT_MESSAGE no_output_message;

    // Dummy parameters.
    bool flag_verbose = false;
    int output_trigger(1);

    return check_mcube_table_orientation
      (no_output_message, isotable, flag_verbose, output_trigger,
       flag_check_all_pairs, error);
  }

  
  // *****************************************************************
  // MCUBE_ISOTABLE_ORIENT_INFO member functions
  // *****************************************************************

  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  Init(const ISOSURFACE_TABLE & isotable)
  {
    const int num_table_entries = isotable.NumTableEntries();
    num_vertices_per_simplex = isotable.NumVerticesPerSimplex();
    num_isosurface_vertices = isotable.NumIsosurfaceVertices();

    entry.resize(num_table_entries);
    for (TABLE_INDEX table_index = 0;
         table_index < NumTableEntries(); table_index++) {
      const ISOSURFACE_VERTEX_INDEX * simplex_vertices =
        isotable.SimplexVertices(table_index);
      const int num_simplices = isotable.NumSimplices(table_index);
      int numc;
      std::vector<_ITYPE> simplex_component;

      entry[table_index].simplex_info.resize(num_simplices);
      
      _FlagVerticesInEachSimplex
        (table_index, simplex_vertices, num_simplices);
      _SetConnectedComponentIndices
        (table_index, simplex_vertices, num_simplices);
      _SetFacetSwapParity(table_index, simplex_vertices, num_simplices);
      _FlagBoundaryFacets(table_index);
    }
  };


  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  _FlagVerticesInEachSimplex
  (const TABLE_INDEX table_index,
   const ISOSURFACE_VERTEX_INDEX * simplex_vertex_list,
   const int num_simplices)
  {
    for (int isimplex = 0; isimplex < num_simplices; isimplex++) {

      entry[table_index].simplex_info[isimplex].in_simplex.reset();

      const ISOSURFACE_VERTEX_INDEX * first_simplex_vertex =
          simplex_vertex_list + isimplex*NumVerticesPerSimplex();

      for (int j = 0; j < NumVerticesPerSimplex(); j++) {
        const ISOSURFACE_VERTEX_INDEX jw = first_simplex_vertex[j];
        entry[table_index].simplex_info[isimplex].in_simplex.set(jw);
      }
    }
  };
  

  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  _SetConnectedComponentIndices
  (const TABLE_INDEX table_index,
   const ISOSURFACE_VERTEX_INDEX * simplex_vertex_list,
   const int num_simplices)
  {
    int num_components;
    std::vector<_ITYPE> simplex_component;

    IJK::get_facet_connected_components_in_simplicial_complex
      (simplex_vertex_list, NumVerticesPerSimplex(), num_simplices,
       simplex_component, num_components);

    entry[table_index].num_connected_components = num_components;

    for (int isimplex = 0; isimplex < simplex_component.size();
         isimplex++) {
      _SetConnectedComponent
        (table_index, isimplex, simplex_component[isimplex]);
    }
  };


  // Set facet swap parity in all facets of simplex isimplex
  //   in table entry table_index.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  _SetSimplexFacetSwapParity
  (const TABLE_INDEX table_index,
   const int isimplex,
   const ISOSURFACE_VERTEX_INDEX * simplex_vertex,
   ISOSURFACE_VERTEX_INDEX temp_simplex_vert[])
  {
    int swap_parity;
    
    entry[table_index].simplex_info[isimplex].facet_swap_parity.reset();

    for (int jloc = 0; jloc < NumFacetsPerSimplex(); jloc++) {
      const ISOSURFACE_VERTEX_INDEX jw = simplex_vertex[jloc];

      IJK::sort_simplex_facet_vertices
        (simplex_vertex, NumVerticesPerSimplex(), jloc,
         temp_simplex_vert, swap_parity);

      if (swap_parity) {
        entry[table_index].simplex_info[isimplex].facet_swap_parity.set(jw);
      }
    }
  };


  // Set facet swap parity in all facets of simplex isimplex
  //   in table entry table_index. (No temporary array.)
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  _SetSimplexFacetSwapParity
  (const TABLE_INDEX table_index,
   const int isimplex,
   const ISOSURFACE_VERTEX_INDEX * simplex_vertex)
  {
    std::vector<ISOSURFACE_VERTEX_INDEX>
      temp_simplex_vert(NumVerticesPerSimplex());

    _SetSimplexFacetSwapParity
      (table_index, isimplex, simplex_vertex, temp_simplex_vert.data());
  };


  // Set facet swap parity of all facets of all simplices
  //   in table entry table_index.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  _SetFacetSwapParity
  (const TABLE_INDEX table_index, 
   const ISOSURFACE_VERTEX_INDEX * simplex_vertex_list,
   const int num_simplices)
  {
    std::vector<ISOSURFACE_VERTEX_INDEX>
      temp_simplex_vert(NumVerticesPerSimplex());
    
    for (int isimplex = 0; isimplex < num_simplices; isimplex++) {

      const ISOSURFACE_VERTEX_INDEX * first_simplex_vertex =
        simplex_vertex_list + isimplex*NumVerticesPerSimplex();

      _SetSimplexFacetSwapParity
        (table_index, isimplex, first_simplex_vertex,
         temp_simplex_vert.data());
    }

  };


  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  _FlagBoundaryFacets(const TABLE_INDEX table_index)
  {
    IJK::PROCEDURE_ERROR
      error("MCUBE_ISOTABLE_ORIENT_INFO::_FlagBoundaryFacets");

    if (!CheckInSimplexIsSet(table_index, error)) {
      error.AddMessage
        ("  Call _FlagVerticesInEachSimplex() before _FlagBoundaryFacets().");
      throw error;
    }

    // Compute boundary simplex facets.
    for (int isimplexA = 0; isimplexA < NumSimplices(table_index);
         isimplexA++) {

      entry[table_index].simplex_info[isimplexA].is_boundary_facet.reset();
      const ISO_VERTEX_BITSET_TYPE in_simplexA =
        entry[table_index].simplex_info[isimplexA].in_simplex;

      // Initialize is_boundary_facet to true for all facets.
      entry[table_index].simplex_info[isimplexA].is_boundary_facet =
        in_simplexA;
      
      for (int isimplexB = 0; isimplexB < NumSimplices(table_index);
           isimplexB++) {
        
        if (isimplexA == isimplexB)
          { continue; }

        const ISO_VERTEX_BITSET_TYPE in_simplexB =
          entry[table_index].simplex_info[isimplexB].in_simplex;

        const ISO_VERTEX_BITSET_TYPE shared_vert =
          in_simplexA & in_simplexB;

        const int num_ones = shared_vert.count();
        
        if (num_ones == NumVerticesPerSimplex()) {
          // Duplicate simplex in list of simplices?!?
          // Ignore simplexB.
          continue;
        }

        if (num_ones == NumVerticesPerSimplexFacet()) {
          // Simplices isimplexA and isimplexB share a facet.
          // Mark facet as not on boundary.

          // Compute bit indicating vertex not in facet.
          ISO_VERTEX_BITSET_TYPE not_in_facet_bitset =
            in_simplexA ^ shared_vert;


          // Flip bit to off.
          not_in_facet_bitset.flip();
          entry[table_index].simplex_info[isimplexA].is_boundary_facet &=
            not_in_facet_bitset;
        }
      }
    }
  };


  // Return true if all table entries are marked as oriented.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  AreAllTableEntriesOriented(TABLE_INDEX & table_index) const
  {
    table_index = 0;
    for (TABLE_INDEX table_indexB = 0; table_indexB < NumTableEntries();
         table_indexB++) {

      if (!AreAllSimplicesOriented(table_indexB)) {
        table_index = table_indexB;
        return false;
      }
    }

    return true;
  }
  

  // Return true if some simplex in connected component
  //   icomponent has a boundary facet.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  ComponentHasBoundaryFacet
  (const TABLE_INDEX table_index, const int icomponent) const
  {
    for (int isimplex = 0; isimplex < NumSimplices(table_index); isimplex++) {
      if (IndexOfConnectedComponent(table_index, isimplex) == icomponent) {
	if (SimplexInfo(table_index, isimplex).HasBoundaryFacet())
	  { return true; }
      }
    }

    return false;
  }


  // Return true if some simplex in table entry has a boundary facet.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  SomeSimplexHasBoundaryFacet
  (const TABLE_INDEX table_index) const
  {
    for (int isimplex = 0; isimplex < NumSimplices(table_index); isimplex++) {
      if (SimplexInfo(table_index, isimplex).HasBoundaryFacet())
	{ return true; }
    }

    return false;
  }


  // Find index of table entry with a single connected component 
  //   that has at least one boundary simplex facet.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  template <typename _TABLE_INDEX>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  FindEntryWithSingleConnectedComponentWithBoundary
  (_TABLE_INDEX & table_index, bool & flag_found) const
  {
    const int ONE(1);

    // Initialize.
    table_index = 0;
    flag_found = false;

    for (TABLE_INDEX table_indexA = 0; table_indexA < NumTableEntries();
	 table_indexA++) {

      if (NumConnectedComponents(table_indexA) == ONE) {

	// Single connected component has index 0.
	if (ComponentHasBoundaryFacet(table_indexA, 0)) {
	  table_index = table_indexA;
	  flag_found = true;
	  return;
	}
      }
    }
  }


  // Convert ISO_VERTEX_BITSET_TYPE to string.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  std::string MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  _ConvertBitsetToString(const ISO_VERTEX_BITSET_TYPE & bitset) const
  {
    const int num_isov = NumIsosurfaceVertices();
    
    std::string bitset_str = bitset.to_string();
    bitset_str = bitset_str.substr(bitset.size()-num_isov);
    return bitset_str;
  };
  

  // Set flag is_connected_component_oriented 
  //   for connected component icomponent.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  SetIsConnectedComponentOriented
    (const TABLE_INDEX table_index, const int icomponent,
     const bool flag)
  {
    entry[table_index].is_connected_component_oriented[icomponent] = flag;
    entry[table_index].num_oriented_connected_components =
      entry[table_index].is_connected_component_oriented.count();
  };


  // Update flag entry[table_index].is_connected_component_oriented. 
  //   - Update by "or" with is connected_component_oriented.
  // - Also modifies entry[table_index].num_oriented_connected_components.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  UpdateIsConnectedComponentOriented
  (const TABLE_INDEX table_index, 
   const ISO_CONNECTED_COMPONENT_BITSET & is_component_oriented)
  {
    entry[table_index].is_connected_component_oriented |=
      is_component_oriented;
    entry[table_index].num_oriented_connected_components =
      entry[table_index].is_connected_component_oriented.count();
  }

  // Return true if two simplices share a facet.
  // - Also returns bitset of shared vertices, number of shared vertices,
  //   comparison of facet swap parities, and bitsets indicating 
  //   unshared vertices.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  DoSimplicesShareFacet
  (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
   const int isimplexA, const int isimplexB,
   ISO_VERTEX_BITSET_TYPE & shared_vert,
   int & num_shared_vertices,
   bool & are_shared_facet_swap_parities_equal,
   ISO_VERTEX_BITSET_TYPE & not_in_facetA,
   ISO_VERTEX_BITSET_TYPE & not_in_facetB) const
  {
    const ISO_VERTEX_BITSET_TYPE in_simplexA =
      SimplexInfo(table_indexA, isimplexA).in_simplex;
    const ISO_VERTEX_BITSET_TYPE in_simplexB =
      SimplexInfo(table_indexB, isimplexB).in_simplex;

    shared_vert = in_simplexA & in_simplexB;
    num_shared_vertices = shared_vert.count();

    if (num_shared_vertices == NumVerticesPerSimplex()) {

      // Duplicate simplex.
      const ISO_VERTEX_BITSET_TYPE facet_swap_parityA =
        SimplexInfo(table_indexA, isimplexA).FacetSwapParity();
      const ISO_VERTEX_BITSET_TYPE facet_swap_parityB =
        SimplexInfo(table_indexB, isimplexB).FacetSwapParity();

      are_shared_facet_swap_parities_equal =
	(facet_swap_parityA == facet_swap_parityB);

      return true;
    }
    else if (num_shared_vertices == NumVerticesPerSimplexFacet()) {
      
      // Simplices share one facet.
      const ISO_VERTEX_BITSET_TYPE facet_swap_parityA =
        SimplexInfo(table_indexA, isimplexA).facet_swap_parity;
      const ISO_VERTEX_BITSET_TYPE facet_swap_parityB =
        SimplexInfo(table_indexB, isimplexB).facet_swap_parity;
      
      not_in_facetA = in_simplexA ^ shared_vert;
      not_in_facetB = in_simplexB ^ shared_vert;

      const bool shared_facet_swap_parityA =
        bool((not_in_facetA & facet_swap_parityA) == 0);
      const bool shared_facet_swap_parityB =
        bool((not_in_facetB & facet_swap_parityB) == 0);

      are_shared_facet_swap_parities_equal =
	(shared_facet_swap_parityA == shared_facet_swap_parityB);

      return true;
    }
    else {
      // Clear arguments. (For code safety.)
      are_shared_facet_swap_parities_equal = false;
      not_in_facetA.reset();
      not_in_facetB.reset();

      // Simplices do not share a facet.
      return false;
    }

  }


  // Return true if two simplices share a facet.
  //  - Version that returns only comparison of swap parities.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  DoSimplicesShareFacet
  (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
   const int isimplexA, const int isimplexB,
   bool & are_shared_facet_swap_parities_equal) const
  {
    ISO_VERTEX_BITSET_TYPE  shared_vert;
    int num_shared_vertices;
    ISO_VERTEX_BITSET_TYPE not_in_facetA;
    ISO_VERTEX_BITSET_TYPE not_in_facetB;

    return DoSimplicesShareFacet
      (table_indexA, table_indexB, isimplexA, isimplexB,
       shared_vert, num_shared_vertices,
       are_shared_facet_swap_parities_equal,
       not_in_facetA, not_in_facetB);
  }


  // Return true if two simplices share a boundary facet.
  // - Also returns bitset of shared vertices, number of shared vertices,
  //   comparison of facet swap parities, and bitsets indicating 
  //   unshared vertices.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  DoSimplicesShareBoundaryFacet
  (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
   const int isimplexA, const int isimplexB,
   ISO_VERTEX_BITSET_TYPE & shared_vert,
   int & num_shared_vertices,
   bool & are_shared_facet_swap_parities_equal,
   ISO_VERTEX_BITSET_TYPE & not_in_facetA,
   ISO_VERTEX_BITSET_TYPE & not_in_facetB) const
  {
    if (DoSimplicesShareFacet
	(table_indexA, table_indexB, isimplexA, isimplexB,
	 shared_vert, num_shared_vertices,
	 are_shared_facet_swap_parities_equal,
	 not_in_facetA, not_in_facetB)) {

      // Simplices A and B have a common facet.

      const ISO_VERTEX_BITSET_TYPE is_boundary_facetA =
	SimplexInfo(table_indexA, isimplexA).IsBoundaryFacet();
      const ISO_VERTEX_BITSET_TYPE is_boundary_facetB =
	SimplexInfo(table_indexB, isimplexB).IsBoundaryFacet();

      if (num_shared_vertices == NumVerticesPerSimplex()) {

	// Duplicate simplex.

	if ((is_boundary_facetA & is_boundary_facetB).none()) {
	  // No facet is a boundary facet in both table entries.
	  return false;
	}

	return true;
      }
      else {
	// Simplices share exactly one facet.
	
	if ((not_in_facetA & is_boundary_facetA).none()) {
	  // Shared facet is not a boundary facet in table_indexA.
	  return false;
	}

	if ((not_in_facetB & is_boundary_facetB).none()) {
	  // Shared facet is not a boundary facet in table_indexB.
	  return false;
	}

	// Shared facet is a boundary facet in both table entries.
	return true;
      }
    }

    return false;
  }


  // Return true if two simplices share a boundary facet.
  // - Version that does not return bitsets indicating 
  //   unshared vertices.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  DoSimplicesShareBoundaryFacet
  (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
   const int isimplexA, const int isimplexB,
   ISO_VERTEX_BITSET_TYPE & shared_vert,
   int & num_shared_vertices,
   bool & are_shared_facet_swap_parities_equal) const
  {
    ISO_VERTEX_BITSET_TYPE not_in_facetA;
    ISO_VERTEX_BITSET_TYPE not_in_facetB;

    return DoSimplicesShareBoundaryFacet
      (table_indexA, table_indexB, isimplexA, isimplexB,
       shared_vert, num_shared_vertices,
       are_shared_facet_swap_parities_equal,
       not_in_facetA, not_in_facetB);
  }

  
  // Orient all simplices in connected component containing isimplexA
  //   consistently with isimplexA.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::    
  OrientConnectedComponent
  (const TABLE_INDEX table_index, const int isimplexA,
   ISOSURFACE_TABLE & isotable)
  {
    const int num_simplices = NumSimplices(table_index);    
    std::vector<bool> is_oriented(num_simplices, false);
    std::stack<int> simplex_stack;

    is_oriented[isimplexA] = true;
    simplex_stack.push(isimplexA);
    while (simplex_stack.size() != 0) {
      const int isimplexB = simplex_stack.top();
      simplex_stack.pop();      

      const int icomponentB = 
	IndexOfConnectedComponent(table_index, isimplexB);

      for (int isimplexC = 0; isimplexC < num_simplices; isimplexC++) {

        if (isimplexB == isimplexC)
          { continue; }
        
        if (is_oriented[isimplexC]) {
          // Simplex C is already oriented.
          continue;
        }

        const int icomponentC =
	  IndexOfConnectedComponent(table_index, isimplexC);

        if (icomponentB != icomponentC) {
          // isimplexB and isimplexC are in different connected components.
          continue;
        }

	bool are_shared_facet_swap_parities_equal;

	if (DoSimplicesShareFacet
	    (table_index, table_index, isimplexB, isimplexC,
	     are_shared_facet_swap_parities_equal)) {

          if (are_shared_facet_swap_parities_equal) {
            // Simplex C has opposite orientation to simplexB.
            isotable.FlipIsoPolyOrientation(table_index, isimplexC);

            // Set facet swap parity.
            const ISOSURFACE_VERTEX_INDEX * simplexC_vertices =
              isotable.SimplexVertices(table_index, isimplexC);

            _SetSimplexFacetSwapParity
              (table_index, isimplexC, simplexC_vertices);
          }

          is_oriented[isimplexC] = true;
          simplex_stack.push(isimplexC);
	}
      }
    }
  }

  
  // Consistently orient all simplices in each
  //   connected component in table entry table_index.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::    
  OrientAllSimplicesInTableEntry
  (const TABLE_INDEX table_index, ISOSURFACE_TABLE & isotable)
  {
    const int num_simplices = NumSimplices(table_index);    
    const int num_connected_components =
      NumConnectedComponents(table_index);
    std::vector<bool> is_oriented(num_connected_components, false);

    for (int isimplex = 0; isimplex < num_simplices; isimplex++) {
      const int icomponent =
        SimplexInfo(table_index, isimplex).IndexOfConnectedComponent();

      if (!is_oriented[icomponent]) {
        OrientConnectedComponent(table_index, isimplex, isotable);
        is_oriented[icomponent] = true;
      }
    }      
  };


  // Flip all simplices in connected component icomponent
  //   in table entry table_index.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  FlipSimplicesInConnectedComponent
    (const TABLE_INDEX table_index, const int icomponent,
     ISOSURFACE_TABLE & isotable)
  {
    for (int isimplex = 0; isimplex < NumSimplices(table_index);
         isimplex++) {
      if (SimplexInfo(table_index,isimplex).IndexOfConnectedComponent() ==
          icomponent) {
        const ISOSURFACE_VERTEX_INDEX * simplex_vertices =
          isotable.SimplexVertices(table_index, isimplex);
        
        isotable.FlipIsoPolyOrientation(table_index,isimplex);
        _SetSimplexFacetSwapParity
          (table_index, isimplex, simplex_vertices);
      }
    }
  };
  

  // Orient simplices in table entry table_indexB
  //   from orientation of some boundary facet of table_indexA.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  OrientTwoTableEntries
  (const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
   ISOSURFACE_TABLE & isotable)
  {
    ISO_VERTEX_BITSET_TYPE shared_vert;
    int num_shared_vertices;

    if (AreAllSimplicesOriented(table_indexB)) {
      // Nothing left to orient.
      return;
    }

    for (int isimplexB = 0; isimplexB < NumSimplices(table_indexB);
         isimplexB++) {

      if (IsSimplexOriented(table_indexB, isimplexB))
        { continue; }

      bool flag_matched = false;
      for (int isimplexA = 0; 
           isimplexA < NumSimplices(table_indexA) && !flag_matched;
           isimplexA++) {

        bool are_shared_facet_swap_parities_equal;
        if (DoSimplicesShareBoundaryFacet
            (table_indexA, table_indexB, isimplexA, isimplexB,
             shared_vert, num_shared_vertices, 
             are_shared_facet_swap_parities_equal)) {

          for (int ifacet = 0; ifacet < facet_info.NumFacets(); ifacet++) {
            if (facet_info.AreVerticesInFacet(shared_vert, ifacet)) {

              if (isotable.AreAllFacetVertexLabelsIdentical
                  (table_indexA, table_indexB, ifacet)) {

                const int icomponentB = 
                  IndexOfConnectedComponent(table_indexB, isimplexB);

                if (!are_shared_facet_swap_parities_equal) {

                  // Simplices have different orientations.

                  // Reverse orientaton of all simplices in icomponentB.
                  FlipSimplicesInConnectedComponent
                    (table_indexB, icomponentB, isotable);

                  // Shared facet swap parities are now equal.
                  are_shared_facet_swap_parities_equal = true;
                }

                // Flag connected component icomponentB as oriented.
                SetIsConnectedComponentOriented
                  (table_indexB, icomponentB, true);

                if (AreAllSimplicesOriented(table_indexB)) {
                  // All simplices are oriented.
                  return;
                }

                flag_matched = true;
                break;
              }
            }
          }
        }
      }
    }
  }
  

  // *****************************************************************
  // MCUBE_ISOTABLE_ORIENT_INFO check routines.
  // *****************************************************************

  // Check against isosurface lookup table.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  template <typename _ISOTABLE_TYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  Check(const _ISOTABLE_TYPE & isotable, IJK::ERROR & error) const
  {
    if (NumIsosurfaceVertices() != isotable.NumIsosurfaceVertices()) {
      error.AddMessage
        ("Programming error. Incorrect value for NumIsosurfaceVertices().");
      return false;
    }
    
    if (NumVerticesPerSimplex() != isotable.NumVerticesPerSimplex()) {
      error.AddMessage
        ("Programming error. Incorrect value for NumVerticesPerSimplex().");
      return false;
    }

    if (NumTableEntries() != isotable.NumTableEntries()) {
      error.AddMessage
        ("Programming error. Incorrect value for NumTableEntries().");
      return false;
    }

    for (TABLE_INDEX table_index = 0; table_index < isotable.NumTableEntries();
       table_index++) {

      if (NumSimplices(table_index) != isotable.NumSimplices(table_index)) {
        error.AddMessage
          ("Programming error. Incorrect number of simplices for table entry ",
           table_index, ".");
        return false;
      }

      if (!CheckInSimplexIsSet(table_index, error))
        { return false; }

      for (int isimplex = 0; isimplex < isotable.NumSimplices(table_index);
           isimplex++) {
        for (int j = 0; j < isotable.NumVerticesPerSimplex(); j++) {
          const ISOSURFACE_VERTEX_INDEX iw =
            isotable.SimplexVertex(table_index, isimplex, j);

          if (!SimplexInfo(table_index, isimplex).InSimplex(iw)) {
            error.AddMessage
              ("Programming error. Missing isosurface vertex ", int(iw),
               " in table entry ", table_index, ", simplex ", isimplex, ".");
            return false;
          }
        }
      }
    }

    return true;
  }


  // Check that bitset entry[table_index].simplex_info[*].in_simplex is set.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckInSimplexIsSet(const TABLE_INDEX table_index, IJK::ERROR & error) const
  {
    for (int isimplex = 0; isimplex < NumSimplices(table_index);
         isimplex++) {
           
      const int num_ones =
        SimplexInfo(table_index,isimplex).in_simplex.count();
    
      if (num_ones != NumVerticesPerSimplex()) {
        if (num_ones == 0) {
          error.AddMessage
            ("Programming error. Bitset in_simplex not set.");
          error.AddMessage
            ("  Table index: ", table_index, "  Simplex: ", isimplex);
        }
        else {
          error.AddMessage
            ("Programming error. Incorrect number of ones bitset in_simplex not set.");
          error.AddMessage
            ("  Table index: ", table_index, "  Simplex: ", isimplex);
          error.AddMessage
            ("  Bitset in_simplex has ", num_ones, " ones.");
          error.AddMessage
            ("  Number of ones should match number of simplex vertices, ",
             NumVerticesPerSimplex(), ".");
        }
      
        return false;
      }
    }

    return true;
  }

  

  // Check that orientations of all pairs of simplices in entry[table_index]
  //   are consistent.
  // - Return true if all orientations in entry[table_index] are consistent.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckOrientationsInTableEntry
  (const TABLE_INDEX table_index, int & isimplexA, int & isimplexB,
   IJK::ERROR & error) const
  {
    isimplexA = 0;
    isimplexB = 0;
    
    for (int jsimplexA = 0; jsimplexA+1 < entry[table_index].NumSimplices();
         jsimplexA++) {
      for (int jsimplexB = jsimplexA+1; 
	   jsimplexB < entry[table_index].NumSimplices();
	   jsimplexB++) {

	bool are_shared_facet_swap_parities_equal;
	if (DoSimplicesShareFacet
	    (table_index, table_index, jsimplexA, jsimplexB,
	     are_shared_facet_swap_parities_equal)) {

	  if (are_shared_facet_swap_parities_equal) {

	    isimplexA = jsimplexA;
	    isimplexB = jsimplexB;

            // Simplex B has opposite orientation to simplex A.
	    error.AddMessage
	      ("Simplices ", isimplexA, " and ", isimplexB,
	       " in table entry ", table_index,
	       " are not consistently oriented.");

	    return false;
	  }
	}
      }
    }
    
    return true;
  }


  // Check that orientations of all pairs of simplices 
  //   in entry[table_index] are consistent. (Doesn't return simplices.)
  // - Version that does not return mismatched simplices.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckOrientationsInTableEntry
  (const TABLE_INDEX table_index, IJK::ERROR & error) const
  {
    int isimplexA, isimplexB;
    
    return CheckOrientationsInTableEntry
      (table_index, isimplexA, isimplexB, error);
  }


  // Check that orientations in every table entry are consistent.
  //  - Does NOT compare orientations of simplices in DIFFERENT table entries.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckOrientationsInEveryTableEntry
  (TABLE_INDEX & table_index, int & isimplexA, int & isimplexB,
   IJK::ERROR & error) const
  {
    const int num_vert_per_simplex = NumVerticesPerSimplex();
    int jsimplexA, jsimplexB;

    // Initialize.
    table_index = 0;
    isimplexA = isimplexB = 0;

    if (num_vert_per_simplex < 2) {
      // Nothing to check
      return true;
    }

    for (TABLE_INDEX table_indexA = 0;
	 table_indexA < NumTableEntries(); table_indexA++) {

      if (!CheckOrientationsInTableEntry
	  (table_indexA, jsimplexA, jsimplexB, error)) {

	table_index = table_indexA;
	isimplexA = jsimplexA;
	isimplexB = jsimplexB;

	return false;
      }
    }

    return true;
  }


  // Check that orientations in every table entry are consistent. 
  //    (Doesn't return table index or simplices.)
  //  - Does NOT compare orientations of simplices in DIFFERENT table entries.
  //  - Version that does not return index of table entry or 
  //    mismatched simplices.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckOrientationsInEveryTableEntry(IJK::ERROR & error) const
  {
    TABLE_INDEX table_index;
    int isimplexA, isimplexB;

    return CheckOrientationsInEveryTableEntry
      (table_index, isimplexA, isimplexB, error);
  }


  // Return true if orientations of simplices in entry[table_indexB] 
  //   are consistent with orientations of simplices in entry[table_indexA].
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  template <typename _ISOTABLE_TYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckOrientationsOfTwoTableEntries
  (const _ISOTABLE_TYPE & isotable,
   const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
   int & isimplexA, int & isimplexB,
   _ISO_CONNECTED_COMPONENT_BITSET_TYPE & componentB_checked,
   IJK::ERROR & error) const
  {
    ISO_VERTEX_BITSET_TYPE shared_vert;
    int num_shared_vertices;

    // Initialize.
    componentB_checked.reset();

    for (int jsimplexB = 0; jsimplexB < NumSimplices(table_indexB);
	 jsimplexB++) {

      const int icomponentB =
        IndexOfConnectedComponent(table_indexB, jsimplexB);

      if (bool(componentB_checked[icomponentB])) {
        // Orientation of some simplex in connected component icomponentB
        //   has already been checked.
        continue;
      }

      bool flag_checked = false;
      for (int jsimplexA = 0; 
           jsimplexA < NumSimplices(table_indexA) && !flag_checked;
           jsimplexA++) {

        bool are_shared_facet_swap_parities_equal;
        if (DoSimplicesShareBoundaryFacet
            (table_indexA, table_indexB, jsimplexA, jsimplexB,
             shared_vert, num_shared_vertices,
             are_shared_facet_swap_parities_equal)) {

          for (int ifacet = 0; ifacet < facet_info.NumFacets(); ifacet++) {
            if (facet_info.AreVerticesInFacet(shared_vert, ifacet)) {

              if (isotable.AreAllFacetVertexLabelsIdentical
                  (table_indexA, table_indexB, ifacet)) {

                if (!are_shared_facet_swap_parities_equal) {

                  // Simplices have different orientations.
                  isimplexA = jsimplexA;
                  isimplexB = jsimplexB;

                  error.AddMessage
                    ("  Simplex ", isimplexA, " in table entry ", table_indexA,
                     " has inconsistent orientation");
                  error.AddMessage
                    ("  with simplex ", isimplexB, " in table entry ",
                     table_indexB, ".");

                  return false;
                }

                // Flag connected component icomponentB as checked.
                componentB_checked[icomponentB] = 1;
                flag_checked = true;

                continue;
              }
            }
          }
        }
      }
    }

    return true;
  }

  
  // Return true if orientations of simplices in entry[table_indexB] 
  //   are consistent with orientations of simplices in entry[table_indexA].
  // - Version that does not return bitset componentB_checked.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  template <typename _ISOTABLE_TYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckOrientationsOfTwoTableEntries
  (const _ISOTABLE_TYPE & isotable,
   const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
   int & isimplexA, int & isimplexB,
   IJK::ERROR & error) const
  {
    ISO_CONNECTED_COMPONENT_BITSET_TYPE  componentB_checked;

    return CheckOrientationsOfTwoTableEntries
      (isotable, table_indexA, table_indexB, isimplexA, isimplexB,
       componentB_checked, error);
  }


  // Check orientation of table entry table_indexA against
  //    all other table entries.
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  template <typename _ISOTABLE_TYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckOrientationOfTableEntryAgainstAllOthers
  (const _ISOTABLE_TYPE & isotable,
   const TABLE_INDEX table_indexA, 
   int & isimplexA, TABLE_INDEX & table_indexB, int & isimplexB,
   IJK::ERROR & error) const
  {
    // Initialize.
    isimplexA = 0;
    table_indexB = 0;
    isimplexB = 0;


    for (int table_index = 0; table_index < NumTableEntries();
	 table_index++) {

      if (table_index == table_indexA)
	{ continue; }

      if (!CheckOrientationsOfTwoTableEntries
	  (isotable, table_indexB, table_index, 
	   isimplexA, isimplexB, error)) {

	    table_indexB = table_index;
	    return false;
      }
    }

    return true;
  }


  // Check orientation of table entry table_indexA against
  //   all other table entries. (Doesn't return table index or simplices.)
  template <typename _ISO_VERTEX_BITSET_TYPE,
            typename _ISO_CONNECTED_COMPONENT_BITSET_TYPE,
            typename _ITYPE>
  template <typename _ISOTABLE_TYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO
  <_ISO_VERTEX_BITSET_TYPE,_ISO_CONNECTED_COMPONENT_BITSET_TYPE,_ITYPE>::
  CheckOrientationOfTableEntryAgainstAllOthers
  (const _ISOTABLE_TYPE & isotable,
   const TABLE_INDEX table_indexA, 
   IJK::ERROR & error) const
  {
    TABLE_INDEX table_indexB;
    int isimplexA, isimplexB;

    return CheckOrientationOfTableEntryAgainstAllOthers
      (isotable, table_indexA, isimplexA, table_indexB,
       isimplexB, error);
  }


  // *****************************************************************
  // FACET_ISO_VERTEX member functions
  // *****************************************************************

  template <typename _ISO_VERTEX_BITSET_TYPE>
  void FACET_ISO_VERTEX<_ISO_VERTEX_BITSET_TYPE>::
  Init(const ISOSURFACE_TABLE & isotable) {

    num_facets = isotable.Polytope().NumFacets();
    num_isosurface_vertices =
      isotable.NumIsosurfaceVertices();

    in_facet.resize(num_facets);

    for (int ifacet = 0; ifacet < NumFacets(); ifacet++) {
      in_facet[ifacet].reset();

      for (int iw = 0; iw < isotable.NumIsosurfaceVertices(); iw++) {
        if (isotable.IsosurfaceVertex(iw).Type() == ISOSURFACE_VERTEX::VERTEX) {
          const int iv = isotable.IsosurfaceVertex(iw).Face();
          if (isotable.Polytope().IsVertexInFacet(ifacet, iv))
            { in_facet[ifacet][iw] = 1; }
        }
        else if (isotable.IsosurfaceVertex(iw).Type() == ISOSURFACE_VERTEX::EDGE) {
          const int ie = isotable.IsosurfaceVertex(iw).Face();
          const int iend0 = isotable.Polytope().EdgeEndpoint(ie, 0);
          const int iend1 = isotable.Polytope().EdgeEndpoint(ie, 1);
          if (isotable.Polytope().IsVertexInFacet(ifacet, iend0) &&
              isotable.Polytope().IsVertexInFacet(ifacet, iend1))
            { in_facet[ifacet][iw] = 1; }
        }
        else if (isotable.IsosurfaceVertex(iw).Type() == ISOSURFACE_VERTEX::FACET) {
          const int jfacet = isotable.IsosurfaceVertex(iw).Face();
          if (ifacet == jfacet)
            { in_facet[ifacet][iw] = 1; }
        }
        else {
          // Vertex is not on any facet.
          continue;
        }
      }
    }
  };


  // Return string representing bitset in_facet[ifacet].
  template <typename _ISO_VERTEX_BITSET_TYPE>
  std::string FACET_ISO_VERTEX<_ISO_VERTEX_BITSET_TYPE>::
  InFacetStr(const int ifacet) const
  {
    const int num_isov = NumIsosurfaceVertices();
    
    std::string in_facet_str = in_facet[ifacet].to_string();
    in_facet_str = in_facet_str.substr(in_facet[ifacet].size()-num_isov);
    return in_facet_str;
  }

}

#endif

