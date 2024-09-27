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

  // *** NEW ***
  
  /*!
   *  @brief Size of bitset used to represent isosurface vertices
   *    in isotable entries.
   *  - Note: Number of isosurface vertices can be 
   *    at most ISO_VERTEX_BITSET_SIZE.
   */
  const int ISO_VERTEX_BITSET_SIZE = 64;

  typedef std::bitset<ISO_VERTEX_BITSET_SIZE> ISO_VERTEX_BITSET;

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

    
    // *** Return string representation. Mainly for debugging. ***

    /// @brief Return string representing bitset in_facet[ifacet].
    /// - Only returns num_isosurface_vertices bits,
    ///   since all other bits are zero.
    std::string InFacetStr(const int ifacet) const;
  };


  // *** NEW ***
  /// @brief Table containing orientation information.
  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  class MCUBE_ISOTABLE_ORIENT_INFO {

  public:
    typedef _ISO_VERTEX_BITSET_TYPE ISO_VERTEX_BITSET_TYPE;
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

      /// @brief is_oriented[js] Flags whether orientation
      ///   of simplex js has been determined.
      std::vector<bool> is_oriented;

      /// @brief Number of connected components.
      _ITYPE num_connected_components;

      /// @brief Set index of connected component containing isimplex.
      void SetConnectedComponent
      (const int isimplex, const _ITYPE icomponent)
      { simplex_info[isimplex].SetConnectedComponent(icomponent); }

      /// Return number of simplices in entry.
      int NumSimplices() const
      { return simplex_info.size(); }
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
    { entry[table_index].SetConnectedComponent(isimplex,icomponent); }

    /*!
     *  @brief Flag vertices in each simplex.
     */
    void _FlagVerticesInEachSimplex
    (const TABLE_INDEX table_index, 
     const ISOSURFACE_VERTEX_INDEX * simplex_vertex_list,
     const int num_simplices);
    
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

    /// @brief Return flag is_oriented for simplex isimplex
    ///   in table entry table_index.
    bool IsOriented(const TABLE_INDEX table_index,
                    const int isimplex)
    { return entry[table_index].is_oriented[isimplex]; }

    /// @brief Return number of connected components
    ///   in table entry table_index
    int NumConnectedComponents
    (const TABLE_INDEX table_index) const
    { return entry[table_index].num_connected_components; }
    
    /// @brief Return Number of facets per simplex.
    int NumFacetsPerSimplex() const
    { return NumVerticesPerSimplex(); }

    /// @brief Return Number of vertices per simplex facet.
    int NumVerticesPerSimplexFacet() const
    { return NumVerticesPerSimplex()-1; }

    
    // *** Set routines ***
    
    /// @brief Set flag is_oriented for simplex isimplex
    ///   in table entry table_index.
    void SetIsOriented(const TABLE_INDEX table_index,
                       const int isimplex,
                       const bool flag)
    { entry[table_index].is_oriented[isimplex] = flag; }


    // *** Check routines ***

    /// @brief check that bitset entry[table_index].simplex_info[*].in_simplex is set.
    bool CheckInSimplexIsSet
    (const TABLE_INDEX, IJK::ERROR & error) const;

    
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

    /// Table containing orientation information.
    template <typename _TABLE_INDEX, typename _VTYPE,
              typename _STYPE, typename _PTYPE, typename _NTYPE>
    class MCUBE_ORIENT_INFO_TABLE{

    public:
      typedef _TABLE_INDEX TABLE_INDEX;
      typedef _VTYPE VERTEX_TYPE;
      typedef _PTYPE PARITY_TYPE;


    protected:
      
      /// Table entry.
      class ORIENT_INFO_ENTRY {

      public:
        std::vector<_VTYPE> boundary_facet_vertex_list;
        std::vector<_STYPE> simplex_containing_boundary_facet;
        std::vector<_PTYPE> boundary_facet_swap_parity;
        std::vector<bool> is_boundary_facet_oriented;
        _NTYPE num_connected_components;
        bool is_oriented;

        // Constructor
        ORIENT_INFO_ENTRY()
        {
          num_connected_components = 0;
          is_oriented = false;
        }

        /// Return number of boundary facets.
        int NumBoundaryFacets() const
        { return (boundary_facet_swap_parity.size()); }

        /// Return number of connected components.
        int NumConnectedComponents() const
        { return num_connected_components; }

        /// Return is_oriented.
        bool IsOriented() const
        { return is_oriented; }
        
        /// Return true if all boundary facets are oriented.
        bool AreAllBoundaryFacetsOriented() const
        {
          for (int i = 0; i < NumBoundaryFacets(); i++) {
            if (!is_boundary_facet_oriented[i])
              { return false; }
          }
          
          return true;
        }
      };

      
    protected:
      int num_vertices_per_facet;
      std::vector<ORIENT_INFO_ENTRY> entry;

      
    public:
      
      /// Constructor.
      MCUBE_ORIENT_INFO_TABLE
      (const int _num_vertices_per_facet,
       const _TABLE_INDEX _num_table_entries)
      {
        num_vertices_per_facet = _num_vertices_per_facet;
        entry.resize(_num_table_entries);       
      }

      // Get functions.

      int NumVerticesPerFacet() const
      { return num_vertices_per_facet; }

      /// Return number of table entries.
      int NumTableEntries() const
      { return entry.size(); }

      /// Return number of boundary facets.
      int NumBoundaryFacets(const _TABLE_INDEX table_index) const
      { return entry[table_index].NumBoundaryFacets(); }

      /// Return number of connected components.
      int NumConnectedComponents(const _TABLE_INDEX table_index) const
      { return entry[table_index].NumConnectedComponents(); }

      /// Return entry[table_index].is_oriented.
      bool IsOriented(const _TABLE_INDEX table_index) const
      { return entry[table_index].IsOriented(); }
      
      /// Return entry[table_index].
      const ORIENT_INFO_ENTRY & operator[]
      (const _TABLE_INDEX table_index) const
      { return entry[table_index]; }
      
      /// @brief Return point to boundary facet vertex list
      ///   for entry[table_index].
      const _VTYPE * BoundaryFacetVertices
      (const _TABLE_INDEX table_index, const int ifacet) const
      { return (entry[table_index].boundary_facet_vertex_list.data() +
                NumVerticesPerFacet()*ifacet); }

      /// Return true if boundary_facet_vertex_list[] contains
      ///   facet facet_vert[].
      template <typename _VTYPEF, typename _ITYPE>
      bool ContainsBoundaryFacet
      (const _TABLE_INDEX table_index, const _VTYPEF * facet_vert,
       _ITYPE & ifacet) const;
      
      /*!
       *  @brief Set boundary facets.
       *  - Boundary facets are facets that are in only one simplex
       *    in simplex_vertex_list[].
       *  @param simplex_vertex_list[] List of simplices.
       *    - Each simplex has (NumVerticesPerFacet()+1) vertices.
       */
      void SetBoundaryFacets
      (const _TABLE_INDEX table_index, 
       const _VTYPE * simplex_vertex_list,
       const int num_simplices);

      /// Set entry[table_index].num_connected_components.
      void SetNumConnectedComponents
      (const _TABLE_INDEX table_index, const int numc)
      { entry[table_index].num_connected_components = numc; }

      /// Set entry[table_index].is_oriented.
      void SetIsOriented
      (const _TABLE_INDEX table_index, const bool flag)
      { entry[table_index].is_oriented = flag; }

      /// Set entry[table_index].is_boundary_facet_oriented.
      void SetIsBoundaryFacetOriented
      (const _TABLE_INDEX table_index, const int ifacet,
       const bool flag)
      { entry[table_index].is_boundary_facet_oriented[ifacet] = flag; }

      /// Set ORIENT_INFO_TABLE from ISOSURFACE_TABLE.
      void Set(const ISOSURFACE_TABLE & isotable);
    };

    
    /// Table indicating which simplices are oriented.
    template <typename _TABLE_INDEX, typename _ITYPEC,
              typename _VTYPE, typename _STYPE, typename _NTYPEC>
    class MCUBE_ORIENT_SIMPLEX_TABLE {

    public:
      typedef _TABLE_INDEX TABLE_INDEX;
      typedef _ITYPEC COMPONENT_INDEX;
      typedef _NTYPEC NUM_TYPE;

    protected:
      
      /// Table entry.
      class ORIENT_SIMPLEX_ENTRY {

      public:

        /// @brief Facet connected component containing simplex i.
        std::vector<COMPONENT_INDEX> connected_component;

        /// @brief Number of connected components.
        NUM_TYPE num_connected_components;

        /// @brief Indicate if simplex i is consistently oriented.
        std::vector<bool> is_simplex_oriented;

        /// @brief If true, all simplices in table entry are consistently oriented.
        bool is_oriented;

        // Constructor
        ORIENT_SIMPLEX_ENTRY()
        {
          is_oriented = false;
          num_connected_components = 0;
        }

        /// @brief Return connected_component[isimplex].
        bool ConnectedComponent(const int isimplex) const
        { return connected_component[isimplex]; }

        /// @brief Return number of connected components.
        NUM_TYPE NumConnectedComponents() const
        { return num_connected_components; }
        
        /// @brief Return is_simplex_oriented[isimplex].
        bool IsSimplexOriented(const int isimplex)
        { return is_simplex_oriented[isimplex]; }
        
        /// @brief Return is_oriented.
        bool IsOriented() const
        { return is_oriented; }

        /// @brief Return number of simplices in table entry.
        int NumSimplices() const
        { return is_simplex_oriented.size(); }
      };

      
    protected:
      std::vector<ORIENT_SIMPLEX_ENTRY> entry;

      
    public:
      
      /// Constructor.
      template <typename _ISOTABLE_TYPE>
      MCUBE_ORIENT_SIMPLEX_TABLE
      (const _ISOTABLE_TYPE & isotable)
      {
        const int num_vertices_per_simplex =
          isotable.NumVerticesPerSimplex();

        entry.resize(isotable.NumTableEntries());

        for (int table_index = 0;
             table_index < isotable.NumTableEntries(); table_index++) {
          const ISOSURFACE_VERTEX_INDEX * simplex_vertices =
            isotable.SimplexVertices(table_index);
          const int num_simplices = isotable.NumSimplices(table_index);

          entry[table_index].connected_component.resize(num_simplices, 0);
          entry[table_index].is_simplex_oriented.resize(num_simplices, false);
          
          IJK::get_facet_connected_components_in_simplicial_complex
            (simplex_vertices, num_vertices_per_simplex, num_simplices,
             entry[table_index].connected_component,
             entry[table_index].num_connected_components);

        }
      };

      bool IsEntryOriented(const _TABLE_INDEX table_index) const
      { return entry[table_index].is_oriented; }

      bool IsSimplexOriented
      (const _TABLE_INDEX table_index, const int isimplex) const
      { return entry[table_index].is_simplex_oriented[isimplex]; }

      COMPONENT_INDEX ConnectedComponent
      (const TABLE_INDEX table_index, const int simplex_index) const
      { return entry[table_index].ConnectedComponent(simplex_index); }

      
      void SetEntryOriented
      (const _TABLE_INDEX table_index, const bool flag)
      { entry[table_index].is_oriented = true; }

      void SetSimplexOriented
      (const _TABLE_INDEX table_index, const int isimplex,
       const bool flag)
      { entry[table_index].is_simplex_oriented[isimplex] = true; }


      /// Set orientation flag of all simplices in connected component
      ///   containing isimplex to flag.
      void SetConnectedComponentOriented
      (const _TABLE_INDEX table_index, const int isimplex,
       const bool flag)
      {
        const COMPONENT_INDEX icomponent =
          ConnectedComponent(table_index, isimplex);
        const int num_simplices =
          entry[table_index].NumSimplices();

        bool flag_entry_oriented = true;
        for (int jsimplex = 0; jsimplex < num_simplices; jsimplex++) {
          const COMPONENT_INDEX jcomponent =
            ConnectedComponent(table_index, jsimplex);
          if (jcomponent == icomponent)
            { SetSimplexOriented(table_index, jsimplex, flag); }

          if (!IsSimplexOriented(table_index, jsimplex))
            { flag_entry_oriented = false; }
        }

        SetEntryOriented(table_index, flag_entry_oriented);
      }

    };

      
    template <typename _ORIENT_INFO_TYPE, typename ITYPEA, typename ITYPEB>
    bool some_boundary_facets_match
    (const _ORIENT_INFO_TYPE & orient_info,
     const TABLE_INDEX table_indexA, const TABLE_INDEX table_indexB,
     ITYPEA & ifacetA, ITYPEB & ifacetB)
    {
      typedef typename _ORIENT_INFO_TYPE::VERTEX_TYPE VERTEX_TYPE;
    
      const int num_facetsA = orient_info[table_indexA].NumBoundaryFacets();
      const int num_facetsB = orient_info[table_indexB].NumBoundaryFacets();

      // Initialize.
      ifacetA = 0;
      ifacetB = 0;
    
      for (int iA = 0; iA < num_facetsA; iA++) {

        const VERTEX_TYPE * facetA_ptr =
          orient_info.BoundaryFacetVertices(table_indexA, iA);
      
        int iB;
        if (orient_info.ContainsBoundaryFacet
            (table_indexB, facetA_ptr, iB)) {
          ifacetA = iA;
          ifacetB = iB;
        
          return true;
        }
      }

      return false;
    }


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

      const int numc = orient_info[istart].num_connected_components;

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
      typedef MCUBE_ORIENT_INFO_TABLE
        <TABLE_INDEX, ISOSURFACE_VERTEX_INDEX, short, unsigned char, short>
        ORIENT_INFO_TABLE_TYPE;
      
      const int num_table_entries = isotable.NumTableEntries();
      const int num_vert_per_simplex = isotable.NumVerticesPerSimplex();

      if (num_vert_per_simplex < 2) {
        // Nothing to orient.
        return;
      }
    
      ORIENT_INFO_TABLE_TYPE orient_info
        (num_vert_per_simplex-1, num_table_entries);
      std::stack<int> stack;
      IJK::PROCEDURE_ERROR error("orient_mcube_table");

      if (!isotable.CheckTableIndex(istart, error))
        { throw error; }

      for (TABLE_INDEX it = 0; it < isotable.NumTableEntries(); it++) {
        int numc;
        isotable.OrientAllSimplicesInTableEntry(it, numc);
        orient_info.SetNumConnectedComponents(it, numc);
      }

      if (!check_orient_simplices_starting_table_entry
          (isotable, orient_info, istart, error))
        { throw error; }
          
      for (TABLE_INDEX it = 0; it < isotable.NumTableEntries(); it++) {
        const ISOSURFACE_VERTEX_INDEX * simplex_vertices =
          isotable.SimplexVertices(it);
      
        orient_info.SetBoundaryFacets
          (it, simplex_vertices, isotable.NumSimplices(it));
      }

      int num_completed = 1;
      orient_info.SetIsOriented(istart, true);
      stack.push(istart);
      while (stack.size() > 0) {
        const int table_indexA = stack.top();
        stack.pop();

        for (int table_indexB = 0;
             table_indexB < num_table_entries; table_indexB++) {

          if (orient_info[table_indexB].is_oriented) { continue; }
          if (table_indexA == table_indexB) { continue; }
          if (isotable.NumSimplices(table_indexB) == 0) { continue; }

          int ifacetA, ifacetB;
          if (some_boundary_facets_match
              (orient_info, table_indexA, table_indexB,
               ifacetA, ifacetB)) {

            if (orient_info[table_indexB].is_boundary_facet_oriented[ifacetB])
              { continue; }
            
            const ISOSURFACE_VERTEX_INDEX * facetA_ptr =
              orient_info.BoundaryFacetVertices(table_indexA, ifacetA);


            bool flag_matched;
            std::vector<bool> is_orientedB
              (isotable.NumSimplices(table_indexB), false);

            isotable.OrientTwoTableEntries
              (table_indexA, table_indexB, flag_matched,
               is_orientedB);

            if (flag_matched) {
              
              orient_info.SetIsBoundaryFacetOriented
                (table_indexB, ifacetB, true);

              if (orient_info[table_indexB].num_connected_components == 1) {
                orient_info.SetIsOriented(table_indexB, true);
                stack.push(table_indexB);
                num_completed++;
              }
              else {
                for (int ifacetB = 0;
                     ifacetB < orient_info[table_indexB].NumBoundaryFacets();
                     ifacetB++) {
                  const int isimplexB =
                    orient_info[table_indexB].simplex_containing_boundary_facet[ifacetB];
              
                  if (is_orientedB[isimplexB]) {
                    orient_info.SetIsBoundaryFacetOriented
                      (table_indexB, ifacetB, true);
                  }
                }

                if (orient_info[table_indexB].AreAllBoundaryFacetsOriented()) {
                  orient_info.SetIsOriented(table_indexB, true);
                  num_completed++;
                }
              }

              if (flag_verbose &&
                  orient_info[table_indexB].is_oriented &&
                  num_completed%output_trigger == 0) {
                out << "  " << num_completed
                    << " out of " << isotable.NumTableEntries()
                    << " isosurface table entries oriented.\n";
                out.flush();
              }
            }
            
          }
        }
      }

      if (flag_verbose && num_completed > output_trigger) {
        out << "  Completed orientation of all " << num_completed
            << " isosurface table entries.\n";
        out.flush();
      }
    }


    /*!
     *  @brief Return true if simplices in each table entry 
     *    are consistently oriented. (Local)
     *  @param isotable Marching Cubes lookup table.
     */
    template <typename OSTREAM_TYPE, typename _ISOTABLE_TYPE,
              typename _NTYPE>
    bool check_mcube_table_consistent_orientation_local
    (OSTREAM_TYPE & out, const _ISOTABLE_TYPE & isotable, 
     const bool flag_verbose, const _NTYPE output_trigger,
     IJK::ERROR & error)
    {
      const int num_vert_per_simplex = isotable.NumVerticesPerSimplex();
      
      if (num_vert_per_simplex < 2) {
        // Nothing to check
        return true;
      }

      int isimplexA, isimplexB;
      for (TABLE_INDEX it = 0; it < isotable.NumTableEntries(); it++) {
        if (!isotable.AreSimplicesConsistentlyOriented
            (it, isimplexA, isimplexB)) {
          error.AddMessage
            ("Simplices ", isimplexA, " and ", isimplexB,
             " in table entry ", it,
             " are not consistently oriented.");
          return false;
        }
      }

      return true;
    }


    /*!
     *  @brief Return true if orientation of isimplexA in table_indexA
     *    matches orientation of isimplexB in table_indexB.
     *  @pre Simplices isimplexA and isimplexB have matching facets.
     *  @param orient_info Orientation information.
     */
    template <typename _ISOTABLE_TYPE,
              typename _ITYPEA, typename _ITYPEB,
              typename _ISIMPLEXA, typename _ISIMPLEXB>
    bool check_mcube_table_two_simplex_orientations
    (const _ISOTABLE_TYPE & isotable,
     const _ITYPEA table_indexA, const _ITYPEB table_indexB,
     _ISIMPLEXA & isimplexA, _ISIMPLEXB isimplexB,
     IJK::ERROR & error)
    {
      if (!isotable.AreTwoTableEntriesConsistentlyOriented
          (table_indexA, table_indexB, isimplexA, isimplexB)) {
        error.AddMessage
          ("  Simplex ", isimplexA, " in table entry ", table_indexA,
           " has inconsistent orientation");
        error.AddMessage
          ("  with simplex ", isimplexB, " in table entry ",
           table_indexB, ".");
        return false;
      }

      return true;
    }

    
    /*!
     *  @brief Return true if orientation of some simplex in table_indexA
     *    matches orientations of some simplex in table_indexB.
     *  @param orient_info Orientation information.
     */    
    template <typename _ISOTABLE_TYPE, typename _ORIENT_INFO,
              typename _ITYPEA, typename _ITYPEB>
    bool check_mcube_table_entry_orientation
    (const _ISOTABLE_TYPE & isotable, const _ORIENT_INFO & orient_info,
     const _ITYPEA table_indexA, const _ITYPEB table_indexB,
     IJK::ERROR & error)
    {
      int isimplexA, isimplexB;
      int ifacetA, ifacetB;
      if (some_boundary_facets_match
          (orient_info, table_indexA, table_indexB,
           ifacetA, ifacetB)) {

        if (!check_mcube_table_two_simplex_orientations
            (isotable, table_indexA, table_indexB,
             isimplexA, isimplexB, error))
          { return false; }
      }

      return true;
    }

    
    /*!
     *  @brief Return true if all simplex lists in Marching Cubes lookup table
     *    are consistently oriented. (Local)
     *  - Slow version: Checks every table entry against every other table entry.
     *  @param isotable Marching Cubes lookup table.
     */
    template <typename OSTREAM_TYPE, typename _ISOTABLE_TYPE,
              typename _NTYPE>
    bool check_mcube_table_orientation_slow_local
    (OSTREAM_TYPE & out, const _ISOTABLE_TYPE & isotable, 
     const bool flag_verbose, const _NTYPE output_trigger,
     IJK::ERROR & error)
    {
      typedef typename _ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
      typedef MCUBE_ORIENT_INFO_TABLE
        <TABLE_INDEX, ISOSURFACE_VERTEX_INDEX, short, unsigned char, short>
        ORIENT_INFO_TABLE_TYPE;
      
      const int num_table_entries = isotable.NumTableEntries();
      const int num_vert_per_simplex = isotable.NumVerticesPerSimplex();

      if (num_vert_per_simplex < 2) {
        // Nothing to check
        return true;
      }

      if (!check_mcube_table_consistent_orientation_local
          (out, isotable, flag_verbose, output_trigger, error))
        { return false; }

      int isimplexA, isimplexB;
      ORIENT_INFO_TABLE_TYPE orient_info
        (num_vert_per_simplex-1, num_table_entries);
      orient_info.Set(isotable);
      
      for (TABLE_INDEX table_indexA = 0; table_indexA < isotable.NumTableEntries();
           table_indexA++) {

        for (TABLE_INDEX table_indexB = table_indexA+1;
             table_indexB < isotable.NumTableEntries(); table_indexB++) {


          if (!check_mcube_table_entry_orientation
              (isotable, orient_info, table_indexA, table_indexB, error))
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
     *  @brief Return true if table_indexA is consistently oriented
     *    with other table entries. (Local)
     *  - Fast version: Returns true once all simplices in table_indexA
     *    have matched some previously oriented table entry.
     *  - Checks table_indexA only against oriented entries 
     *    with exactly one connected component.
     *  @param isotable Marching Cubes lookup table.
     */
    template <typename OSTREAM_TYPE, typename _ISOTABLE_TYPE,
              typename _ORIENT_INFO_TABLE_TYPE,
              typename _ORIENT_SIMPLEX_TABLE_TYPE,
              typename _NTYPE>
    bool check_mcube_table_entry_orientation_fast_local
    (OSTREAM_TYPE & out, const _ISOTABLE_TYPE & isotable,
     const TABLE_INDEX table_indexA,
     const _ORIENT_INFO_TABLE_TYPE & orient_info,
     const bool flag_verbose, const _NTYPE output_trigger,
     _ORIENT_SIMPLEX_TABLE_TYPE & orient_simplex,
     int & num_checked,
     IJK::ERROR & error)
    {
      int isimplexA, isimplexB;
      
      if (orient_simplex.IsEntryOriented(table_indexA)) {
        // Table entry table_indexA is already flagged
        //   as oriented consistently with other table entries.
        return true;
      }

      for (TABLE_INDEX table_indexB = 0;
           table_indexB < isotable.NumTableEntries(); table_indexB++) {

        if ((orient_info.NumConnectedComponents(table_indexB) == 1) &&
            (orient_simplex.IsEntryOriented(table_indexB))) {
              
          int ifacetA, ifacetB;
          if (some_boundary_facets_match
              (orient_info, table_indexA, table_indexB,
               ifacetA, ifacetB)) {

            if (!check_mcube_table_two_simplex_orientations
                (isotable, table_indexA, table_indexB,
                 isimplexA, isimplexB, error))
              { return false; };

            orient_simplex.SetConnectedComponentOriented
              (table_indexA, isimplexA, true);

            if (orient_simplex.IsEntryOriented(table_indexA)) {
              num_checked++;

              if (flag_verbose && (num_checked > 0) &&
                  num_checked%output_trigger == 0) {
                out << "  Checked " << num_checked
                    << " out of " << isotable.NumTableEntries()
                    << " isosurface table entry orientations.\n";
                out.flush();
              }

              // All simplex orientations are consistent
              //   with checked orientations.
              return true;
            }
          }
        }
      }

      return true;
    }

      
    /*!
     *  @brief Return true if all simplex lists in Marching Cubes lookup table
     *    are consistently oriented. (Local)
     *  - Fast version: Skips entries that have already matched
     *    some table entry.
     *  @param isotable Marching Cubes lookup table.
     */
    template <typename OSTREAM_TYPE, typename _ISOTABLE_TYPE,
              typename _NTYPE>
    bool check_mcube_table_orientation_fast_local
    (OSTREAM_TYPE & out, const _ISOTABLE_TYPE & isotable,
     const bool flag_verbose, const _NTYPE output_trigger,
     IJK::ERROR & error)
    {
      typedef typename _ISOTABLE_TYPE::TABLE_INDEX TABLE_INDEX;
      typedef MCUBE_ORIENT_INFO_TABLE
        <TABLE_INDEX, ISOSURFACE_VERTEX_INDEX, short, unsigned char, short>
        ORIENT_INFO_TABLE_TYPE;
      typedef MCUBE_ORIENT_SIMPLEX_TABLE<TABLE_INDEX, short, short, short, short>
        ORIENT_SIMPLEX_TABLE_TYPE;
      const int INIT_CHECK_NUM = 20;
      
      const int num_table_entries = isotable.NumTableEntries();
      const int num_vert_per_simplex = isotable.NumVerticesPerSimplex();

      if (num_vert_per_simplex < 2) {
        // Nothing to check
        return true;
      }

      if (!check_mcube_table_consistent_orientation_local
          (out, isotable, flag_verbose, output_trigger, error))
        { return false; }

      int isimplexA, isimplexB;
      ORIENT_INFO_TABLE_TYPE orient_info
        (num_vert_per_simplex-1, num_table_entries);
      orient_info.Set(isotable);
      ORIENT_SIMPLEX_TABLE_TYPE orient_simplex(isotable);

      // Check first INIT_CHECK_NUM table entries with exactly one connected component.
      int num_checked = 0;
      TABLE_INDEX table_indexA = 0;
      while ((table_indexA < isotable.NumTableEntries()) &&
             (num_checked < INIT_CHECK_NUM)) {

        if (orient_info.NumConnectedComponents(table_indexA) == 1) {
        
          for (TABLE_INDEX table_indexB = table_indexA+1;
               table_indexB < isotable.NumTableEntries(); table_indexB++) {

            if (orient_info.NumConnectedComponents(table_indexB) == 1) {
              
              if (!check_mcube_table_entry_orientation
                  (isotable, orient_info, table_indexA, table_indexB, error))
                { return false; }
            }
          }

          orient_simplex.SetEntryOriented(table_indexA, true);
          num_checked++;
        }

        table_indexA++;
      }
      

      // Next, check table entries with exactly one connected component.
      for (TABLE_INDEX table_indexA = 0; table_indexA < isotable.NumTableEntries();
           table_indexA++) {

        if ((orient_info.NumConnectedComponents(table_indexA) == 1) &&
            (orient_simplex.IsEntryOriented(table_indexA))) {
            
          for (TABLE_INDEX table_indexB = table_indexA+1;
               table_indexB < isotable.NumTableEntries(); table_indexB++) {

            if (table_indexB == table_indexA) { continue; }

            if ((orient_info.NumConnectedComponents(table_indexB) == 1) &&
                (!orient_simplex.IsEntryOriented(table_indexB))) {

              int ifacetA, ifacetB;
              if (some_boundary_facets_match
                  (orient_info, table_indexA, table_indexB,
                   ifacetA, ifacetB)) {

                if (!check_mcube_table_two_simplex_orientations
                    (isotable, table_indexA, table_indexB,
                     isimplexA, isimplexB, error))
                  { return false; }

                orient_simplex.SetEntryOriented(table_indexB, true);
                num_checked++;

                if (flag_verbose && (num_checked > 0) &&
                    num_checked%output_trigger == 0) {
                  out << "  Checked " << num_checked
                      << " out of " << isotable.NumTableEntries()
                      << " isosurface table entry orientations.\n";
                  out.flush();
                }
              }
            }
          }
        }
      }

      // Check table entries with more than one connected component.
      for (TABLE_INDEX table_indexA = 0; table_indexA < isotable.NumTableEntries();
           table_indexA++) {

        if (!orient_simplex.IsEntryOriented(table_indexA)) {

          if (!check_mcube_table_entry_orientation_fast_local
            (out, isotable, table_indexA, orient_info,
             flag_verbose, output_trigger, orient_simplex,
             num_checked, error))
            { return false; }
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
      return check_mcube_table_orientation_slow_local
        (out, isotable, flag_verbose, output_trigger, error);
    }
    else {
      return check_mcube_table_orientation_fast_local
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
  // MCUBE_TABLE_ORIENT_INFO member functions
  // *****************************************************************

  // Return true if boundary_facet_vertex_list[] contains
  //   facet facet_vert[].
  template <typename _TABLE_INDEX, typename _VTYPE,
            typename _STYPE, typename _PTYPE, typename _NTYPE>
  template <typename _VTYPEF, typename _ITYPE>
  bool MCUBE_ORIENT_INFO_TABLE<_TABLE_INDEX,_VTYPE,_STYPE,_PTYPE,_NTYPE>::
  ContainsBoundaryFacet
  (const _TABLE_INDEX table_index, const _VTYPEF * facet_vert,
   _ITYPE & ifacet) const
  {
    for (int ifacetB = 0; ifacetB < NumBoundaryFacets(table_index); ifacetB++) {

      const _VTYPEF * facetB_vert =
        BoundaryFacetVertices(table_index, ifacetB);

      if (std::equal(facet_vert, facet_vert+NumVerticesPerFacet(),
                     facetB_vert)) {
        ifacet = ifacetB;
        return true;
      }
    }

    return false;
  }
  
  
  template <typename _TABLE_INDEX, typename _VTYPE,
            typename _STYPE, typename _PTYPE, typename _NTYPE>
  void MCUBE_ORIENT_INFO_TABLE<_TABLE_INDEX,_VTYPE,_STYPE,_PTYPE,_NTYPE>::
  SetBoundaryFacets
  (const _TABLE_INDEX table_index,
   const _VTYPE * simplex_vertex_list,
   const int num_simplices)
  {
    const int num_vert_per_simplex = NumVerticesPerFacet()+1;
    
    IJK::get_simplex_boundary_facets
      (simplex_vertex_list, num_vert_per_simplex, num_simplices,
       entry[table_index].boundary_facet_vertex_list,
       entry[table_index].simplex_containing_boundary_facet,
       entry[table_index].boundary_facet_swap_parity);

    const int num_boundary_facets =
      entry[table_index].boundary_facet_swap_parity.size();
    
    entry[table_index].is_boundary_facet_oriented.resize
      (num_boundary_facets, false);

  }


  // Set ORIENT_INFO_TABLE from ISOSURFACE_TABLE.
  template <typename _TABLE_INDEX, typename _VTYPE,
            typename _STYPE, typename _PTYPE, typename _NTYPE>
  void MCUBE_ORIENT_INFO_TABLE<_TABLE_INDEX,_VTYPE,_STYPE,_PTYPE,_NTYPE>::
  Set(const ISOSURFACE_TABLE & isotable)

  {
    const int num_vertices_per_simplex = isotable.NumVerticesPerSimplex();
    num_vertices_per_facet = num_vertices_per_simplex-1;
    entry.resize(isotable.NumTableEntries());

    for (TABLE_INDEX it = 0; it < isotable.NumTableEntries(); it++) {
      const ISOSURFACE_VERTEX_INDEX * simplex_vertices =
        isotable.SimplexVertices(it);
      const int num_simplices = isotable.NumSimplices(it);
      int numc;
      std::vector<int> simplex_component;
      
      IJK::get_connected_components_in_simplicial_complex
        (simplex_vertices, num_vertices_per_simplex, num_simplices,
         simplex_component, numc);
      SetNumConnectedComponents(it, numc);
      SetBoundaryFacets(it, simplex_vertices, num_simplices);
    }

  }

  
  // *****************************************************************
  // MCUBE_ISOTABLE_ORIENT_INFO member functions
  // *****************************************************************

  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO<_ISO_VERTEX_BITSET_TYPE,_ITYPE>::
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

      entry[table_index].is_oriented.resize(num_simplices, false);
      entry[table_index].simplex_info.resize(num_simplices);
      
      _FlagVerticesInEachSimplex
        (table_index, simplex_vertices, num_simplices);
      _SetConnectedComponentIndices
        (table_index, simplex_vertices, num_simplices);
      _SetFacetSwapParity(table_index, simplex_vertices, num_simplices);
      _FlagBoundaryFacets(table_index);
    }
  }


  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO<_ISO_VERTEX_BITSET_TYPE,_ITYPE>::
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
  }
  
      
  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO<_ISO_VERTEX_BITSET_TYPE,_ITYPE>::
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


  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO<_ISO_VERTEX_BITSET_TYPE,_ITYPE>::
  _SetFacetSwapParity
  (const TABLE_INDEX table_index, 
   const ISOSURFACE_VERTEX_INDEX * simplex_vertex_list,
   const int num_simplices)
  {
    int swap_parity;
    
    std::vector<ISOSURFACE_VERTEX_INDEX>
      sorted_facet_vert(NumVerticesPerSimplex());
    
    for (int isimplex = 0; isimplex < num_simplices; isimplex++) {

      const ISOSURFACE_VERTEX_INDEX * first_simplex_vertex =
        simplex_vertex_list + isimplex*NumVerticesPerSimplex();

      entry[table_index].simplex_info[isimplex].facet_swap_parity.reset();

      for (int jloc = 0; jloc < NumFacetsPerSimplex(); jloc++) {
        const ISOSURFACE_VERTEX_INDEX jw = first_simplex_vertex[jloc];

        IJK::sort_simplex_facet_vertices
          (first_simplex_vertex, NumVerticesPerSimplex(),
           jloc, sorted_facet_vert.data(), swap_parity);

        if (swap_parity) {
          entry[table_index].simplex_info[isimplex].facet_swap_parity.set(jw);
        }
      }
    }

  };


  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  void MCUBE_ISOTABLE_ORIENT_INFO<_ISO_VERTEX_BITSET_TYPE,_ITYPE>::
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


  // Convert ISO_VERTEX_BITSET_TYPE to string.
  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  std::string MCUBE_ISOTABLE_ORIENT_INFO<_ISO_VERTEX_BITSET_TYPE,_ITYPE>::
  _ConvertBitsetToString(const ISO_VERTEX_BITSET_TYPE & bitset) const
  {
    const int num_isov = NumIsosurfaceVertices();
    
    std::string bitset_str = bitset.to_string();
    bitset_str = bitset_str.substr(bitset.size()-num_isov);
    return bitset_str;
  }
  
  
  template <typename _ISO_VERTEX_BITSET_TYPE, typename _ITYPE>
  bool MCUBE_ISOTABLE_ORIENT_INFO<_ISO_VERTEX_BITSET_TYPE,_ITYPE>::
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

