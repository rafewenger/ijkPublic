/*!
 *  @file ijksimplex.tpp
 *  @brief Template functions for manipulating simplices.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2023 Rephael Wenger

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

#ifndef _IJKSIMPLEX_TPP_
#define _IJKSIMPLEX_TPP_

#include <algorithm>
#include <stack>
#include <vector>

#include "ijk.tpp"
#include "ijkgraph.tpp"

// *** DEBUG ***
#include "ijkprint.tpp"


/// Classes and routines for storing and manipulating simplices.
namespace IJK {

  // *****************************************************************
  //! @name Sorting.
  // *****************************************************************

  ///@{
  
  /*!
   *  @brief Sort vertices of simplex facet.
   *  - Facet is all simplex vertices except iv.
   *  - Move iv to end of list of simplex vertices.
   *  @param simplex_vert[] Simplex vertices.
   *  @param num_simplex_vert Number of simplex vertices.
   *  @param iloc Location (in simplex_vert) of vertex not in facet.
   *    - Move vertex simplex_vert[iloc] to end 
   *      of simplex sorted_facet_vert[].
   *    @pre iloc in range [0..(num_simplex_vert-1)].
   *  @param[out] simplex_sorted_facet_vert[]
   *    Simplex vertices with facet vertices sorted 
   *    in increasing order followed by vertex simplex_vert[iloc].
   *    @pre simplex_sorted_facet_vert[] is pre-allocated
   *    to size num_simplex_vert.
   *  @param[out] swap_parity Parity (0=even, 1=odd) of number 
   *    of swaps to change simplex_vert[] to simplex_sorted_facet_vert[].
   */
  template <typename VTYPE_IN, typename VTYPE_OUT, typename ITYPE,
            typename NTYPEV, typename PTYPES>
  void sort_simplex_facet_vertices
  (const VTYPE_IN simplex_vert[], const NTYPEV num_simplex_vert,
   const ITYPE iloc,
   VTYPE_OUT simplex_sorted_facet_vert[],
   PTYPES & swap_parity)
  {
    swap_parity = 0;

    if (num_simplex_vert <= 1) {
      // Nothing to sort.
      return;
    }

    simplex_sorted_facet_vert[num_simplex_vert-1] =
      simplex_vert[iloc];

    swap_parity = (num_simplex_vert+1 - iloc)%2;
    
    // Insertion sort facet vertices.
    int i1 = 0;
    for (NTYPEV i0 = 0; i0 < num_simplex_vert; i0++) {
      if (i0 != iloc) {
        const int iv0 = simplex_vert[i0];
        simplex_sorted_facet_vert[i1] = iv0;

        int i2 = i1;
        while ((i2 > 0) &&
               (simplex_sorted_facet_vert[i2-1] > iv0)) {
          std:: swap(simplex_sorted_facet_vert[i2],
                     simplex_sorted_facet_vert[i2-1]);
          i2 = i2-1;
          swap_parity = (swap_parity+1)%2;
        }

        i1++;
      }
    }
  }


  /*!
   *  @brief Sort vertices of each simplex.
   *  @param[out] simplex_vertex_list[] Simplex vertices.
   *  @param num_vert_per_simplex Number of simplex vertices.
   *  @param num_simplices Number of simplices.
   */
  template <typename VTYPE, typename NTYPEV, typename NTYPES>
  void sort_simplex_vertices
  (VTYPE simplex_vertex_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices)
  {
    for (NTYPES isimplex = 0; isimplex < num_simplices; isimplex++) {
      VTYPE * simplexA_vert =
        simplex_vertex_list + isimplex*num_vert_per_simplex;

      std::sort(simplexA_vert, simplexA_vert+num_vert_per_simplex);
    }
  }


  /*!
   *  @overload
   *  @brief Sort vertices of each simplex. (C++ STL vector simplex_vertex_list[].)
   *  - Version using C++ STL vector simplex_vertex_list[].
   */
  template <typename VTYPE, typename NTYPEV>
  void sort_simplex_vertices
  (std::vector<VTYPE> & simplex_vertex_list,
   const NTYPEV num_vert_per_simplex)
  {
    const auto num_simplices = simplex_vertex_list.size()/num_vert_per_simplex;

    sort_simplex_vertices
      (simplex_vertex_list.data(), num_vert_per_simplex, num_simplices);
  }

  ///@}

  
  // *****************************************************************
  //! @name Queries.
  // *****************************************************************

  ///@{
  
  /*!
   *  @brief Return true if simplex contains facet.
   *  @param[out] jloc Location of vertex not in facet
   *    in simplex_vert[].
   */
  template <typename VTYPES, typename VTYPEF, typename ITYPE,
            typename NTYPE>
  bool does_simplex_contain_facet
  (const VTYPES simplex_vert[], const NTYPE num_vert_per_simplex,
   const VTYPEF facet_vert[],
   ITYPE & jloc)
  {
    const int ONE(1);
    std::vector<bool> flag_match(num_vert_per_simplex, false);
    int num_match = 0;

    // Initialize.
    jloc = 0;
    
    for (NTYPE i = 0; i+1 < num_vert_per_simplex; i++) {
      for (NTYPE j = 0; j < num_vert_per_simplex; j++) {
        if (facet_vert[i] == simplex_vert[j]) {
          if (!flag_match[j]) {
            num_match++;
            flag_match[j] = true;
            break;
          }
        }
      }
    }

    NTYPE num_mismatch = 0;
    if ((num_match+1) == num_vert_per_simplex) {
      for (NTYPE j = 0; j < num_vert_per_simplex; j++) {
        if (!flag_match[j]) {
          jloc = j;
          num_mismatch++;
        }
      }

      if (num_mismatch != ONE) {
        IJK::PROCEDURE_ERROR error("does_simplex_contain_facet");

        error.AddMessage
          ("Programming error. Incorrect value for num_mismatch.");
        error.AddMessage("  num_mismatch: ", num_mismatch);
        error.AddMessage("  Expected num_mismatch = 1.");
        throw error;
      }

      return true;
    }
    else {
      return false;
    }
  }


  /*!
   *  @brief Return true if some simplex in simplex_vert_list[]
   *    contains facet facet_vert[].
   *  @param[out] jsimplex Index of simplex containing containing facet.
   *  @param[out] jloc Location of vertex not in simplex facet.
   */
  template <typename VTYPES, typename VTYPEF,
            typename NTYPEV, typename NTYPES,
            typename ITYPES, typename ITYPEL>
  bool does_some_simplex_in_list_contain_facet
  (const VTYPES simplex_vert_list[],
   const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   const VTYPEF facet_vert[],
   ITYPES & jsimplex, ITYPEL & jloc)
  {
    // Initialize.
    jsimplex = 0;
    jloc = 0;

    for (NTYPES jsA = 0; jsA < num_simplices; jsA++) {
      const VTYPES * simplexA_vert =
        simplex_vert_list + jsA*num_vert_per_simplex;

      ITYPEL iloc;
      if (does_simplex_contain_facet
          (simplexA_vert, num_vert_per_simplex, facet_vert, iloc)) {
        jsimplex = jsA;
        jloc = iloc;

        return true;
      }
    }

    return false;
  }
  

  /*!
   *  @brief Return true if simplices share a facet.
   *  @pre Simplices have same number of vertices.
   *  @param[out] jlocA Location of vertex not in facet
   *    in simplexA_vert[].
   *  @param[out] jlocB Location of vertex not in facet
   *    in simplexN_vert[].
   */
  template <typename VTYPEA, typename VTYPEB, typename NTYPE,
            typename ITYPEA, typename ITYPEB>
  bool do_simplices_share_a_facet
  (const VTYPEA simplexA_vert[], const VTYPEB simplexB_vert[],
   const NTYPE num_vert_per_simplex,
   ITYPEA & jlocA, ITYPEB & jlocB)
  {
    std::vector<VTYPEA> simplexA_sorted_facet_vert(num_vert_per_simplex);
    int swapA_parity;

    // Initialize
    jlocA = 0;
    jlocB = 0;
    
    for (NTYPE iA = 0; iA < num_vert_per_simplex; iA++) {

      sort_simplex_facet_vertices
        (simplexA_vert, num_vert_per_simplex,
         iA, simplexA_sorted_facet_vert.data(), swapA_parity);

      int iB;
      if (does_simplex_contain_facet
          (simplexB_vert, num_vert_per_simplex,
           simplexA_sorted_facet_vert.data(), iB)) {

        jlocA = iA;
        jlocB = iB;

        return true;
      }
    }

    return false;
  }


  /*!
   *  @brief Return true if simplices share a facet. (Does not return shared facets.)
   *  @pre Simplices have same number of vertices.
   */
  template <typename VTYPEA, typename VTYPEB, typename NTYPE>
  bool do_simplices_share_a_facet
  (const VTYPEA simplexA_vert[], const VTYPEB simplexB_vert[],
   const NTYPE num_vert_per_simplex)
  {
    NTYPE ilocA, ilocB;
    return do_simplices_share_a_facet
      (simplexA_vert, simplexB_vert, num_vert_per_simplex, ilocA, ilocB);
  }

    
  /*!
   *  @brief Return true if simplices intersect.
   *  @pre Simplices have same number of vertices.
   *  @pre Vertices of simplex A and simplex B are sorted.
   *  @param[out] iv_shared Index of some vertex
   *    in simplexA_vert[] and simplexB_vert[].
   */
  template <typename VTYPEA, typename VTYPEB, typename VTYPEC,
            typename NTYPE>
  bool do_simplices_share_a_vertex_sorted
  (const VTYPEA simplexA_sorted_vert[],
   const VTYPEB simplexB_sorted_vert[],
   const NTYPE num_vert_per_simplex,
   VTYPEC & iv_shared)
  {
    // Initialize.
    iv_shared = 0;


    NTYPE iA = 0;
    NTYPE iB = 0;
    while ((iA < num_vert_per_simplex) &&
           (iB < num_vert_per_simplex)) {

      if (simplexA_sorted_vert[iA] < simplexB_sorted_vert[iB]) { iA++; }
      else if (simplexA_sorted_vert[iA] > simplexB_sorted_vert[iB]) { iB++; }
      else {
        // simplexA_sorted_vert[iA] == simplexB_sorted_vert[iB].
        iv_shared = simplexA_sorted_vert[iA];
        return true;
      }
    }

    return false;
  }

    
  /*!
   *  @overload
   *  @brief Return true if simplices intersect.
   *    (Does not return shared vertex.)
   *  - Version that does not return index of shared vertex.
   */
  template <typename VTYPEA, typename VTYPEB, typename NTYPE>
  bool do_simplices_share_a_vertex_sorted
  (const VTYPEA simplexA_sorted_vert[],
   const VTYPEB simplexB_sorted_vert[],
   NTYPE num_vert_per_simplex)
  {
    VTYPEA iv_shared;

    return do_simplices_share_a_vertex_sorted
      (simplexA_sorted_vert, simplexB_sorted_vert,
       num_vert_per_simplex, iv_shared);
  }

  
  /*!
   *  @brief Return number of swaps to sort simplex facet vertes.
   *  @param iloc Location of vertex not in facet
   *    in simplex_vert[].
   */
  template <typename VTYPES, typename ITYPE, typename NTYPE>
  int get_simplex_facet_swap_parity
  (const VTYPES simplex_vert[], const NTYPE num_simplex_vert,
   const ITYPE iloc)
  {
    std::vector<VTYPES> simplex_sorted_facet_vert(num_simplex_vert);
    int swap_parity;

    sort_simplex_facet_vertices
      (simplex_vert, num_simplex_vert, iloc,
       simplex_sorted_facet_vert.data(), swap_parity);

    return swap_parity;
  }


  /*!
   *  @brief Return true if simplex facet is a boundary facet.
   *  - Not optimized to be particularly fast.
   *  @param simplex_vert_list[] Array of simplex vertices.
   *    - Simplex js has vertices:
   *      (simplex_vert_list[js*num_vert_per_simplex],
   *       simplex_vert_list[js*num_vert_per_simplex+1], ...,
   *       simplex_vert_list[(js+1)*num_vert_per_simplex-1]).
   *  @param num_vert_per_simplex Number of vertices
   *    per simplex.
   *    @pre num_vert_per_simplex > 1.
   *  @param num_simplices Number of simplices.
   *  @param isimplexA Index of simplex A.
   *  @param jlocA Location of vertex not in facet of simplex A.
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES,
            typename ITYPES, typename ITYPEL>
  bool is_simplex_facet_a_boundary_facet
  (const VTYPES simplex_vert_list[],
   const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   const ITYPES isimplexA,
   const ITYPEL jlocA)
  {
    const VTYPES * simplexA_vert_ptr =
      simplex_vert_list + isimplexA*num_vert_per_simplex;
    std::vector<VTYPES> simplexA_vert(num_vert_per_simplex);
    IJK::PROCEDURE_ERROR error("is_simplex_facet_a_boundary_facet");

    if (isimplexA >= num_simplices) {
      error.AddMessage
        ("Programming error. Incorrect index value of simplexA.");
      error.AddMessage("  isimplexA: ", isimplexA);
      error.AddMessage("  num_simplices: ", num_simplices);
      error.AddMessage("  isimplexA should be less than num_simplices.");
      throw error;
    }

    if (jlocA >= num_vert_per_simplex) {
      error.AddMessage
        ("Programming error. Incorrect location of simplex vertex.");
      error.AddMessage("  jlocA: ", jlocA);
      error.AddMessage("  num_vert_per_simplex: ", num_vert_per_simplex);
      error.AddMessage("  jlocA should be less than num_vert_per_simplex.");
      throw error;
    }

    int swap_parity;
    sort_simplex_facet_vertices
      (simplexA_vert_ptr, num_vert_per_simplex, jlocA,
       simplexA_vert.data(), swap_parity);

    for (NTYPES isimplexB = 0; isimplexB < num_simplices;
         isimplexB++) {

      if (isimplexA != isimplexB) {
        const VTYPES * simplexB_vert_ptr =
          simplex_vert_list + isimplexB*num_vert_per_simplex;
        int jlocB;
        
        if (does_simplex_contain_facet
            (simplexB_vert_ptr, num_vert_per_simplex,
             simplexA_vert.data(), jlocB))
          { return false; }
      }
    }

    return true;
  }


  /*!
   *  @brief Return true if all pairs of simplices are
   *    consistently oriented.
   *  @param simplex_vert_list[] Array of simplex vertices.
   *    - Simplex js has vertices:
   *      (simplex_vert_list[js*num_vert_per_simplex],
   *       simplex_vert_list[js*num_vert_per_simplex+1], ...,
   *       simplex_vert_list[(js+1)*num_vert_per_simplex-1]).
   *  @param num_vert_per_simplex Number of vertices
   *    per simplex.
   *  @param[out] isimplexA,isimplexB Simplices A and B 
   *    share a facet but have inconsistent orientation 
   *    at that facet.
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES,
            typename ITYPEA, typename ITYPEB>
  bool are_simplices_consistently_oriented
  (const VTYPES simplex_vert_list[],
   const NTYPEV num_vert_per_simplex, const NTYPES num_simplices,
   ITYPEA & isimplexA, ITYPEB & isimplexB)
  {
    const int TWO(2);
    NTYPES ilocA, ilocB;

    // Initialize.
    isimplexA = 0;
    isimplexB = 0;
    
    if (num_vert_per_simplex < TWO) {
      // No facets. Nothing to orient.
      return true;
    }
    
    for (NTYPES jsA = 0; jsA < num_simplices; jsA++) {
      const VTYPES * simplexA_vert =
        simplex_vert_list + jsA*num_vert_per_simplex;

      for (NTYPES jsB = jsA+1; jsB < num_simplices; jsB++) {

        const VTYPES * simplexB_vert =
          simplex_vert_list + jsB*num_vert_per_simplex;

        if (do_simplices_share_a_facet
            (simplexA_vert, simplexB_vert, num_vert_per_simplex,
             ilocA, ilocB)) {

          const int swapA_parity =
            get_simplex_facet_swap_parity
            (simplexA_vert, num_vert_per_simplex, ilocA);

          const int swapB_parity =
            get_simplex_facet_swap_parity
            (simplexB_vert, num_vert_per_simplex, ilocB);

          if (swapA_parity == swapB_parity) {
            // Shared facet has same orientation in both simplices.
            isimplexA = jsA;
            isimplexB = jsB;
            return false;
          }
        }
      }
    }

    return true;
  }
  
  ///@}


  // *****************************************************************
  //! @name Orient.
  // *****************************************************************

  ///@{
  
  /*!
   *  @brief Orient simplices in list consistently with simplex istart.
   *  @pre Each simplex facet is in at most two simplices.
   *  @pre Simplicial complex is orientable.
   *  @param simplex_vert_list[] Array of simplex vertices.
   *    - Simplex js has vertices:
   *      (simplex_vert_list[js*num_vert_per_simplex],
   *       simplex_vert_list[js*num_vert_per_simplex+1], ...,
   *       simplex_vert_list[(js+1)*num_vert_per_simplex-1]).
   *  @param num_vert_per_simplex Number of vertices
   *    per simplex.
   *  @param istart Index of starting simplex.
   *    @pre istart is in range [0..(num_simplices-1)].
   *  @param[out] is_oriented[isimplex] True if simplex isimplex
   *    is oriented.
   *    - Note: orient_simplices ignores any simplices
   *      where is_oriented[isimplex] is true.
   */
  template <typename VTYPES, typename NTYPEV,
            typename NTYPES, typename ITYPE>
  void orient_simplices
  (VTYPES simplex_vert_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices, const ITYPE istart,
   std::vector<bool> & is_oriented)
  {
    const int TWO(2);
    std::stack<NTYPES> simplex_stack;
    NTYPES ilocA, ilocB;

    is_oriented.resize(num_simplices, false);
    
    if (num_vert_per_simplex < TWO) {
      // No facets. Nothing to orient.
      return;
    }

    if (num_simplices <= 0) {
      // No simplices. Nothing to orient.
      return;
    }

    is_oriented[istart] = true;
    simplex_stack.push(istart);

    while (simplex_stack.size() != 0) {

      const NTYPES jsA = simplex_stack.top();
      simplex_stack.pop();
      const VTYPES * simplexA_vert =
        simplex_vert_list + jsA*num_vert_per_simplex;

      for (NTYPES jsB = 0; jsB < num_simplices; jsB++) {

        if (!is_oriented[jsB]) {
          VTYPES * simplexB_vert =
            simplex_vert_list + jsB*num_vert_per_simplex;

          if (do_simplices_share_a_facet
              (simplexA_vert, simplexB_vert, num_vert_per_simplex,
               ilocA, ilocB)) {

            const int swapA_parity =
              get_simplex_facet_swap_parity
              (simplexA_vert, num_vert_per_simplex, ilocA);

            const int swapB_parity =
              get_simplex_facet_swap_parity
              (simplexB_vert, num_vert_per_simplex, ilocB);
          
            if (swapA_parity == swapB_parity) {
              // Shared facet has same orientation in both simplices.
              // Reorient simplexB.
              std::swap
                (simplexB_vert[num_vert_per_simplex-2],
                 simplexB_vert[num_vert_per_simplex-1]);
            }

            // Add simplexB to the stack.
            is_oriented[jsB] = true;
            simplex_stack.push(jsB);
          }
        }
      }
    }
  }

    
  /*!
   *  @brief Orient simplices in list consistently with simplex istart.
   *  @pre Each simplex facet is in at most two simplices.
   *  @pre Simplicial complex is orientable.
   *  @param simplex_vert_list[] Array of simplex vertices.
   *    - Simplex js has vertices:
   *      (simplex_vert_list[js*num_vert_per_simplex],
   *       simplex_vert_list[js*num_vert_per_simplex+1], ...,
   *       simplex_vert_list[(js+1)*num_vert_per_simplex-1]).
   *  @param num_vert_per_simplex Number of vertices
   *    per simplex.
   *  @param istart Index of starting simplex.
   *    @pre istart is in range [0..(num_simplices-1)].
   *  @param[out] num_oriented Number of simplices oriented
   *    starting from istart.
   *    (Equivalent to number of simplices in connected component
   *     containing simplex istart.)
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES,
            typename ITYPE, typename NTYPEO>
  void orient_simplices
  (VTYPES simplex_vert_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices, const ITYPE istart,
   NTYPEO & num_oriented)
  {
    const int TWO(2);
    std::vector<bool> is_oriented(num_simplices, false);

    orient_simplices
      (simplex_vert_list, num_vert_per_simplex, num_simplices,
       istart, is_oriented);

    num_oriented =
      NTYPEO(std::count(is_oriented.begin(), is_oriented.end(), true));
  }


  /*!
   *  @overload
   *  @brief Orient simplices in list consistently with simplex istart.
   *  - Version that does not return number of simplices oriented.
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES,
            typename ITYPE>
  void orient_simplices
  (VTYPES simplex_vert_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices, const ITYPE istart)
  {
    NTYPES num_oriented;

    orient_simplices(simplex_vert_list, num_vert_per_simplex,
                     num_simplices, istart, num_oriented);
    
  }


  /*!
   *  @brief Orient all simplices so that all simplices 
   *    are consistently oriented.
   *    - Orientation of each connected component is arbitrary.
   *  @pre Each simplex facet is in at most two simplices.
   *  @pre Simplicial complex is orientable.
   *  @param simplex_vert_list[] Array of simplex vertices.
   *    - Simplex js has vertices:
   *      (simplex_vert_list[js*num_vert_per_simplex],
   *       simplex_vert_list[js*num_vert_per_simplex+1], ...,
   *       simplex_vert_list[(js+1)*num_vert_per_simplex-1]).
   *  @param num_vert_per_simplex Number of vertices
   *    per simplex.
   *  @param[out] num_components Number of connected components
   *    in simplex_vert_list[].
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES,
            typename NTYPEC>
  void orient_all_simplices
  (VTYPES simplex_vert_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices, NTYPEC & num_components)
  {
    std::vector<bool> is_oriented(num_simplices, false);

    num_components = 0;
    for (NTYPES isimplex = 0; isimplex < num_simplices; isimplex++) {

      if (!is_oriented[isimplex]) {
        num_components++;
        orient_simplices
          (simplex_vert_list, num_vert_per_simplex, num_simplices,
           isimplex, is_oriented);
      }
    }
  }


  /*!
   *  @overload
   *  @brief Orient all simplices so that all simplices 
   *    are consistently oriented. (C++ STL vector simplex_vertex_list[].)
   *  - Version using C++ STL vector simplex_vertex_list[].
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPEC>
  void orient_all_simplices
  (std::vector<VTYPES> & simplex_vertex_list,
   const NTYPEV num_vert_per_simplex,
   NTYPEC & num_components)
  {
    const auto num_simplices =
      simplex_vertex_list.size()/num_vert_per_simplex;

    orient_all_simplices
      (IJK::vector2pointer(simplex_vertex_list),
       num_vert_per_simplex, num_simplices, num_components);
  }


  /*!
   *  @overload
   *  @brief Orient all simplices so that all simplices 
   *    are consistently oriented. (C++ STL vector simplex_vertex_list[].)
   *  - Version using C++ STL vector simplex_vertex_list[].
   */
  template <typename VTYPES, typename NTYPEV>
  void orient_all_simplices
  (std::vector<VTYPES> & simplex_vertex_list,
   const NTYPEV num_vert_per_simplex)
  {
    typedef typename std::vector<VTYPES>::size_type SIZE_TYPE;

    SIZE_TYPE num_components;
    orient_all_simplices
      (simplex_vertex_list, num_vert_per_simplex, num_components);
  }
  
    
  /*!
   *  @brief Orient simplices in list consistently with simplex 0.
   *  @pre Each simplex facet is in at most two simplices.
   *  @pre Simplicial complex is orientable.
   *  @param simplex_vert_list[] Array of simplex vertices.
   *    - Simplex js has vertices:
   *      (simplex_vert_list[js*num_vert_per_simplex],
   *       simplex_vert_list[js*num_vert_per_simplex+1], ...,
   *       simplex_vert_list[(js+1)*num_vert_per_simplex-1]).
   *  @param num_vert_per_simplex Number of vertices
   *    per simplex.
   *  @param istart Index of starting simplex.
   *    @pre istart is in range [0..(num_simplices-1)].
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES,
            typename NTYPEO>
  void orient_simplices_with_simplex0
  (VTYPES simplex_vert_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices, NTYPEO & num_oriented)
  {
    orient_simplices
      (simplex_vert_list, num_vert_per_simplex, num_simplices, 0,
       num_oriented);
  }


  /*!
   *  @brief Orient simplices in list to match facet orientation.
   *  - Return true if found simplex in simplex_vert[] containing
   *    facet facet_vert[].
   *  @pre Only one simplex in simplex_vert_list[] contains 
   *    facet facet_vert[].
   *  @param simplex_vert[] Array of simplex vertices.
   *    - Simplex js has vertices:
   *      (simplex_vert_list[js*num_vert_per_simplex],
   *       simplex_vert_list[js*num_vert_per_simplex+1], ...,
   *       simplex_vert_list[(js+1)*num_vert_per_simplex-1]).
   *  @param num_vert_per_simplex Number of vertices
   *    per simplex.
   *  @param facet_vert[] Facet vertices.
   *  @param facet_swap_parity Facet swap parity.
   *  @param[out] num_oriented Number of simplices oriented
   *    starting from istart.
   *    (Equivalent to number of simplices in connected component
   *     containing boundary facet from list B.)
   */
  template <typename VTYPES, typename VTYPEF,
            typename NTYPEV, typename NTYPES, typename PTYPE,
            typename NTYPEO>
  bool orient_simplices_with_facet
  (VTYPES simplex_vert_list[],
   const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   const VTYPEF facet_vert[],
   const PTYPE facet_swap_parity,
   NTYPEO & num_oriented)
  {
    const int TWO(2);

    // Initialize.
    num_oriented = 0;
    
    if (num_vert_per_simplex < TWO) {
      // No facets. Nothing to orient.
      return false;
    }

    NTYPES isimplexB;
    NTYPES ilocB;
    if (does_some_simplex_in_list_contain_facet
        (simplex_vert_list, num_vert_per_simplex, num_simplices,
         facet_vert, isimplexB, ilocB)) {

      VTYPES * simplexB_vert =
        simplex_vert_list + isimplexB*num_vert_per_simplex;
      
      const int facetB_swap_parity =
        get_simplex_facet_swap_parity
        (simplexB_vert, num_vert_per_simplex, ilocB);

      if (facet_swap_parity != facetB_swap_parity) {
        // facetB has different orientation from facet.
        // Reorient simplexB.
        std::swap
          (simplexB_vert[num_vert_per_simplex-2],
           simplexB_vert[num_vert_per_simplex-1]);
      }

      orient_simplices
        (simplex_vert_list, num_vert_per_simplex,
         num_simplices, isimplexB, num_oriented);

      return true;
    }
    else {
      return false;
    }
  }


  /*!
   *  @brief Reorient simplices based on swap parity.
   *  - Reverse orientation if swap parity is 1.
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES,
            typename PTYPE>
  void reorient_simplices
  (VTYPES simplex_vertex_list[],
   const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   const std::vector<PTYPE> & simplex_swap_parity)
  {
    if (num_vert_per_simplex < 2) {
      // Nothing to orient.
      return;
    }

    const NTYPEV ilast = num_vert_per_simplex - 1;
    
    for (NTYPES isimplex = 0; isimplex < num_simplices; isimplex++) {

      if (simplex_swap_parity[isimplex] == 1) {

        VTYPES * simplexA_ptr =
          simplex_vertex_list + isimplex*num_vert_per_simplex;
        std::swap(simplexA_ptr[ilast-1], simplexA_ptr[ilast]);
      }
    }
  }


  /*!
   *  @brief Reorient simplices based on swap parity. (C++ STL vector.)
   *  - Version using C++ STL vector for simplex_vertex_list[].
   */
  template <typename VTYPES, typename NTYPEV, typename PTYPE>
  void reorient_simplices
  (std::vector<VTYPES> & simplex_vertex_list,
   const NTYPEV num_vert_per_simplex,
   const std::vector<PTYPE> & simplex_swap_parity)
  {
    const auto num_simplices =
      simplex_vertex_list.size()/num_vert_per_simplex;

    reorient_simplices
      (simplex_vertex_list.data(), num_vert_per_simplex,
       num_simplices, simplex_swap_parity);
  }

    
  ///@}


  // *****************************************************************
  //! @name Get boundary facets.
  // *****************************************************************    

  ///@{
    
  /*!
   *  @brief Get boundary facets of simplices.
   */
  template <typename VTYPES, typename VTYPEF, typename STYPE,
            typename NTYPEV, typename NTYPES, typename PTYPE>
  void get_simplex_boundary_facets
  (const VTYPES simplex_vert_list[],
   const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   std::vector<VTYPEF> & boundary_facet_vert_list,
   std::vector<STYPE> & simplex_containing_boundary_facet,
   std::vector<PTYPE> & boundary_facet_swap_parity)
  {
    std::vector<VTYPES>
      simplexA_sorted_facet_vert(num_vert_per_simplex);

    boundary_facet_vert_list.clear();
    boundary_facet_swap_parity.clear();

    if (num_vert_per_simplex < 2) {
      // No facets.
      return;
    }
    
    const NTYPEV num_vert_per_facet = num_vert_per_simplex - 1;
    
    for (NTYPES isimplexA = 0; isimplexA < num_simplices; isimplexA++) {
      for (NTYPEV ilocA = 0; ilocA < num_vert_per_simplex; ilocA++) {

        if (is_simplex_facet_a_boundary_facet
            (simplex_vert_list, num_vert_per_simplex, num_simplices,
             isimplexA, ilocA)) {

          const VTYPES * simplexA_vert_ptr =
            simplex_vert_list + isimplexA*num_vert_per_simplex;

          int swap_parity;
          sort_simplex_facet_vertices
            (simplexA_vert_ptr, num_vert_per_simplex, ilocA,
             simplexA_sorted_facet_vert.data(), swap_parity);

          const NTYPES ifacet = boundary_facet_swap_parity.size();
          boundary_facet_vert_list.resize((ifacet+1)*num_vert_per_facet);
            
          VTYPEF * boundary_facet_ptr =
            boundary_facet_vert_list.data() + (ifacet*num_vert_per_facet);

          std::copy(simplexA_sorted_facet_vert.data(),
                    simplexA_sorted_facet_vert.data()+num_vert_per_facet,
                    boundary_facet_ptr);
          simplex_containing_boundary_facet.push_back(isimplexA);
          boundary_facet_swap_parity.push_back(swap_parity);
        }
      }
    }
  }
  

  /*!
   *  @overload
   *  @brief Get boundary facets of simplices. (C++ STL simplex_vertex_list[].)
   *  - Version using C++ STL vector simplex_vertex_list[].
   */
  template <typename VTYPES, typename VTYPEF, typename STYPE,
            typename NTYPEV, typename PTYPE>
  void get_simplex_boundary_facets
  (const std::vector<VTYPES> & simplex_vertex_list,
   const NTYPEV num_vert_per_simplex,
   std::vector<VTYPEF> & boundary_facet_vert_list,
   std::vector<STYPE> & simplex_containing_boundary_facet,
   std::vector<PTYPE> & boundary_facet_swap_parity)
  {
    const auto num_simplices =
      simplex_vertex_list.size()/num_vert_per_simplex;
    
    get_simplex_boundary_facets
      (simplex_vertex_list.data(), num_vert_per_simplex,
       num_simplices, boundary_facet_vert_list,
       simplex_containing_boundary_facet,
       boundary_facet_swap_parity);
  }

  ///@}


  // *****************************************************************
  //! @name Get connected components.
  // *****************************************************************

  /*!
   *  @brief Get connected components in simplicial complex.
   *  @param simplex_vertex_list[] Array of simplex vertices.
   *  @param num_vert_per_component Number of vertices per simplex.
   *  @param num_simplices Number of simplices.
   *  @param[out] simplex_component[js]
   *     Index of component containing simplex js.
   *    @pre Preallocated to size at least num_simplices.
   *  @param[out] num_components Number of connected components.
   */
  template <typename VTYPE,
            typename NTYPEV, typename NTYPES, typename NTYPEC,
            typename ITYPEC>
  void get_connected_components_in_simplicial_complex
  (const VTYPE simplex_vertex_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   ITYPEC simplex_component[],
   NTYPEC & num_components)
  {
    IJK::PROCEDURE_ERROR
      error("get_connected_components_in_simplicial_complex");
    
    // Initialize;
    num_components = 0;

    if (num_simplices == 0) {
      // Nothing to get.
      return;
    }

    if (!check_array_allocated
        (simplex_vertex_list, "simplex_vertex_list", error))
      { throw error; }
    if (!check_array_allocated
        (simplex_component, "simplex_component", error))
      { throw error; }
    
    const NTYPES simplex_vertex_list_length =
      num_vert_per_simplex*num_simplices;
    std::vector<VTYPE> simplex_sorted_vertex_list
      (simplex_vertex_list_length);

    std::copy(simplex_vertex_list,
              simplex_vertex_list+simplex_vertex_list_length,
              simplex_sorted_vertex_list.begin());

    sort_simplex_vertices
      (simplex_sorted_vertex_list, num_vert_per_simplex);

    // parent[js] Parent of simplex js in union-find tree.
    std::vector<NTYPES> parent(num_simplices);
    IJK::init_union_find_sets(parent.data(), num_simplices);

    for (NTYPES isimplexA = 0; isimplexA < num_simplices;
         isimplexA++) {

      for (NTYPES isimplexB = isimplexA+1;
           isimplexB < num_simplices; isimplexB++) {

        // Use simplices with sorted vertices.
        const VTYPE * simplexA_ptr =
          simplex_sorted_vertex_list.data() +
          isimplexA*num_vert_per_simplex;
        const VTYPE * simplexB_ptr =
          simplex_sorted_vertex_list.data() +
          isimplexB*num_vert_per_simplex;

        if (do_simplices_share_a_vertex_sorted
            (simplexA_ptr, simplexB_ptr, num_vert_per_simplex)) {

          // Union sets containing isimplexA and isimplexB.
          IJK::union_components
            (isimplexA, isimplexB, parent.data());

          // *** DEBUG ***
          /*
          using namespace std;
          cerr << "Union " << isimplexA << " and " << isimplexB << endl;
          cerr << "  parent[" << isimplexA << "]: " << parent[isimplexA] << endl;
          cerr << "  parent[" << isimplexB << "]: " << parent[isimplexB] << endl;
          if (parent[isimplexB] != isimplexB) {
            const NTYPES isimplexC = parent[isimplexB];
            cerr << "  parent[" << isimplexC << "]: " << parent[isimplexC] << endl;
          }
          */
        }
      }
    }

    set_connected_components_from_union_find_tree
      (parent, simplex_component, num_components);
  }


  /*!
   *  @overload
   *  @brief Get connected components in simplicial complex.
   *    (C++ STL vector simplex_component[].)
   *  - Version using C++ STL vector for simplex_component[].
   */
  template <typename VTYPE,
            typename NTYPEV, typename NTYPES, typename NTYPEC,
            typename ITYPEC>
  void get_connected_components_in_simplicial_complex
  (const VTYPE simplex_vertex_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   std::vector<ITYPEC> & simplex_component,
   NTYPEC & num_components)
  {
    simplex_component.resize(num_simplices);
    get_connected_components_in_simplicial_complex
      (simplex_vertex_list, num_vert_per_simplex, num_simplices,
       IJK::vector2pointer(simplex_component), num_components);
  }

    
  /*!
   *  @overload
   *  @brief Get connected components in simplicial complex.
   *    (C++ STL vector simplex_vertex_list[] and simplex_component[].)
   *  - Version using C++ STL vector for simplex_component[].
   *  - Version using C++ STL vector for simplex_vertex_list[].
   */
  template <typename VTYPE, typename NTYPEV, typename NTYPEC,
            typename ITYPEC>
  void get_connected_components_in_simplicial_complex
  (const std::vector<VTYPE> & simplex_vertex_list,
   const NTYPEV num_vert_per_simplex,
   std::vector<ITYPEC> & simplex_component,
   NTYPEC & num_components)
  {
    const auto num_simplices =
      simplex_vertex_list.size()/num_vert_per_simplex;

    get_connected_components_in_simplicial_complex
      (IJK::vector2pointer(simplex_vertex_list),
       num_vert_per_simplex, num_simplices,
       simplex_component, num_components);
  }    


  /*!
   *  @brief Get number of simplices in each connected component.
   *  @param simplex_vertex_list[] Array of simplex vertices.
   *  @param num_vert_per_component Number of vertices per simplex.
   *  @param num_simplices Number of simplices.
   *  @param simplex_component[js]
   *    Index of component containing simplex js.
   *  @param[out] num_simplices_in_component[i]
   *    Number of simplices in component i.
   */
  template <typename ITYPEC,
            typename NTYPES, typename NTYPES2, typename NTYPEC>
  void get_num_simplices_in_each_connected_component
  (const ITYPEC simplex_component[],
   const NTYPES num_simplices,
   std::vector<NTYPES2> & num_simplices_in_component,
   const NTYPEC num_components)
  {
    num_simplices_in_component.resize(num_components);
    
    std::fill(num_simplices_in_component.begin(),
              num_simplices_in_component.end(), 0);
      
    for (NTYPES isimplex = 0; isimplex < num_simplices; isimplex++) {
      const NTYPEC icomponent = simplex_component[isimplex];
      num_simplices_in_component[icomponent]++;
    }
  }


  /*!
   *  @overload
   *  @brief Get number of simplices in connected components.
   *    (C++ STL vector simplex_component[].)
   *  - Version using C++ STL vector simplex_component[].
   *  @param simplex_vertex_list[] Array of simplex vertices.
   *  @param num_vert_per_component Number of vertices per simplex.
   *  @param num_simplices Number of simplices.
   *  @param simplex_component[js]
   *    Index of component containing simplex js.
   *    @pre simplex_component[js] is in range [0..(num_components-1)].
   *  @param[out] num_simplices_in_component[i]
   *    Number of simplices in component i.
   *  @param num_components Number of connected components.
   */
  template <typename ITYPEC,
            typename NTYPES, typename NTYPES2, typename NTYPEC>
  void get_num_simplices_in_each_connected_component
  (const std::vector<ITYPEC> & simplex_component,
   const NTYPES num_simplices,
   std::vector<NTYPES2> & num_simplices_in_component,
   const NTYPEC num_components)
  {
    get_num_simplices_in_each_connected_component
      (IJK::vector2pointer(simplex_component), num_simplices,
       num_simplices_in_component, num_components);
  }

  
  /*!
   *  @brief Get simplices in connected component icomponent.
   *  @param simplex_vertex_list[] Array of simplex vertices.
   *  @param num_vert_per_component Number of vertices per simplex.
   *  @param num_simplices Number of simplices.
   *  @param simplex_component[js]
   *     Index of component containing simplex js.
   *  @param icomponent Index of connected component.
   *  @param[out] component_simplex_vertex_list[]
   *     List of simplex vertices of simplices 
   *     in connected component.
   */
  template <typename VTYPE, typename VTYPEC,
            typename NTYPEV, typename NTYPES,
            typename ITYPEC, typename ITYPEC2>
  void get_simplices_in_connected_component
  (const VTYPE simplex_vertex_list[],
   const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   const ITYPEC simplex_component[],
   const ITYPEC2 icomponent,
   std::vector<VTYPEC> & component_simplex_vertex_list)
  {
    component_simplex_vertex_list.clear();

    for (NTYPES isimplexA = 0; isimplexA < num_simplices;
         isimplexA++) {

      if (simplex_component[isimplexA] == icomponent) {
        const VTYPE * simplexA_ptr =
          simplex_vertex_list + isimplexA*num_vert_per_simplex;

        const auto list_length =
          component_simplex_vertex_list.size();
        component_simplex_vertex_list.resize
          (list_length + num_vert_per_simplex);

        std::copy(simplexA_ptr, simplexA_ptr+num_vert_per_simplex,
                  component_simplex_vertex_list.begin()+list_length);
      }
    }
  }


  /*!
   *  @overload
   *  @brief Get simplices in connected component icomponent. 
   *    (C++ STL vector simplex_vertex_list[] and simplex_component.)
   *  - Version using C++ STL vectors simplex_vertex_list[]
   *    and simplex_component[].
   */
  template <typename VTYPE, typename VTYPEC, typename NTYPEV,
            typename ITYPEC, typename ITYPEC2>
  void get_simplices_in_connected_component
  (const std::vector<VTYPE> & simplex_vertex_list,
   const NTYPEV num_vert_per_simplex,
   const std::vector<ITYPEC> & simplex_component,
   const ITYPEC2 icomponent,
   std::vector<VTYPEC> & component_simplex_vertex_list)
  {
    const auto num_simplices =
      simplex_vertex_list.size()/num_vert_per_simplex;
    
    get_simplices_in_connected_component
      (IJK::vector2pointer(simplex_vertex_list),
       num_vert_per_simplex, num_simplices,
       IJK::vector2pointer(simplex_component), icomponent,
       component_simplex_vertex_list);
  }


  // *****************************************************************
  //! @name Get facet connected components.
  // *****************************************************************

  /*!
   *  @brief Get facet connected components in simplicial complex.
   *  - Simplices in each component are connected across facets.
   *  @param simplex_vertex_list[] Array of simplex vertices.
   *  @param num_vert_per_component Number of vertices per simplex.
   *  @param num_simplices Number of simplices.
   *  @param[out] simplex_component[js]
   *     Index of component containing simplex js.
   *    @pre Preallocated to size at least num_simplices.
   *  @param[out] num_components Number of connected components.
   */
  template <typename VTYPE,
            typename NTYPEV, typename NTYPES, typename NTYPEC,
            typename ITYPEC>
  void get_facet_connected_components_in_simplicial_complex
  (const VTYPE simplex_vertex_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   ITYPEC simplex_component[],
   NTYPEC & num_components)
  {
    IJK::PROCEDURE_ERROR
      error("get_connected_components_in_simplicial_complex");
    
    // Initialize;
    num_components = 0;

    if (num_simplices == 0) {
      // Nothing to get.
      return;
    }

    if (!check_array_allocated
        (simplex_vertex_list, "simplex_vertex_list", error))
      { throw error; }
    if (!check_array_allocated
        (simplex_component, "simplex_component", error))
      { throw error; }
    
    const NTYPES simplex_vertex_list_length =
      num_vert_per_simplex*num_simplices;
    std::vector<VTYPE> simplex_sorted_vertex_list
      (simplex_vertex_list_length);

    std::copy(simplex_vertex_list,
              simplex_vertex_list+simplex_vertex_list_length,
              simplex_sorted_vertex_list.begin());

    sort_simplex_vertices
      (simplex_sorted_vertex_list, num_vert_per_simplex);

    // parent[js] Parent of simplex js in union-find tree.
    std::vector<NTYPES> parent(num_simplices);
    IJK::init_union_find_sets(parent.data(), num_simplices);

    for (NTYPES isimplexA = 0; isimplexA < num_simplices;
         isimplexA++) {

      for (NTYPES isimplexB = isimplexA+1;
           isimplexB < num_simplices; isimplexB++) {

        // Use simplices with sorted vertices.
        const VTYPE * simplexA_ptr =
          simplex_sorted_vertex_list.data() +
          isimplexA*num_vert_per_simplex;
        const VTYPE * simplexB_ptr =
          simplex_sorted_vertex_list.data() +
          isimplexB*num_vert_per_simplex;

        if (do_simplices_share_a_facet
            (simplexA_ptr, simplexB_ptr, num_vert_per_simplex)) {

          // Union sets containing isimplexA and isimplexB.
          IJK::union_components
            (isimplexA, isimplexB, parent.data());
        }
      }
    }

    set_connected_components_from_union_find_tree
      (parent, simplex_component, num_components);
  }


  /*!
   *  @overload
   *  @brief Get facet connected components in simplicial complex.
   *    (C++ STL vector simplex_component[].)
   *  - Simplices in each component are connected across facets.
   *  - Version using C++ STL vector for simplex_component[].
   */
  template <typename VTYPE,
            typename NTYPEV, typename NTYPES, typename NTYPEC,
            typename ITYPEC>
  void get_facet_connected_components_in_simplicial_complex
  (const VTYPE simplex_vertex_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   std::vector<ITYPEC> & simplex_component,
   NTYPEC & num_components)
  {
    simplex_component.resize(num_simplices);
    get_facet_connected_components_in_simplicial_complex
      (simplex_vertex_list, num_vert_per_simplex, num_simplices,
       IJK::vector2pointer(simplex_component), num_components);
  }

    
  /*!
   *  @overload
   *  @brief Get facet connected components in simplicial complex.
   *    (C++ STL vector simplex_vertex_list[] and simplex_component[].)
   *  - Simplices in each component are connected across facets.
   *  - Version using C++ STL vector for simplex_component[].
   *  - Version using C++ STL vector for simplex_vertex_list[].
   */
  template <typename VTYPE, typename NTYPEV, typename NTYPEC,
            typename ITYPEC>
  void get_facet_connected_components_in_simplicial_complex
  (const std::vector<VTYPE> & simplex_vertex_list,
   const NTYPEV num_vert_per_simplex,
   std::vector<ITYPEC> & simplex_component,
   NTYPEC & num_components)
  {
    const auto num_simplices =
      simplex_vertex_list.size()/num_vert_per_simplex;

    get_facet_connected_components_in_simplicial_complex
      (IJK::vector2pointer(simplex_vertex_list),
       num_vert_per_simplex, num_simplices,
       simplex_component, num_components);
  }    

  
  // *****************************************************************
  //! @name Check.
  // *****************************************************************

  ///@{
  
  /*!
   *  @brief Return false and set error message if simplices are NOT
   *    consistently oriented.
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES,
            typename ITYPEA, typename ITYPEB>
  bool check_are_simplices_consistently_oriented
  (const VTYPES simplex_vert_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   ITYPEA & isimplexA, ITYPEB & isimplexB, ERROR & error)
  {
    isimplexA = 0;
    isimplexB = 0;
    
    if (are_simplices_consistently_oriented
        (simplex_vert_list, num_vert_per_simplex, num_simplices,
         isimplexA, isimplexB)) {
      return true;
    }
    else {
      error.AddMessage
        ("Simplices have inconsistent orientations.");
      error.AddMessage("  Simplices ", isimplexA, " and ",
                       isimplexB, " share a facet");
      error.AddMessage
        ("  but have inconsistent orientations on that facet.");
      return false;
    }
  }

    
  /*!
   *  @overload
   *  @brief Return false and set error message if simplices are NOT
   *    consistently oriented.
   *  - Does not return indices of two simplices 
   *    with inconsistent orientations.
   */
  template <typename VTYPES, typename NTYPEV, typename NTYPES>
  bool check_are_simplices_consistently_oriented
  (const VTYPES simplex_vert_list[], const NTYPEV num_vert_per_simplex,
   const NTYPES num_simplices,
   ERROR & error)
  {
    VTYPES isimplexA, isimplexB;

    return check_are_simplices_consistently_oriented
      (simplex_vert_list, num_vert_per_simplex, num_simplices,
       isimplexA, isimplexB, error);
  }

  ///@}
   
}

#endif

