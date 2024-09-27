/*!
 *  @file ijkmesh.tpp
 *  @brief ijk templates for processing polyhedral meshes.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2010-2024 Rephael Wenger

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


#ifndef _IJKMESH_
#define _IJKMESH_

#include "ijk.tpp"
#include "ijklist.tpp"
#include "ijkmap.tpp"

#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>


namespace IJK {

  // **************************************************
  // FORWARD DECLARATIONS
  // **************************************************

  template <int NUM_VERT_PER_POLY, typename VTYPE>
  bool check_vertex_vector_size
  (const std::vector<VTYPE> & vlist, IJK::ERROR & error);


  // **************************************************
  /// @name FUNCTION CLASS FOR TUPLE COMPARISON
  // **************************************************

  //@{

  /// comparison function template
  template <typename DTYPE, typename T> class TUPLE_LESS_THAN {
  protected:
    const DTYPE dimension;
    const T * array;

  public:
    TUPLE_LESS_THAN(const DTYPE d, const T * a):
      dimension(d), array(a) {};
  

    bool operator ()(const int i0, const int i1) const
    { 
      const T * p0 = i0*dimension+array;
      const T * p1 = i1*dimension+array;

      return(std::lexicographical_compare(p0, p0+dimension, p1, p1+dimension)); 
    };
  };

  //@}


  // **************************************************
  /// @name COMPUTE SIMPLEX FACETS
  // **************************************************

  //@{

  /// Compute simplex facets
  template <typename DTYPE, typename SVERT, typename NTYPE,
            typename FVERT, typename FACET_INDEX>
  void compute_simplex_facets
  (const DTYPE simplex_dimension, 
   const SVERT * simplex_vertex, const NTYPE num_simplices,
   std::vector<FVERT> & facet_vertex, 
   std::vector<FACET_INDEX> & simplex_facet)
  {
    const NTYPE numv_per_simplex = simplex_dimension+1;
    const NTYPE num_facets_per_simplex = simplex_dimension+1;
    const NTYPE numv_per_facet =  numv_per_simplex-1;
    const NTYPE total_num_facets = num_simplices*num_facets_per_simplex;
    IJK::PROCEDURE_ERROR error("compute_simplex_facets");

    facet_vertex.clear();
    simplex_facet.clear();

    if (num_simplices == 0) { return; };
    if (num_facets_per_simplex < 1) { return; }

    IJK::ARRAY<FVERT> facet_vertex2(total_num_facets*numv_per_facet);

    // store facets in facet_vertex2[]
    NTYPE m = 0;
    for (NTYPE i = 0; i < num_simplices; i++) {
      const SVERT * firstv = simplex_vertex + i*numv_per_simplex;
      for (NTYPE j = 0; j < num_facets_per_simplex; j++) {
        for (NTYPE k = 0; k < numv_per_simplex; k++) {
          if (j != k) {
            facet_vertex2[m] = firstv[k];
            m++;
          }
        }
      }
    }

    if (m != total_num_facets*numv_per_facet) {
      error.AddMessage
        ("Programming error.  Incorrect number of vertices in facet_vertex2[].");
      throw error;
    }

    // Sort facet vertices
    for (NTYPE i = 0; i < total_num_facets; i++) {
      SVERT * firstv = facet_vertex2.Ptr() + i*numv_per_facet;
      std::sort(firstv, firstv+numv_per_facet);
    }

    // Sort facets
    IJK::ARRAY<NTYPE> index_sorted(total_num_facets);
    for (NTYPE i = 0; i < total_num_facets; i++) 
      { index_sorted[i] = i; }

    TUPLE_LESS_THAN<DTYPE, NTYPE> facet_less_than
      (numv_per_facet, facet_vertex2.PtrConst());
    std::sort(index_sorted.Ptr(), index_sorted.Ptr()+total_num_facets, 
              facet_less_than);

    // Count number of non-duplicate facets
    simplex_facet.resize(total_num_facets);
    NTYPE num_non_duplicates = 1;
    simplex_facet[0] = 0;
    for (NTYPE i = 1; i < total_num_facets; i++) {
      NTYPE j = index_sorted[i-1];
      NTYPE k = index_sorted[i];
      const FVERT * firstv_j = facet_vertex2.PtrConst()+j*numv_per_facet;
      const FVERT * firstv_k = facet_vertex2.PtrConst()+k*numv_per_facet;

      if (!std::equal(firstv_j, firstv_j+numv_per_facet, firstv_k))
        { num_non_duplicates++; };

      simplex_facet[k] = num_non_duplicates-1;
    }

    facet_vertex.resize(num_non_duplicates*numv_per_facet);

    // Copy first facet from facet_vertex2[] into facet_vertex[]
    FVERT * fvert = &(facet_vertex[0]);
    FVERT * fvert2 = facet_vertex2.Ptr()+index_sorted[0]*numv_per_facet;
    std::copy(fvert2, fvert2+numv_per_facet, fvert);
    fvert += numv_per_facet;

    // Copy facet_vertex2[] into facet_vertex[], avoiding duplicates
    for (NTYPE i = 1; i < total_num_facets; i++) {
      NTYPE j = index_sorted[i-1];
      NTYPE k = index_sorted[i];
      if (simplex_facet[k] != simplex_facet[j]) {
        // Copy simplex_facet[k] into facet_vertex[]
        SVERT * fvert2 = facet_vertex2.Ptr() + k*numv_per_facet;
        std::copy(fvert2, fvert2+numv_per_facet, fvert);

        fvert += numv_per_facet;
      }
    }

    if (fvert != &(facet_vertex[0]) + num_non_duplicates*numv_per_facet) {
      error.AddMessage
        ("Programming error.  Incorrect number of vertices in facet_vertex[].");
      throw error;
    }
  }

  //@}


  // **************************************************
  /// @name ADD POLYGON VERTICES
  // **************************************************

  //@{

  /// Add triangle vertices to triangle vertex list.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPEB>
  void add_triangle_vertices
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   std::vector<VTYPEB> & tri_vert_list)
  {
    typedef typename std::vector<VTYPEB>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_TRIANGLE(3);

    const SIZE_TYPE k = tri_vert_list.size();
    tri_vert_list.resize(k+NUM_VERT_PER_TRIANGLE);
    tri_vert_list[k] = v0;
    tri_vert_list[k+1] = v1;
    tri_vert_list[k+2] = v2;
  }

  /// Add triangle vertices to triangle vertex list.
  /// - Version with triangle vertices in array tri_vert.
  template <typename VTYPEA, typename VTYPEB>
  void add_triangle_vertices
  (const VTYPEA tri_vert[],
   std::vector<VTYPEB> & tri_vert_list)
  {
    add_triangle_vertices
      (tri_vert[0], tri_vert[1], tri_vert[2], tri_vert_list);
  }

  /// Add triangle vertices to triangle vertex list.
  /// - Version that includes flag_reverse_orient.
  /// @param flag_reverse_orient If true, add quadrilateral vertices
  ///   in reverse orientation.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPEB>
  void add_triangle_vertices
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2,
   const bool flag_reverse_orient,
   std::vector<VTYPEB> & tri_vert_list)
  {
    if (flag_reverse_orient) {
      add_triangle_vertices(v2, v1, v0, tri_vert_list);
    }
    else {
      add_triangle_vertices(v0, v1, v2, tri_vert_list);
    }
  }


  /// Add triangle vertices to triangle vertex list.
  /// - Version with triangle vertices in array tri_vert.
  /// - Version that includes flag_reverse_orient.
  /// @param flag_reverse_orient If true, add quadrilateral vertices
  ///   in reverse orientation.
  template <typename VTYPEA, typename VTYPEB>
  void add_triangle_vertices
  (const VTYPEA tri_vert[],
   const bool flag_reverse_orient,
   std::vector<VTYPEB> & tri_vert_list)
  {
    add_triangle_vertices
      (tri_vert[0], tri_vert[1], tri_vert[2], flag_reverse_orient, 
       tri_vert_list);
  }

  /// Add quadrilateral vertices to quadrilateral vertex list.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPEB>
  void add_quad_vertices
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3,
   std::vector<VTYPEB> & quad_vert_list)
  {
    typedef typename std::vector<VTYPEB>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD(4);

    const SIZE_TYPE k = quad_vert_list.size();
    quad_vert_list.resize(k+NUM_VERT_PER_QUAD);
    quad_vert_list[k] = v0;
    quad_vert_list[k+1] = v1;
    quad_vert_list[k+2] = v2;
    quad_vert_list[k+3] = v3;

  }

  /// Add quadrilateral vertices to quadrilateral vertex list.
  /// - Version that includes flag_reverse_orient.
  /// @param flag_reverse_orient If true, add quadrilateral vertices
  ///   in reverse orientation.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPEB>
  void add_quad_vertices
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3,
   const bool flag_reverse_orient, std::vector<VTYPEB> & quad_vert_list)
  {
    if (flag_reverse_orient) {
      add_quad_vertices(v3, v2, v1, v0, quad_vert_list);
    }
    else {
      add_quad_vertices(v0, v1, v2, v3, quad_vert_list);
    }
  }

  /// Add quadrilateral vertices to quadrilateral vertex list.
  /// - Version with quad vertices in array quad_vert.
  template <typename VTYPEA, typename VTYPEB>
  void add_quad_vertices
  (const VTYPEA quad_vert[],
   std::vector<VTYPEB> & quad_vert_list)
  {
    add_quad_vertices
      (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3], quad_vert_list);
  }


  /// Add quadrilateral vertices to quadrilateral vertex list.
  /// - Version with triangle vertices in array tri_vert.
  /// - Version that includes flag_reverse_orient.
  /// @param flag_reverse_orient If true, add quadrilateral vertices
  ///   in reverse orientation.
  template <typename VTYPEA, typename VTYPEB>
  void add_quad_vertices
  (const VTYPEA quad_vert[],
   const bool flag_reverse_orient,
   std::vector<VTYPEB> & quad_vert_list)
  {
    add_quad_vertices
      (quad_vert[0], quad_vert[1], quad_vert[2], quad_vert[3], 
       flag_reverse_orient, quad_vert_list);
  }


  /// Add pentagon vertices to pentagon vertex list.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPEB>
  void add_pentagon_vertices
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3,
   const VTYPE4 v4, std::vector<VTYPEB> & pentagon_vert_list)
  {
    typedef typename std::vector<VTYPEB>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_PENTAGON(5);

    const SIZE_TYPE k = pentagon_vert_list.size();
    pentagon_vert_list.resize(k+NUM_VERT_PER_PENTAGON);
    pentagon_vert_list[k] = v0;
    pentagon_vert_list[k+1] = v1;
    pentagon_vert_list[k+2] = v2;
    pentagon_vert_list[k+3] = v3;
    pentagon_vert_list[k+4] = v4;

  }

  /// Add pentagon vertices to pentagon vertex list.
  /// - Version that includes flag_reverse_orient.
  /// @param flag_reverse_orient If true, add pentagon vertices
  ///   in reverse orientation.
  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPE4, typename VTYPEB>
  void add_pentagon_vertices
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, 
   const VTYPE3 v3, const VTYPE4 v4,
   const bool flag_reverse_orient, std::vector<VTYPEB> & pentagon_vert_list)
  {
    if (flag_reverse_orient) {
      add_pentagon_vertices(v4, v3, v2, v1, v0, pentagon_vert_list);
    }
    else {
      add_pentagon_vertices(v0, v1, v2, v3, v4, pentagon_vert_list);
    }
  }

  /// Add pentagon vertices to pentagon vertex list.
  /// - Version that includes flag_reverse_orient.
  /// - Version with pentagon vertices in array pentagon_vert.
  /// @param flag_reverse_orient If true, add pentagon vertices
  ///   in reverse orientation.
  template <typename VTYPE, typename VTYPEB>
  void add_pentagon_vertices
  (const VTYPE pentagon_vert[],
   const bool flag_reverse_orient, std::vector<VTYPEB> & pentagon_vert_list)
  {
    add_pentagon_vertices
      (pentagon_vert[0], pentagon_vert[1], pentagon_vert[2],
       pentagon_vert[3], pentagon_vert[4], flag_reverse_orient,
       pentagon_vert_list);
  }

  //@}


  // **************************************************
  /// @name ADD POLYTOPE VERTICES
  // **************************************************

  //@{

  template <typename VTYPE0, typename VTYPE1, typename VTYPE2,
            typename VTYPE3, typename VTYPEB>
  void add_tetrahedron_vertices
  (const VTYPE0 v0, const VTYPE1 v1, const VTYPE2 v2, const VTYPE3 v3,
   std::vector<VTYPEB> & tet_vert_list)
  {
    typedef typename std::vector<VTYPEB>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_TETRAHEDRON(4);

    const SIZE_TYPE k = tet_vert_list.size();
    tet_vert_list.resize(k+NUM_VERT_PER_TETRAHEDRON);
    tet_vert_list[k] = v0;
    tet_vert_list[k+1] = v1;
    tet_vert_list[k+2] = v2;
    tet_vert_list[k+3] = v3;
  }

  //@}

  // **************************************************
  /// @name DELETE ISOLATED VERTICES
  // **************************************************

  //@{

  /// Set flag[x] to true for every x in list.
  template <bool bool_value, typename ETYPE, typename NTYPE>
  void set_bool_flag
  (const ETYPE * list, const NTYPE list_length, bool * flag)
  {
    for (NTYPE i = 0; i < list_length; i++)
      { flag[list[i]] = bool_value; }
  }


  /*!
   *  @brief Copy coord0 to coord1, removing coordinates for any vertices iv 
   *    where flag[iv] is false.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE0, typename NTYPE1, typename ITYPE>
  void compact_coord_array
  (const DTYPE dimension, const CTYPE0 * coord0, const NTYPE0 numv0,
   const bool * flag, ITYPE * new_index, CTYPE1 * coord1, NTYPE1 & numv1)
  {
    NTYPE0 k = 0;
    for (NTYPE0 iv = 0; iv < numv0; iv++) {
      if (flag[iv]) {
        if (k < iv) {
          const CTYPE0 * coord0_ptr = coord0+dimension*iv;
          CTYPE1 * coord1_ptr = coord1+dimension*k;
          std::copy(coord0_ptr, coord0_ptr+dimension, coord1_ptr);
        }
        new_index[iv] = k;
        k++;
      }
    }

    numv1 = k;
  }


  /*!
   *  @brief Delete vertex coordinates of vertices that are not in vlist.
   *  - Compact coord0[] and relabel remaining vertices.
   *  @param coord0[] Array of vertex coordinates.
   *  @param[out] coord1[] Compacted array of vertex coordinates.
   *     - Unreferenced vertices have been removed.
   *     - coord0[] and coord1[] could be the same array.
   *  @pre Array coord1[] is pre-allocated to length at least numv0*dimension.
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE0, typename NTYPE1,
            typename VTYPE, typename NUMV_TYPE>
  void delete_unreferenced_vertices
  (const DTYPE dimension, const CTYPE0 * coord0, const NTYPE0 numv0,
   VTYPE * vlist, const NUMV_TYPE vlist_length,
   CTYPE1 * coord1, NTYPE1 & numv1)
  {
    // flag_referenced[i] = true if vertex i is referenced.
    IJK::ARRAY<bool> flag_referenced(numv0, false); 

    // new_index[i] = new index of vertex i.
    IJK::ARRAY<VTYPE> new_index(numv0);

    set_bool_flag<true>(vlist, vlist_length, flag_referenced.Ptr());

    compact_coord_array
      (dimension, coord0, numv0, flag_referenced.PtrConst(), new_index.Ptr(),
       coord1, numv1);

    for (NUMV_TYPE j = 0; j < vlist_length; j++)
      { vlist[j] = new_index[vlist[j]]; }
  }


  /*!
   *  @overload
   *  @brief Delete vertex coordinates of vertices that are not in vlist.
   *  - Compact coord[] and relabel remaining vertices.
   *  - Version that alters coord[] and numv.
   *  @param coord0[] Array of vertex coordinates.
  */
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPE, typename NUMV_TYPE>
  void delete_unreferenced_vertices
  (const DTYPE dimension, CTYPE * coord, NTYPE & numv,
   VTYPE * vlist, const NUMV_TYPE vlist_length)
  {
    NTYPE numv1;
    delete_unreferenced_vertices
      (dimension, coord, numv, vlist, vlist_length, coord, numv1);
    numv = numv1;
  }


  /*!
   *  @brief Delete vertex coordinates of vertices that are not 
   *    in vlistA or vlistB.
   *  - Compact coord0[] and relabel remaining vertices.
   *  @param coord0[] Array of vertex coordinates.
   *  @param[out] coord1[] Compacted array of vertex coordinates.
   *       Unreferenced vertices have been deleted.
   *       coord0[] and coord1[] could be the same array.
   *  @pre Array coord1[] is pre-allocated to length at least numv0*dimension.
   *  @param[out] new_index[] Array indicating new index of vertices.
   *       Vertex i is mapped to vertex new_index[i].
   *  @pre Array new_index[] is pre-allocated to length at least numv0.
   *  @param[out] flag_keep[] Array indicating vertices which are kept.
   *  @pre Array flag_keep[] is pre-allocated to length at least numv0.
  */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE0, typename NTYPE1,
            typename VTYPEA, typename VTYPEB,
            typename NUMV_TYPEA, typename NUMV_TYPEB, typename ITYPE>
  void delete_unreferenced_vertices_two_lists
  (const DTYPE dimension, const CTYPE0 * coord0, const NTYPE0 numv0,
   VTYPEA * vlistA, const NUMV_TYPEA vlistA_length,
   VTYPEB * vlistB, const NUMV_TYPEB vlistB_length,
   CTYPE1 * coord1, NTYPE1 & numv1, 
   ITYPE * new_index, bool * flag_keep)
  {
    for (NTYPE0 i = 0; i < numv0; i++)
      { flag_keep[i] = false; }

    set_bool_flag<true>(vlistA, vlistA_length, flag_keep);
    set_bool_flag<true>(vlistB, vlistB_length, flag_keep);

    compact_coord_array
      (dimension, coord0, numv0, flag_keep, new_index, coord1, numv1);

    for (NUMV_TYPEA j = 0; j < vlistA_length; j++)
      { vlistA[j] = new_index[vlistA[j]]; }

    for (NUMV_TYPEB j = 0; j < vlistB_length; j++)
      { vlistB[j] = new_index[vlistB[j]]; }
  }


  /*!
   *  @brief Delete vertex coordinates of vertices which are not in vlistA or vlistB.
   *  - Compact coord0[] and relabel remaining vertices.
   *  @param coord0[] Array of vertex coordinates.
   *  @param[out] coord1[] Compacted array of vertex coordinates.
   *       Unreferenced vertices have been deleted.
   *       coord0[] and coord1[] could be the same array.
   *  @pre Array coord1[] is pre-allocated to length at least numv0*dimension.
  */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename NTYPE0, typename NTYPE1,
            typename VTYPEA, typename VTYPEB,
            typename NUMV_TYPEA, typename NUMV_TYPEB>
  void delete_unreferenced_vertices_two_lists
  (const DTYPE dimension, const CTYPE0 * coord0, const NTYPE0 numv0,
   VTYPEA * vlistA, const NUMV_TYPEA vlistA_length,
   VTYPEB * vlistB, const NUMV_TYPEB vlistB_length,
   CTYPE1 * coord1, NTYPE1 & numv1)
  {
    // new_index[i] = new index of vertex i.
    IJK::ARRAY<NTYPE0> new_index(numv0);

    // flag_keep[i] = true if vertex i is not deleted.
    IJK::ARRAY<bool> flag_keep(numv0);

    delete_unreferenced_vertices_two_lists
      (dimension, coord0, numv0, vlistA, vlistA_length,
       vlistB, vlistB_length, coord1, numv1, 
       new_index.Ptr(), flag_keep.Ptr());
  }


  /*!
   *  @overload
   *  @brief Delete vertex coordinates of vertices that are not 
   *    in vlistA or vlistB.
   *  - Compact coord[] and relabel remaining vertices.
   *  - Version that alters coord[] and numv.
   *  @param coord0[] Array of vertex coordinates.
  */
  template <typename DTYPE, typename CTYPE, typename NTYPE,
            typename VTYPEA, typename VTYPEB,
            typename NUMV_TYPEA, typename NUMV_TYPEB,
            typename ITYPE>
  void delete_unreferenced_vertices_two_lists
  (const DTYPE dimension, CTYPE * coord, NTYPE & numv,
   VTYPEA * vlistA, const NUMV_TYPEA vlistA_length,
   VTYPEB * vlistB, const NUMV_TYPEB vlistB_length,
   ITYPE * new_index, bool * flag_keep)
  {
    NTYPE numv1;
    delete_unreferenced_vertices_two_lists
      (dimension, coord, numv, vlistA, vlistA_length, vlistB, vlistB_length,
       coord, numv1, new_index, flag_keep);
    numv = numv1;
  }


  /*!
   *  @overload
   *  @brief Delete vertex coordinates of vertices that are not 
   *  in vlistA or vlistB.
   *  C++ STL vector format for coord, vlistA and vlistB
  */
  template <typename DTYPE, typename CTYPE,
            typename VTYPEA, typename VTYPEB, typename ITYPE>
  void delete_unreferenced_vertices_two_lists
  (const DTYPE dimension, std::vector<CTYPE> & coord,
   std::vector<VTYPEA> & vlistA, std::vector<VTYPEB> & vlistB,
   std::vector<ITYPE> & new_index, bool * flag_keep)
  {
    typedef typename std::vector<CTYPE>::size_type SIZE_TYPE;

    if (coord.size() == 0) { return; }

    SIZE_TYPE numv = coord.size()/dimension;
    new_index.resize(numv);

    VTYPEA * vlistA_ptr = NULL;
    VTYPEB * vlistB_ptr = NULL;

    if (vlistA.size() > 0) { vlistA_ptr = &(vlistA.front()); }
    if (vlistB.size() > 0) { vlistB_ptr = &(vlistB.front()); }

    delete_unreferenced_vertices_two_lists
      (dimension, &(coord.front()), numv, 
       vlistA_ptr, vlistA.size(), vlistB_ptr, vlistB.size(),
       &(new_index.front()), flag_keep);

    coord.resize(numv*dimension);
  }


  /*!
   *  @overolad
   *  @brief Delete vertex coordinates of vertices that are not 
   *  in vlistA or vlistB. (C++ vector.)
   *  - C++ STL vector format for coord, vlistA and vlistB.
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPEA, typename VTYPEB, typename ITYPE>
  void delete_unreferenced_vertices_two_lists
  (const DTYPE dimension, std::vector<CTYPE> & coord,
   std::vector<VTYPEA> & vlistA, std::vector<VTYPEB> & vlistB,
   std::vector<ITYPE> & new_index, std::vector<bool> & flag_keep)
  {
    typedef typename std::vector<CTYPE>::size_type SIZE_TYPE;

    if (coord.size() == 0) { return; }
    SIZE_TYPE numv = coord.size()/dimension;

    // Create an array of bool since flag_keep may store bool values as bits.
    IJK::ARRAY<bool> flag_keep_array(numv);

    delete_unreferenced_vertices_two_lists
      (dimension, coord, vlistA, vlistB, new_index, flag_keep_array.Ptr());

    // Copy flag_keep_array into flag_keep.
    flag_keep.resize(numv);
    for (SIZE_TYPE i = 0; i < numv; i++) 
      { flag_keep[i] = flag_keep_array[i]; }
  }


  /*!
   *  @overload
   *  @brief Delete vertex coordinates of vertices that are not 
   *    in vlistA or vlistB.
   *  - C++ STL vector format for coord, vlistA and vlistB
   */
  template <typename DTYPE, typename CTYPE,
            typename VTYPEA, typename VTYPEB>
  void delete_unreferenced_vertices_two_lists
  (const DTYPE dimension, std::vector<CTYPE> & coord,
   std::vector<VTYPEA> & vlistA, std::vector<VTYPEB> & vlistB)
  {
    typedef typename std::vector<CTYPE>::size_type SIZE_TYPE;
    const SIZE_TYPE numv = coord.size()/dimension;
    std::vector<SIZE_TYPE> new_index(numv);
    IJK::ARRAY<bool> flag_keep(numv);

    delete_unreferenced_vertices_two_lists
      (dimension, coord, vlistA, vlistB, 
       new_index, flag_keep.Ptr());
  }

  //@}


  // **************************************************
  /// @name COLLAPSE ROUTINES
  // **************************************************

  //@{

  /// Set all vertices of quad to quad_vert[0].
  template <typename VTYPE>
  void collapse_quad_to_vertex(VTYPE quad_vert[])
  {
    quad_vert[1] = quad_vert[0];
    quad_vert[2] = quad_vert[0];
    quad_vert[3] = quad_vert[0];
  }

  
  /// Set all vertices of quad iquad to first vertex in quad.
  template <typename VTYPE, typename ITYPE>
  void collapse_quad_to_vertex
  (std::vector<VTYPE> & quad_vert, const ITYPE iquad)
  {
    const int NUM_VERT_PER_QUAD(4);

    collapse_quad_to_vertex
      (IJK::vector2pointer(quad_vert)+iquad*NUM_VERT_PER_QUAD);
  }

  
  /// Return true if all quad vertices are identical.
  template <typename VTYPE>
  bool is_quad_collapsed_to_one_vertex(const VTYPE quad_vert[])
  {
    if (quad_vert[0] == quad_vert[1]) {
      if (quad_vert[0] == quad_vert[2]) {
        if (quad_vert[0] == quad_vert[3])
          { return true; }
      }
    }

    return false;
  }


  /*!
   *  @overload
   *  @brief Return true if all quad vertices are identical. (C++ vector.)
   *  - Version using C++ STL vector for quad_vert.
   */
  template <typename VTYPE, typename ITYPE>
  bool is_quad_collapsed_to_one_vertex
  (std::vector<VTYPE> & quad_vert, const ITYPE iquad)
  {
    const int NUM_VERT_PER_QUAD(4);

    return is_collapsed_to_one_vertex
      (IJK::vector2pointer(quad_vert)+iquad*NUM_VERT_PER_QUAD);

  }

  //@}


  // **************************************************
  /// @name DELETE POLYTOPES
  // **************************************************

  //@{

  /*!
   *  @brief Delete polytopes with few vertices.
   *  @param min_num_poly_vert Delete all polytopes with fewer
   *    than min_num_poly_vert vertices.
   *  @param vlist[] List of polytope vertices.
   *  @param num_poly_vert[ipoly] Number of vertices in polytope ipoly.
   *  @param first_poly_vert[ipoly] Location of first vertex
   *    in polytope ipoly.
   *  @param num_poly[ipoly] Number of vertices in polytope ipoly.
   */
  template <typename NUMV_TYPE0, typename NUMV_TYPE1,
            typename VTYPE, typename ITYPE, typename NTYPE>
  void delete_poly_few_vert
  (const NUMV_TYPE0 min_num_poly_vert, VTYPE * vlist, 
   NUMV_TYPE1 * num_poly_vert, ITYPE * first_poly_vert, NTYPE & num_poly)
  {
    ITYPE vlist_length = 0;
    NTYPE k = 0;
    for (NTYPE ipoly = 0; ipoly < num_poly; ipoly++) {
      if (num_poly_vert[ipoly] >= min_num_poly_vert) {
        if (k < ipoly) {
          VTYPE * vlist_ptr0 = vlist + first_poly_vert[ipoly];
          VTYPE * vlist_ptr1 = vlist + vlist_length;
          std::copy(vlist_ptr0, vlist_ptr0+num_poly_vert[ipoly], vlist_ptr1);
          num_poly_vert[k] = num_poly_vert[ipoly];
          first_poly_vert[k] = vlist_length;
        }
        vlist_length += num_poly_vert[k];
        k++;
      }
    }
    num_poly = k;
  }


  /*!
   *  @brief Delete polytopes flagged for deletion.
   *  - Returns number of polytopes remaining.
   *  - Removes groups of polytope vertices corresponding 
   *    to deleted polytopes  from vlist[].
   *  @pre Polytopes all have NUM_VERT_PER_POLY vertices.
   *  @param vlist[] Array of polytope vertices.
   *  @param num_poly Number of polytopes.
   *  @param flag_delete[ipoly] True if ipoly is to be deleted.
   */
  template <int NUM_VERT_PER_POLY, typename VTYPE, typename NTYPE>
  NTYPE delete_poly_fixed_numv
  (VTYPE vlist[], const NTYPE num_poly, 
   const std::vector<bool> & flag_delete)
  {
    IJK::PROCEDURE_ERROR error("delete_poly_fixed_numv");

    if (num_poly != flag_delete.size()) {
      error.AddMessage("Programming error. Incorrect size of C++ vector flag_delete[].");
      error.AddMessage("  flag_delete.size() = ", flag_delete.size(), ".");
      error.AddMessage("  Should equal ", num_poly, ", the number of polytopes.");
      throw error;
    }

    if (num_poly > 0) {
      if (vlist == NULL) {
        error.AddMessage("Programming error. Array vlist[] not allocated.");
        throw error;
      }
    }

    NTYPE new_num_poly = 0;    

    for (NTYPE ipoly = 0; ipoly < num_poly; ipoly++) {

      if (!flag_delete[ipoly]) {
        const VTYPE * poly_ptr = vlist + ipoly*NUM_VERT_PER_POLY;
        std::copy(poly_ptr, poly_ptr + NUM_VERT_PER_POLY,
                  vlist + new_num_poly*NUM_VERT_PER_POLY);
        new_num_poly++;
      }
    }

    return(new_num_poly);
  }


  /*!
   *  @overload
   *  @brief Delete polytopes flagged for deletion.
   *  - Returns number of polytopes remaining.
   *  - Version using C++ STL vector for array vlist[].
   *  @pre Polytopes all have NUM_VERT_PER_POLY vertices.
   *  @param[out] vlist[] List of polytope vertices.
   *    - vlist[] is resized to new_num_poly*NUM_VERT_PER_POLY.
   *  @param flag_delete[ipoly] True if ipoly is to be deleted.
   */
  template <int NUM_VERT_PER_POLY, typename VTYPE>
  typename std::vector<VTYPE>::size_type delete_poly_fixed_numv
  (std::vector<VTYPE> & vlist,
   const std::vector<bool> & flag_delete)
  {
    typedef typename std::vector<VTYPE>::size_type NUMBER_TYPE;

    IJK::PROCEDURE_ERROR error("delete_poly_fixed_numv");

    if (!check_vertex_vector_size<NUM_VERT_PER_POLY>(vlist, error))
      { throw error; }

    const NUMBER_TYPE num_poly = vlist.size()/NUM_VERT_PER_POLY;
    const NUMBER_TYPE new_num_poly = 
      delete_poly_fixed_numv<NUM_VERT_PER_POLY>
      (IJK::vector2pointer(vlist), num_poly, flag_delete);

    vlist.resize(new_num_poly*NUM_VERT_PER_POLY);

    return(new_num_poly);
  }


  /*!
   *  @brief Delete quadrilaterals flagged for deletion.
   *  - Returns number of quadrilaterals remaining.
   *  @param[out] vlist[] List of quad vertices.
   *    - vlist[] is resized to new_num_poly*4.
   */
  template <typename VTYPE>
  typename std::vector<VTYPE>::size_type
  delete_quads(std::vector<VTYPE> & vlist, 
               const std::vector<bool> & flag_delete)
  {
    const int NUM_VERT_PER_QUAD(4);

    return(delete_poly_fixed_numv<NUM_VERT_PER_QUAD>
           (vlist, flag_delete));
  }


  /*!
   *  @brief Delete quadrilaterals that have been collapsed
   *    to a single vertex.
   *  @param vlist[] Array of quadrilateral vertices.
   *  @param num_quad Number of quadrilaterals.
   */
  template <typename VTYPE, typename NTYPE>
  NTYPE delete_collapsed_quads
  (VTYPE vlist[], const NTYPE num_quad)
  {
    const int NUM_VERT_PER_QUAD(4);
    IJK::PROCEDURE_ERROR error("delete_collapsed_quads");

    if (num_quad > 0) {
      if (vlist == NULL) {
        error.AddMessage("Programming error. Array vlist[] not allocated.");
        throw error;
      }
    }

    NTYPE new_num_quad = 0;

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      const VTYPE * quad_ptr = vlist + iquad*NUM_VERT_PER_QUAD;

      if (!is_quad_collapsed_to_one_vertex(quad_ptr)) {
        std::copy(quad_ptr, quad_ptr + NUM_VERT_PER_QUAD,
                  vlist + new_num_quad*NUM_VERT_PER_QUAD);
        new_num_quad++;
      }
    }

    return(new_num_quad);
  }

  
  /*!
   *  @overload
   *  @brief Delete quadrilaterals that have been collapsed
   *    to a single vertex.
   *  - Returns number of quadrilaterals remaining.
   */
  template <typename VTYPE>
  typename std::vector<VTYPE>::size_type
  delete_collapsed_quads(std::vector<VTYPE> & vlist)
  {
    typedef typename std::vector<VTYPE>::size_type NUMBER_TYPE;

    const int NUM_VERT_PER_QUAD(4);

    const NUMBER_TYPE num_quad = vlist.size()/NUM_VERT_PER_QUAD;
    const NUMBER_TYPE new_num_quad = 
      delete_collapsed_quads(IJK::vector2pointer(vlist), num_quad);

    vlist.resize(new_num_quad*NUM_VERT_PER_QUAD);

    return(new_num_quad);
  }

  //@}


  // **************************************************
  /// @name DELETE DUPLICATE POLYTOPE VERTICES
  // **************************************************

  //@{

  /// @brief Delete duplicate polytope vertices.
  /// *** ORDER OF vlist and num_poly_vert INCONSISTENT WITH IO ROUTINES ***
  template <typename NUMV_TYPE, typename VTYPE, typename ITYPE, 
            typename NTYPE>
  void delete_poly_vert_duplicate
  (VTYPE * vlist, NUMV_TYPE * num_poly_vert, ITYPE * first_poly_vert, 
   const NTYPE num_poly)
  {
    ITYPE vlist_length = 0;
    for (NTYPE ipoly = 0; ipoly < num_poly; ipoly++) {
      if (num_poly_vert[ipoly] > 0) {
        VTYPE * vlist_ptr0 = vlist + first_poly_vert[ipoly];
        VTYPE * vlist_ptr1 = vlist + vlist_length;
        *vlist_ptr1 = *vlist_ptr0;
        ITYPE k = 1;
        for (ITYPE i = 1; i < num_poly_vert[ipoly]; i++) {
          if (std::find(vlist_ptr1, vlist_ptr1+k, *(vlist_ptr0+i)) ==
              vlist_ptr1+k) {
            *(vlist_ptr1+k) = *(vlist_ptr0+i);
            k++;
          }
        }
        num_poly_vert[ipoly] = k;
        first_poly_vert[ipoly] = vlist_length;
        vlist_length = vlist_length + k;
      }
      else {
        first_poly_vert[ipoly] = vlist_length;
      }
    }
  }


  /// @brief Delete polytope vertices with duplicate coordinates.
  /// *** ORDER OF vlist and num_poly_vert INCONSISTENT WITH IO ROUTINES ***
  template <typename DTYPE, typename CTYPE, typename NUMV_TYPE, typename VTYPE,
            typename ITYPE, typename NTYPE>
  void delete_poly_vert_duplicate_coord
  (const DTYPE dimension, const CTYPE * vertex_coord, VTYPE * vlist, 
   NUMV_TYPE * num_poly_vert, ITYPE * first_poly_vert, const NTYPE num_poly)
  {
    // Duplicate vertices automatically have duplicate coordinates.
    // Delete duplicate vertices.
    delete_poly_vert_duplicate
      (vlist, num_poly_vert, first_poly_vert, num_poly);

    // Delete distinct vertices with duplicate coordinates.
    ITYPE vlist_length = 0;
    for (NTYPE ipoly = 0; ipoly < num_poly; ipoly++) {
      if (num_poly_vert[ipoly] > 0) {
        VTYPE * vlist_ptr0 = vlist + first_poly_vert[ipoly];
        VTYPE * vlist_ptr1 = vlist + vlist_length;
        *vlist_ptr1 = *vlist_ptr0;
        ITYPE k = 1;
        for (ITYPE i = 1; i < num_poly_vert[ipoly]; i++) {
          VTYPE iv0 = *(vlist_ptr0+i);
          bool flag_duplicate = false;
          for (ITYPE j = 0; j < k; j++) {
            VTYPE iv1 = *(vlist_ptr1+j);
            const CTYPE * vcoord0 = vertex_coord+iv0*dimension;
            const CTYPE * vcoord1 = vertex_coord+iv1*dimension;
            if (std::equal(vcoord0, vcoord0+dimension, vcoord1)) {
              flag_duplicate = true;
              break;
            }
          }

          if (!flag_duplicate) {
            *(vlist_ptr1+k) = *(vlist_ptr0+i);
            k++;
          }
        }
        num_poly_vert[ipoly] = k;
        first_poly_vert[ipoly] = vlist_length;
        vlist_length = vlist_length + k;
      }
      else {
        first_poly_vert[ipoly] = vlist_length;
      }
    }
  }

  //@}


  // **************************************************
  /// @name DELETE MESH ELEMENTS OUTSIDE REGION
  // **************************************************

  //@{

  /// Return true if coordinates are is in region.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2>
  bool region_contains_coord
  (const DTYPE dimension, const CTYPE0 * min_coord, const CTYPE1 * max_coord,
   const CTYPE2 * coord)
  {
    for (DTYPE d = 0; d < dimension; d++) {
      if ((coord[d] < min_coord[d]) ||
          (coord[d] > max_coord[d]))
        { return(false); }
    }

    return(true);
  }


  /// Return true if all polytope vertices are is in region.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename VTYPE, typename NUMV_TYPE>
  bool region_contains_poly
  (const DTYPE dimension, const CTYPE1 * min_coord, const CTYPE2 * max_coord,
   const CTYPE0 * vertex_coord, VTYPE * poly_vert, NUMV_TYPE num_poly_vert)
  {
    for (NUMV_TYPE i = 0; i < num_poly_vert; i++) {
      VTYPE iv = poly_vert[i];
      const CTYPE0 * coord = vertex_coord+iv*dimension;
      if (!region_contains_coord(dimension, min_coord, max_coord, coord))
        { return(false); }
    }

    return(true);
  }


  /*!
   *  @brief Delete any polytopes with vertices outside the given region.
   *  - *** ORDER OF vlist and num_poly_vert INCONSISTENT WITH IO ROUTINES ***
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename NUMV_TYPE, typename VTYPE, typename ITYPE, typename NTYPE>
  void delete_poly_outside_region
  (const DTYPE dimension, const CTYPE0 * min_coord, const CTYPE1 * max_coord,
   const CTYPE2 * vertex_coord, 
   VTYPE * vlist, NUMV_TYPE * num_poly_vert, ITYPE * first_poly_vert, 
   NTYPE & num_poly)
  {
    ITYPE vlist_length = 0;
    NTYPE k = 0;
    for (NTYPE ipoly = 0; ipoly < num_poly; ipoly++) {
      VTYPE * polyvert = vlist+first_poly_vert[ipoly];
      if (region_contains_poly(dimension, min_coord, max_coord, vertex_coord, 
                               polyvert, num_poly_vert[ipoly])) {
        if (k < ipoly) {
          VTYPE * vlist_ptr0 = vlist + first_poly_vert[ipoly];
          VTYPE * vlist_ptr1 = vlist + vlist_length;
          std::copy(vlist_ptr0, vlist_ptr0+num_poly_vert[ipoly], vlist_ptr1);
          num_poly_vert[k] = num_poly_vert[ipoly];
          first_poly_vert[k] = vlist_length;
        }
        vlist_length += num_poly_vert[k];
        k++;
      }
    }
    num_poly = k;
  }


  /*!
   *  @brief Delete any vertices outside the given region.
   *  - Also deletes any polytopes incident on those vertices
   *  - *** ORDER OF vlist and num_poly_vert INCONSISTENT WITH IO ROUTINES ***
   */
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, typename CTYPE2,
            typename NTYPE,
            typename NUM_POLYV_TYPE, typename VTYPE, typename ITYPE,
            typename NUMP_TYPE>
  void delete_vert_outside_region
  (const DTYPE dimension, 
   const CTYPE0 * min_coord, const CTYPE1 * max_coord,
   CTYPE2 * vertex_coord, NTYPE & num_vert,
   VTYPE * vlist, NUM_POLYV_TYPE * num_poly_vert, ITYPE * first_poly_vert, 
   NUMP_TYPE & num_poly)
  {
    // new_index[i] = new index of vertex i.
    IJK::ARRAY<VTYPE> new_index(num_vert);

    delete_poly_outside_region
      (dimension, min_coord, max_coord, vertex_coord, 
       vlist, num_poly_vert, first_poly_vert, num_poly);

    NTYPE k = 0;
    for (NTYPE iv = 0; iv < num_vert; iv++) {
      if (region_contains_coord
          (dimension, min_coord, max_coord, vertex_coord+iv*dimension)) {

        if (k < iv) {
          const CTYPE0 * coord0_ptr = vertex_coord+dimension*iv;
          CTYPE1 * coord1_ptr = vertex_coord+dimension*k;
          std::copy(coord0_ptr, coord0_ptr+dimension, coord1_ptr);
        }
        new_index[iv] = k;
        k++;
      }
    }

    num_vert = k;

    ITYPE vlist_length =
      std::accumulate(num_poly_vert, num_poly_vert+num_poly, 0);
    for (ITYPE i = 0; i < vlist_length; i++) 
      { vlist[i] = new_index[vlist[i]]; }
  }

  //@}


  // ***********************************************************************
  /// @name DELETE MESH POLY THAT ARE NOT INCIDENT ON ANY VERTICES IN LIST
  // ***********************************************************************

  //@{


  /*!
   *  @brief Delerte mesh polygons that are not in incident on any vertices
   *    in list.
   */
  template <typename VTYPE0, typename VTYPE1, typename NUMV_TYPE, 
            typename ITYPE, typename NTYPE>
  void delete_poly_not_incident_on_vertices
  (const std::vector<VTYPE0> & vertex_list, VTYPE1 * poly_vlist, 
   NUMV_TYPE * num_poly_vert, ITYPE * first_poly_vert, NTYPE & num_poly)
  {
    ITYPE vlist_length = 0;
    NTYPE k = 0;
    for (NTYPE ipoly = 0; ipoly < num_poly; ipoly++) {
      const VTYPE1 * polyvert = poly_vlist+first_poly_vert[ipoly];

      bool flag_contains = false;
      if (do_lists_intersect
          (polyvert, num_poly_vert[ipoly], vertex_list)) {

        if (k < ipoly) {
          VTYPE1 * poly_vlist_ptr0 = poly_vlist + first_poly_vert[ipoly];
          VTYPE1 * poly_vlist_ptr1 = poly_vlist + vlist_length;
          std::copy(poly_vlist_ptr0, poly_vlist_ptr0+num_poly_vert[ipoly], 
                    poly_vlist_ptr1);
          num_poly_vert[k] = num_poly_vert[ipoly];
          first_poly_vert[k] = vlist_length;
        }
        vlist_length += num_poly_vert[k];
        k++;
      }
    }
    num_poly = k;
  }

  //@}


  // **************************************************
  /// @name REVERSE ORIENTATIONS
  // **************************************************

  //@{

  /// @brief Reverse simplex orientations
  template <typename DTYPE, typename SVERT, typename NTYPE>
  void reverse_simplex_orientations
  (const DTYPE simplex_dimension,
   SVERT * simplex_vertex, const NTYPE num_simplices)
  {
    if (simplex_dimension < 1) { return; }

    const DTYPE num_simplex_vert = simplex_dimension+1;
    const DTYPE inc1 = num_simplex_vert-1;
    const DTYPE inc2 = num_simplex_vert-2;

    SVERT * svert0 = simplex_vertex;
    for (NTYPE i = 0; i < num_simplices; i++) {
      
      std::swap(*(svert0+inc1), *(svert0+inc2));
      svert0 += num_simplex_vert;
    }
  }

  
  /// @brief Reverse simplex orientations
  template <typename DTYPE, typename SVERT>
  void reverse_simplex_orientations
  (const DTYPE simplex_dimension, std::vector<SVERT> & simplex_vertex)
  {
    if (simplex_vertex.size() > 0) {
      reverse_simplex_orientations
        (simplex_dimension, &(simplex_vertex.front()),
         simplex_vertex.size()/(simplex_dimension+1));
    }
  }

  
  /// @brief Reverse quadrilateral orientations
  template <typename QVERT, typename NTYPE>
  void reverse_quad_orientations
  (QVERT * quad_vertex, const NTYPE num_quad)
  {
    const NTYPE NUM_VERT_PER_QUAD = 4;

    QVERT * qvert0 = quad_vertex;
    for (NTYPE i = 0; i < num_quad; i++) {
      std::swap(*(qvert0), *(qvert0+3));
      std::swap(*(qvert0+1), *(qvert0+2));
      qvert0 += NUM_VERT_PER_QUAD;
    }
  }

  
  /// @brief Reverse quadrilateral orientations
  template <typename QVERT>
  void reverse_quad_orientations(std::vector<QVERT> & quad_vertex)
  {
    const int NUM_VERT_PER_QUAD = 4;

    if (quad_vertex.size() > 0) {
      reverse_quad_orientations
        (&(quad_vertex.front()), quad_vertex.size()/NUM_VERT_PER_QUAD);
    }
  }

  
  /// @brief Reverse orientation of a polygon.
  template <typename QVERT, typename NTYPE>
  void reverse_polygon_orientation
  (QVERT * polygon_vertex, const NTYPE num_polygon_vert)
  {
    for (NTYPE i = 0; i < num_polygon_vert/2; i++) {
      NTYPE i2 = num_polygon_vert-i-1;
      std::swap(polygon_vertex[i], polygon_vertex[i2]); 
    }
  }

  
  /// @brief Reverse orientations of polygons in a list.
  template <typename VTYPE, typename NTYPE0, typename NTYPE1, typename ITYPE>
  void reverse_orientations_polygon_list
  (VTYPE * poly_vert, const NTYPE0 * num_poly_vert, 
   const ITYPE * first_poly_vert, const NTYPE1 num_poly)
  {
    for (NTYPE1 ipoly = 0; ipoly < num_poly; ipoly++) {
      reverse_polygon_orientation
        (poly_vert+first_poly_vert[ipoly], num_poly_vert[ipoly]);
    }
  }


  /*!
   *  @overload
   *  @brief Reverse orientations of polygons in a list.
   *  - C++ STL vector format for num_poly_vert[], 
   *    poly_vert[], first_poly_vert[].
   */
  template <typename VTYPE, typename NTYPE, typename ITYPE>
  void reverse_orientations_polygon_list
  (std::vector<VTYPE> & poly_vert, const std::vector<NTYPE> & num_poly_vert, 
   const std::vector<ITYPE> & first_poly_vert)
  {
    if (poly_vert.size() > 0 && num_poly_vert.size() > 0) {
      reverse_orientations_polygon_list
        (&(poly_vert.front()), &(num_poly_vert.front()), 
         &(first_poly_vert.front()), num_poly_vert.size());
    }
  }


  /*!
   *  @brief Reverse cube orientation.
   *  @param cube_vertex[] Array of cube vertices.
   *  @param num_cube_facet_vert Number of vertices per cube facet
   */
  template <typename VTYPE, typename NTYPE>
  void reverse_cube_orientation
  (VTYPE * cube_vertex, const NTYPE num_vert_per_cube_facet)
  {
    // Swap left and right facet vertices to reverse cube orientation.
    for (NTYPE i = 0; i < num_vert_per_cube_facet; i++) {
      std::swap(cube_vertex[i], cube_vertex[i+num_vert_per_cube_facet]);
    }
  }


  /*!
   *  @brief Reverse orientations of cubes in list.
   *  @param cube_vertex[] Array of cube vertices.
   *  @param num_cubes Number of cubes.
   *  @param num_cube_facet_vert Number of vertices per cube facet.
   */
  template <typename VTYPE, typename NTYPE0, typename NTYPE1>
  void reverse_orientations_cube_list
  (VTYPE * cube_vertex, 
   const NTYPE0 num_cubes, const NTYPE1 num_vert_per_cube_facet)
  {
    const NTYPE1 num_vert_per_cube = 2*num_vert_per_cube_facet;

    for (NTYPE0 i = 0; i < num_cubes; i++) {
      reverse_cube_orientation
        (cube_vertex+i*num_vert_per_cube, num_vert_per_cube_facet);
    }
  }


  /*!
   *  @overload
   *  @brief Reverse orientations of cubes in list. (C++ vector.)
   *  - C++ STL vector format for cube_vertex[].
   */
  template <typename VTYPE, typename NTYPE>
  void reverse_orientations_cube_list
  (std::vector<VTYPE> & cube_vertex, const NTYPE num_vert_per_cube_facet)
  {
    typedef typename std::vector<VTYPE>::size_type SIZE_TYPE;

    const NTYPE num_vert_per_cube = 2*num_vert_per_cube_facet;
    const SIZE_TYPE num_cubes = cube_vertex.size()/num_vert_per_cube;

    if (num_cubes > 0) {
      reverse_orientations_cube_list
        (&(cube_vertex.front()), num_cubes, num_vert_per_cube_facet);
    }
  }

  //@}


  // **************************************************
  /// @name REORDER QUADRILATERALS
  // **************************************************

  //@{

  /*!
   *  @brief Reorder quad vertices.
   *  - Swap last two vertices.
   *  - Changes counter-clockwise order to lower-left, lower-right, 
   *    upper-left, upper-right order.
   *  - Changes lower-left, lower-right, upper-left, upper-right order
   *    to counter-clockwise order.
   */
  template <typename VTYPE, typename NTYPE>
  void reorder_quad_vertices
  (VTYPE * quad_vert, const NTYPE num_quad)
  {
    const NTYPE NUM_VERT_PER_QUAD = 4;

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      NTYPE k = iquad*NUM_VERT_PER_QUAD;
      std::swap(quad_vert[k+2], quad_vert[k+3]);
    }
  }

  /*!
   *  \overload
   *  @brief Reorder quad vertices. (C++ vector.)
   *  - C++ vector version of quad_vert.
   */
  template <typename VTYPE>
  void reorder_quad_vertices(std::vector<VTYPE> & quad_vert)
  {
    typedef typename std::vector<VTYPE>::size_type SIZE_TYPE;
    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
    const SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (num_quad == 0) { return; }

    reorder_quad_vertices
      (IJK::vector2pointer(quad_vert), num_quad);
  }
  
  //@}


  // **************************************************
  /// @name GET NON-DEGENERATE QUADRILATERALS
  // **************************************************

  //@{

  /*!
   *  @brief Get a single non-degenerate quadrilateral.
   *  - A quadrilateral with a single edge collapse are reported as triangles.
   *  - Note: Does not clear new_quad_vert or new_tri_vert.
   *  @param quad_vert Four quadrilateral vertices.
   *  @pre quadrilateral vertices are in counter-clockwise order
   *         around the quadrilateral.
   */
  template <typename VTYPE1, typename VTYPE2, typename VTYPE3>
  void get_non_degenerate_quad_ccw
  (const VTYPE1 * quad_vert,
   std::vector<VTYPE2> & new_tri_vert,
   std::vector<VTYPE3> & new_quad_vert)
  {
    typedef typename std::vector<VTYPE1>::size_type SIZE_TYPE;
    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;

    for (SIZE_TYPE j1 = 0; j1 < 2; j1++) {
      SIZE_TYPE j2 = j1+2;
      if (quad_vert[j1] == quad_vert[j2]) {
        // Degenerate quad. Collapses to two edges.
        return;
      }
    }

    for (SIZE_TYPE j1 = 0; j1 < NUM_VERT_PER_QUAD; j1++) {
      SIZE_TYPE j2 = (j1+1)%NUM_VERT_PER_QUAD;
      SIZE_TYPE j3 = (j2+1)%NUM_VERT_PER_QUAD;
      SIZE_TYPE j4 = (j3+1)%NUM_VERT_PER_QUAD;

      if (quad_vert[j1] == quad_vert[j2]) {
        if (quad_vert[j3] != quad_vert[j4]) {
          // Single edge collapse.
          new_tri_vert.push_back(quad_vert[j2]);
          new_tri_vert.push_back(quad_vert[j3]);
          new_tri_vert.push_back(quad_vert[j4]);
          return;
        }
        else {
          // Degenerate quad. Collapses to two edges.
          return;
        }
      }
    }

    // Non-degenerate quadrilateral.
    for (SIZE_TYPE j = 0; j < NUM_VERT_PER_QUAD; j++) {
      new_quad_vert.push_back(quad_vert[j]);
    }
  }


  /*!
   *  @brief Get non-degenerate quadrilaterals.
   *  -Quadrilaterals with a single edge collapse are reported as triangles.
   *  - Note: Does not clear new_quad_vert or new_tri_vert.
   *  @param quad_vert Array of quadrilateral vertices.
   *         - quad_vert[4*i+j] = j'th vertex of quadrilateral i.
   *  @pre Quadrilateral vertices are in counter-clockwise order 
   *         around quadrilateral.
   */
  template <typename VTYPE1, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  void get_non_degenerate_quad_ccw
  (const VTYPE1 * quad_vert, const NTYPE num_quad,
   std::vector<VTYPE2> & new_tri_vert,
   std::vector<VTYPE3> & new_quad_vert)
  {
    const NTYPE NUM_VERT_PER_QUAD = 4;
    NTYPE num_distinct;
    VTYPE1 vlist[NUM_VERT_PER_QUAD];

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      NTYPE k = iquad*NUM_VERT_PER_QUAD;

      IJK::get_non_degenerate_quad_ccw
        (quad_vert+k, new_tri_vert, new_quad_vert);
    }

  }


  /*!
   *  @overload
   *  @brief Get non-degenerate quadrilaterals.
   *  -C++ vector version for quad_vert.
   */
  template <typename VTYPE1, typename VTYPE2, typename VTYPE3>
  void get_non_degenerate_quad_ccw
  (const std::vector<VTYPE1> & quad_vert,
   std::vector<VTYPE2> & new_tri_vert,
   std::vector<VTYPE3> & new_quad_vert)
  {
    typedef typename std::vector<VTYPE1>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;

    SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (quad_vert.size() == 0) { return; }
    get_non_degenerate_quad_ccw
      (IJK::vector2pointer(quad_vert), num_quad,
       new_tri_vert, new_quad_vert);
  }


  /*!
   *  @brief Get a single non-degenerate quadrilateral.
   *  - A quadrilateral with a single edge collapse is reported as a triangle.
   *  -Note: Does not clear new_quad_vert or new_tri_vert.
   *  @param quad_vert Four quadrilateral vertices.
   *  @pre Quadrilateral vertices are ordered bottom-left, bottom-right,
   *         top-left, top-right.
   */
  template <typename VTYPE1, typename VTYPE2, typename VTYPE3>
  void get_non_degenerate_quad_btlr
  (const VTYPE1 * quad_vert,
   std::vector<VTYPE2> & new_tri_vert,
   std::vector<VTYPE3> & new_quad_vert)
  {
    typedef typename std::vector<VTYPE1>::size_type SIZE_TYPE;
    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;
    static const SIZE_TYPE next_vert_ccw[NUM_VERT_PER_QUAD] = { 1, 3, 0, 2};

    for (SIZE_TYPE j1 = 0; j1 < 2; j1++) {
      SIZE_TYPE j2 = 3-j1;
      if (quad_vert[j1] == quad_vert[j2]) {
        // Degenerate quad. Collapses to two edges.
        return;
      }
    }

    for (SIZE_TYPE j1 = 0; j1 < NUM_VERT_PER_QUAD; j1++) {
      SIZE_TYPE j2 = next_vert_ccw[j1];
      SIZE_TYPE j3 = next_vert_ccw[j2];
      SIZE_TYPE j4 = next_vert_ccw[j3];


      if (quad_vert[j1] == quad_vert[j2]) {
        if (quad_vert[j3] != quad_vert[j4]) {
          // Single edge collapse.
          new_tri_vert.push_back(quad_vert[j2]);
          new_tri_vert.push_back(quad_vert[j3]);
          new_tri_vert.push_back(quad_vert[j4]);
          return;
        }
        else {
          // Degenerate quad. Collapses to two edges.
          return;
        }
      }
    }

    // Non-degenerate quadrilateral.
    for (SIZE_TYPE j = 0; j < NUM_VERT_PER_QUAD; j++) {
      new_quad_vert.push_back(quad_vert[j]);
    }
  }


  /*!
   *  @brief Get non-degenerate quadrilaterals.
   *  -Quadrilaterals with a single edge collapse are reported as triangles.
   *  - Note: Does not clear new_quad_vert or new_tri_vert.
   *  @param quad_vert Array of quadrilateral vertices.
   *         quad_vert[4*i+j] = j'th vertex of quadrilateral i.
   *  @pre Quadrilateral vertices are ordered bottom-left, bottom-right,
   *         top-left, top-right.
   */
  template <typename VTYPE1, typename VTYPE2, typename VTYPE3,
            typename NTYPE>
  void get_non_degenerate_quad_btlr
  (const VTYPE1 * quad_vert, const NTYPE num_quad,
   std::vector<VTYPE2> & new_tri_vert,
   std::vector<VTYPE3> & new_quad_vert)
  {
    const NTYPE NUM_VERT_PER_QUAD = 4;
    NTYPE num_distinct;
    VTYPE1 vlist[NUM_VERT_PER_QUAD];

    for (NTYPE iquad = 0; iquad < num_quad; iquad++) {
      NTYPE k = iquad*NUM_VERT_PER_QUAD;

      IJK::get_non_degenerate_quad_btlr
        (quad_vert+k, new_tri_vert, new_quad_vert);
    }

  }

  
  /*!
   *  @overloadu
   *  @brief Get non-degenerate quadrilaterals. (C++ vector.)
   *  - C++ vector version for quad_vert.
   */
  template <typename VTYPE1, typename VTYPE2, typename VTYPE3>
  void get_non_degenerate_quad_btlr
  (const std::vector<VTYPE1> & quad_vert,
   std::vector<VTYPE2> & new_tri_vert,
   std::vector<VTYPE3> & new_quad_vert)
  {
    typedef typename std::vector<VTYPE1>::size_type SIZE_TYPE;

    const SIZE_TYPE NUM_VERT_PER_QUAD = 4;

    SIZE_TYPE num_quad = quad_vert.size()/NUM_VERT_PER_QUAD;

    if (quad_vert.size() == 0) { return; }
    get_non_degenerate_quad_btlr
      (IJK::vector2pointer(quad_vert), num_quad,
       new_tri_vert, new_quad_vert);
  }

  //@}


  // **************************************************
  /// @name GET HEXAHEDRON FACET DIAGONAL
  // **************************************************

  //@{

  /*!
   *  @brief Get diagonal of hexahedron facet jfacet.
   *  @param flag_tri Indicate which diagonal.
   *     - If true, return diagonal containing vertex 0 or 7.
   *     - If false, return diagonal not containing vertex 0 or 7.
   *  @param[out] diagonal_endpoint[] Endpoints of diagonal.
   *    diagonal_endpoint[0] <= diagonal_endpoint[1].
   */
  template <typename VTYPE, typename FTYPE, typename ETYPE>
  void get_hexahedron_facet_diagonal
  (const VTYPE hex_vert[], const FTYPE jfacet, const bool flag_tri,
   ETYPE diagonal_endpoint[2])
  {
    const ETYPE facet_diagA[] = { 0, 6, 0, 5, 0, 3, 1, 7, 2, 7, 4, 7 };
    const ETYPE facet_diagB[] = { 2, 4, 1, 4, 1, 2, 3, 5, 3, 6, 5, 6 };

    if (flag_tri) {
      diagonal_endpoint[0] = hex_vert[facet_diagB[2*jfacet]];
      diagonal_endpoint[1] = hex_vert[facet_diagB[2*jfacet+1]];
    }
    else {
      diagonal_endpoint[0] = hex_vert[facet_diagA[2*jfacet]];
      diagonal_endpoint[1] = hex_vert[facet_diagA[2*jfacet+1]];
    }

    if (diagonal_endpoint[0] > diagonal_endpoint[1])
      { std::swap(diagonal_endpoint[0], diagonal_endpoint[1]); }
  }

  //@}


  // **************************************************
  /// @name SORT COORDINATES
  // **************************************************

  /*!
   *  @brief Sort point coordinates in increasing lexicographic order.
   *  @param dimension Dimension. Number of coordinates per point.
   *  @param coord[] Array of point coordinates.
   *    coord[i*dimension+j] is j'th coordinate of point i.
   *    - Modified by this routine.
   *  @param num_points Number of points.
   *  @param[out] loc_map[i] After sorting, coordinates of point i 
   *    are at coord[loc_map[i]*dimension,loc_map[i]*dimension+1,...].
   *    @pre Array loc_map[] is preallocated to size at least num_points.
   */
  template <typename CTYPE, typename NTYPE, typename LTYPE>
  void sort_point_coordinates
  (const int dimension, CTYPE * coord, const NTYPE num_points,
   LTYPE * loc_map)
  {
    TUPLE_LESS_THAN<int, CTYPE> point_less_than(dimension, coord);

    IJK::ARRAY<NTYPE> index_sorted(num_points);
    std::vector<bool> is_set(num_points, false);

    for (NTYPE i = 0; i < num_points; i++)
      { index_sorted[i] = i; }

    std::sort(index_sorted.Ptr(), index_sorted.Ptr()+num_points,
         point_less_than);

    IJK::invert_permutation(index_sorted.Ptr(), loc_map, num_points);

    // Permute coord using loc_map[].
    for (NTYPE i0 = 0; i0 < num_points; i0++) {
      if (!is_set[i0]) {
        NTYPE i1 = loc_map[i0];
        while (i1 != i0) {
          CTYPE * c0 = coord+i0*dimension;
          CTYPE * c1 = coord+i1*dimension;
          std::swap_ranges(c0, c0+dimension, c1);
          is_set[i1] = true;
          i1 = loc_map[i1];
        }
      }
    }
  }


  // **************************************************
  /// @name SORT SIMPLEX VERTICES AND SIMPLICES
  // **************************************************

  /*!
   *  @brief Sort simplex vertices, preserving simplex orientation.
   *  - Bubble sort vertices of each simplex, in increasing order.
   *  @param numv_per_simplex Number of vertices per simplex.
   *  @param nums Number of simplices
   *  @param simplex_vert[] Array of simplex vertices.
   *    - Simplex[i*numv_per_simplex+j] is the j'th vertex of simplex i.
   */
  template <typename NTYPE0, typename NTYPE1, typename VTYPE>
  void sort_oriented_simplex_vertices
  (const NTYPE0 numv_per_simplex, const NTYPE1 nums, VTYPE simplex_vert[])
  {
    if (numv_per_simplex > 1) {
      for (NTYPE1 js = 0; js < nums; js++) {
        bool orientation = true;
        bool done;
        do {
          done = true;
          for (NTYPE1 i = js*numv_per_simplex; 
               i+1 < (js+1)*numv_per_simplex; i++) {
            VTYPE iv0 = simplex_vert[i];
            VTYPE iv1 = simplex_vert[i+1];
            if (iv0 > iv1) {
              std::swap(simplex_vert[i], simplex_vert[i+1]);
              orientation = !orientation;
              done = false;
            }
          }
        } while (!done);

        // if necessary, swap two last vertices for correct orientation
        if (!orientation) {
          NTYPE1 ilast = js*numv_per_simplex + numv_per_simplex - 1;
          std::swap(simplex_vert[ilast], simplex_vert[ilast-1]);
        }
      }
    }
  }


  /*!
   *  @brief Sort simplices by lexicographic order of vertex indices.
   *  - Calls sort_oriented_simplex_vertices to sort all simplex vertices.
   *  @param numv_per_simplex Number of vertices per simplex.
   *  @param nums Number of simplices
   *  @param simplex_vert[] Array of simplex vertices.
   *    - Simplex[i*numv_per_simplex+j] is the j'th vertex of simplex i.
   */
  template <typename NTYPE0, typename NTYPE1, typename VTYPE>
  void sort_simplices
  (const NTYPE0 numv_per_simplex, const NTYPE1 nums, VTYPE simplex_vert[])
  {
    TUPLE_LESS_THAN<NTYPE0,VTYPE> 
      simplex_vert_less_than(numv_per_simplex, simplex_vert);

    IJK::ARRAY<NTYPE1> index_sorted(nums);
    IJK::ARRAY<NTYPE1> index_sorted2(nums);
  
    sort_oriented_simplex_vertices(numv_per_simplex, nums, simplex_vert);

    for (NTYPE1 i = 0; i < nums; i++)
      { index_sorted[i] = i; }

    // sort simplices by lexicographic order of vertices
    std::sort(index_sorted.Ptr(), index_sorted.Ptr()+nums, 
              simplex_vert_less_than);

    // reverse the permutation
    // Note: Probably could do this in place, but...
    for (NTYPE1 i = 0; i < nums; i++) {
      NTYPE1 k0 = index_sorted[i];
      index_sorted2[k0] = i;
    }

    // permute simplex_vert using index_sorted2
    for (NTYPE1 i0 = 0; i0 < nums; i0++) {
      NTYPE1 i1 = index_sorted2[i0];
      while (i1 != i0) {
        VTYPE * c0 = simplex_vert+i0*numv_per_simplex;
        VTYPE * c1 = simplex_vert+i1*numv_per_simplex;
        std::swap_ranges(c0, c0+numv_per_simplex, c1);
        std::swap(index_sorted2[i0], index_sorted2[i1]);
        i1 = index_sorted2[i0];
      }
    }

  }


  /*!
   *  @overload
   *  @brief Sort simplices by lexicographic order of vertex indices. (C++ vector.)
   *  - Version using C++ STL vector for array simplex_vert[].
   */
  template <typename NTYPE0, typename VTYPE>
  void sort_simplices
  (const NTYPE0 numv_per_simplex, std::vector<VTYPE> & simplex_vert)
  {
    typedef typename std::vector<VTYPE>::size_type SIZE_TYPE;

    SIZE_TYPE nums = simplex_vert.size()/numv_per_simplex;
    sort_simplices(numv_per_simplex, nums, IJK::vector2pointer(simplex_vert));
  }

  // **************************************************
  /// @name CHECK ROUTINES
  // **************************************************

  /*!
   *  @brief Return true if vlist.size() is divisible by NUM_VERT_PER_POLY.
   */
  template <int NUM_VERT_PER_POLY, typename VTYPE>
  bool check_vertex_vector_size
  (const std::vector<VTYPE> & vlist, IJK::ERROR & error)
  {
    if (vlist.size()%NUM_VERT_PER_POLY == 0)
      { return(true); }
    else {
      error.AddMessage
	("Programming error. Incorrect vertex list size.");
      error.AddMessage
	("  Vertex list size is not divisible by ", NUM_VERT_PER_POLY,
	 ", the number of vertices in each polytope.");
      error.AddMessage
	("  vlist.size() = ", vlist.size(), "");

      return(false);
    }
  }

}

#endif
