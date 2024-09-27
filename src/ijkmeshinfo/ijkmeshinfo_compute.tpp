/*!
 *  @file ijkmeshinfo_compute.tpp
 *  @brief Templates for computing angles, edge lengths, etc.
 *  - Version 0.4.0
 */

/*
  IJK: Isosurface Jeneration Code
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


#ifndef _IJKMESHINFO_COMPUTE_TPP_
#define _IJKMESHINFO_COMPUTE_TPP_

#include "ijkmeshinfo.h"

namespace IJKMESHINFO {

  // **************************************************
  //! @name Compute values for each polytope.
  // **************************************************

  /// Compute the minimum and maximum values for each polytope in polymesh.
  template <typename MDATA_TYPE, typename CTYPE,
            typename MINVAL_TYPE0, typename MAXVAL_TYPE0, 
            typename MINVAL_TYPE1, typename MAXVAL_TYPE1, 
            typename PTYPE0, typename PTYPE1,
            typename VTYPE, typename NTYPE0, typename NTYPE1>
  void compute_min_max_plist_values
  (const MDATA_TYPE & mesh_data,
   const COORD_TYPE * vertex_coord, const bool flag_internal,
   const MINVAL_TYPE0 min_init, const MAXVAL_TYPE0 max_init,
   MINVAL_TYPE1 & min_value, MAXVAL_TYPE1 & max_value,
   PTYPE0 & poly_with_min_value, PTYPE1 & poly_with_max_value,
   void compute_min_max_poly_values
   (const MDATA_TYPE &, const VTYPE *, const NTYPE0, const CTYPE *,
    MINVAL_TYPE1 &, MAXVAL_TYPE1 &, NTYPE1 &))
  {
    min_value = min_init;
    max_value = max_init;
    poly_with_min_value = 0;
    poly_with_max_value = 0;

    if (flag_internal) {
      IJK::PROCEDURE_ERROR error("compute_min_max_plist_values");
      if (!check_boundary_facets(mesh_data, error)) { throw error; } 
    }

    bool is_set(false);
    for (int ipoly = 0; ipoly < mesh_data.NumPoly(); ipoly++) {

      if (mesh_data.poly_data[ipoly].IsDegenerate()) { continue; }

      if (flag_internal) {
        if (mesh_data.poly_data[ipoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      MINVAL_TYPE1 min_value_i;
      MAXVAL_TYPE1 max_value_i;
      NTYPE1 num_values;
      compute_min_max_poly_values
        (mesh_data, mesh_data.VertexList(ipoly), 
         mesh_data.NumPolyVert(ipoly),
         vertex_coord, min_value_i, max_value_i, num_values);

      if (num_values > 0) {
        if (!is_set || min_value_i < min_value) { 
          min_value = min_value_i; 
          poly_with_min_value = ipoly;
        }

        if (!is_set || max_value_i > max_value) { 
          max_value = max_value_i; 
          poly_with_max_value = ipoly;
        }

        is_set = true;
      }
    }

  }


  /// Compute the minimum values over all polytopes in polymesh.
  template <typename MDATA_TYPE, typename CTYPE,
            typename MINVAL_TYPE0, typename MINVAL_TYPE1,
            typename PTYPE0,
            typename VTYPE, typename NTYPE0, typename NTYPE1>
  void compute_min_plist_values
  (const MDATA_TYPE & mesh_data,
   const COORD_TYPE * vertex_coord, const bool flag_internal,
   const MINVAL_TYPE0 min_init, MINVAL_TYPE1 & min_value,
   PTYPE0 & poly_with_min_value,
   void compute_poly_value
   (const MDATA_TYPE &, const VTYPE *, const NTYPE0, const CTYPE *,
    MINVAL_TYPE1 &, NTYPE1 &))
  {
    min_value = min_init;
    poly_with_min_value = 0;

    if (flag_internal) {
      IJK::PROCEDURE_ERROR error("compute_min_plist_values");
      if (!check_boundary_facets(mesh_data, error)) { throw error; } 
    }

    bool is_set(false);
    for (int ipoly = 0; ipoly < mesh_data.NumPoly(); ipoly++) {

      if (mesh_data.poly_data[ipoly].IsDegenerate()) { continue; }

      if (flag_internal) {
        if (mesh_data.poly_data[ipoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      MINVAL_TYPE1 value_i;
      NTYPE1 num_values;
      compute_poly_value
        (mesh_data, mesh_data.VertexList(ipoly), 
         mesh_data.NumPolyVert(ipoly),
         vertex_coord, value_i, num_values);

      if (num_values > 0) {
        if (!is_set || value_i < min_value) { 
          min_value = value_i; 
          poly_with_min_value = ipoly;
        }

        is_set = true;
      }
    }

  }


  /// Compute the minimum and maximum values for each polytope in polymesh.
  /// - Return min/max values data.
  template <typename MDATA_TYPE, typename CTYPE,
            typename MINVAL_TYPE0, typename MAXVAL_TYPE0, 
            typename MINVAL_TYPE1, typename MAXVAL_TYPE1, 
            typename PTYPE0, typename PTYPE1,
            typename MINVAL_DATA_TYPE, typename MAXVAL_DATA_TYPE,
            typename VTYPE, typename NTYPE0, typename NTYPE1>
  void compute_min_max_plist_values_data
  (const MDATA_TYPE & mesh_data,
   const COORD_TYPE * vertex_coord, const bool flag_internal,
   const MINVAL_TYPE0 min_init, const MAXVAL_TYPE0 max_init,
   MINVAL_TYPE1 & min_value, MAXVAL_TYPE1 & max_value,
   PTYPE0 & poly_with_min_value, PTYPE1 & poly_with_max_value,
   MINVAL_DATA_TYPE & min_value_data, 
   MAXVAL_DATA_TYPE & max_value_data,
   void compute_min_max_poly_values
   (const MDATA_TYPE &, const VTYPE *, const NTYPE0, const CTYPE *,
    MINVAL_TYPE1 &, MAXVAL_TYPE1 &, 
    MINVAL_DATA_TYPE &, MAXVAL_DATA_TYPE &, NTYPE1 &))
  {
    min_value = min_init;
    max_value = max_init;
    poly_with_min_value = 0;
    poly_with_max_value = 0;

    if (flag_internal) {
      IJK::PROCEDURE_ERROR error("compute_min_max_plist_values_data");
      if (!check_boundary_facets(mesh_data, error)) { throw error; } 
    }

    bool is_set(false);
    for (int ipoly = 0; ipoly < mesh_data.NumPoly(); ipoly++) {

      if (mesh_data.poly_data[ipoly].IsDegenerate()) { continue; }

      if (flag_internal) {
        if (mesh_data.poly_data[ipoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      MINVAL_TYPE1 min_value_i;
      MAXVAL_TYPE1 max_value_i;
      MINVAL_DATA_TYPE min_value_data_i;
      MAXVAL_DATA_TYPE max_value_data_i;
      NTYPE1 num_values;
      compute_min_max_poly_values
        (mesh_data, mesh_data.VertexList(ipoly), 
         mesh_data.NumPolyVert(ipoly),
         vertex_coord, min_value_i, max_value_i, 
         min_value_data_i, max_value_data_i, num_values);

      if (num_values > 0) {
        if (!is_set || min_value_i < min_value) { 
          min_value = min_value_i;
          poly_with_min_value = ipoly;
          min_value_data = min_value_data_i;
        }

        if (!is_set || max_value_i > max_value) { 
          max_value = max_value_i; 
          poly_with_max_value = ipoly;
          max_value_data = max_value_data_i;
        }

        is_set = true;
      }
    }

  }


  /// Compute the minimum and maximum values for each polytope in polymesh
  ///   whose number of vertices equals num_poly_vert.
  template <typename MDATA_TYPE, typename CTYPE,
            typename MINVAL_TYPE0, typename MAXVAL_TYPE0, 
            typename MINVAL_TYPE1, typename MAXVAL_TYPE1, 
            typename PTYPE0, typename PTYPE1,
            typename VTYPE, typename NTYPE0, typename NTYPE1, typename NTYPE2>
  void compute_min_max_plist_values_select_poly_by_numv
  (const MDATA_TYPE & mesh_data,
   const COORD_TYPE * vertex_coord, const bool flag_internal,
   const MINVAL_TYPE0 min_init, const MAXVAL_TYPE0 max_init,
   const NTYPE2 num_poly_vert,
   MINVAL_TYPE1 & min_value, MAXVAL_TYPE1 & max_value,
   PTYPE0 & poly_with_min_value, PTYPE1 & poly_with_max_value,
   void compute_min_max_poly_values
   (const MDATA_TYPE &, const VTYPE *, const NTYPE0, const CTYPE *,
    MINVAL_TYPE1 &, MAXVAL_TYPE1 &, NTYPE1 &))
  {
    min_value = min_init;
    max_value = max_init;
    poly_with_min_value = 0;
    poly_with_max_value = 0;

    if (flag_internal) {
      IJK::PROCEDURE_ERROR error
        ("compute_min_max_plist_values_select_poly_by_numv");
      if (!check_boundary_facets(mesh_data, error)) { throw error; } 
    }

    bool is_set(false);
    for (int ipoly = 0; ipoly < mesh_data.NumPoly(); ipoly++) {

      if (mesh_data.NumPolyVert(ipoly) != num_poly_vert) { continue; }

      if (mesh_data.poly_data[ipoly].IsDegenerate()) { continue; }

      if (flag_internal) {
        if (mesh_data.poly_data[ipoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      MINVAL_TYPE1 min_value_i;
      MAXVAL_TYPE1 max_value_i;
      NTYPE1 num_values;
      compute_min_max_poly_values
        (mesh_data, mesh_data.VertexList(ipoly), 
         mesh_data.NumPolyVert(ipoly),
         vertex_coord, min_value_i, max_value_i, num_values);

      if (num_values > 0) {
        if (!is_set || min_value_i < min_value) { 
          min_value = min_value_i; 
          poly_with_min_value = ipoly;
        }

        if (!is_set || max_value_i > max_value) { 
          max_value = max_value_i; 
          poly_with_max_value = ipoly;
        }

        is_set = true;
      }
    }

  }


  /// Compute the minimum values over all polytopes in polymesh
  ///   whose number of vertices equals num_poly_vert.
  template <typename MDATA_TYPE, typename CTYPE,
            typename MINVAL_TYPE0, typename MINVAL_TYPE1,
            typename PTYPE0,
            typename VTYPE, typename NTYPE0, typename NTYPE1, typename NTYPE2>
  void compute_min_plist_values_select_poly_by_numv
  (const MDATA_TYPE & mesh_data,
   const COORD_TYPE * vertex_coord, const bool flag_internal,
   const MINVAL_TYPE0 min_init, 
   const NTYPE2 num_poly_vert,
   MINVAL_TYPE1 & min_value, 
   PTYPE0 & poly_with_min_value,
   void compute_poly_value
   (const MDATA_TYPE &, const VTYPE *, const NTYPE0, const CTYPE *,
    MINVAL_TYPE1 &, NTYPE1 &))
  {
    min_value = min_init;
    poly_with_min_value = 0;

    if (flag_internal) {
      IJK::PROCEDURE_ERROR error
        ("compute_min_plist_values_select_poly_by_numv");
      if (!check_boundary_facets(mesh_data, error)) { throw error; } 
    }

    bool is_set(false);
    for (int ipoly = 0; ipoly < mesh_data.NumPoly(); ipoly++) {

      if (mesh_data.NumPolyVert(ipoly) != num_poly_vert) { continue; }

      if (mesh_data.poly_data[ipoly].IsDegenerate()) { continue; }

      if (flag_internal) {
        if (mesh_data.poly_data[ipoly].ContainsBoundaryFacet()) 
          { continue; } 
      }

      MINVAL_TYPE1 value_i;
      NTYPE1 num_values;
      compute_poly_value
        (mesh_data, mesh_data.VertexList(ipoly), 
         mesh_data.NumPolyVert(ipoly),
         vertex_coord, value_i, num_values);

      if (num_values > 0) {
        if (!is_set || value_i < min_value) { 
          min_value = value_i; 
          poly_with_min_value = ipoly;
        }

        is_set = true;
      }
    }

  }


  // **************************************************
  //! @name Compute values for each vertex.
  // **************************************************


  /// Compute the minimum and maximum values for each vertex in polymesh.
  template <typename MDATA_TYPE,
            typename VP_INCIDENCE_TYPE, typename CTYPE,
            typename MINVAL_TYPE0, typename MAXVAL_TYPE0, 
            typename MINVAL_TYPE1, typename MAXVAL_TYPE1, 
            typename VTYPE0, typename VTYPE1, typename VTYPE2, 
            typename PTYPE0, typename PTYPE1, 
            typename PTYPE2, typename PTYPE3,
            typename NTYPE>
  void compute_min_max_vlist_values
  (const MDATA_TYPE & mesh_data, 
   const VP_INCIDENCE_TYPE & vertex_poly_incidence,
   const COORD_TYPE * vertex_coord, 
   const bool flag_internal_poly,
   const bool flag_internal_vert,
   const MINVAL_TYPE0 min_init, const MAXVAL_TYPE0 max_init,
   MINVAL_TYPE1 & min_value, MAXVAL_TYPE1 & max_value,
   PTYPE0 & poly_with_min_value, 
   PTYPE1 & poly_with_max_value, 
   VTYPE0 & vertex_with_min_value, 
   VTYPE1 & vertex_with_max_value,
   void compute_min_max_vertex_values
   (const MDATA_TYPE &, const VP_INCIDENCE_TYPE &, 
    const CTYPE *, const VTYPE2, const bool, const bool,
    MINVAL_TYPE1 &, MAXVAL_TYPE1 &, PTYPE2 &, PTYPE3 &, NTYPE &))
  {
    PTYPE2 poly_min;
    PTYPE3 poly_max;

    min_value = min_init;
    max_value = max_init;
    vertex_with_min_value = 0;
    vertex_with_max_value = 0;

    if (flag_internal_poly) {
      IJK::PROCEDURE_ERROR error("compute_min_max_vlist_values");
      if (!check_boundary_facets(mesh_data, error)) { throw error; } 
    }

    bool is_set(false);
    for (int iv = 0; iv < vertex_poly_incidence.NumVertices(); iv++) {

      MINVAL_TYPE1 min_value_i;
      MAXVAL_TYPE1 max_value_i;
      NTYPE num_values;
      compute_min_max_vertex_values
        (mesh_data, vertex_poly_incidence, vertex_coord, 
         iv, flag_internal_poly, flag_internal_vert,
         min_value_i, max_value_i, poly_min, poly_max, num_values);

      if (num_values > 0) {
        if (!is_set || min_value_i < min_value) { 
          min_value = min_value_i; 
          vertex_with_min_value = iv;
          poly_with_min_value = poly_min;
        }

        if (!is_set || max_value_i > max_value) { 
          max_value = max_value_i; 
          vertex_with_max_value = iv;
          poly_with_max_value = poly_max;
        }

        is_set = true;
      }
    }

  }

}

#endif

