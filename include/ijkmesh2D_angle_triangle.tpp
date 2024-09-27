/*!
 *  @file ijkmesh2D_angle_triangle.tpp
 *  @brief ijk template functions for computing triangulations
 *    and triangulation angles when triangles are split/replaced.
 *  - Version 0.4.0
 */


/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2019-2022 Rephael Wenger

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

#ifndef _IJKMESH2D_ANGLE_TRIANGLE_
#define _IJKMESH2D_ANGLE_TRIANGLE_

#include "ijk.tpp"
#include "ijkcoord.tpp"
#include "ijkmesh2D_angle.tpp"
#include "ijkmesh2D_geom.tpp"
#include "ijkmesh2D_split.tpp"
#include "ijktri2D_info.tpp"
#include "ijktri2D_angle.tpp"

#include <algorithm>
#include <numeric>
#include <vector>

// *** DEBUG ***
#include "ijkprint.tpp"

namespace IJK {

  // *****************************************************************
  /// @name Compute triangulations for splitting triangles.
  // *****************************************************************

  ///@{

  /// Compute cosine of the min angle in the triangulation of a triangle.
  /// - Split triangle into 4 subtriangles.
  /// @param splitA_coord[] Split vertex on line segment (vcoord0[], vcoord1[]).
  /// @param splitB_coord[] Split vertex on line segment (vcoord1[], vcoord2[]).
  /// @param splitC_coord[] Split vertex on line segment (vcoord2[], vcoord0[]).
  template <typename DTYPE, typename CTYPE0, 
            typename CTYPEA, typename CTYPEB, typename CTYPEC,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_split_triangle_triangulation_angleX
  (const DTYPE dimension, 
   const CTYPE0 vcoord0[],
   const CTYPE0 vcoord1[],
   const CTYPE0 vcoord2[],
   const CTYPEA splitA_coord[], 
   const CTYPEB splitB_coord[], 
   const CTYPEC splitC_coord[], 
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle,
   bool & flag_zero)
  {
    COS_TYPE cosAB, cosBC, cosAC, cosABC;
    bool flagAB_zero, flagBC_zero, flagAC_zero, flagABC_zero;

    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord0, splitA_coord, splitC_coord, 
       max_small_magnitude, cosAC, flagAC_zero);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord1, splitB_coord, splitA_coord, 
       max_small_magnitude, cosAB, flagAB_zero);
    IJK::compute_cos_min_triangle_angle
      (dimension, vcoord2, splitB_coord, splitC_coord, 
       max_small_magnitude, cosBC, flagBC_zero);
    IJK::compute_cos_min_triangle_angle
      (dimension, splitA_coord, splitB_coord, splitC_coord, 
       max_small_magnitude, cosABC, flagABC_zero);

    flag_zero = (flagAC_zero && flagAB_zero && flagBC_zero && flagABC_zero);
    cos_min_angle = std::max(cosAB, cosBC);
    cos_min_angle = std::max(cos_min_angle, cosAC);
    cos_min_angle = std::max(cos_min_angle, cosABC);
  }


  /// Compute cosine of the min angle in the triangulation of a triangle.
  /// - Split triangle into 4 subtriangles.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename MESH_TYPE, typename ITYPE0,
            typename MTYPE, typename COS_TYPE>
  void compute_cos_min_split_triangle_triangulation_angle
  (const DTYPE dimension, const CTYPE0 vcoord[],
   const CTYPE1 splitA_coord[], 
   const CTYPE1 splitB_coord[], 
   const CTYPE1 splitC_coord[], 
   const MESH_TYPE & mesh,  
   const ITYPE0 itriangle,
   const MTYPE max_small_magnitude,
   COS_TYPE & cos_min_angle,
   bool & flag_zero)
  {
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE 
      VERTEX_INDEX_TYPE;

    const VERTEX_INDEX_TYPE iv0 = mesh.PolygonVertex(itriangle, 0);
    const VERTEX_INDEX_TYPE iv1 = mesh.PolygonVertex(itriangle, 1);
    const VERTEX_INDEX_TYPE iv2 = mesh.PolygonVertex(itriangle, 2);
    const CTYPE0 * vcoord0 = vcoord + dimension*iv0;
    const CTYPE0 * vcoord1 = vcoord + dimension*iv1;
    const CTYPE0 * vcoord2 = vcoord + dimension*iv2;
    COS_TYPE cosAB, cosBC, cosAC, cosABC;
    bool flagAB_zero, flagBC_zero, flagAC_zero, flagABC_zero;

    // Note: splitA_coord[] lies on line segment (iv0, iv1).
    // Note: splitB_coord[] lies on line segment (iv1, iv2).
    // Note: splitC_coord[] lies on line segment (iv2, iv0).

    compute_cos_min_split_triangle_triangulation_angleX
      (dimension, vcoord0, vcoord1, vcoord2, 
       splitA_coord, splitB_coord, splitC_coord,
       max_small_magnitude, cos_min_angle, flag_zero);
  }


  /// Compute cosine of the min angle in the triangulation of a triangle.
  /// - Split triangle into 4 subtriangles.
  /// - Version that also sets triangle_tri_result.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1, 
            typename MESH_TYPE, typename ITYPE0,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_min_split_triangle_triangulation_angle
  (const DTYPE dimension, const CTYPE0 vcoord[],   
   const CTYPE1 splitA_coord[], 
   const CTYPE1 splitB_coord[], 
   const CTYPE1 splitC_coord[], 
   const MESH_TYPE & mesh,  
   const ITYPE0 itriangle,
   const MTYPE max_small_magnitude,
   RESULT_TYPE & triangle_tri_result,
   COS_TYPE & cos_min_angle,
   bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    const NUMBER_TYPE FOUR(4);

    compute_cos_min_split_triangle_triangulation_angle
      (dimension, vcoord, splitA_coord, splitB_coord, splitC_coord,
       mesh, itriangle, max_small_magnitude, cos_min_angle, flag_zero);

    triangle_tri_result.SetNoInterior(cos_min_angle, FOUR, 0);
    triangle_tri_result.flag_zero = flag_zero;
    triangle_tri_result.SetEarFlag(0, true);
    triangle_tri_result.SetEarFlag(1, true);
    triangle_tri_result.SetEarFlag(2, true);
  }

  ///@}


  // *****************************************************************
  /// @name Compute triangulations for replacing triangles.
  // *****************************************************************

  ///@{

  /// Replace triangle with triangle and edge.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 new_triangle_vcoord[],
   const CTYPE2 interior_vcoord0[],
   const CTYPE3 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ihalf_edge0);
    const VERTEX_INDEX_TYPE jv1 = mesh.ToVertexIndex(ihalf_edge0);
    const CTYPE0 * vcoord0 = vcoord + jv0*dimension;
    const CTYPE0 * vcoord1 = vcoord + jv1*dimension;
    COS_TYPE cos_min_triangle_angle;
    bool flag_zeroT;

    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1 = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge2 = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1X = mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge2X = mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge2);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << ".";
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  Triangle ", 
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge0), ".\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  poly0: ",
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge1X), ".\n");
    mesh.PrintPolygonIndexAndVertices
      (cerr, "  poly2: ",
       mesh.IndexOfPolygonContainingHalfEdge(ihalf_edge2X), ".\n");
    */

    
    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, new_triangle_vcoord,
       max_small_magnitude, poly_tri_result[1]);

    compute_split_half_edge_triangulation_to_max_min_angleP
      (dimension, vcoord, new_triangle_vcoord, interior_vcoord0, 
       mesh, ihalf_edge1X, max_small_magnitude, triangulation_settings[0],
       poly_tri_result[0]);

    compute_split_half_edge_triangulation_to_max_min_angleP
      (dimension, vcoord, new_triangle_vcoord, interior_vcoord2,
       mesh, ihalf_edge2X, max_small_magnitude, triangulation_settings[2],
       poly_tri_result[2]);


    cos_min_angle = std::max(poly_tri_result[0].cos_min_triangulation_angle,
                             poly_tri_result[1].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             poly_tri_result[2].cos_min_triangulation_angle);
    flag_zero = (poly_tri_result[0].flag_zero || 
                 poly_tri_result[1].flag_zero ||
                 poly_tri_result[2].flag_zero);
  }


  /// Replace triangle with triangle and edge.
  /// - Version using C++ STL vector vcoord[].
  template <typename DTYPE, typename CTYPE0, 
            typename CTYPE1, typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 new_triangle_vcoord[], 
   const CTYPE2 interior_vcoord0[],
   const CTYPE3 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_replace_triangle_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), new_triangle_vcoord,
       interior_vcoord0, interior_vcoord2, mesh, ihalf_edge0, 
       max_small_magnitude, triangulation_settings, 
       flag_allow_centroidx2,
       poly_tri_result, cos_min_angle, flag_zero);
  }


  /// Replace triangle with triangle and edge.
  /// - Version using C++ STL vector vcoord[].
  /// - Version using C++ STL vector new_triangle_vcoord[].
  template <typename DTYPE, typename CTYPE0,
            typename CTYPE1, typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & new_triangle_vcoord, 
   const std::vector<CTYPE2> & interior_vcoord0,
   const std::vector<CTYPE3> & interior_vcoord2,
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   const bool flag_allow_centroidx2,
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_replace_triangle_triangulation_angle
      (dimension, vcoord, IJK::vector2pointer(new_triangle_vcoord),
       IJK::vector2pointer(interior_vcoord0), 
       IJK::vector2pointer(interior_vcoord2),
       mesh, ihalf_edge0, max_small_magnitude, 
       triangulation_settings, flag_allow_centroidx2,
       poly_tri_result, cos_min_angle, flag_zero);
  }


  /// Replace triangle with triangle and edge.
  /// - Version trying two different interior vertices.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3, typename CTYPE4,
            typename CTYPE5, typename CTYPE6, typename CTYPE7,            
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle_IVx2
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 new_triangle_vcoord[],
   const CTYPE2 interior_vcoord0I[],
   const CTYPE3 interior_vcoord0II[],
   const CTYPE4 interior_vcoord2I[],
   const CTYPE5 interior_vcoord2II[],
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero,
   std::vector<CTYPE6> & selected_vcoord0,
   std::vector<CTYPE7> & selected_vcoord2)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ihalf_edge0);
    const VERTEX_INDEX_TYPE jv1 = mesh.ToVertexIndex(ihalf_edge0);
    const CTYPE0 * vcoord0 = vcoord + jv0*dimension;
    const CTYPE0 * vcoord1 = vcoord + jv1*dimension;
    COS_TYPE cos_min_triangle_angle;
    bool flag_zeroT;

    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1 = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge2 = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1X = mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge2X = mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge2);

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "*** Considering replacing triangle: ", ihalf_edge0, "\n");
    mesh.PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "  Polygon 0: ", ihalf_edge1X, "\n");
    IJK::print_coord3D
      (cerr, "  interior_vcoord0I: ", interior_vcoord0I, "\n");
    IJK::print_coord3D
      (cerr, "  interior_vcoord0II: ", interior_vcoord0II, "\n");
    mesh.PrintIndexAndVerticesOfPolygonContainingHalfEdge
      (cerr, "  Polygon 2: ", ihalf_edge2X, "\n");
    IJK::print_coord3D
      (cerr, "  interior_vcoord2I: ", interior_vcoord2I, "\n");
    IJK::print_coord3D
      (cerr, "  interior_vcoord2II: ", interior_vcoord2II, "\n");
    */


    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, new_triangle_vcoord,
       max_small_magnitude, poly_tri_result[1]);

    compute_cos_max_min_split_half_edge_polygon_triangulation_angleP_IVx2
      (dimension, vcoord, new_triangle_vcoord, 
       interior_vcoord0I, interior_vcoord0II, 
       mesh, ihalf_edge1X, max_small_magnitude, triangulation_settings[0],
       poly_tri_result[0], selected_vcoord0);

    compute_cos_max_min_split_half_edge_polygon_triangulation_angleP_IVx2
      (dimension, vcoord, new_triangle_vcoord,
       interior_vcoord2I, interior_vcoord2II,
       mesh, ihalf_edge2X, max_small_magnitude, triangulation_settings[2],
       poly_tri_result[2], selected_vcoord2);

    cos_min_angle = std::max(poly_tri_result[0].cos_min_triangulation_angle,
                             poly_tri_result[1].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             poly_tri_result[2].cos_min_triangulation_angle);
    flag_zero = (poly_tri_result[0].flag_zero || 
                 poly_tri_result[1].flag_zero ||
                 poly_tri_result[2].flag_zero);
  }


  /// Replace triangle with triangle and edge.
  /// - Version trying two different interior vertices.
  /// - Version using C++ STL vector vcoord[].
  template <typename DTYPE, typename CTYPE0, 
            typename CTYPE1, typename CTYPE2, typename CTYPE3,
            typename CTYPE4, typename CTYPE5, typename CTYPE6,
            typename CTYPE7,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle_IVx2
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 new_triangle_vcoord[], 
   const CTYPE2 interior_vcoord0I[],
   const CTYPE3 interior_vcoord0II[],
   const CTYPE4 interior_vcoord2I[],
   const CTYPE5 interior_vcoord2II[],
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero,
   std::vector<CTYPE6> & selected_vcoord0,
   std::vector<CTYPE7> & selected_vcoord2)
  {
    compute_cos_max_min_replace_triangle_triangulation_angle_IVx2
      (dimension, IJK::vector2pointer(vcoord), new_triangle_vcoord,
       interior_vcoord0I, interior_vcoord0II, 
       interior_vcoord2I, interior_vcoord2II, 
       mesh, ihalf_edge0, max_small_magnitude, triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero,
       selected_vcoord0, selected_vcoord2);
  }


  /// Replace triangle with triangle and edge.
  /// - Version trying two different interior vertices.
  /// - Version using C++ STL vector vcoord[].
  /// - Version using C++ STL vector new_triangle_vcoord[].
  template <typename DTYPE, typename CTYPE0,
            typename CTYPE1, typename CTYPE2, typename CTYPE3,
            typename CTYPE4, typename CTYPE5, typename CTYPE6,
            typename CTYPE7,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle_IVx2
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & new_triangle_vcoord, 
   const std::vector<CTYPE2> & interior_vcoord0I,
   const std::vector<CTYPE3> & interior_vcoord0II,
   const std::vector<CTYPE4> & interior_vcoord2I,
   const std::vector<CTYPE5> & interior_vcoord2II,
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero,
   std::vector<CTYPE6> & selected_vcoord0,
   std::vector<CTYPE7> & selected_vcoord2)
  {
    compute_cos_max_min_replace_triangle_triangulation_angle_IVx2
      (dimension, vcoord, IJK::vector2pointer(new_triangle_vcoord),
       IJK::vector2pointer(interior_vcoord0I), 
       IJK::vector2pointer(interior_vcoord0II),
       IJK::vector2pointer(interior_vcoord2I), 
       IJK::vector2pointer(interior_vcoord2II),
       mesh, ihalf_edge0, max_small_magnitude, 
       triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero,
       selected_vcoord0, selected_vcoord2);
  }


  /// Replace triangle with triangle and edge.
  /// - Allow interior vertex in triangulation.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle_allow_IV
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 new_triangle_vcoord[], 
   const CTYPE2 interior_vcoord0[],
   const CTYPE3 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ihalf_edge0);
    const VERTEX_INDEX_TYPE jv1 = mesh.ToVertexIndex(ihalf_edge0);
    const CTYPE0 * vcoord0 = vcoord + jv0*dimension;
    const CTYPE0 * vcoord1 = vcoord + jv1*dimension;
    COS_TYPE cos_min_triangle_angle;
    bool flag_zeroT;

    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1 = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge0);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge2 = mesh.IndexOfNextHalfEdgeInPolygon(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1X = mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge2X = mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge2);

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, new_triangle_vcoord,
       max_small_magnitude, poly_tri_result[1]);

    compute_cos_max_min_split_half_edge_polygon_triangulationC_angle
      (dimension, vcoord, new_triangle_vcoord, interior_vcoord0,
       mesh, ihalf_edge1X, max_small_magnitude, poly_tri_result[0]);

    compute_cos_max_min_split_half_edge_polygon_triangulationC_angle
      (dimension, vcoord, new_triangle_vcoord, interior_vcoord2,
       mesh, ihalf_edge2X, max_small_magnitude, poly_tri_result[2]);

    cos_min_angle = std::max(poly_tri_result[0].cos_min_triangulation_angle,
                             poly_tri_result[1].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             poly_tri_result[2].cos_min_triangulation_angle);
    flag_zero = (poly_tri_result[0].flag_zero || 
                 poly_tri_result[1].flag_zero ||
                 poly_tri_result[2].flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "  cos_min_triangle_angle: "
         << poly_tri_result[1].cos_min_triangulation_angle << endl;
    cerr << "  poly_tri_result[0].cos_min_triangulation_angle: "
         << poly_tri_result[0].cos_min_triangulation_angle << endl;
    cerr << "  poly_tri_result[0].tri_vertex_index: "
         << poly_tri_result[0].tri_vertex_index << endl;
    cerr << "  poly_tri_result[2].cos_min_triangulation_angle: "
         << poly_tri_result[2].cos_min_triangulation_angle << endl;
    cerr << "  poly_tri_result[2].tri_vertex_index: "
         << poly_tri_result[2].tri_vertex_index << endl;
    */
  }


  /// Replace triangle with triangle and edge.
  /// - Allow interior vertex in triangulation.
  /// - Version using C++ STL vector new_triangle_vcoord[], 
  ///     interior_vcoord0[] and interior_vcoord2[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle_allow_IV
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 new_triangle_vcoord[],
   const CTYPE2 interior_vcoord0[],
   const CTYPE3 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_replace_triangle_triangulation_angle_allow_IV
      (dimension, IJK::vector2pointer(vcoord), new_triangle_vcoord, 
       interior_vcoord0, interior_vcoord2,
       mesh, ihalf_edge0, max_small_magnitude, 
       poly_tri_result, cos_min_angle, flag_zero);
  }


  /// Replace triangle with triangle and edge.
  /// - Allow interior vertex in triangulation.
  /// - Version using C++ STL vector vcoord[].
  /// - Version using C++ STL vector new_triangle_vcoord[], 
  ///     interior_vcoord0[] and interior_vcoord2[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_triangle_triangulation_angle_allow_IV
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & new_triangle_vcoord,
   const std::vector<CTYPE2> & interior_vcoord0,
   const std::vector<CTYPE3> & interior_vcoord2,
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_replace_triangle_triangulation_angle_allow_IV
      (dimension, vcoord, IJK::vector2pointer(new_triangle_vcoord),
       IJK::vector2pointer(interior_vcoord0), 
       IJK::vector2pointer(interior_vcoord2), 
       mesh, ihalf_edge0, max_small_magnitude, 
       poly_tri_result, cos_min_angle, flag_zero);
  }


  /// Replace split edge triangle with triangle and edge.
  /// - Triangle has two split edges.
  /// - Allow interior vertex in triangulation.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_splitE2T_triangulation_angle
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 new_triangle_vcoord[], 
   const CTYPE2 interior_vcoord0[],
   const CTYPE3 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE 
      POLYGON_INDEX_TYPE;

    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ihalf_edge0);
    const VERTEX_INDEX_TYPE jv1 = mesh.ToVertexIndex(ihalf_edge0);
    const CTYPE0 * vcoord0 = vcoord + jv0*dimension;
    const CTYPE0 * vcoord1 = vcoord + jv1*dimension;    
    COS_TYPE cos_min_triangle_angle;
    bool flag_zeroT;

    HALF_EDGE_INDEX_TYPE ihalf_edge1, ihalf_edge2, ihalf_edge3, ihalf_edge4;

    mesh.GetIndicesOfNextFourHalfEdgesInPolygon
      (ihalf_edge0, ihalf_edge1, ihalf_edge2, ihalf_edge3, ihalf_edge4);

    const HALF_EDGE_INDEX_TYPE
      ihalf_edge1X = mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge1);
    const HALF_EDGE_INDEX_TYPE
      ihalf_edge3X = mesh.IndexOfNextHalfEdgeAroundEdge(ihalf_edge3);
    const POLYGON_INDEX_TYPE ipolyA = 
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge1);
    const POLYGON_INDEX_TYPE ipolyC = 
      mesh.IndexOfPolygonContainingNextHalfEdgeAroundEdge(ihalf_edge3);
    const NUMBER_TYPE jlocA_replace =
      mesh.LocationOfHalfEdgeInPolygon(ihalf_edge1X);
    const NUMBER_TYPE jlocC_replace =
      mesh.LocationOfHalfEdgeInPolygon(ihalf_edge3X);

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, new_triangle_vcoord,
       max_small_magnitude, poly_tri_result[1]);

    // *** CHANGE ipolyA, ipolyC to ipoly0, ipoly2 ***
    compute_cos_max_min_replace_vertex_polygon_triangulation_angleP
      (dimension, vcoord, new_triangle_vcoord, interior_vcoord0,
       mesh, ipolyA, jlocA_replace, max_small_magnitude, 
       triangulation_settings[0], poly_tri_result[0]);

    compute_cos_max_min_replace_vertex_polygon_triangulation_angleP
      (dimension, vcoord, new_triangle_vcoord, interior_vcoord2,
       mesh, ipolyC, jlocC_replace, max_small_magnitude, 
       triangulation_settings[2], poly_tri_result[2]);

    cos_min_angle = std::max(poly_tri_result[0].cos_min_triangulation_angle,
                             poly_tri_result[1].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             poly_tri_result[2].cos_min_triangulation_angle);
    flag_zero = (poly_tri_result[0].flag_zero || 
                 poly_tri_result[1].flag_zero ||
                 poly_tri_result[2].flag_zero);

    // *** DEBUG ***
    /*
    using namespace std;
    mesh.PrintPolygonIndexAndVertices(cerr, "  Poly0: ", ipolyA, "\n");
    mesh.PrintPolygonIndexAndVertices(cerr, "  Poly2: ", ipolyC, "\n");
    IJK::print_coord3D(cerr, "  New triangle vcoord: ",
                       new_triangle_vcoord, "\n");
    cerr << "  cos_min_triangle_angle: "
         << poly_tri_result[1].cos_min_triangulation_angle << endl;
    cerr << "  poly_tri_result[0].cos_min_triangulation_angle: "
         << poly_tri_result[0].cos_min_triangulation_angle << endl;
    cerr << "  poly_tri_result[0].tri_vertex_index: "
         << poly_tri_result[0].tri_vertex_index << endl;
    cerr << "  poly_tri_result[2].cos_min_triangulation_angle: "
         << poly_tri_result[2].cos_min_triangulation_angle << endl;
    cerr << "  poly_tri_result[2].tri_vertex_index: "
         << poly_tri_result[2].tri_vertex_index << endl;
    */
  }


  /// Replace split edge triangle with triangle and edge.
  /// - Triangle has two split edges.
  /// - Allow interior vertex in triangulation.
  /// - Version using C++ STL vector new_triangle_vcoord[], 
  ///     interior_vcoord0[] and interior_vcoord2[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_splitE2T_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 new_triangle_vcoord[],
   const CTYPE2 interior_vcoord0[],
   const CTYPE3 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_replace_splitE2T_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), new_triangle_vcoord, 
       interior_vcoord0, interior_vcoord2,
       mesh, ihalf_edge0, max_small_magnitude, triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero);
  }


  /// Replace split edge triangle with triangle and edge.
  /// - Triangle has two split edges.
  /// - Allow interior vertex in triangulation.
  /// - Version using C++ STL vector vcoord[].
  /// - Version using C++ STL vector new_triangle_vcoord[], 
  ///     interior_vcoord0[] and interior_vcoord2[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPE,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_splitE2T_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & new_triangle_vcoord,
   const std::vector<CTYPE2> & interior_vcoord0,
   const std::vector<CTYPE3> & interior_vcoord2,
   const MESH_TYPE & mesh, 
   const ITYPE ihalf_edge0, const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_replace_splitE2T_triangulation_angle
      (dimension, vcoord, IJK::vector2pointer(new_triangle_vcoord),
       IJK::vector2pointer(interior_vcoord0), 
       IJK::vector2pointer(interior_vcoord2), 
       mesh, ihalf_edge0, max_small_magnitude, triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero);
  }


  /// Replace triangle with one split edge by triangle and edge.
  /// @pre ilong_half_edge and ibase_half_edge are adjacent half edges
  ///      in the same polygon.
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPEA, typename ITYPEB,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_splitE1T_triangulation_angle
  (const DTYPE dimension, const CTYPE0 vcoord[], 
   const CTYPE1 new_triangle_vcoord[], 
   const CTYPE2 interior_vcoord0[],
   const CTYPE3 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPEA ibase_half_edge, const ITYPEB ilong_half_edge,
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    typedef typename MESH_TYPE::NUMBER_TYPE NUMBER_TYPE;
    typedef typename MESH_TYPE::VERTEX_INDEX_TYPE VERTEX_INDEX_TYPE;
    typedef typename MESH_TYPE::HALF_EDGE_INDEX_TYPE 
      HALF_EDGE_INDEX_TYPE;
    typedef typename MESH_TYPE::POLYGON_INDEX_TYPE 
      POLYGON_INDEX_TYPE;

    // If true, ilong_half_edge is next half edge of ibase_half_edge.
    const bool flag_is_next_half_edge =
      mesh.IsNextHalfEdge(ibase_half_edge, ilong_half_edge);
    const HALF_EDGE_INDEX_TYPE isplit_half_edge = 
      mesh.GetIndexOfAdjacentHalfEdgeInPolygon
      (ibase_half_edge, !flag_is_next_half_edge);
    const HALF_EDGE_INDEX_TYPE ilong_half_edgeX =
      mesh.IndexOfNextHalfEdgeAroundEdge(ilong_half_edge);
    const HALF_EDGE_INDEX_TYPE isplit_half_edgeX =
      mesh.IndexOfNextHalfEdgeAroundEdge(isplit_half_edge);
    const POLYGON_INDEX_TYPE ipoly0 =
      mesh.IndexOfPolygonContainingHalfEdge(ilong_half_edgeX);
    const POLYGON_INDEX_TYPE itriangle =
      mesh.IndexOfPolygonContainingHalfEdge(ibase_half_edge);
    const POLYGON_INDEX_TYPE ipoly2 =
      mesh.IndexOfPolygonContainingHalfEdge(isplit_half_edgeX);
    const VERTEX_INDEX_TYPE iv_apex =
      mesh.GetIndexOfHalfEdgeEndpoint(ilong_half_edge, flag_is_next_half_edge);
    const VERTEX_INDEX_TYPE iv_split =
      mesh.GetIndexOfHalfEdgeEndpoint(isplit_half_edge, !flag_is_next_half_edge);
    const VERTEX_INDEX_TYPE jv0 = mesh.FromVertexIndex(ibase_half_edge);
    const VERTEX_INDEX_TYPE jv1 = mesh.ToVertexIndex(ibase_half_edge);
    const NUMBER_TYPE jloc2_replace =
      mesh.GetLocationOfHalfEdgeEndpointInPolygon
      (isplit_half_edgeX, flag_is_next_half_edge);

    const CTYPE0 * vcoord0 = vcoord + jv0*dimension;
    const CTYPE0 * vcoord1 = vcoord + jv1*dimension;    

    // *** DEBUG ***
    /*
    using namespace std;
    cerr << "In " << __func__ << endl;
    */

    compute_cos_min_triangle_angle
      (dimension, vcoord0, vcoord1, new_triangle_vcoord,
       max_small_magnitude, poly_tri_result[1]);

    compute_split_half_edge_triangulation_to_max_min_angleP
      (dimension, vcoord, new_triangle_vcoord, interior_vcoord0,
       mesh, ilong_half_edgeX, max_small_magnitude,
       triangulation_settings[0], poly_tri_result[0]);

    // Ignore triangulation_settings[2].
    compute_cos_max_min_replace_vertex_polygon_triangulation_angleP
      (dimension, vcoord, new_triangle_vcoord, interior_vcoord2,
       mesh, ipoly2, jloc2_replace, max_small_magnitude, 
       triangulation_settings[2], poly_tri_result[2]);

    cos_min_angle = std::max(poly_tri_result[0].cos_min_triangulation_angle,
                             poly_tri_result[1].cos_min_triangulation_angle);
    cos_min_angle = std::max(cos_min_angle,
                             poly_tri_result[2].cos_min_triangulation_angle);
    flag_zero = (poly_tri_result[0].flag_zero || 
                 poly_tri_result[1].flag_zero ||
                 poly_tri_result[2].flag_zero);
  }


  /// Replace triangle with one split edge by triangle and edge.
  /// - Version using C++ STL vector new_triangle_vcoord[], 
  ///     interior_vcoord0[] and interior_vcoord2[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPEA, typename ITYPEB,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_splitE1T_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const CTYPE1 new_triangle_vcoord[],
   const CTYPE2 interior_vcoord0[],
   const CTYPE3 interior_vcoord2[],
   const MESH_TYPE & mesh, 
   const ITYPEA ibase_half_edge, ITYPEB ilong_half_edge, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_replace_splitE1T_triangulation_angle
      (dimension, IJK::vector2pointer(vcoord), new_triangle_vcoord, 
       interior_vcoord0, interior_vcoord2,
       mesh, ibase_half_edge, ilong_half_edge, max_small_magnitude,
       triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero);
  }


  /// Replace triangle with one split edge by triangle and edge.
  /// - Version using C++ STL vector vcoord[].
  /// - Version using C++ STL vector new_triangle_vcoord[], 
  ///     interior_vcoord0[] and interior_vcoord2[].
  template <typename DTYPE, typename CTYPE0, typename CTYPE1,
            typename CTYPE2, typename CTYPE3,
            typename MESH_TYPE, typename ITYPEA, typename ITYPEB,
            typename MTYPE, typename RESULT_TYPE, typename COS_TYPE>
  void compute_cos_max_min_replace_splitE1T_triangulation_angle
  (const DTYPE dimension, const std::vector<CTYPE0> & vcoord, 
   const std::vector<CTYPE1> & new_triangle_vcoord,
   const std::vector<CTYPE2> & interior_vcoord0,
   const std::vector<CTYPE3> & interior_vcoord2,
   const MESH_TYPE & mesh, 
   const ITYPEB ibase_half_edge, const ITYPEA ilong_half_edge, 
   const MTYPE max_small_magnitude,
   const POLYGON_TRIANGULATION_SETTINGS triangulation_settings[3],
   RESULT_TYPE poly_tri_result[3],
   COS_TYPE & cos_min_angle, bool & flag_zero)
  {
    compute_cos_max_min_replace_splitE1T_triangulation_angle
      (dimension, vcoord, IJK::vector2pointer(new_triangle_vcoord),
       IJK::vector2pointer(interior_vcoord0), 
       IJK::vector2pointer(interior_vcoord2), 
       mesh, ibase_half_edge, ilong_half_edge, max_small_magnitude,
       triangulation_settings,
       poly_tri_result, cos_min_angle, flag_zero);
  }

  ///@}

}

#endif

